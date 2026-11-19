//! Not from the Java suite, which reads alignments through htsjdk: the port's own BAM / CRAM reading.
//! Covers the tag parser, the reads a region query returns with and without the read cache, and CRAM
//! against BAM.

use super::indel_dedup::test_tmp_dir;
use super::util::*;
use crate::bamio::{for_each_aux, AuxValue, BamSlicer, ReadPool};
use crate::types::*;
use rust_htslib::bam;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Aux, Cigar, CigarString};

const CHR: &str = "1";
const REF_LENGTH: usize = 6000;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self, bound: u64) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 33) % bound
    }
}

fn reference_bases() -> Vec<u8> {
    let mut rng = Lcg(7);
    (0..REF_LENGTH).map(|_| b"ACGT"[rng.next(4) as usize]).collect()
}

/// Writes the reference as an indexed FASTA and returns its path.
fn write_reference(test_name: &str) -> String {
    let path = test_tmp_dir().join(format!("{test_name}.fasta"));
    let _ = std::fs::remove_file(path.with_extension("fasta.fai"));
    let bases = reference_bases();
    let mut fasta = format!(">{CHR}\n").into_bytes();
    for line in bases.chunks(60) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    std::fs::write(&path, fasta).unwrap();
    let path = path.to_str().unwrap().to_string();
    rust_htslib::faidx::Reader::from_path(&path).unwrap();
    path
}

struct TestRecord {
    name: String,
    /// 1-based
    start: i32,
    cigar: Vec<(u32, char)>,
    bases: Vec<u8>,
    quals: Vec<u8>,
    flags: u16,
    mapq: u8,
    mate_start: i32,
    insert_size: i32,
    nm: Option<i32>,
    mate_cigar: Option<String>,
    consensus: Option<&'static str>,
    duplex_index: Option<i32>,
}

/// Coordinate-sorted records with a mix of cigars, flags and tags. The bases follow the reference
/// apart from a few mismatches, and NM is the true edit distance where it is set.
fn test_records() -> Vec<TestRecord> {
    let reference = reference_bases();
    let mut rng = Lcg(11);
    let mut records = Vec::new();
    let mut start = 20;
    for i in 0..900 {
        start += rng.next(12) as i32;
        let cigar: Vec<(u32, char)> = match rng.next(8) {
            0 => vec![(5, 'S'), (40 + rng.next(30) as u32, 'M')],
            1 => vec![(30, 'M'), (3, 'I'), (25, 'M')],
            2 => vec![(28, 'M'), (4, 'D'), (30, 'M'), (6, 'S')],
            3 => vec![(20, 'M'), (50, 'N'), (20, 'M')],
            // covers no reference bases
            4 if i % 5 == 0 => vec![(30, 'S')],
            _ => vec![(30 + rng.next(60) as u32, 'M')],
        };
        let mut bases = Vec::new();
        let mut ref_pos = start as usize - 1;
        let mut mismatches = 0;
        for &(len, op) in &cigar {
            match op {
                'M' => {
                    for k in 0..len as usize {
                        let ref_base = reference[ref_pos + k];
                        if rng.next(40) == 0 {
                            bases.push(if ref_base == b'A' { b'C' } else { b'A' });
                            mismatches += 1;
                        } else {
                            bases.push(ref_base);
                        }
                    }
                    ref_pos += len as usize;
                }
                'I' | 'S' => {
                    for _ in 0..len {
                        bases.push(b"ACGT"[rng.next(4) as usize]);
                    }
                    if op == 'I' {
                        mismatches += len as i32;
                    }
                }
                'D' => {
                    mismatches += len as i32;
                    ref_pos += len as usize;
                }
                _ => ref_pos += len as usize,
            }
        }
        let quals: Vec<u8> = (0..bases.len()).map(|_| 2 + rng.next(38) as u8).collect();
        let mut flags = FLAG_PAIRED | FLAG_PROPER | if rng.next(2) == 0 { FLAG_REVERSE } else { FLAG_MATE_REVERSE };
        flags |= if rng.next(2) == 0 { FLAG_FIRST } else { 0x80 };
        match rng.next(25) {
            0 => flags |= FLAG_DUPLICATE,
            1 => flags |= FLAG_SECONDARY,
            2 | 3 => flags |= FLAG_SUPPLEMENTARY,
            _ => {}
        }
        records.push(TestRecord {
            name: format!("READ_{i:04}"),
            start,
            cigar,
            bases,
            quals,
            flags,
            mapq: rng.next(61) as u8,
            mate_start: start + rng.next(300) as i32,
            insert_size: rng.next(500) as i32 - 250,
            nm: if rng.next(5) == 0 { None } else { Some(mismatches) },
            mate_cigar: if rng.next(4) == 0 { None } else { Some(format!("{}M", 40 + rng.next(60))) },
            consensus: [None, None, Some("DUAL"), Some("SINGLE")][rng.next(4) as usize],
            duplex_index: if rng.next(6) == 0 { Some(rng.next(70000) as i32) } else { None },
        });
    }
    records
}

fn write_alignments(test_name: &str, format: bam::Format, reference: &str) -> String {
    let extension = if matches!(format, bam::Format::Cram) { "cram" } else { "bam" };
    let path = test_tmp_dir().join(format!("{test_name}.{extension}")).to_str().unwrap().to_string();
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", CHR);
    sq.push_tag(b"LN", REF_LENGTH);
    if matches!(format, bam::Format::Cram) {
        // the writer looks for the reference when it writes the header, before one can be set
        sq.push_tag(b"UR", reference);
    }
    header.push_record(&sq);
    {
        let mut writer = bam::Writer::from_path(&path, &header, format).unwrap();
        if matches!(format, bam::Format::Cram) {
            writer.set_reference(reference).unwrap();
        }
        for (i, r) in test_records().iter().enumerate() {
            let cigar = CigarString(
                r.cigar
                    .iter()
                    .map(|&(len, op)| match op {
                        'M' => Cigar::Match(len),
                        'I' => Cigar::Ins(len),
                        'D' => Cigar::Del(len),
                        'N' => Cigar::RefSkip(len),
                        _ => Cigar::SoftClip(len),
                    })
                    .collect(),
            );
            let mut record = bam::Record::new();
            record.set(r.name.as_bytes(), Some(&cigar), &r.bases, &r.quals);
            record.set_tid(0);
            record.set_pos(r.start as i64 - 1);
            record.set_mapq(r.mapq);
            record.set_flags(r.flags);
            record.set_mtid(0);
            record.set_mpos(r.mate_start as i64 - 1);
            record.set_insert_size(r.insert_size as i64);
            // tags SAGE does not read, of every kind, before and between the ones it does
            record.push_aux(b"XA", Aux::Char(b'q')).unwrap();
            record.push_aux(b"XF", Aux::Float(1.5)).unwrap();
            if let Some(nm) = r.nm {
                // every integer width
                let value = match i % 4 {
                    0 => Aux::U8(nm as u8),
                    1 => Aux::I16(nm as i16),
                    2 => Aux::U16(nm as u16),
                    _ => Aux::I32(nm),
                };
                record.push_aux(b"NM", value).unwrap();
            }
            record.push_aux(b"XB", Aux::ArrayI16((&[-3i16, 4, 500][..]).into())).unwrap();
            record.push_aux(b"XS", Aux::String("some text")).unwrap();
            if let Some(mate_cigar) = &r.mate_cigar {
                record.push_aux(b"MC", Aux::String(mate_cigar)).unwrap();
            }
            record.push_aux(b"XH", Aux::HexByteArray("1AE301")).unwrap();
            if let Some(consensus) = r.consensus {
                record.push_aux(b"UT", Aux::String(consensus)).unwrap();
            }
            if !matches!(format, bam::Format::Cram) {
                // CRAM cannot hold a double
                record.push_aux(b"XD", Aux::Double(2.25)).unwrap();
            }
            if let Some(duplex_index) = r.duplex_index {
                record.push_aux(b"YX", Aux::I32(duplex_index)).unwrap();
            }
            record.push_aux(b"XU", Aux::ArrayU8((&[1u8, 2, 3][..]).into())).unwrap();
            writer.write(&record).unwrap();
        }
    }
    bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
    path
}

fn ref_length(cigar: &[(u32, char)]) -> i32 {
    cigar.iter().filter(|(_, op)| matches!(op, 'M' | 'D' | 'N')).map(|(len, _)| *len as i32).sum()
}

/// The records a query of the regions should return, by the rules of an indexed query and the slicer's
/// flag filter.
fn expected_names(regions: &[BaseRegion], keep_supplementaries: bool) -> Vec<String> {
    test_records()
        .iter()
        .filter(|r| r.flags & (FLAG_DUPLICATE | FLAG_SECONDARY) == 0)
        .filter(|r| keep_supplementaries || r.flags & FLAG_SUPPLEMENTARY == 0)
        .filter(|r| {
            // a record covering no reference bases is treated as covering one
            let end = r.start + ref_length(&r.cigar).max(1) - 1;
            regions.iter().any(|region| r.start <= region.end && end >= region.start)
        })
        .map(|r| r.name.clone())
        .collect()
}

fn describe(read: &Read) -> String {
    format!("{read:?}")
}

fn region_sets() -> Vec<Vec<BaseRegion>> {
    vec![
        vec![BaseRegion::new(1, REF_LENGTH as i32)],
        vec![BaseRegion::new(1000, 2000)],
        vec![BaseRegion::new(500, 520), BaseRegion::new(560, 900), BaseRegion::new(3000, 3001)],
        vec![BaseRegion::new(100, 100), BaseRegion::new(101, 101), BaseRegion::new(140, 4000)],
        vec![BaseRegion::new(5990, 6000)],
    ]
}

#[test]
fn test_region_queries_return_the_expected_reads() {
    let reference = write_reference("alignment_reading_queries");
    let path = write_alignments("alignment_reading_queries", bam::Format::Bam, &reference);
    let mut slicer = BamSlicer::open(&path, &reference, 1).unwrap();
    let records = test_records();
    for regions in region_sets() {
        for keep_supplementaries in [true, false] {
            let mut reads: Vec<Read> = Vec::new();
            slicer.slice_ref(CHR, &regions, keep_supplementaries, |read| reads.push(read.clone())).unwrap();
            let names: Vec<String> = reads.iter().map(|r| r.name.clone()).collect();
            assert_eq!(names, expected_names(&regions, keep_supplementaries));
            for read in &reads {
                let expected = records.iter().find(|r| r.name == read.name).unwrap();
                assert_eq!(read.start, expected.start);
                assert_eq!(read.end, expected.start + ref_length(&expected.cigar) - 1);
                assert_eq!(read.bases, expected.bases);
                assert_eq!(read.quals, expected.quals);
                assert_eq!(read.flags, expected.flags);
                assert_eq!(read.mapq, expected.mapq as i32);
                assert_eq!(read.mate_start, expected.mate_start);
                assert_eq!(read.insert_size, expected.insert_size);
                assert_eq!(read.nm, expected.nm);
                assert_eq!(read.mate_cigar, expected.mate_cigar);
                assert_eq!(read.has_consensus_tag, expected.consensus.is_some());
                assert_eq!(read.duplex_index, expected.duplex_index.unwrap_or(-1));
                let cigar: Vec<(u32, char)> = read
                    .cigar
                    .iter()
                    .map(|el| {
                        (el.len as u32, match el.op {
                            Op::M => 'M',
                            Op::I => 'I',
                            Op::D => 'D',
                            Op::N => 'N',
                            _ => 'S',
                        })
                    })
                    .collect();
                assert_eq!(cigar, expected.cigar);
            }
        }
    }
}

/// The candidate stage's caching read gives the consumer what a plain read gives it, and the evidence
/// stage's read gives the same reads from the cache as from the file.
#[test]
fn test_cached_reads_equal_reads_from_the_file() {
    let reference = write_reference("alignment_reading_cache");
    let path = write_alignments("alignment_reading_cache", bam::Format::Bam, &reference);
    let mut caching = BamSlicer::open(&path, &reference, 1).unwrap();
    let mut plain = BamSlicer::open(&path, &reference, 1).unwrap();
    let pool = ReadPool::default();
    for regions in region_sets() {
        for margin in [0, 30] {
            let mut cached: Vec<String> = Vec::new();
            caching.slice_ref_caching(CHR, &regions, margin, |read| cached.push(describe(read))).unwrap();
            let mut direct: Vec<String> = Vec::new();
            plain.slice_ref(CHR, &regions, true, |read| direct.push(describe(read))).unwrap();
            assert_eq!(cached, direct);

            // evidence-stage queries: inside the cached span, reaching into the margin, and outside
            // it (which reads the file)
            let first = regions[0];
            let last = regions[regions.len() - 1];
            let queries = [
                regions.clone(),
                vec![BaseRegion::new((first.start - margin).max(1), first.end + margin)],
                vec![BaseRegion::new(first.start, first.start), BaseRegion::new(last.end, last.end + margin)],
                vec![BaseRegion::new((first.start - 200).max(1), last.end + 200)],
                vec![BaseRegion::new(first.start + (first.end - first.start) / 2, first.end)],
            ];
            for query in queries {
                for keep_supplementaries in [false, true] {
                    let mut from_cache: Vec<String> = Vec::new();
                    caching
                        .slice_pooled(CHR, &query, keep_supplementaries, &pool, |read| {
                            from_cache.push(describe(&read));
                            pool.put(read);
                        })
                        .unwrap();
                    let mut from_file: Vec<String> = Vec::new();
                    plain
                        .slice_pooled(CHR, &query, keep_supplementaries, &pool, |read| {
                            from_file.push(describe(&read));
                            pool.put(read);
                        })
                        .unwrap();
                    assert_eq!(from_cache, from_file);
                }
            }
        }
    }
}

/// A CRAM gives the reads its BAM gives, except that NM is absent where the CRAM left it out: it is
/// not regenerated, as htsjdk does not regenerate it.
#[test]
fn test_cram_reads_equal_bam_reads() {
    let reference = write_reference("alignment_reading_cram");
    let bam_path = write_alignments("alignment_reading_cram", bam::Format::Bam, &reference);
    let cram_path = write_alignments("alignment_reading_cram", bam::Format::Cram, &reference);
    let mut from_bam = BamSlicer::open(&bam_path, &reference, 1).unwrap();
    let mut from_cram = BamSlicer::open(&cram_path, &reference, 1).unwrap();
    let mut without_nm = 0;
    for regions in region_sets() {
        let mut bam_reads: Vec<Read> = Vec::new();
        from_bam.slice_ref(CHR, &regions, true, |read| bam_reads.push(read.clone())).unwrap();
        let mut cram_reads: Vec<Read> = Vec::new();
        from_cram.slice_ref(CHR, &regions, true, |read| cram_reads.push(read.clone())).unwrap();
        assert_eq!(bam_reads.len(), cram_reads.len());
        for (b, c) in bam_reads.iter_mut().zip(cram_reads.iter_mut()) {
            if c.nm.is_none() && b.nm.is_some() {
                without_nm += 1;
                b.nm = None;
            }
            assert_eq!(describe(b), describe(c));
        }
    }
    assert!(without_nm > 0);
}

#[test]
fn test_tag_parser() {
    let mut aux: Vec<u8> = Vec::new();
    aux.extend_from_slice(b"XAAq");
    aux.extend_from_slice(b"NMc\xfe");
    aux.extend_from_slice(b"N2C\xfe");
    aux.extend_from_slice(b"N3s\xfe\xff");
    aux.extend_from_slice(b"N4S\xfe\xff");
    aux.extend_from_slice(b"N5i\xfe\xff\xff\xff");
    aux.extend_from_slice(b"N6I\x10\x00\x01\x00");
    aux.extend_from_slice(b"XFf\x00\x00\xc0\x3f");
    aux.extend_from_slice(b"XDd\x00\x00\x00\x00\x00\x00\x02\x40");
    aux.extend_from_slice(b"MCZ10M2D5M\0");
    aux.extend_from_slice(b"XHH1AE3\0");
    aux.extend_from_slice(b"tpBc\x03\x00\x00\x00\x01\xff\x02");
    aux.extend_from_slice(b"XBBs\x02\x00\x00\x00\x01\x00\x02\x00");
    aux.extend_from_slice(b"XIBI\x01\x00\x00\x00\x01\x00\x02\x00");
    aux.extend_from_slice(b"XEZ\0");
    let parse = |aux: &[u8]| -> Vec<String> {
        let mut seen = Vec::new();
        for_each_aux(aux, |tag, value| {
            let tag = std::str::from_utf8(tag).unwrap();
            seen.push(match value {
                AuxValue::Int(v) => format!("{tag}={v}"),
                AuxValue::Str(s) => format!("{tag}='{}'", std::str::from_utf8(s).unwrap()),
                AuxValue::ByteArray(a) => format!("{tag}={a:?}"),
                AuxValue::Other => format!("{tag}=?"),
            });
        });
        seen
    };
    let expected = [
        "XA=?", "NM=-2", "N2=254", "N3=-2", "N4=65534", "N5=-2", "N6=65552", "XF=?", "XD=?", "MC='10M2D5M'", "XH=?",
        "tp=[1, 255, 2]", "XB=?", "XI=?", "XE=''",
    ];
    assert_eq!(parse(&aux), expected);
    // a truncated or malformed section gives the tags before the damage and stops
    for cut in 0..aux.len() {
        let seen = parse(&aux[..cut]);
        assert!(seen.len() <= expected.len());
        assert_eq!(seen[..], expected[..seen.len()]);
    }
    let mut unknown_type = b"NMc\x01".to_vec();
    unknown_type.extend_from_slice(b"XX?\x01NMc\x02");
    assert_eq!(parse(&unknown_type), ["NM=1"]);
}
