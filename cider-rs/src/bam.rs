//! Reading the BAM (CiderApplication.readBamFile, AsyncBamReader) and writing cider.bam
//! (CiderApplication.writeCiderBam) as htsjdk's BAM writer does.
//!
//! Java queries each region in turn and hands the records to its consumer threads; with one thread the
//! records are screened in BAM order, region by region. The port reads in that order and screens
//! batches of records on its threads, then takes the results in read order, so its candidate lists are
//! those of a one-thread Java run at any thread count (CI9).

use std::collections::HashSet;
use std::io::Write;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use rust_htslib::bam::{self, Read as _};

use crate::compat;
use crate::genes::GeneDatastore;
use crate::read::{CigarOp, CiderRead, ReadRef};
use crate::screener::{adapter_sequence_trim, CiderReadScreener, Screened, VjReadCandidate};

/// a region to query: chromosome, 1-based inclusive start and end
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

/// GenomeRegion.compareTo: the same chromosome by start only, else by HumanChromosome rank
/// (ContigComparator), equal ranks comparing equal
fn java_region_cmp(a: &QueryRegion, b: &QueryRegion) -> std::cmp::Ordering {
    use hmf_common::regions::chromosome_rank_java;
    if a.chromosome == b.chromosome {
        a.start.cmp(&b.start)
    } else {
        chromosome_rank_java(&a.chromosome).cmp(&chromosome_rank_java(&b.chromosome))
    }
}

/// The TreeSet<GenomeRegion> of readBamFile: anchor locations and constant/D regions, each widened by
/// the fragment length, in TreeSet order.
///
/// CI8: GenomeRegion.compareTo looks at the chromosome (its rank) and the start only, so the TreeSet keeps
/// the first of two regions with the same start and drops the other even when it reaches further.
/// Corrected: such a region is kept (ordered after the shorter one); only identical regions are merged.
pub fn query_regions(ds: &GeneDatastore, max_fragment_length: i32) -> Result<Vec<QueryRegion>> {
    let mut candidates = Vec::new();
    for l in &ds.anchor_locations {
        if !l.location.in_primary_assembly() {
            bail!("java.lang.IllegalArgumentException: Failed requirement (anchor location in primary assembly)");
        }
        candidates.push(QueryRegion { chromosome: l.location.chromosome.clone(), start: l.start() - max_fragment_length, end: l.end() + max_fragment_length });
    }
    for r in &ds.constant_diversity_regions {
        if !r.location.in_primary_assembly() {
            bail!("java.lang.IllegalArgumentException: Failed requirement (constant region in primary assembly)");
        }
        candidates.push(QueryRegion {
            chromosome: r.location.chromosome.clone(),
            start: r.location.pos_start - max_fragment_length,
            end: r.location.pos_end + max_fragment_length,
        });
    }
    let java = compat::java("CI8");
    let mut set: Vec<QueryRegion> = Vec::new();
    let mut dropped = 0;
    for c in candidates {
        let cmp = |x: &QueryRegion| {
            if java {
                java_region_cmp(x, &c)
            } else {
                // corrected: rank, name, start, end (the Java order for the human chromosomes)
                use hmf_common::regions::chromosome_rank_java;
                chromosome_rank_java(&x.chromosome)
                    .cmp(&chromosome_rank_java(&c.chromosome))
                    .then_with(|| x.chromosome.cmp(&c.chromosome))
                    .then(x.start.cmp(&c.start))
                    .then(x.end.cmp(&c.end))
            }
        };
        match set.binary_search_by(cmp) {
            Ok(i) => {
                if set[i] != c {
                    dropped += 1;
                }
            }
            Err(i) => set.insert(i, c),
        }
    }
    if dropped > 0 {
        crate::debug!("{dropped} regions with the start of another region dropped (Java TreeSet)");
    }
    Ok(set)
}

/// The reads collected by the screener
pub struct ScreenResult {
    pub vj_read_candidates: Vec<VjReadCandidate>,
    /// the matched reads in screening order, with their htsjdk BAM encoding
    pub all_matched_reads: Vec<(ReadRef, Vec<u8>)>,
    pub vj_mate_candidates: Vec<VjReadCandidate>,
    /// the BAM header, for cider.bam
    pub header_text: Vec<u8>,
    pub targets: Vec<(String, u32)>,
}

/// the default read size of the BAM file, in KB
const FILE_BUFFER_KB: i32 = 1024;

/// ReadRecordKey of the Java screener
#[derive(Hash, PartialEq, Eq)]
struct ReadRecordKey {
    name: Box<[u8]>,
    first_of_pair: bool,
    tid: i32,
    alignment_start: i32,
}

fn mate_cigar(record: &bam::Record) -> Option<String> {
    match record.aux(b"MC") {
        Ok(bam::record::Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

fn to_cider_read(record: &bam::Record, contigs: &[Arc<str>], star: &Arc<str>, order: usize) -> CiderRead {
    let core = &record.inner().core;
    let contig = |tid: i32| if tid >= 0 && (tid as usize) < contigs.len() { contigs[tid as usize].clone() } else { star.clone() };
    let flags = core.flag;
    let cigar = record.raw_cigar().iter().map(|&c| (CigarOp::from_bam(c & 0xf), c >> 4)).collect();
    let quals = record.qual();
    // htsjdk: no qualities (0xff) is an empty array
    let quals = if quals.first() == Some(&0xff) { Vec::new() } else { quals.to_vec() };
    let mut r = CiderRead {
        name: match std::str::from_utf8(record.qname()) {
            Ok(n) => n.to_string(),
            Err(_) => String::from_utf8_lossy(record.qname()).into_owned(),
        },
        flags,
        tid: core.tid,
        contig: contig(core.tid),
        alignment_start: (core.pos + 1) as i32,
        mapq: core.qual,
        mate_tid: core.mtid,
        mate_contig: contig(core.mtid),
        mate_alignment_start: (core.mpos + 1) as i32,
        insert_size: core.isize_ as i32,
        cigar,
        bases: record.seq().as_bytes(),
        quals,
        mate_cigar: mate_cigar(record),
        adapter_trim: (0, 0),
        order,
    };
    r.adapter_trim = adapter_sequence_trim(&r);
    r
}

/// readBamFile: queries the regions in order, drops duplicates and records already seen (the same name,
/// pair end, reference and start), screens the rest, then the candidate mates.
pub fn read_bam(
    bam_path: &str,
    ref_genome: Option<&str>,
    regions: &[QueryRegion],
    screener: &CiderReadScreener,
    threads: usize,
) -> Result<ScreenResult> {
    let mut reader = bam::IndexedReader::from_path(bam_path).map_err(|e| anyhow!("opening {bam_path}: {e}"))?;
    if let Some(r) = ref_genome {
        reader.set_reference(r).map_err(|e| anyhow!("setting reference {r}: {e}"))?;
    }
    if bam_path.ends_with(".cram") {
        // htsjdk does not make up the MD and NM tags of a CRAM record; htslib does unless told not to
        use rust_htslib::htslib;
        unsafe { htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 0 as std::os::raw::c_int) };
    }
    // the read size of the BAM file: the regions lie close together within each locus, so larger reads
    // serve several regions (CIDER_FILE_BUFFER_KB, 0 for the htslib default)
    let buffer_kb: i32 = std::env::var("CIDER_FILE_BUFFER_KB").ok().and_then(|v| v.parse().ok()).unwrap_or(FILE_BUFFER_KB);
    if buffer_kb > 0 {
        use rust_htslib::htslib;
        unsafe { htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_HTS_OPT_BLOCK_SIZE, (buffer_kb << 10) as std::os::raw::c_int) };
    }
    let header = reader.header().clone();
    let header_text = header.as_bytes().to_vec();
    let mut contigs: Vec<Arc<str>> = Vec::new();
    let mut targets = Vec::new();
    for tid in 0..header.target_count() {
        let name = String::from_utf8_lossy(header.tid2name(tid)).into_owned();
        targets.push((name.clone(), header.target_len(tid).unwrap_or(0) as u32));
        contigs.push(Arc::from(name.as_str()));
    }
    let star: Arc<str> = Arc::from("*");
    crate::info!("{} bam record consumer threads started", threads);

    let mut result = ScreenResult { vj_read_candidates: Vec::new(), all_matched_reads: Vec::new(), vj_mate_candidates: Vec::new(), header_text, targets };
    let mut unmatched_mapped: Vec<ReadRef> = Vec::new();
    let mut screen_time = std::time::Duration::ZERO;
    let from_bam = !bam_path.ends_with(".cram");
    let mut take = |batch: Batch, result: &mut ScreenResult, unmatched: &mut Vec<ReadRef>, screen_threads: usize| -> Result<()> {
        let t = std::time::Instant::now();
        let screened = screen_batch(screener, &batch.reads, screen_threads)?;
        for ((s, read), record) in screened.into_iter().zip(batch.reads).zip(batch.records) {
            match s {
                Screened::Candidate(c) => {
                    result.all_matched_reads.push((read, encode_htsjdk_record(&record, from_bam)?));
                    result.vj_read_candidates.push(c);
                }
                Screened::UnmatchedMapped => unmatched.push(read),
                Screened::Unmatched => {}
            }
        }
        screen_time += t.elapsed();
        Ok(())
    };
    let ctx = ReadContext { header: &header, contigs: &contigs, star: &star, regions, bam_path };
    let counts = if threads <= 1 {
        read_regions(&mut reader, &ctx, &mut |b| take(b, &mut result, &mut unmatched_mapped, 1))?
    } else {
        // a reader thread hands the batches, in order, to the screening thread(s)
        std::thread::scope(|scope| -> Result<ReadCounts> {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Batch>(4);
            let ctx = &ctx;
            let handle = scope.spawn(move || read_regions(&mut reader, ctx, &mut |b| tx.send(b).map_err(|_| anyhow!("screening stopped"))));
            let mut failed: Option<anyhow::Error> = None;
            for batch in rx {
                if failed.is_none() {
                    if let Err(e) = take(batch, &mut result, &mut unmatched_mapped, threads - 1) {
                        failed = Some(e);
                    }
                }
            }
            let counts = handle.join().map_err(|_| anyhow!("BAM reader thread panicked"))?;
            if let Some(e) = failed {
                return Err(e);
            }
            counts
        })?
    };
    crate::info!("{} bam reader threads finished", threads);
    crate::debug!(
        "{} records in {} regions: {} duplicates, {} seen before, {} screened ({:.2}s screening), {} unmatched mapped",
        counts.records,
        regions.len(),
        counts.duplicates,
        counts.seen,
        counts.screened,
        screen_time.as_secs_f64(),
        unmatched_mapped.len()
    );
    result.vj_mate_candidates = screener.process_candidate_mates(&result.vj_read_candidates, &unmatched_mapped)?;
    crate::info!("found {} VJ read records with {} mates", result.all_matched_reads.len(), result.vj_mate_candidates.len());
    Ok(result)
}

/// records with their reads, in read order
struct Batch {
    records: Vec<bam::Record>,
    reads: Vec<ReadRef>,
}

struct ReadContext<'a> {
    header: &'a bam::HeaderView,
    contigs: &'a [Arc<str>],
    star: &'a Arc<str>,
    regions: &'a [QueryRegion],
    bam_path: &'a str,
}

#[derive(Default)]
struct ReadCounts {
    records: usize,
    duplicates: usize,
    seen: usize,
    screened: usize,
}

const BATCH: usize = 4096;

/// queries the regions in order; the records not duplicates and not seen before, in batches
fn read_regions(reader: &mut bam::IndexedReader, ctx: &ReadContext, out: &mut dyn FnMut(Batch) -> Result<()>) -> Result<ReadCounts> {
    let mut seen: HashSet<ReadRecordKey> = HashSet::new();
    let mut counts = ReadCounts::default();
    let mut batch = Batch { records: Vec::with_capacity(BATCH), reads: Vec::with_capacity(BATCH) };
    for region in ctx.regions {
        let Some(tid) = ctx.header.tid(region.chromosome.as_bytes()) else { continue };
        let start0 = (region.start as i64 - 1).max(0);
        reader.fetch((tid as i32, start0, region.end as i64)).map_err(|e| anyhow!("querying {}:{}-{}: {e}", region.chromosome, region.start, region.end))?;
        let mut record = bam::Record::new();
        while let Some(r) = reader.read(&mut record) {
            r.map_err(|e| anyhow!("reading {}: {e}", ctx.bam_path))?;
            counts.records += 1;
            if record.is_duplicate() {
                counts.duplicates += 1;
                continue;
            }
            let core = &record.inner().core;
            let paired = core.flag & 0x1 != 0;
            let key = ReadRecordKey { name: record.qname().into(), first_of_pair: !paired || core.flag & 0x40 != 0, tid: core.tid, alignment_start: (core.pos + 1) as i32 };
            if !seen.insert(key) {
                counts.seen += 1;
                continue;
            }
            batch.reads.push(Arc::new(to_cider_read(&record, ctx.contigs, ctx.star, counts.screened)));
            counts.screened += 1;
            batch.records.push(std::mem::replace(&mut record, bam::Record::new()));
            if batch.reads.len() == BATCH {
                out(std::mem::replace(&mut batch, Batch { records: Vec::with_capacity(BATCH), reads: Vec::with_capacity(BATCH) }))?;
            }
        }
    }
    if !batch.reads.is_empty() {
        out(batch)?;
    }
    Ok(counts)
}

/// screens the reads on `threads` threads; the results in read order
fn screen_batch(screener: &CiderReadScreener, reads: &[ReadRef], threads: usize) -> Result<Vec<Screened>> {
    if threads <= 1 || reads.len() < 64 {
        return reads.iter().map(|r| screener.screen(r)).collect();
    }
    // interleaved, so that the expensive unmapped reads spread over the threads
    let n = threads.min(reads.len());
    let parts: Vec<Result<Vec<(usize, Screened)>>> = std::thread::scope(|s| {
        let handles: Vec<_> = (1..n)
            .map(|t| s.spawn(move || reads.iter().enumerate().skip(t).step_by(n).map(|(i, r)| screener.screen(r).map(|x| (i, x))).collect()))
            .collect();
        let mut parts = vec![reads.iter().enumerate().step_by(n).map(|(i, r)| screener.screen(r).map(|x| (i, x))).collect()];
        for h in handles {
            parts.push(h.join().unwrap_or_else(|_| Err(anyhow!("screening thread panicked"))));
        }
        parts
    });
    let mut out: Vec<Option<Screened>> = (0..reads.len()).map(|_| None).collect();
    for p in parts {
        for (i, x) in p? {
            out[i] = Some(x);
        }
    }
    Ok(out.into_iter().map(|x| x.expect("screened")).collect())
}

// ---- cider.bam as htsjdk writes it ----

/// htsjdk BinaryTagCodec.getIntegerType: the smallest type, signed first
fn integer_type(v: i64) -> Result<u8> {
    Ok(if v > u32::MAX as i64 {
        bail!("Integer attribute value too large to be encoded in BAM")
    } else if v > i32::MAX as i64 {
        b'I'
    } else if v > 65535 {
        b'i'
    } else if v > 32767 {
        b'S'
    } else if v > 255 {
        b's'
    } else if v > 127 {
        b'C'
    } else if v >= -128 {
        b'c'
    } else if v >= -32768 {
        b's'
    } else if v >= i32::MIN as i64 {
        b'i'
    } else {
        bail!("Integer attribute value too negative to be encoded in BAM")
    })
}

fn write_integer(out: &mut Vec<u8>, v: i64) -> Result<()> {
    let t = integer_type(v)?;
    out.push(t);
    match t {
        b'c' | b'C' => out.push(v as u8),
        b's' | b'S' => out.extend_from_slice(&(v as u16).to_le_bytes()),
        _ => out.extend_from_slice(&(v as u32).to_le_bytes()),
    }
    Ok(())
}

/// BinaryTagCodec.writeTag of decoded tags: integers in the smallest type, hex strings as byte arrays,
/// the rest unchanged (htsjdk writes a record it did not read from a BAM, as from a CRAM, this way)
fn reencode_aux(aux: &[u8], out: &mut Vec<u8>) -> Result<()> {
    let bad = || anyhow!("malformed BAM tags");
    let mut i = 0;
    while i < aux.len() {
        if i + 3 > aux.len() {
            bail!(bad());
        }
        let t = aux[i + 2];
        out.extend_from_slice(&aux[i..i + 2]);
        i += 3;
        let get = |n: usize, i: usize| -> Result<&[u8]> { aux.get(i..i + n).ok_or_else(bad) };
        match t {
            b'A' => {
                out.push(b'A');
                out.push(get(1, i)?[0]);
                i += 1;
            }
            b'c' | b'C' => {
                let b = get(1, i)?[0];
                write_integer(out, if t == b'c' { b as i8 as i64 } else { b as i64 })?;
                i += 1;
            }
            b's' | b'S' => {
                let v = <[u8; 2]>::try_from(get(2, i)?).unwrap();
                write_integer(out, if t == b's' { i16::from_le_bytes(v) as i64 } else { u16::from_le_bytes(v) as i64 })?;
                i += 2;
            }
            b'i' | b'I' => {
                let v = <[u8; 4]>::try_from(get(4, i)?).unwrap();
                write_integer(out, if t == b'i' { i32::from_le_bytes(v) as i64 } else { u32::from_le_bytes(v) as i64 })?;
                i += 4;
            }
            b'f' => {
                out.push(b'f');
                out.extend_from_slice(get(4, i)?);
                i += 4;
            }
            b'Z' => {
                let end = aux[i..].iter().position(|&b| b == 0).ok_or_else(bad)?;
                out.push(b'Z');
                out.extend_from_slice(&aux[i..=i + end]);
                i += end + 1;
            }
            b'H' => {
                let end = aux[i..].iter().position(|&b| b == 0).ok_or_else(bad)?;
                let bytes: Vec<u8> = aux[i..i + end]
                    .chunks(2)
                    .map(|c| std::str::from_utf8(c).ok().and_then(|x| u8::from_str_radix(x, 16).ok()).ok_or_else(bad))
                    .collect::<Result<_>>()?;
                out.extend_from_slice(b"Bc");
                out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
                out.extend_from_slice(&bytes);
                i += end + 1;
            }
            b'B' => {
                let sub = get(1, i)?[0];
                let n = u32::from_le_bytes(<[u8; 4]>::try_from(get(4, i + 1)?).unwrap()) as usize;
                let size = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => bail!(bad()),
                };
                out.push(b'B');
                out.extend_from_slice(get(5 + n * size, i)?);
                i += 5 + n * size;
            }
            _ => bail!("Unrecognized tag type: {}", t as char),
        }
    }
    Ok(())
}

/// GenomicIndexUtil.regionToBin (reg2bin), 0-based start, end exclusive
fn region_to_bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let b = if beg >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (beg >> 26)
    } else {
        0
    };
    b as u16
}

/// BAMRecordCodec.encode: a record read from a BAM and not changed keeps its variable part (read name,
/// CIGAR, bases, qualities, tags) and bin as they were in the file; one from a CRAM (`from_bam` false) has
/// its tags encoded again and its bin computed (SAMRecord.computeIndexingBin)
pub fn encode_htsjdk_record(record: &bam::Record, from_bam: bool) -> Result<Vec<u8>> {
    let core = &record.inner().core;
    // the variable-length data of the record (read name, CIGAR, bases, qualities, tags)
    let data: &[u8] = unsafe { std::slice::from_raw_parts(record.inner().data as *const u8, record.inner().l_data as usize) };
    let l_qname = core.l_qname as usize;
    let name_len = l_qname - core.l_extranul as usize;
    let n_cigar = core.n_cigar as usize;
    let l_seq = core.l_qseq.max(0) as usize;
    let cigar_start = l_qname;
    let mut body: Vec<u8> = Vec::with_capacity(32 + data.len());
    body.extend_from_slice(&core.tid.to_le_bytes());
    body.extend_from_slice(&(core.pos as i32).to_le_bytes());
    body.push(name_len as u8);
    body.push(core.qual);
    // the bin read from the file (computed for a CRAM record); 0 for no reference
    let bin = if core.tid < 0 {
        0
    } else if from_bam {
        core.bin
    } else {
        let start = core.pos;
        let end = if record.is_unmapped() { 0 } else { start + crate::read::reference_length(&record.raw_cigar().iter().map(|&c| (CigarOp::from_bam(c & 0xf), c >> 4)).collect::<Vec<_>>()) as i64 };
        let end = if end <= start { start + 1 } else { end };
        region_to_bin(start, end)
    };
    body.extend_from_slice(&bin.to_le_bytes());
    body.extend_from_slice(&(n_cigar as u16).to_le_bytes());
    body.extend_from_slice(&core.flag.to_le_bytes());
    body.extend_from_slice(&(l_seq as i32).to_le_bytes());
    body.extend_from_slice(&core.mtid.to_le_bytes());
    body.extend_from_slice(&(core.mpos as i32).to_le_bytes());
    body.extend_from_slice(&(core.isize_ as i32).to_le_bytes());
    body.extend_from_slice(&data[..name_len]);
    if from_bam {
        body.extend_from_slice(&data[cigar_start..]);
    } else {
        let aux_start = cigar_start + n_cigar * 4 + l_seq.div_ceil(2) + l_seq;
        body.extend_from_slice(&data[cigar_start..aux_start]);
        reencode_aux(&data[aux_start..], &mut body)?;
    }
    let mut out = Vec::with_capacity(body.len() + 4);
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(&body);
    Ok(out)
}

/// SAMTextHeaderCodec.encode of the parsed header: @HD, then @SQ (SN, LN, the rest), @RG (ID first),
/// @PG (ID first) and @CO lines, each kind in its order
pub fn htsjdk_header_text(text: &[u8]) -> String {
    let text = String::from_utf8_lossy(text);
    let mut hd: Option<String> = None;
    let (mut sq, mut rg, mut pg, mut co) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let first = |fields: &[&str], key: &str| -> String {
        let mut out = vec![fields[0].to_string()];
        let mut keyed = Vec::new();
        let mut rest = Vec::new();
        for f in &fields[1..] {
            if key.split(',').any(|k| f.starts_with(&format!("{k}:"))) {
                keyed.push(f.to_string());
            } else {
                rest.push(f.to_string());
            }
        }
        // the keys in the given order
        for k in key.split(',') {
            if let Some(f) = keyed.iter().find(|f| f.starts_with(&format!("{k}:"))) {
                out.push(f.clone());
            }
        }
        out.extend(rest);
        out.join("\t")
    };
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "@HD" => hd = Some(line.to_string()),
            "@SQ" => sq.push(first(&fields, "SN,LN")),
            "@RG" => rg.push(first(&fields, "ID")),
            "@PG" => pg.push(first(&fields, "ID")),
            "@CO" => co.push(line.to_string()),
            _ => {}
        }
    }
    let mut out = String::new();
    // htsjdk's decoded header has its own version (VN:1.6), not the input's
    out.push_str(&hd.map(|l| hmf_common::sam_header::htsjdk_hd_line(&l)).unwrap_or_else(|| "@HD\tVN:1.6".to_string()));
    out.push('\n');
    for l in sq.iter().chain(&rg).chain(&pg).chain(&co) {
        out.push_str(l);
        out.push('\n');
    }
    out
}

/// the record fields SAMRecordCoordinateComparator compares
struct SortKey {
    tid: i32,
    pos: i32,
    reverse: bool,
    name: Vec<u8>,
    flags: u16,
    mapq: u8,
    mtid: i32,
    mpos: i32,
    isize: i32,
}

fn sort_key(read: &CiderRead) -> SortKey {
    SortKey {
        tid: read.tid,
        pos: read.alignment_start,
        reverse: read.negative_strand(),
        name: read.name.as_bytes().to_vec(),
        flags: read.flags,
        mapq: read.mapq,
        mtid: read.mate_tid,
        mpos: read.mate_alignment_start,
        isize: read.insert_size,
    }
}

/// SAMRecordCoordinateComparator: reference (unmapped last), start, forward strand first, then name, flags,
/// mapping quality, mate reference, mate start and insert size
fn coordinate_cmp(a: &SortKey, b: &SortKey) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let file_order = match (a.tid, b.tid) {
        (-1, -1) => Ordering::Equal,
        (-1, _) => Ordering::Greater,
        (_, -1) => Ordering::Less,
        _ => a.tid.cmp(&b.tid).then((a.pos as i64 - b.pos as i64).cmp(&0)),
    };
    if file_order != Ordering::Equal {
        return file_order;
    }
    if a.reverse != b.reverse {
        return if a.reverse { Ordering::Greater } else { Ordering::Less };
    }
    a.name
        .cmp(&b.name)
        .then(a.flags.cmp(&b.flags))
        .then(a.mapq.cmp(&b.mapq))
        .then(a.mtid.cmp(&b.mtid))
        .then(a.mpos.cmp(&b.mpos))
        .then(a.isize.cmp(&b.isize))
}

/// htsjdk BlockCompressedOutputStream: 65498-byte blocks deflated at level 5 (stored when that does not
/// fit), then the empty EOF block
struct BgzfWriter<W: Write> {
    out: W,
    buffer: Vec<u8>,
}

const BGZF_BLOCK_SIZE: usize = 65498;
const BGZF_EOF: [u8; 28] = [0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0];

impl<W: Write> BgzfWriter<W> {
    fn new(out: W) -> Self {
        BgzfWriter { out, buffer: Vec::with_capacity(BGZF_BLOCK_SIZE) }
    }
    fn write_all(&mut self, mut data: &[u8]) -> Result<()> {
        while !data.is_empty() {
            let n = (BGZF_BLOCK_SIZE - self.buffer.len()).min(data.len());
            self.buffer.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.buffer.len() == BGZF_BLOCK_SIZE {
                self.deflate_block()?;
            }
        }
        Ok(())
    }
    fn deflate_block(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let mut enc = flate2::write::DeflateEncoder::new(Vec::with_capacity(self.buffer.len()), flate2::Compression::new(5));
        enc.write_all(&self.buffer)?;
        let mut data = enc.finish()?;
        if data.len() + 26 > 65536 {
            let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::none());
            enc.write_all(&self.buffer)?;
            data = enc.finish()?;
        }
        let mut crc = flate2::Crc::new();
        crc.update(&self.buffer);
        let bsize = (data.len() + 25) as u16;
        let mut block = Vec::with_capacity(data.len() + 26);
        block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00]);
        block.extend_from_slice(&bsize.to_le_bytes());
        block.extend_from_slice(&data);
        block.extend_from_slice(&crc.sum().to_le_bytes());
        block.extend_from_slice(&(self.buffer.len() as u32).to_le_bytes());
        self.out.write_all(&block)?;
        self.buffer.clear();
        Ok(())
    }
    fn finish(mut self) -> Result<W> {
        self.deflate_block()?;
        self.out.write_all(&BGZF_EOF)?;
        self.out.flush()?;
        Ok(self.out)
    }
}

/// writeCiderBam: the matched reads under the input's header; sorted by coordinate when the header says
/// so (SAMFileWriterImpl with presorted false), else in screening order
pub fn write_cider_bam(path: &str, header_text: &[u8], targets: &[(String, u32)], reads: &[(ReadRef, Vec<u8>)]) -> Result<()> {
    let text = htsjdk_header_text(header_text);
    let coordinate = text.lines().next().is_some_and(|hd| hd.split('\t').any(|f| f == "SO:coordinate"));
    let mut order: Vec<usize> = (0..reads.len()).collect();
    if coordinate {
        let keys: Vec<SortKey> = reads.iter().map(|(r, _)| sort_key(r)).collect();
        // stable, as Java's sort of the records
        order.sort_by(|&a, &b| coordinate_cmp(&keys[a], &keys[b]));
    }
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let mut w = BgzfWriter::new(std::io::BufWriter::with_capacity(1 << 20, file));
    let mut h: Vec<u8> = Vec::with_capacity(text.len() + 64 * targets.len());
    h.extend_from_slice(b"BAM\x01");
    h.extend_from_slice(&(text.len() as u32).to_le_bytes());
    h.extend_from_slice(text.as_bytes());
    h.extend_from_slice(&(targets.len() as u32).to_le_bytes());
    for (name, len) in targets {
        h.extend_from_slice(&((name.len() + 1) as u32).to_le_bytes());
        h.extend_from_slice(name.as_bytes());
        h.push(0);
        h.extend_from_slice(&len.to_le_bytes());
    }
    w.write_all(&h)?;
    for i in order {
        w.write_all(&reads[i].1)?;
    }
    w.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_types() {
        let t = |v| integer_type(v).unwrap();
        assert_eq!((t(0), t(127), t(128), t(255), t(256), t(40000), t(70000)), (b'c', b'c', b'C', b'C', b's', b'S', b'i'));
        assert_eq!((t(-1), t(-200), t(-40000), t(3_000_000_000)), (b'c', b's', b'i', b'I'));
    }

    #[test]
    fn aux_reencoding() {
        let mut aux = b"NMC\x01ASS".to_vec();
        aux.extend_from_slice(&300u16.to_le_bytes());
        aux.extend_from_slice(b"XAZx\0");
        let mut out = Vec::new();
        reencode_aux(&aux, &mut out).unwrap();
        let mut expected = b"NMc\x01ASs".to_vec();
        expected.extend_from_slice(&300u16.to_le_bytes());
        expected.extend_from_slice(b"XAZx\0");
        assert_eq!(out, expected);
        assert_eq!(region_to_bin(100, 250), 4681);
    }

    #[test]
    fn header_order() {
        let h = b"@HD\tVN:1.6\tSO:coordinate\n@PG\tPN:x\tID:p1\n@SQ\tLN:10\tSN:chr1\tM5:a\n@RG\tID:r\tSM:s\n@CO\tc\n";
        assert_eq!(htsjdk_header_text(h), "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\tM5:a\n@RG\tID:r\tSM:s\n@PG\tID:p1\tPN:x\n@CO\tc\n");
    }

    #[test]
    fn region_set() {
        let r = |c: &str, s, e| QueryRegion { chromosome: c.into(), start: s, end: e };
        let java_cmp = java_region_cmp(&r("chr2", 5, 10), &r("chr10", 1, 2));
        assert_eq!(java_cmp, std::cmp::Ordering::Less);
        assert_eq!(java_region_cmp(&r("chr1", 5, 10), &r("chr1", 5, 20)), std::cmp::Ordering::Equal);
    }
}
