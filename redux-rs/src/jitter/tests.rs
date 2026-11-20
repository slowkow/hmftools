//! Tests of the read analysis and the counting, on synthetic reads. test_microsatellite_read is the Java
//! MsJitterReadTest.

use super::read::{self, ConsensusType, SiteRef};
use super::*;
use rust_htslib::bam::record::{Aux, CigarString};

const DEFAULT_QUAL: u8 = 37;

fn record(start: i32, cigar: &str, bases: &[u8]) -> bam::Record {
    let mut rec = bam::Record::new();
    let cigar = CigarString::try_from(cigar).unwrap();
    rec.set(b"read", Some(&cigar), bases, &vec![DEFAULT_QUAL; bases.len()]);
    rec.set_tid(0);
    rec.set_pos(start as i64 - 1);
    rec.set_mapq(60);
    rec.set_flags(0);
    rec
}

fn set_qual(rec: &mut bam::Record, index: usize, qual: u8) {
    let bases = rec.seq().as_bytes();
    let mut quals = rec.qual().to_vec();
    quals[index] = qual;
    let (name, cigar) = (rec.qname().to_vec(), rec.cigar().take());
    rec.set(&name, Some(&cigar), &bases, &quals);
}

fn analyse(site: &SiteRef, rec: &bam::Record) -> read::MsRead {
    let cigar = rec.raw_cigar();
    let start = rec.pos() as i32 + 1;
    read::analyse(site, rec, start, start + read::reference_length(cigar) - 1, cigar, true)
}

#[test]
fn test_microsatellite_read() {
    let repeat = b"AAAA";
    let site = SiteRef { start: 50, end: 53, unit: repeat, repeat_count: 1 };
    let flank = b"ACGTACCCCCGGGGGACGTA";
    let bases: Vec<u8> = [&flank[..], &repeat[..], &flank[..]].concat();
    let cigar = format!("{}M", bases.len());

    assert!(analyse(&site, &record(31, &cigar, &bases)).valid);

    // a read too close to the flanks on either side is invalid
    assert!(!analyse(&site, &record(46, &cigar, &bases)).valid);
    assert!(!analyse(&site, &record(8, &cigar, &bases)).valid);

    // low quals within the repeat bounds is invalid
    let mut rec = record(31, &cigar, &bases);
    set_qual(&mut rec, 20, 1);
    assert!(!analyse(&site, &rec).valid);
    set_qual(&mut rec, 20, 40);
    set_qual(&mut rec, 23, 1);
    assert!(!analyse(&site, &rec).valid);

    // deletes covering the repeat
    assert!(!analyse(&site, &record(31, "18M2D26M", &bases)).valid);

    // invalidating inserts
    assert!(!analyse(&site, &record(31, "18M2I24M", &bases)).valid);
}

/// Bases for a read: `len` bases cycling CGT so that nothing but the given ranges is an A.
fn bases_with_a(len: usize, a_ranges: &[(usize, usize)]) -> Vec<u8> {
    (0..len).map(|i| if a_ranges.iter().any(|r| i >= r.0 && i < r.1) { b'A' } else { b"CGT"[i % 3] }).collect()
}

// the test site: eight A at 101-108
const SITE: SiteRef = SiteRef { start: 101, end: 108, unit: b"A", repeat_count: 8 };

#[test]
fn jitter_from_indels() {
    // no indel
    let r = analyse(&SITE, &record(81, "40M", &bases_with_a(40, &[(20, 28)])));
    assert!(r.valid);
    assert_eq!((r.repeat_length, r.repeat_units, r.jitter), (8, 8, 0));

    // one base deleted inside the repeat
    let r = analyse(&SITE, &record(81, "22M1D38M", &bases_with_a(60, &[(20, 27)])));
    assert!(r.valid);
    assert_eq!((r.repeat_length, r.jitter), (7, -1));

    // two bases of the unit inserted
    let r = analyse(&SITE, &record(81, "22M2I18M", &bases_with_a(42, &[(20, 30)])));
    assert!(r.valid);
    assert_eq!((r.repeat_length, r.jitter), (10, 2));

    // an insert right after the repeat counts too
    let r = analyse(&SITE, &record(81, "28M1I20M", &bases_with_a(49, &[(20, 29)])));
    assert!(r.valid);
    assert_eq!(r.jitter, 1);

    // an insert that is not the unit
    assert!(!analyse(&SITE, &record(81, "22M2I18M", &bases_with_a(42, &[(20, 22), (24, 30)]))).valid);

    // a delete reaching the base before the repeat, and one running over its end
    assert!(!analyse(&SITE, &record(81, "18M2D40M", &bases_with_a(58, &[(18, 26)]))).valid);
    assert!(!analyse(&SITE, &record(81, "26M3D40M", &bases_with_a(66, &[(20, 26)]))).valid);

    // an insert in the flank
    assert!(!analyse(&SITE, &record(81, "17M1I30M", &bases_with_a(48, &[(21, 29)]))).valid);

    // a delete in the flank leaves too few aligned flank bases
    assert!(!analyse(&SITE, &record(81, "16M2D40M", &bases_with_a(56, &[(18, 26)]))).valid);

    // soft clips do not count as flank
    assert!(analyse(&SITE, &record(96, "10S25M", &bases_with_a(35, &[(15, 23)]))).valid);
    assert!(!analyse(&SITE, &record(97, "10S25M", &bases_with_a(35, &[(14, 22)]))).valid);

    // a unit of two bases: the repeat length is divided by the unit length
    let site = SiteRef { start: 101, end: 108, unit: b"AC", repeat_count: 4 };
    let r = analyse(&site, &record(81, "22M3D38M", &bases_with_a(60, &[])));
    assert!(r.valid);
    assert_eq!((r.repeat_length, r.repeat_units, r.jitter), (5, 2, -2));
}

#[test]
fn base_quality_filter() {
    let bases = bases_with_a(40, &[(20, 28)]);
    let mut rec = record(81, "40M", &bases);
    // the base before the repeat
    set_qual(&mut rec, 19, 1);
    assert!(!analyse(&SITE, &rec).valid);
    // outside the repeat and its one base flanks
    set_qual(&mut rec, 19, 40);
    set_qual(&mut rec, 18, 1);
    assert!(analyse(&SITE, &rec).valid);
    // without the filter
    set_qual(&mut rec, 22, 0);
    assert!(!analyse(&SITE, &rec).valid);
    let cigar = rec.raw_cigar();
    assert!(read::analyse(&SITE, &rec, 81, 120, cigar, false).valid);
    // SBX medium quals are flagged but valid
    set_qual(&mut rec, 22, 27);
    let r = analyse(&SITE, &rec);
    assert!(r.valid && r.has_medium_quals);
    set_qual(&mut rec, 22, 30);
    assert!(!analyse(&SITE, &rec).has_medium_quals);
}

#[test]
fn consensus_types() {
    let bases = bases_with_a(44, &[(19, 23)]);
    let site = SiteRef { start: 50, end: 53, unit: b"A", repeat_count: 4 };
    let mut rec = record(31, "44M", &bases);
    assert_eq!(read::extract_consensus_type(&rec), ConsensusType::None);
    rec.push_aux(b"UT", Aux::String("DUAL")).unwrap();
    assert_eq!(read::extract_consensus_type(&rec), ConsensusType::Dual);
    let cigar = rec.raw_cigar().to_vec();

    assert_eq!(read::microsatellite_boundaries(50, 53, 31, &cigar), (18, 23));
    // no duplex index: simplex
    assert_eq!(read::sbx_consensus_type(ConsensusType::Dual, &site, &rec, 31, &cigar), ConsensusType::Single);
    assert_eq!(read::sbx_consensus_type(ConsensusType::None, &site, &rec, 31, &cigar), ConsensusType::None);

    // forward read: duplex from the index to the end
    rec.push_aux(b"YX", Aux::I32(18)).unwrap();
    assert_eq!(read::sbx_consensus_type(ConsensusType::Dual, &site, &rec, 31, &cigar), ConsensusType::Dual);
    let mut later = record(31, "44M", &bases);
    later.push_aux(b"YX", Aux::U8(19)).unwrap();
    assert_eq!(read::sbx_consensus_type(ConsensusType::Dual, &site, &later, 31, &cigar), ConsensusType::Single);

    // reverse read: duplex from the start to the index
    let mut reverse = record(31, "44M", &bases);
    reverse.set_flags(0x10);
    reverse.push_aux(b"YX", Aux::I16(23)).unwrap();
    assert_eq!(read::sbx_consensus_type(ConsensusType::Dual, &site, &reverse, 31, &cigar), ConsensusType::Dual);
    let mut reverse = record(31, "44M", &bases);
    reverse.set_flags(0x10);
    reverse.push_aux(b"YX", Aux::I16(22)).unwrap();
    assert_eq!(read::sbx_consensus_type(ConsensusType::Dual, &site, &reverse, 31, &cigar), ConsensusType::Single);

    // boundaries with a soft clip, an insert in the repeat and a delete
    let cigar = record(41, "5S11M2I10M", &bases_with_a(28, &[])).raw_cigar().to_vec();
    assert_eq!(read::microsatellite_boundaries(50, 53, 41, &cigar), (13, 20));
    let cigar = record(41, "10M2D20M", &bases_with_a(30, &[])).raw_cigar().to_vec();
    assert_eq!(read::microsatellite_boundaries(50, 53, 41, &cigar), (8, 11));
}

fn test_config(sequencing_type: SequencingType, write_site_file: bool, output_dir: &str) -> JitterConfig {
    JitterConfig {
        sample_id: "SAMPLE".to_string(),
        output_dir: output_dir.to_string(),
        sequencing_type,
        uses_duplex_umis: false,
        apply_bq_filter: true,
        write_site_file,
        coefficients_file: None,
        error_rates_file: None,
    }
}

fn test_sites() -> Sites {
    let s = |chr: &str, start: i32, end: i32, unit: &str| (chr.to_string(), start, end, unit.as_bytes().to_vec());
    // deliberately out of order within chr1
    Sites::build(vec![
        s("chr1", 301, 308, "A"),
        s("chr1", 101, 108, "A"),
        s("chr2", 101, 108, "AC"),
        s("chr1", 501, 512, "CAG"),
        s("chr1", 701, 706, "T"),
    ])
    .unwrap()
}

fn add_reads(p: &mut MsJitterProcessor, chromosome: &str, count: usize, start: i32, cigar: &str, bases: &[u8], consensus: Option<&str>) {
    for _ in 0..count {
        let mut rec = record(start, cigar, bases);
        if let Some(c) = consensus {
            rec.push_aux(b"UT", Aux::String(c)).unwrap();
        }
        p.process_read(&rec, chromosome);
    }
}

#[test]
fn site_order_and_groups() {
    let sites = test_sites();
    assert_eq!(sites.start, vec![101, 301, 501, 701, 101]);
    assert_eq!(sites.file_index, vec![1, 0, 3, 4, 2]);
    assert_eq!(sites.repeat_count, vec![8, 8, 4, 6, 4]);
    assert_eq!(sites.chromosome_range("chr1"), Some((0, 4)));
    assert_eq!(sites.chromosome_range("chr2"), Some((4, 5)));
    assert_eq!(sites.chromosome_range("chr3"), None);
    assert_eq!(sites.chromosome_of(3), "chr1");
    assert_eq!(sites.chromosome_of(4), "chr2");
    assert_eq!(sites.group, vec![0, 0, 6, 0, 3]);
    assert_eq!(unit_group(b"G"), 1);
    assert_eq!(unit_group(b"TA"), 2);
    assert_eq!(unit_group(b"CT"), 4);
    assert_eq!(unit_group(b"GC"), 5);
    assert_eq!(unit_group(b"AAAT"), 7);
    assert_eq!(unit_group(b"AAAAT"), 8);
    assert_eq!(unit_group(b"AA"), NO_GROUP);
    assert_eq!(unit_group(b"N"), NO_GROUP);
    assert_eq!(unit_group(b"AAAAAT"), NO_GROUP);
    let names: Vec<String> = (0..GROUP_COUNT).map(group_name).collect();
    assert_eq!(names, vec!["A/T", "C/G", "AT/TA", "AC/CA/GT/TG", "AG/GA/CT/TC", "CG/GC", "3bp repeat", "4bp repeat", "5bp repeat"]);
}

#[test]
fn counts_tables_and_site_selection() {
    let analyser = MsJitterAnalyser::from_parts(test_config(SequencingType::Illumina, true, ""), test_sites());
    {
        // two workers share the reads of the first site; sorted within each worker
        let mut p1 = analyser.processor();
        let mut p2 = analyser.processor();
        add_reads(&mut p1, "chr1", 20, 81, "40M", &bases_with_a(40, &[(20, 28)]), None);
        add_reads(&mut p2, "chr1", 10, 85, "40M", &bases_with_a(40, &[(16, 24)]), None);
        add_reads(&mut p2, "chr1", 3, 81, "22M1D38M", &bases_with_a(60, &[(20, 27)]), None);
        // a read that contains the site but not its flanks is rejected
        add_reads(&mut p1, "chr1", 2, 98, "40M", &bases_with_a(40, &[(3, 11)]), None);
        // a read that does not contain the site is not seen at all
        add_reads(&mut p1, "chr1", 2, 102, "40M", &bases_with_a(40, &[]), None);
        // second site: 10 of 30 reads with one more unit is over the 30% limit
        add_reads(&mut p1, "chr1", 20, 281, "40M", &bases_with_a(40, &[(20, 28)]), None);
        add_reads(&mut p1, "chr1", 10, 281, "28M1I20M", &bases_with_a(49, &[(20, 29)]), None);
        // third site (3 base unit): 25 reads, 2 with a jitter of -3, counted outside the dense range
        add_reads(&mut p2, "chr1", 23, 481, "60M", &bases_with_a(60, &[]), None);
        add_reads(&mut p2, "chr1", 2, 481, "22M9D38M", &bases_with_a(60, &[]), None);
        // fourth site: too few reads
        add_reads(&mut p2, "chr1", 5, 681, "40M", &bases_with_a(40, &[]), None);
        // chr2 after chr1, then back to an earlier position of chr1 as a new region would
        add_reads(&mut p1, "chr2", 21, 81, "40M", &bases_with_a(40, &[]), Some("DUAL"));
        add_reads(&mut p1, "chr1", 1, 81, "40M", &bases_with_a(40, &[(20, 28)]), None);
        // a chromosome without sites, and an unmapped read
        add_reads(&mut p1, "chr5", 3, 81, "40M", &bases_with_a(40, &[]), None);
        let mut unmapped = record(81, "40M", &bases_with_a(40, &[]));
        unmapped.set_flags(0x4);
        p1.process_read(&unmapped, "chr1");
    }
    let totals = analyser.totals.lock().unwrap();
    assert_eq!(totals.passing, vec![34, 30, 25, 5, 21]);
    assert_eq!(totals.read_count(0, ConsensusType::None), 34);
    assert_eq!(totals.rejected_count(0, ConsensusType::None), 2);
    assert_eq!(totals.jitter_count(0, ConsensusType::None, 0), 31);
    assert_eq!(totals.jitter_count(0, ConsensusType::None, -1), 3);
    assert_eq!(totals.jitter_count(1, ConsensusType::None, 1), 10);
    assert_eq!(totals.jitter_count(2, ConsensusType::None, -3), 2);
    assert_eq!(totals.read_count(4, ConsensusType::None), 0);
    assert_eq!(totals.read_count(4, ConsensusType::Dual), 21);
    assert!(totals.by_type[ConsensusType::Single as usize].is_empty());

    let real_variant: Vec<bool> = (0..5).map(|i| totals.is_real_variant(i)).collect();
    assert_eq!(real_variant, vec![false, true, false, true, false]);

    let tables = analyser.build_group_tables(&totals, &real_variant, ConsensusType::None);
    assert_eq!(tables.len(), 9);
    // A/T: only the first site is kept
    assert_eq!(tables[0].rows.len(), 1);
    assert_eq!(tables[0].rows[0].ref_num_units, 8);
    assert_eq!(tables[0].rows[0].total_read_count, 34);
    assert_eq!(tables[0].rows[0].jitter_read_count(-1), 3);
    // the kept dual-only site adds an empty row to the NONE table
    assert_eq!(tables[3].rows.len(), 1);
    assert_eq!(tables[3].rows[0].total_read_count, 0);
    assert!(tables[3].rows[0].jitter_counts.is_empty());
    assert_eq!(tables[6].rows[0].ref_num_units, 4);
    assert_eq!(tables[6].rows[0].jitter_read_count(-3), 2);
    let dual = analyser.build_group_tables(&totals, &real_variant, ConsensusType::Dual);
    assert_eq!(dual[3].rows[0].total_read_count, 21);

    // the site file: sorted by chromosome, start; NONE only for Illumina without duplex UMIs
    let mut out = Vec::new();
    analyser.write_site_file(&mut out, &totals, &real_variant, &analyser.config.consensus_types()).unwrap();
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 6);
    assert!(lines[0].starts_with("chromosome\tstart\tend\tunit\tconsensusType\tnumReads\tnumReadsRejected\trealVariant\tcount-10\t"));
    assert!(lines[0].ends_with("\tcount+10\treadRepeatLengths"));
    let lengths = format!("{},{}", vec!["7"; 3].join(","), vec!["8"; 31].join(","));
    assert_eq!(lines[1], format!("chr1\t101\t108\tA\tNONE\t34\t2\tfalse\t0\t0\t0\t0\t0\t0\t0\t0\t0\t3\t31\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t{lengths}"));
    assert!(lines[2].starts_with("chr1\t301\t308\tA\tNONE\t30\t0\ttrue\t"));
    assert!(lines[3].starts_with("chr1\t501\t512\tCAG\tNONE\t25\t0\tfalse\t0\t0\t0\t0\t0\t0\t0\t2\t0\t0\t23\t"));
    assert!(lines[5].starts_with("chr2\t101\t108\tAC\tNONE\t0\t0\tfalse\t"));
    assert!(lines[5].ends_with("\t0\t"));
}

#[test]
fn jitter_outside_dense_range_and_window_moves() {
    // many sites so that a worker's window has to move
    let count = 3 * MAX_WINDOW_SITES as i32;
    let records = (0..count).map(|i| ("chr1".to_string(), 1000 + 100 * i, 1007 + 100 * i, b"A".to_vec())).collect();
    let analyser = MsJitterAnalyser::from_parts(test_config(SequencingType::Illumina, false, ""), Sites::build(records).unwrap());
    {
        let mut p = analyser.processor();
        for i in 0..count {
            let n = if i % 1000 == 0 { 3 } else { 1 };
            add_reads(&mut p, "chr1", n, 980 + 100 * i, "40M", &bases_with_a(40, &[(20, 28)]), None);
        }
        // seven bases deleted, then back to the first site with eight units inserted
        add_reads(&mut p, "chr1", 2, 980 + 100 * (count - 1), "21M7D30M", &bases_with_a(51, &[]), None);
        add_reads(&mut p, "chr1", 4, 980, "22M8I18M", &bases_with_a(48, &[(20, 36)]), None);
    }
    let totals = analyser.totals.lock().unwrap();
    assert_eq!(totals.passing.iter().map(|c| *c as usize).sum::<usize>(), count as usize + 2 * 13 + 6);
    assert_eq!(totals.passing[0], 7);
    assert_eq!(totals.passing[1000], 3);
    assert_eq!(totals.passing[count as usize - 1], 3);
    assert_eq!(totals.jitter_count(0, ConsensusType::None, 8), 4);
    assert_eq!(totals.jitter_count(count as usize - 1, ConsensusType::None, -7), 2);
    assert_eq!(totals.jitter_count(count as usize - 1, ConsensusType::None, 0), 1);
    let mut seen = Vec::new();
    totals.for_each_jitter(0, ConsensusType::None, |j, c| seen.push((j, c)));
    assert_eq!(seen, vec![(0, 3), (8, 4)]);
}

#[test]
fn real_variant_thresholds() {
    let records = vec![("chr1".to_string(), 101, 108, b"A".to_vec())];
    let run = |unchanged: usize, plus_one: usize, minus_three: usize, dual_plus_one: usize| {
        let analyser = MsJitterAnalyser::from_parts(test_config(SequencingType::Illumina, false, ""), Sites::build(records.clone()).unwrap());
        {
            let mut p = analyser.processor();
            add_reads(&mut p, "chr1", unchanged, 81, "40M", &bases_with_a(40, &[(20, 28)]), None);
            add_reads(&mut p, "chr1", plus_one, 81, "28M1I20M", &bases_with_a(49, &[(20, 29)]), None);
            add_reads(&mut p, "chr1", minus_three, 81, "22M3D38M", &bases_with_a(60, &[(20, 25)]), None);
            add_reads(&mut p, "chr1", dual_plus_one, 81, "28M1I20M", &bases_with_a(49, &[(20, 29)]), Some("DUAL"));
        }
        let totals = analyser.totals.lock().unwrap();
        totals.is_real_variant(0)
    };
    // under 20 passing reads
    assert!(run(19, 0, 0, 0));
    assert!(!run(20, 0, 0, 0));
    // jitter of one: over 30% of the passing reads
    assert!(!run(70, 30, 0, 0));
    assert!(run(69, 31, 0, 0));
    // jitter of three: over 20%
    assert!(!run(80, 0, 20, 0));
    assert!(run(79, 0, 21, 0));
    // the reads of all consensus types are summed, also types the run does not report
    assert!(run(69, 16, 0, 15));
}

#[test]
fn sbx_counts() {
    let records = vec![("chr1".to_string(), 101, 108, b"A".to_vec())];
    let analyser = MsJitterAnalyser::from_parts(test_config(SequencingType::Sbx, false, ""), Sites::build(records).unwrap());
    {
        let mut p = analyser.processor();
        let bases = bases_with_a(40, &[(20, 28)]);
        // duplex over the repeat
        for _ in 0..5 {
            let mut rec = record(81, "40M", &bases);
            rec.push_aux(b"UT", Aux::String("DUAL")).unwrap();
            rec.push_aux(b"YX", Aux::I32(10)).unwrap();
            p.process_read(&rec, "chr1");
        }
        // duplex read, repeat in the simplex part: passing only
        for _ in 0..4 {
            let mut rec = record(81, "40M", &bases);
            rec.push_aux(b"UT", Aux::String("DUAL")).unwrap();
            rec.push_aux(b"YX", Aux::I32(30)).unwrap();
            p.process_read(&rec, "chr1");
        }
        // medium quality base in the repeat: passing only
        for _ in 0..3 {
            let mut rec = record(81, "40M", &bases);
            set_qual(&mut rec, 24, 28);
            p.process_read(&rec, "chr1");
        }
        // no consensus
        for _ in 0..2 {
            p.process_read(&record(81, "40M", &bases), "chr1");
        }
        // rejected simplex-part read
        let mut rec = record(98, "40M", &bases);
        rec.push_aux(b"UT", Aux::String("DUAL")).unwrap();
        p.process_read(&rec, "chr1");
    }
    let totals = analyser.totals.lock().unwrap();
    assert_eq!(totals.passing[0], 14);
    assert_eq!(totals.read_count(0, ConsensusType::Dual), 5);
    assert_eq!(totals.read_count(0, ConsensusType::Single), 0);
    assert_eq!(totals.rejected_count(0, ConsensusType::Single), 1);
    assert_eq!(totals.read_count(0, ConsensusType::None), 2);
    assert_eq!(analyser.config.consensus_types(), vec![ConsensusType::None, ConsensusType::Dual]);
}

#[test]
fn consensus_types_by_sequencing_type() {
    let mut config = test_config(SequencingType::Illumina, false, "");
    assert_eq!(config.consensus_types(), vec![ConsensusType::None]);
    config.uses_duplex_umis = true;
    assert_eq!(config.consensus_types(), vec![ConsensusType::None, ConsensusType::Single, ConsensusType::Dual]);
    config.sequencing_type = SequencingType::Ultima;
    assert_eq!(config.consensus_types(), vec![ConsensusType::None, ConsensusType::Dual]);
}

#[test]
fn writes_output_files() {
    let dir = std::env::temp_dir().join(format!("redux_rs_jitter_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let output_dir = format!("{}/", dir.to_str().unwrap());
    let coefficients = format!("{output_dir}coef.tsv");
    let error_rates = format!("{output_dir}rates.tsv");
    std::fs::write(&coefficients, "unit\tnumUnits\tcoefficient1\tcoefficient2\nA/T\t10\t100\t2000\n").unwrap();
    std::fs::write(&error_rates, "unit\tnumUnits\terrorRate\nA/T\t10\t0.01\n").unwrap();
    let mut config = test_config(SequencingType::Illumina, true, &output_dir);
    config.coefficients_file = Some(coefficients);
    config.error_rates_file = Some(error_rates);

    let records = vec![("chr1".to_string(), 101, 110, b"A".to_vec()), ("chr1".to_string(), 301, 304, b"T".to_vec())];
    let analyser = MsJitterAnalyser::from_parts(config, Sites::build(records).unwrap());
    {
        let mut p = analyser.processor();
        add_reads(&mut p, "chr1", 540, 81, "40M", &bases_with_a(40, &[(20, 30)]), None);
        add_reads(&mut p, "chr1", 60, 81, "22M1D38M", &bases_with_a(60, &[(20, 29)]), None);
        add_reads(&mut p, "chr1", 25, 281, "40M", &bases_with_a(40, &[]), None);
    }
    analyser.write_analysis_output().unwrap();

    let mut table = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(format!("{output_dir}SAMPLE.redux.ms_table.tsv.gz")).unwrap())
        .read_to_string(&mut table)
        .unwrap();
    let lines: Vec<&str> = table.lines().collect();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[1], "A/T\tNONE\t4\t25\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t25\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");
    assert_eq!(lines[2], "A/T\tNONE\t10\t600\t0\t0\t0\t0\t0\t0\t0\t0\t0\t60\t540\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");

    let params = std::fs::read_to_string(format!("{output_dir}SAMPLE.redux.jitter_params.tsv")).unwrap();
    let lines: Vec<&str> = params.lines().collect();
    assert_eq!(lines.len(), 8);
    assert_eq!(lines[0], "unit\tconsensusType\toptimalScaleRepeat4\toptimalScaleRepeat5\toptimalScaleRepeat6\tscaleFitGradient\tscaleFitIntercept\tmicrosatelliteSkew");
    assert_eq!(lines[1], "A/T\tNONE\t0.1\t0.15\t0.2\t0.06\t-0.16\t1");
    assert_eq!(lines[3], "AT/TA\tNONE\t0.15\t0.2\t0.25\t0.06\t-0.11\t1");
    assert_eq!(lines[7], "3-5bp repeat\tNONE\t0.2\t0.25\t0.3\t0.06\t-0.06\t1");

    // error rate 60 / 600 less the noise 0.01
    let x: f64 = 0.1 - 0.01;
    let expected = 100.0 * x + 2000.0 * fit::java_pow(x, 3.0);
    let prediction = std::fs::read_to_string(format!("{output_dir}SAMPLE.redux.msi_prediction.tsv")).unwrap();
    assert_eq!(prediction, format!("PredictedValue\n{}\n", crate::fmt::java_fixed(expected, 4)));

    assert!(std::path::Path::new(&format!("{output_dir}SAMPLE.repeat.tsv.gz")).exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

// The remaining tests compare with values printed by the Java classes (Math, OLSMultipleLinearRegression,
// JitterModelFitter, MsModelCalculator) for inputs made by the same generator.

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

fn mix(sum: u64, v: f64) -> u64 {
    sum.rotate_left(7).wrapping_add(v.to_bits())
}

#[test]
fn math_functions_match_java() {
    let mut rng = Lcg(42);
    let (mut s_exp, mut s_log, mut s_log_small, mut s_pow) = (0u64, 0u64, 0u64, 0u64);
    for _ in 0..500000 {
        s_exp = mix(s_exp, fit::java_exp(-750.0 * rng.next() + 5.0));
        s_log = mix(s_log, fit::java_log10(rng.next()));
        s_log_small = mix(s_log_small, fit::java_log10(1e-4 + rng.next() * 1e-3));
        s_pow = mix(s_pow, fit::java_pow(rng.next(), 3.0));
    }
    assert_eq!(s_exp, 0xfb8785ad15786ecb);
    assert_eq!(s_log, 0x5ca80d5eccd9752b);
    assert_eq!(s_log_small, 0x059a4c7666979f40);
    assert_eq!(s_pow, 0x0a6cce2af54e1bf5);
}

#[test]
fn regression_matches_java() {
    let mut rng = Lcg(7);
    let expected: [[u64; 2]; 3] = [
        [0xbfe0dfa0a74fba06, 0x40129892a39e8e21],
        [0xbfbd05fdfd597d29, 0x3ffa2bc568f31eb7],
        [0x3f789028a45d0a3f, 0x3fe01454ebbca65d],
    ];
    for (t, want) in expected.iter().enumerate() {
        let n = 3 + t * 3;
        let mut y = Vec::new();
        let mut x = Vec::new();
        for i in 0..n {
            let w = rng.next() * 0.2;
            y.push(w * (0.1 + rng.next()));
            x.push([w * (7 + i) as f64, w]);
        }
        let p = fit::ols_no_intercept(&y, &x).unwrap();
        assert_eq!([p[0].to_bits(), p[1].to_bits()], *want);
    }
}

fn reference_tables() -> Vec<JitterCountsTable> {
    let units = ["A/T", "AC/CA/GT/TG", "3-5bp repeat"];
    let mut rng = Lcg(99);
    let mut tables = Vec::new();
    for (t, unit) in units.iter().enumerate() {
        let mut table = JitterCountsTable::new(unit, ConsensusType::None);
        for num_units in 4..=18 {
            if t == 2 && num_units > 9 {
                break;
            }
            let row = table.get_or_create_row(num_units);
            let total = (200000.0 * fit::java_pow(0.75, (num_units - 4) as f64)) as i32 + (rng.next() * 1000.0) as i32;
            let spread = 0.002 * (num_units - 3) as f64 * (1 + t) as f64;
            for j in -6i32..=6 {
                if j == 0 {
                    continue;
                }
                let c = (total as f64
                    * spread
                    * fit::java_pow(0.2 + 0.03 * num_units as f64, j.abs() as f64)
                    * (if j < 0 { 1.4 } else { 0.8 })
                    * (0.9 + 0.2 * rng.next())) as i32;
                if c > 0 {
                    row.add_reads(j, c);
                }
            }
            row.add_reads(0, total);
        }
        tables.push(table);
    }
    tables
}

#[test]
fn fit_matches_java() {
    let tables = reference_tables();
    let totals: Vec<i32> = tables[0].rows.iter().map(|r| r.total_read_count).collect();
    assert_eq!(totals, vec![200632, 151593, 113454, 86260, 64764, 48946, 37380, 28717, 21471, 16717, 12301, 9611, 8035, 5775, 5179]);
    let totals: Vec<i32> = tables[2].rows.iter().map(|r| r.total_read_count).collect();
    assert_eq!(totals, vec![201380, 152273, 116186, 88146, 67408, 51705]);

    let expected: [[u64; 6]; 3] = [
        [0x3fc1eb851eb851eb, 0x3fc47ae147ae147b, 0x3fc5c28f5c28f5c3, 0x3fa1ad0016f6df74, 0xbfb63c430a881df1, 0x3ff546a147b404a4],
        [0x3fc3333333333333, 0x3fc5c28f5c28f5c3, 0x3fc851eb851eb853, 0x3fb1b78df756a5c3, 0xbfd3a94e36eec1a8, 0x3ff6a74ddc0d1875],
        [0x3fc47ae147ae147b, 0x3fc851eb851eb853, 0x3fc999999999999b, 0x3fb999999999999b, 0xbfe01b4e81b4e81b, 0x3ff2666666666669],
    ];
    let grid = fit::SearchGrid::new();
    for (table, want) in tables.iter().zip(&expected) {
        let p = fit::perform_fit(&grid, table).unwrap();
        let got = [
            p.optimal_scale_repeat4.to_bits(),
            p.optimal_scale_repeat5.to_bits(),
            p.optimal_scale_repeat6.to_bits(),
            p.scale_fit_gradient.to_bits(),
            p.scale_fit_intercept.to_bits(),
            p.microsatellite_skew.to_bits(),
        ];
        assert_eq!(got, *want, "{}", table.repeat_unit);
    }
}

#[test]
fn fit_fails_with_two_rows_for_the_regression() {
    // two repeat counts of 7 or more, both with at least 20,000 reads: Java's regression throws
    let mut table = JitterCountsTable::new("A/T", ConsensusType::None);
    for num_units in [7, 8] {
        let row = table.get_or_create_row(num_units);
        row.add_reads(0, 30000);
        row.add_reads(-1, 600);
        row.add_reads(1, 300);
    }
    let grid = fit::SearchGrid::new();
    assert!(fit::perform_fit(&grid, &table).is_err());
    table.get_or_create_row(9).add_reads(0, 100);
    assert!(fit::perform_fit(&grid, &table).is_ok());
}

#[test]
fn msi_model_matches_java() {
    let mut calc = ms_model::MsModelCalculator::new(ms_model::MsModelParams::for_sequencing_type(SequencingType::Illumina));
    for n in 10..=20 {
        calc.add_coefficients("A/T", n, 1000.5 / (n - 9) as f64, 30000.25 / (n - 9) as f64);
        calc.add_error_rate("A/T", n, 0.0001 * (n - 9) as f64);
    }
    let predicted = calc.calc_predicted_msi_rate(&reference_tables()).unwrap();
    assert_eq!(predicted.to_bits(), 0x404e90f16b714215);
    assert_eq!(crate::fmt::java_fixed(predicted, 4), "61.1324");
}
