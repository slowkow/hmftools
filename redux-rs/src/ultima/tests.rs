//! Port of the Java test consensus/UltimaReadTest, and tests of the attribute strip and of the statistics
//! whose expected values were printed by the Java classes (UltimaRoutines, UltimaStats) for the same input.

use super::*;
use rust_htslib::bam::record::{Aux, AuxArray, Cigar, CigarString};

const LOW_QUAL_BASE: u8 = 12;

/// SeqTechTestUtils.buildUltimaRead: a forward read at chromosome index 0 with NM, tp and t0 (plus the Phred offset)
fn build_ultima_read(read_start: i64, bases: &str, quals: &[u8], tp: &[i8], t0: &[u8]) -> bam::Record {
    let mut rec = bam::Record::new();
    rec.set(b"READ_001", Some(&CigarString(vec![Cigar::Match(quals.len() as u32)])), bases.as_bytes(), quals);
    rec.set_tid(0);
    rec.set_pos(read_start - 1);
    rec.set_mtid(-1);
    rec.set_mpos(-1);
    rec.set_mapq(60);
    rec.set_flags(0);
    rec.push_aux(b"NM", Aux::U8(0)).unwrap();
    rec.push_aux(b"tp", Aux::ArrayI8(AuxArray::from(&tp.to_vec()))).unwrap();
    let t0: Vec<u8> = t0.iter().map(|q| q + 33).collect();
    rec.push_aux(b"t0", Aux::String(std::str::from_utf8(&t0).unwrap())).unwrap();
    rec
}

fn low_qual_tag(rec: &bam::Record) -> Option<String> {
    let t0 = string_attribute(aux_bytes(rec), b"t0").unwrap().unwrap();
    form_low_qual_tag(&rec.seq().as_bytes(), rec.qual(), t0).unwrap()
}

fn set_values(values: &mut [u8], indices: &[usize], value: u8) {
    for &i in indices {
        values[i] = value;
    }
}

/// The attributes of a record as SAM text, in stored order.
fn attributes(rec: &bam::Record) -> String {
    let aux = aux_bytes(rec);
    let mut out = Vec::new();
    let mut i = 0;
    while i + 3 <= aux.len() {
        let n = field_len(&aux[i..]);
        let f = &aux[i + 2..i + n];
        let tag = String::from_utf8_lossy(&aux[i..i + 2]).to_string();
        let text = match f[0] {
            b'Z' => format!("{tag}:Z:{}", String::from_utf8_lossy(&f[1..f.len() - 1])),
            b'B' => {
                let count = u32::from_le_bytes([f[2], f[3], f[4], f[5]]) as usize;
                let values: Vec<String> = (0..count).map(|k| (f[6 + k] as i8).to_string()).collect();
                format!("{tag}:B:{}{}", f[1] as char, values.join(","))
            }
            t => format!("{tag}:{}:{}", t as char, field_int(f).unwrap()),
        };
        out.push(text);
        i += n;
    }
    out.join(" ")
}

/// UltimaReadTest.testUltimaQualTag. The Java test carries @Ignore; its expected values hold for the Java code
/// (checked by running UltimaRoutines.formLowQualTag on the same reads) and for this port.
#[test]
fn test_ultima_qual_tag() {
    let read_start = 100;
    let mut read_bases = "ATTTC";
    let n = read_bases.len();
    let mut base_quals = vec![ULTIMA_MAX_QUAL as u8; n];
    let mut tp_values = vec![0i8; n];
    let mut t0_values = vec![ULTIMA_MAX_QUAL as u8; n];

    // test 0: no quals low enough to trigger creation of the tag
    let read = build_ultima_read(read_start, read_bases, &base_quals, &tp_values, &t0_values);
    assert_eq!(None, low_qual_tag(&read));

    // test 1: quals 35 12 35 12 35, t0 all 35
    set_values(&mut base_quals, &[1, 3], LOW_QUAL_BASE);
    let read = build_ultima_read(read_start, read_bases, &base_quals, &tp_values, &t0_values);
    assert_eq!(Some("3=1-3".to_string()), low_qual_tag(&read));

    // test 2: quals 35 12 35 12 35, t0 35 12 12 12 12
    set_values(&mut t0_values, &[1, 2, 3, 4], LOW_QUAL_BASE);
    let read = build_ultima_read(read_start, read_bases, &base_quals, &tp_values, &t0_values);
    assert_eq!(Some("4=1-4".to_string()), low_qual_tag(&read));

    // test 3: bases A G T T T C, quals 35 35 23 35 23 23, t0 12 12 35 35 35 35
    read_bases = "AGTTTC";
    let n = read_bases.len();
    base_quals = vec![ULTIMA_MAX_QUAL as u8; n];
    tp_values = vec![0i8; n];
    t0_values = vec![ULTIMA_MAX_QUAL as u8; n];
    set_values(&mut base_quals, &[2, 4, 5], 23);
    set_values(&mut t0_values, &[0, 1], LOW_QUAL_BASE);
    let read = build_ultima_read(read_start, read_bases, &base_quals, &tp_values, &t0_values);
    assert_eq!(Some("3=0-1,5".to_string()), low_qual_tag(&read));

    // test 4: successive homopolymers do not affect each other
    read_bases = "ATTTCGGGC";
    let n = read_bases.len();
    base_quals = vec![ULTIMA_MAX_QUAL as u8; n];
    tp_values = vec![0i8; n];
    t0_values = vec![ULTIMA_MAX_QUAL as u8; n];
    set_values(&mut base_quals, &[1, 2, 3], LOW_QUAL_BASE);
    set_values(&mut base_quals, &[5, 6, 7], 23);
    let read = build_ultima_read(read_start, read_bases, &base_quals, &tp_values, &t0_values);
    assert_eq!(Some("3=1-3".to_string()), low_qual_tag(&read));
}

/// A read for the attribute tests: 16 bases at chromosome index 0 (10 bases and unplaced when unmapped).
fn strip_read(flags: u16, bases: &str, quals: &[u8]) -> bam::Record {
    let mut rec = bam::Record::new();
    if flags & FLAG_UNMAPPED != 0 {
        rec.set(b"READ_002", None, bases.as_bytes(), quals);
        rec.set_tid(-1);
        rec.set_pos(-1);
    } else {
        rec.set(b"READ_002", Some(&CigarString(vec![Cigar::Match(quals.len() as u32)])), bases.as_bytes(), quals);
        rec.set_tid(0);
        rec.set_pos(1000);
        rec.set_mapq(60);
    }
    rec.set_mtid(-1);
    rec.set_mpos(-1);
    rec.set_flags(flags);
    rec
}

fn push_t0(rec: &mut bam::Record, t0: &[u8]) {
    let t0: Vec<u8> = t0.iter().map(|q| q + 33).collect();
    rec.push_aux(b"t0", Aux::String(std::str::from_utf8(&t0).unwrap())).unwrap();
}

fn push_tp(rec: &mut bam::Record, value: i8, n: usize) {
    rec.push_aux(b"tp", Aux::ArrayI8(AuxArray::from(&vec![value; n]))).unwrap();
}

/// Strips with UltimaRoutines.stripAttributes(record, stats) and returns the attributes and the counted
/// low-quality bases (None when the read is not counted).
fn strip(rec: &mut bam::Record, drop_qual_tags: bool) -> (String, Option<i32>) {
    let mut routines = UltimaRoutines::new(drop_qual_tags);
    let mut stats = UltimaStats::default();
    routines.strip_attributes(rec, Some(&mut stats)).unwrap();
    (attributes(rec), stats.java_order().first().copied())
}

/// UltimaRoutines.stripAttributes on the reads of a synthetic BAM; the expected attributes are those REDUX 2.1
/// wrote for the same reads, with and without -ultima_drop_qual_tags.
#[test]
fn test_strip_attributes_as_java() {
    let n = 16;
    let max = vec![35u8; n];
    let high_t0 = vec![40u8; n];

    // old ppmSeq lengths balanced at both ends: DUAL; integers take their smallest type; no low-quality bases
    let build = || {
        let mut r = strip_read(0, "CAATGTCTGGAATCTA", &max);
        r.push_aux(b"X1", Aux::I32(1)).unwrap();
        r.push_aux(b"as", Aux::I32(3)).unwrap();
        r.push_aux(b"ts", Aux::U8(4)).unwrap();
        r.push_aux(b"ae", Aux::I32(3)).unwrap();
        r.push_aux(b"te", Aux::I32(3)).unwrap();
        r.push_aux(b"NM", Aux::I32(300)).unwrap();
        r.push_aux(b"AS", Aux::I32(70000)).unwrap();
        r.push_aux(b"MD", Aux::String("16")).unwrap();
        r.push_aux(b"RG", Aux::String("RG1")).unwrap();
        r.push_aux(b"rq", Aux::Float(0.5)).unwrap();
        r.push_aux(b"XS", Aux::I32(-3)).unwrap();
        push_tp(&mut r, 0, n);
        push_t0(&mut r, &high_t0);
        r
    };
    assert_eq!(
        ("t0:Z:IIIIIIIIIIIIIIII RG:Z:RG1 NM:s:300 AS:i:70000 XS:c:-3 UT:Z:DUAL tp:B:c0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0".to_string(), Some(0)),
        strip(&mut build(), false)
    );
    assert_eq!(("RG:Z:RG1 NM:s:300 AS:i:70000 XS:c:-3 UT:Z:DUAL".to_string(), Some(0)), strip(&mut build(), true));

    // st and et both MIXED: DUAL whatever the old lengths say; an unsigned tp array becomes signed; low-quality
    // bases from single bases and a homopolymer
    let build = || {
        let quals = [35, 35, 12, 35, 35, 20, 20, 20, 35, 35, 35, 35, 24, 35, 35, 10];
        let mut r = strip_read(FLAG_REVERSE, "TAACACAAGTGGTGAG", &quals);
        r.push_aux(b"RG", Aux::String("RG1")).unwrap();
        r.push_aux(b"st", Aux::String("MIXED")).unwrap();
        r.push_aux(b"et", Aux::String("MIXED")).unwrap();
        r.push_aux(b"as", Aux::U8(0)).unwrap();
        r.push_aux(b"tp", Aux::ArrayU8(AuxArray::from(&vec![1u8; n]))).unwrap();
        r.push_aux(b"SA", Aux::String("chr21,100,+,8S8M,60,0;")).unwrap();
        r.push_aux(b"NM", Aux::U8(0)).unwrap();
        push_t0(&mut r, &high_t0);
        r
    };
    assert_eq!(
        (
            "t0:Z:IIIIIIIIIIIIIIII SA:Z:chr21,100,+,8S8M,60,0; RG:Z:RG1 NM:c:0 UQ:Z:6=2,5-7,12,15 UT:Z:DUAL tp:B:c1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"
                .to_string(),
            Some(6)
        ),
        strip(&mut build(), false)
    );
    assert_eq!(
        ("SA:Z:chr21,100,+,8S8M,60,0; RG:Z:RG1 NM:c:0 UQ:Z:6=2,5-7,12,15 UT:Z:DUAL".to_string(), Some(6)),
        strip(&mut build(), true)
    );

    // st MIXED and et PLUS: NONE; low t0 values (3 and 19, not 20) make low-quality bases
    let mut t0 = high_t0.clone();
    t0[2] = 3;
    t0[7] = 19;
    t0[8] = 20;
    let mut r = strip_read(0, "GGGGTTGAATTTTAAC", &max);
    r.push_aux(b"st", Aux::String("MIXED")).unwrap();
    r.push_aux(b"et", Aux::String("PLUS")).unwrap();
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    r.push_aux(b"XS", Aux::I32(200)).unwrap();
    r.push_aux(b"AS", Aux::I32(16)).unwrap();
    push_tp(&mut r, -1, n);
    push_t0(&mut r, &t0);
    assert_eq!(
        (
            "t0:Z:II$IIII45IIIIIII RG:Z:RG1 UQ:Z:2=2,7 AS:c:16 XS:C:200 UT:Z:NONE tp:B:c-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1".to_string(),
            Some(2)
        ),
        strip(&mut r, false)
    );

    // only et: the old lengths are not used, NONE; quality 30 everywhere is not low
    let mut r = strip_read(0, "AGGCCCGAAATGTGGT", &vec![30u8; n]);
    r.push_aux(b"et", Aux::String("MIXED")).unwrap();
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    for (tag, v) in [(b"as", 1), (b"ts", 6), (b"ae", 1), (b"te", 6)] {
        r.push_aux(tag, Aux::I32(v)).unwrap();
    }
    push_tp(&mut r, 2, n);
    push_t0(&mut r, &high_t0);
    assert_eq!(("t0:Z:IIIIIIIIIIIIIIII RG:Z:RG1 UT:Z:NONE tp:B:c2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2".to_string(), Some(0)), strip(&mut r, false));

    // a secondary read: no raw qualities, no UT or UQ, not counted, and no need for t0
    let mut r = strip_read(FLAG_SECONDARY, "GTTTGTCTGACTATGA", &max);
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    r.push_aux(b"as", Aux::I32(4)).unwrap();
    r.push_aux(b"NM", Aux::I32(1)).unwrap();
    push_tp(&mut r, 0, n);
    r.push_aux(b"zz", Aux::I32(5)).unwrap();
    assert_eq!(("RG:Z:RG1 NM:c:1".to_string(), None), strip(&mut r, false));

    // a supplementary read: UT and UQ, no raw qualities; the old lengths are unbalanced at the end: NONE
    let mut r = strip_read(FLAG_SUPPLEMENTARY, "CACTTCTGCTTATTCT", &vec![11u8; n]);
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    r.push_aux(b"SA", Aux::String("chr21,5000,+,8S8M,60,0;")).unwrap();
    r.push_aux(b"AS", Aux::I32(40)).unwrap();
    for (tag, v) in [(b"as", 1), (b"ts", 3), (b"ae", 6), (b"te", 2)] {
        r.push_aux(tag, Aux::I32(v)).unwrap();
    }
    push_tp(&mut r, 0, n);
    push_t0(&mut r, &high_t0);
    assert_eq!(("SA:Z:chr21,5000,+,8S8M,60,0; RG:Z:RG1 UQ:Z:16=0-15 AS:c:40 UT:Z:NONE".to_string(), Some(16)), strip(&mut r, false));

    // a sum of old lengths above 8 at the start: NONE; a low-quality run to the last base
    let mut quals = max.clone();
    quals[8..].fill(0);
    let mut r = strip_read(0, "ACTTCTACCTTTTCAG", &quals);
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    for (tag, v) in [(b"as", 8), (b"ts", 1), (b"ae", 2), (b"te", 2)] {
        r.push_aux(tag, Aux::I32(v)).unwrap();
    }
    push_tp(&mut r, 0, n);
    push_t0(&mut r, &high_t0);
    assert_eq!(("t0:Z:IIIIIIIIIIIIIIII RG:Z:RG1 UQ:Z:8=8-15 UT:Z:NONE tp:B:c0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0".to_string(), Some(8)), strip(&mut r, false));

    // homopolymers assessed from their first base below the maximum quality; a ratio at the end below 0.27: NONE
    let quals = [35, 35, 35, 13, 13, 13, 13, 13, 13, 13, 35, 15, 15, 15, 15, 35];
    let mut r = strip_read(0, "CAAGAAATTAAGTATA", &quals);
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    for (tag, v) in [(b"as", 2), (b"ts", 2), (b"ae", 2), (b"te", 6)] {
        r.push_aux(tag, Aux::I32(v)).unwrap();
    }
    push_tp(&mut r, 0, n);
    push_t0(&mut r, &high_t0);
    assert_eq!(("t0:Z:IIIIIIIIIIIIIIII RG:Z:RG1 UQ:Z:12=3-14 UT:Z:NONE tp:B:c0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0".to_string(), Some(12)), strip(&mut r, false));

    // an unmapped read (RegionUnmapper.processFullyUnmappedReads): only the general attributes stay
    let mut r = strip_read(FLAG_UNMAPPED, "ACGTACGTAC", &vec![35u8; 10]);
    push_tp(&mut r, 0, 10);
    r.push_aux(b"XS", Aux::I32(70000)).unwrap();
    r.push_aux(b"as", Aux::I32(3)).unwrap();
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    r.push_aux(b"NM", Aux::I32(256)).unwrap();
    push_t0(&mut r, &vec![40u8; 10]);
    r.push_aux(b"AS", Aux::I32(-200)).unwrap();
    let mut routines = UltimaRoutines::new(false);
    routines.strip_attributes(&mut r, None).unwrap();
    assert_eq!("RG:Z:RG1 NM:s:256 AS:s:-200 XS:i:70000", attributes(&r));

    // a duplicate keeps UT and UQ but not the raw qualities
    let mut r = strip_read(FLAG_DUPLICATE, "CAATGTCTGGAATCTA", &max);
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    push_tp(&mut r, 0, n);
    push_t0(&mut r, &high_t0);
    assert_eq!(("RG:Z:RG1 UT:Z:NONE".to_string(), Some(0)), strip(&mut r, false));
}

/// Where Java throws, the port returns an error: a mapped primary read without t0, or with a short t0.
#[test]
fn test_strip_attributes_java_failures() {
    let mut routines = UltimaRoutines::new(false);
    let mut stats = UltimaStats::default();
    let mut r = strip_read(0, "CAATGTCTGGAATCTA", &[35u8; 16]);
    r.push_aux(b"RG", Aux::String("RG1")).unwrap();
    assert!(routines.strip_attributes(&mut r, Some(&mut stats)).is_err());
    assert_eq!("RG:Z:RG1", attributes(&r));

    let mut r = strip_read(0, "CAATGTCTGGAATCTA", &[35u8; 16]);
    push_t0(&mut r, &[40u8; 15]);
    assert!(routines.strip_attributes(&mut r, Some(&mut stats)).is_err());

    // a longer t0 is read by index (a hard-clipped read)
    let mut r = strip_read(0, "CAATGTCTGGAATCTA", &[35u8; 16]);
    push_t0(&mut r, &[40u8; 20]);
    assert!(routines.strip_attributes(&mut r, Some(&mut stats)).is_ok());
    assert!(stats.java_order() == vec![0]);
}

fn stats_rows(stats: &UltimaStats) -> Vec<(i32, i32)> {
    let path = std::env::temp_dir().join(format!("redux_rs_ultima_stats_test_{}.tsv", std::process::id()));
    stats.write_stats(path.to_str().unwrap()).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let _ = std::fs::remove_file(&path);
    let mut lines = text.lines();
    assert_eq!(Some("LowQualBases\tCount"), lines.next());
    lines.map(|l| l.split_once('\t').unwrap()).map(|(a, b)| (a.parse().unwrap(), b.parse().unwrap())).collect()
}

/// UltimaStats: counts add up, partitions merge, and the rows come in the order of a java.util.HashMap with
/// Integer keys. The expected orders are those of the files REDUX 2.1 wrote for synthetic reads with these counts.
#[test]
fn test_ultima_stats_order_as_java() {
    // 16 shares the first of 16 buckets with 0
    let mut stats = UltimaStats::default();
    for k in [0, 6, 2, 0, 16, 8, 12] {
        stats.add_low_qual_count(k, 1);
    }
    let mut combined = UltimaStats::default();
    combined.merge(&stats);
    assert_eq!(vec![(0, 2), (16, 1), (2, 1), (6, 1), (8, 1), (12, 1)], stats_rows(&combined));

    // 11 counts stay in 16 buckets: by count modulo 16, then by first appearance
    let mut stats = UltimaStats::default();
    for k in [82, 39, 7, 40, 25, 59, 92, 124, 13, 79, 15, 7] {
        stats.add_low_qual_count(k, 1);
    }
    let mut combined = UltimaStats::default();
    combined.merge(&stats);
    let keys: Vec<i32> = stats_rows(&combined).iter().map(|r| r.0).collect();
    assert_eq!(vec![82, 39, 7, 40, 25, 59, 92, 124, 13, 79, 15], keys);
    assert_eq!(2, stats_rows(&combined).iter().find(|r| r.0 == 7).unwrap().1);

    // more than 12 counts double the table: 36 counts are in 64 buckets
    let first_seen = [
        7, 133, 70, 200, 8, 394, 141, 14, 142, 16, 18, 150, 22, 24, 25, 27, 93, 95, 31, 34, 37, 426, 173, 302, 50, 114, 52, 180, 53, 56, 60, 252,
        444, 1085, 126, 383,
    ];
    let mut by_bucket: Vec<i32> = first_seen.to_vec();
    by_bucket.sort_by_key(|k| k % 64);
    let mut stats = UltimaStats::default();
    for k in first_seen {
        stats.add_low_qual_count(k, 3);
    }
    // two partitions that are merged in turn
    let mut other = UltimaStats::default();
    other.add_low_qual_count(383, 1);
    other.add_low_qual_count(7, 1);
    let mut combined = UltimaStats::default();
    combined.merge(&stats);
    combined.merge(&other);
    let rows = stats_rows(&combined);
    assert_eq!(by_bucket, rows.iter().map(|r| r.0).collect::<Vec<i32>>());
    assert_eq!(vec![133, 70, 7, 200, 8, 394, 141, 14, 142, 16, 18, 150], rows.iter().map(|r| r.0).take(12).collect::<Vec<i32>>());
    assert_eq!(4, rows.iter().find(|r| r.0 == 7).unwrap().1);
    assert_eq!(4, rows.iter().find(|r| r.0 == 383).unwrap().1);
    assert_eq!(3, rows.iter().find(|r| r.0 == 1085).unwrap().1);
}

/// R26: with -ultima_drop_qual_tags and the base quality recalibration on, the raw quality attributes stay until
/// the read has been counted and are taken off as it is written; the written read is the one Java writes.
#[test]
fn raw_quals_are_dropped_after_the_recalibration() {
    let n = 16;
    let build = || {
        let mut r = strip_read(0, "CAATGTCTGGAATCTA", &vec![35u8; n]);
        r.push_aux(b"NM", Aux::I32(0)).unwrap();
        push_t0(&mut r, &vec![40u8; n]);
        push_tp(&mut r, 0, n);
        r
    };
    // Java: gone before the recalibration, which then fails
    let mut java = build();
    let (java_attributes, _) = strip(&mut java, true);
    assert!(!java_attributes.contains("t0"));
    assert!(crate::compat::with_java(true, || !UltimaRoutines::for_partition(true, true).drop_after_capture));

    let mut routines = UltimaRoutines::for_partition(true, true);
    assert!(routines.drop_after_capture);
    let mut read = build();
    routines.strip_attributes(&mut read, Some(&mut UltimaStats::default())).unwrap();
    assert!(attributes(&read).contains("t0") && attributes(&read).contains("tp"));
    UltimaRoutines::drop_raw_quals(&mut read);
    assert_eq!(attributes(&read), java_attributes);

    // nothing is held back without the recalibration, or without the option
    assert!(!UltimaRoutines::for_partition(true, false).drop_after_capture);
    assert!(!UltimaRoutines::for_partition(false, true).drop_after_capture);
}
