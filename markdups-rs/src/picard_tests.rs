//! Picard 3.3.0's MarkDuplicates tests (AbstractMarkDuplicatesCommandLineProgramTest, MarkDuplicatesTest),
//! run through the whole port. The tester follows Picard's SamFileTester: records built as
//! SAMRecordSetBuilder builds them (names READ<n> unless given, read length 36 unless set, base quality
//! 10, and MarkDuplicatesTester's TOTAL_MAPPED_REFERENCE_LENGTH scoring), written
//! coordinate sorted, marked, and every record's duplicate flag checked against the test's expectation.

use crate::bam::{self, Rec, SamHeader, SqRecord};
use crate::bgzf::{BgzfReader, BgzfWriter};
use crate::input::Inputs;

const DEFAULT_BASE_QUALITY: u8 = 10;

struct Tester {
    records: Vec<(Vec<u8>, bool)>,
    counter: usize,
    read_length: usize,
    expected_optical: Option<i64>,
    args: Vec<String>,
}

fn cigar_ops(c: &str) -> Vec<u32> {
    let mut ops = Vec::new();
    let mut n = 0u32;
    for ch in c.chars() {
        if let Some(d) = ch.to_digit(10) {
            n = n * 10 + d;
        } else {
            let op = "MIDNSHP=X".find(ch).unwrap() as u32;
            ops.push(n << 4 | op);
            n = 0;
        }
    }
    ops
}

#[allow(clippy::too_many_arguments)]
fn record(name: &str, flag: u16, tid: i32, start: i32, cigar: Option<&str>, mtid: i32, mstart: i32, len: usize, qual: u8) -> Vec<u8> {
    let ops = cigar.map(cigar_ops).unwrap_or_default();
    // the bases a CIGAR covers (SAMRecordSetBuilder makes reads as long as their CIGAR)
    let from_cigar: u32 = ops.iter().filter(|&&o| matches!(o & 0xf, 0 | 1 | 4 | 7 | 8)).map(|&o| o >> 4).sum();
    let len = if ops.is_empty() { len } else { from_cigar as usize };
    let mut b = Vec::new();
    b.extend_from_slice(&tid.to_le_bytes());
    b.extend_from_slice(&(start - 1).to_le_bytes());
    b.push((name.len() + 1) as u8);
    b.push(if flag & 4 != 0 { 0 } else { 255 });
    b.extend_from_slice(&0u16.to_le_bytes());
    b.extend_from_slice(&(ops.len() as u16).to_le_bytes());
    b.extend_from_slice(&flag.to_le_bytes());
    b.extend_from_slice(&(len as i32).to_le_bytes());
    b.extend_from_slice(&mtid.to_le_bytes());
    b.extend_from_slice(&(mstart - 1).to_le_bytes());
    b.extend_from_slice(&0i32.to_le_bytes());
    b.extend_from_slice(name.as_bytes());
    b.push(0);
    for o in &ops {
        b.extend_from_slice(&o.to_le_bytes());
    }
    b.extend(std::iter::repeat_n(0x12u8, len.div_ceil(2)));
    b.extend(std::iter::repeat_n(qual, len));
    b.extend_from_slice(b"RGZ1\0");
    b
}

impl Tester {
    fn new() -> Self {
        // MarkDuplicatesTester's default scoring
        Tester {
            records: Vec::new(),
            counter: 0,
            read_length: 36,
            expected_optical: None,
            args: vec!["--DUPLICATE_SCORING_STRATEGY".into(), "TOTAL_MAPPED_REFERENCE_LENGTH".into()],
        }
    }

    fn name(&mut self) -> String {
        let n = format!("READ{}", self.counter);
        self.counter += 1;
        n
    }

    fn cigar(&self, c: Option<&str>) -> String {
        c.map(str::to_string).unwrap_or_else(|| format!("{}M", self.read_length))
    }

    fn add_unmapped_fragment(&mut self, tid: i32, q: u8) {
        let n = self.name();
        self.records.push((record(&n, 4, tid, 0, None, -1, 0, self.read_length, q), false));
    }

    fn add_mapped_fragment(&mut self, tid: i32, start: i32, dup: bool, q: u8) {
        self.add_fragment(None, tid, start, dup, None, q, false, false, false);
    }

    #[allow(clippy::too_many_arguments)]
    fn add_fragment(&mut self, name: Option<&str>, tid: i32, start: i32, dup: bool, cigar: Option<&str>, q: u8, secondary: bool, supplementary: bool, reverse: bool) {
        let n = name.map(str::to_string).unwrap_or_else(|| self.name());
        let c = self.cigar(cigar);
        let flag = if secondary { 0x100 } else { 0 } | if supplementary { 0x800 } else { 0 } | if reverse { 0x10 } else { 0 };
        self.records.push((record(&n, flag, tid, start, Some(&c), -1, 0, self.read_length, q), dup));
    }

    fn add_unmapped_pair(&mut self, tid: i32, q: u8) {
        let n = self.name();
        let f = 1 | 4 | 8;
        self.records.push((record(&n, f | 0x40, tid, 0, None, tid, 0, self.read_length, q), false));
        self.records.push((record(&n, f | 0x80, tid, 0, None, tid, 0, self.read_length, q), false));
    }

    fn add_mapped_pair(&mut self, tid: i32, s1: i32, s2: i32, d1: bool, d2: bool, q: u8) {
        self.add_mate_pair(None, tid, tid, s1, s2, false, false, d1, d2, None, None, false, true, false, false, false, q);
    }

    #[allow(clippy::too_many_arguments)]
    fn add_pair_cigars(&mut self, tid: i32, s1: i32, s2: i32, d1: bool, d2: bool, c1: &str, c2: &str, strand1: bool, strand2: bool, q: u8) {
        self.add_mate_pair(None, tid, tid, s1, s2, false, false, d1, d2, Some(c1), Some(c2), strand1, strand2, false, false, false, q);
    }

    /// SamFileTester.addMatePair (SAMRecordSetBuilder.addPair): strand true is the negative strand; an
    /// unmapped mate is placed at its mapped mate's position
    #[allow(clippy::too_many_arguments)]
    fn add_mate_pair(
        &mut self,
        name: Option<&str>,
        tid1: i32,
        tid2: i32,
        s1: i32,
        s2: i32,
        u1: bool,
        u2: bool,
        d1: bool,
        d2: bool,
        c1: Option<&str>,
        c2: Option<&str>,
        strand1: bool,
        strand2: bool,
        first_only: bool,
        np1: bool,
        np2: bool,
        q: u8,
    ) {
        let n = name.map(str::to_string).unwrap_or_else(|| self.name());
        let (c1, c2) = (self.cigar(c1), self.cigar(c2));
        let (mut tid1, mut s1, mut tid2, mut s2) = (tid1, s1, tid2, s2);
        if u1 {
            tid1 = tid2;
            s1 = s2;
        }
        if u2 {
            tid2 = tid1;
            s2 = s1;
        }
        let mut f1: u16 = 1 | 0x40;
        let mut f2: u16 = 1 | 0x80;
        if !u1 && !u2 {
            f1 |= 2;
            f2 |= 2;
        }
        if u1 {
            f1 |= 4;
            f2 |= 8;
        }
        if u2 {
            f2 |= 4;
            f1 |= 8;
        }
        if strand1 && !u1 {
            f1 |= 0x10;
            f2 |= 0x20;
        }
        if strand2 && !u2 {
            f2 |= 0x10;
            f1 |= 0x20;
        }
        if np1 {
            f1 |= 0x100;
        }
        if np2 {
            f2 |= 0x100;
        }
        let len = self.read_length;
        self.records.push((record(&n, f1, tid1, s1, if u1 { None } else { Some(&c1) }, tid2, s2, len, q), d1));
        if !first_only {
            self.records.push((record(&n, f2, tid2, s2, if u2 { None } else { Some(&c2) }, tid1, s1, len, q), d2));
        }
    }

    /// Writes the records coordinate sorted, runs MarkDuplicates and checks the flags; returns the
    /// metrics text.
    fn run(mut self) -> String {
        let dir = std::env::temp_dir().join(format!("mdrs_picard_{}_{:?}", std::process::id(), std::thread::current().id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("in.bam").to_string_lossy().to_string();
        let output = dir.join("out.bam").to_string_lossy().to_string();
        let metrics = dir.join("out.metrics").to_string_lossy().to_string();
        let mut h = SamHeader::default();
        h.set_hd(b"SO", Some(b"coordinate"));
        for i in 0..8 {
            h.sq.push(SqRecord { name: format!("chr{}", i + 1).into_bytes(), len: 200_000_000, attrs: Vec::new() });
        }
        h.rg.push(bam::IdRecord { id: b"1".to_vec(), attrs: vec![(b"SM".to_vec(), b"s".to_vec()), (b"LB".to_vec(), b"lib".to_vec())] });
        // SAMFileWriter's coordinate sort: SAMRecordCoordinateComparator, stable
        self.records.sort_by(|a, b| crate::input::coordinate_compare(&Rec(&a.0), &Rec(&b.0)));
        let mut w = BgzfWriter::create(&input, 1, 0).unwrap();
        bam::write_header(&mut w, &h).unwrap();
        for (r, _) in &self.records {
            w.write(&(r.len() as u32).to_le_bytes()).unwrap();
            w.write(r).unwrap();
        }
        w.finish().unwrap();
        drop(w);
        let mut args: Vec<String> = ["--INPUT", &input, "--OUTPUT", &output, "--METRICS_FILE", &metrics].iter().map(|s| s.to_string()).collect();
        args.extend(self.args.iter().cloned());
        let parsed = crate::args::parse(&args, 2).unwrap();
        let cfg = crate::config(&parsed, 1).unwrap();
        let inputs = Inputs { paths: vec![input.clone()], assume_sort_order: None };
        let mut d = crate::markdups::find_duplicates(&cfg, &inputs, 0).unwrap();
        crate::markdups::write_output(&cfg, &inputs, &mut d, 0, 0).unwrap();
        let mut r = BgzfReader::open(&output, 0).unwrap();
        bam::read_header(&mut r).unwrap();
        let mut buf = Vec::new();
        let mut i = 0;
        while r.read_record(&mut buf).unwrap() {
            let rec = Rec(&buf);
            let (_, want) = &self.records[i];
            assert_eq!(rec.flag() & 0x400 != 0, *want, "record {i} ({}) duplicate flag", String::from_utf8_lossy(rec.name()));
            i += 1;
        }
        assert_eq!(i, self.records.len());
        let text = std::fs::read_to_string(&metrics).unwrap();
        if let Some(o) = self.expected_optical {
            let row = text.lines().skip_while(|l| !l.starts_with("LIBRARY")).nth(1).unwrap();
            let got: i64 = row.split('\t').nth(7).unwrap().parse().unwrap();
            assert_eq!(got, o, "optical duplicates");
        }
        std::fs::remove_dir_all(&dir).ok();
        text
    }
}

const Q: u8 = DEFAULT_BASE_QUALITY;

#[test]
fn single_and_two_unmapped_fragments_and_pair() {
    let mut t = Tester::new();
    t.add_unmapped_fragment(-1, Q);
    t.run();
    let mut t = Tester::new();
    t.add_unmapped_fragment(-1, Q);
    t.add_unmapped_fragment(-1, Q);
    t.run();
    let mut t = Tester::new();
    t.add_unmapped_pair(-1, Q);
    t.run();
}

#[test]
fn mapped_fragments() {
    let mut t = Tester::new();
    t.add_mapped_fragment(1, 1, false, Q);
    t.run();
    let mut t = Tester::new();
    t.add_mapped_fragment(0, 1, false, Q);
    t.add_mapped_fragment(0, 1, true, Q);
    t.run();
}

#[test]
fn mapped_pairs() {
    let mut t = Tester::new();
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.run();
    let mut t = Tester::new();
    t.add_mapped_fragment(1, 1, true, Q);
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.run();
    let mut t = Tester::new();
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.run();
    let mut t = Tester::new();
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.run();
    let mut t = Tester::new();
    t.add_mapped_fragment(1, 1, true, Q);
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.run();
}

#[test]
fn terminal_unmapped_records() {
    let mut t = Tester::new();
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.add_unmapped_fragment(-1, Q);
    t.run();
    let mut t = Tester::new();
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.add_unmapped_pair(-1, Q);
    t.run();
}

#[test]
fn optical_duplicate_finding() {
    let mut t = Tester::new();
    t.expected_optical = Some(1);
    t.add_mate_pair(Some("READ0:1:1:1:1"), 1, 1, 1, 100, false, false, false, false, Some("50M"), Some("50M"), false, true, false, false, false, Q);
    t.add_mate_pair(Some("READ1:1:1:1:300"), 1, 1, 1, 100, false, false, true, true, Some("50M"), Some("50M"), false, true, false, false, false, Q);
    t.add_mate_pair(Some("READ2:1:1:1:50"), 1, 1, 1, 100, false, false, true, true, Some("50M"), Some("50M"), false, true, false, false, false, Q);
    t.run();
}

/// testOpticalDuplicateClusterSamePositionNoOpticalDuplicates and the pixel distance variants: two pairs
/// whose names are far apart (or on other tiles) are no optical duplicates; within the distance one is
#[test]
fn optical_cluster_same_position() {
    let pairs = [
        ("RUNID:7:1203:2886:16756", "RUNID:7:1205:3886:16756", 0),
        ("RUNID:7:1203:2886:16756", "RUNID:7:1203:4886:26756", 0),
        ("RUNID:7:1203:2886:34756", "RUNID:7:1203:2886:36756", 0),
        ("RUNID:7:1203:2886:35756", "RUNID:7:1203:2886:35800", 1),
    ];
    for (a, b, optical) in pairs {
        let mut t = Tester::new();
        t.expected_optical = Some(optical);
        t.add_mate_pair(Some(a), 1, 1, 1, 100, false, false, false, false, Some("50M"), Some("50M"), false, true, false, false, false, Q);
        t.add_mate_pair(Some(b), 1, 1, 1, 100, false, false, true, true, Some("50M"), Some("50M"), false, true, false, false, false, Q);
        t.run();
    }
}

#[test]
fn secondary_and_non_primary() {
    let mut t = Tester::new();
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.add_fragment(None, 1, 200, false, None, Q, true, false, false);
    t.run();
    // testMappedFragmentAndMappedPairFirstOfPairNonPrimary
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_mapped_fragment(1, 1, false, Q);
    t.add_mate_pair(None, 1, 1, 200, 0, false, true, false, false, Some("54M22S"), None, false, false, true, true, false, Q);
    t.run();
}

#[test]
fn soft_clipping() {
    // testTwoMappedPairsMatesSoftClipped: no duplicates
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_pair_cigars(1, 10022, 10051, false, false, "76M", "8S68M", false, true, Q);
    t.add_pair_cigars(1, 10022, 10063, false, false, "76M", "5S71M", false, true, Q);
    t.run();
    // testTwoMappedPairsWithSoftClipping: 5' ends 2 and 118 for both... yet no duplicates (different
    // second ends after clipping)
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_pair_cigars(1, 2, 46, false, false, "6S42M28S", "3S73M", false, true, Q);
    t.add_pair_cigars(1, 2, 51, true, true, "6S42M28S", "8S68M", false, true, Q);
    t.run();
}

#[test]
fn unmapped_mates() {
    // testMatePairSecondUnmapped
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_mate_pair(None, 1, 1, 10049, 10049, false, true, false, false, Some("11M2I63M"), None, false, false, false, false, false, Q);
    t.run();
    // testMappedFragmentAndMatePairSecondUnmapped
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_mate_pair(None, 1, 1, 10049, 10049, false, true, false, false, Some("11M2I63M"), None, false, false, false, false, false, Q);
    t.read_length = 36;
    t.add_mapped_fragment(1, 10049, true, Q);
    t.run();
    // testMappedPairAndMatePairSecondUnmapped (MarkDuplicates does not mark unmapped mates)
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_mate_pair(None, 1, 1, 10040, 10040, false, true, true, false, Some("76M"), None, false, false, false, false, false, Q);
    t.add_pair_cigars(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, Q);
    t.run();
    // testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped
    let mut t = Tester::new();
    t.read_length = 101;
    t.add_mate_pair(None, 1, 1, 484071, 484071, false, true, false, false, Some("66S35M"), None, true, false, false, false, false, Q);
    t.add_pair_cigars(1, 484105, 484075, false, false, "35M66S", "30S71M", false, true, Q);
    t.run();
}

#[test]
fn orientations() {
    // testTwoMappedPairsWithOppositeOrientations and NumberTwo: -/+ and +/- of the same fragment
    for second_first in [false, true] {
        let mut t = Tester::new();
        t.read_length = 76;
        if second_first {
            t.add_pair_cigars(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, Q);
            t.add_pair_cigars(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, Q);
        } else {
            t.add_pair_cigars(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, Q);
            t.add_pair_cigars(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, Q);
        }
        t.run();
    }
    // testThreeMappedPairsWithMatchingSecondMate
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_pair_cigars(1, 10181, 10058, false, false, "41S35M", "47M29S", true, false, Q);
    t.add_pair_cigars(1, 10181, 10058, true, true, "37S39M", "44M32S", true, false, Q);
    t.add_pair_cigars(1, 10180, 10058, false, false, "36S40M", "50M26S", true, false, Q);
    t.run();
}

#[test]
fn same_positions() {
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_pair_cigars(1, 4914, 4914, false, false, "37M39S", "73M3S", false, false, Q);
    t.run();
    // testTwoMappedPairWithSamePositionDifferentStrands (+/- and -/+ at one position: FR both)
    let mut t = Tester::new();
    t.read_length = 50;
    t.add_pair_cigars(0, 5604914, 5604914, false, false, "50M", "50M", true, false, Q);
    t.add_pair_cigars(0, 5604914, 5604914, true, true, "50M", "50M", false, true, Q);
    t.run();
    let mut t = Tester::new();
    t.read_length = 50;
    t.add_pair_cigars(0, 5604914, 5604915, false, false, "50M", "50M", true, false, Q);
    t.add_pair_cigars(0, 5604915, 5604914, true, true, "50M", "50M", false, true, Q);
    t.run();
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_pair_cigars(0, 5604914, 5605914, false, false, "37M39S", "73M3S", false, false, Q);
    t.add_pair_cigars(0, 5604914, 5604914, false, false, "37M39S", "73M3S", false, false, Q);
    t.run();
}

#[test]
fn groups_on_different_chromosomes() {
    let mut t = Tester::new();
    for tid in [0, 1] {
        t.add_mapped_pair(tid, 1, 100, false, false, Q);
        t.add_mapped_pair(tid, 1, 100, true, true, Q);
        t.add_mapped_pair(tid, 1, 100, true, true, Q);
    }
    t.run();
    let mut t = Tester::new();
    t.add_mapped_fragment(0, 1, false, Q);
    t.add_mapped_fragment(0, 1, true, Q);
    t.add_mapped_fragment(1, 1, false, Q);
    t.add_mapped_fragment(1, 1, true, Q);
    t.run();
}

#[test]
fn second_end_before_first() {
    let mut t = Tester::new();
    t.read_length = 68;
    t.add_pair_cigars(0, 108855339, 108855323, false, false, "33S35M", "17S51M", true, true, Q);
    t.run();
}

#[test]
fn different_chromosomes_in_opposite_order() {
    let mut t = Tester::new();
    t.read_length = 101;
    t.expected_optical = Some(1);
    t.add_mate_pair(Some("RUNID:6:101:17642:6835"), 0, 1, 123989, 18281, false, false, true, true, Some("37S64M"), Some("52M49S"), false, false, false, false, false, Q);
    t.add_mate_pair(Some("RUNID:6:101:17616:6888"), 1, 0, 18281, 123989, false, false, false, false, Some("52M49S"), Some("37S64M"), false, false, false, false, false, Q);
    t.run();
}

#[test]
fn bulk_fragments_and_pairs() {
    // testBulkFragmentsNoDuplicates / WithDuplicates
    let mut t = Tester::new();
    t.read_length = 100;
    for i in 1..=1000 {
        t.add_mapped_fragment(0, i, false, Q);
    }
    t.run();
    let mut t = Tester::new();
    t.read_length = 100;
    for i in 1..=1000 {
        t.add_mapped_fragment(0, i, false, Q);
        t.add_mapped_fragment(0, i, true, Q);
    }
    t.run();
    let mut t = Tester::new();
    t.read_length = 100;
    for i in 1..=1000 {
        t.add_mapped_pair(0, i, i + 1, false, false, Q);
        t.add_mapped_pair(0, i, i + 1, true, true, Q);
    }
    t.run();
}

/// testOpticalDuplicateBehaviorNullRegex: no optical duplicates without a regex
#[test]
fn optical_null_regex() {
    let mut t = Tester::new();
    t.read_length = 68;
    t.expected_optical = Some(0);
    t.add_mate_pair(Some("RUNID:1:1:15993:13361"), 2, 2, 41212324, 41212310, false, false, false, false, Some("33S35M"), Some("19S49M"), true, true, false, false, false, Q);
    t.add_mate_pair(Some("RUNID:1:1:16020:13352"), 2, 2, 41212324, 41212319, false, false, true, true, Some("33S35M"), Some("28S40M"), true, true, false, false, false, Q);
    t.args.extend(["--READ_NAME_REGEX".to_string(), "null".to_string()]);
    t.run();
}

/// MarkDuplicatesTest.testTwoMappedPairsWithSupplementaryReads: supplementary records are not marked
#[test]
fn supplementary_records() {
    let mut t = Tester::new();
    t.read_length = 76;
    t.add_mapped_pair(1, 1, 100, false, false, Q);
    t.add_mapped_pair(1, 1, 100, true, true, Q);
    t.add_fragment(Some("READ0"), 1, 300, false, Some("40M36S"), Q, false, true, false);
    t.add_fragment(Some("READ1"), 1, 300, false, Some("40M36S"), Q, false, true, false);
    t.run();
}

/// Base qualities decide which pair is kept: the best score stays, whatever its position
#[test]
fn scoring_picks_the_best() {
    let mut t = Tester::new();
    t.args = vec!["--DUPLICATE_SCORING_STRATEGY".into(), "SUM_OF_BASE_QUALITIES".into()];
    t.add_mapped_pair(1, 1, 100, true, true, 20);
    t.add_mapped_pair(1, 1, 100, false, false, 30);
    t.add_mapped_pair(1, 1, 100, true, true, 25);
    t.run();
}
