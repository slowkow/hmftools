//! Ports of the tests of the Java duplicate/UmiDuplicatesTest and duplicate/UmiGroupJitterTest that drive a
//! PartitionReader. The Java tests count the reads passed to a TestBamWriter; here the reader writes unsorted
//! to a temporary BAM, which is read back for the counts.

use super::*;
use crate::bqr::BaseQualRecalibration;
use crate::umi::{DEFAULT_DUPLEX_UMI_DELIM, SINGLE_END_JITTER_COLLAPSE_DISTANCE};
use crate::umi_test::*;
use rust_htslib::bam::Read as _;
use std::sync::atomic::{AtomicUsize, Ordering};

const SUPP_POS_STRAND: char = '+';

/// What a PartitionReader borrows: the test ReduxConfig(refGenome, umiEnabled = true, duplexUmi, ...), the
/// sequence dictionary and a disabled unmapper
struct TestEnv {
    config: Config,
    names: Vec<String>,
    expected_valid_tid: Vec<bool>,
    unmapper: ReadUnmapper,
    bqr: BaseQualRecalibration,
    header: bam::HeaderView,
    out_path: String,
}

impl TestEnv {
    fn new(duplex_umi: bool) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!("redux_rs_umi_tests_{}_{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("ref.fa");
        let bases = generate_random_bases(TEST_CHR_LENGTH);
        std::fs::write(&fasta, format!(">{CHR_1}\n{bases}\n>{CHR_2}\n{bases}\n")).unwrap();
        let fasta = fasta.to_string_lossy().to_string();
        let mut args: Vec<String> =
            ["-sample", "TEST", "-input_bam", "/redux_rs_tests/none.bam", "-ref_genome", &fasta, "-bqr_disable", "-umi_enabled"].iter().map(|s| s.to_string()).collect();
        if duplex_umi {
            args.push("-umi_duplex".to_string());
        }
        let config = Config::from_args(&args).unwrap();
        let genome = crate::refgenome::RefGenome::open(&fasta).unwrap();
        let bqr = BaseQualRecalibration::new(&config, &genome).unwrap();
        let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CHR_1}\tLN:{TEST_CHR_LENGTH}\n@SQ\tSN:{CHR_2}\tLN:{TEST_CHR_LENGTH}\n");
        let header = bam::HeaderView::from_bytes(text.as_bytes());
        let unmapper = ReadUnmapper::disabled(&header);
        let out_path = dir.join("out.bam").to_string_lossy().to_string();
        TestEnv { config, names: names(), expected_valid_tid: vec![true; 2], unmapper, bqr, header, out_path }
    }

    fn new_out(&self) -> RawBamOut {
        let mut out = RawBamOut::create(&self.out_path, 1).unwrap();
        out.write_header(&self.header).unwrap();
        out
    }

    fn partition_reader(&self) -> TestReader<'_> {
        let capture = ReadCapture { jitter: None, bqr: self.bqr.region_reader() };
        let mut reader =
            PartitionReader::new(&self.config, &self.names, &self.expected_valid_tid, &self.unmapper, Some(self.new_out()), None, capture).unwrap();
        // records are written as they are produced
        reader.sorted = None;
        TestReader { reader, env: self, written: Vec::new() }
    }
}

struct TestReader<'a> {
    reader: PartitionReader<'a>,
    env: &'a TestEnv,
    written: Vec<bam::Record>,
}

impl TestReader<'_> {
    fn setup_region(&mut self, chromosome: &str, start: i32, end: i32) {
        self.reader.setup_region(&ChrBaseRegion::new(chromosome, start, end));
    }
    fn process_read(&mut self, read: &bam::Record) {
        let mut read = read.clone();
        self.reader.process_sam_record(&mut read).unwrap();
    }
    /// PartitionReader.flushReadPositions
    fn flush_read_positions(&mut self) {
        let evicted = self.reader.read_cache.evict_all();
        self.reader.process_read_groups(evicted).unwrap();
    }
    fn post_process_region(&mut self) {
        self.reader.post_process_region().unwrap();
    }
    /// every record written so far: the output so far is closed and read back, and a new one started
    fn written_records(&mut self) -> Vec<bam::Record> {
        self.reader.out.take().unwrap().close().unwrap();
        let mut reader = bam::Reader::from_path(&self.env.out_path).unwrap();
        self.written.extend(reader.records().map(|r| r.unwrap()));
        self.reader.out = Some(self.env.new_out());
        self.written.clone()
    }
    /// (BamWriter.nonConsensusWriteCount, BamWriter.consensusWriteCount)
    fn write_counts(&mut self) -> (usize, usize) {
        let records = self.written_records();
        let consensus = records.iter().filter(|r| has_aux(r, b"CR")).count();
        (records.len() - consensus, consensus)
    }
    fn stats(&self) -> &Statistics {
        &self.reader.stats
    }
}

/// a first-in-pair forward read with a reverse mate
fn read_at(read_id: &str, read_pos: i32, mate_pos: i32, supp: Option<&str>) -> bam::Record {
    create_sam_record(read_id, CHR_1, read_pos, &test_read_bases(), TEST_READ_CIGAR, CHR_1, mate_pos, false, false, supp, true, Some(TEST_READ_CIGAR))
}

/// UmiDuplicatesTest.testUmiGroup
#[test]
fn test_umi_group() {
    // 2 primaries in a UMI group, followed by their mates in the same partition and then their supps in a different partition
    let env = TestEnv::new(false);
    let mut reader = env.partition_reader();
    let bases = test_read_bases();
    let umi_id = "TATTAT";
    let (read_pos, mate_pos, supp_pos) = (100, 200, 1800);

    reader.setup_region(CHR_1, 1, 1000);
    let supp = supp_data(CHR_1, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
    let read1 = read_at(&next_read_id(umi_id), read_pos, mate_pos, Some(&supp));
    let read2 = read_at(&next_read_id(umi_id), read_pos, mate_pos, Some(&supp));

    reader.process_read(&read1);
    reader.process_read(&read2);
    reader.flush_read_positions();
    assert_eq!((2, 1), reader.write_counts());

    let mate = |read: &bam::Record| {
        let name = String::from_utf8_lossy(read.qname()).into_owned();
        let mut m = create_sam_record(&name, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false, None, false, Some(TEST_READ_CIGAR));
        set_second_in_pair(&mut m);
        m
    };
    reader.process_read(&mate(&read1));
    assert_eq!((2, 1), reader.write_counts());
    reader.process_read(&mate(&read2));
    reader.flush_read_positions();
    assert_eq!((4, 2), reader.write_counts());

    let supp_of = |read: &bam::Record| {
        let name = String::from_utf8_lossy(read.qname()).into_owned();
        let primary = supp_data(CHR_1, read_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
        create_sam_record(&name, CHR_1, supp_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, false, true, Some(&primary), true, Some(TEST_READ_CIGAR))
    };
    reader.setup_region(CHR_1, 1001, 2000);
    reader.process_read(&supp_of(&read1));
    assert_eq!((4, 2), reader.write_counts());
    reader.process_read(&supp_of(&read2));
    reader.post_process_region();
    assert_eq!((6, 3), reader.write_counts());

    // every read of a group carries the group's UMI, the consensus read has it in its name
    for r in reader.written_records() {
        if has_aux(&r, b"CR") {
            assert!(r.qname().ends_with(b":CNS_TATTAT"));
            assert_eq!(Some(&b"SINGLE"[..]), aux_str(&r, b"UT"));
        } else {
            assert_eq!(Some(&b"TATTAT"[..]), aux_str(&r, b"UI"));
            assert!(has(&r, FLAG_DUPLICATE));
        }
    }
}

/// UmiDuplicatesTest.testSingleUmiGroup
#[test]
fn test_single_umi_group() {
    // 2 UMI groups and a single fragment with the forward orientation, 1 fragment and UMI group on the reverse
    // are matched off with each other to form 3 final groups in total
    let env = TestEnv::new(false);
    let mut reader = env.partition_reader();
    let (read_pos, mate_pos) = (100, 200);
    let (umi1, umi2, umi3, umi4, umi5) = ("AAAAA", "TTTTT", "CCCCC", "GGGGG", "ACGTC");

    reader.setup_region(CHR_1, 1, 1000);
    let mut reads = Vec::new();
    for umi in [umi1, umi1, umi2, umi2, umi2, umi3] {
        let mut r = read_at(&next_read_id(umi), read_pos, mate_pos, None);
        set_second_in_pair(&mut r);
        reads.push(r);
    }
    // and the reverse orientation reads
    for umi in [umi4, umi5, umi5] {
        reads.push(read_at(&next_read_id(umi), read_pos, mate_pos, None));
    }
    for r in &reads {
        reader.process_read(r);
    }
    reader.post_process_region();

    assert_eq!((9, 3), reader.write_counts());
    assert_eq!(3, reader.stats().umi.umi_groups);
    assert_eq!(1, reader.stats().umi.position_fragments.len());
    let counts = &reader.stats().umi.position_fragments[0];
    assert_eq!(2, counts.unique_coord_count);
    assert_eq!(5, counts.unique_fragment_count);
    assert_eq!(1, counts.frequency);
}

/// UmiDuplicatesTest.testSingleUmiGroup2
#[test]
fn test_single_umi_group2() {
    // 2 single fragments with difference coordinates
    let env = TestEnv::new(false);
    let mut reader = env.partition_reader();
    let read1 = read_at(&next_read_id("AAAAA"), 100, 200, None);
    let read2 = read_at(&next_read_id("TTTTT"), 100, 300, None);

    reader.setup_region(CHR_1, 1, 1000);
    reader.process_read(&read1);
    reader.process_read(&read2);
    reader.post_process_region();

    assert_eq!((2, 0), reader.write_counts());
    let counts = &reader.stats().umi.position_fragments[0];
    assert_eq!(2, counts.unique_coord_count);
    assert_eq!(2, counts.unique_fragment_count);
    assert_eq!(1, counts.frequency);
}

/// UmiDuplicatesTest.testDuplexUmiGroup
#[test]
fn test_duplex_umi_group() {
    // 2 fragments have opposite fragments but are linked by their duplex UMIs
    let env = TestEnv::new(true);
    let mut reader = env.partition_reader();
    let d = DEFAULT_DUPLEX_UMI_DELIM;
    let (umi1_part1, umi1_part2) = ("TATTAT", "GCGGCG");
    let (umi2_part1, umi2_part2) = ("TTAATT", "GGCCGG");
    let umi_id1 = format!("{umi1_part1}{d}{umi1_part2}");
    let umi_id1_reversed = format!("{umi1_part2}{d}{}G", &umi1_part1[..5]); // 1 base diff is allowed
    let umi_id2 = format!("{umi2_part1}{d}{umi2_part2}");
    let umi_id2_reversed = format!("{}A{d}{umi2_part1}", &umi2_part2[..5]);

    let read_pos = 100;
    reader.setup_region(CHR_1, 1, 1000);

    // 2 single fragments
    let read1 = read_at(&next_read_id(&umi_id1), read_pos, 200, None);
    let mut read2 = read_at(&next_read_id(&umi_id1_reversed), read_pos, 200, None);
    set_second_in_pair(&mut read2);

    // a pair and a reversed single
    let read3 = read_at(&next_read_id(&umi_id2), read_pos, 300, None);
    let mut read4 = read_at(&next_read_id(&umi_id2_reversed), read_pos, 300, None);
    set_second_in_pair(&mut read4);
    let mut read5 = read_at(&next_read_id(&umi_id2_reversed), read_pos, 300, None);
    set_second_in_pair(&mut read5);

    for r in [&read1, &read2, &read3, &read4, &read5] {
        reader.process_read(r);
    }
    reader.post_process_region();

    assert_eq!((5, 2), reader.write_counts());
    assert_eq!(2, reader.stats().umi.umi_groups);
    assert_eq!(1, reader.stats().umi.position_fragments.len());
    let counts = &reader.stats().umi.position_fragments[0];
    assert_eq!(2, counts.unique_coord_count);
    assert_eq!(2, counts.unique_fragment_count);
    assert_eq!(1, counts.frequency);

    // (rounded size, first size seen, count, dual strand count)
    let frequencies = &reader.stats().frequencies;
    assert_eq!(2, frequencies.len());
    assert_eq!(1, frequencies.iter().find(|f| f.0 == 2).unwrap().3);
    assert_eq!(1, frequencies.iter().find(|f| f.0 == 3).unwrap().3);

    for r in reader.written_records() {
        if has_aux(&r, b"CR") {
            assert_eq!(Some(&b"DUAL"[..]), aux_str(&r, b"UT"));
        }
    }
}

/// UmiDuplicatesTest.testAlternativeFirstInPairGroups
#[test]
fn test_alternative_first_in_pair_groups() {
    // 4 primary reads with alternating first in pair designation for the same lower / upper reads
    let env = TestEnv::new(true);
    let mut reader = env.partition_reader();
    let bases = test_read_bases();
    let (read_pos, mate_pos, supp_pos, supp_pos2) = (100, 1500, 2200, 2400);
    let d = DEFAULT_DUPLEX_UMI_DELIM;
    let umi_id1 = format!("TATTAT{d}GCGGCG");
    let umi_id1_reversed = format!("GCGGCG{d}TATTAT");

    let supp_data1 = supp_data(CHR_1, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
    let supp_data2 = supp_data(CHR_1, supp_pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);

    let read1 = read_at(&next_read_id(&umi_id1), read_pos, mate_pos, Some(&supp_data1));
    let mut read2 = read_at(&next_read_id(&umi_id1_reversed), read_pos, mate_pos, None);
    set_second_in_pair(&mut read2);
    let read3 = read_at(&next_read_id(&umi_id1), read_pos, mate_pos, Some(&supp_data1));
    let mut read4 = read_at(&next_read_id(&umi_id1_reversed), read_pos, mate_pos, Some(&supp_data1));
    set_second_in_pair(&mut read4);
    let name = |r: &bam::Record| String::from_utf8_lossy(r.qname()).into_owned();

    reader.setup_region(CHR_1, 1, 1000);
    for r in [&read1, &read2, &read3, &read4] {
        reader.process_read(r);
    }
    reader.flush_read_positions();
    reader.post_process_region();

    let mate = |id: &str, supp: Option<&str>, second: bool| {
        let mut m = create_sam_record(id, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false, supp, false, Some(TEST_READ_CIGAR));
        if second {
            set_second_in_pair(&mut m);
        }
        m
    };
    let mate1 = mate(&name(&read1), Some(&supp_data2), true);
    let mate2 = mate(&name(&read2), Some(&supp_data2), false);
    let mate3 = mate(&name(&read1), None, true);
    let mate4 = mate(&name(&read2), None, false);

    reader.setup_region(CHR_1, 1001, 2000);
    for r in [&mate1, &mate2, &mate3, &mate4] {
        reader.process_read(r);
    }
    reader.flush_read_positions();
    reader.post_process_region();

    // finally the expected supplementaries
    reader.setup_region(CHR_1, 2001, 3000);
    let supp_data_read1 = supp_data(CHR_1, read_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
    let supp_data_read2 = supp_data(CHR_1, mate_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);
    let supp = |id: &str, pos: i32, other_pos: i32, reversed: bool, primary: &str, mate_reversed: bool, second: bool| {
        let mut s =
            create_sam_record(id, CHR_1, pos, &bases, TEST_READ_CIGAR, CHR_1, other_pos, reversed, true, Some(primary), mate_reversed, Some(TEST_READ_CIGAR));
        if second {
            set_second_in_pair(&mut s);
        }
        s
    };
    let supp1 = supp(&name(&read1), supp_pos, mate_pos, false, &supp_data_read1, true, false);
    // read 2 has no supp
    let supp3 = supp(&name(&read3), supp_pos, mate_pos, false, &supp_data_read1, true, false);
    let supp4 = supp(&name(&read4), supp_pos, mate_pos, false, &supp_data_read1, true, true);
    let mate_supp1 = supp(&name(&read1), supp_pos2, read_pos, true, &supp_data_read2, false, true);
    let mate_supp2 = supp(&name(&read2), supp_pos2, read_pos, true, &supp_data_read2, false, false);

    for r in [&supp1, &supp3, &supp4, &mate_supp1, &mate_supp2] {
        reader.process_read(r);
    }
    reader.post_process_region();

    assert_eq!((13, 4), reader.write_counts());
}

/// the two reads of a jitter test processed in a region of chromosome 1, and the counts written
fn jitter_counts(duplex: bool, read1: &bam::Record, read2: &bam::Record) -> (usize, usize) {
    let env = TestEnv::new(duplex);
    let mut reader = env.partition_reader();
    reader.setup_region(CHR_1, 1, 1000);
    reader.process_read(read1);
    reader.process_read(read2);
    reader.post_process_region();
    reader.write_counts()
}

/// UmiGroupJitterTest.testIlluminaJitterUmiGroupCollapse
#[test]
fn test_illumina_jitter_umi_group_collapse() {
    let umi = format!("AAAAA{DEFAULT_DUPLEX_UMI_DELIM}AAATA");
    let read1 = read_at(&format!("READ_001:{umi}"), 100, 1000, None);
    let read2 = read_at(&format!("READ_002:{umi}"), 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, 1000, None);
    assert_eq!((2, 1), jitter_counts(true, &read1, &read2));
}

/// UmiGroupJitterTest.testIlluminaMateJitterUmiGroupCollapse
#[test]
fn test_illumina_mate_jitter_umi_group_collapse() {
    let umi = "AAAAA";
    let read1 = read_at(&next_read_id(umi), 100, 1000, None);
    let read2 = read_at(&next_read_id(umi), 100, 1000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, None);
    assert_eq!((2, 1), jitter_counts(true, &read1, &read2));
}

/// UmiGroupJitterTest.testIlluminaMateJitterDuplexUmiGroupCollapse
#[test]
fn test_illumina_mate_jitter_duplex_umi_group_collapse() {
    let d = DEFAULT_DUPLEX_UMI_DELIM;
    let umi_id = format!("TATTAT{d}GCGGCG");
    let umi_id_reversed = format!("GCGGCG{d}TATTAT");
    let mut read1 = read_at(&next_read_id(&umi_id), 100, 1000, None);
    let read2 = read_at(&next_read_id(&umi_id_reversed), 100, 1000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, None);
    // as in the Java test, it is the first read that ends up second in pair
    set_first_of_pair(&mut read1, false);
    assert_eq!((2, 1), jitter_counts(true, &read1, &read2));
}

/// UmiGroupJitterTest.testIlluminaJitterUmiGroupNoCollapseDueToNoFixedFragmentEnd
#[test]
fn test_illumina_jitter_umi_group_no_collapse_due_to_no_fixed_fragment_end() {
    let umi = "AAAAA";
    let read1 = read_at(&next_read_id(umi), 100, 1000, None);
    let read2 = read_at(&next_read_id(umi), 101, 1001, None);
    assert_eq!((2, 0), jitter_counts(false, &read1, &read2));
}

/// UmiGroupJitterTest.testIlluminaJitterDuplexUmiGroupNoCollapseDueToUMIMismatch
#[test]
fn test_illumina_jitter_duplex_umi_group_no_collapse_due_to_umi_mismatch() {
    let d = DEFAULT_DUPLEX_UMI_DELIM;
    let umi_id = format!("TATTAT{d}GCGGCG");
    let umi_id_reversed = format!("GCATCG{d}TATTAT");
    let mut read1 = read_at(&next_read_id(&umi_id), 100, 1000, None);
    let read2 = read_at(&next_read_id(&umi_id_reversed), 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, 1000, None);
    set_first_of_pair(&mut read1, false);
    assert_eq!((2, 0), jitter_counts(true, &read1, &read2));
}

/// UmiGroupJitterTest: the two reads of a fragment from its lower and upper positions
fn fragment_records(read_name: &str, lower_pos: i32, upper_pos: i32) -> Vec<bam::Record> {
    let bases = test_read_bases();
    let upper_start = upper_pos - bases.len() as i32 + 1;
    let lower = create_sam_record(read_name, CHR_1, lower_pos, &bases, TEST_READ_CIGAR, CHR_1, upper_start, false, false, None, true, Some(TEST_READ_CIGAR));
    let mut upper = create_sam_record(read_name, CHR_1, upper_start, &bases, TEST_READ_CIGAR, CHR_1, lower_pos, true, false, None, false, Some(TEST_READ_CIGAR));
    flip_first_in_pair(&mut upper);
    vec![lower, upper]
}

/// the fragments' reads processed in order of alignment start in a region of chromosome 1: the counts written and
/// how many read names were not written exactly twice
fn inconsistent_pops_counts(fragments: &[(&str, i32, i32)]) -> ((usize, usize), usize, usize) {
    let env = TestEnv::new(true);
    let mut reader = env.partition_reader();
    reader.setup_region(CHR_1, 1, 2000);
    let mut reads: Vec<bam::Record> = fragments.iter().flat_map(|(name, lower, upper)| fragment_records(name, *lower, *upper)).collect();
    reads.sort_by_key(|r| r.pos());
    for r in &reads {
        reader.process_read(r);
    }
    reader.post_process_region();
    let mut names: Vec<Vec<u8>> = reader.written_records().iter().map(|r| r.qname().to_vec()).collect();
    names.sort();
    let total = names.len();
    let mut not_twice = 0;
    for chunk in names.chunk_by(|a, b| a == b) {
        if chunk.len() != 2 {
            not_twice += 1;
        }
    }
    (reader.write_counts(), total, not_twice)
}

const UMI_ID1: &str = "CATGATC";
const UMI_ID2: &str = "TCCTATC";
const MAX_SOFT_CLIP: i32 = DEFAULT_MAX_SOFT_CLIP;

/// UmiGroupJitterTest.testIlluminaJitterUmiGroupCollapseReadCacheInconsistentPops, ignored in Java
#[test]
#[ignore = "ignored in the Java tests"]
fn test_illumina_jitter_umi_group_collapse_read_cache_inconsistent_pops() {
    let (id1, id2, id3) = (format!("READ_001:{UMI_ID1}"), format!("READ_002:{UMI_ID1}"), format!("READ_003:{UMI_ID2}"));
    let result = inconsistent_pops_counts(&[(&id1, 100, 700), (&id2, 101, 700), (&id3, 100 + MAX_SOFT_CLIP, 1000)]);
    assert_eq!(((6, 2), 8, 0), result);
}

/// UmiGroupJitterTest.testIlluminaJitterUmiGroupCollapseReadCacheInconsistentPopsNotFixedByIncreasingMaxSoftClip
#[test]
fn test_illumina_jitter_umi_group_collapse_read_cache_inconsistent_pops_not_fixed_by_increasing_max_soft_clip() {
    let (id1, id2, id3) = (format!("READ_001:{UMI_ID1}"), format!("READ_002:{UMI_ID1}"), format!("READ_003:{UMI_ID2}"));
    let result =
        inconsistent_pops_counts(&[(&id1, 100, 700), (&id2, 101, 700), (&id3, 100 + MAX_SOFT_CLIP + SINGLE_END_JITTER_COLLAPSE_DISTANCE, 1000)]);
    assert_eq!(((6, 2), 8, 0), result);
}

/// UmiGroupJitterTest.testIlluminaJitterUmiGroupCollapseReadCacheInconsistentPopsNotFixedByIncreasingMaxSoftClipAsAnInitialFilter,
/// ignored in Java
#[test]
#[ignore = "ignored in the Java tests"]
fn test_illumina_jitter_umi_group_collapse_read_cache_inconsistent_pops_not_fixed_as_an_initial_filter() {
    let ids: Vec<String> = (1..=5).map(|i| format!("READ_00{i}:{}", if i < 5 { UMI_ID1 } else { UMI_ID2 })).collect();
    let result = inconsistent_pops_counts(&[
        (&ids[0], 100, 700),
        (&ids[1], 110, 700),
        (&ids[2], 111, 700),
        (&ids[3], 111, 700),
        (&ids[4], 100 + MAX_SOFT_CLIP + SINGLE_END_JITTER_COLLAPSE_DISTANCE, 1000),
    ]);
    assert_eq!(((10, 2), 12, 0), result);
}
