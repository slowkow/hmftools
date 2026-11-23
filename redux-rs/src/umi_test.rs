//! Ports of the Java tests duplicate/UmiGroupsTest, and of the tests of duplicate/UmiDuplicatesTest and
//! duplicate/UmiGroupJitterTest that do not need a PartitionReader (those are in umi_partition_test.rs).
//! The helpers build synthetic reads the way the Java SamRecordTestUtils and TestUtils do.

use crate::consensus::{form_consensus_read_id, select_template_read};
use crate::dup::{DuplicateGroup, FragmentCoords, ReadCache, ReadInfo};
use crate::rec;
use crate::umi::*;
use rust_htslib::bam::{self, record::{Cigar, CigarString}};
use std::sync::atomic::{AtomicI32, Ordering};

pub const CHR_1: &str = "1";
pub const CHR_2: &str = "2";
pub const NO_CHROMOSOME_NAME: &str = "*";
pub const NO_CIGAR: &str = "*";
pub const NO_POSITION: i32 = 0;
pub const TEST_READ_CIGAR: &str = "100M";
pub const FIXED_READ_ID: &str = "123:ABC:1:4455:";
/// length of the two test chromosomes
pub const TEST_CHR_LENGTH: usize = 5000;

/// the test sequence dictionary: chromosomes 1 and 2
pub fn names() -> Vec<String> {
    vec![CHR_1.to_string(), CHR_2.to_string()]
}

/// MockRefGenome.generateRandomBases: cycles through A, C, G, T
pub fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| b"ACGT"[i % 4] as char).collect()
}

pub fn test_read_bases() -> String {
    generate_random_bases(100)
}

/// TestUtils.READ_ID_GEN with UmiDuplicatesTest.nextReadId: a new read id that ends in the UMI
pub fn next_read_id(umi_id: &str) -> String {
    static NEXT: AtomicI32 = AtomicI32::new(1);
    format!("READ_{:03}:{}", NEXT.fetch_add(1, Ordering::Relaxed), umi_id)
}

fn chromosome_index(chr: &str) -> i32 {
    match chr {
        CHR_1 => 0,
        CHR_2 => 1,
        _ => -1,
    }
}

pub fn parse_cigar(cigar: &str) -> CigarString {
    let mut out = Vec::new();
    if cigar == NO_CIGAR {
        return CigarString(out);
    }
    let mut len: u32 = 0;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            len = len * 10 + d;
            continue;
        }
        out.push(match c {
            'M' => Cigar::Match(len),
            'I' => Cigar::Ins(len),
            'D' => Cigar::Del(len),
            'S' => Cigar::SoftClip(len),
            other => panic!("cigar operator {other} is not used in these tests"),
        });
        len = 0;
    }
    CigarString(out)
}

pub fn set_string_attribute(record: &mut bam::Record, tag: &[u8; 2], value: &str) {
    rec::set_attribute(record, tag, Some(&rec::encode_z(value.as_bytes())));
}

/// SupplementaryReadData.asSamTag for chromosome, position, strand, cigar, map quality (NM is 0)
pub fn supp_data(chromosome: &str, position: i32, strand: char, cigar: &str, map_quality: i32) -> String {
    format!("{chromosome},{position},{strand},{cigar},{map_quality},0")
}

/// SamRecordTestUtils.createSamRecord (12 arguments): a paired, first-in-pair read as SAMRecordSetBuilder.addFrag
/// makes it (read group attribute, map quality 255), with base qualities of 37
#[allow(clippy::too_many_arguments)]
pub fn create_sam_record(
    read_id: &str, chr: &str, read_start: i32, read_bases: &str, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool,
    is_supplementary: bool, supp_alignment: Option<&str>, mate_reversed: bool, mate_cigar: Option<&str>,
) -> bam::Record {
    let chr_index = chromosome_index(chr);
    let mut record = bam::Record::new();
    let quals = vec![37u8; read_bases.len()];
    record.set(read_id.as_bytes(), Some(&parse_cigar(cigar)), read_bases.as_bytes(), &quals);
    record.set_tid(chr_index);
    // the builder sets the alignment start only with a contig
    record.set_pos(if chr_index >= 0 { read_start as i64 - 1 } else { -1 });
    record.set_mapq(255);
    record.set_insert_size(0);
    let mut flags = rec::FLAG_PAIRED | rec::FLAG_FIRST;
    if is_reversed {
        flags |= rec::FLAG_REVERSE;
    }
    if chr_index < 0 {
        flags |= rec::FLAG_UNMAPPED;
    }
    if !mate_chr.is_empty() && mate_chr != NO_CHROMOSOME_NAME {
        record.set_mtid(chromosome_index(mate_chr));
        record.set_mpos(mate_start as i64 - 1);
        flags |= rec::FLAG_PROPER_PAIR;
    } else {
        flags |= rec::FLAG_MATE_UNMAPPED;
        record.set_mtid(-1);
        record.set_mpos(read_start as i64 - 1);
    }
    if is_supplementary {
        flags |= rec::FLAG_SUPPLEMENTARY;
    }
    if mate_reversed {
        flags |= rec::FLAG_MATE_REVERSE;
    }
    record.set_flags(flags);
    set_string_attribute(&mut record, b"RG", "1");
    if let Some(supp) = supp_alignment {
        set_string_attribute(&mut record, b"SA", supp);
    }
    if chr == mate_chr {
        record.set_insert_size((read_start - mate_start).abs() as i64);
    }
    if let Some(mc) = mate_cigar {
        set_string_attribute(&mut record, b"MC", mc);
    }
    rec::fix_bin(&mut record);
    record
}

/// setFirstOfPairFlag(first) and setSecondOfPairFlag(!first)
pub fn set_first_of_pair(record: &mut bam::Record, first: bool) {
    rec::set_flag(record, rec::FLAG_FIRST, first);
    rec::set_flag(record, rec::FLAG_SECOND, !first);
}

/// TestUtils.setSecondInPair
pub fn set_second_in_pair(record: &mut bam::Record) {
    set_first_of_pair(record, false);
}

/// SamRecordTestUtils.flipFirstInPair
pub fn flip_first_in_pair(record: &mut bam::Record) {
    let first = rec::has(record, rec::FLAG_FIRST);
    set_first_of_pair(record, !first);
}

/// SamRecordTestUtils.cloneSamRecord
pub fn clone_sam_record(record: &bam::Record, new_read_id: &str) -> bam::Record {
    let mut r = record.clone();
    r.set_qname(new_read_id.as_bytes());
    r
}

/// The Java tests find a read in a group by object identity; here each read carries a label attribute.
fn labelled(mut record: bam::Record, label: &str) -> bam::Record {
    set_string_attribute(&mut record, b"ZT", label);
    record
}

fn label(record: &bam::Record) -> String {
    String::from_utf8_lossy(rec::aux_str(record, b"ZT").unwrap()).into_owned()
}

fn contains(reads: &[bam::Record], wanted: &str) -> bool {
    reads.iter().any(|r| label(r) == wanted)
}

fn all_reads(group: &DuplicateGroup) -> Vec<bam::Record> {
    group.reads.iter().chain(group.non_consensus_reads.iter()).cloned().collect()
}

fn read_name(record: &bam::Record) -> String {
    String::from_utf8_lossy(record.qname()).into_owned()
}

fn umi_config() -> UmiConfig {
    UmiConfig::new(true, false, "", true, false)
}

fn exceeds(first: &str, second: &str) -> bool {
    exceeds_umi_id_diff(first.as_bytes(), second.as_bytes(), DEFAULT_MAX_UMI_BASE_DIFF)
}

fn create_simple_read(read_id: &str, label: &str) -> bam::Record {
    labelled(
        create_sam_record(read_id, CHR_1, 100, &test_read_bases(), TEST_READ_CIGAR, CHR_1, 200, false, false, None, true, Some(TEST_READ_CIGAR)),
        label,
    )
}

fn fixed(umi: &str) -> String {
    format!("{FIXED_READ_ID}{umi}")
}

/// TestUtils.createFragmentCoords
fn create_fragment_coords(record: &bam::Record) -> FragmentCoords {
    FragmentCoords::from_read(record, &names(), false)
}

fn group_with<'a>(groups: &'a [DuplicateGroup], wanted: &str) -> &'a DuplicateGroup {
    groups.iter().find(|g| contains(&g.reads, wanted)).unwrap()
}

/// UmiGroupsTest.testUmiUtils
#[test]
fn test_umi_utils() {
    let umi_id = "TATCGC";
    let read_id = fixed(umi_id);
    assert_eq!(umi_id.as_bytes(), extract_umi_id_from_read_id(read_id.as_bytes()));
    assert_eq!(umi_id.as_bytes(), extract_umi_id_fixed(read_id.as_bytes(), umi_id.len()));
    assert_eq!(umi_id.as_bytes(), umi_config().extract_umi_id(read_id.as_bytes()));

    let read_id1 = fixed("TATCGC");
    let read_id2 = fixed("TATCGG");
    let read_id3 = fixed("TATGGG");
    let read_id4 = fixed("AAACGG");
    let read_id5 = fixed("AAAGGG");
    let read_id6 = fixed("AAAGGGT");

    assert!(!exceeds(&read_id1, &read_id2));
    assert!(!exceeds(&read_id2, &read_id3));
    assert!(exceeds(&read_id1, &read_id3));
    assert!(exceeds(&read_id2, &read_id4));
    assert!(!exceeds(&read_id4, &read_id5));
    assert!(exceeds(&read_id5, &read_id6));
}

/// UmiGroupsTest.testDuplexUmis
#[test]
fn test_duplex_umis() {
    let config = UmiConfig::new(true, true, DEFAULT_DUPLEX_UMI_DELIM, true, false);
    let matches = |a: &str, b: &str| has_duplex_umi_match(a, b, &config.duplex_delim, config.permitted_base_diff).unwrap();

    let umi_id1 = "TATCGC_AAGTCG";
    assert!(!matches(umi_id1, umi_id1));

    assert!(matches(umi_id1, "AAGTCC_TATCGG"));

    // too many different bases
    assert!(!matches(umi_id1, "AAGCCC_TATCGC"));

    // a UMI without the delimiter is an error in Java (R25)
    assert!(crate::compat::with_java(true, || has_duplex_umi_match(umi_id1, "AAGTCCTATCGG", "_", 1)).is_err());
}

#[test]
fn duplex_umi_without_delimiter_is_no_match() {
    // R25 corrected
    assert!(!has_duplex_umi_match("TATCGC_AAGTCG", "AAGTCCTATCGG", "_", 1).unwrap());
    assert!(!has_duplex_umi_match("AAGTCCTATCGG", "TATCGC_AAGTCG", "_", 1).unwrap());
}

#[test]
fn duplex_umi_delimiter_of_two_characters() {
    // R10: Java starts the second half one character after the start of the delimiter
    assert!(has_duplex_umi_match("TATCGC--AAGTCG", "AAGTCG--TATCGC", "--", 0).unwrap());
    assert!(!crate::compat::with_java(true, || has_duplex_umi_match("TATCGC--AAGTCG", "AAGTCG--TATCGC", "--", 0)).unwrap());
}

/// UmiGroupsTest.testUmiGroupAssignment
#[test]
fn test_umi_group_assignment() {
    let config = umi_config();
    let frag1 = create_simple_read(&fixed("TATCGC"), "frag1");
    let frag2 = create_simple_read(&fixed("TATCGC"), "frag2");
    let frag3 = create_simple_read(&fixed("TATCCC"), "frag3");
    let frag4 = create_simple_read(&fixed("AATGGG"), "frag4");
    let frag5 = create_simple_read(&fixed("TATCGC"), "frag5");
    let frag6 = create_simple_read(&fixed("AATGGG"), "frag6");

    let fragment_coords = create_fragment_coords(&frag1);
    let groups = build_umi_groups(&fragment_coords, vec![frag1, frag2, frag3, frag4, frag5, frag6], &config);
    assert_eq!(2, groups.len());

    let group = group_with(&groups, "frag1");
    assert!(contains(&group.reads, "frag2"));
    assert!(contains(&group.reads, "frag3"));
    assert!(contains(&group.reads, "frag5"));

    let group = group_with(&groups, "frag4");
    assert!(contains(&group.reads, "frag6"));

    // check that larger groups aren't merged
    let frag11 = create_simple_read(&fixed("TTTCGT"), "frag11");
    let frag12 = create_simple_read(&fixed("TTCCGT"), "frag12");
    let frag13 = create_simple_read(&fixed("TTACGT"), "frag13");
    let frag14 = create_simple_read(&fixed("TTACAT"), "frag14");
    let frag15 = create_simple_read(&fixed("TAAAAT"), "frag15");
    let frag16 = create_simple_read(&fixed("TTACAG"), "frag16");

    // duplicate some to build bigger UMI groups
    let mut fragments = Vec::new();
    for (frag, count) in [(&frag11, 2), (&frag12, 3), (&frag13, 8), (&frag14, 4), (&frag15, 5), (&frag16, 1)] {
        for _ in 0..count {
            fragments.push(frag.clone());
        }
    }
    let fragment_count = fragments.len();

    let groups = build_umi_groups(&fragment_coords, fragments, &config);
    assert_eq!(2, groups.len());
    assert_eq!(fragment_count, groups.iter().map(|g| g.total_read_count()).sum::<usize>());

    let group = group_with(&groups, "frag11");
    assert_eq!(18, group.total_read_count());
    assert!(contains(&group.reads, "frag12"));
    assert!(contains(&group.reads, "frag13"));
    assert!(contains(&group.reads, "frag14"));
    assert!(contains(&group.reads, "frag16"));

    let group = group_with(&groups, "frag15");
    assert_eq!(5, group.total_read_count());
}

/// UmiGroupsTest.testUmiGroupAssignment2
#[test]
fn test_umi_group_assignment2() {
    // test a final merge of groups with 2 bases difference
    let fragments: Vec<bam::Record> =
        ["TTAAGG", "TTAAGC", "TTAATT", "GGGATT", "GGGAGG", "CCGAGC"].iter().map(|u| create_simple_read(&fixed(u), u)).collect();
    let fragment_coords = create_fragment_coords(&fragments[0]);
    let groups = build_umi_groups(&fragment_coords, fragments, &umi_config());
    assert_eq!(3, groups.len());
}

/// UmiGroupsTest.testUmiGroupAssignment3
#[test]
fn test_umi_group_assignment3() {
    // test a final merge of groups with 4 bases difference if the groups are large enough
    let umi1 = "TTAAGG";
    let umi2 = "TTAGGA"; // 2-base diff from #1
    let umi3 = "CCGGTT"; // unrelated
    let umi4 = "TTCCCC"; // 4-base diff from #1 and #2

    let mut fragments = vec![create_simple_read(&fixed(umi1), "1"), create_simple_read(&fixed(umi2), "2"), create_simple_read(&fixed(umi3), "3")];
    for _ in 0..51 {
        fragments.push(create_simple_read(&fixed(umi4), "4"));
    }
    let fragment_coords = create_fragment_coords(&fragments[0]);
    let groups = build_umi_groups(&fragment_coords, fragments, &umi_config());
    assert_eq!(2, groups.len());
}

/// UmiGroupsTest.testDefinedUmis
#[test]
fn test_defined_umis() {
    let mut config = umi_config();
    let defined_umi1 = "AAAGGG";
    let defined_umi2 = "TTTAAA";
    let defined_umi3 = "CCCAAA";
    config.add_defined_umis(&[defined_umi1, defined_umi2, defined_umi3]);

    let fragments: Vec<bam::Record> =
        [defined_umi1, "AAAGGC", defined_umi1, defined_umi2, "TTTAAC", defined_umi3].iter().map(|u| create_simple_read(&fixed(u), u)).collect();
    let fragment_coords = create_fragment_coords(&fragments[0]);
    let groups = build_umi_groups(&fragment_coords, fragments, &config);
    assert_eq!(3, groups.len());
    // the reads with a base difference take the defined UMI
    for g in &groups {
        assert!([defined_umi1, defined_umi2, defined_umi3].contains(&g.umi.as_deref().unwrap()));
    }
}

/// The order of the defined UMIs is that of a java.util.HashSet: by hash bucket, then insertion
#[test]
fn test_defined_umi_order() {
    let mut config = umi_config();
    // String.hashCode: "Aa" and "BB" collide (2112), "C" is 67; in a table of 16 "C" is in bucket 3, the others in bucket 0
    config.add_defined_umis(&["C", "BB", "Aa"]);
    assert_eq!(Some("BB"), config.match_defined_umi_id(b"BB"));
    // one base from both "BB" and "Aa": the first in iteration order wins
    assert_eq!(Some("BB"), config.match_defined_umi_id(b"Ba"));
    assert_eq!(None, config.match_defined_umi_id(b"CCC"));
}

/// UmiDuplicatesTest.testConsistentReadPairsAfterUmiCollapsing
#[test]
fn test_consistent_read_pairs_after_umi_collapsing() {
    let config = UmiConfig::new(true, true, DEFAULT_DUPLEX_UMI_DELIM, true, false);
    let names = names();
    let bases = test_read_bases();

    let umi_id1 = format!("GTCACTC{DEFAULT_DUPLEX_UMI_DELIM}TAAGATA");
    let umi_id2 = format!("GTCACTC{DEFAULT_DUPLEX_UMI_DELIM}TACGATC");
    let mate_position = 200;
    let id1 = format!("READ_001:{umi_id1}");
    let id2 = format!("READ_002:{umi_id2}");

    let read1 = create_sam_record(&id1, CHR_1, 200, &bases, TEST_READ_CIGAR, CHR_2, mate_position, false, false, None, true, Some(TEST_READ_CIGAR));
    let mate1 = create_sam_record(&id1, CHR_2, 200, &bases, TEST_READ_CIGAR, CHR_1, mate_position, true, false, None, false, Some(TEST_READ_CIGAR));
    let read2 = create_sam_record(&id2, CHR_1, 200, &bases, TEST_READ_CIGAR, CHR_2, mate_position, false, false, None, true, Some(TEST_READ_CIGAR));
    let mate2 = create_sam_record(&id2, CHR_2, 200, &bases, TEST_READ_CIGAR, CHR_1, mate_position, true, false, None, false, Some(TEST_READ_CIGAR));

    let read_coords = FragmentCoords::from_read(&read1, &names, true);
    let mate_coords = FragmentCoords::from_read(&mate1, &names, true);

    let mut read_umi_groups = build_umi_groups(&read_coords, vec![read1, read2], &config);
    let mut mate_umi_groups = build_umi_groups(&mate_coords, vec![mate2, mate1], &config);
    assert_eq!(1, read_umi_groups.len());
    assert_eq!(1, mate_umi_groups.len());

    let consensus_read_name = |group: &mut DuplicateGroup| {
        let template = select_template_read(&mut group.reads, &group.coords);
        form_consensus_read_id(group.reads[template].qname(), group.umi.as_deref())
    };
    assert_eq!(consensus_read_name(&mut read_umi_groups[0]), consensus_read_name(&mut mate_umi_groups[0]));
}

fn form_single_fragments(reads: &[&bam::Record]) -> Vec<ReadInfo> {
    let names = names();
    reads.iter().map(|r| ReadInfo { read: (*r).clone(), coords: FragmentCoords::from_read(r, &names, true) }).collect()
}

/// UmiDuplicatesTest.findMatchingGroup
fn find_matching_group<'a>(groups: &'a [DuplicateGroup], labels: &[&str], check_non_consensus: bool) -> Option<&'a DuplicateGroup> {
    for g in groups {
        if labels.iter().all(|l| contains(&g.reads, l)) {
            return Some(g);
        }
        if check_non_consensus && labels.iter().all(|l| contains(&all_reads(g), l)) {
            return Some(g);
        }
    }
    None
}

/// UmiDuplicatesTest.testPolyGMerging
#[test]
fn test_poly_g_merging() {
    let bases = test_read_bases();
    let umi_id_part1 = "A".repeat(7);
    let umi_id1_part2 = "G".repeat(7);
    let umi_id2_part2 = format!("C{}", "G".repeat(6));

    // test 1: basic case where all mates are unmapped
    let read_name1 = format!("READ_001:{umi_id_part1}{DEFAULT_DUPLEX_UMI_DELIM}{umi_id1_part2}");
    let read_name2 = format!("READ_002:{umi_id_part1}{DEFAULT_DUPLEX_UMI_DELIM}{umi_id2_part2}");
    let read_position = 1000;

    let unmapped_mate_pair = |name: &str| {
        let mut read =
            create_sam_record(name, CHR_1, read_position, &bases, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, None, false, Some(NO_CIGAR));
        rec::set_flag(&mut read, rec::FLAG_MATE_UNMAPPED, true);
        let mut mate =
            create_sam_record(name, NO_CHROMOSOME_NAME, read_position, &bases, NO_CIGAR, CHR_1, read_position, false, false, None, false, Some(TEST_READ_CIGAR));
        flip_first_in_pair(&mut mate);
        (read, mate)
    };
    let (read1, mate1) = unmapped_mate_pair(&read_name1);
    let (read1, mate1) = (labelled(read1, "read1"), labelled(mate1, "mate1"));
    let read2 = labelled(clone_sam_record(&read1, &read_name2), "read2");
    let mate2 = labelled(clone_sam_record(&mate1, &read_name2), "mate2");

    let config = UmiConfig::new(true, true, "_", true, false);

    let (umi_groups, single_fragments) = poly_g_merge(&config, Vec::new(), form_single_fragments(&[&read1, &mate1, &read2, &mate2]));
    assert_eq!(2, umi_groups.len());
    assert!(single_fragments.is_empty());

    let group = find_matching_group(&umi_groups, &["read1", "read2"], false).unwrap();
    assert!(!group.poly_g_unmapped);
    let group = find_matching_group(&umi_groups, &["mate1", "mate2"], true).unwrap();
    assert!(!group.poly_g_unmapped);

    // test 2: UMIs must have poly-G tail
    let name1 = format!("READ_001:{umi_id_part1}{DEFAULT_DUPLEX_UMI_DELIM}{}", "C".repeat(7));
    let name2 = format!("READ_002:{umi_id_part1}{DEFAULT_DUPLEX_UMI_DELIM}{}A", "C".repeat(6));
    let (r1, m1, r2, m2) =
        (clone_sam_record(&read1, &name1), clone_sam_record(&mate1, &name1), clone_sam_record(&read2, &name2), clone_sam_record(&mate2, &name2));
    let (umi_groups, single_fragments) = poly_g_merge(&config, Vec::new(), form_single_fragments(&[&r1, &m1, &r2, &m2]));
    assert_eq!(0, umi_groups.len());
    assert_eq!(4, single_fragments.len());

    // test 3: require unmapped reads to merge
    let read1 = labelled(
        create_sam_record(&read_name1, CHR_1, read_position, &bases, TEST_READ_CIGAR, CHR_2, read_position, false, false, None, false, Some(TEST_READ_CIGAR)),
        "read1",
    );
    let read2 = labelled(clone_sam_record(&read1, &read_name2), "read2");
    let (umi_groups, single_fragments) = poly_g_merge(&config, Vec::new(), form_single_fragments(&[&read1, &read2]));
    assert_eq!(0, umi_groups.len());
    assert_eq!(2, single_fragments.len());

    // test 4: some reads unmapped, some mapped elsewhere, mark unmapped as duplicates
    let (read2, mate2) = unmapped_mate_pair(&read_name2);
    let (read2, mate2) = (labelled(read2, "read2"), labelled(mate2, "mate2"));
    let (umi_groups, single_fragments) = poly_g_merge(&config, Vec::new(), form_single_fragments(&[&read1, &read2, &mate2]));
    assert_eq!(2, umi_groups.len());
    assert_eq!(0, single_fragments.len());

    let group = find_matching_group(&umi_groups, &["read1", "read2"], false).unwrap();
    assert!(!group.poly_g_unmapped);
    let group = find_matching_group(&umi_groups, &["mate2"], true).unwrap();
    assert!(group.poly_g_unmapped);
}

/// UmiGroupJitterTest.testIlluminaJitterUmiGroupCollapsedReadsNotUsedForConsensus
#[test]
fn test_illumina_jitter_umi_group_collapsed_reads_not_used_for_consensus() {
    let config = UmiConfig::new(true, true, "_", true, false);
    let bases = test_read_bases();
    let umi1 = "AAAAA_CCCCC";
    let mate_start = 1000;

    let read1 = create_sam_record(&next_read_id(umi1), CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, mate_start, false, false, None, true, Some(TEST_READ_CIGAR));
    let read2 = create_sam_record(
        &next_read_id(umi1), CHR_1, 100 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, &bases, TEST_READ_CIGAR, CHR_1, mate_start, false, false, None, true,
        Some(TEST_READ_CIGAR),
    );
    let single_fragments = form_single_fragments(&[&read1, &read2]);
    let mut stats = UmiStatistics::default();
    let (umi_groups, single_fragments) = process_umi_groups(&config, true, Vec::new(), single_fragments, false, &mut stats).unwrap();

    assert!(single_fragments.is_empty());
    assert_eq!(1, umi_groups.len());
    assert_eq!(1, umi_groups[0].reads.len());
    assert_eq!(1, umi_groups[0].non_consensus_reads.len());
}

/// UmiGroupJitterTest.testConsistentConsensusReadPairAfterIlluminaJitterUmiGroupCollapse
#[test]
fn test_consistent_consensus_read_pair_after_illumina_jitter_umi_group_collapse() {
    let names = names();
    let bases = test_read_bases();
    let umi_id_part1 = "A".repeat(6);
    let umi_id_part2 = "C".repeat(6);
    let forward_id = format!("READ_001:{umi_id_part1}{DEFAULT_DUPLEX_UMI_DELIM}{umi_id_part2}");
    let reverse_id = format!("READ_002:{umi_id_part2}{DEFAULT_DUPLEX_UMI_DELIM}{umi_id_part1}");

    let read1 = create_sam_record(&forward_id, CHR_1, 200, &bases, TEST_READ_CIGAR, CHR_2, 200, false, false, None, true, Some(TEST_READ_CIGAR));
    let mut read2 = create_sam_record(&reverse_id, CHR_1, 201, &bases, TEST_READ_CIGAR, CHR_2, 200, false, false, None, true, Some(TEST_READ_CIGAR));
    set_first_of_pair(&mut read2, false);

    let mut mate1 = create_sam_record(&forward_id, CHR_2, 200, &bases, TEST_READ_CIGAR, CHR_1, 200, true, false, None, false, Some(TEST_READ_CIGAR));
    set_first_of_pair(&mut mate1, false);
    let mate2 = create_sam_record(&reverse_id, CHR_2, 200, &bases, TEST_READ_CIGAR, CHR_1, 201, true, false, None, false, Some(TEST_READ_CIGAR));

    // duplicate group forming; the reads are in order of alignment start
    let mut read_cache = ReadCache::illumina(true);
    for r in [read1, read2] {
        read_cache.process_read(r, &names);
    }
    let chr1_reads = read_cache.evict_all().unwrap();
    let mut read_cache = ReadCache::illumina(true);
    for r in [mate1, mate2] {
        read_cache.process_read(r, &names);
    }
    let chr2_reads = read_cache.evict_all().unwrap();

    assert_eq!(0, chr1_reads.duplicate_groups.len());
    assert_eq!(2, chr1_reads.single_reads.len());
    assert_eq!(2, chr1_reads.total_read_count());
    assert_eq!(0, chr2_reads.duplicate_groups.len());
    assert_eq!(2, chr2_reads.single_reads.len());
    assert_eq!(2, chr2_reads.total_read_count());

    // umi group forming
    let config = UmiConfig::new(true, true, "_", true, false);
    let mut stats = UmiStatistics::default();
    let (mut chr1_umi_groups, chr1_single_reads) =
        process_umi_groups(&config, true, chr1_reads.duplicate_groups, chr1_reads.single_reads, false, &mut stats).unwrap();
    let (mut chr2_umi_groups, chr2_single_reads) =
        process_umi_groups(&config, true, chr2_reads.duplicate_groups, chr2_reads.single_reads, false, &mut stats).unwrap();

    assert_eq!(1, chr1_umi_groups.len());
    assert_eq!(0, chr1_single_reads.len());
    assert_eq!(1, chr2_umi_groups.len());
    assert_eq!(0, chr2_single_reads.len());

    let template_names = |groups: &mut Vec<DuplicateGroup>| -> Vec<String> {
        let mut v: Vec<String> = groups
            .iter_mut()
            .map(|g| {
                let i = select_template_read(&mut g.reads, &g.coords);
                read_name(&g.reads[i])
            })
            .collect();
        v.sort();
        v
    };
    assert_eq!(template_names(&mut chr1_umi_groups), template_names(&mut chr2_umi_groups));

    let name_sets = |groups: &[DuplicateGroup]| -> Vec<(Vec<String>, Vec<String>)> {
        let sorted_names = |reads: &[bam::Record]| {
            let mut v: Vec<String> = reads.iter().map(read_name).collect();
            v.sort();
            v
        };
        let mut v: Vec<(Vec<String>, Vec<String>)> = groups.iter().map(|g| (sorted_names(&g.reads), sorted_names(&g.non_consensus_reads))).collect();
        v.sort();
        v
    };
    assert_eq!(name_sets(&chr1_umi_groups), name_sets(&chr2_umi_groups));
}

/// The UMI statistics of two partitions merged are those of one pass over both
#[test]
fn test_umi_statistics_merge() {
    let config = UmiConfig::new(true, false, "", true, true);
    let bases = test_read_bases();
    let group = |pos: i32, umi: &str, count: usize| -> DuplicateGroup {
        let reads: Vec<bam::Record> = (0..count)
            .map(|_| create_sam_record(&next_read_id(umi), CHR_1, pos, &bases, TEST_READ_CIGAR, CHR_1, pos + 100, false, false, None, true, Some(TEST_READ_CIGAR)))
            .collect();
        let coords = FragmentCoords::from_read(&reads[0], &names(), true);
        DuplicateGroup::new(reads, coords)
    };
    let batches = || vec![vec![group(100, "AAAAAA", 2)], vec![group(300, "CCCCCC", 3)], vec![group(500, "GGGGGG", 3)], vec![group(700, "TTTTTT", 2)]];

    let mut whole = UmiStatistics::default();
    for b in batches() {
        process_umi_groups(&config, true, b, Vec::new(), true, &mut whole).unwrap();
    }
    let mut first = UmiStatistics::default();
    let mut second = UmiStatistics::default();
    for (i, b) in batches().into_iter().enumerate() {
        process_umi_groups(&config, true, b, Vec::new(), true, if i < 2 { &mut first } else { &mut second }).unwrap();
    }
    first.merge(&second);

    assert_eq!(4, whole.umi_groups);
    assert_eq!(whole.umi_groups, first.umi_groups);
    assert_eq!(1, whole.position_fragments.len());
    let (w, m) = (&whole.position_fragments[0], &first.position_fragments[0]);
    assert_eq!(4, w.frequency);
    assert_eq!(3, w.max_umi_reads_count);
    assert_eq!((w.frequency, w.max_umi_reads_count, w.max_coord_umi_count), (m.frequency, m.max_umi_reads_count, m.max_coord_umi_count));
    // the first group of the largest size gives the details
    assert!(w.umi_group_details.starts_with("1:300_"));
    assert!(m.umi_group_details.starts_with("1:300_") && m.umi_group_details.ends_with(":CCCCCC"));
    assert_eq!(whole.umi_group_frequencies.len(), first.umi_group_frequencies.len());
}

/// R24: with -umi_base_diff_stats a fragment without duplicates is counted as a group of one read. Java casts it
/// to the wrong class and the run ends.
#[test]
fn base_stats_take_a_fragment_without_duplicates() {
    let bases = test_read_bases();
    let read = create_sam_record(&next_read_id("AAAAA_CCCCC"), CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 1000, false, false, None, true, Some(TEST_READ_CIGAR));
    let run = |java: bool| {
        crate::compat::with_java(java, || {
            let config = UmiConfig::new(true, true, "_", true, true);
            let mut stats = UmiStatistics::default();
            process_umi_groups(&config, true, Vec::new(), form_single_fragments(&[&read]), false, &mut stats).map(|(groups, singles)| {
                (groups.len(), singles.len(), stats.umi_group_frequencies.clone())
            })
        })
    };
    assert!(run(true).is_err());
    let (groups, singles, frequencies) = run(false).unwrap();
    assert_eq!((groups, singles), (0, 1));
    // one coordinate group of one read, whose UMI is at distance 0 from itself
    assert_eq!(frequencies.len(), 1);
    assert_eq!((frequencies[0].duplicate_group_count, frequencies[0].read_count, frequencies[0].group_count), (1, 1, 1));
    assert_eq!(frequencies[0].edit_distance_frequency[0], 1);
}
