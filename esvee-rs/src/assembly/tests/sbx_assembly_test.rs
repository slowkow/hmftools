//! SbxAssemblyTest
use super::*;
use crate::assembly::common_utils::is_higher_base_qual_category;
use crate::assembly::constants::{BASE_QUAL_MINIMUM, LOW_BASE_QUAL_THRESHOLD};
use crate::assembly::seq_tech_utils::{pass_sbx_distinct_prime_positions_filter, trim_sbx_uncertain_bases};
use crate::assembly::types::support_read::{SupP, SupportRead};

// hmf-common SbxBamUtils
const SBX_SIMPLEX_QUAL: u8 = 27;
const SBX_DUPLEX_QUAL: u8 = 40;

fn sbx_env() -> TestEnv {
    let env = test_env();
    set_sbx_sequencing();
    env
}

#[test]
fn test_base_qual_comparisons() {
    let _env = sbx_env();
    assert!(is_higher_base_qual_category(SBX_SIMPLEX_QUAL, LOW_BASE_QUAL_THRESHOLD));
    assert!(!is_higher_base_qual_category(SBX_SIMPLEX_QUAL, SBX_SIMPLEX_QUAL));
    assert!(!is_higher_base_qual_category(SBX_SIMPLEX_QUAL, SBX_DUPLEX_QUAL));
    assert!(is_higher_base_qual_category(SBX_DUPLEX_QUAL, SBX_SIMPLEX_QUAL));
}

#[test]
fn test_trim_uncertain_bases() {
    let _env = sbx_env();
    let read_bases = ref_bases_random_100();
    let mut base_qualities = build_default_base_quals(read_bases.len());
    // low qual will be 70% of outer bases
    let n = base_qualities.len();
    for i in 3..10 {
        base_qualities[i] = BASE_QUAL_MINIMUM;
        base_qualities[n - i - 1] = BASE_QUAL_MINIMUM;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    read.m().bam_record_mut().base_qualities = base_qualities;
    trim_sbx_uncertain_bases(read.m());
    assert_eq!(10, read.trim_count_start());
    assert_eq!(10, read.trim_count_end());

    // a read that can trim from both sides past each other
    let mut base_qualities = build_default_base_quals(read_bases.len());
    for q in base_qualities.iter_mut().take(80) {
        *q = BASE_QUAL_MINIMUM;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    read.m().bam_record_mut().base_qualities = base_qualities;
    trim_sbx_uncertain_bases(read.m());
    assert_eq!(80, read.trim_count_start());
    assert_eq!(0, read.trim_count_end());

    // and the other side
    let mut base_qualities = build_default_base_quals(read_bases.len());
    for q in base_qualities.iter_mut().skip(20) {
        *q = BASE_QUAL_MINIMUM;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    read.m().bam_record_mut().base_qualities = base_qualities;
    trim_sbx_uncertain_bases(read.m());
    assert_eq!(0, read.trim_count_start());
    assert_eq!(80, read.trim_count_end());
}

#[test]
fn test_read_prime_range() {
    let _env = sbx_env();
    let read_bases = ref_bases_random_100();
    let cigar = "60M40S";
    let junc_read_distance = 60; // irrelevant for this test
    let support = |read: &ReadP| -> SupP { SupportRead::new(read, SupportType::Junction, junc_read_distance, 0, 0) };
    let mut support_reads: Vec<SupP> = [100, 101, 102].iter().map(|start| support(&create_read(TEST_READ_ID, *start, &read_bases, cigar))).collect();
    assert!(!pass_sbx_distinct_prime_positions_filter(&support_reads));
    support_reads.push(support(&create_read(TEST_READ_ID, 104, &read_bases, cigar)));
    assert!(pass_sbx_distinct_prime_positions_filter(&support_reads));

    // reverse orientation reads
    let cigar = "40S60M";
    let reverse_read = |start: i32| -> ReadP {
        let read = create_read(TEST_READ_ID, start, &read_bases, cigar);
        read.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
        read
    };
    let mut support_reads: Vec<SupP> = [102, 101, 100].iter().map(|start| support(&reverse_read(*start))).collect();
    assert!(!pass_sbx_distinct_prime_positions_filter(&support_reads));
    support_reads.push(support(&reverse_read(98)));
    assert!(pass_sbx_distinct_prime_positions_filter(&support_reads));
}
