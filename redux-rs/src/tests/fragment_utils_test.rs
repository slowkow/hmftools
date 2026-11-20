//! Port of FragmentUtilsTest.

use super::util::*;
use crate::rec::{self, FLAG_PAIRED, FLAG_REVERSE, FLAG_UNMAPPED};

const FORWARD: bool = true;
const REVERSE: bool = false;

/// FragmentUtilsTest.testFragmentCoords
#[test]
fn test_fragment_coords() {
    let bases = test_read_bases();
    let mut read =
        create_sam_record(TEST_READ_ID, CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 200, false, false, None, true, Some(TEST_READ_CIGAR));

    let mut fragment_coords = create_fragment_coords(&read);

    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(299, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert_eq!(REVERSE, fragment_coords.upper_forward);
    assert!(fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:100_1:299:R_L", fragment_coords.key_non_oriented());

    flip_first_in_pair(&mut read);

    fragment_coords = create_fragment_coords(&read);
    assert!(!fragment_coords.frag_forward);

    // a mate upper read
    read = create_sam_record(TEST_READ_ID, CHR_1, 200, &bases, TEST_READ_CIGAR, CHR_1, 100, true, false, None, false, Some(TEST_READ_CIGAR));

    flip_first_in_pair(&mut read);

    fragment_coords = create_fragment_coords(&read);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(299, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert_eq!(REVERSE, fragment_coords.upper_forward);
    assert!(!fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:100_1:299:R_U", fragment_coords.key_non_oriented());

    // reverse orientation fragment with soft-clips at both ends
    read = create_sam_record(TEST_READ_ID, CHR_1, 100, &bases, "80M20S", CHR_1, 400, true, false, None, false, Some("10S80M"));

    fragment_coords = create_fragment_coords(&read);
    assert_eq!(199, fragment_coords.pos_lower);
    assert_eq!(390, fragment_coords.pos_upper);
    assert_eq!(REVERSE, fragment_coords.lower_forward);
    assert_eq!(FORWARD, fragment_coords.upper_forward);
    assert!(fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:199:R_1:390_L", fragment_coords.key_non_oriented());

    // using supplementary data to get primary coords
    let supp_alignment = SuppData::new(CHR_1, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

    read = create_sam_record(
        TEST_READ_ID, CHR_2, 1000, &bases, TEST_READ_CIGAR, CHR_1, 200, false, true, Some(&supp_alignment), true, Some(TEST_READ_CIGAR),
    );

    fragment_coords = create_fragment_coords(&read);
    assert!(fragment_coords.supp.is_some());
    assert_eq!(1000, fragment_coords.supp.unwrap().0);
    assert_eq!(FORWARD, fragment_coords.supp.unwrap().1);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(299, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert_eq!(REVERSE, fragment_coords.upper_forward);
    assert!(fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:100_1:299:R_L_S", fragment_coords.key_non_oriented());

    // test again with an unmapped mate
    read = create_sam_record_basic(TEST_READ_ID, CHR_1, 100, &bases, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, 100, false, false, None);

    fragment_coords = create_fragment_coords(&read);
    assert!(!fragment_coords.unmapped_sourced);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(NO_POSITION, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert!(fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:100", fragment_coords.key_non_oriented());

    rec::set_flag(&mut read, FLAG_REVERSE, true);
    fragment_coords = create_fragment_coords(&read);
    assert_eq!(REVERSE, fragment_coords.lower_forward);
    assert_eq!("1:199:R", fragment_coords.key_non_oriented());

    // the unmapped read
    read = create_sam_record(
        TEST_READ_ID, NO_CHROMOSOME_NAME, 100, &bases, NO_CIGAR, CHR_1, 100, false, false, None, false, Some(TEST_READ_CIGAR),
    );
    rec::set_flag(&mut read, FLAG_UNMAPPED, true);

    fragment_coords = create_fragment_coords(&read);
    assert!(fragment_coords.unmapped_sourced);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(NO_POSITION, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert!(!fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:100_N", fragment_coords.key_non_oriented());

    // an unpaired read
    read = create_sam_record(
        TEST_READ_ID, CHR_1, 100, &bases, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, None, false, Some(NO_CIGAR),
    );
    rec::set_flag(&mut read, FLAG_PAIRED, false);

    fragment_coords = create_fragment_coords(&read);
    assert!(fragment_coords.unpaired);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(199, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert!(fragment_coords.read_is_lower);
    assert!(fragment_coords.frag_forward);
    assert_eq!("1:100_199", fragment_coords.key_non_oriented());

    // matching 5' positions, determine by first-in-pair
    read = create_sam_record(TEST_READ_ID, CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 100, false, false, None, false, Some(TEST_READ_CIGAR));

    fragment_coords = create_fragment_coords(&read);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(100, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert_eq!(FORWARD, fragment_coords.upper_forward);
    assert!(fragment_coords.read_is_lower);
    assert_eq!("1:100_1:100_L", fragment_coords.key_non_oriented());

    flip_first_in_pair(&mut read);

    fragment_coords = create_fragment_coords(&read);
    assert_eq!(100, fragment_coords.pos_lower);
    assert_eq!(100, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert_eq!(FORWARD, fragment_coords.upper_forward);
    assert!(!fragment_coords.read_is_lower);
    assert_eq!("1:100_1:100_U", fragment_coords.key_non_oriented());
}

/// The checks shared by the unpaired tests of FragmentUtilsTest
#[allow(clippy::too_many_arguments)]
fn check_unpaired(
    read_reversed: bool, supp: Option<(char, i32)>, expected_key: &str, expected_chr: &str, position_lower: i32, position_upper: i32,
    read_is_lower: bool, expected_supp: Option<(i32, bool)>, read_position: i32, read_orientation: bool,
) {
    let bases = test_read_bases();
    let read = match supp {
        None => create_sam_record_unpaired(TEST_READ_ID, CHR_1, 100, &bases, TEST_READ_CIGAR, read_reversed, false, None),
        Some((strand, supp_read_pos)) => {
            let primary = SuppData::new(CHR_1, 1_000, strand, TEST_READ_CIGAR, 0);
            create_sam_record_unpaired(TEST_READ_ID, CHR_2, supp_read_pos, &bases, TEST_READ_CIGAR, read_reversed, true, Some(&primary))
        }
    };
    let fragment_coords = create_fragment_coords(&read);

    assert_eq!(expected_chr, fragment_coords.chr_lower);
    assert_eq!(expected_chr, fragment_coords.chr_upper);
    assert_eq!(position_lower, fragment_coords.pos_lower);
    assert_eq!(position_upper, fragment_coords.pos_upper);
    assert_eq!(FORWARD, fragment_coords.lower_forward);
    assert_eq!(REVERSE, fragment_coords.upper_forward);
    assert_eq!(FORWARD, fragment_coords.frag_forward);
    assert_eq!(read_is_lower, fragment_coords.read_is_lower);
    assert_eq!(expected_supp, fragment_coords.supp);
    assert!(!fragment_coords.unmapped_sourced);
    assert!(fragment_coords.unpaired);
    assert_eq!(expected_key, fragment_coords.key);

    assert_eq!(read_position, fragment_coords.read_position());
    assert_eq!(read_orientation, fragment_coords.read_forward());
}

/// FragmentUtilsTest.testFragmentCoordsUnpairedForward
#[test]
fn test_fragment_coords_unpaired_forward() {
    check_unpaired(false, None, "1:100_199", CHR_1, 100, 199, true, None, 100, FORWARD);
}

/// FragmentUtilsTest.testFragmentCoordsUnpairedReverse
#[test]
fn test_fragment_coords_unpaired_reverse() {
    check_unpaired(true, None, "1:100_199_R", CHR_1, 100, 199, false, None, 199, REVERSE);
}

/// FragmentUtilsTest.testFragmentCoordsUnpairedSuppForwardPrimaryForward
#[test]
fn test_fragment_coords_unpaired_supp_forward_primary_forward() {
    check_unpaired(false, Some((SUPP_POS_STRAND, 2_000)), "1:1000_1099_S", CHR_1, 1000, 1099, true, Some((2_000, FORWARD)), 2_000, FORWARD);
}

/// FragmentUtilsTest.testFragmentCoordsUnpairedSuppForwardPrimaryReverse
#[test]
fn test_fragment_coords_unpaired_supp_forward_primary_reverse() {
    check_unpaired(false, Some((SUPP_NEG_STRAND, 2_000)), "1:1000_1099_R_S", CHR_1, 1000, 1099, false, Some((2_000, FORWARD)), 2_000, FORWARD);
}

/// FragmentUtilsTest.testFragmentCoordsUnpairedSuppReversePrimaryForward
#[test]
fn test_fragment_coords_unpaired_supp_reverse_primary_forward() {
    check_unpaired(true, Some((SUPP_POS_STRAND, 2_000)), "1:1000_1099_S", CHR_1, 1000, 1099, true, Some((2_099, REVERSE)), 2_099, REVERSE);
}

/// FragmentUtilsTest.testFragmentCoordsUnpairedSuppReversePrimaryReverse
#[test]
fn test_fragment_coords_unpaired_supp_reverse_primary_reverse() {
    check_unpaired(true, Some((SUPP_NEG_STRAND, 2_000)), "1:1000_1099_R_S", CHR_1, 1000, 1099, false, Some((2_099, REVERSE)), 2_099, REVERSE);
}
