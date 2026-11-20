//! Port of duplicate/DuplicateGroupTest.

use super::util::*;
use crate::dup::{DuplicateGroup, FragmentCoords};
use crate::partition::set_primary_read;

/// DuplicateGroupTest.testNonConsensusPrimaryRead
///
/// DuplicateGroupBuilder.processDuplicateGroups is inlined in the port's PartitionReader; without UMIs and without
/// consensus its effect on a group is setPrimaryRead.
#[test]
fn test_non_consensus_primary_read() {
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();
    let read_pos = 100;
    let mate_pos = 200;

    let mut read1 =
        create_sam_record(&ids.next_id(), CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, None, true, Some(TEST_READ_CIGAR));
    set_base_qualities(&mut read1, 25);
    let fragment_coords = create_fragment_coords(&read1);

    let mut read2 =
        create_sam_record(&ids.next_id(), CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, None, true, Some(TEST_READ_CIGAR));
    set_base_qualities(&mut read2, 25);

    let mut read3 =
        create_sam_record(&ids.next_id(), CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, None, true, Some(TEST_READ_CIGAR));
    set_base_qualities(&mut read3, 11);

    let mut duplicate_group = DuplicateGroup::new(vec![read1, read2, read3], fragment_coords);

    set_primary_read(&mut duplicate_group);

    // isPrimaryRead(read1)
    assert_eq!(Some(0), duplicate_group.primary_read);
}

// CollapseUtils.collapseToNonOrientedKeyWithoutCoordinates serves the UMI and SBX duplicate collapsing, which the port
// does not have.

/// DuplicateGroupTest.testCollapseToKeyWithoutCoordinatesUnmappedRead
#[test]
#[ignore = "not ported: CollapseUtils (UMI / SBX duplicate collapsing) is not in the port"]
fn test_collapse_to_key_without_coordinates_unmapped_read() {}

/// DuplicateGroupTest.testCollapseToKeyWithoutCoordinatesUnmappedMate
#[test]
#[ignore = "not ported: CollapseUtils (UMI / SBX duplicate collapsing) is not in the port"]
fn test_collapse_to_key_without_coordinates_unmapped_mate() {}

/// DuplicateGroupTest.testCollapseToKeyWithoutCoordinatesSupplementaryReadWithUnmappedMate
#[test]
#[ignore = "not ported: CollapseUtils (UMI / SBX duplicate collapsing) is not in the port"]
fn test_collapse_to_key_without_coordinates_supplementary_read_with_unmapped_mate() {}

/// DuplicateGroupTest.testCollapseToKeyWithoutCoordinatesMappedPrimaryPair
#[test]
#[ignore = "not ported: CollapseUtils (UMI / SBX duplicate collapsing) is not in the port"]
fn test_collapse_to_key_without_coordinates_mapped_primary_pair() {}

/// DuplicateGroupTest.testCollapseToKeyWithoutCoordinatesSupplementaryRead
#[test]
#[ignore = "not ported: CollapseUtils (UMI / SBX duplicate collapsing) is not in the port"]
fn test_collapse_to_key_without_coordinates_supplementary_read() {}

/// DuplicateGroupTest.testCollapseToKeyWithoutCoordinatesMateHasSameFivePrimePos: the collapsed key is not ported
/// (see the ignored tests above), but the fragment coordinates the Java test checks first are.
#[test]
fn test_collapse_to_key_without_coordinates_mate_has_same_five_prime_pos_coords_only() {
    const FIVE_PRIME_POS: i32 = 200;
    let bases = test_read_bases();
    let names = names();
    for first_of_pair in [true, false] {
        for read_reversed in [true, false] {
            for mate_reversed in [true, false] {
                let read_start = if !read_reversed { FIVE_PRIME_POS } else { FIVE_PRIME_POS - (bases.len() as i32 - 1) };
                let mate_start = if !mate_reversed { FIVE_PRIME_POS } else { FIVE_PRIME_POS - (bases.len() as i32 - 1) };
                let mut read = create_sam_record(
                    "READ_001", CHR_1, read_start, &bases, TEST_READ_CIGAR, CHR_1, mate_start, read_reversed, false, None, mate_reversed,
                    Some(TEST_READ_CIGAR),
                );
                set_first_of_pair(&mut read, first_of_pair);

                let coords = FragmentCoords::from_read(&read, &names, true);
                assert_eq!(FIVE_PRIME_POS, coords.pos_lower);
                assert_eq!(FIVE_PRIME_POS, coords.pos_upper);
                assert_eq!(first_of_pair, coords.read_is_lower);
            }
        }
    }
}
