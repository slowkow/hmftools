//! Port of unmap/UnmapReadsTest.

use super::util::*;
use crate::rec::{self, FLAG_MATE_UNMAPPED, FLAG_UNMAPPED};
use crate::unmap::{clip_length_from_cigar_str, UnmapStats};
use rust_htslib::bam;

const SOFT_CLIPPED_READ_CIGAR: &str = "21S79M";

fn read_bases() -> String {
    "A".repeat(100)
}

fn read_unmapped(read: &bam::Record) -> bool {
    rec::has(read, FLAG_UNMAPPED)
}
fn mate_unmapped(read: &bam::Record) -> bool {
    rec::has(read, FLAG_MATE_UNMAPPED)
}
fn has_supp_attribute(read: &bam::Record) -> bool {
    rec::has_aux(read, b"SA")
}

#[allow(clippy::too_many_arguments)]
fn create(
    chr: &str, read_start: i32, cigar: &str, mate_chr: &str, mate_start: i32, is_reversed: bool, is_supplementary: bool,
    supp: Option<&SuppData>, mate_reversed: bool, mate_cigar: Option<&str>,
) -> bam::Record {
    create_sam_record(
        TEST_READ_ID, chr, read_start, &read_bases(), cigar, mate_chr, mate_start, is_reversed, is_supplementary, supp, mate_reversed, mate_cigar,
    )
}

/// UnmapReadsTest.testGetSoftClipCountFromCigarStr
#[test]
fn test_get_soft_clip_count_from_cigar_str() {
    assert_eq!(0, clip_length_from_cigar_str(Some(b"151M")));
    assert_eq!(10, clip_length_from_cigar_str(Some(b"10S141M")));
    assert_eq!(10, clip_length_from_cigar_str(Some(b"141M10S")));
    assert_eq!(12, clip_length_from_cigar_str(Some(b"2S138M10S")));
    assert_eq!(12, clip_length_from_cigar_str(Some(b"2H138M10H")));
}

/// UnmapReadsTest.testNoUnmap
#[test]
fn test_no_unmap() {
    // neither read in unmap regions
    let mut read = create(CHR_1, 100, TEST_READ_CIGAR, CHR_2, 600, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(!check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testNoUnmapSoftClipped
#[test]
fn test_no_unmap_soft_clipped() {
    // soft-clipped but not in an unmap region
    let mut read = create(CHR_1, 100, SOFT_CLIPPED_READ_CIGAR, CHR_2, 600, false, false, None, true, Some(SOFT_CLIPPED_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(!check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapReadSoftClipped
#[test]
fn test_unmap_read_soft_clipped() {
    // read unmapped since chimeric, sets read coords to mate
    let mut read = create(CHR_1, 550, SOFT_CLIPPED_READ_CIGAR, CHR_2, 600, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert!(read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert_eq!(600, rec::start(&read));
    assert_eq!(CHR_2, reference_name(&read));
}

/// UnmapReadsTest.testUnmapMateSoftClipped
#[test]
fn test_unmap_mate_soft_clipped() {
    // read unmapped since soft-clipped
    let mut read = create(CHR_1, 100, TEST_READ_CIGAR, CHR_1, 600, false, false, None, true, Some(SOFT_CLIPPED_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(mate_unmapped(&read));
    assert_eq!(100, rec::mate_start(&read));
    assert_eq!(CHR_1, mate_reference_name(&read));
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapMateAndSupplementarySoftClipped
#[test]
fn test_unmap_mate_and_supplementary_soft_clipped() {
    // mate unmapped and supp removed since both are in unmap regions and both have long SCs
    let supp = SuppData::new(CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60);
    let mut read = create(CHR_1, 100, TEST_READ_CIGAR, CHR_1, 600, false, false, Some(&supp), true, Some(SOFT_CLIPPED_READ_CIGAR));

    assert!(has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(mate_unmapped(&read));
    assert_eq!(100, rec::mate_start(&read));
    assert_eq!(CHR_1, mate_reference_name(&read));
    assert!(!has_supp_attribute(&read));
}

fn assert_fully_unplaced(read: &bam::Record) {
    assert!(read_unmapped(read));
    assert!(mate_unmapped(read));
    assert_eq!(NO_POSITION, rec::start(read));
    assert_eq!(NO_CHROMOSOME_NAME, reference_name(read));
    assert_eq!(NO_POSITION, rec::mate_start(read));
    assert_eq!(NO_CHROMOSOME_NAME, mate_reference_name(read));
}

/// UnmapReadsTest.testUnmapReadSoftClippedMateUnmapped
#[test]
fn test_unmap_read_soft_clipped_mate_unmapped() {
    let mut read = create(CHR_1, 550, SOFT_CLIPPED_READ_CIGAR, CHR_1, 550, false, false, None, true, None);
    rec::set_flag(&mut read, FLAG_MATE_UNMAPPED, true);

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert_fully_unplaced(&read);
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapMateSoftClippedReadUnmapped
#[test]
fn test_unmap_mate_soft_clipped_read_unmapped() {
    let mut read = create(CHR_1, 550, NO_CIGAR, CHR_1, 550, false, false, None, true, Some(SOFT_CLIPPED_READ_CIGAR));
    rec::set_flag(&mut read, FLAG_UNMAPPED, true);

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert_fully_unplaced(&read);
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapMateSoftClippedReadBecomesUnmapped
#[test]
fn test_unmap_mate_soft_clipped_read_becomes_unmapped() {
    let mut read = create(CHR_1, 550, TEST_READ_CIGAR, CHR_1, 600, false, false, None, true, Some(SOFT_CLIPPED_READ_CIGAR));
    let mut mate = create(CHR_1, 600, SOFT_CLIPPED_READ_CIGAR, CHR_1, 550, true, false, None, false, Some(TEST_READ_CIGAR));

    assert!(check_transform_read(&mut read, CHR_1));
    assert_fully_unplaced(&read);

    assert!(check_transform_read(&mut mate, CHR_1));
    assert_fully_unplaced(&mate);
}

/// UnmapReadsTest.testUnmapSupplementaryWhenPrimaryIsUnmappedSoftClipped
#[test]
fn test_unmap_supplementary_when_primary_is_unmapped_soft_clipped() {
    // note that the supplementary data of a supplementary read refers to the associated primary read
    let supp_read_data = SuppData::new(CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60);
    let mut read = create(CHR_2, 100, TEST_READ_CIGAR, CHR_2, 1000, false, true, Some(&supp_read_data), true, Some(TEST_READ_CIGAR));

    // when supplementaries are unmapped not all properties are unset since it will be dropped immediately
    assert!(has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_2));
    assert!(read_unmapped(&read));
}

/// UnmapReadsTest.testUnmapReadChimeric
#[test]
fn test_unmap_read_chimeric() {
    let mut read = create(CHR_1, 550, TEST_READ_CIGAR, CHR_2, 600, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert!(read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert_eq!(600, rec::start(&read));
    assert_eq!(CHR_2, reference_name(&read));
}

/// UnmapReadsTest.testUnmapMateChimeric
#[test]
fn test_unmap_mate_chimeric() {
    // the Java test passes the read bases as the mate cigar
    let mut read = create(CHR_2, 100, TEST_READ_CIGAR, CHR_1, 600, false, false, None, true, Some(&read_bases()));

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(mate_unmapped(&read));
    assert_eq!(100, rec::mate_start(&read));
    assert_eq!(CHR_2, mate_reference_name(&read));
}

/// UnmapReadsTest.testSupplementaryNotUnmappedChimeric
#[test]
fn test_supplementary_not_unmapped_chimeric() {
    // note that the supplementary data of a supplementary read refers to the associated primary read
    let supp_read_data = SuppData::new(CHR_2, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let mut read = create(CHR_1, 100, TEST_READ_CIGAR, CHR_2, 600, false, true, Some(&supp_read_data), true, Some(&read_bases()));

    assert!(has_supp_attribute(&read));

    assert!(!check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert!(has_supp_attribute(&read));
}

/// UnmapReadsTest.testSupplementaryNotUnmappedChimeric2
#[test]
fn test_supplementary_not_unmapped_chimeric2() {
    let supp_read_data = SuppData::new(CHR_5, 10000, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let mut read = create(CHR_2, 100, TEST_READ_CIGAR, CHR_2, 200, false, true, Some(&supp_read_data), true, Some(&read_bases()));

    assert!(has_supp_attribute(&read));

    assert!(!check_transform_read(&mut read, CHR_2));
    assert!(!read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert!(has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapSupplementaryWhenPrimaryIsUnmappedChimeric
#[test]
fn test_unmap_supplementary_when_primary_is_unmapped_chimeric() {
    // note that the supplementary data of a supplementary read refers to the associated primary read
    let supp_read_data = SuppData::new(CHR_1, 550, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let mut read = create(CHR_2, 100, TEST_READ_CIGAR, CHR_2, 1000, false, true, Some(&supp_read_data), true, Some(TEST_READ_CIGAR));

    assert!(has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_2));
    assert!(read_unmapped(&read));
}

/// UnmapReadsTest.testUnmapReadHighDepth
#[test]
fn test_unmap_read_high_depth() {
    let mut read = create(CHR_3, 600, TEST_READ_CIGAR, CHR_3, 701, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_3));
    assert!(read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert_eq!(701, rec::start(&read));
    assert_eq!(CHR_3, reference_name(&read));
}

/// UnmapReadsTest.testUnmapMateHighDepth
#[test]
fn test_unmap_mate_high_depth() {
    let mut read = create(CHR_3, 400, TEST_READ_CIGAR, CHR_3, 500, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_3));
    assert!(!read_unmapped(&read));
    assert!(mate_unmapped(&read));
    assert_eq!(400, rec::mate_start(&read));
    assert_eq!(CHR_3, mate_reference_name(&read));
}

/// UnmapReadsTest.testUnmapMateAndSupplementaryHighDepth
#[test]
fn test_unmap_mate_and_supplementary_high_depth() {
    let supp = SuppData::new(CHR_3, 600, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let mut read = create(CHR_3, 400, TEST_READ_CIGAR, CHR_3, 500, false, false, Some(&supp), true, Some(TEST_READ_CIGAR));

    assert!(has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_3));
    assert!(!read_unmapped(&read));
    assert!(mate_unmapped(&read));
    assert_eq!(400, rec::mate_start(&read));
    assert_eq!(CHR_3, mate_reference_name(&read));
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapMultipleSupplementaryAlignments
#[test]
fn test_unmap_multiple_supplementary_alignments() {
    let mut read = create(CHR_1, 100, TEST_READ_CIGAR, CHR_3, 200, false, false, None, true, Some(TEST_READ_CIGAR));

    let supp_data1 = SuppData::new(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let supp_data2 = SuppData::new(CHR_2, 6000, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let supp_data3 = SuppData::new(CHR_3, 600, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    set_string_attribute(&mut read, b"SA", &alignments_to_sam_tag(&[supp_data1, supp_data2, supp_data3]));

    assert!(has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_1));
    assert!(!read_unmapped(&read));
    assert!(!mate_unmapped(&read));
    assert!(!has_supp_attribute(&read));
}

/// UnmapReadsTest.testUnmapSupplementaryWhenUnmappingAssociatedSupplementary
#[test]
fn test_unmap_supplementary_when_unmapping_associated_supplementary() {
    let mut read = create(CHR_1, 100, TEST_READ_CIGAR, CHR_3, 300, false, true, None, true, Some(TEST_READ_CIGAR));

    let supp_data_primary = SuppData::new(CHR_3, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let supp_data_other = SuppData::new(CHR_3, 500, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    set_string_attribute(&mut read, b"SA", &alignments_to_sam_tag(&[supp_data_primary, supp_data_other]));

    assert!(has_supp_attribute(&read));

    assert!(check_transform_read(&mut read, CHR_3));
    assert!(read_unmapped(&read));
}

/// UnmapReadsTest.testUnmapRegionState
#[test]
fn test_unmap_region_state() {
    let unmapper = test_read_unmapper();
    let stats = &mut UnmapStats::default();
    let mut region_state = unmap_region_state(CHR_4, 1, 1000000, chr_location_map().remove(CHR_4).unwrap());

    // primary read
    let mut read = create(CHR_4, 1000, TEST_READ_CIGAR, CHR_5, 100, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(unmapper.check_transform_read(&mut read, &mut region_state, stats));
    assert_eq!(Some(0), region_state.last_matched_index);

    read = create(CHR_4, 2100, TEST_READ_CIGAR, CHR_5, 100, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!unmapper.check_transform_read(&mut read, &mut region_state, stats));
    assert_eq!(Some(0), region_state.last_matched_index);

    read = create(CHR_4, 5100, TEST_READ_CIGAR, CHR_5, 100, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(unmapper.check_transform_read(&mut read, &mut region_state, stats));
    assert_eq!(Some(2), region_state.last_matched_index);

    read = create(CHR_4, 10000, TEST_READ_CIGAR, CHR_5, 100, false, false, None, true, Some(TEST_READ_CIGAR));

    assert!(!unmapper.check_transform_read(&mut read, &mut region_state, stats));
    assert_eq!(Some(3), region_state.last_matched_index);
}
