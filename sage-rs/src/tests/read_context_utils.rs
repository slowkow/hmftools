//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! common/ReadContextUtilsTest.java

use super::util::*;
use crate::constants::DEFAULT_FLANK_LENGTH;
use crate::read_context::{create_context, ReadCigarInfo};
use crate::types::*;

fn build_read_cigar(
    read: &Read, left_flank_index: i32, left_core_index: i32, right_core_index: i32, right_flank_index: i32,
) -> ReadCigarInfo {
    ReadCigarInfo::build(read.start, &read.cigar, left_flank_index, left_core_index, right_core_index, right_flank_index)
}

// ReadContextUtilsTest.testReadBasesCigar
#[test]
fn test_read_bases_cigar() {
    let read_bases = REF_BASES_200[0..40].to_string();
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "40M";
    let read = build_sam_record_quals(100, read_cigar, &read_bases, base_quals.clone());

    let flank_size = DEFAULT_FLANK_LENGTH;

    // test 1: basic core in aligned section
    let mut read_flank_start = 10;
    let mut read_flank_end = 34;
    let mut read_core_start = read_flank_start + flank_size;
    let mut read_core_end = read_flank_end - flank_size;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("25M", cigar_str(&info.cigar));
    assert_eq!(110, info.flank_position_start);
    assert_eq!(120, info.core_position_start);
    assert_eq!(124, info.core_position_end);
    assert_eq!(134, info.flank_position_end);
    assert_eq!(10, info.flank_index_start);
    assert_eq!(34, info.flank_index_end);

    // test 2: flanks ending in soft-clipped bases
    let read_cigar = "20S10M10S";
    let read = build_sam_record_quals(100, read_cigar, &read_bases, base_quals.clone());

    // test 2: flanks both going into the soft-clips
    read_flank_start = 10;
    read_flank_end = 34;
    read_core_start = read_flank_start + flank_size;
    read_core_end = read_flank_end - flank_size;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("10S10M5S", cigar_str(&info.cigar));
    assert_eq!(90, info.flank_position_start);
    assert_eq!(100, info.core_position_start);
    assert_eq!(104, info.core_position_end);
    assert_eq!(114, info.flank_position_end);
    assert_eq!(10, info.flank_index_start);
    assert_eq!(34, info.flank_index_end);

    // test 3: flanks straddling deletes
    let read_cigar = "10M10D20M10D10M";
    let read = build_sam_record_quals(100, read_cigar, &read_bases, base_quals.clone());

    // elements:
    // 100-109  10M   index 0-9
    // 110-119  10D   no change
    // 120-139  20M   index 10-29
    // 140-149  10D   no change
    // 150-159  10M   index 30-39

    read_flank_start = 8;
    read_flank_end = 32;
    read_core_start = read_flank_start + flank_size;
    read_core_end = read_flank_end - flank_size;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("2M10D20M10D3M", cigar_str(&info.cigar));
    assert_eq!(108, info.flank_position_start);
    assert_eq!(128, info.core_position_start);
    assert_eq!(132, info.core_position_end);
    assert_eq!(152, info.flank_position_end);
    assert_eq!(8, info.flank_index_start);
    assert_eq!(32, info.flank_index_end);

    // test 4: cores end in deletes
    let read_cigar = "20M10D20M10D20M";
    let read = build_sam_record_quals(100, read_cigar, &read_bases, base_quals.clone());

    // elements:
    // 100-119  20M   index 0-19
    // 120-129  10D   no change
    // 130-149  20M   index 20-39
    // 150-159  10D   no change
    // 160-179  10M   index 40-59

    read_flank_start = 8;
    read_core_start = 18;
    read_core_end = 42;
    read_flank_end = 52;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("12M10D20M10D13M", cigar_str(&info.cigar));
    assert_eq!(108, info.flank_position_start);
    assert_eq!(118, info.core_position_start);
    assert_eq!(162, info.core_position_end);
    assert_eq!(172, info.flank_position_end);
    assert_eq!(8, info.flank_index_start);
    assert_eq!(52, info.flank_index_end);

    // test 4b: cores start in deletes
    // elements:
    // 100-120    21M     index 0-20
    // 120-121    1D      no change
    // 122-146    25M     index 21-45

    let read_bases =
        format!("{}{}A{}", &REF_BASES_200[100..121], &REF_BASES_200[122..124], &REF_BASES_200[125..147]);
    let read_cigar = "21M1D25M";
    let read = build_sam_record(100, read_cigar, &read_bases);

    let var = SimpleVariant::new(124, b"T", b"A");
    let read_context = create_context(&var, &read, 23, &ref_sequence_200(), DEFAULT_FLANK_LENGTH).unwrap();
    assert_eq!(122, read_context.core_position_start);

    // repeat just testing the cigar building
    read_flank_start = 11;
    read_flank_end = 35;
    read_core_start = read_flank_start + flank_size;
    read_core_end = read_flank_end - flank_size;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!(122, info.core_position_start);

    // test 5: flanks ending in inserts at both ends
    let read_cigar = "10M5I20M5I10M";
    let read = build_sam_record_quals(100, read_cigar, &read_bases, base_quals.clone());

    // elements:
    // 100-109  10M   index 0-9
    // unch     5I    index 10-14
    // 110-129  20M   index 15-34
    // unc      5I    index 35-39
    // 130-139  10M   index 40-49

    read_flank_start = 11;
    read_flank_end = 36;
    read_core_start = read_flank_start + flank_size;
    read_core_end = read_flank_end - flank_size;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("1M5I20M5I1M", cigar_str(&info.cigar));
    assert_eq!(109, info.flank_position_start);
    assert_eq!(116, info.core_position_start);
    assert_eq!(121, info.core_position_end);
    assert_eq!(130, info.flank_position_end);
    assert_eq!(9, info.flank_index_start);
    assert_eq!(40, info.flank_index_end);

    // test 6: core positions ending in inserts
    read_flank_start = 3;
    read_core_start = 13;
    read_core_end = 36;
    read_flank_end = 46;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("7M5I20M5I7M", cigar_str(&info.cigar));
    assert_eq!(103, info.flank_position_start);
    assert_eq!(109, info.core_position_start);
    assert_eq!(130, info.core_position_end);
    assert_eq!(136, info.flank_position_end);
    assert_eq!(3, info.flank_index_start);
    assert_eq!(46, info.flank_index_end);
}

// ReadContextUtilsTest.testReadBasesCigar2
#[test]
fn test_read_bases_cigar2() {
    let read_bases = &REF_BASES_200[0..51];
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "10M2D10M1I30M";
    let read = build_sam_record_quals(100, read_cigar, read_bases, base_quals);

    let flank_size = DEFAULT_FLANK_LENGTH;

    // test 1: basic core in aligned section
    let read_flank_start = 20;
    let read_flank_end = 44;
    let read_core_start = read_flank_start + flank_size;
    let read_core_end = read_flank_end - flank_size;

    let info = build_read_cigar(&read, read_flank_start, read_core_start, read_core_end, read_flank_end);
    assert!(info.is_valid());
    assert_eq!("1M1I24M", cigar_str(&info.cigar));
    assert_eq!(121, info.flank_position_start);
    assert_eq!(131, info.core_position_start);
    assert_eq!(135, info.core_position_end);
    assert_eq!(145, info.flank_position_end);
    assert_eq!(19, info.flank_index_start);
    assert_eq!(44, info.flank_index_end);
}
