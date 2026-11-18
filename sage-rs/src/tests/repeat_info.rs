//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! common/RepeatInfoTest.java

use crate::constants::{MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT};
use crate::repeat::*;

fn repeat(index: i32, bases: &str, count: i32) -> RepeatInfo {
    RepeatInfo::new(index, bases.as_bytes().to_vec(), count)
}

/// RepeatInfo.addIfUnique returns whether the repeat was added; the port returns nothing
fn add_if_unique_checked(repeats: &mut Vec<RepeatInfo>, new: RepeatInfo) -> bool {
    let before = repeats.len();
    add_if_unique(repeats, new);
    repeats.len() == before + 1
}

// RepeatInfoTest.testRepeats
#[test]
fn test_repeats() {
    //           0123456789
    let bases = "AAACCTTTTT";

    // first check limits
    let repeat_info = find_multi_base_repeat(bases.as_bytes(), 6, 1, 4);
    assert!(repeat_info.is_some());
    let repeat_info = repeat_info.unwrap();
    assert_eq!("T", repeat_info.bases_str());
    assert_eq!(4, repeat_info.count);

    let repeat_info = find_multi_base_repeat(bases.as_bytes(), 7, 1, 4);
    assert!(repeat_info.is_none());

    //           01234567890
    let bases = "AACGGTACGGT";

    let repeat_info = find_multi_base_repeat(bases.as_bytes(), 1, 5, 2);
    assert!(repeat_info.is_some());
    let repeat_info = repeat_info.unwrap();
    assert_eq!("ACGGT", repeat_info.bases_str());
    assert_eq!(2, repeat_info.count);
    assert_eq!(10, repeat_info.total_length());
    assert_eq!(5, repeat_info.repeat_length());

    let repeat_info = find_multi_base_repeat(bases.as_bytes(), 2, 5, 2);
    assert!(repeat_info.is_none());
}

// RepeatInfoTest.testMaxRepeats
#[test]
fn test_max_repeats() {
    //           01234567890123456789
    let bases = "ACGTACGTTTTTTGTACGT";

    let max_length = 5;
    let min_count = 3;

    let index_end = bases.len() as i32 - 1;

    let repeat = find_max_repeat(bases.as_bytes(), 0, index_end, max_length, min_count, false, -1);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("T", repeat.bases_str());
    assert_eq!(7, repeat.index);
    assert_eq!(6, repeat.count);

    let repeat = find_max_repeat(bases.as_bytes(), 10, index_end, max_length, min_count, true, -1);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("T", repeat.bases_str());
    assert_eq!(7, repeat.index);
    assert_eq!(6, repeat.count);

    // longer repeats

    //           01234567890123456789
    let bases = "ACGTACGAGAGAGGTACGT";

    let repeat = find_max_repeat(bases.as_bytes(), 0, index_end, max_length, min_count, true, -1);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("GA", repeat.bases_str());
    assert_eq!(6, repeat.index);
    assert_eq!(3, repeat.count);

    //           0123456789012345678901234567890123456789
    let bases = "ACGTACACGGAAAGGAAAGGAAAGGAAAGGAAAGTACGT";

    let repeat = find_max_repeat(bases.as_bytes(), 0, index_end, max_length, min_count, true, -1);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("GGAAA", repeat.bases_str());
    assert_eq!(8, repeat.index);
    assert_eq!(5, repeat.count);

    // finding and extending the start
    let repeat = find_max_repeat(bases.as_bytes(), 17, index_end, max_length, min_count, true, -1);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("GGAAA", repeat.bases_str());
    assert_eq!(8, repeat.index);
    assert_eq!(5, repeat.count);

    // longest when there are lots to choose from

    //                     10        20        30        40        50
    //           012345678901234567890123456789012345678901234567890123456789
    let bases = "AAAACGCGCGCGTTTTTGCTTTGCTTTGCTTTAACGTACGTACGTACGTCCDDCCCDDD";

    let repeat = find_max_repeat(bases.as_bytes(), 0, bases.len() as i32 - 1, max_length, min_count, false, 20);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("TTTGC", repeat.bases_str());
    assert_eq!(14, repeat.index);
    assert_eq!(3, repeat.count);

    // checks that the repeat crosses the required index
    let repeat = find_max_repeat(bases.as_bytes(), 0, bases.len() as i32 - 1, max_length, min_count, false, 40);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("ACGT", repeat.bases_str());
    assert_eq!(33, repeat.index);
    assert_eq!(4, repeat.count);

    let repeat = find_max_repeat(bases.as_bytes(), 0, bases.len() as i32 - 1, max_length, min_count, false, 8);

    assert!(repeat.is_some());
    let repeat = repeat.unwrap();
    assert_eq!("CG", repeat.bases_str());
    assert_eq!(4, repeat.index);
    assert_eq!(4, repeat.count);
}

// RepeatInfoTest.testFindRepeatBoundaries
#[test]
fn test_find_repeat_boundaries() {
    //                0123456789012345678901234567890123456789
    let read_bases = "ACGTAGCTTGCGCGCGCCACACACTACGTAGCT";

    // test 1: a transitioning repeat where the lower repeat is the max

    let read_core_start = 13;
    let read_core_end = 18;

    let repeat_boundaries =
        find_repeat_boundaries(read_bases.as_bytes(), read_core_start, read_core_end, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

    let max_repeat = repeat_boundaries.max_repeat.as_ref().unwrap();
    assert_eq!("GC", max_repeat.bases_str());
    assert_eq!(4, max_repeat.count);
    assert_eq!(8, repeat_boundaries.lower_index);
    assert_eq!(24, repeat_boundaries.upper_index);

    // test 2: now where the upper repeat is the max but the second is still used to set the lower index boundary

    //                0123456789012345678901234567890123456789
    let read_bases = "ACGTAGCTTGCGCGCCACACACACTACGTAGCT";

    let read_core_start = 11;
    let read_core_end = 16;

    let repeat_boundaries =
        find_repeat_boundaries(read_bases.as_bytes(), read_core_start, read_core_end, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

    let max_repeat = repeat_boundaries.max_repeat.as_ref().unwrap();
    assert_eq!("CA", max_repeat.bases_str());
    assert_eq!(4, max_repeat.count);
    assert_eq!(8, repeat_boundaries.lower_index);
    assert_eq!(24, repeat_boundaries.upper_index);

    // test 3: with partial repeats at the end which need to be included
    //                0123456789012345678901234567890123456789
    let read_bases = "AACTTTCTTTCTTTCTTAGGGGCGGGGCGGGGCGGGGTT";

    let read_core_start = 15;
    let read_core_end = 19;

    let repeat_boundaries =
        find_repeat_boundaries(read_bases.as_bytes(), read_core_start, read_core_end, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

    let max_repeat = repeat_boundaries.max_repeat.as_ref().unwrap();
    assert_eq!("G", max_repeat.bases_str());
    assert_eq!(4, max_repeat.count);
    assert_eq!(1, repeat_boundaries.lower_index);
    assert_eq!(37, repeat_boundaries.upper_index);

    // test 4: 3 repeats, the inner one the longest but requiring use of the outer 2

    //                          10        20        30
    //                0123456789012345678901234567890123456789
    let read_bases = "ACTGAAAAAAAACAAAAAAAAAAACTATATATATATTACT";

    let read_core_start = 11;
    let read_core_end = 25;

    let repeat_boundaries =
        find_repeat_boundaries(read_bases.as_bytes(), read_core_start, read_core_end, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

    let max_repeat = repeat_boundaries.max_repeat.as_ref().unwrap();
    assert_eq!("A", max_repeat.bases_str());
    assert_eq!(11, max_repeat.count);
    assert_eq!(3, repeat_boundaries.all_repeats.as_ref().unwrap().len());
    assert_eq!(3, repeat_boundaries.lower_index);
    assert_eq!(36, repeat_boundaries.upper_index);
}

// RepeatInfoTest.testUniqueRepeats
#[test]
fn test_unique_repeats() {
    let mut repeats: Vec<RepeatInfo> = Vec::new();

    let r = repeat(0, "A", 7);
    let r2 = repeat(0, "AAA", 2);
    assert!(r.is_multiple_of(&r2));

    let r2 = repeat(0, "AA", 3);
    assert!(r.is_multiple_of(&r2));

    let r = repeat(0, "ACT", 4);
    let r2 = repeat(0, "ACTACT", 2);
    assert!(r.is_multiple_of(&r2));

    let r = repeat(0, "ACT", 4);
    let r2 = repeat(0, "ACTAGT", 2);
    assert!(!r.is_multiple_of(&r2));

    // too long
    let r = repeat(0, "ACT", 4);
    let r2 = repeat(0, "ACTACT", 3);
    assert!(!r.is_multiple_of(&r2));

    // alternative starts for the same repeat
    let r = repeat(0, "ACT", 4);
    let r2 = repeat(1, "CTA", 3);
    let r3 = repeat(2, "TAC", 3);
    assert!(r.is_alternative_start(&r2));
    assert!(r.is_alternative_start(&r3));

    // not if longer or starts outside first repeat's bounds
    let r2 = repeat(1, "CTA", 5);
    assert!(!r.is_alternative_start(&r2));

    let r2 = repeat(3, "CTA", 3);
    assert!(!r.is_alternative_start(&r2));

    // or different bases
    let r2 = repeat(1, "CTG", 3);
    assert!(!r.is_alternative_start(&r2));

    assert!(add_if_unique_checked(&mut repeats, repeat(0, "ACT", 3)));
    assert!(add_if_unique_checked(&mut repeats, repeat(0, "A", 9)));
    assert!(!add_if_unique_checked(&mut repeats, repeat(0, "AA", 2)));
    assert!(!add_if_unique_checked(&mut repeats, repeat(1, "AA", 2)));
    assert!(!add_if_unique_checked(&mut repeats, repeat(1, "CTA", 3)));
}
