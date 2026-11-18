//! Ported from the Java SAGE test suite: seqtech/UltimaMiscUtilsTest.java

use super::util::*;
use crate::constants::*;
use crate::read_context::{build_context, ReadCigarInfo};
use crate::refseq::RefSequence;
use crate::types::*;
use crate::ultima::*;

/// new ReadCigarInfo(readAlignmentStart, cigar, flankPosStart, flankPosEnd, corePosStart, corePosEnd,
/// flankIndexStart, flankIndexEnd)
#[allow(clippy::too_many_arguments)]
fn read_cigar_info(
    read_alignment_start: i32, cigar: &str, flank_position_start: i32, flank_position_end: i32, core_position_start: i32,
    core_position_end: i32, flank_index_start: i32, flank_index_end: i32,
) -> ReadCigarInfo {
    ReadCigarInfo {
        cigar: parse_cigar(cigar),
        read_alignment_start,
        flank_position_start,
        flank_position_end,
        core_position_start,
        core_position_end,
        flank_index_start,
        flank_index_end,
    }
}

/// extendUltimaCore(readBases, refSequence, readAlignmentStart, cigarElements, readCigarInfo,
/// DEFAULT_FLANK_LENGTH, inAppendMode) over ref bases starting at position 100
fn extend(
    read_bases: &str, ref_bases: &str, read_alignment_start: i32, cigar: &str, info: &ReadCigarInfo, in_append_mode: bool,
) -> Option<ReadCigarInfo> {
    let ref_sequence = RefSequence::from_bases(100, ref_bases.as_bytes().to_vec());
    extend_ultima_core(
        read_bases.as_bytes(), &ref_sequence, read_alignment_start, &parse_cigar(cigar), info, DEFAULT_FLANK_LENGTH,
        in_append_mode,
    )
    .expect("no exception")
}

/// UltimaMiscUtilsTest.checkEqual
fn check_equal(first: &ReadCigarInfo, second: &ReadCigarInfo) -> bool {
    first.flank_index_start == second.flank_index_start
        && first.flank_index_end == second.flank_index_end
        && first.core_position_start == second.core_position_start
        && first.core_position_end == second.core_position_end
        && first.flank_position_start == second.flank_position_start
        && first.flank_position_end == second.flank_position_end
}

/// the assertions each extended case makes: flank index start, flank position start, core position start,
/// core position end, flank index end, flank position end, cigar
fn assert_info(info: &ReadCigarInfo, expected: (i32, i32, i32, i32, i32, i32, &str)) {
    assert_eq!(expected.0, info.flank_index_start);
    assert_eq!(expected.1, info.flank_position_start);
    assert_eq!(expected.2, info.core_position_start);
    assert_eq!(expected.3, info.core_position_end);
    assert_eq!(expected.4, info.flank_index_end);
    assert_eq!(expected.5, info.flank_position_end);
    assert_eq!(expected.6, cigar_str(&info.cigar));
}

// UltimaMiscUtilsTest.testCoreExtensionBasic
#[test]
fn test_core_extension_basic() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    // pos / index                  10             20             30        40
    //                    0123456789012345     67890     12345678901234567890123456789
    let mut ref_bases = "AACCGGTTAACCGGTTTACATTTAACCGGTTAACCGGTT".to_string();
    let mut read_bases = format!("{}G{}", &ref_bases[0..17], &ref_bases[18..]);

    let read_alignment_start = 100;
    let cigar = "39M";
    let info = read_cigar_info(read_alignment_start, cigar, 106, 130, 116, 120, 6, 30);

    let mut new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert_info(&new_info, (3, 103, 113, 123, 33, 133, "31M"));

    // extend further due to non-matching bases from SNVs
    // pos / index         10             20             30        40
    //           0123456789012345     67890     12345678901234567890123456789
    ref_bases = format!("{}{}{}", "AACCGGTTAACCAGGT", "TACAT", "TTTACCGGTTAACCGGTT");
    // SNV at start and longer HP, then a shorter HP and another SNV
    read_bases = format!("{}{}{}", "AACCGGTTAACCGTTT", "TAGAT", "TTAGCCGGTTAACCGGTT");

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert_info(&new_info, (1, 101, 111, 125, 35, 135, "35M"));
    set_test_sequencing_type(None);
}

// UltimaMiscUtilsTest.testCoreExtensionIndels
#[test]
fn test_core_extension_indels() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    // ref pos                    10                       20                      30        40
    //                  012345678901     I     2345     67890     12     3     45678901234567890123456789
    let mut ref_bases = format!("{}{}{}{}{}{}{}", "AACCGGTTAACC", "", "GGTT", "TACAT", "TT", "A", "ACCGGTTAACCGGTT");

    // read index                 10                       20                        30        40
    //                  012345678901     2     3456     78901     23     D     45678901234567890123456789
    // has inserted G before core and deleted A after core
    let mut read_bases = format!("{}{}{}{}{}{}{}", "AACCGGTTAACC", "G", "GGTT", "TAGAT", "TT", "", "ACCGGTTAACCGGTT");

    let read_alignment_start = 100;
    let mut cigar = "12M1I11M1D15M";
    let mut info = read_cigar_info(read_alignment_start, cigar, 107, 130, 116, 120, 17, 30);

    let mut new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert_info(&new_info, (4, 104, 113, 124, 34, 134, "8M1I11M1D11M"));

    // a deleted base immediately before and after the core
    // pos                 10                       20                  30
    //           012345678901234     5     67     89012     3     456789012345678
    ref_bases = format!("{}{}{}{}{}{}", "AACCGGTTAACCGGG", "A", "AA", "TAGAT", "A", "ACCGGTTAACCGGTT");

    //                     10                        20                   30
    // index     012345678901234     D     56     78901           234567890123456
    read_bases = format!("{}{}{}{}{}{}", "AACCGGTTAACCGGG", "", "AA", "TAGAT", "", "ACCGGTTAACCGGTT");
    // flank          FFFFFFFFFF                                  FFFFFFFFFF

    cigar = "15M1D7M1D15M";
    info = read_cigar_info(read_alignment_start, "8M1D7M1D10M", 107, 133, 118, 122, 7, 31);

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert_info(&new_info, (7, 107, 118, 122, 31, 133, "8M1D7M1D10M"));

    // an inserted base immediately before and after the core
    // pos                 10                  20        30
    //           012345678901234     56789     0123456789012345678
    ref_bases = format!("{}{}{}", "AACCGGTTAACCAGG", "ATCTA", "TCCGGTTAACCGGTT");

    //                     10                               20        30
    // pos       012345678901234           56789            01234567890123456
    //                     10                  20                   30
    // index     012345678901234     5     67890      1     234567890123456
    read_bases = format!("{}{}{}{}{}", "AACCGGTTAACCAGG", "T", "ATCTA", "T", "TCCGGTTAACCGGTT");
    // flank          FFFFFFFFFF                            FFFFFFFFFF

    cigar = "15M1I5M1I15M";
    info = read_cigar_info(read_alignment_start, "11M1I5M1I11M", 105, 129, 115, 119, 5, 31);

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert!(check_equal(&info, &new_info));
    set_test_sequencing_type(None);
}

// UltimaMiscUtilsTest.testCoreExtensionMiscExamples
#[test]
fn test_core_extension_misc_examples() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    // test 1: HP starts on ref but switches to read then ref again
    // pos / index                10             20             30        40
    //                  0123456789012345     67890     12345678901234567890123456789
    let mut ref_bases = format!("{}{}{}", "AACCGGTTAACCGGTA", "ACTCA", "ATAAACGGTTAACCGGTT");
    let mut read_bases = format!("{}{}{}", "AACCGGTTAACCGTTT", "ACTCA", "TTTACCGGTTAACCGGTT");
    //                        FFFFFFFFFF               FFFFFFFFFF

    let read_alignment_start = 100;
    let cigar = "39M";
    let info = read_cigar_info(read_alignment_start, "25M", 106, 130, 116, 120, 6, 30);

    let mut new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert_info(&new_info, (1, 101, 111, 127, 37, 137, "37M"));

    // test 2: no HP, but next bases don't match and so require extension
    // pos / index         10             20             30        40
    //           0123456789012345     67890     12345678901234567890123456789
    ref_bases = format!("{}{}{}", "AACCGGTTAACCGATC", "ACTCA", "CTAGACGGTTAACCGGTT");
    read_bases = format!("{}{}{}", "AACCGGTTAACCGATT", "GCTCG", "TTAGACGGTTAACCGGTT");
    //                 FFFFFFFFFF               FFFFFFFFFF

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert_info(&new_info, (4, 104, 114, 122, 32, 132, "29M"));

    // test 3: no HP and core end bases match
    // pos / index         10             20             30        40
    //           0123456789012345     67890     12345678901234567890123456789
    ref_bases = format!("{}{}{}", "AACCGGTTAACCGATC", "ACTCA", "CTAGACGGTTAACCGGTT");
    read_bases = format!("{}{}{}", "AACCGGTTAACCGATT", "ACGCA", "TTAGACGGTTAACCGGTT");
    //                 FFFFFFFFFF               FFFFFFFFFF

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false).expect("read cigar info");
    assert!(check_equal(&info, &new_info));
    set_test_sequencing_type(None);
}

// UltimaMiscUtilsTest.testCoreExtensionSoftClips
#[test]
fn test_core_extension_soft_clips() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    // test 1: left flank is in soft-clip and requires extension, making it invalid
    // pos / index                10             20             30
    //                  0123456789012345     67890     123456789012345678
    let mut ref_bases = format!("{}{}{}", "AACCGGTTAACCGGTT", "TACAT", "TTAACCGGTTAACCGGTT");
    let mut read_bases = format!("{}G{}", &ref_bases[0..18], &ref_bases[19..]);

    let mut read_alignment_start = 108;
    let mut cigar = "8S31M";
    let mut info = read_cigar_info(read_alignment_start, cigar, 106, 130, 116, 120, 6, 30);

    let mut new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false);
    assert!(new_info.is_none());

    // test 2: and now on the right in the core
    read_alignment_start = 100;
    info = read_cigar_info(read_alignment_start, cigar, 106, 130, 116, 120, 6, 30);
    cigar = "22M17S";

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false);
    assert!(new_info.is_none());

    // test 3: accept truncated flank in append mode if we extend from M into S
    cigar = "31M8S";

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, true);
    assert_info(&new_info.expect("read cigar info"), (3, 103, 113, 123, 30, 130, "28M"));

    // test 4: otherwise reject read context
    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false);
    assert!(new_info.is_none());

    // test 5: if core starts in soft-clip but does not require extension, still check/extend the other side
    ref_bases = format!("{}{}{}", "AACCGGTTAACCGGTT", "TACAT", "ATAACCGGTTAACCGGTT");
    read_bases = format!("{}G{}", &ref_bases[0..18], &ref_bases[19..]);
    cigar = "29M10S";

    new_info = extend(&read_bases, &ref_bases, read_alignment_start, cigar, &info, false);
    assert_info(&new_info.expect("read cigar info"), (3, 103, 113, 120, 30, 130, "26M2S"));
    set_test_sequencing_type(None);
}

// UltimaMiscUtilsTest.testBuildContextUltimaNoCoreExtension
#[test]
fn test_build_context_ultima_no_core_extension() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let middle_ref_bases = "CTCGGCCTCCCGGAGGTGCCGG";
    let middle_read_bases = "CTCGGCCTCCCTGAGGTGCCGG";
    let variant_middle_index = 11;

    let ref_padding_size = 20;
    let read_padding_size = 20;

    let ref_padding = "A".repeat(ref_padding_size);
    let read_padding = "A".repeat(read_padding_size);

    let ref_bases = format!("{ref_padding}{middle_ref_bases}{ref_padding}");
    let read_bases = format!("{read_padding}{middle_read_bases}{read_padding}");
    let variant = SimpleVariant::new((1 + ref_padding_size + variant_middle_index) as i32, b"G", b"T");

    let alignment_start = (1 + ref_padding_size - read_padding_size) as i32;
    let cigar = format!("{}M", read_bases.len());
    let read = create_sam_record_unpaired("READ_001", alignment_start, &read_bases, &cigar);
    let var_index_in_read = read_padding_size + variant_middle_index;

    let ref_sequence = RefSequence::from_bases(1, ref_bases.as_bytes().to_vec());

    let read_context =
        build_context(&variant, &read, var_index_in_read as i32, &ref_sequence, DEFAULT_FLANK_LENGTH).expect("read context");

    let expected_core_length = 7;

    assert_eq!(DEFAULT_FLANK_LENGTH, read_context.core_index_start);
    assert_eq!(DEFAULT_FLANK_LENGTH + expected_core_length - 1, read_context.core_index_end);
    assert_eq!(variant.pos - 4, read_context.core_position_start);
    assert_eq!(variant.pos + 2, read_context.core_position_end);
    assert_eq!("TCCCGGA", std::str::from_utf8(&read_context.ref_bases).unwrap());
    assert_eq!("TCCCTGA", read_context.core_str());
    let flank = DEFAULT_FLANK_LENGTH as usize;
    assert_eq!(
        &read_bases[var_index_in_read - 4 - flank..var_index_in_read + 3 + flank],
        std::str::from_utf8(&read_context.read_bases).unwrap()
    );
    set_test_sequencing_type(None);
}

// UltimaMiscUtilsTest.testBuildContextUltimaWithExtension
#[test]
fn test_build_context_ultima_with_extension() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    //                        0123456789012345678
    let middle_ref_bases = "CTTTTTTTTTTAATTGCAA";
    let middle_read_bases = "CTTTTTTTTTAATTGGCAA";
    let variant_middle_index = 14;

    let ref_padding_size = 20;
    let read_padding_size = 20;

    let ref_bases = format!("{}{}{}", "A".repeat(ref_padding_size), middle_ref_bases, "T".repeat(ref_padding_size));
    let read_bases = format!("{}{}{}", "A".repeat(read_padding_size), middle_read_bases, "T".repeat(read_padding_size));
    let variant = SimpleVariant::new((1 + ref_padding_size + variant_middle_index) as i32, b"T", b"G");

    let alignment_start = (1 + ref_padding_size - read_padding_size) as i32;
    let cigar = format!("{}M", read_bases.len());
    let read = create_sam_record_unpaired("READ_001", alignment_start, &read_bases, &cigar);
    let var_index_in_read = read_padding_size + variant_middle_index;

    let ref_sequence = RefSequence::from_bases(1, ref_bases.as_bytes().to_vec());

    let read_context =
        build_context(&variant, &read, var_index_in_read as i32, &ref_sequence, DEFAULT_FLANK_LENGTH).expect("read context");

    assert_eq!(10, read_context.core_index_start);
    assert_eq!(26, read_context.core_index_end);
    assert_eq!("CTTTTTTTTTAATTGGC", read_context.core_str());
    assert_eq!("CTTTTTTTTTTAATTGC", std::str::from_utf8(&read_context.ref_bases).unwrap());
    set_test_sequencing_type(None);
}

// UltimaMiscUtilsTest.testIsMsiIndelTypes
#[test]
fn test_is_msi_indel_types() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    // not an indel (Java passes null units)
    let mut variant = SimpleVariant::new(100, b"A", b"C");
    assert!(!is_msi_indel_of_type(&variant, &[]));

    // no repeat units
    variant = SimpleVariant::new(100, b"A", b"AA");
    assert!(!is_msi_indel_of_type(&variant, &[]));

    // single unit but not MSI
    variant = SimpleVariant::new(100, b"A", b"AAB");
    assert!(!is_msi_indel_of_type(&variant, &["A"]));

    // MSI repeat
    variant = SimpleVariant::new(100, b"A", b"AABABAB");
    assert!(is_msi_indel_of_type(&variant, &["AB"]));
    set_test_sequencing_type(None);
}

/// UltimaMiscUtilsTest.createRead
fn create_read(read_bases: &str) -> Read {
    build_sam_record(100, &format!("{}M", read_bases.len()), read_bases)
}

// UltimaMiscUtilsTest.testIsAdjacentToLongHomopolymers
#[test]
fn test_is_adjacent_to_long_homopolymers() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let long_length = 10usize;
    let ll = long_length as i32;

    // not long enough
    let mut read = create_read(&"A".repeat(2 * long_length - 1));
    assert!(!is_adjacent_to_long_homopolymer(&read, ll - 1, ll));

    // on left
    read = create_read(&format!("{}C{}", "A".repeat(long_length), "T".repeat(long_length - 1)));
    assert!(is_adjacent_to_long_homopolymer(&read, ll, ll));

    // on right
    read = create_read(&format!("{}C{}", "T".repeat(long_length - 1), "A".repeat(long_length)));
    assert!(is_adjacent_to_long_homopolymer(&read, ll - 1, ll));

    // not on left
    read = create_read(&format!("C{}", "A".repeat(2 * long_length - 1)));
    assert!(!is_adjacent_to_long_homopolymer(&read, ll, ll));

    // not on right
    read = create_read(&format!("{}C", "A".repeat(2 * long_length - 1)));
    assert!(!is_adjacent_to_long_homopolymer(&read, ll - 1, ll));

    // invalid bases
    read = create_read(&"N".repeat(2 * long_length));
    assert!(!is_adjacent_to_long_homopolymer(&read, ll, ll));
    set_test_sequencing_type(None);
}
