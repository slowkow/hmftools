//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! common/VariantReadContextTest.java

use super::util::*;
use crate::candidate::soft_clip_insert;
use crate::constants::{BASE_QUAL_MINIMUM, DEFAULT_FLANK_LENGTH};
use crate::matcher::{MatchType, Matcher};
use crate::read_context::{create_context, find_left_homology_shift, has_indel_in_core, VariantReadContext};
use crate::refseq::RefSequence;
use crate::types::*;

const TEST_FLANK_LENGTH: i32 = 5;

fn create_simple_variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

fn ref_sequence(start: i32, bases: &str) -> RefSequence {
    RefSequence::from_bases(start, bases.as_bytes().to_vec())
}

fn s(bases: &[u8]) -> &str {
    std::str::from_utf8(bases).unwrap()
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

/// VariantReadContext.homologyBases
fn homology_bases(read_context: &VariantReadContext) -> &str {
    read_context.homology.as_ref().map_or("", |h| s(&h.bases))
}

// VariantReadContextTest.testSnvMnvVariantReadContext
#[test]
fn test_snv_mnv_variant_read_context() {
    let ref_seq_200 = ref_sequence_200();

    // test 1: most simple case with no repeats or homology
    let var = create_simple_variant(50, "C", "A");

    let read_bases = format!("{}A{}", &REF_BASES_200[30..50], &REF_BASES_200[51..70]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "40M";
    let read = build_sam_record_quals(30, read_cigar, &read_bases, base_quals);

    let flank = TEST_FLANK_LENGTH;

    let read_context = create_context(&var, &read, 20, &ref_seq_200, flank).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(43, read_context.alignment_start);
    assert_eq!(57, read_context.alignment_end);
    assert_eq!(5, read_context.core_index_start);
    assert_eq!(7, read_context.var_index);
    assert_eq!(9, read_context.core_index_end);
    assert_eq!("CGCTG", s(&read_context.ref_bases));
    assert_eq!("GTCACCGATGTCTGT", s(&read_context.read_bases));
    assert_eq!("15M", read_context.read_cigar);
    assert_eq!("GCT", s(read_context.trinucleotide()));
    assert_eq!(5, read_context.core_length());
    assert_eq!(5, read_context.left_flank_length());
    assert_eq!(5, read_context.right_flank_length());
    assert_eq!("CGATG", read_context.core_str());
    assert_eq!("GTCAC", read_context.left_flank_str());
    assert_eq!("TCTGT", read_context.right_flank_str());
    assert!(read_context.homology.is_none());
    assert!(read_context.max_repeat.is_none());

    // test 2:
    // Read bases and cigar:
    // 15M index 0-14, pos 30-44
    // 3I index 15-17, from pos 44
    // 10M index 18-27, pos 45-54, with SNV at index 23
    // 3I index 28-30, from pos 54
    // 15M index 31-45, pos 55-69
    // Initially flank indices are 16-30 of the read bases, but since 16 falls into the insert, it will push back to 14 and 1M at 44
    // and likewise since 30 falls into the last insert it will push back to 31 and pos 55

    let read_bases = format!(
        "{}GGG{}A{}GGG{}",
        &REF_BASES_200[30..45],
        &REF_BASES_200[45..50],
        &REF_BASES_200[51..55],
        &REF_BASES_200[55..70]
    );

    let read_cigar = "15M3I10M3I15M";
    let read = build_sam_record(30, read_cigar, &read_bases);

    let read_context = create_context(&var, &read, 23, &ref_seq_200, flank).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(44, read_context.alignment_start);
    assert_eq!(55, read_context.alignment_end);
    assert_eq!(7, read_context.core_index_start);
    assert_eq!(9, read_context.var_index);
    assert_eq!(11, read_context.core_index_end);
    assert_eq!("CGCTG", s(&read_context.ref_bases));
    assert_eq!("TGGGCACCGATGTCGGGT", s(&read_context.read_bases));
    assert_eq!("1M3I10M3I1M", read_context.read_cigar);
    assert_eq!(5, read_context.core_length());
    assert_eq!(7, read_context.left_flank_length());
    assert_eq!(6, read_context.right_flank_length());
    assert_eq!("CGATG", read_context.core_str());
    assert_eq!("TGGGCAC", read_context.left_flank_str());
    assert_eq!("TCGGGT", read_context.right_flank_str());

    // test 3: with transitional repeats near the core

    // TGCGCGC TA CACACACT -> TGCGCGC GC CACACACT
    let var = create_simple_variant(143, "TA", "GC");

    let read_bases = format!("{}GC{}", &REF_BASES_200[120..143], &REF_BASES_200[145..170]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "50M";
    let read = build_sam_record_quals(120, read_cigar, &read_bases, base_quals);

    let read_context = create_context(&var, &read, 23, &ref_seq_200, flank).unwrap();

    assert!(read_context.is_valid());
    assert!(read_context.max_repeat.is_some());
    let max_repeat = read_context.max_repeat.as_ref().unwrap();
    assert_eq!("GC", max_repeat.bases_str());
    assert_eq!(4, max_repeat.count);
    assert_eq!(131, read_context.alignment_start);
    assert_eq!(157, read_context.alignment_end);
    assert_eq!(5, read_context.core_index_start);
    assert_eq!(12, read_context.var_index);
    assert_eq!(21, read_context.core_index_end);
    assert_eq!("TGCGCGCTACACACACT", s(&read_context.ref_bases));
    assert_eq!("TGCGCGCGCCACACACT", read_context.core_str());
    assert_eq!("27M", read_context.read_cigar);
    assert_eq!(17, read_context.core_length());
    assert_eq!(5, read_context.left_flank_length());
    assert_eq!(5, read_context.right_flank_length());
}

// VariantReadContextTest.testLowBaseQualVariants
#[test]
fn test_low_base_qual_variants() {
    let ref_seq_200 = ref_sequence_200();

    // test 1: most simple case with no repeats or homology
    let var = create_simple_variant(50, "C", "A");

    let read_bases = format!("{}A{}", &REF_BASES_200[30..50], &REF_BASES_200[51..70]);
    let mut base_quals = default_quals(read_bases.len());

    // low qual bases in core
    base_quals[18] = BASE_QUAL_MINIMUM;
    let read_cigar = "40M";
    let read = build_sam_record_quals(30, read_cigar, &read_bases, base_quals);

    let read_context = create_context(&var, &read, 20, &ref_seq_200, DEFAULT_FLANK_LENGTH);
    assert!(read_context.is_none());

    // an indel with low bases in the flank

    let var = create_simple_variant(50, "C", "CTTT");
    let read_bases = format!("{}TTT{}", &REF_BASES_200[30..51], &REF_BASES_200[51..70]);
    let mut base_quals = default_quals(read_bases.len());

    // low qual bases in core
    base_quals[6] = BASE_QUAL_MINIMUM;
    base_quals[31] = BASE_QUAL_MINIMUM;
    let read_cigar = "20M3I20M";
    let read = build_sam_record_quals(30, read_cigar, &read_bases, base_quals);

    let read_context = create_context(&var, &read, 20, &ref_seq_200, DEFAULT_FLANK_LENGTH);
    assert!(read_context.is_none());
}

// VariantReadContextTest.testSnvComplex
#[test]
fn test_snv_complex() {
    let ref_bases = format!(
        "{}{}{}",
        &REF_BASES_200[0..50],
        //       50        60        70
        //       0        90123456789012345        4
        //                             T
        "CTTTTTTTTTATTTTTTTTTTTCTTTTTTTTGAGA",
        &REF_BASES_200[100..150]
    );

    let ref_seq = ref_sequence(0, &ref_bases);

    let read_bases = format!("{}T{}", &ref_bases[40..72], &ref_bases[73..100]);

    let read_cigar = "60M";
    let read = build_sam_record(40, read_cigar, &read_bases);

    let var = create_simple_variant(72, "C", "T");
    let read_context = create_context(&var, &read, 32, &ref_seq, DEFAULT_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(50, read_context.alignment_start);
    assert_eq!(91, read_context.alignment_end);
    assert_eq!(10, read_context.core_index_start);
    assert_eq!(22, read_context.var_index);
    assert_eq!(31, read_context.core_index_end);
    assert_eq!(12, read_context.variant_ref_index());
    assert_eq!("ATTTTTTTTTTTCTTTTTTTTG", s(&read_context.ref_bases));
    assert_eq!("CTTTTTTTTTATTTTTTTTTTTTTTTTTTTTGAGAGAGGCAA", s(&read_context.read_bases));

    // test 1: a core match
    let read_ref_bases = &ref_bases[40..100];
    let ref_read = build_sam_record(40, read_cigar, read_ref_bases);

    // VariantUtils.createReadContextMatcher
    let matcher = Matcher::new(&read_context, true, false);
    assert_eq!(MatchType::Ref, matcher.determine_read_match(&read_context, &ref_read, 32).kind);

    // test 2: repeats starting 1 base outside the basic core
    let ref_bases = "AGGTCGATCGTAAAAAAATGCGTCGCGCGCGCGCGTAGGTCGATCG";
    //               100       110       120       130       140
    //               01234567890123456789012345678901234567890123456789
    //                                   A
    let ref_seq = ref_sequence(100, ref_bases);

    let alt_base = "A";
    let read_bases = format!("{}A{}", &ref_bases[0..20], &ref_bases[21..]);
    let read = build_sam_record(100, &build_cigar_string(read_bases.len()), &read_bases);

    let var = create_simple_variant(120, &ref_bases[20..21], alt_base);
    let read_context = create_context(&var, &read, 20, &ref_seq, DEFAULT_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(100, read_context.alignment_start);
    assert_eq!(145, read_context.alignment_end);
    assert_eq!(10, read_context.core_index_start);
    assert_eq!(20, read_context.var_index);
    assert_eq!(35, read_context.core_index_end);
    assert_eq!(2, read_context.all_repeats.len());
    assert_eq!("TAAAAAAATGAGTCGCGCGCGCGCGT", read_context.core_str());
}

// VariantReadContextTest.testVariantRepeatContexts
#[test]
fn test_variant_repeat_contexts() {
    let ref_bases =
        //             10        20        30        40        50
        //   0123456789012345678901234567890123456789012345678901
            "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

    let ref_seq = ref_sequence(0, ref_bases);

    // test 1: a repeat starting at the variant
    let read_bases = format!("{}C{}", &ref_bases[0..21], &ref_bases[22..51]);

    let read_cigar = build_cigar_string(read_bases.len());
    let read = build_sam_record(1, &read_cigar, &read_bases);

    let var = create_simple_variant(21, "A", "C");

    let read_context = create_context(&var, &read, 21, &ref_seq, 5).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(5, read_context.core_index_start);
    assert!(read_context.max_repeat.is_some());
    let max_repeat = read_context.max_repeat.as_ref().unwrap();
    assert_eq!("C", max_repeat.bases_str());
    assert_eq!(4, max_repeat.count);
    assert_eq!(7, read_context.var_index);
    assert_eq!(11, read_context.core_index_end);

    // test 2 : repeats before and after but only the first is used
    let var = create_simple_variant(30, "A", "C");

    let read_bases = format!("{}{}{}", &ref_bases[0..30], var.alt_str(), &ref_bases[31..51]);

    let read_cigar = build_cigar_string(read_bases.len());
    let read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = create_context(&var, &read, 30, &ref_seq, 5).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(5, read_context.core_index_start);
    assert!(read_context.max_repeat.is_some());
    let max_repeat = read_context.max_repeat.as_ref().unwrap();
    assert_eq!("C", max_repeat.bases_str());
    assert_eq!(4, max_repeat.count);
    assert_eq!(9, read_context.var_index);
    assert_eq!(13, read_context.core_index_end);
}

// VariantReadContextTest.testDeleteVariantReadContext
#[test]
fn test_delete_variant_read_context() {
    let ref_seq_200 = ref_sequence_200();

    // firstly a simple delete will no repeats or homology
    let var = create_simple_variant(50, "CTG", "C");

    let read_bases = format!("{}{}", &REF_BASES_200[30..51], &REF_BASES_200[53..70]);
    let read_cigar = "21M2D17M";
    let read = build_sam_record(30, read_cigar, &read_bases);

    let read_context = create_context(&var, &read, 20, &ref_seq_200, TEST_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(44, read_context.alignment_start);
    assert_eq!(60, read_context.alignment_end);
    assert_eq!(5, read_context.core_index_start);
    assert_eq!(6, read_context.var_index);
    assert_eq!(9, read_context.core_index_end);
    assert_eq!("GCTGTCT", s(&read_context.ref_bases));
    assert_eq!("TCACCGCTCTGTGAC", s(&read_context.read_bases));
    assert_eq!("7M2D8M", read_context.read_cigar);
    assert_eq!(5, read_context.core_length());
    assert_eq!(5, read_context.left_flank_length());
    assert_eq!(5, read_context.right_flank_length());
    assert_eq!("GCTCT", read_context.core_str());
    assert_eq!("TCACC", read_context.left_flank_str());
    assert_eq!("GTGAC", read_context.right_flank_str());
    assert_eq!(49, read_context.core_position_start);
    assert_eq!(55, read_context.core_position_end);

    // test 2: read context spans into the soft-clip bases on the left
    let ref_bases = &REF_BASES_200[57..59];
    let alt_base = &ref_bases[0..1];
    let var = create_simple_variant(57, ref_bases, alt_base);

    let read_bases = format!("{}{}", &REF_BASES_200[1..58], &REF_BASES_200[59..79]);

    // unclipped pos 1-49 = 49S, read index 0-48
    // aligned pos   50-57, read index 49-56
    // del pos       58
    // aligned pos   59-78, read index 57-76
    // so del occurs 49 + 8 = 57 bases into the read
    let read_cigar = "49S8M1D20M";
    let read = build_sam_record(50, read_cigar, &read_bases);

    // full flank so it extends into the soft-clip
    let read_context = create_context(&var, &read, 56, &ref_seq_200, DEFAULT_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(50, read_context.alignment_start);
    assert_eq!(70, read_context.alignment_end);
    assert_eq!(10, read_context.core_index_start);
    assert_eq!(11, read_context.var_index);
    assert_eq!(13, read_context.core_index_end);
    assert_eq!("4S8M1D12M", read_context.read_cigar);
    assert_eq!("GTAC", read_context.core_str());
    assert_eq!("ACCGCTGTCT", read_context.left_flank_str());
    assert_eq!("TCGGAAAAAA", read_context.right_flank_str());

    // test 3: del spanning into the right soft-clip
    let ref_bases = &REF_BASES_200[30..33];
    let alt_base = &ref_bases[0..1];
    let var = create_simple_variant(30, ref_bases, alt_base);

    let read_bases = format!("{}{}", &REF_BASES_200[1..31], &REF_BASES_200[33..53]);
    let read_cigar = "30M2D5M15S";
    let read = build_sam_record(1, read_cigar, &read_bases);

    // full flank so it extends into the soft-clip
    let read_context = create_context(&var, &read, 29, &ref_seq_200, DEFAULT_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(16, read_context.alignment_start);
    assert_eq!(37, read_context.alignment_end);
    assert_eq!(10, read_context.core_index_start);
    assert_eq!(14, read_context.var_index);
    assert_eq!(16, read_context.core_index_end);
    assert_eq!("15M2D5M7S", read_context.read_cigar);
}

// VariantReadContextTest.testInsertVariantReadContext
#[test]
fn test_insert_variant_read_context() {
    let ref_seq_200 = ref_sequence_200();

    let var = create_simple_variant(50, "C", "CTG");

    let read_bases = format!("{}TG{}", &REF_BASES_200[30..51], &REF_BASES_200[51..70]);
    let read_cigar = "21M2I17M";
    let read = build_sam_record(30, read_cigar, &read_bases);

    let read_context = create_context(&var, &read, 20, &ref_seq_200, TEST_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(44, read_context.alignment_start);
    assert_eq!(60, read_context.alignment_end);
    assert_eq!(5, read_context.core_index_start);
    assert_eq!(6, read_context.var_index);
    assert_eq!(13, read_context.core_index_end);
    assert_eq!("GCTGTCT", s(&read_context.ref_bases));
    assert_eq!("TCACCGCTGTGTCTGTGAC", s(&read_context.read_bases));
    assert_eq!("7M2I10M", read_context.read_cigar);
    assert_eq!(9, read_context.core_length());
    assert_eq!(5, read_context.left_flank_length());
    assert_eq!(5, read_context.right_flank_length());
    assert_eq!("TG", homology_bases(&read_context));
    assert_eq!(3, read_context.homology.as_ref().unwrap().length);
    assert_eq!("GCTGTGTCT", read_context.core_str());
    assert_eq!("TCACC", read_context.left_flank_str());
    assert_eq!("GTGAC", read_context.right_flank_str());

    assert_eq!(49, read_context.core_position_start);
    assert_eq!(55, read_context.core_position_end);

    // test 2: insert in left soft-clip with homology

    // a duplicated section leading up to the aligned bases
    // TTCGTCACCGCTGTCTGTGACTCGGAAAAAAAACTCC

    let insert = &REF_BASES_200[40..55];
    let read_bases = format!("{}{}{}", &REF_BASES_200[15..40], insert, &REF_BASES_200[40..110]);
    let read_cigar = "40S70M";
    let read = build_sam_record(40, read_cigar, &read_bases);

    // RefContextConsumer.processSoftClip (static) returns the AltRead's ref and alt
    let (alt_ref, alt_alt) = soft_clip_insert(40, 109, read_bases.as_bytes(), 40, 0, &ref_seq_200, true).unwrap();

    let var = SimpleVariant::new(39, &alt_ref, &alt_alt);

    let read_context = create_context(&var, &read, 39, &ref_seq_200, DEFAULT_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(38, read_context.core_position_start);
    assert_eq!(11, read_context.var_index);
    assert_eq!(58, read_context.core_position_end);
    assert_eq!(16, read_context.homology.as_ref().unwrap().length);

    // test 3: insertion (a duplicated segment) in the right soft clip
    let insert = &REF_BASES_200[85..100];
    let read_bases = format!("{}{}{}", &REF_BASES_200[30..100], insert, &REF_BASES_200[100..125]);
    let read_cigar = "70M40S";
    let read = build_sam_record(30, read_cigar, &read_bases);

    let (alt_ref, alt_alt) = soft_clip_insert(30, 99, read_bases.as_bytes(), 40, 70, &ref_seq_200, false).unwrap();

    let mut read_index = read_bases.len() as i32 - 40 - 1;

    let var = SimpleVariant::new(99, &alt_ref, &alt_alt);

    let left_homology_shift = find_left_homology_shift(&var, &ref_seq_200, read_bases.as_bytes(), read_index);

    read_index -= left_homology_shift;
    let new_position = var.pos - left_homology_shift;
    let new_ref = s(ref_seq_200.position_bases(new_position, new_position).unwrap()).to_string();
    let new_alt =
        format!("{}{}", new_ref, &read_bases[(read_index + 1) as usize..(read_index + var.alt_len()) as usize]);
    let var = create_simple_variant(new_position, &new_ref, &new_alt);

    let read_context = create_context(&var, &read, read_index, &ref_seq_200, DEFAULT_FLANK_LENGTH).unwrap();

    assert!(read_context.is_valid());
    assert_eq!(80, read_context.core_position_start);
    assert_eq!(14, read_context.var_index);
    assert_eq!(116, read_context.core_position_end);
}

// VariantReadContextTest.testIndelInCoreCigarElements
#[test]
fn test_indel_in_core_cigar_elements() {
    let cigar_elements = vec![CigarEl::new(10, Op::M)];

    assert!(!has_indel_in_core(&cigar_elements, 0, 9));

    let cigar_elements = vec![CigarEl::new(5, Op::M), CigarEl::new(1, Op::I), CigarEl::new(5, Op::M)];

    assert!(!has_indel_in_core(&cigar_elements, 0, 4));
    assert!(!has_indel_in_core(&cigar_elements, 6, 10));
    assert!(has_indel_in_core(&cigar_elements, 0, 5));
    assert!(has_indel_in_core(&cigar_elements, 5, 10));

    let cigar_elements = vec![CigarEl::new(5, Op::M), CigarEl::new(1, Op::D), CigarEl::new(5, Op::M)];

    assert!(has_indel_in_core(&cigar_elements, 0, 9));

    let cigar_elements = vec![CigarEl::new(20, Op::M), CigarEl::new(1, Op::D), CigarEl::new(20, Op::M)];

    assert!(!has_indel_in_core(&cigar_elements, 5, 15));
    assert!(!has_indel_in_core(&cigar_elements, 25, 30));
    assert!(has_indel_in_core(&cigar_elements, 19, 20));
}

// VariantReadContextTest.testCoreCigarElements
#[test]
fn test_core_cigar_elements() {
    let core_length = 10;
    check_core_cigar_elements(
        core_length,
        vec![CigarEl::new(2 * TEST_FLANK_LENGTH + core_length, Op::M)],
        vec![CigarEl::new(core_length, Op::M)],
    );

    check_core_cigar_elements(
        core_length,
        vec![
            CigarEl::new(TEST_FLANK_LENGTH, Op::M),
            CigarEl::new(1, Op::D),
            CigarEl::new(core_length + TEST_FLANK_LENGTH, Op::M),
        ],
        vec![CigarEl::new(core_length, Op::M)],
    );

    check_core_cigar_elements(
        core_length,
        vec![
            CigarEl::new(TEST_FLANK_LENGTH + core_length, Op::M),
            CigarEl::new(1, Op::D),
            CigarEl::new(TEST_FLANK_LENGTH, Op::M),
        ],
        vec![CigarEl::new(core_length, Op::M)],
    );

    check_core_cigar_elements(
        core_length,
        vec![
            CigarEl::new(TEST_FLANK_LENGTH, Op::M),
            CigarEl::new(1, Op::I),
            CigarEl::new(core_length - 1 + TEST_FLANK_LENGTH, Op::M),
        ],
        vec![CigarEl::new(1, Op::I), CigarEl::new(core_length - 1, Op::M)],
    );

    check_core_cigar_elements(
        core_length,
        vec![
            CigarEl::new(TEST_FLANK_LENGTH + core_length - 1, Op::M),
            CigarEl::new(1, Op::I),
            CigarEl::new(TEST_FLANK_LENGTH, Op::M),
        ],
        vec![CigarEl::new(core_length - 1, Op::M), CigarEl::new(1, Op::I)],
    );

    check_core_cigar_elements(
        core_length,
        vec![
            CigarEl::new(TEST_FLANK_LENGTH + 1, Op::M),
            CigarEl::new(1, Op::D),
            CigarEl::new(core_length - 1 + TEST_FLANK_LENGTH, Op::M),
        ],
        vec![CigarEl::new(1, Op::M), CigarEl::new(1, Op::D), CigarEl::new(core_length - 1, Op::M)],
    );
}

/// VariantReadContextTest.checkCoreCigarElements. The Java constructor takes nulls and -1 for everything
/// but the cigar and core indices; the port's struct needs values, so placeholders are used.
fn check_core_cigar_elements(core_length: i32, read_cigar: Vec<CigarEl>, expected_core_cigar: Vec<CigarEl>) {
    let core_index_start = TEST_FLANK_LENGTH;
    let core_index_end = core_index_start + core_length - 1;

    let read_context = VariantReadContext {
        variant: SimpleVariant::new(-1, b"", b""),
        alignment_start: -1,
        alignment_end: -1,
        ref_bases: Vec::new(),
        read_bases: Vec::new(),
        core_index_start,
        var_index: -1,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start: -1,
        core_position_end: -1,
        read_cigar: cigar_str(&read_cigar),
        has_indel_in_core: has_indel_in_core(&read_cigar, core_index_start, core_index_end),
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    };

    let actual_core_cigar = read_context_core_cigar(&read_context);
    assert_eq!(expected_core_cigar, actual_core_cigar);
}

/// VariantReadContextTest.readContextCoreCigar
fn read_context_core_cigar(read_context: &VariantReadContext) -> Vec<CigarEl> {
    let mut core_cigar = Vec::new();
    let cigar = parse_cigar(&read_context.read_cigar);
    let mut read_index = 0;
    for el in cigar {
        if read_index > read_context.core_index_end {
            break;
        }

        if !el.op.consumes_read() {
            if read_index > read_context.core_index_start {
                core_cigar.push(el);
            }
            continue;
        }

        let mut read_index_end = read_index + el.len - 1;
        if read_index_end < read_context.core_index_start {
            read_index += el.len;
            continue;
        }

        let read_index_start = read_index.max(read_context.core_index_start);
        read_index_end = read_index_end.min(read_context.core_index_end);
        core_cigar.push(CigarEl::new(read_index_end - read_index_start + 1, el.op));

        read_index += el.len;
    }

    core_cigar
}
