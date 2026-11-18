//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/common/ReadContextMatchingTest.java
//!
//! Java `matcher.determineReadMatch(read, index)` is `Matcher::determine_read_match(ctx, read, index).kind`
//! (the Rust matcher does not own its read context), and `matcher.coversVariant(read, index)` is
//! `Matcher::covers_variant(ctx, read_len, index)`.

use super::util::*;
use crate::constants::DEFAULT_FLANK_LENGTH;
use crate::evidence::average_core_quality;
use crate::matcher::{MatchType, Matcher};
use crate::read_context::{create_context, has_indel_in_core, VariantReadContext};
use crate::types::*;

// VariantUtils.TEST_LEFT_FLANK etc
const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";
const TEST_LEFT_CORE: &str = "AC";
const TEST_RIGHT_CORE: &str = "GT";

/// VariantUtils.createSimpleVariant(position, ref, alt)
fn create_simple_variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

/// VariantUtils.createReadContext(variant, leftCore, rightCore, leftFlank, rightFlank): a context with no
/// repeat or homology
fn create_read_context_flanks(
    variant: &SimpleVariant, left_core: &str, right_core: &str, left_flank: &str, right_flank: &str,
) -> VariantReadContext {
    let core_index_start = left_flank.len() as i32;
    let var_read_index = core_index_start + left_core.len() as i32;
    let core_index_end = var_read_index + variant.alt_len() - 1 + right_core.len() as i32;
    let ref_bases = format!("{}{}{}", left_core, variant.ref_str(), right_core);
    let read_bases = format!("{}{}{}{}{}", left_flank, left_core, variant.alt_str(), right_core, right_flank);

    let alignment_start = variant.pos - var_read_index;

    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + right_core.len() as i32 + right_flank.len() as i32 - 1;

    let read_cigar = vec![CigarEl::new(read_bases.len() as i32, Op::M)];

    let core_position_start = variant.pos - left_core.len() as i32;
    let core_position_end = variant.pos + right_core.len() as i32;

    VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.into_bytes(),
        read_bases: read_bases.into_bytes(),
        core_index_start,
        var_index: var_read_index,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start,
        core_position_end,
        read_cigar: cigar_str(&read_cigar),
        has_indel_in_core: has_indel_in_core(&read_cigar, core_index_start, core_index_end),
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    }
}

/// VariantUtils.createReadContext(variant, leftCore, rightCore)
fn create_read_context_cores(variant: &SimpleVariant, left_core: &str, right_core: &str) -> VariantReadContext {
    create_read_context_flanks(variant, left_core, right_core, TEST_LEFT_FLANK, TEST_RIGHT_FLANK)
}

/// VariantUtils.createReadContext(variant)
fn create_read_context(variant: &SimpleVariant) -> VariantReadContext {
    create_read_context_cores(variant, TEST_LEFT_CORE, TEST_RIGHT_CORE)
}

/// VariantUtils.createReadContextMatcher
fn create_read_context_matcher(read_context: &VariantReadContext) -> Matcher {
    Matcher::new(read_context, true, false)
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

/// ReadContextMatcher.determineReadMatch(read, readVarIndex)
fn determine_read_match(matcher: &Matcher, read_context: &VariantReadContext, read: &Read, read_var_index: i32) -> MatchType {
    matcher.determine_read_match(read_context, read, read_var_index).kind
}

// Java: testReadCoversCore
#[test]
fn test_read_covers_core() {
    let position = 100;
    let variant = create_simple_variant(position, "A", "T");

    let mut read_context = create_read_context(&variant);

    let mut matcher = create_read_context_matcher(&read_context);

    let read_bases = &REF_BASES_200[0..40];
    let read_qualities = default_quals(read_bases.len());
    let cigar = build_cigar_string(read_bases.len());

    let read = build_sam_record_quals(1, &cigar, read_bases, read_qualities);

    assert!(matcher.covers_variant(&read_context, read.len(), 20));
    assert!(matcher.covers_variant(&read_context, read.len(), 0));
    assert!(matcher.covers_variant(&read_context, read.len(), 39));

    // now test with an indel with a wider alt range
    let mut var = create_simple_variant(position, "C", "CCC");
    read_context = create_read_context_cores(&var, "ATGTGTG", "CCA");
    matcher = create_read_context_matcher(&read_context);

    assert!(matcher.covers_variant(&read_context, read.len(), 20));
    assert!(matcher.covers_variant(&read_context, read.len(), 0));
    assert!(matcher.covers_variant(&read_context, read.len(), 36));
    assert!(!matcher.covers_variant(&read_context, read.len(), 37));

    var = create_simple_variant(position, "ACGT", "A");
    read_context = create_read_context_cores(&var, "AA", "GG");
    matcher = create_read_context_matcher(&read_context);

    assert!(matcher.covers_variant(&read_context, read.len(), 20));
    assert!(matcher.covers_variant(&read_context, read.len(), 0));
    assert!(matcher.covers_variant(&read_context, read.len(), 38));
    assert!(!matcher.covers_variant(&read_context, read.len(), 39));
}

// Java: testRefAndCoreMatches
#[test]
fn test_ref_and_core_matches() {
    let left_core = "AC";
    let right_core = "GT";
    let position = 100;
    let ref_ = "A";
    let alt = "T";
    let variant = create_simple_variant(position, ref_, alt);
    let read_context = create_read_context_cores(&variant, left_core, right_core);
    let left_flank = read_context.left_flank_str();
    let right_flank = read_context.right_flank_str();

    let matcher = create_read_context_matcher(&read_context);

    let mut read_bases = format!("{left_flank}{left_core}{ref_}{right_core}{right_flank}");
    let mut read_qualities = default_quals(read_bases.len());
    let mut read_var_index = read_context.left_flank_length() + 2;
    let cigar = build_cigar_string(read_bases.len());

    let mut read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::Ref, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // partial ref core assuming satisfies core-covered
    read_bases = format!("{ref_}{right_core}{right_flank}");
    read = build_sam_record_quals(position, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::Ref, determine_read_match(&matcher, &read_context, &read, 0));

    // same again but ending at the ref
    read_bases = format!("ACGTA{left_flank}{left_core}{ref_}");
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::Ref, determine_read_match(&matcher, &read_context, &read, read_bases.len() as i32 - 1));

    // full match
    read_bases = format!("{left_flank}{left_core}{alt}{right_core}{right_flank}");
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::Full, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // low-qual mismatch but not at variant
    read_bases = format!("{left_flank}AAAGT{right_flank}");
    read_qualities[(read_var_index - 1) as usize] = 10;
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::Ref, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // too many mismatches
    read_bases = format!("{left_flank}AAACT{right_flank}");
    read_qualities[(read_var_index - 1) as usize] = 10;
    read_qualities[(read_var_index + 1) as usize] = 10;
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, read_var_index));
    assert_ne!(MatchType::Ref, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // mismatch at SNV itself cannot be low-qual
    read_bases = format!("{left_flank}ACGGT{right_flank}");
    read_qualities = default_quals(read_bases.len());
    read_qualities[read_var_index as usize] = 10;
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // can be low-qual at the variant base itself
    read_bases = format!("{left_flank}{left_core}{ref_}{right_core}{right_flank}");
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::Ref, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // partial core
    read_bases = format!("AAAAA{left_flank}{left_core}{alt}");
    read_var_index = read_bases.len() as i32 - 1;
    read_qualities = default_quals(read_bases.len());
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::PartialCore, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // same again from the other side
    read_bases = format!("{alt}{right_core}{right_flank}AAAAA");
    read_var_index = 0;
    read_qualities = default_quals(read_bases.len());
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());
    assert_eq!(MatchType::PartialCore, determine_read_match(&matcher, &read_context, &read, read_var_index));
}

// Java: testSimpleAltMatches
#[test]
fn test_simple_alt_matches() {
    let left_core = "AC";
    let right_core = "GT";
    let position = 100;
    let mut ref_ = "AGT";
    let mut alt = "TGC";
    let mut variant = create_simple_variant(position, ref_, alt);
    let mut read_context = create_read_context_cores(&variant, left_core, right_core);

    let mut matcher = Matcher::new(&read_context, true, true);

    let mut read_bases =
        format!("{}TT{}{}{}", read_context.left_flank_str(), alt, right_core, read_context.right_flank_str());
    let mut read_qualities = default_quals(read_bases.len());
    let mut read_var_index = read_context.left_flank_length() + 2;
    let mut cigar = build_cigar_string(read_bases.len());

    let mut read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::SimpleAlt, determine_read_match(&matcher, &read_context, &read, read_var_index));

    read_bases = format!("{}TTTCC{}{}", read_context.left_flank_str(), right_core, read_context.right_flank_str());
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // repeated for an SNV
    ref_ = "A";
    alt = "T";
    variant = create_simple_variant(position, ref_, alt);
    read_context = create_read_context_cores(&variant, left_core, right_core);

    matcher = Matcher::new(&read_context, true, true);

    read_bases = format!("{}{}{}CC{}", read_context.left_flank_str(), left_core, alt, read_context.right_flank_str());
    read_qualities = default_quals(read_bases.len());
    read_var_index = read_context.left_flank_length() + 2;
    cigar = build_cigar_string(read_bases.len());

    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::SimpleAlt, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // test 3: delete
    ref_ = "ATG";
    alt = "A";
    variant = create_simple_variant(position, ref_, alt);
    read_context = create_read_context_cores(&variant, left_core, right_core);

    matcher = Matcher::new(&read_context, true, true);

    read_bases = format!("{}{}{}CC{}", read_context.left_flank_str(), left_core, alt, read_context.right_flank_str());
    read_qualities = default_quals(read_bases.len());
    read_var_index = read_context.left_flank_length() + 2;

    cigar = "13M2D12M".to_string();
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::SimpleAlt, determine_read_match(&matcher, &read_context, &read, read_var_index));

    cigar = "14M2D12M".to_string();
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // test 4: insert
    ref_ = "A";
    alt = "ATTT";
    variant = create_simple_variant(position, ref_, alt);
    read_context = create_read_context_cores(&variant, left_core, right_core);

    matcher = Matcher::new(&read_context, true, true);

    read_bases = format!("{}{}{}CC{}", read_context.left_flank_str(), left_core, alt, read_context.right_flank_str());
    read_qualities = default_quals(read_bases.len());
    read_var_index = read_context.left_flank_length() + 2;
    cigar = "13M3I12M".to_string();

    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::SimpleAlt, determine_read_match(&matcher, &read_context, &read, read_var_index));

    cigar = "13M4I12M".to_string();
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, read_var_index));
}

// Java: testIndelMatches
#[test]
fn test_indel_matches() {
    // basic del ref match
    let position = 50usize;
    let var = create_simple_variant(
        position as i32,
        &REF_BASES_200[position..position + 5],
        &REF_BASES_200[position..position + 1],
    );

    let read_bases = format!("{}{}", &REF_BASES_200[30..position + 1], &REF_BASES_200[position + 5..70]);
    let read_cigar = "21M2D17M";
    let read = build_sam_record(30, read_cigar, &read_bases);

    let ref_sequence = ref_sequence_200();

    // new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH).createContext(var, read, 20, REF_SEQUENCE_200)
    let read_context = create_context(&var, &read, 20, &ref_sequence, DEFAULT_FLANK_LENGTH).unwrap();

    let matcher = create_read_context_matcher(&read_context);

    let ref_bases = &REF_BASES_200[30..70];
    let ref_read = build_sam_record(30, "40M", ref_bases);

    assert_eq!(MatchType::Ref, determine_read_match(&matcher, &read_context, &ref_read, 20));
}

// Java: testIndelHomologyMatches
#[test]
fn test_indel_homology_matches() {
    let position = 64usize;
    let ref_ = &REF_BASES_200[position..position + 1];
    let var = create_simple_variant(position as i32, ref_, &format!("{ref_}AAAA"));

    let mut read_bases = format!("{}AAAA{}", &REF_BASES_200[40..position + 1], &REF_BASES_200[position + 1..90]);
    let mut read_cigar = "25M4I25M";
    let mut read = build_sam_record(40, read_cigar, &read_bases);

    let ref_sequence = ref_sequence_200();

    let read_context = create_context(&var, &read, 24, &ref_sequence, DEFAULT_FLANK_LENGTH).unwrap();

    let matcher = create_read_context_matcher(&read_context);

    assert_eq!(MatchType::Full, determine_read_match(&matcher, &read_context, &read, 24));

    // values: [20, 23, 24, 11, 12, 15]
    // now with a read with low-qual mismatches in non-permitted locations

    read_bases = format!("{}ACAA{}", &REF_BASES_200[40..position + 1], &REF_BASES_200[position + 1..90]);
    read_cigar = "25M4I25M";
    read = build_sam_record(40, read_cigar, &read_bases);
    read.quals[26] = 11;

    assert_eq!(MatchType::Full, determine_read_match(&matcher, &read_context, &read, 24));

    read_bases = format!("{}AAAC{}", &REF_BASES_200[40..position + 1], &REF_BASES_200[position + 1..90]);
    read_cigar = "25M4I25M";
    read = build_sam_record(40, read_cigar, &read_bases);
    read.quals[28] = 11;

    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, 24));

    // 1:64 G>GAAAA read(TCTGTGACTC-GGAAAAAAAAAAAACT-CCCTGACCCC 12M4I20M) pos(53-84) index(10-11-25) repeat(A-12 index(12-23)) homology(AAAA length=8) alt(11-24) ref(GGAAAAAAAACT)
    read_bases = format!(
        "{}AAAA{}T{}",
        &REF_BASES_200[40..position + 1],
        &REF_BASES_200[position + 1..position + 8],
        &REF_BASES_200[position + 9..90]
    );

    read_cigar = "25M4I25M";
    read = build_sam_record(40, read_cigar, &read_bases);
    read.quals[36] = 11;

    assert_eq!(MatchType::None, determine_read_match(&matcher, &read_context, &read, 24));
}

// Java: testFlankMatches
#[test]
fn test_flank_matches() {
    let left_core = "AC";
    let right_core = "GT";
    let position = 100;
    let ref_ = "A";
    let alt = "T";
    let variant = create_simple_variant(position, ref_, alt);
    let read_context = create_read_context_cores(&variant, left_core, right_core);
    let left_flank = read_context.left_flank_str();
    let right_flank = read_context.right_flank_str();

    let matcher = create_read_context_matcher(&read_context);

    let mut flank_mismatch = "AAAAATTTTT".to_string();
    let mut read_bases = format!("{left_flank}{left_core}{alt}{right_core}{flank_mismatch}");
    let mut read_qualities = default_quals(read_bases.len());
    let read_var_index = read_context.left_flank_length() + 2;
    let mut cigar = build_cigar_string(read_bases.len());

    let mut read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::Core, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // low-qual mismatches in the flanks - 3 are permitted
    read_bases = format!("{flank_mismatch}{left_core}{alt}{right_core}{flank_mismatch}");
    read_qualities = default_quals(read_bases.len());

    for i in 0..10 {
        read_qualities[i] = 10;
    }

    for i in read_bases.len() - 11..read_bases.len() {
        read_qualities[i] = 10;
    }

    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::Core, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // now limited to 3
    flank_mismatch = format!("{}TTT", &right_flank[0..7]);
    read_bases = format!("{left_flank}{left_core}{alt}{right_core}{flank_mismatch}");
    read_qualities = default_quals(read_bases.len());

    for i in read_bases.len() - 4..read_bases.len() {
        read_qualities[i] = 10;
    }

    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::Full, determine_read_match(&matcher, &read_context, &read, read_var_index));

    // missing flank still classifies as FULL if there are no mismatches
    read_bases = format!("{left_flank}{left_core}{alt}{right_core}");
    read_qualities = default_quals(read_bases.len());
    cigar = build_cigar_string(read_bases.len());

    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    assert_eq!(MatchType::Full, determine_read_match(&matcher, &read_context, &read, read_var_index));
}

// Java: testAverageCoreQuality
#[test]
fn test_average_core_quality() {
    let left_core = TEST_LEFT_CORE;
    let right_core = TEST_RIGHT_CORE;
    let position = 100;
    let ref_ = "A";
    let alt = "T";
    let variant = create_simple_variant(position, ref_, alt);

    let read_context = create_read_context_cores(&variant, left_core, right_core);
    let _matcher = create_read_context_matcher(&read_context);
    let left_flank = read_context.left_flank_str();
    let right_flank = read_context.right_flank_str();

    let mut read_bases = format!("{left_flank}{left_core}{alt}{right_core}{right_flank}");
    let mut read_qualities = default_quals(read_bases.len());
    let mut read_var_index = read_context.left_flank_length() + left_core.len() as i32;
    let rvi = read_var_index as usize;

    read_qualities[rvi - 2] = 10;
    read_qualities[rvi - 1] = 12;
    read_qualities[rvi] = 14;
    read_qualities[rvi + 1] = 16;
    read_qualities[rvi + 2] = 18;

    let cigar = build_cigar_string(read_bases.len());

    let mut read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    let mut average = average_core_quality(&read_context, &read, read_var_index);
    assert!((14.0 - average).abs() <= 0.01, "average {average}");

    // now with partial cores
    read_bases = format!("{left_flank}{left_core}{alt}");
    read_qualities = default_quals(read_bases.len());

    read_qualities[rvi - 2] = 12;
    read_qualities[rvi - 1] = 13;
    read_qualities[rvi] = 14;

    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    average = average_core_quality(&read_context, &read, read_var_index);
    assert!((13.0 - average).abs() <= 0.01, "average {average}");

    // partial on the left
    read_bases = format!("{alt}{right_core}{right_flank}");
    read_qualities = default_quals(read_bases.len());

    read_qualities[0] = 10;
    read_qualities[1] = 11;
    read_qualities[2] = 12;

    read_var_index = 0;
    read = build_sam_record_quals(position - read_var_index, &cigar, &read_bases, read_qualities.clone());

    average = average_core_quality(&read_context, &read, read_var_index);
    assert!((11.0 - average).abs() <= 0.01, "average {average}");
}
