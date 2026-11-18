//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! evidence/JitterTest.java. JitterMatch is `evidence::Jitter`, JitterMatch.checkJitter / hasJitterMatchType are
//! `evidence::check_jitter` / `has_jitter_match_type`, MsiModelParams + JitterModelParams are
//! `quality::JitterParams`, MsiJitterCalcs (for one sample) is `quality::MsiJitter` and
//! JitterData.setJitterQualFilterState is `quality::set_jitter_qual_filter_state`.

use super::util::*;
use crate::constants::*;
use crate::evidence::{check_jitter, has_jitter_match_type, Jitter, JitterData, ReadContextCounter};
use crate::matcher::Matcher;
use crate::quality::{set_jitter_qual_filter_state, JitterParams, MsiJitter};
use crate::read_context::{create_context, VariantReadContext};
use crate::refseq::RefSequence;
use crate::repeat::{find_max_repeat, RepeatInfo};
use crate::types::*;
use std::sync::Arc;

fn assert_close(expected: f64, actual: f64, delta: f64) {
    assert!((expected - actual).abs() <= delta, "expected {expected} but was {actual} (delta {delta})");
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

fn variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

/// new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH).createContext(...) followed by assertTrue(isValid())
/// where the Java test has it; Java returns null where the Rust port returns None.
fn build_read_context(variant: &SimpleVariant, read: &Read, var_read_index: i32, ref_seq: &RefSequence) -> VariantReadContext {
    create_context(variant, read, var_read_index, ref_seq, DEFAULT_FLANK_LENGTH).expect("read context")
}

fn has_repeat(ctx: &VariantReadContext, bases: &str, index: i32, count: i32) -> bool {
    // RepeatInfo.matches(bases, index, count)
    ctx.all_repeats.iter().any(|x| x.index == index && x.bases == bases.as_bytes() && x.count == count)
}

const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";

/// VariantUtils.createReadContext(variant, leftCore, rightCore): no repeat or homology, default flanks
fn create_read_context(variant: &SimpleVariant, left_core: &str, right_core: &str) -> VariantReadContext {
    let (left_flank, right_flank) = (TEST_LEFT_FLANK, TEST_RIGHT_FLANK);
    let core_index_start = left_flank.len() as i32;
    let var_read_index = core_index_start + left_core.len() as i32;
    let core_index_end = var_read_index + variant.alt_len() - 1 + right_core.len() as i32;
    let ref_bases = format!("{}{}{}", left_core, variant.ref_str(), right_core);
    let read_bases = format!("{}{}{}{}{}", left_flank, left_core, variant.alt_str(), right_core, right_flank);
    let alignment_start = variant.pos - var_read_index;
    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + right_core.len() as i32 + right_flank.len() as i32 - 1;
    VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.into_bytes(),
        read_cigar: format!("{}M", read_bases.len()),
        read_bases: read_bases.into_bytes(),
        core_index_start,
        var_index: var_read_index,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start: variant.pos - left_core.len() as i32,
        core_position_end: variant.pos + right_core.len() as i32,
        has_indel_in_core: false,
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    }
}

/// VariantReadContextBuilder.setMaxRefRepeat
fn set_max_ref_repeat(ctx: &mut VariantReadContext) {
    let ref_repeat_index = if ctx.variant.is_indel() { ctx.variant_ref_index() + 1 } else { ctx.variant_ref_index() };
    ctx.max_ref_repeat = find_max_repeat(
        &ctx.ref_bases, ref_repeat_index, ref_repeat_index, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT + 1, false,
        ref_repeat_index,
    );
}

/// new JitterModelParams(repeatUnit, ConsensusType.NONE, scale4, scale5, scale6, gradient, intercept, skew)
fn jitter_params(repeat_unit: &str, scale4: f64, scale5: f64, scale6: f64, gradient: f64, intercept: f64, skew: f64) -> JitterParams {
    JitterParams { repeat_unit: repeat_unit.to_string(), consensus: "NONE".to_string(), scale4, scale5, scale6, gradient, intercept, skew }
}

/// MsiJitterCalcs.calcErrorRate(readContext, sampleId)
fn calc_error_rate(msi: &MsiJitter, ctx: &VariantReadContext) -> f64 {
    match MsiJitter::find_repeat(ctx) {
        Some(repeat) => msi.calc_error_rate(&ctx.variant, &repeat),
        None => 0.0,
    }
}

/// JitterData.setValues
fn set_values(jitter: &mut JitterData, shortened: i32, lengthened: i32, full_support: i32) {
    jitter.shortened = shortened;
    jitter.lengthened = lengthened;
    jitter.valid_qual_full_support = full_support;
    jitter.qual_boost = 1.0;
    jitter.filter_on_noise = false;
    jitter.hard_filter_on_noise = false;
}

// JitterTest.testJitterMatching
#[test]
fn test_jitter_matching() {
    //                         10        20        30
    //               012345678901234567890123456789012
    let ref_bases = "GTCGTCGTCGTAAAACAAAAACTCGTCGTCGTC";

    let ref_sequence = RefSequence::from_bases(100, ref_bases.as_bytes().to_vec());

    let variant_read_bases = format!("{}A{}", &ref_bases[0..15], &ref_bases[16..]);
    let var = variant(115, "C", "A");
    let read = build_sam_record(100, &build_cigar_string(variant_read_bases.len()), &variant_read_bases);

    let read_context = build_read_context(&var, &read, 15, &ref_sequence);
    assert!(read_context.is_valid());
    assert_eq!(15, read_context.var_index);
    assert_eq!(1, read_context.all_repeats.len());
    let repeat: RepeatInfo = read_context.all_repeats[0].clone();
    assert_eq!(11, repeat.index);
    assert_eq!(10, repeat.count);

    let mut read_quals = default_quals(variant_read_bases.len());

    let matcher = Matcher::new(&read_context, true, true);

    // The Java calls below pass (..., 16, matcher.altIndexLower(), matcher.altIndexUpper(), ...) for the
    // parameters (altIndexLower, altIndexUpper, readVarIndex): altIndexLower() lands in the altIndexUpper
    // slot and altIndexUpper() in the readVarIndex slot. The same values are passed here; the Rust function
    // has no altIndexLower parameter (Java only reads it in a branch that is never reached).

    // low-qual base mismatches are permitted within the core in specific locations (only one) and outside the core

    let read_mismatches = format!(
        "T{}A{}G{}",
        &variant_read_bases[1..11],
        &variant_read_bases[11..22],
        &variant_read_bases[23..]
    );
    read_quals[0] = 11;
    read_quals[23] = 11;

    assert!(has_jitter_match_type(
        &repeat, &read_context, matcher.alt_index_lower, matcher.alt_index_upper, read_mismatches.as_bytes(),
        &read_quals, Jitter::Lengthened, false, 1
    ));

    // 2 core mismatches are not permitted
    let mut read_quals = default_quals(variant_read_bases.len());
    let read_mismatches = format!("{}ATT{}", &variant_read_bases[0..11], &variant_read_bases[13..]);
    read_quals[12] = 11;
    read_quals[13] = 11;

    assert!(!has_jitter_match_type(
        &repeat, &read_context, matcher.alt_index_lower, matcher.alt_index_upper, read_mismatches.as_bytes(),
        &read_quals, Jitter::Lengthened, false, 1
    ));

    // and not within the critical range
    let mut read_quals = default_quals(variant_read_bases.len());
    let read_mismatches =
        format!("{}A{}T{}", &variant_read_bases[0..11], &variant_read_bases[11..15], &variant_read_bases[16..]);
    read_quals[16] = 11;

    assert!(!has_jitter_match_type(
        &repeat, &read_context, matcher.alt_index_lower, matcher.alt_index_upper, read_mismatches.as_bytes(),
        &read_quals, Jitter::Lengthened, false, 1
    ));
}

// JitterTest.testDoubleJitter
#[test]
fn test_double_jitter() {
    //                         10        20        30
    //               0123456789012345678901234567890123456789012345678
    let ref_bases = "GTCGTCGTCGTATAAAAAAAAAACTCGTAGTCGTCATTCGTCATTCATT";

    let var_index = 23usize;
    let req_pos_start = 100;
    let position = req_pos_start + var_index as i32;
    let ref_sequence = RefSequence::from_bases(req_pos_start, ref_bases.as_bytes().to_vec());

    let variant_read_bases = format!("{}G{}", &ref_bases[0..var_index], &ref_bases[var_index + 1..]);
    let var = variant(position, "C", "G");
    let read = build_sam_record(req_pos_start, &build_cigar_string(variant_read_bases.len()), &variant_read_bases);

    let read_context = build_read_context(&var, &read, var_index as i32, &ref_sequence);
    assert!(read_context.is_valid());
    assert_eq!(21, read_context.var_index);
    assert_eq!(1, read_context.all_repeats.len());
    let repeat: RepeatInfo = read_context.all_repeats[0].clone();
    assert_eq!(11, repeat.index);
    assert_eq!(10, repeat.count);

    // VariantUtils.createReadContextMatcher
    let matcher = Matcher::new(&read_context, true, false);

    let read_quals = default_quals(variant_read_bases.len());
    let var_index = var_index as i32;

    // double jitter for variants with long repeats
    let read_bases = format!("{}A{}", &variant_read_bases[0..16], &variant_read_bases[16..]);

    assert!(has_jitter_match_type(
        &repeat, &read_context, matcher.alt_index_upper, var_index + 1, read_bases.as_bytes(), &read_quals,
        Jitter::Lengthened, true, 1
    ));

    let read_bases = format!("{}AA{}", &variant_read_bases[0..16], &variant_read_bases[16..]);

    assert!(has_jitter_match_type(
        &repeat, &read_context, matcher.alt_index_upper, var_index + 2, read_bases.as_bytes(), &read_quals,
        Jitter::Lengthened, true, 2
    ));

    let read_bases = format!("{}{}", &variant_read_bases[0..16], &variant_read_bases[18..]);

    assert!(has_jitter_match_type(
        &repeat, &read_context, matcher.alt_index_upper, var_index - 2, read_bases.as_bytes(), &read_quals,
        Jitter::Shortened, true, 2
    ));
}

// JitterTest.testReadJitter
#[test]
fn test_read_jitter() {
    let var = variant(30, "A", "T");
    let alt = var.alt_str().to_string();

    //                         10        20        30        40        50        60        70        80        90
    //               0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    let ref_bases = "ATGTCATTTTAGCGCGCATTCCTTCCTTCCAAAAAAAAAATGCTGGCACACACACATTAGTCGTAGATTAGTCGTAG";
    let ref_sequence = RefSequence::from_bases(0, ref_bases.as_bytes().to_vec());

    let read_bases = format!("{}{}{}", &ref_bases[1..30], alt, &ref_bases[31..71]);
    let read_cigar = "70M".to_string();
    let read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = build_read_context(&var, &read, 29, &ref_sequence);
    assert!(read_context.is_valid());
    assert_eq!(23, read_context.var_index);
    assert_eq!(2, read_context.all_repeats.len());
    assert!(has_repeat(&read_context, "A", 24, 9));
    assert!(has_repeat(&read_context, "TTCC", 11, 3));

    let matcher = Matcher::new(&read_context, true, true);

    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 29, false);
    assert_eq!(Jitter::None, jitter_match);

    let read_bases = format!("{}TTCC{}{}", &ref_bases[1..30], alt, &ref_bases[31..71]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 33, false);
    assert_eq!(Jitter::Lengthened, jitter_match);

    // as in Java, the cigar from the previous read is reused here
    let read_bases = format!("{}{}{}", &ref_bases[1..26], alt, &ref_bases[31..71]);
    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 25, false);
    assert_eq!(Jitter::Shortened, jitter_match);

    // now test reads where the jitter is after the variant read index
    // 28 bases then the variant at the index = 29 then an extra 'A'
    let read_bases = format!("{}{}A{}", &ref_bases[1..30], alt, &ref_bases[31..71]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 29, false);
    assert_eq!(Jitter::Lengthened, jitter_match);

    let read_bases = format!("{}{}{}", &ref_bases[1..30], alt, &ref_bases[32..71]);
    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 29, false);
    assert_eq!(Jitter::Shortened, jitter_match);

    // test a partial core on each side - still valid for jitter

    // first missing the left core
    let read_bases = format!("{}{}A{}", &ref_bases[13..30], alt, &ref_bases[31..71]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read1 = build_sam_record(13, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 17, false);
    assert_eq!(Jitter::Lengthened, jitter_match);

    // now missing the right core
    let read_bases = format!("{}{}{}", &ref_bases[1..30], alt, &ref_bases[32..46]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 29, false);
    assert_eq!(Jitter::Shortened, jitter_match);

    // test read not covering the repeat so cannot test it correctly
    let read_bases = format!("{}{}{}", &ref_bases[1..30], alt, &ref_bases[32..40]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read1 = build_sam_record(1, &read_cigar, &read_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read1, 29, false);
    assert_eq!(Jitter::None, jitter_match);
}

// JitterTest.testJitterSpanningVariantPos
#[test]
fn test_jitter_spanning_variant_pos() {
    let var = variant(30, "G", "A");
    let alt = var.alt_str().to_string();

    //                         10        20        30        40        50        60        70
    //               01234567890123456789012345678901234567890123456789012345678901234567890
    let ref_bases = "ATGTCATTTTAGCGCGCATTCTACACACACGCACACACACACACACACATTAGTCGTAGATTAGTCGTAGT";
    let ref_sequence = RefSequence::from_bases(0, ref_bases.as_bytes().to_vec());

    let read_bases = format!("{}{}{}", &ref_bases[1..30], alt, &ref_bases[31..71]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = build_read_context(&var, &read, 29, &ref_sequence);

    let matcher = Matcher::new(&read_context, true, false);

    // test 1: insert (ie lengthened) with a matching CIGAR
    let read1_bases = format!("{}AC{}{}", &ref_bases[1..30], alt, &ref_bases[31..71]);
    let read_with_indel = build_sam_record(1, "21M2I49M", &read1_bases);

    let jitter_match = check_jitter(&read_context, &matcher, &read_with_indel, 31, false);
    assert_eq!(Jitter::Lengthened, jitter_match);

    // test 2: insert but without a matching CIGAR
    let read_without_indel = build_sam_record(1, &build_cigar_string(read1_bases.len()), &read1_bases);
    let jitter_match = check_jitter(&read_context, &matcher, &read_without_indel, 29, false);
    assert_eq!(Jitter::Lengthened, jitter_match);

    // test 3: deleted repeat bases matched in the CIGAR
    let read2_bases = format!("{}{}{}", &ref_bases[1..28], alt, &ref_bases[31..71]);
    let read_with_indel = build_sam_record(1, "21M2D47M", &read2_bases);
    let jitter_match = check_jitter(&read_context, &matcher, &read_with_indel, 27, false);
    assert_eq!(Jitter::Shortened, jitter_match);

    // test 4: deleted but all bases appear aligned
    let read_without_indel = build_sam_record(1, &build_cigar_string(read2_bases.len()), &read2_bases);
    let jitter_match = check_jitter(&read_context, &matcher, &read_without_indel, 29, false);
    assert_eq!(Jitter::Shortened, jitter_match);
}

// JitterTest.testJitterEndingOnVariantPos
#[test]
fn test_jitter_ending_on_variant_pos() {
    let var = variant(49, "T", "C");
    let alt = var.alt_str().to_string();

    //                         10        20        30        40        50        60        70
    //               01234567890123456789012345678901234567890123456789012345678901234567890
    let ref_bases = "ATGTCATTTTAGCGCGCATTCTACACACACACACACACACACACACACATTAGTCGTAGATTAGTCGTAGT";
    let ref_sequence = RefSequence::from_bases(0, ref_bases.as_bytes().to_vec());

    let read_bases = format!("{}{}{}", &ref_bases[1..49], alt, &ref_bases[50..71]);
    let read_cigar = build_cigar_string(read_bases.len());
    let read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = build_read_context(&var, &read, 48, &ref_sequence);

    let matcher = Matcher::new(&read_context, true, false);

    let read1_bases = format!("{}CA{}{}", &ref_bases[1..49], alt, &ref_bases[50..71]);
    let read1_with_cigar_elem = build_sam_record(1, "21M2I49M", &read1_bases);
    let jitter_match_with_cigar_elem = check_jitter(&read_context, &matcher, &read1_with_cigar_elem, 50, false);
    assert_eq!(Jitter::Lengthened, jitter_match_with_cigar_elem);

    let read1_no_cigar_elem = build_sam_record(1, &build_cigar_string(read1_bases.len()), &read1_bases);
    let jitter_match_no_cigar_elem = check_jitter(&read_context, &matcher, &read1_no_cigar_elem, 48, false);
    assert_eq!(Jitter::Lengthened, jitter_match_no_cigar_elem);
}

// JitterTest.testMsiJitterCalcs
#[test]
fn test_msi_jitter_calcs() {
    let jitter_params = jitter_params("A/T", 0.05, 0.05, 0.11, 0.0444, -0.1835, 1.0164);

    // MsiModelParams.calcErrorRate
    let calc = jitter_params.calc_error_rate(8, -1, None);
    assert_close(0.002986, calc, 0.0001);

    let calc = jitter_params.calc_error_rate(6, 1, Some(jitter_params.scale6));
    assert_close(0.00011, calc, 0.0001);
}

// JitterTest.testMsiJitterErrorRate
#[test]
fn test_msi_jitter_error_rate() {
    let jitter_params1 = jitter_params("A/T", 0.05, 0.06, 0.07, 0.0444, -0.1835, 1.0164);

    let jitter_params2 = jitter_params("AT/TA", 0.02, 0.03, 0.04, 0.0444, -0.1835, 1.0164);

    // REPEAT_UNIT_3_PLUS_LABEL
    let jitter_params3 = jitter_params("3-5bp repeat", 0.1, 0.11, 0.12, 0.0133, 0.1115, 1.6087);

    // msiJitterCalcs.setSampleParams(TEST_SAMPLE, ...)
    let msi_jitter_calcs = MsiJitter { params: vec![jitter_params1, jitter_params2, jitter_params3], probable_msi: false };

    let var = variant(100, "TA", "T");
    let mut read_context = create_read_context(&var, "GT", "AAAAAAAAAT");
    set_max_ref_repeat(&mut read_context);

    let error_rate = calc_error_rate(&msi_jitter_calcs, &read_context);
    assert_close(2.09e-2, error_rate, 1e-4);

    let var = variant(100, "T", "TAA");
    let read_context = create_read_context(&var, "GT", "AAAAAAAAAAAAT");

    let error_rate = calc_error_rate(&msi_jitter_calcs, &read_context);
    assert_close(2.81e-3, error_rate, 1e-5);

    // too many repeats changing
    let var = variant(100, "T", "TAAAAAAA");
    let read_context = create_read_context(&var, "GT", "AAAAT");

    let error_rate = calc_error_rate(&msi_jitter_calcs, &read_context);
    assert_close(0.0, error_rate, 0.001);

    // 4-base repeat
    let var = variant(100, "ACGTA", "A");
    let read_context = create_read_context(&var, "GT", "CGTACGTACGTACGTACGTGG");

    let error_rate = calc_error_rate(&msi_jitter_calcs, &read_context);
    assert_close(1.99e-4, error_rate, 1e-6);
}

// JitterTest.testJitterFilterQualState
#[test]
fn test_jitter_filter_qual_state() {
    let jitter_params1 = jitter_params("A/C/G/T", 0.05, 0.06, 0.07, 0.05, -0.02, 0.05);

    // msiJitterCalcs.setSampleParams(TEST_SAMPLE, ...)
    let msi_jitter_calcs = MsiJitter { params: vec![jitter_params1], probable_msi: false };

    let var = variant(100, "T", "TA");

    let read_bases = format!("{}{}AAAAAA{}", &REF_BASES_200[80..100], var.alt_str(), &REF_BASES_200[101..121]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = build_cigar_string(read_bases.len());
    let read = build_sam_record_quals(80, &read_cigar, &read_bases, base_quals);

    let fixture = Fixture::new();
    let read_context = build_read_context(&var, &read, 20, &fixture.ref_seq);

    // VariantUtils.createReadCounter(0, readContext)
    let mut read_context_counter = ReadContextCounter::new(
        Arc::new(read_context), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &fixture.qc(), false,
    );

    // lengthened is noise
    set_values(&mut read_context_counter.jitter, 10, 10, 50);

    set_jitter_qual_filter_state(&msi_jitter_calcs, &mut read_context_counter);

    assert_close(1.2, read_context_counter.jitter.qual_boost, 0.01);
    assert!(!read_context_counter.jitter.filter_on_noise);

    // shortened is noise
    set_values(&mut read_context_counter.jitter, 5, 35, 50);

    set_jitter_qual_filter_state(&msi_jitter_calcs, &mut read_context_counter);

    assert_close(1.1, read_context_counter.jitter.qual_boost, 0.01);
    assert!(!read_context_counter.jitter.filter_on_noise);

    // both noise
    set_values(&mut read_context_counter.jitter, 5, 2, 50);

    set_jitter_qual_filter_state(&msi_jitter_calcs, &mut read_context_counter);

    assert_close(1.14, read_context_counter.jitter.qual_boost, 0.01);
    assert!(!read_context_counter.jitter.filter_on_noise);

    // shortened too high
    set_values(&mut read_context_counter.jitter, 21, 1, 10);

    set_jitter_qual_filter_state(&msi_jitter_calcs, &mut read_context_counter);

    assert!(read_context_counter.jitter.filter_on_noise);
    assert_close(1.0, read_context_counter.jitter.qual_boost, 0.01);

    // again for lengthened
    set_values(&mut read_context_counter.jitter, 1, 60, 2);

    set_jitter_qual_filter_state(&msi_jitter_calcs, &mut read_context_counter);

    assert!(read_context_counter.jitter.filter_on_noise);
    assert_close(1.0, read_context_counter.jitter.qual_boost, 0.01);

    // combined
    set_values(&mut read_context_counter.jitter, 126, 126, 35);

    set_jitter_qual_filter_state(&msi_jitter_calcs, &mut read_context_counter);

    assert!(read_context_counter.jitter.filter_on_noise);
    assert_close(1.0, read_context_counter.jitter.qual_boost, 0.01);
}
