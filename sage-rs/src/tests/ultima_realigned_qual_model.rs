//! Ported from the Java SAGE test suite: seqtech/UltimaRealignedQualModelTest.java
//!
//! The Java class shares one static MockRefGenome between its tests; each test here builds its own from
//! the ref homopolymers (populateRefGenomeBases).

use super::util::*;
use crate::constants::*;
use crate::read_context::{has_indel_in_core, VariantReadContext};
use crate::types::*;
use crate::ultima::*;

/// VariantUtils.TEST_LEFT_FLANK etc
const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";
const TEST_LEFT_CORE: &str = "AC";
const TEST_RIGHT_CORE: &str = "GT";

fn hp(base: char, length: i32) -> Homopolymer {
    Homopolymer::new(base as u8, length)
}

fn variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

fn homopolymers_of(bases: &str) -> Vec<Homopolymer> {
    get_homopolymers(bases.as_bytes(), 0, bases.len() as i32 - 1)
}

/// new VariantReadContext(variant, alignmentStart, alignmentEnd, refBases, readBases, readCigar,
/// coreIndexStart, varIndex, coreIndexEnd, null, null, emptyList, corePositionStart, corePositionEnd)
#[allow(clippy::too_many_arguments)]
fn read_context(
    variant: &SimpleVariant, alignment_start: i32, alignment_end: i32, ref_bases: &[u8], read_bases: &[u8], cigar: &str,
    core_index_start: i32, var_index: i32, core_index_end: i32, core_position_start: i32, core_position_end: i32,
) -> VariantReadContext {
    let read_cigar = parse_cigar(cigar);
    VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.to_vec(),
        read_bases: read_bases.to_vec(),
        core_index_start,
        var_index,
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

/// VariantUtils.createReadContext(variant): a context with no repeat or homology
fn create_read_context(variant: &SimpleVariant) -> VariantReadContext {
    let core_index_start = TEST_LEFT_FLANK.len() as i32;
    let var_read_index = core_index_start + TEST_LEFT_CORE.len() as i32;
    let core_index_end = var_read_index + variant.alt_len() - 1 + TEST_RIGHT_CORE.len() as i32;
    let ref_bases = format!("{}{}{}", TEST_LEFT_CORE, variant.ref_str(), TEST_RIGHT_CORE);
    let read_bases = format!("{}{}{}{}{}", TEST_LEFT_FLANK, TEST_LEFT_CORE, variant.alt_str(), TEST_RIGHT_CORE, TEST_RIGHT_FLANK);
    let alignment_start = variant.pos - var_read_index;
    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + (TEST_RIGHT_CORE.len() + TEST_RIGHT_FLANK.len()) as i32 - 1;
    read_context(
        variant, alignment_start, alignment_end, ref_bases.as_bytes(), read_bases.as_bytes(),
        &format!("{}M", read_bases.len()), core_index_start, var_read_index, core_index_end,
        variant.pos - TEST_LEFT_CORE.len() as i32, variant.pos + TEST_RIGHT_CORE.len() as i32,
    )
}

/// UltimaRealignedQualModelTest.createVariantReadContext(variant, readBases, refBases, coreIndexStart, varIndex)
fn create_variant_read_context_at(
    variant: &SimpleVariant, read_bases: &[u8], ref_bases: &[u8], core_index_start: i32, var_index: i32,
) -> VariantReadContext {
    let core_pos_start = variant.pos - (var_index - core_index_start); // assumes no other indels
    let core_pos_end = variant.pos + variant.indel_length() + 2;
    let core_index_end = var_index + variant.indel_length() + 2;
    let align_start = variant.pos - var_index;
    let align_end = variant.pos + (core_index_end - var_index) + DEFAULT_FLANK_LENGTH;
    read_context(
        variant, align_start, align_end, ref_bases, read_bases, "25M", core_index_start, var_index, core_index_end,
        core_pos_start, core_pos_end,
    )
}

/// UltimaRealignedQualModelTest.createVariantReadContext(variant, readBases, refBases)
fn create_variant_read_context(variant: &SimpleVariant, read_bases: &[u8], ref_bases: &[u8]) -> VariantReadContext {
    let core_index_start = DEFAULT_FLANK_LENGTH;
    let var_index = core_index_start + 2;
    create_variant_read_context_at(variant, read_bases, ref_bases, core_index_start, var_index)
}

// UltimaRealignedQualModelTest.testGetHomopolymersInvalidIndices
#[test]
fn test_get_homopolymers_invalid_indices() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let bases = b"AA";
    let mut homopolymers = get_homopolymers(bases, -1, 1);
    assert_eq!(0, homopolymers.len());

    homopolymers = get_homopolymers(bases, 0, 2);
    assert_eq!(0, homopolymers.len());

    homopolymers = get_homopolymers(bases, 1, 0);
    assert_eq!(0, homopolymers.len());
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testMergeSandwichedHomopolymersMultiHomopolymerSkip
#[test]
fn test_merge_sandwiched_homopolymers_multi_homopolymer_skip() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('C', 3), hp('T', 3), hp('C', 2), hp('T', 10), hp('A', 1)];
    let read_homopolymers = vec![hp('C', 1), hp('T', 1), hp('C', 1), hp('T', 3), hp('C', 2), hp('T', 10), hp('A', 1)];

    let var = variant(100, "A", "AC");
    let read_context = create_read_context(&var);
    let merged_homopolymers = merge_sandwiched_homopolymers(&read_context, &ref_homopolymers, &read_homopolymers);

    let expected_homopolymers = ref_homopolymers;
    assert_eq!(expected_homopolymers, merged_homopolymers.ref_homopolymers);
    assert_eq!(expected_homopolymers, merged_homopolymers.read_homopolymers);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testMergeSandwichedHomopolymersValidationCases
#[test]
fn test_merge_sandwiched_homopolymers_validation_cases() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_bases = "TCAAACAAACAAACAAACAAACAAAAAAAAAAAACT";
    let read_bases = "TCAAACAAACAAACAAACAAACAAACAAAAAAACT";

    let ref_homopolymers = homopolymers_of(ref_bases);
    let read_homopolymers = homopolymers_of(read_bases);

    let var = variant(100, "A", "AC");
    let read_context = create_read_context(&var);
    let merged_homopolymers = merge_sandwiched_homopolymers(&read_context, &ref_homopolymers, &read_homopolymers);

    let expected_merged_ref_homopolymers = ref_homopolymers.clone();
    let mut expected_merged_read_homopolymers = ref_homopolymers.clone();
    let n = expected_merged_read_homopolymers.len();
    expected_merged_read_homopolymers[n - 3] = hp('A', 11);

    assert_eq!(expected_merged_ref_homopolymers, merged_homopolymers.ref_homopolymers);
    assert_eq!(expected_merged_read_homopolymers, merged_homopolymers.read_homopolymers);
    set_test_sequencing_type(None);
}

/// SandwichedHomopolymerTestCase: the constructor and check()
fn check_sandwiched_homopolymer_case(ref_bases: &str, read_bases: &str) {
    let ref_homopolymers = homopolymers_of(ref_bases);
    let read_homopolymers = homopolymers_of(read_bases);
    let expected_homopolymers = vec![Homopolymer::new(ref_bases.as_bytes()[0], ref_bases.len() as i32)];

    let var = variant(50, "A", "AA");
    let context_ref_bases = &REF_BASES_200[30..70];
    let context_read_bases = format!("{}A{}", &REF_BASES_200[30..51], &REF_BASES_200[51..70]);
    let read_context = create_variant_read_context(&var, context_read_bases.as_bytes(), context_ref_bases.as_bytes());
    let merged_homopolymers = merge_sandwiched_homopolymers(&read_context, &ref_homopolymers, &read_homopolymers);

    assert_eq!(expected_homopolymers, merged_homopolymers.ref_homopolymers);
    assert_eq!(expected_homopolymers, merged_homopolymers.read_homopolymers);
}

// UltimaRealignedQualModelTest.testSandwichedHomopolymer
#[test]
fn test_sandwiched_homopolymer() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    for bases in ["ATA", "ATGA", "ATTA", "ATCGA", "ATTTA", "ATAGA", "ATATA"] {
        let single_base = bases[0..1].repeat(bases.len());
        // SandwichedHomopolymerTestCase.fromRefBases, then fromReadBases
        check_sandwiched_homopolymer_case(bases, &single_base);
        check_sandwiched_homopolymer_case(&single_base, bases);
    }
    set_test_sequencing_type(None);
}

// common state for the following realigned qual model tests
const CORE_INDEX_START: i32 = 10;
const CORE_POS_START: i32 = 100;
const VAR_CORE_START_OFFSET: i32 = 4;

fn expand(homopolymers: &[Homopolymer]) -> String {
    homopolymers.iter().map(|h| h.expand()).collect()
}

/// UltimaRealignedQualModelTest.populateRefGenomeBases
fn populate_ref_genome_bases(ref_homopolymers: &[Homopolymer]) -> MockRefGenome {
    let ref_base_padding = "C".repeat(99);
    MockRefGenome::new(&format!("{}{}{}", ref_base_padding, expand(ref_homopolymers), ref_base_padding))
}

/// UltimaRealignedQualModelTest.createTestReadContext
fn create_test_read_context(read_homopolymers: &[Homopolymer]) -> VariantReadContext {
    let flank_bases = "C".repeat(DEFAULT_FLANK_LENGTH as usize);
    let var_position = CORE_POS_START + VAR_CORE_START_OFFSET;
    let var_index = CORE_INDEX_START + VAR_CORE_START_OFFSET;
    let read_bases = format!("{}{}{}", flank_bases, expand(read_homopolymers), flank_bases);
    let var = variant(var_position, "AA", "G");
    create_variant_read_context_at(&var, read_bases.as_bytes(), read_bases.as_bytes(), CORE_INDEX_START, var_index)
}

/// populateRefGenomeBases, createTestReadContext and getRealignedVariants with no ref masks
fn realigned_variants(ref_homopolymers: &[Homopolymer], read_homopolymers: &[Homopolymer]) -> Vec<UltimaRealignedQualModel> {
    let ref_genome = populate_ref_genome_bases(ref_homopolymers);
    let read_context = create_test_read_context(read_homopolymers);
    get_realigned_variants(&read_context, &ref_genome, ref_homopolymers, read_homopolymers, &[]).expect("realigned variants")
}

fn matches(realigned_variant: &UltimaRealignedQualModel, position: i32, ref_: &str, alt: &str) -> bool {
    realigned_variant.variant == variant(position, ref_, alt)
}

/// the assertions on an UltimaHomopolymerAdjustment realigned model: refAdjustCount, and the hp start
/// and end indices with the var read index offset added
fn assert_adjustment(realigned_variant: &UltimaRealignedQualModel, ref_adjust_count: i32, hp_start: i32, hp_end: i32) {
    let qual_model = realigned_variant.qual_model.homopolymer_adjustment().expect("homopolymer adjustment");
    assert_eq!(ref_adjust_count, qual_model.ref_adjust_count);
    assert_eq!(hp_start, qual_model.hp_start_index + realigned_variant.var_read_index_offset);
    assert_eq!(hp_end, qual_model.hp_end_index + realigned_variant.var_read_index_offset);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsVariantStartForHomopolymerContraction
#[test]
fn test_get_realigned_variants_variant_start_for_homopolymer_contraction() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 5), hp('A', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 2), hp('A', 5)];

    let core_index_start = 0;
    let var_index = 4;
    let core_position_start = 100;
    let position = 104;
    let var = variant(position, "ATTT", "A");
    let read_context =
        read_context(&var, -1, -1, &[], b"AAAAATTAAAAA", "", core_index_start, var_index, -1, core_position_start, -1);

    // Java uses whatever an earlier test left in the shared ref genome
    let ref_genome = populate_ref_genome_bases(&ref_homopolymers);
    let realigned_variants =
        get_realigned_variants(&read_context, &ref_genome, &ref_homopolymers, &read_homopolymers, &[]).expect("realigned variants");

    assert_eq!(1, realigned_variants.len());
    assert_eq!(var, realigned_variants[0].variant);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsVariantStartForHomopolymerExpansion
#[test]
fn test_get_realigned_variants_variant_start_for_homopolymer_expansion() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 2), hp('A', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 5), hp('A', 5)];

    let core_index_start = 0;
    let var_index = 4;
    let core_position_start = 100;
    let position = 104;
    let var = variant(position, "A", "ATTT");
    let read_context =
        read_context(&var, -1, -1, &[], b"AAAAATTTTTAAAAA", "", core_index_start, var_index, -1, core_position_start, -1);

    // Java uses whatever an earlier test left in the shared ref genome
    let ref_genome = populate_ref_genome_bases(&ref_homopolymers);
    let realigned_variants =
        get_realigned_variants(&read_context, &ref_genome, &ref_homopolymers, &read_homopolymers, &[]).expect("realigned variants");

    assert_eq!(1, realigned_variants.len());
    assert_eq!(var, realigned_variants[0].variant);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testSandwichRefGenome
#[test]
fn test_sandwich_ref_genome() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_bases = "GTAAAAAAAAAAAAGAGAGAAAAAAC";
    let read_core_bases = "GTGAAAAAAAAAAGAGAGAAAAAAC";

    let ref_homopolymers = homopolymers_of(ref_bases);
    let read_homopolymers = homopolymers_of(read_core_bases);

    let core_position_start = 13;
    let var = variant(15, "AA", "G");
    let read_context = create_variant_read_context(&var, read_core_bases.as_bytes(), ref_bases.as_bytes());
    let merged_homopolymers = merge_sandwiched_homopolymers(&read_context, &ref_homopolymers, &read_homopolymers);

    let ref_mask_index = 18;

    assert_eq!(1, merged_homopolymers.ref_masks.len());

    let ref_mask = merged_homopolymers.ref_masks[0];
    assert_eq!(b'A', ref_mask.base_mask);
    assert_eq!(1, ref_mask.pos_end - ref_mask.pos_start + 1);
    assert_eq!(core_position_start + ref_mask_index, ref_mask.pos_start);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsSingleHomopolymerDeletion
#[test]
fn test_get_realigned_variants_single_homopolymer_deletion() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 2), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(1, realigned_variants.len());

    let realigned_variant = &realigned_variants[0];

    assert!(realigned_variant.model_type() == UltimaModelType::HomopolymerDeletion);
    assert_eq!(0, realigned_variant.var_read_index_offset);
    assert!(matches(realigned_variant, core_position_start + 4, "ATT", "A"));
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsSingleHomopolymerContraction
#[test]
fn test_get_realigned_variants_single_homopolymer_contraction() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 5), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 3), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(1, realigned_variants.len());

    let realigned_variant = &realigned_variants[0];

    assert!(realigned_variant.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant.var_read_index_offset);
    assert!(matches(realigned_variant, core_position_start + 4, "ATT", "A"));
    assert_adjustment(realigned_variant, 2, 1, 3);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsTwoHomopolymerDeletions
#[test]
fn test_get_realigned_variants_two_homopolymer_deletions() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 2), hp('C', 1), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerDeletion);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "ATT", "A"));

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerDeletion);
    assert_eq!(0, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 6, "TC", "A"));
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsTwoHomopolymerContractions
#[test]
fn test_get_realigned_variants_two_homopolymer_contractions() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 5), hp('C', 2), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 3), hp('C', 1), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "ATT", "A"));
    assert_adjustment(realigned_variant1, 2, 1, 3);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(3, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 9, "TC", "T"));
    assert_adjustment(realigned_variant2, 1, 4, 4);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsHomopolyerDeletionThenHomopolymerContraction
#[test]
fn test_get_realigned_variants_homopolyer_deletion_then_homopolymer_contraction() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 5), hp('C', 2), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('C', 1), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerDeletion);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "ATTTTT", "A"));

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 9, "TC", "A"));
    assert_adjustment(realigned_variant2, 1, 1, 1);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsHomopolyerContractionThenHomopolymerDeletion
#[test]
fn test_get_realigned_variants_homopolyer_contraction_then_homopolymer_deletion() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 5), hp('C', 2), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 3), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "ATT", "A"));
    assert_adjustment(realigned_variant1, 2, 1, 3);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerDeletion);
    assert_eq!(3, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 9, "TCC", "T"));
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsSingleHomopolymerCreation
#[test]
fn test_get_realigned_variants_single_homopolymer_creation() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 2), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(1, realigned_variants.len());

    let realigned_variant = &realigned_variants[0];

    assert!(realigned_variant.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant.var_read_index_offset);
    assert!(matches(realigned_variant, core_position_start + 4, "A", "ATT"));
    assert_adjustment(realigned_variant, -2, 1, 2);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsSingleHomopolymerExpansion
#[test]
fn test_get_realigned_variants_single_homopolymer_expansion() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 3), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 5), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(1, realigned_variants.len());

    let realigned_variant = &realigned_variants[0];

    assert!(realigned_variant.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant.var_read_index_offset);
    assert!(matches(realigned_variant, core_position_start + 4, "A", "ATT"));
    assert_adjustment(realigned_variant, -2, 1, 5);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsTwoHomopolymerCreations
#[test]
fn test_get_realigned_variants_two_homopolymer_creations() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 2), hp('C', 1), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "A", "ATT"));
    assert_adjustment(realigned_variant1, -2, 1, 2);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(2, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 4, "A", "TC"));
    assert_adjustment(realigned_variant2, -1, 3, 3);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsTwoHomopolymerExpansions
#[test]
fn test_get_realigned_variants_two_homopolymer_expansions() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 3), hp('C', 1), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 5), hp('C', 2), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "A", "ATT"));
    assert_adjustment(realigned_variant1, -2, 1, 5);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(5, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 7, "T", "TC"));
    assert_adjustment(realigned_variant2, -1, 6, 7);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsHomopolyerCreationThenHomopolymerExpansion
#[test]
fn test_get_realigned_variants_homopolyer_creation_then_homopolymer_expansion() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('C', 1), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 5), hp('C', 2), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "A", "ATTTTT"));
    assert_adjustment(realigned_variant1, -5, 1, 5);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(5, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 4, "A", "TC"));
    assert_adjustment(realigned_variant2, -1, 6, 7);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsHomopolyerExpansionThenHomopolymerCreation
#[test]
fn test_get_realigned_variants_homopolyer_expansion_then_homopolymer_creation() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 3), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('T', 5), hp('C', 2), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "A", "ATT"));
    assert_adjustment(realigned_variant1, -2, 1, 5);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(5, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 7, "T", "TCC"));
    assert_adjustment(realigned_variant2, -2, 6, 7);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsHomopolyerCreationThenHomopolymerCreation
#[test]
fn test_get_realigned_variants_homopolyer_creation_then_homopolymer_creation() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 3), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('G', 1), hp('A', 1), hp('T', 3), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    assert_eq!(2, realigned_variants.len());

    let realigned_variant1 = &realigned_variants[0];
    let realigned_variant2 = &realigned_variants[1];

    assert!(realigned_variant1.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, realigned_variant1.var_read_index_offset);
    assert!(matches(realigned_variant1, core_position_start + 4, "A", "AG"));
    assert_adjustment(realigned_variant1, -1, 1, 1);

    assert!(realigned_variant2.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(1, realigned_variant2.var_read_index_offset);
    assert!(matches(realigned_variant2, core_position_start + 4, "A", "GA"));
    assert_adjustment(realigned_variant2, -1, 2, 2);
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testGetRealignedVariantsSNV
#[test]
fn test_get_realigned_variants_snv() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let ref_homopolymers = vec![hp('A', 5), hp('T', 1), hp('G', 5)];
    let read_homopolymers = vec![hp('A', 5), hp('C', 1), hp('G', 5)];
    let core_position_start = CORE_POS_START;

    let realigned_variants = realigned_variants(&ref_homopolymers, &read_homopolymers);

    let del_realigned_variants: Vec<&UltimaRealignedQualModel> =
        realigned_variants.iter().filter(|x| x.model_type() == UltimaModelType::HomopolymerDeletion).collect();

    let adj_realigned_variants: Vec<&UltimaRealignedQualModel> =
        realigned_variants.iter().filter(|x| x.model_type() == UltimaModelType::HomopolymerAdjustment).collect();

    assert_eq!(1, del_realigned_variants.len());
    assert_eq!(1, adj_realigned_variants.len());

    let del_realigned_variant = del_realigned_variants[0];
    let adj_realigned_variant = adj_realigned_variants[0];

    assert!(del_realigned_variant.model_type() == UltimaModelType::HomopolymerDeletion);
    assert_eq!(0, del_realigned_variant.var_read_index_offset);
    assert!(matches(del_realigned_variant, core_position_start + 4, "AT", "A"));

    assert!(adj_realigned_variant.model_type() == UltimaModelType::HomopolymerAdjustment);
    assert_eq!(0, adj_realigned_variant.var_read_index_offset);
    assert!(matches(adj_realigned_variant, core_position_start + 4, "A", "AC"));
    assert_adjustment(adj_realigned_variant, -1, 1, 1);
    set_test_sequencing_type(None);
}

/// UltimaRealignedQualModelTest.createSimpleSnv
fn create_simple_snv() -> SimpleVariant {
    variant(100, "A", "C")
}

// UltimaRealignedQualModelTest.testIsCleanSnvIndel
#[test]
fn test_is_clean_snv_indel() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let var = variant(100, "A", "AC");
    let read_context = create_variant_read_context(&var, b"AGACCT", b"AGACT");
    assert!(!is_clean_snv(&read_context));
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testIsCleanSnvRefCoreLengthMismatch
#[test]
fn test_is_clean_snv_ref_core_length_mismatch() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let read_context = create_variant_read_context(&create_simple_snv(), b"AATAAA", b"AAAAAA");
    assert!(!is_clean_snv(&read_context));
    assert_eq!(5, read_context.core_length());
    assert!(!is_clean_snv(&read_context));
    set_test_sequencing_type(None);
}

/// the read bases (all A) and the core's ref bases with T at each of the read indices given
fn clean_snv_bases(ref_t_indices: &[usize]) -> (String, String) {
    let core_length = 5;
    let flank_length = 5;
    let core_index_start = flank_length;
    let core_index_end = core_index_start + core_length - 1;
    let read_bases = "A".repeat(core_length + 2 * flank_length);
    let mut ref_bases = read_bases.clone().into_bytes();
    for &i in ref_t_indices {
        ref_bases[i] = b'T';
    }
    let ref_bases = String::from_utf8(ref_bases[core_index_start..=core_index_end].to_vec()).unwrap();
    (read_bases, ref_bases)
}

// UltimaRealignedQualModelTest.testIsCleanSnvClean
#[test]
fn test_is_clean_snv_clean() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let var_index = 5 + 1;
    let (read_bases, ref_bases) = clean_snv_bases(&[var_index]);
    let read_context = create_variant_read_context(&create_simple_snv(), read_bases.as_bytes(), ref_bases.as_bytes());
    // as in Java: the context's own core starts at the default flank length, so this is not clean
    assert!(!is_clean_snv(&read_context));
    set_test_sequencing_type(None);
}

// UltimaRealignedQualModelTest.testIsCleanSnvNotClean
#[test]
fn test_is_clean_snv_not_clean() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    let var_index = 5 + 1;
    let (read_bases, ref_bases) = clean_snv_bases(&[var_index, var_index + 1]);
    let read_context = create_variant_read_context(&create_simple_snv(), read_bases.as_bytes(), ref_bases.as_bytes());
    assert!(!is_clean_snv(&read_context));
    set_test_sequencing_type(None);
}
