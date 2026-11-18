//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/dedup/VariantDedupTest.java
//!
//! Java's dedup methods take a list of shared SageVariant objects; here the variants live in a Vec
//! and are checked by index (var1 is `variants[0]`, and so on).

use super::candidate_creation::{TEST_LEFT_FLANK, TEST_RIGHT_FLANK};
use super::util::*;
use crate::candidate::Candidate;
use crate::dedup::{dedup_mnv_overlaps, dedup_mnv_snvs, longer_contains_shorter, MixedGermlineDedup};
use crate::evidence::ReadContextCounter;
use crate::filters::{SageVariant, SoftFilter};
use crate::read_context::{has_indel_in_core, VariantReadContext};
use crate::refdata::{ExonData, TranscriptData};
use crate::types::*;
use std::sync::Arc;

/// VariantUtils.TEST_LEFT_CORE / TEST_RIGHT_CORE
const TEST_LEFT_CORE: &str = "AC";
const TEST_RIGHT_CORE: &str = "GT";

/// VariantUtils.createSimpleVariant(position, ref, alt)
pub(super) fn create_simple_variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

/// VariantUtils.createReadContext(variant, leftCore, rightCore, leftFlank, rightFlank): no repeat or homology
pub(super) fn create_read_context(
    variant: SimpleVariant, left_core: &str, right_core: &str, left_flank: &str, right_flank: &str,
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
        variant,
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

/// VariantUtils.createSageVariant(readContext)
pub(super) fn create_sage_variant_from_context(read_context: VariantReadContext) -> SageVariant {
    let fixture = Fixture::new();
    let ctx = Arc::new(read_context);
    let read_counter = ReadContextCounter::new(ctx.clone(), Tier::LowConfidence, 100, 1, &fixture.qc(), false);
    let candidate =
        Candidate { tier: Tier::HighConfidence, ctx, min_events: 1, full_match_support: 1, core_match_support: 0 };
    SageVariant::new(candidate, Vec::new(), vec![read_counter])
}

/// VariantUtils.createSageVariant(position, ref, alt)
fn create_sage_variant(position: i32, ref_: &str, alt: &str) -> SageVariant {
    let variant = create_simple_variant(position, ref_, alt);
    create_sage_variant_from_context(create_read_context(
        variant, TEST_LEFT_CORE, TEST_RIGHT_CORE, TEST_LEFT_FLANK, TEST_RIGHT_FLANK,
    ))
}

/// TestUtils.setTumorQuality
pub(super) fn set_tumor_quality(variant: &mut SageVariant, count: i32, quality: i32) {
    variant.tumor_counters[0].counts.full = count;
    variant.tumor_counters[0].qualities.full = quality;
}

/// TestUtils.addLocalPhaseSet
pub(super) fn add_local_phase_set(variant: &mut SageVariant, lps: i32, read_count: i32) {
    variant.tumor_counters[0].add_local_phase_set(lps, read_count, 0.0);
}

/// TestUtils.clearFilters
fn clear_filters(variants: &mut [SageVariant]) {
    variants.iter_mut().for_each(|x| x.filters.clear());
}

/// GeneTestUtils.createTransExons(GENE_ID_1, TRANS_ID_1, POS_STRAND, {7, 20}, 9, 7, 25, true, ""):
/// exons 7-16 and 20-29. Coding starts on the first base of exon 1, so its start phase is 0; it holds 10
/// coding bases from phase 1, which gives an end phase of 1 and so a start phase of 1 for exon 2.
fn create_trans_exons() -> TranscriptData {
    TranscriptData {
        strand: 1,
        trans_start: 7,
        trans_end: 29,
        coding_start: Some(7),
        coding_end: Some(25),
        exons: vec![
            ExonData { start: 7, end: 16, rank: 1, phase_start: 0 },
            ExonData { start: 20, end: 29, rank: 2, phase_start: 1 },
        ],
    }
}

// VariantDedupTest.testLongerContainsShorter
#[test]
fn test_longer_contains_shorter() {
    let mnv = create_simple_variant(100, "CAC", "TGT");

    assert!(longer_contains_shorter(&create_simple_variant(100, "C", "T"), &mnv));
    assert!(!longer_contains_shorter(&create_simple_variant(100, "C", "C"), &mnv));

    assert!(longer_contains_shorter(&create_simple_variant(101, "A", "G"), &mnv));
    assert!(longer_contains_shorter(&create_simple_variant(102, "C", "T"), &mnv));

    assert!(longer_contains_shorter(&create_simple_variant(100, "CA", "TG"), &mnv));
    assert!(longer_contains_shorter(&create_simple_variant(101, "AC", "GT"), &mnv));

    assert!(longer_contains_shorter(&create_simple_variant(100, "CAC", "TGT"), &mnv));
}

// VariantDedupTest.testMnvDedup
#[test]
fn test_mnv_dedup() {
    // must overlap, be phased and passing

    // ATG -> CGA
    let mut var1 = create_sage_variant(10, "AT", "CG");
    let mut var2 = create_sage_variant(11, "TG", "GA");

    set_tumor_quality(&mut var1, 5, 1000);
    set_tumor_quality(&mut var2, 4, 900);

    let mut variants = vec![var1, var2];

    dedup_mnv_overlaps(&mut variants);

    assert!(variants[0].is_passing());
    assert!(variants[1].is_passing());

    // must be same LPS
    add_local_phase_set(&mut variants[0], 1, 1);
    add_local_phase_set(&mut variants[1], 1, 1);

    dedup_mnv_overlaps(&mut variants);

    // determined by tumor qual
    assert!(variants[0].is_passing());
    assert!(!variants[1].is_passing());

    clear_filters(&mut variants);

    set_tumor_quality(&mut variants[0], 5, 500);

    dedup_mnv_overlaps(&mut variants);

    // determined by tumor qual
    assert!(!variants[0].is_passing());
    assert!(variants[1].is_passing());

    clear_filters(&mut variants);

    // dedup shorter if has same number of changed bases
    let mut var1 = create_sage_variant(10, "ATG", "CTA");
    let mut var2 = create_sage_variant(12, "GC", "AT");
    add_local_phase_set(&mut var1, 1, 1);
    add_local_phase_set(&mut var2, 1, 1);

    let mut variants = vec![var1, var2];

    dedup_mnv_overlaps(&mut variants);

    assert!(!variants[0].is_passing());
    assert!(variants[1].is_passing());

    clear_filters(&mut variants);

    // dedup longer assuming all bases have changed
    let mut var1 = create_sage_variant(10, "AAG", "CTA");
    let mut var2 = create_sage_variant(12, "GC", "AT");
    add_local_phase_set(&mut var1, 1, 1);
    add_local_phase_set(&mut var2, 1, 1);

    let mut variants = vec![var1, var2];

    dedup_mnv_overlaps(&mut variants);

    assert!(variants[0].is_passing());
    assert!(!variants[1].is_passing());
}

// VariantDedupTest.testMnvSnvDedup
#[test]
fn test_mnv_snv_dedup() {
    // ATG -> CGA
    let mut var1 = create_sage_variant(10, "ATG", "CGA");
    let mut var2 = create_sage_variant(10, "A", "C");
    let mut var3 = create_sage_variant(11, "T", "G");
    let mut var4 = create_sage_variant(12, "G", "A");

    // must be same LPS
    add_local_phase_set(&mut var1, 1, 1);
    add_local_phase_set(&mut var2, 1, 1);
    add_local_phase_set(&mut var3, 1, 1);
    add_local_phase_set(&mut var4, 1, 1);

    let mut variants = vec![var1, var2, var3, var4];

    dedup_mnv_snvs(&mut variants);

    assert!(variants[0].is_passing());
    assert!(!variants[1].is_passing());
    assert!(!variants[2].is_passing());
    assert!(!variants[3].is_passing());

    // more complicated example with variants interspersed
    let var1 = create_sage_variant(10, "A", "G"); // no match
    let mut var2 = create_sage_variant(10, "ATG", "CGA");
    let mut var7 = create_sage_variant(11, "TGA", "CCC");
    let var3 = create_sage_variant(11, "T", "A");
    let mut var4 = create_sage_variant(11, "T", "G"); // dedup with first MNV
    let var5 = create_sage_variant(12, "G", "T");
    let mut var6 = create_sage_variant(12, "G", "A"); // deup with first MNV
    let mut var8 = create_sage_variant(13, "A", "C");

    // must be same LPS
    add_local_phase_set(&mut var2, 1, 1);
    add_local_phase_set(&mut var4, 1, 1);
    add_local_phase_set(&mut var6, 1, 1);
    add_local_phase_set(&mut var7, 2, 1);
    add_local_phase_set(&mut var8, 2, 1);

    // index:                0     1     2     3     4     5     6     7
    let mut variants = vec![var1, var2, var3, var4, var5, var6, var7, var8];

    dedup_mnv_snvs(&mut variants);

    assert!(variants[0].is_passing()); // var1
    assert!(variants[1].is_passing()); // var2
    assert!(variants[2].is_passing()); // var3
    assert!(!variants[3].is_passing()); // var4
    assert!(variants[4].is_passing()); // var5
    assert!(!variants[5].is_passing()); // var6
    assert!(variants[6].is_passing()); // var7
    assert!(!variants[7].is_passing()); // var8
}

// VariantDedupTest.testMixedGermline
#[test]
fn test_mixed_germline() {
    // ATG -> CGA
    let mut var1 = create_sage_variant(10, "ATG", "CGA");
    let mut var2 = create_sage_variant(10, "A", "C");
    let mut var3 = create_sage_variant(12, "G", "A"); // marked as germline

    // must be same LPS
    add_local_phase_set(&mut var1, 1, 1);
    add_local_phase_set(&mut var2, 1, 1);
    add_local_phase_set(&mut var3, 1, 1);

    var3.filters.insert(SoftFilter::MinGermlineDepth);

    let mut variants = vec![var1, var2, var3];

    let transcripts = vec![create_trans_exons()];

    let mut deduper = MixedGermlineDedup::new(&transcripts);
    deduper.dedup(&mut variants);

    assert!(variants[0].is_passing());
    assert!(!variants[1].is_passing());
    assert!(!variants[2].is_passing());

    // test again with the MNV spanning a codon so the SNV is kept
    let mut var1 = create_sage_variant(11, "ATG", "CGA");
    let mut var2 = create_sage_variant(12, "T", "G"); // marked as germline
    let mut var3 = create_sage_variant(13, "G", "A");

    // must be same LPS
    add_local_phase_set(&mut var1, 1, 1);
    add_local_phase_set(&mut var2, 1, 1);
    add_local_phase_set(&mut var3, 1, 1);

    var2.filters.insert(SoftFilter::MinGermlineDepth);

    let mut variants = vec![var1, var2, var3];

    deduper.dedup(&mut variants);

    assert!(!variants[0].is_passing());
    assert!(!variants[1].is_passing());
    assert!(variants[2].is_passing());
}
