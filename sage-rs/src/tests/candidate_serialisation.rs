//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/candidate/CandidateSerialisationTest.java
//!
//! Java's vcf/CandidateSerialisation is `append::to_context` / `append::to_candidate`; a htsjdk
//! VariantContext is `append::VariantContext`.

use super::util::*;
use crate::append::{to_candidate, to_context, VariantContext};
use crate::candidate::Candidate;
use crate::constants::DEFAULT_FLANK_LENGTH;
use crate::read_context::{create_context, VariantReadContext};
use std::sync::Arc;
use crate::refseq::RefSequence;
use crate::types::*;

//             10        20        30        40        50
//   0123456789012345678901234567890123456789012345678901
const REF_BASES: &str = "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

fn ref_sequence() -> RefSequence {
    RefSequence::from_bases(1, REF_BASES.as_bytes().to_vec())
}

/// VariantReadContext.matches
fn matches(a: &VariantReadContext, other: &VariantReadContext) -> bool {
    if a.homology.is_some() != other.homology.is_some() {
        return false;
    }
    if let (Some(h), Some(o)) = (&a.homology, &other.homology) {
        if h.bases != o.bases {
            return false;
        }
    }
    if a.max_repeat.is_some() != other.max_repeat.is_some() {
        return false;
    }
    if let (Some(r), Some(o)) = (&a.max_repeat, &other.max_repeat) {
        // RepeatInfo.matches
        if r.index != o.index || r.bases != o.bases || r.count != o.count {
            return false;
        }
    }
    a.var_index == other.var_index
        && a.alignment_start == other.alignment_start
        && a.alignment_end == other.alignment_end
        && a.ref_bases == other.ref_bases
        && a.read_bases == other.read_bases
        && a.read_cigar == other.read_cigar
        && a.core_position_start == other.core_position_start
        && a.core_position_end == other.core_position_end
}

/// new Candidate(VariantTier.PANEL, readContext, 2, 0) -> toContext -> toCandidate, then
/// assertTrue(recreatedContext.matches(readContext))
fn assert_round_trip(read_context: VariantReadContext, ref_sequence: &RefSequence) {
    let candidate =
        Candidate { tier: Tier::Panel, ctx: Arc::new(read_context), min_events: 2, full_match_support: 0, core_match_support: 0 };

    let variant_context = to_context(CHR_1, &candidate);

    let recreated_candidate = to_candidate(&variant_context, ref_sequence);
    assert!(recreated_candidate.is_ok());

    let recreated_candidate = recreated_candidate.unwrap();
    assert!(matches(&recreated_candidate.ctx, &candidate.ctx));
}

// Java: testCandidateCreation
#[test]
fn test_candidate_creation() {
    let ref_sequence = ref_sequence();

    // test 1: a simple SNV
    let mut position: usize = 21;
    let mut read_bases = format!("{}C{}", &REF_BASES[1..position], &REF_BASES[position + 1..50]);

    let mut read_cigar = format!("{}M", read_bases.len());
    let mut read = build_sam_record(1, &read_cigar, &read_bases);

    let mut var = SimpleVariant::new(position as i32, b"A", b"C");

    let read_context = create_context(&var, &read, position as i32 - 1, &ref_sequence, DEFAULT_FLANK_LENGTH);

    assert!(read_context.as_ref().is_some_and(|x| x.is_valid()));

    assert_round_trip(read_context.unwrap(), &ref_sequence);

    // a delete
    position = 17;
    var = SimpleVariant::new(position as i32, b"ACCC", b"A");

    read_bases = format!("{}{}", &REF_BASES[1..position + 1], &REF_BASES[position + 4..50]);
    read_cigar = "17M3D29M".to_string();
    read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = create_context(&var, &read, position as i32 - 1, &ref_sequence, DEFAULT_FLANK_LENGTH);

    assert!(read_context.as_ref().is_some_and(|x| x.is_valid()));

    assert_round_trip(read_context.unwrap(), &ref_sequence);

    // an insert
    position = 21;
    var = SimpleVariant::new(position as i32, b"A", b"ACCC");

    read_bases = REF_BASES[1..position + 1].to_string(); // up to and including the ref base
    read_bases += &format!("CCC{}", &REF_BASES[position + 2..50]);
    read_cigar = "21M3I27M".to_string();
    read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = create_context(&var, &read, position as i32 - 1, &ref_sequence, DEFAULT_FLANK_LENGTH);

    assert!(read_context.as_ref().is_some_and(|x| x.is_valid()));

    assert_round_trip(read_context.unwrap(), &ref_sequence);
}

// Java: testCreateFromVariantContext
#[test]
fn test_create_from_variant_context() {
    let ref_sequence = ref_sequence();

    let mut genotype = crate::append::Genotype { sample: "SAMPLE".to_string(), ..Default::default() };
    genotype.standard.insert("AD".to_string(), "10,10".to_string());
    genotype.standard.insert("DP".to_string(), "100".to_string());
    genotype.gt = Some("./.".to_string());
    genotype.attributes.insert("AF".to_string(), "0.100".to_string());

    let variant = SimpleVariant::new(17, b"AAAAAT", b"A");

    let mut variant_context = VariantContext {
        chromosome: CHR_1.to_string(),
        position: variant.pos,
        id: ".".to_string(),
        ref_: variant.ref_str().to_string(),
        alts: vec![variant.alt_str().to_string()],
        genotypes: vec![genotype],
        ..Default::default()
    };

    // index                                                 0123456789 012345678901234567890 1234567890
    variant_context
        .info
        .insert("RC_INFO".to_string(), "1-13-TTCAGGTTTG-GAAAAAAATATATATATATAA-GTGCATGAAA-17M5D24M".to_string());

    let candidate = to_candidate(&variant_context, &ref_sequence);
    assert!(candidate.is_ok());
    let candidate = candidate.unwrap();
    assert_eq!(variant, *candidate.variant());
}
