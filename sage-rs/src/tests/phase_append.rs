//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/evidence/PhaseAppendTest.java

use super::util::*;
use crate::append::{to_context, AppendVariantPhaser};
use crate::candidate::Candidate;
use crate::constants::*;
use crate::evidence::ReadContextCounter;
use crate::read_context::{create_context, VariantReadContext};
use crate::refseq::RefSequence;
use crate::types::*;
use std::sync::Arc;

//             10        20        30        40        50        60        70        80        90        100
//   012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
const REF_BASES: &str =
    "CGCAATATTCGAGTGGAAGTGACTCGATTAATCCGGTGCGTTCGTCACCGCTGTCTATTAAGCCGGTGCGTTCGTCACCGCTGTCTATTTACCCGGTGCGTTCGTC";
//                  A                   T                   C                   G

fn alt_bases() -> String {
    format!("{}A{}T{}C{}", &REF_BASES[0..15], &REF_BASES[16..35], &REF_BASES[36..55], &REF_BASES[56..])
}

const TEST_SAMPLE: &str = "SAMPLE";

/// VariantUtils.createReadCounter(readContext)
fn create_read_counter(ctx: &Arc<VariantReadContext>, fixture: &Fixture) -> ReadContextCounter {
    ReadContextCounter::new(ctx.clone(), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &fixture.qc(), false)
}

// Java: testAppendPhasedEvidence
#[test]
fn test_append_phased_evidence() {
    let fixture = Fixture::new();
    let ref_sequence = RefSequence::from_bases(1, REF_BASES.as_bytes().to_vec());
    let alt_bases = alt_bases();

    let mut variant_phaser = AppendVariantPhaser::new();

    variant_phaser.initialise(TEST_SAMPLE);

    let mut position = 15;
    let var1 = SimpleVariant::new(position, b"G", b"A");

    let mut read_bases = &alt_bases[1..35];
    let mut read_cigar = format!("{}M", read_bases.len());
    let mut read = build_sam_record(1, &read_cigar, read_bases);

    let read_context = create_context(&var1, &read, 14, &ref_sequence, DEFAULT_FLANK_LENGTH).unwrap();
    let candidate1 =
        Candidate { tier: Tier::Panel, ctx: Arc::new(read_context), min_events: 1, full_match_support: 0, core_match_support: 0 };

    let mut var_context1 = to_context(CHR_1, &candidate1);

    position = 35;
    let var2 = SimpleVariant::new(position, b"G", b"T");

    read_bases = &alt_bases[16..55];
    read_cigar = format!("{}M", read_bases.len());
    read = build_sam_record(16, &read_cigar, read_bases);

    let read_context = create_context(&var2, &read, 19, &ref_sequence, DEFAULT_FLANK_LENGTH).unwrap();
    let candidate2 =
        Candidate { tier: Tier::Panel, ctx: Arc::new(read_context), min_events: 1, full_match_support: 0, core_match_support: 0 };

    let mut var_context2 = to_context(CHR_1, &candidate2);

    position = 55;
    let var3 = SimpleVariant::new(position, b"T", b"C");

    read_bases = &alt_bases[36..80];
    read_cigar = format!("{}M", read_bases.len());
    read = build_sam_record(36, &read_cigar, read_bases);

    let read_context = create_context(&var3, &read, 19, &ref_sequence, DEFAULT_FLANK_LENGTH).unwrap();
    let candidate3 =
        Candidate { tier: Tier::Panel, ctx: Arc::new(read_context), min_events: 1, full_match_support: 0, core_match_support: 0 };

    let mut var_context3 = to_context(CHR_1, &candidate3);

    position = 75;
    let var4 = SimpleVariant::new(position, b"C", b"G");

    read_bases = &alt_bases[56..99];
    read_cigar = format!("{}M", read_bases.len());
    read = build_sam_record(56, &read_cigar, read_bases);

    let read_context = create_context(&var4, &read, 19, &ref_sequence, DEFAULT_FLANK_LENGTH).unwrap();
    let candidate4 =
        Candidate { tier: Tier::Panel, ctx: Arc::new(read_context), min_events: 1, full_match_support: 0, core_match_support: 0 };

    let mut var_context4 = to_context(CHR_1, &candidate4);

    let lps_id1 = 1; // variants 1 & 2
    let lps_id2 = 2; // variants 1, 2 & 3
    let lps_id3 = 3; // variants 2, 3 & 4

    // the LOCAL_PHASE_SET attribute, as the VCF holds it
    let lps = |ids: &[i32]| ids.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    var_context1.info.insert("LPS".to_string(), lps(&[lps_id1, lps_id2]));
    var_context2.info.insert("LPS".to_string(), lps(&[lps_id1, lps_id2, lps_id3]));
    var_context3.info.insert("LPS".to_string(), lps(&[lps_id2, lps_id3]));
    var_context4.info.insert("LPS".to_string(), lps(&[lps_id3]));

    let rc_counter1 = create_read_counter(&candidate1.ctx, &fixture);
    let rc_counter2 = create_read_counter(&candidate2.ctx, &fixture);
    let rc_counter3 = create_read_counter(&candidate3.ctx, &fixture);
    let rc_counter4 = create_read_counter(&candidate4.ctx, &fixture);

    let candidates = vec![candidate1, candidate2, candidate3, candidate4];
    let variant_contexts = vec![var_context1, var_context2, var_context3, var_context4];

    variant_phaser.register_local_phase_sets(&candidates, &variant_contexts);

    // LPS #1
    // depth support only
    let mut alt_counters: Vec<&ReadContextCounter> = Vec::new();
    let mut ref_counters: Vec<&ReadContextCounter> = vec![&rc_counter1, &rc_counter2];

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);
    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    // alt support
    alt_counters = vec![&rc_counter1, &rc_counter2];
    ref_counters.clear();

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    // LPS #2

    // depth support
    ref_counters = vec![&rc_counter1, &rc_counter2, &rc_counter3];
    alt_counters.clear();

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    // alt support
    alt_counters = vec![&rc_counter1, &rc_counter2, &rc_counter3];
    ref_counters.clear();

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);
    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    // LPS #3

    // alt support
    alt_counters = vec![&rc_counter2, &rc_counter3, &rc_counter4];
    ref_counters.clear();

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);
    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    // subset of LPS 3
    ref_counters = vec![&rc_counter3, &rc_counter4];
    alt_counters.clear();

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    alt_counters = vec![&rc_counter3, &rc_counter4];
    ref_counters.clear();

    variant_phaser.registered_phased_variants(&alt_counters, &ref_counters);

    let lps_read_count_map = variant_phaser.sample_lps_read_counts(TEST_SAMPLE);
    assert!(lps_read_count_map.is_some());
    let lps_read_count_map = lps_read_count_map.unwrap();

    let lps_read_counts = lps_read_count_map.get(&lps_id1);
    assert!(lps_read_counts.is_some());
    let lps_read_counts = lps_read_counts.unwrap();
    assert!(lps_read_counts.subset_variant_counts.is_none());
    assert_eq!(3, lps_read_counts.depth);
    assert_eq!(1, lps_read_counts.alt_support);

    let lps_read_counts = lps_read_count_map.get(&lps_id2);
    assert!(lps_read_counts.is_some());
    let lps_read_counts = lps_read_counts.unwrap();
    assert!(lps_read_counts.subset_variant_counts.is_none());
    assert_eq!(3, lps_read_counts.depth);
    assert_eq!(2, lps_read_counts.alt_support);

    let lps_read_counts = lps_read_count_map.get(&lps_id3);
    assert!(lps_read_counts.is_some());
    let lps_read_counts = lps_read_counts.unwrap();
    assert!(lps_read_counts.subset_variant_counts.is_some());

    assert_eq!(2, lps_read_counts.depth);
    assert_eq!(2, lps_read_counts.alt_support);

    let sub_lps_read_counts = &lps_read_counts.subset_variant_counts.as_ref().unwrap()[0];
    assert_eq!(2, sub_lps_read_counts.depth);
    assert_eq!(1, sub_lps_read_counts.alt_support);
}
