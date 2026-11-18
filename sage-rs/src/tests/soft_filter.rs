//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! filters/SoftFilterTest.java
//!
//! Three Java tests set SageConfig.HighDepthMode, which the Rust port does not implement. With these
//! inputs it changes nothing that is asserted: no read is chimeric, every base quality is high, and
//! its only other effect is to lift the cap on the edge distance penalty in the map qual factor, which
//! feeds minMQF only. The tests are therefore ported as they are, without the mode.
//!
//! A Java SageVariant shares its counters with the test; the Rust one owns them, so where a test goes
//! on using a counter after filtering, the counter is taken back out of the variant.

use super::phasing_groups::{create_read_context_cores, create_read_context_flanks, TEST_LEFT_FLANK, TEST_RIGHT_FLANK};
use super::util::*;
use crate::candidate::Candidate;
use crate::constants::*;
use crate::evidence::{FragmentCoords, FragmentData, ReadContextCounter};
use crate::filters::{SageVariant, SoftFilter, VariantFilters};
use crate::read_context::VariantReadContext;
use crate::types::*;
use std::sync::Arc;

const TEST_READ_ID: &str = "READ_01";

/// MockRefGenome.generateRandomBases: not random, it cycles through the nucleotides
fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect()
}

/// MockRefGenome.getNextBase
fn get_next_base(base: &str) -> String {
    match base {
        "A" => "C",
        "C" => "G",
        "G" => "T",
        "T" => "A",
        other => other,
    }
    .to_string()
}

/// REF_BASES = "X" + generateRandomBases(100)
fn ref_bases() -> String {
    format!("X{}", generate_random_bases(100))
}

/// REF_BASES.substring(start, end)
fn sub(start: usize, end: usize) -> String {
    ref_bases()[start..end].to_string()
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

/// FILTERS.applySoftFilters(variant)
fn apply_soft_filters(fixture: &Fixture, variant: &mut SageVariant) {
    let hla = |_: i32| false;
    VariantFilters { cfg: &fixture.cfg, chr: CHR_1, hla: &hla }.apply_soft_filters(variant);
}

/// VariantUtils.createReadCounter(readContext, TEST_CONFIG, tier)
fn create_read_counter_tier(fixture: &Fixture, ctx: &Arc<VariantReadContext>, tier: Tier) -> ReadContextCounter {
    ReadContextCounter::new(ctx.clone(), tier, DEFAULT_MAX_READ_DEPTH, 1, &fixture.qc(), false)
}

/// VariantUtils.createReadCounter(readContext, isReference)
fn create_read_counter(fixture: &Fixture, ctx: &Arc<VariantReadContext>, is_reference: bool) -> ReadContextCounter {
    ReadContextCounter::new(ctx.clone(), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &fixture.qc(), is_reference)
}

/// VariantUtils.sageVariantFromReadContextCounter
fn sage_variant_from_read_context_counter(counter: ReadContextCounter) -> SageVariant {
    let candidate =
        Candidate { tier: Tier::HighConfidence, ctx: counter.ctx.clone(), min_events: 1, full_match_support: 1, core_match_support: 0 };
    SageVariant::new(candidate, Vec::new(), vec![counter])
}

/// new Candidate(tier, readContext, 1, 1)
fn candidate(tier: Tier, ctx: &Arc<VariantReadContext>) -> Candidate {
    Candidate { tier, ctx: ctx.clone(), min_events: 1, full_match_support: 1, core_match_support: 0 }
}

/// createSnvReadContextCounter
fn create_snv_read_context_counter(fixture: &Fixture, variant_position: usize) -> ReadContextCounter {
    let ref_base = sub(variant_position, variant_position + 1);
    let alt_base = get_next_base(&ref_base);

    let variant = SimpleVariant::new(variant_position as i32, ref_base.as_bytes(), alt_base.as_bytes());

    let left_core = sub(variant_position - MIN_CORE_DISTANCE as usize, variant_position);
    let right_core = sub(variant_position + 1, variant_position + 1 + MIN_CORE_DISTANCE as usize);

    let ctx = create_read_context_cores(variant, &left_core, &right_core);

    create_read_counter_tier(fixture, &ctx, Tier::Panel)
}

/// FragmentCoords.lowerCount / upperCount / atCapacity
fn lower_count(coords: &FragmentCoords) -> usize {
    coords.lower.len()
}
fn upper_count(coords: &FragmentCoords) -> usize {
    coords.upper.len()
}
fn at_capacity(coords: &FragmentCoords) -> bool {
    coords.lower.len() >= coords.max && coords.upper.len() >= coords.max
}

/// new FragmentData(first, second)
fn fragment_data(first: &Read, second: &Read) -> FragmentData {
    FragmentData { first: first.clone(), second: second.clone() }
}

// Java: testTumorQualFilter
#[test]
fn test_tumor_qual_filter() {
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 50;

    let ctx = create_read_context_flanks(SimpleVariant::new(position as i32, b"A", b"T"), &sub(48, 50), &sub(51, 53), &sub(38, 48), &sub(53, 63));

    // SageConfig.HighDepthMode = true: see the module comment
    let mut read_context_counter = create_read_counter_tier(&fixture, &ctx, Tier::Panel);

    let alt_base = read_context_counter.variant().alt_str().to_string();

    let alt_read_bases = format!("{}{}{}", sub(20, position), alt_base, sub(position + 1, 80));
    let read_cigar = build_cigar_string(alt_read_bases.len());

    let alt_read = create_sam_record(&next_read_id(), 20, &alt_read_bases, &read_cigar);

    let ref_read = create_sam_record(&next_read_id(), 20, &sub(20, 80), &read_cigar);

    read_context_counter.process_read(&alt_read, 1, None, &qc);
    read_context_counter.process_read(&alt_read, 1, None, &qc);
    read_context_counter.process_read(&ref_read, 1, None, &qc);
    read_context_counter.process_read(&ref_read, 1, None, &qc);

    let mut variant = sage_variant_from_read_context_counter(read_context_counter);
    apply_soft_filters(&fixture, &mut variant);

    assert!(variant.filters.contains(&SoftFilter::MinTumorQual));

    let mut read_context_counter = variant.tumor_counters.pop().unwrap();
    read_context_counter.process_read(&alt_read, 1, None, &qc);
    read_context_counter.process_read(&alt_read, 1, None, &qc);

    let mut variant = sage_variant_from_read_context_counter(read_context_counter);
    apply_soft_filters(&fixture, &mut variant);

    assert!(!variant.filters.contains(&SoftFilter::MinTumorQual));
}

// Java: testMaxReadEdgeDistanceFilter
#[test]
fn test_max_read_edge_distance_filter() {
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 50;

    // SageConfig.HighDepthMode = true: see the module comment
    let mut read_context_counter = create_snv_read_context_counter(&fixture, position);

    // all reads have the variant near the end of the read, and so fail the max edge distance filter

    let alt_base = read_context_counter.variant().alt_str().to_string();

    let read = create_sam_record(TEST_READ_ID, 45, &format!("{}{}{}", sub(45, position), alt_base, sub(position + 1, 95)), "50M");

    read_context_counter.process_read(&read, 1, None, &qc);
    read_context_counter.process_read(&read, 1, None, &qc);
    read_context_counter.process_read(&read, 1, None, &qc);

    let read2 = create_sam_record(TEST_READ_ID, 45, &format!("{}{}{}", sub(42, position), alt_base, sub(position + 1, 95)), "50M");

    read_context_counter.process_read(&read2, 1, None, &qc);
    read_context_counter.process_read(&read2, 1, None, &qc);

    // factor in soft-clipped bases and use a FULL match rather than CORE
    let read_bases = format!(
        "{}{}{}{}{}",
        TEST_LEFT_FLANK,
        sub(position - 2, position),
        alt_base,
        sub(position + 1, position + 3),
        TEST_RIGHT_FLANK.repeat(5)
    );

    let read3 = create_sam_record(TEST_READ_ID, 46, &read_bases, "8S13M24S");

    read_context_counter.process_read(&read3, 1, None, &qc);

    // a read supporting the ref
    let read4 = create_sam_record(TEST_READ_ID, 1, &sub(1, 99), "98M"); // was 98
    read_context_counter.process_read(&read4, 1, None, &qc);

    assert!((0.06 - read_context_counter.read_edge_distance.max_alt_distance_from_edge).abs() <= 0.01);
    assert!((0.48 - read_context_counter.read_edge_distance.max_distance_from_edge).abs() <= 0.01);

    let mut variant = sage_variant_from_read_context_counter(read_context_counter);

    apply_soft_filters(&fixture, &mut variant);

    assert!(variant.filters.contains(&SoftFilter::MaxEdgeDistance));
}

// Java: testFragmentCoordsFilter
#[test]
fn test_fragment_coords_filter() {
    // variants are filtered if the read coords supporting the alt are not sufficiently unique
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 20;

    // SageConfig.HighDepthMode = true: see the module comment
    let mut read_context_counter = create_snv_read_context_counter(&fixture, position);
    let alt_base = read_context_counter.variant().alt_str().to_string();

    // 2 reads at the start, 3 at the end is insufficient
    let cigar1 = "2S25M3S";

    let mut read1 =
        create_sam_record(TEST_READ_ID, 12, &format!("{}{}{}", sub(10, position), alt_base, sub(position + 1, 40)), cigar1);
    read1.mate_tid = read1.tid;
    read1.mate_start = 100;
    read1.mate_cigar = Some(cigar1.to_string());
    read1.flags |= FLAG_MATE_REVERSE;

    // same frag coordinates processed twice
    read_context_counter.process_read(&read1, 1, None, &qc);
    read_context_counter.process_read(&read1, 1, None, &qc);

    assert_eq!(1, lower_count(&read_context_counter.fragment_coords));
    assert_eq!(1, upper_count(&read_context_counter.fragment_coords));

    // same fragment start (unclipped), different mate end
    let read2 = create_sam_record(TEST_READ_ID, 10, &format!("{}{}{}", sub(10, position), alt_base, sub(position + 1, 40)), "30M");

    let mut mate2 = create_sam_record(TEST_READ_ID, 14, &format!("{}{}{}", sub(14, position), alt_base, sub(position + 1, 44)), "30M");

    mate2.flags |= FLAG_REVERSE;

    read_context_counter.process_read(&read2, 1, Some(&fragment_data(&read2, &mate2)), &qc);
    read_context_counter.process_read(&read2, 1, Some(&fragment_data(&read2, &mate2)), &qc);
    read_context_counter.process_read(&read2, 1, Some(&fragment_data(&read2, &mate2)), &qc);
    read_context_counter.process_read(&read2, 1, Some(&fragment_data(&read2, &mate2)), &qc);

    assert_eq!(1, lower_count(&read_context_counter.fragment_coords));
    assert_eq!(2, upper_count(&read_context_counter.fragment_coords));

    let mut variant = sage_variant_from_read_context_counter(read_context_counter);

    apply_soft_filters(&fixture, &mut variant);

    assert!(variant.filters.contains(&SoftFilter::FragmentCoords));

    let mut read_context_counter = variant.tumor_counters.pop().unwrap();

    // different frag coords pass the filter
    let read3 = create_sam_record(TEST_READ_ID, 10, &format!("{}{}{}", sub(8, position), alt_base, sub(position + 1, 42)), "2S30M");

    let mut mate3 =
        create_sam_record(TEST_READ_ID, 14, &format!("{}{}{}", sub(14, position), alt_base, sub(position + 1, 46)), "30M2S");

    mate3.flags |= FLAG_REVERSE;

    read_context_counter.process_read(&read3, 1, Some(&fragment_data(&read3, &mate3)), &qc);
    read_context_counter.process_read(&read3, 1, Some(&fragment_data(&read3, &mate3)), &qc);

    assert_eq!(2, lower_count(&read_context_counter.fragment_coords));
    assert_eq!(3, upper_count(&read_context_counter.fragment_coords));

    let read4 = create_sam_record(TEST_READ_ID, 11, &format!("{}{}{}", sub(11, position), alt_base, sub(position + 1, 40)), "30M");
    // as in Java, the mate details are set again on read1, not on read4
    read1.mate_tid = read1.tid;
    read1.mate_start = 100;
    read1.mate_cigar = Some(cigar1.to_string());
    read1.flags |= FLAG_MATE_REVERSE;

    read_context_counter.process_read(&read4, 1, None, &qc);

    assert!(at_capacity(&read_context_counter.fragment_coords));
}

// Java: testGermlineRelativeQualFilter
#[test]
fn test_germline_relative_qual_filter() {
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 50;

    let ref_base = sub(position, position + 1);
    let alt_base = get_next_base(&ref_base);
    let variant = SimpleVariant::new(position as i32, ref_base.as_bytes(), alt_base.as_bytes());

    let ctx = create_read_context_flanks(variant, &sub(48, position), &sub(position + 1, 53), &sub(38, 48), &sub(53, 63));

    let mut ref_counter = create_read_counter(&fixture, &ctx, true);
    let mut tumor_counter = create_read_counter(&fixture, &ctx, false);

    let read_pos_start = 20;

    let alt_read_bases = format!("{}{}{}", sub(read_pos_start, position), alt_base, sub(position + 1, 80));
    let read_cigar = build_cigar_string(alt_read_bases.len());

    let alt_read = create_sam_record(&next_read_id(), read_pos_start as i32, &alt_read_bases, &read_cigar);
    tumor_counter.process_read(&alt_read, 1, None, &qc);

    ref_counter.process_read(&alt_read, 1, None, &qc);

    let ref_bases = sub(read_pos_start, 80);
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let ref_read = create_sam_record(&next_read_id(), read_pos_start as i32, &ref_bases, &read_cigar);

    for _ in 0..12 {
        ref_counter.process_read(&ref_read, 1, None, &qc);
    }

    assert_eq!(1, ref_counter.alt_support());
    assert_eq!(13, ref_counter.depth());

    let mut snv = SageVariant::new(candidate(Tier::Panel, &ctx), vec![ref_counter], vec![tumor_counter]);

    apply_soft_filters(&fixture, &mut snv);

    assert!(!snv.filters.contains(&SoftFilter::MaxGermlineRelativeQual));

    let mut ref_counter = snv.ref_counters.pop().unwrap();
    let tumor_counter = snv.tumor_counters.pop().unwrap();

    // now exceeding the limit
    ref_counter.process_read(&alt_read, 1, None, &qc);
    ref_counter.process_read(&alt_read, 1, None, &qc);

    let mut snv = SageVariant::new(candidate(Tier::Panel, &ctx), vec![ref_counter], vec![tumor_counter]);

    apply_soft_filters(&fixture, &mut snv);

    assert!(snv.filters.contains(&SoftFilter::MaxGermlineRelativeQual));
}

// Java: testSimpleAtGermlineVaf
#[test]
fn test_simple_at_germline_vaf() {
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 50;

    let variant = SimpleVariant::new(position as i32, b"ACG", b"TCA");
    let alt_bases = variant.alt_str().to_string();
    let ctx = create_read_context_flanks(variant, &sub(48, position), &sub(53, 55), &sub(38, 48), &sub(55, 63));

    let mut ref_counter = create_read_counter(&fixture, &ctx, true);
    let mut tumor_counter = create_read_counter(&fixture, &ctx, false);

    let read_pos_start = 20;

    // trigger alt in the germline and and full support tumor reads to observe the germline filter still being applied
    let alt_read_bases = format!("{}{}{}", sub(read_pos_start, position), alt_bases, sub(position + 3, 80));
    let read_cigar = build_cigar_string(alt_read_bases.len());

    let alt_read = create_sam_record(&next_read_id(), read_pos_start as i32, &alt_read_bases, &read_cigar);
    tumor_counter.process_read(&alt_read, 1, None, &qc);

    let simple_alt_read_bases = format!("{}G{}{}", sub(read_pos_start, position - 1), alt_bases, sub(position + 3, 80));
    let simple_alt_read = create_sam_record(&next_read_id(), 20, &simple_alt_read_bases, &read_cigar);
    ref_counter.process_read(&simple_alt_read, 1, None, &qc);

    assert_eq!(1, ref_counter.simple_alt_matches);
    assert_eq!(0, ref_counter.alt_support());

    let mut mnv = SageVariant::new(candidate(Tier::HighConfidence, &ctx), vec![ref_counter], vec![tumor_counter]);

    apply_soft_filters(&fixture, &mut mnv);

    assert!(mnv.filters.contains(&SoftFilter::MaxGermlineVaf));
}

// Java: testPanelAdjustedGermlineVaf
#[test]
fn test_panel_adjusted_germline_vaf() {
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let position = 50;

    let ref_base = sub(position, position + 1);
    let alt_base = get_next_base(&ref_base);
    let variant = SimpleVariant::new(position as i32, ref_base.as_bytes(), alt_base.as_bytes());

    let ctx = create_read_context_flanks(variant, &sub(48, position), &sub(position + 1, 53), &sub(38, 48), &sub(53, 63));

    let mut ref_counter = create_read_counter(&fixture, &ctx, true);
    let mut tumor_counter = create_read_counter(&fixture, &ctx, false);

    let read_pos_start = 20;

    let alt_read_bases = format!("{}{}{}", sub(read_pos_start, position), alt_base, sub(position + 1, 80));
    let read_cigar = build_cigar_string(alt_read_bases.len());

    let mut alt_read = create_sam_record(&next_read_id(), read_pos_start as i32, &alt_read_bases, &read_cigar);
    tumor_counter.process_read(&alt_read, 1, None, &qc);

    alt_read.quals[30] = 10;
    ref_counter.process_read(&alt_read, 1, None, &qc);

    let ref_bases = sub(read_pos_start, 80);
    let read_cigar = build_cigar_string(alt_read_bases.len());
    let ref_read = create_sam_record(&next_read_id(), read_pos_start as i32, &ref_bases, &read_cigar);

    for _ in 0..12 {
        ref_counter.process_read(&ref_read, 1, None, &qc);
    }

    assert_eq!(1, ref_counter.alt_support());
    assert_eq!(13, ref_counter.depth());

    let mut snv = SageVariant::new(candidate(Tier::Panel, &ctx), vec![ref_counter], vec![tumor_counter]);

    apply_soft_filters(&fixture, &mut snv);

    assert!(!snv.filters.contains(&SoftFilter::MaxGermlineVaf));
}
