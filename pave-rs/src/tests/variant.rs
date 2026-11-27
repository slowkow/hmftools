//! VariantTest

use super::*;
use hmf_common::vcf::{Genotype, VcfRecord};

use crate::annotation::pon::{PonData, PonFilters};
use crate::constants::{PON_ARTEFACT_FILTER, PON_FILTER, PON_GNOMAD_FILTER, PON_FILTERS_V38};
use crate::task::apply_filters_with;
use crate::variant::MultiPonStatus;

#[test]
fn test_point_mutation() {
    // SNVs and MNVs
    let var = VariantData::new(CHR_1, 100, "A", "G");
    assert!(var.is_base_change());
    assert_eq!(100, var.end_position);
    assert_eq!(1, var.alt_positions().len());
    assert!(var.alt_bases_above(99));
    assert!(!var.alt_bases_above(100));
    assert!(var.alt_bases_below(101));
    assert!(!var.alt_bases_below(100));
    assert!(var.alt_positions_overlap(90, 100));
    assert!(!var.alt_positions_overlap(90, 99));
    assert!(var.alt_positions_overlap(100, 110));
    assert!(!var.alt_positions_overlap(101, 110));
    assert!(var.alt_positions_within(100, 110));
    assert!(var.alt_positions_within(90, 100));

    let var = VariantData::new(CHR_1, 100, "AGA", "GCG");
    assert!(var.is_base_change());
    assert_eq!(102, var.end_position);
    assert_eq!(3, var.alt_positions().len());
    assert!(var.alt_bases_above(99));
    assert!(!var.alt_bases_above(100));
    assert!(!var.alt_bases_above(102));
    assert!(var.alt_bases_below(103));
    assert!(!var.alt_bases_below(102));
    assert!(!var.alt_bases_below(100));
    assert!(var.alt_positions_overlap(90, 100));
    assert!(!var.alt_positions_overlap(103, 110));
    assert!(var.alt_positions_within(100, 110));
    assert!(!var.alt_positions_within(102, 110));
    assert!(var.alt_positions_within(90, 102));
    assert!(!var.alt_positions_within(90, 101));
}

#[test]
fn test_deletion() {
    let var = VariantData::new(CHR_1, 100, "ACGT", "A");
    assert!(var.is_deletion());
    assert_eq!(104, var.end_position);
    assert_eq!(3, var.alt_positions().len());
    assert!(var.alt_bases_above(99));
    assert!(var.alt_bases_above(100));
    assert!(!var.alt_bases_above(101));
    assert!(var.alt_bases_below(104));
    assert!(!var.alt_bases_below(103));
    assert!(!var.alt_positions_overlap(90, 100));
    assert!(var.alt_positions_overlap(101, 110));
    assert!(!var.alt_positions_overlap(104, 110));
    assert!(var.alt_positions_overlap(103, 110));
    assert!(var.alt_positions_within(100, 110));
    assert!(var.alt_positions_within(101, 110));
    assert!(!var.alt_positions_within(102, 110));
    assert!(var.alt_positions_within(90, 104));
    assert!(var.alt_positions_within(90, 103));
    assert!(!var.alt_positions_within(90, 102));
}

#[test]
fn test_insertion() {
    let var = VariantData::new(CHR_1, 100, "A", "AAAAA"); // 4 bases
    assert!(var.is_insert());
    assert_eq!(101, var.end_position);
    assert!(var.alt_positions().is_empty());
    assert!(var.alt_bases_above(100));
    assert!(!var.alt_bases_above(101));
    assert!(var.alt_bases_below(101));
    assert!(!var.alt_bases_below(100));
    assert!(!var.alt_positions_overlap(90, 100));
    assert!(var.alt_positions_overlap(90, 101));
    assert!(!var.alt_positions_overlap(101, 110));
    assert!(var.alt_positions_overlap(100, 110));
    assert!(var.alt_positions_within(100, 110));
    assert!(!var.alt_positions_within(101, 110));
    assert!(var.alt_positions_within(90, 101));
    assert!(!var.alt_positions_within(90, 100));
}

const SAMPLE_ID: &str = "SAMPLE";

/// buildVariantContext: the sample with AF 0.1, the tier and the repeat
fn build_variant_context(variant: &VariantData, tier: &str, repeat_count: i32, repeat_bases: &str) -> VcfRecord {
    let mut rec = VcfRecord::new(&variant.chromosome, variant.position, ".", &variant.ref_, &variant.alt);
    let mut g = Genotype { sample: SAMPLE_ID.to_string(), gt: Some("./.".to_string()), ..Default::default() };
    g.ad = Some("10,10".to_string());
    g.dp = Some("100".to_string());
    g.attrs.push(("AF".to_string(), "0.1".to_string()));
    rec.genotypes.push(g);
    rec.put_text("TIER", tier.to_string());
    rec.put_text("REP_S", repeat_bases.to_string());
    rec.put_text("REP_C", repeat_count.to_string());
    rec
}

fn assert_filters(var: &VariantData, pon_standard: bool, pon_artefact: bool, gnomad: bool) {
    assert_eq!(pon_standard, var.filters.iter().any(|f| f == PON_FILTER));
    assert_eq!(pon_artefact, var.filters.iter().any(|f| f == PON_ARTEFACT_FILTER));
    assert_eq!(gnomad, var.filters.iter().any(|f| f == PON_GNOMAD_FILTER));
}

#[test]
fn test_pon_filters() {
    let mut var = VariantData::new(CHR_1, 100, "A", "T");

    // a hotspot variant first
    var.set_context(build_variant_context(&var, "HOTSPOT", 3, ""));

    let standard_pon = PonFilters::parse(PON_FILTERS_V38).unwrap();
    let mut artefacts = PonData::default();
    artefacts.add_entry(CHR_1, 1, "A", "C", 0, 0, 0, MultiPonStatus::Base);

    apply_filters_with(&mut var, SAMPLE_ID, &standard_pon, Some(&artefacts), true, false, false);
    assert_filters(&var, false, false, false);

    // gnomAD
    var.set_gnomad_frequency(Some(0.02));
    apply_filters_with(&mut var, SAMPLE_ID, &standard_pon, Some(&artefacts), true, false, false);
    assert_filters(&var, false, false, true);

    // standard PON
    var.set_pon_frequency(11, 7, 100, MultiPonStatus::Base);
    apply_filters_with(&mut var, SAMPLE_ID, &standard_pon, Some(&artefacts), true, false, false);
    assert_filters(&var, true, false, true);

    // artefact PON
    artefacts.add_entry(CHR_1, var.position, &var.ref_.clone(), &var.alt.clone(), 11, 11, 100, MultiPonStatus::Base);
    apply_filters_with(&mut var, SAMPLE_ID, &standard_pon, Some(&artefacts), true, false, false);
    assert_filters(&var, true, true, true);

    artefacts.clear();
    artefacts.add_entry(CHR_1, var.position, &var.ref_.clone(), &var.alt.clone(), 5, 11, 20, MultiPonStatus::Base);
    apply_filters_with(&mut var, SAMPLE_ID, &standard_pon, Some(&artefacts), true, false, false);
    assert_filters(&var, true, false, true);

    // a non-hotspot variant
    var.set_context(build_variant_context(&var, "PANEL", 3, ""));
    var.set_gnomad_frequency(Some(0.001));
    apply_filters_with(&mut var, SAMPLE_ID, &standard_pon, Some(&artefacts), true, false, false);
    assert_filters(&var, true, true, true);
}
