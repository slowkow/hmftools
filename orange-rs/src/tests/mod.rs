//! Ported Java unit tests of ORANGE (orange/src/test/java, the classes ported: algo, conversion, util) and tests of
//! the port. Synthetic fixtures; the CUPPA tests read the Java tests' resource files (ORANGE_TEST_RESOURCES, or
//! ../../../../hmftools/orange/src/test/resources) and are skipped without them.

use hmf_common::cytobands::CytoBands;
use hmf_common::linx_files::{LinxBreakend, LinxSvAnnotation};
use hmf_common::purple_files::{ChrArmCopyNumber, GeneCopyNumber, GermlineAmpDel, PurplePurity, PurpleQc};
use hmf_common::rna;
use hmf_common::small_variant::determine_coding_effect;
use hmf_common::tool_results as tr;

use crate::compat::with_java;
use crate::linx::LinxData;
use crate::model::*;
use crate::purple::PurpleData;

const EPSILON: f64 = 1e-2;

fn resources() -> Option<String> {
    let dir = std::env::var("ORANGE_TEST_RESOURCES")
        .unwrap_or_else(|_| concat!(env!("CARGO_MANIFEST_DIR"), "/../../../../hmftools/orange/src/test/resources").to_string());
    if std::path::Path::new(&dir).exists() { Some(dir) } else { None }
}

// ---------------------------------------------------------------------------------------------
// test factories (GermlineDeletionTestFactory, GeneCopyNumberTestFactory, TestPurpleGainDeletionFactory, LinxTestFactory)

fn germline_deletion(gene: &str, reported: bool, status: &str, tumor_cn: f64, start: i32, end: i32) -> GermlineAmpDel {
    GermlineAmpDel {
        gene_name: gene.into(),
        transcript: String::new(),
        chromosome: "1".into(),
        chromosome_band: "p36.13".into(),
        region_start: start,
        region_end: end,
        depth_window_count: 0,
        exon_start: 0,
        exon_end: 0,
        is_partial: true,
        detection_method: "SEGMENT".into(),
        normal_status: status.into(),
        tumor_status: status.into(),
        germline_copy_number: 0.0,
        tumor_copy_number: tumor_cn,
        filter: String::new(),
        cohort_frequency: 0,
        reported: if reported { "REPORTED" } else { "NONE" }.into(),
    }
}

fn gene_copy_number(gene: &str, min: f64, max: f64) -> GeneCopyNumber {
    GeneCopyNumber {
        chromosome: "1".into(),
        start: 0,
        end: 0,
        gene_name: gene.into(),
        trans_name: String::new(),
        is_canonical: true,
        chromosome_band: "p36.13".into(),
        max_copy_number: max,
        min_copy_number: min,
        min_minor_allele_copy_number: 0.0,
        somatic_regions: 1,
        min_regions: 1,
        min_region_start: 0,
        min_region_end: 0,
        depth_window_count: 0,
        gc_content: 1.0,
        min_region_start_support: "NONE".into(),
        min_region_end_support: "NONE".into(),
        min_region_method: "UNKNOWN".into(),
        relative_min_copy_number: min / 2.0,
        reported_status: "NONE".into(),
        driver_type: "UNKNOWN".into(),
    }
}

fn arm(chromosome: &str, arm: &str) -> ChrArmCopyNumber {
    ChrArmCopyNumber {
        chromosome: chromosome.into(),
        arm: arm.into(),
        mean_copy_number: 2.0,
        median_copy_number: 2.0,
        min_copy_number: 2.0,
        max_copy_number: 2.0,
    }
}

fn driver(gene: &str, driver_type: &str) -> PurpleDriver {
    PurpleDriver {
        gene: gene.into(),
        transcript: String::new(),
        driver_type: driver_type.into(),
        driver_likelihood: 0.0,
        likelihood_method: "DEL".into(),
        reported_status: "REPORTED".into(),
        driver_interpretation: "HIGH".into(),
        is_canonical: true,
        category: "TSG".into(),
    }
}

fn gain_del(gene: &str, driver_type: &str, canonical: bool) -> PurpleGainDeletion {
    let mut d = driver(gene, driver_type);
    d.is_canonical = canonical;
    PurpleGainDeletion {
        driver: d,
        germline_amp_del_fields: None,
        chromosome: String::new(),
        chromosome_band: String::new(),
        gene_range: String::new(),
        exon_start: None,
        exon_end: None,
        min_copy_number: 0.0,
        max_copy_number: 0.0,
        relative_copy_number: 0.0,
        arm_copy_number: 0.0,
        min_minor_allele_copies: 0.0,
        tpm: None,
        tpm_percentile: None,
        tpm_fold_change: None,
        plot_filename: None,
    }
}

fn minimal_qc() -> PurpleQc {
    PurpleQc {
        status: vec!["PASS".into()],
        method: "NORMAL".into(),
        copy_number_segments: 0,
        unsupported_copy_number_segments: 0,
        purity: 0.0,
        amber_gender: "FEMALE".into(),
        cobalt_gender: "FEMALE".into(),
        deleted_genes: 0,
        contamination: 0.0,
        germline_aberrations: vec!["NONE".into()],
        amber_mean_depth: 0,
        loh_percent: 0.0,
        tinc_level: 0.0,
    }
}

fn minimal_purity() -> PurplePurity {
    PurplePurity {
        purity: 0.5,
        norm_factor: 1.0,
        score: 0.0,
        diploid_proportion: 1.0,
        ploidy: 2.0,
        sex: "FEMALE".into(),
        fit_method: "NORMAL".into(),
        polyclonal_proportion: 0.0,
        min_purity: 0.5,
        max_purity: 0.5,
        min_ploidy: 2.0,
        max_ploidy: 2.0,
        min_diploid_proportion: 1.0,
        max_diploid_proportion: 1.0,
        somatic_penalty: 0.0,
        whole_genome_duplication: false,
        ms_indels_per_mb: 0.0,
        ms_status: "MSS".into(),
        tml: 0,
        tml_status: "LOW".into(),
        tmb_per_mb: 0.0,
        tmb_status: "LOW".into(),
        sv_tumor_mutational_burden: 0,
        run_mode: "TUMOR_GERMLINE".into(),
        targeted: false,
    }
}

/// PurpleTestFactory.createMinimalTestPurpleData
fn minimal_purple_data() -> PurpleData {
    PurpleData {
        qc: minimal_qc(),
        purity: minimal_purity(),
        somatic_drivers: Vec::new(),
        germline_drivers: None,
        somatic_variants: Vec::new(),
        somatic_transcript_impacts: Vec::new(),
        germline_variants: None,
        germline_transcript_impacts: Vec::new(),
        somatic_gene_copy_numbers: Vec::new(),
        germline_amp_dels: None,
        chr_arm_copy_numbers: Vec::new(),
        variant_plots: Vec::new(),
    }
}

fn minimal_linx_data() -> LinxData {
    LinxData {
        somatic_sv_annotations: Vec::new(),
        somatic_driver_data: Vec::new(),
        somatic_drivers: Vec::new(),
        fusions: Vec::new(),
        somatic_breakends: Vec::new(),
        germline_drivers: None,
        germline_breakends: None,
        germline_sv_annotations: None,
        reportable_event_plots: Vec::new(),
        reportable_event_plot_titles: Vec::new(),
    }
}

fn driver_catalog(gene: &str, driver_type: &str) -> hmf_common::purple_files::DriverCatalog {
    hmf_common::purple_files::DriverCatalog {
        chromosome: "1".into(),
        chromosome_band: "p36.13".into(),
        gene: gene.into(),
        transcript: String::new(),
        is_canonical: true,
        driver: driver_type.into(),
        category: "TSG".into(),
        likelihood_method: "DEL".into(),
        reported_status: "REPORTED".into(),
        driver_likelihood: 1.0,
        missense: 0,
        nonsense: 0,
        splice: 0,
        inframe: 0,
        frameshift: 0,
        biallelic: false,
        min_copy_number: 0.0,
        max_copy_number: 0.0,
    }
}

// ---------------------------------------------------------------------------------------------
// GermlineGainDeletionFactoryTest

#[test]
fn can_transform_reportable_hom_deletions_to_partial() {
    let dels = vec![
        germline_deletion("gene", true, "HOM_DELETION", 0.0, 400, 700),
        germline_deletion("gene", true, "HOM_DELETION", 0.0, 800, 1000),
    ];
    let gcn = gene_copy_number("gene", 1.0, 4.0);
    let g = crate::purple::germline_gains_dels(&dels, &[driver("gene", "GERMLINE_DELETION")], &[gcn], &[arm("1", "P")], None).unwrap();
    assert_eq!(g.len(), 1);
    assert_eq!(g[0].driver.reported_status, "REPORTED");
    assert_eq!(g[0].driver.gene, "gene");
    assert!((g[0].min_copy_number - 0.0).abs() < EPSILON);
    assert!((g[0].max_copy_number - 4.0).abs() < EPSILON);
}

#[test]
fn can_transform_reportable_hom_deletion_to_full() {
    let dels = vec![germline_deletion("gene", true, "HOM_DELETION", 0.0, 100, 1200)];
    let g = crate::purple::germline_gains_dels(&dels, &[driver("gene", "GERMLINE_DELETION")], &[gene_copy_number("gene", 1.0, 1.0)], &[arm("1", "P")], None)
        .unwrap();
    assert_eq!(g.len(), 1);
    assert!((g[0].min_copy_number - 0.0).abs() < EPSILON);
    assert!((g[0].max_copy_number - 1.0).abs() < EPSILON);
}

#[test]
fn can_transform_reportable_partial_hom_deletions_to_full_gene_loss() {
    let dels = vec![
        germline_deletion("gene", true, "HOM_DELETION", 0.1, 200, 300),
        germline_deletion("gene", true, "HOM_DELETION", 0.0, 300, 500),
        germline_deletion("gene", true, "HOM_DELETION", 0.0, 500, 800),
        germline_deletion("gene", true, "HOM_DELETION", 0.2, 700, 2000),
    ];
    let g = crate::purple::germline_gains_dels(&dels, &[driver("gene", "GERMLINE_DELETION")], &[gene_copy_number("gene", 1.0, 4.0)], &[arm("1", "P")], None)
        .unwrap();
    assert_eq!(g.len(), 1);
    assert!((g[0].min_copy_number - 0.0).abs() < EPSILON);
    assert!((g[0].max_copy_number - 4.0).abs() < EPSILON);
}

#[test]
fn germline_gain_deletions_need_a_reported_germline_driver_and_gene_copy_number() {
    let dels = vec![germline_deletion("gene", true, "HOM_DELETION", 0.0, 100, 1200)];
    let mut not_reported = driver("gene", "GERMLINE_DELETION");
    not_reported.reported_status = "NOT_REPORTED".into();
    let gcn = vec![gene_copy_number("gene", 1.0, 1.0)];
    let arms = vec![arm("1", "P")];
    assert!(crate::purple::germline_gains_dels(&dels, &[not_reported], &gcn, &arms, None).unwrap().is_empty());
    assert!(crate::purple::germline_gains_dels(&dels, &[driver("gene", "GERMLINE_MUTATION")], &gcn, &arms, None).unwrap().is_empty());
    assert!(crate::purple::germline_gains_dels(&dels, &[driver("gene", "GERMLINE_DELETION")], &[], &arms, None).unwrap().is_empty());
    // an arm missing from the file is one of copy number 0 (short arms)
    let g = crate::purple::germline_gains_dels(&dels, &[driver("gene", "GERMLINE_DELETION")], &gcn, &[], None).unwrap();
    assert_eq!(g[0].arm_copy_number, 0.0);
    assert_eq!(g[0].gene_range, "PARTIAL");
}

// ---------------------------------------------------------------------------------------------
// PurpleInterpreterTest

#[test]
fn can_interpret_minimal_purple_data() {
    let r = crate::purple::interpret(&minimal_purple_data(), &minimal_linx_data(), None).unwrap();
    assert!(r.germline_drivers.is_none() && r.germline_variants.is_none() && r.germline_gains_dels.is_none());
    assert_eq!(r.fit.qc.status, vec!["PASS"]);
}

fn purple_test_data(all: Vec<GermlineAmpDel>) -> PurpleData {
    let mut d = minimal_purple_data();
    d.somatic_gene_copy_numbers.push(gene_copy_number("gene", 0.0, 0.0));
    d.chr_arm_copy_numbers.push(arm("1", "P"));
    d.germline_amp_dels = Some(all.into_iter().filter(|x| x.reported == "REPORTED").collect());
    d.germline_drivers = Some(vec![driver_catalog("gene", "GERMLINE_DELETION")]);
    d.germline_variants = Some(Vec::new());
    d
}

#[test]
fn can_create_reportable_germline_full_dels() {
    let d = purple_test_data(vec![
        germline_deletion("gene", true, "HET_DELETION", 1.0, 0, 0),
        germline_deletion("gene", true, "HOM_DELETION", 2.0, 0, 0),
    ]);
    let r = crate::purple::interpret(&d, &minimal_linx_data(), None).unwrap();
    assert_eq!(r.germline_gains_dels.unwrap().len(), 1);
}

#[test]
fn can_handle_half_reportable_germline_full_dels() {
    let d = purple_test_data(vec![
        germline_deletion("gene", false, "HET_DELETION", 1.0, 0, 0),
        germline_deletion("gene", true, "HOM_DELETION", 2.0, 0, 0),
    ]);
    let r = crate::purple::interpret(&d, &minimal_linx_data(), None).unwrap();
    assert_eq!(r.germline_gains_dels.unwrap().len(), 1);
}

#[test]
fn arm_aberrations_follow_ploidy_and_gender() {
    let mut d = minimal_purple_data();
    let mut gain = arm("1", "Q");
    gain.mean_copy_number = 3.0; // > 1.25 x 2: GAIN, > 1.4 x 2: HIGH
    gain.median_copy_number = 2.9;
    let mut loss = arm("2", "P");
    loss.mean_copy_number = 1.4; // < 0.75 x 2: LOSS, not < 0.6 x 2: LOW
    let mut y = arm("Y", "Q");
    y.mean_copy_number = 0.0;
    d.chr_arm_copy_numbers = vec![gain, loss, arm("3", "P"), y.clone()];
    let r = crate::purple::interpret(&d, &minimal_linx_data(), None).unwrap();
    let a = &r.arm_copy_number_abberations;
    assert_eq!(a.len(), 2); // Y is left out for a female sample
    assert_eq!((a[0].chromosome.as_str(), a[0].cn_type.as_str(), a[0].driver_interpretation.as_str()), ("chr1", "GAIN", "HIGH"));
    assert_eq!(a[0].copy_number, 2.9); // the median, while type and relative copy number use the mean
    assert_eq!(a[0].relative_copy_number, 1.5);
    assert_eq!((a[1].cn_type.as_str(), a[1].driver_interpretation.as_str()), ("LOSS", "LOW"));
    d.purity.sex = "MALE".into();
    let mut x = arm("X", "P");
    x.mean_copy_number = 1.0; // the expected ploidy of an allosome of a male sample is halved
    d.chr_arm_copy_numbers = vec![x, y];
    let r = crate::purple::interpret(&d, &minimal_linx_data(), None).unwrap();
    assert_eq!(r.arm_copy_number_abberations.len(), 1);
    assert_eq!(r.arm_copy_number_abberations[0].chromosome, "chrY");
}

#[test]
fn somatic_gains_dels_from_amp_and_del_drivers() {
    let mut d = minimal_purple_data();
    let mut amp = driver_catalog("AMPGENE", "AMP");
    amp.min_copy_number = -0.5;
    amp.max_copy_number = 9.0;
    d.somatic_drivers = vec![amp, driver_catalog("DELGENE", "DEL"), driver_catalog("MUTGENE", "MUTATION"), driver_catalog("NOCN", "DEL")];
    let mut g1 = gene_copy_number("AMPGENE", 8.0, 9.0);
    g1.relative_min_copy_number = -1.0; // pre-4.4 PURPLE: min copy number over ploidy
    d.somatic_gene_copy_numbers = vec![g1, gene_copy_number("DELGENE", 0.0, 1.0)];
    let r = crate::purple::interpret(&d, &minimal_linx_data(), None).unwrap();
    let g = &r.somatic_gains_dels;
    assert_eq!(g.len(), 2);
    assert_eq!((g[0].gene_range.as_str(), g[0].min_copy_number, g[0].relative_copy_number), ("FULL", 0.0, 4.0));
    assert_eq!(g[1].gene_range, "PARTIAL");
}

// ---------------------------------------------------------------------------------------------
// CodingEffectDeterminerTest

#[test]
fn coding_effect_determination() {
    assert_eq!(determine_coding_effect(&[]), "NONE");
    assert_eq!(determine_coding_effect(&["START_LOST"]), "NONSENSE_OR_FRAMESHIFT");
    assert_eq!(determine_coding_effect(&["SPLICE_DONOR"]), "SPLICE");
    assert_eq!(determine_coding_effect(&["MISSENSE"]), "MISSENSE");
    assert_eq!(determine_coding_effect(&["SYNONYMOUS"]), "SYNONYMOUS");
    assert_eq!(determine_coding_effect(&["INTRONIC"]), "NONE");
    assert_eq!(determine_coding_effect(&["START_LOST", "SPLICE_DONOR"]), "NONSENSE_OR_FRAMESHIFT");
    assert_eq!(determine_coding_effect(&["MISSENSE", "SPLICE_DONOR"]), "SPLICE");
    assert_eq!(determine_coding_effect(&["MISSENSE", "SYNONYMOUS"]), "MISSENSE");
    assert_eq!(determine_coding_effect(&["INTRONIC", "SYNONYMOUS"]), "SYNONYMOUS");
}

// ---------------------------------------------------------------------------------------------
// LinxBreakendInterpreterTest

fn sv_annotation() -> LinxSvAnnotation {
    LinxSvAnnotation {
        vcf_id_start: "vcfId".into(),
        vcf_id_end: "vcfId".into(),
        sv_id: 1,
        coords_start: "1:1000:1".into(),
        coords_end: "1:2000:-1".into(),
        sv_type: "DEL".into(),
        cluster_id: 0,
        cluster_reason: String::new(),
        fragile_site_start: false,
        fragile_site_end: false,
        is_foldback: false,
        line_type_start: String::new(),
        line_type_end: String::new(),
        junction_copy_number_min: 1.0,
        junction_copy_number_max: 2.0,
        gene_start: String::new(),
        gene_end: String::new(),
        local_topology_id_start: 0,
        local_topology_id_end: 0,
        local_topology_start: String::new(),
        local_topology_end: String::new(),
        local_ti_count_start: 0,
        local_ti_count_end: 0,
    }
}

fn breakend(id: i32, is_start: bool) -> LinxBreakend {
    LinxBreakend {
        id,
        sv_id: 1,
        vcf_id: String::new(),
        coords: String::new(),
        is_start,
        gene: "gene".into(),
        transcript_id: String::new(),
        canonical: true,
        gene_orientation: "Upstream".into(),
        disruptive: false,
        reported_status: "REPORTED".into(),
        undisrupted_copy_number: 0.0,
        region_type: "INTRONIC".into(),
        coding_type: "CODING".into(),
        biotype: String::new(),
        exonic_base_phase: 0,
        next_splice_exon_rank: 0,
        next_splice_exon_phase: 0,
        next_splice_distance: 0,
        total_exon_count: 0,
        exon_up: 0,
        exon_down: 0,
    }
}

#[test]
fn can_interpret_valid_linx_breakends() {
    let mut l = minimal_linx_data();
    l.somatic_breakends = vec![breakend(1, true), breakend(2, false)];
    l.somatic_sv_annotations = vec![sv_annotation()];
    let r = crate::linx::interpret(&l, None, &CytoBands::new(true)).unwrap();
    let (left, right) = (&r.somatic_breakends[0], &r.somatic_breakends[1]);
    assert_eq!((left.chromosome.as_str(), left.chromosome_band.as_str(), left.breakend_type.as_str(), left.orientation), ("1", "p36.33", "DEL", 1));
    assert!((left.junction_copy_number - 1.5).abs() < 1e-6);
    assert_eq!((right.chromosome.as_str(), right.chromosome_band.as_str(), right.orientation), ("1", "p36.33", -1));
    assert!((right.junction_copy_number - 1.5).abs() < 1e-6);
    assert!(r.germline_breakends.is_none());
}

#[test]
fn can_compute_junction_copy_number() {
    assert!((crate::linx::junction_copy_number(&sv_annotation()) - 1.5).abs() < 1e-6);
}

#[test]
fn linx_breakends_take_disruption_drivers_and_reject_duplicate_svs() {
    let mut l = minimal_linx_data();
    let mut unreported = breakend(3, true);
    unreported.reported_status = "NONE".into();
    l.somatic_breakends = vec![breakend(1, true), unreported];
    l.somatic_sv_annotations = vec![sv_annotation()];
    let mut d = driver_catalog("gene", "HOM_DEL_DISRUPTION");
    d.driver_likelihood = 0.75;
    l.somatic_drivers = vec![d];
    let r = crate::linx::interpret(&l, None, &CytoBands::new(true)).unwrap();
    assert_eq!(r.somatic_breakends.len(), 1);
    assert_eq!((r.somatic_breakends[0].driver_type.as_str(), r.somatic_breakends[0].driver_likelihood), ("HOM_DEL_DISRUPTION", 0.75));
    l.somatic_sv_annotations.push(sv_annotation());
    assert!(crate::linx::interpret(&l, None, &CytoBands::new(true)).is_err());
}

// ---------------------------------------------------------------------------------------------
// ImmuneEscapeInterpreterTest (the active test) and the current rules

fn purple_record(drivers: Vec<PurpleDriver>, gains_dels: Vec<PurpleGainDeletion>) -> PurpleRecord {
    let mut r = crate::purple::interpret(&minimal_purple_data(), &minimal_linx_data(), None).unwrap();
    r.somatic_drivers = drivers;
    r.somatic_gains_dels = gains_dels;
    r
}

fn linx_record(breakend_genes: &[(&str, f64)]) -> LinxRecord {
    let mut l = minimal_linx_data();
    for (i, (g, cn)) in breakend_genes.iter().enumerate() {
        let mut b = breakend(i as i32, true);
        b.gene = g.to_string();
        b.undisrupted_copy_number = *cn;
        l.somatic_breakends.push(b);
    }
    l.somatic_sv_annotations = vec![sv_annotation()];
    crate::linx::interpret(&l, None, &CytoBands::new(true)).unwrap()
}

#[test]
fn immune_escape_does_not_crash_on_minimal_data() {
    let e = crate::others::immune_escape(&purple_record(vec![], vec![]), &linx_record(&[]));
    assert!(!e.has_hla_escape && !e.has_antigen_presentation_pathway_escape && !e.has_epigenetic_setdb1_escape);
}

#[test]
fn immune_escape_rules() {
    let none = linx_record(&[]);
    let esc = |p: &PurpleRecord, l: &LinxRecord| crate::others::immune_escape(p, l);
    // HLA: LOH, deletion or a homozygous disruption
    assert!(esc(&purple_record(vec![driver("HLA-A", "LOH")], vec![]), &none).has_hla_escape);
    assert!(esc(&purple_record(vec![driver("HLA-B", "HET_DEL")], vec![]), &none).has_hla_escape);
    assert!(!esc(&purple_record(vec![driver("HLA-A", "AMP")], vec![]), &none).has_hla_escape);
    assert!(esc(&purple_record(vec![], vec![gain_del("HLA-C", "DEL", true)]), &none).has_hla_escape);
    assert!(esc(&purple_record(vec![], vec![]), &linx_record(&[("HLA-A", 0.2)])).has_hla_escape);
    assert!(!esc(&purple_record(vec![], vec![]), &linx_record(&[("HLA-A", 0.5)])).has_hla_escape);
    // the first canonical gain/deletion of the gene decides
    assert!(!esc(&purple_record(vec![], vec![gain_del("B2M", "DEL", false)]), &none).has_antigen_presentation_pathway_escape);
    assert!(esc(&purple_record(vec![], vec![gain_del("B2M", "DEL", true)]), &none).has_antigen_presentation_pathway_escape);
    assert!(!esc(&purple_record(vec![], vec![gain_del("B2M", "HET_DEL", true), gain_del("B2M", "DEL", true)]), &none)
        .has_antigen_presentation_pathway_escape);
    assert!(esc(&purple_record(vec![], vec![gain_del("JAK1", "DEL", true)]), &none).has_ifn_gamma_pathway_escape);
    assert!(esc(&purple_record(vec![], vec![]), &linx_record(&[("CD58", 0.0)])).has_cd58_inactivation_escape);
    // amplifications: AMP only
    assert!(esc(&purple_record(vec![], vec![gain_del("CD274", "AMP", true)]), &none).has_pdl1_overexpression_escape);
    assert!(!esc(&purple_record(vec![], vec![gain_del("CD274", "PARTIAL_AMP", true)]), &none).has_pdl1_overexpression_escape);
    assert!(esc(&purple_record(vec![], vec![gain_del("SETDB1", "AMP", true)]), &none).has_epigenetic_setdb1_escape);
    // OR9: TAPBP is checked only in the corrected mode
    let tapbp = purple_record(vec![], vec![gain_del("TAPBP", "DEL", true)]);
    with_java(false, || assert!(esc(&tapbp, &none).has_antigen_presentation_pathway_escape));
    with_java(true, || assert!(!esc(&tapbp, &none).has_antigen_presentation_pathway_escape));
}

// ---------------------------------------------------------------------------------------------
// Isofox: IsofoxInterpreterTest, IsofoxConversionTest and the selection rules

fn gene_expression(gene: &str, reported: &str, median_cohort: f64, median_cancer: f64) -> rna::GeneExpression {
    rna::GeneExpression {
        gene_name: gene.into(),
        tpm: 10.0,
        spliced_fragments: 0,
        unspliced_fragments: 0,
        median_tpm_cancer: median_cancer,
        percentile_cancer: 0.5,
        median_tpm_cohort: median_cohort,
        percentile_cohort: 0.9,
        reported_status: reported.into(),
    }
}

fn rna_fusion(name: &str, known_type: &str) -> rna::RnaFusion {
    rna::RnaFusion {
        name: name.into(),
        known_type: known_type.into(),
        chromosome_up: "1".into(),
        chromosome_down: "1".into(),
        position_up: 1,
        position_down: 2,
        orientation_up: 1,
        orientation_down: -1,
        junction_type_up: String::new(),
        junction_type_down: String::new(),
        transcript_up: String::new(),
        transcript_down: String::new(),
        exon_up: 0,
        exon_down: 0,
        sv_type: "BND".into(),
        split_fragments: 4,
        realigned_frags: 1,
        discordant_frags: 2,
        depth_up: 10,
        depth_down: 11,
        max_anchor_length_up: 0,
        max_anchor_length_down: 0,
        cohort_frequency: 0,
    }
}

fn minimal_isofox() -> crate::isofox::IsofoxData {
    crate::isofox::IsofoxData {
        summary: rna::RnaStatistics {
            qc_status: vec!["PASS".into()],
            total_fragments: 0,
            duplicate_fragments: 0,
            spliced_fragment_perc: 0.0,
            unspliced_fragment_perc: 0.0,
            alt_fragment_perc: 0.0,
            chimeric_fragment_perc: 0.0,
            spliced_gene_count: 0,
            read_length: 0,
        },
        gene_expressions: Vec::new(),
        fusions: Vec::new(),
        novel_splice_junctions: Vec::new(),
    }
}

#[test]
fn can_interpret_minimal_isofox_data() {
    let r = crate::isofox::interpret(&minimal_isofox(), &[]).unwrap();
    assert!(r.high_expression_genes.is_empty() && r.low_expression_genes.is_empty() && r.fusions.is_empty());
}

#[test]
fn isofox_expression_and_fusion_selection() {
    let mut d = minimal_isofox();
    d.gene_expressions = vec![
        gene_expression("A", "REPORTED", 1.0, 2.0),
        gene_expression("B", "REPORTED", 0.0, 2.0),
        gene_expression("C", "NOT_REPORTED", 1.0, 2.0),
    ];
    d.fusions = vec![rna_fusion("A_B", "KNOWN_PAIR"), rna_fusion("C_D", "KNOWN_PAIR"), rna_fusion("E_F", "NONE"), rna_fusion("G_G", "KNOWN_PAIR"), rna_fusion("H", "KNOWN_PAIR"), rna_fusion("_X", "KNOWN_PAIR")];
    let linx = vec![LinxFusion {
        gene_up: "C".into(),
        cluster_id: 0,
        context_up: String::new(),
        transcript_up: String::new(),
        gene_down: "D".into(),
        context_down: String::new(),
        transcript_down: String::new(),
        reported_type: "KNOWN_PAIR".into(),
        phased: "INFRAME".into(),
        driver_interpretation: "HIGH".into(),
        fused_exon_up: 0,
        fused_exon_down: 0,
        chain_links: 0,
        chain_terminated: false,
        domains_kept: String::new(),
        domains_lost: String::new(),
        junction_copy_number: 0.0,
        rna_support: None,
        plot_filename: None,
    }];
    let r = crate::isofox::interpret(&d, &linx).unwrap();
    assert_eq!(r.high_expression_genes.iter().map(|g| g.gene.as_str()).collect::<Vec<_>>(), vec!["A"]);
    assert_eq!(r.fusions.iter().map(|f| (f.gene_start.as_deref().unwrap(), f.known_type.as_str())).collect::<Vec<_>>(), vec![("A", "KNOWN_PAIR"), ("G", "EXON_DEL_DUP")]);
}

#[test]
fn should_convert_cancer_flag_values_to_null() {
    let good = crate::isofox::convert_gene_expression(&rna::GeneExpression { percentile_cancer: 20.0, ..gene_expression("A", "NONE", 1.0, 10.0) });
    assert_eq!((good.median_tpm_cancer, good.percentile_cancer), (Some(10.0), Some(20.0)));
    let flagged = crate::isofox::convert_gene_expression(&rna::GeneExpression { percentile_cancer: -1.0, ..gene_expression("A", "NONE", 1.0, -1.0) });
    assert_eq!((flagged.median_tpm_cancer, flagged.percentile_cancer), (None, None));
}

// ---------------------------------------------------------------------------------------------
// ViridianInterpreterTest (upstream f5f660cc9f; VirusInterpreterTest and OrangeConversionTest's virus test
// went with VirusInterpreter), and the loading of either tool's file

fn virus_report(name: &str, group: &str, qc: &str, t: Option<&str>, n: i32, cov: f64, depth: f64, copies: Option<f64>, clonality: &str, reported: bool) -> tr::VirusReport {
    tr::VirusReport {
        virus_name: name.into(),
        virus_group: group.into(),
        qc_status: qc.into(),
        reporting_type: t.map(String::from),
        integrations: n,
        coverage_percentage: cov,
        mean_depth: depth,
        copies_per_tumor_cell: copies,
        clonality_status: clonality.into(),
        is_reported: reported,
    }
}

#[test]
fn viridian_interpret() {
    let reported = virus_report("Human papillomavirus type 16", "Human papillomavirus", "AMBIGUOUS_VIRUS", Some("HPV"), 2, 0.9, 20.0, Some(2.5), "CLONAL", true);
    let not_reported = virus_report("JC polyomavirus", "JC polyomavirus", "NO_ABNORMALITIES", None, 0, 0.4, 1.0, None, "PURPLE_FAIL_NO_TUMOR", false);
    let expected_reported = VirusReport {
        name: "Human papillomavirus type 16".into(),
        virus_group: "Human papillomavirus".into(),
        qc_status: "AMBIGUOUS_VIRUS".into(),
        integrations: 2,
        interpretation: Some("HPV".into()),
        percentage_covered: 0.9,
        mean_coverage: 20.0,
        reported: true,
        clonality_status: "CLONAL".into(),
        copies_per_tumor_cell: Some(2.5),
    };
    let expected_not_reported = VirusReport {
        name: "JC polyomavirus".into(),
        virus_group: "JC polyomavirus".into(),
        qc_status: "NO_ABNORMALITIES".into(),
        integrations: 0,
        interpretation: None,
        percentage_covered: 0.4,
        mean_coverage: 1.0,
        reported: false,
        clonality_status: "PURPLE_FAIL_NO_TUMOR".into(),
        copies_per_tumor_cell: None,
    };
    assert_eq!(
        crate::others::interpret_viridian(&[reported, not_reported]),
        ViridianRecord { all_viruses: vec![expected_reported.clone(), expected_not_reported], reportable_viruses: vec![expected_reported] }
    );
    assert_eq!(crate::others::interpret_viridian(&[]), ViridianRecord { all_viruses: vec![], reportable_viruses: vec![] });
}

#[test]
fn virus_types_convert() {
    // ViridianInterpreterTest.convertsEachVirusConstantProperly: VirusType's display names give the
    // VirusInterpretation names; ReportQcStatus and ClonalityStatus have the same names in both
    for (display, name) in [("MCV", "MCV"), ("EBV", "EBV"), ("HPV", "HPV"), ("HBV", "HBV"), ("HHV-8", "HHV8")] {
        assert_eq!(tr::virus_type_from_name(display).unwrap(), name);
    }
    assert_eq!(tr::VIRIDIAN_REPORT_QC_STATUSES, ["NO_ABNORMALITIES", "AMBIGUOUS_VIRUS"]);
    assert_eq!(tr::VIRIDIAN_CLONALITY_STATUSES.len(), 5);
}

#[test]
fn virus_report_json() {
    // as minimally.populated.orange.json writes a VirusReport
    let v = VirusReport {
        name: "JC polyomavirus".into(),
        virus_group: "JC polyomavirus".into(),
        qc_status: "NO_ABNORMALITIES".into(),
        integrations: 0,
        interpretation: None,
        percentage_covered: 0.4,
        mean_coverage: 1.0,
        reported: false,
        clonality_status: "PURPLE_FAIL_NO_TUMOR".into(),
        copies_per_tumor_cell: None,
    };
    let mut w = crate::json::JsonWriter::new();
    crate::json::ToJson::write_json(&ViridianRecord { all_viruses: vec![v], reportable_viruses: vec![] }, &mut w);
    let j = w.finish();
    assert_eq!(
        j,
        "{\n  \"allViruses\": [\n    {\n      \"name\": \"JC polyomavirus\",\n      \"virusGroup\": \"JC polyomavirus\",\n      \"qcStatus\": \"NO_ABNORMALITIES\",\n      \"integrations\": 0,\n      \"interpretation\": null,\n      \"percentageCovered\": 0.4,\n      \"meanCoverage\": 1.0,\n      \"reported\": false,\n      \"clonalityStatus\": \"PURPLE_FAIL_NO_TUMOR\",\n      \"copiesPerTumorCell\": null\n    }\n  ],\n  \"reportableViruses\": []\n}"
    );
}

#[test]
fn virus_files_load() {
    let dir = std::env::temp_dir().join(format!("orange_virus_{}", std::process::id()));
    let (vi, vir, both, none) = (dir.join("vi"), dir.join("viridian"), dir.join("both"), dir.join("none"));
    for d in [&vi, &vir, &both, &none] {
        std::fs::create_dir_all(d).unwrap();
    }
    let vi_text = "taxid\tname\tqcStatus\tintegrations\tinterpretation\tpercentageCovered\tmeanCoverage\texpectedClonalCoverage\treported\tblacklisted\tdriverLikelihood\n\
                   1\tVirus A\tNO_ABNORMALITIES\t3\tHPV\t98.25\t120.5\t44.125\ttrue\tfalse\tHIGH\n\
                   2\tVirus B\tLOW_VIRAL_COVERAGE\t0\tnull\t5.0\t1.25\tnull\tfalse\ttrue\tLOW\n";
    let viridian_text = "name\tqcStatus\tintegrations\tinterpretation\tpercentageCovered\tmeanCoverage\treported\tclonalityStatus\tcopiesPerTumorCell\tvirusGroup\n\
                         Virus C\tNO_ABNORMALITIES\t1\tHHV-8\t0.5\t7.0\ttrue\tCLONAL\t3.0\tGroup C\n";
    std::fs::write(vi.join("T.virus.annotated.tsv"), vi_text).unwrap();
    std::fs::write(vir.join("T.viridian.report.tsv"), viridian_text).unwrap();
    std::fs::write(both.join("T.virus.annotated.tsv"), vi_text).unwrap();
    std::fs::write(both.join("T.viridian.report.tsv"), viridian_text).unwrap();
    let load = |d: &std::path::Path| crate::others::load_virus_dir(&Some(format!("{}/", d.display())), "T");

    // Virus Interpreter's file, converted: blacklisted viruses kept (OR12), coverage as a fraction
    let r = load(&vi).unwrap().unwrap();
    assert_eq!(r.all_viruses.len(), 2);
    assert_eq!(r.reportable_viruses.len(), 1);
    assert_eq!(r.all_viruses[0].percentage_covered, 98.25 * 0.01);
    assert_eq!((r.all_viruses[0].virus_group.as_str(), r.all_viruses[0].clonality_status.as_str()), ("Virus A", "NOT_PROFILED"));
    assert_eq!((r.all_viruses[1].qc_status.as_str(), r.all_viruses[1].copies_per_tumor_cell), ("AMBIGUOUS_VIRUS", None));
    // Viridian's
    let r = load(&vir).unwrap().unwrap();
    assert_eq!((r.all_viruses[0].name.as_str(), r.all_viruses[0].interpretation.as_deref()), ("Virus C", Some("HHV8")));
    assert_eq!(r.reportable_viruses.len(), 1);
    // both files: Virus Interpreter's
    assert_eq!(load(&both).unwrap().unwrap().all_viruses[0].name, "Virus A");
    // a directory without either file fails (Viridian's report is read); no directory, no record
    assert!(load(&none).is_err());
    assert_eq!(crate::others::load_virus_dir(&Some(format!("{}/missing/", dir.display())), "T").unwrap(), None);
    assert_eq!(crate::others::load_virus_dir(&None, "T").unwrap(), None);
    std::fs::remove_dir_all(&dir).unwrap();
}

// ---------------------------------------------------------------------------------------------
// SigsInterpreterTest (with the resource etiologies)

#[test]
fn signatures_are_annotated() {
    let e = tr::snv_signature_etiologies();
    assert_eq!(e.get(&"Sig1".to_string()).map(String::as_str), Some("Deamination of 5-methylcytosine"));
    assert!(e.get(&"MISALLOC".to_string()).is_none());
}

// ---------------------------------------------------------------------------------------------
// CuppaDataFactoryTest

fn cuppa(file: &str) -> Option<CuppaData> {
    let dir = resources()?;
    let entries = tr::read_cuppa_predictions(&format!("{dir}/{file}")).unwrap();
    Some(crate::others::cuppa_data(&entries).unwrap())
}

fn assert_prediction(d: &CuppaData, cancer: &str, expected: Option<[Option<f64>; 6]>) {
    let p = d.predictions.iter().find(|p| p.cancer_type == cancer);
    match expected {
        None => assert!(p.is_none(), "{cancer}"),
        Some(e) => {
            let p = p.unwrap_or_else(|| panic!("{cancer}"));
            let close = |a: Option<f64>, b: Option<f64>| match (a, b) {
                (None, None) => true,
                (Some(a), Some(b)) => (a - b).abs() < 1e-10,
                _ => false,
            };
            assert!(close(Some(p.likelihood), e[0]), "{cancer} likelihood");
            assert!(close(p.genomic_position_classifier, e[1]), "{cancer} gen pos");
            assert!(close(p.snv_pairwise_classifier, e[2]), "{cancer} snv");
            assert!(close(p.feature_classifier, e[3]), "{cancer} feature");
            assert!(close(p.expression_pairwise_classifier, e[4]), "{cancer} expression");
            assert!(close(p.alt_sj_cohort_classifier, e[5]), "{cancer} alt sj");
        }
    }
}

#[test]
fn can_extract_predictions_from_cuppa_including_rna() {
    let Some(d) = cuppa("cuppa/with_rna.vis_data.tsv") else { return };
    assert_eq!(d.predictions.len(), 33);
    assert_eq!(d.best_prediction.cancer_type, "Skin: Melanoma");
    assert_eq!(d.mode, "WGTS");
    assert_eq!(d.max_complex_size, 1000);
    assert_prediction(&d, "Skin: Melanoma", Some([Some(0.9843), Some(0.9643), Some(0.9632), Some(0.9962), Some(0.983), Some(0.9691)]));
    assert_prediction(&d, "Skin: Other", Some([Some(0.0413136), Some(0.0309102), Some(0.0040984), Some(0.0377109), Some(0.0106618), Some(0.0291157)]));
    assert_prediction(&d, "Prostate", Some([Some(0.0325366), Some(0.0229236), Some(0.0008691), Some(0.0169241), Some(0.041203), Some(3.4526e-07)]));
    assert_prediction(&d, "Bone/Soft tissue: Cartilaginous neoplasm", None);
}

#[test]
fn can_extract_predictions_from_cuppa_excluding_rna() {
    let Some(d) = cuppa("test_run/cuppa/tumor_sample.cuppa.vis_data.tsv") else { return };
    assert_eq!(d.predictions.len(), 40);
    assert_eq!(d.best_prediction.cancer_type, "Skin: Melanoma");
    assert_eq!(d.mode, "WGS");
    assert_eq!(d.max_complex_size, 8);
    assert_prediction(&d, "Skin: Melanoma", Some([Some(0.9999710902284965), Some(0.9997621099072876), Some(0.9999910153956342), Some(0.8991388522282147), None, None]));
    assert_prediction(&d, "Skin: Other", Some([Some(0.0), Some(3.159853588544685E-5), Some(5.697607828154083E-6), Some(5.889241560947043E-4), None, None]));
    assert_prediction(&d, "Prostate", Some([Some(0.0), Some(2.8698904215136434E-6), Some(9.922426122823756e-16), Some(3.744164257376693E-5), None, None]));
    assert_prediction(
        &d,
        "Bone/Soft tissue: Cartilaginous neoplasm",
        Some([Some(2.8909771503514605E-5), Some(4.107161442240924e-10), Some(4.414947168645547e-16), Some(7.756820351850527E-10), None, None]),
    );
}

// ---------------------------------------------------------------------------------------------
// PathUtilTest and FileBasedPlotManagerTest

#[test]
fn plot_manager_handles_null_and_missing_files() {
    let dir = std::env::temp_dir().join(format!("orange_plot_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let m = crate::plots::PlotManager::new(&format!("{}/", dir.display()));
    m.create_plot_directory().unwrap();
    assert_eq!(m.process_plot_file(None).unwrap(), None);
    assert_eq!(m.process_plot_file(Some("/no/such/file.png")).unwrap(), None);
    let src = dir.join("source.png");
    std::fs::write(&src, b"png").unwrap();
    with_java(true, || assert_eq!(m.process_plot_file(Some(src.to_str().unwrap())).unwrap().as_deref(), Some("/plot/source.png")));
    with_java(false, || assert_eq!(m.process_plot_file(Some(src.to_str().unwrap())).unwrap().as_deref(), Some("plot/source.png")));
    assert!(dir.join("plot/source.png").exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

// ---------------------------------------------------------------------------------------------
// the JSON as Gson writes it for a whole record

#[test]
fn record_json_layout() {
    let purple = crate::purple::interpret(&minimal_purple_data(), &minimal_linx_data(), None).unwrap();
    let r = OrangeRecord {
        sample_id: "T".into(),
        reference_id: None,
        rna_id: None,
        sampling_date: "2026-01-02".into(),
        experiment_type: "WHOLE_GENOME".into(),
        configured_primary_tumor: vec![OrangeDoidNode { doid: String::new(), doid_term: Some("Skin & 'bone'".into()) }],
        ref_genome_version: "V38".into(),
        pipeline_version: Some("OA_V3_1".into()),
        panel_name: None,
        sequencing_type: "ILLUMINA".into(),
        purple,
        linx: LinxRecord { fusions: Vec::new(), somatic_breakends: Vec::new(), germline_breakends: None },
        isofox: None,
        lilac: None,
        immune_escape: crate::others::immune_escape(&purple_record(vec![], vec![]), &linx_record(&[])),
        viridian: None,
        chord: None,
        cuppa: None,
        peach: Some(Vec::new()),
        sig_allocations: None,
        plots: OrangePlots::default(),
    };
    let j = to_json(&r);
    assert!(j.starts_with("{\n  \"sampleId\": \"T\",\n  \"referenceId\": null,\n  \"rnaId\": null,\n  \"samplingDate\": \"2026-01-02\""));
    assert!(j.contains("\"doidTerm\": \"Skin \\u0026 \\u0027bone\\u0027\""));
    assert!(j.contains("\"germlineBreakends\": null"));
    assert!(j.contains("\"peach\": [],\n  \"sigAllocations\": null,"));
    assert!(j.contains("\"purity\": 0.5,"));
    assert!(j.ends_with("    \"qSeePlot\": null\n  }\n}"));
}
