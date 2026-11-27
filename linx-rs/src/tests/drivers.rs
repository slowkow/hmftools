//! drivers: AmpDriversTest and HomDisruptionsTest

use hmf_common::driver_panel::{DriverCategory, DriverGene, GermlineReporting};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::purple_files::{DriverCatalog, GeneCopyNumber};
use hmf_common::sv::StructuralVariantType::*;

use super::fusion::{set_sv_gene_data_test, FusionEnv};
use super::*;
use crate::drivers::{DriverEventType, DriverGeneAnnotator};

/// GeneTestUtils.createDriverGene: a TSG reported for disruptions
pub fn create_driver_gene(name: &str) -> DriverGene {
    DriverGene {
        gene: name.into(),
        report_missense_and_inframe: false,
        report_nonsense_and_frameshift: false,
        report_splice: false,
        report_deletion: false,
        report_het_deletion: false,
        report_loh: false,
        het_deletion_threshold: hmf_common::driver_panel::DEFAULT_DRIVER_HET_DELETION_THRESHOLD,
        report_disruption: true,
        report_amplification: false,
        amplification_ratio: hmf_common::driver_panel::DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO,
        report_somatic_hotspot: false,
        likelihood_type: DriverCategory::Tsg,
        report_germline_variant: GermlineReporting::None,
        report_germline_hotspot: GermlineReporting::None,
        report_germline_disruption: GermlineReporting::Any,
        report_germline_deletion: GermlineReporting::Any,
        report_germline_amplification: false,
        additional_reported_transcripts: Vec::new(),
        report_pgx: false,
        report_high_expression: false,
        report_low_expression: false,
        report_novel_splice_junction: false,
    }
}

/// SvTestUtils.createGeneCopyNumber
fn create_gene_copy_number(gene: &str, chromosome: &str, min_cn: f64, start: i32, end: i32) -> GeneCopyNumber {
    GeneCopyNumber {
        chromosome: chromosome.into(),
        start,
        end,
        gene_name: gene.into(),
        trans_name: String::new(),
        is_canonical: true,
        chromosome_band: String::new(),
        max_copy_number: min_cn,
        min_copy_number: min_cn,
        min_minor_allele_copy_number: 0.0,
        somatic_regions: 1,
        min_regions: 1,
        min_region_start: 0,
        min_region_end: 0,
        depth_window_count: 0,
        gc_content: 1.0,
        min_region_start_support: "BND".into(),
        min_region_end_support: "BND".into(),
        min_region_method: "BAF_WEIGHTED".into(),
        relative_min_copy_number: min_cn / 2.0,
        reported_status: String::new(),
        driver_type: String::new(),
    }
}

/// SvTestUtils.createDriver
fn create_driver(gene: &str, chromosome: &str, driver: &str, category: &str, biallelic: bool, min_cn: f64) -> DriverCatalog {
    let method = match driver {
        "AMP" | "PARTIAL_AMP" => "AMP",
        "DEL" => "DEL",
        _ => "BIALLELIC",
    };
    DriverCatalog {
        chromosome: chromosome.into(),
        chromosome_band: String::new(),
        gene: gene.into(),
        transcript: String::new(),
        is_canonical: true,
        driver: driver.into(),
        category: category.into(),
        likelihood_method: method.into(),
        reported_status: "REPORTED".into(),
        driver_likelihood: 1.0,
        missense: 0,
        nonsense: 0,
        splice: 0,
        inframe: 0,
        frameshift: 0,
        biallelic,
        min_copy_number: min_cn,
        max_copy_number: min_cn,
    }
}

fn annotator(gc: &'static GeneCache, driver_genes: JavaHashMap<String, DriverGene>) -> DriverGeneAnnotator {
    let dg: &'static JavaHashMap<String, DriverGene> = Box::leak(Box::new(driver_genes));
    let mut a = DriverGeneAnnotator::new(gc, dg, "");
    a.set_sample_purity_data(2.0, false);
    a
}
fn cluster_and_annotate(t: &mut LinxTester, a: &mut DriverGeneAnnotator, vars: &[&Sv]) {
    t.clear_clusters_and_svs();
    t.add(vars);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    a.clear_results(false);
    let map = t.analyser.state.chr_breakend_map.clone();
    a.annotate_svs(&t.sample_id.clone(), &map, &t.cn_loader, &mut t.vis);
}

#[test]
fn test_amp_drivers() {
    let mut t = LinxTester::new();
    t.set_non_cluster_allele_ploidies(1.0, 0.0);
    let mut cache = create_gene_data_cache();
    let (chr, ts, te) = (CHR_1, 100000, 120000);
    let td = TranscriptData {
        trans_id: 1,
        trans_name: "TRAN0001".into(),
        gene_id: GENE_ID_1.into(),
        is_canonical: true,
        strand: 1,
        trans_start: ts,
        trans_end: te,
        coding_start: None,
        coding_end: None,
        bio_type: String::new(),
        ref_seq_id: None,
        exons: Vec::new(),
    };
    add_trans_exon_data(&mut cache, GENE_ID_1, vec![td]);
    add_gene_data(&mut cache, chr, vec![create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, chr, 1, ts, te)]);
    let gc = leak_cache(cache);
    let mut a = annotator(gc, JavaHashMap::new());
    a.add_driver_gene(create_driver(GENE_NAME_1, chr, "AMP", "ONCO", false, 12.0), create_gene_copy_number(GENE_NAME_1, chr, 12.0, ts, te));

    let var1 = create_del(t.next_var_id(), chr, 500, 600);
    let var2 = create_test_sv(t.next_var_id(), chr, chr, ts - 1000, te + 1000, -1, 1, DUP, 10.0);
    cluster_and_annotate(&mut t, &mut a, &[&var1, &var2]);
    let c = t.find_cluster_with_svs(&[&var2]).unwrap();
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(1, d.events.len());
    assert!(d.events[0].cluster() == Some(c));

    let var1 = create_del(t.next_var_id(), chr, 500, 600);
    let var2 = create_test_sv(t.next_var_id(), chr, "2", 10000, 1000, -1, 1, BND, 5.0);
    let var3 = create_del(t.next_var_id(), chr, 20000, 21000);
    let var4 = create_dup(t.next_var_id(), chr, 30000, 31000);
    let var5 = create_test_sv(t.next_var_id(), chr, chr, 40000, 45000, -1, -1, INV, 2.5);
    cluster_and_annotate(&mut t, &mut a, &[&var1, &var2, &var3, &var4, &var5]);
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(2, d.events.len());
    let clusters: Vec<Option<Cluster>> = d.events.iter().map(|e| e.cluster()).collect();
    assert!(clusters.contains(&t.find_cluster_with_svs(&[&var2])));
    assert!(clusters.contains(&t.find_cluster_with_svs(&[&var5])));

    let var1 = create_del(t.next_var_id(), chr, 500, 600);
    let var2 = create_test_sv(t.next_var_id(), chr, "2", 6300000, 1000, 1, 1, BND, 5.0);
    let var3 = create_test_sv(t.next_var_id(), chr, "3", 6250000, 1000, 1, 1, BND, 15.0);
    let var4 = create_test_sv(t.next_var_id(), chr, chr, 200000, 201000, -1, -1, INV, 5.0);
    cluster_and_annotate(&mut t, &mut a, &[&var1, &var2, &var3, &var4]);
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(1, d.events.len());
    assert!(d.events[0].cluster() == t.find_cluster_with_svs(&[&var3]));
}

fn tp53_env(exon_starts: &[i32], ts: i32, te: i32, cs: i32, ce: i32) -> (FusionEnv, &'static TranscriptData) {
    let mut env = FusionEnv::new();
    let (name, id, chr) = ("TP53", "ENSG00000141510", "17");
    let td = create_trans_exons(id, 1, 1, exon_starts, 100, Some(cs), Some(ce), true, "");
    add_trans_exon_data(&mut env.ensembl, id, vec![td.clone()]);
    add_gene_data(&mut env.ensembl, chr, vec![create_ensembl_gene_data(id, name, chr, 1, ts, te)]);
    (env, Box::leak(Box::new(td)))
}

fn driver_genes() -> JavaHashMap<String, DriverGene> {
    let mut m = JavaHashMap::new();
    m.insert("TP53".to_string(), create_driver_gene("TP53"));
    m
}

fn annotate_hom(
    t: &mut LinxTester, gc: &'static GeneCache, fa: &crate::fusion::analyser::FusionDisruptionAnalyser, a: &mut DriverGeneAnnotator, vars: &[&Sv],
    clear_drivers: bool,
) {
    t.clear_clusters_and_svs();
    t.add(vars);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    set_sv_gene_data_test(&t.all_variants, gc, 0);
    fa.annotate_transcripts(&t.all_variants, false);
    a.clear_results(clear_drivers);
    let map = t.analyser.state.chr_breakend_map.clone();
    a.annotate_svs(&t.sample_id.clone(), &map, &t.cn_loader, &mut t.vis);
}

#[test]
fn test_simple_del_disruption() {
    let mut t = LinxTester::new();
    t.set_non_cluster_allele_ploidies(1.0, 0.0);
    let (env, td) = tp53_env(&[10000, 40000, 70000, 100000, 130000], 10000, 130000, 10050, 130050);
    let (gc, fa) = env.analyser();
    let mut a = annotator(gc, driver_genes());
    let chr = "17";
    let del1 = create_del(t.next_var_id(), chr, 12000, 24000);
    let del2 = create_del(t.next_var_id(), chr, 18000, 30000);
    let del3 = create_del(t.next_var_id(), chr, 42000, 85000);
    let del4 = create_del(t.next_var_id(), chr, 75000, 120000);
    annotate_hom(&mut t, gc, &fa, &mut a, &[&del1, &del2, &del3, &del4], false);
    assert_eq!(4, t.clusters().len());
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(1, d.events.len());
    assert_eq!(td, d.trans);
    let c = d.events[0].cluster().unwrap();
    assert_eq!(ResolvedType::DEL, c.resolved_type());
    assert!(contains(&c.svs(), &del4));
    assert_eq!(DriverEventType::HOM_DEL_DISRUPTION, d.events[0].event_type);
}

#[test]
fn test_chained_del_disruption() {
    let mut t = LinxTester::new();
    t.set_non_cluster_allele_ploidies(0.2, 0.0);
    let (env, _) = tp53_env(&[10000, 20000, 30000, 40000, 50000], 10000, 55000, 10002, 50098);
    let (gc, fa) = env.analyser();
    let mut a = annotator(gc, driver_genes());
    let chr = "17";
    let var1 = create_inv(t.next_var_id(), chr, 15000, 25000, 1);
    let var2 = create_dup(t.next_var_id(), chr, 15100, 35000);
    let var3 = create_inv(t.next_var_id(), chr, 25100, 35100, -1);
    annotate_hom(&mut t, gc, &fa, &mut a, &[&var1, &var2, &var3], false);
    assert_eq!(1, t.clusters().len());
    let c = t.clusters()[0].clone();
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(1, d.events.len());
    assert!(d.events[0].cluster() == Some(c));
    assert_eq!(DriverEventType::HOM_DEL_DISRUPTION, d.events[0].event_type);

    let del = create_del(t.next_var_id(), chr, 400, 500);
    let var1 = create_inv(t.next_var_id(), chr, 15000, 25000, 1);
    let var2 = create_inv(t.next_var_id(), chr, 14990, 24990, -1);
    annotate_hom(&mut t, gc, &fa, &mut a, &[&del, &var1, &var2], true);
    assert_eq!(2, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var1, &var2]);
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(1, d.events.len());
    assert!(d.events[0].cluster() == c);
    assert_eq!(DriverEventType::HOM_DEL_DISRUPTION, d.events[0].event_type);
}

#[test]
fn test_dup_disruption() {
    let mut t = LinxTester::new();
    t.set_non_cluster_allele_ploidies(0.2, 0.0);
    let (env, _) = tp53_env(&[10000, 20000, 30000, 40000, 50000], 10000, 55000, 10002, 50098);
    let (gc, fa) = env.analyser();
    let mut a = annotator(gc, driver_genes());
    let var1 = create_dup(t.next_var_id(), "17", 15000, 35000);
    annotate_hom(&mut t, gc, &fa, &mut a, &[&var1], false);
    assert_eq!(1, t.clusters().len());
    let c = t.clusters()[0].clone();
    assert_eq!(1, a.driver_gene_data.len());
    let d = &a.driver_gene_data[0];
    assert_eq!(1, d.events.len());
    assert!(d.events[0].cluster() == Some(c));
    assert_eq!(DriverEventType::HOM_DUP_DISRUPTION, d.events[0].event_type);
}
