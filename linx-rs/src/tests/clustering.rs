//! clustering: SimpleMergeRuleTest and ComplexMergeRuleTest

use hmf_common::sv::StructuralVariantType::*;

use super::*;
use crate::cn::{HomLoss, HomLossEvent};
use ClusteringReason::*;

fn hom_loss(chr: &str, s: i32, e: i32, ts: &str, te: &str, sv1: i32, sv2: i32) -> HomLoss {
    HomLoss::new(HomLossEvent::new(chr, s, e, ts, te, sv1, sv2))
}
fn loh(chr: &str, s: i32, e: i32, ts: &str, te: &str, sv1: i32, sv2: i32) -> Loh {
    Loh::new(LohEvent::new(chr, s, e, ts, te, 1, sv1, sv2))
}

#[test]
fn test_proximity_merge() {
    let mut t = LinxTester::new();
    let var1 = create_del(t.next_var_id(), "1", 1000, 1100);
    let var2 = create_dup(t.next_var_id(), "1", 2000, 2100);
    let var3 = create_ins(t.next_var_id(), "1", 3000, 3100);
    let var4 = create_inv(t.next_var_id(), "1", 4000, 4100, 1);
    let var5 = create_sgl(t.next_var_id(), "1", 5000, -1);
    let var6 = create_bnd(t.next_var_id(), "1", 6000, -1, "5", 1000, 1);
    let var10 = create_ins(t.next_var_id(), "5", 2000, 20000);
    let var11 = create_bnd(t.next_var_id(), "2", 6000, -1, "5", 21000, 1);
    let var12 = create_ins(t.next_var_id(), "2", 7000, 8000);
    let var13 = create_bnd(t.next_var_id(), "3", 6000, -1, "5", 22000, 1);
    t.add(&[&var1, &var2, &var3, &var4, &var5, &var6, &var10, &var11, &var12, &var13]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(1, t.clusters().len());
    assert!(t.has_cluster_with_svs(&[&var1, &var2, &var3, &var4, &var5, &var6, &var10, &var11, &var12, &var13]));
}

#[test]
fn test_simple_svs_demerge() {
    let mut t = LinxTester::new();
    t.clear_clusters_and_svs();
    let del1 = create_del(t.next_var_id(), "2", 1000, 3000);
    let del2 = create_del(t.next_var_id(), "2", 2000, 4000);
    let del3 = create_del(t.next_var_id(), "2", 3500, 8000);
    let var1 = create_del(t.next_var_id(), "1", 20000, 25000);
    let var2 = create_del(t.next_var_id(), "1", 26500, 50000);
    t.add(&[&del1, &del2, &del3, &var1, &var2]);
    set_assembled_link_info(&var1, false, &var2, true);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(3, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&del1, &del2, &del3]).unwrap();
    assert_eq!(ResolvedType::COMPLEX, c.resolved_type());
}

#[test]
fn test_long_del_dup_inv_merge() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 200000, -1);
    let var2 = create_inv(t.next_var_id(), "1", 20000, 220000, -1);
    t.add(&[&var1, &var2]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(2, t.clusters().len());

    t.clear_clusters_and_svs();
    let var3 = create_inv(t.next_var_id(), "1", 40000, 240000, 1);
    t.add(&[&var1, &var2, &var3]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(1, t.clusters().len());
    let c = t.clusters()[0].clone();
    assert!(var1.has_cluster_reason(LONG_DEL_DUP_INV));
    assert!(c.has_cluster_reason(LONG_DEL_DUP_INV));

    t.clear_clusters_and_svs();
    let var1 = create_del(t.next_var_id(), "1", 10000, 30000000);
    let var2 = create_inv(t.next_var_id(), "1", 500000, 10000000, 1);
    let var3 = create_dup(t.next_var_id(), "1", 20000000, 40000000);
    t.add(&[&var1, &var2, &var3]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(3, t.clusters().len());

    t.clear_clusters_and_svs();
    let var1 = create_dup(t.next_var_id(), "1", 10000000, 30001000);
    let var2 = create_dup(t.next_var_id(), "1", 30000000, 50001000);
    let var3 = create_dup(t.next_var_id(), "1", 50000000, 70001000);
    let var4 = create_inv(t.next_var_id(), "1", 70000000, 71000000, -1);
    let var5 = create_del(t.next_var_id(), "1", 1000000, 20000000);
    let var6 = create_del(t.next_var_id(), "1", 20001000, 40000000);
    let var7 = create_del(t.next_var_id(), "1", 40001000, 60000000);
    t.add(&[&var1, &var2, &var3, &var4, &var5, &var6, &var7]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(1, t.clusters().len());
    let c = t.clusters()[0].clone();
    assert!(t.all_variants.iter().any(|x| x.has_cluster_reason(LONG_DEL_DUP_INV)));
    assert!(c.has_cluster_reason(LONG_DEL_DUP_INV));
}

#[test]
fn test_high_facing_jcn_merge() {
    let mut t = LinxTester::new();
    let var1 = create_del(1, "1", 1000, 2000);
    let var2 = create_test_sv(2, "1", "1", 10000, 40000, -1, 1, DUP, 6.0);
    let var3 = create_test_sv(3, "1", "1", 30000, 70000, -1, 1, DUP, 6.0);
    let var4 = create_test_sv(4, "1", "1", 72000, 80000, -1, -1, INV, 2.0);
    t.add(&[&var1, &var2, &var3, &var4]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(2, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var2, &var3, &var4]).unwrap();
    assert!(c.has_cluster_reason(HIGH_JCN));

    t.clear_clusters_and_svs();
    let var1 = create_del(1, "1", 1000, 2000);
    let var2 = create_test_sv(2, "1", "", 10000, 0, -1, 0, SGL, 40.0);
    let var3 = create_test_sv(3, "1", "", 30000, 0, 1, 0, SGL, 6.0);
    let var4 = create_del(4, "1", 40000, 41000);
    let var5 = create_dup(5, "1", 50000, 51000);
    let var6 = create_test_sv(6, "1", "", 60000, 0, -1, 0, SGL, 10.0);
    let var7 = create_test_sv(7, "1", "", 70000, 0, 1, 0, SGL, 36.0);
    let var8 = create_test_sv(8, "1", "", 80000, 0, 1, 0, SGL, 8.0);
    t.add(&[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(4, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var2, &var3, &var6, &var7, &var8]).unwrap();
    assert!(c.has_cluster_reason(HIGH_JCN));
    assert!(c.has_cluster_reason(MAJOR_ALLELE_JCN));

    t.clear_clusters_and_svs();
    let var1 = create_del(1, "1", 1000, 2000);
    let var2 = create_test_sv(2, "1", "", 10000, 0, -1, 0, SGL, 40.0);
    let var3 = create_test_sv(3, "1", "", 30000, 0, 1, 0, SGL, 30.0);
    let var4 = create_test_sv(4, "1", "", 40000, 0, 1, 0, SGL, 8.0);
    let var5 = create_dup(5, "1", 50000, 51000);
    let var6 = create_del(6, "1", 60000, 61000);
    let var7 = create_test_sv(7, "1", "1", 70000, 100000, -1, 1, DUP, 10.0);
    let var8 = create_test_sv(8, "1", "1", 90000, 120000, -1, 1, DUP, 10.0);
    t.add(&[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(5, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var2, &var3, &var4]).unwrap();
    assert!(c.has_cluster_reason(HIGH_JCN));
    let c = t.find_cluster_with_svs(&[&var7, &var8]).unwrap();
    assert!(c.has_cluster_reason(HIGH_JCN));
}

#[test]
fn test_loh_hom_loss_event_merge() {
    let mut t = LinxTester::new();
    // scenario 1: an LOH containing all clustered hom-loss events is clustered too
    let var1 = create_bnd(1, "1", 1000, 1, "2", 100, 1);
    let var2 = create_bnd(2, "1", 100000, -1, "3", 100, 1);
    let var3 = create_del(3, "1", 6500, 6600);
    let var4 = create_bnd(4, "1", 20000, 1, "5", 200, 1);
    let var5 = create_bnd(5, "1", 22000, -1, "5", 100, -1);
    let l = loh("1", 1000, 100000, "BND", "BND", var1.id(), var2.id());
    t.cn_loader.loh_event_data.push(l.clone());
    t.cn_loader.hom_loss_data.push(hom_loss(var3.chromosome(true), var3.position(true), var3.position(false), var3.type_str(), var3.type_str(), var3.id(), var3.id()));
    t.cn_loader.hom_loss_data.push(hom_loss(var4.chromosome(true), var4.position(true), var5.position(true), var4.type_str(), var5.type_str(), var4.id(), var5.id()));
    l.add_hom_loss_events(&t.cn_loader.hom_loss_data);
    t.add(&[&var1, &var2, &var3, &var4, &var5]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(3, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var1, &var2]).unwrap();
    assert!(c.has_cluster_reason(HOM_LOSS));

    // scenario 2: several hom-loss events clustered as their LOH is
    t.clear_clusters_and_svs();
    let var1 = create_del(1, "1", 1000, 100000);
    let var2 = create_bnd(2, "1", 10000, 1, "2", 100, 1);
    let var3 = create_bnd(3, "1", 20000, -1, "3", 100, 1);
    let var4 = create_bnd(4, "1", 30000, 1, "4", 100, 1);
    let var5 = create_bnd(5, "1", 40000, -1, "5", 100, 1);
    t.cn_loader.loh_event_data.clear();
    let l = loh(var1.chromosome(true), var1.position(true), var1.position(false), "DEL", "DEL", var1.id(), var1.id());
    t.cn_loader.loh_event_data.push(l.clone());
    t.cn_loader.hom_loss_data.clear();
    t.cn_loader.hom_loss_data.push(hom_loss(var2.chromosome(true), var2.position(true), var3.position(true), var2.type_str(), var3.type_str(), var2.id(), var3.id()));
    t.cn_loader.hom_loss_data.push(hom_loss(var4.chromosome(true), var4.position(true), var5.position(true), var4.type_str(), var5.type_str(), var4.id(), var5.id()));
    l.add_hom_loss_events(&t.cn_loader.hom_loss_data);
    t.add(&[&var1, &var2, &var3, &var4, &var5]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(2, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var2, &var3, &var4, &var5]).unwrap();
    assert!(c.has_cluster_reason(HOM_LOSS));
    assert!(c.has_cluster_reason(MAJOR_ALLELE_JCN));

    // scenario 3: a hom-loss event overlapping an LOH
    t.clear_clusters_and_svs();
    let var1 = create_del(1, "1", 10000, 30000);
    let var2 = create_bnd(2, "1", 20000, 1, "2", 100, 1);
    let var3 = create_bnd(3, "1", 40000, -1, "3", 100, 1);
    t.cn_loader.loh_event_data.clear();
    let l = loh(var1.chromosome(true), var1.position(true), var3.position(true), var1.type_str(), var3.type_str(), var1.id(), var3.id());
    t.cn_loader.loh_event_data.push(l.clone());
    t.cn_loader.hom_loss_data.clear();
    t.cn_loader.hom_loss_data.push(hom_loss(var2.chromosome(true), var2.position(true), var1.position(false), var2.type_str(), var1.type_str(), var2.id(), var1.id()));
    l.add_hom_loss_events(&t.cn_loader.hom_loss_data);
    t.add(&[&var1, &var2, &var3]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let c = t.find_cluster_with_svs(&[&var2, &var3]).unwrap();
    assert!(c.has_cluster_reason(HOM_LOSS));

    // again with more SVs
    t.clear_clusters_and_svs();
    let var1 = create_del(1, "1", 10000, 30000);
    let var2 = create_del(2, "1", 40000, 60000);
    let var3 = create_bnd(3, "1", 20000, 1, "2", 100, 1);
    let var4 = create_bnd(4, "1", 50000, -1, "3", 100, 1);
    t.cn_loader.loh_event_data.clear();
    let l = loh(var1.chromosome(true), var1.position(true), var2.position(false), var1.type_str(), var2.type_str(), var1.id(), var2.id());
    t.cn_loader.loh_event_data.push(l.clone());
    t.cn_loader.hom_loss_data.clear();
    t.cn_loader.hom_loss_data.push(hom_loss(var2.chromosome(true), var3.position(true), var1.position(false), var3.type_str(), var1.type_str(), var3.id(), var1.id()));
    t.cn_loader.hom_loss_data.push(hom_loss(var2.chromosome(true), var2.position(true), var4.position(true), var2.type_str(), var4.type_str(), var2.id(), var4.id()));
    l.add_hom_loss_events(&t.cn_loader.hom_loss_data);
    t.add(&[&var1, &var2, &var3, &var4]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let c = t.find_cluster_with_svs(&[&var3, &var4]).unwrap();
    assert!(c.has_cluster_reason(HOM_LOSS));
}

// ComplexMergeRuleTest

#[test]
fn test_foldback_merge() {
    let mut t = LinxTester::new();
    let inv1 = create_inv(t.next_var_id(), "1", 100, 200, -1);
    let inv2 = create_inv(t.next_var_id(), "1", 20000, 20100, 1);
    t.add(&[&inv1, &inv2]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(1, t.clusters().len());
    assert!(inv1.has_cluster_reason(FOLDBACKS));
    assert!(inv2.has_cluster_reason(FOLDBACKS));
}

#[test]
fn test_loh_resolving_cluster_merge() {
    let mut t = LinxTester::new();
    let var1 = create_del(1, "1", 10000, 60000);
    let var2 = create_dup(2, "1", 20000, 50000);
    let var3 = create_bnd(3, "1", 30000, 1, "2", 10000, 1);
    let var4 = create_bnd(4, "1", 40000, -1, "2", 20000, 1);
    let var5 = create_del(5, "1", 70000, 120000);
    let var6 = create_dup(6, "1", 80000, 110000);
    let var7 = create_bnd(7, "1", 90000, 1, "2", 80000, 1);
    let var8 = create_bnd(8, "1", 100000, -1, "2", 90000, -1);
    let var9 = create_dup(9, "2", 30000, 60000);
    let var10 = create_del(10, "2", 40000, 50000);
    t.add(&[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8, &var9, &var10]);
    t.cn_loader.loh_event_data.push(loh("1", 10000, 20000, "DEL", "DUP", var1.id(), var2.id()));
    t.cn_loader.loh_event_data.push(loh("1", 50000, 60000, "DUP", "DEL", var2.id(), var1.id()));
    t.cn_loader.loh_event_data.push(loh("1", 110000, 120000, "DUP", "DEL", var6.id(), var5.id()));
    t.cn_loader.loh_event_data.push(loh("2", 20000, 30000, "BND", "DUP", var4.id(), var9.id()));
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(2, t.clusters().len());
    assert_eq!(1, var10.cl().sv_count());
    assert!(var2.has_cluster_reason(LOH_CHAIN));
    assert!(var2.cluster_reason().contains(&var4.id().to_string()));
    assert!(var4.has_cluster_reason(LOH_CHAIN));
    assert!(var4.cluster_reason().contains(&var2.id().to_string()));
}

#[test]
fn test_consistent_breakend_overlap_merge() {
    let mut t = LinxTester::new();
    let consec1 = create_bnd(t.next_var_id(), "1", 100000, 1, "2", 100, -1);
    let consec2 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 200, 1);
    let consec3 = create_inv(t.next_var_id(), "1", 30000, 101000, 1);
    let var1 = create_bnd(t.next_var_id(), "1", 10000, 1, "3", 200, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 10100, -1, "3", 100, -1);
    let var3 = create_del(t.next_var_id(), "1", 60000, 60100);
    let overlap1 = create_bnd(t.next_var_id(), "1", 20000, -1, "4", 200, 1);
    let overlap2 = create_sgl(t.next_var_id(), "1", 40000, -1);
    t.add(&[&consec1, &consec2, &consec3, &var1, &var2, &var3, &overlap1, &overlap2]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(3, t.all_clusters().len());
    assert!(t.find_cluster_with_svs(&[&consec1, &consec2, &consec3, &overlap1, &overlap2]).is_some());
    assert!(overlap1.has_cluster_reason(CONSEC_BREAKS));
    assert!(overlap2.has_cluster_reason(CONSEC_BREAKS));
}

#[test]
fn test_foldbacks_straddling_breakend_merge() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 2000, -1);
    let var2 = create_inv(t.next_var_id(), "1", 100000, 101000, -1);
    let var3 = create_del(t.next_var_id(), "1", 20000, 21000);
    let var4 = create_bnd(t.next_var_id(), "1", 40000, -1, "3", 100, -1);
    t.add(&[&var1, &var2, &var3, &var4]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert!(var1.is_foldback());
    assert!(var2.is_foldback());
    assert_eq!(2, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var1, &var2, &var4]).unwrap();
    assert!(var4.has_cluster_reason(OVERLAP_FOLDBACKS));
    assert!(c.has_cluster_reason(OVERLAP_FOLDBACKS));

    t.clear_clusters_and_svs();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 2000, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 40000, -1, "3", 100, -1);
    t.add(&[&var1, &var2, &var3]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert!(var1.is_foldback());
    assert!(var2.is_foldback());
    assert_eq!(2, t.clusters().len());
}

#[test]
fn test_chained_foldback_merge() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 100, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 2000, -1, "2", 200, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 6101000, 1, "4", 100, -1);
    let var4 = create_bnd(t.next_var_id(), "1", 6102000, 1, "4", 200, 1);
    let var5 = create_dup(t.next_var_id(), "1", 6105000, 6110000);
    let var6 = create_dup(t.next_var_id(), "1", 6108000, 6112000);
    let var7 = create_del(t.next_var_id(), "1", 6111000, 6120000);
    let var8 = create_sgl(t.next_var_id(), "1", 6115000, -1);
    t.add(&[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert!(var1.is_chained_foldback());
    assert!(var2.is_chained_foldback());
    assert!(var3.is_chained_foldback());
    assert!(var4.is_chained_foldback());
    assert_eq!(1, t.clusters().len());
    assert!(var1.has_cluster_reason(FOLDBACKS));
}

#[test]
fn test_distance_ploidy_link_match_merge() {
    let mut t = LinxTester::new();
    t.set_non_cluster_allele_ploidies(4.0, 1.0);
    let var0 = create_test_sv(t.next_var_id(), "1", "0", 500, -1, -1, -1, SGL, 2.0);
    let var1 = create_test_sv(t.next_var_id(), "1", "2", 1000, 100, 1, -1, BND, 2.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 2000, 2000, -1, 1, BND, 2.0);
    let var3 = create_del(t.next_var_id(), "1", 20000, 21000);
    let var4 = create_test_sv(t.next_var_id(), "1", "3", 50000, 100, 1, -1, BND, 1.0);
    let var5 = create_test_sv(t.next_var_id(), "1", "3", 51000, 200, -1, 1, BND, 1.0);
    let var6 = create_test_sv(t.next_var_id(), "1", "1", 100000, 110000, 1, 1, INV, 2.0);
    t.add(&[&var0, &var1, &var2, &var3, &var4, &var5, &var6]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(3, t.clusters().len());
    assert!(t.has_cluster_with_svs(&[&var4, &var5]));
    let c = t.find_cluster_with_svs(&[&var0, &var1, &var2, &var6]).unwrap();
    assert!(var6.has_cluster_reason(TI_JCN_MATCH));
    assert!(c.has_cluster_reason(TI_JCN_MATCH));
}

#[test]
fn test_satellite_sgl_merge() {
    let mut t = LinxTester::new();
    let var0 = create_test_sv_full(t.next_var_id(), "1", "0", 500, -1, 1, 0, SGL, 2.0, 0.0, 1.0, 0.0, 1.0, "", "PASS", SGL_CENTRO_SATELLITE, "SAR");
    let var1 = create_test_sv_full(t.next_var_id(), "1", "0", 50000, -1, -1, 0, SGL, 2.0, 0.0, 1.0, 0.0, 1.0, "", "PASS", "", "SAR");
    t.add(&[&var0, &var1]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(1, t.clusters().len());
    assert!(var0.has_cluster_reason(SATELLITE_SGL));
    assert!(var1.has_cluster_reason(SATELLITE_SGL));
}
