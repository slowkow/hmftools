//! analyser: SyntheticTest, InversionPairTest, TranslocationPairTest, FoldbackTest, FilteringTest,
//! CopyNumberTest, AnnotationTest, ReplicationRepairTest, LineTest and DoubleMinuteTest

use hmf_common::sv::StructuralVariantType::*;

use super::*;
use crate::analysis::classification::{get_synthetic_gap_length, get_synthetic_length, get_synthetic_ti_length};
use crate::analysis::clustering::mark_foldbacks;
use crate::cn::{calc_adjusted_jcn_values, CN_DATA_NO_SV, TELOMERE};
use ArmClusterType as AC;
use ResolvedType as RT;

fn only_cluster(t: &LinxTester) -> Cluster {
    assert_eq!(1, t.clusters().len());
    t.clusters()[0].clone()
}
fn arm_cluster_types(c: &Cluster) -> Vec<AC> {
    c.arm_clusters.borrow().iter().map(|a| a.arm_cluster_type).collect()
}
fn arm_cluster_data(c: &Cluster, t: AC) -> usize {
    arm_cluster_types(c).iter().filter(|x| **x == t).count()
}
fn arm_cluster_ti_count(c: &Cluster, i: usize) -> i32 {
    c.arm_clusters.borrow()[i].ti_count
}
fn run(t: &mut LinxTester, vars: &[&Sv]) {
    t.clear_clusters_and_svs();
    t.add(vars);
    t.pre_clustering_init();
    t.cluster_and_analyse();
}
fn loh(chr: &str, s: i32, e: i32, ts: &str, te: &str, sv1: i32, sv2: i32) -> Loh {
    Loh::new(LohEvent::new(chr, s, e, ts, te, 1, sv1, sv2))
}
fn max_link_distance(chain: &Chain) -> i32 {
    chain.linked_pairs().iter().map(|p| p.position_distance()).max().unwrap()
}

// SyntheticTest

#[test]
fn test_synthetic_del_dup_from_inv_pairs() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 100, 200, 1);
    let var2 = create_inv(t.next_var_id(), "1", 300, 400, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::PAIR_OTHER, c.resolved_type());

    let var1 = create_inv(t.next_var_id(), "1", 100, 300, 1);
    let var2 = create_inv(t.next_var_id(), "1", 200, 400, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL, c.resolved_type());
    assert_eq!(vec![AC::DSB], arm_cluster_types(&c));
    assert_eq!(1, arm_cluster_ti_count(&c, 0));

    let var1 = create_inv(t.next_var_id(), "1", 100, 2400, 1);
    let var2 = create_inv(t.next_var_id(), "1", 90, 3000, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_INV, c.resolved_type());

    let var1 = create_inv(t.next_var_id(), "1", 100, 10400, 1);
    let var2 = create_inv(t.next_var_id(), "1", 250, 10350, -1);
    set_assembled_link_info(&var1, false, &var2, false);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var2.position(true) - var1.position(true));
    assert_eq!(get_synthetic_ti_length(&c), var1.position(false) - var2.position(false) + 1);
    assert_eq!(vec![AC::DSB, AC::TI_ONLY], arm_cluster_types(&c));

    let var1 = create_inv(t.next_var_id(), "1", 500, 4000, 1);
    let var2 = create_inv(t.next_var_id(), "1", 480, 5000, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_INV, c.resolved_type());

    let var1 = create_inv(t.next_var_id(), "1", 200, 10400, 1);
    let var2 = create_inv(t.next_var_id(), "1", 100, 10350, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DUP, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var1.position(true) - var2.position(true));
    assert_eq!(get_synthetic_ti_length(&c), var1.position(false) - var2.position(false) + 1);
    assert_eq!(vec![AC::COMPLEX_OTHER, AC::TI_ONLY], arm_cluster_types(&c));

    let var1 = create_inv(t.next_var_id(), "1", 2300, 2400, 1);
    let var2 = create_inv(t.next_var_id(), "1", 100, 250, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::FB_INV_PAIR, c.resolved_type());
    assert_eq!(vec![AC::COMPLEX_FOLDBACK], arm_cluster_types(&c));
}

#[test]
fn test_synthetic_del_dup_from_bnd_pairs() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 100, 1, "2", 200, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 200, -1, "2", 150, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var2.position(true) - var1.position(true));
    assert_eq!(get_synthetic_ti_length(&c), var1.position(false) - var2.position(false) + 1);

    let var1 = create_bnd(t.next_var_id(), "1", 200, 1, "2", 200, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 100, -1, "2", 150, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DUP, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var1.position(true) - var2.position(true));
    assert_eq!(get_synthetic_ti_length(&c), var1.position(false) - var2.position(false) + 1);
}

#[test]
fn test_synthetic_del_dup_from_dels_and_dups() {
    let mut t = LinxTester::new();
    let var1 = create_del(t.next_var_id(), "1", 300, 350);
    let var2 = create_dup(t.next_var_id(), "1", 100, 400);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DUP, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var1.position(true) - var2.position(true));
    assert_eq!(get_synthetic_ti_length(&c), var2.position(false) - var1.position(false) + 1);

    let var1 = create_del(t.next_var_id(), "1", 100, 350);
    let var2 = create_dup(t.next_var_id(), "1", 200, 400);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var2.position(true) - var1.position(true));
    assert_eq!(get_synthetic_ti_length(&c), var2.position(false) - var1.position(false) + 1);

    let var1 = create_dup(t.next_var_id(), "1", 100, 350);
    let var2 = create_del(t.next_var_id(), "1", 200, 400);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var2.position(false) - var1.position(false));
    assert_eq!(get_synthetic_ti_length(&c), var2.position(true) - var1.position(true) + 1);

    let var1 = create_dup(t.next_var_id(), "1", 100, 2200);
    let var2 = create_dup(t.next_var_id(), "1", 300, 2400);
    set_assembled_link_info(&var1, true, &var2, false);
    t.add_and_cluster(&var1, &var2);
    assert_eq!(2, t.clusters().len());

    let var1 = create_dup(t.next_var_id(), "1", 100, 300);
    let var2 = create_dup(t.next_var_id(), "1", 250, 400);
    set_assembled_link_info(&var1, false, &var2, true);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DUP, c.resolved_type());
}

#[test]
fn test_synthetic_dels_and_dups_with_longer_chains() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 1200, 1, "3", 2000, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 5000, -1, "3", 1600, -1);
    run(&mut t, &[&var1, &var2, &var3]);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL, c.resolved_type());
    assert_eq!(get_synthetic_length(&c), var3.position(true) - var1.position(true));
    assert_eq!(401, get_synthetic_ti_length(&c));

    let var4 = create_dup(t.next_var_id(), "1", 980, 5600);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DUP, c.resolved_type());
    assert!(c.is_synthetic_type());
    assert_eq!(get_synthetic_length(&c), var1.position(true) - var4.position(true));
    assert_eq!(601, get_synthetic_ti_length(&c));
}

#[test]
fn test_synthetic_incompletes() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 1200, 1, "3", 2000, 1);
    run(&mut t, &[&var1, &var2]);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::UNBAL_TRANS, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 1200, 1, "2", 1500, 1);
    run(&mut t, &[&var1, &var2]);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::INV, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_sgl(t.next_var_id(), "2", 1500, 1);
    run(&mut t, &[&var1, &var2]);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::SGL, c.resolved_type());
}

// InversionPairTest

#[test]
fn test_reciprocal_inversions() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 10000, 1);
    let var2 = create_inv(t.next_var_id(), "1", 1500, 10500, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_INV, c.resolved_type());

    let var1 = create_inv(t.next_var_id(), "1", 100, 10010, 1);
    let var2 = create_inv(t.next_var_id(), "1", 90, 10000, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_INV, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 100, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 10000, 1, "2", 200, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 1500, -1, "3", 100, -1);
    let var4 = create_bnd(t.next_var_id(), "1", 10500, -1, "3", 200, 1);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_INV, c.resolved_type());
    assert_eq!(1, c.chains().len());
}

#[test]
fn test_facing_inversions() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 2000, -1);
    let var2 = create_inv(t.next_var_id(), "1", 5000, 6000, 1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::FB_INV_PAIR, c.resolved_type());
    assert_eq!(1, c.chains().len());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 1500, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 2000, -1, "2", 1000, -1);
    let var3 = create_inv(t.next_var_id(), "1", 5000, 6000, 1);
    run(&mut t, &[&var1, &var2, &var3]);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::FB_INV_PAIR, c.resolved_type());
    assert_eq!(1, c.chains().len());
}

#[test]
fn test_inner_overlap_inversions() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 50000, 350000, 1);
    let var2 = create_inv(t.next_var_id(), "1", 1000, 200000, -1);
    t.cn_loader.loh_event_data.push(loh(var2.chromosome(true), 1, var2.position(true), TELOMERE, var2.type_str(), CN_DATA_NO_SV, var2.id()));
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::DUP_TI, c.resolved_type());
    assert_eq!(var1.position(false) - var2.position(false), get_synthetic_length(&c));
    assert_eq!(var1.position(true) - var2.position(true) + 1, get_synthetic_ti_length(&c));
    assert_eq!(var2.position(false) - var1.position(true), get_synthetic_gap_length(&c));

    t.cn_loader.loh_event_data.clear();
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::RECIP_INV_DUPS, c.resolved_type());
    assert_eq!(1, c.chains().len());
    let long_ti = max_link_distance(&c.chains()[0]);
    assert_eq!(long_ti, var1.position(false) - var2.position(true));

    let var1 = create_bnd(t.next_var_id(), "1", 50000, 1, "2", 100, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 350000, 1, "2", 200, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 1000, -1, "3", 100, -1);
    let var4 = create_bnd(t.next_var_id(), "1", 200000, -1, "3", 200, 1);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    // the Java test asserts on the earlier cluster here
    assert!(!c.is_resolved());
    assert_eq!(RT::RECIP_INV_DUPS, c.resolved_type());
    assert_eq!(1, c.chains().len());
    assert_eq!(long_ti, var2.position(true) - var3.position(true));
}

#[test]
fn test_enclosing_inversions() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 200000, 1);
    let var2 = create_inv(t.next_var_id(), "1", 5000, 150000, -1);
    t.cn_loader.loh_event_data.push(loh(var1.chromosome(false), var1.position(false), 1000000, var1.type_str(), TELOMERE, var1.id(), CN_DATA_NO_SV));
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL_TI, c.resolved_type());
    assert_eq!(var2.position(true) - var1.position(true), get_synthetic_length(&c));
    assert_eq!(var1.position(false) - var2.position(false) + 1, get_synthetic_ti_length(&c));
    assert_eq!(var2.position(false) - var2.position(true), get_synthetic_gap_length(&c));

    t.cn_loader.loh_event_data.clear();
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::RECIP_INV_DEL_DUP, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 100, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 200000, -1, "2", 200, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 5000, 1, "3", 100, -1);
    let var4 = create_bnd(t.next_var_id(), "1", 150000, 1, "3", 200, 1);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = t.clusters()[0].clone();
    assert!(!c.is_resolved());
    assert_eq!(RT::RECIP_INV_DEL_DUP, c.resolved_type());
    assert_eq!(1, c.chains().len());
    assert_eq!(max_link_distance(&c.chains()[0]), var4.position(true) - var1.position(true));
}

// TranslocationPairTest

#[test]
fn test_reciprocal_translocations() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 100, 1, "2", 200, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 120, -1, "2", 220, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_TRANS, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 100, 1, "2", 200, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 90, -1, "2", 190, -1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_TRANS, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 1200, 1, "3", 2000, -1);
    let var3 = create_bnd(t.next_var_id(), "1", 4000, -1, "3", 1600, 1);
    run(&mut t, &[&var1, &var2, &var3]);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_TRANS, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 1200, 1, "3", 2000, -1);
    let var3 = create_bnd(t.next_var_id(), "1", 4000, -1, "4", 10000, 1);
    let var4 = create_bnd(t.next_var_id(), "3", 1600, 1, "4", 9400, -1);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_TRANS, c.resolved_type());
    assert_eq!(2, c.chains().len());
    assert!(!c.chains().iter().any(|x| x.linked_pairs().iter().any(|y| y.base_length() > SHORT_TI_LENGTH)));
}

#[test]
fn test_facing_del_translocations() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 50001, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 5000, -1, "2", 150000, 1);
    t.cn_loader.loh_event_data.push(loh(var1.chromosome(false), 1, var1.position(false), TELOMERE, var1.type_str(), CN_DATA_NO_SV, var1.id()));
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::DEL_TI, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 50001, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 5000, -1, "2", 150000, 1);
    t.cn_loader.loh_event_data.clear();
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(c.is_resolved());
    assert_eq!(RT::RECIP_TRANS_DEL_DUP, c.resolved_type());
    assert!(c.chains().is_empty());
}

#[test]
fn test_facing_dup_translocations() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 50000, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 200000, 1, "2", 150000, 1);
    t.cn_loader.loh_event_data.push(loh(var1.chromosome(false), 1, var1.position(false), TELOMERE, var1.type_str(), CN_DATA_NO_SV, var1.id()));
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::DUP_TI, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 50000, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 200000, 1, "2", 150000, 1);
    t.cn_loader.loh_event_data.clear();
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::RECIP_TRANS_DUPS, c.resolved_type());
    assert!(c.chains().is_empty());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "3", 100, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 50000, -1, "3", 200, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 200000, 1, "4", 100, -1);
    let var4 = create_bnd(t.next_var_id(), "2", 150000, 1, "4", 200, 1);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = t.clusters()[0].clone();
    assert!(!c.is_resolved());
    assert_eq!(RT::RECIP_TRANS_DUPS, c.resolved_type());
    assert_eq!(2, c.chains().len());
}

#[test]
fn test_pair_other_translocations() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 10000, 1, "3", 1000, 1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::UNBAL_TRANS_TI, c.resolved_type());

    let var1 = create_bnd(t.next_var_id(), "1", 1000, -1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "2", 5000, -1, "3", 1000, 1);
    t.add_and_cluster(&var1, &var2);
    let c = only_cluster(&t);
    assert!(!c.is_resolved());
    assert_eq!(RT::PAIR_OTHER, c.resolved_type());
}

// FoldbackTest

fn foldback_case(t: &mut LinxTester, cluster_id: i32, vars: &[&Sv]) {
    t.clear_clusters_and_svs();
    let c = Cluster::create(cluster_id);
    for v in vars {
        c.add_variant(v);
    }
    t.analyser.clusters.push(c.clone());
    t.all_variants.extend(c.svs());
    t.pre_clustering_init();
}

#[test]
fn test_inv_foldbacks() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 100, 200, 1);
    foldback_case(&mut t, 0, &[&var1]);
    mark_foldbacks(&t.analyser.state.chr_breakend_map, false);
    assert_eq!(var1.id(), var1.foldback_id(true));
    assert_eq!(var1.id(), var1.foldback_id(false));
    assert_eq!(100, var1.foldback_length(true));

    let dup = create_dup(t.next_var_id(), "1", 190, 400);
    foldback_case(&mut t, 0, &[&var1, &dup]);
    mark_foldbacks(&t.analyser.state.chr_breakend_map, false);
    assert_eq!(var1.id(), var1.foldback_id(true));
    assert_eq!(var1.id(), var1.foldback_id(false));
    assert_eq!(100, var1.foldback_length(true));

    let var2 = create_inv(t.next_var_id(), "1", 100, 200, -1);
    let sgl = create_sgl(t.next_var_id(), "1", 110, 1);
    foldback_case(&mut t, 0, &[&var2, &sgl]);
    mark_foldbacks(&t.analyser.state.chr_breakend_map, false);
    assert_eq!(var2.id(), var2.foldback_id(true));
    assert_eq!(var2.id(), var2.foldback_id(false));
    assert_eq!(100, var2.foldback_length(true));

    let var2 = create_inv(t.next_var_id(), "1", 100, 200, -1);
    let del = create_del(t.next_var_id(), "1", 190, 400);
    foldback_case(&mut t, 1, &[&var2, &del]);
    let db = var2.db_link(false).unwrap();
    assert!(db.has_breakend_of(&var2, false));
    assert!(db.has_breakend_of(&del, true));
    mark_foldbacks(&t.analyser.state.chr_breakend_map, false);
    assert!(var2.foldback_breakend(true).is_none());
    assert!(var2.foldback_breakend(false).is_none());

    let var3 = create_inv(t.next_var_id(), "1", 100, 200, 1);
    let sgl = create_sgl(t.next_var_id(), "1", 90, -1);
    let sgl2 = create_sgl(t.next_var_id(), "1", 190, -1);
    foldback_case(&mut t, 2, &[&var3, &sgl, &sgl2]);
    let db = var3.db_link(true).unwrap();
    assert!(db.has_breakend_of(&var3, true));
    assert!(db.has_breakend_of(&sgl, true));
    let db = var3.db_link(false).unwrap();
    assert!(db.has_breakend_of(&var3, false));
    assert!(db.has_breakend_of(&sgl2, true));
    mark_foldbacks(&t.analyser.state.chr_breakend_map, false);
    assert!(var3.foldback_breakend(true).is_none());
    assert!(var3.foldback_breakend(false).is_none());
}

// FilteringTest

#[test]
fn test_duplicate_breakend_filtering() {
    let mut t = LinxTester::new();
    let var1 = create_del(t.next_var_id(), "1", 1000, 1100);
    let var2 = create_sgl(t.next_var_id(), "1", 2000, -1);
    let var3 = create_sgl(t.next_var_id(), "1", 2001, -1);
    let var4 = create_sgl(t.next_var_id(), "1", 3000, -1);
    let var6 = create_del(t.next_var_id(), "1", 5000, 6000);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var6]);
    assert_eq!(2, t.all_clusters().len());
    assert!(t.has_cluster_with_svs(&[&var1, &var3, &var4, &var6]));
    assert_eq!(RT::DUP_BE, t.find_cluster_with_svs(&[&var2]).unwrap().resolved_type());

    let var1 = create_sv(t.next_var_id(), "1", "1", 1000, 2000, 1, -1, DEL, "1234567890123456789012345678901");
    let var2 = create_bnd(t.next_var_id(), "1", 999, 1, "2", 100, -1);
    let var3 = create_bnd(t.next_var_id(), "1", 2001, -1, "2", 200, 1);
    set_assembled_link_info(&var2, false, &var3, false);
    run(&mut t, &[&var1, &var2, &var3]);
    assert_eq!(2, t.all_clusters().len());
    assert!(t.has_cluster_with_svs(&[&var2, &var3]));
    assert_eq!(RT::DUP_BE, t.find_cluster_with_svs(&[&var1]).unwrap().resolved_type());
}

#[test]
fn test_low_vaf_filtering() {
    let mut t = LinxTester::new();
    let var0 = create_del(t.next_var_id(), "1", 1000, 1100);
    let var1 = create_test_sv_cn(t.next_var_id(), "1", "1", 2000, 2050, -1, -1, INV, 1.0, 1.0, 0.1, 0.1, 0.1, "");
    let var2 = create_test_sv_cn(t.next_var_id(), "1", "2", 10000, 10000, -1, -1, BND, 1.0, 1.0, 0.1, 0.1, 0.1, crate::annotators::POLY_A_MOTIF);
    let var3 = create_test_sv_cn(t.next_var_id(), "1", "2", 20000, 20000, -1, -1, BND, 1.0, 1.0, 0.1, 0.1, 0.1, "");
    let var4 = create_test_sv_cn(t.next_var_id(), "1", "2", 30000, 38000, -1, -1, BND, 1.0, 1.0, 0.1, 0.1, 0.1, "");
    let var5 = create_test_sv_cn(t.next_var_id(), "1", "2", 40000, 40000, -1, -1, BND, 1.0, 1.0, 0.1, 0.1, 0.1, "");
    let var6 = create_test_sv_cn(t.next_var_id(), "1", "", 50000, -1, -1, 0, SGL, 1.0, 0.0, 0.1, 0.0, 0.1, "");
    t.add(&[&var0, &var1, &var2, &var3, &var4, &var5, &var6]);
    for v in &t.all_variants {
        v.set_jcn_recalc_data(0.1, 0.2);
    }
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(5, t.all_clusters().len());
    assert!(t.has_cluster_with_svs(&[&var0]));
    assert!(t.has_cluster_with_svs(&[&var2, &var4, &var5]));
    for v in [&var1, &var3, &var6] {
        assert_eq!(RT::LOW_VAF, t.find_cluster_with_svs(&[v]).unwrap().resolved_type());
    }
}

#[test]
fn test_inferred_pair_filtering() {
    let mut t = LinxTester::new();
    let inf = |t: &mut LinxTester, pos: i32, o: i32, cn_chg: f64, jcn: f64| {
        create_test_sv_cn(t.next_var_id(), "1", "", pos, 0, o, 0, INF, 1.0, 1.0, cn_chg, 0.0, jcn, "")
    };
    let v0 = inf(&mut t, 2000, 1, 1.0, 1.0);
    let v1 = inf(&mut t, 3000, -1, 1.0, 1.0);
    let v2 = inf(&mut t, 4000, -1, 1.0, 1.0);
    let v3 = inf(&mut t, 6000, 1, 1.0, 1.0);
    let v4 = inf(&mut t, 7000, -1, 1.0, 1.0);
    let v5 = inf(&mut t, 8000, 1, 2.0, 2.0);
    run(&mut t, &[&v0, &v1, &v2, &v3, &v4, &v5]);
    assert_eq!(5, t.all_clusters().len());
}

// CopyNumberTest

#[test]
fn test_ploidy_recalcs() {
    let r = calc_adjusted_jcn_values(3.1, 3.6, 58, 2.14, 6.01, 6.66, &[3, 3], Some(&[1, 1]));
    assert!((r.jcn_estimate - 2.93).abs() < 0.01);
    assert!((r.jcn_uncertainty - 0.83).abs() < 0.01);
    let r = calc_adjusted_jcn_values(2.0, 1.3, 69, 1.43, 7.0, 5.0, &[0, 3], Some(&[0, 1]));
    assert!((r.jcn_estimate - 1.48).abs() < 0.01);
    assert!((r.jcn_uncertainty - 1.69).abs() < 0.01);
}

#[test]
fn test_cluster_min_max_copy_number() {
    let mut t = LinxTester::new();
    let var1 = create_del(t.next_var_id(), "1", 100, 200);
    let var2 = create_del(t.next_var_id(), "1", 300, 400);
    let var3 = create_del(t.next_var_id(), "1", 500, 600);
    var1.set_jcn_recalc_data(0.9, 1.6);
    var2.set_jcn_recalc_data(0.1, 1.2);
    var3.set_jcn_recalc_data(0.99, 2.01);
    let c = Cluster::create(0);
    c.add_variant(&var1);
    c.add_variant(&var2);
    c.add_variant(&var3);
    assert!(!c.has_varied_jcn());
    assert!((c.jcn_range()[SE_START] - 1.0).abs() < 0.001);
    assert!((c.jcn_range()[SE_END] - 1.0).abs() < 0.001);
    let var4 = create_del(t.next_var_id(), "1", 700, 800);
    var4.set_jcn_recalc_data(1.1, 1.5);
    c.add_variant(&var4);
    assert!(!c.has_varied_jcn());
    assert!((c.jcn_range()[SE_START] - 1.0).abs() < 0.001);
    assert!((c.jcn_range()[SE_END] - 1.0).abs() < 0.001);
}

fn build_copy_number_data(t: &mut LinxTester) {
    t.analyser.state.reset();
    populate_chromosome_breakend_map(&t.all_variants, &mut t.analyser.state);
    t.populate_copy_number_data(false);
}

fn assert_cn(t: &LinxTester, chr: &str, i: usize, cn: f64, baf: Option<f64>) {
    let list = t.cn_loader.chr_cn_data_map.get(&chr.to_string()).unwrap();
    assert!((list[i].copy_number - cn).abs() < 0.01, "segment {i} copy number {}", list[i].copy_number);
    if let Some(b) = baf {
        assert!((list[i].actual_baf - b).abs() < 0.01, "segment {i} BAF {}", list[i].actual_baf);
    }
}

#[test]
fn test_copy_number_reconstruction() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 1000, 2000, 1, -1, DEL, 1.0);
    t.add(&[&var1]);
    build_copy_number_data(&mut t);
    assert_eq!(4, t.cn_loader.chr_cn_data_map.get(&"1".to_string()).unwrap().len());
    assert_cn(&t, "1", 0, 2.0, Some(0.5));
    assert_cn(&t, "1", 1, 1.0, Some(1.0));
    assert_cn(&t, "1", 2, 2.0, None);

    t.clear_clusters_and_svs();
    let var1 = create_test_sv(1, "1", "1", 1000, 2000, -1, 1, DUP, 2.0);
    t.add(&[&var1]);
    build_copy_number_data(&mut t);
    assert_eq!(4, t.cn_loader.chr_cn_data_map.get(&"1".to_string()).unwrap().len());
    assert_cn(&t, "1", 0, 3.0, Some(0.67));
    assert_cn(&t, "1", 1, 5.0, Some(0.8));
    assert_cn(&t, "1", 2, 3.0, None);

    t.clear_clusters_and_svs();
    let var1 = create_test_sv(1, "1", "1", 1000, 2000, -1, -1, INV, 2.0);
    let var2 = create_test_sv(2, "1", "1", 5000, 6000, 1, 1, INV, 1.0);
    let var3 = create_test_sv(3, "1", "2", 8000, 100, 1, 1, BND, 1.0);
    t.add(&[&var1, &var2, &var3]);
    build_copy_number_data(&mut t);
    assert_eq!(7, t.cn_loader.chr_cn_data_map.get(&"1".to_string()).unwrap().len());
    assert_cn(&t, "1", 0, 1.0, Some(1.0));
    assert_cn(&t, "1", 1, 3.0, Some(0.67));
    assert_cn(&t, "1", 2, 5.0, Some(0.8));
    assert_cn(&t, "1", 3, 4.0, Some(0.75));
}

// AnnotationTest

#[test]
fn test_local_topology() {
    let mut t = LinxTester::new();
    let var1 = create_inv(1, "1", 101000, 107000, 1);
    let var2 = create_bnd(2, "1", 149000, 1, "2", 1000, -1);
    let var3 = create_bnd(3, "1", 250000, 1, "2", 1100, 1);
    let var4 = create_dup(4, "1", 150000, 200500);
    let var5 = create_inv(5, "1", 201000, 204000, -1);
    let var6 = create_inv(6, "1", 301000, 304000, 1);
    let var7 = create_del(7, "1", 305000, 306000);
    let var8 = create_inv(8, "1", 310000, 315000, -1);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8]);
    for v in [&var1, &var5, &var6, &var8] {
        assert!(v.is_foldback());
    }
    let c = only_cluster(&t);
    assert_eq!(6, c.arm_clusters.borrow().len());
    for at in [AC::ISOLATED_BE, AC::TI_ONLY, AC::DSB, AC::FOLDBACK, AC::FOLDBACK_DSB, AC::COMPLEX_FOLDBACK] {
        assert_eq!(1, arm_cluster_data(&c, at), "{at}");
    }
}

#[test]
fn test_multiple_tis() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 1000, 1, "2", 1000, -1);
    let var2 = create_bnd(t.next_var_id(), "1", 1110, -1, "2", 1200, 1);
    let var3 = create_bnd(t.next_var_id(), "1", 1310, 1, "2", 1210, -1);
    let var4 = create_bnd(t.next_var_id(), "1", 1320, -1, "2", 1410, 1);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = only_cluster(&t);
    assert_eq!(2, c.arm_clusters.borrow().len());
    assert_eq!(1, arm_cluster_data(&c, AC::TI_ONLY));
    assert_eq!(1, arm_cluster_data(&c, AC::DSB));
    let tis: i32 = c.arm_clusters.borrow().iter().map(|x| x.ti_count).sum();
    assert_eq!(3, tis);
    assert!(c.arm_clusters.borrow().iter().any(|x| x.ti_count == 2 && x.arm_cluster_type == AC::TI_ONLY));
}

#[test]
fn test_cluster_metrics() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 1000, 49000, 1);
    let var2 = create_inv(t.next_var_id(), "1", 10000, 50000, -1);
    let del1 = create_del(t.next_var_id(), "1", 20000, 20200);
    let del2 = create_del(t.next_var_id(), "1", 30000, 30200);
    let bnd1 = create_bnd(t.next_var_id(), "1", 52000, 1, "2", 100, -1);
    let bnd2 = create_bnd(t.next_var_id(), "1", 52050, -1, "2", 200, 1);
    let del3 = create_del(t.next_var_id(), "1", 60000, 60100);
    let bnd3 = create_bnd(t.next_var_id(), "1", 74000, 1, "2", 10000, -1);
    let bnd4 = create_bnd(t.next_var_id(), "1", 75000, -1, "2", 11000, 1);
    let var3 = create_del(t.next_var_id(), "1", 76000, 90000);
    let var4 = create_dup(t.next_var_id(), "1", 77000, 90100);
    run(&mut t, &[&var1, &var2, &var3, &var4, &del1, &del2, &del3, &bnd1, &bnd2, &bnd3, &bnd4]);
    assert_eq!(4, t.clusters().len());
    let c = t.clusters()[0].clone();
    let mut m = c.metrics.borrow().clone();
    m.populate(&c, &t.analyser.state.chr_breakend_map);
    assert_eq!(5, m.db_count);
    assert_eq!(1, m.short_db_count);
    assert_eq!(76000, m.total_range);
    assert_eq!(3, m.traversed_del_count);
    assert_eq!(1, m.traversed_short_del_count);
    assert_eq!(5, m.implied_ti_count);
    assert_eq!(500, m.traversed_del_length);
}

// ReplicationRepairTest

#[test]
fn test_replication_repair_one_break_dup_inversions() {
    let mut t = LinxTester::new();
    let var1 = create_inv(1, "1", 1000, 10000, 1);
    let var2 = create_dup(2, "1", 1100, 9500);
    let var3 = create_inv(3, "1", 1500, 10100, -1);
    run(&mut t, &[&var1, &var2, &var3]);
    let c = only_cluster(&t);
    t.analyser.annotate_clusters();
    assert_eq!(1, c.chains().len());
    assert!(c.has_annotation(CLUSTER_ANNOT_REP_REPAIR));
}

#[test]
fn test_replication_repair_one_break_foldback_inversion() {
    let mut t = LinxTester::new();
    let var1 = create_inv(1, "1", 1000, 10000, 1);
    let var2 = create_inv(2, "1", 1100, 1500, -1);
    let var3 = create_del(3, "1", 9500, 10100);
    run(&mut t, &[&var1, &var2, &var3]);
    assert!(var2.is_foldback());
    let c = only_cluster(&t);
    t.analyser.annotate_clusters();
    assert_eq!(1, c.chains().len());
    assert!(c.has_annotation(CLUSTER_ANNOT_REP_REPAIR));
}

// LineTest

use crate::annotators::{contains_motif, has_line_poly_a_or_t_motif, POLY_A_MOTIF, POLY_T_MOTIF};
use LineElementType::{KNOWN, SUSPECT};

const SHORT_BASES: &str = "GCTG";
const RANDOM_BASES: &str = "GCTGTGAGCTGATCG";

fn sgl_mapping(chr: &str, pos: i32, o: i8) -> SglMapping {
    SglMapping { chromosome: chr.into(), position: pos, orientation: o, cigar: String::new(), qual_score: 1 }
}

#[test]
fn test_known_line_marking() {
    let mut t = LinxTester::new();
    let bnd1 = create_bnd(t.next_var_id(), "1", 100, -1, "2", 100, 1);
    let bnd2 = create_bnd(t.next_var_id(), "1", 200, 1, "2", 200, -1);
    let sgl = create_sgl(t.next_var_id(), "1", 500, 1);
    let inf = create_inf(t.next_var_id(), "1", 1500, 1);
    bnd1.add_line_element(KNOWN, true);
    bnd2.add_line_element(KNOWN, true);
    let c = Cluster::create(0);
    for v in [&bnd1, &bnd2, &sgl, &inf] {
        c.add_variant(v);
    }
    t.analyser.clusters.push(c.clone());
    t.all_variants.extend(c.svs());
    LineElementAnnotator::new(5000).mark_line_cluster(&c, None);
    assert!(c.has_linking_line_elements());
}

#[test]
fn test_line_motifs() {
    let a8 = &POLY_A_MOTIF[..8];
    assert!(contains_motif(&format!("{a8}{SHORT_BASES}{a8}"), POLY_A_MOTIF));
    assert!(!contains_motif(&format!("{}{SHORT_BASES}{a8}", &POLY_A_MOTIF[..7]), POLY_A_MOTIF));
    assert!(!has_line_poly_a_or_t_motif(SHORT_BASES, 1, true, false));
    assert!(!has_line_poly_a_or_t_motif(RANDOM_BASES, 1, true, false));

    let ins = format!("{SHORT_BASES}{POLY_A_MOTIF}{RANDOM_BASES}");
    assert!(has_line_poly_a_or_t_motif(&ins, 1, true, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, false, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, true, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, false, false));

    let ins = format!("{SHORT_BASES}{POLY_T_MOTIF}{RANDOM_BASES}");
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, true, false));
    assert!(has_line_poly_a_or_t_motif(&ins, 1, false, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, true, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, false, false));

    let ins = format!("{RANDOM_BASES}{POLY_T_MOTIF}{SHORT_BASES}");
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, true, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, false, false));
    assert!(has_line_poly_a_or_t_motif(&ins, -1, true, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, false, false));

    let ins = format!("{RANDOM_BASES}{POLY_A_MOTIF}{SHORT_BASES}");
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, true, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, false, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, true, false));
    assert!(has_line_poly_a_or_t_motif(&ins, -1, false, false));
    assert!(!has_line_poly_a_or_t_motif(&ins, 1, true, false));
    assert!(has_line_poly_a_or_t_motif(&ins, 1, true, true));

    let ins = format!("{SHORT_BASES}{POLY_T_MOTIF}{RANDOM_BASES}");
    assert!(!has_line_poly_a_or_t_motif(&ins, -1, true, false));
    assert!(has_line_poly_a_or_t_motif(&ins, -1, true, true));
}

#[test]
fn test_suspect_bnd_line_marking() {
    let mut t = LinxTester::new();
    let poly_a = format!("{SHORT_BASES}{POLY_A_MOTIF}{RANDOM_BASES}");
    let poly_t = format!("{RANDOM_BASES}{POLY_T_MOTIF}{SHORT_BASES}");
    let insert_poly_a = format!("{RANDOM_BASES}{POLY_A_MOTIF}{SHORT_BASES}");

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 150, -1, 1, BND, &poly_t);
    let bnd2 = create_sv(t.next_var_id(), "1", "2", 200, 170, 1, -1, BND, "");
    t.add_and_cluster(&bnd1, &bnd2);
    assert!(bnd1.has_line_element(SUSPECT, true));
    assert!(bnd2.has_line_element(SUSPECT, true));
    assert!(!bnd1.is_line_element(false));
    assert!(!bnd2.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 100, -1, 1, BND, &poly_t);
    let bnd2 = create_sv(t.next_var_id(), "1", "3", 200, 110, 1, -1, BND, &poly_a);
    t.add_and_cluster(&bnd1, &bnd2);
    assert!(bnd1.has_line_element(SUSPECT, true));
    assert!(bnd2.has_line_element(SUSPECT, true));
    assert!(!bnd1.is_line_element(false));
    assert!(!bnd2.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 100, -1, 1, BND, &poly_t);
    let sgl = create_sv(t.next_var_id(), "2", "", 90, -1, -1, 0, SGL, &insert_poly_a);
    t.add_and_cluster(&bnd1, &sgl);
    assert!(bnd1.has_line_element(SUSPECT, true));
    assert!(!sgl.is_line_element(true));
    assert!(!bnd1.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 100, -1, 1, BND, &poly_t);
    let bnd2 = create_sv(t.next_var_id(), "1", "3", 200, 110, 1, -1, BND, "");
    t.add_and_cluster(&bnd1, &bnd2);
    assert!(bnd1.has_line_element(SUSPECT, true));
    assert!(bnd2.has_line_element(SUSPECT, true));
    assert!(!bnd1.is_line_element(false));
    assert!(!bnd2.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 100, 1, 1, BND, "");
    let bnd2 = create_sv(t.next_var_id(), "1", "2", 110, 110, -1, -1, BND, &poly_t);
    t.add_and_cluster(&bnd1, &bnd2);
    for v in [&bnd1, &bnd2] {
        assert!(!v.is_line_element(true));
        assert!(!v.is_line_element(false));
    }
    assert!(!t.clusters()[0].has_linking_line_elements());

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 100, -1, 1, BND, "");
    let bnd2 = create_sv(t.next_var_id(), "1", "3", 150, 200, -1, -1, BND, "");
    let del = create_sv(t.next_var_id(), "1", "1", 1000, 11000, 1, -1, DEL, &poly_a);
    run(&mut t, &[&bnd1, &bnd2, &del]);
    assert!(bnd1.has_line_element(SUSPECT, true));
    assert!(bnd2.has_line_element(SUSPECT, true));
    assert!(del.has_line_element(SUSPECT, true));
    assert!(!del.is_line_element(false));
    assert!(!bnd1.is_line_element(false));
    assert!(!bnd2.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let insert_poly_t = format!("{SHORT_BASES}{POLY_T_MOTIF}{RANDOM_BASES}");
    let sgl = create_sv(t.next_var_id(), "1", "", 100, 0, 1, 0, SGL, &insert_poly_t);
    let bnd1 = create_sv(t.next_var_id(), "1", "2", 110, 100, -1, 1, BND, "");
    let bnd2 = create_sv(t.next_var_id(), "2", "3", 200, 1000, 1, -1, BND, "");
    run(&mut t, &[&bnd1, &bnd2, &sgl]);
    assert!(bnd1.has_line_element(SUSPECT, false));
    assert!(bnd2.has_line_element(SUSPECT, true));
    assert!(!sgl.in_line_element());
    assert!(!bnd1.is_line_element(true));
    assert!(!bnd2.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let bnd1 = create_sv(t.next_var_id(), "1", "2", 100, 100, -1, 1, BND, "");
    let bnd2 = create_sv(t.next_var_id(), "1", "2", 200, 110, 1, -1, BND, &poly_a);
    bnd1.add_line_element(KNOWN, false);
    t.add_and_cluster(&bnd1, &bnd2);
    assert!(!bnd1.has_line_element(SUSPECT, true));
    assert!(!bnd2.has_line_element(SUSPECT, true));
    assert!(!bnd1.is_line_element(true));
    assert!(!bnd2.is_line_element(false));
    assert!(!t.clusters()[0].has_linking_line_elements());

    let inv = create_sv(t.next_var_id(), "1", "1", 100, 1050000, 1, 1, INV, &poly_a);
    let sgl = create_sv(t.next_var_id(), "1", "0", 1050020, 0, -1, 0, SGL, &poly_a);
    t.add_and_cluster(&inv, &sgl);
    assert!(inv.has_line_element(SUSPECT, true));
    assert!(!inv.has_line_element(SUSPECT, false));
    assert!(t.clusters()[0].has_linking_line_elements());
}

#[test]
fn test_suspect_sgl_line_marking() {
    let mut t = LinxTester::new();
    let insert_poly_a = format!("{RANDOM_BASES}{POLY_A_MOTIF}{SHORT_BASES}");
    let sgl1 = create_sv(t.next_var_id(), "1", "", 100, 0, 1, 0, SGL, "");
    let sgl2 = create_sv(t.next_var_id(), "1", "", 120, 0, -1, 0, SGL, &insert_poly_a);
    t.add_and_cluster(&sgl1, &sgl2);
    assert!(!sgl1.is_line_element(false));
    assert!(!sgl2.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());

    let sgl1 = create_sv(t.next_var_id(), "1", "", 120, 0, -1, 0, SGL, &insert_poly_a);
    run(&mut t, &[&sgl1]);
    assert!(!sgl1.is_line_element(false));
    assert!(!t.clusters()[0].has_linking_line_elements());

    let sgl1 = create_test_sv_cn(t.next_var_id(), "1", "", 120, 0, -1, 0, SGL, 2.0, 0.0, 0.1, 0.0, 0.1, &insert_poly_a);
    run(&mut t, &[&sgl1]);
    assert!(!sgl1.is_line_element(false));
    assert!(t.clusters()[0].has_linking_line_elements());
}

#[test]
fn test_line_chaining() {
    let mut t = LinxTester::new();
    let var1 = create_bnd(t.next_var_id(), "1", 100, -1, "2", 100, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 200, 1, "2", 200, -1);
    let var3 = create_bnd(t.next_var_id(), "1", 1100, -1, "3", 100, 1);
    let var4 = create_bnd(t.next_var_id(), "1", 1200, 1, "3", 200, -1);
    for v in [&var1, &var2, &var3, &var4] {
        v.add_line_element(KNOWN, true);
    }
    set_assembled_link_info(&var1, true, &var2, true);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    let c = only_cluster(&t);
    assert_eq!(2, c.chains().len());

    let var1 = create_bnd(t.next_var_id(), "1", 100, -1, "3", 100, 1);
    let var2 = create_bnd(t.next_var_id(), "1", 200, 1, "3", 200, -1);
    let var3 = create_bnd(t.next_var_id(), "1", 300, -1, "2", 100, 1);
    let var4 = create_inv(t.next_var_id(), "1", 400, 600, 1);
    let var5 = create_bnd(t.next_var_id(), "1", 500, -1, "2", 200, -1);
    set_assembled_link_info(&var3, true, &var4, true);
    let var6 = create_bnd(t.next_var_id(), "3", 1000, -1, "4", 100, 1);
    let var7 = create_bnd(t.next_var_id(), "3", 1100, 1, "4", 200, -1);
    let var8 = create_bnd(t.next_var_id(), "1", 2000, -1, "5", 100, 1);
    let var9 = create_bnd(t.next_var_id(), "1", 2100, 1, "6", 200, -1);
    let all = [&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8, &var9];
    for v in all {
        v.add_line_element(KNOWN, true);
    }
    run(&mut t, &all);
    let c = only_cluster(&t);
    assert_eq!(3, c.chains().len());
    assert!(t.find_chain_with_svs(&c, &[&var1, &var2]).is_some());
    assert!(t.find_chain_with_svs(&c, &[&var3, &var4, &var5]).is_some());
    assert!(t.find_chain_with_svs(&c, &[&var6, &var7]).is_some());
}

#[test]
fn test_line_sgl_chaining() {
    let mut t = LinxTester::new();
    let insert_poly_t = format!("{SHORT_BASES}{POLY_T_MOTIF}{RANDOM_BASES}");
    let sgl1 = create_sv(t.next_var_id(), "1", "", 100, 0, 1, 0, SGL, &insert_poly_t);
    let sgl2 = create_sv(t.next_var_id(), "1", "", 120, 0, -1, 0, SGL, "");
    sgl1.add_sgl_mapping(sgl_mapping(CHR_2, 1000, -1));
    sgl1.add_sgl_mapping(sgl_mapping("5", 10000, -1));
    sgl2.add_sgl_mapping(sgl_mapping(CHR_2, 1100, 1));
    sgl2.add_sgl_mapping(sgl_mapping("3", 2000, -1));
    run(&mut t, &[&sgl1, &sgl2]);
    assert_eq!(1, only_cluster(&t).chains().len());

    let sgl1 = create_sv(t.next_var_id(), "1", "", 100, 0, 1, 0, SGL, &insert_poly_t);
    let bnd = create_sv(t.next_var_id(), "1", "2", 120, 1100, -1, 1, BND, "");
    bnd.add_line_element(KNOWN, false);
    sgl1.add_sgl_mapping(sgl_mapping(CHR_2, 1000, -1));
    sgl1.add_sgl_mapping(sgl_mapping("5", 10000, -1));
    run(&mut t, &[&sgl1, &bnd]);
    let c = only_cluster(&t);
    assert!(c.has_linking_line_elements());
    assert_eq!(1, c.chains().len());
}

// DoubleMinuteTest

use crate::analysis::double_minute::SEG_DATA_COUNT;

#[test]
fn test_simple_dup_dm() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 500, 600, 1, -1, DEL, 1.0);
    let dup = create_test_sv(2, "1", "1", 50000, 55000, -1, 1, DUP, 10.0);
    dup.set_jcn_recalc_data(8.0, 12.0);
    run(&mut t, &[&var1, &dup]);
    assert_eq!(2, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&dup]).unwrap();
    assert_eq!(1, c.sv_count());
    assert_eq!(RT::DOUBLE_MINUTE, c.resolved_type());
    assert_eq!(1, c.chains().len());
    assert_eq!(1, c.chains()[0].link_count());
    assert!(c.chains()[0].is_closed_loop());
}

#[test]
fn test_chained_dm() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 1000, 6000, -1, -1, INV, 8.0);
    let var2 = create_test_sv(2, "1", "1", 3000, 8000, 1, 1, INV, 8.0);
    run(&mut t, &[&var1, &var2]);
    let c = t.find_cluster_with_svs(&[&var1, &var2]).unwrap();
    assert_eq!(RT::DOUBLE_MINUTE, c.resolved_type());
    assert!(contains(&c.double_minute_svs(), &var1));
    assert!(contains(&c.double_minute_svs(), &var2));
    assert_eq!(1, c.chains().len());
    let chain = c.chains()[0].clone();
    assert_eq!(2, chain.link_count());
    assert_eq!(2, chain.sv_count());
    assert!(chain.is_closed_loop());
}

#[test]
fn test_invalid_open_breakend_dms() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(t.next_var_id(), "1", "2", 1001, 100, -1, 1, BND, 8.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 20000, 200, 1, -1, BND, 8.0);
    let var3 = create_test_sv(t.next_var_id(), "1", "1", 2000, 10001, 1, -1, DEL, 8.0);
    let var4 = create_test_sv(t.next_var_id(), "1", "1", 11000, 16000, 1, -1, DEL, 8.0);
    run(&mut t, &[&var1, &var2, &var3, &var4]);
    assert_eq!(1, t.analyser.dm_finder.double_minutes.len());
    let dm = t.analyser.dm_finder.double_minutes.get(&0).unwrap();
    assert_eq!(2, dm.simple_dels);
    assert_eq!(2, dm.open_breakends);
    assert_eq!(2, dm.closed_breakends);
    assert_eq!(RT::RECIP_TRANS_DEL_DUP, only_cluster(&t).resolved_type());
}

#[test]
fn test_internal_external_breakend_dms() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(t.next_var_id(), "1", "2", 1000, 1000, -1, -1, BND, 8.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 10000, 10000, 1, 1, BND, 8.0);
    run(&mut t, &[&var1, &var2]);
    let c = t.find_cluster_with_svs(&[&var1, &var2]).unwrap();
    assert_eq!(RT::DOUBLE_MINUTE, c.resolved_type());
    assert!(contains(&c.double_minute_svs(), &var1));
    assert!(contains(&c.double_minute_svs(), &var2));

    let var1 = create_test_sv(t.next_var_id(), "1", "2", 1000, 1000, -1, -1, BND, 8.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 10000, 10000, 1, 1, BND, 8.0);
    let var3 = create_test_sv(t.next_var_id(), "1", "3", 2000, 100, 1, -1, BND, 1.5);
    let var4 = create_test_sv(t.next_var_id(), "1", "3", 2100, 200, -1, 1, BND, 1.5);
    let var5 = create_test_sv(t.next_var_id(), "1", "4", 3000, 100, -1, 1, BND, 1.5);
    let var6 = create_test_sv(t.next_var_id(), "1", "4", 3100, 200, 1, -1, BND, 1.5);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6]);
    assert_eq!(1, t.analyser.dm_finder.double_minutes.len());
    let dm = t.analyser.dm_finder.double_minutes.values().next().unwrap();
    assert_eq!(1, dm.chain_int_ext_data.len());
    assert_eq!(4.0, dm.chain_int_ext_data.get(&dm.chains[0].id()).unwrap()[SEG_DATA_COUNT]);
    assert_eq!(RT::COMPLEX, only_cluster(&t).resolved_type());

    let var1 = create_test_sv(t.next_var_id(), "1", "2", 1000, 1000, -1, -1, BND, 8.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 10000, 10000, 1, 1, BND, 8.0);
    let var3 = create_test_sv(t.next_var_id(), "1", "3", 2000, 100, 1, -1, BND, 1.5);
    let var4 = create_test_sv(t.next_var_id(), "1", "3", 2100, 200, -1, 1, BND, 1.5);
    set_assembled_link_info(&var3, false, &var4, false);
    let var5 = create_test_sv(t.next_var_id(), "1", "4", 3000, 100, -1, 1, BND, 1.5);
    let var6 = create_test_sv(t.next_var_id(), "1", "4", 3100, 200, 1, -1, BND, 1.5);
    set_assembled_link_info(&var5, true, &var6, true);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6]);
    assert_eq!(1, t.analyser.dm_finder.double_minutes.len());
    let dm = t.analyser.dm_finder.double_minutes.values().next().unwrap();
    assert!(dm.chain_int_ext_data.is_empty());
    let c = only_cluster(&t);
    assert_eq!(RT::COMPLEX, c.resolved_type());
    assert!(c.has_annotation(CLUSTER_ANNOT_DM));
    assert_eq!(2, c.double_minute_svs().len());
}

#[test]
fn test_invalid_chained_foldback_dms() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(t.next_var_id(), "1", "1", 1000, 1100, -1, -1, INV, 8.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 5000, 1000, 1, -1, BND, 16.0);
    let var3 = create_test_sv(t.next_var_id(), "2", "2", 5100, 5000, 1, 1, INV, 8.0);
    let var4 = create_test_sv(t.next_var_id(), "1", "3", 2000, 1000, 1, 1, BND, 2.0);
    let var5 = create_test_sv(t.next_var_id(), "1", "4", 3000, 1000, -1, 1, BND, 2.0);
    let var6 = create_test_sv(t.next_var_id(), "1", "5", 4000, 1000, 1, 1, BND, 2.0);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6]);
    assert_eq!(RT::COMPLEX, only_cluster(&t).resolved_type());
}

#[test]
fn test_multiple_dms_in_cluster() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 1000, 10000, -1, 1, DUP, 8.0);
    let var2 = create_test_sv(2, "1", "1", 2000, 3000, 1, -1, DEL, 1.0);
    set_assembled_link_info(&var1, true, &var2, true);
    let var3 = create_test_sv(3, "1", "1", 12000, 12500, -1, -1, INV, 8.0);
    let var4 = create_test_sv(4, "1", "1", 20000, 20500, 1, 1, INV, 8.0);
    let var5 = create_test_sv(5, "1", "1", 13000, 14000, 1, -1, DEL, 1.0);
    let var6 = create_test_sv(6, "1", "2", 15000, 100, 1, -1, BND, 1.0);
    let var7 = create_test_sv(7, "1", "2", 16000, 200, -1, 1, BND, 1.0);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6, &var7]);
    let c = only_cluster(&t);
    assert_eq!(RT::COMPLEX, c.resolved_type());
    assert!(c.has_annotation(CLUSTER_ANNOT_DM));
    for v in [&var1, &var3, &var4] {
        assert!(contains(&c.double_minute_svs(), v));
    }
    assert_eq!(2, c.chains().iter().filter(|x| x.is_double_minute()).count());

    let var1 = create_test_sv(t.next_var_id(), "1", "1", 25000, 35000, -1, 1, DUP, 8.0);
    let var2 = create_test_sv(t.next_var_id(), "1", "2", 100, 1100, -1, 1, BND, 8.0);
    let var3 = create_test_sv(t.next_var_id(), "1", "2", 10000, 1000, 1, -1, BND, 8.0);
    let var4 = create_test_sv(t.next_var_id(), "1", "1", 1000, 20000, 1, 1, INV, 1.5);
    let var5 = create_test_sv(t.next_var_id(), "1", "1", 2000, 21000, -1, -1, INV, 1.5);
    let var6 = create_test_sv(t.next_var_id(), "1", "1", 3000, 22000, 1, 1, INV, 1.5);
    let var7 = create_test_sv(t.next_var_id(), "1", "1", 36000, 36500, -1, 1, DUP, 8.0);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6, &var7]);
    let c = only_cluster(&t);
    assert_eq!(RT::COMPLEX, c.resolved_type());
    assert!(c.has_annotation(CLUSTER_ANNOT_DM));
    assert_eq!(1, c.double_minute_svs().len());
    assert!(contains(&c.double_minute_svs(), &var1));
    assert_eq!(1, c.chains().iter().filter(|x| x.is_double_minute()).count());
}
