//! chaining: ChainingSimpleTest, LinkingTest, ChainRoutinesTest and ChainingComplexTest (ChainingActualTest
//! is left out: its fixtures are real sample data)

use hmf_common::sv::StructuralVariantType::*;

use super::*;
use crate::chaining::utils::*;

fn first_cluster(t: &LinxTester) -> Cluster {
    assert_eq!(1, t.clusters().len());
    t.clusters()[0].clone()
}
fn only_chain(c: &Cluster) -> Chain {
    assert_eq!(1, c.chains().len());
    c.chains()[0].clone()
}
fn pair(b1: Be, b2: Be) -> Pair {
    Pair::from_breakends(&b1, &b2)
}

// ChainingSimpleTest

#[test]
fn test_fully_assembled_chain() {
    let mut t = LinxTester::new();
    let var1 = create_del(0, "1", 100, 200);
    let var2 = create_del(1, "1", 300, 400);
    let var3 = create_del(2, "1", 500, 600);
    let var4 = create_del(3, "1", 700, 800);
    set_assembled_link_info(&var1, false, &var2, true);
    set_assembled_link_info(&var2, false, &var3, true);
    set_assembled_link_info(&var3, false, &var4, true);
    t.add(&[&var1, &var2, &var3, &var4]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let c = first_cluster(&t);
    let chain = only_chain(&c);
    assert_eq!(3, chain.link_count());
    assert_eq!(3, c.assembly_linked_pairs().len());
}

#[test]
fn test_chain_not_closed() {
    let mut t = LinxTester::new();
    let var1 = create_inv(1, "1", 100, 200, -1);
    let var2 = create_del(2, "1", 300, 400);
    let var3 = create_del(3, "1", 500, 600);
    let var4 = create_inv(4, "1", 700, 800, 1);
    t.add(&[&var1, &var2, &var3, &var4]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let c = first_cluster(&t);
    assert_eq!(3, only_chain(&c).link_count());
}

#[test]
fn test_partially_assembled_chain() {
    let mut t = LinxTester::new();
    let var0 = create_del(0, "1", 100, 200);
    let var1 = create_del(1, "1", 300, 400);
    let var2 = create_del(2, "1", 500, 600);
    let var3 = create_del(3, "1", 700, 800);
    set_assembled_link_info(&var1, false, &var2, true);
    t.add(&[&var0, &var1, &var2, &var3]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let c = first_cluster(&t);
    let chain = only_chain(&c);
    assert_eq!(3, chain.link_count());
    assert_eq!(4, chain.sv_count());
}

#[test]
fn test_skip_spanning_assemblies() {
    let mut t = LinxTester::new();
    let var0 = create_del(0, "1", 100, 200);
    let var1 = create_del(1, "1", 300, 400);
    let var2 = create_del(2, "1", 500, 600);
    let var3 = create_del(3, "1", 700, 800);
    let var4 = create_del(4, "1", 900, 1000);
    set_assembled_link_info(&var0, false, &var1, true);
    set_assembled_link_info(&var1, false, &var2, true);
    set_assembled_link_info(&var2, false, &var3, true);
    set_assembled_link_info(&var3, false, &var4, true);
    set_assembled_link_info(&var0, false, &var2, true);
    set_assembled_link_info(&var2, false, &var4, true);
    t.add(&[&var0, &var1, &var2, &var3, &var4]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let c = first_cluster(&t);
    assert!(!c.requires_replication());
    let chain = only_chain(&c);
    assert_eq!(1, var0.assembled_linked_pairs(false).len());
    assert_eq!(1, var2.assembled_linked_pairs(true).len());
    assert_eq!(1, var2.assembled_linked_pairs(false).len());
    assert_eq!(1, var4.assembled_linked_pairs(true).len());
    assert_eq!(4, chain.link_count());
    assert_eq!(5, chain.sv_count());
}

#[test]
fn test_zero_base_inferred_chains() {
    let mut t = LinxTester::new();
    let var0 = create_test_sv(1, "1", "1", 1000, 2000, 1, -1, DEL, 1.0);
    let var1 = create_test_sv(1, "1", "0", 21000, -1, 1, 0, INF, 1.0);
    let var2 = create_test_sv(2, "1", "0", 21000, -1, -1, 0, SGL, 1.0);
    let var3 = create_test_sv(3, "1", "2", 40000, 100, 1, -1, BND, 4.0);
    let var4 = create_test_sv(4, "1", "0", 40000, -1, -1, 1, INF, 4.0);
    let var5 = create_test_sv(5, "1", "2", 44000, 200, -1, 1, BND, 4.0);
    let var6 = create_test_sv(6, "1", "0", 45000, -1, 1, 1, SGL, 4.0);
    t.add(&[&var0, &var1, &var2, &var3, &var4, &var5, &var6]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    assert_eq!(3, t.clusters().len());
    let c = t.find_cluster_with_svs(&[&var1, &var2]).unwrap();
    assert!(c.chains().is_empty());
    let c = t.find_cluster_with_svs(&[&var3, &var4, &var5, &var6]).unwrap();
    assert_eq!(1, c.chains().len());
    assert_eq!(3, c.chains()[0].link_count());
    let p = var4.linked_pair(true).unwrap();
    assert!(p.has_variant(&var3));
}

// LinkingTest

#[test]
fn test_linked_pairs() {
    let var1 = create_dup(1, "1", 100, 200);
    let var2 = create_dup(2, "1", 300, 400);
    let lp1 = Pair::from_svs(&var1, &var2, false, true);
    assert!(lp1.first() == var1);
    assert!(lp1.second() == var2);
    assert!(!lp1.first_link_on_start());
    assert!(lp1.second_link_on_start());
    assert!(lp1.first_unlinked_on_start());
    assert!(!lp1.second_unlinked_on_start());
    lp1.switch_svs();
    assert!(lp1.first() == var2);
    assert!(lp1.second() == var1);
    assert!(lp1.first_link_on_start());
    assert!(!lp1.second_link_on_start());

    let var3 = create_del(3, "1", 100, 200);
    let var4 = create_del(4, "1", 210, 400);
    let lp2 = Pair::from_svs(&var3, &var4, false, true);
    let var5 = create_del(4, "1", 250, 400);
    let lp3 = Pair::from_svs(&var3, &var5, false, true);
    let _ = lp3.same_variants(&lp2);
    let _ = lp3.has_link_clash(&lp2);
}

#[test]
fn test_deletion_bridges() {
    let mut t = LinxTester::new();
    let var1 = create_inv(t.next_var_id(), "1", 100, 1000, 1);
    let var2 = create_inf(t.next_var_id(), "1", 1001, -1);
    t.add(&[&var1, &var2]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let db = var1.db_link(false).unwrap();
    assert!(Some(db.clone()) == var2.db_link(true));
    assert_eq!(0, db.length());

    t.clear_clusters_and_svs();
    let var1 = create_inv(t.next_var_id(), "1", 100, 1000, 1);
    let var2 = create_inf(t.next_var_id(), "1", 1000, -1);
    t.add(&[&var1, &var2]);
    t.pre_clustering_init();
    t.cluster_and_analyse();
    let db = var1.db_link(false).unwrap();
    assert!(Some(db.clone()) == var2.db_link(true));
    assert_eq!(-1, db.length());
}

// ChainRoutinesTest

#[test]
fn test_chain_routines() {
    let var1 = create_del(1, "1", 1100, 1200);
    let var2 = create_del(2, "1", 1300, 1400);
    let lp1 = Pair::from_svs(&var1, &var2, false, true);
    lp1.set_is_assembled();
    let chain = Chain::create(0);
    chain.add_link(&lp1, true);
    assert!(chain.first_link_open_on_start());
    assert!(!chain.last_link_open_on_start());

    let var3 = create_del(3, "1", 1500, 1600);
    let lp2 = Pair::from_svs(&var3, &var2, true, false);
    lp2.set_is_assembled();
    assert!(!chain.can_add_linked_pair_to_start(&lp2));
    assert!(chain.can_add_linked_pair_to_end(&lp2));
    chain.add_link(&lp2, false);
    assert!(chain.first_sv() == var1);
    assert!(chain.last_sv() == var3);
    assert!(chain.sv_list()[1] == var2);

    let var4 = create_del(4, "1", 900, 1000);
    let lp3 = Pair::from_svs(&var1, &var4, true, false);
    assert!(chain.can_add_linked_pair_to_start(&lp3));
    assert!(!chain.can_add_linked_pair_to_end(&lp3));
    chain.add_link(&lp3, true);
    assert!(chain_is_valid(&chain));
    assert!(chain.first_sv() == var4);
    assert!(chain.last_sv() == var3);

    let var5 = create_dup(5, "1", 800, 1700);
    let lp4 = Pair::from_svs(&var5, &var4, true, true);
    assert!(chain.can_add_linked_pair_to_start(&lp4));
    chain.add_link(&lp4, true);
    assert_eq!(2, chain.assembly_link_count());
    assert!(chain_is_valid(&chain));

    let lp5 = Pair::from_svs(&var5, &var3, false, false);
    assert!(chain.can_add_linked_pair_to_end(&lp5));
    assert!(chain.link_would_close_chain(&lp5));
    chain.close_chain_with("CLOSE", 0);
    assert!(chain.is_closed_loop());
    assert_eq!(5, chain.link_count());

    let d = breakends_are_chained(&chain, &var4, false, &var3, true);
    assert_eq!(3, d[CHAIN_LINK_COUNT]);
    assert_eq!(2, d[CHAIN_ASSEMBLY_LINK_COUNT]);
    let d = breakends_are_chained(&chain, &var3, true, &var4, false);
    assert_eq!(3, d[CHAIN_LINK_COUNT]);
    assert_eq!(2, d[CHAIN_ASSEMBLY_LINK_COUNT]);
    let d = breakends_are_chained(&chain, &var1, false, &var2, true);
    assert_eq!(1, d[CHAIN_LINK_COUNT]);
    assert_eq!(1, d[CHAIN_ASSEMBLY_LINK_COUNT]);
    assert_eq!(0, breakends_are_chained(&chain, &var1, false, &var2, false)[CHAIN_LINK_COUNT]);
    assert_eq!(0, breakends_are_chained(&chain, &var1, true, &var2, false)[CHAIN_LINK_COUNT]);
    assert_eq!(0, breakends_are_chained(&chain, &var1, true, &var2, true)[CHAIN_LINK_COUNT]);
    assert_eq!(0, breakends_are_chained(&chain, &var5, false, &var1, true)[CHAIN_LINK_COUNT]);
}

fn two_link_chain(id: i32, a: &Sv, b: &Sv, c: &Sv) -> Chain {
    let chain = Chain::create(id);
    chain.add_link(&pair(a.be(false), b.be(true)), true);
    chain.add_link(&pair(b.be(false), c.be(true)), false);
    chain
}

#[test]
fn test_chain_reconciliation() {
    let v: Vec<Sv> = (1..=7).map(|i| create_del(i, "1", i * 100, i * 100 + 50)).collect();
    let chain1 = two_link_chain(1, &v[0], &v[1], &v[2]);
    let chain2 = two_link_chain(2, &v[2], &v[3], &v[4]);
    let mut list = vec![chain1.clone(), chain2.clone()];
    reconcile_chains(&mut list, false, 0, false);
    assert_eq!(1, list.len());
    assert_eq!(4, chain1.link_count());
    assert_eq!(5, chain1.sv_count());
    assert!(chain_is_valid(&chain1));

    let chain1 = two_link_chain(1, &v[0], &v[1], &v[2]);
    let chain3 = two_link_chain(3, &v[4], &v[5], &v[6]);
    let mut list = vec![chain3.clone(), chain1.clone(), chain2.clone()];
    reconcile_chains(&mut list, false, 0, false);
    assert_eq!(1, list.len());
    assert_eq!(6, chain3.link_count());
    assert_eq!(7, chain3.sv_count());
    assert!(chain_is_valid(&chain3));

    let chain1 = two_link_chain(1, &v[0], &v[1], &v[2]);
    let chain2 = Chain::create(2);
    chain2.add_link(&pair(v[1].be(false), v[3].be(true)), true);
    chain2.add_link(&pair(v[3].be(false), v[4].be(true)), false);
    let mut list = vec![chain1.clone(), chain2.clone()];
    reconcile_chains(&mut list, false, 0, false);
    assert_eq!(2, list.len());
    let mut list = vec![chain2, chain1];
    reconcile_chains(&mut list, false, 0, false);
    assert_eq!(2, list.len());
}

#[test]
fn test_chain_splitting_by_foldback() {
    let var1 = create_del(1, "1", 300, 400);
    let var2 = create_del(2, "1", 500, 600);
    let var3 = create_del(3, "1", 700, 800);
    let chain = two_link_chain(0, &var1, &var2, &var3);
    assert!(chain_is_valid(&chain));
    let var4 = create_inv(4, "1", 900, 1000, 1);
    foldback_chain_on_link(&chain, &pair(var4.be(true), var3.be(false)), &pair(var4.be(false), var3.be(false)));
    assert_eq!(6, chain.link_count());
    assert_eq!(4, chain.sv_count());
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(true) == chain.open_breakend(false));

    let chain = two_link_chain(0, &var1, &var2, &var3);
    assert!(chain_is_valid(&chain));
    let var4 = create_inv(4, "1", 100, 200, -1);
    foldback_chain_on_link(&chain, &pair(var1.be(true), var4.be(true)), &pair(var4.be(false), var1.be(true)));
    assert_eq!(6, chain.link_count());
    assert_eq!(4, chain.sv_count());
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(true) == chain.open_breakend(false));
}

#[test]
fn test_chain_splitting_by_foldback_chain() {
    let var1 = create_del(1, "1", 300, 400);
    let var2 = create_del(2, "1", 500, 600);
    let var3 = create_del(3, "1", 700, 800);
    let chain = two_link_chain(0, &var1, &var2, &var3);
    assert!(chain_is_valid(&chain));
    let var4 = create_bnd(4, "1", 900, 1, "2", 100, -1);
    let var5 = create_del(5, "2", 500, 600);
    let var6 = create_bnd(6, "1", 1000, -1, "2", 1000, 1);
    let fb = Chain::create(1);
    fb.add_link(&pair(var4.be(false), var5.be(true)), true);
    fb.add_link(&pair(var5.be(false), var6.be(false)), false);
    foldback_chain_on_chain(&chain, &fb, &pair(var4.be(true), var3.be(false)), &pair(var6.be(true), var3.be(false)));
    assert_eq!(8, chain.link_count());
    assert_eq!(6, chain.sv_count());
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(true) == chain.open_breakend(false));

    let chain = two_link_chain(0, &var1, &var2, &var3);
    assert!(chain_is_valid(&chain));
    let var4 = create_bnd(4, "1", 100, 1, "2", 100, -1);
    let var5 = create_del(5, "2", 500, 600);
    let var6 = create_bnd(6, "1", 200, -1, "2", 1000, 1);
    let fb = Chain::create(1);
    fb.add_link(&pair(var4.be(false), var5.be(true)), true);
    fb.add_link(&pair(var5.be(false), var6.be(false)), false);
    foldback_chain_on_chain(&chain, &fb, &pair(var1.be(true), var4.be(true)), &pair(var1.be(true), var6.be(true)));
    assert_eq!(8, chain.link_count());
    assert_eq!(6, chain.sv_count());
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(true) == chain.open_breakend(false));
}

#[test]
fn test_chain_splitting_by_complex_dup() {
    let var1 = create_del(1, "1", 300, 400);
    let var2 = create_del(2, "1", 500, 600);
    let var3 = create_del(3, "1", 700, 800);
    let chain = two_link_chain(0, &var1, &var2, &var3);
    assert!(chain_is_valid(&chain));
    let var4 = create_dup(4, "1", 200, 900);
    duplicate_chain_on_link(&chain, &pair(var4.be(true), var1.be(true)), &pair(var4.be(false), var3.be(false)));
    assert_eq!(6, chain.link_count());
    assert_eq!(4, chain.sv_count());
    assert!(chain_is_valid(&chain));
}

#[test]
fn test_chain_section_reversal() {
    let var1 = create_inv(1, "1", 100, 200, -1);
    let var2 = create_inv(2, "1", 300, 400, 1);
    let chain = Chain::create(0);
    chain.add_link(&pair(var1.be(false), var2.be(true)), true);
    reverse_section_on_breakend(&chain, &var1.be(false));
    assert!(chain.open_breakend(true) == Some(var1.be(false)));
    assert!(chain.open_breakend(false) == Some(var2.be(false)));
    reverse_section_on_breakend(&chain, &var2.be(true));
    assert!(chain.open_breakend(true) == Some(var1.be(false)));
    assert!(chain.open_breakend(false) == Some(var2.be(true)));

    let v: Vec<Sv> = (1..=5).map(|i| create_del(i, "1", i * 200 - 100, i * 200)).collect();
    let chain = Chain::create(0);
    chain.add_link(&pair(v[0].be(false), v[1].be(true)), true);
    for i in 1..4 {
        chain.add_link(&pair(v[i].be(false), v[i + 1].be(true)), false);
    }
    reverse_section_on_breakend(&chain, &v[2].be(false));
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(true) == Some(v[2].be(false)));
    reverse_section_on_breakend(&chain, &v[0].be(true));
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(true) == Some(v[0].be(true)));
    reverse_section_on_breakend(&chain, &v[1].be(true));
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(false) == Some(v[1].be(true)));
    reverse_section_on_breakend(&chain, &v[4].be(false));
    assert!(chain_is_valid(&chain));
    assert!(chain.open_breakend(false) == Some(v[4].be(false)));
}

// ChainingComplexTest

fn run(t: &mut LinxTester, vars: &[&Sv], ploidy_calcs: bool) {
    t.add(vars);
    t.pre_clustering_init_ploidy(ploidy_calcs);
    t.cluster_and_analyse();
}

#[test]
fn test_bfb_chain1() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 1000, 2000, -1, -1, INV, 8.0);
    let var2 = create_test_sv(2, "1", "1", 9000, 10000, 1, 1, INV, 4.0);
    let var3 = create_test_sv(3, "1", "1", 6000, 7000, 1, 1, INV, 2.0);
    let var4 = create_test_sv(4, "1", "1", 3000, 4000, 1, 1, INV, 1.0);
    let var5 = create_test_sv(5, "1", "2", 12000, 100, 1, 1, BND, 1.0);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5], false);
    for v in [&var1, &var2, &var3, &var4] {
        assert_eq!(v.id(), v.foldback_id(true));
    }
    let chain = only_chain(&first_cluster(&t));
    assert_eq!(5, chain.sv_count());
    assert_eq!(15, chain.link_count());
}

#[test]
fn test_bfb_chain2() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 2000, 3000, -1, -1, INV, 3.0);
    let var2 = create_test_sv(2, "1", "1", 9000, 10000, 1, 1, INV, 1.0);
    let var3 = create_test_sv(3, "1", "1", 5000, 6000, 1, 1, INV, 1.0);
    let var4 = create_test_sv(4, "1", "1", 1000, 8000, 1, 1, INV, 1.0);
    run(&mut t, &[&var1, &var2, &var3, &var4], false);
    for v in [&var1, &var2, &var3] {
        assert_eq!(v.id(), v.foldback_id(true));
    }
    assert_eq!(5, only_chain(&first_cluster(&t)).link_count());
}

#[test]
fn test_bfb_chain_with_chained_foldbacks() {
    let mut t = LinxTester::new();
    let a1 = create_test_sv(0, "1", "1", 2000, 5500, -1, -1, INV, 4.0);
    let a2 = create_test_sv(1, "1", "1", 3000, 5600, -1, 1, DUP, 4.0);
    set_assembled_link_info(&a1, false, &a2, false);
    let b = create_test_sv(2, "1", "1", 9000, 10000, 1, 1, INV, 2.0);
    let c = create_test_sv(3, "1", "1", 5000, 6000, 1, 1, INV, 1.0);
    let r = create_test_sv(4, "1", "1", 1000, 8000, 1, 1, INV, 1.0);
    run(&mut t, &[&a1, &a2, &b, &c, &r], false);
    assert!(b.is_foldback());
    assert!(c.is_foldback());
    assert!(a1.is_chained_foldback());
    assert!(a2.is_chained_foldback());
    assert_eq!(a2.id(), a1.foldback_id(true));
    assert_eq!(a1.id(), a2.foldback_id(true));
    assert_eq!(b.id(), b.foldback_id(true));
    assert_eq!(2 * 4 + 2 + 2 - 1, only_chain(&first_cluster(&t)).link_count());
}

#[test]
fn test_complex_dup_assembly_chain() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 2000, 3000, 1, -1, DEL, 2.0);
    let var2 = create_test_sv(2, "1", "2", 4000, 100, 1, -1, BND, 2.0);
    let var3 = create_test_sv(3, "1", "2", 5000, 200, -1, 1, BND, 2.0);
    let dup = create_test_sv(4, "1", "1", 1000, 6000, -1, 1, DUP, 1.0);
    set_assembled_link_info(&var1, false, &var2, true);
    set_assembled_link_info(&var2, false, &var3, false);
    run(&mut t, &[&var1, &var2, &var3, &dup], true);
    let c = first_cluster(&t);
    assert!(c.foldbacks().is_empty());
    assert_eq!(6, only_chain(&c).link_count());
}

#[test]
fn test_complex_multi_assembly_chain() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 2000, 1000, 1, -1, DEL, 5.0);
    let var2 = create_test_sv(2, "1", "1", 1500, 9000, 1, 1, INV, 5.0);
    let var3 = create_test_sv(3, "1", "1", 1200, 6000, 1, -1, DEL, 4.0);
    let var4 = create_test_sv(4, "1", "2", 1300, 200, -1, 1, BND, 2.0);
    let var5 = create_test_sv(5, "1", "2", 1400, 500, -1, -1, BND, 1.0);
    let var6 = create_test_sv(6, "1", "1", 100, 900, 1, -1, DEL, 2.0);
    let var7 = create_test_sv(7, "1", "1", 800, 2000, -1, -1, INV, 1.0);
    let var8 = create_test_sv(8, "1", "1", 1100, 10000, 1, 1, INV, 2.0);
    set_assembled_link_info(&var1, false, &var2, true);
    set_assembled_link_info(&var1, false, &var8, true);
    set_assembled_link_info(&var1, false, &var3, true);
    set_assembled_link_info(&var2, true, &var4, true);
    set_assembled_link_info(&var2, true, &var5, true);
    set_assembled_link_info(&var3, true, &var6, false);
    set_assembled_link_info(&var3, true, &var7, true);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6, &var7, &var8], true);
    let c = first_cluster(&t);
    assert_eq!(7, c.assembly_linked_pairs().len());
    let n: i32 = c.chains().iter().map(|x| x.assembly_link_count() as i32).sum();
    assert_eq!(7, n);
}

#[test]
fn test_bfb_chain_with_2_chained_foldbacks() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 1000, 2000, -1, -1, INV, 2.0);
    let var2 = create_test_sv(2, "1", "1", 3000, 4000, 1, 1, INV, 1.0);
    let var3 = create_test_sv(3, "1", "1", 6000, 7000, 1, -1, DEL, 4.0);
    let var4 = create_test_sv(4, "1", "1", 9000, 10000, 1, 1, INV, 2.0);
    let var5 = create_test_sv(5, "1", "1", 5000, 15000, -1, -1, INV, 1.0);
    let var6 = create_test_sv(6, "1", "2", 15500, 1000, 1, 1, BND, 1.0);
    set_assembled_link_info(&var1, true, &var2, true);
    set_assembled_link_info(&var1, true, &var2, false);
    set_assembled_link_info(&var3, false, &var4, true);
    set_assembled_link_info(&var3, false, &var4, false);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5, &var6], false);
    assert!(var1.is_foldback());
    assert!(var3.is_foldback());
    let chain = only_chain(&first_cluster(&t));
    assert_eq!(10, chain.link_count());
    assert_eq!(6, chain.sv_count());
    assert!((chain.jcn() - 1.0).abs() < 0.01);
}

#[test]
fn test_bfb_chain_to_complex_chain() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "2", 1000, 1000, -1, -1, BND, 1.0);
    let var2 = create_test_sv(2, "1", "2", 6000, 3000, -1, 1, BND, 1.0);
    let var3 = create_test_sv(3, "1", "2", 10000, 20000, 1, -1, BND, 2.0);
    let var4 = create_test_sv(4, "1", "2", 11000, 25000, -1, 1, BND, 2.0);
    let var5 = create_test_sv(5, "1", "3", 28000, 2000, 1, 1, BND, 1.0);
    set_assembled_link_info(&var1, false, &var2, false);
    set_assembled_link_info(&var3, false, &var4, false);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5], true);
    assert!(var1.is_foldback());
    assert!(var2.is_foldback());
    let chain = only_chain(&first_cluster(&t));
    assert_eq!(6, chain.link_count());
    assert_eq!(5, chain.sv_count());
    assert!((chain.jcn() - 1.0).abs() < 0.01);
}

#[test]
fn test_complex_dup_simple() {
    let mut t = LinxTester::new();
    let var1 = create_test_sv(1, "1", "1", 4000, 5000, 1, -1, DEL, 2.0);
    let var2 = create_test_sv(2, "1", "1", 2000, 6000, -1, 1, DUP, 1.0);
    let var3 = create_test_sv(3, "1", "1", 7000, 8000, 1, -1, DEL, 1.0);
    let var4 = create_test_sv(4, "1", "2", 500, 100, 1, -1, BND, 1.0);
    let var5 = create_test_sv(5, "1", "2", 1000, 200, -1, 1, BND, 1.0);
    run(&mut t, &[&var1, &var2, &var3, &var4, &var5], false);
    let c = first_cluster(&t);
    assert!(c.foldbacks().is_empty());
    let chain = only_chain(&c);
    assert_eq!(5, chain.link_count());
    assert_eq!(5, chain.sv_count());
}

#[test]
fn test_complex_dup_chain() {
    let mut t = LinxTester::new();
    let a = create_test_sv(1, "1", "1", 1000, 5000, 1, 1, INV, 1.0);
    let b = create_test_sv(2, "1", "1", 4000, 9000, -1, 1, DUP, 2.0);
    let c = create_test_sv(3, "1", "1", 7000, 8000, 1, -1, DEL, 2.0);
    let d = create_test_sv(4, "1", "1", 3000, 6000, 1, -1, DEL, 2.0);
    let e = create_test_sv(5, "1", "1", 2000, 10000, -1, -1, INV, 1.0);
    let dup = create_test_sv(6, "1", "1", 2500, 4500, -1, 1, DUP, 1.0);
    run(&mut t, &[&a, &b, &c, &d, &e, &dup], false);
    let cl = first_cluster(&t);
    assert!(cl.foldbacks().is_empty());
    assert_eq!(8, only_chain(&cl).link_count());
}
