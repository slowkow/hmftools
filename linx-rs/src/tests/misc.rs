//! misc: MiscTest and JcnCalcsTest, with tests of the corrected Java faults (LX2 to LX5)

use super::*;
use crate::cn::{calc_poisson_observed_given_prob, CnDataLoader, READ_COUNT_PROB_HIGH, READ_COUNT_PROB_LOW};
use crate::utils::{calc_consistency_sv, find_centromere_breakend_index, get_chromosomal_arm_length};

#[test]
fn test_consistency() {
    assert_eq!(0, calc_consistency_sv(&create_del(1, "1", 100, 200)));
    assert_eq!(0, calc_consistency_sv(&create_ins(1, "1", 100, 200)));
    assert_eq!(0, calc_consistency_sv(&create_dup(1, "1", 100, 200)));
    assert_eq!(2, calc_consistency_sv(&create_inv(1, "1", 100, 200, 1)));
    assert_eq!(0, calc_consistency_sv(&create_bnd(1, "1", 100, 1, "2", 100, -1)));
    assert_eq!(1, calc_consistency_sv(&create_sgl(1, "1", 100, 1)));
}

#[test]
fn test_breakend_lists() {
    let mut t = LinxTester::new();
    let centromere = get_chromosomal_arm_length("1", Arm::P);
    let q = centromere + 10000000;
    let var1 = create_del(t.next_var_id(), "1", 100, 200);
    let var2 = create_del(t.next_var_id(), "1", 300, 400);
    let var3 = create_del(t.next_var_id(), "1", q + 1000, q + 2000);
    let var4 = create_del(t.next_var_id(), "1", q + 10000, q + 20000);
    for (vars, p, qi) in [(vec![&var1, &var2, &var3, &var4], 3, 4), (vec![&var1, &var2], 3, -1), (vec![&var3, &var4], -1, 0)] {
        t.clear_clusters_and_svs();
        t.add(&vars);
        t.pre_clustering_init();
        let list = t.analyser.state.chr_breakend_map.get(&"1".to_string()).cloned();
        assert_eq!(p, find_centromere_breakend_index(list.as_ref(), Arm::P));
        assert_eq!(qi, find_centromere_breakend_index(list.as_ref(), Arm::Q));
    }
}

#[test]
fn test_expected_read_count() {
    const LOW: f64 = READ_COUNT_PROB_LOW;
    const HIGH: f64 = READ_COUNT_PROB_HIGH;
    let first: Vec<(i32, i32)> = (0..=2000).map(|i| (calc_poisson_observed_given_prob(i, LOW), calc_poisson_observed_given_prob(i, HIGH))).collect();
    for _ in 0..=10 {
        for i in 0..=2000 {
            assert_eq!(first[i as usize], (calc_poisson_observed_given_prob(i, LOW), calc_poisson_observed_given_prob(i, HIGH)));
        }
    }
}

// the corrected Java faults

#[test]
fn test_lx2_expect_single_chromosome() {
    assert!(CnDataLoader::expect_single_chromosome(true, "X"));
    assert!(CnDataLoader::expect_single_chromosome(true, "chrY"));
    assert!(!CnDataLoader::expect_single_chromosome(false, "chrX"));
    // Java: GRCh38's names never match
    crate::compat::with_java(true, || assert!(!CnDataLoader::expect_single_chromosome(true, "chrX")));
    crate::compat::with_java(true, || assert!(CnDataLoader::expect_single_chromosome(true, "X")));
}

#[test]
fn test_lx5_has_short_arm() {
    use crate::drivers::has_short_arm;
    assert!(has_short_arm("13"));
    assert!(has_short_arm("chr21"));
    assert!(!has_short_arm("chr1"));
    crate::compat::with_java(true, || assert!(!has_short_arm("chr13")));
    crate::compat::with_java(true, || assert!(has_short_arm("22")));
}
