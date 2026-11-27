//! ClusterMetrics and ClusterAnnotations' cluster features (reportClusterFeatures): the chain arm
//! consistency (CHROMO_THRIP), the deletion metrics, replication before repair (REP_REPAIR) and ChainMetrics.
//! Of these only TotalDeleted reaches a sample file (germline disruption reportability); the rest feeds
//! LINX's cohort files and unit tests.

use std::collections::HashMap;

use hmf_common::java_hash::JavaHashMap;
use hmf_common::sv::StructuralVariantType::*;

use super::{find_end_index, find_start_index};
use crate::types::*;
use crate::utils::{calc_consistency_be, copy_numbers_equal};

/// ClusterMetrics
#[derive(Clone, Debug)]
pub struct ClusterMetrics {
    pub origin_arms: i32,
    pub fragment_arms: i32,
    pub consistent_arms: i32,
    pub complex_arms: i32,
    pub db_count: i32,
    pub short_db_count: i32,
    pub total_db_length: i64,
    pub traversed_del_length: i32,
    pub traversed_del_count: i32,
    pub traversed_short_del_count: i32,
    pub implied_ti_count: i32,
    pub chained_length: i64,
    pub total_range: i64,
    pub traversed_range: i64,
    pub total_deleted: i64,
    pub valid_allele_jcn_segment_perc: f64,
}

impl Default for ClusterMetrics {
    fn default() -> Self {
        ClusterMetrics {
            origin_arms: 0,
            fragment_arms: 0,
            consistent_arms: 0,
            complex_arms: 0,
            db_count: 0,
            short_db_count: 0,
            total_db_length: 0,
            traversed_del_length: 0,
            traversed_del_count: 0,
            traversed_short_del_count: 0,
            implied_ti_count: 0,
            chained_length: 0,
            total_range: 0,
            traversed_range: 0,
            total_deleted: 0,
            valid_allele_jcn_segment_perc: 1.0,
        }
    }
}

impl ClusterMetrics {
    /// populate: deletion bridges, simple DELs and traversed DELs of other clusters
    pub fn populate(&mut self, cluster: &SvCluster, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
        if cluster.sv_count() == 0 {
            return;
        }
        let mut cross_centromere: Vec<String> = Vec::new();
        for p in cluster.linked_pairs() {
            if p.first_breakend().arm() != p.second_breakend().arm() {
                let c = p.chromosome();
                if !cross_centromere.contains(&c) {
                    cross_centromere.push(c);
                }
            }
        }
        for (chromosome, list) in cluster.chr_breakend_entries() {
            let chr_list = chr_breakend_map.get(&chromosome);
            let start = find_start_index(&list);
            let end = find_end_index(&list);
            if start < end {
                self.total_range += (list[end as usize].position() - list[start as usize].position()) as i64;
            }
            let mut i = start;
            while i < end {
                let b = &list[i as usize];
                let var = b.sv().clone();
                let next = &list[i as usize + 1];
                let distance = next.position() - b.position();
                let is_db = b.db_link().is_some() && b.db_link() == next.db_link();
                let is_simple_del = !is_db && var.sv_type() == DEL && b.forward_orient() && next.sv() == &var;
                let is_traversed_db = !is_db && !is_simple_del && b.forward_orient() && next.reverse_orient();
                if is_db || is_simple_del || is_traversed_db {
                    self.db_count += 1;
                    let del_length = if is_db { b.db_link().unwrap().length().max(0) } else if is_simple_del { var.length() } else { distance };
                    self.total_db_length += del_length as i64;
                    if del_length <= SHORT_DB_LENGTH {
                        self.short_db_count += 1;
                    }
                } else if next.chr_pos_index() > b.chr_pos_index() + 2 && b.reverse_orient() && next.forward_orient() {
                    let chr_list = chr_list.expect("chromosome breakends");
                    let mut implied = 0;
                    let mut j = b.chr_pos_index() + 1;
                    while j < next.chr_pos_index() - 1 {
                        let nc = &chr_list[j as usize];
                        let nc_var = nc.sv().clone();
                        if nc_var.sv_type() == DEL && chr_list[j as usize + 1].sv() == &nc_var {
                            self.traversed_del_length += nc_var.length();
                            self.traversed_del_count += 1;
                            if nc_var.length() <= SHORT_DB_LENGTH {
                                self.traversed_short_del_count += 1;
                                implied += 1;
                            } else {
                                implied += 2;
                            }
                            j += 1;
                        }
                        j += 1;
                    }
                    self.implied_ti_count += implied;
                }
                if b.arm() == Arm::P && next.arm() == Arm::Q && b.forward_orient() && next.reverse_orient() && !cross_centromere.contains(&chromosome) {
                    self.total_range -= distance as i64;
                }
                i += 1;
            }
        }
    }
}

/// ChainMetrics
#[derive(Default, Debug)]
pub struct ChainMetrics {
    pub internal_tis: i32,
    pub internal_ti_cn_gain: i32,
    pub internal_short_tis: i32,
    pub external_tis: i32,
    pub external_ti_cn_gain: i32,
    pub external_short_tis: i32,
    pub overlapping_tis: i32,
    pub chain_ends_face: i32,
    pub chain_ends_away: i32,
}

pub fn extract_chain_metrics(chain: &SvChain) -> ChainMetrics {
    let mut m = ChainMetrics::default();
    let (Some(cs), Some(ce)) = (chain.open_breakend(true), chain.open_breakend(false)) else { return m };
    let lower = if cs.position() < ce.position() { cs.clone() } else { ce.clone() };
    let upper = if cs == lower { ce.clone() } else { cs.clone() };
    if ce.chr_arm() == cs.chr_arm() {
        if lower.forward_orient() && upper.reverse_orient() {
            m.chain_ends_away += 1;
        } else if lower.reverse_orient() && upper.forward_orient() {
            m.chain_ends_face += 1;
        }
    }
    for p in chain.linked_pairs() {
        if p.first().sv_type() == SGL || p.second().sv_type() == SGL {
            continue;
        }
        let lt = p.location_type();
        if lt == LOCATION_TYPE_INTERNAL {
            m.internal_tis += 1;
            if p.base_length() <= SHORT_TI_LENGTH {
                m.internal_short_tis += 1;
            }
            if p.has_copy_number_gain() {
                m.internal_ti_cn_gain += 1;
            }
        } else if lt == LOCATION_TYPE_REMOTE || lt == LOCATION_TYPE_EXTERNAL {
            m.external_tis += 1;
            if p.base_length() <= SHORT_TI_LENGTH {
                m.external_short_tis += 1;
            }
            if p.has_copy_number_gain() {
                m.external_ti_cn_gain += 1;
            }
        }
        if p.overlap_count() > 0 {
            m.overlapping_tis += 1;
        }
    }
    m
}

const INCONSISTENT_ARM: i32 = -2;
const CONSISTENT_ARM: i32 = 0;

/// annotateClusterChains: arm consistency of chain ends and unlinked SVs, and the chromothripsis annotation
pub fn annotate_cluster_chains(cluster: &Cluster) {
    if cluster.sv_count() == 1 || (cluster.arm_count() == 1 && cluster.chains().is_empty()) {
        return;
    }
    let mut fragment_arms: Vec<String> = Vec::new();
    // a HashMap in Java, read for its size and values only
    let mut origin_arms: HashMap<String, i32> = HashMap::new();
    let mut inconsistent_chains = 0;
    let mut long_ti_count = 0;
    for chain in cluster.chains() {
        let consistent = chain.is_consistent();
        cluster.metrics.borrow_mut().chained_length += chain.length(false);
        for p in chain.linked_pairs() {
            if p.first().is_sgl_breakend() || p.second().is_sgl_breakend() {
                continue;
            }
            let chr_arm = p.first().be(p.first_link_on_start()).chr_arm().to_string();
            if p.base_length() <= SHORT_TI_LENGTH {
                if !fragment_arms.contains(&chr_arm) {
                    fragment_arms.push(chr_arm);
                }
            } else {
                long_ti_count += 1;
            }
        }
        if !consistent {
            inconsistent_chains += 1;
        }
        for s in [SE_START, SE_END] {
            if let Some(b) = chain.open_breakend(is_start(s)) {
                add_breakend_to_arm_consistency(&mut origin_arms, &b);
            }
        }
    }
    for var in cluster.unlinked_svs() {
        for s in [SE_START, SE_END] {
            if s == SE_END && var.is_sgl_breakend() {
                continue;
            }
            add_breakend_to_arm_consistency(&mut origin_arms, &var.be(is_start(s)));
        }
    }
    fragment_arms.retain(|x| !origin_arms.contains_key(x));
    let consistent_arms = origin_arms.values().filter(|x| **x >= -1 && **x <= 1).count() as i32;
    let mut arm_cluster_counts: HashMap<String, i32> = HashMap::new();
    for ac in cluster.arm_clusters.borrow().iter() {
        if ac.arm_cluster_type == ArmClusterType::TI_ONLY {
            continue;
        }
        *arm_cluster_counts.entry(ac.chr_arm()).or_insert(0) += 1;
    }
    let complex_arms = arm_cluster_counts.values().filter(|x| **x > 1).count() as i32;
    {
        let mut m = cluster.metrics.borrow_mut();
        m.origin_arms = origin_arms.len() as i32;
        m.fragment_arms = fragment_arms.len() as i32;
        m.consistent_arms = consistent_arms;
        m.complex_arms = complex_arms;
    }
    if cluster.sv_count() > 2 {
        let is_complex = cluster.requires_replication() || !cluster.foldbacks().is_empty();
        let is_complete = inconsistent_chains == 0 && consistent_arms == origin_arms.len() as i32;
        if is_complete && !is_complex && long_ti_count > 0 && cluster.resolved_type() == ResolvedType::COMPLEX {
            cluster.add_annotation(CLUSTER_ANNOT_SHATTERING);
        }
    }
}

fn add_breakend_to_arm_consistency(map: &mut HashMap<String, i32>, breakend: &Be) {
    let chr_arm = breakend.chr_arm().to_string();
    let Some(&c) = map.get(&chr_arm) else {
        map.insert(chr_arm, calc_consistency_be(breakend));
        return;
    };
    if c == INCONSISTENT_ARM {
        return;
    }
    if c == CONSISTENT_ARM || calc_consistency_be(breakend) == c {
        map.insert(chr_arm, INCONSISTENT_ARM);
    } else {
        map.insert(chr_arm, CONSISTENT_ARM);
    }
}

/// annotateClusterDeletions
pub fn annotate_cluster_deletions(cluster: &Cluster, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
    if cluster.sv_count() == 1 || cluster.resolved_type() == ResolvedType::LINE {
        return;
    }
    let mut m = cluster.metrics.borrow().clone();
    m.populate(cluster, chr_breakend_map);
    *cluster.metrics.borrow_mut() = m;
}

/// annotateReplicationBeforeRepair
pub fn annotate_replication_before_repair(cluster: &Cluster) {
    if cluster.sv_count() == 1 || cluster.sv_count() > 100 {
        return;
    }
    if cluster.resolved_type() != ResolvedType::COMPLEX && cluster.resolved_type() != ResolvedType::SIMPLE_GRP {
        return;
    }
    for chain in cluster.chains() {
        if !chain.is_consistent() || chain.is_closed_loop() {
            continue;
        }
        let (Some(cs), Some(ce)) = (chain.open_breakend(true), chain.open_breakend(false)) else { continue };
        if ce.chr_arm() != cs.chr_arm() {
            continue;
        }
        if chain.has_repeated_sv() {
            continue;
        }
        let cm = extract_chain_metrics(&chain);
        if cm.chain_ends_away != 1 || cm.internal_ti_cn_gain == 0 || cm.overlapping_tis == 0 {
            continue;
        }
        let gain_length = copy_number_gain_length(cluster, &chain);
        if gain_length == 0 {
            continue;
        }
        let chain_range = (ce.position() - cs.position()).abs() as i64;
        let gain_percent = gain_length as f64 / chain_range as f64;
        if gain_percent >= 0.05 && chain_range >= 1000 {
            cluster.add_annotation(CLUSTER_ANNOT_REP_REPAIR);
        }
    }
}

fn copy_number_gain_length(cluster: &Cluster, chain: &SvChain) -> i64 {
    let cs = chain.open_breakend(true).unwrap();
    let ce = chain.open_breakend(false).unwrap();
    let lower = cs.cluster_chr_pos_index().min(ce.cluster_chr_pos_index());
    let upper = cs.cluster_chr_pos_index().max(ce.cluster_chr_pos_index());
    let ends_cn = cs.copy_number().max(ce.copy_number());
    let mut gain: i64 = 0;
    let map = cluster.chr_breakend_map();
    let list = map.get(&cs.chromosome().to_string()).expect("cluster chromosome breakends");
    let svs = chain.sv_list();
    let mut prev_pos = 0;
    let mut prev_cn = 0.0;
    let mut in_segment = false;
    let mut net = 0.0;
    let mut i = lower + 1;
    while i < upper - 1 {
        let b = &list[i as usize];
        i += 1;
        if !contains(&svs, b.sv()) {
            return 0;
        }
        if !in_segment {
            if b.forward_orient() {
                if prev_pos > 0 && b.copy_number() > ends_cn && !copy_numbers_equal(b.copy_number(), ends_cn) {
                    gain += (b.position() - prev_pos) as i64;
                }
                continue;
            }
            in_segment = true;
            prev_pos = b.position();
            prev_cn = b.copy_number();
            net = b.jcn();
        } else {
            let seg = b.position() - prev_pos;
            if b.reverse_orient() {
                if prev_cn > ends_cn && !copy_numbers_equal(prev_cn, ends_cn) {
                    gain += seg as i64;
                }
                prev_pos = b.position();
                prev_cn = b.copy_number();
                net += b.jcn();
            } else {
                if b.copy_number() > ends_cn && !copy_numbers_equal(b.copy_number(), ends_cn) {
                    gain += seg as i64;
                }
                net -= b.jcn();
                prev_pos = b.position();
                if copy_numbers_equal(net, 0.0) {
                    in_segment = false;
                }
            }
        }
    }
    gain
}

/// reportClusterFeatures (DoubleMinuteFinder.reportCluster only logs)
pub fn report_cluster_features(cluster: &Cluster, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
    annotate_cluster_chains(cluster);
    annotate_cluster_deletions(cluster, chr_breakend_map);
    annotate_replication_before_repair(cluster);
    let mut m = cluster.metrics.borrow_mut();
    if m.total_deleted == 0 {
        m.total_deleted = m.total_db_length;
    }
}
