//! ClusterClassification and PairResolution: the resolved type of a cluster, and the classification of
//! two-break clusters (reciprocal translocations and inversions, TI types, synthetic DELs and DUPs).

use hmf_common::java_hash::JavaHashMap;
use hmf_common::sv::StructuralVariantType as SvType;

use crate::chaining::jcn_limits::{jcn_match, jcn_match_svs};
use crate::chaining::utils::{get_min_templated_insertion_length, reverse_section_on_breakend};
use crate::types::ResolvedType::*;
use crate::types::*;
use crate::utils::calc_consistency_be;

pub const NO_LENGTH: i32 = -1;

// LinxCommonTypes super types
pub const SUPER_TYPE_SIMPLE: &str = "SIMPLE";
pub const SUPER_TYPE_INSERTION: &str = "INSERTION";
pub const SUPER_TYPE_ARTIFACT: &str = "ARTIFACT";
pub const SUPER_TYPE_RECIPROCAL: &str = "RECIPROCAL";
pub const SUPER_TYPE_TEMPLATED_INSERTION: &str = "TEMPLATED_INSERTION";
pub const SUPER_TYPE_DOUBLE_MINUTE: &str = "DOUBLE_MINUTE";
pub const SUPER_TYPE_INCOMPLETE: &str = "INCOMPLETE";
pub const SUPER_TYPE_COMPLEX: &str = "COMPLEX";

pub fn get_cluster_category(cluster: &SvCluster) -> &'static str {
    let rt = cluster.resolved_type();
    if rt == LINE {
        return SUPER_TYPE_INSERTION;
    }
    if rt.is_simple() || rt == UNBAL_TRANS {
        return SUPER_TYPE_SIMPLE;
    }
    if is_filtered_resolved_type(rt) {
        return SUPER_TYPE_ARTIFACT;
    }
    if is_cluster_reciprocal_type(rt) {
        return SUPER_TYPE_RECIPROCAL;
    }
    if is_cluster_templated_insertion_type(rt) {
        return SUPER_TYPE_TEMPLATED_INSERTION;
    }
    if rt == DOUBLE_MINUTE {
        return SUPER_TYPE_DOUBLE_MINUTE;
    }
    if is_incomplete_type(rt) || rt == PAIR_OTHER {
        return SUPER_TYPE_INCOMPLETE;
    }
    SUPER_TYPE_COMPLEX
}

pub fn is_simple_single_sv(cluster: &SvCluster) -> bool {
    cluster.sv_count() == 1 && cluster.sv0().is_simple_type()
}

pub fn is_incomplete_type(rt: ResolvedType) -> bool {
    rt == ResolvedType::INV || rt == ResolvedType::SGL || rt == ResolvedType::INF
}

pub fn is_synthetic_type(cluster: &SvCluster) -> bool {
    let rt = cluster.resolved_type();
    if rt.is_simple() {
        return if cluster.sgl_breakend_count() == 2 { false } else { cluster.sv_count() > 1 };
    }
    if is_cluster_pair_type(rt) {
        return cluster.sv_count() > 2;
    }
    if rt == UNBAL_TRANS {
        return cluster.sv_count() > 1;
    }
    if is_incomplete_type(rt) {
        return cluster.sv_count() > 1;
    }
    false
}

pub fn is_filtered_resolved_type(rt: ResolvedType) -> bool {
    rt == DUP_BE || rt == LOW_VAF || rt == PAIR_INF || rt == SGL_MAPPED_INF
}

pub fn set_cluster_resolved_state(
    cluster: &Cluster, is_final: bool, is_germline: bool, long_del_threshold: i32, long_dup_threshold: i32,
    chr_breakend_map: &JavaHashMap<String, Vec<Be>>,
) {
    if cluster.resolved_type() != NONE {
        return;
    }
    if cluster.has_linking_line_elements() {
        cluster.set_resolved(true, LINE);
        return;
    }
    if cluster.sv_count() <= 2 && !cluster.svs().iter().any(|x| !x.is_simple_type()) {
        if cluster.chain_count() == 1 {
            classify_pair_clusters(cluster, long_del_threshold as i64, long_dup_threshold as i64, chr_breakend_map, is_germline);
            if cluster.resolved_type() != NONE {
                return;
            }
        }
        let has_long_svs = cluster.svs().iter().any(|v| {
            (v.sv_type() == SvType::DEL && v.length() >= long_del_threshold) || (v.sv_type() == SvType::DUP && v.length() >= long_dup_threshold)
        });
        if cluster.sv_count() == 1 {
            match cluster.sv0().sv_type() {
                SvType::DEL => cluster.set_resolved(!has_long_svs, ResolvedType::DEL),
                SvType::DUP => cluster.set_resolved(!has_long_svs, ResolvedType::DUP),
                SvType::INS => cluster.set_resolved(!has_long_svs, ResolvedType::INS),
                _ => {}
            }
        } else {
            cluster.set_resolved(!has_long_svs, SIMPLE_GRP);
        }
        return;
    }
    if cluster.sgl_breakend_count() == 0 {
        classify_pair_clusters(cluster, long_del_threshold as i64, long_dup_threshold as i64, chr_breakend_map, is_germline);
    } else if cluster.sv_count() == 2 && cluster.is_consistent() && cluster.sgl_breakend_count() == 2 {
        let rt = classify_single_pair_resolved_type(&cluster.sv0(), &cluster.sv(1).unwrap());
        if rt != NONE {
            cluster.set_resolved(true, rt);
        }
    }
    if cluster.resolved_type() != NONE {
        return;
    }
    if is_final {
        if cluster.sv_count() == 1 {
            match cluster.sv0().sv_type() {
                SvType::BND => cluster.set_resolved(false, UNBAL_TRANS),
                SvType::INV => cluster.set_resolved(false, ResolvedType::INV),
                SvType::INF => cluster.set_resolved(false, ResolvedType::INF),
                SvType::SGL => cluster.set_resolved(false, ResolvedType::SGL),
                _ => {}
            }
            return;
        }
        mark_synthetic_incompletes(cluster);
        if cluster.resolved_type() != NONE {
            return;
        }
        if cluster.sv_count() == 2 {
            cluster.set_resolved(false, PAIR_OTHER);
            return;
        }
        cluster.set_resolved(false, COMPLEX);
    }
}

pub fn get_synthetic_length(cluster: &SvCluster) -> i32 {
    if cluster.resolved_type() == SGL_PAIR_DEL || cluster.resolved_type() == SGL_PAIR_DUP {
        (cluster.sv0().position(true) - cluster.sv(1).unwrap().position(true)).abs()
    } else {
        if cluster.chain_count() != 1 {
            return NO_LENGTH;
        }
        let chain = &cluster.chains()[0];
        (chain.open_breakend(true).unwrap().position() - chain.open_breakend(false).unwrap().position()).abs()
    }
}

pub fn get_synthetic_ti_length(cluster: &SvCluster) -> i32 {
    if cluster.chain_count() != 1 {
        return NO_LENGTH;
    }
    // Java's getAsInt throws on an empty chain
    cluster.chains()[0].linked_pairs().iter().map(|x| x.base_length()).max().expect("chain links")
}

pub fn get_synthetic_gap_length(cluster: &SvCluster) -> i32 {
    if cluster.chain_count() != 1 {
        return NO_LENGTH;
    }
    let chain = &cluster.chains()[0];
    let cs = chain.open_breakend(true).unwrap();
    let ce = chain.open_breakend(false).unwrap();
    if cs.chr_arm() != ce.chr_arm() {
        return NO_LENGTH;
    }
    let mut min_distance = NO_LENGTH;
    for pair in chain.linked_pairs() {
        if pair.base_length() > SHORT_TI_LENGTH {
            if pair.chromosome() != cs.chromosome() {
                return NO_LENGTH;
            }
        } else if pair.chromosome() != cs.chromosome() {
            continue;
        }
        let (f, s) = (pair.first_breakend().position(), pair.second_breakend().position());
        let ds = (cs.position() - f).abs().min((cs.position() - s).abs());
        let de = (ce.position() - f).abs().min((ce.position() - s).abs());
        let gap = ds.min(de);
        if min_distance == NO_LENGTH || gap < min_distance {
            min_distance = gap;
        }
    }
    min_distance
}

pub fn mark_synthetic_incompletes(cluster: &Cluster) {
    if cluster.sv_count() > 5 || !cluster.is_fully_chained(false) || cluster.chain_count() != 1 {
        return;
    }
    let chain = cluster.chains()[0].clone();
    let mut total = 0;
    let mut longest = 0;
    for pair in chain.linked_pairs() {
        if pair.base_length() > SHORT_TI_LENGTH {
            return;
        }
        longest = longest.max(pair.base_length());
        total += pair.base_length();
    }
    let rt = if chain.last_sv().sv_type() == SvType::INF || chain.first_sv().sv_type() == SvType::INF {
        ResolvedType::INF
    } else if chain.last_sv().sv_type() == SvType::SGL || chain.first_sv().sv_type() == SvType::SGL {
        ResolvedType::SGL
    } else {
        let s = chain.open_breakend(true).unwrap();
        let e = chain.open_breakend(false).unwrap();
        if s.chromosome() != e.chromosome() || s.arm() != e.arm() { UNBAL_TRANS } else { ResolvedType::INV }
    };
    hmf_common::debug!(
        "cluster({}) chain(links=({} len={} tiLen({}) marked as {}",
        cluster.id(),
        chain.link_count(),
        total,
        longest,
        rt
    );
    cluster.set_resolved(false, rt);
}

pub fn classify_single_pair_resolved_type(sgl1: &Sv, sgl2: &Sv) -> ResolvedType {
    if sgl1.sgl_to_centromere_or_telomere() || sgl2.sgl_to_centromere_or_telomere() {
        return NONE;
    }
    if sgl1.is_inferred_sgl() && sgl2.is_inferred_sgl() {
        return NONE;
    }
    let b1 = sgl1.be(true);
    let b2 = sgl2.be(true);
    if b1.orientation() == b2.orientation() {
        return NONE;
    }
    if !jcn_match_svs(sgl1, sgl2) {
        return NONE;
    }
    let faces = (b1.position() < b2.position() && b1.reverse_orient()) || (b2.position() < b1.position() && b2.reverse_orient());
    let length = (b1.position() - b2.position()).abs();
    if faces {
        let min_ti = get_min_templated_insertion_length(&b1, &b2);
        if length + 1 >= min_ti { SGL_PAIR_DUP } else { SGL_PAIR_INS }
    } else if length >= MIN_DEL_LENGTH {
        SGL_PAIR_DEL
    } else {
        SGL_PAIR_INS
    }
}

// PairResolution

pub fn is_cluster_pair_type(rt: ResolvedType) -> bool {
    is_cluster_reciprocal_type(rt) || is_cluster_templated_insertion_type(rt)
}

pub fn is_cluster_reciprocal_type(rt: ResolvedType) -> bool {
    matches!(rt, RECIP_TRANS | RECIP_TRANS_DEL_DUP | RECIP_TRANS_DUPS | RECIP_INV | RECIP_INV_DEL_DUP | RECIP_INV_DUPS | RESOLVED_FOLDBACK | FB_INV_PAIR)
}

pub fn is_cluster_templated_insertion_type(rt: ResolvedType) -> bool {
    matches!(rt, DEL_TI | DUP_TI | UNBAL_TRANS_TI)
}

pub fn classify_pair_clusters(
    cluster: &Cluster, long_del_threshold: i64, long_dup_threshold: i64, chr_breakend_map: &JavaHashMap<String, Vec<Be>>, is_germline: bool,
) {
    if cluster.chain_count() > 2 || cluster.sgl_breakend_count() > 0 {
        return;
    }
    if cluster.has_annotation(CLUSTER_ANNOT_DM) {
        return;
    }
    let mut long_ti_link: Option<Pair> = None;
    let mut cluster_chains: Vec<Chain> = cluster.chains();
    let unlinked = cluster.unlinked_svs();
    let mut unchained_sv: Option<Sv> = unlinked.first().cloned();
    let mut unchained_sv_count = unlinked.len();
    let mut long_ti_links: Vec<Pair> = Vec::new();
    let mut existing_chain_modified = false;
    if !cluster_chains.is_empty() {
        let chains = cluster.chains();
        long_ti_links.extend(chains[0].linked_pairs().into_iter().filter(|x| x.base_length() > SHORT_TI_LENGTH));
        if chains.len() == 2 {
            long_ti_links.extend(chains[1].linked_pairs().into_iter().filter(|x| x.base_length() > SHORT_TI_LENGTH));
        }
        if long_ti_links.len() > 1 {
            return;
        }
        if long_ti_links.len() == 1 {
            long_ti_link = Some(long_ti_links[0].clone());
        }
    }
    let is_single_chain = cluster.is_fully_chained(false) && cluster_chains.len() == 1;
    if is_single_chain && long_ti_link.is_none() {
        classify_synthetic_del_dups(cluster, long_del_threshold, long_dup_threshold);
        return;
    }
    if !is_germline && long_ti_links.iter().any(|x| is_interrupted_ti(x, chr_breakend_map)) {
        return;
    }
    if is_single_chain {
        if let Some(lt) = &long_ti_link {
            let chain = cluster_chains[0].clone();
            cluster_chains.clear();
            existing_chain_modified = true;
            let pairs = chain.linked_pairs();
            if pairs.len() > 1 {
                let mut new_chain = Chain::create(chain.id());
                if &pairs[0] == lt || &pairs[pairs.len() - 1] == lt {
                    unchained_sv_count = 1;
                    unchained_sv = Some(if &pairs[0] == lt { pairs[0].first() } else { pairs[pairs.len() - 1].second() });
                    for pair in &pairs {
                        if pair != lt {
                            new_chain.add_link(pair, false);
                        }
                    }
                    cluster_chains.push(new_chain);
                } else {
                    for pair in &pairs {
                        if pair == lt {
                            cluster_chains.push(new_chain);
                            new_chain = Chain::create(chain.id() + 1);
                            continue;
                        }
                        new_chain.add_link(pair, false);
                    }
                    cluster_chains.push(new_chain);
                }
            }
        }
    }
    let (start_be1, end_be1, start_be2, end_be2, uniform_jcn);
    if cluster_chains.is_empty() {
        if cluster.sv_count() != 2 {
            return;
        }
        let var1 = cluster.sv0();
        let var2 = cluster.sv(1).unwrap();
        start_be1 = var1.breakend(true);
        end_be1 = var1.breakend(false);
        start_be2 = var2.breakend(true);
        end_be2 = var2.breakend(false);
        uniform_jcn = jcn_match_svs(&var1, &var2);
    } else if cluster.is_fully_chained(false) && cluster_chains.len() == 2 {
        let c1 = &cluster_chains[0];
        let c2 = &cluster_chains[1];
        start_be1 = c1.open_breakend(true);
        end_be1 = c1.open_breakend(false);
        start_be2 = c2.open_breakend(true);
        end_be2 = c2.open_breakend(false);
        uniform_jcn = jcn_match(c1.jcn(), c1.jcn_uncertainty(), c2.jcn(), c2.jcn_uncertainty());
    } else if cluster_chains.len() == 1 && unchained_sv_count == 1 {
        let var = unchained_sv.clone().unwrap();
        let chain = &cluster_chains[0];
        start_be1 = chain.open_breakend(true);
        end_be1 = chain.open_breakend(false);
        start_be2 = var.breakend(true);
        end_be2 = var.breakend(false);
        uniform_jcn = jcn_match(var.jcn(), var.jcn_uncertainty(), chain.jcn(), chain.jcn_uncertainty());
    } else {
        return;
    }
    // Java dereferences each of these (a null throws)
    let (s1, e1, s2, e2) = (start_be1.unwrap(), end_be1.unwrap(), start_be2.unwrap(), end_be2.unwrap());
    if !is_germline && are_interrupted_breakends(&s1, &e1, &s2, &e2, chr_breakend_map) {
        return;
    }
    if s1.chromosome() == e1.chromosome() && s2.chromosome() == e2.chromosome() && s1.chromosome() == s2.chromosome() {
        if !is_single_chain {
            return;
        }
        if s1.orientation() == e1.orientation() && s2.orientation() == e2.orientation() && s1.orientation() != s2.orientation() {
            classify_inversion_pair_clusters(
                cluster,
                long_del_threshold,
                long_dup_threshold,
                &s1,
                &e1,
                &s2,
                &e2,
                uniform_jcn,
                long_ti_link.as_ref().expect("long TI link"),
            );
        }
    } else {
        classify_translocation_pair_clusters(
            cluster,
            long_del_threshold,
            long_dup_threshold,
            &s1,
            &e1,
            &s2,
            &e2,
            uniform_jcn,
            long_ti_link.as_ref(),
            if existing_chain_modified { Some(&cluster_chains) } else { None },
        );
    }
}

fn classify_synthetic_del_dups(cluster: &Cluster, long_del_threshold: i64, long_dup_threshold: i64) {
    if !cluster.is_fully_chained(true) || cluster.chain_count() != 1 || cluster.sgl_breakend_count() > 0 {
        return;
    }
    let chain = cluster.chains()[0].clone();
    let cs = chain.open_breakend(true).unwrap();
    let ce = chain.open_breakend(false).unwrap();
    if cs.chromosome() != ce.chromosome() || cs.arm() != ce.arm() {
        return;
    }
    if cs.orientation() == ce.orientation() {
        return;
    }
    let face_away = (cs.position() < ce.position()) == cs.forward_orient();
    let total = chain.length(false);
    let synthetic_length = (cs.position() - ce.position()).abs() as i64;
    let rt = if face_away { ResolvedType::DEL } else { ResolvedType::DUP };
    hmf_common::debug!(
        "cluster({}) chain(links=({} len={} synLen({}) marked as {}",
        cluster.id(),
        chain.link_count(),
        total,
        synthetic_length,
        rt
    );
    let within = if rt == ResolvedType::DEL { synthetic_length < long_del_threshold } else { synthetic_length < long_dup_threshold };
    cluster.set_resolved(within, rt);
    set_pair_length_data(cluster);
}

fn set_pair_length_data(cluster: &Cluster) {
    if cluster.chain_count() != 1 {
        return;
    }
    let a = get_synthetic_length(cluster);
    let b = get_synthetic_ti_length(cluster);
    let c = get_synthetic_gap_length(cluster);
    cluster.add_annotation(&format!("PairLen={a};{b};{c}"));
}

fn is_loh_bounded_ti(pair: &Pair) -> bool {
    let (f, s) = (Some(pair.first_breakend()), Some(pair.second_breakend()));
    for loh in pair.first().cl().loh_events() {
        let (ls, le) = (loh.breakend(true), loh.breakend(false));
        if ls == f || ls == s || le == f || le == s {
            return true;
        }
    }
    false
}

#[allow(clippy::too_many_arguments)]
fn classify_translocation_pair_clusters(
    cluster: &Cluster, long_del_threshold: i64, long_dup_threshold: i64, start_be1: &Be, end_be1: &Be, start_be2: &Be, end_be2: &Be,
    uniform_jcn: bool, longest_ti_pair: Option<&Pair>, rearranged_chains: Option<&Vec<Chain>>,
) {
    let consistency = calc_consistency_be(start_be1) + calc_consistency_be(end_be1) + calc_consistency_be(start_be2) + calc_consistency_be(end_be2);
    if consistency != 0 {
        return;
    }
    let (mut arm1_be1, mut arm2_be1, mut arm1_be2, mut arm2_be2): (Option<&Be>, Option<&Be>, Option<&Be>, Option<&Be>) = (None, None, None, None);
    if start_be1.chr_arm() != end_be1.chr_arm() && start_be2.chr_arm() != end_be2.chr_arm() {
        arm1_be1 = Some(start_be1);
        arm2_be1 = Some(end_be1);
        if start_be2.chr_arm() == start_be1.chr_arm() && end_be2.chr_arm() == end_be1.chr_arm() {
            arm1_be2 = Some(start_be2);
            arm2_be2 = Some(end_be2);
        } else if end_be2.chr_arm() == start_be1.chr_arm() && start_be2.chr_arm() == end_be1.chr_arm() {
            arm1_be2 = Some(end_be2);
            arm2_be2 = Some(start_be2);
        }
    }
    let (Some(arm1_be1), Some(arm2_be1), Some(arm1_be2), Some(arm2_be2)) = (arm1_be1, arm2_be1, arm1_be2, arm2_be2) else {
        if longest_ti_pair.is_some() && cluster.chain_count() == 1 && cluster.unlinked_sv_count() == 0 {
            let chain = &cluster.chains()[0];
            if chain.open_breakend(true).unwrap().chr_arm() != chain.open_breakend(false).unwrap().chr_arm() {
                cluster.set_resolved(false, UNBAL_TRANS_TI);
            }
        }
        return;
    };
    let arm1_db1 = arm1_be1.db_link();
    let arm2_db1 = arm2_be1.db_link();
    let arm1_matching_db = arm1_db1.is_some() && arm1_db1 == arm1_be2.db_link();
    let arm2_matching_db = arm2_db1.is_some() && arm2_db1 == arm2_be2.db_link();
    let rt;
    let mut is_resolved = true;
    if arm1_matching_db && arm2_matching_db {
        rt = RECIP_TRANS;
        let (l1, l2) = (arm1_db1.as_ref().unwrap().length(), arm2_db1.as_ref().unwrap().length());
        is_resolved = (l1 as i64) <= long_del_threshold && (l2 as i64) <= long_del_threshold;
        cluster.add_annotation(&format!("PairLen={l1};{l2};{NO_LENGTH}"));
    } else {
        let Some(lt) = longest_ti_pair else { return };
        set_pair_length_data(cluster);
        let loh_bounded = is_loh_bounded_ti(lt);
        rt = if !arm1_matching_db && !arm2_matching_db {
            if loh_bounded && uniform_jcn { DUP_TI } else { RECIP_TRANS_DUPS }
        } else if loh_bounded && uniform_jcn {
            DEL_TI
        } else {
            RECIP_TRANS_DEL_DUP
        };
        if lt.base_length() as i64 > long_dup_threshold {
            is_resolved = false;
        } else {
            let non1 = if lt.has_breakend(start_be1) { end_be1 } else { start_be1 };
            let non2 = if lt.has_breakend(start_be2) { end_be2 } else { start_be2 };
            let distance = (non1.position() - non2.position()).abs() as i64;
            if rt == DEL_TI || rt == RECIP_TRANS_DEL_DUP {
                if distance > long_del_threshold {
                    is_resolved = false;
                }
            } else if distance > long_dup_threshold {
                is_resolved = false;
            }
        }
    }
    hmf_common::debug!(
        "cluster({}) longestTI({}) resolvedType({}) isResolved({})",
        cluster.id(),
        longest_ti_pair.map(|p| p.base_length().to_string()).unwrap_or_else(|| "none".into()),
        rt,
        is_resolved
    );
    if let Some(chains) = rearranged_chains {
        let existing = cluster.chains()[0].clone();
        cluster.dissolve_links_and_chains();
        for nc in chains {
            nc.set_jcn_data(existing.jcn(), existing.jcn_uncertainty());
            cluster.add_chain(nc, true);
        }
    }
    cluster.set_resolved(is_resolved, rt);
}

#[allow(clippy::too_many_arguments)]
fn classify_inversion_pair_clusters(
    cluster: &Cluster, long_del_threshold: i64, long_dup_threshold: i64, start_be1: &Be, end_be1: &Be, start_be2: &Be, end_be2: &Be,
    uniform_jcn: bool, longest_ti_pair: &Pair,
) {
    let lower_be1 = if start_be1.position() < end_be1.position() { start_be1 } else { end_be1 };
    let upper_be1 = if lower_be1 == start_be1 { end_be1 } else { start_be1 };
    let lower_be2 = if start_be2.position() < end_be2.position() { start_be2 } else { end_be2 };
    let upper_be2 = if lower_be2 == start_be2 { end_be2 } else { start_be2 };
    if upper_be1.position() < lower_be2.position() || upper_be2.position() < lower_be1.position() {
        let facing = (upper_be1.position() < lower_be2.position()) == upper_be1.reverse_orient();
        if !facing {
            return;
        }
        cluster.set_resolved(false, FB_INV_PAIR);
        return;
    }
    let lower_db1 = lower_be1.db_link();
    let upper_db1 = upper_be1.db_link();
    let lower_matching = lower_db1.is_some() && lower_db1 == lower_be2.db_link();
    let upper_matching = upper_db1.is_some() && upper_db1 == upper_be2.db_link();
    if lower_matching && upper_matching {
        let (l1, l2) = (lower_db1.unwrap().length(), upper_db1.unwrap().length());
        let is_resolved = (l1 as i64) <= long_del_threshold && (l2 as i64) <= long_del_threshold;
        cluster.set_resolved(is_resolved, RECIP_INV);
        cluster.add_annotation(&format!("PairLen={l1};{l2};{NO_LENGTH}"));
        return;
    } else if (lower_be1.forward_orient() && lower_be1.position() < lower_be2.position() && upper_be1.position() < upper_be2.position())
        || (lower_be2.forward_orient() && lower_be2.position() < lower_be1.position() && upper_be2.position() < upper_be1.position())
    {
        return;
    }
    set_pair_length_data(cluster);
    let loh_bounded = is_loh_bounded_ti(longest_ti_pair);
    let rt;
    if (lower_be1.position() < lower_be2.position() && upper_be1.position() > upper_be2.position())
        || (lower_be2.position() < lower_be1.position() && upper_be2.position() > upper_be1.position())
    {
        if loh_bounded && uniform_jcn {
            rt = DEL_TI;
        } else {
            let inner = (upper_be1.position() - lower_be1.position()).min(upper_be2.position() - lower_be2.position());
            if inner < MIN_SIMPLE_DUP_DEL_CUTOFF {
                rt = RESOLVED_FOLDBACK;
            } else {
                let chain = cluster.chains()[0].clone();
                let to_switch = if lower_be1.position() < lower_be2.position() {
                    if longest_ti_pair.has_breakend(lower_be2) { lower_be2 } else { upper_be2 }
                } else if longest_ti_pair.has_breakend(lower_be1) {
                    lower_be1
                } else {
                    upper_be1
                };
                if reverse_section_on_breakend(&chain, to_switch) {
                    hmf_common::debug!("cluster({}) reconfigured chain:", cluster.id());
                }
                rt = RECIP_INV_DEL_DUP;
            }
        }
    } else if loh_bounded && uniform_jcn {
        rt = DUP_TI;
    } else {
        let chain = cluster.chains()[0].clone();
        // Java uses lowerBe1's partner in both branches
        let to_switch = if lower_be1.position() > lower_be2.position() {
            if longest_ti_pair.has_breakend(lower_be1) { lower_be1.clone() } else { longest_ti_pair.other_breakend(lower_be1) }
        } else if longest_ti_pair.has_breakend(lower_be2) {
            lower_be2.clone()
        } else {
            longest_ti_pair.other_breakend(lower_be1)
        };
        if reverse_section_on_breakend(&chain, &to_switch) {
            hmf_common::debug!("cluster({}) reconfigured chain:", cluster.id());
        }
        rt = RECIP_INV_DUPS;
    }
    let chain = cluster.chains()[0].clone();
    let longest_ti_length = chain.linked_pairs().iter().map(|x| x.position_distance() as i64).max().expect("chain links");
    let cs = chain.open_breakend(true).unwrap();
    let ce = chain.open_breakend(false).unwrap();
    let non_ti_distance = (cs.position() - ce.position()).abs() as i64;
    let is_resolved = if rt == RECIP_INV_DUPS || rt == DUP_TI {
        longest_ti_length <= long_dup_threshold && non_ti_distance <= long_dup_threshold
    } else if rt == RECIP_INV_DEL_DUP || rt == DEL_TI {
        longest_ti_length <= long_dup_threshold && non_ti_distance <= long_del_threshold
    } else {
        false
    };
    hmf_common::debug!(
        "cluster({}) longestTI({}) resolvedType({}) isResolved({})",
        cluster.id(),
        longest_ti_pair.base_length(),
        rt,
        is_resolved
    );
    cluster.set_resolved(is_resolved, rt);
}

fn is_interrupted_ti(pair: &Pair, map: &JavaHashMap<String, Vec<Be>>) -> bool {
    is_interrupted_breakend_pair(&pair.breakend(true), &pair.breakend(false), map)
}

fn are_interrupted_breakends(s1: &Be, e1: &Be, s2: &Be, e2: &Be, map: &JavaHashMap<String, Vec<Be>>) -> bool {
    (s1.chr_arm() == s2.chr_arm() && is_interrupted_breakend_pair(s1, s2, map))
        || (s1.chr_arm() == e2.chr_arm() && is_interrupted_breakend_pair(s1, e2, map))
        || (e1.chr_arm() == s2.chr_arm() && is_interrupted_breakend_pair(e1, s2, map))
        || (e1.chr_arm() == e2.chr_arm() && is_interrupted_breakend_pair(e1, e2, map))
}

fn is_interrupted_breakend_pair(be1: &Be, be2: &Be, map: &JavaHashMap<String, Vec<Be>>) -> bool {
    let list = map.get(&be1.chromosome().to_string()).expect("breakends of chromosome");
    let lower = if be1.position() < be2.position() { be1 } else { be2 };
    let upper = if lower == be1 { be2 } else { be1 };
    let start_index = lower.chr_pos_index();
    let end_index = upper.chr_pos_index();
    for i in (start_index + 1)..end_index {
        let var = list[i as usize].sv();
        if var == be1.sv() || var == be2.sv() || var.cluster() == be1.sv().cluster() {
            continue;
        }
        if var.cl().sv_count() > 1 {
            return true;
        }
        if !var.is_simple_type() {
            return true;
        }
        if var.be(true).chr_pos_index() < start_index || var.be(false).chr_pos_index() > end_index {
            return true;
        }
    }
    false
}
