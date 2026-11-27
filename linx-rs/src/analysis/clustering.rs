//! ClusteringState, ClusteringPrep, SvFiltering, SimpleClustering, ComplexClustering, FoldbackFinder and
//! the arm clusters (ArmCluster.buildArmClusters / setFeatures).

use std::cell::Cell;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

use hmf_common::java_hash::JavaHashMap;
use hmf_common::sv::StructuralVariantData;
use hmf_common::sv::StructuralVariantType::{self as SvType, *};

use crate::annotators::{has_line_repeat_class, has_poly_a_or_t_motif};
use crate::analysis::classification::{get_synthetic_length, is_simple_single_sv};
use crate::chaining::jcn_limits::{jcn_match_breakends, jcn_match_svs};
use crate::chaining::utils::{breakends_are_chained, breakends_cannot_link, get_min_templated_insertion_length, have_linked_assemblies, CHAIN_ASSEMBLY_LINK_COUNT, CHAIN_LENGTH, CHAIN_LINK_COUNT};
use crate::cn::{HomLoss, Loh, CN_DATA_NO_SV};
use crate::idmap::IdMap;
use crate::types::ClusteringReason::*;
use crate::types::*;
use crate::utils::*;

/// ClusteringState
#[derive(Default)]
pub struct ClusteringState {
    next_cluster_id: Cell<i32>,
    /// every breakend on a chromosome, by ascending position
    pub chr_breakend_map: JavaHashMap<String, Vec<Be>>,
    pub loh_event_list: Vec<Loh>,
    pub hom_loss_list: Vec<HomLoss>,
    pub del_cutoff_length: i32,
    pub dup_cutoff_length: i32,
}

impl ClusteringState {
    pub fn next_cluster_id(&self) -> i32 {
        let id = self.next_cluster_id.get();
        self.next_cluster_id.set(id + 1);
        id
    }
    pub fn reset(&mut self) {
        self.next_cluster_id.set(0);
        self.del_cutoff_length = 0;
        self.dup_cutoff_length = 0;
        self.chr_breakend_map.clear();
    }
    pub fn breakend_list(&self, chr: &str) -> Option<&Vec<Be>> {
        self.chr_breakend_map.get(&chr.to_string())
    }
}

// ClusteringPrep

const DEL_DUP_LENGTH_TRIM_COUNT: i32 = 5;
const MAX_ARM_COUNT: i32 = 41;

pub fn populate_chromosome_breakend_map(all: &[Sv], state: &mut ClusteringState) {
    for var in all {
        add_sv_to_chr_breakend_map(var, &mut state.chr_breakend_map);
    }
    for (_, list) in state.chr_breakend_map.iter() {
        for (i, b) in list.iter().enumerate() {
            b.set_chr_pos_index(i as i32);
        }
    }
}

/// ClusteringPrep.linkSglMappedInferreds: an SGL whose alignment maps to an INF becomes one SV
pub fn link_sgl_mapped_inferreds(all: &mut Vec<Sv>) {
    let sgls: Vec<Sv> = all.iter().filter(|x| x.sv_type() == SGL && !x.sgl_mappings().is_empty()).cloned().collect();
    let infs: Vec<Sv> = all.iter().filter(|x| x.sv_type() == INF).cloned().collect();
    for sgl in &sgls {
        let mut matched = false;
        for mapping in sgl.sgl_mappings() {
            for inf in &infs {
                if mapping.chromosome != inf.chromosome(true) {
                    continue;
                }
                if mapping.orientation != inf.orientation(true) {
                    continue;
                }
                if (mapping.position - inf.position(true)).abs() > 1000 {
                    continue;
                }
                let new_var = merge_sgl_mapped_inferred(sgl, inf, &mapping, sgl.id());
                hmf_common::debug!("new SV({}) from sgl({}) with mapping to inf({})", new_var, sgl.pos_id(), inf.pos_id());
                all.push(new_var.clone());
                new_var.set_linked_svs(sgl, inf);
                new_var.add_cluster_reason(SGL_MAPPING_INF, inf.id());
                sgl.set_linked_svs(&new_var, inf);
                inf.set_linked_svs(&new_var, sgl);
                matched = true;
                break;
            }
            if matched {
                break;
            }
        }
    }
}

pub fn merge_sgl_mapped_inferred(sgl: &Sv, inf: &Sv, mapping: &SglMapping, sv_id: i32) -> Sv {
    let s = sgl.sv_data();
    let i = inf.sv_data();
    let new_type = if s.start_chromosome == i.start_chromosome {
        if s.start_orientation == i.start_orientation {
            INV
        } else if (s.start_orientation == 1) == (s.start_position < i.start_position) {
            DEL
        } else {
            DUP
        }
    } else {
        BND
    };
    let d = StructuralVariantData {
        id: sv_id,
        vcf_id_start: s.vcf_id_start.clone(),
        vcf_id_end: String::new(),
        start_chromosome: s.start_chromosome.clone(),
        end_chromosome: i.start_chromosome.clone(),
        start_position: s.start_position,
        end_position: mapping.position,
        start_orientation: s.start_orientation,
        end_orientation: mapping.orientation,
        start_homology_sequence: s.start_homology_sequence.clone(),
        end_homology_sequence: i.start_homology_sequence.clone(),
        junction_copy_number: s.junction_copy_number,
        start_af: s.start_af,
        end_af: i.start_af,
        adjusted_start_af: s.adjusted_start_af,
        adjusted_end_af: i.adjusted_start_af,
        adjusted_start_copy_number: s.adjusted_start_copy_number,
        adjusted_end_copy_number: i.adjusted_start_copy_number,
        adjusted_start_copy_number_change: s.adjusted_start_copy_number_change,
        adjusted_end_copy_number_change: i.adjusted_start_copy_number_change,
        insert_sequence: s.insert_sequence.clone(),
        sv_type: Some(new_type),
        filter: s.filter.clone(),
        quality_score: s.quality_score,
        event: s.event.clone(),
        start_tumor_variant_fragment_count: s.start_tumor_variant_fragment_count,
        start_tumor_reference_fragment_count: s.start_tumor_reference_fragment_count,
        start_normal_variant_fragment_count: s.start_normal_variant_fragment_count,
        start_normal_reference_fragment_count: s.start_normal_reference_fragment_count,
        end_tumor_variant_fragment_count: i.start_tumor_variant_fragment_count,
        end_tumor_reference_fragment_count: i.start_tumor_reference_fragment_count,
        end_normal_variant_fragment_count: i.start_normal_variant_fragment_count,
        end_normal_reference_fragment_count: i.start_normal_reference_fragment_count,
        start_interval_offset_start: s.start_interval_offset_start,
        start_interval_offset_end: s.start_interval_offset_end,
        end_interval_offset_start: i.start_interval_offset_start,
        end_interval_offset_end: i.start_interval_offset_end,
        inexact_homology_offset_start: s.inexact_homology_offset_start,
        inexact_homology_offset_end: i.inexact_homology_offset_start,
        start_linked_by: s.start_linked_by.clone(),
        end_linked_by: i.start_linked_by.clone(),
        insert_sequence_alignments: s.insert_sequence_alignments.clone(),
        insert_sequence_repeat_class: s.insert_sequence_repeat_class.clone(),
        insert_sequence_repeat_type: s.insert_sequence_repeat_type.clone(),
        insert_sequence_repeat_orientation: s.insert_sequence_repeat_orientation,
        insert_sequence_repeat_coverage: s.insert_sequence_repeat_coverage,
        start_anchoring_support_distance: s.start_anchoring_support_distance,
        end_anchoring_support_distance: i.start_anchoring_support_distance,
        pon_count: 0,
    };
    Sv::from_data(Rc::new(d))
}

pub fn set_simple_variant_lengths(state: &mut ClusteringState) {
    let mut del_cutoff = 0;
    let mut dup_cutoff = 0;
    let mut del_lengths: Vec<i32> = Vec::new();
    let mut dup_lengths: Vec<i32> = Vec::new();
    let mut simple_arm_count = 0;
    for (_, list) in state.chr_breakend_map.iter() {
        let mut p_inv = false;
        let mut q_inv = false;
        for b in list {
            let var = b.sv();
            if var.sv_type() != INV {
                continue;
            }
            if !p_inv && (var.arm(true) == Arm::P || var.arm(false) == Arm::P) {
                p_inv = true;
            }
            if !q_inv && (var.arm(true) == Arm::Q || var.arm(false) == Arm::Q) {
                q_inv = true;
            }
            if p_inv && q_inv {
                break;
            }
        }
        if p_inv && q_inv {
            continue;
        }
        if !p_inv {
            simple_arm_count += 1;
        }
        if !q_inv {
            simple_arm_count += 1;
        }
        for b in list {
            if !b.uses_start() || !(b.sv_type() == DEL || b.sv_type() == DUP) {
                continue;
            }
            let var = b.sv();
            if p_inv && (var.arm(true) == Arm::P || var.arm(false) == Arm::P) {
                continue;
            }
            if q_inv && (var.arm(true) == Arm::Q || var.arm(false) == Arm::Q) {
                continue;
            }
            if var.sv_type() == DEL {
                del_lengths.push(var.length());
            } else if var.sv_type() == DUP {
                dup_lengths.push(var.length());
            }
        }
    }
    let trim_count = jround(simple_arm_count as f64 / MAX_ARM_COUNT as f64 * DEL_DUP_LENGTH_TRIM_COUNT as f64) as i32;
    if del_lengths.len() as i32 > trim_count {
        del_lengths.sort();
        del_cutoff = del_lengths[del_lengths.len() - trim_count as usize - 1];
    }
    if dup_lengths.len() as i32 > trim_count {
        dup_lengths.sort();
        dup_cutoff = dup_lengths[dup_lengths.len() - trim_count as usize - 1];
    }
    del_cutoff = del_cutoff.max(MIN_SIMPLE_DUP_DEL_CUTOFF).min(MAX_SIMPLE_DUP_DEL_CUTOFF);
    dup_cutoff = dup_cutoff.max(MIN_SIMPLE_DUP_DEL_CUTOFF).min(MAX_SIMPLE_DUP_DEL_CUTOFF);
    hmf_common::debug!(
        "simple dels count({}) cutoff-length({}), dups count({}) cutoff-length({}) simpleArms({}) trimCount({})",
        del_lengths.len(),
        del_cutoff,
        dup_lengths.len(),
        dup_cutoff,
        simple_arm_count,
        trim_count
    );
    state.del_cutoff_length = del_cutoff;
    state.dup_cutoff_length = dup_cutoff;
}

pub fn annotate_nearest_sv_data(map: &JavaHashMap<String, Vec<Be>>) {
    for (_, list) in map.iter() {
        let n = list.len();
        for i in 0..n {
            let b = &list[i];
            let var = b.sv();
            let prev = if i > 0 { Some(&list[i - 1]) } else { None };
            let next = if i < n - 1 { Some(&list[i + 1]) } else { None };
            let mut closest = -1;
            if let Some(p) = prev.filter(|p| p.sv() != var) {
                closest = b.position() - p.position();
            }
            if let Some(nx) = next.filter(|nx| nx.sv() != var) {
                let d = nx.position() - b.position();
                if closest < 0 || d < closest {
                    closest = d;
                }
            }
            if closest < 0 {
                continue;
            }
            if var.nearest_sv_distance() >= 0 && closest >= var.nearest_sv_distance() {
                continue;
            }
            var.set_nearest_sv_distance(closest);
            if var.sv_type() == BND || var.is_sgl_breakend() || var.sv_type() == INS {
                var.set_nearest_sv_relation(RELATION_TYPE_NEIGHBOUR);
            } else if next.is_some_and(|nx| nx.sv() == var) || prev.is_some_and(|p| p.sv() == var) {
                var.set_nearest_sv_relation(RELATION_TYPE_NEIGHBOUR);
            } else {
                var.set_nearest_sv_relation(RELATION_TYPE_OVERLAP);
            }
        }
    }
}

pub fn associate_breakend_cn_events(sample_id: &str, state: &ClusteringState) {
    let mut missed = 0;
    let mut current_chr = String::new();
    let mut list: Option<&Vec<Be>> = None;
    let mut list_set = false;
    for loh in &state.loh_event_list {
        if !(loh.is_sv_event(true) || loh.is_sv_event(false)) {
            continue;
        }
        if !list_set || list.is_none() || current_chr != loh.chromosome {
            list = state.breakend_list(&loh.chromosome);
            current_chr = loh.chromosome.clone();
            list_set = true;
        }
        let Some(bl) = list else { continue };
        for b in bl {
            let var = b.sv();
            if b.forward_orient() && var.id() == loh.start_sv {
                let skip_inv = var.sv_type() == INV
                    && (b.position() - loh.pos_start).abs() > (b.other_breakend().unwrap().position() - loh.pos_start).abs();
                if !skip_inv {
                    loh.set_breakend(Some(b.clone()), true);
                    var.cl().add_loh_event(loh);
                }
            }
            if b.reverse_orient() && var.id() == loh.end_sv {
                let skip_inv =
                    var.sv_type() == INV && (b.position() - loh.pos_end).abs() > (b.other_breakend().unwrap().position() - loh.pos_end).abs();
                if !skip_inv {
                    loh.set_breakend(Some(b.clone()), false);
                    var.cl().add_loh_event(loh);
                }
            }
            if loh.matched_both_svs() {
                break;
            }
        }
        if loh.start_sv != CN_DATA_NO_SV && loh.breakend(true).is_none() {
            missed += 1;
        }
        if loh.end_sv != CN_DATA_NO_SV && loh.breakend(false).is_none() {
            missed += 1;
        }
    }
    for hl in &state.hom_loss_list {
        if hl.start_sv == CN_DATA_NO_SV && hl.end_sv == CN_DATA_NO_SV {
            continue;
        }
        let Some(bl) = state.breakend_list(&hl.chromosome) else { continue };
        for b in bl {
            if b.forward_orient() && b.sv().id() == hl.start_sv {
                hl.set_breakend(Some(b.clone()), true);
            }
            if b.reverse_orient() && b.sv().id() == hl.end_sv {
                hl.set_breakend(Some(b.clone()), false);
            }
            if hl.matched_both_svs() {
                break;
            }
        }
        if hl.start_sv != CN_DATA_NO_SV && hl.breakend(true).is_none() {
            missed += 1;
        }
        if hl.end_sv != CN_DATA_NO_SV && hl.breakend(false).is_none() {
            missed += 1;
        }
    }
    if missed > 0 {
        hmf_common::debug!("sample({}) missed {} links to LOH and hom-loss events", sample_id, missed);
    }
}

// SvFiltering

const LOW_VAF_THRESHOLD: f64 = 0.05;
const SHORT_INV_DISTANCE: i32 = 100;
const ISOLATED_BND_DISTANCE: i32 = 5000;
const PERMITTED_SGL_DUP_BE_DISTANCE: i32 = 1;

/// SvFiltering
#[derive(Default)]
pub struct SvFiltering {
    /// SV and exclusion reason: a HashMap keyed by SvVarData (LX1)
    excluded_svs: IdMap<Sv, ResolvedType>,
    non_sgl_duplicate_breakends: HashMap<Be, Be>,
}

impl SvFiltering {
    pub fn apply_filters(&mut self, state: &mut ClusteringState) {
        self.excluded_svs = IdMap::new();
        self.non_sgl_duplicate_breakends.clear();
        let mut removal: HashMap<String, HashSet<Be>> = HashMap::new();
        let entries: Vec<(String, Vec<Be>)> = state.chr_breakend_map.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        for (chromosome, list) in &entries {
            removal.entry(chromosome.clone()).or_default();
            let n = list.len();
            for i in 0..n {
                let b = &list[i];
                let var = b.sv().clone();
                if var.is_sgl_breakend() && var.linked_svs().is_some() {
                    self.excluded_svs.insert(var.clone(), ResolvedType::SGL_MAPPED_INF);
                    removal.get_mut(chromosome).unwrap().insert(b.clone());
                    continue;
                }
                if (var.sv_type() == BND || var.sv_type() == SGL) && is_isolated_low_vaf_bnd(&var, state) {
                    removal.get_mut(chromosome).unwrap().insert(b.clone());
                    if var.sv_type() == BND {
                        let o = b.other_breakend().unwrap();
                        removal.entry(o.chromosome().to_string()).or_default().insert(o);
                    }
                    self.excluded_svs.insert(var.clone(), ResolvedType::LOW_VAF);
                    continue;
                }
                if i >= n - 1 {
                    break;
                }
                let next = &list[i + 1];
                let next_var = next.sv().clone();
                if var.sv_type() == INV && next_var == var && is_low_vaf_inversion(b, next) {
                    let r = removal.get_mut(chromosome).unwrap();
                    r.insert(b.clone());
                    r.insert(next.clone());
                    self.excluded_svs.insert(var.clone(), ResolvedType::LOW_VAF);
                    continue;
                }
                if var.sv_type() == INF
                    && next_var.sv_type() == INF
                    && b.orientation() != next.orientation()
                    && jcn_match_svs(&var, &next_var)
                    && !self.excluded_svs.contains_key(&var)
                    && !self.excluded_svs.contains_key(&next_var)
                {
                    let r = removal.get_mut(chromosome).unwrap();
                    r.insert(b.clone());
                    r.insert(next.clone());
                    self.excluded_svs.insert(var.clone(), ResolvedType::PAIR_INF);
                    self.excluded_svs.insert(next_var.clone(), ResolvedType::PAIR_INF);
                    continue;
                }
                if next_var == var {
                    continue;
                }
                if b.orientation() != next.orientation() {
                    continue;
                }
                let distance = next.position() - b.position();
                if distance <= PERMITTED_SGL_DUP_BE_DISTANCE && (var.sv_type() == SGL || next_var.sv_type() == SGL) {
                    if var.sv_type() == SGL {
                        self.excluded_svs.insert(var.clone(), ResolvedType::DUP_BE);
                        removal.get_mut(chromosome).unwrap().insert(b.clone());
                    } else if next_var.sv_type() == SGL {
                        self.excluded_svs.insert(next_var.clone(), ResolvedType::DUP_BE);
                        removal.get_mut(chromosome).unwrap().insert(next.clone());
                    }
                } else if !var.is_sgl_breakend() && !next_var.is_sgl_breakend() {
                    self.check_candidate_spanning_variant(b, next, chromosome, &mut removal);
                    self.check_candidate_spanning_variant(next, b, chromosome, &mut removal);
                }
            }
        }
        for (chr, set) in &removal {
            if set.is_empty() {
                continue;
            }
            let list = state.chr_breakend_map.get_mut(chr).expect("breakends of chromosome");
            list.retain(|b| !set.contains(b));
            for (i, b) in list.iter().enumerate() {
                b.set_chr_pos_index(i as i32);
            }
        }
    }

    fn check_candidate_spanning_variant(&mut self, b: &Be, next: &Be, chromosome: &str, removal: &mut HashMap<String, HashSet<Be>>) {
        if b.sv().sv_data().insert_sequence.len() as i32 >= MIN_TEMPLATED_INSERTION_LENGTH && !next.sv().assembly_links(!next.uses_start()).is_empty() {
            let lo = next.position() - 1;
            let hi = next.position() + 1;
            if b.position() >= lo && b.position() <= hi {
                let other = b.other_breakend().unwrap();
                let Some(next2) = self.non_sgl_duplicate_breakends.get(&other).cloned() else {
                    self.non_sgl_duplicate_breakends.insert(b.clone(), next.clone());
                    return;
                };
                if have_linked_assemblies(next.sv(), next2.sv(), !next.uses_start(), !next2.uses_start()) {
                    let var = b.sv().clone();
                    removal.get_mut(chromosome).unwrap().insert(b.clone());
                    if var.sv_type() == BND {
                        removal.entry(other.chromosome().to_string()).or_default().insert(other.clone());
                    } else {
                        removal.get_mut(chromosome).unwrap().insert(other.clone());
                    }
                    self.excluded_svs.insert(var, ResolvedType::DUP_BE);
                }
            }
        }
    }

    pub fn cluster_excluded_variants(&self, clusters: &mut Vec<Cluster>, state: &ClusteringState) {
        for (var, reason) in self.excluded_svs.iter() {
            if let Some(c) = var.cluster() {
                if !contains(clusters, &c) {
                    clusters.push(c);
                }
                continue;
            }
            let c = Cluster::create(state.next_cluster_id());
            c.add_variant(var);
            c.set_resolved(true, *reason);
            clusters.push(c);
        }
    }
}

fn is_isolated_low_vaf_bnd(var: &Sv, state: &ClusteringState) -> bool {
    if !has_low_jcn(var) {
        return false;
    }
    if has_poly_a_or_t_motif(var) || has_line_repeat_class(var) {
        return false;
    }
    for s in [SE_START, SE_END] {
        if s == SE_END && var.is_sgl_breakend() {
            continue;
        }
        let b = var.breakend_se(s).unwrap();
        let list = state.breakend_list(b.chromosome()).unwrap();
        let i = b.chr_pos_index();
        let down = i == 0 || b.position() - list[(i - 1) as usize].position() >= ISOLATED_BND_DISTANCE;
        let up = i >= list.len() as i32 - 1 || list[(i + 1) as usize].position() - b.position() >= ISOLATED_BND_DISTANCE;
        if !down || !up {
            return false;
        }
    }
    true
}

fn is_low_vaf_inversion(b: &Be, next: &Be) -> bool {
    let var = b.sv();
    if next.sv() != var || next.position() - b.position() > SHORT_INV_DISTANCE {
        return false;
    }
    has_low_jcn(var) || var.calc_vaf(true) < LOW_VAF_THRESHOLD
}

// SimpleClustering

pub fn has_low_jcn(var: &SvVarData) -> bool {
    if var.sv_type() == INS {
        return false;
    }
    var.jcn_max() < LOW_JCN_THRESHOLD
}

pub fn add_cluster_reasons(var1: &Sv, var2: &Sv, reason: ClusteringReason) {
    var1.add_cluster_reason(reason, var2.id());
    var2.add_cluster_reason(reason, var1.id());
}

pub fn variants_have_different_jcn(var1: &SvVarData, var2: &SvVarData) -> bool {
    let low1 = var1.jcn() < LOW_JCN_THRESHOLD;
    let low2 = var2.jcn() < LOW_JCN_THRESHOLD;
    if low1 == low2 {
        return false;
    }
    if (var1.jcn() - var2.jcn()).abs() < MAX_COPY_NUM_DIFF {
        return false;
    }
    (low1 && var1.jcn_max() < var2.jcn_min()) || (low2 && var2.jcn_max() < var1.jcn_min())
}

pub fn breakends_have_different_jcn(b1: &SvBreakend, b2: &SvBreakend) -> bool {
    variants_have_different_jcn(b1.sv(), b2.sv())
}

fn breakends_violate_loh_hom_loss(b: &Be, other: &Be) -> bool {
    if b.chromosome() != other.chromosome() {
        return false;
    }
    let events: Vec<Loh> = b.cl().loh_events().into_iter().filter(|x| x.breakend(true).as_ref() == Some(b) || x.breakend(false).as_ref() == Some(b)).collect();
    if events.is_empty() {
        return false;
    }
    events.iter().any(|x| other.position() > x.pos_start && other.position() < x.pos_end)
}

pub fn variants_violate_loh_hom_loss(var1: &Sv, var2: &Sv) -> bool {
    for s1 in [SE_START, SE_END] {
        if s1 == SE_END && var1.is_sgl_breakend() {
            continue;
        }
        for s2 in [SE_START, SE_END] {
            if s2 == SE_END && var2.is_sgl_breakend() {
                continue;
            }
            let (b1, b2) = (var1.breakend_se(s1).unwrap(), var2.breakend_se(s2).unwrap());
            if breakends_violate_loh_hom_loss(&b1, &b2) || breakends_violate_loh_hom_loss(&b2, &b1) {
                return true;
            }
        }
    }
    false
}

fn exceeds_dup_del_cutoff_length(t: SvType, length: i32, state: &ClusteringState) -> bool {
    match t {
        DEL => length > state.del_cutoff_length,
        DUP => length > state.dup_cutoff_length,
        _ => false,
    }
}

fn breakends_in_close_link(var1: &Sv, var2: &Sv) -> [bool; 2] {
    let mut data = [false, false];
    for s1 in [SE_START, SE_END] {
        if s1 == SE_END && var1.is_sgl_breakend() {
            continue;
        }
        let b1 = var1.breakend_se(s1).unwrap();
        for s2 in [SE_START, SE_END] {
            // Java: getBreakend(se2) of an SGL is null and its getChrArm throws
            let b2 = var2.breakend_se(s2).expect("breakend");
            if b1.chr_arm() != b2.chr_arm() {
                continue;
            }
            if b1.orientation() != b2.orientation() {
                data[0] = true;
                if (b1.position() - b2.position()).abs() < MAX_MERGE_DISTANCE {
                    data[1] = true;
                    return data;
                }
            }
        }
    }
    data
}

pub fn skip_cluster_type(cluster: &SvCluster) -> bool {
    if cluster.resolved_type() == ResolvedType::LINE {
        return true;
    }
    if is_simple_single_sv(cluster) && cluster.sv0().nearest_sv_relation() == RELATION_TYPE_NEIGHBOUR {
        return true;
    }
    cluster.is_resolved()
}

pub struct SimpleClustering {
    proximity_distance: i32,
    sample_id: String,
}

impl SimpleClustering {
    pub fn new(proximity_distance: i32) -> Self {
        SimpleClustering { proximity_distance, sample_id: String::new() }
    }
    pub fn initialise(&mut self, sample_id: &str) {
        self.sample_id = sample_id.to_string();
    }

    pub fn cluster_by_proximity(&self, clusters: &mut Vec<Cluster>, state: &ClusteringState) {
        for (_, list) in state.chr_breakend_map.iter() {
            let mut current = 0;
            while current < list.len() {
                let b = &list[current];
                let var = b.sv().clone();
                let next_index = current + 1;
                let Some(next) = list.get(next_index) else {
                    if var.cluster().is_none() {
                        let c = Cluster::create(state.next_cluster_id());
                        c.add_variant(&var);
                        clusters.push(c);
                    }
                    break;
                };
                let cluster = var.cluster();
                let next_var = next.sv().clone();
                let next_cluster = next_var.cluster();
                if cluster.is_some() && cluster == next_cluster {
                    // already clustered
                } else if (next.position() - b.position()).abs() > self.proximity_distance {
                    if cluster.is_none() {
                        let c = Cluster::create(state.next_cluster_id());
                        c.add_variant(&var);
                        clusters.push(c);
                    }
                } else if var == next_var {
                    if cluster.is_none() {
                        let c = Cluster::create(state.next_cluster_id());
                        c.add_variant(&var);
                        clusters.push(c);
                    }
                } else {
                    match (&cluster, &next_cluster) {
                        (None, None) => {
                            let c = Cluster::create(state.next_cluster_id());
                            c.add_variant(&var);
                            c.add_variant(&next_var);
                            c.add_cluster_reason(PROXIMITY);
                            clusters.push(c);
                        }
                        (Some(c), Some(nc)) => {
                            c.merge_other_cluster_log(nc, false);
                            c.add_cluster_reason(PROXIMITY);
                            remove_item(clusters, nc);
                        }
                        (None, Some(nc)) => {
                            nc.add_variant(&var);
                            nc.add_cluster_reason(PROXIMITY);
                        }
                        (Some(c), None) => {
                            c.add_variant(&next_var);
                            c.add_cluster_reason(PROXIMITY);
                        }
                    }
                    if !var.has_cluster_reason(PROXIMITY) {
                        var.add_cluster_reason(PROXIMITY, next_var.id());
                    }
                    if !next_var.has_cluster_reason(PROXIMITY) {
                        next_var.add_cluster_reason(PROXIMITY, var.id());
                    }
                }
                current = next_index;
            }
        }
    }

    pub fn merge_clusters(&self, clusters: &mut Vec<Cluster>, state: &ClusteringState) {
        for cluster in clusters.iter() {
            if cluster.is_resolved() {
                continue;
            }
            mark_cluster_long_del_dups(cluster, state);
            mark_cluster_inversions(cluster);
        }
        let init = clusters.len();
        self.merge_on_loh_events(clusters, state);
        merge_on_high_facing_jcn(clusters, state);
        merge_on_major_allele_jcn_bounds(clusters, state);
        merge_loh_resolving_clusters(clusters, state);
        let mut found = true;
        let mut iterations = 0;
        while found {
            found = merge_on_overlapping_inv_dup_dels(clusters, false);
            iterations += 1;
            if iterations >= 10 {
                if found {
                    hmf_common::warn!("sample({}) exiting simple merge loop after {} iterations with merge just found", self.sample_id, iterations);
                }
                break;
            }
        }
        if clusters.len() < init {
            hmf_common::debug!("reduced cluster count({} -> {}) iterations({})", init, clusters.len(), iterations);
        }
    }

    pub fn merge_long_del_dup_clusters(&self, clusters: &mut Vec<Cluster>) {
        let mut found = true;
        let mut iterations = 0;
        while found {
            found = merge_on_overlapping_inv_dup_dels(clusters, true);
            iterations += 1;
            if iterations >= 5 {
                if found {
                    hmf_common::warn!("sample({}) exiting simple merge loop after {} iterations with merge just found", self.sample_id, iterations);
                }
                break;
            }
        }
    }

    fn merge_on_loh_events(&self, clusters: &mut Vec<Cluster>, state: &ClusteringState) {
        if state.loh_event_list.is_empty() && state.hom_loss_list.is_empty() {
            return;
        }
        for loh in &state.loh_event_list {
            if !loh.matched_both_svs() || loh.has_incomplete_hom_loss_events() {
                continue;
            }
            let bs = loh.breakend(true).unwrap();
            let be = loh.breakend(false).unwrap();
            let cs = bs.cl();
            let ce = be.cl();
            if cs != ce {
                if cs.has_linking_line_elements() || ce.has_linking_line_elements() {
                    continue;
                }
                if breakends_have_different_jcn(&bs, &be) {
                    continue;
                }
                add_cluster_reasons(bs.sv(), be.sv(), LOH);
                cs.add_cluster_reason(LOH);
                cs.merge_other_cluster(&ce);
                remove_item(clusters, &ce);
            }
        }
        for loh in &state.loh_event_list {
            if !loh.has_incomplete_hom_loss_events() {
                continue;
            }
            if loh.clustered() || loh.whole_arm_loss() {
                let unclustered: Vec<HomLoss> = loh
                    .hom_loss_events()
                    .into_iter()
                    .filter(|x| x.matched_both_svs())
                    .filter(|x| !x.same_sv())
                    .filter(|x| x.pos_start > loh.pos_start)
                    .filter(|x| x.pos_end < loh.pos_end)
                    .filter(|x| !x.clustered())
                    .collect();
                for hl in &unclustered {
                    let hs = hl.breakend(true).unwrap();
                    let he = hl.breakend(false).unwrap();
                    let c = hs.cl();
                    let oc = he.cl();
                    if c == oc {
                        continue;
                    }
                    if breakends_have_different_jcn(&hs, &he) {
                        continue;
                    }
                    add_cluster_reasons(hs.sv(), he.sv(), HOM_LOSS);
                    c.add_cluster_reason(HOM_LOSS);
                    c.merge_other_cluster(&oc);
                    remove_item(clusters, &oc);
                }
                continue;
            }
            if !loh.matched_both_svs() || loh.clustered() {
                continue;
            }
            if loh.breakend(true).unwrap().cl().has_linking_line_elements() || loh.breakend(false).unwrap().cl().has_linking_line_elements() {
                continue;
            }
            let mut incomplete = false;
            for hl in loh.hom_loss_events() {
                if !(hl.pos_start > loh.pos_start && hl.pos_end < loh.pos_end) || !hl.matched_both_svs() || !hl.clustered() {
                    incomplete = true;
                    break;
                }
            }
            if !incomplete {
                let bs = loh.breakend(true).unwrap();
                let be = loh.breakend(false).unwrap();
                if breakends_have_different_jcn(&bs, &be) {
                    continue;
                }
                let c = bs.cl();
                let oc = be.cl();
                loh.set_is_valid(true);
                add_cluster_reasons(bs.sv(), be.sv(), HOM_LOSS);
                c.add_cluster_reason(HOM_LOSS);
                c.merge_other_cluster(&oc);
                remove_item(clusters, &oc);
            }
            let mut unclustered: Vec<Be> = Vec::new();
            let mut all_valid = true;
            for s in [SE_START, SE_END] {
                let is = is_start(s);
                unclustered.push(loh.breakend(is).unwrap());
                for hl in loh.hom_loss_events() {
                    if !hl.matched_both_svs() {
                        all_valid = false;
                        break;
                    }
                    unclustered.push(hl.breakend(is).unwrap());
                }
                if !all_valid {
                    break;
                }
            }
            if all_valid {
                let mut i = 0;
                while i < unclustered.len() {
                    let b1 = unclustered[i].clone();
                    let mut found = false;
                    for j in i + 1..unclustered.len() {
                        if b1.cluster() == unclustered[j].cluster() {
                            unclustered.remove(j);
                            unclustered.remove(i);
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        i += 1;
                    }
                }
            }
            if unclustered.len() == 2 {
                let (b1, b2) = (unclustered[0].clone(), unclustered[1].clone());
                let c = b1.cl();
                let oc = b2.cl();
                loh.set_is_valid(true);
                if c == oc {
                    continue;
                }
                if breakends_have_different_jcn(&b1, &b2) {
                    continue;
                }
                add_cluster_reasons(b1.sv(), b2.sv(), HOM_LOSS);
                c.add_cluster_reason(HOM_LOSS);
                c.merge_other_cluster(&oc);
                remove_item(clusters, &oc);
            }
        }
    }
}

fn mark_cluster_inversions(cluster: &Cluster) {
    if cluster.type_count(INV) == 0 {
        return;
    }
    if cluster.is_resolved() && cluster.resolved_type().is_simple() {
        return;
    }
    for var in cluster.svs() {
        if var.sv_type() == INV && !var.is_cross_arm() {
            cluster.register_inversion(&var);
        }
    }
}

fn mark_cluster_long_del_dups(cluster: &Cluster, state: &ClusteringState) {
    if cluster.is_synthetic_type() && cluster.resolved_type().is_simple() && !cluster.is_resolved() {
        let syn = get_synthetic_length(cluster);
        if (cluster.resolved_type() == ResolvedType::DEL && syn >= state.del_cutoff_length)
            || (cluster.resolved_type() == ResolvedType::DUP && syn >= state.dup_cutoff_length)
        {
            for var in cluster.svs() {
                cluster.register_long_del_dup(&var);
            }
        }
    }
    if cluster.type_count(DEL) > 0 || cluster.type_count(DUP) > 0 {
        for var in cluster.svs() {
            if var.is_cross_arm() {
                continue;
            }
            if exceeds_dup_del_cutoff_length(var.sv_type(), var.length(), state) {
                cluster.register_long_del_dup(&var);
            }
        }
    }
}

fn merge_on_overlapping_inv_dup_dels(clusters: &mut Vec<Cluster>, allow_del_dup_overlaps: bool) -> bool {
    let ddi: Vec<Cluster> = clusters
        .iter()
        .filter(|x| !x.inversions().is_empty() || !x.long_del_dups().is_empty())
        .filter(|x| !x.has_linking_line_elements())
        .cloned()
        .collect();
    if ddi.len() <= 1 {
        return false;
    }
    let mut merged: Vec<Cluster> = Vec::new();
    let mut index1 = 0;
    while index1 < ddi.len() {
        let c1 = &ddi[index1];
        if contains(&merged, c1) {
            index1 += 1;
            continue;
        }
        let mut merged_other = false;
        let mut reason = ClusteringReason::NONE;
        let mut c1_svs = c1.long_del_dups();
        c1_svs.extend(c1.inversions());
        let mut index2 = index1 + 1;
        while index2 < ddi.len() {
            let c2 = &ddi[index2];
            if contains(&merged, c2) {
                index2 += 1;
                continue;
            }
            let mut c2_svs = c2.long_del_dups();
            c2_svs.extend(c2.inversions());
            let mut del_dup_overlap_count = 0;
            let mut close_link_pairs = 0;
            'outer: for var1 in &c1_svs {
                for var2 in &c2_svs {
                    let contains_inv = var1.sv_type() == INV || var2.sv_type() == INV;
                    if !allow_del_dup_overlaps && !contains_inv {
                        continue;
                    }
                    if var1.chromosome(true) != var2.chromosome(true) {
                        continue;
                    }
                    if var1.position(false) < var2.position(true) || var1.position(true) > var2.position(false) {
                        continue;
                    }
                    let enclosed = (var1.position(true) < var2.position(true) && var1.position(false) > var2.position(false))
                        || (var2.position(true) < var1.position(true) && var2.position(false) > var1.position(false));
                    if allow_del_dup_overlaps && enclosed {
                        continue;
                    }
                    if variants_violate_loh_hom_loss(var1, var2) {
                        continue;
                    }
                    if variants_have_different_jcn(var1, var2) {
                        continue;
                    }
                    if allow_del_dup_overlaps
                        && (!copy_numbers_equal(var1.copy_number(true), var2.copy_number(true))
                            || !copy_numbers_equal(var1.copy_number(false), var2.copy_number(false)))
                    {
                        continue;
                    }
                    if !contains_inv {
                        del_dup_overlap_count += 1;
                    }
                    let cl = breakends_in_close_link(var1, var2);
                    if cl[0] {
                        close_link_pairs += 1;
                    }
                    let close_link = cl[1];
                    if (close_link && contains_inv) || (del_dup_overlap_count >= 3 && close_link_pairs >= 3) {
                        reason = LONG_DEL_DUP_INV;
                        add_cluster_reasons(var1, var2, reason);
                        merged_other = true;
                        break 'outer;
                    }
                }
            }
            if merged_other {
                c1.merge_other_cluster(c2);
                c1.add_cluster_reason(reason);
                merged.push(c2.clone());
                break;
            }
            index2 += 1;
        }
        if !merged_other {
            index1 += 1;
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(clusters, x);
    }
    true
}

fn merge_on_major_allele_jcn_bounds(clusters: &mut Vec<Cluster>, state: &ClusteringState) -> bool {
    let mut merged: Vec<Cluster> = Vec::new();
    let mut ci = 0;
    while ci < clusters.len() {
        let cluster = clusters[ci].clone();
        if contains(&merged, &cluster) || skip_cluster_type(&cluster) {
            ci += 1;
            continue;
        }
        let mut merged_other = false;
        'chr: for (chr, list) in cluster.chr_breakend_entries() {
            let full = state.breakend_list(&chr).unwrap();
            for i in 0..=1 {
                let traverse_up = i == 0;
                let mut index: i32 = if traverse_up { -1 } else { list.len() as i32 };
                loop {
                    index += if traverse_up { 1 } else { -1 };
                    if index < 0 || index as usize >= list.len() {
                        break;
                    }
                    let b = &list[index as usize];
                    if b.reverse_orient() != traverse_up {
                        continue;
                    }
                    let mut be_jcn = b.jcn();
                    let mut chr_index = b.chr_pos_index();
                    let mut opposing: Vec<Be> = Vec::new();
                    loop {
                        chr_index += if traverse_up { 1 } else { -1 };
                        if chr_index < 0 || chr_index as usize >= full.len() {
                            break;
                        }
                        let nb = &full[chr_index as usize];
                        if (nb.position() - b.position()).abs() > MAX_MERGE_DISTANCE {
                            break;
                        }
                        if nb.arm() != b.arm() {
                            break;
                        }
                        if nb.orientation() == b.orientation() {
                            continue;
                        }
                        if skip_cluster_type(&nb.cl()) {
                            continue;
                        }
                        if nb.cl() == cluster {
                            be_jcn -= nb.jcn();
                            if be_jcn <= 0.0 {
                                break;
                            }
                            continue;
                        }
                        if nb.is_assembled_link() {
                            continue;
                        }
                        if breakends_have_different_jcn(b, nb) {
                            continue;
                        }
                        opposing.push(nb.clone());
                        let following_major = nb.major_allele_jcn(!traverse_up);
                        if !copy_numbers_equal(be_jcn, following_major) && be_jcn > following_major {
                            let max_opposing = opposing.iter().map(|x| x.jcn()).fold(f64::NEG_INFINITY, jmax);
                            let ob = opposing.iter().find(|x| copy_numbers_equal(x.jcn(), max_opposing)).unwrap().clone();
                            let other_cluster = ob.cl();
                            add_cluster_reasons(b.sv(), ob.sv(), MAJOR_ALLELE_JCN);
                            other_cluster.add_cluster_reason(MAJOR_ALLELE_JCN);
                            cluster.add_cluster_reason(MAJOR_ALLELE_JCN);
                            cluster.merge_other_cluster(&other_cluster);
                            merged.push(other_cluster);
                            merged_other = true;
                            break;
                        }
                    }
                    if merged_other {
                        break;
                    }
                }
                if merged_other {
                    break;
                }
            }
            if merged_other {
                break 'chr;
            }
        }
        if !merged_other {
            ci += 1;
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(clusters, x);
    }
    true
}

fn merge_on_high_facing_jcn(clusters: &mut Vec<Cluster>, state: &ClusteringState) -> bool {
    let mut merged: Vec<Cluster> = Vec::new();
    for (_, list) in state.chr_breakend_map.iter() {
        for b in list {
            if b.forward_orient() {
                continue;
            }
            if b.cl().resolved_type() == ResolvedType::LINE {
                continue;
            }
            if b.jcn() < HIGH_JCN_THRESHOLD {
                continue;
            }
            let adj = b.major_allele_jcn(true);
            if b.jcn() / jmax(adj, 0.01) < ADJACENT_JCN_RATIO {
                continue;
            }
            for nb in list.iter().skip((b.chr_pos_index() + 1) as usize) {
                if nb.reverse_orient() {
                    continue;
                }
                let next_adj = nb.major_allele_jcn(false);
                let high = nb.jcn() >= HIGH_JCN_THRESHOLD && nb.jcn() / jmax(next_adj, 0.01) >= ADJACENT_JCN_RATIO;
                if high {
                    if nb.cl() == b.cl() || nb.cl().resolved_type() == ResolvedType::LINE {
                        continue;
                    }
                    let c = b.cl();
                    let oc = nb.cl();
                    add_cluster_reasons(b.sv(), nb.sv(), HIGH_JCN);
                    oc.add_cluster_reason(HIGH_JCN);
                    c.add_cluster_reason(HIGH_JCN);
                    c.merge_other_cluster(&oc);
                    merged.push(oc);
                }
                if nb.major_allele_jcn(false) < HIGH_JCN_THRESHOLD {
                    break;
                }
            }
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(clusters, x);
    }
    true
}

fn merge_loh_resolving_clusters(clusters: &mut Vec<Cluster>, state: &ClusteringState) -> bool {
    let with_loh: Vec<Cluster> = clusters.iter().filter(|x| !x.loh_events().is_empty()).filter(|x| !x.has_linking_line_elements()).cloned().collect();
    let mut merged: Vec<Cluster> = Vec::new();
    for loh_cluster in &with_loh {
        if contains(&merged, loh_cluster) {
            continue;
        }
        let mut loh_index = 0;
        while loh_index < loh_cluster.loh_events().len() {
            let loh = loh_cluster.loh_events()[loh_index].clone();
            if !loh.is_valid() {
                loh_index += 1;
                continue;
            }
            for s in [SE_START, SE_END] {
                let Some(loh_b) = loh.breakend(is_start(s)) else { continue };
                if loh_b.sv_type() != DUP {
                    continue;
                }
                if loh_b.cluster().as_ref() != Some(loh_cluster) {
                    continue;
                }
                let full = state.breakend_list(loh_b.chromosome()).unwrap();
                let other = loh_b.other_breakend().unwrap();
                let mut index = other.chr_pos_index();
                let traverse_up = other.reverse_orient();
                let mut resolving_cluster: Option<Cluster> = None;
                let mut resolving_breakend: Option<Be> = None;
                loop {
                    index += if traverse_up { 1 } else { -1 };
                    if index < 0 || index as usize >= full.len() {
                        break;
                    }
                    let nb = &full[index as usize];
                    if nb == &loh_b {
                        let Some(rc) = &resolving_cluster else { break };
                        let rb = resolving_breakend.as_ref().unwrap();
                        add_cluster_reasons(loh_b.sv(), rb.sv(), LOH_CHAIN);
                        rc.add_cluster_reason(LOH_CHAIN);
                        loh_cluster.add_cluster_reason(LOH_CHAIN);
                        loh_cluster.merge_other_cluster(rc);
                        merged.push(rc.clone());
                        break;
                    }
                    if nb.orientation() == other.orientation() {
                        continue;
                    }
                    if nb.is_assembled_link() {
                        continue;
                    }
                    if breakends_have_different_jcn(&loh_b, nb) {
                        continue;
                    }
                    let oc = nb.cl();
                    if &oc == loh_cluster {
                        break;
                    }
                    if contains(&merged, &oc) || oc.is_resolved() {
                        continue;
                    }
                    if resolving_cluster.is_some() {
                        break;
                    }
                    resolving_breakend = Some(nb.clone());
                    resolving_cluster = Some(oc);
                }
            }
            loh_index += 1;
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(clusters, x);
    }
    true
}

// ComplexClustering

const MAX_COMPLEX_CLUSTER_MERGE_SIZE: usize = 5;

pub fn apply_complex_rules(sample_id: &str, clusters: &mut Vec<Cluster>, state: &ClusteringState, foldbacks_only: bool) {
    let mut complex: Vec<Cluster> = clusters.iter().filter(|x| !x.is_resolved()).cloned().collect();
    let mut iterations = 1;
    let mut found = true;
    while found {
        found = false;
        let mut index1 = 0;
        while index1 < complex.len() {
            let c1 = complex[index1].clone();
            if c1.is_resolved() {
                index1 += 1;
                continue;
            }
            let mut index2 = index1 + 1;
            while index2 < complex.len() {
                let c2 = complex[index2].clone();
                if c2.is_resolved() {
                    index2 += 1;
                    continue;
                }
                let mut can_merge = can_merge_clusters_on_foldbacks(&c1, &c2);
                if !can_merge && !foldbacks_only {
                    can_merge = can_merge_clusters_on_common_arms(&c1, &c2);
                }
                if !can_merge {
                    index2 += 1;
                    continue;
                }
                found = true;
                c1.merge_other_cluster_log(&c2, false);
                complex.remove(index2);
                remove_item(clusters, &c2);
            }
            index1 += 1;
        }
        if !foldbacks_only {
            if merge_breakend_straddled_clusters(&mut complex, clusters, state) {
                found = true;
            }
            if merge_facing_jcn_linked_clusters(&mut complex, clusters, state) {
                found = true;
            }
            if merge_single_satellite_repeats(&mut complex, clusters) {
                found = true;
            }
        }
        iterations += 1;
        if iterations > 20 {
            hmf_common::warn!("sample({}) reached {} iterations of clustering merging", sample_id, iterations);
            break;
        }
    }
}

fn can_merge_clusters_on_foldbacks(c1: &Cluster, c2: &Cluster) -> bool {
    let f1 = c1.foldbacks();
    let f2 = c2.foldbacks();
    if f1.is_empty() || f2.is_empty() {
        return false;
    }
    for var1 in &f1 {
        for s1 in [SE_START, SE_END] {
            let v1s = is_start(s1);
            if s1 == SE_END && var1.is_sgl_breakend() {
                continue;
            }
            if var1.foldback_breakend(v1s).is_none() {
                continue;
            }
            for var2 in &f2 {
                for s2 in [SE_START, SE_END] {
                    let v2s = is_start(s2);
                    if s2 == SE_END && var2.is_sgl_breakend() {
                        continue;
                    }
                    if var2.foldback_breakend(v2s).is_none() {
                        continue;
                    }
                    if var1.chromosome(v1s) != var2.chromosome(v2s) || var1.arm(v1s) != var2.arm(v2s) {
                        continue;
                    }
                    if variants_have_different_jcn(var1, var2) {
                        continue;
                    }
                    add_cluster_reasons(var1, var2, FOLDBACKS);
                    c1.add_cluster_reason(FOLDBACKS);
                    c2.add_cluster_reason(FOLDBACKS);
                    return true;
                }
            }
        }
    }
    false
}

fn can_merge_clusters_on_common_arms(c1: &Cluster, c2: &Cluster) -> bool {
    c1.set_arm_links();
    c2.set_arm_links();
    let l1 = c1.unlinked_remote_svs();
    let l2 = c2.unlinked_remote_svs();
    for var1 in &l1 {
        for var2 in &l2 {
            if !have_same_chr_arms(var1, var2) {
                continue;
            }
            if variants_violate_loh_hom_loss(var1, var2) {
                continue;
            }
            if variants_have_different_jcn(var1, var2) {
                continue;
            }
            add_cluster_reasons(var1, var2, COMMON_ARMS);
            c1.add_cluster_reason(COMMON_ARMS);
            c2.add_cluster_reason(COMMON_ARMS);
            return true;
        }
    }
    false
}

fn merge_breakend_straddled_clusters(complex: &mut Vec<Cluster>, all: &mut Vec<Cluster>, state: &ClusteringState) -> bool {
    let mut merged: Vec<Cluster> = Vec::new();
    let mut ci = 0;
    while ci < complex.len() {
        let cluster = complex[ci].clone();
        if contains(&merged, &cluster) || cluster.is_resolved() {
            ci += 1;
            continue;
        }
        let mut merged_other = false;
        for (chr, list) in cluster.chr_breakend_entries() {
            let full = state.breakend_list(&chr).unwrap();
            for i in 0..list.len().saturating_sub(1) {
                let lower = list[i].clone();
                let mut upper = list[i + 1].clone();
                if lower.arm() != upper.arm() {
                    continue;
                }
                let is_foldback_pair = lower.is_foldback()
                    && upper.is_foldback()
                    && lower.foldback_breakend().as_ref() != Some(&upper)
                    && !(lower.forward_orient() && upper.reverse_orient());
                if !is_foldback_pair {
                    if lower.orientation() == upper.orientation() {
                        let front = if lower.forward_orient() { &lower } else { &upper };
                        if front.db_link().is_some_and(|d| d.length() < 0) {
                            continue;
                        }
                    } else {
                        if i + 2 >= list.len() || list[i + 2].orientation() != lower.orientation() {
                            continue;
                        }
                        let back = if lower.reverse_orient() { &lower } else { &upper };
                        if back.db_link().is_some_and(|d| d.length() < 0) {
                            upper = list[i + 2].clone();
                        } else {
                            continue;
                        }
                    }
                    if !jcn_match_breakends(&lower, &upper) {
                        continue;
                    }
                }
                let ci_lower = lower.chr_pos_index();
                let ci_upper = upper.chr_pos_index();
                let mut other_b: Option<Be> = None;
                let mut unclustered = 0;
                for j in (ci_lower + 1)..ci_upper {
                    let b = &full[j as usize];
                    let oc = b.cl();
                    if oc.is_resolved() || oc == cluster {
                        continue;
                    }
                    unclustered += 1;
                    other_b = Some(b.clone());
                }
                if unclustered != 1 {
                    continue;
                }
                let ob = other_b.unwrap();
                let oc = ob.cl();
                if oc == cluster || oc.is_resolved() || contains(&merged, &oc) {
                    continue;
                }
                if cluster.sv_count() > MAX_COMPLEX_CLUSTER_MERGE_SIZE && oc.sv_count() > MAX_COMPLEX_CLUSTER_MERGE_SIZE {
                    continue;
                }
                if !is_foldback_pair && ob.orientation() == lower.orientation() {
                    continue;
                }
                let reason = if is_foldback_pair { OVERLAP_FOLDBACKS } else { CONSEC_BREAKS };
                add_cluster_reasons(ob.sv(), lower.sv(), reason);
                oc.add_cluster_reason(reason);
                cluster.add_cluster_reason(reason);
                cluster.merge_other_cluster(&oc);
                merged.push(oc);
                merged_other = true;
                break;
            }
            if merged_other {
                break;
            }
        }
        if !merged_other {
            ci += 1;
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(complex, x);
    }
    for x in &merged {
        remove_item(all, x);
    }
    true
}

fn merge_facing_jcn_linked_clusters(complex: &mut Vec<Cluster>, all: &mut Vec<Cluster>, state: &ClusteringState) -> bool {
    let mut merged: Vec<Cluster> = Vec::new();
    let mut ci = 0;
    while ci < complex.len() {
        let cluster = complex[ci].clone();
        if contains(&merged, &cluster) || cluster.is_resolved() {
            ci += 1;
            continue;
        }
        let mut merged_other = false;
        for (chr, list) in cluster.chr_breakend_entries() {
            let full = state.breakend_list(&chr).unwrap();
            for i in 0..=1 {
                let traverse_up = i == 0;
                let boundary = if traverse_up { list[list.len() - 1].clone() } else { list[0].clone() };
                if traverse_up == boundary.forward_orient() {
                    continue;
                }
                let mut index = boundary.chr_pos_index();
                loop {
                    index += if traverse_up { 1 } else { -1 };
                    if index < 0 || index as usize >= full.len() {
                        break;
                    }
                    let nb = &full[index as usize];
                    if breakends_have_different_jcn(&boundary, nb) {
                        continue;
                    }
                    if (nb.position() - boundary.position()).abs() > MAX_MERGE_DISTANCE {
                        break;
                    }
                    let oc = nb.cl();
                    if oc == cluster || skip_cluster_type(&oc) {
                        continue;
                    }
                    if cluster.sv_count() > MAX_COMPLEX_CLUSTER_MERGE_SIZE && oc.sv_count() > MAX_COMPLEX_CLUSTER_MERGE_SIZE {
                        continue;
                    }
                    if nb.orientation() == boundary.orientation() {
                        break;
                    }
                    if copy_numbers_equal(boundary.jcn(), nb.jcn()) {
                        add_cluster_reasons(boundary.sv(), nb.sv(), TI_JCN_MATCH);
                        oc.add_cluster_reason(TI_JCN_MATCH);
                        cluster.add_cluster_reason(TI_JCN_MATCH);
                        cluster.merge_other_cluster(&oc);
                        merged.push(oc);
                        merged_other = true;
                    }
                    break;
                }
                if merged_other {
                    break;
                }
            }
            if merged_other {
                break;
            }
        }
        if !merged_other {
            ci += 1;
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(complex, x);
    }
    for x in &merged {
        remove_item(all, x);
    }
    true
}

fn merge_single_satellite_repeats(complex: &mut Vec<Cluster>, all: &mut Vec<Cluster>) -> bool {
    let with_sr: Vec<Cluster> = complex.iter().filter(|x| x.svs().iter().any(|y| y.sgl_to_satellite_repeats())).cloned().collect();
    if with_sr.is_empty() {
        return false;
    }
    let mut one_or_no_sgls: Vec<Cluster> = with_sr.iter().filter(|x| x.type_count(SGL) <= 1).cloned().collect();
    let mut merged: Vec<Cluster> = Vec::new();
    for sr in &with_sr {
        if contains(&merged, sr) {
            continue;
        }
        let arms1: Vec<String> =
            sr.svs().iter().filter(|x| x.sgl_to_satellite_repeats()).map(|x| make_chr_arm_str(x.chromosome(true), x.arm(true))).collect();
        let mut index = 0;
        while index < one_or_no_sgls.len() {
            let sgl_cluster = one_or_no_sgls[index].clone();
            if &sgl_cluster == sr {
                index += 1;
                continue;
            }
            let mut did_merge = false;
            for var in sgl_cluster.svs() {
                if var.sgl_to_satellite_repeats() && arms1.contains(&make_chr_arm_str(var.chromosome(true), var.arm(true))) {
                    let chromosome = var.chromosome(true).to_string();
                    let other = sr.svs().into_iter().find(|x| x.sgl_to_satellite_repeats() && x.chromosome(true) == chromosome).unwrap();
                    if variants_have_different_jcn(&var, &other) {
                        continue;
                    }
                    add_cluster_reasons(&other, &var, SATELLITE_SGL);
                    sr.add_cluster_reason(SATELLITE_SGL);
                    sr.merge_other_cluster(&sgl_cluster);
                    merged.push(sgl_cluster.clone());
                    did_merge = true;
                    break;
                }
            }
            if did_merge {
                one_or_no_sgls.remove(index);
            } else {
                index += 1;
            }
        }
    }
    if merged.is_empty() {
        return false;
    }
    for x in &merged {
        remove_item(complex, x);
    }
    for x in &merged {
        remove_item(all, x);
    }
    true
}

// FoldbackFinder

pub fn mark_foldbacks(map: &JavaHashMap<String, Vec<Be>>, recheck_chained_svs: bool) -> bool {
    let mut found_any = false;
    for (_, list) in map.iter() {
        if list.len() < 2 {
            continue;
        }
        for i in 0..list.len() - 1 {
            let b = &list[i];
            if b.is_assembled_link() {
                continue;
            }
            if recheck_chained_svs && (b.is_foldback() || b.cl().sv_count() == 1) {
                continue;
            }
            let mut be_front: Option<Be> = None;
            let mut be_back: Option<Be> = None;
            let mut j = i + 1;
            while j < list.len() {
                let mut nb = list[j].clone();
                if j == i + 1 && b.reverse_orient() && nb.forward_orient() && breakends_cannot_link(b, &nb) {
                    j += 1;
                    continue;
                }
                if j + 1 < list.len() && nb.is_assembled_link() {
                    if let Some(asmb) = nb.sv().linked_pair(nb.uses_start()) {
                        let nn = &list[j + 1];
                        if &asmb.other_breakend(&nb) == nn {
                            j += 2;
                            continue;
                        }
                    }
                }
                if b.forward_orient() && nb.reverse_orient() && j < list.len() - 1 {
                    let nn = list[j + 1].clone();
                    if nn.orientation() == b.orientation() && nn.position() - nb.position() < get_min_templated_insertion_length(&nb, &nn) {
                        nb = nn;
                    }
                }
                if nb.orientation() == b.orientation() && !nb.is_assembled_link() {
                    be_front = Some(if b.forward_orient() { b.clone() } else { nb.clone() });
                    be_back = Some(if b.forward_orient() { nb.clone() } else { b.clone() });
                }
                break;
            }
            let mut found = false;
            if let (Some(f), Some(bk)) = (&be_front, &be_back) {
                let db = f.sv().db_link(f.uses_start());
                if db.is_none_or(|d| d.length() > 0)
                    && (!recheck_chained_svs || (f.cluster() == bk.cluster() && !f.is_foldback() && !bk.is_foldback()))
                {
                    found = check_foldback_breakends(f, bk);
                }
            }
            if !found && !recheck_chained_svs {
                check_replicated_breakend_foldback(b);
            }
            if found {
                found_any = true;
            }
        }
    }
    found_any
}

fn check_foldback_breakends(be_start: &Be, be_end: &Be) -> bool {
    let var_end = be_end.sv().clone();
    let var_start = be_start.sv().clone();
    if var_end.sv_type() == INS || var_start.sv_type() == INS {
        return false;
    }
    let cluster = var_end.cl();
    let single_sv = var_end == var_start;
    if single_sv {
        if var_start.sv_type() != INV {
            return false;
        }
    } else if var_start.cluster().as_ref() != Some(&cluster) {
        return false;
    }
    if cluster.is_resolved() {
        return false;
    }
    let be_end_uses_start = be_end.uses_start();
    let be_start_uses_start = be_start.uses_start();
    let mut chain_info = "0;0;0".to_string();
    if single_sv {
        if let Some(chain) = cluster.find_chain(&var_end) {
            let mut bnd_links = 0;
            for pair in chain.linked_pairs() {
                if pair.first() == var_end && pair.second().sv_type() == BND {
                    bnd_links += 1;
                } else if pair.second() == var_end && pair.first().sv_type() == BND {
                    bnd_links += 1;
                }
                if bnd_links == 2 {
                    return false;
                }
            }
        }
    } else {
        let c1 = cluster.find_chain(&var_end);
        let c2 = cluster.find_chain(&var_start);
        let (Some(c1), Some(c2)) = (c1, c2) else { return false };
        if c1 != c2 {
            return false;
        }
        let cd = breakends_are_chained(&c1, &var_end, !be_end_uses_start, &var_start, !be_start_uses_start);
        if cd[CHAIN_LINK_COUNT] == 0 {
            return false;
        }
        if cd[CHAIN_LENGTH] > MAX_FOLDBACK_CHAIN_LENGTH {
            return false;
        }
        chain_info = format!("{};{};{}", cd[CHAIN_LINK_COUNT], cd[CHAIN_ASSEMBLY_LINK_COUNT], cd[CHAIN_LENGTH]);
    }
    let length = (be_end.position() - be_start.position()).abs();
    if var_end.foldback_breakend(be_end_uses_start).is_some() && !single_sv && var_end.foldback_length(be_end_uses_start) < length {
        return false;
    } else if var_start.foldback_breakend(be_start_uses_start).is_some() && !single_sv && var_start.foldback_length(be_start_uses_start) < length {
        return false;
    }
    if let Some(o) = var_end.foldback_breakend(be_end_uses_start) {
        o.sv().set_foldback_link(o.uses_start(), None, -1, "");
    }
    if let Some(o) = var_start.foldback_breakend(be_start_uses_start) {
        o.sv().set_foldback_link(o.uses_start(), None, -1, "");
    }
    var_end.set_foldback_link(be_end_uses_start, Some(be_start.clone()), length, &chain_info);
    var_start.set_foldback_link(be_start_uses_start, Some(be_end.clone()), length, &chain_info);
    if single_sv {
        hmf_common::debug!("cluster({}) foldback inversion SV({}) length({})", cluster.id(), var_end.pos_id(), length);
    } else {
        hmf_common::debug!("cluster({}) foldback be1({}) be2({}) length({})", cluster.id(), be_end, be_start, length);
    }
    true
}

fn check_replicated_breakend_foldback(b: &Be) {
    let var = b.sv();
    if var.assembled_linked_pairs(!b.uses_start()).len() < 2 {
        return;
    }
    let Some(chain) = var.cl().find_chain(var) else { return };
    if chain.open_breakend(true).as_ref() == Some(b) && chain.open_breakend(false).as_ref() == Some(b) {
        let chain_info = format!("{};{};{}", 1, 1, chain.length(false));
        var.set_foldback_link(b.uses_start(), Some(b.clone()), 0, &chain_info);
    }
}

// ArmCluster

pub fn build_arm_clusters(cluster: &Cluster) {
    let entries = cluster.chr_breakend_entries();
    let mut acs = std::mem::take(&mut *cluster.arm_clusters.borrow_mut());
    let find = |acs: &Vec<ArmCluster>, b: &Be| acs.iter().position(|ac| contains(&ac.breakends, b));
    for (_, list) in &entries {
        let mut prev: Option<usize> = None;
        for b in list {
            let var = b.sv();
            if var.is_foldback() {
                if let Some(ofb) = var.foldback_breakend(b.uses_start()) {
                    if let Some(i) = find(&acs, &ofb) {
                        acs[i].add_breakend(b);
                        continue;
                    }
                }
            }
            if let Some(p) = prev {
                if b.arm() == acs[p].arm && b.position() - acs[p].end_pos <= DEFAULT_PROXIMITY_DISTANCE {
                    acs[p].add_breakend(b);
                    continue;
                }
            }
            let mut found = false;
            for (ai, ac) in acs.iter_mut().enumerate() {
                if b.chromosome() != ac.chromosome || b.arm() != ac.arm {
                    continue;
                }
                if b.position() >= ac.start_pos - DEFAULT_PROXIMITY_DISTANCE && b.position() <= ac.end_pos + DEFAULT_PROXIMITY_DISTANCE {
                    ac.add_breakend(b);
                    found = true;
                    prev = Some(ai);
                    break;
                }
            }
            if !found {
                let mut ac = ArmCluster::new(acs.len() as i32, cluster.id(), b.chromosome(), b.arm());
                ac.add_breakend(b);
                acs.push(ac);
                prev = Some(acs.len() - 1);
            }
        }
    }
    for ac in acs.iter_mut() {
        set_arm_cluster_features(ac);
    }
    *cluster.arm_clusters.borrow_mut() = acs;
}

fn set_arm_cluster_features(ac: &mut ArmCluster) {
    use ArmClusterType::*;
    ac.ti_count = 0;
    let bes = ac.breakends.clone();
    if bes.len() == 1 {
        ac.arm_cluster_type = ISOLATED_BE;
        return;
    }
    if bes.len() == 2 {
        let (be1, be2) = (&bes[0], &bes[1]);
        let (var1, var2) = (be1.sv(), be2.sv());
        if var1 == var2 {
            if var1.sv_type() == DEL {
                ac.arm_cluster_type = DSB;
                return;
            } else if var1.sv_type() == DUP {
                ac.arm_cluster_type = SIMPLE_DUP;
                return;
            }
        }
        if var1.foldback_id(be1.uses_start()) == var2.id() {
            ac.arm_cluster_type = FOLDBACK;
            return;
        }
        if let Some(ti) = var1.linked_pair(be1.uses_start()) {
            if ti.has_breakend(be2) {
                ac.arm_cluster_type = TI_ONLY;
                ac.ti_count = 1;
                return;
            }
        }
        let db = be1.sv().db_link(be1.uses_start());
        if db.is_some() && db == be2.sv().db_link(be2.uses_start()) {
            ac.arm_cluster_type = DSB;
            return;
        }
    }
    let mut dsb_count = 0;
    let mut foldback_count = 0;
    let mut suspect_line = 0;
    let mut ti_count = 0;
    let mut ti_breakends: Vec<Be> = Vec::new();
    for b in &bes {
        if b.sv().is_line_element(b.uses_start()) {
            suspect_line += 1;
        }
        if b.sv().is_foldback() {
            foldback_count += 1;
        }
        if b.reverse_orient() {
            if let Some(ti) = b.sv().linked_pair(b.uses_start()) {
                let o = ti.other_breakend(b);
                if contains(&bes, &o) {
                    ti_breakends.push(b.clone());
                    ti_breakends.push(o);
                    ti_count += 1;
                }
            }
        }
    }
    // Java: (foldbackCount % 2) == 2 is never true
    foldback_count /= 2;
    ac.ti_count = ti_count;
    if ti_breakends.len() == bes.len() {
        ac.arm_cluster_type = TI_ONLY;
        return;
    }
    if suspect_line > 0 {
        ac.arm_cluster_type = COMPLEX_LINE;
        return;
    } else if foldback_count > 1 {
        ac.arm_cluster_type = COMPLEX_FOLDBACK;
        return;
    }
    let non_ti_count = bes.len() as i32 - ti_count * 2;
    let unlinked: Vec<Be> = bes.iter().filter(|x| !contains(&ti_breakends, x)).cloned().collect();
    for i in 0..unlinked.len().saturating_sub(1) {
        let (b1, b2) = (&unlinked[i], &unlinked[i + 1]);
        if b1.forward_orient() && b2.reverse_orient() {
            dsb_count += 1;
        } else if b1.reverse_orient() && b2.forward_orient() && breakends_cannot_link(b1, b2) {
            dsb_count += 1;
        }
    }
    ac.arm_cluster_type = if non_ti_count == 1 {
        ISOLATED_BE
    } else if foldback_count == 1 && dsb_count == 0 {
        FOLDBACK
    } else if non_ti_count == 3 && foldback_count == 1 && dsb_count == 1 {
        FOLDBACK_DSB
    } else if foldback_count == 0 && dsb_count == 1 && non_ti_count == 2 {
        DSB
    } else if unlinked.len() == 2 && unlinked[0].orientation() == unlinked[1].orientation() {
        SAME_ORIENT
    } else {
        COMPLEX_OTHER
    };
}
