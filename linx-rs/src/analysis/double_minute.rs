//! DoubleMinuteFinder and DoubleMinuteData: high-JCN clusters forming closed (ecDNA) chains. The
//! DOUBLE_MINUTES annotation extension's ecDNA file (reportCluster) is not ported.

use std::collections::HashSet;

use hmf_common::java_hash::JavaHashMap;
use hmf_common::regions::ChrBaseRegion;
use hmf_common::sv::StructuralVariantType::*;

use crate::analysis::classification::classify_single_pair_resolved_type;
use crate::chaining::finder::{ChainFinder, LR_METHOD_DM_CLOSE};
use crate::types::*;
use crate::utils::*;

pub const JCN_UPPER_THRESHOLD: f64 = 8.0;
const JCN_THRESHOLD: f64 = 5.0;
const LOWER_ADJACENT_JCN_RATIO: f64 = 2.0;
const MIN_PERC_OF_MAX_JCN: f64 = 0.25;
pub const MIN_SEGMENT_DEPTH_WINDOW_COUNT: i32 = 6;

pub struct DoubleMinuteFinder {
    chain_finder: ChainFinder,
    /// getDoubleMinutes: by cluster id, the DMs (and every candidate when logging)
    pub double_minutes: JavaHashMap<i32, DoubleMinuteData>,
    processed_clusters: Vec<i32>,
    show_candidates: bool,
    log_candidates: bool,
}

impl DoubleMinuteFinder {
    pub fn new(log_candidates: bool, show_candidates: bool) -> Self {
        DoubleMinuteFinder { chain_finder: ChainFinder::new(), double_minutes: JavaHashMap::new(), processed_clusters: Vec::new(), show_candidates, log_candidates }
    }
    pub fn clear(&mut self) {
        self.double_minutes.clear();
        self.processed_clusters.clear();
    }

    pub fn analyse_cluster(&mut self, cluster: &Cluster, reassess: bool, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
        if self.processed_clusters.contains(&cluster.id()) {
            if !reassess {
                return;
            }
            cluster.clear_double_minute_data();
        } else {
            self.processed_clusters.push(cluster.id());
        }
        let mut has_valid_sv = false;
        let mut candidates: Vec<Sv> = Vec::new();
        for var in cluster.svs() {
            if var.jcn() < JCN_THRESHOLD {
                continue;
            }
            let r = get_adjacent_major_ap_ratio(&var);
            if r >= ADJACENT_JCN_RATIO {
                has_valid_sv = true;
            }
            if r >= LOWER_ADJACENT_JCN_RATIO {
                candidates.push(var.clone());
            }
        }
        if !has_valid_sv {
            return;
        }
        let max_dm_jcn = candidates.iter().map(|x| x.jcn()).fold(None, |a: Option<f64>, v| Some(a.map_or(v, |a| jmax(a, v)))).unwrap_or(0.0);
        let min_jcn = jmax(JCN_THRESHOLD, MIN_PERC_OF_MAX_JCN * max_dm_jcn);
        let dm_svs: Vec<Sv> = candidates.iter().filter(|x| x.jcn() >= min_jcn).cloned().collect();
        if dm_svs.len() == 1 {
            let var = &dm_svs[0];
            if var.sv_type() != DUP {
                return;
            }
            if !variant_exceeds_both_adjacent_jcn(var) {
                return;
            }
        }
        let mut dm_chains = self.create_dm_chains(cluster, &dm_svs, false);
        if dm_svs.len() > 2 && cluster.requires_replication() {
            let varied = self.create_dm_chains(cluster, &dm_svs, true);
            let uniform_chained = count_chained_svs(&dm_svs, &dm_chains);
            let uniform_closed = count_closeable_chains(&dm_chains);
            let varied_chained = count_chained_svs(&dm_svs, &varied);
            let varied_closed = count_closeable_chains(&varied);
            if varied_closed == varied.len() && uniform_closed < dm_chains.len() && varied_chained >= uniform_chained {
                dm_chains = varied;
            }
        }
        let chained_count = count_chained_svs(&dm_svs, &dm_chains);
        let unchained_count = dm_svs.len() - chained_count;
        for c in &dm_chains {
            if !c.is_closed_loop() && c.could_close_chain() {
                let max_index = c.linked_pairs().iter().map(|x| x.link_index()).max().unwrap_or(0);
                c.close_chain_with(LR_METHOD_DM_CLOSE, max_index + 1);
            }
        }
        let mut fully_chained = !dm_chains.is_empty();
        if unchained_count > 0 {
            fully_chained = false;
        }
        for c in &dm_chains {
            if c.is_closed_loop() {
                continue;
            }
            if !c.chain_end_sv(true).unwrap().is_sgl_breakend() && !c.chain_end_sv(false).unwrap().is_sgl_breakend() {
                fully_chained = false;
                break;
            }
        }
        let mut dm = DoubleMinuteData::new(cluster, dm_svs.clone());
        dm.chains = dm_chains;
        dm.fully_chained = fully_chained;
        dm.annotate(chr_breakend_map);
        let is_dm = dm.is_double_minute;
        if !self.log_candidates && !self.show_candidates && !is_dm {
            return;
        }
        hmf_common::debug!(
            "cluster({}) dmSVs({}) chains({}) unchainedSVs({}) {}",
            cluster.id(),
            dm_svs.len(),
            dm.chains.len(),
            unchained_count,
            if fully_chained { "fully chained" } else { "invalid chain" }
        );
        if is_dm {
            cluster.set_double_minute_data(dm.valid_svs.clone(), dm.valid_chains.clone());
            for c in &dm.valid_chains {
                add_single_dup_chain(cluster, c);
                c.set_double_minute(true);
            }
            if cluster.sv_count() <= 2 && cluster.sv_count() == dm_svs.len() {
                cluster.set_resolved(false, ResolvedType::DOUBLE_MINUTE);
            }
            cluster.add_annotation(CLUSTER_ANNOT_DM);
        } else if self.show_candidates {
            cluster.set_double_minute_data(dm.svs.clone(), dm.chains.clone());
            for c in &dm.chains {
                add_single_dup_chain(cluster, c);
                c.set_double_minute(true);
            }
        }
        self.double_minutes.insert(cluster.id(), dm);
    }

    fn create_dm_chains(&mut self, cluster: &Cluster, dm_svs: &[Sv], apply_replication: bool) -> Vec<Chain> {
        let mut dm_chains: Vec<Chain> = Vec::new();
        let mut chain_id = 0;
        let mut chain_sv_list: Vec<Sv> = Vec::new();
        for var in dm_svs {
            if var.sv_type() != DUP {
                chain_sv_list.push(var.clone());
                continue;
            }
            let has_overlap = dm_svs.iter().filter(|x| *x != var).any(|x| {
                position_within(x.position(true), var.position(true), var.position(false))
                    || position_within(x.position(false), var.position(true), var.position(false))
            });
            if !has_overlap && variant_exceeds_both_adjacent_jcn(var) {
                let chain = Chain::create(chain_id);
                chain_id += 1;
                let pair = Pair::from_svs(var, var, true, false);
                chain.add_link(&pair, true);
                chain.set_jcn_data(var.jcn(), var.jcn_uncertainty());
                chain.close_chain();
                chain.set_double_minute(true);
                dm_chains.push(chain);
            } else {
                chain_sv_list.push(var.clone());
            }
        }
        self.chain_finder.initialise_subset(cluster, &chain_sv_list, apply_replication);
        self.chain_finder.form_chains(false);
        for c in self.chain_finder.unique_chains() {
            c.set_id(chain_id);
            chain_id += 1;
            dm_chains.push(c.clone());
        }
        self.chain_finder.clear();
        dm_chains
    }
}

fn add_single_dup_chain(cluster: &Cluster, c: &Chain) {
    if c.sv_count() == 1 && c.sv_list()[0].sv_type() == DUP {
        let dup = c.sv_list()[0].clone();
        if !cluster.chains().iter().any(|x| x.sv_count() == 1 && x.has_sv(&dup)) {
            cluster.add_chain(c, false);
        }
    }
}

fn count_chained_svs(svs: &[Sv], chains: &[Chain]) -> usize {
    svs.iter().filter(|x| chains.iter().any(|y| y.has_sv(x))).count()
}
fn count_closeable_chains(chains: &[Chain]) -> usize {
    chains.iter().filter(|x| x.could_close_chain()).count()
}

pub fn get_adjacent_major_ap_ratio(var: &SvVarData) -> f64 {
    let mut r = get_major_allele_jcn_ratio(&var.be(true));
    if !var.is_sgl_breakend() {
        r = jmax(get_major_allele_jcn_ratio(&var.be(false)), r);
    }
    r
}

pub fn variant_exceeds_both_adjacent_jcn(var: &SvVarData) -> bool {
    if var.is_sgl_breakend() {
        return false;
    }
    jmin(get_major_allele_jcn_ratio(&var.be(true)), get_major_allele_jcn_ratio(&var.be(false))) >= ADJACENT_JCN_RATIO
}

pub fn get_major_allele_jcn_ratio(b: &SvBreakend) -> f64 {
    b.jcn() / jmax(b.major_allele_jcn(b.reverse_orient()), 0.01)
}

pub const SEG_DATA_COUNT: usize = 0;
const SEG_DATA_MAX: usize = 1;
const SEG_DATA_SUM: usize = 2;
const OPEN_CHAIN_ID: i32 = -1;
const MIN_CLOSED_SEG_LENGTH: i64 = 1500;
const MIN_CLOSED_SEG_RATIO: f64 = 0.66;
const LOW_JCN_BUFFER: f64 = 4.0;

pub struct DoubleMinuteData {
    cluster: Cluster,
    pub svs: Vec<Sv>,
    pub unchained_svs: Vec<Sv>,
    pub max_jcn: f64,
    pub min_adj_ma_jcn_ratio: f64,
    pub chains: Vec<Chain>,
    pub fully_chained: bool,
    pub valid_chains: Vec<Chain>,
    pub valid_svs: Vec<Sv>,
    pub closed_segment_length: i64,
    pub closed_breakends: i32,
    pub closed_jcn_total: f64,
    pub open_breakends: i32,
    pub open_jcn_total: f64,
    pub open_jcn_max: f64,
    pub simple_dels: i32,
    pub non_segment_foldbacks: i32,
    pub non_segment_foldback_jcn_total: f64,
    pub chains_centromere: bool,
    pub low_depth_count: bool,
    fb_internal_data: [f64; 3],
    pub chain_int_ext_data: JavaHashMap<i32, [f64; 3]>,
    chain_sgl_internal_data: JavaHashMap<i32, [f64; 3]>,
    chain_inf_internal_data: JavaHashMap<i32, [f64; 3]>,
    dm_regions: Vec<ChrBaseRegion>,
    pub is_double_minute: bool,
}

#[derive(Clone, Copy)]
enum SegMap {
    IntExt,
    Sgl,
    Inf,
}

impl DoubleMinuteData {
    pub fn new(cluster: &Cluster, svs: Vec<Sv>) -> Self {
        DoubleMinuteData {
            cluster: cluster.clone(),
            svs,
            unchained_svs: Vec::new(),
            max_jcn: 0.0,
            min_adj_ma_jcn_ratio: 0.0,
            chains: Vec::new(),
            fully_chained: false,
            valid_chains: Vec::new(),
            valid_svs: Vec::new(),
            closed_segment_length: 0,
            closed_breakends: 0,
            closed_jcn_total: 0.0,
            open_breakends: 0,
            open_jcn_total: 0.0,
            open_jcn_max: 0.0,
            simple_dels: 0,
            non_segment_foldbacks: 0,
            non_segment_foldback_jcn_total: 0.0,
            chains_centromere: false,
            low_depth_count: false,
            fb_internal_data: [0.0; 3],
            chain_int_ext_data: JavaHashMap::new(),
            chain_sgl_internal_data: JavaHashMap::new(),
            chain_inf_internal_data: JavaHashMap::new(),
            dm_regions: Vec::new(),
            is_double_minute: false,
        }
    }

    pub fn annotate(&mut self, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
        for var in self.svs.clone() {
            if !self.chains.iter().any(|x| x.has_sv(&var)) {
                self.unchained_svs.push(var.clone());
            }
            let r = get_adjacent_major_ap_ratio(&var);
            self.min_adj_ma_jcn_ratio = if self.min_adj_ma_jcn_ratio == 0.0 { r } else { jmin(r, self.min_adj_ma_jcn_ratio) };
            self.max_jcn = jmax(self.max_jcn, var.jcn());
        }
        self.build_dm_regions();
        let all_pairs: Vec<Pair> = self.chains.iter().flat_map(|x| x.linked_pairs()).collect();
        self.set_chain_characteristics(chr_breakend_map, &all_pairs);
        for var in self.cluster.foldbacks() {
            if var.is_chained_foldback() {
                continue;
            }
            if all_pairs.iter().any(|x| {
                x.chromosome() == var.chromosome(true)
                    && hmf_common::regions::positions_overlap(
                        var.position(true),
                        var.position(false),
                        x.breakend(true).position(),
                        x.breakend(false).position(),
                    )
            }) {
                continue;
            }
            self.non_segment_foldbacks += 1;
            self.non_segment_foldback_jcn_total += var.jcn();
        }
        self.is_double_minute = self.check_criteria(&all_pairs);
    }

    fn set_chain_characteristics(&mut self, chr_breakend_map: &JavaHashMap<String, Vec<Be>>, all_pairs: &[Pair]) {
        let mut observed_breakends: HashSet<Be> = HashSet::new();
        let mut observed_svs: HashSet<Sv> = HashSet::new();
        for chain in self.chains.clone() {
            let simple_dels: Vec<Sv> =
                chain.sv_list().into_iter().filter(|x| x.sv_type() == DEL).filter(|x| x.nearest_sv_relation() == RELATION_TYPE_NEIGHBOUR).collect();
            self.simple_dels += simple_dels.len() as i32;
            for pair in chain.linked_pairs() {
                self.closed_segment_length += pair.base_length() as i64;
                if !contains(&simple_dels, &pair.first()) {
                    self.closed_breakends += 1;
                    self.closed_jcn_total += pair.first().jcn();
                }
                if !contains(&simple_dels, &pair.second()) {
                    self.closed_breakends += 1;
                    self.closed_jcn_total += pair.second().jcn();
                }
                if pair.first_breakend().arm() != pair.second_breakend().arm() {
                    self.chains_centromere = true;
                }
                self.set_pair_characteristics(&chain, &pair, chr_breakend_map, all_pairs, &mut observed_breakends, &mut observed_svs);
            }
            if !chain.is_closed_loop() {
                if !chain.first_sv().is_sgl_breakend() {
                    self.open_jcn_total += chain.first_sv().jcn();
                    self.open_breakends += 1;
                    self.open_jcn_max = jmax(self.open_jcn_max, chain.first_sv().jcn());
                }
                if !chain.last_sv().is_sgl_breakend() {
                    self.open_jcn_total += chain.last_sv().jcn();
                    self.open_breakends += 1;
                    self.open_jcn_max = jmax(self.open_jcn_max, chain.last_sv().jcn());
                }
            }
        }
        for var in &self.unchained_svs {
            self.open_jcn_total += var.jcn() * 2.0;
            self.open_breakends += 2;
            self.open_jcn_max = jmax(self.open_jcn_max, var.jcn());
        }
    }

    fn seg_data(&mut self, chain: &Chain, which: SegMap) -> &mut [f64; 3] {
        let id = if chain.is_closed_loop() { chain.id() } else { OPEN_CHAIN_ID };
        let map = match which {
            SegMap::IntExt => &mut self.chain_int_ext_data,
            SegMap::Sgl => &mut self.chain_sgl_internal_data,
            SegMap::Inf => &mut self.chain_inf_internal_data,
        };
        if !map.contains_key(&id) {
            map.insert(id, [0.0; 3]);
        }
        map.get_mut(&id).unwrap()
    }

    fn set_pair_characteristics(
        &mut self, chain: &Chain, pair: &Pair, chr_breakend_map: &JavaHashMap<String, Vec<Be>>, all_pairs: &[Pair],
        observed_breakends: &mut HashSet<Be>, observed_svs: &mut HashSet<Sv>,
    ) {
        let list = chr_breakend_map.get(&pair.chromosome()).expect("breakends of chromosome").clone();
        let start_index = pair.breakend(true).chr_pos_index();
        let end_index = pair.breakend(false).chr_pos_index();
        let mut assembled: HashSet<Be> = HashSet::new();
        let mut index = start_index + 1;
        while index < end_index {
            let b = list[index as usize].clone();
            if contains(&self.svs, b.sv()) {
                index += 1;
                continue;
            }
            if observed_breakends.contains(&b) {
                index += 1;
                continue;
            }
            observed_breakends.insert(b.clone());
            if observed_svs.contains(b.sv()) {
                index += 1;
                continue;
            }
            observed_svs.insert(b.sv().clone());
            if !b.sv().is_sgl_breakend() {
                let other = b.other_breakend().unwrap();
                let mut other_in_segment = all_pairs.iter().any(|x| {
                    x.chromosome() == other.chromosome() && position_within(other.position(), x.breakend(true).position(), x.breakend(false).position())
                });
                if !other_in_segment {
                    other_in_segment = self.dm_regions.iter().any(|x| x.chromosome == other.chromosome() && x.contains_position(other.position()));
                }
                if !other_in_segment {
                    let in_short_link = assembled.contains(&b)
                        || (if index < end_index - 1 { in_short_assembled_link(&b, &list[(index + 1) as usize]) } else { false });
                    if in_short_link {
                        assembled.insert(b.clone());
                        assembled.insert(list[(index + 1) as usize].clone());
                    } else {
                        let d = self.seg_data(chain, SegMap::IntExt);
                        d[SEG_DATA_COUNT] += 1.0;
                        d[SEG_DATA_SUM] += b.jcn();
                        d[SEG_DATA_MAX] = jmax(d[SEG_DATA_MAX], b.jcn());
                    }
                }
            }
            let mut sv_type = b.sv_type();
            if sv_type != SGL && sv_type != INF {
                if b.is_foldback() && !b.sv().is_chained_foldback() {
                    sv_type = INV;
                } else {
                    index += 1;
                    continue;
                }
            }
            let data: &mut [f64; 3] = if sv_type == INV {
                &mut self.fb_internal_data
            } else {
                if sv_type == SGL && index < end_index - 1 {
                    let nb = list[(index + 1) as usize].clone();
                    if in_sgl_pair_cluster(&b, &nb) {
                        observed_breakends.insert(nb.clone());
                        observed_svs.insert(nb.sv().clone());
                        index += 1;
                        continue;
                    }
                }
                self.seg_data(chain, if sv_type == SGL { SegMap::Sgl } else { SegMap::Inf })
            };
            data[SEG_DATA_COUNT] += 1.0;
            data[SEG_DATA_SUM] += b.jcn();
            data[SEG_DATA_MAX] = jmax(data[SEG_DATA_MAX], b.jcn());
            index += 1;
        }
    }

    fn build_dm_regions(&mut self) {
        if self.svs.len() == 1 {
            let var = &self.svs[0];
            if var.sv_type() == DUP {
                self.dm_regions.push(ChrBaseRegion::new(var.chromosome(true), var.position(true), var.position(false)));
            }
            return;
        }
        for (chromosome, list) in self.cluster.chr_breakend_entries() {
            if !self.svs.iter().any(|x| x.chromosome(true) == chromosome || x.chromosome(false) == chromosome) {
                continue;
            }
            let mut current: Option<ChrBaseRegion> = None;
            let mut region_start_cn = 0.0;
            for b in &list {
                if !contains(&self.svs, b.sv()) {
                    continue;
                }
                if b.reverse_orient() {
                    if current.is_some() {
                        continue;
                    }
                    if get_major_allele_jcn_ratio(b) >= ADJACENT_JCN_RATIO {
                        region_start_cn = b.copy_number_low_side();
                        current = Some(ChrBaseRegion::new(&chromosome, b.position(), 0));
                    }
                } else {
                    let Some(cr) = current.as_mut() else {
                        if let Some(last) = self.dm_regions.last_mut() {
                            last.end = b.position();
                        }
                        continue;
                    };
                    let low = b.copy_number_low_side();
                    if low < region_start_cn || copy_numbers_equal(region_start_cn, low) || get_major_allele_jcn_ratio(b) >= ADJACENT_JCN_RATIO {
                        region_start_cn = 0.0;
                        cr.end = b.position();
                        self.dm_regions.push(current.take().unwrap());
                    }
                }
            }
        }
    }

    fn check_criteria(&mut self, all_pairs: &[Pair]) -> bool {
        let ratio = self.closed_breakends as f64 / (self.closed_breakends + self.open_breakends) as f64;
        if ratio < MIN_CLOSED_SEG_RATIO {
            return false;
        }
        if self.svs.len() == 1 {
            let var = &self.svs[0];
            if var.sv_type() != DUP {
                return false;
            }
            if !variant_exceeds_both_adjacent_jcn(var) {
                return false;
            }
        }
        if self.max_jcn < JCN_UPPER_THRESHOLD {
            if self.svs.iter().any(|x| x.is_sgl_breakend()) {
                return false;
            }
            if self.chains.is_empty() || self.chains.iter().any(|x| !x.is_closed_loop()) {
                return false;
            }
            if self.svs.len() > 1 && !self.svs.iter().any(|x| variant_exceeds_both_adjacent_jcn(x)) {
                return false;
            }
        }
        if all_pairs.iter().all(have_low_adjacent_depth_window_count) {
            self.low_depth_count = true;
            return false;
        }
        self.set_valid_chains()
    }

    fn set_valid_chains(&mut self) -> bool {
        let mut has_valid = false;
        let foldback_jcn = hmf_common::java_math::double_stream_sum(
            self.cluster.foldbacks().iter().filter(|x| !contains(&self.svs, x)).map(|x| if x.is_chained_foldback() { x.jcn() * 0.5 } else { x.jcn() }),
        );
        let mut chain_ids: Vec<i32> = vec![OPEN_CHAIN_ID];
        chain_ids.extend(self.chains.iter().filter(|x| x.is_closed_loop()).map(|x| x.id()));
        let mut non_closed_svs: Vec<Sv> = self.unchained_svs.clone();
        let open_chains: Vec<Chain> = self.chains.iter().filter(|x| !x.is_closed_loop()).cloned().collect();
        for c in &open_chains {
            non_closed_svs.extend(c.sv_list());
        }
        let non_closed_length: i64 = self.chains.iter().filter(|x| !x.is_closed_loop()).map(|x| x.length(false)).sum();
        for chain_id in chain_ids {
            let chain = self.chains.iter().find(|x| x.id() == chain_id).cloned();
            if non_closed_svs.is_empty() && chain_id == OPEN_CHAIN_ID {
                continue;
            }
            if let Some(c) = &chain {
                if c.length(false) < MIN_CLOSED_SEG_LENGTH {
                    continue;
                }
            } else if non_closed_length < MIN_CLOSED_SEG_LENGTH {
                continue;
            }
            let mut max_chain_jcn = 0.0;
            let mut max_chain_foldback_jcn = 0.0;
            let svs: Vec<Sv> = match &chain {
                Some(c) => c.sv_list(),
                None => non_closed_svs.clone(),
            };
            for var in &svs {
                max_chain_jcn = jmax(max_chain_jcn, var.jcn());
                if var.is_foldback() {
                    max_chain_foldback_jcn = if max_chain_foldback_jcn == 0.0 { var.jcn() } else { jmin(max_chain_foldback_jcn, var.jcn()) };
                }
            }
            let max_jcn = if max_chain_foldback_jcn > 0.0 { max_chain_foldback_jcn } else { max_chain_jcn };
            let sgl = self.chain_sgl_internal_data.get(&chain_id).map(|x| x[SEG_DATA_MAX]).unwrap_or(0.0);
            let inf = self.chain_inf_internal_data.get(&chain_id).map(|x| x[SEG_DATA_MAX]).unwrap_or(0.0);
            let sgl_inf_max = jmax(sgl, inf);
            let int_ext_total = self.chain_int_ext_data.get(&chain_id).map(|x| x[SEG_DATA_SUM]).unwrap_or(0.0);
            let opposing = int_ext_total + sgl_inf_max + foldback_jcn;
            if opposing + LOW_JCN_BUFFER <= max_jcn {
                has_valid = true;
                if let Some(c) = &chain {
                    self.valid_chains.push(c.clone());
                    self.valid_svs.extend(c.sv_list());
                } else {
                    self.valid_svs.extend(non_closed_svs.iter().cloned());
                    self.valid_chains.extend(open_chains.iter().cloned());
                }
            }
        }
        has_valid
    }
}

fn in_sgl_pair_cluster(b: &Be, nb: &Be) -> bool {
    if nb.sv_type() != SGL || nb.cluster() != b.cluster() {
        return false;
    }
    b.cl().resolved_type().is_simple() || classify_single_pair_resolved_type(b.sv(), nb.sv()).is_simple()
}

fn in_short_assembled_link(b: &Be, nb: &Be) -> bool {
    if b.linked_pairs().iter().filter(|x| x.is_assembled()).any(|x| x.has_breakend(nb)) {
        return true;
    }
    let (Some(o), Some(no)) = (b.other_breakend(), nb.other_breakend()) else { return false };
    o.linked_pairs().iter().filter(|x| x.is_assembled()).any(|x| x.has_breakend(&no))
}

fn have_low_adjacent_depth_window_count(pair: &Pair) -> bool {
    if !pair.first().is_sgl_breakend() || !pair.second().is_sgl_breakend() {
        return false;
    }
    if (pair.first_breakend().chr_pos_index() - pair.second_breakend().chr_pos_index()).abs() != 1 {
        return false;
    }
    let first = pair.first_breakend();
    let second = pair.second_breakend();
    let d1 = first.sv().copy_number_data(true, first.forward_orient()).expect("CN data").depth_window_count;
    let d2 = second.sv().copy_number_data(true, second.forward_orient()).expect("CN data").depth_window_count;
    d1.min(d2) < MIN_SEGMENT_DEPTH_WINDOW_COUNT
}
