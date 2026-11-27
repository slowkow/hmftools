//! Package `analysis`: clustering, chaining calls, classification and cluster annotation (ClusterAnalyser).
//!
//! Left out, since it only feeds a diagnostic file: reportUnderclustering (UNDER_CLUSTERING).

pub mod annotations;
pub mod classification;
pub mod clustering;
pub mod double_minute;

use hmf_common::java_hash::JavaHashMap;

use crate::analysis::classification::{is_filtered_resolved_type, is_simple_single_sv, set_cluster_resolved_state};
use crate::analysis::clustering::*;
use crate::analysis::double_minute::DoubleMinuteFinder;
use crate::annotators::{LineElementAnnotator, PseudoGeneFinder};
use crate::chaining::finder::ChainFinder;
use crate::chaining::utils::{create_assembly_linked_pairs, find_deletion_bridges, get_min_templated_insertion_length};
use crate::cn::{CnDataLoader, TelomereCentromereCnData};
use crate::types::*;
use crate::utils::*;

pub use annotations::ClusterMetrics;

const SMALL_CLUSTER_SIZE: usize = 3;

/// ClusterAnalyser
pub struct ClusterAnalyser {
    pub is_germline: bool,
    pub chaining_sv_limit: i32,
    pub state: ClusteringState,
    filters: SvFiltering,
    simple_clustering: SimpleClustering,
    pub dm_finder: DoubleMinuteFinder,
    sample_id: String,
    pub clusters: Vec<Cluster>,
    pub artifact_clusters: Vec<Cluster>,
    all_variants: Vec<Sv>,
    chain_finder: ChainFinder,
    pub run_validation_checks: bool,
}

impl ClusterAnalyser {
    pub fn new(proximity_distance: i32, is_germline: bool, chaining_sv_limit: i32, log_dm_candidates: bool, show_dm_candidates: bool) -> Self {
        let mut chain_finder = ChainFinder::new();
        chain_finder.set_use_allele_ploidies(true);
        ClusterAnalyser {
            is_germline,
            chaining_sv_limit,
            state: ClusteringState::default(),
            filters: SvFiltering::default(),
            simple_clustering: SimpleClustering::new(proximity_distance),
            dm_finder: DoubleMinuteFinder::new(log_dm_candidates, show_dm_candidates),
            sample_id: String::new(),
            clusters: Vec::new(),
            artifact_clusters: Vec::new(),
            all_variants: Vec::new(),
            chain_finder,
            run_validation_checks: false,
        }
    }

    pub fn set_sample_data(&mut self, sample_id: &str, all: &[Sv], loader: &CnDataLoader) {
        self.sample_id = sample_id.to_string();
        self.all_variants = all.to_vec();
        self.clusters.clear();
        self.artifact_clusters.clear();
        self.simple_clustering.initialise(sample_id);
        self.chain_finder.set_sample_id(sample_id);
        self.state.loh_event_list = loader.loh_event_data.clone();
        self.state.hom_loss_list = loader.hom_loss_data.clone();
    }

    /// getAllClusters: the clusters, then the artifact clusters
    pub fn all_clusters(&self) -> Vec<Cluster> {
        let mut v = self.clusters.clone();
        v.extend(self.artifact_clusters.iter().cloned());
        v
    }

    pub fn pre_clustering_preparation(&mut self) {
        self.state.reset();
        populate_chromosome_breakend_map(&self.all_variants, &mut self.state);
        self.filters.apply_filters(&mut self.state);
        annotate_nearest_sv_data(&self.state.chr_breakend_map);
        find_deletion_bridges(&self.state.chr_breakend_map);
        set_simple_variant_lengths(&mut self.state);
    }

    pub fn cluster_and_analyse(&mut self, line_annotator: &LineElementAnnotator, pseudo: Option<&PseudoGeneFinder>, tc: &JavaHashMap<String, TelomereCentromereCnData>) -> bool {
        if self.is_germline {
            return self.cluster_and_analyse_germline(line_annotator, pseudo);
        }
        self.clusters.clear();
        self.artifact_clusters.clear();
        self.dm_finder.clear();
        self.filters.cluster_excluded_variants(&mut self.artifact_clusters, &self.state);
        for c in &self.artifact_clusters {
            build_arm_clusters(c);
        }
        self.simple_clustering.cluster_by_proximity(&mut self.clusters, &self.state);
        for c in self.clusters.iter().filter(|x| x.sv_count() > 1) {
            c.update_cluster_details();
        }
        for c in &self.clusters {
            line_annotator.mark_line_cluster(c, pseudo);
        }
        associate_breakend_cn_events(&self.sample_id, &self.state);
        self.find_limited_chains();
        let mut clusters = std::mem::take(&mut self.clusters);
        self.simple_clustering.merge_clusters(&mut clusters, &self.state);
        self.clusters = clusters;
        // logDetails, which updates the cluster details
        for c in self.clusters.iter().filter(|x| x.sv_count() > 1) {
            c.update_cluster_details();
        }
        mark_foldbacks(&self.state.chr_breakend_map, false);
        let mut clusters = std::mem::take(&mut self.clusters);
        apply_complex_rules(&self.sample_id, &mut clusters, &self.state, false);
        self.simple_clustering.merge_long_del_dup_clusters(&mut clusters);
        self.clusters = clusters;
        self.find_links_and_chains();
        self.dissolve_simple_groups();
        let dm_cluster_ids: Vec<i32> = self.clusters.iter().filter(|x| x.has_annotation(CLUSTER_ANNOT_DM)).map(|x| x.id()).collect();
        if mark_foldbacks(&self.state.chr_breakend_map, true) {
            let mut clusters = std::mem::take(&mut self.clusters);
            apply_complex_rules(&self.sample_id, &mut clusters, &self.state, true);
            self.clusters = clusters;
        }
        for cluster in self.clusters.clone() {
            if !cluster.double_minute_svs().is_empty() || dm_cluster_ids.contains(&cluster.id()) {
                self.dm_finder.analyse_cluster(&cluster, true, &self.state.chr_breakend_map);
            }
            if !cluster.has_annotation(CLUSTER_ANNOT_DM) && !cluster.is_resolved() {
                bfb_analyse_cluster(&cluster, tc);
            }
            if !cluster.is_resolved() {
                if cluster.resolved_type() != ResolvedType::NONE {
                    if cluster.sv_count() == 1 && cluster.resolved_type().is_simple() {
                        cluster.set_resolved(true, cluster.resolved_type());
                    }
                } else {
                    self.set_cluster_resolved_state(&cluster, true);
                    cluster.update_cluster_details();
                }
            }
            cluster.cache_linked_pairs();
            build_arm_clusters(&cluster);
        }
        true
    }

    fn cluster_and_analyse_germline(&mut self, line_annotator: &LineElementAnnotator, pseudo: Option<&PseudoGeneFinder>) -> bool {
        self.clusters.clear();
        self.artifact_clusters.clear();
        self.filters.cluster_excluded_variants(&mut self.artifact_clusters, &self.state);
        for c in &self.artifact_clusters {
            build_arm_clusters(c);
        }
        self.simple_clustering.cluster_by_proximity(&mut self.clusters, &self.state);
        for c in self.clusters.iter().filter(|x| x.sv_count() > 1) {
            c.update_cluster_details();
        }
        for c in &self.clusters {
            line_annotator.mark_line_cluster(c, pseudo);
        }
        self.find_limited_chains();
        for c in self.clusters.iter().filter(|x| x.sv_count() > 1) {
            c.update_cluster_details();
        }
        mark_foldbacks(&self.state.chr_breakend_map, false);
        self.dissolve_simple_groups();
        for cluster in self.clusters.clone() {
            if !cluster.is_resolved() {
                if cluster.resolved_type() != ResolvedType::NONE {
                    if cluster.sv_count() == 1 && cluster.resolved_type().is_simple() {
                        cluster.set_resolved(true, cluster.resolved_type());
                    }
                } else {
                    self.set_cluster_resolved_state(&cluster, true);
                    cluster.update_cluster_details();
                }
            }
            cluster.cache_linked_pairs();
            build_arm_clusters(&cluster);
        }
        true
    }

    fn find_limited_chains(&mut self) {
        let check_dms = !self.is_germline;
        for cluster in self.clusters.clone() {
            if check_dms && is_simple_single_sv(&cluster) {
                self.dm_finder.analyse_cluster(&cluster, false, &self.state.chr_breakend_map);
                self.set_cluster_resolved_state(&cluster, false);
                continue;
            }
            let is_simple = cluster.sv_count() <= SMALL_CLUSTER_SIZE && cluster.is_consistent() && !cluster.has_varied_jcn();
            cluster.set_assembly_linked_pairs(create_assembly_linked_pairs(&cluster));
            cluster.determine_requires_replication();
            if check_dms {
                self.dm_finder.analyse_cluster(&cluster, false, &self.state.chr_breakend_map);
            }
            self.find_chains(&cluster, !is_simple);
            if is_simple {
                self.set_cluster_resolved_state(&cluster, false);
            }
        }
    }

    fn find_links_and_chains(&mut self) {
        for cluster in self.clusters.clone() {
            if cluster.resolved_type() == ResolvedType::LINE {
                continue;
            }
            if cluster.resolved_type() != ResolvedType::NONE {
                continue;
            }
            if is_simple_single_sv(&cluster) || cluster.is_fully_chained(false) || cluster.sv_count() < 2 {
                self.set_cluster_resolved_state(&cluster, true);
                continue;
            }
            cluster.dissolve_links_and_chains();
            self.dm_finder.analyse_cluster(&cluster, true, &self.state.chr_breakend_map);
            cluster.determine_requires_replication();
            self.find_chains(&cluster, false);
            self.set_cluster_resolved_state(&cluster, true);
            cluster.update_cluster_details();
        }
    }

    fn dissolve_simple_groups(&mut self) {
        let simple_groups: Vec<Cluster> = self.clusters.iter().filter(|x| x.resolved_type() == ResolvedType::SIMPLE_GRP).cloned().collect();
        let mut discard: Vec<Sv> = Vec::new();
        for c in &simple_groups {
            remove_item(&mut self.clusters, c);
            discard.extend(c.svs());
        }
        for var in &discard {
            let c = Cluster::create(self.state.next_cluster_id());
            var.clear_clustering_data();
            c.add_variant(var);
            self.dm_finder.analyse_cluster(&c, false, &self.state.chr_breakend_map);
            self.set_cluster_resolved_state(&c, true);
            self.clusters.push(c);
        }
    }

    fn set_cluster_resolved_state(&self, cluster: &Cluster, is_final: bool) {
        set_cluster_resolved_state(
            cluster,
            is_final,
            self.is_germline,
            self.state.del_cutoff_length,
            self.state.dup_cutoff_length,
            &self.state.chr_breakend_map,
        );
    }

    fn find_chains(&mut self, cluster: &Cluster, assembled_links_only: bool) {
        if self.chaining_sv_limit > 0 && cluster.sv_count() as i32 > self.chaining_sv_limit {
            hmf_common::debug!(
                "sample({}) skipping chaining large cluster({}) with SV count({})",
                self.sample_id,
                cluster.id(),
                cluster.sv_count()
            );
            return;
        }
        cluster.chains_mut().clear();
        self.chain_finder.initialise(cluster);
        self.chain_finder.form_chains(assembled_links_only);
        self.chain_finder.add_chains(cluster);
        if let Some(range) = self.chain_finder.calc_range_data() {
            let mut m = cluster.metrics.borrow_mut();
            m.valid_allele_jcn_segment_perc = self.chain_finder.valid_allele_copy_number_segment_perc();
            m.traversed_range = range[0];
            m.total_deleted = range[1];
        }
        self.chain_finder.clear();
    }

    /// annotateClusters: the templated insertion annotations, then the cluster features
    pub fn annotate_clusters(&self) {
        for c in &self.clusters {
            annotate_templated_insertions(c, &self.state.chr_breakend_map);
        }
        for c in &self.clusters {
            annotations::report_cluster_features(c, &self.state.chr_breakend_map);
        }
    }
}

/// BfbFinder.analyseCluster
fn bfb_analyse_cluster(cluster: &Cluster, tc: &JavaHashMap<String, TelomereCentromereCnData>) {
    if cluster.foldbacks().is_empty() || cluster.sv_count() == 1 {
        return;
    }
    let mut cluster_max_jcn = 0.0;
    for var in cluster.svs() {
        cluster_max_jcn = if var.is_sgl_breakend() { jmax(cluster_max_jcn, var.jcn_min() * 0.5) } else { jmax(cluster_max_jcn, var.jcn_min()) };
    }
    let mut sum_fb = 0.0;
    let mut max_fb = 0.0;
    let mut fb_count = 0.0;
    let mut max_sgl = 0.0;
    for var in cluster.svs() {
        if var.is_foldback() {
            sum_fb += if var.is_chained_foldback() { var.jcn() * 0.5 } else { var.jcn() };
            fb_count += if var.is_chained_foldback() { 0.5 } else { 1.0 };
            max_fb = jmax(max_fb, var.jcn());
        } else if var.is_sgl_breakend() {
            max_sgl = jmax(var.jcn(), max_sgl);
        }
    }
    let mut max_arm_end_cn = 0.0;
    for g in cluster.arm_groups().iter() {
        if let Some(t) = tc.get(&g.chromosome) {
            let c = if g.arm == Arm::P { t.centromere_p_arm } else { t.centromere_q_arm };
            let te = if g.arm == Arm::P { t.telomere_p_arm } else { t.telomere_q_arm };
            max_arm_end_cn = jmax(jmax(te, c), max_arm_end_cn);
        }
    }
    let fb_sgl_count = fb_count + if max_sgl > 0.1 * cluster_max_jcn { 1.0 } else { 0.0 };
    let arm_bfb_cn = max_arm_end_cn * 2f64.powf(fb_sgl_count);
    let max_bfb = jmin(2.0 * sum_fb + max_sgl, jmin(6.0 * max_fb, arm_bfb_cn));
    if max_bfb > cluster_max_jcn && fb_count > 0.0 {
        cluster.add_annotation(CLUSTER_ANNOT_BFB);
    }
}

/// ClusterAnnotations.annotateTemplatedInsertions
pub fn annotate_templated_insertions(cluster: &Cluster, chr_breakend_map: &JavaHashMap<String, Vec<Be>>) {
    if cluster.chain_count() == 0 {
        return;
    }
    for chain in cluster.chains() {
        let cs = chain.open_breakend(true);
        let ce = chain.open_breakend(false);
        let mut chain_lower: Option<Be> = None;
        let mut chain_upper: Option<Be> = None;
        let mut start_end_same_arm = false;
        let mut ends_cn_match = false;
        let mut ends_cn = 0.0;
        if let (Some(s), Some(e)) = (&cs, &ce) {
            let lower = if s.position() < e.position() { s.clone() } else { e.clone() };
            let upper = if &lower == s { e.clone() } else { s.clone() };
            start_end_same_arm = e.chr_arm() == s.chr_arm();
            let (scn, ecn) = (s.copy_number(), e.copy_number());
            ends_cn_match = copy_numbers_equal(scn, ecn);
            if ends_cn_match {
                ends_cn = jmax(scn, ecn);
            }
            chain_lower = Some(lower);
            chain_upper = Some(upper);
        }
        let pairs = chain.linked_pairs();
        for pair in &pairs {
            if pair.first().is_sgl_breakend() || pair.second().is_sgl_breakend() {
                continue;
            }
            let pl = pair.breakend(true);
            let pu = pair.breakend(false);
            let pair_cn = link_pair_copy_number(pair);
            if ends_cn_match && pair_cn > ends_cn && !copy_numbers_equal(pair_cn, ends_cn) {
                pair.set_has_copy_number_gain(true);
            }
            if start_end_same_arm {
                let cl = chain_lower.as_ref().unwrap();
                let cu = chain_upper.as_ref().unwrap();
                if pl.chromosome() == cl.chromosome() {
                    let lb = get_min_templated_insertion_length(&pl, cl);
                    let ub = get_min_templated_insertion_length(&pu, cu);
                    if pl.position() >= cl.position() - lb && pu.position() <= cu.position() + ub {
                        pair.set_location_type(LOCATION_TYPE_INTERNAL);
                    } else {
                        pair.set_location_type(LOCATION_TYPE_EXTERNAL);
                    }
                } else {
                    pair.set_location_type(LOCATION_TYPE_REMOTE);
                }
            } else if cs.as_ref().is_some_and(|s| pl.chr_arm() == s.chr_arm()) || ce.as_ref().is_some_and(|e| pl.chr_arm() == e.chr_arm()) {
                pair.set_location_type(LOCATION_TYPE_EXTERNAL);
            } else {
                pair.set_location_type(LOCATION_TYPE_REMOTE);
            }
            let mut unique: Vec<Pair> = Vec::new();
            let mut overlap_count = 0;
            for other in &pairs {
                if other == pair {
                    continue;
                }
                if unique.iter().any(|x| x.matches(other)) {
                    continue;
                }
                if pair.chromosome() != other.chromosome() || pair.position_distance() <= other.position_distance() {
                    continue;
                }
                unique.push(other.clone());
                if pair.breakend(true).position() > other.breakend(false).position() || pair.breakend(false).position() < other.breakend(true).position() {
                    continue;
                }
                let os = pair.breakend(true).position().max(other.breakend(true).position());
                let oe = pair.breakend(false).position().min(other.breakend(false).position());
                if oe - os >= SHORT_TI_LENGTH {
                    overlap_count += 1;
                }
            }
            pair.set_overlap_count(overlap_count);
            let list = chr_breakend_map.get(&pair.chromosome()).expect("breakends of chromosome");
            let next = next_cluster_sv_data(cluster, list, pair);
            pair.set_next_sv_data(next[0], next[1]);
            pair.set_traversed_sv_count(traversed_sv_count(cluster, list, pair.breakend(true).chr_pos_index(), pair.breakend(false).chr_pos_index()));
        }
    }
}

fn link_pair_copy_number(pair: &Pair) -> f64 {
    let mut cn = [0.0, 0.0];
    for s in [SE_START, SE_END] {
        let b = pair.breakend_se(s);
        cn[s] = match b.db_link() {
            Some(db) if db.length() < 1 => jmax(b.copy_number() - db.other_breakend(&b).jcn(), 0.0),
            _ => b.copy_number(),
        };
    }
    (cn[0] + cn[1]) * 0.5
}

fn traversed_sv_count(cluster: &Cluster, list: &[Be], lower: i32, upper: i32) -> i32 {
    if lower >= upper - 1 {
        return 0;
    }
    if lower < 0 || upper as usize >= list.len() {
        hmf_common::error!("invalid indices({} & {}) vs breakend list size({})", lower, upper, list.len());
        return 0;
    }
    let mut count = 0;
    for i in (lower + 1)..=(upper - 1) {
        let oc = list[i as usize].cl();
        if &oc == cluster || oc.is_resolved() {
            continue;
        }
        count += 1;
    }
    count
}

fn next_cluster_sv_data(cluster: &Cluster, list: &[Be], pair: &Pair) -> [i32; 2] {
    let mut d = [-1, -1];
    let lower = pair.breakend(true);
    let upper = pair.breakend(false);
    let li = lower.chr_pos_index();
    let ui = upper.chr_pos_index();
    if li > 0 {
        for i in (0..li).rev() {
            let b = &list[i as usize];
            let dist = lower.position() - b.position();
            if &b.cl() == cluster {
                d[1] = dist;
                break;
            } else if !is_filtered_resolved_type(b.cl().resolved_type()) && d[0] == -1 {
                d[0] = dist;
            }
        }
    }
    if ui < list.len() as i32 - 1 {
        for b in list.iter().skip((ui + 1) as usize) {
            let dist = b.position() - upper.position();
            if &b.cl() == cluster {
                if d[1] == -1 || dist < d[1] {
                    d[1] = dist;
                }
                break;
            } else if !is_filtered_resolved_type(b.cl().resolved_type()) && (d[0] == -1 || dist < d[0]) {
                d[0] = dist;
            }
        }
    }
    if d[0] == -1 {
        d[0] = d[1];
    }
    d
}

/// ClusterMetrics.findStartIndex: the first breakend of a cluster on a chromosome, past any short TIs
pub fn find_start_index(list: &[Be]) -> i32 {
    let mut start = 0usize;
    while start + 1 < list.len() {
        let lower = &list[start];
        let upper = &list[start + 1];
        if upper.position() - lower.position() <= SHORT_TI_LENGTH
            && lower.reverse_orient()
            && upper.forward_orient()
            && lower.linked_pairs().iter().any(|x| upper.linked_pairs().contains(x))
        {
            start += 2;
            continue;
        }
        break;
    }
    start as i32
}

/// ClusterMetrics.findEndIndex
pub fn find_end_index(list: &[Be]) -> i32 {
    let mut end = list.len() as i32 - 1;
    while end >= 1 {
        let lower = &list[end as usize - 1];
        let upper = &list[end as usize];
        if upper.position() - lower.position() <= SHORT_TI_LENGTH
            && lower.reverse_orient()
            && upper.forward_orient()
            && lower.linked_pairs().iter().any(|x| upper.linked_pairs().contains(x))
        {
            end -= 2;
            continue;
        }
        break;
    }
    end
}
