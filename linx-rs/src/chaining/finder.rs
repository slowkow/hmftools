//! ChainFinder, with ChainLinkAllocator, ChainRuleSelector, ProposedLinks, ChainState, BreakendJcn and
//! FoldbackBreakendPair: forms chains from a cluster's SVs. The Java classes share their lists by
//! reference; here they are one struct, the allocator's and rule selector's methods in their own impl
//! blocks. ChainDiagnostics only logs and writes cohort files, and is left out.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use hmf_common::java_hash::JavaHashMap;
use hmf_common::sv::StructuralVariantType::*;

use crate::chaining::jcn_limits::*;
use crate::chaining::line::LineChainer;
use crate::chaining::utils::*;
use crate::cn::JcnCalcData;
use crate::idmap::IdMap;
use crate::types::*;
use crate::utils::{copy_numbers_equal, doubles_equal, format_jcn, jmax, jmin};

pub const MIN_CHAINING_JCN_LEVEL: f64 = 0.05;
pub const LR_METHOD_DM_CLOSE: &str = "DM_CLOSE";

/// ChainingRule
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainingRule {
    ASSEMBLY,
    CA_JCN_SUPPORT,
    ONLY,
    FOLDBACK_SPLIT,
    COMP_DUP_SPLIT,
    FOLDBACK,
    JCN_MATCH,
    JCN_OVERLAP,
    ADJACENT,
    ADJACENT_MATCH,
    JCN_MAX,
    NEAREST,
}

impl ChainingRule {
    pub fn priority(self) -> i32 {
        use ChainingRule::*;
        match self {
            CA_JCN_SUPPORT => 9,
            ONLY => 8,
            FOLDBACK_SPLIT | COMP_DUP_SPLIT => 7,
            FOLDBACK => 6,
            JCN_MATCH => 5,
            JCN_OVERLAP => 4,
            ADJACENT | ADJACENT_MATCH => 3,
            JCN_MAX => 2,
            NEAREST => 1,
            ASSEMBLY => 0,
        }
    }
}

impl std::fmt::Display for ChainingRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

fn calc_rule_priority(rules: &[ChainingRule]) -> i32 {
    rules.iter().map(|r| 2f64.powi(r.priority()) as i32).sum()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinkSkipType {
    CLOSING,
    JCN_MISMATCH,
}

pub const PM_MATCHED: &str = "Matched";
pub const PM_NONE: &str = "None";
pub const PM_OVERLAP: &str = "Overlap";

/// ProposedLinks
pub struct ProposedLinks {
    pub links: Vec<Pair>,
    rules: Vec<ChainingRule>,
    priority: i32,
    chain_target: Option<Chain>,
    foldback_chain: Option<Chain>,
    jcn: f64,
    breakend_jcn: HashMap<Be, f64>,
    exhaust_breakend: HashMap<Be, bool>,
    jcn_match_type: &'static str,
    shortest_distance: i32,
}

impl ProposedLinks {
    pub fn single(link: &Pair, rule: ChainingRule) -> Self {
        let rules = vec![rule];
        ProposedLinks {
            links: vec![link.clone()],
            priority: calc_rule_priority(&rules),
            rules,
            chain_target: None,
            foldback_chain: None,
            jcn: 0.0,
            breakend_jcn: HashMap::new(),
            exhaust_breakend: HashMap::new(),
            jcn_match_type: PM_NONE,
            shortest_distance: link.base_length(),
        }
    }
    pub fn multi(links: Vec<Pair>, rule: ChainingRule, target: Option<Chain>, foldback_chain: Option<Chain>) -> Self {
        let rules = vec![rule];
        let shortest = links.iter().map(|x| x.position_distance()).min().expect("proposed links");
        ProposedLinks {
            links,
            priority: calc_rule_priority(&rules),
            rules,
            chain_target: target,
            foldback_chain,
            jcn: 0.0,
            breakend_jcn: HashMap::new(),
            exhaust_breakend: HashMap::new(),
            jcn_match_type: PM_NONE,
            shortest_distance: shortest,
        }
    }
    pub fn target_chain(&self) -> Option<&Chain> {
        self.chain_target.as_ref()
    }
    pub fn foldback_chain(&self) -> Option<&Chain> {
        self.foldback_chain.as_ref()
    }
    pub fn splitting_rule(&self) -> ChainingRule {
        if self.rules.contains(&ChainingRule::FOLDBACK_SPLIT) { ChainingRule::FOLDBACK_SPLIT } else { ChainingRule::COMP_DUP_SPLIT }
    }
    pub fn multi_connection(&self) -> bool {
        self.links.len() > 1
    }
    pub fn jcn(&self) -> f64 {
        self.jcn
    }
    pub fn set_lower_jcn(&mut self, jcn: f64) {
        if jcn >= self.jcn || self.multi_connection() {
            return;
        }
        self.jcn = jcn;
        let entries: Vec<(Be, f64)> = self.breakend_jcn.iter().map(|(k, v)| (k.clone(), *v)).collect();
        for (b, v) in entries {
            self.exhaust_breakend.insert(b, doubles_equal(v, self.jcn));
        }
    }
    pub fn add_comp_dup_breakends(&mut self, cd_start: &Be, cd_end: &Be, cd_jcn: f64, other_start: &Be, other_end: &Be, other_jcn: f64) {
        if cd_jcn == 0.0 || other_jcn == 0.0 {
            return;
        }
        self.jcn = cd_jcn;
        self.breakend_jcn.insert(cd_start.clone(), cd_jcn);
        self.breakend_jcn.insert(cd_end.clone(), cd_jcn);
        self.exhaust_breakend.insert(cd_start.clone(), true);
        self.exhaust_breakend.insert(cd_end.clone(), true);
        self.breakend_jcn.insert(other_start.clone(), other_jcn);
        self.breakend_jcn.insert(other_end.clone(), other_jcn);
    }
    pub fn add_foldback_breakends(&mut self, fb_start: &Be, fb_end: &Be, fb_jcn: f64, other: &Be, other_jcn: f64, other_unc: f64) {
        if fb_jcn == 0.0 || other_jcn == 0.0 {
            return;
        }
        if copy_numbers_equal(fb_jcn * 2.0, other_jcn) {
            self.jcn_match_type = PM_MATCHED;
        } else if jcn_overlap(fb_jcn * 2.0, fb_start.jcn_uncertainty(), other_jcn, other_unc) {
            self.jcn_match_type = PM_OVERLAP;
        }
        self.jcn = if self.jcn_match_type == PM_NONE { fb_jcn } else { (fb_jcn + other_jcn / 2.0) * 0.5 };
        self.breakend_jcn.insert(fb_start.clone(), fb_jcn);
        self.breakend_jcn.insert(fb_end.clone(), fb_jcn);
        if fb_start.sv() == fb_end.sv() && copy_numbers_equal(fb_start.jcn(), fb_jcn) {
            self.exhaust_breakend.insert(fb_start.clone(), true);
            self.exhaust_breakend.insert(fb_end.clone(), true);
        }
        self.breakend_jcn.insert(other.clone(), other_jcn);
    }
    pub fn add_breakend_ploidies(&mut self, b1: &Be, jcn1: f64, b2: &Be, jcn2: f64) {
        if jcn1 == 0.0 || jcn2 == 0.0 {
            return;
        }
        self.breakend_jcn.insert(b1.clone(), jcn1);
        self.breakend_jcn.insert(b2.clone(), jcn2);
        if copy_numbers_equal(jcn1, jcn2) {
            self.jcn_match_type = PM_MATCHED;
        } else if jcn_overlap(jcn1, b1.jcn_uncertainty(), jcn2, b2.jcn_uncertainty()) {
            self.jcn_match_type = PM_OVERLAP;
        }
        self.jcn = if self.jcn_match_type == PM_NONE { if jcn1 > jcn2 { jcn2 } else { jcn1 } } else { (jcn1 + jcn2) * 0.5 };
    }
    pub fn breakend_jcn(&self, b: &Be) -> f64 {
        // Java unboxes a missing entry: NullPointerException
        *self.breakend_jcn.get(b).expect("breakend JCN")
    }
    pub fn exhaust_breakend(&self, b: &Be) -> bool {
        self.exhaust_breakend.get(b).copied().unwrap_or(false)
    }
    pub fn override_breakend_jcn_matched(&mut self, b: &Be, matched: bool) {
        self.exhaust_breakend.insert(b.clone(), matched);
    }
    pub fn jcn_match_type(&self) -> &'static str {
        self.jcn_match_type
    }
    pub fn link_jcn_match(&self) -> bool {
        self.jcn_match_type != PM_NONE
    }
    pub fn add_rule(&mut self, rule: ChainingRule) {
        if self.rules.contains(&rule) {
            return;
        }
        self.rules.push(rule);
        self.priority = calc_rule_priority(&self.rules);
    }
    pub fn priority(&self) -> i32 {
        self.priority
    }
    pub fn top_rule(&self) -> ChainingRule {
        self.rules[0]
    }
    pub fn has_rule(&self, rule: ChainingRule) -> bool {
        self.rules.contains(&rule)
    }
    pub fn shortest_link_distance(&self) -> i32 {
        self.shortest_distance
    }
    pub fn is_valid(&self) -> bool {
        if self.links.is_empty() || self.links.len() > 2 {
            return false;
        }
        if self.breakend_jcn.is_empty() {
            return false;
        }
        if self.links.len() == 1 && self.breakend_jcn.len() != 2 {
            return false;
        }
        if self.links.len() == 2 && self.breakend_jcn.len() < 3 && !self.links[0].matches(&self.links[1]) {
            return false;
        }
        self.jcn > 0.0
    }
}

impl std::fmt::Display for ProposedLinks {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "pair({}) rule({}) jcn({}) match({}) length({}) priority({})",
            self.links[0],
            self.rules[0],
            format_jcn(self.jcn),
            self.jcn_match_type,
            self.shortest_distance,
            self.priority
        )
    }
}

/// ChainState: the JCN allocated to each breakend of an SV during chaining
pub struct ChainState {
    pub sv: Sv,
    pub min_jcn: f64,
    pub jcn: f64,
    pub max_jcn: f64,
    exhausted_level: f64,
    breakend_count: [f64; 2],
    connections: [Vec<Be>; 2],
}

pub type ChainStateRef = Rc<RefCell<ChainState>>;

impl ChainState {
    pub fn new(var: &Sv, single_jcn: Option<f64>) -> Self {
        let (jcn, min_jcn, max_jcn) = match single_jcn {
            Some(j) => (j, j, j),
            None => (var.jcn(), var.jcn_min(), jmax(var.jcn_max(), var.jcn())),
        };
        ChainState {
            sv: var.clone(),
            min_jcn,
            jcn,
            max_jcn,
            exhausted_level: jmin(EXHAUSTED_JCN_PERC * jcn, EXHAUSTED_JCN_ABS),
            breakend_count: [0.0, 0.0],
            connections: [Vec::new(), Vec::new()],
        }
    }
    pub fn add(&mut self, is_start: bool, link: f64) {
        self.breakend_count[se(is_start)] += link;
    }
    pub fn set(&mut self, is_start: bool, link: f64) {
        self.breakend_count[se(is_start)] = link;
    }
    pub fn current_count(&self) -> f64 {
        jmin(self.breakend_count[0], self.breakend_count[1])
    }
    pub fn breakend_count(&self, is_start: bool) -> f64 {
        self.breakend_count[se(is_start)]
    }
    pub fn breakend_exhausted(&self, is_start: bool) -> bool {
        if self.breakend_count(is_start) == 0.0 {
            return false;
        }
        self.unlinked(is_start) <= self.exhausted_level
    }
    pub fn breakend_exhausted_vs_max(&self, is_start: bool) -> bool {
        self.max_unlinked(is_start) <= self.exhausted_level
    }
    pub fn unlinked_total(&self) -> f64 {
        jmax(self.jcn - self.current_count(), 0.0)
    }
    pub fn unlinked(&self, is_start: bool) -> f64 {
        self.unlinked_se(se(is_start))
    }
    pub fn unlinked_se(&self, s: usize) -> f64 {
        jmax(self.jcn - self.breakend_count[s], 0.0)
    }
    pub fn max_unlinked(&self, is_start: bool) -> f64 {
        jmax(self.max_jcn - self.breakend_count[se(is_start)], 0.0)
    }
    pub fn unique_connections(&self, is_start: bool) -> usize {
        self.connections[se(is_start)].len()
    }
    pub fn has_connections(&self) -> bool {
        !self.connections[0].is_empty() || !self.connections[1].is_empty()
    }
    pub fn add_connection(&mut self, b: &Be, is_start: bool) {
        let l = &mut self.connections[se(is_start)];
        if !contains(l, b) {
            l.push(b.clone());
        }
    }
}

/// BreakendJcn
pub struct BreakendJcn {
    pub unchained_jcn: f64,
    pub total_chain_jcn: f64,
    pub max_jcn_chain: Option<Chain>,
    pub chains: Vec<Chain>,
    pub other_breakend_exhausted: bool,
}

impl BreakendJcn {
    fn simple(unchained: f64, other_exhausted: bool) -> Self {
        BreakendJcn { unchained_jcn: unchained, total_chain_jcn: 0.0, max_jcn_chain: None, chains: Vec::new(), other_breakend_exhausted: other_exhausted }
    }
    pub fn unlinked_jcn(&self) -> f64 {
        if self.total_chain_jcn > 0.0 { jmax(self.max_jcn_chain.as_ref().unwrap().jcn(), self.unchained_jcn) } else { self.unchained_jcn }
    }
    pub fn multi_connections(&self) -> bool {
        self.chains.len() > 1 || (!self.chains.is_empty() && self.unchained_jcn > 0.0)
    }
    pub fn exhausted_in_chain(&self) -> bool {
        self.chains.len() == 1 && doubles_equal(self.unchained_jcn, 0.0)
    }
}

/// FoldbackBreakendPair
#[derive(Clone)]
pub struct FoldbackBreakendPair {
    pub breakend_start: Be,
    pub breakend_end: Be,
    pub ploidy: f64,
    pub chain: Option<Chain>,
}

impl FoldbackBreakendPair {
    fn is_chained(&self) -> bool {
        self.chain.is_some()
    }
    fn matches(&self, o: &FoldbackBreakendPair) -> bool {
        (o.breakend_end == self.breakend_end && o.breakend_start == self.breakend_start)
            || (o.breakend_start == self.breakend_end && o.breakend_end == self.breakend_start)
    }
}

impl std::fmt::Display for FoldbackBreakendPair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} & {} ploidy({})", self.breakend_start, self.breakend_end, format_jcn(self.ploidy))
    }
}

fn add_by_ploidy(pairs: &mut Vec<FoldbackBreakendPair>, new_pair: FoldbackBreakendPair) {
    let mut index = 0;
    while index < pairs.len() {
        if pairs[index].ploidy < new_pair.ploidy {
            break;
        }
        index += 1;
    }
    pairs.insert(index, new_pair);
}

/// ChainFinder with its allocator and rule selector
pub struct ChainFinder {
    cluster_id: i32,
    sample_id: String,
    has_replication: bool,
    is_line_cluster: bool,
    sv_list: Vec<Sv>,
    foldbacks: Vec<Sv>,
    double_minute_svs: Vec<Sv>,
    assembled_links: Vec<Pair>,
    /// the breakend lists, as the cluster's (or the subset's) map iterates them
    chr_breakend_map: Vec<(String, Vec<Be>)>,
    chr_breakend_lookup: HashMap<String, usize>,
    is_cluster_subset: bool,
    subset_breakend_cluster_index_map: HashMap<Be, i32>,
    complex_dup_candidates: IdMap<Sv, Vec<Pair>>,
    adjacent_matching_pairs: Vec<Pair>,
    adjacent_pairs: Vec<Pair>,
    pub chains: Vec<Chain>,
    unique_chains: Vec<Chain>,
    pub jcn_limits: ChainJcnLimits,
    sv_breakend_possible_links: HashMap<Be, Vec<Pair>>,
    finder_is_valid: bool,
    run_validation: bool,
    use_allele_jcns: bool,
    line_chainer: LineChainer,
    // ChainLinkAllocator
    skipped_pairs: HashMap<Pair, LinkSkipType>,
    link_index: i32,
    allocator_is_valid: bool,
    pair_skipped: bool,
    chains_split: bool,
    unique_pairs: Vec<Pair>,
    next_chain_id: i32,
    sv_connections: JavaHashMap<i32, ChainStateRef>,
    sv_completed_connections: Vec<ChainStateRef>,
    // ChainRuleSelector
    rules_to_apply: Vec<ChainingRule>,
    foldback_breakend_pairs: Vec<FoldbackBreakendPair>,
    foldbacks_initialised: bool,
}

impl Default for ChainFinder {
    fn default() -> Self {
        Self::new()
    }
}

impl ChainFinder {
    pub fn new() -> Self {
        ChainFinder {
            cluster_id: -1,
            sample_id: String::new(),
            has_replication: false,
            is_line_cluster: false,
            sv_list: Vec::new(),
            foldbacks: Vec::new(),
            double_minute_svs: Vec::new(),
            assembled_links: Vec::new(),
            chr_breakend_map: Vec::new(),
            chr_breakend_lookup: HashMap::new(),
            is_cluster_subset: false,
            subset_breakend_cluster_index_map: HashMap::new(),
            complex_dup_candidates: IdMap::new(),
            adjacent_matching_pairs: Vec::new(),
            adjacent_pairs: Vec::new(),
            chains: Vec::new(),
            unique_chains: Vec::new(),
            jcn_limits: ChainJcnLimits::new(),
            sv_breakend_possible_links: HashMap::new(),
            finder_is_valid: true,
            run_validation: false,
            use_allele_jcns: false,
            line_chainer: LineChainer::new(),
            skipped_pairs: HashMap::new(),
            link_index: 0,
            allocator_is_valid: true,
            pair_skipped: false,
            chains_split: false,
            unique_pairs: Vec::new(),
            next_chain_id: 0,
            sv_connections: JavaHashMap::new(),
            sv_completed_connections: Vec::new(),
            rules_to_apply: Vec::new(),
            foldback_breakend_pairs: Vec::new(),
            foldbacks_initialised: false,
        }
    }

    pub fn set_sample_id(&mut self, sample_id: &str) {
        self.sample_id = sample_id.to_string();
    }
    pub fn set_use_allele_ploidies(&mut self, t: bool) {
        self.use_allele_jcns = t;
    }
    pub fn set_run_validation(&mut self, t: bool) {
        self.run_validation = t;
    }

    pub fn clear(&mut self) {
        self.cluster_id = -1;
        self.is_line_cluster = false;
        self.sv_list.clear();
        self.foldbacks.clear();
        self.double_minute_svs.clear();
        self.is_cluster_subset = false;
        self.assembled_links.clear();
        self.chr_breakend_map.clear();
        self.chr_breakend_lookup.clear();
        self.adjacent_matching_pairs.clear();
        self.adjacent_pairs.clear();
        self.subset_breakend_cluster_index_map.clear();
        self.complex_dup_candidates.clear();
        self.chains.clear();
        self.unique_chains.clear();
        self.sv_breakend_possible_links.clear();
        self.finder_is_valid = true;
        self.line_chainer.clear();
    }

    fn set_breakend_map(&mut self, entries: Vec<(String, Vec<Be>)>) {
        self.chr_breakend_lookup = entries.iter().enumerate().map(|(i, (k, _))| (k.clone(), i)).collect();
        self.chr_breakend_map = entries;
    }
    fn breakend_list(&self, chr: &str) -> Option<&Vec<Be>> {
        self.chr_breakend_lookup.get(chr).map(|i| &self.chr_breakend_map[*i].1)
    }

    /// initialise(cluster): chain all the SVs of a cluster
    pub fn initialise(&mut self, cluster: &Cluster) {
        self.clear();
        self.cluster_id = cluster.id();
        self.is_line_cluster = cluster.resolved_type() == ResolvedType::LINE;
        if self.is_line_cluster {
            self.line_chainer.initialise(&self.sample_id, cluster);
        }
        self.sv_list = cluster.svs();
        self.foldbacks = cluster.foldbacks();
        self.double_minute_svs = cluster.double_minute_svs();
        self.assembled_links = cluster.assembly_linked_pairs();
        self.set_breakend_map(cluster.chr_breakend_entries());
        self.has_replication = cluster.requires_replication();
        self.is_cluster_subset = false;
        self.allocator_initialise(self.cluster_id);
        self.rule_selector_initialise(self.has_replication);
    }

    /// initialise(cluster, svList, requiresReplication): chain a subset of a cluster's SVs
    pub fn initialise_subset(&mut self, cluster: &Cluster, sv_list: &[Sv], requires_replication: bool) {
        self.clear();
        self.is_cluster_subset = true;
        self.cluster_id = cluster.id();
        self.is_line_cluster = false;
        let mut map: JavaHashMap<String, Vec<Be>> = JavaHashMap::new();
        for (chr, list) in cluster.chr_breakend_entries() {
            let mut bl: Vec<Be> = Vec::new();
            for b in &list {
                if contains(sv_list, b.sv()) {
                    self.subset_breakend_cluster_index_map.insert(b.clone(), bl.len() as i32);
                    bl.push(b.clone());
                }
            }
            if !bl.is_empty() {
                map.insert(chr.clone(), bl);
            }
        }
        self.set_breakend_map(map.iter().map(|(k, v)| (k.clone(), v.clone())).collect());
        self.has_replication = requires_replication;
        for var in sv_list {
            self.sv_list.push(var.clone());
            // Java: `if(var.isFoldback() && mFoldbacks.contains(var))`, never true for an empty list
            if var.is_foldback() && contains(&self.foldbacks, var) {
                self.foldbacks.push(var.clone());
            }
            for s in [SE_START, SE_END] {
                for link in var.assembled_linked_pairs(is_start(s)) {
                    let other = link.other_sv(var);
                    if !contains(sv_list, &other) {
                        continue;
                    }
                    self.assembled_links.push(link);
                }
            }
        }
        self.allocator_initialise(self.cluster_id);
        self.rule_selector_initialise(self.has_replication);
        self.double_minute_svs.extend(sv_list.iter().cloned());
    }

    pub fn unique_chains(&self) -> &[Chain] {
        &self.unique_chains
    }
    pub fn valid_allele_copy_number_segment_perc(&self) -> f64 {
        self.jcn_limits.valid_allele_jcn_segment_perc()
    }
    pub fn calc_range_data(&self) -> Option<[i64; 2]> {
        self.jcn_limits.calc_range_data()
    }

    pub fn form_chains(&mut self, assembled_links_only: bool) {
        if self.is_line_cluster {
            self.line_chainer.form_chains();
            self.unique_chains.extend(self.line_chainer.chains().iter().cloned());
            return;
        }
        if self.sv_list.len() < 2 {
            return;
        }
        if self.sv_list.len() >= 4 {
            hmf_common::debug!(
                "cluster({}) starting chaining with assemblyLinks({}) svCount({}) replication({})",
                self.cluster_id,
                self.assembled_links.len(),
                self.sv_list.len(),
                self.has_replication
            );
        }
        let entries = self.chr_breakend_map.clone();
        self.jcn_limits.initialise(self.cluster_id, entries);
        self.build_chains(assembled_links_only);
        self.check_chains();
        self.remove_identical_chains();
        if !self.is_valid() {
            hmf_common::warn!("cluster({}) chain finding failed", self.cluster_id);
        }
    }

    fn is_valid(&self) -> bool {
        self.finder_is_valid && self.allocator_is_valid
    }

    fn remove_identical_chains(&mut self) {
        if !self.has_replication {
            self.unique_chains.extend(self.chains.iter().cloned());
            return;
        }
        for new_chain in self.chains.clone() {
            let mut matched = false;
            for chain in &self.unique_chains {
                if identical_chain(chain, &new_chain, false, false) {
                    hmf_common::debug!(
                        "cluster({}) skipping duplicate chain({}) ploidy({}) vs origChain({}) ploidy({})",
                        self.cluster_id,
                        new_chain.id(),
                        format_jcn(new_chain.jcn()),
                        chain.id(),
                        format_jcn(chain.jcn())
                    );
                    chain.set_jcn_data(chain.jcn() + new_chain.jcn(), chain.jcn_uncertainty());
                    matched = true;
                    break;
                }
            }
            if !matched {
                self.unique_chains.push(new_chain);
            }
        }
    }

    pub fn add_chains(&self, cluster: &Cluster) {
        for chain in &self.unique_chains {
            // checkAddNewChain: with replication, only those in the unique list, which all are
            cluster.add_chain(chain, false);
        }
        for (i, chain) in cluster.chains().iter().enumerate() {
            hmf_common::debug!(
                "cluster({}) added chain({}) ploidy({}) with {} linked pairs:",
                self.cluster_id,
                chain.id(),
                format_jcn(chain.jcn()),
                chain.link_count()
            );
            chain.set_id(i as i32);
        }
    }

    fn build_chains(&mut self, assembled_links_only: bool) {
        let svs = self.sv_list.clone();
        self.populate_sv_jcn_map(&svs, self.has_replication);
        let links = self.assembled_links.clone();
        self.add_assembly_links_to_chains(&links, self.has_replication);
        if assembled_links_only {
            return;
        }
        if self.use_allele_jcns && self.has_replication {
            self.jcn_limits.determine_breakend_jcns();
        }
        self.determine_possible_links();
        let mut iterations_without_new_links = 0;
        loop {
            self.clear_skipped_state();
            let last_added_index = self.link_index;
            let mut proposed = self.find_proposed_links();
            if proposed.is_empty() {
                if !self.pair_skipped {
                    break;
                }
            } else {
                self.process_proposed_links(&mut proposed);
                if self.run_validation {
                    self.check_chains();
                }
                if !self.is_valid() {
                    return;
                }
            }
            if last_added_index == self.link_index {
                iterations_without_new_links += 1;
                if iterations_without_new_links >= 50 {
                    hmf_common::info!(
                        "cluster({}) {} iterations without adding a link, skipped pairs: closing({}) ploidyMismatch({})",
                        self.cluster_id,
                        iterations_without_new_links,
                        self.skipped_pair_count(LinkSkipType::CLOSING),
                        self.skipped_pair_count(LinkSkipType::JCN_MISMATCH)
                    );
                    break;
                }
            } else {
                iterations_without_new_links = 0;
            }
        }
        if self.skipped_pair_count(LinkSkipType::JCN_MISMATCH) > 0 {
            let mut mismatches = 0;
            let mut overlaps = 0;
            for (pair, t) in &self.skipped_pairs {
                if *t != LinkSkipType::JCN_MISMATCH {
                    continue;
                }
                mismatches += 1;
                if jcn_match_breakends(&pair.first_breakend(), &pair.second_breakend()) {
                    overlaps += 1;
                }
            }
            if overlaps > 0 {
                hmf_common::info!(
                    "sample({}) cluster({}) ploidy skips({}) with overlaps({}) chainCount({})",
                    self.sample_id,
                    self.cluster_id,
                    mismatches,
                    overlaps,
                    self.chains.len()
                );
            }
        }
        if self.chains.len() < 50 {
            let next_id = self.next_chain_id;
            reconcile_chains(&mut self.chains, false, next_id, true);
        }
        self.check_double_minute_chains();
    }

    fn cluster_chr_breakend_index(&self, breakend: &Be) -> i32 {
        if !self.is_cluster_subset {
            return breakend.cluster_chr_pos_index();
        }
        self.subset_breakend_cluster_index_map.get(breakend).copied().unwrap_or(-1)
    }

    fn determine_possible_links(&mut self) {
        let entries = self.chr_breakend_map.clone();
        for (chromosome, list) in &entries {
            if list.len() < 2 {
                continue;
            }
            for i in 0..list.len() - 1 {
                let lower = &list[i];
                if below_jcn_threshold(lower.sv()) {
                    continue;
                }
                if !lower.reverse_orient() {
                    continue;
                }
                if self.already_linked_breakend(lower) {
                    continue;
                }
                let lower_sv = lower.sv().clone();
                let lower_valid_ap = self.use_allele_jcns
                    && has_valid_allele_jcn_data(self.cluster_chr_breakend_index(lower), self.jcn_limits.chr_allele_jcns.get(chromosome));
                let lower_jcn = self.unlinked_breakend_count(lower);
                let mut skipped_non_assembled_index: i32 = -1;
                for j in (i + 1)..list.len() {
                    let upper = &list[j];
                    if below_jcn_threshold(upper.sv()) {
                        continue;
                    }
                    if skipped_non_assembled_index == -1 && !upper.is_assembled_link() {
                        skipped_non_assembled_index = j as i32;
                    }
                    if !upper.forward_orient() {
                        continue;
                    }
                    if upper.sv() == lower.sv() {
                        continue;
                    }
                    if self.already_linked_breakend(upper) {
                        continue;
                    }
                    let distance = (upper.position() - lower.position() + 1) as i64;
                    let min_ti = get_min_templated_insertion_length(lower, upper);
                    if distance < min_ti as i64 {
                        continue;
                    }
                    let upper_sv = upper.sv().clone();
                    if (lower_sv.sv_type() == INF || upper_sv.sv_type() == INF)
                        && distance <= 1
                        && (!jcn_exceeds_major_allele_jcn(lower) || !jcn_exceeds_major_allele_jcn(upper))
                    {
                        continue;
                    }
                    let new_pair = Pair::from_breakends(lower, upper);
                    let mut adjacent = false;
                    if j == i + 1 {
                        if lower.db_link().is_none() || lower.db_link() != upper.db_link() {
                            adjacent = true;
                        }
                    } else if j == i + 2 {
                        let middle = &list[i + 1];
                        if (lower.db_link().is_some() && lower.db_link() == middle.db_link())
                            || (upper.db_link().is_some() && upper.db_link() == middle.db_link())
                        {
                            adjacent = true;
                        }
                    }
                    if adjacent {
                        self.adjacent_pairs.push(new_pair.clone());
                        if copy_numbers_equal(lower_jcn, self.unlinked_breakend_count(upper)) {
                            self.adjacent_matching_pairs.push(new_pair.clone());
                        }
                    }
                    self.sv_breakend_possible_links.entry(lower.clone()).or_default().push(new_pair.clone());
                    self.sv_breakend_possible_links.entry(upper.clone()).or_default().insert(0, new_pair.clone());
                    if skipped_non_assembled_index == -1 || skipped_non_assembled_index == j as i32 {
                        if !lower.sv().is_foldback() {
                            self.check_is_complex_dup_sv(lower, upper);
                        }
                        if !upper.sv().is_foldback() {
                            self.check_is_complex_dup_sv(upper, lower);
                        }
                    }
                    if lower_valid_ap {
                        let allele = self.jcn_limits.chr_allele_jcns.get(chromosome);
                        let bi = self.cluster_chr_breakend_index(upper);
                        if has_valid_allele_jcn_data(bi, allele) {
                            let cluster_ap = allele.unwrap()[bi as usize].cluster_jcn();
                            if cluster_ap < CLUSTER_ALLELE_JCN_MIN {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    fn check_is_complex_dup_sv(&mut self, lower_jcn_breakend: &Be, higher_jcn_breakend: &Be) {
        let var = lower_jcn_breakend.sv().clone();
        if var.jcn() < 1.0 {
            return;
        }
        if var.is_sgl_breakend() || var.sv_type() == DEL || higher_jcn_breakend.sv().is_sgl_breakend() {
            return;
        }
        if self.complex_dup_candidates.contains_key(&var) {
            return;
        }
        if !self.sv_connections.contains_key(&var.id()) {
            return;
        }
        let other_sv = higher_jcn_breakend.sv().clone();
        if var.jcn() > other_sv.jcn() || var.jcn_min() * 2.0 > other_sv.jcn_max() {
            return;
        }
        if !jcn_match_for_splits(var.jcn(), var.jcn_uncertainty(), other_sv.jcn(), other_sv.jcn_uncertainty()) {
            return;
        }
        let other_breakend = var.be(!lower_jcn_breakend.uses_start());
        let Some(list) = self.breakend_list(other_breakend.chromosome()).cloned() else { return };
        let traverse_up = other_breakend.reverse_orient();
        let mut index = self.cluster_chr_breakend_index(&other_breakend);
        loop {
            index += if traverse_up { 1 } else { -1 };
            if index < 0 || index as usize >= list.len() {
                break;
            }
            let breakend = &list[index as usize];
            if breakend == lower_jcn_breakend {
                break;
            }
            if breakend.is_assembled_link() {
                continue;
            }
            if breakend.orientation() == other_breakend.orientation() {
                break;
            }
            let d = (breakend.position() - other_breakend.position()).abs() + 1;
            if d < get_min_templated_insertion_length(breakend, &other_breakend) {
                continue;
            }
            let other_sv2 = breakend.sv().clone();
            if !jcn_match_for_splits(var.jcn(), var.jcn_uncertainty(), other_sv2.jcn(), other_sv2.jcn_uncertainty()) {
                return;
            }
            let links = vec![Pair::from_breakends(lower_jcn_breakend, higher_jcn_breakend), Pair::from_breakends(&other_breakend, breakend)];
            hmf_common::debug!("identified complex dup({} {}) ploidy({})", var.pos_id(), var.sv_type(), format_jcn(var.jcn()));
            self.complex_dup_candidates.insert(var.clone(), links);
            break;
        }
    }

    fn already_linked_breakend(&self, breakend: &Be) -> bool {
        breakend.is_assembled_link() && self.unlinked_breakend_count(breakend) == 0.0
    }

    fn check_double_minute_chains(&mut self) {
        if self.double_minute_svs.is_empty() {
            return;
        }
        let dm_chains = self.double_minute_svs[0].cl().double_minute_chains();
        let mut unmatched = dm_chains.clone();
        for chain in self.chains.clone() {
            let matches_dm = dm_chains.iter().any(|d| chain.has_matching_svs(d));
            if matches_dm && chain.could_close_chain() {
                let cs = chain.open_breakend(true).unwrap();
                let ce = chain.open_breakend(false).unwrap();
                if ce.sv() == cs.sv() && cs != ce {
                    chain.close_chain();
                } else {
                    chain.close_chain_with(LR_METHOD_DM_CLOSE, self.link_index);
                    hmf_common::debug!("cluster({}) closed DM chain({})", self.cluster_id, chain.id());
                }
            }
            for d in &dm_chains {
                if identical_chain(&chain, d, false, true) {
                    chain.set_double_minute(true);
                    remove_item(&mut unmatched, d);
                    break;
                }
            }
        }
        if !unmatched.is_empty() {
            hmf_common::debug!("cluster({}) added pre-formed DM chain", self.cluster_id);
            self.chains.extend(unmatched);
        }
    }

    fn check_chains(&mut self) {
        for chain in &self.chains {
            if !chain_is_valid(chain) {
                hmf_common::error!("cluster({}) has invalid chain({})", self.cluster_id, chain.id());
                self.finder_is_valid = false;
            }
        }
        for i in 0..self.chains.len().saturating_sub(1) {
            let c1 = self.chains[i].linked_pairs();
            for j in i + 1..self.chains.len() {
                for pair in self.chains[j].linked_pairs() {
                    if c1.iter().any(|x| x == &pair) {
                        hmf_common::error!("cluster({}) chain({}) and chain({}) share pair({})", self.cluster_id, self.chains[i].id(), self.chains[j].id(), pair);
                        self.finder_is_valid = false;
                    }
                }
            }
        }
    }
}

pub fn below_jcn_threshold(var: &SvVarData) -> bool {
    var.jcn() <= MIN_CHAINING_JCN_LEVEL
}

// ChainLinkAllocator
impl ChainFinder {
    fn allocator_initialise(&mut self, cluster_id: i32) {
        self.cluster_id = cluster_id;
        self.allocator_is_valid = true;
        self.link_index = 0;
        self.pair_skipped = false;
        self.chains_split = false;
        self.next_chain_id = 0;
        self.unique_pairs.clear();
        self.skipped_pairs.clear();
        self.sv_connections.clear();
        self.sv_completed_connections.clear();
    }

    fn clear_skipped_state(&mut self) {
        self.pair_skipped = false;
        self.chains_split = false;
    }

    fn conn(&self, var: &SvVarData) -> Option<ChainStateRef> {
        self.sv_connections.get(&var.id()).cloned()
    }

    fn populate_sv_jcn_map(&mut self, sv_list: &[Sv], has_replication: bool) {
        let uniform = if !has_replication {
            let s = hmf_common::java_math::double_stream_sum(sv_list.iter().map(|x| x.jcn()));
            Some(s / sv_list.len() as f64)
        } else {
            None
        };
        for var in sv_list {
            if below_jcn_threshold(var) {
                hmf_common::debug!(
                    "cluster({}) skipping SV({}) with low JCN({} min={} max={})",
                    self.cluster_id,
                    var.id(),
                    format_jcn(var.jcn()),
                    format_jcn(var.jcn_min()),
                    format_jcn(var.jcn_max())
                );
                continue;
            }
            self.sv_connections.insert(var.id(), Rc::new(RefCell::new(ChainState::new(var, uniform))));
        }
    }

    fn add_assembly_links_to_chains(&mut self, assembled_links: &[Pair], has_replication: bool) {
        if assembled_links.is_empty() {
            return;
        }
        if !has_replication {
            for pair in assembled_links {
                let mut pl = ProposedLinks::single(pair, ChainingRule::ASSEMBLY);
                let (b1, b2) = (pair.breakend(true), pair.breakend(false));
                pl.add_breakend_ploidies(&b1, self.unlinked_breakend_count(&b1), &b2, self.unlinked_breakend_count(&b2));
                if !pl.is_valid() {
                    hmf_common::debug!("cluster({}) skipping assembled link({}) with low JCN", self.cluster_id, pl);
                    continue;
                }
                self.add_links(&mut pl);
            }
            return;
        }
        let mut single_link_breakends: Vec<Be> = Vec::new();
        let mut both_multi_pairs: Vec<Pair> = Vec::new();
        let mut assembly_links: Vec<Pair> = assembled_links.to_vec();
        let mut index = 0;
        while index < assembly_links.len() {
            let pair = assembly_links[index].clone();
            let (fb, sb) = (pair.first_breakend(), pair.second_breakend());
            let first_single = fb.sv().max_assembled_breakend() <= 1;
            let second_single = sb.sv().max_assembled_breakend() <= 1;
            let first_jcn = self.unlinked_breakend_count(&fb);
            let second_jcn = self.unlinked_breakend_count(&sb);
            if first_jcn == 0.0 || second_jcn == 0.0 {
                hmf_common::debug!(
                    "cluster({}) skipping assembled pair({}) with low ploidy({} & {})",
                    self.cluster_id,
                    pair,
                    format_jcn(first_jcn),
                    format_jcn(second_jcn)
                );
                assembly_links.remove(index);
                continue;
            }
            if first_single && second_single {
                let mut pl = ProposedLinks::single(&pair, ChainingRule::ASSEMBLY);
                pl.add_breakend_ploidies(&fb, first_jcn, &sb, second_jcn);
                self.add_links(&mut pl);
                assembly_links.remove(index);
                continue;
            }
            index += 1;
            if first_single {
                single_link_breakends.push(fb.clone());
            } else if second_single {
                single_link_breakends.push(sb.clone());
            }
            if !first_single && !second_single {
                both_multi_pairs.push(pair.clone());
            }
        }
        index = 0;
        while index < assembly_links.len() {
            let pair = assembly_links[index].clone();
            if !contains(&both_multi_pairs, &pair) {
                let (fb, sb) = (pair.first_breakend(), pair.second_breakend());
                let first_single = contains(&single_link_breakends, &fb);
                let second_single = contains(&single_link_breakends, &sb);
                let first_jcn = if first_single { self.unlinked_breakend_count(&fb) } else { self.max_unlinked_breakend_count(&fb) };
                let second_jcn = if second_single { self.unlinked_breakend_count(&sb) } else { self.max_unlinked_breakend_count(&sb) };
                if first_jcn == 0.0 || second_jcn == 0.0 {
                    assembly_links.remove(index);
                    continue;
                }
                let mut pl = ProposedLinks::single(&pair, ChainingRule::ASSEMBLY);
                pl.add_breakend_ploidies(&fb, first_jcn, &sb, second_jcn);
                if !first_single && pl.exhaust_breakend(&fb) {
                    pl.override_breakend_jcn_matched(&fb, false);
                } else if !second_single && pl.exhaust_breakend(&sb) {
                    pl.override_breakend_jcn_matched(&sb, false);
                }
                self.add_links(&mut pl);
                assembly_links.remove(index);
                continue;
            }
            index += 1;
        }
        index = 0;
        let mut linked_pair = true;
        let mut iterations = 0;
        while index < both_multi_pairs.len() && !both_multi_pairs.is_empty() {
            iterations += 1;
            let pair = both_multi_pairs[index].clone();
            let (fb, sb) = (pair.first_breakend(), pair.second_breakend());
            let first_remaining =
                fb.sv().assembled_linked_pairs(fb.uses_start()).iter().filter(|x| contains(&both_multi_pairs, x)).count();
            let second_remaining =
                sb.sv().assembled_linked_pairs(sb.uses_start()).iter().filter(|x| contains(&both_multi_pairs, x)).count();
            if first_remaining == 0 || second_remaining == 0 {
                hmf_common::error!(
                    "cluster({}) pair({}) unexpected remaining assembly link count: first({}) second({})",
                    self.cluster_id,
                    pair,
                    first_remaining,
                    second_remaining
                );
                break;
            }
            let first_jcn = self.max_unlinked_breakend_count(&fb);
            let second_jcn = self.max_unlinked_breakend_count(&sb);
            if first_jcn == 0.0 || second_jcn == 0.0 {
                both_multi_pairs.remove(index);
                continue;
            }
            let mut pl = ProposedLinks::single(&pair, ChainingRule::ASSEMBLY);
            if first_remaining == 1 || second_remaining == 1 {
                pl.add_breakend_ploidies(&fb, first_jcn, &sb, second_jcn);
            } else if !linked_pair {
                pl.add_breakend_ploidies(&fb, first_jcn / first_remaining as f64, &sb, second_jcn / second_remaining as f64);
            } else {
                index += 1;
                if index >= both_multi_pairs.len() {
                    index = 0;
                }
                if iterations > both_multi_pairs.len() * 3 {
                    hmf_common::debug!("cluster({}) assembly multi-connection breakends missed", self.cluster_id);
                    break;
                }
                continue;
            }
            if first_remaining > 1 && pl.exhaust_breakend(&fb) {
                pl.override_breakend_jcn_matched(&fb, false);
            }
            if second_remaining > 1 && pl.exhaust_breakend(&sb) {
                pl.override_breakend_jcn_matched(&sb, false);
            }
            self.add_links(&mut pl);
            linked_pair = true;
            both_multi_pairs.remove(index);
        }
        let _ = linked_pair;
        if !self.chains.is_empty() {
            hmf_common::debug!("created {} partial chains from {} assembly links", self.chains.len(), assembled_links.len());
        }
    }

    fn process_proposed_links(&mut self, list: &mut Vec<ProposedLinks>) {
        let mut link_added = false;
        while !list.is_empty() {
            if let Some(t) = list[0].target_chain() {
                if !contains(&self.chains, t) {
                    break;
                }
            }
            let mut pl = list.remove(0);
            if !pl.is_valid() {
                hmf_common::debug!("cluster({}) skipping invalid proposed links: {}", self.cluster_id, pl);
                continue;
            }
            link_added |= self.add_links(&mut pl);
            if !self.allocator_is_valid {
                return;
            }
            if pl.multi_connection() {
                break;
            }
        }
        if link_added {
            if self.chains_split {
                self.skipped_pairs.clear();
            } else {
                self.skipped_pairs.retain(|_, v| *v == LinkSkipType::JCN_MISMATCH);
            }
        }
    }

    fn add_links(&mut self, pl: &mut ProposedLinks) -> bool {
        let new_pair = pl.links[0].clone();
        let top_rule = pl.top_rule().to_string();
        for l in &pl.links {
            l.set_link_reason(&top_rule, self.link_index);
        }
        let mut add_links_to_new_chain = true;
        let mut reconcile = false;
        if pl.target_chain().is_some() {
            add_links_to_new_chain = false;
            self.add_complex_links_to_existing_chain(pl);
            reconcile = true;
        } else if pl.multi_connection() {
            // managed by the rule-finder
        } else {
            let mut target_chain: Option<Chain> = None;
            let mut pair_linked_on_first = false;
            let mut add_to_start = false;
            let mut matches_chain_jcn = false;
            let mut new_sv_jcn = 0.0;
            let allow_chain_splits = self.allow_dm_chain_splits(&pl.links[0]);
            let pair_jcn_matched = jcn_match_breakends(&new_pair.first_breakend(), &new_pair.second_breakend());
            let mut required_chains: [Option<Chain>; 2] = [None, None];
            let mut all_chains: Vec<Chain> = Vec::new();
            for s in [SE_START, SE_END] {
                let pb = new_pair.breakend_se(s);
                let bj = self.breakend_jcn_data(&pb);
                if bj.exhausted_in_chain() {
                    let mc = bj.max_jcn_chain.clone();
                    all_chains.push(mc.clone().expect("max JCN chain"));
                    required_chains[s] = mc;
                    continue;
                }
                for c in &bj.chains {
                    if !contains(&all_chains, c) {
                        all_chains.push(c.clone());
                    }
                }
            }
            if required_chains[0].is_some() || required_chains[1].is_some() {
                if let (Some(rs), Some(re)) = (&required_chains[0], &required_chains[1]) {
                    if !jcn_match(rs.jcn(), rs.jcn_uncertainty(), re.jcn(), re.jcn_uncertainty()) && !pair_jcn_matched {
                        if !allow_chain_splits {
                            self.add_skipped_pair(&new_pair, LinkSkipType::JCN_MISMATCH);
                            return false;
                        }
                    } else {
                        matches_chain_jcn = true;
                    }
                    if rs == re {
                        self.add_skipped_pair(&new_pair, LinkSkipType::CLOSING);
                        return false;
                    }
                }
                let tc = required_chains[0].clone().or_else(|| required_chains[1].clone()).unwrap();
                let (fbk, sbk) = (Some(new_pair.first_breakend()), Some(new_pair.second_breakend()));
                if tc.open_breakend(true) == fbk {
                    pair_linked_on_first = true;
                    add_to_start = true;
                } else if tc.open_breakend(true) == sbk {
                    pair_linked_on_first = false;
                    add_to_start = true;
                } else if tc.open_breakend(false) == fbk {
                    pair_linked_on_first = true;
                    add_to_start = false;
                } else if tc.open_breakend(false) == sbk {
                    pair_linked_on_first = false;
                    add_to_start = false;
                }
                let new_unc = if pair_linked_on_first { new_pair.second().jcn_uncertainty() } else { new_pair.first().jcn_uncertainty() };
                if !jcn_match(tc.jcn(), tc.jcn_uncertainty(), pl.jcn(), new_unc) && !pair_jcn_matched {
                    if !allow_chain_splits {
                        self.add_skipped_pair(&new_pair, LinkSkipType::JCN_MISMATCH);
                        return false;
                    } else {
                        matches_chain_jcn = false;
                        new_sv_jcn = pl.jcn();
                    }
                } else {
                    matches_chain_jcn = true;
                    new_sv_jcn = tc.jcn();
                }
                if matches_chain_jcn {
                    for s in [SE_START, SE_END] {
                        if let Some(rc) = &required_chains[s] {
                            let b = new_pair.breakend_se(s);
                            if !(rc.open_breakend(true).as_ref() == Some(&b) && rc.open_breakend(false).as_ref() == Some(&b)) {
                                pl.override_breakend_jcn_matched(&b, true);
                            }
                        }
                    }
                }
                target_chain = Some(tc);
            } else {
                for chain in &all_chains {
                    let mut can_add_to_start = [false, false];
                    let mut links_to_first = false;
                    let ploidy_matched = jcn_match(pl.jcn(), chain.jcn_uncertainty(), chain.jcn(), chain.jcn_uncertainty());
                    for s in [SE_START, SE_END] {
                        let Some(cb) = chain.open_breakend(is_start(s)) else { continue };
                        if cb == new_pair.first_breakend() {
                            links_to_first = true;
                        } else if cb == new_pair.second_breakend() {
                            links_to_first = false;
                        } else {
                            continue;
                        }
                        can_add_to_start[s] = true;
                    }
                    if !can_add_to_start[0] && !can_add_to_start[1] {
                        continue;
                    }
                    let could_close = if can_add_to_start[0] && can_add_to_start[1] { chain.link_would_close_chain(&new_pair) } else { false };
                    if could_close {
                        self.add_skipped_pair(&new_pair, LinkSkipType::CLOSING);
                        return false;
                    }
                    let new_breakend = if links_to_first { new_pair.second_breakend() } else { new_pair.first_breakend() };
                    if !ploidy_matched {
                        if pl.link_jcn_match() {
                            continue;
                        }
                        if target_chain.as_ref().is_some_and(|t| t.jcn() > chain.jcn()) {
                            continue;
                        }
                    } else if target_chain.is_some() && matches_chain_jcn {
                        continue;
                    }
                    target_chain = Some(chain.clone());
                    add_to_start = can_add_to_start[SE_START];
                    pair_linked_on_first = links_to_first;
                    new_sv_jcn = pl.breakend_jcn(&new_breakend);
                    if ploidy_matched {
                        matches_chain_jcn = true;
                    }
                }
            }
            if !allow_chain_splits && target_chain.as_ref().is_some_and(|t| !matches_chain_jcn && t.jcn() > pl.jcn()) {
                hmf_common::debug!(
                    "skipping targetChain({} ploidy={}) for proposedLink({}) on ploidy mismatch",
                    target_chain.as_ref().unwrap().id(),
                    format_jcn(target_chain.as_ref().unwrap().jcn()),
                    pl
                );
                target_chain = None;
            }
            if let Some(tc) = target_chain {
                add_links_to_new_chain = false;
                reconcile = true;
                self.add_links_to_existing_chain(pl, &tc, add_to_start, pair_linked_on_first, matches_chain_jcn, new_sv_jcn);
            }
        }
        if add_links_to_new_chain {
            reconcile = self.add_links_to_new_chain(pl);
        }
        if reconcile {
            reconcile_chains(&mut self.chains, false, 0, false);
        }
        self.register_new_link(pl);
        self.link_index += 1;
        true
    }

    fn add_complex_links_to_existing_chain(&mut self, pl: &ProposedLinks) {
        let target = pl.target_chain().unwrap().clone();
        let matches_chain_jcn = pl.link_jcn_match();
        let new_sv_jcn = pl.jcn();
        if !matches_chain_jcn && target.jcn() > new_sv_jcn * 2.0 {
            let new_chain = Chain::create(self.next_chain_id);
            self.next_chain_id += 1;
            self.chains.push(new_chain.clone());
            new_chain.copy_from(&target);
            if target.jcn() > new_sv_jcn * 2.0 {
                new_chain.set_jcn_data(target.jcn() - new_sv_jcn * 2.0, target.jcn_uncertainty());
                target.set_jcn_data(new_sv_jcn * 2.0, target.jcn_uncertainty());
            }
        }
        if pl.splitting_rule() == ChainingRule::FOLDBACK_SPLIT {
            if let Some(fc) = pl.foldback_chain() {
                foldback_chain_on_chain(&target, fc, &pl.links[0], &pl.links[1]);
                remove_item(&mut self.chains, fc);
            } else {
                foldback_chain_on_link(&target, &pl.links[0], &pl.links[1]);
            }
        } else {
            duplicate_chain_on_link(&target, &pl.links[0], &pl.links[1]);
        }
        let new_jcn = target.jcn() * 0.5;
        let new_unc = target.jcn_uncertainty() / UNCERTAINTY_SCALE_FACTOR;
        target.set_jcn_data(new_jcn, new_unc);
        self.chains_split = true;
    }

    fn add_links_to_existing_chain(
        &mut self, pl: &mut ProposedLinks, target: &Chain, add_to_start: bool, pair_linked_on_first: bool, matches_chain_jcn: bool,
        new_sv_jcn: f64,
    ) {
        let new_pair = pl.links[0].clone();
        let new_sv_breakend = if pair_linked_on_first { new_pair.second_breakend() } else { new_pair.first_breakend() };
        if !matches_chain_jcn && target.jcn() > pl.jcn() {
            let new_chain = Chain::create(self.next_chain_id);
            self.next_chain_id += 1;
            self.chains.push(new_chain.clone());
            new_chain.copy_from(target);
            self.chains_split = true;
            new_chain.set_jcn_data(target.jcn() - new_sv_jcn, target.jcn_uncertainty());
            target.set_jcn_data(new_sv_jcn, target.jcn_uncertainty());
        }
        if matches_chain_jcn || target.jcn() > new_sv_jcn {
            let d1 = if !pl.link_jcn_match() {
                JcnCalcData { jcn_estimate: pl.jcn(), jcn_uncertainty: new_sv_breakend.jcn_uncertainty(), valid: true }
            } else {
                JcnCalcData { jcn_estimate: pl.breakend_jcn(&new_sv_breakend), jcn_uncertainty: new_sv_breakend.jcn_uncertainty(), valid: true }
            };
            let d2 = JcnCalcData { jcn_estimate: target.jcn(), jcn_uncertainty: target.jcn_uncertainty(), valid: true };
            let p = calc_jcn_uncertainty(&d1, &d2);
            target.set_jcn_data(p.jcn_estimate, p.jcn_uncertainty);
        } else {
            pl.set_lower_jcn(target.jcn());
        }
        target.add_link(&pl.links[0], add_to_start);
    }

    fn add_links_to_new_chain(&mut self, pl: &ProposedLinks) -> bool {
        let mut reconcile = false;
        let new_pair = pl.links[0].clone();
        let mut new_chain: Option<Chain> = None;
        for pair in &pl.links {
            if let Some(nc) = &new_chain {
                if nc.can_add_linked_pair_to_start(pair) {
                    nc.add_link(pair, true);
                } else if nc.can_add_linked_pair_to_end(pair) {
                    nc.add_link(pair, false);
                } else {
                    new_chain = None;
                    reconcile = true;
                }
            }
            if new_chain.is_none() {
                let nc = Chain::create(self.next_chain_id);
                self.next_chain_id += 1;
                self.chains.push(nc.clone());
                nc.add_link(pair, true);
                let (d1, d2) = if !pl.link_jcn_match() || pl.multi_connection() {
                    (
                        JcnCalcData { jcn_estimate: pl.jcn(), jcn_uncertainty: new_pair.first().jcn_uncertainty(), valid: true },
                        JcnCalcData { jcn_estimate: pl.jcn(), jcn_uncertainty: new_pair.second().jcn_uncertainty(), valid: true },
                    )
                } else {
                    (
                        JcnCalcData {
                            jcn_estimate: pl.breakend_jcn(&new_pair.first_breakend()),
                            jcn_uncertainty: new_pair.first().jcn_uncertainty(),
                            valid: true,
                        },
                        JcnCalcData {
                            jcn_estimate: pl.breakend_jcn(&new_pair.second_breakend()),
                            jcn_uncertainty: new_pair.second().jcn_uncertainty(),
                            valid: true,
                        },
                    )
                };
                let p = calc_jcn_uncertainty(&d1, &d2);
                nc.set_jcn_data(p.jcn_estimate, p.jcn_uncertainty);
                new_chain = Some(nc);
            }
        }
        reconcile
    }

    fn register_new_link(&mut self, pl: &ProposedLinks) {
        let mut exhausted: Vec<Be> = Vec::new();
        let can_use_max_jcn = pl.top_rule() == ChainingRule::ASSEMBLY;
        for new_pair in &pl.links {
            self.jcn_limits.assign_link_jcn(new_pair, pl.jcn());
            self.skipped_pairs.remove(new_pair);
            self.remove_opposite_links(new_pair);
            for s in [SE_START, SE_END] {
                let breakend = new_pair.breakend_se(s);
                if contains(&exhausted, &breakend) {
                    continue;
                }
                let other_pair_breakend = new_pair.other_breakend(&breakend);
                let var = breakend.sv().clone();
                let conn = self.conn(&var);
                let be_is_start = breakend.uses_start();
                let Some(conn) = conn.filter(|c| !c.borrow().breakend_exhausted_vs_max(be_is_start)) else {
                    hmf_common::error!("breakend({}) breakend already exhausted with proposedLink({})", breakend, pl);
                    self.allocator_is_valid = false;
                    return;
                };
                if pl.exhaust_breakend(&breakend) {
                    let j = conn.borrow().jcn;
                    conn.borrow_mut().set(be_is_start, j);
                } else if !conn.borrow().has_connections() {
                    if pl.link_jcn_match() {
                        let j = conn.borrow().jcn;
                        conn.borrow_mut().set(be_is_start, j);
                    } else {
                        conn.borrow_mut().add(be_is_start, pl.jcn());
                    }
                } else {
                    self.update_breakend_allocated_ploidy(&conn, &breakend);
                }
                let breakend_exhausted =
                    if can_use_max_jcn { conn.borrow().breakend_exhausted_vs_max(be_is_start) } else { conn.borrow().breakend_exhausted(be_is_start) };
                if breakend_exhausted {
                    exhausted.push(breakend.clone());
                }
                conn.borrow_mut().add_connection(&other_pair_breakend, be_is_start);
            }
            if !self.matches_existing_pair(new_pair) {
                self.unique_pairs.push(new_pair.clone());
            }
        }
        for breakend in &exhausted {
            if let Some(conn) = self.conn(breakend.sv()) {
                let other_exhausted = if can_use_max_jcn {
                    conn.borrow().breakend_exhausted_vs_max(!breakend.uses_start())
                } else {
                    conn.borrow().breakend_exhausted(!breakend.uses_start())
                };
                if other_exhausted {
                    self.check_sv_complete(&conn);
                }
            }
            self.remove_possible_links(breakend);
        }
    }

    fn update_breakend_allocated_ploidy(&self, conn: &ChainStateRef, breakend: &Be) {
        let chains: Vec<&Chain> = self.chains.iter().filter(|x| x.has_sv(breakend.sv())).collect();
        let mut open_chained = 0.0;
        let mut contained = 0.0;
        for chain in chains {
            let cj = chain.jcn();
            if chain.open_breakend(true).as_ref() == Some(breakend) {
                open_chained += cj;
            }
            if chain.open_breakend(false).as_ref() == Some(breakend) {
                open_chained += cj;
            }
            contained += chain.linked_pairs().iter().filter(|x| x.has_breakend(breakend)).count() as f64 * cj;
        }
        let other_exhausted = conn.borrow().breakend_exhausted(!breakend.uses_start());
        if open_chained == 0.0 && other_exhausted {
            let j = conn.borrow().jcn;
            conn.borrow_mut().set(breakend.uses_start(), j);
        } else {
            conn.borrow_mut().set(breakend.uses_start(), contained);
        }
    }

    fn remove_possible_links(&mut self, breakend: &Be) {
        let Some(possible) = self.sv_breakend_possible_links.remove(breakend) else { return };
        for pair in &possible {
            let other = pair.other_breakend(breakend);
            let Some(others) = self.sv_breakend_possible_links.get_mut(&other) else { continue };
            remove_item(others, pair);
            if others.is_empty() {
                self.sv_breakend_possible_links.remove(&other);
            }
        }
    }

    fn remove_opposite_links(&mut self, pair: &Pair) {
        for s in [SE_START, SE_END] {
            let Some(other) = pair.breakend_se(s).other_breakend() else { continue };
            let Some(possible) = self.sv_breakend_possible_links.get_mut(&other) else { continue };
            if possible.is_empty() {
                self.sv_breakend_possible_links.remove(&other);
                continue;
            }
            let Some(other_pair_breakend) = pair.breakend(!is_start(s)).other_breakend() else { continue };
            if let Some(i) = possible.iter().position(|x| x.has_breakend(&other) && x.has_breakend(&other_pair_breakend)) {
                possible.remove(i);
                if possible.is_empty() {
                    self.sv_breakend_possible_links.remove(&other);
                }
            }
        }
    }

    fn check_sv_complete(&mut self, conn: &ChainStateRef) {
        let (done, id) = {
            let c = conn.borrow();
            (c.breakend_exhausted(true) && (c.sv.is_sgl_breakend() || c.breakend_exhausted(false)), c.sv.id())
        };
        if done {
            self.sv_connections.remove(&id);
            self.sv_completed_connections.push(conn.clone());
        }
    }

    pub fn unlinked_breakend_count(&self, breakend: &SvBreakend) -> f64 {
        let Some(conn) = self.conn(breakend.sv()) else { return 0.0 };
        let c = conn.borrow();
        if !c.breakend_exhausted(breakend.uses_start()) { c.unlinked(breakend.uses_start()) } else { 0.0 }
    }

    fn unlinked_breakend_count_chains(&self, breakend: &Be, limit_by_chains: bool) -> f64 {
        if limit_by_chains { self.breakend_jcn_data(breakend).unlinked_jcn() } else { self.unlinked_breakend_count(breakend) }
    }

    fn breakend_jcn_data(&self, breakend: &Be) -> BreakendJcn {
        let Some(conn) = self.conn(breakend.sv()) else { return BreakendJcn::simple(0.0, true) };
        let c = conn.borrow();
        let unlinked = if !c.breakend_exhausted(breakend.uses_start()) { c.unlinked(breakend.uses_start()) } else { 0.0 };
        let other_exhausted = c.breakend_exhausted(!breakend.uses_start());
        if unlinked == 0.0 {
            return BreakendJcn::simple(0.0, other_exhausted);
        }
        if !c.has_connections() {
            return BreakendJcn::simple(unlinked, other_exhausted);
        }
        let chains = self.chains_with_open_breakend(breakend);
        if chains.is_empty() {
            return BreakendJcn::simple(unlinked, other_exhausted);
        }
        if chains.len() == 1 && other_exhausted {
            let chain = chains[0].clone();
            return BreakendJcn { unchained_jcn: 0.0, total_chain_jcn: chain.jcn(), max_jcn_chain: Some(chain), chains, other_breakend_exhausted: other_exhausted };
        }
        let mut total = 0.0;
        let mut max_chain: Option<Chain> = None;
        for chain in &chains {
            if chain.open_breakend(true).as_ref() == Some(breakend) {
                total += chain.jcn();
                if max_chain.as_ref().is_none_or(|m| m.jcn() < chain.jcn()) {
                    max_chain = Some(chain.clone());
                }
            }
            if chain.open_breakend(false).as_ref() == Some(breakend) {
                total += chain.jcn();
                if max_chain.as_ref().is_none_or(|m| m.jcn() < chain.jcn()) {
                    max_chain = Some(chain.clone());
                }
            }
        }
        let unchained = jmax(unlinked - total, 0.0);
        BreakendJcn { unchained_jcn: unchained, total_chain_jcn: total, max_jcn_chain: max_chain, chains, other_breakend_exhausted: false }
    }

    fn chains_with_open_breakend(&self, breakend: &Be) -> Vec<Chain> {
        self.chains
            .iter()
            .filter(|x| x.open_breakend(true).as_ref() == Some(breakend) || x.open_breakend(false).as_ref() == Some(breakend))
            .cloned()
            .collect()
    }

    fn max_unlinked_breakend_count(&self, breakend: &Be) -> f64 {
        let Some(conn) = self.conn(breakend.sv()) else { return 0.0 };
        let c = conn.borrow();
        if !c.breakend_exhausted(breakend.uses_start()) {
            c.unlinked(breakend.uses_start())
        } else if !c.breakend_exhausted_vs_max(breakend.uses_start()) {
            c.max_unlinked(breakend.uses_start())
        } else {
            0.0
        }
    }

    fn matches_existing_pair(&self, pair: &Pair) -> bool {
        self.unique_pairs.iter().any(|x| pair.matches(x))
    }
    fn has_skipped_pairs(&self, pair: &Pair) -> bool {
        self.skipped_pairs.keys().any(|x| x.matches(pair))
    }
    fn skipped_pair_count(&self, t: LinkSkipType) -> usize {
        self.skipped_pairs.values().filter(|x| **x == t).count()
    }
    fn add_skipped_pair(&mut self, pair: &Pair, t: LinkSkipType) {
        if self.has_skipped_pairs(pair) {
            return;
        }
        self.pair_skipped = true;
        self.skipped_pairs.insert(pair.clone(), t);
    }
    fn remove_skipped_pairs(&self, proposed: &mut Vec<ProposedLinks>) {
        if proposed.is_empty() || self.skipped_pairs.is_empty() {
            return;
        }
        proposed.retain(|p| !p.links.iter().any(|x| self.has_skipped_pairs(x)));
    }
    fn allow_dm_chain_splits(&self, pair: &Pair) -> bool {
        if self.double_minute_svs.is_empty() {
            return false;
        }
        let f = contains(&self.double_minute_svs, &pair.first());
        let s = contains(&self.double_minute_svs, &pair.second());
        f != s
    }
}

// ChainRuleSelector
impl ChainFinder {
    fn rule_selector_initialise(&mut self, has_replication: bool) {
        use ChainingRule::*;
        self.rules_to_apply.clear();
        if has_replication {
            self.rules_to_apply.extend([FOLDBACK_SPLIT, COMP_DUP_SPLIT, ONLY, JCN_MATCH, ADJACENT, JCN_MAX]);
        } else {
            self.rules_to_apply.extend([ONLY, ADJACENT_MATCH]);
        }
        self.rules_to_apply.push(NEAREST);
        self.foldbacks_initialised = false;
        self.foldback_breakend_pairs.clear();
    }

    fn find_proposed_links(&mut self) -> Vec<ProposedLinks> {
        use ChainingRule::*;
        let mut proposed: Vec<ProposedLinks> = Vec::new();
        let rules = self.rules_to_apply.clone();
        for rule in rules {
            proposed = match rule {
                ONLY => self.find_single_option_pairs(proposed),
                FOLDBACK_SPLIT => self.find_foldback_links(proposed),
                COMP_DUP_SPLIT => self.find_complex_dups(proposed),
                JCN_MATCH => self.find_jcn_match_pairs(proposed),
                ADJACENT => self.find_adjacent_pairs(proposed),
                ADJACENT_MATCH => self.find_adjacent_matching_pairs(proposed),
                JCN_MAX => self.find_highest_jcn(proposed),
                NEAREST => self.find_nearest(proposed),
                _ => proposed,
            };
            if proposed.is_empty() {
                continue;
            }
            self.remove_skipped_pairs(&mut proposed);
            let has_support = self.check_cluster_jcn_support(&mut proposed);
            cull_by_priority(&mut proposed);
            if proposed.len() == 1 && has_support {
                return proposed;
            }
        }
        proposed
    }

    fn find_single_option_pairs(&self, mut proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        if self.sv_breakend_possible_links.is_empty() {
            return proposed;
        }
        for pl in proposed.iter_mut() {
            let mut single = false;
            for pair in pl.links.clone() {
                for s in [SE_START, SE_END] {
                    let b = pair.breakend(is_start(s));
                    let Some(possible) = self.sv_breakend_possible_links.get(&b) else {
                        hmf_common::error!("breakend({}) has no possible pairs, from proposedLink: {}", b, pl);
                        continue;
                    };
                    if possible.len() == 1 {
                        single = true;
                        break;
                    }
                }
                if single {
                    pl.add_rule(ChainingRule::ONLY);
                    break;
                }
            }
        }
        // BreakendComparator: by SV id, the start breakend first
        let mut sorted: Vec<Be> = self.sv_breakend_possible_links.keys().cloned().collect();
        sorted.sort_by(|a, b| (a.sv().id(), !a.uses_start()).cmp(&(b.sv().id(), !b.uses_start())));
        for limiting in &sorted {
            let pairs = &self.sv_breakend_possible_links[limiting];
            if pairs.is_empty() {
                hmf_common::warn!("breakend({}) has no possibles left, should be purged", limiting);
                continue;
            }
            if pairs.len() >= 2 {
                continue;
            }
            let new_pair = pairs[0].clone();
            if self.has_skipped_pairs(&new_pair) {
                continue;
            }
            if proposed.iter().any(|x| x.links[0] == new_pair) {
                continue;
            }
            let jcn_first = self.unlinked_breakend_count_chains(&new_pair.first_breakend(), true);
            let jcn_second = self.unlinked_breakend_count_chains(&new_pair.second_breakend(), true);
            if jcn_first == 0.0 || jcn_second == 0.0 {
                continue;
            }
            let mut pl = ProposedLinks::single(&new_pair, ChainingRule::ONLY);
            pl.add_breakend_ploidies(&new_pair.first_breakend(), jcn_first, &new_pair.second_breakend(), jcn_second);
            let mut add_new = true;
            let mut index = 0;
            while index < proposed.len() {
                let other_pair = proposed[index].links[0].clone();
                if !other_pair.has_link_clash(&new_pair) && !other_pair.opposite_match(&new_pair) {
                    index += 1;
                    continue;
                }
                if self.has_replication {
                    if copy_numbers_equal(proposed[index].jcn(), pl.jcn()) {
                        if pl.shortest_link_distance() > proposed[index].shortest_link_distance() {
                            add_new = false;
                            break;
                        }
                        proposed.remove(index);
                    }
                    // Java increments the index after removing too, skipping the next entry
                    index += 1;
                } else {
                    if pl.shortest_link_distance() > proposed[index].shortest_link_distance() {
                        add_new = false;
                        break;
                    }
                    proposed.remove(index);
                }
            }
            if add_new {
                proposed.push(pl);
            }
        }
        proposed
    }

    fn update_foldback_breakends(&mut self) {
        let mut existing_chained: Vec<FoldbackBreakendPair> = Vec::new();
        if !self.foldbacks_initialised {
            self.foldbacks_initialised = true;
            self.foldback_breakend_pairs.clear();
            for foldback in self.foldbacks.clone() {
                if foldback.is_single_breakend_foldback() || foldback.is_chained_foldback() {
                    continue;
                }
                let fs = foldback.be(true);
                let fe = foldback.be(false);
                let jcn = jmin(self.unlinked_breakend_count(&fs), self.unlinked_breakend_count(&fe));
                if jcn == 0.0 {
                    continue;
                }
                add_by_ploidy(&mut self.foldback_breakend_pairs, FoldbackBreakendPair { breakend_start: fs, breakend_end: fe, ploidy: jcn, chain: None });
            }
        } else {
            let mut index = 0;
            while index < self.foldback_breakend_pairs.len() {
                let fb = self.foldback_breakend_pairs[index].clone();
                let bp_start = self.breakend_jcn_data(&fb.breakend_start);
                let bp_end = self.breakend_jcn_data(&fb.breakend_end);
                let min_unlinked;
                let mut chain_has_simple_foldback = false;
                if fb.is_chained() {
                    min_unlinked = jmin(bp_start.unlinked_jcn(), bp_end.unlinked_jcn());
                } else {
                    min_unlinked = jmin(bp_start.unchained_jcn, bp_end.unchained_jcn);
                    if !bp_start.chains.is_empty() && bp_end.chains.is_empty() {
                        chain_has_simple_foldback = bp_start.chains.iter().any(|x| chain_has_foldback(x, fb.breakend_start.sv()));
                    }
                }
                if min_unlinked < MIN_CHAINING_JCN_LEVEL || chain_has_simple_foldback {
                    self.foldback_breakend_pairs.remove(index);
                    continue;
                }
                if fb.is_chained() {
                    existing_chained.push(fb);
                } else {
                    self.foldback_breakend_pairs[index].ploidy = min_unlinked;
                }
                index += 1;
            }
        }
        for chain in self.chains.clone() {
            let (Some(cs), Some(ce)) = (chain.open_breakend(true), chain.open_breakend(false)) else { continue };
            if cs.chromosome() != ce.chromosome() || cs.orientation() != ce.orientation() {
                continue;
            }
            if !self.sv_breakend_possible_links.contains_key(&cs) || !self.sv_breakend_possible_links.contains_key(&ce) {
                continue;
            }
            let si = cs.cluster_chr_pos_index();
            let ei = ce.cluster_chr_pos_index();
            if cs != ce && (si - ei).abs() != 1 {
                if (si - ei).abs() > 2 {
                    continue;
                }
                let start_is_front = (cs.position() < ce.position()) == cs.forward_orient();
                let front = if start_is_front { &cs } else { &ce };
                if front.sv().db_link(front.uses_start()).is_some_and(|d| d.length() > 0) {
                    continue;
                }
            }
            let fb = FoldbackBreakendPair { breakend_start: cs, breakend_end: ce, ploidy: chain.jcn(), chain: Some(chain.clone()) };
            let already = self.foldback_breakend_pairs.iter().any(|x| x.matches(&fb));
            if let Some(i) = existing_chained.iter().position(|x| x.matches(&fb)) {
                existing_chained.remove(i);
            }
            if !already {
                add_by_ploidy(&mut self.foldback_breakend_pairs, fb);
            } else if let Some(i) = self.foldback_breakend_pairs.iter().position(|x| x.matches(&fb)) {
                self.foldback_breakend_pairs[i] = fb;
            }
        }
        for x in &existing_chained {
            if let Some(i) = self.foldback_breakend_pairs.iter().position(|p| p.matches(x)) {
                self.foldback_breakend_pairs.remove(i);
            }
        }
    }

    fn find_foldback_links(&mut self, proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        const A: i32 = 5;
        const B: i32 = 3;
        const C: i32 = 2;
        const D: i32 = 1;
        self.update_foldback_breakends();
        if self.foldback_breakend_pairs.is_empty() {
            return proposed;
        }
        let mut new_proposed: Vec<ProposedLinks> = Vec::new();
        let mut link_score = 0;
        let mut last_foldback_ploidy = 0.0;
        let fb_pairs = self.foldback_breakend_pairs.clone();
        for (fbi, fb) in fb_pairs.iter().enumerate() {
            if !new_proposed.is_empty() && !copy_numbers_equal(fb.ploidy, last_foldback_ploidy) {
                break;
            }
            last_foldback_ploidy = fb.ploidy;
            let fs = &fb.breakend_start;
            let fe = &fb.breakend_end;
            let foldback = fs.sv().clone();
            let foldback_chain = fb.chain.clone();
            let foldback_ploidy = fb.ploidy;
            let foldback_uncertainty = if let Some(c) = &foldback_chain { c.jcn_uncertainty() } else { 0.0 };
            let (Some(on_start), Some(on_end)) = (self.sv_breakend_possible_links.get(fs), self.sv_breakend_possible_links.get(fe)) else {
                continue;
            };
            let mut pairs_on_start = on_start.clone();
            let pairs_on_end = on_end.clone();
            self.cull_dual_option_pairs(fs, &mut pairs_on_start);
            for pair_start in &pairs_on_start {
                let non_fb_var = pair_start.other_sv(&foldback);
                let other = pair_start.other_breakend(fs);
                let Some(pair_end) = pairs_on_end.iter().find(|x| x.other_breakend(fe) == other).cloned() else { continue };
                let non_fb = self.breakend_jcn_data(&other);
                let non_fb_ploidy = non_fb.unlinked_jcn();
                let allow_splits = !(copy_numbers_equal(foldback_ploidy, non_fb_ploidy) && foldback_ploidy < 1.0 && non_fb_ploidy < 1.0);
                if allow_splits
                    && jcn_match_for_splits(foldback_ploidy, foldback_uncertainty, non_fb_ploidy, non_fb_var.jcn_uncertainty())
                    && !non_fb.multi_connections()
                {
                    if link_score < A {
                        new_proposed.clear();
                        link_score = A;
                    }
                    let mut pl = ProposedLinks::multi(
                        vec![pair_start.clone(), pair_end.clone()],
                        ChainingRule::FOLDBACK_SPLIT,
                        non_fb.max_jcn_chain.clone(),
                        foldback_chain.clone(),
                    );
                    pl.add_foldback_breakends(fs, fe, foldback_ploidy, &other, non_fb_ploidy, non_fb_var.jcn_uncertainty());
                    new_proposed.push(pl);
                    continue;
                }
                if link_score >= A {
                    continue;
                }
                let other_fb = fb_pairs
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != fbi)
                    .map(|(_, x)| x)
                    .find(|x| x.breakend_start == other || x.breakend_end == other)
                    .cloned();
                if !non_fb.multi_connections() && jcn_match(foldback_ploidy, foldback_uncertainty, non_fb_ploidy, non_fb_var.jcn_uncertainty()) {
                    if link_score < B {
                        new_proposed.clear();
                        link_score = B;
                    }
                    // Java: getSvConnections().get(foldback), which exists for an unexhausted foldback
                    let fb_conn = self.conn(&foldback).expect("foldback chain state");
                    let fb_breakend = {
                        let c = fb_conn.borrow();
                        if c.unlinked(fs.uses_start()) > c.unlinked(fe.uses_start()) { fs.clone() } else { fe.clone() }
                    };
                    let mut pl = ProposedLinks::single(&Pair::from_breakends(&fb_breakend, &other), ChainingRule::FOLDBACK);
                    pl.add_breakend_ploidies(&fb_breakend, foldback_ploidy, &other, non_fb_ploidy);
                    new_proposed.push(pl);
                    continue;
                }
                if link_score >= B {
                    continue;
                }
                if let Some(ofb) = &other_fb {
                    if !jcn_match_for_splits(foldback_ploidy, foldback_uncertainty, non_fb_ploidy, non_fb_var.jcn_uncertainty())
                        && non_fb_ploidy > foldback_ploidy * 2.0
                    {
                        if link_score < C {
                            new_proposed.clear();
                            link_score = C;
                        }
                        let target = non_fb.max_jcn_chain.clone().filter(|m| non_fb_ploidy == m.jcn());
                        let mut pl = ProposedLinks::multi(
                            vec![pair_start.clone(), pair_end.clone()],
                            ChainingRule::FOLDBACK_SPLIT,
                            target,
                            foldback_chain.clone(),
                        );
                        pl.add_foldback_breakends(fs, fe, foldback_ploidy, &other, non_fb_ploidy, non_fb_var.jcn_uncertainty());
                        new_proposed.push(pl);
                        continue;
                    }
                    if allow_splits && jcn_match_for_splits(non_fb_ploidy, non_fb_var.jcn_uncertainty(), foldback_ploidy, foldback_uncertainty) {
                        link_score = D;
                        let mut pl = ProposedLinks::multi(
                            vec![Pair::from_breakends(&ofb.breakend_start, fs), Pair::from_breakends(&ofb.breakend_end, fs)],
                            ChainingRule::FOLDBACK_SPLIT,
                            foldback_chain.clone(),
                            ofb.chain.clone(),
                        );
                        pl.add_foldback_breakends(&ofb.breakend_start, &ofb.breakend_end, ofb.ploidy, fs, foldback_ploidy, foldback_uncertainty);
                        new_proposed.push(pl);
                        continue;
                    }
                }
            }
        }
        let _ = link_score;
        new_proposed
    }

    fn cull_dual_option_pairs(&self, source: &Be, pairs: &mut Vec<Pair>) {
        let mut i = 0;
        while i < pairs.len() {
            let breakend = pairs[i].other_breakend(source);
            let mut removed = false;
            for j in i + 1..pairs.len() {
                let breakend2 = pairs[j].other_breakend(source);
                if breakend.sv() != breakend2.sv() {
                    continue;
                }
                let Some(conn) = self.conn(breakend.sv()) else { continue };
                let c = conn.borrow();
                let (u1, u2) = (c.unlinked(breakend.uses_start()), c.unlinked(breakend2.uses_start()));
                if copy_numbers_equal(u1, u2) && pairs[i].position_distance() > pairs[j].position_distance() {
                    pairs.remove(i);
                    removed = true;
                } else if u1 < u2 {
                    pairs.remove(i);
                    removed = true;
                } else {
                    pairs.remove(j);
                }
                break;
            }
            if !removed {
                i += 1;
            }
        }
    }

    fn find_complex_dups(&mut self, proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        if self.complex_dup_candidates.is_empty() {
            return proposed;
        }
        if !proposed.is_empty() {
            return proposed;
        }
        let mut new_proposed: Vec<ProposedLinks> = Vec::new();
        let mut invalid: Vec<Sv> = Vec::new();
        let entries: Vec<(Sv, Vec<Pair>)> = self.complex_dup_candidates.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        for (comp_dup, links) in &entries {
            let conn = self.conn(comp_dup);
            let Some(conn) = conn.filter(|c| !c.borrow().has_connections()) else {
                invalid.push(comp_dup.clone());
                continue;
            };
            let cd_ploidy = conn.borrow().jcn;
            if cd_ploidy == 0.0 {
                continue;
            }
            let cd_start = comp_dup.be(true);
            let cd_end = comp_dup.be(false);
            if links.len() != 2 {
                invalid.push(comp_dup.clone());
                continue;
            }
            let pair1 = &links[0];
            let pair2 = &links[1];
            let other_sv1 = pair1.other_sv(comp_dup);
            let other_sv2 = pair2.other_sv(comp_dup);
            if other_sv1 == other_sv2 {
                let other_conn = self.conn(&other_sv1);
                let ob1 = other_sv1.be(true);
                let ob2 = other_sv1.be(false);
                let Some(other_conn) = other_conn.filter(|c| !c.borrow().has_connections()) else {
                    invalid.push(comp_dup.clone());
                    continue;
                };
                let mut pl = ProposedLinks::multi(vec![pair1.clone(), pair2.clone()], ChainingRule::COMP_DUP_SPLIT, None, None);
                let oj = other_conn.borrow().jcn;
                pl.add_comp_dup_breakends(&cd_start, &cd_end, cd_ploidy, &ob1, &ob2, oj);
                new_proposed.push(pl);
            } else {
                let ob1 = if pair1.has_breakend(&cd_start) { pair1.other_breakend(&cd_start) } else { pair1.other_breakend(&cd_end) };
                let ob2 = if pair2.has_breakend(&cd_start) { pair2.other_breakend(&cd_start) } else { pair2.other_breakend(&cd_end) };
                for chain in &self.chains {
                    let cbs = chain.open_breakend(true);
                    let cbe = chain.open_breakend(false);
                    if !jcn_match_for_splits(cd_ploidy, comp_dup.jcn_uncertainty(), chain.jcn(), chain.jcn_uncertainty()) {
                        continue;
                    }
                    let (o1, o2) = (Some(ob1.clone()), Some(ob2.clone()));
                    if (cbs == o1 && cbe == o2) || (cbs == o2 && cbe == o1) {
                        let mut pl =
                            ProposedLinks::multi(vec![pair1.clone(), pair2.clone()], ChainingRule::COMP_DUP_SPLIT, Some(chain.clone()), None);
                        pl.add_comp_dup_breakends(&cd_start, &cd_end, cd_ploidy, cbs.as_ref().unwrap(), cbe.as_ref().unwrap(), chain.jcn());
                        new_proposed.push(pl);
                        break;
                    }
                }
            }
        }
        for x in &invalid {
            self.complex_dup_candidates.remove(x);
        }
        new_proposed
    }

    fn conn_values(&self) -> Vec<ChainStateRef> {
        self.sv_connections.values().cloned().collect()
    }

    fn find_jcn_match_pairs(&mut self, mut proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        let mut new_proposed: Vec<ProposedLinks> = Vec::new();
        let has_support = any_links_have_jcn_support(&proposed);
        if !proposed.is_empty() {
            for p in proposed.iter_mut() {
                if p.jcn_match_type() == PM_MATCHED {
                    p.add_rule(ChainingRule::JCN_MATCH);
                }
            }
            for p in proposed.iter_mut() {
                if p.jcn_match_type() == PM_OVERLAP {
                    p.add_rule(ChainingRule::JCN_OVERLAP);
                }
            }
            if has_support {
                return proposed;
            }
        }
        let mut added: Vec<Pair> = Vec::new();
        for conn in self.conn_values() {
            let var = conn.borrow().sv.clone();
            for be in [SE_START, SE_END] {
                if var.is_sgl_breakend() && be == SE_END {
                    continue;
                }
                let s = is_start(be);
                let breakend_ploidy = conn.borrow().unlinked_se(be);
                if breakend_ploidy == 0.0 {
                    continue;
                }
                let breakend = var.be(s);
                let Some(sv_links) = self.sv_breakend_possible_links.get(&breakend).cloned() else { continue };
                let mut best: Option<ProposedLinks> = None;
                for pair in &sv_links {
                    if self.has_skipped_pairs(pair) {
                        continue;
                    }
                    if contains(&added, pair) {
                        continue;
                    }
                    let other = pair.other_breakend(&breakend);
                    let other_ploidy = self.unlinked_breakend_count_chains(&other, true);
                    if other_ploidy == 0.0 {
                        continue;
                    }
                    // Java passes the SV's uncertainty and the breakend's JCN to jcnOverlap in swapped order (LX4)
                    let overlap = if crate::compat::java("LX4") {
                        jcn_overlap(var.jcn_uncertainty(), breakend_ploidy, other_ploidy, other.jcn_uncertainty())
                    } else {
                        jcn_overlap(breakend_ploidy, var.jcn_uncertainty(), other_ploidy, other.jcn_uncertainty())
                    };
                    let match_type = if copy_numbers_equal(other_ploidy, breakend_ploidy) {
                        PM_MATCHED
                    } else if overlap {
                        PM_OVERLAP
                    } else {
                        continue;
                    };
                    if best.as_ref().is_some_and(|b| b.top_rule() == ChainingRule::JCN_OVERLAP) && match_type != PM_MATCHED {
                        continue;
                    }
                    let mut pl = ProposedLinks::single(
                        pair,
                        if match_type == PM_MATCHED { ChainingRule::JCN_MATCH } else { ChainingRule::JCN_OVERLAP },
                    );
                    pl.add_breakend_ploidies(&breakend, breakend_ploidy, &other, other_ploidy);
                    best = Some(pl);
                    added.push(pair.clone());
                    if match_type == PM_MATCHED {
                        break;
                    }
                }
                if let Some(b) = best {
                    new_proposed.push(b);
                }
            }
        }
        self.check_cluster_jcn_support(&mut new_proposed);
        if proposed.is_empty() || any_links_have_jcn_support(&new_proposed) { new_proposed } else { proposed }
    }

    fn find_adjacent_matching_pairs(&mut self, mut proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        if self.adjacent_matching_pairs.is_empty() {
            return proposed;
        }
        if !proposed.is_empty() {
            for p in proposed.iter_mut() {
                if p.links.iter().any(|x| contains(&self.adjacent_matching_pairs, x)) {
                    p.add_rule(ChainingRule::ADJACENT);
                }
            }
            return proposed;
        }
        let mut new_proposed: Vec<ProposedLinks> = Vec::new();
        while !self.adjacent_matching_pairs.is_empty() {
            let next = self.adjacent_matching_pairs.remove(0);
            if self.matches_existing_pair(&next) {
                continue;
            }
            let pf = self.unlinked_breakend_count_chains(&next.first_breakend(), true);
            let ps = self.unlinked_breakend_count_chains(&next.second_breakend(), true);
            if pf == 0.0 || ps == 0.0 {
                continue;
            }
            if !copy_numbers_equal(pf, ps) {
                continue;
            }
            let mut pl = ProposedLinks::single(&next, ChainingRule::ADJACENT);
            pl.add_rule(ChainingRule::JCN_MATCH);
            pl.add_breakend_ploidies(&next.first_breakend(), pf, &next.second_breakend(), ps);
            new_proposed.push(pl);
        }
        restrict_proposed_links(proposed, new_proposed, ChainingRule::ADJACENT)
    }

    fn find_adjacent_pairs(&mut self, mut proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        if self.adjacent_pairs.is_empty() {
            return proposed;
        }
        let has_support = any_links_have_jcn_support(&proposed);
        if !proposed.is_empty() {
            for p in proposed.iter_mut() {
                if p.links.iter().any(|y| contains(&self.adjacent_pairs, y)) {
                    p.add_rule(ChainingRule::ADJACENT);
                }
            }
            if has_support {
                return proposed;
            }
        }
        let mut new_proposed: Vec<ProposedLinks> = Vec::new();
        while !self.adjacent_pairs.is_empty() {
            let next = self.adjacent_pairs.remove(0);
            if self.matches_existing_pair(&next) {
                continue;
            }
            let pf = self.unlinked_breakend_count_chains(&next.first_breakend(), true);
            let ps = self.unlinked_breakend_count_chains(&next.second_breakend(), true);
            if pf == 0.0 || ps == 0.0 {
                continue;
            }
            let mut pl = ProposedLinks::single(&next, ChainingRule::ADJACENT);
            pl.add_breakend_ploidies(&next.first_breakend(), pf, &next.second_breakend(), ps);
            new_proposed.push(pl);
        }
        self.check_cluster_jcn_support(&mut new_proposed);
        if proposed.is_empty() || any_links_have_jcn_support(&new_proposed) { new_proposed } else { proposed }
    }

    fn find_highest_jcn(&mut self, mut proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        let mut new_proposed: Vec<ProposedLinks> = Vec::new();
        let has_support = any_links_have_jcn_support(&proposed);
        if !proposed.is_empty() {
            let max_ploidy = proposed.iter().map(|x| x.jcn()).fold(f64::NEG_INFINITY, jmax);
            for p in proposed.iter_mut() {
                if copy_numbers_equal(max_ploidy, p.jcn()) {
                    p.add_rule(ChainingRule::JCN_MAX);
                }
            }
            if has_support {
                return proposed;
            }
        }
        let mut current_max = 0.0;
        let added: Vec<Pair> = Vec::new();
        for conn in self.conn_values() {
            let var = conn.borrow().sv.clone();
            for be in [SE_START, SE_END] {
                if var.is_sgl_breakend() && be == SE_END {
                    continue;
                }
                let s = is_start(be);
                let mut breakend_ploidy = conn.borrow().unlinked_se(be);
                if breakend_ploidy == 0.0 {
                    continue;
                }
                if !copy_numbers_equal(breakend_ploidy, current_max) && breakend_ploidy < current_max {
                    continue;
                }
                let breakend = var.be(s);
                let sv_links = self.sv_breakend_possible_links.get(&breakend).cloned();
                breakend_ploidy = self.unlinked_breakend_count_chains(&breakend, true);
                let Some(sv_links) = sv_links else { continue };
                for pair in &sv_links {
                    if contains(&added, pair) {
                        continue;
                    }
                    if self.has_skipped_pairs(pair) {
                        continue;
                    }
                    let other = pair.other_breakend(&breakend);
                    let other_ploidy = self.unlinked_breakend_count_chains(&other, true);
                    if other_ploidy == 0.0 {
                        continue;
                    }
                    let min_pair = jmin(other_ploidy, breakend_ploidy);
                    if !copy_numbers_equal(current_max, min_pair) && min_pair < current_max {
                        continue;
                    }
                    current_max = jmax(min_pair, current_max);
                    let mut pl = ProposedLinks::single(pair, ChainingRule::JCN_MAX);
                    pl.add_breakend_ploidies(&breakend, breakend_ploidy, &other, other_ploidy);
                    new_proposed.push(pl);
                }
            }
        }
        self.check_cluster_jcn_support(&mut new_proposed);
        if proposed.is_empty() || any_links_have_jcn_support(&new_proposed) { new_proposed } else { proposed }
    }

    fn find_nearest(&mut self, mut proposed: Vec<ProposedLinks>) -> Vec<ProposedLinks> {
        if proposed.is_empty() {
            for conn in self.conn_values() {
                let var = conn.borrow().sv.clone();
                for be in [SE_START, SE_END] {
                    if var.is_sgl_breakend() && be == SE_END {
                        continue;
                    }
                    let s = is_start(be);
                    let breakend_ploidy = conn.borrow().unlinked_se(be);
                    if breakend_ploidy == 0.0 {
                        continue;
                    }
                    let breakend = var.be(s);
                    let breakend_ploidy = self.unlinked_breakend_count_chains(&breakend, true);
                    let Some(sv_links) = self.sv_breakend_possible_links.get(&breakend).cloned() else { continue };
                    for pair in &sv_links {
                        if self.has_skipped_pairs(pair) {
                            continue;
                        }
                        let other = pair.other_breakend(&breakend);
                        let other_ploidy = self.unlinked_breakend_count_chains(&other, true);
                        if other_ploidy == 0.0 {
                            continue;
                        }
                        let mut pl = ProposedLinks::single(pair, ChainingRule::NEAREST);
                        pl.add_breakend_ploidies(&breakend, breakend_ploidy, &other, other_ploidy);
                        proposed.push(pl);
                    }
                }
            }
        }
        let mut shortest: Vec<ProposedLinks> = Vec::new();
        for pl in proposed {
            let mut index = 0;
            let mut add_new = true;
            while index < shortest.len() {
                if links_have_clash(&pl.links, &shortest[index].links) {
                    if pl.shortest_link_distance() < shortest[index].shortest_link_distance() {
                        shortest.remove(index);
                    } else {
                        add_new = false;
                        break;
                    }
                } else {
                    index += 1;
                }
            }
            if !add_new {
                continue;
            }
            index = 0;
            while index < shortest.len() {
                if pl.shortest_link_distance() < shortest[index].shortest_link_distance() {
                    break;
                }
                index += 1;
            }
            shortest.insert(index, pl);
        }
        shortest
    }

    fn check_cluster_jcn_support(&self, proposed: &mut [ProposedLinks]) -> bool {
        let mut any = false;
        for pl in proposed.iter_mut() {
            if pl.has_rule(ChainingRule::CA_JCN_SUPPORT) {
                any = true;
                continue;
            }
            if pl.links.iter().any(|x| !self.jcn_limits.link_has_jcn_support(x, pl.jcn())) {
                continue;
            }
            pl.add_rule(ChainingRule::CA_JCN_SUPPORT);
            any = true;
        }
        any
    }
}

fn chain_has_foldback(chain: &SvChain, var: &Sv) -> bool {
    let lp = chain.linked_pairs();
    for i in 0..lp.len().saturating_sub(1) {
        let (pair, next) = (&lp[i], &lp[i + 1]);
        if &pair.second() == var && &next.first() == var && pair.first_breakend() == next.second_breakend() {
            return true;
        }
    }
    false
}

fn any_links_have_jcn_support(proposed: &[ProposedLinks]) -> bool {
    proposed.iter().any(|x| x.has_rule(ChainingRule::CA_JCN_SUPPORT))
}

fn restrict_proposed_links(mut proposed: Vec<ProposedLinks>, new_proposed: Vec<ProposedLinks>, new_rule: ChainingRule) -> Vec<ProposedLinks> {
    if proposed.is_empty() {
        return new_proposed;
    }
    if new_proposed.is_empty() {
        return proposed;
    }
    for pl in proposed.iter_mut() {
        if new_proposed.iter().any(|n| pl.links.iter().any(|p| contains(&n.links, p))) {
            pl.add_rule(new_rule);
        }
    }
    proposed
}

fn cull_by_priority(proposed: &mut Vec<ProposedLinks>) {
    if proposed.len() <= 1 {
        return;
    }
    let max = proposed.iter().map(|x| x.priority()).max().unwrap();
    proposed.retain(|x| x.priority() >= max);
}
