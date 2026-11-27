//! ChimericReadTracker, ChimericUtils, SupplementaryJunctionData and HardFilteredCache's per-chromosome part.

use crate::common::*;
use crate::compat;
use crate::config::{IsofoxConfig, MAX_NOVEL_SJ_DISTANCE};
use crate::fusion_frag::*;
use crate::gene::GeneCollection;
use crate::read::{clipped_side, parse_cigar, Read, OP_D, OP_H, OP_M, OP_N, OP_S};
use hmf_common::ensembl::{EnsemblDataCache, TranscriptData};
use hmf_common::java_seq_hash::SeqHashMap;
use hmf_common::known_fusion::KnownFusionType;
use std::collections::{HashMap, HashSet};

const CHIMERIC_SHORT_INV_MIN_LENGTH: i32 = 2000;

#[derive(Clone, Debug, Default)]
pub struct ChimericStats {
    pub chimeric_junctions: i32,
    pub local_inter_gene_frags: i32,
    pub candidate_realign_frags: i32,
    pub translocations: i32,
    pub inversions: i32,
    pub excluded: i32,
}

impl ChimericStats {
    #[allow(dead_code)]
    pub fn merge(&mut self, o: &ChimericStats) {
        self.chimeric_junctions += o.chimeric_junctions;
        self.local_inter_gene_frags += o.local_inter_gene_frags;
        self.candidate_realign_frags += o.candidate_realign_frags;
        self.translocations += o.translocations;
        self.inversions += o.inversions;
        self.excluded += o.excluded;
    }
}

/// ChimericReadGroup: reads by index into the gene collection's read store
#[derive(Clone, Debug)]
pub struct ChimericReadGroup {
    pub reads: Vec<usize>,
    pub complete: bool,
}

pub fn read_group_complete(reads: &[usize], store: &[Read]) -> bool {
    let mut supp = 0;
    let mut non_supp = 0;
    let mut exp_supp = 0;
    let mut exp_non_supp = 1;
    for i in reads {
        let r = &store[*i];
        if r.is_read_paired() && !r.is_mate_unmapped() {
            exp_non_supp = 2;
        }
        if r.is_supplementary() {
            supp += 1;
        } else {
            non_supp += 1;
            if r.has_supp_alignment() {
                exp_supp += 1;
            }
        }
    }
    exp_non_supp == non_supp && exp_supp == supp
}

impl ChimericReadGroup {
    pub fn new(reads: Vec<usize>, store: &[Read]) -> Self {
        let complete = read_group_complete(&reads, store);
        ChimericReadGroup { reads, complete }
    }
    fn add_read(&mut self, r: usize, store: &[Read]) {
        self.reads.push(r);
        self.complete = read_group_complete(&self.reads, store);
    }
}

/// Read.getJunctionMatchingTransRefs
fn junction_matching_trans_refs(read: &Read, gc: &GeneCollection, pos: i32, is_start: bool) -> Vec<TransExonRef> {
    let mut out = Vec::new();
    for (ri, mt) in read.mapped_regions.entries() {
        if !mt.exon_boundary() {
            continue;
        }
        let r = &gc.exon_regions[ri];
        if (is_start && r.end == pos) || (!is_start && r.start == pos) {
            out.extend(r.trans_exon_refs.iter().copied());
        }
    }
    out
}

fn has_candidate_junction_soft_clips_read(read: &Read, se: usize) -> bool {
    if !read.is_soft_clipped_no_region_match(se) {
        return false;
    }
    let sc = if se == SE_START { read.left_clip_length() } else { read.right_clip_length() };
    above_threshold(sc)
}

/// ChimericUtils.isRealignedFragmentCandidate
pub fn is_realigned_fragment_candidate(read: &Read) -> bool {
    has_candidate_junction_soft_clips_read(read, SE_START) || has_candidate_junction_soft_clips_read(read, SE_END)
}

/// ChimericUtils.setHasMultipleKnownSpliceGenes
pub fn set_has_multiple_known_splice_genes(reads: &[usize], store: &mut [Read], gc: &GeneCollection, known_pairs: &[(String, String)]) -> bool {
    let mut split_read: Option<usize> = None;
    let mut refs: [Vec<TransExonRef>; 2] = [Vec::new(), Vec::new()];
    for &ri in reads {
        let read = &store[ri];
        if read.has_split {
            if let Some(j) = crate::fusion_frag::read_split_junction(read) {
                refs[0].extend(junction_matching_trans_refs(read, gc, j[0], true));
                refs[1].extend(junction_matching_trans_refs(read, gc, j[1], false));
                split_read = Some(ri);
            }
            break;
        }
        if has_candidate_junction_soft_clips_read(read, SE_START) {
            refs[0].extend(junction_matching_trans_refs(read, gc, read.coords_boundary(SE_START), false));
        }
        if has_candidate_junction_soft_clips_read(read, SE_END) {
            refs[1].extend(junction_matching_trans_refs(read, gc, read.coords_boundary(SE_END), true));
        }
    }
    if refs[0].is_empty() || refs[1].is_empty() {
        return false;
    }
    let mut known_pair = false;
    let mut gene_match = false;
    for s in &refs[0] {
        for e in &refs[1] {
            if s.trans_id == e.trans_id {
                gene_match = true;
                if s.exon_rank == e.exon_rank - 1 {
                    return false;
                }
                continue;
            }
            if s.gene_id == e.gene_id {
                gene_match = true;
            } else if !known_pair && known_pairs.iter().any(|(u, d)| u == s.gene_id && d == e.gene_id) {
                known_pair = true;
            }
        }
    }
    if !gene_match || known_pair {
        if let Some(sr) = split_read {
            store[sr].has_inter_gene_split = true;
        }
        return true;
    }
    false
}

/// ChimericUtils.isInversion
fn is_inversion(reads: &[usize], store: &[Read]) -> bool {
    if reads.len() == 2 && !reads.iter().any(|r| store[*r].is_supplementary()) {
        let (r1, r2) = (&store[reads[0]], &store[reads[1]]);
        return r1.chromosome == r2.chromosome && r1.orient() == r2.orient() && (r1.pos_start - r2.pos_start).abs() >= CHIMERIC_SHORT_INV_MIN_LENGTH;
    }
    if !reads.iter().any(|r| store[*r].has_supp_alignment()) || reads.iter().filter(|r| store[**r].is_supplementary()).count() != 1 {
        return false;
    }
    let mut existing_orient = 0i8;
    let mut existing_chr = "";
    for &ri in reads {
        let r = &store[ri];
        if !r.has_supp_alignment() {
            continue;
        }
        if existing_chr.is_empty() {
            existing_chr = r.chromosome;
        } else if existing_chr != r.chromosome {
            return false;
        }
        let l = if r.is_soft_clipped_no_region_match(SE_START) { r.left_clip_length() } else { 0 };
        let rr = if r.is_soft_clipped_no_region_match(SE_END) { r.right_clip_length() } else { 0 };
        if l == 0 && rr == 0 {
            return false;
        }
        let o = if l >= rr { ORIENT_FWD } else { ORIENT_REV };
        if existing_orient == 0 {
            existing_orient = o;
        } else if existing_orient == o {
            return true;
        }
    }
    false
}

/// HardFilteredCache.formChromosomePairString: lower-ranked chromosome first
pub fn form_chromosome_pair_string(c1: &str, c2: &str) -> String {
    let (r1, r2) = (hmf_common::regions::chromosome_rank_java(c1), hmf_common::regions::chromosome_rank_java(c2));
    if r1 <= r2 { format!("{c1}_{c2}") } else { format!("{c2}_{c1}") }
}

#[derive(Clone, Debug)]
struct SuppJunction {
    read_ids: Vec<String>,
    local_pos: i32,
    local_orient: i8,
    remote_chr: String,
    remote_pos: i32,
    match_count: i32,
}

impl SuppJunction {
    fn matches(&self, o: &SuppJunction) -> bool {
        self.local_pos == o.local_pos && self.remote_chr == o.remote_chr && self.remote_pos == o.remote_pos && self.local_orient == o.local_orient
    }
    fn from_reads(reads: &[usize], store: &[Read]) -> Option<SuppJunction> {
        let r = reads.iter().map(|i| &store[*i]).find(|r| r.has_supp_alignment())?;
        let supp = r.supp_data.as_ref()?;
        let left = r.longest_soft_clip_is_left() == Some(true);
        let (local_pos, local_orient) = if left { (r.coords_boundary(SE_START), ORIENT_REV) } else { (r.coords_boundary(SE_END), ORIENT_FWD) };
        let cigar = parse_cigar(&supp.cigar).unwrap_or_default();
        // ClippedSide.fromCigar(cigar, true)
        let clipped = |e: Option<&crate::read::CigarElem>| e.map(|e| e.op == OP_S || e.op == OP_H).unwrap_or(false);
        let sc_left = if clipped(cigar.first()) { cigar[0].len as i32 } else { 0 };
        let sc_right = if clipped(cigar.last()) { cigar[cigar.len() - 1].len as i32 } else { 0 };
        if sc_left == 0 && sc_right == 0 {
            return None;
        }
        let is_left = if sc_left > 0 && sc_right > 0 { sc_left >= sc_right } else { sc_left > 0 };
        let remote_pos = if is_left {
            supp.position
        } else {
            let skipped: i32 = cigar.iter().filter(|e| e.op == OP_N || e.op == OP_M || e.op == OP_D).map(|e| e.len as i32).sum();
            supp.position + skipped - 1
        };
        Some(SuppJunction { read_ids: vec![r.id.to_string()], local_pos, local_orient, remote_chr: supp.chromosome.clone(), remote_pos, match_count: 0 })
    }
}

/// the known fusion pairs' and promiscuous genes' ids (registerKnownFusionPairs, the same for every chromosome)
#[derive(Clone, Debug, Default)]
pub struct KnownPairs {
    pub pairs: Vec<(String, String)>,
    pub genes: HashSet<String>,
}

impl KnownPairs {
    pub fn compute(config: &IsofoxConfig, cache: &EnsemblDataCache) -> KnownPairs {
        let mut k = KnownPairs::default();
        for t in [KnownFusionType::KNOWN_PAIR, KnownFusionType::PROMISCUOUS_5, KnownFusionType::PROMISCUOUS_3] {
            for d in config.fusions.known_fusions.data_by_type(t) {
                let up = if !d.five_gene.is_empty() { cache.gene_data_by_name(&d.five_gene) } else { None };
                let down = if !d.three_gene.is_empty() { cache.gene_data_by_name(&d.three_gene) } else { None };
                if t == KnownFusionType::KNOWN_PAIR {
                    if let (Some(u), Some(dn)) = (up, down) {
                        k.pairs.push((u.gene_id.clone(), dn.gene_id.clone()));
                    }
                }
                if let Some(u) = up {
                    k.genes.insert(u.gene_id.clone());
                }
                if let Some(dn) = down {
                    k.genes.insert(dn.gene_id.clone());
                }
            }
        }
        k
    }
}

/// ChimericReadTracker: the chimeric reads of a chromosome's gene collections
pub struct ChimericReadTracker {
    pub run_fusions: bool,
    pub enabled: bool,
    known_pair_gene_ids: Vec<(String, String)>,
    reference_known_gene_ids: HashSet<String>,
    known_gene_ids: HashSet<String>,
    pub chimeric_read_map: SeqHashMap<String, ChimericReadGroup>,
    pub fusion_read_group_map: SeqHashMap<String, FusionReadGroup>,
    pub junction_rac: JunctionRacFragments,
    pub local_chimeric_reads: Vec<Vec<usize>>,
    candidate_realigned_groups: Vec<Vec<usize>>,
    supplementary_junctions: HashMap<i32, Vec<SuppJunction>>,
    pub hard_filtered_read_ids: HashMap<String, HashSet<String>>,
    post_gene_read_map: HashMap<String, Vec<(Box<str>, u16)>>,
    previous_post_gene_read_map: HashMap<String, Vec<(Box<str>, u16)>>,
    pub stats: ChimericStats,
    gc_region_bounds: [i32; 2],
    gc_end_of_chromosome: bool,
    gc_chromosome: &'static str,
    min_hard_filter_frags: i32,
    skip_non_genic: bool,
    write_chimeric_only: bool,
    run_validations: bool,
}

impl ChimericReadTracker {
    pub fn new(config: &IsofoxConfig) -> Self {
        let run_fusions = config.run_function(crate::config::IsofoxFunction::Fusions);
        ChimericReadTracker {
            run_fusions,
            enabled: run_fusions || config.run_function(crate::config::IsofoxFunction::AltSpliceJunctions),
            known_pair_gene_ids: Vec::new(),
            reference_known_gene_ids: HashSet::new(),
            known_gene_ids: HashSet::new(),
            chimeric_read_map: SeqHashMap::new(),
            fusion_read_group_map: SeqHashMap::new(),
            junction_rac: JunctionRacFragments::default(),
            local_chimeric_reads: Vec::new(),
            candidate_realigned_groups: Vec::new(),
            supplementary_junctions: HashMap::new(),
            hard_filtered_read_ids: HashMap::new(),
            post_gene_read_map: HashMap::new(),
            previous_post_gene_read_map: HashMap::new(),
            stats: ChimericStats::default(),
            gc_region_bounds: [0, 0],
            gc_end_of_chromosome: false,
            gc_chromosome: "",
            min_hard_filter_frags: config.fusions.min_hard_filter_frags,
            skip_non_genic: config.fusions.skip_non_genic,
            write_chimeric_only: config.fusions.write_chimeric_only,
            run_validations: config.run_validations,
        }
    }

    /// registerKnownFusionPairs, from the pairs worked out once for the run
    pub fn register_known_fusion_pairs(&mut self, known: &KnownPairs) {
        self.known_pair_gene_ids = known.pairs.clone();
        self.reference_known_gene_ids = known.genes.clone();
    }

    pub fn initialise(&mut self, gc: &GeneCollection) {
        self.known_gene_ids = gc.gene_ids.iter().filter(|g| self.reference_known_gene_ids.contains(**g)).map(|g| g.to_string()).collect();
        self.junction_rac = JunctionRacFragments::default();
        self.previous_post_gene_read_map = std::mem::take(&mut self.post_gene_read_map);
        self.gc_region_bounds = gc.region_bounds;
        self.gc_end_of_chromosome = gc.end_of_chromosome;
        self.gc_chromosome = gc.chromosome;
    }

    pub fn clear_data(&mut self) {
        self.chimeric_read_map.clear();
        self.fusion_read_group_map.clear();
        self.candidate_realigned_groups.clear();
        self.stats = ChimericStats::default();
        self.local_chimeric_reads.clear();
        self.supplementary_junctions.clear();
    }

    /// extractJunctionRacFragments
    pub fn extract_junction_rac_fragments(&mut self) -> JunctionRacFragments {
        self.junction_rac.purge_empty_junctions();
        self.stats.chimeric_junctions += self.junction_rac.junction_count() as i32;
        self.stats.candidate_realign_frags += self.junction_rac.fragment_count() as i32;
        std::mem::take(&mut self.junction_rac)
    }

    pub fn is_chimeric(&self, reads: &[usize], store: &mut [Read], gc: &GeneCollection, multi_mapped: bool) -> bool {
        if reads.iter().any(|r| store[*r].is_chimeric() || !store[*r].within_gene_collection()) {
            return true;
        }
        if !multi_mapped && self.enabled && reads.iter().any(|r| store[*r].has_split) {
            return set_has_multiple_known_splice_genes(reads, store, gc, &self.known_pair_gene_ids);
        }
        false
    }

    fn in_immune_region(&self, read: &Read, immune: &[hmf_common::regions::ChrBaseRegion]) -> bool {
        let in_region = immune.iter().any(|x| x.chromosome == read.chromosome && positions_overlap(read.pos_start, read.pos_end, x.start, x.end));
        if in_region && immune.iter().any(|x| x.chromosome == read.mate_chromosome && position_within(read.mate_pos, x.start, x.end)) {
            return true;
        }
        if let Some(s) = &read.supp_data {
            if in_region && immune.iter().any(|x| x.chromosome == s.chromosome && position_within(s.position, x.start, x.end)) {
                return true;
            }
        }
        false
    }

    fn add_intronic_transcript_data(read: &mut Read, transcripts: &[&'static TranscriptData]) {
        if !(read.overlaps_gene_collection() && read.mapped_regions.is_empty()) {
            return;
        }
        let mut refs: Vec<TransExonRef> = Vec::new();
        for t in transcripts {
            if !read.mapped_coords.alignments_within(t.trans_start, t.trans_end) {
                continue;
            }
            for i in 0..t.exons.len().saturating_sub(1) {
                let (e, n) = (&t.exons[i], &t.exons[i + 1]);
                if read.mapped_coords.alignments_within(e.end, n.start) {
                    refs.push(TransExonRef::new(t, e.rank.min(n.rank)));
                    break;
                }
            }
        }
        if !refs.is_empty() {
            read.intron_trans_refs = Some(refs);
        }
    }

    pub fn add_realignment_candidates(&mut self, reads: &[usize], store: &[Read]) {
        if reads.iter().any(|r| store[*r].is_duplicate()) {
            return;
        }
        self.candidate_realigned_groups.push(reads.to_vec());
    }

    pub fn add_chimeric_fragment(&mut self, reads: &[usize], store: &mut [Read], gc: &GeneCollection, immune: &[hmf_common::regions::ChrBaseRegion]) {
        if reads.iter().any(|r| store[*r].is_duplicate() || self.in_immune_region(&store[*r], immune)) {
            return;
        }
        for r in reads {
            Self::add_intronic_transcript_data(&mut store[*r], &gc.transcripts);
        }
        let id = store[reads[0]].id.to_string();
        if self.run_validations && self.chimeric_read_map.contains_key(&id) {
            hmf_common::error!("overriding chimeric read({})", id);
            let g = self.chimeric_read_map.get_mut(&id).unwrap();
            for r in reads {
                g.add_read(*r, store);
            }
        } else {
            let g = ChimericReadGroup::new(reads.to_vec(), store);
            if !self.write_chimeric_only {
                self.chimeric_read_map.insert(id, g);
            }
        }
    }

    pub fn add_supplementary_read(&mut self, ri: usize, store: &mut [Read], gc: &GeneCollection, immune: &[hmf_common::regions::ChrBaseRegion]) {
        if self.in_immune_region(&store[ri], immune) {
            return;
        }
        let id = store[ri].id.to_string();
        match self.chimeric_read_map.get_mut(&id) {
            None => {
                let g = ChimericReadGroup::new(vec![ri], store);
                self.chimeric_read_map.insert(id, g);
            }
            Some(g) => g.add_read(ri, store),
        }
        Self::add_intronic_transcript_data(&mut store[ri], &gc.transcripts);
    }

    /// postProcessChimericReads: `lone_reads` are the reads left in the fragment tracker, in Java order
    #[allow(clippy::too_many_arguments)]
    pub fn post_process_chimeric_reads(
        &mut self,
        depth: &mut BaseDepth,
        lone_reads: &[usize],
        store: &mut [Read],
        gc: &mut GeneCollection,
        immune: &[hmf_common::regions::ChrBaseRegion],
        known_splice_sites: &HashMap<String, HashSet<i32>>,
    ) {
        for &ri in lone_reads {
            let r = &store[ri];
            if r.is_duplicate() || r.is_mate_unmapped() || self.in_immune_region(r, immune) {
                continue;
            }
            // IF12: these reads' depth was counted already (processIncompleteReads); Java counts it twice
            if compat::java("IF12") {
                let coords = store[ri].mapped_coords.alignments.clone();
                depth.process_read(&coords);
            }
            Self::add_intronic_transcript_data(&mut store[ri], &gc.transcripts);
            let id = store[ri].id.to_string();
            match self.chimeric_read_map.get_mut(&id) {
                None => {
                    let g = ChimericReadGroup::new(vec![ri], store);
                    self.chimeric_read_map.insert(id, g);
                }
                Some(g) => g.add_read(ri, store),
            }
        }

        let mut to_remove: Vec<String> = Vec::new();
        for id in self.chimeric_read_map.keys_vec() {
            let (mut reads, complete) = {
                let g = self.chimeric_read_map.get(&id).unwrap();
                (g.reads.clone(), g.complete)
            };
            let read_id = store[reads[0]].id.to_string();
            let count = reads.len();
            if self.run_fusions && self.skip_non_genic_reads(&mut reads, store) {
                self.chimeric_read_map.get_mut(&id).unwrap().reads = reads;
                to_remove.push(read_id);
                continue;
            }
            let removed = reads.len() < count;
            self.chimeric_read_map.get_mut(&id).unwrap().reads = reads.clone();
            if !removed && !self.keep_chimeric_group(&reads, complete, store, gc) {
                to_remove.push(read_id);
                continue;
            }
            if self.run_fusions {
                self.collect_candidate_junctions(&reads, store, gc);
                self.cache_supplementary_junction_candidate(&reads, store, gc);
            }
        }
        for id in to_remove {
            self.chimeric_read_map.remove(&id);
        }
        gc.add_count(FragmentType::Chimeric, self.chimeric_read_map.len() as i64);
        self.apply_hard_filter(known_splice_sites);

        if self.run_fusions {
            let candidates = std::mem::take(&mut self.candidate_realigned_groups);
            for reads in &candidates {
                let fusion_reads: Vec<FusionRead> = reads.iter().map(|r| FusionRead::new(&store[*r], &gc.exon_regions)).collect();
                let id = store[reads[0]].id.to_string();
                self.junction_rac.check_add_candidate_group(&id, &fusion_reads);
            }
            for (id, g) in self.chimeric_read_map.iter() {
                if self.skip_non_genic && !g.reads.iter().any(|r| store[*r].within_gene_collection()) {
                    continue;
                }
                if self.write_chimeric_only {
                    continue;
                }
                let mut fusion_reads: Vec<FusionRead> = g.reads.iter().map(|r| FusionRead::new(&store[*r], &gc.exon_regions)).collect();
                for fr in &mut fusion_reads {
                    fr.set_read_junction_depth(depth);
                }
                let gid = store[g.reads[0]].id.to_string();
                let _ = id;
                self.fusion_read_group_map.insert(gid.clone(), FusionReadGroup::new(&gid, fusion_reads));
            }
            if self.write_chimeric_only {
                self.chimeric_read_map.clear();
            }
            self.junction_rac.set_read_junction_depth(depth);
        } else {
            self.chimeric_read_map.clear();
            self.candidate_realigned_groups.clear();
            self.stats = ChimericStats::default();
        }
    }

    fn keep_chimeric_group(&mut self, reads: &[usize], complete: bool, store: &mut [Read], gc: &GeneCollection) -> bool {
        if reads.iter().any(|r| store[*r].is_translocation()) {
            self.stats.translocations += 1;
            return true;
        }
        if complete && is_inversion(reads, store) {
            self.stats.inversions += 1;
            return true;
        }
        let span_gene_splice_sites = if reads.iter().any(|r| store[*r].has_inter_gene_split) {
            true
        } else {
            set_has_multiple_known_splice_genes(reads, store, gc, &self.known_pair_gene_ids)
        };
        if reads.iter().any(|r| store[*r].spans_gene_collections()) {
            self.stats.local_inter_gene_frags += 1;
            return true;
        }
        if !complete {
            return true;
        }
        if reads.iter().any(|r| !store[*r].within_gene_collection()) {
            if reads.iter().filter(|r| store[**r].fully_non_genic()).count() == reads.len() {
                return false;
            }
            let min_pos = reads.iter().map(|r| store[*r].coords_boundary(SE_START)).min().unwrap_or(0);
            let max_pos = reads.iter().map(|r| store[*r].coords_boundary(SE_END)).max().unwrap_or(0);
            if self.gc_region_bounds[0] - min_pos > MAX_NOVEL_SJ_DISTANCE || max_pos - self.gc_region_bounds[1] > MAX_NOVEL_SJ_DISTANCE {
                return true;
            }
        }
        if span_gene_splice_sites {
            self.stats.local_inter_gene_frags += 1;
            return true;
        }
        self.local_chimeric_reads.push(reads.to_vec());
        false
    }

    fn skip_non_genic_reads(&mut self, reads: &mut Vec<usize>, store: &[Read]) -> bool {
        let post: Vec<usize> = if !self.gc_end_of_chromosome { reads.iter().copied().filter(|r| store[*r].pos_start > self.gc_region_bounds[1]).collect() } else { Vec::new() };
        if post.len() == reads.len() {
            return true;
        }
        let pre: Vec<usize> = reads.iter().copied().filter(|r| store[*r].pos_start < self.gc_region_bounds[0]).collect();
        if !pre.is_empty() {
            let id = store[pre[0]].id.to_string();
            if let Some(prev) = self.previous_post_gene_read_map.get(&id) {
                for p in &pre {
                    if prev.iter().any(|(pid, pflags)| store[*p].matches(pid, *pflags)) {
                        reads.retain(|x| x != p);
                    }
                }
                if reads.is_empty() {
                    return true;
                }
            }
        }
        if !post.is_empty() {
            let id = store[reads[0]].id.to_string();
            self.post_gene_read_map.insert(id, post.iter().map(|r| (store[*r].id.clone(), store[*r].flags)).collect());
        }
        false
    }

    fn add_junction(&mut self, ri: usize, store: &mut [Read], se: usize, pos: i32, orient: i8) {
        if pos > 0 {
            store[ri].set_junction_position(se, pos);
            self.junction_rac.add_junction(pos, orient);
        }
    }

    fn collect_candidate_junctions(&mut self, reads: &[usize], store: &mut [Read], gc: &GeneCollection) {
        let split = reads.iter().copied().find(|r| store[*r].has_split && (store[*r].spans_gene_collections() || store[*r].has_inter_gene_split));
        if let Some(sr) = split {
            if let Some(j) = read_split_junction(&store[sr]) {
                self.add_junction(sr, store, SE_START, j[0], ORIENT_FWD);
                self.add_junction(sr, store, SE_END, j[1], ORIENT_REV);
            }
            return;
        }
        if let Some(supp) = reads.iter().copied().find(|r| store[*r].has_supp_alignment()) {
            let (side, len) = clipped_side(&store[supp]);
            if len >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH {
                let pos = store[supp].coords_boundary(side);
                self.add_junction(supp, store, side, pos, if side == SE_START { ORIENT_REV } else { ORIENT_FWD });
            }
            return;
        }
        if reads.len() > 1 || self.known_gene_ids.is_empty() || !self.read_in_known_gene_gc(reads, store, gc) {
            return;
        }
        let r = reads[0];
        let (side, _) = clipped_side(&store[r]);
        if side == SE_START {
            let pos = store[r].coords_boundary(SE_START);
            self.add_junction(r, store, SE_START, pos, ORIENT_REV);
        } else {
            let pos = store[r].coords_boundary(SE_END);
            self.add_junction(r, store, SE_END, pos, ORIENT_FWD);
        }
    }

    fn read_in_known_gene_gc(&self, reads: &[usize], store: &[Read], gc: &GeneCollection) -> bool {
        if self.known_gene_ids.is_empty() {
            return false;
        }
        reads.iter().any(|r| {
            store[*r].mapped_regions.keys().iter().any(|ri| gc.exon_regions[*ri].trans_exon_refs.iter().any(|t| self.known_gene_ids.contains(t.gene_id)))
        })
    }

    fn apply_hard_filter(&mut self, known_splice_sites: &HashMap<String, HashSet<i32>>) {
        if self.min_hard_filter_frags <= 1 || !self.known_gene_ids.is_empty() {
            return;
        }
        let local = self.gc_chromosome;
        let matches = |chr: &str, pos: i32| known_splice_sites.get(chr).map(|s| s.contains(&pos)).unwrap_or(false);
        let juncs = std::mem::take(&mut self.supplementary_junctions);
        for list in juncs.values() {
            for j in list {
                if j.match_count + 1 >= self.min_hard_filter_frags {
                    continue;
                }
                if matches(local, j.local_pos) || matches(&j.remote_chr, j.remote_pos) {
                    continue;
                }
                let pair = form_chromosome_pair_string(local, &j.remote_chr);
                let set = self.hard_filtered_read_ids.entry(pair).or_default();
                for id in &j.read_ids {
                    set.insert(id.clone());
                    self.chimeric_read_map.remove(id);
                }
            }
        }
    }

    fn cache_supplementary_junction_candidate(&mut self, reads: &[usize], store: &[Read], gc: &GeneCollection) {
        if self.min_hard_filter_frags <= 1 || self.read_in_known_gene_gc(reads, store, gc) {
            return;
        }
        let Some(j) = SuppJunction::from_reads(reads, store) else { return };
        let list = self.supplementary_junctions.entry(j.local_pos).or_default();
        if let Some(m) = list.iter_mut().find(|x| x.matches(&j)) {
            m.match_count += 1;
        } else {
            list.push(j);
        }
    }
}
