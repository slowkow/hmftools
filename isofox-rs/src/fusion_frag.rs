//! FusionRead, FusionReadGroup, FusionTransExon, FusionFragment and FusionFragmentBuilder, and the
//! realignment-candidate fragments of junctions (JunctionRacFragments).

use crate::common::*;
use crate::gene::RegionReadData;
use crate::read::{Read, SuppData, FLAG_FIRST, FLAG_MATE_UNMAPPED, FLAG_PAIRED, FLAG_SUPPLEMENTARY, OP_N};
use hmf_common::ensembl::TranscriptData;
use hmf_common::java_hash::JavaHashMap;
use hmf_common::regions::BaseRegion;
use std::collections::HashMap;

pub const REALIGN_MIN_SOFT_CLIP_BASE_LENGTH: i32 = 3;
pub const SOFT_CLIP_JUNC_BUFFER: i32 = 3;
pub const JUNCTION_BASE_LENGTH: usize = 10;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FusionTransExon {
    pub trans_id: i32,
    pub exon_rank: i32,
}

pub fn merge_unique(existing: &mut Vec<FusionTransExon>, new: &[FusionTransExon]) {
    for t in new {
        if !existing.iter().any(|x| x == t) {
            existing.push(*t);
        }
    }
}

fn matches_within_range(te: &TransExonRef, f: &FusionTransExon, max_diff: i32) -> bool {
    if te.trans_id != f.trans_id {
        return false;
    }
    if max_diff >= 0 {
        f.exon_rank >= te.exon_rank && f.exon_rank <= te.exon_rank + max_diff
    } else {
        f.exon_rank >= te.exon_rank + max_diff && f.exon_rank <= te.exon_rank
    }
}

pub fn has_match_within_range(te_refs: &[TransExonRef], f_refs: &[FusionTransExon], max_diff: i32) -> bool {
    te_refs.iter().any(|x| f_refs.iter().any(|y| matches_within_range(x, y, max_diff)))
}

/// FusionRead: what the fusion routines keep of a read
#[derive(Clone, Debug)]
pub struct FusionRead {
    pub chromosome: &'static str,
    pub positions: [i32; 2],
    pub orientation: i8,
    pub mate_chromosome: &'static str,
    pub mate_pos_start: i32,
    pub flags: u16,
    pub map_quality: i16,
    pub soft_clip_lengths: [i32; 2],
    pub boundary_bases: [Vec<u8>; 2],
    pub read_base_length: i32,
    pub gene_collections: [i32; 2],
    pub is_genic: [bool; 2],
    pub has_inter_gene_split: bool,
    pub has_supp_alignment: bool,
    pub supp_data: Option<SuppData>,
    pub contains_split: bool,
    pub mapped_coords: Vec<BaseRegion>,
    pub junction_positions: [i32; 2],
    boundary_depth: Option<HashMap<i32, i32>>,
    junction_depth: Option<[i32; 2]>,
    trans_exon_refs: Vec<FusionTransExon>,
    region_match_type: RegionMatchType,
    upper_trans_exon_refs: Option<Vec<FusionTransExon>>,
    upper_region_match_type: RegionMatchType,
}

fn possibly_exonic(m: RegionMatchType) -> bool {
    matches!(m, RegionMatchType::ExonBoundary | RegionMatchType::ExonIntron | RegionMatchType::WithinExon)
}

impl FusionRead {
    pub fn new(read: &Read, regions: &[RegionReadData]) -> FusionRead {
        let sc = [
            if read.is_soft_clipped_no_region_match(SE_START) { read.left_clip_length() } else { 0 },
            if read.is_soft_clipped_no_region_match(SE_END) { read.right_clip_length() } else { 0 },
        ];
        let bases = read.read_bases();
        let len = bases.len() as i32;
        let start_bases = (sc[0] + 5).min(len) as usize;
        let end_bases = (sc[1] + 5).min(len) as usize;
        let mut f = FusionRead {
            chromosome: read.chromosome,
            positions: [read.pos_start, read.pos_end],
            orientation: read.orient(),
            mate_chromosome: read.mate_chromosome,
            mate_pos_start: read.mate_pos,
            flags: read.flags,
            map_quality: read.mapq as i16,
            soft_clip_lengths: sc,
            boundary_bases: [bases[..start_bases].to_vec(), bases[bases.len() - end_bases..].to_vec()],
            read_base_length: len,
            gene_collections: read.gene_collections,
            is_genic: read.is_genic,
            has_inter_gene_split: read.has_inter_gene_split,
            has_supp_alignment: read.has_supp_alignment(),
            supp_data: read.supp_data.clone(),
            contains_split: read.has_split,
            mapped_coords: read.mapped_coords.without_inferred(),
            junction_positions: [0, 0],
            boundary_depth: None,
            junction_depth: None,
            trans_exon_refs: Vec::new(),
            region_match_type: RegionMatchType::None,
            upper_trans_exon_refs: None,
            upper_region_match_type: RegionMatchType::None,
        };
        f.extract_read_trans_exon_refs(read, regions);
        if let Some(j) = read.junction_positions {
            f.junction_positions = j;
            f.junction_depth = Some([0, 0]);
        }
        f
    }

    fn extract_read_trans_exon_refs(&mut self, read: &Read, regions: &[RegionReadData]) {
        if !read.mapped_regions.is_empty() {
            for (ri, mt) in read.mapped_regions.entries() {
                if !possibly_exonic(mt) && possibly_exonic(self.region_match_type) {
                    continue;
                }
                if possibly_exonic(mt) && !possibly_exonic(self.region_match_type) {
                    self.region_match_type = mt;
                    self.trans_exon_refs.clear();
                }
                let refs: Vec<FusionTransExon> = regions[ri].trans_exon_refs.iter().map(|t| FusionTransExon { trans_id: t.trans_id, exon_rank: t.exon_rank }).collect();
                merge_unique(&mut self.trans_exon_refs, &refs);
            }
        } else if let Some(intron) = &read.intron_trans_refs {
            // extractTopTransExonRefs over the single INTRON entry
            let existing = self.region_match_type;
            if RegionMatchType::Intron.match_rank() >= existing.match_rank() {
                if RegionMatchType::Intron.match_rank() > existing.match_rank() {
                    self.trans_exon_refs.clear();
                }
                let refs: Vec<FusionTransExon> = intron.iter().map(|t| FusionTransExon { trans_id: t.trans_id, exon_rank: t.exon_rank }).collect();
                merge_unique(&mut self.trans_exon_refs, &refs);
                self.region_match_type = RegionMatchType::Intron;
            }
        }
    }

    pub fn coords_boundary(&self, se: usize) -> i32 {
        if se == SE_START { self.mapped_coords[0].start } else { self.mapped_coords[self.mapped_coords.len() - 1].end }
    }
    pub fn spans_gene_collections(&self) -> bool {
        self.gene_collections[0] != self.gene_collections[1]
    }
    pub fn is_soft_clipped(&self, se: usize) -> bool {
        self.soft_clip_lengths[se] > 0
    }
    pub fn is_longest_soft_clip(&self, se: usize) -> bool {
        self.soft_clip_lengths[se] > self.soft_clip_lengths[switch_index(se)]
    }
    pub fn is_supplementary(&self) -> bool {
        self.flags & FLAG_SUPPLEMENTARY != 0
    }
    pub fn is_first_of_pair(&self) -> bool {
        self.flags & FLAG_FIRST != 0
    }
    pub fn is_read_paired(&self) -> bool {
        self.flags & FLAG_PAIRED != 0
    }
    pub fn is_mate_unmapped(&self) -> bool {
        self.flags & FLAG_MATE_UNMAPPED != 0
    }
    pub fn trans_exon_refs(&self, se: usize) -> Option<&Vec<FusionTransExon>> {
        if self.spans_gene_collections() {
            if se == SE_START { Some(&self.trans_exon_refs) } else { self.upper_trans_exon_refs.as_ref() }
        } else {
            Some(&self.trans_exon_refs)
        }
    }
    pub fn region_match_type(&self, se: usize) -> RegionMatchType {
        if se == SE_END && self.upper_trans_exon_refs.is_some() { self.upper_region_match_type } else { self.region_match_type }
    }
    pub fn set_upper_trans_exon_refs(&mut self, refs: Vec<FusionTransExon>, mt: RegionMatchType) {
        if self.upper_trans_exon_refs.is_none() {
            self.upper_trans_exon_refs = Some(refs);
        }
        self.upper_region_match_type = mt;
    }
    pub fn set_read_junction_depth(&mut self, depth: &BaseDepth) {
        if self.junction_positions[0] > 0 || self.junction_positions[1] > 0 {
            let jd = self.junction_depth.get_or_insert([0, 0]);
            if self.junction_positions[0] > 0 {
                jd[0] = depth.depth_at_base(self.junction_positions[0]);
            }
            if self.junction_positions[1] > 0 {
                jd[1] = depth.depth_at_base(self.junction_positions[1]);
            }
            return;
        }
        let bd = self.boundary_depth.get_or_insert_with(HashMap::new);
        for c in &self.mapped_coords {
            bd.insert(c.start, depth.depth_at_base(c.start));
            bd.insert(c.end, depth.depth_at_base(c.end));
        }
    }
    pub fn junction_depth(&self, se: usize, pos: i32) -> i32 {
        if let Some(jd) = &self.junction_depth {
            return if self.junction_positions[se] == pos { jd[se] } else { 0 };
        }
        if let Some(bd) = &self.boundary_depth {
            return bd.get(&pos).copied().unwrap_or(0);
        }
        0
    }
    /// FusionUtils.checkMissingGeneData: the upper section's transcript-exons of a read spanning into the
    /// next gene collection
    pub fn check_missing_gene_data(&mut self, transcripts: &[&TranscriptData]) {
        if !self.is_genic[SE_END] || !self.spans_gene_collections() {
            return;
        }
        let upper = self.mapped_coords[self.mapped_coords.len() - 1];
        let mut refs: Vec<FusionTransExon> = Vec::new();
        let mut top = RegionMatchType::None;
        for t in transcripts {
            if !positions_within(upper.start, upper.end, t.trans_start, t.trans_end) {
                continue;
            }
            for e in &t.exons {
                if !positions_overlap(upper.start, upper.end, e.start, e.end) {
                    continue;
                }
                if e.start > upper.end {
                    break;
                }
                let mt = if upper.start == e.start || upper.end == e.end {
                    RegionMatchType::ExonBoundary
                } else if positions_within(upper.start, upper.end, e.start, e.end) {
                    RegionMatchType::WithinExon
                } else {
                    RegionMatchType::ExonIntron
                };
                if mt.match_rank() < top.match_rank() {
                    continue;
                }
                if mt.match_rank() > top.match_rank() {
                    refs.clear();
                    top = mt;
                }
                refs.push(FusionTransExon { trans_id: t.trans_id, exon_rank: e.rank });
            }
        }
        self.set_upper_trans_exon_refs(refs, top);
    }
}

/// FusionUtils.findSplitReadJunction
pub fn fusion_split_junction(read: &FusionRead) -> Option<[i32; 2]> {
    if !read.contains_split {
        return None;
    }
    Some(read.junction_positions)
}

/// ChimericUtils.findSplitReadJunction: the longest N-split of a read
pub fn read_split_junction(read: &Read) -> Option<[i32; 2]> {
    if !read.has_split {
        return None;
    }
    let mut pos = read.pos_start;
    let mut max: Option<u32> = None;
    let mut max_start = 0;
    let n = read.cigar.len();
    for e in &read.cigar[..n.saturating_sub(1)] {
        if e.op == OP_N && (max.is_none() || e.len > max.unwrap()) {
            max = Some(e.len);
            max_start = pos - 1;
        }
        if crate::read::consumes_ref(e.op) {
            pos += e.len as i32;
        }
    }
    max.map(|l| [max_start, max_start + l as i32 + 1])
}

pub fn above_threshold(sc: i32) -> bool {
    sc >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH
}

pub fn has_candidate_junction_soft_clips(read: &FusionRead) -> bool {
    above_threshold(read.soft_clip_lengths[0]) || above_threshold(read.soft_clip_lengths[1])
}

/// FusionReadGroup
#[derive(Clone, Debug)]
pub struct FusionReadGroup {
    pub read_id: Box<str>,
    pub reads: Vec<FusionRead>,
    pub complete: bool,
}

impl FusionReadGroup {
    pub fn new(read_id: &str, reads: Vec<FusionRead>) -> Self {
        let mut g = FusionReadGroup { read_id: read_id.into(), reads, complete: false };
        g.complete = g.group_complete();
        g
    }
    pub fn add_read(&mut self, r: FusionRead) {
        self.reads.push(r);
        self.complete = self.group_complete();
    }
    pub fn merge(&mut self, other: FusionReadGroup) {
        for r in other.reads {
            self.add_read(r);
        }
    }
    pub fn has_supp_alignment(&self) -> bool {
        self.reads.iter().any(|r| r.has_supp_alignment)
    }
    pub fn find_other_chromosome(&self, chromosome: &str) -> Option<String> {
        for r in &self.reads {
            if r.is_read_paired() && r.mate_chromosome != chromosome {
                return Some(r.mate_chromosome.to_string());
            }
            if let Some(s) = &r.supp_data {
                return Some(s.chromosome.clone());
            }
        }
        None
    }
    fn group_complete(&self) -> bool {
        let mut supp = 0;
        let mut non_supp = 0;
        let mut exp_supp = 0;
        let mut exp_non_supp = 1;
        for r in &self.reads {
            if r.is_read_paired() && !r.is_mate_unmapped() {
                exp_non_supp = 2;
            }
            if r.is_supplementary() {
                supp += 1;
            } else {
                non_supp += 1;
                if r.has_supp_alignment {
                    exp_supp += 1;
                }
            }
        }
        exp_non_supp == non_supp && exp_supp == supp
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragType {
    Unknown,
    MatchedJunction,
    RealignCandidate,
    Realigned,
    Discordant,
    DiscordantJunction,
}

impl FragType {
    pub fn is_junction_type(self) -> bool {
        self == FragType::MatchedJunction || self == FragType::DiscordantJunction
    }
    #[allow(dead_code)]
    pub fn name(self) -> &'static str {
        match self {
            FragType::Unknown => "UNKNOWN",
            FragType::MatchedJunction => "MATCHED_JUNCTION",
            FragType::RealignCandidate => "REALIGN_CANDIDATE",
            FragType::Realigned => "REALIGNED",
            FragType::Discordant => "DISCORDANT",
            FragType::DiscordantJunction => "DISCORDANT_JUNCTION",
        }
    }
}

pub fn form_location(chromosome: &str, gc: i32) -> String {
    format!("{chromosome}:{gc}")
}

/// FusionFragment
#[derive(Clone, Debug)]
pub struct FusionFragment {
    pub read_id: Box<str>,
    pub reads: Vec<FusionRead>,
    pub has_supp: bool,
    pub gene_collections: [i32; 2],
    pub chromosomes: [&'static str; 2],
    pub orientations: [i8; 2],
    pub junction_positions: [i32; 2],
    pub junction_orientations: [i8; 2],
    pub frag_type: FragType,
    pub location_ids: [String; 2],
    pub region_match_types: [RegionMatchType; 2],
    pub trans_exon_refs: [Vec<FusionTransExon>; 2],
}

impl FusionFragment {
    pub fn new(group: FusionReadGroup) -> FusionFragment {
        let has_supp = group.has_supp_alignment();
        let mut f = FusionFragment {
            read_id: group.read_id,
            reads: group.reads,
            has_supp,
            gene_collections: [NO_GENE_ID, NO_GENE_ID],
            chromosomes: ["", ""],
            orientations: [0, 0],
            junction_positions: [-1, -1],
            junction_orientations: [0, 0],
            frag_type: FragType::Unknown,
            location_ids: [String::new(), String::new()],
            region_match_types: [RegionMatchType::None, RegionMatchType::None],
            trans_exon_refs: [Vec::new(), Vec::new()],
        };
        set_fragment_properties(&mut f);
        if f.frag_type != FragType::Unknown {
            for se in [SE_START, SE_END] {
                f.location_ids[se] = form_location(f.chromosomes[se], f.gene_collections[se]);
            }
        }
        f.extract_transcript_exon_data();
        f
    }

    pub fn is_single_chromosome(&self) -> bool {
        self.chromosomes[0] == self.chromosomes[1]
    }
    pub fn is_single_gene_collection(&self) -> bool {
        self.chromosomes[0] == self.chromosomes[1] && self.gene_collections[0] == self.gene_collections[1]
    }
    pub fn is_unspliced(&self) -> bool {
        self.region_match_types[0] == RegionMatchType::Intron && self.region_match_types[1] == RegionMatchType::Intron
    }
    pub fn is_spliced(&self) -> bool {
        self.region_match_types[0].exon_boundary() && self.region_match_types[1].exon_boundary()
    }
    pub fn location_pair(&self) -> String {
        if self.is_single_gene_collection() { self.location_ids[0].clone() } else { format!("{}_{}", self.location_ids[0], self.location_ids[1]) }
    }
    pub fn position_hash(&self) -> String {
        format!("{}_{}_{}_{}", self.junction_positions[0], self.junction_orientations[0], self.junction_positions[1], self.junction_orientations[1])
    }
    #[allow(dead_code)]
    pub fn implied_sv_type(&self) -> &'static str {
        implied_sv_type(&self.chromosomes, self.junction_orientations)
    }
    /// readsByLocation: indices of the reads in the location of one junction
    pub fn reads_by_location(&self, se: usize) -> Vec<usize> {
        if self.frag_type == FragType::Unknown {
            return Vec::new();
        }
        if self.is_single_gene_collection() {
            return (0..self.reads.len()).collect();
        }
        (0..self.reads.len())
            .filter(|i| {
                let r = &self.reads[*i];
                r.chromosome == self.chromosomes[se] && (r.gene_collections[0] == self.gene_collections[se] || r.gene_collections[1] == self.gene_collections[se])
            })
            .collect()
    }

    fn extract_transcript_exon_data(&mut self) {
        for se in [SE_START, SE_END] {
            for ri in 0..self.reads.len() {
                let read = &self.reads[ri];
                if !self.is_single_gene_collection() {
                    if read.chromosome != self.chromosomes[se] {
                        continue;
                    }
                    if read.gene_collections[0] != self.gene_collections[se] && read.gene_collections[1] != self.gene_collections[se] {
                        continue;
                    }
                }
                if self.junction_positions[se] > 0 && !position_within(self.junction_positions[se], read.positions[0], read.positions[1]) {
                    continue;
                }
                let Some(refs) = read.trans_exon_refs(se) else { continue };
                let mt = read.region_match_type(se);
                if mt.match_rank() < self.region_match_types[se].match_rank() {
                    continue;
                }
                if mt.match_rank() > self.region_match_types[se].match_rank() {
                    self.region_match_types[se] = mt;
                    self.trans_exon_refs[se].clear();
                }
                let refs = refs.clone();
                merge_unique(&mut self.trans_exon_refs[se], &refs);
            }
        }
    }

    /// validateTranscriptExons
    pub fn validate_transcript_exons(&mut self, trans: &[&TranscriptData], se: usize) {
        if self.junction_positions[se] < 0 || !self.region_match_types[se].exon_boundary() {
            return;
        }
        let pos = self.junction_positions[se];
        let orient = self.junction_orientations[se];
        let mut i = 0;
        while i < self.trans_exon_refs[se].len() {
            let te = self.trans_exon_refs[se][i];
            let Some(t) = trans.iter().find(|t| t.trans_id == te.trans_id) else {
                self.trans_exon_refs[se].remove(i);
                continue;
            };
            let exon = t.exons.iter().find(|x| (orient == 1 && x.end == pos) || (orient == -1 && x.start == pos));
            match exon {
                Some(e) if e.rank == te.exon_rank => i += 1,
                _ => {
                    self.trans_exon_refs[se].remove(i);
                }
            }
        }
    }

    pub fn soft_clip_bases(&self, pos: i32, orient: i8) -> String {
        let se = if orient == ORIENT_FWD { SE_END } else { SE_START };
        for r in &self.reads {
            if r.soft_clip_lengths[se] == 0 {
                continue;
            }
            if r.positions[se] == pos {
                let sc = r.soft_clip_lengths[se] as usize;
                let b = &r.boundary_bases[se];
                let s = if orient == ORIENT_FWD { &b[b.len() - sc..] } else { &b[..sc] };
                return String::from_utf8_lossy(s).to_string();
            }
        }
        String::new()
    }
}

// ---- FusionFragmentBuilder ----

fn set_fragment_properties(f: &mut FusionFragment) {
    let mut chromosomes: Vec<&'static str> = Vec::with_capacity(2);
    for r in &f.reads {
        if !chromosomes.contains(&r.chromosome) {
            chromosomes.push(r.chromosome);
        }
    }
    if chromosomes.len() == 1 {
        f.chromosomes = [chromosomes[0], chromosomes[0]];
    } else if lower_chromosome(chromosomes[0], chromosomes[1]) {
        f.chromosomes = [chromosomes[0], chromosomes[1]];
    } else {
        f.chromosomes = [chromosomes[1], chromosomes[0]];
    }
    if f.has_supp {
        set_supp_align_junction_data(f);
        return;
    }
    if f.is_single_chromosome() {
        let split = f.reads.iter().position(|r| r.contains_split && (r.spans_gene_collections() || r.has_inter_gene_split));
        if let Some(si) = split {
            set_split_read_junction_data(f, si);
            return;
        }
        if f.reads.iter().any(has_candidate_junction_soft_clips) && !f.reads.iter().any(|r| r.spans_gene_collections()) {
            if f.reads.len() == 1 && !f.reads[0].is_read_paired() {
                set_single_soft_clip_junction_data(f);
                return;
            }
            if f.reads.len() == 2 && f.reads[0].gene_collections[0] == f.reads[1].gene_collections[1] {
                set_single_soft_clip_junction_data(f);
                return;
            }
        }
    }
    set_non_junction_data(f);
}

fn set_supp_align_junction_data(f: &mut FusionFragment) {
    let mut jp = [-1, -1];
    let mut jo: [i8; 2] = [0, 0];
    let mut gcs = [-1, -1];
    let mut chrs: [&'static str; 2] = ["", ""];
    let mut pi = 0;
    let mut longest_first = false;
    let mut longest = 0;
    for r in &f.reads {
        if !r.is_supplementary() && r.has_supp_alignment {
            let m = r.soft_clip_lengths[0].max(r.soft_clip_lengths[1]);
            if m > longest {
                longest = m;
                longest_first = r.is_first_of_pair();
            }
        }
    }
    for r in &f.reads {
        if !r.has_supp_alignment || r.is_first_of_pair() != longest_first {
            continue;
        }
        chrs[pi] = r.chromosome;
        if r.is_longest_soft_clip(SE_START) {
            jp[pi] = r.coords_boundary(SE_START);
            jo[pi] = ORIENT_REV;
            gcs[pi] = r.gene_collections[SE_START];
        } else {
            jp[pi] = r.coords_boundary(SE_END);
            jo[pi] = ORIENT_FWD;
            gcs[pi] = r.gene_collections[SE_END];
        }
        pi += 1;
        if pi >= 2 {
            break;
        }
    }
    if jp[0] < 0 || jp[1] < 0 {
        return;
    }
    let lower = if chrs[0] == chrs[1] {
        if jp[0] <= jp[1] { 0 } else { 1 }
    } else if lower_chromosome(chrs[0], chrs[1]) {
        0
    } else {
        1
    };
    for se in [SE_START, SE_END] {
        let i = if se == SE_START { lower } else { switch_index(lower) };
        f.junction_positions[se] = jp[i];
        f.junction_orientations[se] = jo[i];
        f.orientations[se] = jo[i];
        f.gene_collections[se] = gcs[i];
    }
    f.frag_type = FragType::MatchedJunction;
}

fn set_split_read_junction_data(f: &mut FusionFragment, si: usize) {
    if let Some(j) = fusion_split_junction(&f.reads[si]) {
        for se in [SE_START, SE_END] {
            f.gene_collections[se] = f.reads[si].gene_collections[se];
            f.junction_positions[se] = j[se];
        }
        f.junction_orientations = [ORIENT_FWD, ORIENT_REV];
        f.orientations = [ORIENT_FWD, ORIENT_REV];
        f.frag_type = FragType::MatchedJunction;
    }
}

fn set_single_soft_clip_junction_data(f: &mut FusionFragment) {
    let mut realign: Option<usize> = None;
    let mut max_sc = 0;
    let mut side = 0;
    for (i, r) in f.reads.iter().enumerate() {
        for se in [SE_START, SE_END] {
            if !above_threshold(r.soft_clip_lengths[se]) {
                continue;
            }
            if r.soft_clip_lengths[se] > max_sc {
                max_sc = r.soft_clip_lengths[se];
                side = se;
                realign = Some(i);
            }
        }
    }
    let Some(ri) = realign else { return };
    f.junction_positions[0] = f.reads[ri].coords_boundary(side);
    f.junction_orientations[0] = if side == SE_START { ORIENT_REV } else { ORIENT_FWD };
    f.orientations[0] = f.junction_orientations[0];
    let gc = f.reads[ri].gene_collections[side];
    f.gene_collections = [gc, gc];
    f.frag_type = FragType::RealignCandidate;
}

fn set_non_junction_data(f: &mut FusionFragment) {
    let mut chr_gcs: Vec<String> = Vec::with_capacity(2);
    let mut chromosomes: Vec<&'static str> = Vec::with_capacity(2);
    let mut positions: HashMap<String, i32> = HashMap::new();
    let mut gene_collections: Vec<i32> = Vec::with_capacity(2);
    let mut reads: HashMap<String, usize> = HashMap::new();
    let mut groups: JavaHashMap<String, Vec<usize>> = JavaHashMap::with_capacity(hmf_common::java_hash::guava_capacity(2));
    for (ri, r) in f.reads.iter().enumerate() {
        for se in [SE_START, SE_END] {
            if !r.spans_gene_collections() && se == SE_END {
                continue;
            }
            let id = form_location(r.chromosome, r.gene_collections[se]);
            match groups.get_mut(&id) {
                None => {
                    groups.insert(id.clone(), vec![ri]);
                    chr_gcs.push(id.clone());
                    chromosomes.push(r.chromosome);
                    gene_collections.push(r.gene_collections[se]);
                    positions.insert(id.clone(), r.coords_boundary(se));
                    reads.insert(id, ri);
                }
                Some(g) => {
                    if r.coords_boundary(se) < positions[&id] {
                        let p = positions[&id].min(r.coords_boundary(se));
                        positions.insert(id.clone(), p);
                        reads.insert(id.clone(), ri);
                    }
                    g.push(ri);
                }
            }
        }
    }
    if groups.len() == 1 {
        f.gene_collections = [gene_collections[0], gene_collections[0]];
        if f.reads.len() == 2 && f.reads[0].orientation == f.reads[1].orientation {
            f.frag_type = FragType::Discordant;
            f.orientations = [f.reads[0].orientation, f.reads[0].orientation];
        } else {
            f.orientations = [ORIENT_FWD, ORIENT_REV];
        }
        return;
    }
    let lower = if chromosomes[0] == chromosomes[1] {
        if positions[&chr_gcs[0]] <= positions[&chr_gcs[1]] { 0 } else { 1 }
    } else if lower_chromosome(chromosomes[0], chromosomes[1]) {
        0
    } else {
        1
    };
    for se in [SE_START, SE_END] {
        let i = if se == SE_START { lower } else { switch_index(lower) };
        f.gene_collections[se] = gene_collections[i];
        f.orientations[se] = f.reads[reads[&chr_gcs[i]]].orientation;
    }
    f.frag_type = FragType::Discordant;
    if groups.len() != 2 || f.reads.len() != 2 || !f.is_single_chromosome() {
        return;
    }
    let mut sc_pos = [-1, -1];
    let mut sc_orient: [i8; 2] = [0, 0];
    for se in [SE_START, SE_END] {
        let i = if se == SE_START { lower } else { switch_index(lower) };
        let r = &f.reads[reads[&chr_gcs[i]]];
        let required = switch_index(se);
        if !r.is_soft_clipped(required) {
            break;
        }
        if required == SE_START {
            sc_pos[se] = r.coords_boundary(SE_START);
            sc_orient[se] = ORIENT_REV;
        } else {
            sc_pos[se] = r.coords_boundary(SE_END);
            sc_orient[se] = ORIENT_FWD;
        }
    }
    if sc_pos[0] > 0 && sc_pos[0] < sc_pos[1] {
        f.frag_type = FragType::DiscordantJunction;
        for se in [SE_START, SE_END] {
            // IF9: Java indexes the soft-clip positions, already stored by junction, by read again, so
            // with the reads in upper-lower order the two junctions are swapped
            let i = if crate::compat::java("IF9") {
                if se == SE_START { lower } else { switch_index(lower) }
            } else {
                se
            };
            f.junction_positions[se] = sc_pos[i];
            f.junction_orientations[se] = sc_orient[i];
        }
    }
}

/// FusionReadData.softClippedReadSupportsJunction
pub fn soft_clipped_read_supports_junction(read: &FusionRead, junc_se: usize, pos: i32, orient: i8, junction_bases: Option<&[String; 2]>) -> bool {
    if orient == ORIENT_FWD {
        if !read.is_soft_clipped(SE_END) {
            return false;
        }
        let boundary = read.coords_boundary(SE_END);
        if !position_within(boundary, pos, pos + SOFT_CLIP_JUNC_BUFFER) {
            return false;
        }
        let sc = read.soft_clip_lengths[SE_END];
        if sc < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH {
            return false;
        }
        let Some(jb) = junction_bases else { return true };
        let adj = if boundary > pos { boundary - pos } else { 0 };
        let b = &read.boundary_bases[SE_END];
        let extra_len = ((sc + adj) as usize).min(b.len());
        let mut extra = &b[b.len() - extra_len..];
        if extra.len() > JUNCTION_BASE_LENGTH {
            extra = &extra[..JUNCTION_BASE_LENGTH];
        }
        jb[switch_index(junc_se)].as_bytes().starts_with(extra)
    } else {
        if !read.is_soft_clipped(SE_START) {
            return false;
        }
        let boundary = read.coords_boundary(SE_START);
        if !position_within(boundary, pos - SOFT_CLIP_JUNC_BUFFER, pos) {
            return false;
        }
        let sc = read.soft_clip_lengths[SE_START];
        if sc < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH {
            return false;
        }
        let Some(jb) = junction_bases else { return true };
        let adj = if boundary < pos { pos - boundary } else { 0 };
        let b = &read.boundary_bases[SE_START];
        let extra_len = ((sc + adj) as usize).min(b.len());
        let mut extra = &b[..extra_len];
        if extra.len() > JUNCTION_BASE_LENGTH {
            extra = &extra[extra.len() - JUNCTION_BASE_LENGTH..];
        }
        jb[switch_index(junc_se)].as_bytes().ends_with(extra)
    }
}

/// JunctionRacFragments: realignment-candidate fragments by the junction they may support
#[derive(Default)]
pub struct JunctionRacFragments {
    pub fragments: Vec<FusionFragment>,
    pos_junctions: JavaHashMap<i32, Vec<usize>>,
    neg_junctions: JavaHashMap<i32, Vec<usize>>,
    /// one fragment per read group (Java's mRacFragments)
    pub rac_fragments: Vec<usize>,
}

impl JunctionRacFragments {
    pub fn junction_count(&self) -> usize {
        self.pos_junctions.len() + self.neg_junctions.len()
    }
    pub fn fragment_count(&self) -> usize {
        self.rac_fragments.len()
    }
    pub fn add_junction(&mut self, pos: i32, orient: i8) {
        let m = if orient == ORIENT_FWD { &mut self.pos_junctions } else { &mut self.neg_junctions };
        if !m.contains_key(&pos) {
            m.insert(pos, Vec::new());
        }
    }
    pub fn junction_fragments(&self, orient: i8, pos: i32) -> Option<&Vec<usize>> {
        if orient == ORIENT_FWD { self.pos_junctions.get(&pos) } else { self.neg_junctions.get(&pos) }
    }
    /// checkAddCandidateGroup
    pub fn check_add_candidate_group(&mut self, read_id: &str, reads: &[FusionRead]) -> bool {
        let mut added = false;
        for i in 0..2 {
            let orient = if i == 0 { ORIENT_FWD } else { ORIENT_REV };
            let se = if i == 0 { SE_START } else { SE_END };
            let positions: Vec<i32> = if i == 0 { self.pos_junctions.keys().copied().collect() } else { self.neg_junctions.keys().copied().collect() };
            for pos in positions {
                if orient == ORIENT_FWD {
                    if reads.iter().any(|r| r.coords_boundary(SE_END) > pos + SOFT_CLIP_JUNC_BUFFER) {
                        continue;
                    }
                } else if reads.iter().any(|r| r.coords_boundary(SE_START) < pos - SOFT_CLIP_JUNC_BUFFER) {
                    continue;
                }
                if reads.iter().any(|r| soft_clipped_read_supports_junction(r, se, pos, orient, None)) {
                    let frag = FusionFragment::new(FusionReadGroup::new(read_id, reads.to_vec()));
                    if frag.frag_type != FragType::RealignCandidate {
                        continue;
                    }
                    self.fragments.push(frag);
                    let fi = self.fragments.len() - 1;
                    let m = if i == 0 { &mut self.pos_junctions } else { &mut self.neg_junctions };
                    m.get_mut(&pos).unwrap().push(fi);
                    if !added {
                        added = true;
                        self.rac_fragments.push(fi);
                    }
                }
            }
        }
        added
    }
    pub fn purge_empty_junctions(&mut self) {
        self.pos_junctions.retain(|_, v| !v.is_empty());
        self.neg_junctions.retain(|_, v| !v.is_empty());
    }
    pub fn set_read_junction_depth(&mut self, depth: &BaseDepth) {
        for fi in &self.rac_fragments {
            for r in &mut self.fragments[*fi].reads {
                r.set_read_junction_depth(depth);
            }
        }
    }
}
