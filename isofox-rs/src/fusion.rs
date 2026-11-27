//! FusionReadData, FusionData, PassingFusions, FusionUtils' junction adjustments and FusionFinder.

use crate::common::*;
use crate::compat;
use crate::config::{GeneRegionFilters, IsofoxConfig};
use crate::fusion_frag::*;
use crate::refgen::RefSource;
use hmf_common::ensembl::{EnsemblDataCache, GeneData, TranscriptData};
use hmf_common::java_hash::{JavaHashMap, JavaHashSet};
use hmf_common::java_seq_hash::SeqHashMap;
use hmf_common::known_fusion::{KnownFusionCache, KnownFusionType};
use hmf_common::regions::BaseRegion;
use std::collections::HashMap;

const FILTER_COHORT_LIMIT_KNOWN: i32 = 5;
const FILTER_COHORT_LIMIT_NOT_KNOWN: i32 = 2;
const FILTER_MIN_MAP_QUAL: i16 = 10;
const CLOSE_MATCH_BASES: i32 = 10;
const CLOSE_MATCH_FRAG_RATIO: f64 = 0.2;
const POSITION_REALIGN_DISTANCE: i32 = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JunctionType {
    Known,
    Canonical,
    Unknown,
}

impl JunctionType {
    pub fn name(self) -> &'static str {
        match self {
            JunctionType::Known => "KNOWN",
            JunctionType::Canonical => "CANONICAL",
            JunctionType::Unknown => "UNKNOWN",
        }
    }
    fn ordinal(self) -> i32 {
        match self {
            JunctionType::Known => 0,
            JunctionType::Canonical => 1,
            JunctionType::Unknown => 2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixedKnownType {
    KnownPair,
    KnownProm3,
    KnownOther,
    Prom5Known,
    Prom5Prom3,
    Prom5Other,
    OtherProm3,
    Other,
}

impl MixedKnownType {
    pub fn name(self) -> &'static str {
        match self {
            MixedKnownType::KnownPair => "KNOWN_PAIR",
            MixedKnownType::KnownProm3 => "KNOWN_PROM3",
            MixedKnownType::KnownOther => "KNOWN_OTHER",
            MixedKnownType::Prom5Known => "PROM5_KNOWN",
            MixedKnownType::Prom5Prom3 => "PROM5_PROM3",
            MixedKnownType::Prom5Other => "PROM5_OTHER",
            MixedKnownType::OtherProm3 => "OTHER_PROM3",
            MixedKnownType::Other => "OTHER",
        }
    }
    fn has_known_pair_gene(self) -> bool {
        matches!(self, MixedKnownType::KnownPair | MixedKnownType::KnownProm3 | MixedKnownType::KnownOther | MixedKnownType::Prom5Known)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FilterType {
    NotSet,
    Pass,
    Cohort,
    FragmentCount,
    AlleleFrequency,
    AnchorDistance,
}

impl FilterType {
    pub fn name(self) -> &'static str {
        match self {
            FilterType::NotSet => "NOT_SET",
            FilterType::Pass => "PASS",
            FilterType::Cohort => "COHORT",
            FilterType::FragmentCount => "FRAGMENT_COUNT",
            FilterType::AlleleFrequency => "ALLELE_FREQUENCY",
            FilterType::AnchorDistance => "ANCHOR_DISTANCE",
        }
    }
}

/// what the fusion routines need of the Ensembl data and the run's settings
pub struct FusionContext<'a> {
    pub cache: &'static EnsemblDataCache,
    pub trans_by_id: &'a HashMap<i32, &'static TranscriptData>,
    pub known: &'a KnownFusionCache,
    pub filters: &'a GeneRegionFilters,
    pub max_fragment_length: i32,
    pub cohort: &'a CohortFusions,
    pub min_hard_filter_frags: i32,
}

/// FusionReadData: a fusion candidate
pub struct FusionReadData {
    pub id: i32,
    pub location_id: String,
    /// the batch index of the fragment that established it
    pub initial: usize,
    frag_counts: Vec<(FragType, i32)>,
    incomplete: bool,
    related_proximate: Vec<i32>,
    related_spliced: Vec<i32>,
    pub chromosomes: [&'static str; 2],
    pub junction_positions: [i32; 2],
    pub junction_orientations: [i8; 2],
    soft_clip_bases: [String; 2],
    pub junction_ref_bases: [String; 2],
    pub post_junction_ref_bases: [String; 2],
    pub junction_splice_bases: [String; 2],
    split_junction_overlap: i32,
    split_junction_overlap_applied: bool,
    pub junction_types: [JunctionType; 2],
    pub trans_exon_refs: [Vec<TransExonRef>; 2],
    read_depth: [i32; 2],
    fusion_genes: [Option<&'static GeneData>; 2],
    stream_indices: [usize; 2],
    max_split_lengths: [i32; 2],
    junction_homology: [i32; 2],
}

impl FusionReadData {
    fn new(id: i32, frags: &[FusionFragment], fi: usize) -> FusionReadData {
        let f = &frags[fi];
        let mut depth = [0, 0];
        for se in [SE_START, SE_END] {
            for r in &f.reads {
                let d = r.junction_depth(se, f.junction_positions[se]);
                if d > 0 {
                    depth[se] = d;
                    break;
                }
            }
        }
        let mut fusion = FusionReadData {
            id,
            location_id: f.location_pair(),
            initial: fi,
            frag_counts: Vec::new(),
            incomplete: false,
            related_proximate: Vec::new(),
            related_spliced: Vec::new(),
            chromosomes: f.chromosomes,
            junction_positions: f.junction_positions,
            junction_orientations: f.junction_orientations,
            soft_clip_bases: [String::new(), String::new()],
            junction_ref_bases: [String::new(), String::new()],
            post_junction_ref_bases: [String::new(), String::new()],
            junction_splice_bases: [String::new(), String::new()],
            split_junction_overlap: 0,
            split_junction_overlap_applied: false,
            junction_types: [JunctionType::Unknown, JunctionType::Unknown],
            trans_exon_refs: [Vec::new(), Vec::new()],
            read_depth: depth,
            fusion_genes: [None, None],
            stream_indices: [SE_START, SE_END],
            max_split_lengths: [0, 0],
            junction_homology: [0, 0],
        };
        fusion.add_fusion_fragment(frags, fi);
        fusion
    }

    pub fn frag_type_count(&self, t: FragType) -> i32 {
        self.frag_counts.iter().find(|x| x.0 == t).map(|x| x.1).unwrap_or(0)
    }
    fn add_frag_type_count(&mut self, t: FragType, n: i32) {
        match self.frag_counts.iter_mut().find(|x| x.0 == t) {
            Some(e) => e.1 += n,
            None => self.frag_counts.push((t, n)),
        }
    }
    pub fn supporting_fragment_count(&self) -> i32 {
        self.frag_counts.iter().filter(|x| x.0.is_junction_type() || x.0 == FragType::Discordant).map(|x| x.1).sum()
    }

    /// addFusionFragment (fragments are not cached: CacheFragments is off with the pipeline's options)
    fn add_fusion_fragment(&mut self, frags: &[FusionFragment], fi: usize) {
        let f = &frags[fi];
        self.add_frag_type_count(f.frag_type, 1);
        if f.frag_type.is_junction_type() || fi == self.initial {
            for se in [SE_START, SE_END] {
                let reads: Vec<&FusionRead> = f.reads_by_location(se).into_iter().map(|i| &f.reads[i]).collect();
                set_max_split_mapped_length(se, &reads, self.junction_positions, self.junction_orientations, &mut self.max_split_lengths);
            }
        }
    }
    /// addFusionFragment for a realignment-candidate fragment outside the batch
    fn add_rac_fragment(&mut self, f: &FusionFragment) {
        self.add_frag_type_count(f.frag_type, 1);
        if f.frag_type.is_junction_type() {
            for se in [SE_START, SE_END] {
                let reads: Vec<&FusionRead> = f.reads_by_location(se).into_iter().map(|i| &f.reads[i]).collect();
                set_max_split_mapped_length(se, &reads, self.junction_positions, self.junction_orientations, &mut self.max_split_lengths);
            }
        }
    }

    fn is_known_spliced(&self, frags: &[FusionFragment]) -> bool {
        frags[self.initial].is_spliced()
    }
    fn is_unspliced(&self, frags: &[FusionFragment]) -> bool {
        let f = &frags[self.initial];
        f.is_unspliced() && f.frag_type.is_junction_type()
    }
    fn has_viable_genes(&self) -> bool {
        self.fusion_genes[0].is_some() && self.fusion_genes[1].is_some()
    }

    fn set_junction_bases_all(&mut self, frags: &[FusionFragment], refg: &mut RefSource) {
        for se in [SE_START, SE_END] {
            self.set_junction_bases(frags, refg, se);
        }
    }

    fn set_junction_bases(&mut self, frags: &[FusionFragment], refg: &mut RefSource, se: usize) {
        let f = &frags[self.initial];
        if !f.frag_type.is_junction_type() {
            return;
        }
        if crate::version::v201() {
            // ISOFOX 2.0.1: 10 bases each side of the junction, and no soft-clip bases
            let jbl = JUNCTION_BASE_LENGTH as i32;
            let pos = self.junction_positions[se];
            if self.junction_orientations[se] == ORIENT_FWD {
                let b = refg.base_string(self.chromosomes[se], pos - jbl + 1, pos + jbl);
                if b.len() != 2 * JUNCTION_BASE_LENGTH {
                    return;
                }
                self.junction_ref_bases[se] = b[..JUNCTION_BASE_LENGTH].to_string();
                self.post_junction_ref_bases[se] = b[JUNCTION_BASE_LENGTH..].to_string();
                self.junction_splice_bases[se] = b[JUNCTION_BASE_LENGTH..JUNCTION_BASE_LENGTH + 2].to_string();
            } else {
                let b = refg.base_string(self.chromosomes[se], pos - jbl, pos + jbl - 1);
                if b.len() != 2 * JUNCTION_BASE_LENGTH {
                    return;
                }
                self.junction_ref_bases[se] = b[JUNCTION_BASE_LENGTH..].to_string();
                self.junction_splice_bases[se] = b[JUNCTION_BASE_LENGTH - 2..JUNCTION_BASE_LENGTH].to_string();
                self.post_junction_ref_bases[se] = b[..JUNCTION_BASE_LENGTH].to_string();
            }
            return;
        }
        let rbl = (JUNCTION_BASE_LENGTH * 2) as i32;
        let pos = self.junction_positions[se];
        let sc = if f.frag_type == FragType::MatchedJunction { f.soft_clip_bases(pos, self.junction_orientations[se]) } else { String::new() };
        let jbl = JUNCTION_BASE_LENGTH;
        if self.junction_orientations[se] == ORIENT_FWD {
            let b = refg.base_string(self.chromosomes[se], pos - rbl + 1, pos + jbl as i32);
            // an exception (a junction near a contig end) leaves the bases as they were
            if b.len() != 3 * jbl {
                return;
            }
            self.junction_ref_bases[se] = b[..rbl as usize].to_string();
            self.post_junction_ref_bases[se] = b[rbl as usize..].to_string();
            self.junction_splice_bases[se] = self.post_junction_ref_bases[se][..2].to_string();
            self.soft_clip_bases[se] = if sc.len() > jbl { sc[..jbl].to_string() } else { sc };
        } else {
            let b = refg.base_string(self.chromosomes[se], pos - jbl as i32, pos + rbl - 1);
            if b.len() != 3 * jbl {
                return;
            }
            self.junction_ref_bases[se] = b[jbl..].to_string();
            self.post_junction_ref_bases[se] = b[..jbl].to_string();
            self.junction_splice_bases[se] = self.post_junction_ref_bases[se][jbl - 2..].to_string();
            let n = sc.len();
            self.soft_clip_bases[se] = if n > jbl { sc[n - jbl..].to_string() } else { sc };
        }
    }

    fn check_homology_position_adjustment(&mut self, frags: &[FusionFragment]) {
        if !frags[self.initial].frag_type.is_junction_type() {
            return;
        }
        let mut other_orient = self.junction_orientations[SE_END];
        let same = self.junction_orientations[SE_START] == other_orient;
        let sc = &self.soft_clip_bases[SE_START];
        if sc.len() < JUNCTION_BASE_LENGTH {
            return;
        }
        let mut other = self.junction_ref_bases[SE_END].clone();
        if same {
            other = String::from_utf8(hmf_common::codon::reverse_complement_bases(other.as_bytes())).unwrap_or_default();
            other_orient = -other_orient;
        }
        let overlap = if other_orient == ORIENT_FWD {
            match other.rfind(sc.as_str()) {
                None => return,
                Some(i) => other.len() as i32 - (i + sc.len()) as i32,
            }
        } else {
            match other.find(sc.as_str()) {
                None => return,
                Some(i) => i as i32,
            }
        };
        self.split_junction_overlap = overlap.max(0);
    }

    fn set_homology_offsets(&mut self) {
        let start_pos_orient = self.junction_orientations[SE_START] == ORIENT_FWD;
        for se in [SE_START, SE_END] {
            let ps = switch_index(se);
            for i in 1..JUNCTION_BASE_LENGTH {
                if self.junction_ref_bases[se].len() < i || self.post_junction_ref_bases[ps].len() < i {
                    break;
                }
                let jr = &self.junction_ref_bases[se];
                let junction_str = if self.junction_orientations[se] == ORIENT_FWD { &jr[jr.len() - i..] } else { &jr[..i] };
                let pj = &self.post_junction_ref_bases[ps];
                let adjacent = if self.junction_orientations[ps] == ORIENT_FWD { &pj[..i] } else { &pj[pj.len() - i..] };
                // IF10 (kept as Java): for an inversion Java reverses the adjacent bases but discards the result. A
                // reversal would not be right either (the other side reads the opposite strand, and the two
                // breakends shift in opposite directions), so the comparison stays as Java makes it.
                if junction_str != adjacent {
                    break;
                }
                let retreating = (se == SE_START && start_pos_orient) || (se == SE_END && !start_pos_orient);
                if retreating {
                    self.junction_homology[SE_START] = -(i as i32);
                } else {
                    self.junction_homology[SE_END] = i as i32;
                }
            }
        }
    }

    fn match_within_homology(&self, o: &FusionReadData) -> bool {
        if self.junction_positions[0] - o.junction_positions[0] != self.junction_positions[1] - o.junction_positions[1] {
            return false;
        }
        for se in [SE_START, SE_END] {
            if self.junction_orientations[se] != o.junction_orientations[se] {
                return false;
            }
            if self.junction_positions[se] == o.junction_positions[se] {
                continue;
            }
            if self.junction_positions[se] < o.junction_positions[se] {
                if self.junction_positions[se] + self.junction_homology[SE_END] < o.junction_positions[se] {
                    return false;
                }
            } else if self.junction_positions[se] + self.junction_homology[SE_START] > o.junction_positions[se] {
                return false;
            }
        }
        true
    }

    fn is_close_match(&self, o: &FusionReadData) -> bool {
        for se in [SE_START, SE_END] {
            if (self.junction_positions[se] - o.junction_positions[se]).abs() > CLOSE_MATCH_BASES {
                return false;
            }
        }
        let ratio = self.frag_type_count(FragType::MatchedJunction) as f64 / o.frag_type_count(FragType::MatchedJunction) as f64;
        ratio >= 1.0 / CLOSE_MATCH_FRAG_RATIO || ratio <= CLOSE_MATCH_FRAG_RATIO
    }

    fn set_stream_data(&mut self, up: &[&'static GeneData], down: &[&'static GeneData], start_is_upstream: bool) {
        if !up.is_empty() {
            self.stream_indices[FS_UP] = if start_is_upstream { SE_START } else { SE_END };
            if down.is_empty() {
                self.stream_indices[FS_DOWN] = switch_index(self.stream_indices[FS_UP]);
            }
            self.fusion_genes[FS_UP] = if is_null_gene(up[0]) { None } else { Some(up[0]) };
        }
        if !down.is_empty() {
            self.stream_indices[FS_DOWN] = if start_is_upstream { SE_END } else { SE_START };
            if up.is_empty() {
                self.stream_indices[FS_UP] = switch_index(self.stream_indices[FS_DOWN]);
            }
            self.fusion_genes[FS_DOWN] = if is_null_gene(down[0]) { None } else { Some(down[0]) };
        }
    }

    fn gene_strand_by_position(&self, se: usize) -> Option<i8> {
        self.fusion_genes[se].map(|g| g.strand)
    }
    pub fn gene_name(&self, stream: usize) -> &str {
        self.fusion_genes[stream].map(|g| g.gene_name.as_str()).unwrap_or("")
    }

    fn trans_exon_refs_by_stream(&self, fs: usize) -> &Vec<TransExonRef> {
        if self.has_viable_genes() { &self.trans_exon_refs[self.stream_indices[fs]] } else { &self.trans_exon_refs[fs] }
    }

    fn add_related(&mut self, id: i32, spliced: bool) {
        let l = if spliced { &mut self.related_spliced } else { &mut self.related_proximate };
        if !l.contains(&id) {
            l.push(id);
        }
    }

    fn can_add_discordant_fragment(&self, frags: &[FusionFragment], f: &FusionFragment, max_dist: i32) -> bool {
        let mut implied = f.reads[0].read_base_length * 2;
        let mut del_dup = true;
        // ISOFOX 2.0.1 checks every discordant fragment as a DEL or DUP one
        if f.reads.len() >= 2 && !crate::version::v201() {
            let (r1, r2) = (&f.reads[0], &f.reads[1]);
            if r1.chromosome != r2.chromosome || r1.orientation == r2.orientation {
                del_dup = false;
            }
        }
        let initial = &frags[self.initial];
        for se in [SE_START, SE_END] {
            if f.orientations[se] != self.junction_orientations[se] {
                return false;
            }
            let is_up = self.stream_indices[FS_UP] == se;
            let diff = if is_up {
                -2
            } else if f.region_match_types[se] == RegionMatchType::Intron {
                1
            } else {
                2
            };
            if !has_match_within_range(&self.trans_exon_refs[se], &f.trans_exon_refs[se], diff) {
                return false;
            }
            let reads = f.reads_by_location(se);
            if reads.is_empty() {
                return false;
            }
            let read = &f.reads[reads[0]];
            let fpos = read.coords_boundary(switch_index(se));
            if (self.junction_orientations[se] == ORIENT_FWD) == (fpos > self.junction_positions[se])
                && (del_dup || read.is_soft_clipped(se))
                && !soft_clipped_read_supports_junction(read, se, self.junction_positions[se], self.junction_orientations[se], Some(&self.junction_ref_bases))
            {
                return false;
            }
            if f.region_match_types[se] == RegionMatchType::Intron || initial.region_match_types[se] == RegionMatchType::Intron {
                if fpos < self.junction_positions[se] {
                    implied += self.junction_positions[se] - fpos;
                } else {
                    implied += fpos - self.junction_positions[se];
                }
            }
        }
        !(del_dup && implied > max_dist)
    }

    fn can_realign_fragment_to_junction(&self, frags: &[FusionFragment], f: &FusionFragment) -> bool {
        let initial = &frags[self.initial];
        let mut supporting = false;
        for se in [SE_START, SE_END] {
            for r in &f.reads {
                if r.chromosome != self.chromosomes[se] || initial.gene_collections[se] != r.gene_collections[se] {
                    continue;
                }
                if self.junction_orientations[se] == 1 && r.coords_boundary(SE_END) > self.junction_positions[se] + SOFT_CLIP_JUNC_BUFFER {
                    return false;
                } else if self.junction_orientations[se] == -1 && r.coords_boundary(SE_START) < self.junction_positions[se] - SOFT_CLIP_JUNC_BUFFER {
                    return false;
                }
                if !supporting && soft_clipped_read_supports_junction(r, se, self.junction_positions[se], self.junction_orientations[se], Some(&self.junction_ref_bases)) {
                    supporting = true;
                }
            }
        }
        supporting
    }

    fn to_fusion_data(&self, frags: &[FusionFragment]) -> FusionData {
        let valid = self.has_viable_genes() && !self.incomplete;
        let initial = &frags[self.initial];
        let mut split = 0;
        let mut realigned = 0;
        let mut discordant = 0;
        // mFragmentCounts: a HashMap keyed by enum (identity hashes), so in insertion order under hashCode=2
        let mut disc_sum = 0;
        for (t, n) in &self.frag_counts {
            match t {
                FragType::MatchedJunction => split = *n,
                FragType::Discordant | FragType::DiscordantJunction => {
                    discordant = *n;
                    disc_sum += *n;
                }
                FragType::Realigned => realigned = *n,
                _ => {}
            }
        }
        // IF11: Java keeps only the last of the DISCORDANT and DISCORDANT_JUNCTION counts
        if !compat::java("IF11") {
            discordant = disc_sum;
        }
        let total = split + realigned + discordant;
        let mut d = FusionData {
            id: self.id,
            valid_genes: valid,
            chromosomes: [""; 2],
            junction_positions: [0; 2],
            junction_orientations: [0; 2],
            junction_types: [JunctionType::Unknown; 2],
            sv_type: implied_sv_type(&self.chromosomes, self.junction_orientations),
            read_type: if initial.has_supp {
                "SuppAlign"
            } else if initial.reads.iter().any(|r| r.contains_split) {
                "Split"
            } else if initial.frag_type == FragType::DiscordantJunction {
                "Discordant"
            } else {
                "Other"
            },
            gene_ids: [String::new(), String::new()],
            gene_names: [String::new(), String::new()],
            strands: [0; 2],
            coverage: [0; 2],
            anchor_distance: [0; 2],
            transcripts: [String::new(), String::new()],
            exons: [0; 2],
            total_frags: total,
            split_frags: split,
            realigned_frags: realigned,
            discordant_frags: discordant,
            related_spliced: self.related_spliced.clone(),
            init_read_id: initial.read_id.to_string(),
            all_trans_data: [String::new(), String::new()],
            mixed_known_type: MixedKnownType::Other,
            known_type: KnownFusionType::NONE,
            cohort_frequency: 0,
            filter: FilterType::NotSet,
        };
        for fs in [FS_UP, FS_DOWN] {
            let si = self.stream_indices[fs];
            d.chromosomes[fs] = self.chromosomes[si];
            d.junction_positions[fs] = self.junction_positions[si];
            d.junction_orientations[fs] = self.junction_orientations[si];
            // ISOFOX 2.0.1 keeps the junction types on the initial fragment
            d.junction_types[fs] = if crate::version::v201() { initial.junction_types[si] } else { self.junction_types[si] };
            if let Some(g) = self.fusion_genes[fs] {
                d.gene_ids[fs] = g.gene_id.clone();
                d.gene_names[fs] = g.gene_name.clone();
                d.strands[fs] = g.strand;
            }
            d.coverage[fs] = self.read_depth[si].max(split);
            d.anchor_distance[fs] = self.max_split_lengths[si];
            let refs = self.trans_exon_refs_by_stream(fs);
            if refs.is_empty() {
                d.all_trans_data[fs] = "NONE".to_string();
            } else {
                let mut td: Vec<String> = Vec::new();
                for r in refs {
                    td.push(format!("{}-{}", r.trans_name, r.exon_rank));
                    if r.is_canonical() {
                        d.transcripts[fs] = r.trans_name.to_string();
                        d.exons[fs] = r.exon_rank;
                    }
                }
                if d.transcripts[fs].is_empty() {
                    d.transcripts[fs] = refs[0].trans_name.to_string();
                    d.exons[fs] = refs[0].exon_rank;
                }
                d.all_trans_data[fs] = td.join(";");
            }
        }
        d
    }
}

/// FusionUtils.setMaxSplitMappedLength
fn set_max_split_mapped_length(se: usize, reads: &[&FusionRead], jp: [i32; 2], jo: [i8; 2], max: &mut [i32; 2]) {
    let matching: Vec<&&FusionRead> = reads.iter().filter(|r| position_within(jp[se], r.positions[0], r.positions[1])).collect();
    if matching.is_empty() {
        return;
    }
    let coords: Vec<BaseRegion> = if matching.len() == 1 { matching[0].mapped_coords.clone() } else { derive_common_regions(&matching[0].mapped_coords, &matching[1].mapped_coords) };
    let mut mapped = 0;
    for c in &coords {
        if jo[se] == ORIENT_REV {
            if c.end < jp[se] {
                continue;
            }
            mapped += c.end - jp[se].max(c.start) + 1;
        } else {
            if c.start > jp[se] {
                break;
            }
            mapped += jp[se].min(c.end) - c.start + 1;
        }
    }
    max[se] = mapped.max(max[se]);
}

/// FusionData: a fusion as written
#[derive(Clone, Debug)]
pub struct FusionData {
    pub id: i32,
    pub valid_genes: bool,
    pub chromosomes: [&'static str; 2],
    pub junction_positions: [i32; 2],
    pub junction_orientations: [i8; 2],
    pub junction_types: [JunctionType; 2],
    pub sv_type: &'static str,
    pub read_type: &'static str,
    pub gene_ids: [String; 2],
    pub gene_names: [String; 2],
    pub strands: [i8; 2],
    pub coverage: [i32; 2],
    pub anchor_distance: [i32; 2],
    pub transcripts: [String; 2],
    pub exons: [i32; 2],
    pub total_frags: i32,
    pub split_frags: i32,
    pub realigned_frags: i32,
    pub discordant_frags: i32,
    /// the related spliced fusions' ids in insertion order (a HashSet<Integer> in Java)
    pub related_spliced: Vec<i32>,
    pub init_read_id: String,
    pub all_trans_data: [String; 2],
    pub mixed_known_type: MixedKnownType,
    pub known_type: KnownFusionType,
    pub cohort_frequency: i32,
    pub filter: FilterType,
}

impl FusionData {
    pub fn name(&self) -> String {
        format!("{}_{}", self.gene_names[0], self.gene_names[1])
    }
    fn is_related(&self, o: &FusionData) -> bool {
        if !self.related_spliced.contains(&o.id) {
            return false;
        }
        for se in [SE_START, SE_END] {
            if self.gene_names[se].is_empty() || o.gene_names[se] != self.gene_names[se] {
                return false;
            }
        }
        true
    }
    fn has_known_splice_sites(&self) -> bool {
        self.junction_types[0] == JunctionType::Known && self.junction_types[1] == JunctionType::Known
    }
    fn allele_frequency(&self) -> f64 {
        let c = self.coverage[0].max(self.coverage[1]) as f64;
        if c > 0.0 { (self.split_frags + self.realigned_frags) as f64 / c } else { 0.0 }
    }
    fn total_fragments(&self) -> i32 {
        self.split_frags + self.realigned_frags + self.discordant_frags
    }
    fn supporting_fragments(&self) -> i32 {
        self.split_frags + self.discordant_frags
    }

    pub fn header() -> String {
        let mut h: Vec<String> = vec!["FusionId".into(), "KnownFusionType".into()];
        for s in ["Up", "Down"] {
            for f in ["Chr", "Pos", "Orient", "GeneId", "GeneName", "Strand", "JuncType", "Transcript", "Exon", "Coverage", "MaxAnchorLength", "TransData"] {
                h.push(format!("{f}{s}"));
            }
        }
        for f in [
            "SVType",
            "RnaKnownFusionType",
            "TotalFragments",
            "SplitFrags",
            "RealignedFrags",
            "DiscordantFrags",
            "RelatedSplicedIds",
            "CohortCount",
            "ValidGenes",
            "Filter",
            "ReadType",
            "InitReadId",
        ] {
            h.push(f.into());
        }
        h.join("\t")
    }

    /// toTsv, with the final ids
    pub fn to_tsv(&self, final_id: &dyn Fn(i32) -> i32) -> String {
        let mut v: Vec<String> = vec![format!("Id_{}", final_id(self.id)), self.known_type.name()];
        for fs in [FS_UP, FS_DOWN] {
            v.push(self.chromosomes[fs].to_string());
            v.push(self.junction_positions[fs].to_string());
            v.push(self.junction_orientations[fs].to_string());
            v.push(self.gene_ids[fs].clone());
            v.push(self.gene_names[fs].clone());
            v.push(self.strands[fs].to_string());
            v.push(self.junction_types[fs].name().to_string());
            v.push(self.transcripts[fs].clone());
            v.push(self.exons[fs].to_string());
            v.push(self.coverage[fs].to_string());
            v.push(self.anchor_distance[fs].to_string());
            v.push(self.all_trans_data[fs].clone());
        }
        v.push(self.sv_type.to_string());
        v.push(self.mixed_known_type.name().to_string());
        v.push(self.total_frags.to_string());
        v.push(self.split_frags.to_string());
        v.push(self.realigned_frags.to_string());
        v.push(self.discordant_frags.to_string());
        if self.related_spliced.is_empty() {
            v.push("NONE".into());
        } else {
            let mut set: JavaHashSet<i32> = JavaHashSet::new();
            for id in &self.related_spliced {
                set.insert(final_id(*id));
            }
            v.push(set.iter().map(|x| format!("Id_{x}")).collect::<Vec<_>>().join(";"));
        }
        v.push(self.cohort_frequency.to_string());
        v.push(self.valid_genes.to_string());
        v.push(self.filter.name().to_string());
        v.push(self.read_type.to_string());
        v.push(self.init_read_id.clone());
        v.join("\t")
    }

    pub fn pass_header() -> &'static str {
        "Name\tKnownType\tChromosomeUp\tChromosomeDown\tPositionUp\tPositionDown\tOrientationUp\tOrientationDown\tJunctionTypeUp\tJunctionTypeDown\tTranscriptUp\tTranscriptDown\tExonUp\tExonDown\tSvType\tSplitFrags\tRealignedFrags\tDiscordantFrags\tDepthUp\tDepthDown\tMaxAnchorLengthUp\tMaxAnchorLengthDown\tCohortFrequency"
    }

    /// RnaFusionFile.write of buildRnaFusion
    pub fn pass_line(&self) -> String {
        [
            format!("{}_{}", self.gene_names[FS_UP], self.gene_names[FS_DOWN]),
            self.known_type.name(),
            self.chromosomes[0].to_string(),
            self.chromosomes[1].to_string(),
            self.junction_positions[0].to_string(),
            self.junction_positions[1].to_string(),
            self.junction_orientations[0].to_string(),
            self.junction_orientations[1].to_string(),
            self.junction_types[0].name().to_string(),
            self.junction_types[1].name().to_string(),
            self.transcripts[0].clone(),
            self.transcripts[1].clone(),
            self.exons[0].to_string(),
            self.exons[1].to_string(),
            self.sv_type.to_string(),
            self.split_frags.to_string(),
            self.realigned_frags.to_string(),
            self.discordant_frags.to_string(),
            self.coverage[0].to_string(),
            self.coverage[1].to_string(),
            self.anchor_distance[0].to_string(),
            self.anchor_distance[1].to_string(),
            self.cohort_frequency.to_string(),
        ]
        .join("\t")
    }
}

/// the fusion cohort file (PassingFusions' cohort fusions)
#[derive(Default)]
pub struct CohortFusions {
    by_pair: HashMap<String, HashMap<i32, Vec<([String; 2], [i32; 2], [i8; 2], i32)>>>,
}

impl CohortFusions {
    pub fn load(filename: Option<&str>) -> anyhow::Result<CohortFusions> {
        let mut c = CohortFusions::default();
        let Some(f) = filename else { return Ok(c) };
        let delim = if f.ends_with(".tsv") || f.ends_with(".tsv.gz") { '\t' } else { ',' };
        let mut reader = hmf_common::gz::LineReader::open(f)?;
        let header = match reader.next_line()? {
            Some(h) => std::str::from_utf8(h)?.to_string(),
            None => return Ok(c),
        };
        let fi = hmf_common::java_text::fields_index_map(&header, delim);
        let g = |n: &str| fi.get(n).copied().ok_or_else(|| anyhow::anyhow!("fusion cohort file has no {n}"));
        let (cu, cd, pu, pd, ou, od, cc) = (g("ChrUp")?, g("ChrDown")?, g("PosUp")?, g("PosDown")?, g("OrientUp")?, g("OrientDown")?, g("CohortCount")?);
        let mut count = 0;
        while let Some(line) = reader.next_line()? {
            let line = std::str::from_utf8(line)?;
            let v: Vec<&str> = hmf_common::java_text::split_all(line, delim);
            let chrs = [v[cu].to_string(), v[cd].to_string()];
            let pos = [v[pu].parse()?, v[pd].parse()?];
            let orient = [v[ou].parse()?, v[od].parse()?];
            let n: i32 = v[cc].parse()?;
            count += 1;
            let pair = format!("{}_{}", chrs[0], chrs[1]);
            c.by_pair.entry(pair).or_default().entry(pos[0]).or_default().push((chrs, pos, orient, n));
        }
        hmf_common::info!("loaded {} cohort fusions from file({})", count, f);
        Ok(c)
    }
    fn find(&self, d: &FusionData) -> Option<i32> {
        let pair = format!("{}_{}", d.chromosomes[0], d.chromosomes[1]);
        let l = self.by_pair.get(&pair)?.get(&d.junction_positions[0])?;
        l.iter().find(|(c, p, o, _)| c[0] == d.chromosomes[0] && c[1] == d.chromosomes[1] && *p == d.junction_positions && *o == d.junction_orientations).map(|x| x.3)
    }
}

const KNOWN_PAIR_KNOWN_SITE_REQ_FRAGS: i32 = 2;
const KNOWN_PAIR_NON_KNOWN_SITE_REQ_FRAGS: i32 = 4;
const MIN_ANCHOR_DISTANCE: i32 = 20;
const AF_KNOWN_SPLICED_TIER: f64 = 0.002;
const AF_KNOWN_TIER: f64 = 0.005;
const AF_UNKNOWN_TIER: f64 = 0.05;
const LOCAL_FUSION_THRESHOLD: i32 = 1_000_000;

fn is_short_local_fusion(f: &FusionData) -> bool {
    // ISOFOX 2.0.1 does not require opposite orientations
    f.chromosomes[0] == f.chromosomes[1]
        && (crate::version::v201() || f.junction_orientations[0] != f.junction_orientations[1])
        && (f.junction_positions[1] - f.junction_positions[0]).abs() < LOCAL_FUSION_THRESHOLD
}

fn has_known_splice_site(f: &FusionData) -> bool {
    use JunctionType::*;
    f.has_known_splice_sites() || (f.junction_types[0] == Canonical && f.junction_types[1] == Known) || (f.junction_types[0] == Known && f.junction_types[1] == Canonical)
}

fn is_passing_fusion(f: &mut FusionData, known_type: MixedKnownType, has_known_sites: bool) -> bool {
    use JunctionType::*;
    if known_type == MixedKnownType::KnownPair && !is_short_local_fusion(f) {
        let req = if has_known_sites || has_known_splice_site(f) { KNOWN_PAIR_KNOWN_SITE_REQ_FRAGS } else { KNOWN_PAIR_NON_KNOWN_SITE_REQ_FRAGS };
        if f.total_fragments() < req {
            f.filter = FilterType::FragmentCount;
            return false;
        }
        return true;
    }
    if f.anchor_distance[0].min(f.anchor_distance[1]) < MIN_ANCHOR_DISTANCE && f.discordant_frags == 0 {
        f.filter = FilterType::AnchorDistance;
        return false;
    }
    let (af, frags) = if has_known_sites {
        (AF_KNOWN_SPLICED_TIER, 2)
    } else if (f.junction_types[0] == Canonical && f.junction_types[1] == Known) || (f.junction_types[0] == Known && f.junction_types[1] == Canonical) {
        (AF_KNOWN_TIER, 3)
    } else if f.junction_types[0] == Canonical && f.junction_types[1] == Canonical {
        (AF_KNOWN_TIER, 4)
    } else {
        (AF_UNKNOWN_TIER, 10)
    };
    if f.allele_frequency() < af {
        f.filter = FilterType::AlleleFrequency;
        return false;
    }
    // ISOFOX 2.0.1 counts the realigned fragments too, not the discordant only
    let support = if crate::version::v201() { f.total_fragments() } else { f.supporting_fragments() };
    if support < frags {
        f.filter = FilterType::FragmentCount;
        return false;
    }
    let limit = if known_type.has_known_pair_gene() { FILTER_COHORT_LIMIT_KNOWN } else { FILTER_COHORT_LIMIT_NOT_KNOWN };
    if f.cohort_frequency >= limit {
        f.filter = FilterType::Cohort;
        return false;
    }
    true
}

fn mark_known_type(known: &KnownFusionCache, f: &mut FusionData) {
    if f.gene_names[0].is_empty() || f.gene_names[1].is_empty() {
        return;
    }
    for k in known.data_by_type(KnownFusionType::KNOWN_PAIR) {
        if k.five_gene == f.gene_names[0] && k.three_gene == f.gene_names[1] {
            f.known_type = KnownFusionType::KNOWN_PAIR;
            return;
        }
    }
    for k in known.data_by_type(KnownFusionType::PROMISCUOUS_5) {
        if k.five_gene == f.gene_names[0] {
            f.known_type = KnownFusionType::PROMISCUOUS_5;
            return;
        }
    }
    for k in known.data_by_type(KnownFusionType::PROMISCUOUS_3) {
        if k.three_gene == f.gene_names[1] {
            f.known_type = KnownFusionType::PROMISCUOUS_3;
            return;
        }
    }
}

fn mark_known_gene_types(known: &KnownFusionCache, f: &mut FusionData) {
    let mut is_known = [false, false];
    for k in known.data_by_type(KnownFusionType::KNOWN_PAIR) {
        if k.five_gene == f.gene_names[0] && k.three_gene == f.gene_names[1] {
            f.mixed_known_type = MixedKnownType::KnownPair;
            return;
        }
        if k.five_gene == f.gene_names[0] {
            is_known[0] = true;
        }
        if k.three_gene == f.gene_names[1] {
            is_known[1] = true;
        }
    }
    let prom = [!is_known[0] && known.has_promiscuous_five_gene(&f.gene_names[0]), !is_known[1] && known.has_promiscuous_three_gene(&f.gene_names[1])];
    f.mixed_known_type = if is_known[0] {
        if prom[1] { MixedKnownType::KnownProm3 } else { MixedKnownType::KnownOther }
    } else if is_known[1] {
        if prom[0] { MixedKnownType::Prom5Known } else { MixedKnownType::KnownOther }
    } else if prom[0] {
        if prom[1] { MixedKnownType::Prom5Prom3 } else { MixedKnownType::Prom5Other }
    } else if prom[1] {
        MixedKnownType::OtherProm3
    } else {
        MixedKnownType::Other
    };
}

/// PassingFusions.findPassingFusions: the indices of the passing fusions, in Java's order
pub fn find_passing_fusions(known: &KnownFusionCache, cohort: &CohortFusions, all: &mut [FusionData]) -> Vec<usize> {
    let mut passing: Vec<usize> = Vec::new();
    let mut non_passing_related: Vec<usize> = Vec::new();
    let mut low_support: SeqHashMap<String, Vec<usize>> = SeqHashMap::new();
    for i in 0..all.len() {
        let f = &mut all[i];
        mark_known_type(known, f);
        mark_known_gene_types(known, f);
        if let Some(c) = cohort.find(f) {
            f.cohort_frequency = c;
        }
        let (mk, ks) = (f.mixed_known_type, f.has_known_splice_sites());
        if is_passing_fusion(f, mk, ks) {
            passing.push(i);
        } else {
            let short_local = is_short_local_fusion(f);
            if !short_local && !f.related_spliced.is_empty() {
                non_passing_related.push(i);
            }
            if f.known_type == KnownFusionType::KNOWN_PAIR && f.filter == FilterType::FragmentCount && !short_local {
                let name = f.name();
                match low_support.get_mut(&name) {
                    Some(l) => l.push(i),
                    None => {
                        low_support.insert(name, vec![i]);
                    }
                }
            }
        }
    }
    for (_, l) in low_support.iter() {
        let known_sites: i32 = l.iter().filter(|i| has_known_splice_site(&all[**i])).map(|i| all[*i].total_fragments()).sum();
        let total: i32 = l.iter().map(|i| all[*i].total_fragments()).sum();
        if known_sites >= KNOWN_PAIR_KNOWN_SITE_REQ_FRAGS || total >= KNOWN_PAIR_NON_KNOWN_SITE_REQ_FRAGS {
            passing.extend(l.iter().copied());
        }
    }
    for i in non_passing_related {
        let mut matched = false;
        for pi in passing.clone() {
            if all[i].is_related(&all[pi]) && all[pi].has_known_splice_sites() {
                let (mk, ks) = (all[pi].mixed_known_type, all[pi].has_known_splice_sites());
                if is_passing_fusion(&mut all[i], mk, ks) {
                    matched = true;
                    break;
                }
            }
        }
        if matched {
            passing.push(i);
        }
    }
    for i in &passing {
        all[*i].filter = FilterType::Pass;
    }
    passing
}

// ---- FusionUtils: junction position adjustments ----

fn find_exon_boundary_match(t: &'static TranscriptData, pos: i32, orient: i8, buffer: i32) -> Option<(TransExonRef, i32, bool)> {
    let mut lo = pos;
    let mut hi = pos;
    if buffer > 0 {
        if orient == ORIENT_FWD {
            lo -= buffer;
        } else {
            hi += buffer;
        }
    }
    for e in &t.exons {
        let ep = if orient == ORIENT_FWD { e.end } else { e.start };
        if position_within(ep, lo, hi) {
            return Some((TransExonRef::new(t, e.rank), ep, pos == ep));
        }
    }
    None
}

fn check_fusion_position_adjustments_vs_known_exons(fusion: &mut FusionReadData, frags: &[FusionFragment], lists: &[Vec<&'static TranscriptData>; 2], refg: &mut RefSource) {
    let mut buffer = fusion.split_junction_overlap;
    let mut applied = false;
    let mut adjusts = [0, 0];
    for se in [SE_START, SE_END] {
        let mut matched: Vec<TransExonRef> = Vec::new();
        let mut top: Option<(TransExonRef, i32, bool)> = None;
        let mut pos = fusion.junction_positions[se];
        let orient = fusion.junction_orientations[se];
        if adjusts[se] != 0 {
            pos -= adjusts[se] * orient as i32;
        }
        for t in &lists[se] {
            let Some(m) = find_exon_boundary_match(t, pos, orient, buffer) else { continue };
            if top.is_none() || (m.2 && !top.unwrap().2) {
                top = Some(m);
                matched.clear();
                matched.push(m.0);
            } else if m.2 == top.unwrap().2 {
                matched.push(m.0);
            }
        }
        if !matched.is_empty() {
            fusion.trans_exon_refs[se] = matched;
        }
        if let Some(t) = top {
            fusion.junction_types[se] = JunctionType::Known;
            if !applied {
                if t.1 != pos {
                    let a = (t.1 - pos).abs();
                    adjusts[se] = a;
                    adjusts[switch_index(se)] = buffer - a;
                } else {
                    adjusts[switch_index(se)] = buffer;
                }
                buffer = 0;
                applied = true;
            }
        }
    }
    if applied {
        fusion.split_junction_overlap_applied = true;
        for se in [SE_START, SE_END] {
            if adjusts[se] != 0 {
                fusion.junction_positions[se] -= adjusts[se] * fusion.junction_orientations[se] as i32;
                fusion.set_junction_bases(frags, refg, se);
            }
        }
    }
}

fn matches_canonical_splice_junction(orient: i8, bases: &str, strand: Option<i8>) -> bool {
    match strand {
        Some(s) => {
            let donor = orient == s;
            (donor && canonical_donor(bases, s)) || (!donor && canonical_acceptor(bases, s))
        }
        None => canonical_donor(bases, orient) || canonical_acceptor(bases, -orient),
    }
}

fn check_fusion_position_adjustments_vs_canonical(fusion: &mut FusionReadData, frags: &[FusionFragment], refg: &mut RefSource) {
    for se in [SE_START, SE_END] {
        if fusion.junction_types[se] == JunctionType::Known {
            continue;
        }
        if matches_canonical_splice_junction(fusion.junction_orientations[se], &fusion.junction_splice_bases[se], fusion.gene_strand_by_position(se)) {
            fusion.junction_types[se] = JunctionType::Canonical;
        }
    }
    if fusion.split_junction_overlap == 0 || fusion.split_junction_overlap_applied {
        return;
    }
    let adj = fusion.split_junction_overlap;
    let (s, e) = (fusion.junction_types[0].ordinal(), fusion.junction_types[1].ordinal());
    let mut adjusts = [0, 0];
    if s < e {
        adjusts[SE_END] = adj;
    } else if e < s {
        adjusts[SE_START] = adj;
    } else {
        let half = adj / 2;
        adjusts[SE_START] = if adj % 2 == 0 { half } else { half + 1 };
        adjusts[SE_END] = adj - adjusts[SE_START];
    }
    for se in [SE_START, SE_END] {
        if adjusts[se] != 0 {
            fusion.junction_positions[se] -= adjusts[se] * fusion.junction_orientations[se] as i32;
            fusion.set_junction_bases(frags, refg, se);
            if fusion.junction_types[se] != JunctionType::Known {
                fusion.junction_types[se] = if matches_canonical_splice_junction(fusion.junction_orientations[se], &fusion.junction_splice_bases[se], fusion.gene_strand_by_position(se)) {
                    JunctionType::Canonical
                } else {
                    JunctionType::Unknown
                };
            }
        }
    }
}

/// the fusions written by one batch, with their ids still local to the chromosome
pub struct FusionBatchOutput {
    pub fusions: Vec<FusionData>,
    pub passing: Vec<usize>,
}

/// FusionFinder: the fusion candidates of one chromosome
#[allow(dead_code)]
pub struct FusionFinder {
    pub chromosome: &'static str,
    /// the next id, local to this chromosome
    next_id: i32,
    pub local_fusion_positions: std::collections::HashSet<String>,
    pub hard_filtered_count: i32,
    pub excluded_filtered_count: i32,
    pub duplicate_filtered_count: i32,
    pub outputs: Vec<FusionBatchOutput>,
}

/// the fusion candidates of one call to processReadGroups
struct Batch {
    frags: Vec<FusionFragment>,
    all: Vec<usize>,
    candidates: SeqHashMap<String, Vec<usize>>,
    fusions: Vec<Option<FusionReadData>>,
    by_location: HashMap<String, HashMap<String, usize>>,
    discordant: SeqHashMap<String, Vec<usize>>,
}

impl FusionFinder {
    pub fn new(chromosome: &'static str) -> Self {
        FusionFinder {
            chromosome,
            next_id: 0,
            local_fusion_positions: Default::default(),
            hard_filtered_count: 0,
            excluded_filtered_count: 0,
            duplicate_filtered_count: 0,
            outputs: Vec::new(),
        }
    }

    pub fn ids_used(&self) -> i32 {
        self.next_id
    }

    /// processReadGroups
    pub fn process_read_groups(&mut self, ctx: &FusionContext, refg: &mut RefSource, groups: Vec<FusionReadGroup>, rac: &mut dyn FnMut(&str, i32) -> Option<*mut JunctionRacFragments>) {
        if groups.is_empty() {
            return;
        }
        let mut b = Batch { frags: Vec::new(), all: Vec::new(), candidates: SeqHashMap::new(), fusions: Vec::new(), by_location: HashMap::new(), discordant: SeqHashMap::new() };
        let v201 = crate::version::v201();
        for g in groups {
            // ISOFOX 2.0.1: a group with a duplicate read is dropped once complete (not all reads are marked)
            if v201 && g.has_duplicate_read() {
                self.duplicate_filtered_count += 1;
                continue;
            }
            if g.reads.iter().any(|r| (v201 || r.is_read_paired()) && ctx.filters.skip_read_start(r.mate_chromosome, r.mate_pos_start)) {
                self.excluded_filtered_count += 1;
                continue;
            }
            let f = FusionFragment::new(g);
            if f.frag_type == FragType::Unknown {
                continue;
            }
            b.frags.push(f);
            b.all.push(b.frags.len() - 1);
        }
        self.form_initial_fusions(ctx, refg, &mut b);
        self.reconcile_fusions(&mut b);
        self.assign_discordant_fragments(ctx, &mut b);
        self.create_local_fusions(ctx, refg, &mut b);
        self.assign_realign_candidate_fragments(&mut b, rac);
        self.hard_filter_fusions(ctx, &mut b);
        self.check_local_duplicates(&mut b);
        self.mark_related_fusions(&mut b);
        self.write_data(ctx, &mut b);
    }

    fn form_initial_fusions(&mut self, ctx: &FusionContext, refg: &mut RefSource, b: &mut Batch) {
        let all = std::mem::take(&mut b.all);
        for fi in all {
            match b.frags[fi].frag_type {
                FragType::MatchedJunction => {
                    if b.frags[fi].has_supp {
                        self.create_or_update_fusion(ctx, refg, b, fi);
                    } else {
                        cache_discordant(b, fi);
                    }
                }
                FragType::RealignCandidate => {}
                _ => cache_discordant(b, fi),
            }
        }
        b.by_location.clear();
    }

    fn find_existing_fusion(b: &Batch, fi: usize) -> Option<usize> {
        let f = &b.frags[fi];
        let pair = format!("{}_{}", f.chromosomes[0], f.chromosomes[1]);
        b.by_location.get(&pair)?.get(&f.position_hash()).copied()
    }

    fn create_or_update_fusion(&mut self, ctx: &FusionContext, refg: &mut RefSource, b: &mut Batch, fi: usize) -> Option<usize> {
        if let Some(existing) = Self::find_existing_fusion(b, fi) {
            let mut fusion = b.fusions[existing].take().unwrap();
            fusion.add_fusion_fragment(&b.frags, fi);
            b.fusions[existing] = Some(fusion);
            return None;
        }
        let loc = b.frags[fi].location_pair();
        if !b.candidates.contains_key(&loc) {
            b.candidates.insert(loc.clone(), Vec::new());
        }
        let id = self.next_id;
        self.next_id += 1;
        let mut fusion = FusionReadData::new(id, &b.frags, fi);
        fusion.set_junction_bases_all(&b.frags, refg);
        if crate::version::v201() {
            fusion.set_homology_offsets();
            self.set_gene_data_201(ctx, &mut fusion, &mut b.frags);
        } else {
            fusion.check_homology_position_adjustment(&b.frags);
            self.set_gene_data(ctx, refg, &mut fusion, &mut b.frags);
            fusion.set_homology_offsets();
        }
        b.fusions.push(Some(fusion));
        let idx = b.fusions.len() - 1;
        b.candidates.get_mut(&loc).unwrap().push(idx);
        let f = &b.frags[fi];
        let pair = format!("{}_{}", f.chromosomes[0], f.chromosomes[1]);
        b.by_location.entry(pair).or_default().insert(f.position_hash(), idx);
        Some(idx)
    }

    fn can_create_discordant_fusion(&mut self, ctx: &FusionContext, refg: &mut RefSource, b: &mut Batch, fi: usize) -> bool {
        if b.frags[fi].frag_type != FragType::DiscordantJunction {
            return false;
        }
        if Self::find_existing_fusion(b, fi).is_some() {
            return true;
        }
        let mut fusion = FusionReadData::new(0, &b.frags, fi);
        fusion.set_junction_bases_all(&b.frags, refg);
        if crate::version::v201() {
            fusion.set_homology_offsets();
            self.set_gene_data_201(ctx, &mut fusion, &mut b.frags);
        } else {
            self.set_gene_data(ctx, refg, &mut fusion, &mut b.frags);
            fusion.set_homology_offsets();
        }
        if !fusion.has_viable_genes() {
            return false;
        }
        if !ctx.known.has_known_fusion(fusion.gene_name(FS_UP), fusion.gene_name(FS_DOWN)) {
            return false;
        }
        if crate::version::v201() {
            // the fragment's own junction types, set by setGeneData
            let jt = b.frags[fi].junction_types;
            return jt[FS_UP] == JunctionType::Known || jt[FS_DOWN] == JunctionType::Known;
        }
        fusion.junction_types[FS_UP] == JunctionType::Known || fusion.junction_types[FS_DOWN] == JunctionType::Known
    }

    fn set_gene_data(&mut self, ctx: &FusionContext, refg: &mut RefSource, fusion: &mut FusionReadData, frags: &mut [FusionFragment]) {
        let fi = fusion.initial;
        let mut genes_by_pos: [Vec<&'static GeneData>; 2] = [Vec::new(), Vec::new()];
        let mut lists: [Vec<&'static TranscriptData>; 2] = [Vec::new(), Vec::new()];
        for se in [SE_START, SE_END] {
            let mut splice_gene_ids: JavaHashSet<String> = JavaHashSet::new();
            let refs = frags[fi].trans_exon_refs[se].clone();
            for te in refs {
                let Some(t) = ctx.trans_by_id.get(&te.trans_id).copied() else { continue };
                if !lists[se].iter().any(|x| std::ptr::eq(*x, t)) {
                    lists[se].push(t);
                }
                fusion.trans_exon_refs[se].push(TransExonRef::new(t, te.exon_rank));
                splice_gene_ids.insert(t.gene_id.clone());
            }
            frags[fi].validate_transcript_exons(&lists[se], se);
            if !splice_gene_ids.is_empty() {
                genes_by_pos[se] = splice_gene_ids.iter().filter_map(|g| ctx.cache.gene_data_by_id(g)).collect();
            }
        }
        check_fusion_position_adjustments_vs_known_exons(fusion, frags, &lists, refg);
        let mut valid_lists: [Vec<&'static TranscriptData>; 2] = [Vec::new(), Vec::new()];
        for se in [SE_START, SE_END] {
            for t in &lists[se] {
                if fusion.trans_exon_refs[se].iter().any(|x| x.trans_id == t.trans_id) {
                    valid_lists[se].push(t);
                }
            }
        }
        if genes_by_pos[0].len() > 1 || genes_by_pos[1].len() > 1 {
            let mut matched = false;
            'outer: for g1 in genes_by_pos[0].clone() {
                for g2 in genes_by_pos[1].clone() {
                    if ctx.known.has_known_fusion(&g1.gene_name, &g2.gene_name) || ctx.known.has_known_fusion(&g2.gene_name, &g1.gene_name) {
                        matched = true;
                        genes_by_pos[0] = vec![g1];
                        genes_by_pos[1] = vec![g2];
                        break 'outer;
                    }
                }
            }
            if !matched {
                for se in [SE_START, SE_END] {
                    if genes_by_pos[se].len() > 1 {
                        prioritise_known_fusion_gene(ctx.known, &mut genes_by_pos[se]);
                    }
                }
            }
        }
        let orients = fusion.junction_orientations;
        let mut found_both = false;
        let mut found_one = false;
        for se in [SE_START, SE_END] {
            let up_i = se;
            let down_i = switch_index(se);
            let mut up: Vec<&'static GeneData> = genes_by_pos[up_i].iter().copied().filter(|g| g.strand == orients[up_i]).collect();
            let mut down: Vec<&'static GeneData> = genes_by_pos[down_i].iter().copied().filter(|g| g.strand == -orients[down_i]).collect();
            prioritise_longest_coding_fusion_gene(&mut up, &valid_lists[up_i]);
            prioritise_longest_coding_fusion_gene(&mut down, &valid_lists[down_i]);
            if !up.is_empty() && !down.is_empty() {
                if found_both {
                    fusion.incomplete = true;
                    break;
                }
                found_both = true;
                fusion.set_stream_data(&up, &down, se == SE_START);
            } else if !found_both && !found_one && (!up.is_empty() || !down.is_empty()) {
                found_one = true;
                fusion.set_stream_data(&up, &down, se == SE_START);
            }
        }
        check_fusion_position_adjustments_vs_canonical(fusion, frags, refg);
    }

    /// setGeneData in ISOFOX 2.0.1: no junction position adjustments; the fusion keeps every transcript-exon
    /// of the initial fragment; the fragment's junction types are set (KNOWN by the exon validation, then
    /// CANONICAL by the splice motif)
    fn set_gene_data_201(&mut self, ctx: &FusionContext, fusion: &mut FusionReadData, frags: &mut [FusionFragment]) {
        let fi = fusion.initial;
        let mut genes_by_pos: [Vec<&'static GeneData>; 2] = [Vec::new(), Vec::new()];
        let mut valid_lists: [Vec<&'static TranscriptData>; 2] = [Vec::new(), Vec::new()];
        for se in [SE_START, SE_END] {
            let mut list: Vec<&'static TranscriptData> = Vec::new();
            let mut splice_gene_ids: JavaHashSet<String> = JavaHashSet::new();
            let refs = frags[fi].trans_exon_refs[se].clone();
            for te in refs {
                let Some(t) = ctx.trans_by_id.get(&te.trans_id).copied() else { continue };
                if !list.iter().any(|x| std::ptr::eq(*x, t)) {
                    list.push(t);
                }
                fusion.trans_exon_refs[se].push(TransExonRef::new(t, te.exon_rank));
                splice_gene_ids.insert(t.gene_id.clone());
            }
            frags[fi].validate_transcript_exons(&list, se);
            if !splice_gene_ids.is_empty() {
                genes_by_pos[se] = splice_gene_ids.iter().filter_map(|g| ctx.cache.gene_data_by_id(g)).collect();
                valid_lists[se] = list.iter().copied().filter(|t| frags[fi].trans_exon_refs[se].iter().any(|y| y.trans_id == t.trans_id)).collect();
            }
        }
        if genes_by_pos[0].len() > 1 || genes_by_pos[1].len() > 1 {
            let mut matched = false;
            'outer: for g1 in genes_by_pos[0].clone() {
                for g2 in genes_by_pos[1].clone() {
                    if ctx.known.has_known_fusion(&g1.gene_name, &g2.gene_name) || ctx.known.has_known_fusion(&g2.gene_name, &g1.gene_name) {
                        matched = true;
                        genes_by_pos[0] = vec![g1];
                        genes_by_pos[1] = vec![g2];
                        break 'outer;
                    }
                }
            }
            if !matched {
                for se in [SE_START, SE_END] {
                    if genes_by_pos[se].len() > 1 {
                        prioritise_known_fusion_gene(ctx.known, &mut genes_by_pos[se]);
                    }
                }
            }
        }
        let orients = fusion.junction_orientations;
        let mut found_both = false;
        let mut found_one = false;
        for se in [SE_START, SE_END] {
            let up_i = se;
            let down_i = switch_index(se);
            let mut up: Vec<&'static GeneData> = genes_by_pos[up_i].iter().copied().filter(|g| g.strand == orients[up_i]).collect();
            let mut down: Vec<&'static GeneData> = genes_by_pos[down_i].iter().copied().filter(|g| g.strand == -orients[down_i]).collect();
            prioritise_longest_coding_fusion_gene(&mut up, &valid_lists[up_i]);
            prioritise_longest_coding_fusion_gene(&mut down, &valid_lists[down_i]);
            if !up.is_empty() && !down.is_empty() {
                if found_both {
                    fusion.incomplete = true;
                    break;
                }
                found_both = true;
                fusion.set_stream_data(&up, &down, se == SE_START);
            } else if !found_both && !found_one && (!up.is_empty() || !down.is_empty()) {
                found_one = true;
                fusion.set_stream_data(&up, &down, se == SE_START);
            }
        }
        // getGeneStrands: the genes' strands by position, when both genes are known
        let strands = if fusion.has_viable_genes() {
            let (up, down) = (fusion.fusion_genes[FS_UP].unwrap().strand, fusion.fusion_genes[FS_DOWN].unwrap().strand);
            Some(if fusion.stream_indices[FS_UP] == SE_START { [up, down] } else { [down, up] })
        } else {
            None
        };
        let splice = fusion.junction_splice_bases.clone();
        frags[fi].set_junction_types_201(strands, &splice);
    }

    fn reconcile_fusions(&mut self, b: &mut Batch) {
        let keys = b.candidates.keys_vec();
        for k in keys {
            let mut list = b.candidates.get(&k).unwrap().clone();
            if list.len() == 1 {
                continue;
            }
            let mut i = 0;
            while i + 1 < list.len() {
                let mut f1 = list[i];
                let known_spliced = b.fusions[f1].as_ref().unwrap().is_known_spliced(&b.frags);
                let mut j = i + 1;
                while j < list.len() {
                    let f2 = list[j];
                    let (hm, close, f2_known) = {
                        let a = b.fusions[f1].as_ref().unwrap();
                        let c = b.fusions[f2].as_ref().unwrap();
                        (a.match_within_homology(c), a.is_close_match(c), c.is_known_spliced(&b.frags))
                    };
                    if hm {
                        let n = b.fusions[f2].as_ref().unwrap().frag_type_count(FragType::MatchedJunction);
                        b.fusions[f1].as_mut().unwrap().add_frag_type_count(FragType::MatchedJunction, n);
                        list.remove(j);
                        continue;
                    }
                    if !(known_spliced && f2_known) && close {
                        let c1 = b.fusions[f1].as_ref().unwrap().frag_type_count(FragType::MatchedJunction);
                        let c2 = b.fusions[f2].as_ref().unwrap().frag_type_count(FragType::MatchedJunction);
                        if c1 < c2 {
                            list[i] = f2;
                            f1 = f2;
                        }
                        list.remove(j);
                    } else {
                        j += 1;
                    }
                }
                i += 1;
            }
            *b.candidates.get_mut(&k).unwrap() = list;
        }
    }

    fn assign_discordant_to(&self, ctx: &FusionContext, b: &mut Batch, key: &str, fusions: &[usize]) {
        let frags = b.discordant.get(&key.to_string()).cloned().unwrap_or_default();
        let mut allocated: Vec<usize> = Vec::new();
        for fi in &frags {
            if b.frags[*fi].trans_exon_refs[0].is_empty() || b.frags[*fi].trans_exon_refs[1].is_empty() {
                continue;
            }
            for &fu in fusions {
                let (can_add, can_realign) = {
                    let fusion = b.fusions[fu].as_ref().unwrap();
                    let f = &b.frags[*fi];
                    let ca = fusion.can_add_discordant_fragment(&b.frags, f, ctx.max_fragment_length);
                    (ca, ca && fusion.can_realign_fragment_to_junction(&b.frags, f))
                };
                if can_add {
                    b.frags[*fi].frag_type = if can_realign { FragType::Realigned } else { FragType::Discordant };
                    let mut fusion = b.fusions[fu].take().unwrap();
                    fusion.add_fusion_fragment(&b.frags, *fi);
                    b.fusions[fu] = Some(fusion);
                    if !allocated.contains(fi) {
                        allocated.push(*fi);
                    }
                }
            }
        }
        if let Some(l) = b.discordant.get_mut(&key.to_string()) {
            l.retain(|x| !allocated.contains(x));
        }
    }

    fn assign_discordant_fragments(&mut self, ctx: &FusionContext, b: &mut Batch) {
        for key in b.discordant.keys_vec() {
            let Some(fusions) = b.candidates.get(&key).cloned() else { continue };
            self.assign_discordant_to(ctx, b, &key, &fusions);
        }
    }

    fn create_local_fusions(&mut self, ctx: &FusionContext, refg: &mut RefSource, b: &mut Batch) {
        let mut new_fusions: Vec<usize> = Vec::new();
        for key in b.discordant.keys_vec() {
            let frags = b.discordant.get(&key).cloned().unwrap_or_default();
            let mut allocated: Vec<usize> = Vec::new();
            for fi in frags {
                let t = b.frags[fi].frag_type;
                if t != FragType::MatchedJunction && t != FragType::DiscordantJunction {
                    continue;
                }
                if t == FragType::DiscordantJunction && !self.can_create_discordant_fusion(ctx, refg, b, fi) {
                    continue;
                }
                let Some(fu) = self.create_or_update_fusion(ctx, refg, b, fi) else { continue };
                new_fusions.push(fu);
                allocated.push(fi);
            }
            if let Some(l) = b.discordant.get_mut(&key) {
                l.retain(|x| !allocated.contains(x));
            }
        }
        for key in b.discordant.keys_vec() {
            let fusions: Vec<usize> = new_fusions.iter().copied().filter(|f| b.fusions[*f].as_ref().unwrap().location_id == key).collect();
            if fusions.is_empty() {
                continue;
            }
            self.assign_discordant_to(ctx, b, &key, &fusions);
        }
    }

    fn assign_realign_candidate_fragments(&mut self, b: &mut Batch, rac: &mut dyn FnMut(&str, i32) -> Option<*mut JunctionRacFragments>) {
        for key in b.candidates.keys_vec() {
            let list = b.candidates.get(&key).cloned().unwrap();
            for fu in list {
                for se in [SE_START, SE_END] {
                    let (chr, gc, pos, orient) = {
                        let f = b.fusions[fu].as_ref().unwrap();
                        (f.chromosomes[se], b.frags[f.initial].gene_collections[se], f.junction_positions[se], f.junction_orientations[se])
                    };
                    let Some(ptr) = rac(chr, gc) else { continue };
                    // SAFETY: the pointer is to a JunctionRacFragments the caller owns for the duration of
                    // this call, and no other reference to it is live here
                    let racs: &mut JunctionRacFragments = unsafe { &mut *ptr };
                    let Some(list) = racs.junction_fragments(orient, pos).cloned() else { continue };
                    for ri in list {
                        let can = b.fusions[fu].as_ref().unwrap().can_realign_fragment_to_junction(&b.frags, &racs.fragments[ri]);
                        if can {
                            racs.fragments[ri].frag_type = FragType::Realigned;
                            b.fusions[fu].as_mut().unwrap().add_rac_fragment(&racs.fragments[ri]);
                        }
                    }
                }
            }
        }
    }

    fn hard_filter_fusions(&mut self, ctx: &FusionContext, b: &mut Batch) {
        if crate::version::v201() && ctx.min_hard_filter_frags <= 1 {
            return;
        }
        for key in b.candidates.keys_vec() {
            let mut list = b.candidates.get(&key).cloned().unwrap();
            let mut i = 0;
            while i < list.len() {
                if self.hard_filter_fusion(ctx, b, list[i]) {
                    self.hard_filtered_count += 1;
                    list.remove(i);
                } else {
                    i += 1;
                }
            }
            *b.candidates.get_mut(&key).unwrap() = list;
        }
    }

    fn hard_filter_fusion(&self, ctx: &FusionContext, b: &Batch, fu: usize) -> bool {
        let f = b.fusions[fu].as_ref().unwrap();
        let (up, down) = (f.gene_name(FS_UP), f.gene_name(FS_DOWN));
        if crate::version::v201() {
            // ISOFOX 2.0.1: no map quality filter; all fragments count; the initial fragment's junction types
            if ctx.min_hard_filter_frags <= 1 {
                return false;
            }
            let total: i32 = f.frag_counts.iter().map(|x| x.1).sum();
            if total >= ctx.min_hard_filter_frags {
                return false;
            }
            if ctx.known.has_known_fusion(up, down) {
                return false;
            }
            let jt = b.frags[f.initial].junction_types;
            return jt[SE_START] != JunctionType::Known && jt[SE_END] != JunctionType::Known;
        }
        let known_pair = ctx.known.has_known_fusion(up, down);
        let promiscuous = !known_pair && (ctx.known.has_promiscuous_five_gene(up) || ctx.known.has_promiscuous_three_gene(down));
        if !known_pair && !promiscuous {
            // getFragments(MATCHED_JUNCTION) with fragments not cached is the initial fragment where it is a
            // split one
            let initial = &b.frags[f.initial];
            if initial.frag_type == FragType::MatchedJunction {
                let mut valid_mq = false;
                let mut has_supps = false;
                for r in &initial.reads {
                    if r.is_supplementary() {
                        has_supps = true;
                        if r.map_quality >= FILTER_MIN_MAP_QUAL {
                            valid_mq = true;
                            break;
                        }
                    }
                }
                if has_supps && !valid_mq {
                    return true;
                }
            }
        }
        if ctx.min_hard_filter_frags <= 1 {
            return false;
        }
        if f.supporting_fragment_count() >= ctx.min_hard_filter_frags {
            return false;
        }
        if known_pair {
            return false;
        }
        f.junction_types[0] != JunctionType::Known && f.junction_types[1] != JunctionType::Known
    }

    fn check_local_duplicates(&mut self, b: &mut Batch) {
        for key in b.candidates.keys_vec() {
            let mut list = b.candidates.get(&key).cloned().unwrap();
            let mut i = 0;
            while i < list.len() {
                let f = b.fusions[list[i]].as_ref().unwrap();
                let initial = &b.frags[f.initial];
                if !initial.has_supp {
                    let pair = initial.position_hash();
                    if self.local_fusion_positions.contains(&pair) {
                        list.remove(i);
                        continue;
                    }
                    self.local_fusion_positions.insert(pair);
                }
                i += 1;
            }
            *b.candidates.get_mut(&key).unwrap() = list;
        }
    }

    fn mark_related_fusions(&mut self, b: &mut Batch) {
        for key in b.candidates.keys_vec() {
            let list = b.candidates.get(&key).cloned().unwrap();
            if list.len() <= 1 {
                continue;
            }
            // IF8: Java's inner loop stops one short, so the last fusion of a location is never compared
            let end = if compat::java("IF8") { list.len().saturating_sub(1) } else { list.len() };
            for i in 0..list.len() - 1 {
                let (spliced, unspliced, up1, down1, id1) = {
                    let f1 = b.fusions[list[i]].as_ref().unwrap();
                    (f1.is_known_spliced(&b.frags), f1.is_unspliced(&b.frags), f1.trans_exon_refs_by_stream(FS_UP).clone(), f1.trans_exon_refs_by_stream(FS_DOWN).clone(), f1.id)
                };
                for j in i + 1..end {
                    let (f2_unspliced, f2_spliced, up2, down2, id2, jp2) = {
                        let f2 = b.fusions[list[j]].as_ref().unwrap();
                        (f2.is_unspliced(&b.frags), f2.is_known_spliced(&b.frags), f2.trans_exon_refs_by_stream(FS_UP).clone(), f2.trans_exon_refs_by_stream(FS_DOWN).clone(), f2.id, f2.junction_positions)
                    };
                    if spliced && f2_unspliced {
                        if TransExonRef::has_match(&up1, &up2) && TransExonRef::has_match_within_range(&down1, &down2, -1) {
                            b.fusions[list[j]].as_mut().unwrap().add_related(id1, true);
                            continue;
                        }
                    } else if unspliced && f2_spliced && TransExonRef::has_match(&up1, &up2) && TransExonRef::has_match_within_range(&down2, &down1, -1) {
                        b.fusions[list[i]].as_mut().unwrap().add_related(id2, true);
                        continue;
                    }
                    let jp1 = b.fusions[list[i]].as_ref().unwrap().junction_positions;
                    let similar = (0..2).any(|se| (jp1[se] - jp2[se]).abs() <= POSITION_REALIGN_DISTANCE);
                    if similar {
                        b.fusions[list[i]].as_mut().unwrap().add_related(id2, false);
                        b.fusions[list[j]].as_mut().unwrap().add_related(id1, false);
                    }
                }
            }
        }
    }

    fn write_data(&mut self, ctx: &FusionContext, b: &mut Batch) {
        if b.candidates.is_empty() {
            return;
        }
        let mut all: Vec<FusionData> = Vec::new();
        for (_, list) in b.candidates.iter() {
            for fu in list {
                all.push(b.fusions[*fu].as_ref().unwrap().to_fusion_data(&b.frags));
            }
        }
        let passing = find_passing_fusions(ctx.known, ctx.cohort, &mut all);
        self.outputs.push(FusionBatchOutput { fusions: all, passing });
    }
}

fn cache_discordant(b: &mut Batch, fi: usize) {
    let loc = b.frags[fi].location_pair();
    match b.discordant.get_mut(&loc) {
        Some(l) => l.push(fi),
        None => {
            b.discordant.insert(loc, vec![fi]);
        }
    }
}

fn prioritise_longest_coding_fusion_gene(genes: &mut Vec<&'static GeneData>, trans: &[&'static TranscriptData]) {
    if genes.len() <= 1 {
        return;
    }
    let mut longest: Option<&'static GeneData> = None;
    let mut longest_len = 0;
    for g in genes.iter() {
        let m = trans
            .iter()
            .filter(|t| t.gene_id == g.gene_id)
            .filter_map(|t| match (t.coding_start, t.coding_end) {
                (Some(s), Some(e)) => Some(e - s),
                _ => None,
            })
            .max()
            .unwrap_or(0);
        if m > longest_len {
            longest_len = m;
            longest = Some(g);
        }
    }
    genes.clear();
    // Java adds a null where no gene has coding bases; the stream data then takes it as no gene
    if let Some(g) = longest {
        genes.push(g);
    } else {
        genes.push(NULL_GENE.get_or_init(null_gene));
    }
}

/// stands in for the null GeneData Java puts in the list (prioritiseLongestCodingFusionGene)
static NULL_GENE: std::sync::OnceLock<GeneData> = std::sync::OnceLock::new();

fn null_gene() -> GeneData {
    GeneData {
        gene_id: String::new(),
        gene_name: String::new(),
        chromosome: String::new(),
        strand: 0,
        gene_start: 0,
        gene_end: 0,
        karyotype_band: String::new(),
        synonyms: String::new(),
    }
}

fn prioritise_known_fusion_gene(known: &KnownFusionCache, genes: &mut Vec<&'static GeneData>) {
    if genes.len() <= 1 {
        return;
    }
    let mut culled: Vec<&'static GeneData> = genes.iter().copied().filter(|g| known.has_known_pair_gene(&g.gene_name)).collect();
    if culled.is_empty() {
        culled = genes.iter().copied().filter(|g| known.has_promiscuous_five_gene(&g.gene_name) || known.has_promiscuous_three_gene(&g.gene_name)).collect();
    }
    if culled.len() == 1 {
        *genes = culled;
    }
}

pub fn is_null_gene(g: &GeneData) -> bool {
    g.gene_id.is_empty() && g.strand == 0
}

#[allow(dead_code)]
fn unused(_: &IsofoxConfig, _: &JavaHashMap<String, i32>) {}
