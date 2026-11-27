//! Package `annotators`: fragile sites, known and suspected LINE elements (LineElementAnnotator,
//! LineClusterState) and pseudogene insertions (PseudoGeneFinder).

use hmf_common::java_hash::JavaHashMap;
use hmf_common::java_id_hash::IdHashMap;
use hmf_common::regions::ChrBaseRegion;
use hmf_common::sv::StructuralVariantType::*;

use crate::gene::{BreakendGene, GeneCache};
use crate::types::*;
use crate::utils::*;
use crate::vis::{VisGeneAnnotationType, VisGeneData, VisSampleData};

const FRAGILE_SITES_37: &str = include_str!("../resources/fragile_sites.37.csv");
const FRAGILE_SITES_38: &str = include_str!("../resources/fragile_sites.38.csv");
const LINE_ELEMENTS_37: &str = include_str!("../resources/line_elements.37.csv");
const LINE_ELEMENTS_38: &str = include_str!("../resources/line_elements.38.csv");

fn resource_or_file(filename: Option<&str>, v37: &'static str, v38: &'static str, what: &str) -> Vec<String> {
    match filename {
        Some(f) => match hmf_common::java_text::read_lines(f) {
            Ok(l) => l,
            Err(e) => {
                hmf_common::error!("failed to read {} file: {}", what, e);
                Vec::new()
            }
        },
        None => hmf_common::java_text::lines_of(if ref_genome_is_v37() { v37 } else { v38 }),
    }
}

/// FragileSiteAnnotator
#[derive(Default)]
pub struct FragileSiteAnnotator {
    fragile_sites: Vec<ChrBaseRegion>,
}

impl FragileSiteAnnotator {
    pub fn load_fragile_sites_file(&mut self, filename: Option<&str>) {
        let lines = resource_or_file(filename, FRAGILE_SITES_37, FRAGILE_SITES_38, "fragile site");
        if let Some(r) = load_config_file(&lines, Some(ref_genome_is_v37()), ',') {
            self.fragile_sites.extend(r);
        }
        hmf_common::debug!("loaded {} known fragile sites", self.fragile_sites.len());
    }
    pub fn is_fragile_site(&self, sv: &SvVarData, use_start: bool) -> bool {
        self.fragile_sites.iter().any(|r| r.chromosome == sv.chromosome(use_start) && r.contains_position(sv.position(use_start)))
    }
}

pub const LINE_ELEMENT_PROXIMITY_DISTANCE: i32 = 5000;
pub const POLY_A_MOTIF: &str = "AAAAAAAAAA";
pub const POLY_T_MOTIF: &str = "TTTTTTTTTT";

/// LineElementAnnotator
pub struct LineElementAnnotator {
    known_line_elements: Vec<ChrBaseRegion>,
    proximity_distance: i32,
}

impl LineElementAnnotator {
    pub fn new(proximity_distance: i32) -> Self {
        LineElementAnnotator { known_line_elements: Vec::new(), proximity_distance }
    }

    pub fn load_line_elements_file(&mut self, filename: Option<&str>) {
        let lines = resource_or_file(filename, LINE_ELEMENTS_37, LINE_ELEMENTS_38, "line element");
        if let Some(r) = load_config_file(&lines, Some(ref_genome_is_v37()), ',') {
            self.known_line_elements.extend(r);
        }
        hmf_common::debug!("loaded {} known line elements", self.known_line_elements.len());
    }

    pub fn set_known_line_elements(&self, sv: &SvVarData) {
        if self.known_line_elements.is_empty() {
            return;
        }
        for se in [SE_START, SE_END] {
            let s = is_start(se);
            for r in &self.known_line_elements {
                if r.chromosome != sv.chromosome(s) {
                    continue;
                }
                if position_within(sv.position(s), r.start - LINE_ELEMENT_PROXIMITY_DISTANCE, r.end + LINE_ELEMENT_PROXIMITY_DISTANCE) {
                    hmf_common::debug!("var({}) found in known line element({} -> {})", sv.pos_id(), r.chromosome, r.start);
                    sv.add_line_element(LineElementType::KNOWN, s);
                }
            }
        }
    }

    /// markLineCluster
    pub fn mark_line_cluster(&self, cluster: &Cluster, pseudo_gene_finder: Option<&PseudoGeneFinder>) {
        if crate::analysis::classification::is_filtered_resolved_type(cluster.resolved_type()) {
            return;
        }
        let mut state = LineClusterState::new(self.proximity_distance);
        let mut has_suspected = false;
        for (_, list) in cluster.chr_breakend_entries() {
            state.clear();
            for i in 0..list.len() {
                let breakend = &list[i];
                if breakend.has_line_element(LineElementType::SUSPECT) {
                    continue;
                }
                state.add_breakend(breakend);
                if !state.is_suspected() {
                    continue;
                }
                let mut proximate: Vec<Be> = vec![breakend.clone()];
                for next in list.iter().skip(i + 1) {
                    if (next.position() - breakend.position()).abs() > self.proximity_distance {
                        break;
                    }
                    proximate.push(next.clone());
                }
                for j in (0..i).rev() {
                    let prev = &list[j];
                    if (breakend.position() - prev.position()).abs() > self.proximity_distance {
                        break;
                    }
                    proximate.push(prev.clone());
                }
                if let Some(pgf) = pseudo_gene_finder {
                    if let Some(pd) = proximate.iter().filter(|x| x.sv_type() == DEL).find(|x| pgf.variant_matches_pseudogene_exons(x.sv())) {
                        hmf_common::debug!("cluster({}) proximate DEL({}) matches pseudogene exon boundaries", cluster.id(), pd.sv().pos_id());
                        continue;
                    }
                }
                hmf_common::debug!(
                    "cluster({}) chromosome({}) marking {} breakends as suspect, state({})",
                    cluster.id(),
                    breakend.chromosome(),
                    proximate.len(),
                    state.to_str()
                );
                for p in &proximate {
                    if p.has_line_element(LineElementType::SUSPECT) {
                        continue;
                    }
                    p.sv().add_line_element(LineElementType::SUSPECT, p.uses_start());
                }
                has_suspected = true;
            }
        }
        self.check_is_line_cluster(cluster, has_suspected, state.has_insert_breakends());
    }

    fn check_is_line_cluster(&self, cluster: &Cluster, has_suspected: bool, has_insert_sites: bool) {
        let svs = cluster.svs();
        let known = |x: &Sv| x.has_line_element(LineElementType::KNOWN, true) || x.has_line_element(LineElementType::KNOWN, false);
        let all_non_sgl_in_known = !svs.iter().filter(|x| !x.is_sgl_breakend()).any(|x| !known(x));
        if all_non_sgl_in_known && cluster.type_count(BND) >= 1 {
            hmf_common::debug!("cluster({}) marked as line with all known({})", cluster.id(), svs.iter().filter(|x| known(x)).count());
            cluster.mark_as_line();
            return;
        }
        if has_suspected {
            if cluster.sv_count() <= 10 {
                cluster.mark_as_line();
                return;
            }
            let non_sgl_count = cluster.sv_count() as i32 - cluster.sgl_breakend_count();
            let non_sgl_line_count = svs.iter().filter(|x| !x.is_sgl_breakend()).filter(|x| x.in_line_element()).count() as i32;
            if non_sgl_line_count * 2 >= non_sgl_count {
                cluster.mark_as_line();
                return;
            }
        }
        if cluster.sv_count() <= 2 && cluster.sgl_breakend_count() == cluster.sv_count() as i32 {
            let has_sgl_line_repeat_class = svs.iter().filter(|x| x.sv_type() == SGL).any(|x| has_line_repeat_class(x));
            if cluster.type_count(SGL) == 2 && (has_sgl_line_repeat_class || has_insert_sites) {
                let var1 = &svs[0];
                let var2 = &svs[1];
                let db = var1.db_link(true);
                let in_short_db = db.as_ref().is_some_and(|d| Some(d) == var2.db_link(true).as_ref() && d.length() <= MIN_DEL_LENGTH);
                if in_short_db {
                    hmf_common::debug!("cluster({}) marked as line with SGL pair in DB and with insert poly A/T", cluster.id());
                    cluster.mark_as_line();
                    return;
                }
            }
            if cluster.type_count(SGL) == 1 && (has_sgl_line_repeat_class || has_insert_sites) {
                let var = svs.iter().find(|x| x.sv_type() == SGL).unwrap();
                if cluster.sv_count() == 1 {
                    let cn_low_side = var.be(true).copy_number_low_side();
                    let cn_change = var.copy_number_change(true);
                    if (cn_low_side <= 0.0 && cn_change > MAX_COPY_NUM_DIFF) || !copy_numbers_equal(cn_low_side, cn_low_side + cn_change) {
                        return;
                    }
                }
                hmf_common::debug!(
                    "cluster({}) marked as line with SGL with LINE repeat-class({}) insertPolyAT({})",
                    cluster.id(),
                    has_sgl_line_repeat_class,
                    has_insert_sites
                );
                cluster.mark_as_line();
            }
        }
    }
}

/// LineElementAnnotator.hasPolyAorTMotif
pub fn has_poly_a_or_t_motif(var: &SvVarData) -> bool {
    [SE_START, SE_END].iter().any(|se| var.breakend_se(*se).is_some_and(|b| has_line_source_motif(&b) || has_line_insert_motif(&b)))
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LineSuspectReason {
    NONE,
    BND_PAIR_POLY_AT,
    BND_PAIR_NO_DB,
    BND_SGL_REMOTE_DB_POLY_AT,
}

const SPANNING_LOCAL_MIN_DISTANCE: i32 = 1000000;
const POLY_AT_BASES_SEARCH: usize = 20;
const REPEAT_CLASS_LINE: &str = "LINE/L1";

/// LineClusterState
pub struct LineClusterState {
    proximity_distance: i32,
    breakends: Vec<Be>,
    source_breakends: Vec<Be>,
    insert_breakends: Vec<Be>,
    has_remote_db: bool,
    suspect_reason: LineSuspectReason,
}

impl LineClusterState {
    pub fn new(proximity_distance: i32) -> Self {
        LineClusterState {
            proximity_distance,
            breakends: Vec::new(),
            source_breakends: Vec::new(),
            insert_breakends: Vec::new(),
            has_remote_db: false,
            suspect_reason: LineSuspectReason::NONE,
        }
    }
    pub fn clear(&mut self) {
        self.breakends.clear();
        self.source_breakends.clear();
        self.insert_breakends.clear();
        self.suspect_reason = LineSuspectReason::NONE;
        self.has_remote_db = false;
    }
    pub fn is_suspected(&self) -> bool {
        self.suspect_reason != LineSuspectReason::NONE
    }
    pub fn add_breakend(&mut self, breakend: &Be) {
        if let Some(last) = self.breakends.last() {
            if breakend.position() - last.position() > self.proximity_distance {
                self.clear();
            }
        }
        self.breakends.push(breakend.clone());
        if has_line_source_motif(breakend) {
            self.source_breakends.push(breakend.clone());
        } else if has_line_insert_motif(breakend) {
            self.insert_breakends.push(breakend.clone());
        }
        self.suspect_reason = self.check_suspected_criteria();
    }
    pub fn has_insert_breakends(&self) -> bool {
        !self.insert_breakends.is_empty()
    }
    pub fn to_str(&self) -> String {
        format!(
            "{:?}: breakends({} src={} insert={}) remoteDB({})",
            self.suspect_reason,
            self.breakends.len(),
            self.source_breakends.len(),
            self.insert_breakends.len(),
            self.has_remote_db
        )
    }

    fn check_suspected_criteria(&mut self) -> LineSuspectReason {
        if self.suspect_reason != LineSuspectReason::NONE {
            return self.suspect_reason;
        }
        if self.source_breakends.len() >= 2 {
            return LineSuspectReason::BND_PAIR_POLY_AT;
        }
        let spanning: Vec<Be> = self
            .breakends
            .iter()
            .filter(|x| x.sv_type() == BND || (!x.sv().is_sgl_breakend() && x.sv().length() > SPANNING_LOCAL_MIN_DISTANCE))
            .filter(|x| !x.other_breakend().unwrap().has_line_element(LineElementType::KNOWN))
            .cloned()
            .collect();
        if spanning.len() >= 2
            && !self.source_breakends.is_empty()
            && spanning.iter().any(|x| x.db_link().is_none_or(|d| d.length() > MIN_DEL_LENGTH))
        {
            return LineSuspectReason::BND_PAIR_NO_DB;
        }
        for breakend in &spanning {
            let other = breakend.other_breakend().unwrap();
            let Some(db) = other.db_link() else { continue };
            let remote_other = db.other_breakend(&other);
            if remote_other.sv_type() != SGL || !self.is_remote_isolated_deletion_bridge(&db) {
                continue;
            }
            self.has_remote_db = true;
            if has_line_source_motif(breakend) {
                self.source_breakends.push(breakend.clone());
                return LineSuspectReason::BND_SGL_REMOTE_DB_POLY_AT;
            } else if has_line_insert_motif(&remote_other) {
                self.insert_breakends.push(remote_other.clone());
                return LineSuspectReason::BND_SGL_REMOTE_DB_POLY_AT;
            }
        }
        LineSuspectReason::NONE
    }

    fn is_remote_isolated_deletion_bridge(&self, db: &Db) -> bool {
        if db.length() > MIN_DEL_LENGTH {
            return false;
        }
        if db.lower().has_line_element(LineElementType::KNOWN) || db.upper().has_line_element(LineElementType::KNOWN) {
            return false;
        }
        let cluster = db.lower_sv().cl();
        let map = cluster.chr_breakend_map();
        let list = map.get(&db.chromosome().to_string()).expect("cluster breakends of chromosome");
        let lower = db.breakend(true);
        let start_index = lower.cluster_chr_pos_index();
        if start_index > 0 && lower.position() - list[(start_index - 1) as usize].position() < self.proximity_distance {
            return false;
        }
        let upper = db.breakend(false);
        let end_index = upper.cluster_chr_pos_index();
        if end_index < list.len() as i32 - 1 && list[(end_index + 1) as usize].position() - upper.position() < self.proximity_distance {
            return false;
        }
        true
    }
}

pub fn has_line_repeat_class(var: &SvVarData) -> bool {
    var.sv_data().insert_sequence_repeat_class == REPEAT_CLASS_LINE
}

pub fn has_line_source_motif(breakend: &SvBreakend) -> bool {
    let flip = flip_expected_sequence(breakend);
    let either_end = !breakend.sv().is_sgl_breakend();
    has_line_poly_a_or_t_motif(&breakend.sv().sv_data().insert_sequence, breakend.orientation(), !flip, either_end)
}

pub fn has_line_insert_motif(breakend: &SvBreakend) -> bool {
    let flip = flip_expected_sequence(breakend);
    let either_end = !breakend.sv().is_sgl_breakend();
    has_line_poly_a_or_t_motif(&breakend.sv().sv_data().insert_sequence, breakend.orientation(), flip, either_end)
}

fn flip_expected_sequence(breakend: &SvBreakend) -> bool {
    // Java: getOtherBreakend() of a start breakend is only dereferenced when it is not the start
    !breakend.uses_start() && breakend.other_breakend().is_some_and(|o| breakend.orientation() == o.orientation())
}

/// LineClusterState.hasLinePolyAorTMotif
pub fn has_line_poly_a_or_t_motif(ins: &str, orientation: i8, is_source: bool, either_end: bool) -> bool {
    let len = ins.len();
    if len < POLY_A_MOTIF.len() {
        return false;
    }
    let required = if (orientation == 1) == is_source { POLY_A_MOTIF } else { POLY_T_MOTIF };
    if (orientation == 1 || either_end) && contains_motif(&ins[..len.min(POLY_AT_BASES_SEARCH)], required) {
        return true;
    }
    if (orientation == -1 || either_end) && contains_motif(&ins[len.saturating_sub(POLY_AT_BASES_SEARCH)..], required) {
        return true;
    }
    false
}

pub fn contains_motif(sequence: &str, motif: &str) -> bool {
    if sequence.contains(motif) {
        return true;
    }
    let base = motif.as_bytes()[0];
    sequence.bytes().filter(|b| *b == base).count() >= 16
}

struct PseudoGeneMatch {
    gene: String,
    trans_id: i32,
    trans_name: String,
    exon_rank: i32,
    exon_length: i32,
    homology_offset: [i32; 2],
    position_mismatch: [i32; 2],
}

impl PseudoGeneMatch {
    fn is_homology_match(&self, is_start: bool) -> bool {
        self.position_mismatch[se(is_start)] == 0
    }
}

/// PseudoGeneFinder
pub struct PseudoGeneFinder {
    cache: Option<&'static GeneCache>,
}

impl PseudoGeneFinder {
    pub fn new(cache: Option<&'static GeneCache>) -> Self {
        PseudoGeneFinder { cache }
    }

    pub fn check_pseudo_gene_annotations(&self, clusters: &[Cluster], vis: &mut VisSampleData) {
        let Some(cache) = self.cache else { return };
        for cluster in clusters {
            // a HashMap keyed by LinkedPair: identity-hash order
            let mut pair_matches: IdHashMap<Pair, Vec<PseudoGeneMatch>> = IdHashMap::new();
            let mut matched_genes: Vec<BreakendGene> = Vec::new();
            for pair in cluster.linked_pairs() {
                if pair.base_length() > SHORT_TI_LENGTH * 8 {
                    continue;
                }
                let lower = pair.breakend(true);
                let upper = pair.breakend(false);
                let lower_genes = lower.sv().genes_list(lower.uses_start());
                let upper_genes = upper.sv().genes_list(upper.uses_start());
                if lower_genes.is_empty() || upper_genes.is_empty() {
                    continue;
                }
                let ld = lower.sv().sv_data();
                let ud = upper.sv().sv_data();
                let lower_homology = if lower.uses_start() { &ld.start_homology_sequence } else { &ld.end_homology_sequence };
                // Java takes the upper breakend's homology by the lower breakend's side
                let upper_homology = if lower.uses_start() { &ud.start_homology_sequence } else { &ud.end_homology_sequence };
                for gene1 in &lower_genes {
                    for gene2 in &upper_genes {
                        if gene1.gene_name() != gene2.gene_name() {
                            continue;
                        }
                        let matches = self.find_pseudo_gene_exon_matches(
                            cache,
                            gene1,
                            lower.position(),
                            upper.position(),
                            lower_homology.len() as i32,
                            upper_homology.len() as i32,
                        );
                        if !matches.is_empty() {
                            pair_matches.insert(pair.clone(), matches);
                            if !matched_genes.iter().any(|x| x.gene_name() == gene1.gene_name()) {
                                matched_genes.push(gene1.clone());
                            }
                        }
                    }
                }
            }
            if pair_matches.is_empty() {
                continue;
            }
            for gene in &matched_genes {
                let Some((max_trans_id, max_trans_name)) = find_most_common_transcript(gene.gene_name(), &pair_matches) else { continue };
                let matched_pairs: Vec<Pair> =
                    pair_matches.keys().filter(|x| is_relevant_match(x, &pair_matches, max_trans_id)).cloned().collect();
                if matched_pairs.is_empty() {
                    continue;
                }
                let mut gene_data = VisGeneData::new(
                    cluster.id(),
                    gene.gene_id(),
                    gene.gene_name(),
                    &max_trans_name,
                    max_trans_id,
                    &gene.chromosome(),
                    VisGeneAnnotationType::PSEUDOGENE,
                );
                for pair in &matched_pairs {
                    if let Some(m) = pair_matches.get(pair).unwrap().iter().find(|x| x.trans_id == max_trans_id) {
                        let exon_match_data = format!("{};{};{};{}", gene.gene_name(), max_trans_name, m.exon_rank, m.exon_length);
                        pair.set_exon_match_data(&exon_match_data);
                        gene_data.exon_position_offsets.insert(m.exon_rank, [0, 0]);
                    }
                }
                let Some(trans) = cache.transcript_data(&gene_data.gene_id, &gene_data.trans_name) else { continue };
                if trans.exons.is_empty() {
                    continue;
                }
                for exon in &trans.exons {
                    let has_pos_offsets = gene_data.exon_position_offsets.contains_key(&exon.rank);
                    let lost = gene_data.exons_lost_offsets.contains_key(&exon.rank);
                    if has_pos_offsets && !lost {
                        continue;
                    }
                    if !has_pos_offsets {
                        gene_data.exon_position_offsets.insert(exon.rank, [0, 0]);
                    }
                    if !lost {
                        gene_data.exons_lost_offsets.insert(exon.rank, [0, 0]);
                    }
                }
                vis.add_gene_exon_data(gene_data);
            }
        }
    }

    fn find_pseudo_gene_exon_matches(
        &self, cache: &GeneCache, gene: &BreakendGene, pos_start: i32, pos_end: i32, start_hom_len: i32, end_hom_len: i32,
    ) -> Vec<PseudoGeneMatch> {
        let mut out = Vec::new();
        let Some(trans_list) = cache.transcripts(gene.gene_id()) else { return out };
        for t in trans_list {
            for exon in &t.exons {
                let start_within = (exon.start - pos_start).abs() <= start_hom_len;
                let end_within = (exon.end - pos_end).abs() <= end_hom_len;
                if !start_within && !end_within {
                    continue;
                }
                if !start_within && pos_start < t.trans_start {
                    continue;
                }
                if !end_within && pos_end > t.trans_end {
                    continue;
                }
                let mut m = PseudoGeneMatch {
                    gene: gene.gene_name().to_string(),
                    trans_id: t.trans_id,
                    trans_name: t.trans_name.clone(),
                    exon_rank: exon.rank,
                    exon_length: exon.length(),
                    homology_offset: [0, 0],
                    position_mismatch: [0, 0],
                };
                if start_within {
                    m.homology_offset[SE_START] = pos_start - exon.start;
                } else {
                    m.position_mismatch[SE_START] = exon.start - pos_start;
                }
                if end_within {
                    m.homology_offset[SE_END] = pos_end - exon.end;
                } else {
                    m.position_mismatch[SE_END] = pos_end - exon.end;
                }
                out.push(m);
            }
        }
        out
    }

    /// variantMatchesPseudogeneExons. Java takes both gene lists from the start breakend (LX3).
    pub fn variant_matches_pseudogene_exons(&self, var: &SvVarData) -> bool {
        let Some(cache) = self.cache else { return false };
        if var.is_sgl_breakend() {
            return false;
        }
        let genes_start = var.genes_list(true);
        let genes_end = var.genes_list(crate::compat::java("LX3"));
        if genes_start.is_empty() || genes_end.is_empty() {
            return false;
        }
        let pos_start = var.position(true);
        let pos_end = var.position(false);
        let shl = var.sv_data().start_homology_sequence.len() as i32;
        let ehl = var.sv_data().end_homology_sequence.len() as i32;
        for gene in &genes_start {
            if !genes_end.iter().any(|x| x.gene_name() == gene.gene_name()) {
                continue;
            }
            let Some(trans_list) = cache.transcripts(gene.gene_id()) else { continue };
            for t in trans_list {
                if !t.exons.iter().any(|x| (x.start - pos_start).abs() <= shl || (x.end - pos_start).abs() <= shl) {
                    continue;
                }
                if t.exons.iter().any(|x| (x.start - pos_end).abs() <= ehl || (x.end - pos_end).abs() <= ehl) {
                    return true;
                }
            }
        }
        false
    }
}

/// PseudoGeneFinder.isPseudogeneDeletion
pub fn is_pseudogene_deletion(var: &SvVarData, del_start: i32, del_end: i32, trans: &hmf_common::ensembl::TranscriptData) -> bool {
    let shl = var.sv_data().start_homology_sequence.len() as i32;
    let ehl = var.sv_data().end_homology_sequence.len() as i32;
    for i in 0..trans.exons.len().saturating_sub(1) {
        let exon = &trans.exons[i];
        let next = &trans.exons[i + 1];
        if (exon.end - del_start).abs() <= shl && (next.start - del_end).abs() <= ehl {
            return true;
        }
    }
    false
}

fn find_pseudo_match<'a>(breakend: &Be, map: &'a IdHashMap<Pair, Vec<PseudoGeneMatch>>, trans_id: i32) -> Option<&'a PseudoGeneMatch> {
    for (pair, matches) in map.iter() {
        if pair.has_breakend(breakend) {
            return matches.iter().find(|x| x.trans_id == trans_id);
        }
    }
    None
}

fn is_relevant_match(pair: &Pair, map: &IdHashMap<Pair, Vec<PseudoGeneMatch>>, trans_id: i32) -> bool {
    let Some(m) = map.get(pair).unwrap().iter().find(|x| x.trans_id == trans_id) else { return false };
    if m.is_homology_match(true) && m.is_homology_match(false) {
        return true;
    }
    for se in [SE_START, SE_END] {
        let breakend = pair.breakend_se(se);
        let Some(other) = breakend.other_breakend() else { continue };
        if find_pseudo_match(&other, map, trans_id).is_some() {
            return true;
        }
    }
    false
}

/// findMostCommonTranscript: (trans id, name) of the first transcript to reach the highest count
fn find_most_common_transcript(gene_name: &str, map: &IdHashMap<Pair, Vec<PseudoGeneMatch>>) -> Option<(i32, String)> {
    let mut counts: JavaHashMap<i32, i32> = JavaHashMap::new();
    let mut max: Option<(i32, String)> = None;
    let mut max_count = 0;
    for (_, matches) in map.iter() {
        for m in matches {
            if m.gene != gene_name {
                continue;
            }
            let c = counts.get(&m.trans_id).copied().unwrap_or(0) + 1;
            counts.insert(m.trans_id, c);
            if c > max_count {
                max_count = c;
                max = Some((m.trans_id, m.trans_name.clone()));
            }
        }
    }
    max
}
