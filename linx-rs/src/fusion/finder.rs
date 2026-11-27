//! FusionFinder (fusions between the genes of two breakends) and SpecialFusions (promiscuous enhancer
//! targets, and the known fusion genes LINX adds to the gene cache).

use std::cell::Cell;

use hmf_common::ensembl::{GeneData, TranscriptData};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::known_fusion::KnownFusionType::*;
use hmf_common::known_fusion::{KnownFusionCache, KnownFusionData, FS_DOWN, FS_UP, NEG_STRAND, POS_STRAND};
use hmf_common::sv::StructuralVariantType::BND;
use hmf_common::transcript_utils::{tick_phase_forward, PHASE_NONE};

use super::*;
use crate::gene::{BreakendGene, BreakendTrans, GeneCache, TranscriptCodingType, TranscriptRegionType};
use crate::types::*;

pub struct FusionFinder {
    pub config: FusionConfig,
    pub gene_cache: &'static GeneCache,
    pub known_fusion_cache: &'static KnownFusionCache,
    next_fusion_id: Cell<i32>,
}

impl FusionFinder {
    pub fn new(config: FusionConfig, gene_cache: &'static GeneCache, known_fusion_cache: &'static KnownFusionCache) -> Self {
        FusionFinder { config, gene_cache, known_fusion_cache, next_fusion_id: Cell::new(0) }
    }
    pub fn has_valid_config_data(&self) -> bool {
        self.known_fusion_cache.has_valid_data()
    }
    pub fn reset(&self) {
        self.next_fusion_id.set(0)
    }
    pub fn next_fusion_id(&self) -> i32 {
        let id = self.next_fusion_id.get();
        self.next_fusion_id.set(id + 1);
        id
    }

    /// findFusions: every valid transcript pair between the two breakends' genes
    pub fn find_fusions(&self, genes1: &[BreakendGene], genes2: &[BreakendGene]) -> Vec<GeneFusion> {
        let mut potential: Vec<GeneFusion> = Vec::new();
        if !self.has_valid_config_data() {
            return potential;
        }
        let kfc = self.known_fusion_cache;
        for start_gene in genes1 {
            if !self.config.process_sv_gene_debug(start_gene.gene_name(), start_gene.var_id()) {
                continue;
            }
            let start_upstream = start_gene.is_upstream();
            let start_enhancer = kfc.within_enhancer_region(&start_gene.chromosome(), start_gene.position());
            for end_gene in genes2 {
                if !self.config.process_sv_gene_debug(end_gene.gene_name(), end_gene.var_id()) {
                    continue;
                }
                let end_upstream = end_gene.is_upstream();
                let end_enhancer = kfc.within_enhancer_region(&end_gene.chromosome(), end_gene.position());
                if start_enhancer || end_enhancer {
                    if start_enhancer && end_enhancer {
                        continue;
                    }
                    self.check_enhancer_fusion(start_gene, end_gene, &mut potential);
                    continue;
                }
                if start_upstream == end_upstream {
                    continue;
                }
                let up_gene = if start_upstream { start_gene } else { end_gene };
                let down_gene = if !start_upstream { start_gene } else { end_gene };
                let known_pair = kfc.has_known_fusion(up_gene.gene_name(), down_gene.gene_name());
                let known_unmappable3 = kfc.has_known_unmappable3_fusion(up_gene.gene_name(), down_gene.gene_name());
                for up_trans in up_gene.transcripts() {
                    if !is_valid_upstream_transcript(&up_trans, !known_pair, self.config.require_upstream_biotypes) {
                        continue;
                    }
                    for down_trans in down_gene.transcripts() {
                        let fusion = if known_unmappable3 {
                            let f = GeneFusion::create(&up_trans, &down_trans, true);
                            f.set_id(self.next_fusion_id());
                            Some(f)
                        } else {
                            if !is_valid_downstream_transcript(&down_trans) {
                                continue;
                            }
                            let exon_del_dup_candidate =
                                up_trans.trans_name() == down_trans.trans_name() && kfc.is_exon_del_dup_trans(up_trans.trans_name());
                            transcript_pair_has_fusion(&up_trans, &down_trans, &self.config, known_pair, exon_del_dup_candidate)
                        };
                        let Some(fusion) = fusion else { continue };
                        fusion.set_id(self.next_fusion_id());
                        self.set_known_fusion_type(&fusion);
                        potential.push(fusion);
                    }
                }
            }
        }
        potential
    }

    fn check_enhancer_fusion(&self, start_gene: &BreakendGene, end_gene: &BreakendGene, potential: &mut Vec<GeneFusion>) {
        let kfc = self.known_fusion_cache;
        let start_is = kfc.within_enhancer_region(&start_gene.chromosome(), start_gene.position());
        let end_is = kfc.within_enhancer_region(&end_gene.chromosome(), end_gene.position());
        if start_is == end_is {
            return;
        }
        let enhancer_gene = if start_is { start_gene } else { end_gene };
        let down_gene = if start_is { end_gene } else { start_gene };
        let ech = enhancer_gene.chromosome();
        let epos = enhancer_gene.position();
        let mut kf_data: Option<&KnownFusionData> =
            kfc.data_by_type(ENHANCER_KNOWN_PAIR).find(|x| x.three_gene == down_gene.gene_name() && x.within_gene_region(&ech, epos));
        if kf_data.is_none() {
            kf_data = kfc.data_by_type(ENHANCER_PROMISCUOUS).find(|x| x.within_gene_region(&ech, epos));
            if kf_data.is_none() {
                return;
            }
            if ech == down_gene.chromosome() {
                let distance = (epos - down_gene.position()).abs();
                if distance < MIN_ENHANCER_PARTNER_GENE_DISTANCE {
                    return;
                }
            }
        }
        let kf_data = kf_data.unwrap();
        let known_type = kf_data.fusion_type;
        let mut down_transcripts: Vec<BreakendTrans> = Vec::new();
        let mut max_downstream_distance = -1;
        let mut allow_coding = false;
        if kf_data.gene_strand().is_none() {
            max_downstream_distance = MAX_ENHANCER_PARTNER_GENE_DISTANCE_NO_ORIENT;
            allow_coding = true;
        } else if kf_data.downstream_distance(FS_DOWN) > 0 {
            max_downstream_distance = kf_data.downstream_distance(FS_DOWN);
        }
        for t in down_gene.transcripts() {
            if t.coding_bases == 0 && !t.is_upstream() {
                if t.distance_upstream() <= MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN {
                    down_transcripts.push(t);
                }
            } else if t.coding_type() == TranscriptCodingType::CODING && allow_coding {
                down_transcripts.push(t);
            } else if t.coding_type() == TranscriptCodingType::UTR_3P && t.is_upstream() && max_downstream_distance > 0 {
                if t.distance_downstream() <= max_downstream_distance {
                    down_transcripts.push(t);
                }
            }
        }
        let up_trans = generate_enhancer_transcript(enhancer_gene, kf_data);
        for down_trans in down_transcripts {
            if has_matching_enhancer_fusion(potential, known_type, &up_trans, &down_trans) {
                continue;
            }
            let fusion = GeneFusion::create(&up_trans, &down_trans, true);
            fusion.set_id(self.next_fusion_id());
            fusion.set_known_type(known_type);
            potential.push(fusion);
        }
    }

    /// findTopReportableFusion
    pub fn find_top_reportable_fusion(&self, fusions: &[GeneFusion]) -> Option<GeneFusion> {
        let mut candidates: Vec<GeneFusion> = Vec::new();
        for fusion in fusions {
            let reasons = determine_reportability(fusion);
            fusion.set_reportable_reasons(&reasons);
            if !reasons.is_empty() {
                continue;
            }
            candidates.push(fusion.clone());
        }
        if candidates.is_empty() {
            return find_top_priority_fusion(fusions);
        }
        let reportable = find_top_priority_fusion(&candidates);
        if let Some(r) = &reportable {
            if check_protein_domains(r.known_type()) {
                self.set_fusion_protein_features(r);
                if valid_protein_domains(r) {
                    r.set_reportable(true);
                } else {
                    r.add_reportable_reason(FusionReportableReason::PROTEIN_DOMAINS);
                }
            } else {
                r.set_reportable(true);
            }
        }
        reportable
    }

    /// setFusionProteinFeatures
    pub fn set_fusion_protein_features(&self, fusion: &GeneFusion) {
        if fusion.protein_features_set() {
            return;
        }
        fusion.set_protein_features_set();
        let down = fusion.downstream_trans();
        if down.non_coding() {
            return;
        }
        if (fusion.known_type() == ENHANCER_KNOWN_PAIR || fusion.known_type() == ENHANCER_PROMISCUOUS) && down.post_coding() {
            return;
        }
        let Some(protein_data) = self.gene_cache.ensembl.transcript_protein_data_map().get(&down.trans_id()) else { return };
        if protein_data.is_empty() {
            return;
        }
        let mut processed: Vec<&str> = Vec::new();
        for i in 0..protein_data.len() {
            let pf = &protein_data[i];
            let feature = pf.hit_description.as_str();
            if processed.contains(&feature) {
                continue;
            }
            let feature_start = pf.seq_start;
            let mut feature_end = pf.seq_end;
            for other in &protein_data[i + 1..] {
                if other.hit_description == feature {
                    feature_end = other.seq_end.max(feature_end);
                }
            }
            let mut preserved = protein_feature_preserved(down, true, feature_start, feature_end);
            if !preserved && down.gene().gene_id() == fusion.upstream_trans().gene().gene_id() {
                preserved = protein_feature_preserved(fusion.upstream_trans(), false, feature_start, feature_end);
            }
            down.add_protein_feature(feature, preserved);
            processed.push(feature);
        }
    }

    fn set_known_fusion_type(&self, fusion: &GeneFusion) {
        let kfc = self.known_fusion_cache;
        let up_gene = fusion.transcripts()[FS_UP].gene_name();
        let down_gene = fusion.transcripts()[FS_DOWN].gene_name();
        if kfc.has_known_fusion(up_gene, down_gene) {
            fusion.set_known_type(KNOWN_PAIR);
            return;
        }
        if up_gene == down_gene {
            if kfc.within_known_exon_ranges(
                EXON_DEL_DUP,
                fusion.transcripts()[FS_UP].trans_name(),
                fusion.breakend_exon(true),
                fusion.fused_exon(true),
                fusion.breakend_exon(false),
                fusion.fused_exon(false),
            ) {
                fusion.set_known_type(EXON_DEL_DUP);
                fusion.set_known_exons();
                return;
            }
            fusion.set_known_type(NONE);
            return;
        }
        let mut requires_exon_range = false;
        if kfc.has_promiscuous_three_gene(down_gene) {
            fusion.set_known_type(PROMISCUOUS_3);
            fusion.set_promiscuous(FS_DOWN);
            let trans_name = fusion.transcripts()[FS_DOWN].trans_name();
            requires_exon_range = kfc.requires_promiscuous_exon_range(PROMISCUOUS_3, trans_name);
            if kfc.within_promiscuous_exon_range(PROMISCUOUS_3, trans_name, fusion.breakend_exon(false), fusion.fused_exon(false)) {
                fusion.set_known_exons();
            }
        }
        if kfc.has_promiscuous_five_gene(up_gene) {
            fusion.set_known_type(PROMISCUOUS_5);
            fusion.set_promiscuous(FS_UP);
            let trans_name = fusion.transcripts()[FS_UP].trans_name();
            requires_exon_range = kfc.requires_promiscuous_exon_range(PROMISCUOUS_5, trans_name);
            if kfc.within_promiscuous_exon_range(PROMISCUOUS_5, trans_name, fusion.breakend_exon(true), fusion.fused_exon(true)) {
                fusion.set_known_exons();
            }
        }
        let kt = fusion.known_type();
        if (kt == PROMISCUOUS_5 || kt == PROMISCUOUS_3)
            && kfc.is_high_impact_promiscuous(kt, &fusion.gene_name(FS_UP), &fusion.gene_name(FS_DOWN))
            && (!requires_exon_range || fusion.known_exons())
        {
            fusion.set_high_impact_promiscuous();
        }
    }
}

/// validFusionTranscript
pub fn valid_fusion_transcript(t: &BreakendTrans) -> bool {
    if t.is_upstream() { is_valid_upstream_transcript(t, true, false) } else { is_valid_downstream_transcript(t) }
}

fn is_valid_upstream_transcript(t: &BreakendTrans, require_disruptive: bool, require_biotypes: bool) -> bool {
    if t.is_promoter() {
        return false;
    }
    if require_disruptive && !t.is_disruptive() {
        return false;
    }
    if require_biotypes && !REQUIRED_BIOTYPES.contains(&t.bio_type()) {
        return false;
    }
    true
}

fn is_valid_downstream_transcript(t: &BreakendTrans) -> bool {
    !(t.post_coding() || t.non_coding() || t.exon_count() == 1)
}

/// checkFusionLogic
pub fn check_fusion_logic(up: &BreakendTrans, down: &BreakendTrans, config: &FusionConfig) -> Option<GeneFusion> {
    if !is_valid_upstream_transcript(up, true, config.require_upstream_biotypes) || !is_valid_downstream_transcript(down) {
        return None;
    }
    transcript_pair_has_fusion(up, down, config, false, false)
}

fn transcript_pair_has_fusion(
    up: &BreakendTrans, down: &BreakendTrans, config: &FusionConfig, is_known_pair: bool, exon_del_dup_candidate: bool,
) -> Option<GeneFusion> {
    let mut check_exact_match = false;
    if up.pre_coding() {
        if up.is_exonic() && !down.is_exonic() {
            return None;
        } else if down.is_coding() && !is_known_pair {
            return None;
        } else if down.pre_coding() && up.gene().gene_id() == down.gene().gene_id() {
            return None;
        }
    } else if up.is_coding() {
        if down.non_coding() {
            return None;
        }
        if down.pre_coding() && !is_known_pair {
            return None;
        }
        if up.is_exonic() {
            if !down.is_exonic() && (!config.allow_exon_skipping || !down.is_intronic()) {
                return None;
            }
            if up.gene().var_id() != down.gene().var_id() {
                return None;
            }
            check_exact_match = down.is_exonic();
        }
    } else if up.non_coding() {
        if up.is_exonic() && !down.is_exonic() && !is_known_pair {
            return None;
        } else if down.is_coding() && !is_known_pair {
            return None;
        }
    }
    if is_irrelevant_same_gene(up, down) {
        return None;
    }
    let mut phase_matched;
    let mut skipped_up = 0;
    let mut skipped_down = 0;
    if !check_exact_match && down.is_exonic() && down.exon_downstream == down.exon_count() && !down.pre_coding() {
        return None;
    }
    if check_exact_match {
        phase_matched = exon_to_exon_in_phase(up, down);
        if phase_matched || !config.require_phase_match {
            return Some(GeneFusion::create(up, down, phase_matched));
        }
    } else {
        if up.is_exonic() && down.is_exonic() {
            phase_matched = ((up.pre_coding() || up.non_coding()) && down.pre_coding()) || (up.post_coding() && down.post_coding());
        } else {
            phase_matched = up.phase == down.phase;
        }
        if !phase_matched && config.allow_exon_skipping && (up.gene().gene_id() != down.gene().gene_id() || exon_del_dup_candidate) {
            let up_alt: Vec<(i32, i32)> = up.alternative_phasing().iter().map(|(k, v)| (*k, *v)).collect();
            let down_alt: Vec<(i32, i32)> = down.alternative_phasing().iter().map(|(k, v)| (*k, *v)).collect();
            if let Some((_, v)) = up_alt.iter().find(|(k, _)| *k == down.phase) {
                phase_matched = true;
                skipped_up = *v;
            }
            if !phase_matched {
                if let Some((_, v)) = down_alt.iter().find(|(k, _)| *k == up.phase) {
                    phase_matched = true;
                    skipped_down = *v;
                }
                if !phase_matched {
                    // the inner break leaves the outer loop running: the last upstream phasing with a match wins
                    for (ku, vu) in &up_alt {
                        if let Some((_, vd)) = down_alt.iter().find(|(kd, _)| kd == ku) {
                            phase_matched = true;
                            skipped_up = *vu;
                            skipped_down = *vd;
                        }
                    }
                }
            }
        }
        if up.is_exonic() && skipped_up == 1 {
            skipped_up = 0;
        }
        if down.is_exonic() && skipped_down == 1 {
            skipped_down = 0;
        }
        if phase_matched || !config.require_phase_match {
            let fusion = GeneFusion::create(up, down, phase_matched);
            fusion.set_exons_skipped(skipped_up, skipped_down);
            return Some(fusion);
        }
    }
    None
}

fn exon_to_exon_in_phase(up: &BreakendTrans, down: &BreakendTrans) -> bool {
    let up_phase = up.exonic_base_phase;
    let down_phase = down.exonic_base_phase;
    if up_phase == PHASE_NONE && down_phase == PHASE_NONE {
        return true;
    }
    tick_phase_forward(up_phase, 1) == down_phase
}

fn has_matching_enhancer_fusion(
    potential: &[GeneFusion], known_type: hmf_common::known_fusion::KnownFusionType, up: &BreakendTrans, down: &BreakendTrans,
) -> bool {
    for fusion in potential {
        if fusion.known_type() != known_type {
            continue;
        }
        if fusion.sv_id(true) != up.gene().var_id() || fusion.sv_id(false) != down.gene().var_id() {
            continue;
        }
        if fusion.upstream_trans().trans_name() != up.trans_name() {
            continue;
        }
        if fusion.downstream_trans() == down {
            return true;
        }
    }
    false
}

fn leak_transcript(t: TranscriptData) -> crate::gene::Trans {
    Box::leak(Box::new(t))
}
fn leak_gene(g: GeneData) -> crate::gene::Gene {
    Box::leak(Box::new(g))
}

#[allow(clippy::too_many_arguments)]
fn new_transcript(
    trans_name: &str, gene_id: &str, strand: i8, trans_start: i32, trans_end: i32,
) -> crate::gene::Trans {
    leak_transcript(TranscriptData {
        trans_id: 0,
        trans_name: trans_name.to_string(),
        gene_id: gene_id.to_string(),
        is_canonical: false,
        strand,
        trans_start,
        trans_end,
        coding_start: None,
        coding_end: None,
        bio_type: String::new(),
        ref_seq_id: None,
        exons: Vec::new(),
    })
}

fn new_gene(gene_id: &str, gene_name: &str, chromosome: &str, strand: i8, start: i32, end: i32) -> crate::gene::Gene {
    leak_gene(GeneData {
        gene_id: gene_id.to_string(),
        gene_name: gene_name.to_string(),
        chromosome: chromosome.to_string(),
        strand,
        gene_start: start,
        gene_end: end,
        karyotype_band: String::new(),
        synonyms: String::new(),
    })
}

fn generate_enhancer_transcript(gene: &BreakendGene, kf: &KnownFusionData) -> BreakendTrans {
    let region = kf.gene_region().expect("enhancer gene region");
    let trans = new_transcript(&format!("@{}", kf.five_gene), gene.gene_id(), gene.strand(), region.start, region.end);
    let t = BreakendTrans::create(gene, trans, -1, -1, PHASE_NONE, PHASE_NONE, 0, 0);
    t.set_coding_type(TranscriptCodingType::ENHANCER);
    t.set_region_type(TranscriptRegionType::IG);
    t.set_is_disruptive(true);
    t
}

/// isIrrelevantSameGene
pub fn is_irrelevant_same_gene(up: &BreakendTrans, down: &BreakendTrans) -> bool {
    if up.gene_name() != down.gene_name() {
        return false;
    }
    if up.trans_name() != down.trans_name() {
        return true;
    }
    if up.non_coding() {
        return true;
    }
    up.is_intronic() && down.is_intronic() && up.exon_upstream == down.exon_upstream
}

fn protein_feature_preserved(t: &BreakendTrans, is_downstream: bool, feature_start: i32, feature_end: i32) -> bool {
    if t.pre_coding() {
        is_downstream
    } else if t.post_coding() {
        !is_downstream
    } else if is_downstream {
        feature_start * 3 - 2 >= t.coding_bases
    } else {
        feature_end * 3 <= t.coding_bases
    }
}

// SpecialFusions

/// getSpecificPreGeneDistances: the ENHANCER_KNOWN_PAIR 3' genes
pub fn specific_pre_gene_distances(kfc: &KnownFusionCache) -> JavaHashMap<String, i32> {
    let mut m = JavaHashMap::new();
    for x in kfc.data_by_type(ENHANCER_KNOWN_PAIR) {
        m.insert(x.three_gene.clone(), MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN);
    }
    m
}

/// cacheSpecialFusionGenes: downstream gene distances and alternative 3' gene regions
pub fn cache_special_fusion_genes(cache: &GeneCache, kfc: &KnownFusionCache) {
    for kf in kfc.data() {
        if kf.downstream_distance(FS_UP) > 0 {
            if let Some(g) = cache.gene_data_by_name(&kf.five_gene) {
                cache.add_downstream_gene_annotations(g, kf.downstream_distance(FS_UP));
            }
        }
        if kf.downstream_distance(FS_DOWN) > 0 {
            if let Some(g) = cache.gene_data_by_name(&kf.three_gene) {
                cache.add_downstream_gene_annotations(g, kf.downstream_distance(FS_DOWN));
            }
        }
        if !kf.three_gene_alt_regions().is_empty() {
            if let Some(g) = cache.gene_data_by_name(&kf.three_gene) {
                if cache.alternative_gene_data.borrow().iter().any(|x| x.gene_id == g.gene_id) {
                    continue;
                }
                for r in kf.three_gene_alt_regions() {
                    let alt = new_gene(&g.gene_id, &g.gene_name, &r.chromosome, g.strand, r.start, r.end);
                    cache.alternative_gene_data.borrow_mut().push(alt);
                }
            }
        }
    }
}

/// SpecialFusions.findFusions: promiscuous enhancer targets
pub fn find_special_fusions(finder: &FusionFinder, sv_list: &[Sv]) -> Vec<GeneFusion> {
    let mut fusions = Vec::new();
    let targets: Vec<&KnownFusionData> = finder.known_fusion_cache.data_by_type(PROMISCUOUS_ENHANCER_TARGET).collect();
    if targets.is_empty() {
        return fusions;
    }
    for var in sv_list {
        if var.is_sgl_breakend() {
            continue;
        }
        for s in [SE_START, SE_END] {
            let breakend = var.be(is_start(s));
            for kf in &targets {
                if let Some(f) = form_enhancer_target_fusion(finder, kf, var, &breakend, s) {
                    fusions.push(f);
                }
            }
        }
    }
    fusions
}

fn form_enhancer_target_fusion(finder: &FusionFinder, kf: &KnownFusionData, var: &Sv, breakend: &Be, se_index: usize) -> Option<GeneFusion> {
    if let Some(o) = kf.gene_strand() {
        if breakend.orientation() != o {
            return None;
        }
    }
    if !kf.within_gene_region(breakend.chromosome(), breakend.position()) {
        return None;
    }
    if var.sv_type() != BND && var.length() < ENHANCER_PROMISCUOUS_MIN_DISTANCE {
        return None;
    }
    let other = var.be(!is_start(se_index));
    let other_is_start = se_index == SE_END;
    let mut up_trans: Option<BreakendTrans> = var.genes_list(other_is_start).iter().find_map(|x| x.canonical());
    if up_trans.is_none() {
        let g = new_gene("", "", other.chromosome(), POS_STRAND, other.position(), other.position());
        let gene = BreakendGene::create(var.id(), other_is_start, g);
        gene.set_sv_data(var.sv_data(), var.jcn());
        let trans = new_transcript("", "", POS_STRAND, 0, 0);
        let t = BreakendTrans::create(&gene, trans, 0, 0, PHASE_NONE, PHASE_NONE, 0, 0);
        t.set_coding_type(TranscriptCodingType::ENHANCER);
        t.set_region_type(TranscriptRegionType::UNKNOWN);
        up_trans = Some(t);
    }
    let down_trans: BreakendTrans;
    let downstream_gene = var.genes_list(se_index == SE_START).into_iter().find(|x| x.gene_name() == kf.three_gene);
    if let Some(dg) = downstream_gene {
        // Java goes on with a null transcript where the gene has no canonical one
        down_trans = dg.canonical()?;
    } else {
        let region = kf.gene_region().expect("known fusion gene region");
        let (g, trans) = match finder.gene_cache.gene_data_by_name(&kf.three_gene) {
            Some(g) => (g, finder.gene_cache.ensembl.canonical_transcript_data(&g.gene_id)?),
            None => {
                let strand = if kf.gene_strand().is_none_or(|o| o == 1) { NEG_STRAND } else { POS_STRAND };
                let g = new_gene("", &kf.three_gene, &region.chromosome, strand, region.start, region.end);
                (g, new_transcript("", &kf.three_gene, strand, region.start, region.end))
            }
        };
        let gene = BreakendGene::create(var.id(), se_index == SE_START, g);
        gene.set_sv_data(var.sv_data(), var.jcn());
        let t = BreakendTrans::create(&gene, trans, 1, 1, PHASE_NONE, PHASE_NONE, 0, 0);
        t.set_coding_type(TranscriptCodingType::UTR_5P);
        t.set_region_type(TranscriptRegionType::UPSTREAM);
        down_trans = t;
    }
    let fusion = GeneFusion::create(&up_trans.unwrap(), &down_trans, false);
    fusion.set_id(finder.next_fusion_id());
    fusion.set_known_type(PROMISCUOUS_ENHANCER_TARGET);
    let cluster = var.cl();
    fusion.set_annotations(FusionAnnotations {
        cluster_id: cluster.id(),
        cluster_count: cluster.sv_count() as i32,
        resolved_type: cluster.resolved_type().to_string(),
        chain_info: None,
        terminated_up: false,
        terminated_down: false,
    });
    Some(fusion)
}
