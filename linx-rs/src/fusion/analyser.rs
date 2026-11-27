//! FusionDisruptionAnalyser: marks disrupted transcripts, finds fusions from single SVs and along chains,
//! the reportable disruptions, and writes the breakend and fusion files.

use anyhow::Result;
use hmf_common::driver_panel::DriverGene;
use hmf_common::java_hash::JavaHashMap;
use hmf_common::known_fusion::KnownFusionType::*;
use hmf_common::known_fusion::{FS_DOWN, FS_UP};

use super::disruption::{is_disruptive_in_transcript, undisrupted_copy_number, DisruptionFinder};
use super::finder::{find_special_fusions, valid_fusion_transcript, FusionFinder};
use super::neo::NeoEpitopeWriter;
use super::*;
use crate::gene::{BreakendGene, BreakendTrans};
use crate::output::{
    breakend_lines, fusion_lines, LinxBreakend, LinxFusion, BREAKEND_ORIENTATION_DOWNSTREAM, BREAKEND_ORIENTATION_UPSTREAM,
};
use crate::types::*;
use crate::vis::{VisFusion, VisGeneAnnotationType, VisGeneData, VisSampleData};

struct ValidTraversalData {
    link: Pair,
    is_valid: bool,
    fusion_direction: i32,
    is_precoding_upstream: bool,
}

pub struct FusionDisruptionAnalyser {
    pub finder: FusionFinder,
    pub disruption_finder: DisruptionFinder,
    neo: Option<NeoEpitopeWriter>,
    sample_id: String,
    output_dir: String,
    run_fusions: bool,
    is_germline: bool,
    write_vis: bool,
    driver_genes: &'static JavaHashMap<String, DriverGene>,
    /// getFusions: every candidate fusion
    pub fusions: Vec<GeneFusion>,
    pub unique_fusions: Vec<GeneFusion>,
    /// getInvalidFusions (a HashMap keyed by fusion in Java, for the RNA matching): one per fusion name
    pub invalid_fusions: Vec<(GeneFusion, String)>,
}

impl FusionDisruptionAnalyser {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        finder: FusionFinder, disruption_finder: DisruptionFinder, output_dir: &str, run_fusions: bool, is_germline: bool, write_vis: bool,
        write_neo_epitopes: bool, driver_genes: &'static JavaHashMap<String, DriverGene>,
    ) -> Self {
        let neo = if write_neo_epitopes {
            Some(NeoEpitopeWriter::new(output_dir, finder.gene_cache, finder.known_fusion_cache))
        } else {
            None
        };
        FusionDisruptionAnalyser {
            finder,
            disruption_finder,
            neo,
            sample_id: String::new(),
            output_dir: output_dir.to_string(),
            run_fusions,
            is_germline,
            write_vis,
            driver_genes,
            fusions: Vec::new(),
            unique_fusions: Vec::new(),
            invalid_fusions: Vec::new(),
        }
    }

    /// the constructor's cacheSpecialFusionGenes (when running fusions): the known fusion genes' downstream
    /// distances and alternative regions, added to the gene cache
    pub fn cache_special_fusion_genes(&self) {
        super::finder::cache_special_fusion_genes(self.finder.gene_cache, self.finder.known_fusion_cache);
    }

    /// annotateTranscripts: marks transcripts disruptive, then drops those of no further use
    pub fn annotate_transcripts(&self, sv_list: &[Sv], purge_invalid: bool) {
        self.disruption_finder.mark_transcripts_disruptive(sv_list);
        if !purge_invalid {
            return;
        }
        for var in sv_list {
            for be in [SE_START, SE_END] {
                if be == SE_END && var.is_sgl_breakend() {
                    continue;
                }
                for gene in var.genes_list(is_start(be)) {
                    gene.retain_transcripts(|t| t.is_disruptive() || valid_fusion_transcript(t) || t.is_canonical());
                }
            }
        }
    }

    pub fn run(&mut self, sample_id: &str, full_sv_list: &[Sv], clusters: &[Cluster], vis: &mut VisSampleData) -> Result<()> {
        self.sample_id = sample_id.to_string();
        self.unique_fusions.clear();
        self.finder.reset();
        let sv_list: Vec<Sv> = full_sv_list.iter().filter(|x| x.linked_svs().is_none() || !x.is_sgl_breakend()).cloned().collect();
        if self.is_germline {
            self.disruption_finder.find_reportable_disruptions(&sv_list, clusters, Some(vis));
            if self.output_dir.is_empty() {
                return Ok(());
            }
            return self.disruption_finder.write_germline_disruptions(sample_id, &self.output_dir);
        }
        if self.run_fusions && self.finder.has_valid_config_data() {
            self.find_fusions(&sv_list, clusters)?;
        }
        self.disruption_finder.find_reportable_disruptions(&sv_list, clusters, Some(vis));
        if self.run_fusions {
            let unique = self.extract_unique_fusions();
            self.unique_fusions.extend(unique);
            for f in &self.unique_fusions {
                if f.known_type() != NONE {
                    self.finder.set_fusion_protein_features(f);
                }
            }
        }
        let transcripts = self.build_transcript_list(&sv_list);
        let (fusions, breakends) = convert_breakends_and_fusions(&self.unique_fusions, self.driver_genes, &transcripts);
        // no output directory (Java's null): nothing written
        if !self.output_dir.is_empty() {
            hmf_common::purple_files::write_lines(&format!("{}/{}.linx.breakend.tsv", self.output_dir, sample_id), &breakend_lines(&breakends))?;
            hmf_common::purple_files::write_lines(&format!("{}/{}.linx.fusion.tsv", self.output_dir, sample_id), &fusion_lines(&fusions))?;
        }
        if self.run_fusions {
            self.add_visualisation_data(vis);
        }
        Ok(())
    }

    pub fn close(&mut self) -> Result<()> {
        if let Some(n) = self.neo.as_mut() {
            n.close()?;
        }
        Ok(())
    }

    fn find_fusions(&mut self, sv_list: &[Sv], clusters: &[Cluster]) -> Result<()> {
        if self.sample_id.is_empty() {
            return Ok(());
        }
        self.fusions.clear();
        self.invalid_fusions.clear();
        if let Some(n) = self.neo.as_mut() {
            n.initialise_sample(&self.sample_id)?;
        }
        self.find_single_sv_fusions(sv_list)?;
        self.find_chained_fusions(clusters)?;
        let special = find_special_fusions(&self.finder, sv_list);
        self.fusions.extend(special);
        Ok(())
    }

    fn breakend_gene_list(&self, var: &Sv, is_start: bool) -> Vec<BreakendGene> {
        if var.is_sgl_breakend() && !is_start {
            let genes = var.genes_list(false);
            let kfc = self.finder.known_fusion_cache;
            if var.sgl_mappings().iter().any(|m| kfc.within_enhancer_region(&m.chromosome, m.position)) {
                return genes;
            }
            return genes.into_iter().filter(|x| kfc.is_single_breakend_candidate(x.gene_name(), x.is_upstream())).collect();
        }
        var.genes_list(is_start)
    }

    fn find_single_sv_fusions(&mut self, sv_list: &[Sv]) -> Result<()> {
        let kfc = self.finder.known_fusion_cache;
        for var in sv_list {
            if var.is_sgl_breakend() && var.sgl_mappings().is_empty() {
                continue;
            }
            let cluster = var.cl();
            if cluster.sv_count() > 1 && cluster.find_chain(var).is_some() {
                let ig_candidate = !var.is_sgl_breakend()
                    && (kfc.within_enhancer_region(var.chromosome(true), var.position(true))
                        != kfc.within_enhancer_region(var.chromosome(false), var.position(false)));
                if !ig_candidate {
                    continue;
                }
            }
            let gs = self.breakend_gene_list(var, true);
            let ge = self.breakend_gene_list(var, false);
            if gs.is_empty() || ge.is_empty() {
                continue;
            }
            let mut fusions = self.finder.find_fusions(&gs, &ge);
            if let Some(n) = self.neo.as_mut() {
                let cn = if var.is_sgl_breakend() { var.copy_number(true) } else { (var.copy_number(true) + var.copy_number(false)) * 0.5 };
                n.process_fusion_candidate(&gs, &ge, None, None, None, None, cn)?;
            }
            if fusions.is_empty() {
                continue;
            }
            if self.finder.config.log_reportable_only {
                fusions.retain(|x| is_reportable(x));
            }
            for f in fusions {
                f.set_annotations(FusionAnnotations {
                    cluster_id: cluster.id(),
                    cluster_count: cluster.sv_count() as i32,
                    resolved_type: cluster.resolved_type().to_string(),
                    chain_info: None,
                    terminated_up: false,
                    terminated_down: false,
                });
                self.fusions.push(f);
            }
        }
        Ok(())
    }

    fn find_chained_fusions(&mut self, clusters: &[Cluster]) -> Result<()> {
        for cluster in clusters {
            if cluster.sv_count() == 1 || cluster.chains().is_empty() {
                continue;
            }
            let mut chain_fusions: Vec<GeneFusion> = Vec::new();
            let mut valid_pairs: Vec<ValidTraversalData> = Vec::new();
            for chain in cluster.chains() {
                self.find_chain_fusions(cluster, &chain, &mut chain_fusions, &mut valid_pairs)?;
            }
            let log_reportable_only = self.finder.config.log_reportable_only;
            self.fusions.extend(chain_fusions.into_iter().filter(|x| !log_reportable_only || x.reportable()).filter(|x| x.has_annotations()));
        }
        Ok(())
    }

    fn find_chain_fusions(
        &mut self, cluster: &Cluster, chain: &Chain, chain_fusions: &mut Vec<GeneFusion>, valid_pairs: &mut Vec<ValidTraversalData>,
    ) -> Result<()> {
        let links = chain.linked_pairs();
        let n = links.len();
        for lp1 in 0..=n {
            let (lower_sv, lower_breakend) = if lp1 < n {
                (links[lp1].first(), links[lp1].first_breakend().other_breakend())
            } else {
                (links[lp1 - 1].second(), Some(links[lp1 - 1].second_breakend()))
            };
            if lower_sv.is_sgl_breakend() && lower_sv.sgl_mappings().is_empty() {
                continue;
            }
            let lower_pos;
            let genes_lower;
            if let Some(lb) = &lower_breakend {
                lower_pos = lb.position();
                genes_lower = self.breakend_gene_list(&lower_sv, lb.uses_start());
            } else {
                genes_lower = self.breakend_gene_list(&lower_sv, false);
                if genes_lower.is_empty() {
                    continue;
                }
                lower_pos = genes_lower[0].position();
            }
            if genes_lower.is_empty() {
                continue;
            }
            let mut traversed: Vec<Pair> = Vec::new();
            for lp2 in lp1..=n {
                let (upper_sv, upper_breakend): (Sv, Option<Be>) = if lp2 < n {
                    (links[lp2].first(), Some(links[lp2].first_breakend()))
                } else if chain.is_closed_loop() {
                    (links[lp2 - 1].second(), links[lp2 - 1].second_breakend().other_breakend())
                } else {
                    (chain.chain_end_sv(false).expect("chain end SV"), chain.open_breakend(false))
                };
                let genes_upper = if let Some(ub) = &upper_breakend {
                    self.breakend_gene_list(&upper_sv, ub.uses_start())
                } else {
                    let g = self.breakend_gene_list(&upper_sv, false);
                    if g.is_empty() {
                        continue;
                    }
                    g
                };
                if lp2 > lp1 {
                    traversed.push(links[lp2 - 1].clone());
                }
                if genes_upper.is_empty() {
                    continue;
                }
                let mut fusions = self.finder.find_fusions(&genes_lower, &genes_upper);
                if let Some(neo) = self.neo.as_mut() {
                    let lower_link = if lp1 == 0 { None } else { Some(&links[lp1 - 1]) };
                    let upper_link = if lp2 < n { Some(&links[lp2]) } else { None };
                    let cn = match (&upper_breakend, &lower_breakend) {
                        (Some(u), Some(l)) => (l.copy_number() + u.copy_number()) * 0.5,
                        _ => lower_sv.copy_number(true),
                    };
                    neo.process_fusion_candidate(
                        &genes_lower,
                        &genes_upper,
                        Some(&traversed),
                        Some(&self.disruption_finder),
                        lower_link,
                        upper_link,
                        cn,
                    )?;
                }
                if fusions.is_empty() {
                    continue;
                }
                if lp2 > lp1 {
                    fusions.retain(|x| !x.is_exonic());
                }
                if self.finder.config.log_reportable_only {
                    fusions.retain(|x| is_reportable(x));
                    if fusions.is_empty() {
                        continue;
                    }
                }
                let mut valid_traversal_count = 0;
                for fusion in fusions {
                    let up = fusion.upstream_trans().clone();
                    let up_gene_strand = up.gene().strand() as i32;
                    let is_precoding_upstream = up.pre_coding();
                    let fusion_lower_to_upper = up.gene().position() == lower_pos;
                    let downstream_genes =
                        if lower_pos == fusion.downstream_trans().gene().position() { &genes_lower } else { &genes_upper };
                    let non_disruptive_chain = if !traversed.is_empty() && !downstream_genes.is_empty() {
                        self.disruption_finder.is_non_disruptive_chained_pair(&up, downstream_genes)
                    } else {
                        false
                    };
                    let mut total_link_length: i64 = 0;
                    let mut valid_traversal = true;
                    let allow_invalid = allow_suspect_chains(fusion.known_type());
                    let mut all_assembled = true;
                    for pair in &traversed {
                        total_link_length += pair.base_length() as i64;
                        if pair.is_inferred() {
                            all_assembled = false;
                        }
                        let orient = if fusion_lower_to_upper { pair.first_breakend().orientation() } else { pair.second_breakend().orientation() }
                            as i32;
                        let fusion_direction = if orient != up_gene_strand { up_gene_strand } else { -up_gene_strand };
                        if valid_traversal && !self.has_valid_traversal(valid_pairs, pair, fusion_direction, is_precoding_upstream) {
                            valid_traversal = false;
                            if !allow_invalid {
                                break;
                            }
                        }
                    }
                    if !valid_traversal && !allow_invalid {
                        self.record_invalid_fusion(&fusion, "InvalidTraversal");
                        continue;
                    }
                    if valid_traversal {
                        valid_traversal_count += 1;
                    }
                    let mut terminated = [false, false];
                    for fs in [FS_UP, FS_DOWN] {
                        let transcript = if fs == FS_UP { fusion.upstream_trans() } else { fusion.downstream_trans() };
                        let is_lower = lower_pos == transcript.gene().position();
                        let is_chain_end = (is_lower && lp1 == 0) || (!is_lower && lp2 == n);
                        if is_chain_end {
                            terminated[fs] = false;
                        } else {
                            let breakend = if is_lower { &lower_breakend } else { &upper_breakend };
                            match breakend {
                                None => terminated[fs] = false,
                                Some(b) => {
                                    let link_index = if is_lower { lp1 as i32 - 1 } else { lp2 as i32 };
                                    if link_index as usize >= n {
                                        hmf_common::error!(
                                            "cluster({}) link index({}) exceeds chain link size({})",
                                            cluster.id(),
                                            link_index,
                                            n
                                        );
                                        break;
                                    }
                                    terminated[fs] = check_transcript_disruption_info(b, transcript, chain, link_index);
                                }
                            }
                        }
                    }
                    let links_count = (lp2 - lp1) as i32;
                    fusion.set_annotations(FusionAnnotations {
                        cluster_id: cluster.id(),
                        cluster_count: cluster.sv_count() as i32,
                        resolved_type: cluster.resolved_type().to_string(),
                        chain_info: Some(FusionChainInfo {
                            chain_id: chain.id(),
                            chain_links: links_count,
                            chain_length: total_link_length as i32,
                            valid_traversal,
                            traversal_assembled: all_assembled,
                            non_disruptive: non_disruptive_chain,
                        }),
                        terminated_up: terminated[SE_START],
                        terminated_down: terminated[SE_END],
                    });
                    let chain_length_ok = total_link_length <= FUSION_MAX_CHAIN_LENGTH;
                    let traversal_ok = valid_traversal || allow_invalid;
                    if chain_length_ok && traversal_ok && !non_disruptive_chain {
                        if !has_identical_fusion(&fusion, chain_fusions) {
                            chain_fusions.push(fusion);
                        }
                    } else {
                        let reason = if !valid_traversal {
                            "TraversesSPA"
                        } else if !chain_length_ok {
                            "LongChain"
                        } else {
                            "NonDisruptiveChain"
                        };
                        self.record_invalid_fusion(&fusion, reason);
                    }
                }
                if valid_traversal_count == 0 && lp2 > lp1 {
                    break;
                }
            }
        }
        Ok(())
    }

    fn record_invalid_fusion(&mut self, fusion: &GeneFusion, reason: &str) {
        let name = fusion.name();
        if self.invalid_fusions.iter().any(|(f, _)| f.name() == name) {
            return;
        }
        self.invalid_fusions.push((fusion.clone(), reason.to_string()));
    }

    fn has_valid_traversal(&self, valid_pairs: &mut Vec<ValidTraversalData>, pair: &Pair, fusion_direction: i32, precoding: bool) -> bool {
        if let Some(e) = valid_pairs
            .iter()
            .find(|x| x.fusion_direction == fusion_direction && x.is_precoding_upstream == precoding && x.link.matches(pair))
        {
            return e.is_valid;
        }
        let valid = !self.disruption_finder.pair_traverses_gene(pair, fusion_direction, precoding);
        valid_pairs.push(ValidTraversalData { link: pair.clone(), is_valid: valid, fusion_direction, is_precoding_upstream: precoding });
        valid
    }

    fn extract_unique_fusions(&self) -> Vec<GeneFusion> {
        let mut sv_id_pair_fusions: JavaHashMap<String, Vec<GeneFusion>> = JavaHashMap::new();
        for f in &self.fusions {
            if !persist_fusion(f) {
                continue;
            }
            let key = f.sv_id_pair();
            if let Some(l) = sv_id_pair_fusions.get_mut(&key) {
                l.push(f.clone());
            } else {
                sv_id_pair_fusions.insert(key, vec![f.clone()]);
            }
        }
        let mut gene_pair_fusions: JavaHashMap<String, Vec<GeneFusion>> = JavaHashMap::new();
        for (_, similar) in sv_id_pair_fusions.iter() {
            let mut candidates: Vec<GeneFusion> = similar
                .iter()
                .filter(|x| x.known_type() == ENHANCER_PROMISCUOUS && x.downstream_trans().trans_data.is_canonical)
                .cloned()
                .collect();
            let top = self.finder.find_top_reportable_fusion(similar);
            if let Some(t) = top {
                if !contains(&candidates, &t) {
                    candidates.push(t);
                }
            }
            for f in candidates {
                let name = f.name();
                if let Some(l) = gene_pair_fusions.get_mut(&name) {
                    l.push(f);
                } else {
                    gene_pair_fusions.insert(name, vec![f]);
                }
            }
        }
        let mut unique = Vec::new();
        for (_, fusions) in gene_pair_fusions.iter() {
            if fusions.len() == 1 {
                unique.push(fusions[0].clone());
            } else if let Some(t) = find_top_priority_fusion(fusions) {
                unique.push(t);
            }
        }
        unique
    }

    fn build_transcript_list(&self, sv_list: &[Sv]) -> Vec<BreakendTrans> {
        let mut transcripts: Vec<BreakendTrans> = Vec::new();
        for var in sv_list {
            for be in [SE_START, SE_END] {
                if var.is_sgl_breakend() && be == SE_END {
                    continue;
                }
                for gene in var.genes_list(is_start(be)) {
                    for t in gene.transcripts() {
                        if t.is_canonical() || self.disruption_finder.matches_disruption_transcript(gene.gene_id(), t.trans_data) {
                            transcripts.push(t);
                        }
                    }
                }
            }
        }
        for f in &self.unique_fusions {
            if !contains(&transcripts, f.upstream_trans()) {
                transcripts.push(f.upstream_trans().clone());
            }
            if !contains(&transcripts, f.downstream_trans()) {
                transcripts.push(f.downstream_trans().clone());
            }
        }
        for t in &transcripts {
            if t.undisrupted_copy_number_set() {
                continue;
            }
            if let Some(var) = sv_list.iter().find(|x| x.id() == t.gene().var_id()) {
                if let Some(b) = var.breakend(t.gene().is_start()) {
                    t.set_undisrupted_copy_number(undisrupted_copy_number(&b));
                }
            }
        }
        transcripts
    }

    fn add_visualisation_data(&self, vis: &mut VisSampleData) {
        if !self.write_vis {
            return;
        }
        let mut vis_fusions = Vec::new();
        for f in &self.unique_fusions {
            if !(f.reportable() || self.finder.config.write_all_vis_fusions) {
                continue;
            }
            let cluster_id = f.annotations().map(|a| a.cluster_id).unwrap_or(-1);
            let up = f.upstream_trans();
            let down = f.downstream_trans();
            for t in [up, down] {
                vis.add_gene_exon_data(VisGeneData::new(
                    cluster_id,
                    t.gene().gene_id(),
                    t.gene_name(),
                    t.trans_name(),
                    t.trans_id(),
                    &t.gene().chromosome(),
                    VisGeneAnnotationType::FUSION,
                ));
            }
            vis_fusions.push(VisFusion {
                cluster_id,
                reportable: f.reportable(),
                gene_name_up: up.gene_name().to_string(),
                transcript_up: up.trans_name().to_string(),
                chr_up: up.gene().chromosome(),
                pos_up: up.gene().position(),
                strand_up: up.gene().strand() as i32,
                region_type_up: up.region_type().to_string(),
                fused_exon_up: f.fused_exon(true),
                gene_name_down: down.gene_name().to_string(),
                transcript_down: down.trans_name().to_string(),
                chr_down: down.gene().chromosome(),
                pos_down: down.gene().position(),
                strand_down: down.gene().strand() as i32,
                region_type_down: down.region_type().to_string(),
                fused_exon_down: f.fused_exon(false),
            });
        }
        vis.fusions.extend(vis_fusions);
    }
}

fn persist_fusion(f: &GeneFusionData) -> bool {
    if !f.valid_chain_traversal() && !allow_suspect_chains(f.known_type()) {
        return false;
    }
    !(!f.is_enhancer() && f.downstream_trans().has_negative_prev_splice_acceptor_distance())
}

fn has_identical_fusion(new: &GeneFusionData, fusions: &[GeneFusion]) -> bool {
    fusions.iter().any(|f| {
        new.upstream_trans().gene().var_id() == f.upstream_trans().gene().var_id()
            && new.downstream_trans().gene().var_id() == f.downstream_trans().gene().var_id()
            && new.upstream_trans().trans_name() == f.upstream_trans().trans_name()
            && new.downstream_trans().trans_name() == f.downstream_trans().trans_name()
            && new.downstream_trans().has_negative_prev_splice_acceptor_distance()
                == f.downstream_trans().has_negative_prev_splice_acceptor_distance()
    })
}

fn check_transcript_disruption_info(breakend: &Be, transcript: &BreakendTrans, chain: &Chain, mut link_index: i32) -> bool {
    let links = chain.linked_pairs();
    let traverse_up = &links[link_index as usize].first_breakend() == breakend;
    let is_upstream = transcript.is_upstream();
    while link_index >= 0 && (link_index as usize) < links.len() {
        let pair = &links[link_index as usize];
        let next = if traverse_up { pair.second_breakend() } else { pair.first_breakend() };
        if next.forward_orient() {
            if next.position() > transcript.trans_end() {
                break;
            } else if !is_upstream && transcript.coding_end() > 0 && next.position() > transcript.coding_end() {
                break;
            }
        } else if next.position() < transcript.trans_start() {
            break;
        } else if !is_upstream && transcript.coding_start() > 0 && next.position() < transcript.coding_start() {
            break;
        }
        if is_disruptive_in_transcript(next.sv(), transcript.trans_data) {
            return true;
        }
        link_index += if traverse_up { 1 } else { -1 };
    }
    false
}

/// FusionWriter.convertBreakendsAndFusions
pub fn convert_breakends_and_fusions(
    gene_fusions: &[GeneFusion], driver_genes: &JavaHashMap<String, DriverGene>, transcripts: &[BreakendTrans],
) -> (Vec<LinxFusion>, Vec<LinxBreakend>) {
    let mut breakends = Vec::new();
    let mut fusions = Vec::new();
    for (id, t) in transcripts.iter().enumerate() {
        let g = t.gene();
        let mut reported = "NONE";
        if let Some(dg) = driver_genes.get(&g.gene_data.gene_name) {
            if t.reportable_disruption() {
                reported = if dg.report_disruption { "REPORTED" } else { "NOT_REPORTED" };
            }
        }
        breakends.push(LinxBreakend {
            id: id as i32,
            sv_id: g.var_id(),
            vcf_id: g.vcf_id(),
            coords: g.coords_str(),
            is_start: g.is_start(),
            gene: t.gene_name().to_string(),
            transcript_id: t.trans_name().to_string(),
            canonical: t.is_canonical(),
            gene_orientation: if t.is_upstream() { BREAKEND_ORIENTATION_UPSTREAM } else { BREAKEND_ORIENTATION_DOWNSTREAM },
            disruptive: t.is_disruptive(),
            reported_status: reported.to_string(),
            undisrupted_copy_number: t.undisrupted_copy_number(),
            region_type: t.region_type(),
            coding_type: t.coding_type(),
            biotype: t.bio_type().to_string(),
            exonic_base_phase: t.phase,
            next_splice_exon_rank: t.next_splice_exon_rank(),
            next_splice_exon_phase: t.phase,
            next_splice_distance: if t.is_upstream() { t.prev_splice_acceptor_distance() } else { t.next_splice_acceptor_distance() },
            total_exon_count: t.trans_data.exons.len() as i32,
            exon_up: t.exon_upstream,
            exon_down: t.exon_downstream,
        });
    }
    // transIdMap: the last breakend of each transcript
    let breakend_id = |t: &BreakendTrans| -> i32 { transcripts.iter().rposition(|x| x == t).expect("fusion transcript breakend") as i32 };
    for f in gene_fusions {
        let up = f.upstream_trans().gene();
        let down = f.downstream_trans().gene();
        fusions.push(LinxFusion {
            five_prime_breakend_id: breakend_id(f.upstream_trans()),
            three_prime_breakend_id: breakend_id(f.downstream_trans()),
            name: f.name(),
            reported: f.reportable(),
            reported_type: f.known_type_str(),
            reportable_reasons: f.reportable_reasons().iter().map(|r| r.to_string()).collect::<Vec<_>>().join(";"),
            phased: format!("{:?}", f.phase_type()),
            likelihood: format!("{:?}", f.likelihood_type()),
            five_prime_vcf_id: up.vcf_id(),
            three_prime_vcf_id: down.vcf_id(),
            five_prime_coords: up.coords_str(),
            three_prime_coords: down.coords_str(),
            chain_length: f.chain_length(),
            chain_links: f.chain_links(),
            chain_terminated: f.is_terminated(),
            domains_kept: f.downstream_trans().protein_features_kept(),
            domains_lost: f.downstream_trans().protein_features_lost(),
            skipped_exons_up: f.exons_skipped(true),
            skipped_exons_down: f.exons_skipped(false),
            fused_exon_up: f.fused_exon(true),
            fused_exon_down: f.fused_exon(false),
        });
    }
    (fusions, breakends)
}
