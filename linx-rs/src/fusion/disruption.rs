//! DisruptionFinder: which breakend transcripts are disrupted, and the reportable disruptions of driver
//! genes; with GermlineDisruptions (package `germline`) for germline mode.

use std::collections::HashMap;

use hmf_common::driver_panel::{DriverGene, GermlineReporting};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::purple_files::{write_driver_catalog, DriverCatalog};
use hmf_common::sv::StructuralVariantType::*;

use super::SvDisruptionData;
use crate::analysis::find_end_index;
use crate::analysis::find_start_index;
use crate::annotators::is_pseudogene_deletion;
use crate::gene::{BreakendGene, BreakendTrans, Gene, GeneCache, Trans, TranscriptCodingType, TranscriptRegionType};
use crate::idmap::IdMap;
use crate::output::{breakend_lines, LinxBreakend, BREAKEND_ORIENTATION_DOWNSTREAM, BREAKEND_ORIENTATION_UPSTREAM};
use crate::types::*;
use crate::vis::{VisGeneAnnotationType, VisGeneData, VisSampleData};

pub const MAX_NON_DISRUPTED_CHAIN_LENGTH: i32 = 5000;

pub struct DisruptionFinder {
    gene_cache: &'static GeneCache,
    /// gene id to its additional reported transcripts
    disruption_gene_transcripts: HashMap<String, Vec<String>>,
    driver_genes: &'static JavaHashMap<String, DriverGene>,
    pub disruptions: Vec<SvDisruptionData>,
    is_germline: bool,
    pub germline: Option<GermlineDisruptions>,
}

impl DisruptionFinder {
    pub fn new(is_germline: bool, gene_cache: &'static GeneCache, driver_genes: &'static JavaHashMap<String, DriverGene>) -> Self {
        let mut map = HashMap::new();
        for dg in driver_genes.values() {
            if !is_germline && !dg.report_disruption {
                continue;
            }
            if let Some(g) = gene_cache.gene_data_by_name(&dg.gene) {
                map.insert(g.gene_id.clone(), dg.additional_reported_transcripts.clone());
            }
        }
        DisruptionFinder {
            gene_cache,
            disruption_gene_transcripts: map,
            driver_genes,
            disruptions: Vec::new(),
            is_germline,
            germline: if is_germline { Some(GermlineDisruptions::new(gene_cache, driver_genes)) } else { None },
        }
    }

    /// addDisruptionGene
    pub fn add_disruption_gene(&mut self, gene_id: &str) {
        self.disruption_gene_transcripts.insert(gene_id.to_string(), Vec::new());
    }

    fn matches_disruption_gene(&self, gene_id: &str) -> bool {
        self.disruption_gene_transcripts.contains_key(gene_id)
    }

    pub fn matches_disruption_transcript(&self, gene_id: &str, trans: &hmf_common::ensembl::TranscriptData) -> bool {
        let Some(other) = self.disruption_gene_transcripts.get(gene_id) else { return false };
        if trans.is_canonical {
            return true;
        }
        other.contains(&trans.trans_name)
    }

    /// markTranscriptsDisruptive
    pub fn mark_transcripts_disruptive(&self, sv_list: &[Sv]) {
        let mut clusters: Vec<Cluster> = Vec::new();
        for var in sv_list {
            self.mark_sv_transcripts_disruptive(var);
            if var.is_inferred_sgl() {
                for g in var.genes_list(true) {
                    for t in g.transcripts() {
                        t.set_is_disruptive(false);
                    }
                }
            }
            let cluster = var.cl();
            if !cluster.chains().is_empty() && !contains(&clusters, &cluster) {
                clusters.push(cluster);
            }
        }
        // a HashSet in Java: each cluster's checks touch only its own breakends, so the order does not matter
        for c in &clusters {
            for chain in c.chains() {
                self.check_chain_through_introns(&chain);
                self.check_traverse_to_same_intron(&chain);
            }
        }
        for c in &clusters {
            for chain in c.chains() {
                self.check_other_ends_non_genic(&chain);
            }
        }
    }

    fn mark_sv_transcripts_disruptive(&self, var: &Sv) {
        let genes_start = var.genes_list(true);
        let genes_end = var.genes_list(false);
        if genes_start.is_empty() && genes_end.is_empty() {
            return;
        }
        let cluster = var.cl();
        for s in [SE_START, SE_END] {
            if s == SE_END && var.is_sgl_breakend() {
                continue;
            }
            let genes = if s == SE_START { &genes_start } else { &genes_end };
            for gene in genes {
                if cluster.resolved_type() == ResolvedType::LINE {
                    for t in gene.transcripts() {
                        if t.is_intronic() && t.is_disruptive() {
                            mark_non_disruptive(&t);
                        }
                    }
                }
            }
        }
        if var.is_simple_type() {
            self.mark_non_disruptive_gene_transcripts(&genes_start, &genes_end);
            if var.sv_type() == DUP {
                self.mark_non_disruptive_dup(var);
            }
        }
    }

    fn mark_non_disruptive_gene_transcripts(&self, genes_start: &[BreakendGene], genes_end: &[BreakendGene]) {
        for gs in genes_start {
            let Some(ge) = genes_end.iter().find(|x| x.gene_id() == gs.gene_id()) else { continue };
            self.mark_non_disruptive_transcripts(&gs.transcripts(), &ge.transcripts());
        }
    }

    fn mark_non_disruptive_transcripts(&self, list1: &[BreakendTrans], list2: &[BreakendTrans]) -> bool {
        let mut found = false;
        for t1 in list1 {
            let Some(t2) = list2.iter().find(|x| x.trans_id() == t1.trans_id()) else { continue };
            if self.mark_non_disruptive_same_transcript(t1, t2) {
                found = true;
            }
        }
        found
    }

    fn mark_non_disruptive_same_transcript(&self, t1: &BreakendTrans, t2: &BreakendTrans) -> bool {
        let mut mark = false;
        if t1.exon_upstream == t2.exon_upstream && !t1.is_exonic() && !t2.is_exonic() {
            mark = true;
        } else if t1.gene().var_id() == t2.gene().var_id() && self.is_germline {
            use TranscriptCodingType::UTR_3P;
            if t1.coding_type() == UTR_3P && t2.coding_type() == UTR_3P {
                mark = true;
            } else if (t1.coding_type() == UTR_3P || t2.coding_type() == UTR_3P) && t1.gene().sv_type() == Some(DUP) {
                mark = true;
            }
        }
        if mark {
            mark_non_disruptive(t1);
            mark_non_disruptive(t2);
            return true;
        }
        false
    }

    fn mark_non_disruptive_dup(&self, var: &Sv) {
        if var.cl().sv_count() > 1 {
            return;
        }
        for s in [SE_START, SE_END] {
            let genes = var.genes_list(is_start(s));
            if genes.is_empty() {
                continue;
            }
            let others = var.genes_list(!is_start(s));
            if !others.is_empty() && others.iter().any(|x| genes.iter().any(|y| y.gene_name() != x.gene_name())) {
                continue;
            }
            for gene in &genes {
                if others.is_empty() {
                    gene.transcripts().iter().for_each(mark_non_disruptive);
                    continue;
                }
                if let Some(other) = others.iter().find(|x| x.gene_id() == gene.gene_id()) {
                    for t in gene.transcripts() {
                        let matching = other.transcripts().into_iter().find(|x| x.trans_id() == t.trans_id());
                        match matching {
                            None => mark_non_disruptive(&t),
                            Some(m) => {
                                use TranscriptRegionType::*;
                                let within = t.region_type() == EXONIC || t.region_type() == INTRONIC;
                                let within_other = m.region_type() == EXONIC || m.region_type() == INTRONIC;
                                if !within || !within_other {
                                    mark_non_disruptive(&t);
                                    mark_non_disruptive(&m);
                                }
                            }
                        }
                    }
                    if gene.is_upstream() {
                        for up in gene.transcripts() {
                            let down = other.transcripts().into_iter().find(|x| x.trans_id() == up.trans_id());
                            if let Some(down) = down {
                                if up.exon_upstream == 1 && down.exon_downstream <= 2 && !up.is_exonic() {
                                    mark_non_disruptive(&up);
                                    mark_non_disruptive(&down);
                                }
                            }
                        }
                    }
                } else {
                    gene.transcripts().iter().for_each(mark_non_disruptive);
                }
            }
        }
    }

    fn check_chain_through_introns(&self, chain: &Chain) {
        let (Some(cs), Some(ce)) = (chain.open_breakend(true), chain.open_breakend(false)) else { return };
        if cs.sv().sv_type() == DEL || ce.sv().sv_type() == DEL {
            return;
        }
        let gs = cs.genes_list();
        let ge = ce.genes_list();
        if !gs.is_empty() || !ge.is_empty() {
            if gs.iter().any(|x| x.transcripts().iter().any(|y| y.is_disruptive())) {
                return;
            }
            if ge.iter().any(|x| x.transcripts().iter().any(|y| y.is_disruptive())) {
                return;
            }
        }
        for pair in chain.linked_pairs() {
            if !pair.is_assembled() {
                continue;
            }
            let b1 = pair.first_breakend();
            let b2 = pair.second_breakend();
            if b1.sv().sv_type() == DEL || b2.sv().sv_type() == DEL {
                continue;
            }
            let genes1 = b1.genes_list();
            let genes2 = b2.genes_list();
            if genes1.is_empty() || genes2.is_empty() {
                continue;
            }
            if genes1.iter().any(|x| in_deletion_bridge_within_gene(chain, &b1, x.gene_name())) {
                continue;
            } else if genes2.iter().any(|x| in_deletion_bridge_within_gene(chain, &b2, x.gene_name())) {
                continue;
            }
            for g1 in &genes1 {
                let Some(g2) = genes2.iter().find(|x| x.gene_id() == g1.gene_id()) else { continue };
                for t1 in g1.transcripts() {
                    if t1.region_type() != TranscriptRegionType::INTRONIC {
                        continue;
                    }
                    let t2 = g2.transcripts().into_iter().find(|x| {
                        x.trans_id() == t1.trans_id() && x.region_type() == TranscriptRegionType::INTRONIC && x.exon_downstream == t1.exon_downstream
                    });
                    let Some(t2) = t2 else { continue };
                    if !t1.is_disruptive() && !t2.is_disruptive() {
                        continue;
                    }
                    self.mark_non_disruptive_same_transcript(&t1, &t2);
                }
            }
        }
    }

    fn check_traverse_to_same_intron(&self, chain: &Chain) {
        if chain.open_breakend(true).is_none() || chain.open_breakend(false).is_none() {
            return;
        }
        let links = chain.linked_pairs();
        let mut link_index: i32 = -1;
        loop {
            link_index += 1;
            if link_index as usize >= links.len() {
                break;
            }
            let pair = &links[link_index as usize];
            // Java dereferences the other breakend unchecked: an SGL cannot start a link's lower side here
            let Some(breakend) = pair.first_breakend().other_breakend() else { continue };
            let genes = breakend.genes_list();
            if genes.is_empty() {
                continue;
            }
            let mut trans_list = disrupted_transcripts(&genes);
            if trans_list.is_empty() {
                continue;
            }
            let mut next_link_index = link_index - 1;
            let mut found = false;
            let mut traversed: Vec<Pair> = Vec::new();
            let mut chain_length = 0;
            loop {
                next_link_index += 1;
                if next_link_index as usize >= links.len() {
                    break;
                }
                let next_pair = &links[next_link_index as usize];
                if self.pair_traverses_gene(next_pair, 0, false) {
                    break;
                }
                chain_length += next_pair.base_length();
                if chain_length > MAX_NON_DISRUPTED_CHAIN_LENGTH {
                    break;
                }
                traversed.push(next_pair.clone());
                let Some(next_breakend) = next_pair.second_breakend().other_breakend() else { continue };
                if next_breakend.orientation() == breakend.orientation() || next_breakend.chromosome() != breakend.chromosome() {
                    continue;
                }
                let other_genes = next_breakend.genes_list();
                if other_genes.is_empty() || !other_genes.iter().any(|x| genes.iter().any(|y| x.gene_name() == y.gene_name())) {
                    continue;
                }
                let other_trans = disrupted_transcripts(&other_genes);
                if !other_trans.is_empty() && self.mark_non_disruptive_transcripts(&trans_list, &other_trans) {
                    found = true;
                    for tp in &traversed {
                        let tt1 = disrupted_transcripts(&tp.first_breakend().genes_list());
                        let tt2 = disrupted_transcripts(&tp.second_breakend().genes_list());
                        if !tt1.is_empty() && !tt2.is_empty() {
                            self.mark_non_disruptive_transcripts(&tt1, &tt2);
                        }
                    }
                    trans_list.retain(|t| t.is_disruptive());
                }
                break;
            }
            if found {
                link_index = next_link_index;
            }
        }
    }

    fn check_other_ends_non_genic(&self, chain: &Chain) {
        for pair in chain.linked_pairs() {
            let b1 = pair.first_breakend();
            let b2 = pair.second_breakend();
            if b1.sv().sv_type() == DEL || b2.sv().sv_type() == DEL {
                continue;
            }
            let genes1 = b1.genes_list();
            let genes2 = b2.genes_list();
            if genes1.is_empty() || genes2.is_empty() {
                continue;
            }
            let tl1 = disrupted_transcripts(&genes1);
            let tl2 = disrupted_transcripts(&genes2);
            if tl1.is_empty() || tl2.is_empty() {
                continue;
            }
            if genes1.iter().any(|x| in_deletion_bridge_within_gene(chain, &b1, x.gene_name())) {
                continue;
            } else if genes2.iter().any(|x| in_deletion_bridge_within_gene(chain, &b2, x.gene_name())) {
                continue;
            }
            let o1 = b1.other_breakend().is_some_and(|o| o.genes_list().is_empty());
            let o2 = b2.other_breakend().is_some_and(|o| o.genes_list().is_empty());
            if o1 && o2 {
                self.mark_non_disruptive_transcripts(&tl1, &tl2);
            }
        }
    }

    /// pairTraversesGene: whether a templated insertion crosses a splice acceptor
    pub fn pair_traverses_gene(&self, pair: &Pair, fusion_direction: i32, is_precoding_upstream: bool) -> bool {
        let lower = pair.breakend(true).position();
        let upper = pair.breakend(false).position();
        let Some(genes) = self.gene_cache.ensembl.chr_gene_data_map().get(&pair.chromosome()) else { return false };
        for g in genes {
            if lower > g.gene_end {
                continue;
            }
            if upper < g.gene_start {
                break;
            }
            if fusion_direction == 0 || g.strand as i32 == fusion_direction {
                let Some(trans_list) = self.gene_cache.transcripts(&g.gene_id) else { continue };
                for t in trans_list {
                    for e in &t.exons {
                        if e.rank == 1 {
                            continue;
                        }
                        if (g.strand == 1 && lower <= e.start && upper >= e.start) || (g.strand == -1 && lower <= e.end && upper >= e.end) {
                            if is_precoding_upstream && lower <= e.start && upper >= e.end {
                                if g.strand == 1 && t.coding_start.is_none_or(|cs| upper < cs) {
                                    continue;
                                } else if g.strand == -1 && t.coding_end.is_none_or(|ce| lower > ce) {
                                    continue;
                                }
                            }
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// findReportableDisruptions
    pub fn find_reportable_disruptions(&mut self, sv_list: &[Sv], clusters: &[Cluster], vis: Option<&mut VisSampleData>) {
        self.disruptions.clear();
        let mut vis = vis;
        for var in sv_list {
            for be in [SE_START, SE_END] {
                if be == SE_END && var.is_sgl_breakend() {
                    continue;
                }
                let genes: Vec<BreakendGene> =
                    var.genes_list(is_start(be)).into_iter().filter(|x| self.matches_disruption_gene(x.gene_id())).collect();
                for gene in &genes {
                    let disrupted: Vec<BreakendTrans> = gene.transcripts().into_iter().filter(|x| x.is_disruptive()).collect();
                    for t in disrupted {
                        if !self.matches_disruption_transcript(gene.gene_id(), t.trans_data) {
                            continue;
                        }
                        let breakend = var.be(is_start(be));
                        if !t.undisrupted_copy_number_set() {
                            t.set_undisrupted_copy_number(undisrupted_copy_number(&breakend));
                        }
                        t.set_reportable_disruption(true);
                        let d = SvDisruptionData::new(
                            var,
                            gene.is_start(),
                            gene.gene_data,
                            t.trans_data,
                            [t.exon_upstream, t.exon_downstream],
                            t.coding_type(),
                            t.region_type(),
                            t.undisrupted_copy_number(),
                            breakend.copy_number(),
                        );
                        d.set_reportable(true);
                        if var.sv_type() == DEL && is_pseudogene_deletion(var, var.position(true), var.position(false), t.trans_data) {
                            d.mark_pseudogene_deletion();
                        }
                        self.disruptions.push(d);
                        if let Some(v) = vis.as_deref_mut() {
                            v.add_gene_exon_data(VisGeneData::new(
                                var.cl().id(),
                                gene.gene_id(),
                                gene.gene_name(),
                                "",
                                0,
                                &gene.chromosome(),
                                VisGeneAnnotationType::DISRUPTION,
                            ));
                        }
                    }
                }
            }
        }
        if let Some(g) = self.germline.as_mut() {
            g.find_gene_deletions(clusters);
        }
    }

    /// addReportableDisruptions
    pub fn add_reportable_disruptions(&self, catalog: &mut Vec<DriverCatalog>) {
        for d in &self.disruptions {
            if !d.reportable() {
                continue;
            }
            if catalog.iter().filter(|x| x.gene == d.gene.gene_name).any(|x| {
                x.driver == "HOM_DUP_DISRUPTION" || x.driver == "HOM_DEL_DISRUPTION" || x.driver == "DISRUPTION"
            }) {
                continue;
            }
            let Some(dg) = self.driver_genes.get(&d.gene.gene_name) else { continue };
            let reported = if dg.report_disruption { "REPORTED" } else { "NOT_REPORTED" };
            let biallelic = d.undisrupted_copy_number < DISRUPTION_CN_REPORTABLE_THRESHOLD;
            catalog.push(DriverCatalog {
                chromosome: d.gene.chromosome.clone(),
                chromosome_band: d.gene.karyotype_band.clone(),
                gene: d.gene.gene_name.clone(),
                transcript: d.transcript.trans_name.clone(),
                is_canonical: d.transcript.is_canonical,
                driver: "DISRUPTION".into(),
                category: dg.likelihood_type.name().into(),
                likelihood_method: "DISRUPTION".into(),
                reported_status: reported.into(),
                driver_likelihood: if biallelic { 1.0 } else { 0.0 },
                missense: 0,
                nonsense: 0,
                splice: 0,
                inframe: 0,
                frameshift: 0,
                biallelic,
                min_copy_number: d.undisrupted_copy_number,
                max_copy_number: d.max_copy_number,
            });
        }
    }

    /// writeGermlineDisruptions
    pub fn write_germline_disruptions(&self, sample_id: &str, output_dir: &str) -> anyhow::Result<()> {
        self.germline.as_ref().expect("germline disruptions").write_germline_svs(&self.disruptions, sample_id, output_dir)
    }

    /// isNonDisruptiveChainedPair
    pub fn is_non_disruptive_chained_pair(&self, up: &BreakendTrans, down_genes: &[BreakendGene]) -> bool {
        if up.region_type() == TranscriptRegionType::EXONIC {
            return false;
        }
        let Some(down_gene) = down_genes.iter().find(|x| x.gene_id() == up.gene().gene_id()) else { return false };
        for dt in down_gene.transcripts() {
            if dt.trans_id() != up.trans_id() {
                continue;
            }
            if dt.region_type() == up.region_type() && dt.exon_upstream == up.exon_upstream {
                return true;
            }
        }
        false
    }
}

fn mark_non_disruptive(t: &BreakendTrans) {
    if t.is_disruptive() {
        t.set_is_disruptive(false);
    }
}

fn disrupted_transcripts(genes: &[BreakendGene]) -> Vec<BreakendTrans> {
    let mut out = Vec::new();
    for g in genes {
        out.extend(g.transcripts().into_iter().filter(|x| x.is_disruptive()));
    }
    out
}

fn in_deletion_bridge_within_gene(chain: &Chain, breakend: &Be, gene_name: &str) -> bool {
    let Some(db) = breakend.db_link() else { return false };
    let other = db.other_breakend(breakend);
    if other.cluster() != breakend.cluster() {
        return false;
    }
    if other.cl().find_chain(other.sv()).as_ref() == Some(chain) {
        return false;
    }
    other.genes_list().iter().filter(|x| x.gene_name() == gene_name).any(|x| x.transcripts().iter().any(|y| y.is_disruptive()))
}

/// getUndisruptedCopyNumber
pub fn undisrupted_copy_number(breakend: &Be) -> f64 {
    let mut cn = breakend.copy_number_low_side();
    if let Some(db) = breakend.db_link() {
        if db.length() < 0 {
            cn -= db.other_breakend(breakend).jcn();
        }
    }
    cn
}

/// isDisruptiveInTranscript
pub fn is_disruptive_in_transcript(var: &Sv, trans: &hmf_common::ensembl::TranscriptData) -> bool {
    if var.is_sgl_breakend() {
        return true;
    }
    for s in [SE_START, SE_END] {
        let matching = var
            .genes_list(is_start(s))
            .iter()
            .filter(|x| x.gene_id() == trans.gene_id)
            .filter_map(|x| x.transcripts().into_iter().find(|y| y.trans_id() == trans.trans_id))
            .next();
        if matching.is_some_and(|m| m.is_disruptive()) {
            return true;
        }
    }
    false
}

// GermlineDisruptions

const MAX_DELETE_LENGTH: i32 = 3000000;
const MAX_PON_COUNT: i32 = 10;
const MAX_SGL_MAPPED_LENGTH: i32 = 500000;
const FILTER_PSEUDOGENE: &str = "PSEUDOGENE";

pub struct GermlineDisruptions {
    gene_cache: &'static GeneCache,
    driver_gene_data: Vec<Gene>,
    driver_genes: &'static JavaHashMap<String, DriverGene>,
    disruptions: Vec<SvDisruptionData>,
    reportable_sgls: Vec<Sv>,
}

impl GermlineDisruptions {
    pub fn new(gene_cache: &'static GeneCache, driver_genes: &'static JavaHashMap<String, DriverGene>) -> Self {
        let driver_gene_data = driver_genes.keys().filter_map(|n| gene_cache.gene_data_by_name(n)).collect();
        GermlineDisruptions { gene_cache, driver_gene_data, driver_genes, disruptions: Vec::new(), reportable_sgls: Vec::new() }
    }

    pub fn find_gene_deletions(&mut self, clusters: &[Cluster]) {
        self.disruptions.clear();
        self.reportable_sgls.clear();
        for cluster in clusters {
            if cluster.type_count(SGL) > 0 {
                for var in cluster.svs().iter().filter(|x| x.sv_type() == SGL) {
                    self.check_sgl_mappings(var);
                }
            }
            if cluster.sv_count() == 1 && cluster.sv0().sv_type() != DEL {
                continue;
            }
            self.check_cluster_gene_deletions(cluster);
        }
    }

    fn canonical(&self, g: Gene) -> Option<Trans> {
        let t = self.gene_cache.transcript_data(&g.gene_id, "");
        if t.is_none() {
            hmf_common::error!("gene({}:{}) missing canonical transcript", g.gene_id, g.gene_name);
        }
        t
    }

    fn check_cluster_gene_deletions(&mut self, cluster: &Cluster) {
        if cluster.resolved_type() == ResolvedType::LINE {
            return;
        }
        if cluster.sv_count() > 1 && !cluster.is_fully_chained(true) {
            return;
        }
        for (chromosome, list) in cluster.chr_breakend_entries() {
            let start_index = find_start_index(&list);
            let end_index = find_end_index(&list);
            let mut i = start_index;
            while i < end_index {
                let breakend = &list[i as usize];
                let var = breakend.sv().clone();
                let next = &list[i as usize + 1];
                let is_db = breakend.db_link().is_some() && breakend.db_link() == next.db_link();
                let is_simple_del = !is_db && var.sv_type() == DEL && breakend.forward_orient() && next.sv() == &var;
                if !is_db && !is_simple_del {
                    return;
                }
                let del_start = breakend.position();
                let del_end = next.position();
                if del_end - del_start > MAX_DELETE_LENGTH {
                    return;
                }
                for g in self.driver_gene_data.clone() {
                    if g.chromosome != chromosome {
                        continue;
                    }
                    if g.gene_start >= del_start && g.gene_end <= del_end {
                        let Some(t) = self.canonical(g) else { continue };
                        let exons = [1, t.exons.len() as i32 + 1];
                        self.disruptions.push(SvDisruptionData::new(
                            breakend.sv(),
                            breakend.uses_start(),
                            g,
                            t,
                            exons,
                            TranscriptCodingType::UNKNOWN,
                            TranscriptRegionType::UPSTREAM,
                            breakend.copy_number_low_side(),
                            breakend.copy_number(),
                        ));
                        self.disruptions.push(SvDisruptionData::new(
                            next.sv(),
                            next.uses_start(),
                            g,
                            t,
                            exons,
                            TranscriptCodingType::UNKNOWN,
                            TranscriptRegionType::DOWNSTREAM,
                            breakend.copy_number_low_side(),
                            breakend.copy_number(),
                        ));
                    }
                }
                i += 1;
            }
        }
    }

    fn check_sgl_mappings(&mut self, var: &Sv) {
        let bs = var.be(true);
        for mapping in var.sgl_mappings() {
            if mapping.chromosome != bs.chromosome() {
                continue;
            }
            if mapping.orientation == bs.orientation() {
                continue;
            }
            let (pos_start, pos_end, implied) = if mapping.position < bs.position() {
                (mapping.position, bs.position(), if mapping.orientation == 1 { DEL } else { DUP })
            } else {
                (bs.position(), mapping.position, if bs.forward_orient() { DEL } else { DUP })
            };
            if (mapping.position - bs.position()).abs() > MAX_SGL_MAPPED_LENGTH {
                continue;
            }
            for g in self.driver_gene_data.clone() {
                if g.chromosome != bs.chromosome() {
                    continue;
                }
                let Some(t) = self.canonical(g) else { continue };
                let mut disruptive = false;
                if implied == DEL && g.gene_start >= pos_start && g.gene_end <= pos_end {
                    disruptive = true;
                } else {
                    let within = |p: i32| p >= t.trans_start && p <= t.trans_end;
                    if !within(pos_start) && !within(pos_end) {
                        continue;
                    }
                    let mut rank_start = -1;
                    let mut rank_end = -1;
                    for i in 0..t.exons.len() {
                        let e = &t.exons[i];
                        if (pos_start >= e.start && pos_start <= e.end) || (pos_end >= e.start && pos_end <= e.end) {
                            disruptive = true;
                            break;
                        }
                        let next = t.exons.get(i + 1);
                        if let Some(n) = next {
                            if pos_start > e.end && pos_start < n.start {
                                rank_start = e.rank.min(n.rank);
                            }
                            if pos_end > e.end && pos_end < n.start {
                                rank_end = e.rank.min(n.rank);
                            }
                        }
                    }
                    if rank_start != rank_end {
                        if implied == DUP {
                            if (t.pos_strand() && rank_start >= 1) || (!t.pos_strand() && rank_end >= 1) {
                                disruptive = true;
                            }
                        } else {
                            disruptive = true;
                        }
                    }
                }
                if disruptive {
                    let d = SvDisruptionData::new(
                        bs.sv(),
                        bs.uses_start(),
                        g,
                        t,
                        [1, t.exons.len() as i32 + 1],
                        TranscriptCodingType::UNKNOWN,
                        TranscriptRegionType::UPSTREAM,
                        bs.copy_number_low_side(),
                        bs.copy_number(),
                    );
                    if implied == DEL && is_pseudogene_deletion(var, pos_start, pos_end, t) {
                        d.mark_pseudogene_deletion();
                    }
                    self.disruptions.push(d);
                    if !contains(&self.reportable_sgls, var) {
                        self.reportable_sgls.push(var.clone());
                    }
                }
            }
        }
    }

    /// writeGermlineSVs: the germline disruption, breakend and driver catalog files
    pub fn write_germline_svs(&self, standard: &[SvDisruptionData], sample_id: &str, output_dir: &str) -> anyhow::Result<()> {
        let mut germline_svs: Vec<String> = Vec::new();
        let mut breakends: Vec<LinxBreakend> = Vec::new();
        let mut drivers: Vec<DriverCatalog> = Vec::new();
        self.populate_germline_svs(standard, &mut germline_svs, &mut breakends, &mut drivers);
        let mut lines = vec![germline_disruption_header()];
        lines.extend(germline_svs);
        hmf_common::purple_files::write_lines(&format!("{output_dir}/{sample_id}.linx.germline.disruption.tsv"), &lines)?;
        hmf_common::purple_files::write_lines(&format!("{output_dir}/{sample_id}.linx.germline.breakend.tsv"), &breakend_lines(&breakends))?;
        write_driver_catalog(&format!("{output_dir}/{sample_id}.linx.germline.driver.catalog.tsv"), &drivers)?;
        Ok(())
    }

    pub fn populate_germline_svs(
        &self, standard: &[SvDisruptionData], germline_svs: &mut Vec<String>, breakends: &mut Vec<LinxBreakend>,
        drivers: &mut Vec<DriverCatalog>,
    ) {
        // a HashMap keyed by SV (identity) in Java (LX1)
        let mut sv_map: IdMap<Sv, Vec<&SvDisruptionData>> = IdMap::new();
        for d in self.disruptions.iter().chain(standard.iter()) {
            if let Some(l) = sv_map.get_mut(&d.var) {
                l.push(d);
            } else {
                sv_map.insert(d.var.clone(), vec![d]);
            }
        }
        let mut breakend_id = 0;
        for (var, list) in sv_map.iter() {
            let sv_data = var.sv_data();
            let cluster = var.cl();
            let pon_count = sv_data.pon_count;
            let mut all_filters = JavaHashSetStr::from_filter(&sv_data.filter);
            let mut gene_name = String::new();
            for d in list {
                let gene = d.gene;
                let transcript = d.transcript;
                let mut reportable = !self.disruptions.iter().any(|x| std::ptr::eq(x, *d)) && self.is_reportable(d);
                let mut reported_status = "NONE";
                if reportable {
                    if let Some(dg) = self.driver_genes.get(&d.gene.gene_name) {
                        if dg.report_germline_disruption == GermlineReporting::VariantNotLost && sv_data.junction_copy_number < 0.1 {
                            reportable = false;
                        } else if dg.report_germline_disruption == GermlineReporting::None {
                            reported_status = "NOT_REPORTED";
                        } else {
                            reported_status = "REPORTED";
                        }
                    }
                }
                if d.is_pseudogene_deletion() {
                    all_filters.remove("PASS");
                    all_filters.add(FILTER_PSEUDOGENE);
                }
                let is_upstream = (var.orientation(d.is_start) == 1) == gene.forward_strand();
                let mut b = LinxBreakend {
                    id: breakend_id,
                    sv_id: var.id(),
                    vcf_id: sv_data.vcf_id_start.clone(),
                    coords: var.coords_str(d.is_start),
                    is_start: d.is_start,
                    gene: gene.gene_name.clone(),
                    transcript_id: transcript.trans_name.clone(),
                    canonical: transcript.is_canonical,
                    gene_orientation: if is_upstream { BREAKEND_ORIENTATION_UPSTREAM } else { BREAKEND_ORIENTATION_DOWNSTREAM },
                    disruptive: false,
                    reported_status: reported_status.to_string(),
                    undisrupted_copy_number: d.undisrupted_copy_number,
                    region_type: d.region_type,
                    coding_type: d.coding_type,
                    biotype: transcript.bio_type.clone(),
                    exonic_base_phase: 0,
                    next_splice_exon_rank: d.exons[1],
                    next_splice_exon_phase: 0,
                    next_splice_distance: 0,
                    total_exon_count: transcript.exons.len() as i32,
                    exon_up: d.exons[0],
                    exon_down: d.exons[1],
                };
                breakend_id += 1;
                let bg = var.genes_list(d.is_start).into_iter().find(|x| x.gene_name() == d.gene.gene_name);
                if let Some(bt) = bg.and_then(|g| g.canonical()) {
                    b.region_type = bt.region_type();
                    b.coding_type = bt.coding_type();
                    b.disruptive = bt.is_disruptive();
                    b.exonic_base_phase = bt.phase;
                    b.next_splice_exon_rank = bt.next_splice_exon_rank();
                    b.next_splice_exon_phase = bt.phase;
                    b.next_splice_distance =
                        if bt.is_upstream() { bt.prev_splice_acceptor_distance() } else { bt.next_splice_acceptor_distance() };
                    b.exon_up = bt.exon_upstream;
                    b.exon_down = bt.exon_downstream;
                }
                breakends.push(b);
                if reportable && !drivers.iter().any(|x| x.gene == gene.gene_name) {
                    let dg = self.driver_genes.get(&d.gene.gene_name).expect("driver gene");
                    drivers.push(DriverCatalog {
                        chromosome: gene.chromosome.clone(),
                        chromosome_band: gene.karyotype_band.clone(),
                        gene: gene.gene_name.clone(),
                        transcript: d.transcript.trans_name.clone(),
                        is_canonical: d.transcript.is_canonical,
                        driver: "GERMLINE_DISRUPTION".into(),
                        category: dg.likelihood_type.name().into(),
                        likelihood_method: "GERMLINE".into(),
                        reported_status: reported_status.into(),
                        driver_likelihood: 1.0,
                        missense: 0,
                        nonsense: 0,
                        splice: 0,
                        inframe: 0,
                        frameshift: 0,
                        biallelic: false,
                        min_copy_number: 0.0,
                        max_copy_number: 0.0,
                    });
                }
                gene_name = gene.gene_name.clone();
            }
            germline_svs.push(germline_disruption_line(var, &all_filters.join(";"), &gene_name, &cluster, pon_count));
        }
    }

    fn is_reportable(&self, d: &SvDisruptionData) -> bool {
        let var = &d.var;
        if var.sv_data().filter != "PASS" {
            return false;
        }
        if var.sv_data().pon_count > MAX_PON_COUNT {
            return false;
        }
        let is_del_or_dup = var.sv_type() == DEL || var.sv_type() == DUP;
        if d.coding_type == TranscriptCodingType::UTR_3P && !is_del_or_dup {
            return false;
        }
        if d.is_pseudogene_deletion() {
            return false;
        }
        let cluster = var.cl();
        if cluster.sv_count() == 1 {
            if var.sv_type() == DEL && var.length() > MAX_DELETE_LENGTH {
                return false;
            }
        } else if cluster.metrics.borrow().total_deleted > MAX_DELETE_LENGTH as i64 {
            return false;
        }
        if contains(&self.reportable_sgls, var) {
            return true;
        }
        let rt = cluster.resolved_type();
        if matches!(rt, ResolvedType::DEL | ResolvedType::DUP | ResolvedType::INS | ResolvedType::RECIP_INV) {
            return true;
        }
        is_del_or_dup
    }
}

/// Sets.newHashSet(filter) of one string, iterated as Java's HashSet<String> does
struct JavaHashSetStr(hmf_common::java_hash::JavaHashSet<String>);
impl JavaHashSetStr {
    fn from_filter(filter: &str) -> Self {
        let mut s = hmf_common::java_hash::JavaHashSet::new();
        s.insert(filter.to_string());
        JavaHashSetStr(s)
    }
    fn remove(&mut self, v: &str) {
        self.0.remove(&v.to_string());
    }
    fn add(&mut self, v: &str) {
        self.0.insert(v.to_string());
    }
    fn join(&self, d: &str) -> String {
        self.0.iter().cloned().collect::<Vec<_>>().join(d)
    }
}

fn germline_disruption_header() -> String {
    crate::output::germline_disruption_header()
}

fn germline_disruption_line(var: &Sv, filters: &str, gene_name: &str, cluster: &Cluster, pon_count: i32) -> String {
    crate::output::germline_disruption_line(var, filters, gene_name, cluster, pon_count)
}
