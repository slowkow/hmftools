//! The `novel` package: alternative splice junctions (AltSpliceJunctionFinder, AltSpliceJunction, the cohort
//! frequencies) and retained introns (RetainedIntronFinder, RetainedIntron).

use crate::common::*;
use crate::compat;
use crate::gene::GeneCollection;
use crate::read::Read;
use crate::refgen::RefSource;
use anyhow::Result;
use hmf_common::known_fusion::{KnownFusionCache, KnownFusionType};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AsjContext {
    Unknown,
    SpliceJunc,
    Exonic,
    Intronic,
    Mixed,
}

impl AsjContext {
    pub fn name(self) -> &'static str {
        match self {
            AsjContext::Unknown => "UNKNOWN",
            AsjContext::SpliceJunc => "SPLICE_JUNC",
            AsjContext::Exonic => "EXONIC",
            AsjContext::Intronic => "INTRONIC",
            AsjContext::Mixed => "MIXED",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AsjType {
    SkippedExons,
    Novel5Prime,
    Novel3Prime,
    NovelExon,
    NovelIntron,
    MixedTrans,
    Intronic,
    ExonIntron,
    Circular,
}

impl AsjType {
    pub fn name(self) -> &'static str {
        match self {
            AsjType::SkippedExons => "SKIPPED_EXONS",
            AsjType::Novel5Prime => "NOVEL_5_PRIME",
            AsjType::Novel3Prime => "NOVEL_3_PRIME",
            AsjType::NovelExon => "NOVEL_EXON",
            AsjType::NovelIntron => "NOVEL_INTRON",
            AsjType::MixedTrans => "MIXED_TRANS",
            AsjType::Intronic => "INTRONIC",
            AsjType::ExonIntron => "EXON_INTRON",
            AsjType::Circular => "CIRCULAR",
        }
    }
}

/// AltSjCohortCache
pub struct AltSjCohortCache {
    freq: HashMap<String, i32>,
}

impl AltSjCohortCache {
    pub fn load(filename: &str) -> Result<AltSjCohortCache> {
        let mut c = AltSjCohortCache { freq: HashMap::new() };
        if !std::path::Path::new(filename).exists() {
            hmf_common::error!("missing alt-SJ cohort file");
            return Ok(c);
        }
        let delim = if filename.ends_with(".tsv") || filename.ends_with(".tsv.gz") { '\t' } else { ',' };
        let mut reader = hmf_common::gz::LineReader::open(filename)?;
        let header = match reader.next_line()? {
            Some(h) => std::str::from_utf8(h)?.to_string(),
            None => return Ok(c),
        };
        let f = hmf_common::java_text::fields_index_map(&header, delim);
        let get = |n: &str| f.get(n).copied().ok_or_else(|| anyhow::anyhow!("alt-SJ cohort file has no {n} column"));
        let (ci, si, ei, ni) = (get("Chromosome")?, get("SjStart")?, get("SjEnd")?, get("SampleCount")?);
        while let Some(line) = reader.next_line()? {
            let line = std::str::from_utf8(line)?;
            let v: Vec<&str> = hmf_common::java_text::split_all(line, delim);
            let n: i32 = v[ni].parse()?;
            let s: i32 = v[si].parse()?;
            let e: i32 = v[ei].parse()?;
            c.freq.insert(format!("{};{};{}", v[ci], s, e), n);
        }
        hmf_common::info!("loaded {} alt-SJ sites from cohort file({})", c.freq.len(), filename);
        Ok(c)
    }
    pub fn frequency(&self, chr: &str, s: i32, e: i32) -> i32 {
        self.freq.get(&format!("{chr};{s};{e}")).copied().unwrap_or(0)
    }
}

/// AltSpliceJunction
#[derive(Clone, Debug)]
pub struct AltSpliceJunction {
    pub chromosome: &'static str,
    pub sj: [i32; 2],
    pub sj_start_regions: Vec<usize>,
    pub sj_end_regions: Vec<usize>,
    pub contexts: [AsjContext; 2],
    pub sj_type: AsjType,
    pub fragment_count: i32,
    pub position_counts: [i32; 2],
    pub candidate_trans_ids: Vec<i32>,
    pub gene_id: Option<&'static str>,
    pub gene_name: Option<&'static str>,
    pub selected_transcripts: [Option<&'static str>; 2],
    pub selected_exons: [i32; 2],
    pub transcript_names: [String; 2],
    pub base_context: [String; 2],
    pub donor_acceptor_bases: String,
    pub nearest_exon_distance: [i32; 2],
    pub cohort_frequency: i32,
    pub initial_read_id: String,
    pub not_known: bool,
}

impl AltSpliceJunction {
    fn new(chromosome: &'static str, sj: [i32; 2], t: AsjType, read_id: &str, contexts: [AsjContext; 2], starts: Vec<usize>, ends: Vec<usize>) -> Self {
        AltSpliceJunction {
            chromosome,
            sj,
            sj_start_regions: starts,
            sj_end_regions: ends,
            contexts,
            sj_type: t,
            fragment_count: 0,
            position_counts: [0, 0],
            candidate_trans_ids: Vec::new(),
            gene_id: None,
            gene_name: None,
            selected_transcripts: [None, None],
            selected_exons: [-1, -1],
            transcript_names: [String::new(), String::new()],
            base_context: [String::new(), String::new()],
            donor_acceptor_bases: String::new(),
            nearest_exon_distance: [0, 0],
            cohort_frequency: -1,
            initial_read_id: read_id.to_string(),
            not_known: false,
        }
    }
    fn matches(&self, o: &AltSpliceJunction) -> bool {
        self.chromosome == o.chromosome && self.sj == o.sj
    }

    fn set_candidate_transcripts(&mut self, gc: &GeneCollection, candidate_regions: &[usize]) {
        let regions = &gc.exon_regions;
        let mut valid: Vec<i32> = Vec::new();
        if !self.sj_start_regions.is_empty() || !self.sj_end_regions.is_empty() {
            for se in [SE_START, SE_END] {
                let l = if se == SE_START { &self.sj_start_regions } else { &self.sj_end_regions };
                for r in l {
                    for t in &regions[*r].trans_exon_refs {
                        valid.push(t.trans_id);
                    }
                }
            }
        } else {
            for &ri in candidate_regions {
                let r = &regions[ri];
                if position_within(self.sj[0], r.start, r.end) && position_within(self.sj[1], r.start, r.end) {
                    valid.extend(r.trans_exon_refs.iter().map(|x| x.trans_id));
                    continue;
                }
                if position_within(self.sj[0], r.start, r.end) {
                    for x in &r.trans_exon_refs {
                        if r.post_regions.iter().any(|p| regions[*p].has_trans_id(x.trans_id)) {
                            valid.push(x.trans_id);
                        }
                    }
                }
                if position_within(self.sj[1], r.start, r.end) {
                    for x in &r.trans_exon_refs {
                        if r.pre_regions.iter().any(|p| regions[*p].has_trans_id(x.trans_id)) {
                            valid.push(x.trans_id);
                        }
                    }
                }
            }
        }
        for t in valid {
            if !self.candidate_trans_ids.contains(&t) {
                self.candidate_trans_ids.push(t);
            }
        }
    }

    fn cull_non_matched_transcripts(&mut self, gc: &GeneCollection, valid: &[i32]) {
        let regions = &gc.exon_regions;
        self.sj_start_regions.retain(|r| valid.iter().any(|t| regions[*r].has_trans_id(*t)));
        self.sj_end_regions.retain(|r| valid.iter().any(|t| regions[*r].has_trans_id(*t)));
        self.candidate_trans_ids.retain(|t| valid.contains(t));
    }

    fn set_base_context(&mut self, refg: &mut RefSource, chromosome: &str) {
        for se in [SE_START, SE_END] {
            let start_in_exon = self.sj_type != AsjType::Circular;
            let start_offset = if (se == SE_START) == start_in_exon { 1 } else { 10 };
            let end_offset = if start_offset == 1 { 10 } else { 1 };
            self.base_context[se] = refg.base_string(chromosome, self.sj[se] - start_offset, self.sj[se] + end_offset);
        }
        // IF6: Java never sets the donor-acceptor bases, so the splice-motif strand never guides the choice of gene
        if !compat::java("IF6") {
            self.donor_acceptor_bases = donor_acceptor_bases(&self.base_context);
        }
    }

    fn known_splice_base_strand(&self) -> i8 {
        match self.donor_acceptor_bases.as_str() {
            "GT-AG" | "GC-AG" => 1,
            "CT-AC" | "CT-GC" => -1,
            _ => 0,
        }
    }

    fn set_transcript_and_exon(&mut self, gc: &GeneCollection, se: usize, regions: &[usize]) {
        if regions.is_empty() {
            return;
        }
        for r in regions {
            for t in &gc.exon_regions[*r].trans_exon_refs {
                if t.is_canonical() {
                    self.selected_transcripts[se] = Some(t.trans_name);
                    self.selected_exons[se] = t.exon_rank;
                    return;
                }
            }
        }
        if let Some(t) = gc.exon_regions[regions[0]].trans_exon_refs.first() {
            self.selected_transcripts[se] = Some(t.trans_name);
            self.selected_exons[se] = t.exon_rank;
        }
    }

    pub(crate) fn calc_nearest_exon_boundary(&mut self, gc: &GeneCollection, se: usize, gene_idx: usize) -> i32 {
        if (se == SE_START && !self.sj_start_regions.is_empty()) || (se == SE_END && !self.sj_end_regions.is_empty()) {
            let regions = if se == SE_START { self.sj_start_regions.clone() } else { self.sj_end_regions.clone() };
            self.set_transcript_and_exon(gc, se, &regions);
            return 0;
        }
        let gene = &gc.genes[gene_idx];
        let pos = self.sj[se];
        let forward = gene.gene.strand == 1;
        let five_prime = (se == SE_START) == forward;
        let exonic = self.contexts[se] == AsjContext::Exonic || self.contexts[se] == AsjContext::Mixed;
        let search_forwards = (five_prime && exonic) || (!five_prime && !exonic);
        let mut nearest = 0;
        let mut nearest_region: Option<usize> = None;
        for &ri in &gene.exon_regions {
            let r = &gc.exon_regions[ri];
            if exonic && !position_within(pos, r.start, r.end) {
                continue;
            }
            if position_within(pos, r.start, r.end) {
                if !r.trans_exon_refs.iter().any(|x| self.candidate_trans_ids.contains(&x.trans_id)) {
                    continue;
                }
                let d = if search_forwards { pos - r.end } else { r.start - pos };
                if nearest == 0 || d > nearest {
                    nearest = d;
                    nearest_region = Some(ri);
                }
            } else {
                if exonic {
                    continue;
                }
                if (search_forwards && r.start < pos) || (!search_forwards && pos < r.end) {
                    continue;
                }
                let d = if search_forwards { r.start - pos } else { pos - r.end };
                if nearest == 0 || d < nearest {
                    nearest = d;
                    nearest_region = Some(ri);
                }
            }
        }
        if let Some(r) = nearest_region {
            self.set_transcript_and_exon(gc, se, &[r]);
        }
        nearest
    }

    fn generate_transcript_names(&self, gc: &GeneCollection, regions: &[usize]) -> String {
        let mut names: Vec<&str> = Vec::new();
        for r in regions {
            for t in &gc.exon_regions[*r].trans_exon_refs {
                if self.candidate_trans_ids.contains(&t.trans_id) {
                    names.push(t.trans_name);
                }
            }
        }
        names.join(";")
    }

    fn calc_summary_data(&mut self, gc: &GeneCollection, gene_idx: Option<usize>) {
        for se in [SE_START, SE_END] {
            let regions = if se == SE_START { self.sj_start_regions.clone() } else { self.sj_end_regions.clone() };
            self.transcript_names[se] = if regions.is_empty() { "NONE".to_string() } else { self.generate_transcript_names(gc, &regions) };
            // Java dereferences the gene, which is always found
            self.nearest_exon_distance[se] = self.calc_nearest_exon_boundary(gc, se, gene_idx.unwrap_or(0));
        }
    }

    /// NovelSpliceJunctionFile.write
    pub fn novel_line(&self) -> String {
        [
            self.gene_name.unwrap_or("null").to_string(),
            self.chromosome.to_string(),
            self.sj[0].to_string(),
            self.sj[1].to_string(),
            self.sj_type.name().to_string(),
            self.selected_transcripts[0].unwrap_or("").to_string(),
            self.selected_transcripts[1].unwrap_or("").to_string(),
            self.selected_exons[0].to_string(),
            self.selected_exons[1].to_string(),
            self.fragment_count.to_string(),
            self.position_counts[0].to_string(),
            self.position_counts[1].to_string(),
            self.contexts[0].name().to_string(),
            self.contexts[1].name().to_string(),
            self.base_context[0].clone(),
            self.base_context[1].clone(),
            self.cohort_frequency.to_string(),
        ]
        .join("\t")
    }

    /// AltSpliceJunction.write: the unfiltered file's line
    pub fn unfiltered_line(&self) -> String {
        [
            self.gene_id.unwrap_or("null").to_string(),
            self.novel_line(),
            if self.not_known { "NOT_KNOWN".to_string() } else { "PASS".to_string() },
            self.transcript_names[0].clone(),
            self.transcript_names[1].clone(),
            self.nearest_exon_distance[0].to_string(),
            self.nearest_exon_distance[1].to_string(),
            self.initial_read_id.clone(),
        ]
        .join("\t")
    }
}

pub fn donor_acceptor_bases(ctx: &[String; 2]) -> String {
    let s = if ctx[0].len() == 12 { &ctx[0][2..4] } else { "" };
    let e = if ctx[1].len() == 12 { &ctx[1][8..10] } else { "" };
    format!("{s}-{e}")
}

pub const ALT_SJ_HEADER: &str = "GeneName\tChromosome\tSjStart\tSjEnd\tType\tTranscriptStart\tTranscriptEnd\tExonStart\tExonEnd\tFragCount\tDepthStart\tDepthEnd\tRegionStart\tRegionEnd\tBaseStart\tBaseEnd\tCohortFrequency";

pub fn alt_sj_unfiltered_header() -> String {
    format!("GeneId\t{ALT_SJ_HEADER}\tFilter\tTransStart\tTransEnd\tNearestStartExon\tNearestEndExon\tInitialReadId")
}

const MAX_HOMOLOGY_LENGTH: i32 = 5;

/// AltSpliceJunctionFinder: the novel junctions of one gene collection
pub struct AltSpliceJunctionFinder<'a> {
    pub enabled: bool,
    pub junctions: Vec<AltSpliceJunction>,
    cohort: Option<&'a AltSjCohortCache>,
    known_exon_del_dup_genes: HashSet<String>,
    known_fusions: &'a KnownFusionCache,
}

impl<'a> AltSpliceJunctionFinder<'a> {
    pub fn new(enabled: bool, cohort: Option<&'a AltSjCohortCache>, known_fusions: &'a KnownFusionCache) -> Self {
        let known_exon_del_dup_genes = known_fusions.data_by_type(KnownFusionType::EXON_DEL_DUP).map(|x| x.five_gene.clone()).collect();
        AltSpliceJunctionFinder { enabled, junctions: Vec::new(), cohort, known_exon_del_dup_genes, known_fusions }
    }

    pub fn clear(&mut self) {
        self.junctions.clear();
    }

    /// evaluateFragmentReads
    pub fn evaluate_fragment_reads(&mut self, gc: &GeneCollection, refg: &mut RefSource, genes: &[usize], reads: &[&Read], related: &[i32]) {
        if !self.enabled {
            return;
        }
        if genes.is_empty() || reads.iter().any(|r| r.is_duplicate() || r.is_multi_mapped()) {
            return;
        }
        let min_start = genes.iter().map(|g| gc.genes[*g].gene.gene_start).min().unwrap_or(0);
        // IF3: Java bounds the upper side by the genes' latest start, not their end, so a junction more than
        // 500 kb past the start of a long gene is dropped
        let max_upper = if compat::java("IF3") {
            genes.iter().map(|g| gc.genes[*g].gene.gene_start).max().unwrap_or(0)
        } else {
            genes.iter().map(|g| gc.genes[*g].gene.gene_end).max().unwrap_or(0)
        };
        let bounds = [min_start - crate::config::MAX_NOVEL_SJ_DISTANCE, max_upper + crate::config::MAX_NOVEL_SJ_DISTANCE];
        if reads.iter().any(|r| !positions_within(r.pos_start, r.pos_end, bounds[0], bounds[1])) {
            return;
        }
        let candidate_genes: Vec<usize> = genes
            .iter()
            .copied()
            .filter(|g| {
                let gd = gc.genes[*g].gene;
                reads.iter().any(|r| positions_within(r.pos_start, r.pos_end, gd.gene_start, gd.gene_end))
            })
            .collect();
        if candidate_genes.is_empty() {
            return;
        }
        if reads.len() == 1 {
            if is_candidate(reads[0]) {
                self.register_from_read(gc, refg, &candidate_genes, reads[0], related);
            }
            return;
        }
        let (r1, r2) = (reads[0], reads[1]);
        if r1.has_supp_alignment() && r2.has_supp_alignment() {
            self.register_from_reads(gc, &candidate_genes, r1, r2, related);
        } else {
            let mut first: Option<usize> = None;
            if is_candidate(r1) {
                first = self.register_from_read(gc, refg, &candidate_genes, r1, related);
            }
            let mut second: Option<usize> = None;
            if is_candidate(r2) {
                if first.is_some() && positions_overlap(r1.pos_start, r1.pos_end, r2.pos_start, r2.pos_end) {
                    return;
                }
                second = self.register_from_read(gc, refg, &candidate_genes, r2, related);
            }
            if let (Some(a), Some(b)) = (first, second) {
                self.check_novel_exon(gc, a, b);
            }
        }
    }

    pub(crate) fn junction_matches_known(gc: &GeneCollection, genes: &[usize], sj: [i32; 2], trans_ids: &[i32]) -> bool {
        for g in genes {
            for t in &gc.genes[*g].transcripts {
                if !trans_ids.contains(&t.trans_id) {
                    continue;
                }
                for i in 0..t.exons.len().saturating_sub(1) {
                    if t.exons[i].end == sj[0] && t.exons[i + 1].start == sj[1] {
                        return true;
                    }
                }
            }
        }
        false
    }

    pub(crate) fn create_from_read(&self, gc: &GeneCollection, refg: &mut RefSource, read: &Read, related: &[i32]) -> Option<AltSpliceJunction> {
        let mut sj = [0, 0];
        let coords = &read.mapped_coords.alignments;
        if read.mapped_coords.inferred(SE_START) {
            sj = [coords[1].end, coords[2].start];
        } else {
            let entries = read.mapped_regions.entries();
            for i in 0..coords.len().saturating_sub(1) {
                let (js, je) = (coords[i].end, coords[i + 1].start);
                if matches_known_splice_site(gc, &entries, js, je) {
                    continue;
                }
                sj = [js, je];
                break;
            }
            if sj[0] == 0 || sj[1] == 0 {
                return None;
            }
        }
        self.check_junction_homology(gc, refg, &mut sj, read);
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut ctx = [AsjContext::Unknown, AsjContext::Unknown];
        classify_regions(gc, read, sj, &mut starts, &mut ends, &mut ctx);
        let t = classify_splice_junction(gc, related, &starts, &ends, ctx);
        let mut asj = AltSpliceJunction::new(gc.chromosome, sj, t, &read.id, ctx, starts, ends);
        let keys = read.mapped_regions.keys();
        asj.set_candidate_transcripts(gc, &keys);
        Some(asj)
    }

    fn check_junction_homology(&self, gc: &GeneCollection, refg: &mut RefSource, sj: &mut [i32; 2], read: &Read) {
        let regions = &gc.exon_regions;
        let keys = read.mapped_regions.keys();
        let mut matched: Option<Vec<usize>> = None;
        let mut matches_both = false;
        'outer: for se in [SE_START, SE_END] {
            let p = sj[se];
            for &ri in &keys {
                let r = &regions[ri];
                let mb = (se == SE_START && p == r.end) || (se == SE_END && p == r.start);
                if mb {
                    if se == SE_START {
                        matched.get_or_insert_with(Vec::new).push(ri);
                    } else if matched.as_ref().map(|m| m.contains(&ri)).unwrap_or(false) {
                        matches_both = true;
                        break 'outer;
                    }
                }
            }
        }
        if matches_both {
            return;
        }
        let mut bases = [String::new(), String::new()];
        let mut nearest_exon = [0, 0];
        for se in [SE_START, SE_END] {
            bases[se] = refg.base_string(gc.chromosome, sj[se] - MAX_HOMOLOGY_LENGTH, sj[se] + MAX_HOMOLOGY_LENGTH);
            let mut nearest = MAX_HOMOLOGY_LENGTH + 1;
            for &ri in &keys {
                let r = &regions[ri];
                let d = if (sj[se] - r.start).abs() < (sj[se] - r.end).abs() { sj[se] - r.start } else { sj[se] - r.end };
                if d.abs() < nearest.abs() {
                    nearest = d;
                }
            }
            nearest_exon[se] = nearest;
        }
        if nearest_exon[0].abs() > MAX_HOMOLOGY_LENGTH && nearest_exon[1].abs() > MAX_HOMOLOGY_LENGTH {
            return;
        }
        let min_distance = nearest_exon[0].min(nearest_exon[1]);
        let jb = MAX_HOMOLOGY_LENGTH as usize;
        let d = min_distance.unsigned_abs() as usize;
        // IF5: Java compares the start's bases with themselves, so any junction within 5 bases of an exon
        // boundary is shifted to it; corrected, the bases moved must match at both ends
        let other = if compat::java("IF5") { 0 } else { 1 };
        let sub = |s: &str, a: usize, b: usize| -> String { s.get(a..b).unwrap_or("").to_string() };
        if min_distance > 0 {
            let (start_bases, end_bases) = if other == 0 {
                (sub(&bases[0], jb - d, jb), sub(&bases[0], jb - d, jb))
            } else {
                // the last d bases of the first section move to the end of the junction's gap
                (sub(&bases[0], jb + 1 - d, jb + 1), sub(&bases[1], jb - d, jb))
            };
            if start_bases != end_bases {
                return;
            }
        } else if min_distance < 0 {
            let (start_bases, end_bases) = if other == 0 {
                (sub(&bases[0], jb, jb + d), sub(&bases[0], jb, jb + d))
            } else {
                // the first d bases of the second section move to the start of the gap
                (sub(&bases[0], jb + 1, jb + 1 + d), sub(&bases[1], jb, jb + d))
            };
            if start_bases != end_bases {
                return;
            }
        }
        sj[0] -= min_distance;
        sj[1] -= min_distance;
    }

    pub(crate) fn create_from_reads(&self, gc: &GeneCollection, r1: &Read, r2: &Read, related: &[i32]) -> Option<AltSpliceJunction> {
        let mut sj = [0, 0];
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        let mut candidates: Vec<usize> = Vec::new();
        for r in [r1, r2] {
            let use_left = r.longest_soft_clip_is_left() == Some(true);
            if use_left {
                sj[SE_START] = r.coords_boundary(SE_START);
            } else {
                sj[SE_END] = r.coords_boundary(SE_END);
            }
            for k in r.mapped_regions.keys() {
                if !candidates.contains(&k) {
                    candidates.push(k);
                }
            }
        }
        if sj[0] == 0 || sj[1] == 0 {
            return None;
        }
        let mut asj = if sj[0] < sj[1] {
            let flipped = [sj[1], sj[0]];
            let mut fctx = [AsjContext::Unknown, AsjContext::Unknown];
            classify_regions(gc, r1, flipped, &mut ends, &mut starts, &mut fctx);
            classify_regions(gc, r2, flipped, &mut ends, &mut starts, &mut fctx);
            let ctx = [fctx[1], fctx[0]];
            AltSpliceJunction::new(gc.chromosome, sj, AsjType::Circular, &r1.id, ctx, starts, ends)
        } else {
            sj = [sj[1], sj[0]];
            let mut ctx = [AsjContext::Unknown, AsjContext::Unknown];
            classify_regions(gc, r1, sj, &mut ends, &mut starts, &mut ctx);
            classify_regions(gc, r2, sj, &mut ends, &mut starts, &mut ctx);
            let t = classify_splice_junction(gc, related, &starts, &ends, ctx);
            AltSpliceJunction::new(gc.chromosome, sj, t, &r1.id, ctx, starts, ends)
        };
        asj.set_candidate_transcripts(gc, &candidates);
        Some(asj)
    }

    pub(crate) fn check_novel_exon(&mut self, gc: &GeneCollection, a: usize, b: usize) {
        let (ta, tb) = (self.junctions[a].sj_type, self.junctions[b].sj_type);
        if !(ta == AsjType::Novel3Prime && tb == AsjType::Novel5Prime) && !(ta == AsjType::Novel5Prime && tb == AsjType::Novel3Prime) {
            return;
        }
        let mut r1 = self.junctions[a].sj_start_regions.clone();
        r1.extend(self.junctions[a].sj_end_regions.iter().copied());
        let mut r2 = self.junctions[b].sj_start_regions.clone();
        r2.extend(self.junctions[b].sj_end_regions.iter().copied());
        let regions = &gc.exon_regions;
        let mut common: Vec<i32> = Vec::new();
        for x in &r1 {
            for t1 in regions[*x].trans_exon_refs.iter().map(|t| t.trans_id) {
                for y in &r2 {
                    if regions[*y].trans_exon_refs.iter().any(|t| t.trans_id == t1) {
                        if !common.contains(&t1) {
                            common.push(t1);
                        }
                        break;
                    }
                }
            }
        }
        if common.is_empty() {
            return;
        }
        self.junctions[a].sj_type = AsjType::NovelExon;
        self.junctions[b].sj_type = AsjType::NovelExon;
        self.junctions[a].cull_non_matched_transcripts(gc, &common);
        self.junctions[b].cull_non_matched_transcripts(gc, &common);
    }

    fn register(&mut self, gc: &GeneCollection, genes: &[usize], asj: AltSpliceJunction, related: &[i32]) -> Option<usize> {
        if let Some(i) = self.junctions.iter().position(|x| x.matches(&asj)) {
            self.junctions[i].fragment_count += 1;
            return Some(i);
        }
        if Self::junction_matches_known(gc, genes, asj.sj, related) {
            return None;
        }
        let mut asj = asj;
        asj.fragment_count += 1;
        self.junctions.push(asj);
        Some(self.junctions.len() - 1)
    }

    fn register_from_read(&mut self, gc: &GeneCollection, refg: &mut RefSource, genes: &[usize], read: &Read, related: &[i32]) -> Option<usize> {
        let asj = self.create_from_read(gc, refg, read, related)?;
        self.register(gc, genes, asj, related)
    }

    fn register_from_reads(&mut self, gc: &GeneCollection, genes: &[usize], r1: &Read, r2: &Read, related: &[i32]) -> Option<usize> {
        let asj = self.create_from_reads(gc, r1, r2, related)?;
        self.register(gc, genes, asj, related)
    }

    pub fn set_position_depth(&mut self, depth: &BaseDepth) {
        for a in &mut self.junctions {
            for se in [SE_START, SE_END] {
                a.position_counts[se] += depth.depth_at_base(a.sj[se]);
            }
        }
    }

    pub fn prioritise_genes(&mut self, gc: &GeneCollection, refg: &mut RefSource) {
        for a in &mut self.junctions {
            a.set_base_context(refg, gc.chromosome);
        }
        for a in &mut self.junctions {
            let strand = a.known_splice_base_strand();
            let mut candidates: Vec<usize> = (0..gc.genes.len())
                .filter(|g| gc.genes[*g].exon_regions.iter().any(|r| a.sj_start_regions.contains(r) || a.sj_end_regions.contains(r)))
                .collect();
            if candidates.is_empty() {
                candidates = (0..gc.genes.len()).collect();
            }
            let mut top: Option<usize> = None;
            let mut top_match = 0;
            for g in candidates {
                if strand != 0 && gc.genes[g].gene.strand != strand {
                    continue;
                }
                let m = gc.genes[g].transcripts.iter().filter(|t| a.candidate_trans_ids.contains(&t.trans_id)).count();
                if m > top_match {
                    top_match = m;
                    top = Some(g);
                }
            }
            let g = gc.genes[top.unwrap_or(0)].gene;
            a.gene_id = Some(&g.gene_id);
            a.gene_name = Some(&g.gene_name);
        }
    }

    pub fn finalise(&mut self, gc: &GeneCollection) {
        for i in 0..self.junctions.len() {
            let gid = self.junctions[i].gene_id;
            let gene_idx = gc.genes.iter().position(|g| Some(g.gene.gene_id.as_str()) == gid);
            self.junctions[i].calc_summary_data(gc, gene_idx);
            if let Some(c) = self.cohort {
                let a = &self.junctions[i];
                let f = c.frequency(a.chromosome, a.sj[0], a.sj[1]);
                self.junctions[i].cohort_frequency = f;
            }
            if !self.matches_known_exon_del_dup(&self.junctions[i]) {
                self.junctions[i].not_known = true;
            }
        }
    }

    fn matches_known_exon_del_dup(&self, a: &AltSpliceJunction) -> bool {
        if !a.gene_name.map(|g| self.known_exon_del_dup_genes.contains(g)).unwrap_or(false) {
            return false;
        }
        if a.selected_transcripts[0].is_none() || a.selected_transcripts[0] != a.selected_transcripts[1] {
            return false;
        }
        if a.selected_exons[0] < 0 || a.selected_exons[1] < 0 {
            return false;
        }
        let trans = a.selected_transcripts[0].unwrap();
        let mut up = a.selected_exons[0];
        let mut down = a.selected_exons[1];
        for k in self.known_fusions.data_by_type(KnownFusionType::EXON_DEL_DUP) {
            if k.specific_exons_trans_name() != trans {
                continue;
            }
            let five = k.five_gene_exon_range();
            let three = k.three_gene_exon_range();
            let is_dup = five[0] > three[0];
            if is_dup {
                if a.sj_type != AsjType::Circular {
                    continue;
                }
                std::mem::swap(&mut up, &mut down);
            } else if a.sj_type != AsjType::SkippedExons {
                continue;
            }
            if up >= five[0] && up <= five[1] && down >= three[0] && down <= three[1] {
                return true;
            }
        }
        false
    }

    /// writeAltSpliceJunctions: the unfiltered and passing lines
    pub fn write(&self, unfiltered: &mut Vec<String>, passing: &mut Vec<String>) {
        for a in &self.junctions {
            unfiltered.push(a.unfiltered_line());
            if !a.not_known {
                passing.push(a.novel_line());
            }
        }
    }
}

fn is_candidate(read: &Read) -> bool {
    if !read.has_split {
        return false;
    }
    if read.trans_class.values_unordered().any(|t| *t == TransMatchType::SpliceJunction) {
        return false;
    }
    read.mapped_coords.original_alignment_count() != 1
}

fn matches_known_splice_site(gc: &GeneCollection, entries: &[(usize, RegionMatchType)], donor: i32, acceptor: i32) -> bool {
    let regions = &gc.exon_regions;
    let mut donors: Vec<usize> = Vec::new();
    let mut acceptors: Vec<usize> = Vec::new();
    for (ri, mt) in entries {
        if *mt != RegionMatchType::ExonBoundary {
            continue;
        }
        if donor == regions[*ri].end {
            donors.push(*ri);
        } else if acceptor == regions[*ri].start {
            acceptors.push(*ri);
        }
    }
    if donors.is_empty() || acceptors.is_empty() {
        return false;
    }
    for d in &donors {
        for dr in &regions[*d].trans_exon_refs {
            for a in &acceptors {
                for ar in &regions[*a].trans_exon_refs {
                    if dr.trans_id == ar.trans_id && (dr.exon_rank - ar.exon_rank).abs() == 1 {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn classify_regions(gc: &GeneCollection, read: &Read, sj: [i32; 2], starts: &mut Vec<usize>, ends: &mut Vec<usize>, ctx: &mut [AsjContext; 2]) {
    let regions = &gc.exon_regions;
    let entries = read.mapped_regions.entries();
    let mut matched_trans: Vec<i32> = Vec::new();
    for (ri, mt) in &entries {
        if *mt == RegionMatchType::None {
            continue;
        }
        let r = &regions[*ri];
        if r.end == sj[SE_START] {
            if r.post_regions.is_empty() {
                continue;
            }
            ctx[SE_START] = AsjContext::SpliceJunc;
            starts.push(*ri);
            matched_trans.extend(r.trans_exon_refs.iter().map(|x| x.trans_id));
        }
        if r.start == sj[SE_END] {
            if r.pre_regions.is_empty() {
                continue;
            }
            ctx[SE_END] = AsjContext::SpliceJunc;
            ends.push(*ri);
            matched_trans.extend(r.trans_exon_refs.iter().map(|x| x.trans_id));
        }
    }
    for se in [SE_START, SE_END] {
        if ctx[se] == AsjContext::SpliceJunc {
            continue;
        }
        for (ri, mt) in &entries {
            if *mt == RegionMatchType::None {
                continue;
            }
            let r = &regions[*ri];
            if !matched_trans.is_empty() && !r.trans_exon_refs.iter().any(|x| matched_trans.contains(&x.trans_id)) {
                continue;
            }
            if position_within(sj[se], r.start, r.end) {
                ctx[se] = AsjContext::Exonic;
            }
        }
    }
    for c in ctx.iter_mut() {
        if *c == AsjContext::Unknown {
            *c = AsjContext::Intronic;
        }
    }
}

fn classify_splice_junction(gc: &GeneCollection, trans_ids: &[i32], starts: &[usize], ends: &[usize], ctx: [AsjContext; 2]) -> AsjType {
    let regions = &gc.exon_regions;
    if !starts.is_empty() || !ends.is_empty() {
        let mut has_start = false;
        let mut has_end = false;
        let mut start_trans = -1;
        let mut end_trans = -1;
        for t in trans_ids {
            let ms = starts.iter().any(|r| regions[*r].has_trans_id(*t));
            let me = ends.iter().any(|r| regions[*r].has_trans_id(*t));
            if ms && me {
                return AsjType::SkippedExons;
            }
            if ms {
                has_start = true;
                start_trans = *t;
            } else if me {
                has_end = true;
                end_trans = *t;
            }
        }
        if has_start && has_end {
            return AsjType::MixedTrans;
        }
        if has_start {
            return if gc.strand(start_trans) == 1 { AsjType::Novel3Prime } else { AsjType::Novel5Prime };
        } else if has_end {
            return if gc.strand(end_trans) == 1 { AsjType::Novel5Prime } else { AsjType::Novel3Prime };
        }
    }
    if ctx[0] == AsjContext::Intronic && ctx[1] == AsjContext::Intronic {
        return AsjType::Intronic;
    }
    if ctx[0] == AsjContext::Exonic && ctx[1] == AsjContext::Exonic {
        return AsjType::NovelIntron;
    }
    AsjType::ExonIntron
}

// ---- retained introns ----

const RI_MIN_FRAG_COUNT: i32 = 3;
const RI_MIN_SPLICED_FRAG_COUNT: i32 = 1;

#[derive(Clone, Debug)]
pub struct RetainedIntron {
    pub regions: Vec<usize>,
    pub is_start: bool,
    pub position: i32,
    pub fragment_count: i32,
    pub spliced_fragment_count: i32,
    pub depth: i32,
    pub transcript_info: String,
}

impl RetainedIntron {
    fn new(gc: &GeneCollection, regions: Vec<usize>, is_start: bool) -> Self {
        let position = if is_start { gc.exon_regions[regions[0]].start } else { gc.exon_regions[regions[0]].end };
        let mut info: Vec<String> = Vec::new();
        for r in &regions {
            for t in &gc.exon_regions[*r].trans_exon_refs {
                info.push(format!("{}-{}", t.trans_name, t.exon_rank));
            }
        }
        RetainedIntron { regions, is_start, position, fragment_count: 0, spliced_fragment_count: 0, depth: 0, transcript_info: info.join(";") }
    }
    fn matches(&self, o: &RetainedIntron) -> bool {
        o.position == self.position && o.is_start == self.is_start
    }
    fn add_fragment_count(&mut self, spliced: bool) {
        self.fragment_count += 1;
        if spliced {
            self.spliced_fragment_count += 1;
        }
    }
    pub fn type_name(&self, forward: bool) -> &'static str {
        if forward == self.is_start { "ACCEPTOR" } else { "DONOR" }
    }
}

pub const RETAINED_INTRON_HEADER: &str = "GeneId\tGeneName\tChromosome\tStrand\tPosition\tType\tFragCount\tSplicedFragCount\tTotalDepth\tTranscriptInfo";

#[derive(Default)]
pub struct RetainedIntronFinder {
    pub enabled: bool,
    pub introns: Vec<RetainedIntron>,
}

impl RetainedIntronFinder {
    pub fn evaluate_fragment_reads(&mut self, gc: &GeneCollection, reads: &[&Read]) {
        if !self.enabled {
            return;
        }
        if reads.iter().any(|r| r.is_duplicate() || r.is_multi_mapped()) {
            return;
        }
        let mut spliced_trans: Vec<i32> = Vec::new();
        let mut found: Vec<RetainedIntron> = Vec::new();
        for r in reads {
            if r.has_split {
                for (t, m) in r.trans_class.iter() {
                    if *m == TransMatchType::SpliceJunction {
                        spliced_trans.push(t);
                    }
                }
                continue;
            }
            if let Some(ri) = evaluate_read(gc, r) {
                found.push(ri);
            }
        }
        if found.is_empty() {
            return;
        }
        if found.len() == 2 {
            if found[0].is_start != found[1].is_start {
                if found[0].regions.iter().any(|x| found[1].regions.contains(x)) {
                    return;
                }
            } else if found[0].matches(&found[1]) {
                found.remove(1);
            }
        }
        for ri in found {
            let spliced = ri.regions.iter().any(|x| spliced_trans.iter().any(|t| gc.exon_regions[*x].has_trans_id(*t)));
            if let Some(e) = self.introns.iter_mut().find(|x| x.matches(&ri)) {
                e.add_fragment_count(spliced);
            } else {
                let mut ri = ri;
                ri.add_fragment_count(spliced);
                self.introns.push(ri);
            }
        }
    }

    pub fn set_position_depth(&mut self, depth: &BaseDepth) {
        for r in &mut self.introns {
            r.depth = depth.depth_at_base(r.position);
        }
    }

    pub fn write(&self, gc: &GeneCollection, out: &mut Vec<String>) {
        for ri in &self.introns {
            if ri.fragment_count < RI_MIN_FRAG_COUNT && ri.spliced_fragment_count < RI_MIN_SPLICED_FRAG_COUNT {
                continue;
            }
            for g in &gc.genes {
                if !g.transcripts.iter().any(|t| ri.regions.iter().any(|r| gc.exon_regions[*r].has_trans_id(t.trans_id))) {
                    continue;
                }
                let depth = ri.depth.max(ri.fragment_count);
                out.push(format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    g.gene.gene_id,
                    g.gene.gene_name,
                    g.gene.chromosome,
                    g.gene.strand,
                    ri.position,
                    ri.type_name(g.gene.forward_strand()),
                    ri.fragment_count,
                    ri.spliced_fragment_count,
                    depth,
                    ri.transcript_info
                ));
            }
        }
    }
}

fn evaluate_read(gc: &GeneCollection, read: &Read) -> Option<RetainedIntron> {
    let entries = read.mapped_regions.entries();
    if entries.iter().any(|(_, m)| *m != RegionMatchType::ExonIntron) {
        return None;
    }
    let mut spanned_pos = 0;
    let mut spanned_is_start = false;
    let mut candidates: Vec<usize> = Vec::new();
    for (ri, _) in &entries {
        let r = &gc.exon_regions[*ri];
        for se in [SE_START, SE_END] {
            let uses_start = se == SE_START;
            let pos = if uses_start { r.start } else { r.end };
            if !read.mapped_coords.alignments.iter().any(|x| position_within(pos, x.start, x.end)) {
                continue;
            }
            if (uses_start && r.pre_regions.is_empty()) || (!uses_start && r.post_regions.is_empty()) {
                continue;
            }
            spanned_is_start = uses_start;
            if !candidates.is_empty() {
                if uses_start {
                    if pos > spanned_pos {
                        continue;
                    }
                    if spanned_pos > pos {
                        candidates.clear();
                        spanned_pos = pos;
                    }
                } else {
                    if pos < spanned_pos {
                        continue;
                    }
                    if spanned_pos < pos {
                        candidates.clear();
                        spanned_pos = pos;
                    }
                }
            } else {
                spanned_pos = pos;
            }
            candidates.push(*ri);
        }
    }
    if candidates.is_empty() {
        return None;
    }
    Some(RetainedIntron::new(gc, candidates, spanned_is_start))
}
