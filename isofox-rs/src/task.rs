//! ChromosomeTaskExecutor, ChromosomeFusions and the gene-collection part of FusionFinder (the read groups
//! completed across gene collections).

use crate::adjusts::GcRatioCounts;
use crate::allocator::FragmentAllocator;
use crate::common::*;
use crate::config::{IsofoxConfig, IsofoxFunction};
use crate::expression::{ExpectedCountsCache, GeneCollectionSummary, GeneResult, TranscriptExpression, TranscriptResult};
use crate::fusion::{FusionContext, FusionFinder};
use crate::fusion_frag::{FusionRead, FusionReadGroup, JunctionRacFragments};
use crate::gene::{create_gene_read_data, find_next_overlapping_genes, find_unique_bases, GeneCollection};
use crate::novel::AltSjCohortCache;
use crate::refgen::RefSource;
use anyhow::Result;
use hmf_common::ensembl::{EnsemblDataCache, GeneData};
use hmf_common::java_seq_hash::SeqHashMap;
use std::collections::{HashMap, HashSet};

/// the per-chromosome output lines, written in the chromosome order of a single-threaded Java run
#[derive(Default)]
pub struct ChrOutput {
    pub gene_collections: Vec<String>,
    pub alt_sj_unfiltered: Vec<String>,
    pub alt_sj_pass: Vec<String>,
    pub retained_introns: Vec<String>,
}

/// a fusion read group in the finder's arena; merged groups share their reads, as Java's do
#[derive(Clone, Debug)]
pub struct ArenaGroup {
    pub read_id: Box<str>,
    pub reads: Vec<usize>,
    pub complete: bool,
}

/// FusionFinder's read groups that span gene collections
#[derive(Default)]
pub struct GroupArena {
    pub reads: Vec<FusionRead>,
    pub groups: Vec<ArenaGroup>,
    /// mChimericPartialReadGroups
    pub partials: SeqHashMap<String, usize>,
    /// mSpanningReadGroups
    pub spanning: Vec<usize>,
}

impl GroupArena {
    fn add_group(&mut self, g: FusionReadGroup) -> usize {
        let mut ids = Vec::with_capacity(g.reads.len());
        for r in g.reads {
            self.reads.push(r);
            ids.push(self.reads.len() - 1);
        }
        self.groups.push(ArenaGroup { read_id: g.read_id, reads: ids, complete: g.complete });
        self.groups.len() - 1
    }
    fn recompute_complete(&mut self, gid: usize) {
        let reads: Vec<FusionRead> = self.groups[gid].reads.iter().map(|r| self.reads[*r].clone()).collect();
        self.groups[gid].complete = FusionReadGroup::new("", reads).complete;
    }
    /// existing.merge(src): the source's reads are added to (and shared with) the existing group
    fn merge(&mut self, existing: usize, src: usize) {
        let src_reads = self.groups[src].reads.clone();
        for r in src_reads {
            self.groups[existing].reads.push(r);
        }
        self.recompute_complete(existing);
    }
    pub fn snapshot(&self, gid: usize) -> FusionReadGroup {
        let reads: Vec<FusionRead> = self.groups[gid].reads.iter().map(|r| self.reads[*r].clone()).collect();
        let complete = self.groups[gid].complete;
        FusionReadGroup { read_id: self.groups[gid].read_id.clone(), reads, complete }
    }
}

/// FusionReadGroup.mergeChimericReadMaps over the arena
fn merge_chimeric_read_maps(arena: &mut GroupArena, complete: &mut Vec<usize>, source: &[(String, usize)]) {
    for (read_id, gid) in source {
        if arena.groups[*gid].complete {
            complete.push(*gid);
        } else {
            match arena.partials.get(read_id).copied() {
                None => {
                    arena.partials.insert(read_id.clone(), *gid);
                }
                Some(e) => {
                    arena.merge(e, *gid);
                    if arena.groups[e].complete {
                        arena.partials.remove(read_id);
                        complete.push(e);
                    }
                }
            }
        }
    }
}

pub struct ChromosomeResult {
    pub chromosome: &'static str,
    pub summaries: Vec<GeneCollectionSummary>,
    pub combined_counts: FragmentTypeCounts,
    pub non_enriched_gc: GcRatioCounts,
    pub gc_ratio_counts: Option<GcRatioCounts>,
    pub total_reads: i32,
    pub out: ChrOutput,
    pub finder: Option<FusionFinder>,
    pub arena: GroupArena,
    pub rac: HashMap<i32, JunctionRacFragments>,
    pub hard_filtered_read_ids: HashMap<String, HashSet<String>>,
}

pub struct Shared<'a> {
    pub config: &'a IsofoxConfig,
    pub cache: &'static EnsemblDataCache,
    pub expected: Option<&'a ExpectedCountsCache>,
    pub alt_sj_cohort: Option<&'a AltSjCohortCache>,
    pub known_splice_sites: &'a HashMap<String, HashSet<i32>>,
    pub fusion_ctx: Option<&'a FusionContext<'a>>,
    pub known_pairs: &'a crate::chimeric::KnownPairs,
}

/// ChromosomeTaskExecutor.assignTranscriptCounts for one chromosome
pub fn run_chromosome(shared: &Shared, chromosome: &'static str, genes: &[&'static GeneData]) -> Result<ChromosomeResult> {
    let config = shared.config;
    let mut refg = RefSource::open(config.ref_genome_file.as_deref().unwrap_or(""))?;
    let mut alloc = FragmentAllocator::new(config, shared.alt_sj_cohort, shared.known_splice_sites)?;
    alloc.chimeric.register_known_fusion_pairs(shared.known_pairs);
    let expression = shared.expected.filter(|_| config.exp_counts_file.is_some()).map(|c| TranscriptExpression::new(config, c));
    let run_fusions = config.run_function(IsofoxFunction::Fusions);

    let mut result = ChromosomeResult {
        chromosome,
        summaries: Vec::new(),
        combined_counts: FragmentTypeCounts::default(),
        non_enriched_gc: GcRatioCounts::new(),
        gc_ratio_counts: None,
        total_reads: 0,
        out: ChrOutput::default(),
        finder: if run_fusions { Some(FusionFinder::new(chromosome)) } else { None },
        arena: GroupArena::default(),
        rac: HashMap::new(),
        hard_filtered_read_ids: HashMap::new(),
    };

    if genes.len() > 10 {
        hmf_common::info!("chr({}) processing {} genes", chromosome, genes.len());
    }
    let genes_filtered = config.filters.has_specific_regions();
    let mut collection_id = 0;
    let mut idx = 0;
    let mut overlapping: Vec<&'static GeneData> = Vec::new();
    let mut last_end = 1;
    let mut genes_processed = 0;
    let mut next_log = 100;

    while idx < genes.len() {
        idx = find_next_overlapping_genes(genes, idx, &mut overlapping);
        let gene_list = create_gene_read_data(&overlapping, shared.cache);
        if gene_list.is_empty() {
            // Java fails on a collection of genes without transcripts (an index out of bounds)
            anyhow::bail!("genes({}) have no transcripts", overlapping.iter().map(|g| g.gene_name.as_str()).collect::<Vec<_>>().join(";"));
        }
        let mut gc = GeneCollection::new(collection_id, gene_list);
        collection_id += 1;
        if !genes_filtered {
            gc.non_genic_positions[SE_START] = last_end;
            if idx < genes.len() {
                gc.non_genic_positions[SE_END] = genes[idx].gene_start - 1;
            } else {
                let chr_end = config.ref_genome_version.chromosome_length(chromosome);
                gc.non_genic_positions[SE_END] = (gc.non_genic_positions[SE_START] + 1).max(chr_end - 1000);
                gc.end_of_chromosome = true;
            }
        } else {
            gc.non_genic_positions[SE_START] = if last_end == 1 { gc.region_bounds[SE_START] - 10000 } else { last_end };
            gc.non_genic_positions[SE_END] = if idx < genes.len() { genes[idx].gene_start - 1 } else { gc.region_bounds[SE_END] + 10000 };
        }

        analyse_bam_reads(shared, &mut alloc, &mut gc, &mut refg, expression.as_ref(), &mut result)?;

        genes_processed += gc.genes.len();
        result.total_reads = alloc.total_bam_read_count;
        last_end = gc.region_bounds[SE_END] + 1;
        if genes_processed >= next_log {
            next_log += 100;
            hmf_common::info!("chr({}) processed {} of {} genes", chromosome, genes_processed, genes.len());
        }
    }
    if genes.len() > 10 {
        hmf_common::info!("chr({}) processing complete", chromosome);
    }
    result.gc_ratio_counts = alloc.gc_ratio_counts.take();
    result.hard_filtered_read_ids = std::mem::take(&mut alloc.chimeric.hard_filtered_read_ids);
    Ok(result)
}

fn analyse_bam_reads(
    shared: &Shared,
    alloc: &mut FragmentAllocator,
    gc: &mut GeneCollection,
    refg: &mut RefSource,
    expression: Option<&TranscriptExpression>,
    result: &mut ChromosomeResult,
) -> Result<()> {
    let config = shared.config;
    // the reference over the gene collection, read when first needed (its exons, its reads' GC ratios)
    refg.set_window_lazy(gc.chromosome, gc.region_bounds[SE_START] - 1000, gc.region_bounds[SE_END] + 1000);
    if config.ref_genome_file.is_some() {
        for r in &mut gc.exon_regions {
            let b = refg.bases(gc.chromosome, r.start, r.end);
            r.set_ref_bases(b);
        }
        find_unique_bases(&mut gc.exon_regions);
    }
    let region = if config.run_function(IsofoxFunction::Fusions) {
        gc.non_genic_positions
    } else {
        [gc.region_bounds[SE_START] - 100, gc.region_bounds[SE_END] + 100]
    };
    if region[0] >= region[1] {
        hmf_common::warn!("invalid geneCollection({}) region({} -> {})", gc.gene_names(10), region[0], region[1]);
        return Ok(());
    }
    alloc.process_bam(gc, region, refg)?;

    post_bam_read_transcript_counts(shared, alloc, gc, expression, result);
    if config.run_function(IsofoxFunction::AltSpliceJunctions) || config.run_function(IsofoxFunction::RetainedIntrons) {
        let out = &mut result.out;
        alloc.annotate_novel_locations(gc, refg, &mut out.alt_sj_unfiltered, &mut out.alt_sj_pass, &mut out.retained_introns);
    }
    if let (Some(finder), Some(ctx)) = (result.finder.as_mut(), shared.fusion_ctx) {
        on_gene_collection_complete(ctx, refg, finder, &mut result.arena, &mut result.rac, alloc, gc);
    }
    alloc.clear_cache();
    Ok(())
}

fn post_bam_read_transcript_counts(shared: &Shared, alloc: &mut FragmentAllocator, gc: &GeneCollection, expression: Option<&TranscriptExpression>, result: &mut ChromosomeResult) {
    let config = shared.config;
    if config.run_statistics_only() {
        result.combined_counts.combine(&gc.fragment_counts);
        return;
    }
    if !config.run_function(IsofoxFunction::TranscriptCounts) {
        return;
    }
    let mut summary = GeneCollectionSummary::new(gc.chr_id(), gc.gene_ids.clone(), gc.gene_names(10), alloc.trans_combo_data.clone());
    summary.trans_id_names = gc.transcripts.iter().map(|t| (t.trans_id, t.trans_name.as_str())).collect();
    if let Some(e) = expression {
        e.run_transcript_estimation(&mut summary, false);
    }
    for gi in 0..gc.genes.len() {
        for t in &gc.genes[gi].transcripts {
            summary.transcript_results.push(TranscriptResult::new(gc, gi, t, &config.fragment_size_data));
        }
        summary.gene_results.push(GeneResult::new(gc, gi));
        summary.set_fit_allocations();
        summary.assign_multi_mapped_fragments();
        for t in &mut summary.transcript_results {
            t.raw_fit_allocation = t.fit_allocation;
        }
    }
    if let Some(g) = &alloc.gene_gc_ratio_counts {
        result.non_enriched_gc.merge_ratio_counts(&g.counts);
    }
    result.combined_counts.combine(&gc.fragment_counts);
    summary.allocate_residuals_to_genes();
    result.out.gene_collections.push(gene_collection_line(gc));
    if !config.apply_gc_bias_adjust() {
        summary.trans_category_counts.clear();
    }
    result.summaries.push(summary);
}

fn gene_collection_line(gc: &GeneCollection) -> String {
    let c = &gc.fragment_counts;
    [
        gc.chr_id(),
        gc.genes.len().to_string(),
        gc.chromosome.to_string(),
        gc.region_bounds[0].to_string(),
        gc.region_bounds[1].to_string(),
        c.type_count(FragmentType::Total).to_string(),
        c.type_count(FragmentType::Duplicate).to_string(),
        c.type_count(FragmentType::TransSupporting).to_string(),
        c.type_count(FragmentType::Unspliced).to_string(),
        c.type_count(FragmentType::Alt).to_string(),
        c.type_count(FragmentType::Chimeric).to_string(),
        c.type_count(FragmentType::MultiMapped).to_string(),
        c.type_count(FragmentType::ForwardStrand).to_string(),
        c.type_count(FragmentType::ReverseStrand).to_string(),
        gc.gene_names(gc.genes.len()),
    ]
    .join("\t")
}

pub const GENE_COLLECTION_HEADER: &str =
    "GeneSetId\tGeneCount\tChromosome\tRangeStart\tRangeEnd\tTotalFragments\tDuplicates\tSupportingTrans\tUnspliced\tAltSJ\tChimeric\tMultiMapped\tForwardStrand\tReverseStrand\tGenes";

/// ChromosomeFusions.onGeneCollectionComplete
fn on_gene_collection_complete(
    ctx: &FusionContext,
    refg: &mut RefSource,
    finder: &mut FusionFinder,
    arena: &mut GroupArena,
    rac: &mut HashMap<i32, JunctionRacFragments>,
    alloc: &mut FragmentAllocator,
    gc: &GeneCollection,
) {
    let new_groups: Vec<(String, FusionReadGroup)> = alloc.chimeric.fusion_read_group_map.drain_ordered();
    let complete = process_new_chimeric_read_groups(arena, gc, &alloc.base_depth, new_groups);

    let racs = alloc.chimeric.extract_junction_rac_fragments();
    if racs.fragment_count() > 0 {
        rac.insert(gc.id, racs);
    }

    let groups: Vec<FusionReadGroup> = complete.iter().map(|g| arena.snapshot(*g)).collect();
    let chromosome = gc.chromosome;
    let mut lookup = |chr: &str, gcid: i32| -> Option<*mut JunctionRacFragments> {
        if chr != chromosome {
            return None;
        }
        rac.get_mut(&gcid).map(|r| r as *mut JunctionRacFragments)
    };
    finder.process_read_groups(ctx, refg, groups, &mut lookup);
}

/// FusionFinder.processNewChimericReadGroups: the groups now complete, in Java's order
pub fn process_new_chimeric_read_groups(arena: &mut GroupArena, gc: &GeneCollection, depth: &BaseDepth, new_groups: Vec<(String, FusionReadGroup)>) -> Vec<usize> {
    let mut source: Vec<(String, usize)> = Vec::with_capacity(new_groups.len());
    for (id, g) in new_groups {
        let gid = arena.add_group(g);
        source.push((id, gid));
    }
    let mut complete: Vec<usize> = Vec::new();
    merge_chimeric_read_maps(arena, &mut complete, &source);
    let spanning: Vec<usize> = source.iter().map(|(_, g)| *g).filter(|g| arena.groups[*g].reads.iter().any(|r| arena.reads[*r].gene_collections[SE_END] == NO_GENE_ID)).collect();
    let gene_completed = reconcile_spanning_read_groups(arena, gc, &spanning, depth);
    complete.retain(|g| !spanning.contains(g));
    for g in gene_completed {
        if !complete.contains(&g) {
            complete.push(g);
        }
    }
    complete
}

/// FusionFinder.reconcileSpanningReadGroups
fn reconcile_spanning_read_groups(arena: &mut GroupArena, gc: &GeneCollection, spanning: &[usize], depth: &BaseDepth) -> Vec<usize> {
    let mut complete: Vec<usize> = Vec::new();
    let mut i = 0;
    while i < arena.spanning.len() {
        let gid = arena.spanning[i];
        let mut missing = false;
        let reads = arena.groups[gid].reads.clone();
        for ri in reads {
            let read = &mut arena.reads[ri];
            if read.gene_collections[SE_END] != NO_GENE_ID {
                continue;
            }
            let b = read.coords_boundary(SE_END);
            if !position_within(b, gc.non_genic_positions[0], gc.non_genic_positions[1]) {
                missing = true;
                continue;
            }
            if position_within(b, gc.region_bounds[0], gc.region_bounds[1]) {
                read.gene_collections[SE_END] = gc.id;
                read.is_genic[SE_END] = true;
                read.check_missing_gene_data(&gc.transcripts);
            } else {
                read.gene_collections[SE_END] = gc.id;
                read.is_genic[SE_END] = false;
            }
            read.set_read_junction_depth(depth);
        }
        if !missing {
            arena.spanning.remove(i);
            if arena.groups[gid].complete {
                complete.push(gid);
            }
        } else {
            i += 1;
        }
    }
    arena.spanning.extend_from_slice(spanning);
    complete
}

/// FusionFinder.extractIncompleteReadGroups: the partial groups by the chromosome they link to
pub fn extract_incomplete_read_groups(arena: &mut GroupArena, chromosome: &str, hard_filtered: &mut HashMap<String, HashSet<String>>) -> SeqHashMap<String, SeqHashMap<String, FusionReadGroup>> {
    let pair = format!("{chromosome}_{chromosome}");
    if let Some(set) = hard_filtered.remove(&pair) {
        let mut set = set;
        let matched: Vec<String> = set.iter().filter(|x| arena.partials.contains_key(*x)).cloned().collect();
        for m in &matched {
            set.remove(m);
            arena.partials.remove(m);
        }
    }
    let mut out: SeqHashMap<String, SeqHashMap<String, FusionReadGroup>> = SeqHashMap::new();
    for (_, gid) in arena.partials.iter() {
        let g = arena.snapshot(*gid);
        let Some(other) = g.find_other_chromosome(chromosome) else { continue };
        if !hmf_common::segmentation::HumanChromosome::contains(&other) {
            continue;
        }
        if !out.contains_key(&other) {
            out.insert(other.clone(), SeqHashMap::new());
        }
        let id = g.read_id.to_string();
        out.get_mut(&other).unwrap().insert(id, g);
    }
    arena.partials.clear();
    out
}
