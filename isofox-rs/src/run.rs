//! Isofox: the run as a whole. Chromosomes are processed on `-threads` threads; whatever Java writes from
//! its chromosome threads, or orders by them (fusion ids, the inter-chromosomal fusions), is written here
//! in the order of a single-threaded Java run, so the output does not depend on the thread count (IF7).

use crate::adjusts::{calc_fragment_sizes, calc_percentile_data, merge_fragment_lengths, set_config_fragment_length_data, GcRatioCounts, GcTranscriptCalculator};
use crate::chimeric::form_chromosome_pair_string;
use crate::common::*;
use crate::config::*;
use crate::expression::*;
use crate::fusion::{CohortFusions, FusionContext, FusionData};
use crate::fusion_frag::{FusionReadGroup, JunctionRacFragments};
use crate::novel::{alt_sj_unfiltered_header, AltSjCohortCache, ALT_SJ_HEADER, RETAINED_INTRON_HEADER};
use crate::refgen::RefSource;
use crate::task::{extract_incomplete_read_groups, run_chromosome, ChromosomeResult, Shared, GENE_COLLECTION_HEADER};
use anyhow::{Context, Result};
use hmf_common::ensembl::{EnsemblDataCache, GeneData, TranscriptData};
use hmf_common::java_hash::JavaHashMap;
use hmf_common::java_seq_hash::SeqHashMap;
use std::collections::{HashMap, HashSet};
use std::io::Write;

/// An output file. Lines are buffered until `close`; a writer dropped unclosed (a failed run) discards its
/// buffer, as Java's exit leaves its BufferedWriters unflushed, so a run failing before any output leaves
/// the files empty in both.
pub struct Writer {
    path: String,
    out: Option<std::io::BufWriter<std::fs::File>>,
}

impl Writer {
    pub fn create(path: &str) -> Result<Writer> {
        let f = std::fs::File::create(path).with_context(|| format!("failed to create {path}"))?;
        Ok(Writer { path: path.to_string(), out: Some(std::io::BufWriter::with_capacity(1 << 20, f)) })
    }
    pub fn line(&mut self, s: &str) -> Result<()> {
        let out = self.out.as_mut().expect("writer open");
        out.write_all(s.as_bytes())?;
        out.write_all(b"\n")?;
        Ok(())
    }
    pub fn close(mut self) -> Result<()> {
        let mut out = self.out.take().expect("writer open");
        out.flush().with_context(|| format!("failed to write {}", self.path))?;
        Ok(())
    }
}

impl Drop for Writer {
    fn drop(&mut self) {
        if let Some(out) = self.out.take() {
            let _ = out.into_parts();
        }
    }
}

fn chromosome_order(config: &IsofoxConfig, chr_genes: &JavaHashMap<String, Vec<&'static GeneData>>) -> Vec<String> {
    let mut order: Vec<String> = Vec::new();
    for c in PRIORITISED_CHROMOSOMES {
        order.push(config.ref_genome_version.versioned_chromosome(c));
    }
    for c in hmf_common::regions::HUMAN_CHROMOSOMES {
        let v = config.ref_genome_version.versioned_chromosome(c);
        if !order.contains(&v) {
            order.push(v);
        }
    }
    order.into_iter().filter(|c| chr_genes.contains_key(c)).collect()
}

/// getChromosomeGeneLists
fn chromosome_gene_lists(config: &IsofoxConfig, cache: &'static EnsemblDataCache) -> JavaHashMap<String, Vec<&'static GeneData>> {
    let all = cache.chr_gene_data_map();
    let f = &config.filters;
    let mut out: JavaHashMap<String, Vec<&'static GeneData>> = JavaHashMap::new();
    if !f.has_specific_regions() {
        for (k, v) in all.iter() {
            out.insert(k.clone(), v.iter().collect());
        }
        return out;
    }
    if f.specific_regions.regions.is_empty() {
        for (k, v) in all.iter() {
            if f.specific_regions.include_chromosome(k) && !v.is_empty() {
                out.insert(k.clone(), v.iter().collect());
            }
        }
    } else {
        for r in &f.specific_regions.regions {
            // Java dereferences the chromosome's genes unchecked
            let Some(list) = all.get(&r.chromosome) else { continue };
            let in_region: Vec<&'static GeneData> = list.iter().filter(|g| positions_overlap(r.start, r.end, g.gene_start, g.gene_end)).collect();
            if !in_region.is_empty() {
                if !out.contains_key(&r.chromosome) {
                    out.insert(r.chromosome.clone(), Vec::new());
                }
                let l = out.get_mut(&r.chromosome).unwrap();
                for g in in_region {
                    if !l.iter().any(|x| std::ptr::eq(*x, g)) {
                        l.push(g);
                    }
                }
            }
        }
    }
    out
}

/// FusionTaskManager and HardFilteredCache: the inter-chromosomal stage
struct FusionTaskManager {
    incomplete: SeqHashMap<String, SeqHashMap<String, FusionReadGroup>>,
    pair_filtered: HashMap<String, HashSet<String>>,
    hard_filtered_count: i64,
}

impl FusionTaskManager {
    /// addIncompleteReadGroup
    fn add_incomplete_read_group(
        &mut self,
        chromosome: &str,
        mut chr_incomplete: SeqHashMap<String, SeqHashMap<String, FusionReadGroup>>,
        chr_hard_filtered: &mut HashMap<String, HashSet<String>>,
    ) -> Vec<FusionReadGroup> {
        let mut complete: Vec<FusionReadGroup> = Vec::new();
        let entries = chr_incomplete.keys_vec();
        for other in entries {
            if !self.incomplete.contains_key(&other) {
                // putAll: the groups stay in chr_incomplete too, where later hard-filter passes look
                let new_groups: Vec<(String, FusionReadGroup)> =
                    chr_incomplete.get(&other).map(|m| m.iter().map(|(k, v)| (k.clone(), v.clone())).collect()).unwrap_or_default();
                if !self.incomplete.contains_key(&chromosome.to_string()) {
                    self.incomplete.insert(chromosome.to_string(), SeqHashMap::new());
                }
                let target = self.incomplete.get_mut(&chromosome.to_string()).unwrap();
                for (k, v) in new_groups {
                    target.insert(k, v);
                }
            } else {
                // removePartialGroupsWithHardFilteredMatch(existingGroups, chrHardFilteredReadIds)
                {
                    let existing = self.incomplete.get_mut(&other).unwrap();
                    for set in chr_hard_filtered.values_mut() {
                        let matched: Vec<String> = set.iter().filter(|x| existing.contains_key(*x)).cloned().collect();
                        for m in matched {
                            set.remove(&m);
                            existing.remove(&m);
                        }
                    }
                }
                self.remove_hard_filtered_reads(chromosome, &mut chr_incomplete, chr_hard_filtered);
                let new_groups: Vec<(String, FusionReadGroup)> =
                    chr_incomplete.get(&other).map(|m| m.iter().map(|(k, v)| (k.clone(), v.clone())).collect()).unwrap_or_default();
                let existing = self.incomplete.get_mut(&other).unwrap();
                for (read_id, src) in new_groups {
                    if src.complete {
                        complete.push(src);
                    } else {
                        match existing.get_mut(&read_id) {
                            None => {
                                existing.insert(read_id, src);
                            }
                            Some(e) => {
                                e.merge(src);
                                if e.complete {
                                    let g = existing.remove(&read_id).unwrap();
                                    complete.push(g);
                                }
                            }
                        }
                    }
                }
                self.pair_filtered.remove(&form_chromosome_pair_string(chromosome, &other));
            }
        }
        // addHardFilteredReads
        for (pair, set) in chr_hard_filtered.iter() {
            if set.is_empty() {
                continue;
            }
            let s = self.pair_filtered.entry(pair.clone()).or_default();
            for id in set {
                s.insert(id.clone());
            }
            self.hard_filtered_count += set.len() as i64;
        }
        hmf_common::info!("chr({}) complete({}) inter-chromosomal fusion groups", chromosome, complete.len());
        complete
    }

    /// HardFilteredCache.removeHardFilteredReads
    fn remove_hard_filtered_reads(&mut self, chromosome: &str, chr_incomplete: &mut SeqHashMap<String, SeqHashMap<String, FusionReadGroup>>, chr_hard_filtered: &mut HashMap<String, HashSet<String>>) {
        for other in chr_incomplete.keys_vec() {
            let pair = form_chromosome_pair_string(chromosome, &other);
            let groups = chr_incomplete.get_mut(&other).unwrap();
            let filtered: Vec<String> = groups.keys_vec().into_iter().filter(|x| self.pair_filtered.get(&pair).map(|s| s.contains(x)).unwrap_or(false)).collect();
            for f in &filtered {
                groups.remove(f);
            }
            if let (Some(new), Some(existing)) = (chr_hard_filtered.get_mut(&pair), self.pair_filtered.get_mut(&pair)) {
                let matched: Vec<String> = new.iter().filter(|x| existing.contains(*x)).cloned().collect();
                for m in matched {
                    new.remove(&m);
                    existing.remove(&m);
                }
            }
            if let Some(new) = chr_hard_filtered.get_mut(&pair) {
                for f in &filtered {
                    new.remove(f);
                }
            }
        }
    }
}

pub fn run(cb: &hmf_common::config_builder::ConfigBuilder) -> Result<i32> {
    let started = std::time::Instant::now();
    let mut config = IsofoxConfig::from(cb)?;
    if !config.is_valid() {
        hmf_common::error!("missing config options, exiting");
        return Ok(1);
    }

    // output files Java creates at the start
    let gene_collection_writer = if config.run_function(IsofoxFunction::TranscriptCounts) {
        let mut w = Writer::create(&config.form_output_file("gene_collection.tsv"))?;
        w.line(GENE_COLLECTION_HEADER)?;
        Some(w)
    } else {
        None
    };
    let not_fusions_only = !config.run_fusions_only();
    let alt_sj_writers = if not_fusions_only && config.run_function(IsofoxFunction::AltSpliceJunctions) {
        let mut u = Writer::create(&config.form_output_file("alt_splice_junc_unfiltered.tsv"))?;
        u.line(&alt_sj_unfiltered_header())?;
        let mut p = Writer::create(&config.form_output_file("alt_splice_junc.tsv"))?;
        p.line(ALT_SJ_HEADER)?;
        Some((u, p))
    } else {
        None
    };
    let retained_writer = if not_fusions_only && config.run_function(IsofoxFunction::RetainedIntrons) {
        let mut w = Writer::create(&config.form_output_file("retained_intron.tsv"))?;
        w.line(RETAINED_INTRON_HEADER)?;
        Some(w)
    } else {
        None
    };

    // the Ensembl data, kept for the life of the program
    let ensembl_dir = config.ensembl_data_dir.clone().unwrap_or_default();
    let mut ens = EnsemblDataCache::new(&ensembl_dir, config.ref_genome_version);
    ens.set_required_data(true, false, false, config.canonical_only);
    ens.load(false)?;
    let cache: &'static EnsemblDataCache = Box::leak(Box::new(ens));

    let expected = if config.exp_counts_file.is_some() || config.apply_gc_bias_adjust() { Some(ExpectedCountsCache::new(&mut config)?) } else { None };
    let mut gc_calc = if config.apply_gc_bias_adjust() { Some(GcTranscriptCalculator::new(&config)?) } else { None };

    let run_fusions = config.run_function(IsofoxFunction::Fusions);
    let mut trans_by_id: HashMap<i32, &'static TranscriptData> = HashMap::new();
    let mut known_splice_sites: HashMap<String, HashSet<i32>> = HashMap::new();
    let cohort_fusions = if run_fusions { CohortFusions::load(config.fusions.cohort_file.as_deref())? } else { CohortFusions::default() };
    let fusion_writers = if run_fusions {
        for (_, list) in cache.transcript_data_map().iter() {
            for t in list {
                trans_by_id.insert(t.trans_id, t);
            }
        }
        for (chr, genes) in cache.chr_gene_data_map().iter() {
            let set = known_splice_sites.entry(chr.clone()).or_default();
            for g in genes {
                if let Some(ts) = cache.transcripts(&g.gene_id) {
                    for t in ts {
                        for e in &t.exons {
                            set.insert(e.start);
                            set.insert(e.end);
                        }
                    }
                }
            }
        }
        let mut f = Writer::create(&config.form_output_file("fusions.tsv"))?;
        f.line(&FusionData::header())?;
        let mut p = Writer::create(&config.form_output_file("pass_fusions.tsv"))?;
        p.line(FusionData::pass_header())?;
        Some((f, p))
    } else {
        None
    };

    let alt_sj_cohort = match &config.alt_sj_cohort_file {
        Some(f) if config.run_function(IsofoxFunction::AltSpliceJunctions) => Some(AltSjCohortCache::load(f)?),
        _ => None,
    };
    let gene_distribution = match &config.gene_distribution_file {
        Some(f) if config.run_function(IsofoxFunction::TranscriptCounts) => Some(CohortGenePercentiles::load(f, config.cancer_type.as_deref())?),
        _ => None,
    };

    let chr_genes = chromosome_gene_lists(&config, cache);
    if chr_genes.is_empty() {
        hmf_common::error!("no chromosome tasks created");
        return Ok(1);
    }

    if config.run_function(IsofoxFunction::NeoEpitopes) {
        crate::neo::calc_fragment_support(&config, cache)?;
        return Ok(0);
    }

    hmf_common::info!("sample({}) running RNA analysis", config.sample_id);
    if let Some(e) = &expected {
        if !e.valid {
            hmf_common::warn!("invalid expected counts cache");
            return Ok(1);
        }
    }

    // fragment lengths
    let mut max_read_length = 0;
    let mut fragment_lengths: Vec<FragmentSize> = Vec::new();
    if config.require_fragment_length_calcs() {
        let required = config.fragment_length_sampling_count / chr_genes.len() as i32;
        let chrs: Vec<(String, Vec<&'static GeneData>)> = chr_genes.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        let results = parallel_map(config.threads, &chrs, |(chr, genes)| calc_fragment_sizes(&config, cache, chr, genes, required))?;
        let mut gene_lines: Vec<String> = Vec::new();
        for r in &results {
            max_read_length = max_read_length.max(r.max_read_length);
            merge_fragment_lengths(&mut fragment_lengths, &r.fragment_lengths);
            gene_lines.extend(r.gene_lines.iter().cloned());
        }
        if config.apply_fragment_length_adjust {
            hmf_common::info!("max observed read length({}) set", max_read_length);
            set_config_fragment_length_data(&mut config, &fragment_lengths);
        }
        if config.write_type(WriteType::FragLength) && !fragment_lengths.is_empty() {
            let mut w = Writer::create(&config.form_output_file("frag_length.csv"))?;
            w.line("FragmentLength,Count")?;
            for f in &fragment_lengths {
                w.line(&format!("{},{}", f.length, f.frequency))?;
            }
            w.close()?;
        }
        if config.write_type(WriteType::FragLengthByGene) {
            let mut w = Writer::create(&config.form_output_file("frag_length_by_gene.csv"))?;
            w.line("GeneNames,GeneCount,Chromosome,GenesStart,GenesEnd,FragmentLength,Count")?;
            for l in &gene_lines {
                w.line(l)?;
            }
            w.close()?;
            close_all(gene_collection_writer, alt_sj_writers, retained_writer, fusion_writers)?;
            return Ok(0);
        }
    }

    // chromosome tasks
    let order = chromosome_order(&config, &chr_genes);
    let fusion_ctx = FusionContext {
        cache,
        trans_by_id: &trans_by_id,
        known: &config.fusions.known_fusions,
        filters: &config.filters,
        max_fragment_length: config.max_fragment_length,
        cohort: &cohort_fusions,
        min_hard_filter_frags: config.fusions.min_hard_filter_frags,
    };
    let known_pairs = crate::chimeric::KnownPairs::compute(&config, cache);
    let shared = Shared {
        known_pairs: &known_pairs,
        config: &config,
        cache,
        expected: expected.as_ref(),
        alt_sj_cohort: alt_sj_cohort.as_ref(),
        known_splice_sites: &known_splice_sites,
        fusion_ctx: if run_fusions { Some(&fusion_ctx) } else { None },
    };
    let tasks: Vec<(&'static str, Vec<&'static GeneData>)> = order
        .iter()
        .map(|c| {
            let genes = chr_genes.get(c).unwrap().clone();
            let chr: &'static str = &genes[0].chromosome;
            (chr, genes)
        })
        .collect();
    let mut results: Vec<ChromosomeResult> = parallel_map(config.threads, &tasks, |(chr, genes)| run_chromosome(&shared, chr, genes))?;
    let total_reads: i64 = results.iter().map(|r| r.total_reads as i64).sum();
    hmf_common::info!("read {} total BAM records", total_reads);

    // the inter-chromosomal fusions, chromosome by chromosome as a single Java thread completes them
    let mut fusion_lines: Vec<(usize, Vec<String>, Vec<String>)> = Vec::new();
    if run_fusions {
        let mut rac_all: HashMap<&'static str, HashMap<i32, JunctionRacFragments>> = HashMap::new();
        for r in &mut results {
            rac_all.insert(r.chromosome, std::mem::take(&mut r.rac));
        }
        let mut manager = FusionTaskManager { incomplete: SeqHashMap::new(), pair_filtered: HashMap::new(), hard_filtered_count: 0 };
        let mut refg = RefSource::open(config.ref_genome_file.as_deref().unwrap_or(""))?;
        for r in results.iter_mut() {
            let mut hard = std::mem::take(&mut r.hard_filtered_read_ids);
            let incomplete = extract_incomplete_read_groups(&mut r.arena, r.chromosome, &mut hard);
            let groups = manager.add_incomplete_read_group(r.chromosome, incomplete, &mut hard);
            if let Some(finder) = r.finder.as_mut() {
                let mut lookup = |chr: &str, gcid: i32| -> Option<*mut JunctionRacFragments> {
                    rac_all.get_mut(chr).and_then(|m| m.get_mut(&gcid)).map(|x| x as *mut JunctionRacFragments)
                };
                finder.process_read_groups(&fusion_ctx, &mut refg, groups, &mut lookup);
                if finder.hard_filtered_count > 0 {
                    hmf_common::info!(
                        "chr({}) fusion processing complete, filtered(hard={} excluded={})",
                        r.chromosome,
                        finder.hard_filtered_count,
                        finder.excluded_filtered_count
                    );
                }
            }
            r.arena = Default::default();
        }
        // final fusion ids: each chromosome's ids follow those of the chromosomes before it
        let mut offset = 0;
        for (ci, r) in results.iter().enumerate() {
            let Some(finder) = &r.finder else { continue };
            let base = offset;
            let id_of = move |id: i32| base + id;
            let mut lines = Vec::new();
            let mut pass = Vec::new();
            for b in &finder.outputs {
                for f in &b.fusions {
                    lines.push(f.to_tsv(&id_of));
                }
                for p in &b.passing {
                    pass.push(b.fusions[*p].pass_line());
                }
            }
            fusion_lines.push((ci, lines, pass));
            offset += finder.ids_used();
        }
        hmf_common::info!("all fusion tasks complete: hardFiltered({})", manager.hard_filtered_count);
    }

    // write the per-chromosome lines in the single-thread order
    let mut gene_collection_writer = gene_collection_writer;
    let mut alt_sj_writers = alt_sj_writers;
    let mut retained_writer = retained_writer;
    for r in &results {
        if let Some(w) = gene_collection_writer.as_mut() {
            for l in &r.out.gene_collections {
                w.line(l)?;
            }
        }
        if let Some((u, p)) = alt_sj_writers.as_mut() {
            for l in &r.out.alt_sj_unfiltered {
                u.line(l)?;
            }
            for l in &r.out.alt_sj_pass {
                p.line(l)?;
            }
        }
        if let Some(w) = retained_writer.as_mut() {
            for l in &r.out.retained_introns {
                w.line(l)?;
            }
        }
    }
    let mut fusion_writers = fusion_writers;
    if let Some((f, p)) = fusion_writers.as_mut() {
        for (_, lines, pass) in &fusion_lines {
            for l in lines {
                f.line(l)?;
            }
            for l in pass {
                p.line(l)?;
            }
        }
    }

    if !config.run_fusions_only() {
        process_bam_fragments(&config, &mut results, gc_calc.as_mut(), gene_distribution.as_ref(), cache, &fragment_lengths, max_read_length, expected.as_ref())?;
    }
    close_all(gene_collection_writer, alt_sj_writers, retained_writer, fusion_writers)?;
    hmf_common::info!("Isofox complete, mins({:.3})", started.elapsed().as_secs_f64() / 60.0);
    Ok(0)
}

fn close_all(a: Option<Writer>, b: Option<(Writer, Writer)>, c: Option<Writer>, d: Option<(Writer, Writer)>) -> Result<()> {
    if let Some(w) = a {
        w.close()?;
    }
    if let Some((x, y)) = b {
        x.close()?;
        y.close()?;
    }
    if let Some(w) = c {
        w.close()?;
    }
    if let Some((x, y)) = d {
        x.close()?;
        y.close()?;
    }
    Ok(())
}

/// runs `f` over the items on `threads` threads, in the order of a thread pool taking them in turn; the
/// results are in item order
pub fn parallel_map<T: Sync, R: Send, F: Fn(&T) -> Result<R> + Sync>(threads: usize, items: &[T], f: F) -> Result<Vec<R>> {
    let n = items.len();
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<Result<R>>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads.max(1).min(n.max(1)) {
            s.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if i >= n {
                    break;
                }
                let r = f(&items[i]);
                *slots[i].lock().unwrap() = Some(r);
            });
        }
    });
    let mut out = Vec::with_capacity(n);
    for s in slots {
        out.push(s.into_inner().unwrap().unwrap()?);
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn process_bam_fragments(
    config: &IsofoxConfig,
    results: &mut [ChromosomeResult],
    gc_calc: Option<&mut GcTranscriptCalculator>,
    gene_distribution: Option<&CohortGenePercentiles>,
    cache: &'static EnsemblDataCache,
    fragment_lengths: &[FragmentSize],
    max_read_length: i32,
    expected: Option<&ExpectedCountsCache>,
) -> Result<()> {
    let mut total = FragmentTypeCounts::default();
    let mut gc_counts = GcRatioCounts::new();
    for r in results.iter() {
        total.combine(&r.combined_counts);
        gc_counts.merge_ratio_counts(&r.non_enriched_gc.counts);
    }
    if let Some(calc) = gc_calc {
        hmf_common::info!("applying GC adjustments and transcript re-fit");
        for r in results.iter() {
            let allocs: Vec<(&str, &JavaHashMap<String, f64>)> = r.summaries.iter().map(|s| (s.gene_names.as_str(), &s.fit_allocations)).collect();
            calc.generate_gc_counts_from_fit(&allocs);
        }
        calc.calc_gc_ratio_adjustments(&gc_counts);
        let adj = calc.gc_ratio_adjustments.clone();
        if let Some(exp) = expected {
            let te = TranscriptExpression::new(config, exp);
            let te = &te;
            let adj = &adj;
            let mut chunks: Vec<&mut ChromosomeResult> = results.iter_mut().collect();
            parallel_for_each(config.threads, &mut chunks, |r| {
                for s in &mut r.summaries {
                    s.apply_gc_adjustments(adj);
                    te.run_transcript_estimation(s, true);
                    s.set_fit_allocations();
                    s.allocate_residuals_to_genes();
                }
            });
        }
    }

    let all: Vec<&GeneCollectionSummary> = results.iter().flat_map(|r| r.summaries.iter()).collect();
    let factors = calc_tpm_factors(&all);
    let panel = PanelTpmNormaliser::new(config.panel_tpm_norm_file.as_deref())?;
    let mut splice_gene_count = 0;
    let mut gene_writer: Option<Writer> = None;
    let mut trans_writer: Option<Writer> = None;
    for r in results.iter_mut() {
        for s in &mut r.summaries {
            set_transcripts_per_million(s, factors);
            panel.apply(s);
            set_cohort_distribution_values(s, gene_distribution, config.cancer_type.as_deref(), &config.driver_genes);
            splice_gene_count += s.splice_genes_count();
        }
        // writeExpressionResults
        for s in &r.summaries {
            for g in &s.gene_results {
                if !panel.gene_ids.is_empty() && !panel.gene_ids.contains(&g.gene.gene_id) {
                    continue;
                }
                if gene_writer.is_none() {
                    let mut w = Writer::create(&config.form_output_file("gene_data.tsv"))?;
                    w.line(&GeneResult::header())?;
                    gene_writer = Some(w);
                }
                gene_writer.as_mut().unwrap().line(&g.to_line())?;
            }
            for t in &s.transcript_results {
                if !panel.gene_ids.is_empty() && !panel.gene_ids.contains(&t.trans.gene_id) {
                    continue;
                }
                let gene = s.gene_results.iter().find(|g| g.gene.gene_id == t.trans.gene_id).map(|g| g.gene);
                let Some(gene) = gene else { continue };
                if trans_writer.is_none() {
                    let mut w = Writer::create(&config.form_output_file("transcript_data.tsv"))?;
                    w.line(&TranscriptResult::header())?;
                    trans_writer = Some(w);
                }
                trans_writer.as_mut().unwrap().line(&t.to_line(gene))?;
            }
        }
    }
    if let Some(w) = gene_writer {
        w.close()?;
    }
    if let Some(w) = trans_writer {
        w.close()?;
    }

    if config.run_function(IsofoxFunction::TranscriptCounts) || config.run_function(IsofoxFunction::Statistics) {
        let median_gc = gc_counts.percentile_ratio(0.5);
        let mut low_cov = 2_500_000;
        let mut spliced_thr = 17000;
        if panel.enabled() {
            let total_genes: usize = cache.chr_gene_data_map().values().map(|v| v.len()).sum();
            let cov = panel.gene_count() as f64 / total_genes as f64;
            low_cov = (low_cov as f64 * cov * PANEL_LOW_COVERAGE_FACTOR) as i32;
            spliced_thr = (cov * 17000.0) as i32;
        }
        let read_length = if max_read_length > 0 { max_read_length } else { config.read_length };
        write_summary(config, &total, splice_gene_count, median_gc, fragment_lengths, read_length, low_cov, spliced_thr)?;
    }
    Ok(())
}

fn parallel_for_each<T: Send, F: Fn(&mut T) + Sync>(threads: usize, items: &mut [T], f: F) {
    let chunk = items.len().div_ceil(threads.max(1)).max(1);
    std::thread::scope(|s| {
        for part in items.chunks_mut(chunk) {
            let f = &f;
            s.spawn(move || {
                for x in part {
                    f(x);
                }
            });
        }
    });
}

#[allow(clippy::too_many_arguments)]
fn write_summary(config: &IsofoxConfig, c: &FragmentTypeCounts, splice_genes: i32, median_gc: f64, fl: &[FragmentSize], read_length: i32, low_cov: i32, spliced_thr: i32) -> Result<()> {
    use hmf_common::fmt::java_fixed;
    let total = c.type_count(FragmentType::Total);
    let dups = c.type_count(FragmentType::Duplicate);
    let denom = total as f64;
    let fwd = c.type_count(FragmentType::ForwardStrand);
    let strand_total = (fwd + c.type_count(FragmentType::ReverseStrand)) as f64;
    let fwd_perc = if strand_total > 0.0 { fwd as f64 / strand_total } else { 0.0 };
    let lengths = calc_percentile_data(fl, &[0.05, 0.5, 0.95]);
    let mut qc: Vec<&str> = Vec::new();
    if total - dups < low_cov as i64 {
        qc.push("FAIL_LOW_COVERAGE");
    } else {
        let dup_perc = if total > 0 { dups as f64 / total as f64 } else { 0.0 };
        if dup_perc > 0.9 {
            qc.push("WARN_DUPLICATE_RATE");
        }
        if splice_genes >= 0 && splice_genes < spliced_thr {
            qc.push("WARN_SPLICED_GENE_COVERAGE");
        }
        if qc.is_empty() {
            qc.push("PASS");
        }
    }
    let pick = |i: usize| if lengths.is_empty() { 0.0 } else { lengths[i] };
    let line = [
        config.sample_id.clone(),
        qc.join(";"),
        total.to_string(),
        dups.to_string(),
        java_fixed(c.type_count(FragmentType::TransSupporting) as f64 / denom, 3),
        java_fixed(c.type_count(FragmentType::Unspliced) as f64 / denom, 3),
        java_fixed(c.type_count(FragmentType::Alt) as f64 / denom, 3),
        java_fixed(c.type_count(FragmentType::Chimeric) as f64 / denom, 3),
        splice_genes.to_string(),
        read_length.to_string(),
        java_fixed(pick(0), 0),
        java_fixed(pick(1), 0),
        java_fixed(pick(2), 0),
        java_fixed(median_gc, 3),
        java_fixed(fwd_perc, 3),
    ]
    .join("\t");
    let mut w = Writer::create(&config.form_output_file("summary.tsv"))?;
    w.line("SampleId\tQcStatus\tTotalFragments\tDuplicateFragments\tSplicedFragmentPerc\tUnsplicedFragmentPerc\tAltFragmentPerc\tChimericFragmentPerc\tSplicedGeneCount\tReadLength\tFragLength5th\tFragLength50th\tFragLength95th\tMedianGCRatio\tForwardStrandPercent")?;
    w.line(&line)?;
    w.close()
}
