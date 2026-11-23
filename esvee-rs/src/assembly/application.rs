//! AssemblyApplication: the assembly step from junctions and prep reads to assemblies, alignments and breakends.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::{bail, Context, Result};

use crate::common::compat;
use crate::assembly::alignment::alignment_checker::{self, AlignmentChecker};
use crate::assembly::alignment::bwa::{BwaMemAlignParams, BwaMemAligner, IBwaMemAligner};
use crate::assembly::alignment::assembly_aligner::run_alignment;
use crate::assembly::alignment::breakend::BreakendP;
use crate::assembly::config::{form_prep_input_filename, AssemblyConfig, WriteType};
use crate::assembly::output::assembly_read_writer::AssemblyReadWriter;
use crate::assembly::output::breakend_writer::BreakendWriter;
use crate::assembly::output::vcf_writer::VcfWriter;
use crate::assembly::constants::*;
use crate::assembly::globals::*;
use crate::assembly::junction_assembler::MIN_HOTSPOT_JUNCTION_SUPPORT;
use crate::assembly::junction_group_assembler::JunctionGroupAssembler;
use crate::assembly::mem_stats;
use crate::assembly::alignment::assembly_alignment::{AsmAlignP, AssemblyAlignment};
use crate::assembly::assembly_utils::set_assembly_outcome;
use crate::assembly::indel_builder::is_weak_indel_based_unlinked_assembly;
use crate::assembly::output::assembly_writer::AssemblyWriter;
use crate::assembly::phase::phase_group_builder::{PhaseGroupBuildWriter, PhaseGroupBuilder};
use crate::assembly::phase::phase_set_builder::PhaseSetBuilder;
use crate::assembly::phase::remote_read_extractor::RemoteReadExtractor;
use crate::assembly::perf::PerfCounter;
use crate::assembly::read::bam_reader::BamReader;
use crate::assembly::saga::SagaMatcherFactory;
use crate::assembly::types::phase_group::PhaseGroupP;
use crate::assembly::types::{AssemblyOutcome, junction_assembly::JunctionAssembly};
use crate::common::java_id_hash::IdHashSet;
use crate::assembly::read::read_stats::ReadStats;
use crate::assembly::types::junction::{JuncP, Junction};
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::junction_group::{JuncGroupP, JunctionGroup};
use crate::common::fmt::java_fixed;
use crate::common::java_sort::java_sort;

// PrepConstants
const MIN_JUNCTION_SUPPORT: i32 = 2;
const DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT: i32 = 5;
const FRAG_LENGTH_DIST_PERCENTILE: f64 = 0.9975;
const FRAG_LENGTH_1_STD_DEV_PERCENTILE: f64 = 0.16;

/// ResultsWriter
pub struct ResultsWriter {
    pub assembly_writer: AssemblyWriter,
    pub read_writer: AssemblyReadWriter,
    pub breakend_writer: BreakendWriter,
    pub phase_group_build_writer: PhaseGroupBuildWriter,
}

pub struct AssemblyApplication {
    config: AssemblyConfig,
    chr_junctions_map: HashMap<String, Vec<JuncP>>,
    junction_group_map: BTreeMap<String, Vec<JuncGroupP>>,
    task_count: usize,
    results: ResultsWriter,
    saga_matcher_factory: Option<SagaMatcherFactory>,
    perf_counters: Vec<PerfCounter>,
}

pub fn run(argv: &[String]) -> Result<()> {
    if crate::version_log_due() {
        crate::info!("Esvee version {}", crate::VERSION);
    }
    let config = AssemblyConfig::from_args(argv)?;
    let results = ResultsWriter { assembly_writer: AssemblyWriter::new(&config), read_writer: AssemblyReadWriter::new(&config), breakend_writer: BreakendWriter::new(&config), phase_group_build_writer: PhaseGroupBuildWriter::new(&config) };
    let saga_matcher_factory = match &config.saga_fasta_file {
        Some(f) => match SagaMatcherFactory::new(f) {
            Ok(factory) => Some(factory),
            Err(e) => {
                // Java ends here with an uncaught exception, leaving the output files it has opened empty
                crate::error!("thread exception: {}", e);
                std::mem::forget(results);
                std::process::exit(1);
            }
        },
        None => None,
    };
    let mut app = AssemblyApplication { config, chr_junctions_map: HashMap::new(), junction_group_map: BTreeMap::new(), task_count: 1, results, saga_matcher_factory, perf_counters: Vec::new() };
    app.run()
}

/// FragmentSizeDistribution.loadFragmentLengthBounds: (lower, upper) bounds
fn load_fragment_length_bounds(filename: &str) -> Result<(i32, i32)> {
    let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read fragment length file: {filename}"))?;
    let mut length_frequencies: Vec<(i32, i32)> = Vec::new();
    for line in text.lines().skip(1) {
        let mut values = line.split('\t');
        let length: i32 = values.next().unwrap_or("").parse()?;
        let frequency: i32 = values.next().unwrap_or("").parse()?;
        length_frequencies.push((length, frequency));
    }
    if length_frequencies.is_empty() {
        return Ok((0, 0));
    }
    let total_fragments: i32 = length_frequencies.iter().map(|x| x.1).sum();
    let mut cumulative_total = 0;
    let required_min_total = (total_fragments as f64 * (1.0 - FRAG_LENGTH_DIST_PERCENTILE)).floor() as i32;
    let required_max_total = (total_fragments as f64 * FRAG_LENGTH_DIST_PERCENTILE).floor() as i32;
    let mut lower_bound = 0;
    let mut upper_bound = 0;
    for (length, frequency) in length_frequencies {
        if lower_bound == 0 && frequency + cumulative_total >= required_min_total {
            lower_bound = length;
        }
        if upper_bound == 0 && frequency + cumulative_total >= required_max_total {
            upper_bound = length;
        } else {
            cumulative_total += frequency;
        }
    }
    let _ = FRAG_LENGTH_1_STD_DEV_PERCENTILE;
    Ok((lower_bound, upper_bound))
}

/// DiscordantStats.loadDiscordantStats(...).discordantRate()
fn load_discordant_rate(filename: &str) -> Result<f64> {
    let text = std::fs::read_to_string(filename).with_context(|| format!("failed to load Esvee discordant stat file({filename})"))?;
    let line = text.lines().nth(1).context("discordant stats file has no data line")?;
    let values: Vec<i64> = line.split('\t').map(|v| v.parse::<i64>()).collect::<std::result::Result<_, _>>()?;
    let total_reads = values[0];
    let total_discordant: i64 = values[2..].iter().sum();
    Ok(if total_reads > 0 { total_discordant as f64 / total_reads as f64 } else { 0.0 })
}

impl AssemblyApplication {
    fn run(&mut self) -> Result<()> {
        let start_time = Instant::now();
        crate::info!(
            "writing to output directory({}){}",
            self.config.output_dir.as_deref().unwrap_or("null"),
            self.config.output_id.as_ref().map(|id| format!(" outputId({id})")).unwrap_or_default()
        );
        self.load_fragment_length_bounds()?;
        let disc_stats_filename = form_prep_input_filename(&self.config.prep_dir, self.config.sample_id(), "disc_stats.tsv", self.config.output_id.as_deref());
        set_sample_discordant_rate(load_discordant_rate(&disc_stats_filename)?);
        self.load_junction_files()?;

        for (chr, junctions) in &self.chr_junctions_map {
            self.junction_group_map.insert(chr.clone(), JunctionGroup::build_junction_groups(junctions, BAM_READ_JUNCTION_BUFFER));
        }
        let junction_count: usize = self.chr_junctions_map.values().map(|x| x.len()).sum();
        self.task_count = junction_count.min(self.config.threads.max(0) as usize);

        mem_stats::init_from_env();
        mem_stats::log_stage("start");
        self.run_primary_assembly()?;
        mem_stats::log_stage("primary assembly done");
        self.form_phase_groups()?;
        mem_stats::log_stage("phase sets done");
        self.release_retained_reads();
        mem_stats::log_stage("retained reads released");
        mem_stats::log_assembly_breakdown(&self.junction_group_map.values().flatten().copied().collect::<Vec<_>>());
        let (final_assemblies, assembly_alignments) = self.gather_assemblies();
        self.run_alignment(&assembly_alignments)?;
        mem_stats::log_stage("alignment done");

        let mut breakends: Vec<BreakendP> = Vec::new();
        for x in &assembly_alignments {
            breakends.extend(x.breakends().iter().copied());
        }
        java_sort(&mut breakends, |a, b| a.compare_to(b));
        for (i, breakend) in breakends.iter().enumerate() {
            breakend.m().set_id(i as i32);
        }
        self.write_assembly_output(&final_assemblies);
        self.write_variants(&assembly_alignments, &breakends)?;
        if self.config.run_visualiser {
            crate::assembly::vis::run_visualiser(&self.config, &assembly_alignments);
        }
        self.results.assembly_writer.close();
        self.results.read_writer.close();
        self.results.breakend_writer.close();
        self.results.phase_group_build_writer.close();
        mem_stats::log_stage("output written");

        if perf_debug() {
            for perf_counter in &self.perf_counters {
                perf_counter.log_stats();
            }
        }
        let secs = start_time.elapsed().as_secs_f64();
        crate::info!("Esvee assembly complete, mins({})", java_fixed(secs / 60.0, 3));
        Ok(())
    }

    fn load_fragment_length_bounds(&self) -> Result<()> {
        let fragment_length_file = form_prep_input_filename(&self.config.prep_dir, self.config.sample_id(), "fragment_length.tsv", self.config.output_id.as_deref());
        if !Path::new(&fragment_length_file).exists() {
            bail!("missing fragment length file: {fragment_length_file}");
        }
        let (lower, upper) = load_fragment_length_bounds(&fragment_length_file)?;
        if lower > 0 && upper > lower {
            crate::info!("fragment length bounds(min={} max={})", lower, upper);
            set_max_observed_concordant_frag_length(upper);
        }
        Ok(())
    }

    fn load_junction_files(&mut self) -> Result<()> {
        if !self.config.specific_junctions.is_empty() {
            for junction in &self.config.specific_junctions {
                self.chr_junctions_map.entry(junction.chromosome.to_string()).or_default().push(*junction);
            }
            return Ok(());
        }
        let Some(junction_file) = self.config.junction_file.clone().filter(|f| Path::new(f).exists()) else {
            bail!("invalid junction file({})", self.config.junction_file.as_deref().unwrap_or("null"));
        };
        let mut min_junction_frags = MIN_JUNCTION_SUPPORT;
        let mut min_hotspot_frags = MIN_HOTSPOT_JUNCTION_SUPPORT;
        let mut min_discordant_frags = DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
        let rate = sample_discordant_rate();
        if rate >= self.config.discordant_rate_increment {
            let discordant_rate_factor = (rate / self.config.discordant_rate_increment).floor() as i32;
            min_junction_frags += discordant_rate_factor * DISC_RATE_JUNC_INCREMENT;
            min_hotspot_frags += discordant_rate_factor * DISC_RATE_JUNC_INCREMENT;
            min_discordant_frags += discordant_rate_factor * DISC_RATE_DISC_ONLY_INCREMENT;
            crate::info!("raised min fragments(hotspot={} junction={} disc-only={}) for discordantRate({})", min_hotspot_frags, min_junction_frags, min_discordant_frags, java_fixed(rate, 3));
        }
        let loaded = Junction::load_junctions(&junction_file, &self.config.specific_chr_regions, min_junction_frags, min_hotspot_frags, min_discordant_frags)?;
        self.chr_junctions_map.extend(loaded);
        Ok(())
    }

    fn run_primary_assembly(&mut self) -> Result<()> {
        let task_count = self.task_count;
        let mut junction_groups: Vec<JuncGroupP> = Vec::new();
        for groups in self.junction_group_map.values() {
            junction_groups.extend(groups.iter().copied());
        }
        java_sort(&mut junction_groups, |a, b| (-a.range()).cmp(&(-b.range())));

        let config = &self.config;
        // one reference aligner with default parameters shared by the tasks (Java opens one per task)
        let aligner: Option<Arc<dyn IBwaMemAligner>> = if config.assembly_map_qual_threshold > 0 {
            let image = config.ref_genome_image_file.as_deref().context("no ref genome image file")?;
            Some(Arc::new(BwaMemAligner::new(image, &BwaMemAlignParams::DEFAULT)?))
        } else {
            None
        };
        let decoy_aligner: Option<Arc<dyn IBwaMemAligner>> = match &config.decoy_genome {
            Some(d) => Some(Arc::new(BwaMemAligner::new(d, &BwaMemAlignParams::DEFAULT)?)),
            None => None,
        };
        let decoy_writer = alignment_checker::initialise_writer(config);

        let next = AtomicUsize::new(0);
        let results: Mutex<Vec<(usize, ReadStats)>> = Mutex::new(Vec::new());
        let perf: Mutex<PerfCounter> = Mutex::new(PerfCounter::new("PrimaryAssembly"));
        // decoy assemblies and decoy match rows by junction group index
        let decoy_results: Mutex<Vec<(usize, Vec<AsmP>, Vec<String>)>> = Mutex::new(Vec::new());
        let groups = &junction_groups;
        std::thread::scope(|scope| -> Result<()> {
            let mut handles = Vec::new();
            for task in 0..task_count {
                let bam_reader = BamReader::new(config)?;
                let checker = AlignmentChecker::new(config, aligner.clone(), decoy_aligner.clone(), decoy_writer.clone());
                let next = &next;
                let results = &results;
                let decoy_results = &decoy_results;
                let perf = &perf;
                let saga_matcher = self.saga_matcher_factory.as_ref().map(|f| f.create_sequence_matcher());
                handles.push(scope.spawn(move || {
                    let mut assembler = JunctionGroupAssembler::new(config, bam_reader, checker, saga_matcher);
                    let mut local_decoys: Vec<(usize, Vec<AsmP>, Vec<String>)> = Vec::new();
                    let mut perf_counter = PerfCounter::new("PrimaryAssembly");
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= groups.len() {
                            break;
                        }
                        let start = Instant::now();
                        assembler.process_junction_group(groups[i]);
                        perf_counter.stop_check_log(start, || format!("juncGroup({})", *groups[i]));
                        let decoys = assembler.take_decoy_assemblies();
                        let rows = assembler.take_decoy_match_rows();
                        if !decoys.is_empty() || !rows.is_empty() {
                            local_decoys.push((i, decoys, rows));
                        }
                    }
                    results.lock().unwrap().push((task, assembler.read_stats().clone()));
                    decoy_results.lock().unwrap().extend(local_decoys);
                    perf.lock().unwrap().merge(&perf_counter);
                }));
            }
            for h in handles {
                if h.join().is_err() {
                    std::process::exit(1);
                }
            }
            Ok(())
        })?;
        let mut results = results.into_inner().unwrap();
        results.sort_by_key(|r| r.0);
        self.perf_counters.push(perf.into_inner().unwrap());
        // written in junction group order: the order of a single-threaded Java run
        let mut decoy_results = decoy_results.into_inner().unwrap();
        decoy_results.sort_by_key(|r| r.0);
        for (_, _, rows) in &decoy_results {
            alignment_checker::write_match_rows(decoy_writer.as_ref(), rows);
        }

        if self.results.assembly_writer.enabled() {
            for (_, decoys, _) in &decoy_results {
                for decoy_assembly in decoys {
                    decoy_assembly.m().set_outcome_override(AssemblyOutcome::Decoy, true);
                    self.results.assembly_writer.write_assembly(*decoy_assembly);
                }
            }
        }

        let mut assembly_count = 0;
        let mut junction_read_count = 0;
        let mut candidate_read_count = 0;
        for junction_group in &junction_groups {
            for assembly in junction_group.junction_assemblies() {
                assembly_count += 1;
                junction_read_count += assembly.support_count();
                candidate_read_count += assembly.candidate_support().len();
                candidate_read_count += assembly.unmapped_reads().len();
            }
        }
        let mut combined_read_stats = ReadStats::default();
        for (_, stats) in &results {
            combined_read_stats.merge(stats);
        }
        crate::info!("created {} junction assemblies from reads(junc={} candidate={})", assembly_count, junction_read_count, candidate_read_count);
        crate::info!("extracted read stats: {}", combined_read_stats);
        if dev_debug() && combined_read_stats.identical_supplementaries > 0 {
            crate::info!("filtered identical supplementaries({})", combined_read_stats.identical_supplementaries);
        }
        Ok(())
    }

    fn form_phase_groups(&mut self) -> Result<()> {
        let config = &self.config;
        let mut phase_group_builder = PhaseGroupBuilder::new(config, &self.junction_group_map);
        phase_group_builder.build_groups(&mut self.results.phase_group_build_writer);
        self.perf_counters.extend(phase_group_builder.perf_counters().iter().cloned());
        let phase_groups: Vec<PhaseGroupP> = phase_group_builder.phase_groups().clone();
        if phase_groups.is_empty() {
            return Ok(());
        }
        crate::info!("building phase sets from {} phase groups", phase_groups.len());
        let task_count = self.task_count.max(1);
        let next = AtomicUsize::new(0);
        // remote reads matched per phase group, by index
        let matched: Mutex<Vec<(usize, i32)>> = Mutex::new(Vec::new());
        let perf: Mutex<PerfCounter> = Mutex::new(PerfCounter::new("PhaseSets"));
        let groups = &phase_groups;
        std::thread::scope(|scope| -> Result<()> {
            let mut handles = Vec::new();
            for _ in 0..task_count {
                let bam_reader = BamReader::new(config)?;
                let next = &next;
                let matched = &matched;
                let perf = &perf;
                handles.push(scope.spawn(move || {
                    let mut remote_read_extractor = RemoteReadExtractor::new(Some(bam_reader), Some(config));
                    let limit = phase_processing_limit();
                    let mut perf_counter = PerfCounter::new("PhaseSets");
                    let mut last_memory_release = Instant::now();
                    loop {
                        // junction groups free their reads as their assemblies finish here; hand the freed
                        // memory back now and then rather than at the end of the stage
                        if last_memory_release.elapsed().as_millis() >= 500 {
                            crate::assembly::read::read_heap::release_free_memory();
                            last_memory_release = Instant::now();
                        }
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= groups.len() {
                            perf.lock().unwrap().merge(&perf_counter);
                            break;
                        }
                        let phase_group = groups[i];
                        let start = Instant::now();
                        if limit > 0 && phase_group.assembly_count() as i32 > limit {
                            crate::debug!("phase group({}) above processing limited, skipped", *phase_group);
                            continue;
                        }
                        let mut phase_set_builder = PhaseSetBuilder::new(config.ref_genome(), &mut remote_read_extractor, phase_group);
                        phase_set_builder.run();
                        drop(phase_set_builder);
                        matched.lock().unwrap().push((i, remote_read_extractor.remote_reads_matched()));
                        perf_counter.stop_check_log(start, || format!("phaseGroupId({}) assemblies({})", phase_group.id(), phase_group.assembly_count()));
                    }
                }));
            }
            for h in handles {
                if h.join().is_err() {
                    std::process::exit(1);
                }
            }
            Ok(())
        })?;
        self.perf_counters.push(perf.into_inner().unwrap());
        let matched = matched.into_inner().unwrap();
        let total_remote_read_matched = if compat::java("A3") {
            // Java: sums a counter per task that is reset for every phase group, which on one thread is the
            // count of the last phase group processed
            matched.into_iter().max_by_key(|x| x.0).map(|x| x.1).unwrap_or(0)
        } else {
            matched.into_iter().map(|x| x.1).sum()
        };
        crate::info!("created {} phase sets, remote reads extracted({})", phase_groups.iter().map(|x| x.phase_sets().len()).sum::<usize>(), total_remote_read_matched);
        Ok(())
    }

    /// Phase set building has dropped the candidate reads of the assemblies it handled; the reads each
    /// junction group kept for them can now be freed unless something still refers to them.
    fn release_retained_reads(&mut self) {
        // assemblies that phase set building did not clean up (phase groups over the processing limit, the
        // assemblies left out of very large groups, branched assemblies) still hold reads; with these options
        // nothing uses those reads from here on
        if !write_candidate_reads() && clear_cached_read() {
            let mut phase_groups: IdHashSet<PhaseGroupP> = IdHashSet::new();
            for junction_group in self.junction_group_map.values().flatten() {
                for assembly in junction_group.junction_assemblies() {
                    assembly.m().drop_reads_after_phasing();
                    if let Some(phase_group) = assembly.phase_group() {
                        phase_groups.insert(phase_group);
                    }
                }
            }
            for phase_group in phase_groups.iter() {
                for assembly in phase_group.assemblies().iter().chain(phase_group.derived_assemblies().iter()) {
                    assembly.m().drop_reads_after_phasing();
                }
            }
        }
        let mut still_referenced = 0;
        for junction_groups in self.junction_group_map.values() {
            for junction_group in junction_groups {
                still_referenced += junction_group.m().release_retained_reads();
            }
        }
        crate::assembly::read::read_heap::release_free_memory();
        if mem_stats::enabled() {
            crate::info!("memory: reads still referenced after phase sets({})", still_referenced);
        }
    }

    fn gather_assemblies(&mut self) -> (Vec<AsmP>, Vec<AsmAlignP>) {
        let mut all_assemblies: Vec<AsmP> = Vec::new();
        let mut assembly_alignments: Vec<AsmAlignP> = Vec::new();
        let mut phase_groups: IdHashSet<PhaseGroupP> = IdHashSet::new();
        for junction_groups in self.junction_group_map.values() {
            for junction_group in junction_groups {
                for assembly in junction_group.junction_assemblies() {
                    set_assembly_outcome(*assembly);
                    if let Some(pg) = assembly.phase_group() {
                        phase_groups.insert(pg);
                    }
                    all_assemblies.push(*assembly);
                }
            }
        }
        if let Some(factory) = &self.saga_matcher_factory {
            // SAGA matching for assemblies created during phasing, which have not been matched yet
            let mut saga_matcher = factory.create_sequence_matcher();
            for phase_group in phase_groups.iter() {
                for assembly in phase_group.derived_assemblies() {
                    assembly.m().match_to_saga(&mut saga_matcher);
                }
            }
        }
        for phase_group in phase_groups.iter() {
            all_assemblies.extend(phase_group.derived_assemblies().iter().copied());
            for phase_set in phase_group.phase_sets() {
                if phase_set.has_valid_assembly_alignment() {
                    assembly_alignments.push(phase_set.assembly_alignment().unwrap());
                }
            }
            for assembly in phase_group.assemblies() {
                if JunctionAssembly::phase_set(*assembly).is_none() && !skip_unlinked_junction_assembly(*assembly) {
                    assembly_alignments.push(AssemblyAlignment::from_assembly(*assembly));
                }
            }
        }
        java_sort(&mut all_assemblies, |a, b| a.junction().compare_to(&b.junction()));
        for (assembly_id, assembly) in all_assemblies.iter().enumerate() {
            assembly.m().set_id(assembly_id as i32);
            if JunctionAssembly::phase_set(*assembly).is_none() {
                let id = assembly.phase_group().expect("assembly phase group").m().next_phase_set_id();
                assembly.m().set_non_phase_set_id(id);
            }
        }
        for (id, assembly_alignment) in assembly_alignments.iter().enumerate() {
            assembly_alignment.m().set_id(id as i32);
        }
        (all_assemblies, assembly_alignments)
    }

    fn run_alignment(&mut self, assembly_alignments: &[AsmAlignP]) -> Result<()> {
        let config = &self.config;
        let image = config.ref_genome_image_file.as_deref().context("no ref genome image file")?;
        let base = BwaMemAlignParams::DEFAULT;
        let params = BwaMemAlignParams { mismatch_penalty: base.mismatch_penalty + BWA_PENALTY_ADJUST, gap_open_penalty: base.gap_open_penalty + BWA_PENALTY_ADJUST, ..base };
        let aligner = BwaMemAligner::new(image, &params)?;
        let perf_counter = run_alignment(config, &aligner, self.saga_matcher_factory.as_ref(), assembly_alignments);
        self.perf_counters.push(perf_counter);
        Ok(())
    }

    fn write_assembly_output(&mut self, assemblies: &[AsmP]) {
        if self.results.assembly_writer.enabled() || self.results.read_writer.enabled() {
            crate::debug!("writing assembly TSV output");
            for assembly in assemblies {
                self.results.assembly_writer.write_assembly(*assembly);
                self.results.read_writer.write_assembly_reads(*assembly);
            }
        }
    }

    fn write_variants(&mut self, assembly_alignments: &[AsmAlignP], breakends: &[BreakendP]) -> Result<()> {
        if self.config.write_types.contains(&WriteType::Vcf) {
            crate::debug!("writing variant VCF with {} breakends", breakends.len());
            let mut vcf_writer = VcfWriter::new(&self.config)?;
            for breakend in breakends {
                vcf_writer.add_breakend(*breakend);
            }
            vcf_writer.close()?;
        }
        if self.config.write_types.contains(&WriteType::Breakend) {
            for x in assembly_alignments {
                self.results.breakend_writer.write_breakends(*x);
            }
        }
        Ok(())
    }
}

/// Alignment.skipUnlinkedJunctionAssembly
pub fn skip_unlinked_junction_assembly(assembly: AsmP) -> bool {
    if matches!(assembly.outcome(), AssemblyOutcome::DupBranched | AssemblyOutcome::Secondary | AssemblyOutcome::SuppOnly) {
        return true;
    }
    !assembly.is_saga_matched() && is_weak_indel_based_unlinked_assembly(&assembly)
}
