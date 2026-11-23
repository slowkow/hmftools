//! ESVEE depth annotation: a port of the Java package esvee.depth (DepthAnnotator).
//!
//! Reads the breakends of the assembly step's VCF, counts for every breakend and sample the fragments
//! which support the reference (REF: a read over the breakend; REFPAIR: a read on each side), and writes
//! the VCF again with those counts and the allele fraction (AF) added.

use crate::common::compat;
use anyhow::{bail, Result};
use std::sync::atomic::{AtomicBool, AtomicU8, Ordering};
use std::sync::Mutex;

// Log levels as log4j orders them. The crate's shared log macros only know about debug on or off, so
// the step has its own, which honour -log_level.
pub(crate) const LEVEL_NAMES: [&str; 8] = ["OFF", "FATAL", "ERROR", "WARN", "INFO", "DEBUG", "TRACE", "ALL"];
pub(crate) const LEVEL_ERROR: u8 = 2;
pub(crate) const LEVEL_WARN: u8 = 3;
pub(crate) const LEVEL_INFO: u8 = 4;
pub(crate) const LEVEL_DEBUG: u8 = 5;
pub(crate) const LEVEL_TRACE: u8 = 6;
static LOG_LEVEL: AtomicU8 = AtomicU8::new(LEVEL_INFO);

pub(crate) fn set_log_level(level: u8) {
    LOG_LEVEL.store(level, Ordering::Relaxed);
    crate::common::log::DEBUG.store(level >= LEVEL_DEBUG, Ordering::Relaxed);
}

pub(crate) fn log_enabled(level: u8) -> bool {
    level <= LOG_LEVEL.load(Ordering::Relaxed)
}

macro_rules! info {
    ($($arg:tt)*) => { if $crate::depth::log_enabled($crate::depth::LEVEL_INFO) { $crate::common::log::write("INFO ", format_args!($($arg)*)) } };
}
macro_rules! warn {
    ($($arg:tt)*) => { if $crate::depth::log_enabled($crate::depth::LEVEL_WARN) { $crate::common::log::write("WARN ", format_args!($($arg)*)) } };
}
macro_rules! error {
    ($($arg:tt)*) => { if $crate::depth::log_enabled($crate::depth::LEVEL_ERROR) { $crate::common::log::write("ERROR", format_args!($($arg)*)) } };
}
macro_rules! debug {
    ($($arg:tt)*) => { if $crate::depth::log_enabled($crate::depth::LEVEL_DEBUG) { $crate::common::log::write("DEBUG", format_args!($($arg)*)) } };
}
macro_rules! trace {
    ($($arg:tt)*) => { if $crate::depth::log_enabled($crate::depth::LEVEL_TRACE) { $crate::common::log::write("TRACE", format_args!($($arg)*)) } };
}

mod bam;
mod config;
mod fmt;
mod task;
#[cfg(test)]
mod tests;
mod unmapped;
mod vcf;

use config::{form_output_file, DepthConfig};
// for the combined mode (crate::esvee)
pub(crate) use config::parse_java_double;
pub(crate) use vcf::parse_java_int;
#[cfg(test)]
pub(crate) use config::accepted_options;
use task::{DepthTask, ReadSource, ALLELE_FRACTION, REF_DEPTH, REF_DEPTH_PAIR};
use vcf::{VariantContext, VcfHeader, VcfWriter};

const HUMAN_CHROMOSOMES: [&str; 24] =
    ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"];

/// ExcludedRegions.getPolyGRegions
const POLY_G_REGIONS_V37: [(&str, i32, i32); 3] = [("2", 33141260, 33141700), ("4", 41218427, 41218467), ("17", 42646418, 42646458)];
const POLY_G_REGIONS_V38: [(&str, i32, i32); 3] = [("chr2", 32916190, 32916630), ("chr4", 41216410, 41216450), ("chr17", 44569050, 44569090)];

const ALLELE_FRACTION_DESC: &str = "Allele frequency of the breakend";
const REF_DEPTH_DESC: &str = "Count of reads mapping across this breakend";
const REF_DEPTH_PAIR_DESC: &str = "Count of reference read pairs spanning this breakend supporting the reference allele";

/// DepthAnnotator.excludeVariant
fn exclude_variant(config: &DepthConfig, chromosome: &str, position: i32) -> bool {
    let poly_g = if config.is_38 { &POLY_G_REGIONS_V38 } else { &POLY_G_REGIONS_V37 };
    // D3: Java tests the position alone against the poly-G regions, so a variant at those positions is
    // dropped on any chromosome
    let any_chromosome = compat::java("D3");
    if poly_g.iter().any(|(chr, start, end)| (any_chromosome || *chr == chromosome) && position >= *start && position <= *end) {
        return true;
    }
    if !config.specific_regions.is_empty() {
        return !config.specific_regions.iter().any(|(chr, start, end)| chr == chromosome && position >= *start && position <= *end);
    }
    // -specific_chr (D14; the list is empty with Java behaviour)
    !config.specific_chromosomes.is_empty() && !config.specific_chromosomes.iter().any(|chr| chr == chromosome)
}

/// The records to annotate, by chromosome, in the order the chromosomes are processed and written: the
/// reading loop of DepthAnnotator.run with the chromosome loops that follow it. Also whether any record
/// was left after the exclusions (in Java, whether its map of chromosomes has an entry).
fn group_by_chromosome(config: &DepthConfig, records: impl IntoIterator<Item = VariantContext>) -> (Vec<(String, Vec<VariantContext>)>, bool) {
    // D2: the copy Java makes of each record takes getFilters(), which is an empty set and not null for a
    // record whose FILTER is missing, so every unfiltered record is written as PASS
    let missing_filter_is_pass = compat::java("D2");
    // D5: Java starts a new list whenever the chromosome changes and puts it in place of an earlier one,
    // so of a chromosome which appears in two blocks only the last is kept
    let last_block_only = compat::java("D5");
    // D4: Java only goes over the names of the human chromosomes; records on any other contig are dropped
    let human_only = compat::java("D4");

    // each chromosome's blocks of consecutive records, chromosomes in order of first appearance
    let mut blocks: Vec<(String, Vec<Vec<VariantContext>>)> = Vec::new();
    let mut index_of: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut current: Option<usize> = None;
    for mut variant in records {
        if exclude_variant(config, &variant.chromosome, variant.position) {
            continue;
        }
        if current.is_none_or(|i| blocks[i].0 != variant.chromosome) {
            let index = *index_of.entry(variant.chromosome.clone()).or_insert_with(|| {
                blocks.push((variant.chromosome.clone(), Vec::new()));
                blocks.len() - 1
            });
            if last_block_only {
                blocks[index].1.clear();
            }
            blocks[index].1.push(Vec::new());
            current = Some(index);
        }
        if missing_filter_is_pass {
            variant.filters.get_or_insert_with(Vec::new);
        }
        blocks[current.unwrap()].1.last_mut().unwrap().push(variant);
    }

    // the human chromosomes in their own order, then (corrected) the other contigs as they first appeared
    let mut by_chromosome: Vec<Option<(String, Vec<VariantContext>)>> = blocks.into_iter().map(|(chr, blocks)| Some((chr, merge_blocks(blocks)))).collect();
    let mut ordered = Vec::with_capacity(by_chromosome.len());
    for name in HUMAN_CHROMOSOMES {
        if let Some(&index) = index_of.get(&config.versioned_chromosome(name)) {
            ordered.extend(by_chromosome[index].take());
        }
    }
    if !human_only {
        ordered.extend(by_chromosome.into_iter().flatten());
    }
    (ordered, !index_of.is_empty())
}

/// One list from a chromosome's blocks of records. The variants of a chromosome are sliced, and the
/// output is indexed, in order of position, so the blocks are merged by position: each block keeps its
/// own order, and of two records at one position the one from the earlier block comes first.
fn merge_blocks(mut blocks: Vec<Vec<VariantContext>>) -> Vec<VariantContext> {
    if blocks.len() == 1 {
        return blocks.pop().unwrap();
    }
    let mut merged = Vec::with_capacity(blocks.iter().map(|b| b.len()).sum());
    let mut heads: Vec<std::iter::Peekable<std::vec::IntoIter<VariantContext>>> = blocks.into_iter().map(|b| b.into_iter().peekable()).collect();
    loop {
        let mut next: Option<(usize, i32)> = None;
        for (index, head) in heads.iter_mut().enumerate() {
            if let Some(variant) = head.peek() {
                if next.is_none_or(|(_, position)| variant.position < position) {
                    next = Some((index, variant.position));
                }
            }
        }
        let Some((index, _)) = next else { break };
        merged.extend(heads[index].next());
    }
    merged
}

/// The text of the warning when no record is left to annotate. D15: Java's message has one placeholder
/// and two arguments, so it shows the number of records where the file name was meant.
fn all_filtered_message(vcf_count: i32, input_vcf: &str) -> String {
    if compat::java("D15") { format!("all variants filtered from vcf({vcf_count})") } else { format!("all variants filtered from vcf({input_vcf})") }
}

/// The header lines of the tags this step writes: the first part of DepthAnnotator.writeVcf.
fn add_header_lines(config: &DepthConfig, header: &mut VcfHeader) {
    if !header.defs.has_info(ALLELE_FRACTION) {
        header.add_info_line(ALLELE_FRACTION, 1, "Float", ALLELE_FRACTION_DESC);
    }
    if !header.defs.has_format(ALLELE_FRACTION) {
        header.add_format_line(ALLELE_FRACTION, 1, "Float", ALLELE_FRACTION_DESC);
    }
    // D11: with -vcf_tag_prefix Java only adds the lines of the prefixed tags when the input already has
    // the lines of the plain ones, so on a VCF from the assembly step the writer stops at the first record
    let java_prefix_lines = compat::java("D11");
    for (tag, description) in [(REF_DEPTH, REF_DEPTH_DESC), (REF_DEPTH_PAIR, REF_DEPTH_PAIR_DESC)] {
        let had_plain = header.defs.has_format(tag);
        if !had_plain {
            header.add_format_line(tag, 1, "Integer", description);
            header.add_info_line(tag, 1, "Integer", description);
        }
        if config.vcf_tag_prefix.is_none() {
            continue;
        }
        let prefixed = config.vcf_tag(tag);
        if java_prefix_lines {
            if had_plain {
                header.add_format_line(&prefixed, 1, "Integer", description);
                header.add_info_line(&prefixed, 1, "Integer", description);
            }
        } else {
            if !header.defs.has_format(&prefixed) {
                header.add_format_line(&prefixed, 1, "Integer", description);
            }
            if !header.defs.has_info(&prefixed) {
                header.add_info_line(&prefixed, 1, "Integer", description);
            }
        }
    }
}

/// The keys of a java.util.HashMap<Integer, ?> in iteration order, given the order they were first put.
fn java_hash_map_order(keys_in_insertion_order: &[i32]) -> Vec<i32> {
    let mut capacity = 16usize;
    let mut buckets: Vec<Vec<i32>> = vec![Vec::new(); capacity];
    let bucket_of = |key: i32, capacity: usize| ((key ^ ((key as u32) >> 16) as i32) as u32 as usize) & (capacity - 1);
    let mut size = 0usize;
    for &key in keys_in_insertion_order {
        buckets[bucket_of(key, capacity)].push(key);
        size += 1;
        if size > capacity * 3 / 4 {
            capacity *= 2;
            let mut resized: Vec<Vec<i32>> = vec![Vec::new(); capacity];
            for bucket in &buckets {
                for &k in bucket {
                    resized[bucket_of(k, capacity)].push(k);
                }
            }
            buckets = resized;
        }
    }
    buckets.into_iter().flatten().collect()
}

/// DepthAnnotator.analyseVariantDistribution: debug figures on how the variants group
fn analyse_variant_distribution(config: &DepthConfig, chr_variants: &[(String, Vec<VariantContext>)]) {
    let mut frequency_keys: Vec<i32> = Vec::new();
    let mut frequencies: std::collections::HashMap<i32, i32> = std::collections::HashMap::new();
    let mut total_groups = 0;
    let mut total_variants = 0i32;
    let mut estimated_reads = 0i64;

    for (chromosome, variants) in chr_variants {
        let mut solo_variants = 0;
        let mut group_count = 0;
        let mut index = 0;
        while index < variants.len() {
            let mut variant_count = 1i32;
            let pos_start = variants[index].position;
            let mut pos_end = pos_start;
            let mut next_index = index + 1;
            while next_index < variants.len() {
                let next_position = variants[next_index].position;
                if next_position.wrapping_sub(pos_end) > config.proximity_distance {
                    break;
                }
                pos_end = next_position;
                variant_count += 1;
                next_index += 1;
            }
            if !frequencies.contains_key(&variant_count) {
                frequency_keys.push(variant_count);
            }
            *frequencies.entry(variant_count).or_insert(0) += 1;
            group_count += 1;
            total_groups += 1;
            total_variants += variant_count;
            if variant_count == 1 {
                solo_variants += 1;
            }
            estimated_reads += (pos_end - pos_start + 2000) as i64;
            index += variant_count as usize;
        }
        debug!(
            "chr({}) variants({}) group({}) soloVariants({} pct={})",
            chromosome,
            variants.len(),
            group_count,
            solo_variants,
            fmt::java_fixed(solo_variants as f64 / variants.len() as f64, 3)
        );
    }

    let mut large_group_count = 0;
    let mut large_variants_count = 0;
    for key in java_hash_map_order(&frequency_keys) {
        let frequency = frequencies[&key];
        if key >= 25 {
            large_group_count += 1;
            large_variants_count += key * frequency;
        } else {
            debug!("group count({}) frequency({}) total variants({})", key, frequency, key * frequency);
        }
    }
    debug!("large group count({}) total variants({})", large_group_count, large_variants_count);
    debug!("total variants({}) groups({}) estimated reads({})", total_variants, total_groups, estimated_reads);
}

/// DepthAnnotator.writeVcf
fn write_vcf(config: &DepthConfig, header: &mut VcfHeader, tasks: &[DepthTask]) -> Result<()> {
    let Some(sample_id) = config.sample_ids.first() else { bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0") };
    let output_vcf = form_output_file(config.output_dir.as_deref(), sample_id, "ref_depth.vcf.gz", config.output_id.as_deref());
    info!("writing VCF: {}", output_vcf);

    add_header_lines(config, header);

    // the tasks are in the order of the chromosomes (see group_by_chromosome)
    let mut writer = VcfWriter::create(&output_vcf, header)?;
    for task in tasks {
        for variant in &task.variants {
            writer.write(&variant.encode(header)?)?;
        }
    }
    writer.close()
}

struct Job<'a> {
    task: usize,
    group: &'a mut [task::VariantInfo],
}

#[derive(Default)]
struct Progress {
    started: bool,
    remaining_groups: usize,
    total_reads: i32,
    perf: task::PerfCounter,
}

/// TaskExecutor.executeTasks over the DepthTask of each chromosome. Java gives each chromosome to one
/// thread; here the unit handed to a worker is one group of variants (one slice of every BAM), which
/// gives the same counts and lets the threads share a chromosome that has most of the variants. Each
/// worker has its own readers. Returns the first failure.
fn execute_tasks(config: &DepthConfig, tasks: &mut [DepthTask]) -> Result<()> {
    let open_sources = || -> Result<Vec<Box<dyn ReadSource>>> {
        let mut sources: Vec<Box<dyn ReadSource>> = Vec::new();
        for path in &config.bam_files {
            sources.push(Box::new(bam::BamSource::open(path, &config.ref_genome)?));
        }
        Ok(sources)
    };

    let chromosomes: Vec<(String, usize)> = tasks.iter().map(|t| (t.chromosome().to_string(), t.variant_count())).collect();
    let group_sizes: Vec<Vec<usize>> = tasks.iter().map(|t| t.group_ranges().iter().map(|r| r.len()).collect()).collect();
    let progress: Vec<Mutex<Progress>> = group_sizes.iter().map(|sizes| Mutex::new(Progress { remaining_groups: sizes.len(), ..Default::default() })).collect();

    let mut jobs: Vec<Job> = Vec::new();
    for (index, task) in tasks.iter_mut().enumerate() {
        let mut rest = task.infos_mut();
        for &size in &group_sizes[index] {
            let (group, tail) = rest.split_at_mut(size);
            jobs.push(Job { task: index, group });
            rest = tail;
        }
    }
    let job_count = jobs.len();
    let queue = Mutex::new(jobs.into_iter());
    let failed = AtomicBool::new(false);

    let work = || -> Result<()> {
        let mut sources = open_sources()?;
        let mut slicer = DepthTask::new("", config, Vec::new());
        loop {
            if failed.load(Ordering::Relaxed) {
                return Ok(());
            }
            let Some(job) = queue.lock().unwrap().next() else { return Ok(()) };
            let (chromosome, variant_count) = &chromosomes[job.task];
            {
                let mut state = progress[job.task].lock().unwrap();
                if !state.started {
                    state.started = true;
                    info!("chr({}) processing {} variants", chromosome, variant_count);
                }
            }
            let (read_count, seconds) = slicer.slice_group(chromosome, job.group, &mut sources)?;

            let mut state = progress[job.task].lock().unwrap();
            state.total_reads = state.total_reads.wrapping_add(read_count);
            state.perf.add(seconds);
            state.remaining_groups -= 1;
            if state.remaining_groups == 0 {
                let mut processed = 0;
                for size in &group_sizes[job.task] {
                    processed += size;
                    if processed % 1000 == 0 {
                        debug!("chr({}) processed {} variants", chromosome, processed);
                    }
                }
                info!("chr({}) complete for {} variants, total reads({})", chromosome, processed, state.total_reads);
            }
        }
    };

    let workers = (config.threads.max(1) as usize).min(job_count.max(1));
    if workers <= 1 {
        work()?;
    } else {
        let failure: Mutex<Option<anyhow::Error>> = Mutex::new(None);
        std::thread::scope(|scope| {
            for _ in 0..workers {
                scope.spawn(|| {
                    if let Err(e) = work() {
                        failed.store(true, Ordering::Relaxed);
                        failure.lock().unwrap().get_or_insert(e);
                    }
                });
            }
        });
        if let Some(e) = failure.into_inner().unwrap() {
            return Err(e);
        }
    }

    for (task, state) in tasks.iter_mut().zip(progress) {
        task.perf = state.into_inner().unwrap().perf;
        task.set_depth_values()?;
    }
    Ok(())
}

/// DepthAnnotator.main. A failure which is an uncaught exception in Java is reported as Java's default
/// handler reports it, and ends the run with status 1.
pub fn run(argv: &[String]) -> Result<()> {
    if let Err(e) = annotate(argv) {
        error!("thread exception: {:#}", e);
        std::process::exit(1);
    }
    Ok(())
}

fn annotate(argv: &[String]) -> Result<()> {
    let config = DepthConfig::from_args(argv)?;

    if !std::path::Path::new(&config.input_vcf).exists() {
        error!("missing VCF config or file");
        std::process::exit(1);
    }
    if config.bam_files.len() != config.sample_ids.len() {
        error!("inconsistent samples and BAM files");
        std::process::exit(1);
    }

    info!("running depth annotation for samples: [{}]", config.sample_ids.join(", "));
    let start_time = std::time::Instant::now();

    let (mut header, lines) = vcf::read_vcf(&config.input_vcf)?;

    // establishGenotypeIds
    let mut genotype_ids: Vec<usize> = Vec::new();
    for sample_id in &config.sample_ids {
        match header.samples.iter().position(|x| x == sample_id) {
            Some(index) => genotype_ids.push(index),
            None => {
                error!("sample({}) not found in VCF genotype names", sample_id);
                std::process::exit(1);
            }
        }
    }

    let mut vcf_count = 0;
    let mut records: Vec<VariantContext> = Vec::with_capacity(lines.len());
    for line in &lines {
        vcf_count += 1;
        match VariantContext::parse(line, &header) {
            Ok(variant) => records.push(variant),
            Err(e) => {
                error!("error reading vcf({}): {}", config.input_vcf, e);
                std::process::exit(1);
            }
        }
    }
    drop(lines);
    let (chr_variants, any_record) = group_by_chromosome(&config, records);

    info!("loaded {} variants from vcf({})", vcf_count, config.input_vcf);

    if !any_record {
        warn!("{}", all_filtered_message(vcf_count, &config.input_vcf));
        return write_vcf(&config, &mut header, &[]);
    }

    if config.perf_log_time > 0.0 {
        analyse_variant_distribution(&config, &chr_variants);
    }

    let mut tasks: Vec<DepthTask> = Vec::new();
    for (chromosome, variants) in chr_variants {
        let mut task = DepthTask::new(&chromosome, &config, genotype_ids.clone());
        task.add_variants(variants)?;
        tasks.push(task);
    }

    if let Err(e) = execute_tasks(&config, &mut tasks) {
        error!("task execution error: {:#}", e);
        std::process::exit(1);
    }

    if let Some(unmap_regions_file) = &config.unmap_regions_file {
        let unmapped_region_depth = unmapped::UnmappedRegionDepth::load(unmap_regions_file)?;
        let mut variants: Vec<&mut VariantContext> = tasks.iter_mut().flat_map(|t| t.variants.iter_mut()).collect();
        unmapped_region_depth.set_unmapped_regions_depth(&genotype_ids, &config.vcf_tag(REF_DEPTH), &config.vcf_tag(REF_DEPTH_PAIR), &mut variants)?;
    }

    write_vcf(&config, &mut header, &tasks)?;

    info!("Esvee depth annotation complete, mins({})", fmt::java_fixed(start_time.elapsed().as_millis() as f64 / 60000.0, 3));

    if config.perf_log_time > 0.0 {
        let Some((first, rest)) = tasks.split_first_mut() else { bail!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0") };
        for task in rest.iter() {
            first.perf.merge(&task.perf);
        }
        let perf = &first.perf;
        if perf.interval_count > 0 {
            info!(
                "PerfStat(Slice) intervals({}) total({} avg={} max={})",
                perf.interval_count,
                fmt::java_fixed(perf.total_time, 3),
                fmt::java_fixed(perf.total_time / perf.interval_count as f64, 3),
                fmt::java_fixed(perf.max_time, 3)
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod order_tests {
    use super::java_hash_map_order;

    #[test]
    fn hash_map_iteration_order() {
        assert_eq!(java_hash_map_order(&[5, 1, 3]), vec![1, 3, 5]);
        // 17 shares a bucket with 1 in a table of 16
        assert_eq!(java_hash_map_order(&[17, 1, 2]), vec![17, 1, 2]);
        let keys: Vec<i32> = (0..13).rev().chain([40]).collect();
        // the 13th entry grows the table to 32, where 40 falls in bucket 8 after 8
        assert_eq!(java_hash_map_order(&keys), vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 40, 9, 10, 11, 12]);
    }
}
