//! BamMetrics.run: the partitions read by worker threads, the unmapped reads, and the output files
//! (MetricsWriter, TargetRegionStats, OffTargetAnalyser and GeneCoverage writers).

use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{self, Read as _};
use rust_htslib::htslib;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use hmf_common::metrics::{BamFlagCounts, BamMetricCoverage, BamMetricFragmentLength, BamMetricSummary, FlagStatType};
use hmf_common::regions::{compare_chr_regions, partition_chromosome, ChrBaseRegion};
use hmf_common::segmentation::HumanChromosome;

use crate::compat;
use crate::config::MetricsConfig;
use crate::coverage::{coverage_rows, DUPLICATE, LOW_BASE_QUAL, LOW_MAP_QUAL, MATE_UNMAPPED, MAX_COVERAGE, OVERLAPPED};
use crate::genes::GeneCoverage;
use crate::partition::{PartitionReader, Task, WorkerStats};
use crate::record::Rec;
use crate::targets::{ENRICHED_HEADER, TARGET_REGIONS_HEADER};

pub const COVERAGE_LEVELS: [i32; 14] = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100];

fn is_cram(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".cram")
}

/// an indexed reader of the BAM or CRAM, with the reference for a CRAM
fn open_reader(config: &MetricsConfig) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(&config.bam_file).with_context(|| format!("opening {}", config.bam_file))?;
    if is_cram(&config.bam_file) {
        reader.set_reference(&config.ref_genome_file).with_context(|| format!("setting CRAM reference {}", config.ref_genome_file))?;
    }
    if config.file_buffer_kb > 0 {
        // SAFETY: the option takes one int argument
        let ret = unsafe { htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_HTS_OPT_BLOCK_SIZE, (config.file_buffer_kb << 10) as std::os::raw::c_int) };
        if ret != 0 {
            anyhow::bail!("setting the file buffer size for {}", config.bam_file);
        }
    }
    Ok(reader)
}

/// The regions BamMetrics.run reads: each human chromosome's partitions (or those of the specific
/// regions), or with -only_target the target regions, sorted
pub fn build_tasks(config: &MetricsConfig) -> Vec<Task> {
    let specific = &config.specific_chr_regions;
    let mut regions: Vec<ChrBaseRegion> = Vec::new();
    if !config.only_target_regions {
        for c in HumanChromosome::ALL {
            let chr = config.ref_gen_version.versioned_chromosome(c.short_name());
            if specific.exclude_chromosome(&chr) {
                continue;
            }
            regions.extend(partition_chromosome(&chr, config.ref_gen_version, &specific.regions, config.partition_size));
        }
    } else {
        for (chr, list) in config.target_regions.iter() {
            if specific.exclude_chromosome(chr) {
                continue;
            }
            for r in list {
                if specific.include_positions(r.start, r.end) {
                    regions.push(ChrBaseRegion::new(chr, r.start, r.end));
                }
            }
        }
        regions.sort_by(compare_chr_regions);
    }
    let mut tasks: Vec<Task> = Vec::with_capacity(regions.len());
    for (i, r) in regions.iter().enumerate() {
        let previous = if config.only_target_regions && i > 0 && regions[i - 1].chromosome == r.chromosome {
            Some((regions[i - 1].start, regions[i - 1].end))
        } else {
            None
        };
        tasks.push(Task { id: i, region: r.clone(), previous, unmapped: false });
    }
    tasks
}

/// PartitionThread.run: partitions taken in order until none are left
fn run_worker(config: &MetricsConfig, genes: &GeneCoverage, tasks: &[Task], next: &AtomicUsize) -> Result<WorkerStats> {
    let mut reader = open_reader(config)?;
    let mut record = bam::Record::new();
    let mut out = WorkerStats::new(config, genes.exon_count);
    let mut pr = PartitionReader::new(config, genes);
    let mut warned: Vec<String> = Vec::new();
    // BT3: Java takes the QC-fail flag as passing QC, so passing unmapped reads count as failed.
    let java_bt3 = compat::java("BT3");
    loop {
        let i = next.fetch_add(1, Ordering::Relaxed);
        let Some(task) = tasks.get(i) else { break };
        if task.unmapped {
            // the fully unmapped reads, to the flag stats' total
            reader.fetch(bam::FetchDefinition::Unmapped)?;
            while let Some(r) = reader.read(&mut record) {
                r.with_context(|| format!("reading {}", config.bam_file))?;
                let fails = record.flags() & crate::record::FLAG_QC_FAIL != 0;
                out.flag_stats.increment(FlagStatType::Total, if java_bt3 { fails } else { !fails });
            }
            continue;
        }
        if task.id > 0 && task.id % 10 == 0 {
            crate::info!("processing partition({}), remaining({})", task.id, tasks.len().saturating_sub(i + 1));
        }
        crate::debug!("processing region({})", task.region);
        let started = Instant::now();
        let tid = reader.header().tid(task.region.chromosome.as_bytes());
        pr.begin(task, tid.map(|t| t as i32).unwrap_or(-1))?;
        match tid {
            None => {
                if !warned.contains(&task.region.chromosome) {
                    warned.push(task.region.chromosome.clone());
                }
                crate::warn!("cannot find sequence index for chromosome {} in bam header", task.region.chromosome);
            }
            Some(tid) => {
                reader.fetch((tid, (task.region.start - 1) as i64, task.region.end as i64))?;
                while let Some(r) = reader.read(&mut record) {
                    r.with_context(|| format!("reading {}", config.bam_file))?;
                    pr.process(&Rec::of(&record), &mut out).with_context(|| format!("partition {}", task.region))?;
                }
            }
        }
        let elapsed = started.elapsed().as_secs_f64();
        pr.end(&mut out, task.id, elapsed);
    }
    Ok(out)
}

/// BamMetrics.run; the exit status
pub fn run(config: &MetricsConfig) -> Result<i32> {
    crate::info!("sample({}) starting bam metrics", config.sample_id);
    let started = Instant::now();

    // BT19: Java does not create the output directory; without it every writer logs an error and the
    // run ends with status 0 and no files. Corrected: the directory is created.
    if !std::path::Path::new(&config.output_dir).exists() && !compat::java("BT19") {
        std::fs::create_dir_all(&config.output_dir).with_context(|| format!("creating {}", config.output_dir))?;
    }

    let genes = GeneCoverage::load(config)?;

    if !config.target_regions.is_empty() {
        let count: usize = config.target_regions.iter().map(|(_, v)| v.len()).sum();
        let bases: i64 = config.target_regions.iter().map(|(_, v)| v.iter().map(|r| r.base_length() as i64).sum::<i64>()).sum();
        crate::info!("capturing data for {count} target regions, total base count({bases})");
    }
    // BT17: with a fragment overlap threshold but no -write_off_target, Java opens the enriched-region
    // file but has no analyser to take the sites from: a NullPointerException ends the run
    if !config.target_regions.is_empty() && config.high_fragment_overlap_threshold > 0 && !config.write_off_target && compat::java("BT17") {
        return Err(anyhow!("java.lang.NullPointerException: no off-target analyser for the enriched-region file"));
    }

    let mut tasks = build_tasks(config);
    let threads = if tasks.len() <= 1 || config.threads <= 1 { 1 } else { tasks.len().min(config.threads as usize) };
    if threads > 1 {
        crate::info!("splitting {} regions across {} threads", tasks.len(), config.threads);
    }
    // the reads without a position (Java reads them after the partitions, on one thread): a task of
    // their own, taken first; their counts are sums, so the order makes no difference
    if !config.specific_chr_regions.has_filters() {
        let unmapped = Task { id: usize::MAX, region: ChrBaseRegion::new("*", 0, 0), previous: None, unmapped: true };
        tasks.insert(0, unmapped);
    }
    let next = AtomicUsize::new(0);
    let results: Vec<Result<WorkerStats>> = if tasks.is_empty() {
        Vec::new()
    } else {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..threads).map(|_| scope.spawn(|| run_worker(config, &genes, &tasks, &next))).collect();
            handles.into_iter().map(|h| h.join().unwrap_or_else(|_| Err(anyhow!("worker thread panicked")))).collect()
        })
    };
    let mut combined = WorkerStats::new(config, genes.exon_count);
    for r in results {
        combined.merge(r?);
    }
    crate::info!("all regions complete");

    write_per_partition_files(config, &mut combined)?;

    combined.coverage.finalise(config.exclude_zero_coverage);
    write_results(config, &mut combined);

    crate::info!("totalReads({}) stats: {}", combined.read_counts.total, combined.coverage.describe());

    let depths = std::mem::take(&mut combined.exon_depths);
    if let Err(e) = genes.write_results(&config.output_dir, &config.sample_id, &depths) {
        crate::error!("failed to write gene coverage: {e:#}");
    }
    if config.perf_debug {
        crate::info!("PerfStats: {} partitions, {} reads, {:.3} seconds", tasks.len(), combined.reads, started.elapsed().as_secs_f64());
    }
    crate::info!("BamMetrics complete, mins({})", hmf_common::fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(0)
}

fn write_text(path: &str, header: &str, rows: impl Iterator<Item = String>) -> Result<()> {
    let mut s = String::new();
    s.push_str(header);
    s.push('\n');
    for r in rows {
        s.push_str(&r);
        s.push('\n');
    }
    std::fs::write(path, s).with_context(|| format!("writing {path}"))
}

/// The files written partition by partition in Java (in the order the threads finish them, BT9); here
/// in partition order
fn write_per_partition_files(config: &MetricsConfig, combined: &mut WorkerStats) -> Result<()> {
    combined.rows.sort_by_key(|r| r.0);
    if !config.only_target_regions {
        let header = "Chromosome\tPositionStart\tPositionEnd\tTotalReads\tDuplicateReads\tChimericReads\tInterPartition\tUnmappedReads\tProcessTime";
        let path = config.form_filename("partition_stats");
        if let Err(e) = write_text(&path, header, combined.rows.iter().filter_map(|r| r.1.clone())) {
            crate::error!("failed to initialise partition stats file: {e:#}");
        }
    }
    if !config.target_regions.is_empty() {
        let path = format!("{}{}{}.target_regions.tsv", config.output_dir, config.sample_id, hmf_common::metrics::BAM_METRICS_FILE_ID);
        if let Err(e) = write_text(&path, TARGET_REGIONS_HEADER, combined.rows.iter().flat_map(|r| r.2.clone())) {
            crate::error!("failed to initialise target regions writer: {e:#}");
        }
        if config.high_fragment_overlap_threshold > 0 {
            let path = config.form_filename("off_target_enriched");
            if let Err(e) = write_text(&path, ENRICHED_HEADER, combined.rows.iter().flat_map(|r| r.3.clone())) {
                crate::error!("failed to initialise high fragment overlaps writer: {e:#}");
            }
        }
    }
    Ok(())
}

/// MetricsWriter.writeResults: summary, coverage, fragment lengths, flag counts, off-target overlaps
fn write_results(config: &MetricsConfig, combined: &mut WorkerStats) {
    let metrics = &mut combined.coverage;
    let territory = metrics.genome_territory();
    let stats = metrics.statistics();
    let percents: Vec<f64> = COVERAGE_LEVELS.iter().map(|l| metrics.calc_coverage_frequency(*l, territory as f64)).collect();
    let summary = BamMetricSummary {
        total_region_bases: territory,
        total_reads: combined.read_counts.total,
        duplicate_reads: combined.read_counts.duplicates,
        dual_strand_reads: combined.read_counts.dual_strand,
        off_target_reads: combined.read_counts.off_target,
        mean_coverage: stats.mean,
        sd_coverage: stats.std_deviation,
        // (int)Math.round(median)
        median_coverage: hmf_common::java_math::java_round(stats.median) as i32,
        mad_coverage: stats.mad,
        low_map_qual_percent: metrics.calc_filtered_percentage(LOW_MAP_QUAL),
        duplicate_percent: metrics.calc_filtered_percentage(DUPLICATE),
        unmapped_percent: metrics.calc_filtered_percentage(MATE_UNMAPPED),
        low_base_qual_percent: metrics.calc_filtered_percentage(LOW_BASE_QUAL),
        overlapping_read_percent: metrics.calc_filtered_percentage(OVERLAPPED),
        capped_coverage_percent: metrics.calc_filtered_percentage(MAX_COVERAGE),
        coverage_levels: COVERAGE_LEVELS.to_vec(),
        coverage_percents: percents,
    };
    if let Err(e) = summary.write(&BamMetricSummary::generate_filename(&config.output_dir, &config.sample_id)) {
        crate::error!("failed to write metrics file: {e:#}");
    }
    if let Err(e) = BamMetricCoverage::write(&BamMetricCoverage::generate_filename(&config.output_dir, &config.sample_id), &coverage_rows(&metrics.coverage_frequency)) {
        crate::error!("failed to write coverage frequency file: {e:#}");
    }
    if let Err(e) = BamMetricFragmentLength::write(
        &BamMetricFragmentLength::generate_filename(&config.output_dir, &config.sample_id),
        &combined.fragment_lengths.rows(),
    ) {
        crate::error!("failed to write frag lengths file: {e:#}");
    }
    if let Err(e) = BamFlagCounts::write(&BamFlagCounts::generate_filename(&config.output_dir, &config.sample_id), &combined.flag_stats.stats) {
        crate::error!("failed to write flag stats file: {e:#}");
    }
    if config.write_off_target {
        let mut keys: Vec<i32> = combined.off_target_overlaps.keys().copied().collect();
        keys.sort_unstable();
        let path = config.form_filename("off_target_overlaps");
        if let Err(e) = write_text(&path, "FragmentOverlaps\tFrequency", keys.iter().map(|k| format!("{k}\t{}", combined.off_target_overlaps[k]))) {
            crate::error!("failed to writer off-target fragment overlaps: {e:#}");
        }
    }
}
