//! COBALT in Rust: a port of hmftools COBALT 3.1 (see ../../hmftools/cobalt).
//!
//! Usage: `cobalt <the Java options>`: reads the same files and writes the same outputs as the Java
//! CobaltApplication, in its three modes (tumor/normal, tumor-only, germline-only) and for targeted panels.
//! `-java_compat` reproduces the Java faults the port corrects (NOTES.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod calc;
mod chromosome;
mod compat;
mod config;
mod constants;
mod count;
mod gc;
#[cfg(test)]
mod e2e_tests;
#[cfg(test)]
mod java_tests;
mod output;
mod stats;
mod targeted;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
// the paths the shared code's macros use (crate::log)
#[allow(unused_imports)]
use hmf_common::{fmt, log};

use anyhow::{bail, Result};

use crate::calc::{CalcConfig, Calculation, ChrRatios, Consolidator, GcMedianReadDepth, MedianRatio, Sample};
use crate::config::{CobaltConfig, ConfigBuilder, Parsed};
use crate::count::{BamJob, ReadFilter, SampleCounts};
use crate::targeted::Scope;

pub const VERSION: &str = "3.1";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(cobalt_main(&argv));
}

/// CobaltApplication.main: the exit status
pub fn cobalt_main(argv: &[String]) -> i32 {
    // the JVM options of the bioconda wrapper (cobalt -Xmx16g ...) are accepted and ignored
    let argv: Vec<String> = argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect();
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let cb = match ConfigBuilder::parse(&argv) {
        Ok(Parsed::Run(cb)) => cb,
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("Cobalt version {VERSION}");
    let config = match CobaltConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            // thrown by Java's CobaltConfig constructor, outside any handler
            error!("{e}");
            return 1;
        }
    };
    if let Err(e) = config.validate() {
        error!("config loading failed: {e}");
        return 1;
    }
    // CO5: Executors.newFixedThreadPool(threads) throws for -threads 0 or less, outside the handler.
    // Corrected: 0 (or less) is taken as 1, as the option's description says.
    if config.threads < 1 && compat::java("CO5") {
        error!("Exception in thread \"main\" java.lang.IllegalArgumentException (Executors.newFixedThreadPool({}))", config.threads);
        return 1;
    }
    let started = std::time::Instant::now();
    info!("reading GC Profile from {}", config.gc_profile_path);
    if let Err(e) = run(&config) {
        error!("error running Cobalt: {e:#}");
        return 1;
    }
    info!("Cobalt complete, mins({})", fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    0
}

fn threads(config: &CobaltConfig) -> usize {
    config.threads.max(1) as usize
}

/// CobaltApplication.run, inside its exception handler
fn run(config: &CobaltConfig) -> Result<()> {
    let Some(specific) = &config.specific_regions else {
        bail!("java.lang.NullPointerException (invalid specific regions)");
    };
    // the BAMs: tumor first, as Java submits its tasks
    let mut jobs: Vec<BamJob> = Vec::new();
    let ref_genome = config.ref_genome_path.as_deref();
    let tumor_job = match &config.tumor_bam {
        Some(path) => {
            jobs.push(BamJob { path, chromosomes: count::load_chromosomes(path, ref_genome, specific)? });
            Some(jobs.len() - 1)
        }
        None => None,
    };
    let reference_job = match &config.reference_bam {
        Some(path) => {
            jobs.push(BamJob { path, chromosomes: count::load_chromosomes(path, ref_genome, specific)? });
            Some(jobs.len() - 1)
        }
        None => None,
    };
    let filter = ReadFilter { min_mapping_quality: config.min_mapping_quality, include_duplicates: config.include_duplicates };
    let phase = std::time::Instant::now();
    let mut counts = count::count_reads(&jobs, ref_genome, &filter, threads(config))?;
    debug!("PerfStats: read counting {:.1}s", phase.elapsed().as_secs_f64());
    let empty = || SampleCounts { chromosomes: Vec::new() };
    let reference_counts = match reference_job {
        Some(i) => std::mem::replace(&mut counts[i], empty()),
        None => empty(),
    };
    let tumor_counts = match tumor_job {
        Some(i) => std::mem::replace(&mut counts[i], empty()),
        None => empty(),
    };
    if tumor_job.is_some() {
        info!("tumor depths({}) collected", tumor_counts.window_count());
    }
    if reference_job.is_some() {
        info!("reference depths({}) collected", reference_counts.window_count());
    }

    let phase = std::time::Instant::now();
    let result = calculate(config, &tumor_counts, &reference_counts)?;
    debug!("PerfStats: ratio calculation {:.1}s", phase.elapsed().as_secs_f64());
    drop(tumor_counts);
    drop(reference_counts);

    let phase = std::time::Instant::now();
    let collated = calc::collate(&result.tumor, &result.reference)?;
    output::write_ratios(&config.ratios_file_name(), config.v38, &collated)?;
    if let Some(tumor) = &config.tumor_id {
        output::write_gc_median(&config.gc_median_file_name(tumor), &result.tumor_stats)?;
    }
    if let Some(reference) = &config.reference_id {
        output::write_median_ratios(&config.median_ratios_file_name(reference), &result.median_ratios)?;
        output::write_gc_median(&config.gc_median_file_name(reference), &result.reference_stats)?;
    }
    debug!("PerfStats: writing ratio files {:.1}s", phase.elapsed().as_secs_f64());
    if !config.skip_pcf_calc {
        let segment = |sample: &str, value: &dyn Fn(&calc::CobaltRatio) -> f64| -> Result<()> {
            let mut points = Vec::with_capacity(collated.iter().map(|c| c.len()).sum());
            for c in &collated {
                let chromosome = c.chr.versioned(config.v38);
                for i in 0..c.len() {
                    let row = c.row(i);
                    points.push(hmf_common::segmentation::SegmentPoint { chromosome: chromosome.clone(), position: row.position, value: value(&row) });
                }
            }
            let phase = std::time::Instant::now();
            hmf_common::segmentation::write_windowed_segments(&points, config.v38, config.pcf_gamma as f64, threads(config), &config.pcf_file_name(sample))?;
            debug!("PerfStats: segmentation of {sample} {:.1}s", phase.elapsed().as_secs_f64());
            Ok(())
        };
        if let Some(tumor) = &config.tumor_id {
            segment(tumor, &|r| r.tumor_gc_ratio)?;
        }
        if let Some(reference) = &config.reference_id {
            segment(reference, &|r| r.reference_gc_diploid_ratio)?;
        }
    }
    output::write_version(config.out_dir())
}

struct CalculatorResult {
    tumor: Vec<ChrRatios>,
    reference: Vec<ChrRatios>,
    median_ratios: Vec<MedianRatio>,
    tumor_stats: GcMedianReadDepth,
    reference_stats: GcMedianReadDepth,
}

fn log_read_depth_info(c: &Calculation, id: &str) {
    info!("{id} sample median({}) mean({})", fmt::java_fixed(c.gc_median.median, 2), fmt::java_fixed(c.gc_median.mean, 2));
}

/// CobaltCalculator
fn calculate(config: &CobaltConfig, tumor_counts: &SampleCounts, reference_counts: &SampleCounts) -> Result<CalculatorResult> {
    if tumor_counts.window_count() == 0 && reference_counts.window_count() == 0 {
        bail!("java.lang.IllegalArgumentException (CobaltCalculator: no read depths)");
    }
    // WindowStatuses: the GC profile, excluded regions and the diploid regions of tumor-only mode
    let phase = std::time::Instant::now();
    let profiles = gc::load_gc_profile(&config.gc_profile_path)?;
    debug!("PerfStats: GC profile {:.1}s", phase.elapsed().as_secs_f64());
    let excluded = gc::excluded_windows(&gc::excluded_regions(config.v38)?, &profiles)?;
    let diploid = match &config.tumor_only_diploid_bed {
        Some(path) => Some(gc::load_diploid_regions(path)?),
        None => None,
    };
    // Java checks the regions only when there are some (an empty file is no file)
    let diploid = diploid.filter(|d| d.0.iter().any(|l| !l.is_empty()));
    let statuses = gc::WindowStatuses::new(profiles, &excluded, diploid.as_ref());
    let (gc_bucket_min, gc_bucket_max) = config.gc_bucket_limits();
    let (gc_ratio_min, gc_ratio_max) = config.gc_ratio_limits();
    let cfg = CalcConfig { gc_bucket_min, gc_bucket_max, gc_ratio_min, gc_ratio_max };
    let has_reference = config.has_reference_id();

    let mut result = CalculatorResult {
        tumor: Vec::new(),
        reference: Vec::new(),
        median_ratios: Vec::new(),
        tumor_stats: GcMedianReadDepth::no_results(),
        reference_stats: GcMedianReadDepth::no_results(),
    };
    if let Some(path) = &config.target_region_norm_file {
        let panels = targeted::load_normalisation_file(path, config.v38)?;
        if panels.is_empty() {
            error!("invalid target region norm file scope state");
            std::process::exit(1);
        }
        let mut consolidator: Option<Consolidator> = None;
        let mut tumor_stats: Option<GcMedianReadDepth> = None;
        let mut median_ratios: Option<Vec<MedianRatio>> = None;
        for (i, panel) in panels.iter().enumerate() {
            let scope = Scope::Targeted(panel);
            let tumor = calc::calculate(tumor_counts, &statuses, &scope, &cfg, Sample::Tumor)?;
            log_read_depth_info(&tumor, &format!("tumor panel {i}"));
            if tumor_stats.is_none() {
                tumor_stats = Some(tumor.gc_median.clone());
            }
            if consolidator.is_none() {
                consolidator = tumor.consolidator;
            }
            calc::merge_panel_ratios(&mut result.tumor, tumor.ratios)?;
            if has_reference {
                let reference = calc::calculate(
                    reference_counts,
                    &statuses,
                    &scope,
                    &cfg,
                    Sample::Reference { consolidator: consolidator.as_mut().unwrap(), diploid_normalise: false, v38: config.v38 },
                )?;
                log_read_depth_info(&reference, &format!("reference panel {i}"));
                calc::merge_panel_ratios(&mut result.reference, reference.ratios)?;
                if median_ratios.is_none() {
                    median_ratios = Some(calc::panel_median_ratios(&mut result.reference, config.v38));
                    result.reference_stats = reference.gc_median;
                }
            }
        }
        result.tumor_stats = tumor_stats.unwrap();
        result.median_ratios = median_ratios.unwrap_or_default();
    } else {
        let scope = Scope::WholeGenome;
        let tumor = calc::calculate(tumor_counts, &statuses, &scope, &cfg, Sample::Tumor)?;
        log_read_depth_info(&tumor, "tumor");
        result.tumor_stats = tumor.gc_median;
        result.tumor = tumor.ratios;
        let mut consolidator = tumor.consolidator.unwrap();
        if has_reference {
            let reference = calc::calculate(
                reference_counts,
                &statuses,
                &scope,
                &cfg,
                Sample::Reference { consolidator: &mut consolidator, diploid_normalise: true, v38: config.v38 },
            )?;
            log_read_depth_info(&reference, "reference");
            result.median_ratios = reference.median_ratios.unwrap_or_default();
            result.reference_stats = reference.gc_median;
            result.reference = reference.ratios;
        }
    }
    Ok(result)
}
