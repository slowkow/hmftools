//! TEAL in Rust: a port of hmftools TEAL 1.4 (see ../teal).
//!
//! Usage: `teal [<class> | prep | pipeline | app] <the Java options>`, where the class is
//! com.hartwig.hmftools.teal.TealPipelineTelbamApp (`prep`: the telbam of each sample from its BAM),
//! com.hartwig.hmftools.teal.TealPipelineApp (`pipeline`: telomere lengths and break ends from the
//! telbams and the COBALT, PURPLE and BamTools outputs) or com.hartwig.hmftools.teal.TealApplication
//! (`app`, also without a class as the jar's main class: both, with the depths on the command line).
//! JVM options (-Xmx..., -XX:..., -D...=...) are dropped, as the bioconda wrapper passes them to java.
//! `-java_compat` reproduces the Java faults the port corrects (REPORT.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod aligner;
mod bgzf;
mod breakend;
mod compat;
mod config;
mod inputs;
mod prep;
mod read;
mod telbam_reader;
mod tellength;
mod telo;
#[cfg(test)]
mod e2e_tests;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
// the paths the shared code's macros use (crate::log)
#[allow(unused_imports)]
use hmf_common::{fmt, log};

use anyhow::{anyhow, bail, Context, Result};

use crate::config::{App, CommonParams, ConfigBuilder, Parsed, PipelineParams, TealParams, APP_NAME, VERSION};
use crate::tellength::{SampleType, TelomereLengthCalc};

const PREP_CLASS: &str = "com.hartwig.hmftools.teal.TealPipelineTelbamApp";
const PIPELINE_CLASS: &str = "com.hartwig.hmftools.teal.TealPipelineApp";
const APP_CLASS: &str = "com.hartwig.hmftools.teal.TealApplication";

fn is_jvm_option(a: &str) -> bool {
    a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-Xss") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).filter(|a| !is_jvm_option(a)).collect();
    std::process::exit(teal_main(&argv));
}

/// The entry points' main: the exit status
pub fn teal_main(argv: &[String]) -> i32 {
    let (app, rest): (App, Vec<String>) = match argv.first().map(|s| s.as_str()) {
        Some(PREP_CLASS) | Some("prep") => (App::Prep, argv[1..].to_vec()),
        Some(PIPELINE_CLASS) | Some("pipeline") => (App::Pipeline, argv[1..].to_vec()),
        Some(APP_CLASS) | Some("app") => (App::Standalone, argv[1..].to_vec()),
        Some(c) if c.starts_with("com.hartwig.") => {
            eprintln!("teal (Rust port): main class {c} is not ported; only {PREP_CLASS}, {PIPELINE_CLASS} and {APP_CLASS}");
            return 2;
        }
        _ => (App::Standalone, argv.to_vec()),
    };
    let rest = match compat::take_options(&rest) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let cb = match ConfigBuilder::parse(app, &rest) {
        Ok(Parsed::Run(cb)) => cb,
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_value("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("{APP_NAME} version {VERSION}");
    let result = match app {
        App::Prep => run_prep(&cb),
        App::Pipeline => run_pipeline(&cb),
        App::Standalone => TealParams::from(&cb).and_then(|p| run_teal(&p)),
    };
    match result {
        Ok(code) => code,
        Err(e) => {
            error!("{e:#}");
            1
        }
    }
}

fn check_create_output_dir(common: &CommonParams) -> Result<bool> {
    let Some(dir) = &common.output_dir else { bail!("java.lang.NullPointerException (no -output_dir)") };
    if std::path::Path::new(dir).exists() {
        return Ok(true);
    }
    Ok(std::fs::create_dir_all(dir).is_ok())
}

fn threads(common: &CommonParams) -> i32 {
    common.threads.max(1)
}

fn telbam(bam: &str, common: &CommonParams, out: &str) -> Result<()> {
    crate::info!("params: TelbamParams(bamFile={bam}, refGenomeFile={}, telbamFile={out}, threadCount={})", common.ref_genome.as_deref().unwrap_or("null"), common.threads);
    prep::process_bam(&prep::TelbamParams { bam_file: bam, ref_genome: common.ref_genome.as_deref(), telbam_file: out, threads: threads(common) })
}

/// TealPipelineTelbamApp.run
fn run_prep(cb: &ConfigBuilder) -> Result<i32> {
    let common = CommonParams::from(cb)?;
    common.validate()?;
    if !check_create_output_dir(&common)? {
        error!("failed to create output directory({})", common.output_dir.as_deref().unwrap_or("null"));
        return Ok(1);
    }
    info!("Teal version: {VERSION}");
    info!("starting telomeric analysis");
    info!("{common:?}");
    let started = std::time::Instant::now();
    info!("creating telbam files");
    if !common.tumor_only() {
        let bam = common.reference_bam.clone().ok_or_else(|| anyhow!("java.lang.NullPointerException (no -reference_bam)"))?;
        telbam(&bam, &common, &common.germline_telbam_path())?;
    }
    if !common.reference_only() {
        let bam = common.tumor_bam.clone().ok_or_else(|| anyhow!("java.lang.NullPointerException (no -tumor_bam)"))?;
        telbam(&bam, &common, &common.tumor_telbam_path())?;
    }
    let secs = started.elapsed().as_secs();
    info!("Teal run complete, time taken: {}m {}s", secs / 60, secs % 60);
    Ok(0)
}

/// TealPipelineApp.run: the values from the pipeline files, then TealApplication.run
fn run_pipeline(cb: &ConfigBuilder) -> Result<i32> {
    let pp = PipelineParams::from(cb)?;
    pp.validate()?;
    let mut p = TealParams::with_common(pp.common.clone());
    if let (Some(purple), Some(tumor)) = (&pp.purple, &pp.common.tumor_sample_id) {
        let (purity, ploidy) = inputs::read_purity(purple, tumor)?;
        p.tumor_purity = purity;
        p.tumor_ploidy = ploidy;
    }
    if let Some(cobalt) = &pp.cobalt {
        if let Some(tumor) = &pp.common.tumor_sample_id {
            let (mean, gc50) = inputs::read_gc_median(cobalt, tumor)?;
            p.tumor_mean_read_depth = Some(mean);
            p.tumor_gc50_read_depth = Some(gc50);
        }
        if let Some(reference) = &pp.common.reference_sample_id {
            let (mean, gc50) = inputs::read_gc_median(cobalt, reference)?;
            p.germline_mean_read_depth = Some(mean);
            p.germline_gc50_read_depth = Some(gc50);
        }
    }
    if let Some(path) = &pp.tumor_wgs_metrics {
        p.tumor_duplicate_proportion = inputs::read_duplicate_percent(path)?;
        info!("Loaded tumor WGS metrics from {path}");
    }
    if let Some(path) = &pp.reference_wgs_metrics {
        p.germline_duplicate_proportion = inputs::read_duplicate_percent(path)?;
        info!("Loaded reference WGS metrics from {path}");
    }
    info!("loaded teal params {p:?} from pipeline files");
    run_teal(&p)
}

/// The GC 50 read depth to use: Java uses COBALT's -1 for a missing bucket (TE9) and has none when the
/// option is not given (TE8); corrected, both fall back to the mean read depth, as TelLengthApp does
fn gc50(value: Option<f64>, mean: f64) -> Result<f64> {
    match value {
        Some(v) if v >= 0.0 || crate::compat::java("TE9") => Ok(v),
        Some(_) => {
            info!("no gc50 read depth in the COBALT file, using meanDepth({mean})");
            Ok(mean)
        }
        None if crate::compat::java("TE8") => bail!("java.lang.NullPointerException (no gc50 read depth)"),
        None => {
            info!("no gc50_read_depth provided, using meanDepth({mean})");
            Ok(mean)
        }
    }
}

/// TelLengthApp.calcTelomereLength: returns the raw (sample mix) length
#[allow(clippy::too_many_arguments)]
fn tel_length(
    common: &CommonParams,
    sample: &str,
    sample_type: SampleType,
    telbam_file: &str,
    output: &str,
    purity: f64,
    ploidy: f64,
    duplicate_proportion: f64,
    mean: f64,
    gc50_depth: Option<f64>,
    germline_length: Option<f64>,
) -> Result<f64> {
    let gc50_depth = gc50(gc50_depth, mean)?;
    info!("processing telbam: {telbam_file}");
    let groups = telbam_reader::read_telbam(telbam_file, &[], None)?;
    let mut calc = TelomereLengthCalc::new(common.sequencing_type, purity, ploidy, duplicate_proportion, mean, gc50_depth, germline_length);
    calc.process_read_groups(groups.iter());
    calc.write_tsv(output, sample, sample_type)?;
    Ok(calc.telomere_length())
}

/// TealApplication.run
fn run_teal(p: &TealParams) -> Result<i32> {
    p.validate()?;
    let c = &p.common;
    if !check_create_output_dir(c)? {
        error!("failed to create output directory({})", c.output_dir.as_deref().unwrap_or("null"));
        return Ok(1);
    }
    info!("starting telomeric analysis");
    info!("{p:?}");
    let started = std::time::Instant::now();
    let tumor_only = c.tumor_only();
    let germline_only = c.reference_only();
    info!("creating telbam files");
    if !tumor_only {
        let bam = c.reference_bam.clone().ok_or_else(|| anyhow!("java.lang.NullPointerException (no -reference_bam)"))?;
        telbam(&bam, c, &c.germline_telbam_path())?;
    }
    if !germline_only {
        let bam = c.tumor_bam.clone().ok_or_else(|| anyhow!("java.lang.NullPointerException (no -tumor_bam)"))?;
        telbam(&bam, c, &c.tumor_telbam_path())?;
    }
    let mut germline_length: Option<f64> = None;
    if !tumor_only {
        let sample = c.reference_sample_id.as_deref().unwrap_or("null");
        let mean = p.germline_mean_read_depth.context("java.lang.NullPointerException (no reference mean read depth)")?;
        germline_length = Some(tel_length(
            c,
            sample,
            SampleType::Ref,
            &c.germline_telbam_path(),
            &c.germline_tellength_path(),
            1.0,
            2.0,
            p.germline_duplicate_proportion,
            mean,
            p.germline_gc50_read_depth,
            None,
        )?);
    }
    if !germline_only {
        let sample = c.tumor_sample_id.as_deref().unwrap_or("null");
        let mean = p.tumor_mean_read_depth.context("java.lang.NullPointerException (no tumor mean read depth)")?;
        tel_length(
            c,
            sample,
            SampleType::Tumor,
            &c.tumor_telbam_path(),
            &c.tumor_tellength_path(),
            p.tumor_purity,
            p.tumor_ploidy,
            p.tumor_duplicate_proportion,
            mean,
            p.tumor_gc50_read_depth,
            germline_length,
        )?;
    }
    if !germline_only && !tumor_only {
        let mut be = breakend::BreakEndParams::new(c.v37);
        be.sample_id = c.tumor_sample_id.clone().unwrap_or_default();
        be.tumor_telbam = c.tumor_telbam_path();
        be.germline_telbam = c.germline_telbam_path();
        be.output_file = c.breakend_path();
        breakend::find_break_ends(&be)?;
    }
    info!("Teal run complete, mins({})", fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(0)
}
