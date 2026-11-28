//! qsee in Rust: a port of hmftools qsee 1.0 (see ../qsee).
//!
//! Usage: `qsee-rs <the Java options>`: QseeApplication, which collects QC features from PURPLE,
//! BamTools, COBALT, ESVEE and REDUX outputs, writes `<tumor>.qsee.vis.data.tsv.gz`,
//! `<tumor>.qsee.status.tsv.gz` and `<tumor>.qsee.exon_coverage_qc.tsv.gz` as Java does, and draws the visual
//! report (as the R script would) unless `-skip_plot` is given. The jar's other entry points are taken by their
//! class name (full or simple) as the first argument: `CohortPercentilesTrainer`, `DefaultThresholdsWriter`,
//! `ExonCoverageQcFileBuilder`, `ExonExclusionsFileBuilder` and `QseePlot`. `-java_compat` reproduces the Java
//! faults the port corrects (REPORT.md).

mod app;
mod cohort;
mod common;
mod compat;
mod config;
mod coverage;
mod feature;
mod plot;
mod prep;
mod status;
mod vis;
#[cfg(test)]
mod tests;

use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::error;

pub const VERSION: &str = "1.0";
pub const APP_NAME: &str = "Qsee";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(qsee_main(&argv));
}

/// The JVM options of the bioconda wrapper (`qsee -Xmx12240656794 ...`), accepted and ignored
pub fn strip_jvm_options(argv: &[String]) -> Vec<String> {
    argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Entry {
    App,
    Trainer,
    Thresholds,
    /// ExonCoverageQcFileBuilder (GeneLodCalculator before upstream bdceb2302a)
    ExonCoverageQc,
    /// ExonExclusionsFileBuilder (upstream bdceb2302a)
    ExonExclusions,
    Plot,
}

/// the entry point named by the first argument, as `java -cp qsee.jar <class>` selects it
fn entry_of(first: Option<&String>) -> Option<Entry> {
    let name = first?.as_str();
    let simple = name.rsplit('.').next().unwrap_or(name);
    match simple {
        "QseeApplication" => Some(Entry::App),
        "CohortPercentilesTrainer" => Some(Entry::Trainer),
        "DefaultThresholdsWriter" => Some(Entry::Thresholds),
        "ExonCoverageQcFileBuilder" => Some(Entry::ExonCoverageQc),
        "ExonExclusionsFileBuilder" => Some(Entry::ExonExclusions),
        "QseePlot" => Some(Entry::Plot),
        _ => None,
    }
}

/// the exit status: 0, or 1 where Java exits with 1 (System.exit(1) or an uncaught exception)
fn finish(result: anyhow::Result<()>) -> i32 {
    match result {
        Ok(()) => 0,
        Err(e) => {
            if e.downcast_ref::<app::Exit>().is_none() {
                error!("{e:#}");
            }
            1
        }
    }
}

pub fn qsee_main(argv: &[String]) -> i32 {
    let argv = strip_jvm_options(argv);
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let (entry, args) = match entry_of(argv.first()) {
        Some(e) => (e, argv[1..].to_vec()),
        None => (Entry::App, argv),
    };
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    match entry {
        Entry::App => config::register(&mut cb, true, false),
        Entry::Trainer => config::register(&mut cb, false, true),
        Entry::ExonCoverageQc => config::register(&mut cb, false, false),
        Entry::ExonExclusions => config::register_with(&mut cb, false, false, true),
        Entry::Plot => config::register_plot(&mut cb),
        Entry::Thresholds => {
            cb.add_path("output_dir", false, "output_dir");
            cb.add_flag("no_overridable_only", "Also write thresholds that are not overridable");
            cb.add_flag("targeted_mode", "Run in targeted mode (use different threshold defaults for targeted mode)");
        }
    }
    match cb.parse(&args) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    }
    if let Err(e) = config::check_log_level(&cb) {
        error!("{e}");
        return 1;
    }
    finish(match entry {
        Entry::App => config::QseeConfig::from(&cb).and_then(|c| app::run(&c, &cb)),
        Entry::ExonCoverageQc => config::QseeConfig::from(&cb).and_then(|c| coverage::run(&c).map(|_| ())),
        Entry::ExonExclusions => config::QseeConfig::from(&cb).and_then(|c| {
            coverage::run_exclusions_builder(&c, cb.decimal("min_low_coverage_lod_vaf"), cb.decimal("min_low_coverage_cohort_fraction"))
        }),
        Entry::Plot => config::QseePlotConfig::from(&cb).and_then(|pc| vis::run_plot(&pc, &config::log_level_name(&cb))),
        Entry::Trainer => (|| {
            let c = config::QseeConfig::from(&cb)?;
            let t = cohort::TrainConfig {
                num_percentiles: cb.integer("num_percentiles"),
                percentile_increment: cb.decimal("percentile_increment"),
                write_cohort_features: cb.has_flag("write_cohort_features"),
            };
            cohort::train(&c, &t).map_err(app::prep_error)
        })(),
        Entry::Thresholds => {
            let overridable_only = !cb.has_flag("no_overridable_only");
            // QS2 (Java wrote the targeted thresholds without -targeted_mode) is fixed upstream (bdceb2302a), so
            // both modes follow the flag
            let targeted = cb.has_flag("targeted_mode");
            status::write_default_thresholds(&cb.value_or("output_dir", ""), overridable_only, targeted)
        }
    })
}
