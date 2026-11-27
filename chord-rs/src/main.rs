//! CHORD in Rust: a port of hmftools CHORD 2.1.2 (see ../chord), with the R prediction
//! (chord_predict.R and its randomForest model CHORD.rds) done in Rust, so that no R is needed.
//!
//! Usage: `chord-rs <the Java options>` as ChordApplication: counts each sample's mutation contexts
//! from PURPLE's somatic small-variant and SV VCFs, writes `<sample>.chord.mutation_contexts.tsv`, and
//! predicts homologous recombination deficiency into `<sample>.chord.prediction.tsv`.
//! `chord-rs prep <options>` (ChordDataPrep) only counts; `chord-rs predict -mut_contexts_file <tsv>
//! -output_file <tsv> [-chord_model_file <rds>]` (ChordModel) only predicts. The Java class names work
//! in place of `prep` and `predict`. `-java_compat` reproduces the Java faults the port corrects
//! (REPORT.md, CH1-CH10).

mod compat;
mod config;
mod indel;
mod output;
mod predict;
mod prep;
mod refgenome;
mod snv;
mod sv;
#[cfg(test)]
mod tests;
mod vcf_file;

use anyhow::{bail, Result};

use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::{error, info};

use crate::config::{ChordConfig, PredictConfig, APP_NAME, VERSION};

/// log4j's standard levels
const LOG4J_LEVELS: [&str; 8] = ["OFF", "FATAL", "ERROR", "WARN", "INFO", "DEBUG", "TRACE", "ALL"];
/// the levels chord_predict.R knows (LOG_LEVEL); Java passes it the logger's level
const R_LOG_LEVELS: [&str; 6] = ["TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"];

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let status = std::panic::catch_unwind(|| chord_main(&argv)).unwrap_or(1);
    std::process::exit(status);
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// ChordApplication: prep and predict
    Full,
    /// ChordDataPrep.main
    Prep,
    /// ChordModel.main (PredictConfig)
    Predict,
}

/// The exit status of a run with these arguments
pub fn chord_main(argv: &[String]) -> i32 {
    // the JVM options of the bioconda wrapper (chord -Xmx12g ...) are accepted and ignored
    let mut argv: Vec<String> = argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect();
    let mode = match argv.first().map(String::as_str) {
        Some("prep") | Some("com.hartwig.hmftools.chord.prep.ChordDataPrep") => Mode::Prep,
        Some("predict") | Some("com.hartwig.hmftools.chord.predict.ChordModel") => Mode::Predict,
        Some("com.hartwig.hmftools.chord.ChordApplication") => Mode::Full,
        _ => {
            if argv.first().is_some_and(|a| !a.starts_with('-')) {
                eprintln!("unknown command {}: expected prep, predict or the CHORD options", argv[0]);
                return 1;
            }
            Mode::Full
        }
    };
    if mode != Mode::Full || argv.first().is_some_and(|a| a == "com.hartwig.hmftools.chord.ChordApplication") {
        argv.remove(0);
    }
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    if mode == Mode::Predict {
        config::register_predict_config(&mut cb);
    } else {
        config::register_config(&mut cb);
    }
    match cb.parse(&argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    }
    // ConfigUtils.setLogLevel: log4j's Level.valueOf rejects other names
    if mode != Mode::Predict && cb.has_value("log_level") {
        let level = cb.value("log_level").unwrap_or_default().trim().to_ascii_uppercase();
        if !LOG4J_LEVELS.contains(&level.as_str()) {
            error!("java.lang.IllegalArgumentException: Unknown level constant [{level}].");
            return 1;
        }
    }
    let result = if mode == Mode::Predict {
        run_predict(&PredictConfig::from(&cb))
    } else {
        match ChordConfig::from(&cb) {
            Ok(c) => run(&c, mode == Mode::Prep),
            Err(e) => {
                error!("{e:#}");
                return 1;
            }
        }
    };
    match result {
        Ok(()) => 0,
        Err(e) => {
            error!("{e:#}");
            1
        }
    }
}

/// ChordApplication.run (and ChordDataPrep.run with `prep_only`)
fn run(config: &ChordConfig, prep_only: bool) -> Result<()> {
    if config.sample_ids.is_empty() {
        // ChordOutput takes the first sample id for the file names before anything else
        bail!("Failed to run CHORD: java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0");
    }
    let rows = prep::run(config)?;
    if prep_only {
        return Ok(());
    }
    if compat::java("CH2") && !R_LOG_LEVELS.contains(&config.log_level.as_str()) {
        // chord_predict.R looks the level up in its LOG_LEVEL list: ALL and OFF are not there, and its
        // first log message stops the script ("argument is of length zero")
        info!("Starting CHORD predict");
        error!("Error in if (log_level$severity >= LOG_LEVEL[[GLOBAL_LOG_LEVEL]]$severity) message(log_message) : argument is of length zero");
        error!("Failed to get predictions from CHORD model: java.io.IOException: R execution failed for script: chord_predict.R");
        bail!("prediction failed");
    }
    let table = predict::ContextTable::from_counts(&rows);
    let predictions_file = output::predictions_file(config);
    predict::predict_to_file(&table, config.chord_model_file.as_deref(), &predictions_file).map_err(|e| {
        error!("Failed to get predictions from CHORD model: {e:#}");
        e
    })
}

/// ChordModel.main
fn run_predict(config: &PredictConfig) -> Result<()> {
    if config.chord_model_file.is_none() && compat::java("CH7") {
        // Java copies the embedded model next to the output file: new File(output).getParent() is null
        // for a file name without a directory
        let parent = std::path::Path::new(&config.output_file).parent().filter(|p| !p.as_os_str().is_empty());
        if parent.is_none() {
            info!("Starting CHORD predict");
            error!("Failed to get predictions from CHORD model: java.nio.file.NoSuchFileException: null/CHORD.tmp.rds");
            bail!("prediction failed");
        }
    }
    let table = predict::ContextTable::read(&config.mut_contexts_file)?;
    predict::predict_to_file(&table, config.chord_model_file.as_deref(), &config.output_file)
}
