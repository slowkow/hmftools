//! Virus Interpreter in Rust: a port of hmftools Virus Interpreter 1.7.1 (see ../virus-interpreter).
//!
//! Usage: `virus-interpreter-rs <the Java options>`: annotates VIRUSBreakend's summary of the viruses found
//! in a tumor with their names, interpretation and expected clonal coverage, and decides which are reported,
//! writing `<sample>.virus.annotated.tsv` as the Java VirusInterpreterApplication does.
//! `-java_compat` reproduces the Java faults the port corrects (REPORT.md, VI1-VI3).

mod algo;
mod compat;
mod files;
#[cfg(test)]
mod tests;

use std::collections::HashSet;

use anyhow::{bail, Context, Result};
use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::fmt::java_double_to_string;
use hmf_common::metrics::BamMetricSummary;
use hmf_common::{error, info};

pub const VERSION: &str = "1.7.1";
pub const APP_NAME: &str = "VirusInterpreter";

const SAMPLE: &str = "sample";
const PURPLE_DIR: &str = "purple_dir";
const VIRUS_BREAKEND_TSV: &str = "virus_breakend_tsv";
const TAXONOMY_DB_TSV: &str = "taxonomy_db_tsv";
const VIRUS_REPORTING_DB_TSV: &str = "virus_reporting_db_tsv";
const VIRUS_BLACKLISTING_DB_TSV: &str = "virus_blacklisting_db_tsv";
const BAM_METRICS_TUMOR_DIR: &str = "bam_metrics_tumor_dir";
/// the name oncoanalyser passes (Hartwig's virus-interpreter 1.7.2-beta.1 image registers it)
const TUMOR_METRICS_DIR: &str = "tumor_metrics_dir";
const OUTPUT_DIR: &str = "output_dir";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(virus_interpreter_main(&argv));
}

/// VirusInterpreterConfig.registerConfig
fn register(cb: &mut ConfigBuilder) {
    cb.add_str(SAMPLE, true, "Sample ID");
    cb.add_path(PURPLE_DIR, true, "Path to Purple pipeline files");
    cb.add_path(VIRUS_BREAKEND_TSV, true, "Path towards the virus breakend TSV");
    cb.add_path(TAXONOMY_DB_TSV, true, "Path towards a TSV containing a mapping from taxid to taxonomy name");
    cb.add_path(VIRUS_REPORTING_DB_TSV, true, "Path towards a TSV containing reported of viruses");
    cb.add_path(VIRUS_BLACKLISTING_DB_TSV, true, "Path towards a TSV containing viruses to blacklist");
    cb.add_path(BAM_METRICS_TUMOR_DIR, true, "Path to BamMetrics tumor pipeline files");
    cb.add_str(OUTPUT_DIR, false, "Output directory");
    cb.add_logging_options();
}

/// The JVM options of the bioconda wrapper are dropped, and the pipeline's `-tumor_metrics_dir` is read as
/// the jar's `-bam_metrics_tumor_dir`
fn prepare_args(argv: &[String]) -> Vec<String> {
    argv.iter()
        .filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:")))
        .map(|a| if a == &format!("-{TUMOR_METRICS_DIR}") { format!("-{BAM_METRICS_TUMOR_DIR}") } else { a.clone() })
        .collect()
}

pub fn virus_interpreter_main(argv: &[String]) -> i32 {
    let argv = prepare_args(argv);
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    register(&mut cb);
    match cb.parse(&argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    }
    match run(&cb) {
        Ok(()) => 0,
        Err(e) => {
            error!("{e:#}");
            1
        }
    }
}

/// FileWriterUtils.checkAddDirSeparator
fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// File.getParent: the path up to its last separator, or "null" if it has none
fn file_parent(path: &str) -> String {
    match path.rfind('/') {
        Some(0) => "/".to_string(),
        Some(i) => path[..i].to_string(),
        None => "null".to_string(),
    }
}

fn run(cb: &ConfigBuilder) -> Result<()> {
    let sample = cb.value(SAMPLE).unwrap();
    let purple_dir = cb.value(PURPLE_DIR).unwrap();
    let breakend_tsv = cb.value(VIRUS_BREAKEND_TSV).unwrap();
    let taxonomy_tsv = cb.value(TAXONOMY_DB_TSV).unwrap();
    let reporting_tsv = cb.value(VIRUS_REPORTING_DB_TSV).unwrap();
    let blacklisting_tsv = cb.value(VIRUS_BLACKLISTING_DB_TSV).unwrap();
    let metrics_dir = cb.value(BAM_METRICS_TUMOR_DIR).unwrap();
    // parseOutputDir; Java writes to "null<sample>..." in the working directory without -output_dir (VI1)
    let output_prefix = match cb.value(OUTPUT_DIR) {
        Some(d) => check_add_dir_separator(&d),
        None if compat::java("VI1") => "null".to_string(),
        None => String::new(),
    };

    // Java loads the whole taxonomy database first; the port reads the summary first so as to keep only the
    // names it looks up (VI3), then reads the files in Java's order of checks
    let breakends = files::read_virus_breakends(&breakend_tsv);
    let wanted: Option<HashSet<i32>> = breakends.as_ref().ok().map(|b| b.iter().map(|v| v.reference_taxid).collect());

    info!("Loading taxonomy db from {}", taxonomy_tsv);
    let taxonomy = files::load_taxonomy_db(&taxonomy_tsv, wanted.as_ref())?;

    info!("Loading virus blacklisting db from {}", blacklisting_tsv);
    let blacklisted = files::load_blacklisting_db(&blacklisting_tsv)?;

    info!("Building virus reporting db model from {}", reporting_tsv);
    let reporting = files::load_reporting_db(&reporting_tsv)?;

    info!("Loading virus breakends from {}", file_parent(&breakend_tsv));
    let breakends = breakends?;
    info!("Loaded {} virus breakends from {}", breakends.len(), breakend_tsv);

    info!("Loading purity context from Purple dir {}", purple_dir);
    let qc = hmf_common::purple_files::read_purple_qc(&hmf_common::purple_files::purple_qc_filename(&purple_dir, &sample))?;
    let purity = hmf_common::purple_files::read_purity_context(&purple_dir, &sample)?;

    info!("Running coverage analysis based on tumor WGS metrics file {}", metrics_dir);
    let metrics_file = BamMetricSummary::generate_filename(&metrics_dir, &sample);
    let metrics = BamMetricSummary::read(&metrics_file).with_context(|| format!("java.nio.file.NoSuchFileException: {metrics_file}"))?;
    let expected = algo::expected_clonal_coverage(purity.purity, purity.ploidy, metrics.mean_coverage);
    info!("Determined the expected clonal coverage to be {}", java_double_to_string(expected));

    let algo = algo::VirusInterpreterAlgo {
        taxonomy_db: &taxonomy,
        blacklisted_taxids: &blacklisted,
        reporting_db: &reporting,
        expected_clonal_coverage: expected,
    };
    let viruses = algo.analyze(&breakends, &qc.status)?;
    info!("Interpreter classified {} viruses as reportable", viruses.iter().filter(|v| v.reported).count());

    let filename = files::annotated_virus_filename(&output_prefix, &sample);
    if let Err(e) = std::fs::write(&filename, files::annotated_virus_text(&viruses)) {
        bail!("java.nio.file.NoSuchFileException: {filename} ({e})");
    }
    info!("Written {} annotated viruses to {}", viruses.len(), filename);
    Ok(())
}
