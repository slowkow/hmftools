//! PEACH in Rust: a port of hmftools PEACH 2.0 (see ../peach).
//!
//! Usage: `peach-rs <the Java options>`: calls the pharmacogenomic haplotypes of the panel's genes from the
//! germline variants of a sample, and writes the events, the events per gene, every haplotype combination,
//! the best combination with its function and linked drugs, and the QC status per gene, as the Java
//! PeachApplication does. `-java_compat` reproduces the Java faults the port corrects (REPORT.md, PE1-PE2).

mod algo;
mod compat;
mod loaders;
mod output;
mod panel;
#[cfg(test)]
mod tests;

use anyhow::{bail, Context, Result};
use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::{error, info};

pub const VERSION: &str = "2.0";
pub const APP_NAME: &str = "Peach";

const VCF_FILE: &str = "vcf_file";
const SAMPLE_NAME: &str = "sample_name";
const HAPLOTYPES_FILE: &str = "haplotypes_file";
const DRUGS_FILE: &str = "drugs_file";
const FUNCTION_FILE: &str = "function_file";
const OUTPUT_DIR: &str = "output_dir";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(peach_main(&argv));
}

/// PeachConfig.addOptions (then addLoggingOptions again, which registers nothing new)
pub fn register(cb: &mut ConfigBuilder) {
    cb.add_path(VCF_FILE, true, "VCF input file");
    cb.add_path(HAPLOTYPES_FILE, true, "Haplotype config file");
    cb.add_str(SAMPLE_NAME, true, "Name of sample in VCF to call haplotypes for");
    cb.add_path(DRUGS_FILE, false, "Config file of relevant drugs");
    cb.add_path(FUNCTION_FILE, false, "Config file for haplotype function");
    cb.add_str(OUTPUT_DIR, false, "Output directory");
    cb.add_logging_options();
}

/// PeachConfig
pub struct PeachConfig {
    pub vcf_file: String,
    pub sample_name: String,
    pub haplotypes_file: String,
    /// with a trailing separator
    pub output_dir: String,
    pub drugs_file: Option<String>,
    pub function_file: Option<String>,
}

impl PeachConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<PeachConfig> {
        let (Some(vcf_file), Some(haplotypes_file), Some(sample_name), Some(output_dir)) =
            (cb.value(VCF_FILE), cb.value(HAPLOTYPES_FILE), cb.value(SAMPLE_NAME), cb.value(OUTPUT_DIR))
        else {
            bail!("java.lang.IllegalArgumentException: invalid config");
        };
        Ok(PeachConfig {
            vcf_file,
            sample_name,
            haplotypes_file,
            output_dir: check_add_dir_separator(&output_dir),
            drugs_file: if cb.has_value(DRUGS_FILE) { cb.value(DRUGS_FILE) } else { None },
            function_file: if cb.has_value(FUNCTION_FILE) { cb.value(FUNCTION_FILE) } else { None },
        })
    }
    fn path(&self, suffix: &str) -> String {
        format!("{}{}{}", check_add_dir_separator(&self.output_dir), self.sample_name, suffix)
    }
    pub fn events_output_path(&self) -> String {
        self.path(".peach.events.tsv")
    }
    pub fn events_per_gene_output_path(&self) -> String {
        self.path(".peach.gene.events.tsv")
    }
    pub fn all_haplotype_combinations_output_path(&self) -> String {
        self.path(".peach.haplotypes.all.tsv")
    }
    pub fn best_haplotype_combinations_output_path(&self) -> String {
        self.path(".peach.haplotypes.best.tsv")
    }
    pub fn qc_status_output_path(&self) -> String {
        self.path(".peach.qc.tsv")
    }
}

/// FileWriterUtils.checkAddDirSeparator
pub fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

pub fn peach_main(argv: &[String]) -> i32 {
    let argv: Vec<String> = argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect();
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
    let config = match PeachConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e:#}");
            return 1;
        }
    };
    match run(&config) {
        Ok(()) => {
            info!("finished running PEACH");
            0
        }
        Err(e) => {
            error!("unrecoverable error encountered: {e:#}");
            1
        }
    }
}

fn run(config: &PeachConfig) -> Result<()> {
    info!("running PEACH");
    info!("creating output directory");
    std::fs::create_dir_all(&config.output_dir)
        .with_context(|| format!("java.lang.RuntimeException: could not create output directory: {}", config.output_dir))?;

    info!("load config files");
    info!("load haplotype config");
    let panel = panel::load_haplotype_panel(&config.haplotypes_file)?;

    let drugs = match &config.drugs_file {
        None => {
            info!("skip loading drugs since input file not provided");
            None
        }
        Some(f) => {
            info!("load drugs");
            Some(loaders::DrugInfoStore { drug_infos: loaders::load_drug_infos(f)? })
        }
    };
    let functions = match &config.function_file {
        None => {
            info!("skip loading haplotype functions since input file not provided");
            None
        }
        Some(f) => {
            info!("load haplotype functions");
            Some(loaders::HaplotypeFunctionStore { functions: loaders::load_functions(f)? })
        }
    };

    info!("load events");
    let events = loaders::load_relevant_events(&config.vcf_file, &config.sample_name, &panel.relevant_variant_positions())
        .with_context(|| format!("java.lang.RuntimeException: failed to read VCF file: {}", config.vcf_file))?;

    info!("call haplotypes");
    let analyses = algo::gene_to_haplotype_analysis(&panel, &events)?;

    // the files in Java's order, so that a failure leaves the same files written
    let r = (|| -> Result<()> {
        info!("write events output file");
        output::write_lines(&config.events_output_path(), &output::events_lines(&events))?;
        info!("write events per gene output file");
        output::write_lines(&config.events_per_gene_output_path(), &output::events_per_gene_lines(&analyses))?;
        info!("write all haplotype combinations output file");
        output::write_lines(&config.all_haplotype_combinations_output_path(), &output::all_combinations_lines(&analyses))?;
        info!("write best haplotype combination output file");
        let genotypes = output::extract_genotypes(&analyses, drugs.as_ref(), functions.as_ref())?;
        output::write_lines(&config.best_haplotype_combinations_output_path(), &output::genotype_lines(&genotypes))?;
        info!("write qc status output file");
        output::write_lines(&config.qc_status_output_path(), &output::qc_lines(&analyses))?;
        Ok(())
    })();
    r
}
