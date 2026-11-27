//! SIGS in Rust: a port of hmftools SIGS 1.3 (see ../sigs).
//!
//! Usage: `sigs-rs <the Java options>`: counts a sample's passing SNVs of a somatic VCF in the 96
//! trinucleotide buckets, fits the counts to the signature definitions by non-negative least squares, and
//! writes `<sample>.sig.snv_counts.csv` and `<sample>.sig.allocation.tsv` as the Java SignatureFitter does.
//! `-java_compat` reproduces the Java faults the port corrects (REPORT.md, SI1-SI7).

mod compat;
mod config;
mod fitter;
mod lsq;
mod matrix;
mod position_freq;
mod snv_loader;
#[cfg(test)]
mod tests;

use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::error;

pub const VERSION: &str = "1.3";
pub const APP_NAME: &str = "Sigs";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(sigs_main(&argv));
}

/// The JVM options of the bioconda wrapper (`sigs -Xmx9663676416 ...`), accepted and ignored
pub fn strip_jvm_options(argv: &[String]) -> Vec<String> {
    argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect()
}

/// SignatureFitter.main: the exit status
pub fn sigs_main(argv: &[String]) -> i32 {
    let argv = strip_jvm_options(argv);
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let mut cb = ConfigBuilder::new(APP_NAME, VERSION);
    config::register(&mut cb);
    match cb.parse(&argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    }
    let result = config::SigsConfig::from(&cb).and_then(|c| fitter::run(&c, &config::VariantFilters::from(&cb)));
    match result {
        Ok(()) => 0,
        Err(e) => {
            error!("{e:#}");
            1
        }
    }
}
