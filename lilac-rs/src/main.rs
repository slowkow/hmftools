//! LILAC in Rust: a port of hmftools LILAC 2.1 (HLA typing; see ../../hmftools/lilac).
//!
//! Usage: `lilac-rs <the Java options>`: reads the same files and writes the same outputs as the Java
//! LilacApplication: the HLA class I and class II gene groups, tumor/normal, tumor-only and germline-only,
//! with the PURPLE copy numbers and somatic variants and an RNA BAM. `-java_compat` reproduces the Java
//! faults the port corrects (REPORT.md).

// The Java data model is ported whole: some fields and helpers are only read by the tests or the debug output.
#![allow(dead_code)]

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod app;
mod bam;
mod compat;
mod complexes;
mod config;
mod constants;
mod coverage;
mod evidence;
mod fragment;
mod hla;
mod hlay;
mod indel;
#[cfg(test)]
mod java_tests;
mod math;
mod output;
mod qc;
mod reader;
mod refdata;
mod seq;
mod tok;
mod variant;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};
// the paths the shared code's macros use (crate::log)
#[allow(unused_imports)]
use hmf_common::{fmt, log};

pub const VERSION: &str = "2.1";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(lilac_main(&argv));
}

/// LilacApplication.main: the exit status
pub fn lilac_main(argv: &[String]) -> i32 {
    // the JVM options of the bioconda wrapper (lilac -Xmx16g ...) are accepted and ignored
    let argv: Vec<String> = argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect();
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let cb = match config::ConfigBuilder::parse(&argv) {
        Ok(config::Parsed::Run(cb)) => cb,
        Ok(config::Parsed::Exit) => return 0,
        Err(_) => return 1,
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("Lilac version {VERSION}");
    let config = match config::LilacConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    if std::fs::create_dir_all(&config.output_dir).is_err() && !std::path::Path::new(&config.output_dir).exists() {
        error!("failed to create output directory: {}", config.output_dir);
        return 1;
    }
    info!("writing output to {}", config.output_dir);
    if let Err(e) = app::run(&config) {
        error!("{e:#}");
        return 1;
    }
    0
}
