//! ISOFOX in Rust: a port of hmftools ISOFOX 2.1 (see ../isofox).
//!
//! Usage: `isofox-rs <the Java options>`: reads an RNA BAM and writes the same files as the Java Isofox
//! (gene and transcript expression, summary statistics, gene collections, novel splice junctions, retained
//! introns, fusions). `-java_compat` reproduces the Java faults the port corrects (REPORT.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod adjusts;
mod allocator;
mod bam;
mod chimeric;
mod common;
mod compat;
mod config;
mod expression;
mod fragment;
mod fusion;
mod fusion_frag;
mod gene;
mod neo;
mod novel;
mod read;
mod read_trans;
mod refgen;
mod run;
mod task;
#[cfg(test)]
mod tests;

use hmf_common::config_builder::{ConfigBuilder, Parsed};
use hmf_common::error;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    // the bioconda wrapper passes JVM options; they are ignored
    let argv: Vec<String> = argv
        .into_iter()
        .filter(|a| !(a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))))
        .collect();
    let mut cb = ConfigBuilder::new(config::APP_NAME, config::VERSION);
    config::register(&mut cb);
    match cb.parse(&argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    }
    match run::run(&cb) {
        Ok(code) => {
            if code != 0 {
                hmf_common::info!("Isofox RNA analysis failed");
            }
            std::process::exit(code)
        }
        Err(e) => {
            error!("Isofox RNA analysis failed: {e:#}");
            std::process::exit(1);
        }
    }
}
