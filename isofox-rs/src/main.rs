//! ISOFOX in Rust: a port of hmftools ISOFOX 2.1 (see ../isofox), and of ISOFOX 2.0.1, the version for
//! STAR-aligned RNA (`-isofox_version 2.0.1`, see version.rs).
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
mod version;
#[cfg(test)]
mod tests201;
mod allocator201;

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
    // the ISOFOX version to reproduce: the option or ISOFOX_VERSION, else the BAM's aligner (STAR: 2.0.1)
    let (argv, chosen) = match version::take_option(&argv) {
        Ok(x) => x,
        Err(e) => {
            error!("{e}");
            std::process::exit(1);
        }
    };
    let selected = match chosen {
        Some(v) => v,
        None => match version::bam_file_arg(&argv).and_then(bam::header_text) {
            Some(h) => version::detect_from_header_text(&h).unwrap_or(version::IsofoxVersion::V2_1),
            None => version::IsofoxVersion::V2_1,
        },
    };
    version::set(selected);
    let mut cb = ConfigBuilder::new(config::APP_NAME, selected.java_version());
    config::register(&mut cb);
    match cb.parse(&argv) {
        Ok(Parsed::Run) => {}
        Ok(Parsed::Exit) => return,
        Err(_) => std::process::exit(1),
    }
    if version::v201() {
        let why = if chosen.is_some() { "as requested" } else { "BAM aligned by STAR" };
        hmf_common::info!("reproducing ISOFOX 2.0.1 ({})", why);
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
