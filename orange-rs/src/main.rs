//! orange-rs: hmftools ORANGE 5.1 in Rust (the JSON report). See README.md and NOTES.md.

mod algo;
mod compat;
mod config;
mod doid;
mod isofox;
mod json;
mod linx;
mod model;
mod others;
mod plots;
mod purple;

#[cfg(test)]
mod tests;

use anyhow::Result;
use hmf_common::config_builder::{ConfigBuilder, Parsed};

/// The JVM options of the bioconda wrapper (`orange -Xmx16g ...`) are accepted and ignored
fn is_jvm_option(a: &str) -> bool {
    a.starts_with("-Xmx") || a.starts_with("-Xms") || a.starts_with("-Xss") || a.starts_with("-XX:") || (a.starts_with("-D") && a.contains('='))
}

fn run() -> Result<i32> {
    let argv: Vec<String> = std::env::args().skip(1).filter(|a| !is_jvm_option(a)).collect();
    let argv = hmf_common::compat::take_options(&argv, "ORANGE", compat::valid_id)?;
    let mut cb = ConfigBuilder::new(config::APP_NAME, config::VERSION);
    config::register(&mut cb);
    match cb.parse(&argv) {
        Ok(Parsed::Exit) => return Ok(0),
        Ok(Parsed::Run) => {}
        Err(_) => return Ok(1),
    }
    let cfg = config::OrangeConfig::from_builder(&cb)?;
    let record = algo::run(&cfg)?;
    // ReportWriter.write: the PDF is not written by the port (out of scope; see NOTES.md). OR10: Java writes the
    // PDF first, and its front page reads the CHORD record of every run with a reference sample, so a tumor/normal
    // run without CHORD (targeted, or no -chord_dir) fails there and writes no JSON.
    if compat::java("OR10") && record.reference_id.is_some() && record.chord.is_none() {
        anyhow::bail!(
            "java.lang.NullPointerException: Cannot invoke \"com.hartwig.hmftools.datamodel.chord.ChordRecord.hrStatus()\" because \"chord\" is null (FrontPageTables.hrDeficiencyString)"
        );
    }
    // OR11: the PDF's fusion table shows RNA support once any fusion has it, and reads it for every fusion, so
    // a run with RNA where one reported fusion has RNA support and another has none fails there
    let fusions = &record.linx.fusions;
    if compat::java("OR11") && fusions.iter().any(|f| f.rna_support.is_some()) && fusions.iter().any(|f| f.rna_support.is_none()) {
        anyhow::bail!(
            "java.lang.NullPointerException: Cannot invoke \"AllelicDepth.alleleReadCount()\" because the return value of \"LinxFusion.rnaSupport()\" is null (DnaFusionTable.build)"
        );
    }
    let mut name = format!("{}{}.orange.", cfg.output_dir, record.sample_id);
    if let Some(id) = &cfg.output_id {
        name.push_str(id);
        name.push('.');
    }
    name.push_str("json");
    hmf_common::info!("writing JSON report to {name} ");
    std::fs::write(&name, model::to_json(&record))?;
    hmf_common::info!("Orange complete");
    Ok(0)
}

fn main() {
    match run() {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            hmf_common::error!("{e:#}");
            std::process::exit(1);
        }
    }
}
