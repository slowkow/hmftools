//! The Java unit tests of CHORD, ported (one module per Java test class), and tests of the port.

mod chord_application;
mod indel_prep;
mod indel_variant;
mod port;
mod snv_prep;
mod structural_variant;
mod sv_prep;

use crate::config::ChordConfig;
use crate::prep::{LazyRefGenome, MutContextCount};

/// ChordTestDataPaths
pub const MINIMAL_SAMPLE: &str = "MINIMAL_SAMPLE";
pub const EMPTY_SAMPLE: &str = "EMPTY_SAMPLE";
pub const NON_STANDARD_NUC_GENOME_SAMPLE: &str = "NON_STANDARD_NUC_GENOME_SAMPLE";

pub fn testdata(rel: &str) -> String {
    format!("{}/testdata/{rel}", env!("CARGO_MANIFEST_DIR"))
}

pub fn input_vcf_dir() -> String {
    testdata("vcf/")
}

pub fn dummy_genome_fasta() -> String {
    testdata("fasta/dummy_genome.fasta")
}

pub fn non_standard_nuc_genome_fasta() -> String {
    testdata("fasta/non_standard_nuc_genome.fasta")
}

/// A fresh directory for one test's outputs
pub fn tmp_output_dir(name: &str) -> String {
    let dir = std::env::temp_dir().join(format!("chord_rs_test_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    format!("{}/", dir.display())
}

/// ChordConfig.Builder's defaults
pub fn config() -> ChordConfig {
    ChordConfig { threads: 1, output_id: Some(String::new()), log_level: "INFO".into(), ..Default::default() }
}

pub fn mc(name: &str, count: i32) -> MutContextCount {
    MutContextCount::new(name, count)
}

/// SnvPrep.countMutationContexts
pub fn snv_counts(config: &ChordConfig, sample: &str) -> anyhow::Result<Vec<MutContextCount>> {
    let Some(ref_path) = &config.ref_genome_file else {
        anyhow::bail!("java.lang.IllegalStateException: Ref genome must be provided to run SnvPrep");
    };
    let variants = crate::vcf_file::load_variants(&config.snv_indel_vcf_file(sample)?, config.include_non_pass, "")?;
    let snvs = crate::snv::load_snvs(&variants, "");
    crate::snv::count_mutation_contexts(&snvs, &LazyRefGenome::new(ref_path), false, "", |_| Ok(()))
}

/// IndelPrep.countMutationContexts
pub fn indel_counts(config: &ChordConfig, sample: &str) -> anyhow::Result<Vec<MutContextCount>> {
    let Some(ref_path) = &config.ref_genome_file else {
        anyhow::bail!("java.lang.IllegalStateException: Ref genome must be provided to run IndelPrep");
    };
    let variants = crate::vcf_file::load_variants(&config.snv_indel_vcf_file(sample)?, config.include_non_pass, "")?;
    let indels = crate::indel::load_indels(&variants, "");
    crate::indel::count_mutation_contexts(&indels, &LazyRefGenome::new(ref_path), false, "", |_| Ok(()))
}

/// SvPrep.countMutationContexts
pub fn sv_counts(config: &ChordConfig, sample: &str) -> anyhow::Result<Vec<MutContextCount>> {
    let variants = crate::vcf_file::load_variants(&config.sv_vcf_file(sample)?, config.include_non_pass, "")?;
    let svs = crate::sv::load_svs(&variants, "")?;
    crate::sv::count_mutation_contexts(&svs, "")
}
