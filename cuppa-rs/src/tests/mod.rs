//! Tests: the Java unit tests of CUPPA's prep, ported one module per class (`java_*`), pycuppa's tests of
//! the prediction components (`pycuppa`), and tests of the port (`port`: numpy's arithmetic against a
//! fixture written by numpy, the model resource, HG008's prediction, the findings in both modes).

mod java_category_prep;
mod java_cuppa_data_prep;
mod java_data_item;
mod java_prediction_runner;
mod java_prep_config;
mod java_viral_insertion_type;
mod port;
mod pycuppa;

use crate::prep::config::{CategoryType, PrepConfig};

pub fn testdata(rel: &str) -> String {
    format!("{}/testdata/{rel}", env!("CARGO_MANIFEST_DIR"))
}

/// TestPrepConfigBuilder's defaults: V37, DNA categories, the test resources' alt-SJ sites and liftover
/// directory, write by category, 1 thread
pub fn test_config(samples: &[&str], rna: Option<&[&str]>) -> PrepConfig {
    PrepConfig {
        sample_ids: samples.iter().map(|s| s.to_string()).collect(),
        rna_sample_ids: rna.unwrap_or(samples).iter().map(|s| s.to_string()).collect(),
        sample_data_dir: String::new(),
        linx_dir: String::new(),
        purple_dir: String::new(),
        viridian_dir: String::new(),
        virus_dir: None,
        isofox_dir: String::new(),
        somatic_variants_dir: testdata("java/liftover_output/"),
        categories: CategoryType::ALL[..4].to_vec(),
        v38: false,
        alt_sj_sites: Some(testdata("java/alt_sj.selected_loci.minimal.tsv")),
        output_dir: String::new(),
        output_id: None,
        write_by_category: true,
        threads: 1,
    }
}

/// TestPrepConfigBuilder.sampleDataDir: every tool directory
pub fn with_sample_data_dir(mut c: PrepConfig, dir: &str) -> PrepConfig {
    c.sample_data_dir = dir.to_string();
    c.purple_dir = dir.to_string();
    c.linx_dir = dir.to_string();
    c.viridian_dir = dir.to_string();
    c.isofox_dir = dir.to_string();
    c.somatic_variants_dir = dir.to_string();
    c
}

/// a fresh directory under the system's temporary directory
pub fn temp_dir(name: &str) -> String {
    let d = std::env::temp_dir().join(format!("cuppa-rs-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    format!("{}/", d.display())
}
