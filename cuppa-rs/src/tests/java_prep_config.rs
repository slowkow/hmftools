//! PrepConfigTest

use hmf_common::config_builder::{ConfigBuilder, Parsed};

use super::testdata;
use crate::prep::config::{self, wildcard, CategoryType, PrepConfig};

const SELECTED: &str = "SAMPLE_1";

fn from_args(args: &[String]) -> PrepConfig {
    let mut cb = ConfigBuilder::new(config::APP_NAME, config::VERSION);
    config::register(&mut cb, false);
    assert!(matches!(cb.parse(args), Ok(Parsed::Run)));
    PrepConfig::from(&cb).unwrap()
}

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

#[test]
fn can_parse_command_line_args_multi_sample() {
    let data = testdata("java/pipeline_output/");
    let sites = testdata("java/alt_sj.selected_loci.minimal.tsv");
    let c = from_args(&args(&[
        "-sample_id_file",
        &testdata("java/pipeline_output/sample_ids.tsv"),
        "-categories",
        "ALL",
        "-ref_genome_version",
        "V37",
        "-sample_data_dir",
        &format!("{data}*"),
        "-purple_dir",
        &format!("{data}*/purple/"),
        "-linx_dir",
        &format!("{data}*/linx/"),
        "-viridian_dir",
        &format!("{data}*/viridian/"),
        "-isofox_dir",
        &format!("{data}*/isofox/"),
        "-ref_alt_sj_sites",
        &sites,
        "-output_dir",
        "/tmp/",
        "-threads",
        "1",
        "-write_by_category",
    ]));
    assert_eq!(c.sample_ids.len(), 5);
    assert_eq!(c.rna_sample_ids.len(), 5);
    assert!(!c.v38);
    assert!(c.write_by_category);
    assert_eq!(c.threads, 1);
    assert_eq!(c.categories, CategoryType::ALL.to_vec());
    assert_eq!(c.alt_sj_sites.as_deref(), Some(sites.as_str()));
    assert_eq!(c.purple_dir(SELECTED), format!("{data}{SELECTED}/purple/"));
    assert_eq!(c.linx_dir(SELECTED), format!("{data}{SELECTED}/linx/"));
    assert_eq!(wildcard(&c.viridian_dir, SELECTED), format!("{data}{SELECTED}/viridian/"));
    assert_eq!(wildcard(&c.isofox_dir, SELECTED), format!("{data}{SELECTED}/isofox/"));
    // SampleIdsLoader: a blank RnaSampleId is the SampleId, NO_RNA kept
    assert_eq!(c.rna_sample_ids[0], "MINIMAL_SAMPLE");
    assert_eq!(c.rna_sample_ids[3], "NO_RNA");
    assert_eq!(c.rna_sample_ids[4], "TUMOR_SAMPLE_RNA");
}

#[test]
fn sample_dir_used_when_tool_dir_not_specified() {
    let data = testdata("java/pipeline_output/");
    let c = from_args(&args(&["-sample_id_file", &testdata("java/pipeline_output/sample_ids.tsv"), "-sample_data_dir", &format!("{data}*"), "-output_dir", "/tmp/"]));
    assert_eq!(c.purple_dir(SELECTED), format!("{data}{SELECTED}"));
    // the defaults: DNA categories, V37, 1 thread
    assert_eq!(c.categories, CategoryType::ALL[..4].to_vec());
    assert!(!c.v38);
    assert_eq!(c.threads, 1);
}

#[test]
fn can_build_from_test_prep_config_builder() {
    let data = testdata("java/pipeline_output/");
    let mut c = super::with_sample_data_dir(super::test_config(&["sample_1", "sample_2"], None), &format!("{data}*"));
    c.linx_dir = format!("{data}*/linx/");
    assert_eq!(c.purple_dir(SELECTED), format!("{data}{SELECTED}"));
    assert_eq!(c.linx_dir(SELECTED), format!("{data}{SELECTED}/linx/"));
}
