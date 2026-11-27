//! CategoryPrepTest

use std::collections::HashMap;

use super::{test_config, testdata, with_sample_data_dir};
use crate::prep::config::{CategoryType, PrepConfig};
use crate::prep::items::{DataItem, Index, ItemType, Source};
use crate::prep::{drivers, rna, somatic, sv, traits};

const TUMOR: &str = "TUMOR_SAMPLE";
const TUMOR_RNA: &str = "TUMOR_SAMPLE_RNA";

fn config() -> PrepConfig {
    let mut c = with_sample_data_dir(test_config(&[TUMOR], Some(&[TUMOR_RNA])), &format!("{}*", testdata("java/pipeline_output/")));
    c.categories = CategoryType::ALL.to_vec();
    c.output_dir = "/tmp/".to_string();
    c
}

fn map(items: &[DataItem]) -> HashMap<String, String> {
    items.iter().map(|d| (d.index.key.clone(), d.value.clone())).collect()
}

fn item(s: Source, t: ItemType, k: &str, v: &str) -> DataItem {
    DataItem { index: Index { source: s, item_type: t, key: k.to_string() }, value: v.to_string() }
}

#[test]
fn can_extract_snv_features() {
    let items = somatic::extract(&config(), TUMOR).unwrap();
    for t in [ItemType::Snv96, ItemType::TumorMutationalBurden, ItemType::GenPos, ItemType::Signature] {
        assert!(items.iter().any(|d| d.index.item_type == t));
    }
    let m = map(&items);
    assert_eq!(m["C>T_TCC"], "2");
    assert_eq!(m["snv_count"], "8");
    assert_eq!(m["SIG_7_UV"], "6.4");
    assert_eq!(m["1_500000"], "1");
    assert_eq!(m["1_1000000"], "1");
    assert_eq!(m["7_55000000"], "1");
}

#[test]
fn can_extract_snv_features_from_liftover_generic_variants_files() {
    let mut c = test_config(&[TUMOR], None);
    c.categories = vec![CategoryType::Snv];
    c.v38 = true;
    let m = map(&somatic::extract(&c, TUMOR).unwrap());
    assert_eq!(m["C>T_TCC"], "2");
    assert_eq!(m["snv_count"], "8");
    assert_eq!(m["SIG_7_UV"], "6.4");
    assert_eq!(m["chr1_500000"], "1");
    assert_eq!(m["chr1_1000000"], "1");
    assert_eq!(m["chr7_55000000"], "1");
}

#[test]
fn can_extract_sv_features() {
    let items = sv::extract(&config(), TUMOR).unwrap();
    let expected = [("LINE", "3"), ("SIMPLE_DEL_20KB_1MB", "25"), ("SIMPLE_DUP_32B_200B", "4"), ("SIMPLE_DUP_100KB_5MB", "2"), ("MAX_COMPLEX_SIZE", "8"), ("TELOMERIC_SGL", "0")];
    assert_eq!(items.len(), 6);
    for (d, (k, v)) in items.iter().zip(expected) {
        assert_eq!(*d, item(Source::Dna, ItemType::SvCount, k, v));
    }
}

#[test]
fn can_extract_trait_features() {
    let items = traits::extract(&config(), TUMOR).unwrap();
    assert_eq!(
        items,
        vec![
            item(Source::Dna, ItemType::SampleTrait, "is_male", "1"),
            item(Source::Dna, ItemType::TumorMutationalBurden, "indels_per_mb", "0.1207"),
            item(Source::Dna, ItemType::SampleTrait, "whole_genome_duplication", "1"),
        ]
    );
}

#[test]
fn can_extract_event_features() {
    let items = drivers::extract(&config(), TUMOR).unwrap();
    assert_eq!(
        items,
        vec![
            item(Source::Dna, ItemType::Driver, "BRAF.mut", "1.0000"),
            item(Source::Dna, ItemType::Driver, "MYC.amp", "1.0000"),
            item(Source::Dna, ItemType::Driver, "CDKN2A.mut", "1.0000"),
            item(Source::Dna, ItemType::Driver, "SFTPB.indel", "1.0000"),
            item(Source::Dna, ItemType::Fusion, "TMPRSS2_ERG", "1.0000"),
            item(Source::Dna, ItemType::Fusion, "BRAF_PROM3", "1.0000"),
            item(Source::Dna, ItemType::Fusion, "FGFR2_PROM5", "1.0000"),
            item(Source::Dna, ItemType::Virus, "MERKEL", "1.0000"),
        ]
    );
}

#[test]
fn can_extract_gene_exp_features() {
    let items = rna::extract_gene_expression(&config(), TUMOR_RNA).unwrap();
    assert_eq!(
        items,
        vec![item(Source::Rna, ItemType::Expression, "BRAF", "3.434e+00"), item(Source::Rna, ItemType::Expression, "TP53", "3.870e+00")]
    );
}

#[test]
fn can_extract_alt_sj_features() {
    let items = rna::extract_alt_sj(&config(), TUMOR_RNA).unwrap().unwrap();
    assert_eq!(
        items,
        vec![
            item(Source::Rna, ItemType::AltSj, "7;140426316;140439612", "2"),
            item(Source::Rna, ItemType::AltSj, "10;89623836;89623905", "1"),
            item(Source::Rna, ItemType::AltSj, "13;32901736;32901958", "1"),
            item(Source::Rna, ItemType::AltSj, "20;46474688;46476232", "1"),
        ]
    );
}
