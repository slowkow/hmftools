//! CuppaDataPrepTest (the items kept in memory by Java are read back from the files written here)

use super::{temp_dir, test_config, testdata, with_sample_data_dir};
use crate::prep::config::CategoryType;
use crate::prep::items::{read_item_list, read_matrix, DataItem, Index, ItemType, Source};
use crate::prep::{output_path, run};

fn item(s: Source, t: ItemType, k: &str, v: &str) -> DataItem {
    DataItem { index: Index { source: s, item_type: t, key: k.to_string() }, value: v.to_string() }
}

#[test]
fn can_run_single_sample_prep() {
    let out = temp_dir("single");
    let mut c = with_sample_data_dir(test_config(&["TUMOR_SAMPLE"], Some(&["TUMOR_SAMPLE_RNA"])), &format!("{}*", testdata("java/pipeline_output/")));
    c.categories = CategoryType::ALL.to_vec();
    c.output_dir = out.clone();
    assert!(c.is_single_sample());
    run(&c).unwrap();
    let path = output_path(&c, None);
    assert!(std::path::Path::new(&path).exists());
    let items = read_item_list(&path).unwrap();
    assert_eq!(items.len(), 6229);
    assert_eq!(items[45], item(Source::Dna, ItemType::Snv96, "C>T_TCC", "2"));
    assert_eq!(items[items.len() - 1], item(Source::Rna, ItemType::AltSj, "20;46474688;46476232", "1"));
    let _ = std::fs::remove_dir_all(&out);
}

#[test]
fn can_run_multi_sample_prep() {
    let out = temp_dir("multi");
    let ids = ["MINIMAL_SAMPLE", "MINIMAL_SAMPLE", "MINIMAL_SAMPLE", "MINIMAL_SAMPLE_NO_RNA", "TUMOR_SAMPLE"];
    let rna = ["MINIMAL_SAMPLE", "MINIMAL_SAMPLE", "MINIMAL_SAMPLE", "NO_RNA", "TUMOR_SAMPLE_RNA"];
    let mut c = with_sample_data_dir(test_config(&ids, Some(&rna)), &format!("{}*", testdata("java/pipeline_output/")));
    c.categories = CategoryType::ALL.to_vec();
    c.output_dir = out.clone();
    c.threads = 2;
    run(&c).unwrap();
    for cat in CategoryType::ALL {
        let path = output_path(&c, Some(cat));
        assert!(std::path::Path::new(&path).exists(), "{path}");
        let m = read_matrix(&path).unwrap();
        let column = |s: usize| -> Vec<Option<String>> { m.indexes.iter().map(|i| m.values[i][s].clone()).collect() };
        // the duplicated sample has the same values in every column (no thread-unsafe state)
        assert_eq!(column(0), column(1), "{}", cat.name());
        assert_eq!(column(0), column(2), "{}", cat.name());
        assert_ne!(column(0), column(4), "{}", cat.name());
        if cat.is_rna() {
            assert!(column(3).iter().all(|v| v.is_none()), "NO_RNA sample has no RNA values");
        }
    }
    let _ = std::fs::remove_dir_all(&out);
}

#[test]
fn can_read_single_sample_output_file() {
    let items = read_item_list(&testdata("java/prep_output/TUMOR_SAMPLE.cuppa_data.modified.tsv")).unwrap();
    let expected = [
        item(Source::Dna, ItemType::Snv96, "C>T_ACG", "1"),
        item(Source::Dna, ItemType::TumorMutationalBurden, "snv_count", "8"),
        item(Source::Dna, ItemType::Signature, "SIG_7_UV", "6.4"),
        item(Source::Dna, ItemType::GenPos, "1_500000", "1"),
        item(Source::Dna, ItemType::SvCount, "LINE", "3"),
        item(Source::Dna, ItemType::SampleTrait, "is_male", "1"),
        item(Source::Dna, ItemType::Driver, "BRAF.mut", "1"),
        item(Source::Dna, ItemType::Fusion, "TMPRSS2_ERG", "1"),
        item(Source::Dna, ItemType::Virus, "MERKEL", "1"),
        item(Source::Rna, ItemType::Expression, "BRAF", "3.43E+00"),
        item(Source::Rna, ItemType::AltSj, "7;140426316;140439612", "2"),
    ];
    for (i, e) in expected.iter().enumerate() {
        assert_eq!(items[i], *e);
    }
}

#[test]
fn can_read_multi_sample_output_file() {
    let m = read_matrix(&testdata("java/prep_output/cuppa_data.cohort.gene_exp.modified.tsv")).unwrap();
    let braf = Index { source: Source::Rna, item_type: ItemType::Expression, key: "BRAF".into() };
    let tp53 = Index { source: Source::Rna, item_type: ItemType::Expression, key: "TP53".into() };
    assert_eq!(m.indexes[0], braf);
    assert_eq!(m.values[&braf], vec![None, Some("3.43E+00".to_string())]);
    assert_eq!(m.indexes[1], tp53);
    assert_eq!(m.values[&tp53], vec![None, Some("3.87E+00".to_string())]);
}
