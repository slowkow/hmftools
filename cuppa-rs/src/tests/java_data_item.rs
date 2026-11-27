//! DataItemTest

use std::collections::HashMap;

use crate::prep::items::{DataItem, Index, ItemType, Source};

fn idx(s: Source, t: ItemType, k: &str) -> Index {
    Index { source: s, item_type: t, key: k.to_string() }
}

#[test]
fn can_sort_index() {
    let mut indexes = vec![
        idx(Source::Dna, ItemType::Signature, "SIG_2_13"),
        idx(Source::Dna, ItemType::Signature, "SIG_1"),
        idx(Source::Dna, ItemType::GenPos, "X_0"),
        idx(Source::Dna, ItemType::GenPos, "10_0"),
        idx(Source::Dna, ItemType::GenPos, "2_0"),
        idx(Source::Rna, ItemType::AltSj, "X;0;1"),
        idx(Source::Rna, ItemType::AltSj, "10;0;10"),
        idx(Source::Rna, ItemType::AltSj, "10;0;1"),
        idx(Source::Rna, ItemType::AltSj, "2;0;1"),
    ];
    indexes.sort_by(|a, b| a.compare(b));
    let keys: Vec<&str> = indexes.iter().map(|i| i.key.as_str()).collect();
    assert_eq!(keys, ["SIG_1", "SIG_2_13", "2_0", "10_0", "X_0", "2;0;1", "10;0;1", "10;0;10", "X;0;1"]);
}

#[test]
fn can_create_data_item_matrix_with_unique_keys() {
    let item = |k: &str, v: &str| DataItem { index: idx(Source::Dna, ItemType::Driver, k), value: v.to_string() };
    let per_sample = [vec![item("BRAF.mut", "1"), item("TP53.mut", "0.5"), item("TERT.mut", "1")], vec![item("BRAF.mut", "1")]];
    let mut m: HashMap<Index, Vec<Option<String>>> = HashMap::new();
    for (s, items) in per_sample.iter().enumerate() {
        for d in items {
            m.entry(d.index.clone()).or_insert_with(|| vec![None; 2])[s] = Some(d.value.clone());
        }
    }
    assert_eq!(m.len(), 3);
    assert_eq!(m[&idx(Source::Dna, ItemType::Driver, "BRAF.mut")], vec![Some("1".to_string()), Some("1".to_string())]);
    assert_eq!(m[&idx(Source::Dna, ItemType::Driver, "TP53.mut")], vec![Some("0.5".to_string()), None]);
}
