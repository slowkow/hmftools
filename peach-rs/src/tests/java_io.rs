//! DrugInfoLoaderTest, HaplotypeFunctionLoaderTest, PanelLoaderTest, HaplotypeEventLoaderTest, DrugInfoStoreTest,
//! HaplotypeCombinationTest, PeachConfigTest, and the output tests (EventsFileTest, EventsPerGeneFileTest,
//! AllHaplotypeCombinationsFileTest, QcStatusFileTest, PeachGenotypeExtractorTest)

use std::collections::{BTreeMap, HashSet};

use super::testdata;
use crate::algo::{EventCounts, HaplotypeAnalysis, HaplotypeCombination, PeachQcStatus};
use crate::loaders::*;
use crate::output::*;
use crate::panel::*;

fn chr(c: &str) -> Chromosome {
    Chromosome::from_string(c).unwrap()
}

#[test]
fn drug_info_loader() {
    let d = load_drug_infos(&testdata("drugs.tsv")).unwrap();
    assert_eq!(d.len(), 4);
    let info = |drug: &str, gene: &str, url: &str| DrugInfo { drug_name: drug.into(), gene_name: gene.into(), prescription_info_url: url.into() };
    assert_eq!(d[0], info("5-Fluorouracil", "DPYD", "https://www.pharmgkb.org/guidelineAnnotation/PA166104939"));
    assert_eq!(d[1], info("Capecitabine", "DPYD", "https://www.pharmgkb.org/guidelineAnnotation/PA166104963"));
    assert_eq!(d[2], info("Tegafur", "DPYD", "https://www.pharmgkb.org/guidelineAnnotation/PA166104944"));
    assert_eq!(d[3], info("Irinotecan", "UGT1A1", "https://www.pharmgkb.org/guidelineAnnotation/PA166104951"));
}

#[test]
fn haplotype_function_loader() {
    let f = load_functions(&testdata("function.tsv")).unwrap();
    assert_eq!(f.len(), 11);
    assert_eq!(f[0], HaplotypeFunction { gene_name: "DPYD".into(), haplotype_name: "*1".into(), function: "Normal Function".into() });
    assert_eq!(f[10], HaplotypeFunction { gene_name: "UGT1A1".into(), haplotype_name: "*6".into(), function: "Reduced Function".into() });
}

#[test]
fn panel_loader() {
    let p = load_haplotype_panel(&testdata("haplotypes.complicated.37.tsv")).unwrap();
    assert_eq!(p.genes.len(), 2);
    let d = &p.genes["DPYD"].default_haplotype;
    assert_eq!((d.name.as_str(), d.is_wild_type, d.events_to_ignore.len()), ("*9A", false, 0));
    let u = &p.genes["UGT1A1"].default_haplotype;
    assert_eq!((u.name.as_str(), u.is_wild_type, u.events_to_ignore.len()), ("*1", true, 0));
    let ugt = &p.genes["UGT1A1"].non_default_haplotypes;
    assert_eq!(ugt.len(), 4);
    assert_eq!((ugt[0].name.as_str(), ugt[0].is_wild_type), ("*28", false));
    assert_eq!(ugt[0].events.iter().map(|e| e.id()).collect::<Vec<_>>(), vec!["VAR_chr2_234668879_C_CAT"]);
    let dpyd = &p.genes["DPYD"].non_default_haplotypes;
    assert_eq!(dpyd.len(), 6);
    assert_eq!((dpyd[1].name.as_str(), dpyd[1].is_wild_type), ("*2A", false));
    assert_eq!(dpyd[1].events.iter().map(|e| e.id()).collect::<Vec<_>>(), vec!["VAR_chr1_97915614_C_T", "VAR_chr1_98348885_G_A"]);
}

fn expected_events() -> EventCounts {
    [
        ("VAR_chr1_98205966_GATGA_G", Some(1)),
        ("VAR_chr1_98205966_G_C", Some(1)),
        ("VAR_chr1_98348885_G_A", Some(2)),
        ("VAR_chr2_234668879_C_CAT", None),
        ("VAR_chr2_234668880_A_G", Some(0)),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_string(), v))
    .collect()
}

fn test_positions() -> BTreeMap<Chromosome, HashSet<i32>> {
    let mut m = BTreeMap::new();
    m.insert(chr("1"), [98205966, 98348885].into_iter().collect());
    m.insert(chr("2"), [234668879, 234668880].into_iter().collect());
    m.insert(chr("4"), [1733823].into_iter().collect());
    m
}

#[test]
fn haplotype_event_loader() {
    for java in [true, false] {
        crate::compat::with_java(java, || {
            assert_eq!(load_relevant_events(&testdata("variants.vcf"), "FAKER", &test_positions()).unwrap(), expected_events());
            assert_eq!(load_relevant_events(&testdata("variants.vcf.gz"), "FAKER", &test_positions()).unwrap(), expected_events());
            let panel = super::java_algo::default_test_panel();
            assert_eq!(load_relevant_events(&testdata("variants.vcf"), "FAKER", &panel.relevant_variant_positions()).unwrap(), expected_events());
        });
    }
}

#[test]
fn drug_info_store_missing_info() {
    let store = DrugInfoStore { drug_infos: vec![] };
    assert_eq!(store.prescription_info_url("FAKE_GENE", "FAKE_DRUG").unwrap(), None);
}

#[test]
fn haplotype_combination_equals() {
    let combo = HaplotypeCombination::of(&[("*1", 2)]);
    assert_eq!(combo, combo.clone());
    assert_eq!(HaplotypeCombination::of(&[("*1", 2)]), combo);
    assert_ne!(HaplotypeCombination::of(&[("*1", 3)]), combo);
    assert_ne!(HaplotypeCombination::of(&[("*1", 2), ("*7", 1)]), combo);
}

#[test]
fn peach_config() {
    for maximal in [false, true] {
        let mut cb = hmf_common::config_builder::ConfigBuilder::new("Peach", "2.0");
        crate::register(&mut cb);
        let mut args: Vec<String> = vec![
            "-vcf_file".into(),
            testdata("variants.vcf.gz"),
            "-haplotypes_file".into(),
            testdata("haplotypes.complicated.37.tsv"),
            "-sample_name".into(),
            "FAKENAME".into(),
            "-output_dir".into(),
            "/path/to/output".into(),
        ];
        if maximal {
            args.extend(["-drugs_file".into(), testdata("drugs.tsv"), "-function_file".into(), testdata("function.tsv")]);
        }
        assert!(matches!(cb.parse(&args), Ok(hmf_common::config_builder::Parsed::Run)));
        let c = crate::PeachConfig::from(&cb).unwrap();
        assert_eq!(c.vcf_file, testdata("variants.vcf.gz"));
        assert_eq!(c.sample_name, "FAKENAME");
        assert_eq!(c.output_dir, "/path/to/output/");
        assert_eq!(c.drugs_file.is_some(), maximal);
        assert_eq!(c.function_file.is_some(), maximal);
        assert_eq!(c.events_output_path(), "/path/to/output/FAKENAME.peach.events.tsv");
        assert_eq!(c.events_per_gene_output_path(), "/path/to/output/FAKENAME.peach.gene.events.tsv");
        assert_eq!(c.all_haplotype_combinations_output_path(), "/path/to/output/FAKENAME.peach.haplotypes.all.tsv");
        assert_eq!(c.best_haplotype_combinations_output_path(), "/path/to/output/FAKENAME.peach.haplotypes.best.tsv");
        assert_eq!(c.qc_status_output_path(), "/path/to/output/FAKENAME.peach.qc.tsv");
    }
}

fn analysis(events: &[(&str, Option<i32>)], combos: Vec<HaplotypeCombination>, default: &str, qc: PeachQcStatus, best: Option<HaplotypeCombination>) -> HaplotypeAnalysis {
    HaplotypeAnalysis {
        event_counts: events.iter().map(|(e, c)| (e.to_string(), *c)).collect(),
        haplotype_combinations: combos,
        default_haplotype_name: default.into(),
        wild_type_haplotype_name: "*1".into(),
        qc_status: qc,
        best_haplotype_combination: best,
    }
}

fn c(pairs: &[(&str, i32)]) -> HaplotypeCombination {
    HaplotypeCombination::of(pairs)
}

#[test]
fn events_file() {
    assert_eq!(events_lines(&EventCounts::new()), vec!["event\tcount"]);
    let e: EventCounts = [("EVENT_3", Some(1)), ("EVENT_1", Some(2)), ("EVENT_2", None)].into_iter().map(|(k, v)| (k.to_string(), v)).collect();
    assert_eq!(events_lines(&e), vec!["event\tcount", "EVENT_1\t2", "EVENT_2\tUNKNOWN", "EVENT_3\t1"]);
}

#[test]
fn events_per_gene_file() {
    assert_eq!(events_per_gene_lines(&BTreeMap::new()), vec!["gene\tevent\tcount"]);
    let mut m = BTreeMap::new();
    m.insert("FAKE4".to_string(), analysis(&[("EVENT_1", None)], vec![], "*9", PeachQcStatus::FailNoCombinationFound, None));
    m.insert("FAKE3".to_string(), analysis(&[("EVENT_1", Some(2))], vec![c(&[("*4", 1), ("*3", 1)])], "*9", PeachQcStatus::Pass, Some(c(&[("*4", 1), ("*3", 1)]))));
    m.insert(
        "FAKE2".to_string(),
        analysis(&[("EVENT_3", Some(1)), ("EVENT_1", Some(2))], vec![c(&[("2373C>T", 1), ("*9", 1)]), c(&[("*2", 1), ("*5", 1)])], "*9", PeachQcStatus::FailNoUniqueBestCombinationFound, None),
    );
    m.insert("FAKE1".to_string(), analysis(&[], vec![c(&[("*1", 2)])], "*1", PeachQcStatus::Pass, Some(c(&[("*1", 2)]))));
    assert_eq!(
        events_per_gene_lines(&m),
        vec!["gene\tevent\tcount", "FAKE2\tEVENT_1\t2", "FAKE2\tEVENT_3\t1", "FAKE3\tEVENT_1\t2", "FAKE4\tEVENT_1\tUNKNOWN"]
    );
}

#[test]
fn all_haplotype_combinations_file() {
    assert_eq!(all_combinations_lines(&BTreeMap::new()), vec!["gene\tcombination\tcount\tnonWildTypeCount"]);
    let mut m = BTreeMap::new();
    m.insert("FAKE3".to_string(), analysis(&[("EVENT_1", Some(2))], vec![c(&[("*9", 2)])], "*9", PeachQcStatus::Pass, Some(c(&[("*9", 2)]))));
    m.insert(
        "FAKE2".to_string(),
        analysis(
            &[("EVENT_1", Some(1)), ("EVENT_3", Some(2))],
            vec![c(&[("2373C>T", 1), ("*9", 1)]), c(&[("*2", 1), ("*1", 2)]), c(&[("*3", 1), ("*1", 1)])],
            "*9",
            PeachQcStatus::Pass,
            Some(c(&[("*3", 1), ("*1", 1)])),
        ),
    );
    m.insert("FAKE1".to_string(), analysis(&[], vec![c(&[("*1", 2)])], "*1", PeachQcStatus::Pass, Some(c(&[("*1", 2)]))));
    assert_eq!(
        all_combinations_lines(&m),
        vec![
            "gene\tcombination\tcount\tnonWildTypeCount",
            "FAKE1\t(*1,2)\t2\t0",
            "FAKE2\t(*1,1);(*3,1)\t2\t1",
            "FAKE2\t(*9,1);(2373C>T,1)\t2\t2",
            "FAKE2\t(*1,2);(*2,1)\t3\t1",
            "FAKE3\t(*9,2)\t2\t2",
        ]
    );
}

#[test]
fn qc_status_file() {
    assert_eq!(qc_lines(&BTreeMap::new()), vec!["gene\tstatus"]);
    let mut m = BTreeMap::new();
    m.insert("FAKE3".to_string(), analysis(&[("EVENT_1", Some(2))], vec![], "*9", PeachQcStatus::FailNoCombinationFound, None));
    m.insert("FAKE2".to_string(), analysis(&[], vec![], "*9", PeachQcStatus::FailNoUniqueBestCombinationFound, None));
    m.insert("FAKE1".to_string(), analysis(&[], vec![c(&[("*1", 2)])], "*1", PeachQcStatus::Pass, Some(c(&[("*1", 2)]))));
    m.insert("FAKE4".to_string(), analysis(&[], vec![], "*9", PeachQcStatus::WarnTooManyAllelesFound, Some(c(&[("*3", 2), ("*2", 1)]))));
    m.insert("FAKE5".to_string(), analysis(&[("EVENT_6", None)], vec![], "*9", PeachQcStatus::FailEventWithUnknownCount, None));
    assert_eq!(
        qc_lines(&m),
        vec![
            "gene\tstatus",
            "FAKE1\tPASS",
            "FAKE2\tFAIL_NO_UNIQUE_BEST_COMBINATION_FOUND",
            "FAKE3\tFAIL_NO_COMBINATION_FOUND",
            "FAKE4\tWARN_TOO_MANY_ALLELES_FOUND",
            "FAKE5\tFAIL_EVENT_WITH_UNKNOWN_COUNT",
        ]
    );
}

fn fake_analyses() -> BTreeMap<String, HaplotypeAnalysis> {
    let mut m = BTreeMap::new();
    m.insert("FAKE1".to_string(), analysis(&[], vec![c(&[("*1", 2)])], "*1", PeachQcStatus::Pass, Some(c(&[("*1", 2)]))));
    m.insert(
        "FAKE2".to_string(),
        analysis(&[("EVENT_1", Some(1)), ("EVENT_3", Some(2))], vec![c(&[("2373C>T", 1), ("*9", 1)]), c(&[("*2", 1), ("*5", 1)])], "*9", PeachQcStatus::FailNoUniqueBestCombinationFound, None),
    );
    m.insert("FAKE3".to_string(), analysis(&[("EVENT_1", Some(2))], vec![c(&[("*4", 1), ("*3", 1)])], "*9", PeachQcStatus::Pass, Some(c(&[("*4", 1), ("*3", 1)]))));
    m
}

fn g(gene: &str, allele: &str, count: i32, function: &str, drugs: &str, urls: &str) -> PeachGenotype {
    PeachGenotype {
        gene: gene.into(),
        allele: allele.into(),
        allele_count: count,
        function: function.into(),
        linked_drugs: drugs.into(),
        url_prescription_info: urls.into(),
    }
}

fn fake_genotypes() -> Vec<PeachGenotype> {
    vec![
        g("FAKE1", "*1", 2, "", "", ""),
        g("FAKE2", "Unresolved Haplotype", 2, "", "", ""),
        g("FAKE3", "*3", 1, "", "", ""),
        g("FAKE3", "*4", 1, "", "", ""),
    ]
}

#[test]
fn peach_genotype_extractor() {
    let drugs = DrugInfoStore { drug_infos: load_drug_infos(&testdata("drugs.tsv")).unwrap() };
    let functions = HaplotypeFunctionStore { functions: load_functions(&testdata("function.tsv")).unwrap() };
    // testEmpty
    for (d, f) in [(None, None), (None, Some(&functions)), (Some(&drugs), None), (Some(&drugs), Some(&functions))] {
        assert!(extract_genotypes(&BTreeMap::new(), d, f).unwrap().is_empty());
    }
    // testNonEmptyNoDrugOrFunctionInfo, testNonEmptyNoMatchingDrugOrFunction
    assert_eq!(extract_genotypes(&fake_analyses(), None, None).unwrap(), fake_genotypes());
    assert_eq!(extract_genotypes(&fake_analyses(), Some(&drugs), Some(&functions)).unwrap(), fake_genotypes());
    // testNonEmptyWithMatchingDrugAndFunction
    let mut m = fake_analyses();
    m.insert("DPYD".to_string(), analysis(&[("EVENT_DPYD_1", Some(1))], vec![c(&[("*2A", 2)])], "*1", PeachQcStatus::Pass, Some(c(&[("*2A", 2)]))));
    m.insert(
        "UGT1A1".to_string(),
        analysis(&[("EVENT_UGT1A1_1", Some(1))], vec![c(&[("*28", 1), ("*1", 1)])], "*1", PeachQcStatus::Pass, Some(c(&[("*28", 1), ("*1", 1)]))),
    );
    let mut expected = vec![g(
        "DPYD",
        "*2A",
        2,
        "No Function",
        "5-Fluorouracil;Capecitabine;Tegafur",
        "https://www.pharmgkb.org/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/guidelineAnnotation/PA166104963;https://www.pharmgkb.org/guidelineAnnotation/PA166104944",
    )];
    expected.extend(fake_genotypes());
    expected.push(g("UGT1A1", "*1", 1, "Normal Function", "Irinotecan", "https://www.pharmgkb.org/guidelineAnnotation/PA166104951"));
    expected.push(g("UGT1A1", "*28", 1, "Reduced Function", "Irinotecan", "https://www.pharmgkb.org/guidelineAnnotation/PA166104951"));
    assert_eq!(extract_genotypes(&m, Some(&drugs), Some(&functions)).unwrap(), expected);
}
