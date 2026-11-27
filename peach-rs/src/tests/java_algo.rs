//! PeachAlgoTest, with TestHaplotypePanelFactory

use std::collections::{BTreeMap, BTreeSet};

use super::testdata;
use crate::algo::{gene_to_haplotype_analysis, EventCounts, HaplotypeAnalysis, HaplotypeCombination, PeachQcStatus};
use crate::panel::*;

fn ev(chr: &str, pos: i32, r: &str, a: &str) -> VariantHaplotypeEvent {
    VariantHaplotypeEvent::new(Chromosome::from_string(chr).unwrap(), pos, r, a)
}

fn nd(name: &str, wild: bool, events: Vec<VariantHaplotypeEvent>) -> NonDefaultHaplotype {
    NonDefaultHaplotype::new(name, wild, events).unwrap()
}

fn dpyd_panel() -> GeneHaplotypePanel {
    let wt = || ev("1", 98348885, "G", "A");
    GeneHaplotypePanel {
        default_haplotype: DefaultHaplotype { name: "*9A".into(), is_wild_type: false, events_to_ignore: vec![] },
        non_default_haplotypes: vec![
            nd("*1", true, vec![wt()]),
            nd("*2A", false, vec![ev("1", 97915614, "C", "T"), wt()]),
            nd("*7", false, vec![ev("1", 98205966, "GATGA", "G"), wt()]),
            nd("*13", false, vec![ev("1", 97981343, "A", "C"), wt()]),
            nd("*B3", false, vec![ev("1", 98039419, "C", "T"), ev("1", 98045449, "G", "C"), wt()]),
            nd("2846A>T", false, vec![ev("1", 97547947, "T", "A"), wt()]),
        ],
        wild_type_haplotype_name: "*1".into(),
    }
}

fn ugt1a1_panel() -> GeneHaplotypePanel {
    GeneHaplotypePanel {
        default_haplotype: DefaultHaplotype { name: "*1".into(), is_wild_type: true, events_to_ignore: vec![] },
        non_default_haplotypes: vec![
            nd("*28", false, vec![ev("2", 234668879, "C", "CAT")]),
            nd("*37", false, vec![ev("2", 234668879, "C", "CATAT")]),
            nd("*36", false, vec![ev("2", 234668879, "CAT", "C")]),
            nd("*6", false, vec![ev("2", 234669144, "G", "A")]),
        ],
        wild_type_haplotype_name: "*1".into(),
    }
}

pub fn default_test_panel() -> HaplotypePanel {
    let mut p = HaplotypePanel::default();
    p.genes.insert("DPYD".into(), dpyd_panel());
    p.genes.insert("UGT1A1".into(), ugt1a1_panel());
    p
}

fn single(gene: &str, panel: GeneHaplotypePanel) -> HaplotypePanel {
    let mut p = HaplotypePanel::default();
    p.genes.insert(gene.into(), panel);
    p
}

fn counts(pairs: &[(&str, Option<i32>)]) -> EventCounts {
    pairs.iter().map(|(e, c)| (e.to_string(), *c)).collect()
}

fn c(pairs: &[(&str, i32)]) -> HaplotypeCombination {
    HaplotypeCombination::of(pairs)
}

struct Expected<'a> {
    events: EventCounts,
    combinations: Vec<HaplotypeCombination>,
    default: &'a str,
    wild_type: &'a str,
    qc: PeachQcStatus,
    best: Option<HaplotypeCombination>,
}

fn assert_analysis(e: &Expected, a: &HaplotypeAnalysis) {
    assert_eq!(e.default, a.default_haplotype_name);
    assert_eq!(e.wild_type, a.wild_type_haplotype_name);
    assert_eq!(e.events, a.event_counts);
    assert_eq!(e.combinations.len(), a.haplotype_combinations.len());
    for ec in &e.combinations {
        assert_eq!(a.haplotype_combinations.iter().filter(|x| *x == ec).count(), 1, "combination {ec:?}");
    }
    assert_eq!(e.qc, a.qc_status);
    assert_eq!(e.best, a.best_haplotype_combination);
}

fn pass(default: &'static str, wild: &'static str, combo: HaplotypeCombination, events: EventCounts) -> Expected<'static> {
    Expected { events, combinations: vec![combo.clone()], default, wild_type: wild, qc: PeachQcStatus::Pass, best: Some(combo) }
}

fn ugt_unchanged() -> Expected<'static> {
    pass("*1", "*1", c(&[("*1", 2)]), EventCounts::new())
}

fn run(panel: &HaplotypePanel, events: &EventCounts) -> BTreeMap<String, HaplotypeAnalysis> {
    gene_to_haplotype_analysis(panel, events).unwrap()
}

#[test]
fn test_no_panel_or_events() {
    assert!(run(&HaplotypePanel::default(), &EventCounts::new()).is_empty());
}

#[test]
fn test_no_events() {
    let r = run(&default_test_panel(), &EventCounts::new());
    assert_eq!(r.len(), 2);
    assert_analysis(&pass("*9A", "*1", c(&[("*9A", 2)]), EventCounts::new()), &r["DPYD"]);
    assert_analysis(&ugt_unchanged(), &r["UGT1A1"]);
}

#[test]
fn test_homozygous_ref_events() {
    let events = counts(&[("VAR_chr1_98348885_G_A", Some(0)), ("VAR_chr1_98348885_G_T", Some(0))]);
    let r = run(&default_test_panel(), &events);
    assert_analysis(&pass("*9A", "*1", c(&[("*9A", 2)]), events.clone()), &r["DPYD"]);
    assert_analysis(&ugt_unchanged(), &r["UGT1A1"]);
}

#[test]
fn test_no_call_relevant_and_overlapping_event() {
    // testNoCallRelevantEvent and testNoCallOverlappingEvent (the same in Java)
    let events = counts(&[("VAR_chr1_98348885_G_T", None)]);
    let r = run(&default_test_panel(), &events);
    let e = Expected { events, combinations: vec![], default: "*9A", wild_type: "*1", qc: PeachQcStatus::FailEventWithUnknownCount, best: None };
    assert_analysis(&e, &r["DPYD"]);
    assert_analysis(&ugt_unchanged(), &r["UGT1A1"]);
}

#[test]
fn test_wild_types() {
    let events = counts(&[("VAR_chr1_98348885_G_A", Some(2))]);
    let r = run(&default_test_panel(), &events);
    assert_analysis(&pass("*9A", "*1", c(&[("*1", 2)]), events.clone()), &r["DPYD"]);
    assert_analysis(&ugt_unchanged(), &r["UGT1A1"]);
}

#[test]
fn test_heterozygous_non_wild_type() {
    let events = counts(&[
        ("VAR_chr1_98348885_G_A", Some(2)),
        ("VAR_chr1_98205966_GATGA_G", Some(1)),
        ("VAR_chr1_98039419_C_T", Some(1)),
        ("VAR_chr1_98045449_G_C", Some(1)),
        ("VAR_chr2_234668879_C_CAT", Some(1)),
        ("VAR_chr2_234668879_CAT_C", Some(1)),
    ]);
    let r = run(&default_test_panel(), &events);
    let dpyd_events = counts(&[
        ("VAR_chr1_98348885_G_A", Some(2)),
        ("VAR_chr1_98205966_GATGA_G", Some(1)),
        ("VAR_chr1_98039419_C_T", Some(1)),
        ("VAR_chr1_98045449_G_C", Some(1)),
    ]);
    assert_analysis(&pass("*9A", "*1", c(&[("*7", 1), ("*B3", 1)]), dpyd_events), &r["DPYD"]);
    let ugt_events = counts(&[("VAR_chr2_234668879_C_CAT", Some(1)), ("VAR_chr2_234668879_CAT_C", Some(1))]);
    assert_analysis(&pass("*1", "*1", c(&[("*28", 1), ("*36", 1)]), ugt_events), &r["UGT1A1"]);
}

#[test]
fn test_heterozygous_wild_type() {
    let events = counts(&[("VAR_chr1_98348885_G_A", Some(1)), ("VAR_chr2_234669144_G_A", Some(1))]);
    let r = run(&default_test_panel(), &events);
    assert_analysis(&pass("*9A", "*1", c(&[("*9A", 1), ("*1", 1)]), counts(&[("VAR_chr1_98348885_G_A", Some(1))])), &r["DPYD"]);
    assert_analysis(&pass("*1", "*1", c(&[("*1", 1), ("*6", 1)]), counts(&[("VAR_chr2_234669144_G_A", Some(1))])), &r["UGT1A1"]);
}

#[test]
fn test_heterozygous_default_haplotype() {
    let events = counts(&[("VAR_chr1_98348885_G_A", Some(1)), ("VAR_chr1_97915614_C_T", Some(1)), ("VAR_chr2_234668879_C_CATAT", Some(1))]);
    let r = run(&default_test_panel(), &events);
    let dpyd_events = counts(&[("VAR_chr1_98348885_G_A", Some(1)), ("VAR_chr1_97915614_C_T", Some(1))]);
    assert_analysis(&pass("*9A", "*1", c(&[("*2A", 1), ("*9A", 1)]), dpyd_events), &r["DPYD"]);
    assert_analysis(&pass("*1", "*1", c(&[("*1", 1), ("*37", 1)]), counts(&[("VAR_chr2_234668879_C_CATAT", Some(1))])), &r["UGT1A1"]);
}

#[test]
fn test_homozygous_non_wild_type() {
    let events = counts(&[("VAR_chr1_97981343_A_C", Some(2)), ("VAR_chr1_98348885_G_A", Some(2)), ("VAR_chr2_234668879_C_CAT", Some(2))]);
    let r = run(&default_test_panel(), &events);
    let dpyd_events = counts(&[("VAR_chr1_98348885_G_A", Some(2)), ("VAR_chr1_97981343_A_C", Some(2))]);
    assert_analysis(&pass("*9A", "*1", c(&[("*13", 2)]), dpyd_events), &r["DPYD"]);
    assert_analysis(&pass("*1", "*1", c(&[("*28", 2)]), counts(&[("VAR_chr2_234668879_C_CAT", Some(2))])), &r["UGT1A1"]);
}

#[test]
fn test_too_many_alleles_ugt1a1() {
    let events = counts(&[("VAR_chr2_234668879_C_CAT", Some(2)), ("VAR_chr2_234669144_G_A", Some(2))]);
    let r = run(&default_test_panel(), &events);
    assert_analysis(&pass("*9A", "*1", c(&[("*9A", 2)]), EventCounts::new()), &r["DPYD"]);
    let combo = c(&[("*28", 2), ("*6", 2)]);
    let e = Expected { events, combinations: vec![combo.clone()], default: "*1", wild_type: "*1", qc: PeachQcStatus::WarnTooManyAllelesFound, best: Some(combo) };
    assert_analysis(&e, &r["UGT1A1"]);
}

#[test]
fn test_overlapping_genes() {
    let panel = load_haplotype_panel(&testdata("haplotypes.overlap.tsv")).unwrap();
    let events = counts(&[("VAR_chr3_1234567_CAT_C", Some(1)), ("VAR_chr3_1234567_C_CAT", Some(1))]);
    let r = run(&panel, &events);
    assert_eq!(r.len(), 3);
    assert_analysis(&pass("*1", "*1", c(&[("*2", 1), ("*4", 1)]), events.clone()), &r["FAKE1"]);
    assert_analysis(&pass("*1", "*1", c(&[("*2", 1), ("*1", 1)]), counts(&[("VAR_chr3_1234567_C_CAT", Some(1))])), &r["FAKE2"]);
    let e = Expected { events, combinations: vec![], default: "*1", wild_type: "*1", qc: PeachQcStatus::FailNoCombinationFound, best: None };
    assert_analysis(&e, &r["FAKE3"]);
}

#[test]
fn test_ambiguous_calling_prefer_fewer_non_wild_type_haplotypes() {
    let panel = load_haplotype_panel(&testdata("haplotypes.ambiguous.tsv")).unwrap();
    let events = counts(&[("VAR_chr3_1234567_C_CAT", Some(2)), ("VAR_chr3_1234765_G_C", Some(1))]);
    let r = run(&panel, &events);
    let e = Expected {
        events,
        combinations: vec![c(&[("*2", 1), ("*3", 1)]), c(&[("*3", 2), ("*4", 1)])],
        default: "*1",
        wild_type: "*1",
        qc: PeachQcStatus::Pass,
        best: Some(c(&[("*2", 1), ("*3", 1)])),
    };
    assert_analysis(&e, &r["FAKE1"]);
}

#[test]
fn test_ambiguous_calling_prefer_fewer_haplotypes() {
    let e1 = ev("1", 100, "A", "C");
    let e2 = ev("1", 200, "A", "ATT");
    let gene = GeneHaplotypePanel {
        default_haplotype: DefaultHaplotype { name: "*2".into(), is_wild_type: false, events_to_ignore: vec![] },
        non_default_haplotypes: vec![nd("*1", true, vec![e1.clone()]), nd("*3", true, vec![e2.clone()]), nd("*4", true, vec![e1.clone(), e2.clone()])],
        wild_type_haplotype_name: "*1".into(),
    };
    let events = counts(&[(&e1.id(), Some(2)), (&e2.id(), Some(2))]);
    let r = run(&single("FAKE1", gene), &events);
    let e = Expected {
        events: counts(&[("VAR_chr1_100_A_C", Some(2)), ("VAR_chr1_200_A_ATT", Some(2))]),
        combinations: vec![c(&[("*1", 2), ("*3", 2)]), c(&[("*1", 1), ("*3", 1), ("*4", 1)]), c(&[("*4", 2)])],
        default: "*2",
        wild_type: "*1",
        qc: PeachQcStatus::Pass,
        best: Some(c(&[("*4", 2)])),
    };
    assert_analysis(&e, &r["FAKE1"]);
}

#[test]
fn test_ambiguous_calling_no_best_haplotype() {
    let panel = load_haplotype_panel(&testdata("haplotypes.ambiguous.tsv")).unwrap();
    let events = counts(&[("VAR_chr3_1234567_C_CAT", Some(1)), ("VAR_chr3_1234765_G_C", Some(1)), ("VAR_chr3_1234777_AA_T", Some(1))]);
    let r = run(&panel, &events);
    let e = Expected {
        events,
        combinations: vec![c(&[("*2", 1), ("*5", 1)]), c(&[("*3", 1), ("*6", 1)]), c(&[("*3", 1), ("*4", 1), ("*5", 1)])],
        default: "*1",
        wild_type: "*1",
        qc: PeachQcStatus::FailNoUniqueBestCombinationFound,
        best: None,
    };
    assert_analysis(&e, &r["FAKE1"]);
}

#[test]
fn test_ambiguous_calling_prefer_wild_type() {
    let panel = load_haplotype_panel(&testdata("haplotypes.ambiguous.tsv")).unwrap();
    let events = counts(&[("VAR_chr3_1234567_C_CAT", Some(1)), ("VAR_chr3_1234765_G_C", Some(1))]);
    let r = run(&panel, &events);
    let e = Expected {
        events,
        combinations: vec![c(&[("*2", 1), ("*1", 1)]), c(&[("*3", 1), ("*4", 1)])],
        default: "*1",
        wild_type: "*1",
        qc: PeachQcStatus::Pass,
        best: Some(c(&[("*2", 1), ("*1", 1)])),
    };
    assert_analysis(&e, &r["FAKE1"]);
}

/// haplotype name to its event ids, the default haplotype with none
fn haplotype_events(gene: &str, panel: &HaplotypePanel) -> BTreeMap<String, Vec<String>> {
    let g = &panel.genes[gene];
    let mut m: BTreeMap<String, Vec<String>> =
        g.non_default_haplotypes.iter().map(|h| (h.name.clone(), h.events.iter().map(|e| e.id()).collect())).collect();
    m.insert(g.default_haplotype.name.clone(), vec![]);
    m
}

fn all_homozygous(gene: &str, panel: HaplotypePanel) {
    for (name, ids) in haplotype_events(gene, &panel) {
        let events: EventCounts = ids.iter().map(|e| (e.clone(), Some(2))).collect();
        let r = run(&panel, &events);
        assert_eq!(r.len(), 1);
        assert_eq!(r[gene].best_haplotype_combination, Some(c(&[(&name, 2)])), "{name}");
    }
}

fn all_heterozygous(gene: &str, panel: HaplotypePanel) {
    let m = haplotype_events(gene, &panel);
    let names: Vec<&String> = m.keys().collect();
    for i in 0..names.len() {
        for j in i + 1..names.len() {
            let mut n: BTreeMap<String, i32> = BTreeMap::new();
            for id in m[names[i]].iter().chain(m[names[j]].iter()) {
                *n.entry(id.clone()).or_insert(0) += 1;
            }
            let events: EventCounts = n.into_iter().map(|(k, v)| (k, Some(v))).collect();
            let r = run(&panel, &events);
            assert_eq!(r[gene].best_haplotype_combination, Some(c(&[(names[i], 1), (names[j], 1)])), "{} {}", names[i], names[j]);
        }
    }
}

#[test]
fn test_simple_homozygous_dpyd_haplotype_combinations() {
    all_homozygous("DPYD", single("DPYD", dpyd_panel()));
}

#[test]
fn test_simple_heterozygous_dpyd_haplotype_combinations() {
    all_heterozygous("DPYD", single("DPYD", dpyd_panel()));
}

#[test]
fn test_simple_homozygous_ugt1a1_haplotype_combinations() {
    all_homozygous("UGT1A1", single("UGT1A1", ugt1a1_panel()));
}

#[test]
fn test_simple_heterozygous_ugt1a1_haplotype_combinations() {
    all_heterozygous("UGT1A1", single("UGT1A1", ugt1a1_panel()));
}

#[test]
fn relevant_positions_of_the_test_panel() {
    // HaplotypePanelTest.testGetRelevantVariantPositionsSameChromosome
    let mut p = HaplotypePanel::default();
    for (gene, e) in [("FAKE1", ev("1", 1000, "A", "T")), ("FAKE2", ev("1", 3000, "AG", "A"))] {
        p.genes.insert(
            gene.into(),
            GeneHaplotypePanel {
                default_haplotype: DefaultHaplotype { name: "*1".into(), is_wild_type: true, events_to_ignore: vec![] },
                non_default_haplotypes: vec![nd("*2", false, vec![e])],
                wild_type_haplotype_name: "*1".into(),
            },
        );
    }
    let m = p.relevant_variant_positions();
    assert_eq!(m.len(), 1);
    let positions: BTreeSet<i32> = m[&Chromosome::from_string("1").unwrap()].iter().copied().collect();
    assert_eq!(positions, [1000, 3000, 3001].into_iter().collect());
}
