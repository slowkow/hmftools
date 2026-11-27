//! VirusInterpreterAlgoTest, VirusReportingDbModelTest, VirusBlacklistingDbFileTest, VirusReportingDbFileTest,
//! CoveragesAnalyzerTest, TaxonomyDbTest and TaxonomyDbFileTest

use std::collections::HashMap;

use super::testdata;
use crate::algo::{expected_clonal_coverage, has_acceptable_purple_quality, VirusInterpreterAlgo};
use crate::files::*;

fn statuses(s: &[&str]) -> Vec<String> {
    s.iter().map(|x| x.to_string()).collect()
}

fn reporting(interp: &'static str, integrated: Option<i32>, non_integrated: Option<i32>, l: VirusLikelihoodType) -> VirusReportingDb {
    VirusReportingDb {
        virus_interpretation: interp,
        integrated_minimal_coverage: integrated,
        non_integrated_minimal_coverage: non_integrated,
        virus_driver_likelihood_type: l,
    }
}

struct TestAlgo {
    taxonomy: TaxonomyDb,
    blacklist: Vec<i32>,
    reporting: VirusReportingDbModel,
}

impl TestAlgo {
    fn new(name: &str) -> TestAlgo {
        let mut species = HashMap::new();
        species.insert(1, reporting("EBV", None, Some(90), VirusLikelihoodType::High));
        species.insert(2, reporting("EBV", None, None, VirusLikelihoodType::High));
        species.insert(3, reporting("EBV", None, None, VirusLikelihoodType::High));
        let mut names = HashMap::new();
        names.insert(1, name.to_string());
        TestAlgo { taxonomy: TaxonomyDb { names, count: 1 }, blacklist: vec![11646], reporting: VirusReportingDbModel { species } }
    }
    fn algo(&self) -> VirusInterpreterAlgo<'_> {
        VirusInterpreterAlgo { taxonomy_db: &self.taxonomy, blacklisted_taxids: &self.blacklist, reporting_db: &self.reporting, expected_clonal_coverage: 34.5 }
    }
}

fn breakend(genus: i32, species: i32, integrations: i32, coverage: f64) -> VirusBreakend {
    VirusBreakend { reference_taxid: 1, taxid_genus: genus, taxid_species: species, integrations, coverage, ..Default::default() }
}

#[test]
fn can_test_is_purple_qc_pass() {
    assert!(has_acceptable_purple_quality(&statuses(&["WARN_DELETED_GENES"])));
    assert!(has_acceptable_purple_quality(&statuses(&["PASS"])));
    assert!(!has_acceptable_purple_quality(&statuses(&["FAIL_CONTAMINATION"])));
    assert!(!has_acceptable_purple_quality(&statuses(&["FAIL_NO_TUMOR"])));
}

#[test]
fn can_analyze_virus_breakends() {
    let name = "Human gammaherpesvirus 4";
    let t = TestAlgo::new(name);
    let algo = t.algo();
    let mut low_cov = breakend(2, 1, 2, 0.0);
    low_cov.qc_status = "LOW_VIRAL_COVERAGE";
    let mut deep = breakend(2, 1, 0, 91.0);
    deep.mean_depth = 35.0;
    let breakends = vec![
        breakend(2, 1, 2, 0.0),
        breakend(1, 4, 2, 0.0),
        low_cov,
        deep,
        breakend(2, 1, 0, 90.0),
        breakend(2, 2, 0, 91.0),
        breakend(2, 2, 0, 90.0),
    ];
    let qc = hmf_common::purple_files::read_purple_qc(&testdata("genomic/sample.purple.qc")).unwrap();
    let viruses = algo.analyze(&breakends, &qc.status).unwrap();
    assert_eq!(viruses.len(), 7);
    let reported: Vec<&AnnotatedVirus> = viruses.iter().filter(|v| v.reported).collect();
    assert_eq!(reported.len(), 4);
    assert_eq!(reported[0].name, name);
    assert_eq!(reported[0].integrations, 2);
    assert_eq!(reported[0].interpretation, Some("EBV"));
    for r in &reported[1..] {
        assert_eq!(r.name, name);
        assert_eq!(r.integrations, 0);
        assert_eq!(r.interpretation, Some("EBV"));
    }
    assert_eq!(algo.determine_minimal_coverage_virus(0, 1), Some(90));
    assert_eq!(algo.determine_minimal_coverage_virus(1, 1), None);
}

#[test]
fn can_determine_high_risk_virus() {
    let t = TestAlgo::new("");
    let algo = t.algo();
    assert_eq!(algo.virus_likelihood_type(&breakend(2, 1, 2, 0.0), true), VirusLikelihoodType::High);
    assert_eq!(algo.virus_likelihood_type(&breakend(2, 2, 2, 0.0), true), VirusLikelihoodType::High);
    assert_eq!(algo.virus_likelihood_type(&breakend(2, 3, 2, 0.0), false), VirusLikelihoodType::Unknown);
}

#[test]
fn can_report_test() {
    let t = TestAlgo::new("");
    let algo = t.algo();
    assert!(algo.report(&breakend(2, 1, 1, 0.0), 1.0, &statuses(&["FAIL_NO_TUMOR"])));
    assert!(algo.report(&breakend(2, 1, 1, 0.0), 1.0, &statuses(&["FAIL_CONTAMINATION"])));
    assert!(!algo.report(&breakend(2, 1, 0, 0.0), 1.0, &statuses(&["FAIL_NO_TUMOR"])));
    assert!(!algo.report(&breakend(2, 1, 0, 0.0), 1.0, &statuses(&["FAIL_CONTAMINATION"])));
}

#[test]
fn can_determine_blacklisted_virus() {
    let t = TestAlgo::new("");
    let algo = t.algo();
    assert!(algo.blacklist(11646, 0));
    assert!(!algo.blacklist(10000, 0));
}

#[test]
fn can_interpret_virus() {
    let mut species = HashMap::new();
    species.insert(1, reporting("EBV", None, None, VirusLikelihoodType::High));
    let model = VirusReportingDbModel { species };
    assert_eq!(model.interpret_virus_species(1), Some("EBV"));
    assert!(model.has_interpretation(1));
    assert!(!model.has_interpretation(3));
}

#[test]
fn can_read_virus_blacklisting_db_tsv() {
    assert_eq!(load_blacklisting_db(&testdata("virus_interpreter/virus_blacklisting_db.tsv")).unwrap().len(), 3);
}

#[test]
fn can_read_virus_reporting_db_tsv() {
    let model = load_reporting_db(&testdata("virus_interpreter/virus_reporting_db.tsv")).unwrap();
    assert_eq!(model.count(), 3);
    for (taxid, likelihood) in [(1, VirusLikelihoodType::Low), (2, VirusLikelihoodType::High), (3, VirusLikelihoodType::High)] {
        assert!(model.has_interpretation(taxid));
        assert_eq!(model.interpret_virus_species(taxid), Some("MCV"));
        assert_eq!(model.non_integrated_minimal_coverage(taxid), Some(90));
        assert_eq!(model.integrated_minimal_coverage(taxid), None);
        assert_eq!(model.virus_likelihood_type(taxid), likelihood);
    }
    assert!(!model.has_interpretation(4));
    assert_ne!(model.interpret_virus_species(2), Some("HPV"));
    assert_eq!(model.interpret_virus_species(5), None);
}

#[test]
fn can_read_real_virus_reporting_db_tsv() {
    assert_eq!(load_reporting_db(&testdata("virus_interpreter/real_virus_reporting_db.tsv")).unwrap().count(), 14);
}

#[test]
fn can_calculate_expected_clonal_coverage() {
    let dir = testdata("genomic");
    let purity = hmf_common::purple_files::read_purity_context(&dir, "sample").unwrap();
    let metrics = hmf_common::metrics::BamMetricSummary::read(&format!("{dir}/sample.bam_metric.summary.tsv")).unwrap();
    let e = expected_clonal_coverage(purity.purity, purity.ploidy, metrics.mean_coverage);
    assert!((e - 31.94).abs() < 0.1, "{e}");
}

#[test]
fn can_match_taxid_to_name() {
    let mut names = HashMap::new();
    names.insert(1, "species1".to_string());
    let db = TaxonomyDb { names, count: 1 };
    assert_eq!(db.count, 1);
    assert_eq!(db.lookup_name(1), "species1");
    assert_ne!(db.lookup_name(1), "species2");
    assert!(db.taxid_exists(1));
    assert!(!db.taxid_exists(3));
}

#[test]
fn can_read_taxonomy_db_tsv() {
    assert_eq!(load_taxonomy_db(&testdata("virus_interpreter/taxonomy_db.tsv"), None).unwrap().count, 1);
}
