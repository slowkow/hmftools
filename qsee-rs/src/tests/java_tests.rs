//! The Java unit tests of qsee 1.0 (qsee/src/test/java), ported. Each test names its Java class and
//! method. Java's `Map.of` depth counts have no defined order; the port keeps them as lists in the order
//! written here, which changes only the summation order (the tests' tolerances cover it).

use hmf_common::driver_panel::DriverGene;
use hmf_common::metrics::{BamFlagStats, BamMetricSummary, DepthRange, ExonDepth, GeneDepth, ValueFrequency};
use hmf_common::purple_files::{PurplePurity, PurpleQc};
use hmf_common::redux_files::{BqrRecord, ConsensusType, JitterTableRow};

use crate::cohort::{FeatureMatrix, PercentileTransformer};
use crate::common::{BinnedFrequencies, SampleIds, SampleType};
use crate::config::SequencingType;
use crate::feature::{BaseQualBin, DiscordantFragGroup, Feature, FeatureKey, FeatureType, NumberFormat, SourceTool, SummaryTableFeature};
use crate::coverage::{
    calc_exon_exclusions, coverage_qc_per_exon, min_sample_vaf_to_detect_variant, ExonCoverageQc, ExonExclusions, VariantDetectabilityStatus, DEFAULT_LOD_MIN_ALT_READS,
    LOD_VAF_LOWER_BOUND,
};
use crate::prep::{self, BamMetricsData, BaseQualBinner, PurityContext, DELETIONS_KEY, INSERTIONS_KEY};
use crate::status::{purple_to_qsee_status, try_form_display_string, ComparisonOperator, QcThreshold, QcStatusType, ThresholdKey, ThresholdRegistry, PURPLE_QC_STATUSES};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol || (a.is_nan() && b.is_nan()) || a == b
}

fn dr(lower: i32, upper: i32) -> DepthRange {
    DepthRange { lower, upper }
}

// ---- coverage.ExonCoverageQcTest (upstream bdceb2302a; GeneLodCalculatorTest before) ----

fn exon(gene: &str, rank: i32, median: f64, mean: f64, counts: Vec<(DepthRange, i32)>) -> ExonDepth {
    ExonDepth { depth: GeneDepth { gene: gene.into(), chromosome: "1".into(), pos_start: 100, pos_end: 199, depth_counts: counts }, exon_rank: rank, median_depth: median, mean_depth: mean }
}

fn exclusions(entries: &[(&str, i32)]) -> ExonExclusions {
    let mut e = ExonExclusions::empty();
    for (g, r) in entries {
        e.add(g, *r);
    }
    e
}

#[test]
fn exon_coverage_qc_can_calc_coverage_qc_per_exon() {
    let mixed = vec![(dr(20, 20), 10), (dr(100, 100), 90)];
    let low = vec![(dr(2, 2), 10), (dr(100, 100), 90)];
    let high = vec![(dr(500, 500), 100)];
    let exons = vec![
        exon("REPORTABLE_GENE", 1, 100.0, 92.0, mixed.clone()),
        exon("REPORTABLE_GENE", 2, 100.0, 90.2, low),
        exon("REPORTABLE_GENE", 3, 100.0, 500.0, high),
        exon("REPORTABLE_GENE", 4, 100.0, 92.0, mixed.clone()),
        exon("NON_REPORTABLE_GENE", 1, 100.0, 92.0, mixed),
    ];
    let expected: [(i32, f64, VariantDetectabilityStatus); 4] = [
        (1, 0.132, VariantDetectabilityStatus::Warn),
        (2, 1.0, VariantDetectabilityStatus::Warn), // depth too low
        (3, 0.0125, VariantDetectabilityStatus::Pass),
        (4, 0.132, VariantDetectabilityStatus::Excluded), // excluded despite a high LoD VAF
    ];
    let qc = coverage_qc_per_exon(
        "TUMOR_SAMPLE", SampleType::Tumor, &exons, &exclusions(&[("REPORTABLE_GENE", 4)]), &["REPORTABLE_GENE".to_string()], &ThresholdRegistry::create_default(false),
        DEFAULT_LOD_MIN_ALT_READS,
    )
    .unwrap();
    assert_eq!(qc.len(), 4); // the non-reportable gene left out
    for (rank, lod, status) in expected {
        let g = qc.iter().find(|g| g.exon_rank == rank).unwrap();
        assert_eq!(g.gene_name, "REPORTABLE_GENE");
        assert_eq!(g.median_depth, 100.0);
        assert!(close(g.lod_vaf, lod, 1e-3), "{} vs {lod}", g.lod_vaf);
        assert_eq!(g.status, status);
    }
}

// ---- coverage.ExonExclusionsFileBuilderTest (upstream bdceb2302a) ----

#[test]
fn exon_exclusions_from_cohort() {
    let sample_count = 20;
    // exon rank, tumor samples with a high LoD VAF, the high LoD VAF, normal samples with it, excluded
    let rows: [(i32, usize, f64, usize, bool); 4] = [(1, 6, 0.07, 0, true), (2, 5, 0.07, 0, false), (3, 20, 0.069, 0, false), (4, 0, 0.07, 20, false)];
    let mut qc = Vec::new();
    for (rank, tumor_high, high_vaf, normal_high, _) in rows {
        for (st, n_high) in [(SampleType::Tumor, tumor_high), (SampleType::Normal, normal_high)] {
            for i in 0..sample_count {
                qc.push(ExonCoverageQc {
                    sample_id: format!("{}_{i}", st.name()),
                    sample_type: st,
                    gene_name: "GENE".into(),
                    exon_rank: rank,
                    median_depth: 100.0,
                    lod_vaf: if i < n_high { high_vaf } else { 0.01 },
                    status: VariantDetectabilityStatus::Pass,
                    depth_counts: Vec::new(),
                });
            }
        }
    }
    let ex = calc_exon_exclusions(&qc, 0.07, 0.3);
    for (rank, _, _, _, excluded) in rows {
        assert_eq!(ex.is_excluded("GENE", rank), excluded, "exon {rank}");
    }
}

// ---- coverage.LimitOfDetectionCalcsTest ----

#[test]
fn limit_of_detection_can_calc_lod_vaf() {
    let counts = vec![(dr(20, 20), 10), (dr(100, 100), 90)];
    assert!(close(min_sample_vaf_to_detect_variant(&counts, DEFAULT_LOD_MIN_ALT_READS), 0.132, 1e-3));
    let low = vec![(dr(2, 2), 10), (dr(100, 100), 90)];
    assert_eq!(min_sample_vaf_to_detect_variant(&low, DEFAULT_LOD_MIN_ALT_READS), 1.0);
    // a depth high enough for a LoD VAF below the floor
    let high = vec![(dr(1000, 1000), 100)];
    assert_eq!(min_sample_vaf_to_detect_variant(&high, DEFAULT_LOD_MIN_ALT_READS), LOD_VAF_LOWER_BOUND);
}

// ---- status.PurpleQCStatusConverterTest ----

#[test]
fn purple_qc_status_converter_maps_every_status() {
    let t = ThresholdRegistry::create_without_thresholds();
    let statuses: Vec<_> = PURPLE_QC_STATUSES.iter().map(|s| purple_to_qsee_status(&t, s).unwrap()).collect();
    assert_eq!(statuses.len(), PURPLE_QC_STATUSES.len());
}

// ---- status.QcStatusTest ----

#[test]
fn qc_status_formats_percent_thresholds() {
    let f = |v: f64, nf: NumberFormat| try_form_display_string(QcStatusType::Fail, Some(ComparisonOperator::LessThan), v, nf).unwrap();
    assert_eq!(f(0.10, NumberFormat::Number), "<0.1");
    assert_eq!(f(0.10, NumberFormat::Percent), "<10%");
    assert_eq!(f(0.101, NumberFormat::Percent), "<10.1%");
}

#[test]
fn qc_status_formats_no_thresholds() {
    let f = |nf: NumberFormat| try_form_display_string(QcStatusType::Fail, Some(ComparisonOperator::LessThan), f64::NAN, nf).unwrap();
    assert_eq!(f(NumberFormat::Percent), "");
    assert_eq!(f(NumberFormat::Number), "");
}

// ---- status.ThresholdRegistryTest ----

#[test]
fn threshold_registry_can_override_thresholds_from_file() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/threshold_overrides.tsv");
    let t = crate::status::read_threshold_overrides(path, false).unwrap();
    let key = ThresholdKey::new(SampleType::Tumor, FeatureType::SummaryTable, SummaryTableFeature::MappedProportion.name(), QcStatusType::Fail);
    let expected = QcThreshold { key: key.clone(), operator: Some(ComparisonOperator::LessThan), threshold: f64::NEG_INFINITY, determined_elsewhere: false };
    assert_eq!(*t.get(&key).unwrap(), expected);
    let key = ThresholdKey::new(SampleType::Tumor, FeatureType::SummaryTable, SummaryTableFeature::LowMapQual.name(), QcStatusType::Warn);
    let expected = QcThreshold { key: key.clone(), operator: Some(ComparisonOperator::GreaterThan), threshold: f64::NAN, determined_elsewhere: false };
    assert_eq!(*t.get(&key).unwrap(), expected);
}

// ---- common.BinnedFrequenciesTest ----

fn binned() -> BinnedFrequencies {
    BinnedFrequencies::new(vec![0, 10, 25, 50, 100], vec![200.0; 5]).unwrap()
}

#[test]
fn binned_frequencies_bin_ends() {
    assert_eq!(binned().bin_ends, vec![10, 25, 50, 100, i64::MAX]);
}

#[test]
fn binned_frequencies_densities() {
    let expected = [200.0 / 10.0, 200.0 / 15.0, 200.0 / 25.0, 200.0 / 50.0, 200.0];
    for (a, e) in binned().frequency_densities().iter().zip(expected) {
        assert!(close(*a, e, 0.01));
    }
}

#[test]
fn binned_frequencies_proportional_densities() {
    let expected = [0.02, 0.01334, 0.008, 0.004, 0.2];
    for (a, e) in binned().proportional_densities().iter().zip(expected) {
        assert!(close(*a, e, 0.001));
    }
}

// ---- common.SampleIdsLoaderTest ----

fn lines(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn sample_ids_from_strings() {
    let s = SampleIds::from_strings(Some("TUMOR001,TUMOR002"), Some("REF001,REF002")).unwrap();
    assert_eq!((s.tumor_ids.len(), s.reference_ids.len()), (2, 2));
    let s = SampleIds::from_strings(Some("TUMOR001,TUMOR002"), None).unwrap();
    assert_eq!((s.tumor_ids.len(), s.reference_ids.len()), (2, 0));
    assert!(SampleIds::from_strings(Some("TUMOR001,TUMOR002"), Some("REF001")).is_err());
    assert!(SampleIds::from_strings(None, Some("REF001,REF002")).is_err());
}

#[test]
fn sample_ids_from_lines() {
    let s = SampleIds::from_lines(&lines(&["TumorId\tReferenceId", "TUMOR001\tREF001", "TUMOR002\tREF002"])).unwrap();
    assert_eq!((s.tumor_ids.len(), s.reference_ids.len()), (2, 2));
    let s = SampleIds::from_lines(&lines(&["TumorId", "TUMOR001", "TUMOR002"])).unwrap();
    assert_eq!((s.tumor_ids.len(), s.reference_ids.len()), (2, 0));
    assert!(SampleIds::from_lines(&lines(&["TumorId\tReferenceId", "TUMOR001\tREF001", "TUMOR002\t"])).is_err());
    assert!(SampleIds::from_lines(&lines(&["TumorId\tReferenceId", "TUMOR001\t", "\tREF002"])).is_err());
}

// ---- cohort.FeatureMatrixTest ----

fn test_feature(name: &str, value: f64) -> Feature {
    Feature::new(FeatureKey { name: name.into(), feature_type: None, source_tool: None }, value)
}

#[test]
fn feature_matrix_adds_rows_incrementally() {
    let mut m = FeatureMatrix::with_num_samples(4);
    m.add_row("sample1", vec![test_feature("feature1", 1.1), test_feature("feature2", 1.2)]).unwrap();
    m.add_row("sample2", vec![test_feature("feature1", 2.1), test_feature("feature2", 2.2)]).unwrap();
    m.add_row("sample3", vec![test_feature("feature1", 3.1), test_feature("feature2", 3.2), test_feature("feature3", 3.3)]).unwrap();
    m.add_row("sample4", vec![test_feature("feature1", 4.1), test_feature("feature2", 4.2), test_feature("feature3", 4.3)]).unwrap();
    assert_eq!(m.sample_ids, vec!["sample1", "sample2", "sample3", "sample4"]);
    assert_eq!(m.keys.iter().map(|k| k.name.as_str()).collect::<Vec<_>>(), vec!["feature1", "feature2", "feature3"]);
    let v = m.feature_values();
    let expected = [[1.1, 1.2, f64::NAN], [2.1, 2.2, f64::NAN], [3.1, 3.2, 3.3], [4.1, 4.2, 4.3]];
    for (row, e) in v.iter().zip(expected) {
        for (a, b) in row.iter().zip(e) {
            assert!(close(*a, b, 0.0));
        }
    }
}

/// canAddRowsConcurrently: Java adds the rows from 1000 threads; the port adds rows in sample order, and
/// this checks the result for rows added in a shuffled order
#[test]
fn feature_matrix_adds_rows_in_any_order() {
    let ids: Vec<String> = (0..1000).map(|i| format!("sample{i}")).collect();
    let mut m = FeatureMatrix::with_samples(&ids);
    let mut order: Vec<usize> = (0..1000).collect();
    let mut seed = 7u64;
    for i in (1..order.len()).rev() {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        order.swap(i, (seed >> 33) as usize % (i + 1));
    }
    for s in order {
        let features: Vec<Feature> = (0..500).map(|f| test_feature(&format!("feature{f}"), s as f64 + 0.1 * f as f64)).collect();
        m.add_row(&ids[s], features).unwrap();
    }
    assert_eq!(m.sample_ids, ids);
    assert_eq!(m.keys.len(), 500);
    assert!(m.keys.iter().enumerate().all(|(i, k)| k.name == format!("feature{i}")));
    let v = m.feature_values();
    for s in 0..1000 {
        for f in 0..500 {
            assert_eq!(v[s][f], s as f64 + 0.1 * f as f64);
        }
    }
}

// ---- cohort.PercentileTransformerTest ----

const COHORT_VALUES: [f64; 12] = [0.0, 0.0, 5.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, f64::NAN, f64::INFINITY];

fn default_transformer() -> PercentileTransformer {
    let mut t = PercentileTransformer::with_increment(10.0).unwrap();
    t.fit(&COHORT_VALUES).unwrap();
    t
}

fn arrays_close(a: &[f64], b: &[f64], tol: f64) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| close(*x, *y, tol))
}

#[test]
fn percentile_transformer_from_num_or_increment() {
    let e = [0.0, 25.0, 50.0, 75.0, 100.0];
    assert!(arrays_close(&PercentileTransformer::with_num_percentiles(5).unwrap().percentiles, &e, 0.1));
    assert!(arrays_close(&PercentileTransformer::with_increment(25.0).unwrap().percentiles, &e, 0.1));
}

#[test]
fn percentile_transformer_nan_and_inf_ignored() {
    let with_nan = default_transformer();
    let mut without = PercentileTransformer::with_increment(10.0).unwrap();
    let values: Vec<f64> = COHORT_VALUES.iter().copied().filter(|v| !v.is_nan()).collect();
    without.fit(&values).unwrap();
    assert!(arrays_close(&with_nan.percentiles, &without.percentiles, 0.001));
    assert!(arrays_close(&with_nan.ref_values, &without.ref_values, 0.001));
}

#[test]
fn percentile_transformer_dedup() {
    let t = default_transformer();
    assert!(arrays_close(&t.percentiles_deduped, &[5.0, 20.0, 30.0, 70.0], 0.001));
    assert!(arrays_close(&t.ref_values_deduped, &[0.0, 4.0, 8.5, 10.0], 0.001));
}

#[test]
fn percentile_transformer_percentile_to_value() {
    let t = default_transformer();
    assert!(close(t.percentile_to_feature_value(0.0).unwrap(), 0.0, 0.001));
    assert!(close(t.percentile_to_feature_value(25.0).unwrap(), (4.0 + 8.5) / 2.0, 0.001));
}

#[test]
fn percentile_transformer_value_to_percentile() {
    let t = default_transformer();
    for (v, p) in [(0.0, 5.0), (4.0, 20.0), (8.5, 30.0), (10.0, 70.0), (2.0, 12.5), (7.0, 26.67)] {
        assert!(close(t.feature_value_to_percentile(v).unwrap(), p, 0.01), "{v}");
    }
    assert_eq!(t.feature_value_to_percentile(1000.0).unwrap(), f64::INFINITY);
    assert_eq!(t.feature_value_to_percentile(-1000.0).unwrap(), f64::NEG_INFINITY);
}

#[test]
fn percentile_transformer_one_cohort_value() {
    let mut t = PercentileTransformer::with_num_percentiles(3).unwrap();
    t.fit(&[5.0]).unwrap();
    assert!(arrays_close(&t.ref_values_deduped, &[5.0], 0.1));
    assert!(arrays_close(&t.percentiles_deduped, &[50.0], 0.1));
}

#[test]
fn percentile_transformer_all_nan() {
    let mut t = PercentileTransformer::with_num_percentiles(3).unwrap();
    t.fit(&[f64::NAN, f64::NAN]).unwrap();
    assert!(t.feature_value_to_percentile(0.0).unwrap().is_nan());
}

#[test]
fn percentile_transformer_from_prefit() {
    let t = PercentileTransformer::from_prefit(&[0.0, 25.0, 50.0, 75.0, 100.0], &[0.0, 0.0, 10.0, 10.0, 10.0]).unwrap();
    assert!(arrays_close(&t.percentiles_deduped, &[12.5, 75.0], 0.01));
    assert!(arrays_close(&t.ref_values_deduped, &[0.0, 10.0], 0.01));
}

// ---- feature.FeatureKeyTest ----

#[test]
fn feature_key_names_sorted() {
    use FeatureType::*;
    use SourceTool::*;
    let k = |n: &str, t: FeatureType, s: SourceTool| FeatureKey::new(n, t, s);
    let mut keys = vec![
        k("PURITY", SummaryTable, Purple),
        k("COVERAGE_ABOVE_100", SummaryTable, BamMetrics),
        k("COVERAGE_ABOVE_20", SummaryTable, BamMetrics),
        k("DiscordantFragType:INV_SHORT", DiscordantFragFreq, Esvee),
        k("DiscordantFragType:DEL_DUP_MEDIUM", DiscordantFragFreq, Esvee),
        k("DiscordantFragType:DEL_DUP_SHORT", DiscordantFragFreq, Esvee),
        k("ReadCount:>=100", DuplicateFreq, Redux),
        k("ReadCount:10", DuplicateFreq, Redux),
        k("ReadCount:1", DuplicateFreq, Redux),
        k("ConsensusType:NONE;RepeatUnitType:>=3bp repeat;RefNumUnits:4", MsIndelErrorRates, Redux),
        k("ConsensusType:NONE;RepeatUnitType:2bp repeat;RefNumUnits:4", MsIndelErrorRates, Redux),
        k("ConsensusType:NONE;RepeatUnitType:C/G repeat;RefNumUnits:4", MsIndelErrorRates, Redux),
        k("ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:4", MsIndelErrorRates, Redux),
        k("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;OriginalQualBin:HIGH (30+)", BqrPerOrigQual, Redux),
        k("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;OriginalQualBin:LOW (0-29)", BqrPerOrigQual, Redux),
    ];
    keys.sort_by(|a, b| a.compare(b));
    let names: Vec<&str> = keys.iter().map(|k| k.name.as_str()).collect();
    assert_eq!(
        names,
        vec![
            "COVERAGE_ABOVE_20",
            "COVERAGE_ABOVE_100",
            "PURITY",
            "DiscordantFragType:DEL_DUP_SHORT",
            "DiscordantFragType:DEL_DUP_MEDIUM",
            "DiscordantFragType:INV_SHORT",
            "ReadCount:1",
            "ReadCount:10",
            "ReadCount:>=100",
            "ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;OriginalQualBin:LOW (0-29)",
            "ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;OriginalQualBin:HIGH (30+)",
            "ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:4",
            "ConsensusType:NONE;RepeatUnitType:C/G repeat;RefNumUnits:4",
            "ConsensusType:NONE;RepeatUnitType:2bp repeat;RefNumUnits:4",
            "ConsensusType:NONE;RepeatUnitType:>=3bp repeat;RefNumUnits:4",
        ]
    );
}

// ---- feature.category.SummaryTablePrepTest ----

fn test_coverage() -> Vec<ValueFrequency> {
    [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000].iter().map(|v| ValueFrequency::new(*v, 1000)).collect()
}

pub fn test_purity_context(statuses: &[&str]) -> PurityContext {
    let qc = PurpleQc {
        status: statuses.iter().map(|s| s.to_string()).collect(),
        method: "NORMAL".into(),
        copy_number_segments: 100,
        unsupported_copy_number_segments: 1,
        purity: 0.85,
        amber_gender: "FEMALE".into(),
        cobalt_gender: "FEMALE".into(),
        deleted_genes: 10,
        contamination: 0.01,
        germline_aberrations: vec!["NONE".into()],
        amber_mean_depth: 100,
        loh_percent: 0.01,
        tinc_level: 0.1,
    };
    let purity = PurplePurity {
        purity: 0.85,
        norm_factor: 1.0,
        score: 1.0,
        diploid_proportion: 0.95,
        ploidy: 2.0,
        sex: "FEMALE".into(),
        fit_method: "NORMAL".into(),
        polyclonal_proportion: 0.1,
        min_purity: 0.8,
        max_purity: 0.9,
        min_ploidy: 2.0,
        max_ploidy: 2.2,
        min_diploid_proportion: 0.75,
        max_diploid_proportion: 0.85,
        somatic_penalty: 0.0,
        whole_genome_duplication: false,
        ms_indels_per_mb: 1.0,
        ms_status: "MSS".into(),
        tml: 1000,
        tml_status: "HIGH".into(),
        tmb_per_mb: 10.0,
        tmb_status: "HIGH".into(),
        sv_tumor_mutational_burden: 10,
        run_mode: "TUMOR_GERMLINE".into(),
        targeted: false,
    };
    PurityContext { qc, purity }
}

#[test]
fn summary_table_prop_bases_above_coverage() {
    let t = ThresholdRegistry::create_without_thresholds();
    let f = prep::coverage_above_features(&test_coverage(), SampleType::Tumor, &t).unwrap();
    let get = |x: SummaryTableFeature| f.iter().find(|(k, _)| *k == x).unwrap().1.value;
    assert!(close(get(SummaryTableFeature::CoverageAbove10), 0.9, 0.01));
    assert!(close(get(SummaryTableFeature::CoverageAbove30), 0.7, 0.01));
    assert!(close(get(SummaryTableFeature::CoverageAbove100), 0.6, 0.01));
    assert!(close(get(SummaryTableFeature::CoverageAbove250), 0.4, 0.01));
}

#[test]
fn summary_table_highest_tinc_status_selected() {
    let ctx = test_purity_context(&["WARN_LOW_PURITY", "WARN_TINC", "FAIL_TINC"]);
    let t = ThresholdRegistry::create_without_thresholds();
    let m = prep::purple_qc_status_mappings(&ctx.qc.status, &t).unwrap();
    assert_eq!(prep::tinc_qc_status(&m).status_type, QcStatusType::Fail);
}

#[test]
fn summary_table_purple_features() {
    let ctx = test_purity_context(&["PASS"]);
    let t = ThresholdRegistry::create_without_thresholds();
    let f = prep::purple_features(&ctx, &t).unwrap();
    let expected = SummaryTableFeature::ALL.iter().filter(|f| f.source_tool() == SourceTool::Purple).count();
    assert_eq!(f.len(), expected);
}

#[test]
fn summary_table_bam_metrics_features() {
    let summary = BamMetricSummary {
        total_region_bases: 1_000_000_000,
        total_reads: 10_000_000,
        duplicate_reads: 1_000_000,
        dual_strand_reads: 1_000_000,
        off_target_reads: 0,
        mean_coverage: 100.0,
        sd_coverage: 10.0,
        median_coverage: 100,
        mad_coverage: 10,
        low_map_qual_percent: 0.01,
        duplicate_percent: 0.01,
        unmapped_percent: 0.01,
        low_base_qual_percent: 0.01,
        overlapping_read_percent: 0.01,
        capped_coverage_percent: 0.01,
        coverage_levels: vec![30, 60, 100],
        coverage_percents: vec![0.95, 0.9, 0.85],
    };
    let flag = BamFlagStats {
        unique_read_count: 0,
        secondary_count: 0,
        supplementary_count: 0,
        duplicate_proportion: 0.0,
        mapped_proportion: 0.0,
        paired_in_sequencing_proportion: 0.0,
        properly_paired_proportion: 0.0,
        with_itself_and_mate_mapped_proportion: 0.0,
        singleton_proportion: 0.0,
    };
    let data = BamMetricsData { summary: Some(summary), coverage: Some(test_coverage()), flag_stats: Some(flag), exon_depths: Some(Vec::new()), missing: Vec::new() };
    let f = prep::bam_metrics_features(&data, Some(&[]), SampleType::Tumor, &ThresholdRegistry::create_without_thresholds()).unwrap();
    let expected = SummaryTableFeature::ALL.iter().filter(|f| f.source_tool() == SourceTool::BamMetrics).count();
    assert_eq!(f.len(), expected);
}

// ---- feature.category.MissedGeneVariantPrepTest ----

fn driver_gene(name: &str, somatic: bool, germline: bool) -> DriverGene {
    use hmf_common::driver_panel::GermlineReporting;
    let g = if germline { GermlineReporting::Any } else { GermlineReporting::None };
    DriverGene {
        gene: name.into(),
        report_missense_and_inframe: somatic,
        report_nonsense_and_frameshift: somatic,
        report_splice: somatic,
        report_deletion: false,
        report_het_deletion: false,
        report_loh: false,
        het_deletion_threshold: 0.0,
        report_disruption: false,
        report_amplification: false,
        amplification_ratio: 0.0,
        report_somatic_hotspot: somatic,
        likelihood_type: hmf_common::driver_panel::DriverCategory::Tsg,
        report_germline_variant: g,
        report_germline_hotspot: g,
        report_germline_disruption: g,
        report_germline_deletion: g,
        report_germline_amplification: false,
        additional_reported_transcripts: Vec::new(),
        report_pgx: false,
        report_high_expression: false,
        report_low_expression: false,
        report_novel_splice_junction: false,
    }
}

/// MissedGeneVariantPrepTest's EXON_DEPTHS (upstream bdceb2302a): GERMLINE_GENE's two exons sum to the mixed counts
fn exon_depths() -> Vec<ExonDepth> {
    let mixed = vec![(dr(0, 0), 0), (dr(15, 15), 100), (dr(20, 20), 100), (dr(40, 49), 800)];
    let part1 = vec![(dr(0, 0), 0), (dr(15, 15), 100), (dr(20, 20), 100)];
    let part2 = vec![(dr(40, 49), 800)];
    let none = vec![(dr(0, 0), 1)];
    let e = |n: &str, c: &str, rank: i32, counts: Vec<(DepthRange, i32)>| ExonDepth {
        depth: GeneDepth { gene: n.into(), chromosome: c.into(), pos_start: 100, pos_end: 200, depth_counts: counts },
        exon_rank: rank,
        median_depth: 0.0,
        mean_depth: 0.0,
    };
    vec![
        e("GERMLINE_GENE", "1", 1, part1),
        e("GERMLINE_GENE", "1", 2, part2),
        e("GERMLINE_CHR_X_GENE", "X", 1, mixed.clone()),
        e("SOMATIC_GENE", "1", 1, mixed.clone()),
        e("SOMATIC_GENE", "1", 2, none.clone()), // excluded
        e("SOMATIC_CHR_X_GENE", "X", 1, mixed),
        e("SOMATIC_NO_COVERAGE_GENE", "1", 1, none.clone()),
        e("SOMATIC_EXCLUDED_GENE", "1", 1, none), // every exon excluded
    ]
}

fn missed_gene_exclusions() -> ExonExclusions {
    exclusions(&[("SOMATIC_GENE", 2), ("SOMATIC_EXCLUDED_GENE", 1)])
}

fn panel_config() -> crate::config::QseeConfig {
    let mut c = crate::tests::port::test_config();
    c.driver_genes = vec![
        driver_gene("GERMLINE_GENE", false, true),
        driver_gene("GERMLINE_CHR_X_GENE", false, true),
        driver_gene("SOMATIC_GENE", true, false),
        driver_gene("SOMATIC_CHR_X_GENE", true, false),
        driver_gene("SOMATIC_NO_COVERAGE_GENE", true, false),
        driver_gene("SOMATIC_EXCLUDED_GENE", true, false),
    ];
    c
}

#[test]
fn missed_gene_variant_germline_likelihoods() {
    let c = panel_config();
    let genes = prep::reportable_genes(&c, SampleType::Normal);
    let f = prep::missed_variant_likelihoods(&exon_depths(), &genes, &missed_gene_exclusions(), 1.0, "FEMALE", DEFAULT_LOD_MIN_ALT_READS).unwrap();
    let expected = [("Gene:GERMLINE_GENE", 0.000389385), ("Gene:GERMLINE_CHR_X_GENE", 0.000389385)];
    assert_eq!(f.len(), expected.len());
    for (a, (n, v)) in f.iter().zip(expected) {
        assert_eq!(a.key.name, n);
        assert!(close(a.value, v, 1e-9), "{} {}", a.value, v);
    }
}

#[test]
fn missed_gene_variant_somatic_likelihoods() {
    let c = panel_config();
    let genes = prep::reportable_genes(&c, SampleType::Tumor);
    let f = prep::missed_variant_likelihoods(&exon_depths(), &genes, &missed_gene_exclusions(), 0.5, "MALE", DEFAULT_LOD_MIN_ALT_READS).unwrap();
    let expected = [("Gene:SOMATIC_GENE", 0.032975381), ("Gene:SOMATIC_CHR_X_GENE", 0.000389385), ("Gene:SOMATIC_NO_COVERAGE_GENE", 1.0)];
    assert_eq!(f.len(), expected.len());
    for (a, (n, v)) in f.iter().zip(expected) {
        assert_eq!(a.key.name, n);
        assert!(close(a.value, v, 1e-9), "{} {}", a.value, v);
    }
}

// ---- feature.category.MsIndelErrorPrepTest ----

fn jitter_row(ref_units: i32, unit: &str, ct: ConsensusType, del: i32, ins: i32, total: i32) -> JitterTableRow {
    let mut r = JitterTableRow::new(ref_units, unit, ct);
    r.set_jitter_read_count(DELETIONS_KEY, del);
    r.set_jitter_read_count(INSERTIONS_KEY, ins);
    r.total_read_count = total;
    r
}

fn sorted_counts(r: &JitterTableRow) -> Vec<(i32, i32)> {
    let mut c = r.jitter_counts.clone();
    c.sort();
    c
}

#[test]
fn ms_indel_aggregate_by_jitter() {
    let mut r = JitterTableRow::new(0, "", ConsensusType::None);
    let (jitters, counts) = ([-3, -2, -1, 0, 1, 2, 3], [1, 2, 5, 0, 6, 3, 1]);
    r.total_read_count = counts.iter().sum();
    for (j, c) in jitters.iter().zip(counts) {
        r.set_jitter_read_count(*j, c);
    }
    let out = &prep::aggregate_by_jitter(&[r])[0];
    assert_eq!(out.jitter_read_count(DELETIONS_KEY), 8);
    assert_eq!(out.jitter_read_count(INSERTIONS_KEY), 10);
}

#[test]
fn ms_indel_aggregate_by_repeat_type() {
    use ConsensusType::{None as N, Single as S};
    let table = vec![
        jitter_row(0, "A/T", S, 1, 1, 2),
        jitter_row(0, "A/T", S, 1, 1, 2),
        jitter_row(1, "A/T", S, 2, 2, 4),
        jitter_row(1, "A/T", S, 2, 2, 4),
        jitter_row(1, "A/T", N, 3, 3, 6),
        jitter_row(1, "A/T", N, 3, 3, 6),
        jitter_row(1, "C/G", N, 4, 4, 8),
        jitter_row(1, "C/G", N, 4, 4, 8),
        jitter_row(1, "NN/NN", N, 5, 5, 10),
        jitter_row(1, "MM/MM/MM/MM", N, 5, 5, 10),
        jitter_row(1, "1bp repeat", N, 6, 6, 12),
        jitter_row(1, "2bp repeat", N, 6, 6, 12),
        jitter_row(1, "3bp repeat", N, 6, 6, 12),
    ];
    let expected = vec![
        jitter_row(0, "A/T repeat", S, 2, 2, 4),
        jitter_row(1, "A/T repeat", S, 4, 4, 8),
        jitter_row(1, "A/T repeat", N, 6, 6, 12),
        jitter_row(1, "C/G repeat", N, 8, 8, 16),
        jitter_row(1, "2bp repeat", N, 10, 10, 20),
        jitter_row(1, ">=3bp repeat", N, 18, 18, 36),
    ];
    let out = prep::aggregate_by_repeat_type(&table).unwrap();
    assert_eq!(out.len(), expected.len());
    for (a, e) in out.iter().zip(&expected) {
        assert_eq!((a.ref_num_units, &a.repeat_unit, a.consensus_type, a.total_read_count), (e.ref_num_units, &e.repeat_unit, e.consensus_type, e.total_read_count));
        assert_eq!(sorted_counts(a), sorted_counts(e));
    }
}

#[test]
fn ms_indel_phred_scores() {
    let f = prep::ms_indel_phred_scores(&[jitter_row(4, "A/T repeat", ConsensusType::None, 6, 4, 100)]);
    assert_eq!(f[0].key.feature_type, Some(FeatureType::MsIndelErrorRates));
    assert!(close(f[0].value, 10.0, 0.01));
    assert_eq!(f[1].key.feature_type, Some(FeatureType::MsIndelErrorBias));
    assert!(close(f[1].value, -1.76, 0.01));
}

// ---- feature.category.BaseQualRecalibrationPrepTest ----

fn bqr(rt: ConsensusType, r: char, a: char, ctx: &str, count: i32, q: i32, recal: f64) -> prep::StandardBqrRecord {
    prep::StandardBqrRecord::from(&BqrRecord { ref_base: r as u8, alt: a as u8, trinucleotide_context: ctx.as_bytes().to_vec(), quality: q as i8, read_type: rt, count, recalibrated_quality: recal })
}

fn bqr_records() -> Vec<prep::StandardBqrRecord> {
    use ConsensusType::{Dual as D, None as N};
    vec![
        bqr(D, 'C', 'G', "CCC", 90, 11, 13.0),
        bqr(D, 'C', 'G', "CCC", 90, 17, 19.0),
        bqr(D, 'G', 'C', "GGG", 10, 11, 10.5),
        bqr(D, 'G', 'C', "GGG", 10, 17, 16.5),
        bqr(N, 'C', 'A', "CCC", 70, 37, 39.0),
        bqr(N, 'C', 'A', "CCC", 70, 42, 44.0),
        bqr(N, 'G', 'T', "GGG", 30, 37, 36.5),
        bqr(N, 'G', 'T', "GGG", 30, 42, 41.5),
        bqr(N, 'C', 'A', "ACA", 90, 37, 41.0),
        bqr(N, 'C', 'A', "ACA", 90, 42, 46.0),
        bqr(N, 'G', 'T', "TGT", 10, 37, 36.5),
        bqr(N, 'G', 'T', "TGT", 10, 42, 41.5),
        // T>A_ATA: a count below the minimum
        bqr(N, 'T', 'A', "ATA", 2, 37, 39.0),
    ]
}

// ---- feature.category.bqr.StandardTrinucleotideContextTest (upstream 883ab11bba) ----

#[test]
fn standard_context_t_to_g() {
    use prep::StandardTrinucleotideContext as S;
    assert_eq!(S::from_nucleotides('T', 'G', "GTA").mutation_context_string().unwrap(), "G[T>G]A +");
    assert_eq!(S::from_nucleotides('A', 'C', "TAC").mutation_context_string().unwrap(), "G[T>G]A -");
}

#[test]
fn standard_context_c_to_t() {
    use prep::StandardTrinucleotideContext as S;
    assert_eq!(S::from_nucleotides('C', 'T', "GCA").mutation_context_string().unwrap(), "G[C>T]A +");
    assert_eq!(S::from_nucleotides('G', 'A', "TGC").mutation_context_string().unwrap(), "G[C>T]A -");
}

#[test]
fn standard_context_mutation_strings() {
    use prep::StandardTrinucleotideContext as S;
    let forward_cpg = S::from_nucleotides('C', 'T', "ACG");
    let reverse_cpg = S::from_nucleotides('G', 'A', "CGT"); // the reverse of ACG
    let non_cpg = S::from_nucleotides('C', 'T', "ACA");
    let non_c_to_t = S::from_nucleotides('C', 'A', "ACG");
    assert_eq!(forward_cpg.mutation_string(), "C>T");
    assert_eq!(non_c_to_t.mutation_string(), "C>A");
    assert_eq!(forward_cpg.mutation_string_split_cpg().unwrap(), "C>T CpG");
    assert_eq!(reverse_cpg.mutation_string_split_cpg().unwrap(), "C>T CpG");
    assert_eq!(non_cpg.mutation_string_split_cpg().unwrap(), "C>T other");
    assert_eq!(non_c_to_t.mutation_string_split_cpg().unwrap(), "C>A");
}

#[test]
fn bqr_bin_ranges_for_illumina() {
    let b = BaseQualBinner::new(SequencingType::Illumina);
    assert_eq!((b.range(BaseQualBin::Low).unwrap().lower, b.range(BaseQualBin::Low).unwrap().upper), (0, 29));
    assert!(b.range(BaseQualBin::Medium).is_none());
    assert_eq!((b.range(BaseQualBin::High).unwrap().lower, b.range(BaseQualBin::High).unwrap().upper), (30, i8::MAX));
}

#[test]
fn bqr_change_in_qual_per_trinuc_context() {
    let f = prep::bqr_per_trinuc_context(&bqr_records(), &BaseQualBinner::new(SequencingType::Illumina)).unwrap();
    let expected = [
        ("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:CCC;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", 2.0),
        ("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:CCC;ReadOrientation:REVERSE;OriginalQualBin:HIGH (30+)", -0.5),
        ("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ACA;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", 4.0),
        ("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ACA;ReadOrientation:REVERSE;OriginalQualBin:HIGH (30+)", -0.5),
        ("ReadType:NONE;StandardMutation:T>A;StandardTrinucContext:ATA;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", f64::NAN),
    ];
    assert_eq!(f.len(), expected.len());
    for (a, (n, v)) in f.iter().zip(expected) {
        assert_eq!(a.key.name, n);
        assert!(if v.is_nan() { a.value.is_nan() } else { close(a.value, v, 0.001) }, "{n}: {}", a.value);
    }
}

#[test]
fn bqr_change_in_qual_per_original_qual() {
    let f = prep::bqr_per_original_qual(&bqr_records(), &BaseQualBinner::new(SequencingType::Illumina)).unwrap();
    let expected = [
        ("ReadType:DUAL;StandardMutation:C>G;StandardTrinucContext:ANY;ReadOrientation:FORWARD;OriginalQualBin:LOW (0-29)", 2.0),
        ("ReadType:DUAL;StandardMutation:C>G;StandardTrinucContext:ANY;ReadOrientation:REVERSE;OriginalQualBin:LOW (0-29)", -0.5),
        ("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", 3.118),
        ("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;ReadOrientation:REVERSE;OriginalQualBin:HIGH (30+)", -0.5),
        ("ReadType:NONE;StandardMutation:T>A;StandardTrinucContext:ANY;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", f64::NAN),
    ];
    assert_eq!(f.len(), expected.len());
    for (a, (n, v)) in f.iter().zip(expected) {
        assert_eq!(a.key.name, n);
        assert!(if v.is_nan() { a.value.is_nan() } else { close(a.value, v, 0.001) }, "{n}: {}", a.value);
    }
}

// ---- feature.category.DiscordantFragFreqPrepTest ----

#[test]
fn discordant_proportions() {
    let stats = hmf_common::esvee_files::EsveeDiscordantStats { total_reads: 100, prep_reads: i32::MIN as i64, type_counts: [1; 11] };
    let f = prep::discordant_features(&stats);
    let per_group = |g: DiscordantFragGroup| (0..11).filter(|i| DiscordantFragGroup::from_type_index(*i) == g).count() as f64;
    assert_eq!(f.len(), DiscordantFragGroup::ALL.len());
    for (a, g) in f.iter().zip(DiscordantFragGroup::ALL) {
        assert_eq!(a.key.name, format!("DiscordantFragType:{}", g.name()));
        assert!(close(a.value, per_group(g) * 0.01, 0.001));
    }
}

// ---- feature.category.GcBiasPrepTest ----

#[test]
fn gc_bias_normalises_median_read_depths() {
    let g = hmf_common::cobalt_files::GcMedianReadDepth { mean: 100.0, median: 100.0, median_per_bucket: vec![(26, 105.0), (25, 100.0)] };
    let f = prep::gc_bias_features(&g);
    assert_eq!(f[0].key.name, "GCBucket:25");
    assert!(close(f[0].value, 1.0, 0.001));
    assert_eq!(f[1].key.name, "GCBucket:26");
    assert!(close(f[1].value, 1.05, 0.001));
}
