//! Tests of the port: number formats, the threshold display, the key order, the quality bins of each
//! sequencing type, input parsing details, and each corrected finding in both modes.

use crate::cohort::{CohortPercentiles, FeaturePercentiles};
use crate::common::{decimal_4, decimal_8, java_binary_search, parse_multi_field_as_array, parse_multi_field_as_pairs, SampleIds, SampleType};
use crate::compat::with_java;
use crate::config::{convert_wildcard_sample_path, parse_categories, PrepCategory, QseeConfig, SequencingType};
use crate::feature::{BaseQualBin, FeatureKey, FeatureType, NumberFormat, RepeatType, SourceTool};
use crate::prep::{self, is_diploid, BaseQualBinner};
use crate::status::{try_form_display_string, ComparisonOperator, QcStatusType, ThresholdKey, ThresholdRegistry};

/// a configuration with Java's defaults and no inputs
pub fn test_config() -> QseeConfig {
    QseeConfig {
        samples: SampleIds { tumor_ids: vec!["T".into()], reference_ids: vec!["N".into()] },
        sample_data_dir: String::new(),
        bam_metrics_dir_tumor: String::new(),
        bam_metrics_dir_ref: String::new(),
        redux_dir_tumor: String::new(),
        redux_dir_ref: String::new(),
        cobalt_dir: String::new(),
        esvee_dir: String::new(),
        purple_dir: String::new(),
        driver_genes: Vec::new(),
        cohort_percentiles_file: None,
        thresholds_file: None,
        thresholds: ThresholdRegistry::create_default(false),
        sequencing_type: SequencingType::Illumina,
        categories: PrepCategory::ALL.to_vec(),
        targeted_mode: false,
        allow_missing_input: false,
        skip_plot: true,
        write_lod_details: false,
        output_dir: String::new(),
        output_id: None,
        threads: 1,
    }
}

#[test]
fn decimal_formats() {
    // DecimalFormat("0.########") with "Inf"; DecimalFormat("#.####")
    assert_eq!(decimal_8(0.96615902345), "0.96615902");
    assert_eq!(decimal_8(71.004), "71.004");
    assert_eq!(decimal_8(0.0), "0");
    assert_eq!(decimal_8(-0.0), "-0");
    assert_eq!(decimal_8(f64::INFINITY), "Inf");
    assert_eq!(decimal_8(f64::NEG_INFINITY), "-Inf");
    assert_eq!(decimal_8(f64::NAN), "NaN");
    assert_eq!(decimal_8(1e-9), "0");
    assert_eq!(decimal_4(0.96615902), "0.9662");
    assert_eq!(decimal_4(f64::INFINITY), "∞");
}

#[test]
fn threshold_display_as_big_decimal() {
    use ComparisonOperator::*;
    let f = |op, v: f64, nf| try_form_display_string(QcStatusType::Warn, Some(op), v, nf).unwrap();
    assert_eq!(f(LessThan, 70.0, NumberFormat::Number), "<70");
    assert_eq!(f(LessThan, 0.95, NumberFormat::Percent), "<95%");
    assert_eq!(f(GreaterThanOrEqual, 0.25, NumberFormat::Percent), ">=25%");
    assert_eq!(f(GreaterThan, 280.0, NumberFormat::Log10), ">280");
    assert_eq!(f(GreaterThan, 1e-7, NumberFormat::Number), ">0.0000001");
    assert_eq!(f(GreaterThan, 1e-7, NumberFormat::Percent), ">0.00001%");
    assert_eq!(f(LessThanOrEqual, 29.501, NumberFormat::Number), "<=29.501");
    assert_eq!(f(LessThan, 1e21, NumberFormat::Number), "<1000000000000000000000");
    assert_eq!(f(LessThan, -0.5, NumberFormat::Percent), "<-50%");
    assert_eq!(f(LessThan, 0.0, NumberFormat::Percent), "<0%");
    assert_eq!(try_form_display_string(QcStatusType::None, Some(LessThan), 1.0, NumberFormat::Number).unwrap(), "");
    // QS6: an infinite threshold
    with_java(true, || assert!(try_form_display_string(QcStatusType::Fail, Some(LessThan), f64::NEG_INFINITY, NumberFormat::Percent).is_err()));
    with_java(false, || assert_eq!(try_form_display_string(QcStatusType::Fail, Some(LessThan), f64::NEG_INFINITY, NumberFormat::Percent).unwrap(), "<-Inf%"));
}

#[test]
fn default_thresholds() {
    let t = ThresholdRegistry::create_default(false);
    // 5 common (both types), 9 decided by PURPLE, 7 tumor and 7 normal coverage thresholds
    assert_eq!(t.get_all().len(), 10 + 9 + 14);
    let tg = ThresholdRegistry::create_default(true);
    let key = ThresholdKey::new(SampleType::Tumor, FeatureType::SummaryTable, "DUPLICATE_READS", QcStatusType::Warn);
    assert_eq!((t.get(&key).unwrap().threshold, tg.get(&key).unwrap().threshold), (0.3, 0.99));
    let key = ThresholdKey::new(SampleType::Normal, FeatureType::SummaryTable, "DUAL_STRAND_READS", QcStatusType::Warn);
    assert!(tg.get(&key).unwrap().threshold.is_nan());
    // a NaN threshold never fails
    let key = ThresholdKey::new(SampleType::Normal, FeatureType::SummaryTable, "COVERAGE_ABOVE_60", QcStatusType::Warn);
    assert_eq!(t.get(&key).unwrap().sample_qc_status(0.0).status_type, QcStatusType::Pass);
}

#[test]
fn contamination_operators_qs4() {
    let key = ThresholdKey::new(SampleType::Tumor, FeatureType::SummaryTable, "CONTAMINATION", QcStatusType::Warn);
    with_java(true, || assert_eq!(ThresholdRegistry::create_default(false).get(&key).unwrap().operator, Some(ComparisonOperator::GreaterThan)));
    with_java(false, || assert_eq!(ThresholdRegistry::create_default(false).get(&key).unwrap().operator, Some(ComparisonOperator::GreaterThanOrEqual)));
}

#[test]
fn feature_key_order() {
    let k = |n: &str, t, s| FeatureKey::new(n, t, s);
    let mut keys = vec![
        k("Gene:TP53", FeatureType::MissedVariantLikelihood, SourceTool::BamMetrics),
        k("GCBucket:30", FeatureType::GcBias, SourceTool::Cobalt),
        k("GCBucket:4", FeatureType::GcBias, SourceTool::Cobalt),
        k("Gene:BRCA1", FeatureType::MissedVariantLikelihood, SourceTool::BamMetrics),
        k("ReadDepth:10", FeatureType::CoverageDistribution, SourceTool::BamMetrics),
        k("ReadDepth:9", FeatureType::CoverageDistribution, SourceTool::BamMetrics),
        k("MEAN_COVERAGE", FeatureType::SummaryTable, SourceTool::BamMetrics),
        FeatureKey { name: "x".into(), feature_type: None, source_tool: None },
        k("ConsensusType:DUAL;RepeatUnitType:A/T repeat;RefNumUnits:4", FeatureType::MsIndelErrorBias, SourceTool::Redux),
        k("ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:10", FeatureType::MsIndelErrorBias, SourceTool::Redux),
        k("ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:9", FeatureType::MsIndelErrorBias, SourceTool::Redux),
    ];
    keys.sort_by(|a, b| a.compare(b));
    let names: Vec<&str> = keys.iter().map(|k| k.name.as_str()).collect();
    assert_eq!(
        names,
        vec![
            "MEAN_COVERAGE",
            "ReadDepth:9",
            "ReadDepth:10",
            "GCBucket:4",
            "GCBucket:30",
            "Gene:BRCA1",
            "Gene:TP53",
            "ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:9",
            "ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:10",
            "ConsensusType:DUAL;RepeatUnitType:A/T repeat;RefNumUnits:4",
            "x",
        ]
    );
    // equal by name and type, whatever the source tool
    assert_eq!(k("a", FeatureType::GcBias, SourceTool::Cobalt), k("a", FeatureType::GcBias, SourceTool::Redux));
}

#[test]
fn multi_field_parsing() {
    assert_eq!(parse_multi_field_as_array("A:1;B:2"), vec!["A", "1", "B", "2"]);
    assert_eq!(parse_multi_field_as_array("A:"), vec!["A"]);
    assert_eq!(parse_multi_field_as_pairs("A:1;B:x y").unwrap(), vec![("A".to_string(), "1".to_string()), ("B".to_string(), "x y".to_string())]);
    assert!(parse_multi_field_as_pairs("A").is_err());
}

#[test]
fn binary_search_like_java() {
    let a = [-1.0, 0.0, 2.5, 7.0];
    assert_eq!(java_binary_search(&a, 2.5), 2);
    assert_eq!(java_binary_search(&a, 3.0), -4);
    assert_eq!(java_binary_search(&a, -0.0), -2);
    assert_eq!(java_binary_search(&a, -5.0), -1);
    assert_eq!(java_binary_search(&a, 9.0), -5);
}

#[test]
fn quality_bins_per_sequencing_type() {
    let r = |t, bin| BaseQualBinner::new(t).range(bin).map(|r| r.display());
    assert_eq!(r(SequencingType::Illumina, BaseQualBin::Low).as_deref(), Some("LOW (0-29)"));
    assert_eq!(r(SequencingType::Illumina, BaseQualBin::High).as_deref(), Some("HIGH (30+)"));
    assert_eq!(r(SequencingType::Sbx, BaseQualBin::Low).as_deref(), Some("LOW (0-26)"));
    assert_eq!(r(SequencingType::Sbx, BaseQualBin::Medium).as_deref(), Some("MEDIUM (27-29)"));
    assert_eq!(r(SequencingType::Sbx, BaseQualBin::High).as_deref(), Some("HIGH (30+)"));
    assert_eq!(r(SequencingType::Ultima, BaseQualBin::Low).as_deref(), Some("LOW (0-34)"));
    assert_eq!(r(SequencingType::Ultima, BaseQualBin::High).as_deref(), Some("HIGH (35+)"));
    assert!(BaseQualBinner::new(SequencingType::Illumina).bin_from(-56).is_err());
    assert_eq!(BaseQualBinner::new(SequencingType::Illumina).bin_from(127).unwrap().bin, BaseQualBin::High);
}

#[test]
fn repeat_units() {
    assert_eq!(RepeatType::from_repeat_unit("A/T").unwrap(), RepeatType::RepeatAT);
    assert_eq!(RepeatType::from_repeat_unit("AC/CA/GT/TG").unwrap(), RepeatType::Repeat2Bp);
    assert_eq!(RepeatType::from_repeat_unit("A_/x").unwrap(), RepeatType::Repeat2Bp);
    assert_eq!(RepeatType::from_repeat_unit("12bp repeat").unwrap(), RepeatType::Repeat3BpOrLonger);
    assert!(RepeatType::from_repeat_unit("A/G").is_err());
    assert!(RepeatType::from_repeat_unit("bp repeat").is_err());
    assert!(RepeatType::from_repeat_unit("3bp repeats").is_err());
    assert!(RepeatType::from_repeat_unit("AC/\n").is_err());
}

#[test]
fn chromosomes() {
    assert!(is_diploid("chr1", "MALE").unwrap());
    assert!(is_diploid("22", "MALE").unwrap());
    assert!(is_diploid("chrX", "FEMALE").unwrap());
    assert!(!is_diploid("chrX", "MALE").unwrap());
    assert!(!is_diploid("CHRY", "FEMALE").unwrap());
    assert!(is_diploid("chrM", "FEMALE").is_err());
    assert!(is_diploid("chr01", "FEMALE").is_err());
    assert!(is_diploid("chrx", "FEMALE").is_err());
    assert!(is_diploid("chr23", "FEMALE").is_err());
}

#[test]
fn categories() {
    assert_eq!(parse_categories(None).unwrap().len(), 10);
    let c = parse_categories(Some("gc_Bias,MS_INDEL_ERROR")).unwrap();
    assert!(!c.contains(&PrepCategory::GcBias) && !c.contains(&PrepCategory::MsIndelError) && c.len() == 8);
    assert!(parse_categories(Some("SUMMARY_TABLE_PURPLE")).is_err());
    assert!(parse_categories(Some("NOPE")).is_err());
}

#[test]
fn missing_cohort_sample_type_qs3() {
    let mut c = CohortPercentiles::default();
    let key = FeatureKey::new("GCBucket:30", FeatureType::GcBias, SourceTool::Cobalt);
    c.add(FeaturePercentiles { sample_type: SampleType::Tumor, key: key.clone(), percentiles: vec![0.0, 100.0], ref_values: vec![0.0, 1.0] });
    with_java(true, || {
        assert!(c.warn_if_missing(SampleType::Tumor, &key).is_ok());
        assert!(c.warn_if_missing(SampleType::Normal, &key).is_err());
    });
    with_java(false, || assert!(c.warn_if_missing(SampleType::Normal, &key).is_ok()));
}

#[test]
fn tool_directories_qs7() {
    let c = test_config();
    with_java(true, || {
        assert_eq!(c.purple_dir("T").unwrap(), "");
        assert_eq!(hmf_common::purple_files::purple_purity_filename(&c.purple_dir("T").unwrap(), "T"), "/T.purple.purity.tsv");
    });
    with_java(false, || assert_eq!(hmf_common::purple_files::purple_purity_filename(&c.purple_dir("T").unwrap(), "T"), "./T.purple.purity.tsv"));
    // BamMetrics and the duplicate frequencies add a separator only where one is missing, in Java too
    assert_eq!(c.bam_metrics_dir("T", SampleType::Tumor).unwrap(), "");
    assert_eq!(hmf_common::redux_files::duplicate_freq_filename(&c.redux_dir_unjoined("T", SampleType::Tumor).unwrap(), "T"), "T.redux.duplicate_freq.tsv");
}

#[test]
fn wildcards_qs11() {
    assert_eq!(convert_wildcard_sample_path("dir/*/x", "S1").unwrap(), "dir/S1/x");
    assert_eq!(convert_wildcard_sample_path("dir/", "S$1").unwrap(), "dir/");
    with_java(true, || {
        assert!(convert_wildcard_sample_path("dir/*/", "S$1").is_err());
        assert_eq!(convert_wildcard_sample_path("dir/*/", "S$0").unwrap(), "dir/S*/");
        assert_eq!(convert_wildcard_sample_path("dir/*/", "S\\x").unwrap(), "dir/Sx/");
    });
    with_java(false, || assert_eq!(convert_wildcard_sample_path("dir/*/", "S$1").unwrap(), "dir/S$1/"));
}

#[test]
fn contamination_warning_qs5() {
    let mut ctx = crate::tests::java_tests_support::purity_context(&["WARN_CONTAMINATION"]);
    ctx.qc.contamination = 0.05;
    let t = ThresholdRegistry::create_default(false);
    let status = |java: bool| with_java(java, || prep::purple_features(&ctx, &t).unwrap().into_iter().find(|f| f.key.name == "CONTAMINATION").unwrap().metadata.qc_status);
    let j = status(true);
    assert_eq!((j.status_type, j.threshold), (QcStatusType::Pass, 0.1));
    let d = status(false);
    assert_eq!((d.status_type, d.threshold), (QcStatusType::Warn, 0.02));
}

#[test]
fn purple_threshold_overrides_qs9() {
    let dir = std::env::temp_dir().join(format!("qsee_qs9_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("o.tsv");
    std::fs::write(&p, "SampleType\tFeatureType\tFeatureName\tQcStatusType\tComparisonOperator\tThreshold\nTUMOR\tSUMMARY_TABLE\tPURITY\tWARN\tLESS_THAN\t0.3\n").unwrap();
    let key = ThresholdKey::new(SampleType::Tumor, FeatureType::SummaryTable, "PURITY", QcStatusType::Warn);
    let java = with_java(true, || crate::status::read_threshold_overrides(p.to_str().unwrap(), false).unwrap());
    assert_eq!(java.get(&key).unwrap().threshold, 0.3);
    let fixed = with_java(false, || crate::status::read_threshold_overrides(p.to_str().unwrap(), false).unwrap());
    assert_eq!(fixed.get(&key).unwrap().threshold, 0.2);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn default_thresholds_writer_qs2() {
    let dir = std::env::temp_dir().join(format!("qsee_qs2_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let d = dir.to_str().unwrap();
    crate::status::write_default_thresholds(d, true, false).unwrap();
    let text = std::fs::read_to_string(dir.join("qsee.thresholds.default.tsv")).unwrap();
    assert!(text.starts_with("SampleType\tFeatureType\tFeatureName\tQcStatusType\tComparisonOperator\tThreshold\tOverridable\n"));
    assert!(text.contains("TUMOR\tSUMMARY_TABLE\tDUPLICATE_READS\tWARN\tGREATER_THAN\t0.3\ttrue\n"));
    assert!(text.contains("TUMOR\tSUMMARY_TABLE\tCOVERAGE_ABOVE_250\tWARN\tLESS_THAN\tNaN\ttrue\n"));
    assert!(!text.contains("PURITY"));
    // the tumor rows first
    let first_normal = text.find("\nNORMAL").unwrap();
    assert!(!text[first_normal..].contains("\nTUMOR"));
    crate::status::write_default_thresholds(d, false, true).unwrap();
    let text = std::fs::read_to_string(dir.join("qsee.thresholds.default.tsv")).unwrap();
    assert!(text.contains("TUMOR\tSUMMARY_TABLE\tPURITY\tWARN\tNA\t0.2\tfalse\n"));
    assert!(text.contains("TUMOR\tSUMMARY_TABLE\tDUPLICATE_READS\tWARN\tGREATER_THAN\t0.99\ttrue\n"));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn bqr_missing_file_qs1() {
    let mut c = test_config();
    c.redux_dir_tumor = "/nonexistent-qsee-dir".into();
    let r = with_java(true, || prep::extract_sample_data(&c, PrepCategory::BaseQualRecalibration, "T", SampleType::Tumor));
    assert!(matches!(r, Err(prep::PrepError::Fatal(_))));
    let r = with_java(false, || prep::extract_sample_data(&c, PrepCategory::BaseQualRecalibration, "T", SampleType::Tumor));
    assert!(matches!(r, Err(prep::PrepError::Missing(_))));
}

#[test]
fn weighted_mean_like_commons_math() {
    assert_eq!(prep::weighted_mean(&[1.0, 3.0], &[1.0, 3.0]).unwrap(), 2.5);
    assert!(prep::weighted_mean(&[1.0], &[-1.0]).is_err());
    assert!(prep::weighted_mean(&[1.0], &[0.0]).is_err());
}

#[test]
fn sample_ids_details() {
    // String.split(","): trailing empty ids dropped
    let s = SampleIds::from_strings(Some("A,B,"), None).unwrap();
    assert_eq!(s.tumor_ids, vec!["A", "B"]);
    assert!(SampleIds::from_lines(&["TumorId".to_string()]).is_err());
    assert!(SampleIds::from_lines(&["ReferenceId".to_string(), "R".to_string()]).is_err());
    assert!(s.is_single_patient() == false);
}

/// timing of the 8-place DecimalFormat (run with --ignored --nocapture)
#[test]
#[ignore]
fn decimal_8_timing() {
    let mut x = 0.123456789f64;
    let t = std::time::Instant::now();
    let mut n = 0usize;
    for i in 0..400_000 {
        x = (x * 1.000123 + i as f64 * 1e-7) % 3.0;
        n += decimal_8(x).len();
    }
    eprintln!("400000 values in {:?} ({n} chars)", t.elapsed());
}

#[test]
fn percentile_increment_qs12() {
    use crate::cohort::PercentileTransformer;
    with_java(true, || {
        assert!(PercentileTransformer::with_increment(0.1).is_err());
        assert_eq!(PercentileTransformer::with_increment(12.5).unwrap().percentiles.len(), 9);
    });
    with_java(false, || {
        let t = PercentileTransformer::with_increment(0.1).unwrap();
        assert_eq!(t.percentiles.len(), 1001);
        assert!(PercentileTransformer::with_increment(30.0).is_err());
        assert_eq!(PercentileTransformer::with_increment(0.2).unwrap().percentiles.len(), 501);
    });
}
