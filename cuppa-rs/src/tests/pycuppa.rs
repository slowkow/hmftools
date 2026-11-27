//! pycuppa's tests of the components cuppa-rs ports (tests/components, tests/sample_data,
//! tests/classifier/test_cuppa_prediction.py); those that train a classifier are not ported (training is
//! not part of the port).

use crate::predict::classifier::{combine, cos_sim, max_scale, non_best_scale, sex_filter};
use crate::predict::numpy;

fn round3(v: f64) -> f64 {
    numpy::round(v, 3)
}

// test_prob_combiner.py -----------------------------------------------------------------------------------------
#[test]
fn can_multiply_probs_without_zeroes() {
    // dna Breast 0.9, Lung 0.1; rna Breast 0.8, Lung 0.2
    let c = combine(&[vec![0.9, 0.1]], &[vec![0.8, 0.2]], 0.01);
    assert_eq!(round3(c[0][0]), 0.973);
}

#[test]
fn can_multiply_probs_with_zeroes() {
    let c = combine(&[vec![1.0, 0.0]], &[vec![0.8, 0.2]], 0.01);
    assert_eq!(round3(c[0][0]), 0.998);
}

// test_preprocessing.py ------------------------------------------------------------------------------------------
#[test]
fn log1p_transformer() {
    assert_eq!(round3((2.718f64 - 1.0).ln_1p()), 1.0);
}

#[test]
fn max_scaler_with_clipping() {
    // maxima of the training columns: 10 and 11
    assert_eq!(max_scale(100.0, 10.0), 1.0);
    assert_eq!(max_scale(-100.0, 11.0), 0.0);
    // 0 / 0 is 0, x / 0 the largest value clipped to 1
    assert_eq!(max_scale(0.0, 0.0), 0.0);
    assert_eq!(max_scale(5.0, 0.0), 1.0);
}

// test_profile_similarity.py -------------------------------------------------------------------------------------
fn norms(rows: &[&[f64]]) -> Vec<f64> {
    rows.iter().map(|r| numpy::pairwise_sum(&r.iter().map(|v| v * v).collect::<Vec<f64>>()).sqrt()).collect()
}

#[test]
fn array_vs_array_cos_sim_is_correct() {
    let counts = [1.0, 0.5, 0.25, 0.1];
    let profile: &[f64] = &[0.75, 0.65, 0.51, 0.4];
    let r = cos_sim(&[&counts], &[profile], &norms(&[profile]));
    assert_eq!(round3(r[0][0]), 0.911);
}

#[test]
fn matrix_vs_matrix_cos_sim_values_are_correct() {
    let counts: [&[f64]; 2] = [&[105.0, 8.0, 21.0], &[54.0, 294.0, 3.0]];
    let profiles: [&[f64]; 4] = [&[0.9, 0.05, 0.05], &[0.2, 0.6, 0.2], &[0.33, 0.33, 0.33], &[0.1, 0.3, 0.6]];
    let r = cos_sim(&counts, &profiles, &norms(&profiles));
    assert_eq!((r.len(), r[0].len()), (2, 4));
    assert_eq!(numpy::round(r[0][0], 2), 0.99);
    assert_eq!(numpy::round(r[1][1], 2), 0.95);
}

#[test]
fn non_best_scaling_gives_correct_cos_sim_values() {
    let mut row = vec![0.95, 0.92, 0.87, 0.86, 0.75];
    non_best_scale(&mut row, 5);
    let got: Vec<f64> = row.iter().map(|v| numpy::round(*v, 2)).collect();
    assert_eq!(got, vec![0.95, 0.79, 0.57, 0.54, 0.25]);
}

// test_prob_overriders.py: SexProbFilter --------------------------------------------------------------------------
fn kw(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn prostate_prob_for_male_sample_is_increased_after_override() {
    let classes = kw(&["Prostate", "Ovarian", "Skin"]);
    let mut rows = vec![vec![0.8, 0.2, f64::NAN]];
    sex_filter(&mut rows, &classes, &kw(&["prostate"]), &kw(&["ovarian"]), &[Some(true)]);
    assert!(rows[0][0] > 0.8);
    assert_eq!(rows[0][1], 0.0);
}

#[test]
fn ovarian_prob_for_female_sample_is_increased_after_override() {
    let classes = kw(&["Prostate", "Ovarian", "Skin"]);
    let mut rows = vec![vec![0.1, 0.9, f64::NAN]];
    sex_filter(&mut rows, &classes, &kw(&["prostate"]), &kw(&["ovarian"]), &[Some(false)]);
    assert!(rows[0][1] > 0.9);
    assert_eq!(rows[0][0], 0.0);
}

#[test]
fn unknown_sex_leaves_probs_unmodified_before_renormalising() {
    let classes = kw(&["Prostate", "Ovarian"]);
    let mut rows = vec![vec![0.25, 0.25]];
    sex_filter(&mut rows, &classes, &kw(&["prostate"]), &kw(&["ovarian"]), &[None]);
    assert_eq!(rows[0], vec![0.5, 0.5]);
}

// test_cuppa_prediction.py: CuppaPredSummaryBuilder on dummy predictions -------------------------------------------
#[test]
fn can_build_summary_from_dummy_predictions() {
    use crate::predict::output::{summarize, PredRow, Prediction};
    let p = Prediction {
        samples: vec!["sample_1".to_string()],
        classes: kw(&["cancer_type_1", "cancer_type_2"]),
        supertypes: kw(&["cancer_type_1", "cancer_type_2"]),
        rows: vec![
            PredRow { sample: 0, data_type: "prob", clf_group: Some("dna"), clf_name: Some("dna_combined"), feat_name: None, feat_value: Some(0.0), values: vec![0.0, f64::NAN] },
            PredRow {
                sample: 0,
                data_type: "feat_contrib",
                clf_group: Some("dna"),
                clf_name: Some("event"),
                feat_name: Some("event.tmb.snv_count".into()),
                feat_value: Some(0.0),
                values: vec![0.0, 0.0],
            },
            PredRow { sample: 0, data_type: "sig_quantile", clf_group: None, clf_name: None, feat_name: Some("SBS7".into()), feat_value: Some(0.0), values: vec![0.0, 0.0] },
        ],
    };
    let s = summarize(&p, false).unwrap();
    assert_eq!(s.len(), 1);
    assert_eq!(s[0].classes[0].as_deref(), Some("cancer_type_1"));
    // NaN probabilities rank last and give no class
    assert_eq!(s[0].classes[1], None);
}
