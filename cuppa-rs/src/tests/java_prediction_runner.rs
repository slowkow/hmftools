//! PredictionRunnerTest (ignored in Java, as it needs Python and a classifier): the jar's main class with the
//! embedded classifier's exported file, on the test resources' tumour sample and its RNA, all categories.

use super::{temp_dir, testdata};

#[test]
fn can_predict_from_pipeline_output() {
    let out = temp_dir("runner");
    let args: Vec<String> = [
        "-sample",
        "TUMOR_SAMPLE",
        "-rna_sample",
        "TUMOR_SAMPLE_RNA",
        "-categories",
        "ALL",
        "-sample_data_dir",
        &format!("{}*", testdata("java/pipeline_output/")),
        "-ref_alt_sj_sites",
        &testdata("java/alt_sj.selected_loci.minimal.tsv"),
        "-classifier_path",
        &format!("{}/resources/cuppa_model.38.bin.gz", env!("CARGO_MANIFEST_DIR")),
        "-output_dir",
        &out,
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(crate::prep::runner_main(&args), 0);
    let summary = std::fs::read_to_string(format!("{out}TUMOR_SAMPLE.cuppa.pred_summ.tsv")).unwrap();
    let rows: Vec<&str> = summary.lines().collect();
    // combined, the 4 DNA rows and the 3 RNA rows
    assert_eq!(rows.len(), 9);
    assert!(rows[1].starts_with("TUMOR_SAMPLE\tcombined\tcombined\t"));
    assert!(std::path::Path::new(&format!("{out}TUMOR_SAMPLE.cuppa.vis_data.tsv")).exists());
    assert!(std::path::Path::new(&format!("{out}TUMOR_SAMPLE.cuppa_data.tsv.gz")).exists());
    let _ = std::fs::remove_dir_all(&out);
}

#[test]
fn can_predict_from_input_features() {
    // the features of the public HG008 sample, as the pipeline wrote them
    let out = temp_dir("runner-features");
    let args: Vec<String> = [
        "-sample",
        "HG008-T",
        "-features_path",
        &testdata("hg008/HG008-T.cuppa_data.tsv.gz"),
        "-classifier_path",
        &format!("{}/resources/cuppa_model.38.bin.gz", env!("CARGO_MANIFEST_DIR")),
        "-output_dir",
        &out,
        "-clf_group",
        "DNA",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    assert_eq!(crate::prep::runner_main(&args), 0);
    let got = std::fs::read_to_string(format!("{out}HG008-T.cuppa.pred_summ.tsv")).unwrap();
    assert_eq!(got, std::fs::read_to_string(testdata("hg008/HG008-T.cuppa.pred_summ.tsv")).unwrap());
    let _ = std::fs::remove_dir_all(&out);
}
