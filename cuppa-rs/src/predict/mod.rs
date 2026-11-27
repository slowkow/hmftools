//! pycuppa's prediction (`python -m cuppa.predict`, cuppa/runners/prediction_runner.py) in Rust: the
//! features file in, `<sample>.cuppa.pred_summ.tsv` and `<sample>.cuppa.vis_data.tsv` out, then the plot
//! `<sample>.cuppa.vis.png` (or `cuppa.vis.pdf` for several samples), which pycuppa draws with the R script
//! cuppa/visualization/plot_predictions.R (ggplot2, ggh4x, patchwork) and the port draws natively (`crate::plot`).

pub mod classifier;
pub mod cv;
pub mod features;
pub mod model;
pub mod numpy;
pub mod output;
pub mod pyfmt;
pub mod pylog;

use anyhow::{anyhow, bail, Result};

/// The options of cuppa.predict (RunnerArgParser.get_kwargs_predict)
#[derive(Clone, Debug, Default)]
pub struct PredictArgs {
    pub classifier_path: String,
    pub features_path: String,
    pub output_dir: String,
    pub sample_id: Option<String>,
    pub clf_group: String,
    pub cv_predictions_path: Option<String>,
    pub compress_tsv_files: bool,
    pub force_plot: bool,
    pub rename_classes: Option<Vec<(String, String)>>,
    pub log_to_file: bool,
    pub log_path: Option<String>,
    pub log_format: String,
}

const OPTIONS: [(&str, bool); 12] = [
    ("classifier_path", true),
    ("features_path", true),
    ("output_dir", true),
    ("sample_id", true),
    ("clf_group", true),
    ("cv_predictions_path", true),
    ("compress_tsv_files", false),
    ("force_plot", false),
    ("rename_classes", true),
    ("log_to_file", false),
    ("log_path", true),
    ("log_format", true),
];

const USAGE: &str = "usage: cuppa-rs predict [-h] --classifier_path CLASSIFIER_PATH --features_path FEATURES_PATH --output_dir OUTPUT_DIR [--sample_id SAMPLE_ID] [--clf_group CLF_GROUP] [--cv_predictions_path CV_PREDICTIONS_PATH] [--compress_tsv_files] [--force_plot] [--rename_classes RENAME_CLASSES] [--log_to_file] [--log_path LOG_PATH] [--log_format LOG_FORMAT]";

/// Parse errors carry argparse's exit status 2
pub struct UsageError(pub String);

/// argparse: `--name value`, `--name=value`, unique prefixes of the names, flags without values
pub fn parse_args(argv: &[String]) -> std::result::Result<Option<PredictArgs>, UsageError> {
    let mut values: Vec<(String, String)> = Vec::new();
    let mut i = 0;
    while i < argv.len() {
        let a = &argv[i];
        if a == "-h" || a == "--help" {
            println!("{USAGE}");
            return Ok(None);
        }
        let Some(body) = a.strip_prefix("--") else {
            return Err(UsageError(format!("unrecognized arguments: {a}")));
        };
        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (body, None),
        };
        let exact: Vec<&(&str, bool)> = OPTIONS.iter().filter(|(n, _)| *n == name).collect();
        let cands: Vec<&(&str, bool)> = if exact.is_empty() { OPTIONS.iter().filter(|(n, _)| n.starts_with(name)).collect() } else { exact };
        if cands.is_empty() {
            return Err(UsageError(format!("unrecognized arguments: {a}")));
        }
        if cands.len() > 1 {
            let names: Vec<String> = cands.iter().map(|(n, _)| format!("--{n}")).collect();
            return Err(UsageError(format!("ambiguous option: --{name} could match {}", names.join(", "))));
        }
        let (opt, takes) = *cands[0];
        if takes {
            let v = match inline {
                Some(v) => v,
                None => {
                    i += 1;
                    match argv.get(i) {
                        Some(v) if !v.starts_with("--") || v == "--" => v.clone(),
                        _ => return Err(UsageError(format!("argument --{opt}: expected one argument"))),
                    }
                }
            };
            values.push((opt.to_string(), v));
        } else {
            if inline.is_some() {
                return Err(UsageError(format!("argument --{opt}: ignored explicit argument '{}'", inline.unwrap())));
            }
            values.push((opt.to_string(), "true".to_string()));
        }
        i += 1;
    }
    let get = |n: &str| values.iter().rev().find(|(k, _)| k == n).map(|(_, v)| v.clone());
    let mut missing: Vec<String> = Vec::new();
    for req in ["classifier_path", "features_path", "output_dir"] {
        if get(req).is_none() {
            missing.push(format!("--{req}"));
        }
    }
    if !missing.is_empty() {
        return Err(UsageError(format!("the following arguments are required: {}", missing.join(", "))));
    }
    let rename_classes = match get("rename_classes") {
        None => None,
        Some(s) => {
            let mut pairs = Vec::new();
            for pair in s.split(',') {
                let parts: Vec<&str> = pair.split('=').collect();
                if parts.len() != 2 {
                    return Err(UsageError(format!("argument --rename_classes: invalid comma_and_equals_sep_str_to_dict value: '{s}'")));
                }
                pairs.push((parts[0].to_string(), parts[1].to_string()));
            }
            Some(pairs)
        }
    };
    Ok(Some(PredictArgs {
        classifier_path: get("classifier_path").unwrap(),
        features_path: get("features_path").unwrap(),
        output_dir: get("output_dir").unwrap(),
        sample_id: get("sample_id"),
        clf_group: get("clf_group").unwrap_or_else(|| "all".to_string()),
        cv_predictions_path: get("cv_predictions_path"),
        compress_tsv_files: get("compress_tsv_files").is_some(),
        force_plot: get("force_plot").is_some(),
        rename_classes,
        log_to_file: get("log_to_file").is_some(),
        log_path: get("log_path"),
        log_format: get("log_format").unwrap_or_else(|| pylog::DEFAULT_FORMAT.to_string()),
    }))
}

/// `python -m cuppa.predict <argv>`: the exit status
pub fn main(argv: &[String]) -> i32 {
    let args = match parse_args(argv) {
        Ok(Some(a)) => a,
        Ok(None) => return 0,
        Err(UsageError(msg)) => {
            eprintln!("{USAGE}");
            eprintln!("cuppa-rs predict: error: {msg}");
            return 2;
        }
    };
    match run(&args) {
        Ok(()) => 0,
        Err(e) => {
            pylog::error("PredictionRunner", "run", &format!("{e:#}"));
            1
        }
    }
}

fn join(dir: &str, file: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { format!("{dir}{file}") } else { format!("{dir}/{file}") }
}

/// PredictionRunner.run
pub fn run(a: &PredictArgs) -> Result<()> {
    let log_path = if a.log_to_file && a.log_path.is_none() { Some(join(&a.output_dir, "predict.log")) } else { a.log_path.clone() };
    pylog::init(&a.log_format, log_path.as_deref()).map_err(|e| anyhow!("cannot open log file: {e}"))?;
    let clf_group = a.clf_group.to_lowercase();

    pylog::info("CuppaClassifier", "from_file", &format!("Loading classifier from: {}", a.classifier_path));
    if !a.classifier_path.ends_with(".pickle") && !a.classifier_path.ends_with(".pickle.gz") && !a.classifier_path.ends_with(".bin.gz") {
        bail!("ValueError: Path must end with .pickle or .pickle.gz");
    }
    let model = model::load_for_path(&a.classifier_path)?;

    let feats = features::load(&a.features_path, a.sample_id.as_deref())?;
    let x = features::fill_missing(&feats, &model);
    let n_samples = x.samples.len();
    if n_samples == 1 && a.sample_id.is_none() {
        pylog::info("PredictionRunner", "run", "`sample_id` must be provided when predicting on one sample");
        bail!("ValueError");
    }

    // get_predictions
    let mut pred = match &a.cv_predictions_path {
        None => {
            let results = classifier::predict(&model, &x)?;
            output::build(&model, &x, &results)
        }
        Some(path) => cv::predict_with_cv(&model, &x, path)?,
    };
    match clf_group.as_str() {
        "all" => {}
        "dna" | "rna" => {
            pylog::info("CuppaPrediction", "subset_probs_by_clf_groups", &format!("Subsetting probabilties for clf_groups: {clf_group}"));
            pred.rows.retain(|r| r.data_type != "prob" || r.clf_group == Some(if clf_group == "dna" { "dna" } else { "rna" }));
        }
        _ => {
            pylog::error("PredictionRunner", "get_predictions", "`clf_group` must be one of: 'all', 'dna', 'rna'");
            bail!("ValueError");
        }
    }

    // do_rename_classes
    let mut model = model;
    if let Some(renames) = &a.rename_classes {
        let shown: Vec<String> = renames.iter().map(|(k, v)| format!("'{k}': '{v}'")).collect();
        pylog::info("PredictionRunner", "do_rename_classes", &format!("Renaming classes {{old: new}}: {{{}}}", shown.join(", ")));
        let invalid: Vec<&str> = renames.iter().map(|(k, _)| k.as_str()).filter(|k| !pred.classes.iter().any(|c| c == k)).collect();
        if !invalid.is_empty() {
            pylog::error("PredictionRunner", "do_rename_classes", &format!("Failed to rename non-existing classes: {}", invalid.join(", ")));
            bail!("ValueError");
        }
        for c in pred.classes.iter_mut() {
            if let Some((_, new)) = renames.iter().find(|(k, _)| k == c) {
                *c = new.clone();
            }
        }
        for row in model.cv.iter_mut() {
            if let Some((_, new)) = renames.iter().find(|(k, _)| *k == row.class) {
                row.class = new.clone();
            }
        }
    }

    // get_pred_summ, get_vis_data (built before anything is written)
    let show_extra_info = clf_group != "rna";
    pylog::info("CuppaPredSummaryBuilder", "build", "Summarizing predictions");
    let summary = output::summarize(&pred, show_extra_info)?;
    let vis_samples: Vec<usize> = match &a.sample_id {
        Some(id) if n_samples > 1 => {
            let found: Vec<usize> = (0..n_samples).filter(|s| x.samples[*s] == *id).collect();
            if found.is_empty() {
                pylog::error("CuppaPrediction", "get_rows", "Some `keys` were not found in index level `sample_id`");
                bail!("LookupError");
            }
            pylog::info("CuppaVisDataBuilder", "build", &format!("Building vis data for sample '{id}' from {n_samples} samples"));
            found
        }
        _ => {
            if n_samples == 1 {
                pylog::info("CuppaVisDataBuilder", "build", &format!("Building vis data for sample '{}'", x.samples[0]));
            } else {
                pylog::info("CuppaVisDataBuilder", "build", &format!("Building vis data for {n_samples} samples"));
            }
            (0..n_samples).collect()
        }
    };
    // the vis data is built (and checked) in every mode, written only outside RNA mode
    let vis_tmp = output::vis_rows(&model, &pred, &vis_samples)?;

    let out_dir = std::path::Path::new(&a.output_dir);
    if !out_dir.exists() {
        std::fs::create_dir(out_dir).map_err(|e| anyhow!("FileNotFoundError: [Errno 2] No such file or directory: '{}' ({e})", a.output_dir))?;
        pylog::info("PredictionRunner", "run", &format!("Created output dir: {}", a.output_dir));
    }
    let affix = |name: &str| -> String {
        let mut f = match &a.sample_id {
            Some(id) => format!("{id}.{name}"),
            None => name.to_string(),
        };
        if a.compress_tsv_files && f.ends_with(".tsv") {
            f.push_str(".gz");
        }
        join(&a.output_dir, &f)
    };
    let pred_summ_path = affix("cuppa.pred_summ.tsv");
    pylog::info("CuppaPredSummary", "to_tsv", &format!("Writing prediction summary to: {pred_summ_path}"));
    output::write_summary(&pred_summ_path, &pred, &summary, show_extra_info)?;

    if clf_group == "rna" {
        pylog::warning("PredictionRunner", "run", "Skipping plotting predictions as it is unsupported in RNA-only mode");
        return Ok(());
    }
    let vis_path = affix("cuppa.vis_data.tsv");
    pylog::info("CuppaVisData", "to_tsv", &format!("Writing vis_data to: {vis_path}"));
    output::write_vis_rows(&vis_path, &pred, &vis_tmp)?;

    let plot_name = if n_samples == 1 || a.sample_id.is_some() { "cuppa.vis.png" } else { "cuppa.vis.pdf" };
    let plot_path = affix(plot_name);
    let plotted_samples = vis_samples.len();
    const PLOT_MAX_SAMPLES: usize = 10;
    if plotted_samples > 1 && !plot_path.ends_with(".pdf") {
        pylog::error("CuppaVisPlotter", "_check_plot_path_extension", "`plot_path` must end with .pdf for multi-sample plotting");
        bail!("Exception");
    }
    if plotted_samples > PLOT_MAX_SAMPLES && !a.force_plot {
        pylog::warning(
            "CuppaVisPlotter",
            "_check_should_plot",
            &format!("Skipping plotting predictions for {plotted_samples} (>{PLOT_MAX_SAMPLES}) samples. Please use arg --force_plot to plot anyway"),
        );
        return Ok(());
    }
    // pycuppa runs `Rscript --vanilla plot_predictions.R <vis data> <plot>`; the port draws it natively. A plot that
    // R could not draw (no snv_count feature, no CV performance) fails the run as pycuppa's CalledProcessError does,
    // after the data files are written.
    pylog::info("RscriptExecutor", "run", &format!("Drawing the plot natively: {plot_path}"));
    if let Err(e) = crate::plot::draw(&vis_path, &plot_path) {
        pylog::error("RscriptExecutor", "raise_if_error", &format!("Plotting failed: {e:#}"));
        bail!("CalledProcessError");
    }
    Ok(())
}
