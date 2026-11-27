//! `--cv_predictions_path`: the predictions of samples found in a CuppaPrediction file (as
//! CuppaPrediction.to_tsv writes it during training: two header rows, "index" over the index columns and
//! the supertypes over the class columns) are taken from it; the other samples are predicted
//! (PredictionRunner.get_predictions).

use std::collections::HashSet;

use anyhow::{anyhow, bail, Result};

use super::classifier;
use super::features::{read_text_lines, Prepared};
use super::model::Model;
use super::output::{self, PredRow, Prediction};
use super::pylog;

fn static_str(s: &str, allowed: &'static [&'static str], what: &str) -> Result<&'static str> {
    allowed.iter().find(|a| **a == s).copied().ok_or_else(|| anyhow!("cross-validation predictions: unexpected {what} '{s}'"))
}

const DATA_TYPES: &[&str] = &["prob", "feat_contrib", "sig_quantile"];
const CLF_GROUPS: &[&str] = &["combined", "dna", "rna"];
const CLF_NAMES: &[&str] = &["combined", "dna_combined", "gen_pos", "snv96", "event", "rna_combined", "gene_exp", "alt_sj"];

struct CvRow {
    sample: String,
    data_type: &'static str,
    clf_group: Option<&'static str>,
    clf_name: Option<&'static str>,
    feat_name: Option<String>,
    feat_value: Option<f64>,
    values: Vec<f64>,
}

fn parse_float(s: &str) -> Result<f64> {
    if s.is_empty() || s == "nan" || s == "NaN" {
        return Ok(f64::NAN);
    }
    s.parse::<f64>().map_err(|_| anyhow!("cross-validation predictions: not a number: '{s}'"))
}

/// CuppaPrediction.from_tsv: (class columns, rows)
fn read_cv(path: &str) -> Result<(Vec<String>, Vec<CvRow>)> {
    let lines = read_text_lines(path)?;
    if lines.len() < 2 {
        bail!("cross-validation predictions {path}: no header");
    }
    let top: Vec<&str> = lines[0].split('\t').collect();
    let names: Vec<&str> = lines[1].split('\t').collect();
    if top.len() != names.len() {
        bail!("cross-validation predictions {path}: header rows of different lengths");
    }
    let col = |n: &str| names.iter().position(|c| *c == n).filter(|i| top[*i] == "index");
    let need = |n: &str| col(n).ok_or_else(|| anyhow!("KeyError: \"Missing required columns: {n}\""));
    let (si, di, gi, ni, fi, vi) = (need("sample_id")?, need("data_type")?, need("clf_group")?, need("clf_name")?, need("feat_name")?, need("feat_value")?);
    let class_cols: Vec<usize> = (0..names.len()).filter(|i| top[*i] != "index").collect();
    let classes: Vec<String> = class_cols.iter().map(|i| names[*i].to_string()).collect();
    let mut rows = Vec::new();
    for line in &lines[2..] {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let opt = |s: &str| if s.is_empty() { None } else { Some(s.to_string()) };
        rows.push(CvRow {
            sample: get(si).to_string(),
            data_type: static_str(get(di), DATA_TYPES, "data_type")?,
            clf_group: match get(gi) {
                "" => None,
                g => Some(static_str(g, CLF_GROUPS, "clf_group")?),
            },
            clf_name: match get(ni) {
                "" => None,
                n => Some(static_str(n, CLF_NAMES, "clf_name")?),
            },
            feat_name: opt(get(fi)),
            feat_value: match get(vi) {
                "" => None,
                v => Some(parse_float(v)?),
            },
            values: class_cols.iter().map(|i| parse_float(get(*i))).collect::<Result<Vec<f64>>>()?,
        });
    }
    Ok((classes, rows))
}

pub fn predict_with_cv(m: &Model, x: &Prepared, path: &str) -> Result<Prediction> {
    pylog::info("PredictionRunner", "get_predictions", &format!("Loading cross-validation predictions from: {path}"));
    let (cv_classes, cv_rows) = read_cv(path)?;
    let cv_samples: HashSet<&str> = cv_rows.iter().map(|r| r.sample.as_str()).collect();
    let n = x.samples.len();
    let new_idx: Vec<usize> = (0..n).filter(|s| !cv_samples.contains(x.samples[*s].as_str())).collect();
    let cv_idx: Vec<usize> = (0..n).filter(|s| cv_samples.contains(x.samples[*s].as_str())).collect();

    let (classes, _) = output::class_order(&m.classes);
    let mut cv_sorted = cv_classes.clone();
    cv_sorted.sort();
    let mut model_sorted = classes.clone();
    model_sorted.sort();
    if cv_sorted != model_sorted {
        bail!("cross-validation predictions {path}: their classes differ from the classifier's (not supported)");
    }
    let col_map: Vec<usize> = classes.iter().map(|c| cv_classes.iter().position(|k| k == c).unwrap()).collect();

    let mut per_sample: Vec<Vec<PredRow>> = vec![Vec::new(); n];
    if !new_idx.is_empty() {
        pylog::info("PredictionRunner", "get_predictions", &format!("Computing predictions for {} / {} samples", new_idx.len(), n));
        let sub = subset(x, &new_idx);
        let results = classifier::predict(m, &sub)?;
        let p = output::build(m, &sub, &results);
        for r in p.rows {
            let s = new_idx[r.sample];
            per_sample[s].push(PredRow { sample: s, ..r });
        }
    }
    if !cv_idx.is_empty() {
        for r in &cv_rows {
            if let Some(s) = (0..n).find(|s| x.samples[*s] == r.sample) {
                per_sample[s].push(PredRow {
                    sample: s,
                    data_type: r.data_type,
                    clf_group: r.clf_group,
                    clf_name: r.clf_name,
                    feat_name: r.feat_name.clone(),
                    feat_value: r.feat_value,
                    values: col_map.iter().map(|i| r.values[*i]).collect(),
                });
            }
        }
        pylog::info(
            "PredictionRunner",
            "get_predictions",
            &format!("Using pre-computed cross-validation predictions for {} / {} samples", cv_idx.len(), n),
        );
    }
    let (classes, supertypes) = output::class_order(&m.classes);
    Ok(Prediction { samples: x.samples.clone(), classes, supertypes, rows: per_sample.into_iter().flatten().collect() })
}

/// the rows of the given samples
fn subset(x: &Prepared, idx: &[usize]) -> Prepared {
    Prepared {
        samples: idx.iter().map(|s| x.samples[*s].clone()).collect(),
        dna: idx.iter().map(|s| x.dna[*s].clone()).collect(),
        rna: idx.iter().map(|s| x.rna[*s].clone()).collect(),
        sex_missing: idx.iter().map(|s| x.sex_missing[*s]).collect(),
        n_gen_pos: x.n_gen_pos,
        n_snv96: x.n_snv96,
        n_event: x.n_event,
        n_sig: x.n_sig,
        n_gene_exp: x.n_gene_exp,
    }
}
