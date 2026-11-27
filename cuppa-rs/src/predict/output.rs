//! pycuppa's prediction table (CuppaPrediction), its summary (CuppaPredSummaryBuilder, written to
//! `<sample>.cuppa.pred_summ.tsv`) and the visualisation data (CuppaVisDataBuilder, written to
//! `<sample>.cuppa.vis_data.tsv`), as pandas orders, rounds and writes them.

use std::io::Write;

use anyhow::{bail, Context, Result};

use super::classifier::SampleResult;
use super::features::{Prepared, PREDICT_NA_FILL_VALUE};
use super::model::Model;
use super::numpy::{argsort, round};
use super::pyfmt::{csv_float, repr};

/// One row of the prediction table
#[derive(Clone)]
pub struct PredRow {
    pub sample: usize,
    pub data_type: &'static str,
    pub clf_group: Option<&'static str>,
    pub clf_name: Option<&'static str>,
    pub feat_name: Option<String>,
    pub feat_value: Option<f64>,
    /// one value per class column (in `Prediction::classes` order), NaN where missing
    pub values: Vec<f64>,
}

pub struct Prediction {
    pub samples: Vec<String>,
    /// the class columns, ordered as CuppaPrediction.from_data_frame orders them
    pub classes: Vec<String>,
    /// the supertype of each class column (the first header row of CuppaPrediction.to_tsv)
    #[allow(dead_code)]
    pub supertypes: Vec<String>,
    pub rows: Vec<PredRow>,
}

/// `str.extract("(^[A-Za-z/ ]+)")`
fn supertype(c: &str) -> Option<String> {
    let n = c.chars().take_while(|ch| ch.is_ascii_alphabetic() || *ch == '/' || *ch == ' ').count();
    if n == 0 { None } else { Some(c.chars().take(n).collect()) }
}

/// `str.replace("^[A-Za-z/ ]+: ", "", regex=True)`
fn subtype(c: &str) -> String {
    let n = c.chars().take_while(|ch| ch.is_ascii_alphabetic() || *ch == '/' || *ch == ' ').count();
    // the greedy class may give back characters for ": " to follow; the class holds no ':' so it cannot
    let rest: String = c.chars().skip(n).collect();
    if n > 0 && rest.starts_with(": ") { rest[2..].to_string() } else { c.to_string() }
}

/// CuppaPrediction._columns_to_class_metadata: sorted by supertype, "Other" subtypes last, subtype
pub fn class_order(classes: &[String]) -> (Vec<String>, Vec<String>) {
    let mut meta: Vec<(String, bool, String, String)> = classes
        .iter()
        .map(|c| {
            let sup = supertype(c).unwrap_or_default();
            let sub = subtype(c);
            (sup, sub.ends_with("Other"), sub, c.clone())
        })
        .collect();
    meta.sort_by(|a, b| (&a.0, a.1, &a.2).cmp(&(&b.0, b.1, &b.2)));
    (meta.iter().map(|m| m.3.clone()).collect(), meta.iter().map(|m| m.0.clone()).collect())
}

fn clf_group_of(name: &str) -> &'static str {
    match name {
        "combined" => "combined",
        "dna_combined" | "gen_pos" | "snv96" | "event" => "dna",
        _ => "rna",
    }
}

/// CuppaPredictionBuilder.build: per sample, the probabilities (by classifier group and name), the event
/// feature contributions and the signature quantiles
pub fn build(m: &Model, x: &Prepared, results: &[SampleResult]) -> Prediction {
    let (classes, supertypes) = class_order(&m.classes);
    let col_of = |names: &[String]| -> Vec<Option<usize>> { classes.iter().map(|c| names.iter().position(|n| n == c)).collect() };
    let dna_map = col_of(&m.dna_combined.lr.classes);
    let rna_map = col_of(&m.rna_combined.lr.classes);
    let gp_map = col_of(&m.gen_pos.lr.classes);
    let sn_map = col_of(&m.snv96.lr.classes);
    let ev_map = col_of(&m.event.lr.classes);
    let ge_map = col_of(&m.gene_exp.lr.classes);
    let sj_map = col_of(&m.alt_sj.lr.classes);
    let sig_map = col_of(&m.sig.classes);
    let spread = |v: &[f64], map: &[Option<usize>]| -> Vec<f64> { map.iter().map(|i| i.map_or(f64::NAN, |i| v[i])).collect() };

    let mut rows = Vec::new();
    for (s, r) in results.iter().enumerate() {
        let mut probs: Vec<(&'static str, Vec<f64>)> = Vec::new();
        if let Some(c) = &r.combined {
            // the combined classes are the rna_combined ones
            probs.push(("combined", spread(c, &rna_map)));
        }
        probs.push(("dna_combined", spread(&r.dna_combined, &dna_map)));
        probs.push(("gen_pos", spread(&r.gen_pos, &gp_map)));
        probs.push(("snv96", spread(&r.snv96, &sn_map)));
        probs.push(("event", spread(&r.event, &ev_map)));
        if let Some(v) = &r.rna_combined {
            probs.push(("rna_combined", spread(v, &rna_map)));
        }
        if let Some(v) = &r.gene_exp {
            probs.push(("gene_exp", spread(v, &ge_map)));
        }
        if let Some(v) = &r.alt_sj {
            probs.push(("alt_sj", spread(v, &sj_map)));
        }
        for (name, values) in probs {
            rows.push(PredRow { sample: s, data_type: "prob", clf_group: Some(clf_group_of(name)), clf_name: Some(name), feat_name: None, feat_value: None, values });
        }
        let ev = x.event(s);
        for (f, fname) in m.event.lr.features.iter().enumerate() {
            let per_class: Vec<f64> = (0..m.event.lr.n_classes()).map(|c| r.contribs[c][f]).collect();
            rows.push(PredRow {
                sample: s,
                data_type: "feat_contrib",
                clf_group: Some("dna"),
                clf_name: Some("event"),
                feat_name: Some(fname.clone()),
                feat_value: Some(ev[m.req_event.iter().position(|n| n == fname).unwrap()]),
                values: spread(&per_class, &ev_map),
            });
        }
        let sig = x.sig(s);
        for (f, fname) in m.sig.features.iter().enumerate() {
            // feat_name.replace("^\w+[.]", "", regex=True)
            let stripped = match fname.find('.') {
                Some(p) if p > 0 && fname[..p].chars().all(|c| c.is_alphanumeric() || c == '_') => fname[p + 1..].to_string(),
                _ => fname.clone(),
            };
            rows.push(PredRow {
                sample: s,
                data_type: "sig_quantile",
                clf_group: None,
                clf_name: None,
                feat_name: Some(stripped),
                feat_value: Some(sig[f]),
                values: spread(&r.sig_quantiles[f], &sig_map),
            });
        }
    }
    Prediction { samples: x.samples.clone(), classes, supertypes, rows }
}

/// a text cell of DataFrame.to_csv (QUOTE_MINIMAL), NaN as ""
fn csv_text(s: Option<&str>) -> String {
    match s {
        None => String::new(),
        Some(t) if t.contains(['\t', '"', '\n', '\r']) => format!("\"{}\"", t.replace('"', "\"\"")),
        Some(t) => t.to_string(),
    }
}

/// get_top_cols for one row: names and values in descending order (np.fliplr(np.argsort(row)), NaN as
/// -inf), the -inf ones as NaN
fn top_cols(row: &[f64], top_n: usize) -> Vec<Option<usize>> {
    let filled: Vec<f64> = row.iter().map(|v| if v.is_nan() { f64::NEG_INFINITY } else { *v }).collect();
    let order = argsort(&filled);
    order.iter().rev().take(top_n).map(|i| if filled[*i] == f64::NEG_INFINITY { None } else { Some(*i) }).collect()
}

/// `str(x)` of a float cell after DataFrame.astype(str) (numpy's repr; "nan" for NaN)
fn astype_str(x: f64) -> String {
    repr(x)
}

/// CuppaPredSummaryBuilder._signif(x, precision)
fn signif(x: f64, precision: i32) -> f64 {
    let xp = if x.is_finite() && x != 0.0 { x.abs() } else { 10f64.powf((precision - 1) as f64) };
    let mags = 10f64.powf((precision - 1) as f64 - xp.log10().floor());
    round(x * mags, 0) / mags
}

/// `_collapse_columns_to_strings`: None stands for NaN, which spreads through `+`
fn collapse(cols: &[Option<String>], sep: &str) -> Option<String> {
    let mut s = cols.first().cloned().flatten();
    for c in cols.iter().skip(1) {
        let with_sep = match c {
            None => None,
            Some(t) if t.is_empty() => Some(String::new()),
            Some(t) => Some(format!("{sep}{t}")),
        };
        s = match (s, with_sep) {
            (Some(a), Some(b)) => Some(a + &b),
            _ => None,
        };
    }
    s
}

pub struct SummaryRow {
    pub sample: usize,
    pub clf_group: &'static str,
    pub clf_name: &'static str,
    pub classes: Vec<Option<String>>,
    pub probs: Vec<f64>,
    pub extra_info: Option<Option<String>>,
    pub extra_format: Option<String>,
}

/// CuppaPredSummaryBuilder.build (top 3 classes and their probabilities, with the extra information of
/// the snv96 and event rows)
pub fn summarize(p: &Prediction, show_extra_info: bool) -> Result<Vec<SummaryRow>> {
    const TOP_N: usize = 3;
    let mut out: Vec<SummaryRow> = Vec::new();
    for r in p.rows.iter().filter(|r| r.data_type == "prob") {
        let top = top_cols(&r.values, TOP_N);
        let mut classes: Vec<Option<String>> = top.iter().map(|i| i.map(|i| p.classes[i].clone())).collect();
        let mut probs: Vec<f64> = top.iter().map(|i| i.map_or(f64::NAN, |i| r.values[i])).collect();
        while classes.len() < TOP_N {
            classes.push(None);
            probs.push(f64::NAN);
        }
        out.push(SummaryRow {
            sample: r.sample,
            clf_group: r.clf_group.unwrap(),
            clf_name: r.clf_name.unwrap(),
            classes,
            probs,
            extra_info: None,
            extra_format: None,
        });
    }
    if !show_extra_info {
        return Ok(out);
    }

    for s in 0..p.samples.len() {
        // signatures (extra_info_sigs) -------------------------------------------------------------------
        let tmb = p
            .rows
            .iter()
            .find(|r| r.sample == s && r.data_type == "feat_contrib" && r.feat_name.as_deref().is_some_and(|n| n.ends_with("snv_count")))
            .and_then(|r| r.feat_value);
        let sigs: Vec<&PredRow> = p
            .rows
            .iter()
            .filter(|r| r.sample == s && r.data_type == "sig_quantile" && r.feat_value != Some(PREDICT_NA_FILL_VALUE))
            .collect();
        let sig_strings: Option<(Option<String>, String)> = if sigs.is_empty() {
            None
        } else {
            let Some(tmb) = tmb else { bail!("KeyError: no event.tmb.snv_count contribution for {}", p.samples[s]) };
            if tmb.is_nan() {
                bail!("pandas.errors.IntCastingNaNError: Cannot convert non-finite values (NA or inf) to integer");
            }
            let tmb = tmb.trunc() as i64 as f64;
            // rank(method="first", ascending=False): by value descending, ties in order of appearance
            let mut order: Vec<usize> = (0..sigs.len()).collect();
            order.sort_by(|a, b| {
                let (va, vb) = (sigs[*a].feat_value.unwrap(), sigs[*b].feat_value.unwrap());
                vb.partial_cmp(&va).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(b))
            });
            let labels: Vec<Option<String>> = order
                .iter()
                .map(|i| {
                    let r = sigs[*i];
                    let v = r.feat_value.unwrap();
                    let value = round(v, 1);
                    let perc = round((v / tmb) * 100.0, 1);
                    if value == 0.0 {
                        Some(String::new())
                    } else {
                        Some(format!("{}={},{}%", r.feat_name.as_deref().unwrap(), astype_str(value), astype_str(perc)))
                    }
                })
                .collect();
            let strings = collapse(&labels, "; ");
            let fmt = match &strings {
                Some(t) if t.is_empty() => String::new(),
                _ => "{sig_name}={count},{percent}".to_string(),
            };
            Some((strings, fmt))
        };
        if let Some((strings, fmt)) = sig_strings {
            if let Some(row) = out.iter_mut().find(|o| o.sample == s && o.clf_name == "snv96") {
                row.extra_info = Some(strings);
                row.extra_format = Some(fmt);
            }
        }

        // event features (extra_info_feat) ----------------------------------------------------------------
        let Some(event_row) = out.iter().position(|o| o.sample == s && o.clf_name == "event") else { continue };
        let contribs: Vec<&PredRow> = p
            .rows
            .iter()
            .filter(|r| r.sample == s && r.data_type == "feat_contrib" && r.feat_name.as_deref() != Some("_prior"))
            .collect();
        let mut per_class: Vec<Option<String>> = Vec::new();
        for class in out[event_row].classes.clone().into_iter().flatten() {
            let ci = p.classes.iter().position(|c| *c == class).unwrap();
            let row: Vec<f64> = contribs.iter().map(|r| r.values[ci]).collect();
            let top = top_cols(&row, TOP_N);
            let cells: Vec<Option<String>> = top
                .iter()
                .map(|t| {
                    t.map(|i| {
                        let r = contribs[i];
                        let name = r.feat_name.as_deref().unwrap();
                        let name = name.strip_prefix("event.").unwrap_or(name);
                        let value = signif(r.feat_value.unwrap(), 2);
                        format!("{}={},{}", name, astype_str(value), astype_str(round(row[i], 1)))
                    })
                })
                .collect();
            per_class.push(collapse(&cells, "; "));
        }
        let wide: Vec<Option<String>> = per_class.iter().enumerate().map(|(k, s)| s.as_ref().map(|t| format!("{}: {}", k + 1, t))).collect();
        out[event_row].extra_info = Some(collapse(&wide, " || "));
        out[event_row].extra_format = Some("{pred_class_n}: {feat_name}={value},{log_odds}".to_string());
    }
    Ok(out)
}

pub fn write_summary(path: &str, p: &Prediction, rows: &[SummaryRow], show_extra_info: bool) -> Result<()> {
    let mut text = String::new();
    text.push_str("sample_id\tclf_group\tclf_name\tpred_class_1\tpred_class_2\tpred_class_3\tpred_prob_1\tpred_prob_2\tpred_prob_3");
    if show_extra_info {
        text.push_str("\textra_info\textra_info_format");
    }
    text.push('\n');
    for r in rows {
        let mut f: Vec<String> = vec![csv_text(Some(&p.samples[r.sample])), r.clf_group.to_string(), r.clf_name.to_string()];
        f.extend(r.classes.iter().map(|c| csv_text(c.as_deref())));
        f.extend(r.probs.iter().map(|v| csv_float(*v)));
        if show_extra_info {
            f.push(csv_text(r.extra_info.clone().flatten().as_deref()));
            f.push(csv_text(r.extra_format.as_deref()));
        }
        text.push_str(&f.join("\t"));
        text.push('\n');
    }
    write_text(path, &text)
}

/// writes a text file, gzip-compressed if the name ends with .gz (pandas' compression='infer')
pub fn write_text(path: &str, text: &str) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("cannot create {path}"))?;
    if path.ends_with(".gz") {
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(text.as_bytes())?;
        gz.finish()?;
    } else {
        let mut w = std::io::BufWriter::new(file);
        w.write_all(text.as_bytes())?;
        w.flush()?;
    }
    Ok(())
}

struct VisRow {
    sample: Option<usize>,
    data_type: &'static str,
    clf_group: Option<&'static str>,
    clf_name: Option<String>,
    feat_name: Option<String>,
    feat_value: Option<f64>,
    cancer_type: String,
    data_value: f64,
}

pub struct VisData {
    rows: Vec<VisRow>,
    ranks: Vec<usize>,
    groups: Vec<usize>,
}

/// CuppaVisDataBuilder.build for the given samples
pub fn vis_rows(m: &Model, p: &Prediction, samples: &[usize]) -> Result<VisData> {
    const MIN_DRIVER_LIKELIHOOD: f64 = 0.2;
    const MIN_DRIVER_CONTRIB: f64 = 0.2;
    if !p.rows.iter().any(|r| r.feat_name.as_deref() == Some("event.tmb.snv_count")) {
        bail!("LookupError: feature contributions / values for 'event.tmb.snv_count' do not exist");
    }
    let mut vis: Vec<VisRow> = Vec::new();
    let long = |vis: &mut Vec<VisRow>, r: &PredRow, odds: bool, strip: bool| {
        for (c, v) in p.classes.iter().zip(&r.values) {
            let feat_name = if strip {
                // feat_name.str.split(".", n=1).map(lambda x: x[1])
                r.feat_name.as_deref().map(|n| n.split_once('.').map(|(_, b)| b.to_string()).unwrap_or_default())
            } else {
                r.feat_name.clone()
            };
            vis.push(VisRow {
                sample: Some(r.sample),
                data_type: r.data_type,
                clf_group: r.clf_group,
                clf_name: r.clf_name.map(String::from),
                feat_name,
                feat_value: r.feat_value,
                cancer_type: c.clone(),
                data_value: if odds { std::f64::consts::E.powf(*v) } else { *v },
            });
        }
    };
    let rows: Vec<&PredRow> = p.rows.iter().filter(|r| samples.contains(&r.sample)).collect();
    for r in rows.iter().filter(|r| r.data_type == "prob") {
        long(&mut vis, r, false, false);
    }
    for r in rows.iter().filter(|r| r.data_type == "sig_quantile") {
        long(&mut vis, r, false, false);
    }
    let contribs: Vec<&PredRow> = rows.iter().filter(|r| r.data_type == "feat_contrib").copied().collect();
    let subset = |pattern: &str| -> Result<Vec<&PredRow>> {
        let v: Vec<&PredRow> = contribs.iter().filter(|r| r.feat_name.as_deref().is_some_and(|n| n.contains(pattern))).copied().collect();
        if v.is_empty() {
            bail!("LookupError: No features found with regex '{pattern}'");
        }
        Ok(v)
    };
    let row_max = |r: &PredRow| -> f64 {
        // DataFrame.max(axis=1), skipping NaN
        let mut mx = f64::NAN;
        for v in &r.values {
            if !v.is_nan() && (mx.is_nan() || *v > mx) {
                mx = *v;
            }
        }
        mx
    };
    let mut selected: Vec<&PredRow> = Vec::new();
    selected.extend(subset("tmb")?);
    selected.extend(subset("trait")?);
    for pattern in ["fusion", "virus"] {
        selected.extend(subset(pattern)?.into_iter().filter(|r| row_max(r) > 0.0 && r.feat_value.unwrap_or(f64::NAN) > 0.0));
    }
    selected.extend(subset("driver")?.into_iter().filter(|r| row_max(r) >= MIN_DRIVER_CONTRIB || r.feat_value.unwrap_or(f64::NAN) >= MIN_DRIVER_LIKELIHOOD));
    selected.extend(subset("sv")?);
    for r in selected {
        long(&mut vis, r, true, true);
    }

    // cross-validation performance of dna_combined: n_total, recall and precision of the predicted classes
    let cv_rows: Vec<&super::model::CvRow> = m.cv.iter().filter(|c| c.clf_name == "dna_combined" && p.classes.contains(&c.class)).collect();
    let mut uniq: Vec<&str> = cv_rows.iter().map(|c| c.class.as_str()).collect();
    uniq.sort();
    uniq.dedup();
    if uniq.len() != p.classes.len() {
        bail!("LookupError: `self.cv_performance` is missing some classes found in `self.predictions`");
    }
    for (name, get) in [
        ("n_total", (|c: &super::model::CvRow| c.n_total as f64) as fn(&super::model::CvRow) -> f64),
        ("recall", |c: &super::model::CvRow| c.recall),
        ("precision", |c: &super::model::CvRow| c.precision),
    ] {
        for c in &cv_rows {
            vis.push(VisRow {
                sample: None,
                data_type: "cv_performance",
                clf_group: Some("dna"),
                clf_name: Some("dna_combined".to_string()),
                feat_name: Some(name.to_string()),
                feat_value: None,
                cancer_type: c.class.clone(),
                data_value: get(c),
            });
        }
    }

    // add_ranks: groups by (sample_id, data_type, clf_name, feat_name), NaN keys included; rank by
    // data_value descending (method "first", NaN last); group numbers in order of first appearance
    let mut group_of: Vec<usize> = Vec::with_capacity(vis.len());
    let mut keys: Vec<(Option<usize>, &str, Option<&str>, Option<&str>)> = Vec::new();
    for v in &vis {
        let k = (v.sample, v.data_type, v.clf_name.as_deref(), v.feat_name.as_deref());
        let g = match keys.iter().position(|x| *x == k) {
            Some(g) => g,
            None => {
                keys.push(k);
                keys.len() - 1
            }
        };
        group_of.push(g);
    }
    let mut ranks = vec![0usize; vis.len()];
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); keys.len()];
    for (i, g) in group_of.iter().enumerate() {
        members[*g].push(i);
    }
    for mem in &members {
        let mut order = mem.clone();
        order.sort_by(|a, b| {
            let (va, vb) = (vis[*a].data_value, vis[*b].data_value);
            match (va.is_nan(), vb.is_nan()) {
                (true, true) => a.cmp(b),
                (true, false) => std::cmp::Ordering::Greater,
                (false, true) => std::cmp::Ordering::Less,
                _ => vb.partial_cmp(&va).unwrap().then(a.cmp(b)),
            }
        });
        for (rank, i) in order.iter().enumerate() {
            ranks[*i] = rank + 1;
        }
    }

    Ok(VisData { rows: vis, ranks, groups: group_of })
}

/// CuppaVisData.to_tsv
pub fn write_vis_rows(path: &str, p: &Prediction, d: &VisData) -> Result<()> {
    let (vis, ranks, group_of) = (&d.rows, &d.ranks, &d.groups);
    let mut text = String::from("sample_id\tdata_type\tclf_group\tclf_name\tfeat_name\tfeat_value\tcancer_type\tdata_value\trank\trank_group\n");
    for (i, v) in vis.iter().enumerate() {
        let f = [
            csv_text(v.sample.map(|s| p.samples[s].as_str())),
            v.data_type.to_string(),
            csv_text(v.clf_group),
            csv_text(v.clf_name.as_deref()),
            csv_text(v.feat_name.as_deref()),
            v.feat_value.map_or(String::new(), csv_float),
            csv_text(Some(&v.cancer_type)),
            csv_float(v.data_value),
            ranks[i].to_string(),
            group_of[i].to_string(),
        ];
        text.push_str(&f.join("\t"));
        text.push('\n');
    }
    write_text(path, &text)
}
