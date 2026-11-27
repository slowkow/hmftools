//! The CUPPA summary plot, `<sample>.cuppa.vis.png` (one sample) or `cuppa.vis.pdf` (several, a page each), drawn
//! natively. pycuppa draws it with the R script `cuppa/visualization/plot_predictions.R` (ggplot2 3.5.1, ggh4x
//! 0.2.8, patchwork 1.2.0 in the pipeline's image) from the vis data file; this is a port of that script on
//! hmf-plot-rs's ggplot2 and patchwork layout: the same data steps (labels, rounding, number formats, supertype
//! sums, signature quantile classes, odds ratio labels, the training set performance), the same four heat maps
//! (facet_grid2 by row group and cancer supertype with free scales and space, tiles, values, the first-row line,
//! colour bars and the quantile legend), stacked with heights by their numbers of rows, at the script's size.

use std::collections::HashMap;

use anyhow::{Context, Result, bail};
use hmf_plot::color::{GradientN, Rgba, r_colour};
use hmf_plot::ggp::facet::{self, Grid, Strip};
use hmf_plot::ggp::geom::{self, Panel};
use hmf_plot::ggp::guides::{self, AxisSpec, GuideBox, KeyGlyph, LegendSpec, LegendTheme, Side};
use hmf_plot::ggp::plot::{PlotParts, plot_table};
use hmf_plot::ggp::{Ctx, LineEl, RectEl, TextEl, TitleGrob, Version, rect_grob};
use hmf_plot::gg::rfmt;
use hmf_plot::gtable::{Bx, Canvas, Grob, Gtable};
use hmf_plot::patchwork::{self, Annotation, Layout};
use hmf_plot::rdev::{Device, Face};
use hmf_plot::{Figure, png_out::PngMeta};

/// One row of the vis data as R's read.delim reads it (empty numbers are NA).
#[derive(Clone, Debug)]
struct VisRow {
    sample_id: String,
    data_type: String,
    clf_group: String,
    clf_name: String,
    feat_name: String,
    feat_value: Option<f64>,
    cancer_type: String,
    data_value: Option<f64>,
}

fn read_vis(path: &str) -> Result<Vec<VisRow>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).with_context(|| format!("no column {name} in {path}"));
    let (cs, cdt, cg, cn, cf, cfv, cct, cdv) =
        (col("sample_id")?, col("data_type")?, col("clf_group")?, col("clf_name")?, col("feat_name")?, col("feat_value")?, col("cancer_type")?, col("data_value")?);
    let num = |s: &str| -> Option<f64> {
        match s {
            "" | "NA" => None,
            "inf" | "Inf" => Some(f64::INFINITY),
            "-inf" | "-Inf" => Some(f64::NEG_INFINITY),
            "nan" | "NaN" => None,
            _ => s.parse().ok(),
        }
    };
    let mut rows = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let g = |i: usize| f.get(i).copied().unwrap_or("").to_string();
        rows.push(VisRow {
            sample_id: g(cs),
            data_type: g(cdt),
            clf_group: g(cg),
            clf_name: g(cn),
            feat_name: g(cf),
            feat_value: num(&g(cfv)),
            cancer_type: g(cct),
            data_value: num(&g(cdv)),
        });
    }
    Ok(rows)
}

/// CANCER_TYPE_METADATA: supertype (the leading run of letters, '/' and spaces) and subtype (what follows its
/// ": ") of each cancer type, in order of appearance, and the number of subtypes per supertype.
struct CancerTypes {
    #[allow(dead_code)] // read by the tests
    types: Vec<String>,
    supertype: Vec<usize>,
    subtype: Vec<String>,
    supertypes: Vec<String>,
    subtypes_per_supertype: Vec<usize>,
    index: HashMap<String, usize>,
}

fn leading_run(s: &str) -> &str {
    let n = s.bytes().take_while(|b| b.is_ascii_alphabetic() || *b == b'/' || *b == b' ').count();
    &s[..n]
}

impl CancerTypes {
    fn new(rows: &[VisRow]) -> CancerTypes {
        let mut types: Vec<String> = Vec::new();
        let mut index = HashMap::new();
        for r in rows {
            if !index.contains_key(&r.cancer_type) {
                index.insert(r.cancer_type.clone(), types.len());
                types.push(r.cancer_type.clone());
            }
        }
        let mut supertypes: Vec<String> = Vec::new();
        let mut supertype = Vec::new();
        let mut subtype = Vec::new();
        for t in &types {
            let prefix = leading_run(t);
            let sup = match supertypes.iter().position(|s| s == prefix) {
                Some(i) => i,
                None => {
                    supertypes.push(prefix.to_string());
                    supertypes.len() - 1
                }
            };
            supertype.push(sup);
            // str_replace("^[A-Za-z/ ]+: ", ""): the run then ": "
            let sub = if !prefix.is_empty() && t[prefix.len()..].starts_with(": ") { t[prefix.len() + 2..].to_string() } else { t.clone() };
            subtype.push(sub);
        }
        let mut subtypes_per_supertype = vec![0; supertypes.len()];
        for &s in &supertype {
            subtypes_per_supertype[s] += 1;
        }
        CancerTypes { types, supertype, subtype, supertypes, subtypes_per_supertype, index }
    }
}

// -------------------------------------------------------------------------------------------------------------
// R's number formatting

/// as.character of a number (15 significant digits), as paste0 and sprintf("%s") give it.
fn r_chr(x: f64) -> String {
    if x.is_nan() {
        return "NaN".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    rfmt::r_as_character(if x == 0.0 { 0.0 } else { x })
}

fn r_round(x: f64, d: i32) -> f64 {
    rfmt::r_round(x, d)
}

/// signif(x, 1)
fn r_signif1(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let e = x.abs().log10().floor() as i32;
    r_round(x, -e)
}

/// formatC(x, format = "e", digits = d): C's %.de
fn format_c_e(x: f64, d: usize) -> String {
    let s = format!("{:.*e}", d, x);
    // Rust writes 1.2e4; C writes 1.2e+04
    let (m, e) = s.split_once('e').unwrap();
    let ev: i32 = e.parse().unwrap();
    format!("{m}e{}{:02}", if ev < 0 { '-' } else { '+' }, ev.abs())
}

/// A formatted number: R's format_numbers returns numbers (rounded) or strings (formatC); sapply makes them all
/// strings when any is one.
#[derive(Clone, Debug)]
enum Fmt {
    Num(f64),
    Str(String),
}

impl Fmt {
    fn text(&self) -> String {
        match self {
            Fmt::Num(x) => r_chr(*x),
            Fmt::Str(s) => s.clone(),
        }
    }
}

/// NUMBER_FORMATTING_SIGNATURE_QUANTILES
fn format_sig_quantile(x: f64) -> Fmt {
    if (10000.0..f64::INFINITY).contains(&x) {
        Fmt::Str(format_c_e(x, 0))
    } else if (10.0..10000.0).contains(&x) {
        Fmt::Num(r_round(x, 0))
    } else if (1.0..10.0).contains(&x) {
        Fmt::Num(r_round(x, 1))
    } else if (-1.0..1.0).contains(&x) {
        Fmt::Num(r_round(x, 2))
    } else if (-10.0..-1.0).contains(&x) {
        Fmt::Num(r_round(x, 1))
    } else if (-10000.0..-10.0).contains(&x) {
        Fmt::Num(r_round(x, 0))
    } else if x < -10000.0 {
        // formatC(-Inf, format = "e") is "-Inf"
        if x == f64::NEG_INFINITY { Fmt::Str("-Inf".into()) } else { Fmt::Str(format_c_e(x, 0)) }
    } else {
        Fmt::Num(x)
    }
}

/// NUMBER_FORMATTING_ODDS_RATIOS
fn format_odds(x: f64) -> Fmt {
    if (10000.0..f64::INFINITY).contains(&x) {
        Fmt::Str(format_c_e(x, 1))
    } else if (10.0..10000.0).contains(&x) {
        Fmt::Num(r_round(x, 0))
    } else if (1.0..10.0).contains(&x) {
        Fmt::Num(r_round(x, 1))
    } else if (0.01..1.0).contains(&x) {
        Fmt::Num(r_signif1(x))
    } else if x < 0.01 && x > f64::NEG_INFINITY {
        Fmt::Str(format_c_e(x, 0))
    } else {
        Fmt::Num(x)
    }
}

/// stringr::str_wrap(x, width, whitespace_only = FALSE) as stringi 1.8 wraps (ICU line breaks: after spaces, which
/// are dropped, and after '/' and '-' before a letter; lines chosen by minimum raggedness, cost exponent 2, the
/// last line free; a word longer than the width on a line of its own).
pub fn str_wrap(s: &str, width: usize) -> String {
    // tokens: (text, separator after it: " " when the break is a space, "" after '/' or '-')
    let chars: Vec<char> = s.chars().collect();
    let mut tokens: Vec<(String, &str)> = Vec::new();
    let mut cur = String::new();
    let mut i = 0;
    while i < chars.len() {
        let ch = chars[i];
        if ch == ' ' {
            if !cur.is_empty() {
                tokens.push((std::mem::take(&mut cur), " "));
            }
            while i < chars.len() && chars[i] == ' ' {
                i += 1;
            }
            continue;
        }
        cur.push(ch);
        if (ch == '/' || ch == '-') && chars.get(i + 1).is_some_and(|c| c.is_alphabetic()) {
            tokens.push((std::mem::take(&mut cur), ""));
        }
        i += 1;
    }
    if !cur.is_empty() {
        tokens.push((cur, ""));
    }
    let n = tokens.len();
    if n == 0 {
        return String::new();
    }
    let len = |a: usize, b: usize| -> usize {
        // tokens a..=b on one line
        let mut l = 0;
        for k in a..=b {
            l += tokens[k].0.chars().count();
            if k < b {
                l += tokens[k].1.len();
            }
        }
        l
    };
    // f[j]: the least cost of lines ending with token j - 1
    let mut f = vec![f64::INFINITY; n + 1];
    let mut back = vec![0usize; n + 1];
    f[0] = 0.0;
    for j in 1..=n {
        for i in (1..=j).rev() {
            let l = len(i - 1, j - 1);
            if l > width && i < j {
                break;
            }
            let cost = if j == n { 0.0 } else if l > width { 0.0 } else { ((width - l) as f64).powi(2) };
            if f[i - 1] + cost < f[j] {
                f[j] = f[i - 1] + cost;
                back[j] = i - 1;
            }
        }
    }
    let mut lines = Vec::new();
    let mut j = n;
    while j > 0 {
        let i = back[j];
        let mut line = String::new();
        for k in i..j {
            line.push_str(&tokens[k].0);
            if k + 1 < j {
                line.push_str(tokens[k].1);
            }
        }
        lines.push(line);
        j = i;
    }
    lines.reverse();
    lines.join("\n")
}

// -------------------------------------------------------------------------------------------------------------
// the heat map data

/// One tile of a heat map.
#[derive(Clone, Debug)]
struct Cell {
    row_label: String,
    row_group: String,
    cancer_type: usize,
    /// the fill colour (from the value and the scale)
    fill: Rgba,
    label: String,
    nudge_x: f64,
    outline: Option<Rgba>,
}

enum XAxis {
    Top,
    None,
}

/// A heat map to draw (plot_heatmap's arguments and its scales and labels).
struct Heatmap {
    cells: Vec<Cell>,
    x_axis: XAxis,
    show_first_row_hline: bool,
    y_labels: Box<dyn Fn(&str) -> String>,
    title: String,
    subtitle: Option<String>,
    legend: LegendSpec,
    legend_just: (f64, f64),
}

fn col(name: &str) -> Rgba {
    r_colour(name).unwrap_or_else(|| panic!("colour {name}"))
}

const PROB_DECIMALS: i32 = 2;

fn revalue_clf(name: &str) -> String {
    match name {
        "combined" => "COMBINED",
        "dna_combined" => "DNA COMBINED",
        "gen_pos" => "GENOMIC POSITION",
        "snv96" => "SNV96",
        "event" => "EVENT",
        "rna_combined" => "RNA COMBINED",
        "gene_exp" => "GENE EXPRESSION",
        "alt_sj" => "ALT SJ",
        other => other,
    }
    .to_string()
}

/// scale_fill_gradient(low, high, limits = c(0, 1), na.value = "grey95"): the colour of a value.
fn gradient_fill(ramp: &GradientN, v: Option<f64>, limits: (f64, f64), na: Rgba) -> Rgba {
    match v {
        Some(x) if x.is_finite() && x >= limits.0 && x <= limits.1 => ramp.at((x - limits.0) / (limits.1 - limits.0)).unwrap_or(na),
        _ => na,
    }
}

fn snv_count(rows: &[&VisRow]) -> Result<f64> {
    match rows.iter().find(|r| r.feat_name.ends_with("snv_count")).and_then(|r| r.feat_value) {
        Some(v) => Ok(v),
        None => bail!("`snv_count` was not found as a feature in `VIS_DATA`"),
    }
}

/// get_plot_data_probs and plot_probs
fn probs(rows: &[&VisRow], ct: &CancerTypes) -> Result<Heatmap> {
    let data: Vec<&VisRow> = rows.iter().copied().filter(|r| r.data_type == "prob").collect();
    // classifiers without data are removed
    let mut has_data: HashMap<&str, bool> = HashMap::new();
    for r in &data {
        *has_data.entry(r.clf_name.as_str()).or_insert(false) |= r.data_value.is_some();
    }
    let data: Vec<&VisRow> = data.into_iter().filter(|r| has_data[r.clf_name.as_str()]).collect();
    struct P {
        row_label: String,
        row_group: String,
        cancer_type: usize,
        value: Option<f64>,
        label: String,
        nudge: f64,
    }
    let mut plain: Vec<P> = data
        .iter()
        .map(|r| P {
            row_label: revalue_clf(&r.clf_name),
            row_group: r.clf_group.to_uppercase(),
            cancer_type: ct.index[&r.cancer_type],
            value: r.data_value,
            label: r.data_value.map(|v| r_chr(r_round(v, PROB_DECIMALS))).unwrap_or_default(),
            nudge: 0.0,
        })
        .collect();
    // supertype sums of the combined classifiers: split by interaction(row_label, cancer_supertype), the row label
    // (sorted levels) varying fastest
    let mut combined_labels: Vec<String> = plain.iter().filter(|p| p.row_label.to_lowercase().contains("combined")).map(|p| p.row_label.clone()).collect();
    combined_labels.sort();
    combined_labels.dedup();
    let mut groups: Vec<P> = Vec::new();
    for sup in 0..ct.supertypes.len() {
        for lab in &combined_labels {
            let idx: Vec<usize> = (0..plain.len()).filter(|&i| plain[i].row_label == *lab && ct.supertype[plain[i].cancer_type] == sup).collect();
            if idx.is_empty() {
                continue;
            }
            if idx.len() == 1 {
                let p = &plain[idx[0]];
                groups.push(P { row_label: format!("[GROUP] {}", p.row_label), row_group: p.row_group.clone(), cancer_type: p.cancer_type, value: p.value, label: p.label.clone(), nudge: 0.0 });
                continue;
            }
            let sum = r_round(idx.iter().filter_map(|&i| plain[i].value).sum::<f64>(), PROB_DECIMALS);
            let middle = idx.len().div_ceil(2);
            let nudge = if idx.len() % 2 == 0 { 0.5 } else { 0.0 };
            for (k, &i) in idx.iter().enumerate() {
                let p = &plain[i];
                groups.push(P {
                    row_label: format!("[GROUP] {}", p.row_label),
                    row_group: p.row_group.clone(),
                    cancer_type: p.cancer_type,
                    value: Some(sum),
                    label: if k + 1 == middle { r_chr(sum) } else { String::new() },
                    nudge,
                });
            }
        }
    }
    let group_labels: Vec<String> = groups.iter().map(|g| g.row_label.clone()).collect();
    let ramp = GradientN::new(&[col("white"), col("mediumseagreen")], None);
    let na = col("grey95");
    let mut cells = Vec::new();
    for p in groups.into_iter().chain(plain.drain(..)) {
        let outline = if group_labels.contains(&p.row_label) { None } else { Some(col("darkgrey")) };
        cells.push(Cell { fill: gradient_fill(&ramp, p.value, (0.0, 1.0), na), row_label: p.row_label, row_group: p.row_group, cancer_type: p.cancer_type, label: p.label, nudge_x: p.nudge, outline });
    }
    let snv = snv_count(rows)?;
    let warning = if snv < 50.0 { " | WARNING: Sample has <50 SNVs. Predictions may be unreliable".to_string() } else { String::new() };
    Ok(Heatmap {
        cells,
        x_axis: XAxis::Top,
        show_first_row_hline: true,
        y_labels: Box::new(|s| s.to_string()),
        title: format!("Probabilities by classifier{warning}"),
        subtitle: Some("Cancer group (strips) and subtype (label)".into()),
        legend: {
            let (b, l) = unit_labels();
            bar_legend("Probability", &|v| gradient_fill(&ramp, Some(v), (0.0, 1.0), na), (0.0, 1.0), b, l, 3.5)
        },
        legend_just: (0.0, 1.0),
    })
}

/// A colour bar legend as the script asks for it (frame and ticks black, the frame 0.25 mm, barheight in lines):
/// the colours of 300 values from limit to limit (transformed space) and the breaks' places on the bar.
fn bar_legend(title: &str, map: &dyn Fn(f64) -> Rgba, limits: (f64, f64), breaks: Vec<f64>, labels: Vec<String>, barheight: f64) -> LegendSpec {
    let nbin = 300;
    let (l0, l1) = limits;
    let by = (l1 - l0) / (nbin - 1) as f64;
    let colours = (0..nbin).map(|k| map(l0 + k as f64 * by)).collect();
    let breaks = breaks.iter().map(|&b| (0.5 + (b - l0) / (l1 - l0) * (nbin - 1) as f64) / nbin as f64).collect();
    LegendSpec::Bar { title: title.into(), colours, breaks, labels, bar_width: 0.0, bar_height: barheight, frame: Some((Rgba::BLACK, 0.25)), ticks: Some((Rgba::BLACK, 0.5 / hmf_plot::ggp::PT)), tick_len: 0.0 }
}

fn unit_labels() -> (Vec<f64>, Vec<String>) {
    (vec![0.0, 0.25, 0.5, 0.75, 1.0], ["0.00", "0.25", "0.50", "0.75", "1.00"].map(String::from).to_vec())
}

/// plot_cv_performance
fn cv_performance(rows: &[&VisRow], ct: &CancerTypes) -> Heatmap {
    let data: Vec<&VisRow> = rows.iter().copied().filter(|r| r.data_type == "cv_performance").collect();
    let max_n = data.iter().filter(|r| r.feat_name == "n_total").filter_map(|r| r.data_value).fold(f64::NEG_INFINITY, f64::max);
    let ramp = GradientN::new(&[col("white"), col("#fffedc"), col("#76c99b")], None);
    let na = col("grey95");
    let cells = data
        .iter()
        .map(|r| {
            let n_row = r.feat_name == "n_total";
            let label = match r.data_value {
                Some(v) if n_row => r_chr(v),
                Some(v) => r_chr(r_round(v, 2)),
                None => "NA".into(),
            };
            let value = if n_row { r.data_value.map(|v| v / max_n) } else { r.data_value };
            Cell { row_label: r.feat_name.clone(), row_group: "training set".into(), cancer_type: ct.index[&r.cancer_type], fill: gradient_fill(&ramp, value, (0.0, 1.0), na), label, nudge_x: 0.0, outline: Some(col("lightgrey")) }
        })
        .collect();
    Heatmap {
        cells,
        x_axis: XAxis::None,
        show_first_row_hline: true,
        y_labels: Box::new(|s| match s {
            "n_total" => "Total no. of samples".into(),
            "recall" => "Recall (prop. of total correct)".into(),
            "precision" => "Precision (prop. of predicted correct)".into(),
            other => other.into(),
        }),
        title: "DNA COMBINED: Training set performance".into(),
        subtitle: None,
        legend: {
            let (b, l) = unit_labels();
            bar_legend("Metric value", &|v| gradient_fill(&ramp, Some(v), (0.0, 1.0), na), (0.0, 1.0), b, l, 3.5)
        },
        legend_just: (0.0, 0.0),
    }
}

/// plot_signatures
fn signatures(rows: &[&VisRow], ct: &CancerTypes) -> Result<Heatmap> {
    let data: Vec<&VisRow> = rows.iter().copied().filter(|r| r.data_type == "sig_quantile").collect();
    let snv = snv_count(rows)?;
    let labels = ["<0 (below expected range)", "0", ">0 - 0.95 (in expected range)", "0.95 - 1.2 (above expected range)", ">1.2 (well above expected range)"];
    let colours = ["lightcyan2", "grey95", "#C9E7CD", "#f5dfdf", "#e8b6b6"];
    // sapply: all strings when one is
    let fmts: Vec<Fmt> = data.iter().map(|r| format_sig_quantile(r.data_value.unwrap_or(f64::NAN))).collect();
    let mut present = [false; 5];
    let mut any_na = false;
    let mut cells = Vec::new();
    for (r, f) in data.iter().zip(&fmts) {
        let v = r.data_value.unwrap_or(f64::NAN);
        // cut(breaks = c(-Inf, -0.004, 0.004, 0.95, 1.2, Inf)), right-closed: -Inf lies in no interval (NA), so R
        // draws it in the colour of "0" with an "NA" key; corrected, it is below the expected range (CU9)
        let class = if v.is_nan() || (v == f64::NEG_INFINITY && crate::compat::java("CU9")) {
            None
        } else if v <= -0.004 {
            Some(0)
        } else if v <= 0.004 {
            Some(1)
        } else if v <= 0.95 {
            Some(2)
        } else if v <= 1.2 {
            Some(3)
        } else {
            Some(4)
        };
        match class {
            Some(c) => present[c] = true,
            None => any_na = true,
        }
        // scale_fill_manual with named values: "0" is not a name, so it (and NA) take na.value
        let fill = match class {
            Some(c) if c != 1 => col(colours[c]),
            _ => col("grey95"),
        };
        let fv = r.feat_value.unwrap_or(f64::NAN);
        let perc = r_round(fv / snv * 100.0, 2);
        let feat_value = r_round(fv, 0);
        let name = r.feat_name.replacen("SBS", "Sig", 1).replace('_', " ");
        let row_label = format!("{} = {} ({}%)", name, r_chr(feat_value), r_chr(perc));
        cells.push(Cell { row_label, row_group: "signature".into(), cancer_type: ct.index[&r.cancer_type], fill, label: f.text(), nudge_x: 0.0, outline: Some(col("grey50")) });
    }
    let mut entries: Vec<(String, KeyGlyph)> = (0..5)
        .filter(|&c| present[c] && c != 1)
        .map(|c| (labels[c].to_string(), KeyGlyph::Polygon { fill: col(colours[c]), colour: Some(Rgba::BLACK), linewidth: 0.25 }))
        .collect();
    // na.translate: missing values get a key of their own, in na.value
    if any_na {
        entries.push(("NA".to_string(), KeyGlyph::Polygon { fill: col("grey95"), colour: Some(Rgba::BLACK), linewidth: 0.25 }));
    }
    Ok(Heatmap {
        cells,
        x_axis: XAxis::None,
        show_first_row_hline: false,
        y_labels: Box::new(|s| s.to_string()),
        title: "SNV96: Mutational signatures".into(),
        subtitle: None,
        legend: LegendSpec::Keys { title: "Quantile in subtype cohort".into(), entries, key_width: 0.0, key_height: 0.0 },
        legend_just: (0.0, 1.0),
    })
}

/// plot_feat_contrib
fn feat_contrib(rows: &[&VisRow], ct: &CancerTypes) -> Heatmap {
    let data: Vec<&VisRow> = rows.iter().copied().filter(|r| r.data_type == "feat_contrib").collect();
    let (types, names): (Vec<String>, Vec<String>) = data
        .iter()
        .map(|r| match r.feat_name.split_once('.') {
            Some((a, b)) => (a.to_string(), b.to_string()),
            None => (r.feat_name.clone(), String::new()),
        })
        .unzip();
    // odds ratio labels
    let fmts: Vec<Fmt> = data.iter().map(|r| format_odds(r.data_value.unwrap_or(f64::NAN))).collect();
    let any_str = fmts.iter().any(|f| matches!(f, Fmt::Str(_)));
    let odds: Vec<String> = fmts
        .iter()
        .map(|f| {
            let s = f.text();
            // odds_ratios[odds_ratios == 1] <- "": a string compares with "1", a number with 1
            let is_one = match f {
                Fmt::Num(x) => *x == 1.0,
                Fmt::Str(_) => false,
            } || (any_str && s == "1");
            if is_one { String::new() } else { s.replacen("-0", "-", 1).replacen("+0", "", 1) }
        })
        .collect();
    // feature names
    let names: Vec<String> = types
        .iter()
        .zip(&names)
        .map(|(t, n)| match t.as_str() {
            "driver" => n.replacen('.', " ", 1),
            "fusion" => n.replacen('_', "-", 1),
            "trait" | "tmb" | "sv" => n.replace('_', " "),
            _ => n.clone(),
        })
        .collect();
    // feature values: round to 1 decimal; "0", "1" and "absent" assigned make the vector character
    let raw: Vec<f64> = data.iter().map(|r| r.feat_value.unwrap_or(f64::NAN)).collect();
    let as_char = raw.iter().any(|&v| v == 0.0 || v == 1.0 || v == -0.00000001);
    let values: Vec<String> = raw
        .iter()
        .zip(&types)
        .map(|(&v, t)| {
            let s = if v == 0.0 {
                "0".to_string()
            } else if v == 1.0 {
                "1".to_string()
            } else if v == -0.00000001 {
                "absent".to_string()
            } else {
                r_chr(r_round(v, 1))
            };
            if t == "trait" {
                // feat_values == 1: a string compares with "1", a number with 1; TRUE / FALSE
                let one = if as_char { s == "1" } else { r_round(v, 1) == 1.0 };
                if as_char {
                    (if one { "TRUE" } else { "FALSE" }).to_string()
                } else {
                    (if one { "1" } else { "0" }).to_string()
                }
            } else {
                s
            }
        })
        .collect();
    // the colour scale: log odds, white at 1 (rescale_mid), indianred below, mediumseagreen above
    let logs: Vec<f64> = data.iter().map(|r| r.data_value.unwrap_or(f64::NAN).ln()).collect();
    let (lmin, lmax) = logs.iter().filter(|v| v.is_finite()).fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
    let ramp = GradientN::new(&[col("indianred"), col("white"), col("mediumseagreen")], None);
    let span = (lmin - 0.0).abs().max((lmax - 0.0).abs());
    let fill_of = |l: f64| -> Rgba {
        if !l.is_finite() || l < lmin || l > lmax {
            return col("grey50");
        }
        let x = if span == 0.0 { 0.5 } else { (l - 0.0) / span / 2.0 + 0.5 };
        ramp.at(x).unwrap_or(col("grey50"))
    };
    let cells = data
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let row_label = if types[i] == "driver" { format!("{} ({})", names[i], values[i]) } else { format!("{} = {}", names[i], values[i]) };
            let row_group = match types[i].as_str() {
                "driver" => "driver (DL)",
                "virus" => "virus (DL)",
                "fusion" => "fusion (DL)",
                t => t,
            }
            .to_string();
            Cell { row_label, row_group, cancer_type: ct.index[&r.cancer_type], fill: fill_of(logs[i]), label: odds[i].clone(), nudge_x: 0.0, outline: Some(col("lightgrey")) }
        })
        .collect();
    // legend breaks: 10^seq(-10, 10, by = 1 or 2) within the limits; labels sprintf("%g")
    let dmax = data.iter().filter_map(|r| r.data_value).fold(f64::NEG_INFINITY, f64::max);
    let dmin = data.iter().filter_map(|r| r.data_value).fold(f64::INFINITY, f64::min);
    let range = dmax.log10() - dmin.log10();
    let by = if range < 5.0 { 1 } else { 2 };
    let breaks: Vec<f64> = (0..).map(|k| -10 + k * by).take_while(|&e| e <= 10).map(|e| 10f64.powi(e)).filter(|b| b.ln() >= lmin - 1e-12 && b.ln() <= lmax + 1e-12).collect();
    let labels: Vec<String> = breaks.iter().map(|&b| sprintf_g(b)).collect();
    let log_breaks: Vec<f64> = breaks.iter().map(|b| b.ln()).collect();
    let legend = bar_legend("Odds ratio\n(subtype / not subtype)", &fill_of, (lmin, lmax), log_breaks, labels, 5.0);
    Heatmap {
        cells,
        x_axis: XAxis::None,
        show_first_row_hline: false,
        y_labels: Box::new(|s| s.to_string()),
        title: "EVENT: Feature contributions".into(),
        subtitle: None,
        legend,
        legend_just: (0.0, 1.0),
    }
}

/// sprintf("%g", x) for the powers of ten used as breaks
fn sprintf_g(x: f64) -> String {
    let e = x.log10().round() as i32;
    if (-4..6).contains(&e) {
        let s = format!("{x}");
        s.trim_end_matches(".0").to_string()
    } else {
        format!("1e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    }
}

// -------------------------------------------------------------------------------------------------------------
// drawing

/// theme_bw() with the script's theme, for one device.
struct Theme {
    ctx: Ctx,
}

impl Theme {
    fn half_line(&self) -> f64 {
        5.5
    }
    fn axis_text(&self) -> TextEl {
        let mut t = TextEl::new(8.8);
        t.colour = col("grey30");
        t
    }
    fn strip_text(&self) -> TextEl {
        let mut t = TextEl::new(8.8);
        t.colour = col("grey10");
        let m = self.ctx.pt(0.8 * self.half_line());
        t.margin = [m, m, m, m];
        t
    }
    fn ticks(&self) -> LineEl {
        LineEl::solid(col("grey20"), 0.5)
    }
    fn tick_len(&self) -> f64 {
        self.ctx.pt(self.half_line() / 2.0)
    }
}

/// plot_heatmap: the plot's table for patchwork.
fn heatmap_table(th: &Theme, ct: &CancerTypes, hm: Heatmap) -> (Gtable, usize, usize) {
    let ctx = th.ctx;
    // factor levels in order of appearance
    let mut row_labels: Vec<String> = Vec::new();
    let mut row_groups: Vec<String> = Vec::new();
    for c in &hm.cells {
        if !row_labels.contains(&c.row_label) {
            row_labels.push(c.row_label.clone());
        }
        if !row_groups.contains(&c.row_group) {
            row_groups.push(c.row_group.clone());
        }
    }
    // the facet columns: supertypes present; the cancer types present in each (factor order)
    let mut types_present: Vec<usize> = hm.cells.iter().map(|c| c.cancer_type).collect();
    types_present.sort();
    types_present.dedup();
    let mut sups: Vec<usize> = types_present.iter().map(|&t| ct.supertype[t]).collect();
    sups.sort();
    sups.dedup();
    let col_types: Vec<Vec<usize>> = sups.iter().map(|&s| types_present.iter().copied().filter(|&t| ct.supertype[t] == s).collect()).collect();
    // the facet rows: groups; the row labels present in each (levels order), reversed by limits = rev
    let row_levels: Vec<Vec<String>> = row_groups.iter().map(|g| row_labels.iter().filter(|l| hm.cells.iter().any(|c| &c.row_group == g && &c.row_label == *l)).cloned().collect()).collect();
    let n_rows = row_labels.len();
    let n_cols = types_present.len();

    let ypos = |gi: usize, label: &str| -> f64 {
        let levels = &row_levels[gi];
        let k = levels.iter().position(|l| l == label).unwrap();
        (levels.len() - k) as f64
    };
    let xpos = |ci: usize, t: usize| -> f64 { (col_types[ci].iter().position(|&x| x == t).unwrap() + 1) as f64 };

    // panels
    let mut panels: Vec<Vec<Grob>> = Vec::new();
    for (gi, g) in row_groups.iter().enumerate() {
        let mut row = Vec::new();
        for (ci, &s) in sups.iter().enumerate() {
            let cells: Vec<Cell> = hm.cells.iter().filter(|c| &c.row_group == g && ct.supertype[c.cancer_type] == s).cloned().collect();
            let nx = col_types[ci].len() as f64;
            let ny = row_levels[gi].len() as f64;
            let pos: Vec<(f64, f64)> = cells.iter().map(|c| (xpos(ci, c.cancer_type), ypos(gi, &c.row_label))).collect();
            // geom_hline: the script gives a row per (cancer type, row label), but GeomHline draws unique0() of its
            // aesthetics, so one line per panel, at the first row's bottom
            let n_hlines = usize::from(hm.show_first_row_hline && !cells.is_empty());
            row.push(Grob::draw_fn(move |cv, b| {
                let p = Panel { b, xr: (0.5, nx + 0.5), yr: (0.5, ny + 0.5), flip: false };
                // panel.background: white
                hmf_plot::ggp::rect(cv, &ctx, b.x0, b.y0, b.x1, b.y1, Some(Rgba::WHITE), None, Vec::new());
                for (c, &(x, y)) in cells.iter().zip(&pos) {
                    geom::data_rect(cv, &ctx, &p, x - 0.5, x + 0.5, y - 0.5, y + 0.5, Some(c.fill), c.outline.map(|o| (o, 0.1)));
                }
                for (c, &(x, y)) in cells.iter().zip(&pos) {
                    geom::text(cv, &ctx, &p, x + c.nudge_x, y, &c.label, 2.7, Rgba::BLACK, 0.5, 0.5);
                }
                for _ in 0..n_hlines {
                    geom::hline(cv, &ctx, &p, ny - 0.5, Rgba::BLACK, 0.5, "A1");
                }
                // panel.border: 0.75 mm grey20, unfilled
                hmf_plot::ggp::rect(cv, &ctx, b.x0, b.y0, b.x1, b.y1, None, Some((col("grey20"), 0.75 * hmf_plot::ggp::PT)), Vec::new());
            }));
        }
        panels.push(row);
    }

    // axes: x on top (labels rotated), y on the left
    let top_axis = matches!(hm.x_axis, XAxis::Top);
    let mut x_text = th.axis_text();
    x_text.angle = 90.0;
    x_text.hjust = 0.0;
    x_text.vjust = 0.5;
    x_text.margin = [0.0, 0.0, ctx.pt(0.8 * th.half_line() / 2.0), 0.0];
    let axes_top: Vec<guides::Axis> = col_types
        .iter()
        .map(|types| {
            let n = types.len() as f64;
            guides::axis(
                ctx,
                AxisSpec {
                    side: Side::Top,
                    breaks: (0..types.len()).map(|k| (k as f64 + 0.5) / n).collect(),
                    labels: types.iter().map(|&t| str_wrap(&ct.subtype[t], 20)).collect(),
                    text: top_axis.then(|| x_text.clone()),
                    ticks: top_axis.then(|| th.ticks()),
                    tick_len: th.tick_len(),
                    major_ticks: true,
                    extra_ticks: Vec::new(),
                },
            )
        })
        .collect();
    let axes_bottom: Vec<guides::Axis> = col_types.iter().map(|_| guides::Axis::zero()).collect();
    let mut y_text = th.axis_text();
    y_text.hjust = 1.0;
    y_text.margin = [0.0, ctx.pt(0.8 * th.half_line() / 2.0), 0.0, 0.0];
    let axes_left: Vec<guides::Axis> = row_levels
        .iter()
        .map(|levels| {
            let n = levels.len() as f64;
            // limits = rev: the last level at position 1
            let rev: Vec<&String> = levels.iter().rev().collect();
            guides::axis(
                ctx,
                AxisSpec {
                    side: Side::Left,
                    breaks: (0..rev.len()).map(|k| (k as f64 + 0.5) / n).collect(),
                    labels: rev.iter().map(|l| (hm.y_labels)(l)).collect(),
                    text: Some(y_text.clone()),
                    ticks: Some(th.ticks()),
                    tick_len: th.tick_len(),
                    major_ticks: true,
                    extra_ticks: Vec::new(),
                },
            )
        })
        .collect();
    let axes_right: Vec<guides::Axis> = row_levels.iter().map(|_| guides::Axis::zero()).collect();

    // strips: themed (ggh4x strip_themed) with an x axis, vanilla otherwise
    let strip_bg_lw = 0.5;
    let strips_x: Vec<Strip> = sups
        .iter()
        .map(|&s| {
            if !top_axis {
                return Strip { grob: Grob::Zero, size: 0.0 };
            }
            let has_sub = ct.subtypes_per_supertype[s] > 1;
            let mut el = th.strip_text();
            el.colour = if has_sub { Rgba::BLACK } else { facet::INVISIBLE };
            let bg = RectEl { fill: Some(if has_sub { col("lightgrey") } else { Rgba::WHITE }), colour: Some(if has_sub { Rgba::BLACK } else { facet::INVISIBLE }), linewidth: strip_bg_lw };
            facet::strip(ctx, &ct.supertypes[s], Some(&el), Some(bg), false)
        })
        .collect();
    let strips_y: Vec<Strip> = row_groups
        .iter()
        .map(|g| {
            let ws = " ".repeat(21);
            let label = format!("{ws}\n{g}\n{ws}");
            let mut el = th.strip_text();
            el.hjust = 1.0;
            let bg = if top_axis {
                RectEl { fill: Some(col("lightgrey")), colour: Some(Rgba::BLACK), linewidth: strip_bg_lw }
            } else {
                RectEl { fill: Some(col("grey85")), colour: Some(col("grey20")), linewidth: strip_bg_lw }
            };
            facet::strip(ctx, &label, Some(&el), Some(bg), true)
        })
        .collect();

    let grid = Grid {
        widths: col_types.iter().map(|t| t.len() as f64).collect(),
        heights: row_levels.iter().map(|l| l.len() as f64).collect(),
        panels,
        spacing_x: ctx.pt(3.0),
        spacing_y: ctx.pt(3.0),
        clip: true,
        axes_top,
        axes_bottom,
        axes_left,
        axes_right,
        strips_x: Some(strips_x),
        strips_y: Some(strips_y),
        switch_x: false,
        switch_y: true,
        outside_x: true,
        outside_y: true,
        padding: ctx.pt(th.half_line() / 2.0),
        pad_always: true,
    };
    let table = grid.build();

    // the legend
    let lines = |n: f64| ctx.dev.bigpt(n * 12.0 * 1.2);
    let key = lines(1.2);
    let legend_spec = match hm.legend {
        LegendSpec::Bar { title, colours, breaks, labels, bar_height, frame, ticks, .. } => {
            LegendSpec::Bar { title, colours, breaks, labels, bar_width: key, bar_height: lines(bar_height), frame, ticks, tick_len: 0.2 * key }
        }
        LegendSpec::Keys { title, entries, .. } => LegendSpec::Keys { title, entries, key_width: key, key_height: key },
    };
    let mut ltext = TextEl::new(8.8);
    ltext.hjust = 0.0;
    let mut ltitle = TextEl::new(11.0);
    ltitle.hjust = 0.0;
    let lth = LegendTheme {
        title: ltitle,
        text: ltext,
        key_spacing: ctx.pt(5.5),
        margin: [0.0; 4],
        background: Some(RectEl { fill: Some(Rgba::WHITE), colour: None, linewidth: 0.5 }),
        key_background: Some(RectEl { fill: Some(Rgba::WHITE), colour: None, linewidth: 0.5 }),
    };
    let leg = guides::legend(ctx, &legend_spec, &lth);
    let gbox = guides::guide_box(ctx, GuideBox { legends: vec![leg], horizontal: false, just: hm.legend_just, box_just: (0.0, 1.0), spacing: ctx.pt(11.0), box_margin: [0.0; 4] });

    let mut title_el = TextEl::new(13.2);
    title_el.face = Face::Bold;
    title_el.hjust = 0.0;
    title_el.vjust = 1.0;
    title_el.margin = [0.0, 0.0, ctx.pt(5.5), 0.0];
    let mut sub_el = TextEl::new(11.0);
    sub_el.hjust = 0.0;
    sub_el.vjust = 1.0;
    sub_el.margin = [0.0, 0.0, ctx.pt(5.5), 0.0];
    let parts = PlotParts {
        panels: table,
        xlab_t: None,
        xlab_b: None,
        ylab_l: None,
        ylab_l_custom: None,
        ylab_r: None,
        legend_right: Some(gbox),
        legend_left: None,
        legend_top: None,
        legend_bottom: None,
        legend_box_spacing: ctx.pt(11.0),
        title: Some(TitleGrob::new(&hm.title, &title_el, true, true)),
        subtitle: hm.subtitle.as_ref().map(|s| TitleGrob::new(s, &sub_el, true, true)),
        caption: None,
        title_on_plot: false,
        margin: [0.0, 0.0, ctx.pt(20.0), 0.0],
        background: Some(RectEl { fill: Some(Rgba::WHITE), colour: Some(Rgba::WHITE), linewidth: 0.5 }),
    };
    let _ = Bx::new;
    (plot_table(ctx, parts), n_rows, n_cols)
}

/// plot_one_sample: the four heat maps stacked; returns the patchwork's table, the heights (rows) and widths (types).
fn one_sample(th: &Theme, all: &[VisRow], ct: &CancerTypes, sample: &str, annotation_title: Option<&str>) -> Result<(Gtable, Vec<usize>, Vec<usize>)> {
    let rows: Vec<&VisRow> = all.iter().filter(|r| r.sample_id == sample || r.data_type == "cv_performance").collect();
    let p1 = probs(&rows, ct)?;
    let p2 = signatures(&rows, ct)?;
    let p3 = feat_contrib(&rows, ct);
    if !rows.iter().any(|r| r.data_type == "cv_performance") {
        bail!("CV performance data is missing from vis data");
    }
    let p4 = cv_performance(&rows, ct);
    let mut patches = Vec::new();
    let mut heights = Vec::new();
    let mut widths = Vec::new();
    for hm in [p1, p2, p3, p4] {
        let (gt, nr, nc) = heatmap_table(th, ct, hm);
        if nr > 0 {
            heights.push(nr);
        }
        widths.push(nc);
        patches.push(patchwork::plot(gt));
    }
    let layout = Layout::column(patches.len(), heights.iter().map(|&h| Some(h as f64)).collect());
    let gt = patchwork::build(patches, &layout);
    // plot_annotation: theme_grey's margins (5.5 pt), and the sample's title on cohort pages
    let ctx = th.ctx;
    let m = ctx.pt(5.5);
    let title = annotation_title.map(|t| {
        let mut el = TextEl::new(13.2);
        el.face = Face::Bold;
        el.hjust = 0.0;
        el.vjust = 1.0;
        el.margin = [0.0, 0.0, ctx.pt(5.5), 0.0];
        let g = TitleGrob::new(t, &el, true, true);
        let h = g.height(&ctx).abs();
        (g.grob(ctx), h)
    });
    let table = patchwork::annotate(gt, Annotation { margin: [m, m, m, m], title, subtitle: None, background: Some(rect_grob(ctx, RectEl { fill: Some(Rgba::WHITE), colour: Some(Rgba::WHITE), linewidth: 0.5 })) });
    Ok((table, heights, widths))
}

/// Draws the plot of the vis data file: a PNG for one sample, a PDF of a page per sample for several.
pub fn draw(vis_path: &str, plot_path: &str) -> Result<()> {
    let rows = read_vis(vis_path)?;
    let ct = CancerTypes::new(&rows);
    let mut samples: Vec<String> = Vec::new();
    for r in rows.iter().filter(|r| r.data_type == "prob") {
        if !samples.contains(&r.sample_id) {
            samples.push(r.sample_id.clone());
        }
    }
    if samples.len() == 1 {
        if !plot_path.ends_with(".png") && !plot_path.ends_with(".pdf") {
            bail!("Plot output file must have .pdf or .png as the extension");
        }
        let pdf = plot_path.ends_with(".pdf");
        let dev = if pdf { Device::Pdf } else { Device::Png { dpi: 300.0 } };
        let th = Theme { ctx: Ctx { dev, v: Version::V351 } };
        let (table, heights, widths) = one_sample(&th, &rows, &ct, &samples[0], None)?;
        // height <- sum(heights) * 0.25 + 3; width <- max(widths) * 0.35 + 6 (inches)
        let h_in = heights.iter().sum::<usize>() as f64 * 0.25 + 3.0;
        let w_in = *widths.iter().max().unwrap_or(&0) as f64 * 0.35 + 6.0;
        if pdf {
            return write_pdf(plot_path, &[table], w_in, h_in);
        }
        // png(width, height, units = "in", res = 300): R takes the integer part of res * size
        let (wpx, hpx) = ((300.0 * w_in) as u32, (300.0 * h_in) as u32);
        let mut fig = Figure::new(wpx, hpx, Rgba::WHITE);
        fig.meta = PngMeta { dpi: Some(300.0), transparent_key: None };
        let mut c = Canvas::new(&mut fig);
        table.draw(&mut c, Bx::new(0.0, 0.0, wpx as f64, hpx as f64));
        dump_boxes(&table, wpx as f64, hpx as f64, 1);
        fig.write_png(plot_path).with_context(|| format!("writing {plot_path}"))?;
        return Ok(());
    }
    if !plot_path.ends_with(".pdf") {
        bail!("Plot output file extension must be .pdf for multi-sample plotting");
    }
    let th = Theme { ctx: Ctx { dev: Device::Pdf, v: Version::V351 } };
    let mut pages = Vec::new();
    for s in &samples {
        let (table, _, _) = one_sample(&th, &rows, &ct, s, Some(&format!("Sample: {s}")))?;
        pages.push(table);
    }
    write_pdf(plot_path, &pages, 20.0, 12.5)
}

/// With CUPPA_PLOT_BOXES=<file>, every grob's box after layout (for comparing the layout with R's viewports).
fn dump_boxes(table: &Gtable, w: f64, h: f64, page: usize) {
    if let Ok(path) = std::env::var("CUPPA_PLOT_BOXES") {
        let path = if page > 1 { format!("{path}.{page}") } else { path };
        let mut out = Vec::new();
        table.boxes(Bx::new(0.0, 0.0, w, h), "", &mut out);
        let text: String = out.iter().map(|(n, b)| format!("{:.4}\t{:.4}\t{:.4}\t{:.4}\t{n}\n", b.x0, b.y0, b.x1, b.y1)).collect();
        let _ = std::fs::write(path, text);
    }
}

fn write_pdf(path: &str, pages: &[Gtable], w_in: f64, h_in: f64) -> Result<()> {
    let (w, h) = (w_in * 72.0, h_in * 72.0);
    let mut figs = Vec::new();
    for (i, t) in pages.iter().enumerate() {
        let mut fig = Figure::new(w.round() as u32, h.round() as u32, Rgba::WHITE);
        let mut c = Canvas::new(&mut fig);
        t.draw(&mut c, Bx::new(0.0, 0.0, w, h));
        dump_boxes(t, w, h, i + 1);
        figs.push(fig);
    }
    hmf_plot::pdf_figure::write_r_pdf(path, &figs, w, h).map_err(|e| anyhow::anyhow!("writing {path}: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn str_wrap_as_stringi() {
        // stringr::str_wrap(x, width = 20, whitespace_only = FALSE) in the pipeline's CUPPA image (stringi 1.8.4)
        let cases = [
            ("Cartilaginous neoplasm", "Cartilaginous\nneoplasm"),
            ("Undiff. sarcoma", "Undiff. sarcoma"),
            ("Pilocytic astrocytoma", "Pilocytic\nastrocytoma"),
            ("Colorectum/Small intestine/Appendix", "Colorectum/Small\nintestine/Appendix"),
            ("Ovary/Fallopian tube", "Ovary/Fallopian tube"),
            ("Bile duct/Gallbladder", "Bile duct/\nGallbladder"),
            ("Non-small cell: LUAD", "Non-small cell: LUAD"),
            ("Acute myeloid leukemia", "Acute myeloid\nleukemia"),
            ("Myeloproliferative neoplasm", "Myeloproliferative\nneoplasm"),
            ("Colorectum/Small intestine", "Colorectum/Small\nintestine"),
            ("A very long cancer type name without breaks", "A very long cancer\ntype name without\nbreaks"),
            ("Head-and-neck/other-tissue type", "Head-and-neck/other-\ntissue type"),
            ("x/y/z/aaaaaaaaaaaaaaaaaaaaaaa", "x/y/z/\naaaaaaaaaaaaaaaaaaaaaaa"),
            ("word word word word word word", "word word word word\nword word"),
        ];
        for (s, w) in cases {
            assert_eq!(str_wrap(s, 20), w, "{s}");
        }
    }

    #[test]
    fn number_formats_as_r() {
        assert_eq!(format_c_e(27120.0, 0), "3e+04");
        assert_eq!(format_c_e(12345.0, 1), "1.2e+04");
        assert_eq!(format_c_e(0.00012, 0), "1e-04");
        assert_eq!(format_sig_quantile(150.10119047619045).text(), "150");
        assert_eq!(format_sig_quantile(2.1540104211155717).text(), "2.2");
        assert_eq!(format_sig_quantile(0.42066701121819233).text(), "0.42");
        assert_eq!(format_sig_quantile(f64::INFINITY).text(), "Inf");
        assert_eq!(format_odds(0.8965181405087232).text(), "0.9");
        assert_eq!(format_odds(0.0123).text(), "0.01");
        assert_eq!(sprintf_g(0.1), "0.1");
        assert_eq!(sprintf_g(10.0), "10");
        assert_eq!(sprintf_g(1e-6), "1e-06");
        let ct = CancerTypes::new(&[]);
        assert!(ct.types.is_empty());
        assert_eq!(leading_run("Lung: Non-small cell: LUAD"), "Lung");
    }
}
