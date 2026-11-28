//! qsee's visual report, `<tumor>.qsee.vis.report.pdf` (two pages) and (for a single patient) `.page_1.png` and
//! `.page_2.png`, drawn natively: a port of upstream 883ab11bba's R script `plot_qc.R` (ggplot2 4.0.3, patchwork 1.3.2
//! in the qsee image) on hmf-plot-rs's ggplot2 and patchwork layout. The same data steps (cohort merge, multi-field
//! split, factor orders, the peak marks, the missed-variant selection, the summary tables' formats and QC statuses, the
//! BQR values of too few counts hidden), the same plots (distributions with cohort ribbons, box plots and point ranges
//! against the cohort, facet grids, log ticks, the BQR plots' read strands as dodged triangles, the stacked plots with
//! their axes and axis titles collected), the custom legend, the report title, page 1's design of 4 x 9 cells and page
//! 2's BQR by SNV96 context, at 17 x 10 inches.

pub mod data;
pub mod gg;

use std::rc::Rc;

use anyhow::{Context, Result};
use hmf_plot::color::Rgba;
use hmf_plot::display::{Item, LineCap, Stroke};
use hmf_plot::gg::rfmt;
use hmf_plot::ggp::guides::{self, GuideBox, KeyGlyph, LegendSpec, LegendTheme};
use hmf_plot::ggp::scales;
use hmf_plot::ggp::{Ctx, PT, RectEl, STROKE, TextEl, TitleGrob, Version, dash_pattern_on, rect_grob, stroke};
use hmf_plot::gtable::{Bx, Canvas, Grob, Gtable};
use hmf_plot::patchwork::{self, Annotation, Layout, Patch};
use hmf_plot::rdev::{Device, Face};
use hmf_plot::{Figure, png_out::PngMeta};

use data::{Inputs, PlotData, Row};
use gg::{BASE, HALF, LINE, Label, Layers, Map, Pos, Spec, Theme, col};

const TUMOR_COL: &str = "#D1392C";
const NORMAL_COL: &str = "#4A7DB4";

fn type_colour(st: &str) -> Rgba {
    col(if st == "NORMAL" { NORMAL_COL } else { TUMOR_COL })
}

fn qc_colour(status: &str) -> Rgba {
    match status {
        "FAIL" => col("#F3C27B"),
        "WARN" => col("#FCFFC6"),
        _ => Rgba::TRANSPARENT,
    }
}

// -------------------------------------------------------------------------------------------------------------
// scales

/// How a continuous scale's breaks are chosen.
#[derive(Clone, Copy)]
enum Breaks {
    /// scales::extended_breaks (identity) or breaks_log (log10)
    Auto,
    /// scales::breaks_width(w)
    Width(f64),
}

/// A continuous position scale trained on `values` (data space), with scale limits (NA: from the data), coord
/// limits, 5% expansion, breaks on the expanded range, labels of the breaks inside it, and log ticks.
#[allow(clippy::too_many_arguments)]
fn continuous(values: &[f64], scale_limits: (Option<f64>, Option<f64>), coord_limits: (Option<f64>, Option<f64>), log: bool, breaks: Breaks, labels: &dyn Fn(&[f64]) -> Vec<String>, logticks: bool) -> Pos {
    let tr = |v: f64| if log { v.log10() } else { v };
    let inv = |v: f64| if log { 10f64.powf(v) } else { v };
    let finite: Vec<f64> = values.iter().map(|&v| tr(v)).filter(|v| v.is_finite()).collect();
    let (mut lo, mut hi) = if finite.is_empty() { (0.0, 1.0) } else { (finite.iter().copied().fold(f64::INFINITY, f64::min), finite.iter().copied().fold(f64::NEG_INFINITY, f64::max)) };
    if let Some(l) = scale_limits.0 {
        lo = tr(l);
    }
    if let Some(h) = scale_limits.1 {
        hi = tr(h);
    }
    let limits = (lo, hi);
    let (mut clo, mut chi) = limits;
    if let Some(l) = coord_limits.0 {
        clo = tr(l);
    }
    if let Some(h) = coord_limits.1 {
        chi = tr(h);
    }
    let range = scales::expand((clo, chi), 0.05, 0.0);
    let data_range = (inv(range.0), inv(range.1));
    let raw: Vec<f64> = match breaks {
        Breaks::Auto => {
            if log {
                scales::breaks_log(data_range)
            } else {
                rfmt::extended_breaks(data_range.0, data_range.1, 5)
            }
        }
        Breaks::Width(w) => scales::breaks_width(data_range, w),
    };
    let kept: Vec<f64> = raw.into_iter().filter(|&b| {
        let t = tr(b);
        t >= range.0 && t <= range.1
    }).collect();
    let labs = labels(&kept);
    let ticks = if logticks {
        scales::log_ticks((inv(limits.0), inv(limits.1))).into_iter().map(|(v, k)| (tr(v), k)).filter(|(t, _)| *t >= range.0 && *t <= range.1).collect()
    } else {
        Vec::new()
    };
    Pos { range, log, breaks: kept.iter().map(|&b| tr(b)).collect(), labels: labs, logticks: ticks }
}

/// A discrete position scale of n levels (positions 1..n) with labels, expanded by 0.6 on each side; continuous
/// positions (dodged points, box edges, lines between levels) widen it only where they reach past that
/// (ggplot2 4.0.3's expand_limits_discrete_trans).
fn discrete(n: usize, labels: Vec<String>, cont: &[f64]) -> Pos {
    let (mut lo, mut hi) = if n == 0 { (f64::INFINITY, f64::NEG_INFINITY) } else { (1.0 - 0.6, n as f64 + 0.6) };
    for &v in cont.iter().filter(|v| v.is_finite()) {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if !lo.is_finite() {
        // nothing at all: expand_limits_continuous_trans(c(0, 1))
        lo = -0.6;
        hi = 1.6;
    }
    Pos { range: (lo, hi), log: false, breaks: (1..=n).map(|i| i as f64).collect(), labels, logticks: Vec::new() }
}

fn default_labels(xs: &[f64]) -> Vec<String> {
    rfmt::r_format(xs, true)
}

/// R's signif(x, 12)
fn signif12(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let e = 11 - x.abs().log10().floor() as i32;
    let m = 10f64.powi(e.abs());
    if e >= 0 { (x * m).round() / m } else { (x / m).round() * m }
}

fn percent(xs: &[f64]) -> Vec<String> {
    scales::label_percent(xs, None, 100.0, false)
}

fn percent_drop0(xs: &[f64]) -> Vec<String> {
    scales::label_percent(xs, None, 100.0, true)
}

/// labels = function(x) ifelse(x > 0, paste0("+", x), x)
fn plus_labels(xs: &[f64]) -> Vec<String> {
    xs.iter().map(|&x| if x > 0.0 { format!("+{}", rfmt::r_as_character(x)) } else { rfmt::r_as_character(x) }).collect()
}

/// format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
fn plain_labels(xs: &[f64]) -> Vec<String> {
    scales::format_fixed(xs, 0, "", true)
}

// -------------------------------------------------------------------------------------------------------------
// geoms

fn path(c: &mut Canvas, ctx: &Ctx, pts: Vec<(f64, f64)>, colour: Rgba, linewidth: f64, lty: &str) {
    let lwd = linewidth * PT;
    let dash = if lty.is_empty() { Vec::new() } else { dash_pattern_on(ctx, lty, lwd) };
    c.fig.polyline(pts, stroke(ctx, colour, lwd, dash, LineCap::Butt));
}

/// geom_point with shape 21 (fill and outline) or 15 (filled square); size and stroke as ggplot2 4.0.3
fn point21(c: &mut Canvas, ctx: &Ctx, at: (f64, f64), colour: Rgba, fill: Rgba, size: f64, stroke_w: f64) {
    hmf_plot::ggp::geom::point(c, ctx, at.0, at.1, 21, colour, fill, size, stroke_w);
}

/// geom_label (ggplot2 4.0.3's labelGrob, just "centre"): the text centred on its point (vjust 0.5 of the "M" height);
/// the box padded by `padding` left and right, max(padding, descent) on top and max(padding - descent, 0) below,
/// centred on the point (`y + (0.5 - vjust) * height`; only the text is raised, by `(1 - vjust) * descent`); corners
/// of radius 0.15 lines; the border `border_mm` wide.
#[allow(clippy::too_many_arguments)]
fn label(c: &mut Canvas, ctx: &Ctx, at: (f64, f64), s: &str, size_mm: f64, text_col: Rgba, fill: Rgba, border: Rgba, border_mm: f64, padding: f64) {
    let dev = ctx.dev;
    let fs = size_mm * PT;
    let w = dev.str_width(s, fs, Face::Plain);
    let asc = dev.m_ascent(fs, Face::Plain);
    let desc = dev.font_descent(fs, Face::Plain);
    let (mt, mb) = (padding.max(desc), (padding - desc).max(0.0));
    let (bw, bh) = (w + 2.0 * padding, asc + desc + mt + mb);
    // the box centred on the point
    let (x0, y0) = (at.0 - bw / 2.0, at.1 - bh / 2.0);
    let r = dev.lines(0.15, 12.0, 1.2).min(bw / 2.0).min(bh / 2.0);
    hmf_plot::ggp::drawn_box("rrvp", Bx::new(x0, y0, x0 + bw, y0 + bh));
    let mut pts = Vec::new();
    let corner = |pts: &mut Vec<(f64, f64)>, cx: f64, cy: f64, a0: f64| {
        for i in 0..=8 {
            let a = (a0 + 90.0 * i as f64 / 8.0).to_radians();
            pts.push((cx + r * a.cos(), cy + r * a.sin()));
        }
    };
    corner(&mut pts, x0 + bw - r, y0 + r, -90.0);
    corner(&mut pts, x0 + bw - r, y0 + bh - r, 0.0);
    corner(&mut pts, x0 + r, y0 + bh - r, 90.0);
    corner(&mut pts, x0 + r, y0 + r, 180.0);
    let st = (!border.is_transparent()).then(|| stroke(ctx, border, border_mm * PT, Vec::new(), LineCap::Round));
    c.fig.polygon(pts, (!fill.is_transparent()).then_some(fill), st);
    gg::text(c, ctx, at.0, at.1, s, fs, text_col, 0.5, 0.5);
}

// -------------------------------------------------------------------------------------------------------------
// the plots

/// The theme's base plot spec (one panel).
fn base_spec(th: &Theme, title: Option<&str>, xlab: Option<&str>, ylab: Option<&str>) -> Spec {
    Spec {
        title: title.map(String::from),
        subtitle: None,
        xlab: xlab.map(|s| Label::Text(s.into())),
        ylab: ylab.map(|s| Label::Text(s.into())),
        ylab_r: None,
        h: Vec::new(),
        v: Vec::new(),
        flip: false,
        strips_top: None,
        strips_right: None,
        layers: Vec::new(),
        h_text: Some(th.axis_text_h()),
        v_text: Some(th.axis_text_v()),
        h_ticks: true,
        v_ticks: true,
        right_axis: false,
        grid_h: false,
        grid_v: false,
        spacing_x: th.pt(HALF),
        spacing_y: th.pt(HALF),
        margin: th.margin(),
        title_blank: false,
    }
}

/// plot_missing_data: an empty plot with the labels and "Missing sample data" in the middle.
fn missing(th: &Theme, title: &str, xlab: &str, ylab: &str) -> Gtable {
    gg::table(th, missing_spec(th, Some(title), xlab, ylab))
}

/// plot_missing_data's spec
fn missing_spec(th: &Theme, title: Option<&str>, xlab: &str, ylab: &str) -> Spec {
    let mut spec = base_spec(th, title, Some(xlab), Some(ylab));
    let pos = Pos { range: (-0.5, 0.5), log: false, breaks: Vec::new(), labels: Vec::new(), logticks: Vec::new() };
    spec.h = vec![pos.clone()];
    spec.v = vec![pos];
    spec.h_text = None;
    spec.v_text = None;
    spec.h_ticks = false;
    spec.v_ticks = false;
    spec.layers = vec![vec![Rc::new(|c: &mut Canvas, ctx: &Ctx, m: &Map| {
        let (x, y) = m.at(0.0, 0.0);
        gg::text(c, ctx, x, y, "Missing sample data", BASE, Rgba::BLACK, 0.5, 0.5);
    }) as Layers]];
    spec
}

/// plot_distribution: cohort ribbons (5th to 95th percentile), the sample's lines, the peak marks, a line at 0.
fn distribution(th: &Theme, d: &PlotData, x: &str, invert_normal: bool, mark_peak: bool, hline0: bool) -> (Vec<f64>, Vec<f64>, Layers) {
    let st_levels = d.levels("SampleType");
    let invert = invert_normal && st_levels.len() > 1;
    struct P {
        st: String,
        x: f64,
        y: Option<f64>,
        lo: Option<f64>,
        hi: Option<f64>,
    }
    let pts: Vec<P> = d
        .rows
        .iter()
        .map(|r| {
            let s = if invert && r.sample_type == "NORMAL" { -1.0 } else { 1.0 };
            P { st: r.sample_type.clone(), x: r.num(x), y: r.value.map(|v| s * v), lo: r.pct[0].map(|v| s * v), hi: r.pct[4].map(|v| s * v) }
        })
        .collect();
    // the peak marks: per sample type, the rows floor(n * 0.05):ceiling(n * 0.95), the largest |value|
    let mut peaks: Vec<(String, f64, f64)> = Vec::new();
    if mark_peak {
        for st in &st_levels {
            let g: Vec<&P> = pts.iter().filter(|p| &p.st == st).collect();
            let n = g.len() as f64;
            let (a, b) = ((n * 0.05).floor() as usize, (n * 0.95).ceil() as usize);
            let slice: Vec<&&P> = g.iter().enumerate().filter(|(i, _)| i + 1 >= a.max(1) && i + 1 <= b).map(|(_, p)| p).collect();
            let mut best: Option<&&P> = None;
            for p in &slice {
                let v = p.y.map(f64::abs).unwrap_or(f64::NAN);
                if best.map(|q| v > q.y.map(f64::abs).unwrap_or(f64::NAN)).unwrap_or(!v.is_nan()) {
                    best = Some(p);
                }
            }
            if let Some(p) = best {
                peaks.push((st.clone(), p.x, p.y.unwrap()));
            }
        }
    }
    let mut xs: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let mut ys: Vec<f64> = pts.iter().filter_map(|p| p.y).collect();
    if th_has_cohort(d) {
        ys.extend(pts.iter().filter_map(|p| p.lo));
        ys.extend(pts.iter().filter_map(|p| p.hi));
    }
    if hline0 {
        ys.push(0.0);
    }
    for (_, px, h) in &peaks {
        xs.push(*px);
        ys.push(0.0);
        ys.push(*h);
    }
    let cohort = th_has_cohort(d);
    let _ = th;
    // drawing: per group (sample type, in level order) sorted by x
    let mut groups: Vec<(String, Vec<(f64, Option<f64>, Option<f64>, Option<f64>)>)> = Vec::new();
    for st in &st_levels {
        let mut g: Vec<(f64, Option<f64>, Option<f64>, Option<f64>)> = pts.iter().filter(|p| &p.st == st).map(|p| (p.x, p.y, p.lo, p.hi)).collect();
        g.sort_by(|a, b| a.0.total_cmp(&b.0));
        groups.push((st.clone(), g));
    }
    let layers: Layers = Rc::new(move |c: &mut Canvas, ctx: &Ctx, m: &Map| {
        if hline0 {
            let a = m.at(m.h.range.0, 0.0);
            let b = m.at(m.h.range.1, 0.0);
            path(c, ctx, vec![(m.b.x0, a.1), (m.b.x1, b.1)], col("grey90"), 0.21, "");
        }
        if cohort {
            for (st, g) in &groups {
                let fill = type_colour(st).with_alpha(0.1);
                // runs without missing bounds
                let mut run: Vec<(f64, f64, f64)> = Vec::new();
                let flush = |run: &mut Vec<(f64, f64, f64)>, c: &mut Canvas| {
                    if run.len() >= 2 {
                        let mut poly: Vec<(f64, f64)> = run.iter().map(|&(x, _, hi)| m.at(x, hi)).collect();
                        poly.extend(run.iter().rev().map(|&(x, lo, _)| m.at(x, lo)));
                        c.fig.polygon(poly, Some(fill), None);
                    }
                    run.clear();
                };
                for &(x, _, lo, hi) in g {
                    match (lo, hi) {
                        (Some(l), Some(h)) => run.push((x, l, h)),
                        _ => flush(&mut run, c),
                    }
                }
                flush(&mut run, c);
            }
        }
        for (st, g) in &groups {
            let mut seg: Vec<(f64, f64)> = Vec::new();
            for &(x, y, _, _) in g {
                match y {
                    Some(v) if v.is_finite() || v.is_infinite() => seg.push(m.at(x, v)),
                    _ => {
                        if seg.len() >= 2 {
                            path(c, ctx, std::mem::take(&mut seg), type_colour(st), LINE, "");
                        }
                        seg.clear();
                    }
                }
            }
            if seg.len() >= 2 {
                path(c, ctx, seg, type_colour(st), LINE, "");
            }
        }
        for (st, px, h) in &peaks {
            path(c, ctx, vec![m.at(*px, 0.0), m.at(*px, *h)], type_colour(st), 0.2, "");
        }
    });
    (xs, ys, layers)
}

fn th_has_cohort(d: &PlotData) -> bool {
    d.rows.iter().any(|r| r.pct.iter().any(|p| p.is_some()))
}

/// Coverage, fragment length, duplicate frequency and GC bias.
fn distribution_plot(th: &Theme, inp: &Inputs, ft: &str, title: &str, xlab: &str, ylab: &str) -> Gtable {
    let Some(d) = inp.plot_data(ft) else {
        return missing(th, title, xlab, ylab);
    };
    let (x, invert, peak, hline) = match ft {
        "COVERAGE_DISTRIBUTION" => ("ReadDepth", true, true, true),
        "FRAG_LENGTH_DISTRIBUTION" => ("FragLength", true, true, true),
        "DUPLICATE_FREQ" => ("ReadCount", false, false, false),
        _ => ("GCBucket", false, false, false),
    };
    let (xs, ys, layers) = distribution(th, &d, x, invert, peak, hline);
    let mut spec = base_spec(th, Some(title), Some(xlab), Some(ylab));
    let none = (None, None);
    match ft {
        "COVERAGE_DISTRIBUTION" | "FRAG_LENGTH_DISTRIBUTION" => {
            spec.h = vec![continuous(&xs, none, none, false, Breaks::Auto, &default_labels, false)];
            // label_percent()(abs(x)): the breaks' rounding noise makes |b| and |-b| differ in the last bits, so
            // scales::precision sees a gap below sqrt(eps) and labels at accuracy 1 ("2%" and "3%" for 2.5%, "0%" for
            // 0.25%); corrected, the absolute values are taken at 12 significant digits first (QS13)
            let abs_percent = |v: &[f64]| percent(&v.iter().map(|x| if crate::compat::java("QS13") { x.abs() } else { signif12(x.abs()) }).collect::<Vec<_>>());
            spec.v = vec![continuous(&ys, none, none, false, Breaks::Auto, &abs_percent, false)];
        }
        "DUPLICATE_FREQ" => {
            spec.h = vec![continuous(&xs, none, none, true, Breaks::Auto, &default_labels, true)];
            spec.v = vec![continuous(&ys, none, none, false, Breaks::Auto, &percent, false)];
            spec.grid_h = true;
        }
        _ => {
            let pct1 = |v: &[f64]| scales::label_percent(v, None, 1.0, false);
            spec.h = vec![continuous(&xs, none, none, false, Breaks::Auto, &pct1, false)];
            spec.v = vec![continuous(&ys, none, none, false, Breaks::Auto, &default_labels, false)];
        }
    }
    spec.layers = vec![vec![layers]];
    gg::table(th, spec)
}

/// The pairwise comparison (plot_pairwise_comparison): the sample's points against the cohort's box plots (BOX) or
/// dashed ranges (POINT_RANGE), dodged by sample type (width 0.5, preserve "single").
#[derive(Clone, Copy, PartialEq)]
enum Pairwise {
    Box,
    PointRange,
}

/// One element of a pairwise plot: discrete position, group (sample type level index), value, percentiles.
#[derive(Clone)]
struct PwRow {
    x: f64,
    st: String,
    group: usize,
    y: Option<f64>,
    pct: [Option<f64>; 5],
}

/// position_dodge(width = 0.5, preserve = "single") of rows at positions: the shift of each row.
fn dodge(rows: &[PwRow], keys: &[usize]) -> Vec<f64> {
    // n: the most groups at any one position (per panel; the caller passes one panel's rows)
    let mut n = 1;
    let mut by_x: Vec<(f64, Vec<usize>)> = Vec::new();
    for (i, r) in rows.iter().enumerate() {
        match by_x.iter_mut().find(|(x, _)| *x == r.x) {
            Some((_, v)) => v.push(i),
            None => by_x.push((r.x, vec![i])),
        }
    }
    for (_, v) in &by_x {
        let mut gs: Vec<usize> = v.iter().map(|&i| keys[i]).collect();
        gs.sort();
        gs.dedup();
        n = n.max(gs.len());
    }
    let mut shift = vec![0.0; rows.len()];
    if n == 1 {
        return shift;
    }
    for (_, v) in &by_x {
        let mut gs: Vec<usize> = v.iter().map(|&i| keys[i]).collect();
        gs.sort();
        gs.dedup();
        for &i in v {
            let idx = gs.iter().position(|&g| g == keys[i]).unwrap() as f64 + 1.0;
            shift[i] = 0.5 * ((idx - 0.5) / n as f64 - 0.5);
        }
    }
    shift
}

/// The layers and continuous extents of a pairwise plot of one panel: (x positions trained, y values trained, layers).
fn pairwise_layers(rows: Vec<PwRow>, kind: Pairwise, cohort: bool, box_width: f64, hlines: Vec<f64>, vlines: Vec<f64>, n_x: usize, extra: Vec<(f64, f64)>) -> (Vec<f64>, Vec<f64>, Layers) {
    let keys: Vec<usize> = rows.iter().map(|r| r.group).collect();
    let n_groups = {
        let mut g = keys.clone();
        g.sort();
        g.dedup();
        g.len()
    };
    let _ = n_groups;
    let shift = dodge(&rows, &keys);
    // dodged positions; boxes 0.4 * scale wide, split by the number of groups
    let n_dodge = {
        let mut n = 1;
        let mut xs: Vec<f64> = rows.iter().map(|r| r.x).collect();
        xs.sort_by(f64::total_cmp);
        xs.dedup();
        for x in xs {
            let mut gs: Vec<usize> = rows.iter().filter(|r| r.x == x).map(|r| r.group).collect();
            gs.sort();
            gs.dedup();
            n = n.max(gs.len());
        }
        n as f64
    };
    let half_box = if n_dodge > 1.0 { 0.4 * box_width / n_dodge / 2.0 } else { 0.4 * box_width / 2.0 };
    let mut xc: Vec<f64> = Vec::new();
    let mut yc: Vec<f64> = Vec::new();
    // collide() gives points and ranges xmin and xmax of x -+ width / 2 (0.25), which pos_dodge narrows to
    // x' -+ width / n / 2; both train the x scale
    let half_point = if n_dodge > 1.0 { 0.5 / n_dodge / 2.0 } else { 0.25 };
    for (r, s) in rows.iter().zip(&shift) {
        let x = r.x + s;
        if let Some(y) = r.y {
            xc.push(x - half_point);
            xc.push(x + half_point);
            yc.push(y);
        }
        if cohort {
            let has = r.pct.iter().any(|p| p.is_some());
            if has {
                match kind {
                    Pairwise::Box => {
                        xc.push(x - half_box);
                        xc.push(x + half_box);
                        yc.extend(r.pct.iter().filter_map(|p| *p));
                    }
                    Pairwise::PointRange => {
                        xc.push(x - half_point);
                        xc.push(x + half_point);
                        yc.extend([r.pct[0], r.pct[4]].iter().filter_map(|p| *p));
                    }
                }
            }
        }
    }
    yc.extend(hlines.iter().copied());
    xc.extend(vlines.iter().copied());
    for &(x, y) in &extra {
        xc.push(x);
        yc.push(y);
    }
    let _ = n_x;
    let layers: Layers = Rc::new(move |c: &mut Canvas, ctx: &Ctx, m: &Map| {
        let grey90 = col("grey90");
        for &h in &hlines {
            let (a, b) = if m.flip { (m.at(h, m.v.range.0), m.at(h, m.v.range.1)) } else { (m.at(m.h.range.0, h), m.at(m.h.range.1, h)) };
            let _ = (a, b);
            // geom_hline: across the panel at y = h (vertical in a flipped plot)
            if m.flip {
                let x = m.b.npc_x(m.h.npc(h));
                path(c, ctx, vec![(x, m.b.y0), (x, m.b.y1)], grey90, 0.21, "");
            } else {
                let y = m.b.npc_y(m.v.npc(h));
                path(c, ctx, vec![(m.b.x0, y), (m.b.x1, y)], grey90, 0.21, "");
            }
        }
        for &v in &vlines {
            if m.flip {
                let y = m.b.npc_y(m.v.npc(v));
                path(c, ctx, vec![(m.b.x0, y), (m.b.x1, y)], grey90, 0.21, "");
            } else {
                let x = m.b.npc_x(m.h.npc(v));
                path(c, ctx, vec![(x, m.b.y0), (x, m.b.y1)], grey90, 0.21, "");
            }
        }
        if cohort {
            match kind {
                Pairwise::Box => {
                    for (r, s) in rows.iter().zip(&shift) {
                        let [Some(mn), Some(lo), Some(md), Some(up), Some(mx)] = r.pct else { continue };
                        let x = r.x + s;
                        let outline = col("grey70");
                        let lw = BASE * 0.025;
                        // whiskers, then the box (crossbar) and its median (2 x the line width)
                        path(c, ctx, vec![m.at(x, up), m.at(x, mx)], outline, lw, "");
                        path(c, ctx, vec![m.at(x, lo), m.at(x, mn)], outline, lw, "");
                        let (a, b) = (m.at(x - half_box, lo), m.at(x + half_box, up));
                        hmf_plot::ggp::rect(c, ctx, a.0, a.1, b.0, b.1, Some(type_colour(&r.st).with_alpha(0.3)), Some((outline, lw * PT)), Vec::new());
                        path(c, ctx, vec![m.at(x - half_box, md), m.at(x + half_box, md)], outline, 2.0 * lw, "");
                    }
                }
                Pairwise::PointRange => {
                    for (r, s) in rows.iter().zip(&shift) {
                        let (Some(mn), Some(mx)) = (r.pct[0], r.pct[4]) else { continue };
                        let x = r.x + s;
                        path(c, ctx, vec![m.at(x, mn), m.at(x, mx)], type_colour(&r.st), 0.3, "11");
                    }
                }
            }
        }
        for (r, s) in rows.iter().zip(&shift) {
            let Some(y) = r.y else { continue };
            let at = m.at(r.x + s, y);
            match kind {
                Pairwise::Box => point21(c, ctx, at, Rgba::BLACK, type_colour(&r.st), BASE * 0.15, LINE),
                Pairwise::PointRange => point21(c, ctx, at, type_colour(&r.st), type_colour(&r.st), BASE * 0.07, LINE),
            }
        }
    });
    (xc, yc, layers)
}

/// The level index of a sample type among the levels.
fn group_of(levels: &[String], st: &str) -> usize {
    levels.iter().position(|l| l == st).unwrap_or(0)
}

/// Discordant fragment frequency (F).
fn discordant(th: &Theme, inp: &Inputs) -> Gtable {
    let (title, xlab, ylab) = ("Discordant fragment frequency", "Discordant fragment type", "Reads");
    let Some(d) = inp.plot_data("DISCORDANT_FRAG_FREQ") else {
        return missing(th, title, xlab, ylab);
    };
    // DisplayName = reverse_levels(DisplayName)
    let mut names = d.levels("DisplayName");
    names.reverse();
    let st_levels = d.levels("SampleType");
    let rows: Vec<PwRow> = d
        .rows
        .iter()
        .map(|r| PwRow { x: (names.iter().position(|n| n == r.field("DisplayName")).unwrap() + 1) as f64, st: r.sample_type.clone(), group: group_of(&st_levels, &r.sample_type), y: r.value, pct: r.pct })
        .collect();
    let all_zero = d.rows.iter().all(|r| r.value == Some(0.0));
    let (xc, yc, layers) = pairwise_layers(rows, Pairwise::Box, inp.has_cohort, 1.0, Vec::new(), Vec::new(), names.len(), Vec::new());
    let mut spec = base_spec(th, Some(title), Some(ylab), Some(xlab));
    spec.flip = true;
    spec.v = vec![discrete(names.len(), names.clone(), &xc)];
    spec.h = vec![if all_zero {
        continuous(&yc, (None, None), (None, None), false, Breaks::Auto, &percent_drop0, false)
    } else {
        continuous(&yc, (None, None), (None, None), true, Breaks::Auto, &percent_drop0, true)
    }];
    spec.grid_h = true;
    spec.spacing_x = th.pt(-0.5);
    spec.layers = vec![vec![layers]];
    gg::table(th, spec)
}

/// Genes with potential missed variants (G).
fn missed_variants(th: &Theme, inp: &Inputs) -> Gtable {
    let (title, xlab, ylab) = ("Genes with potential missed variants", "Genes (top 20 per sample type)", "Missed variant likelihood (>2%)");
    let Some(d) = inp.plot_data("MISSED_VARIANT_LIKELIHOOD") else {
        return missing(th, title, xlab, ylab);
    };
    // filter >= 0.02; per sample type (in level order) the 20 largest (arrange(-value), stable); then Gene and
    // SampleType factors in order of appearance, reversed
    let st_levels_all = d.levels("SampleType");
    let mut kept: Vec<&Row> = Vec::new();
    for st in &st_levels_all {
        let mut g: Vec<&Row> = d.rows.iter().filter(|r| &r.sample_type == st && r.value.is_some_and(|v| v >= 0.02)).collect();
        g.sort_by(|a, b| b.value.unwrap().total_cmp(&a.value.unwrap()));
        kept.extend(g.into_iter().take(20));
    }
    // dplyr::arrange inside group_by sorts the whole table by -value (groups are not kept together)
    kept.sort_by(|a, b| b.value.unwrap().total_cmp(&a.value.unwrap()));
    let mut genes: Vec<String> = Vec::new();
    let mut sts: Vec<String> = Vec::new();
    for r in &kept {
        if !genes.iter().any(|g| g == r.field("Gene")) {
            genes.push(r.field("Gene").to_string());
        }
        if !sts.contains(&r.sample_type) {
            sts.push(r.sample_type.clone());
        }
    }
    let unique_genes = genes.clone();
    genes.reverse();
    sts.reverse();
    let gene_count = unique_genes.len();
    let split = gene_count > 10;
    let mut spec = base_spec(th, Some(title), Some(ylab), Some(xlab));
    spec.flip = true;
    spec.grid_h = true;
    let to_row = |r: &Row, levels: &[String]| PwRow { x: (levels.iter().position(|g| g == r.field("Gene")).unwrap() + 1) as f64, st: r.sample_type.clone(), group: group_of(&sts, &r.sample_type), y: r.value, pct: r.pct };
    if !split {
        let rows: Vec<PwRow> = kept.iter().map(|r| to_row(r, &genes)).collect();
        let (xc, yc, layers) = pairwise_layers(rows, Pairwise::Box, inp.has_cohort, 1.0, Vec::new(), Vec::new(), genes.len(), Vec::new());
        spec.v = vec![discrete(genes.len(), genes.clone(), &xc)];
        spec.h = vec![continuous(&yc, (None, None), (Some(0.0), None), false, Breaks::Auto, &percent_drop0, false)];
        spec.layers = vec![vec![layers]];
        spec.spacing_x = th.pt(3.0);
        return gg::table(th, spec);
    }
    // facet_wrap(~ GeneGroup, scales = "free_y"): two panels side by side; the genes past the midpoint in the second
    let midpoint = rfmt::r_round(gene_count as f64 / 2.0, 0) as usize;
    let group_of_gene = |g: &str| usize::from(unique_genes.iter().position(|x| x == g).unwrap() + 1 > midpoint);
    let mut spec_h = Vec::new();
    let mut layers_row = Vec::new();
    let mut vs = Vec::new();
    let mut all_y = Vec::new();
    let mut panels = Vec::new();
    for grp in 0..2 {
        let levels: Vec<String> = genes.iter().filter(|g| group_of_gene(g) == grp).cloned().collect();
        let rows: Vec<PwRow> = kept.iter().filter(|r| group_of_gene(r.field("Gene")) == grp).map(|r| to_row(r, &levels)).collect();
        let (xc, yc, layers) = pairwise_layers(rows, Pairwise::Box, inp.has_cohort, 1.0, Vec::new(), Vec::new(), levels.len(), Vec::new());
        all_y.extend(yc);
        vs.push(discrete(levels.len(), levels.clone(), &xc));
        panels.push(layers);
    }
    let hpos = continuous(&all_y, (None, None), (Some(0.0), None), false, Breaks::Auto, &percent_drop0, false);
    for l in panels {
        spec_h.push(hpos.clone());
        layers_row.push(l);
    }
    // the wrap's panels are drawn as a grid of one row with their own vertical axes (free y)
    let mut h_text = th.axis_text_h();
    h_text.angle = 30.0;
    h_text.hjust = 1.0;
    h_text.vjust = 1.0;
    let mut v_text = th.axis_text_v();
    v_text.size = BASE * 0.5;
    spec.h_text = Some(h_text);
    spec.v_text = Some(v_text);
    spec.h = spec_h;
    spec.v = vec![vs[0].clone()];
    spec.layers = vec![layers_row];
    spec.spacing_x = th.pt(3.0);
    let mut s2 = spec;
    s2.v = vs;
    gg::table_wrap(th, s2)
}

/// BQR by original base quality (H) and by SNV96 context (I).
/// What a report plot is in patchwork terms: a ggplot's table, or a built patchwork's (nested where it is used).
enum Built {
    Plot(Gtable),
    Patchwork(Gtable),
}

/// The patch of a ggplot or of a patchwork placed in a design as is (not freed).
fn patch_of(b: Built) -> Patch {
    match b {
        Built::Plot(gt) => patchwork::plot(gt),
        Built::Patchwork(gt) => patchwork::nested(gt),
    }
}

/// free(p, "label"): a ggplot's labels freed; a patchwork's (free.patchwork) through merge
fn free_label_of(b: Built) -> Patch {
    match b {
        Built::Plot(gt) => {
            let mut p = patchwork::plot(gt);
            patchwork::free_label(&mut p.gt, [true; 4]);
            p
        }
        Built::Patchwork(gt) => {
            let mut p = patchwork::nested(gt);
            patchwork::free_label(&mut p.gt, [true; 4]);
            patchwork::merge(p)
        }
    }
}

/// The keys of a subplot's axes and axis titles, for patchwork's collect_axes and collect_axis_titles (grobs
/// compared by what they draw).
fn set_collect_keys(gt: &mut Gtable, h: &[Pos], h_text: &Option<TextEl>, h_ticks: bool, xlab: Option<&str>, ylab: Option<&str>, ylab_r: Option<&str>) {
    let hk: Vec<String> = h.iter().map(|p| format!("{:?}|{:?}|{:?}", p.range, p.breaks, p.labels)).collect();
    let tk = h_text.as_ref().map(|t| format!("{}:{}:{}:{}:{}", t.size, t.angle, t.hjust, t.vjust, t.margin[0])).unwrap_or_default();
    gt.set_keys("axis-b", &format!("axis-b[{}][{tk}][{h_ticks}]", hk.join(";")));
    if let Some(x) = xlab {
        gt.set_keys("xlab-b", &format!("xlab-b[{x}]"));
    }
    if let Some(y) = ylab {
        gt.set_keys("ylab-l", &format!("ylab-l[{y}]"));
    }
    if let Some(y) = ylab_r {
        gt.set_keys("ylab-r", &format!("ylab-r[{y}]"));
    }
}

/// READ_ORIENTATION_SHAPES: FORWARD 24 (triangle up), REVERSE 25 (down)
fn orientation_shape(o: &str) -> u8 {
    if o == "REVERSE" { 25 } else { 24 }
}

/// One element of a BQR panel: discrete position, dodge group order (the read strand's level), sample type,
/// value (NaN as NA), percentiles (NA where the value is).
#[derive(Clone)]
struct BqrRow {
    x: f64,
    group: usize,
    st: String,
    orientation: String,
    y: Option<f64>,
    pct: [Option<f64>; 5],
}

/// position_dodge(width = 0.6) (preserve "total"): at each position the groups present, in group order
fn dodge_total(rows: &[BqrRow]) -> (Vec<f64>, Vec<usize>) {
    let mut shift = vec![0.0; rows.len()];
    let mut ns = vec![1usize; rows.len()];
    let mut xs: Vec<f64> = rows.iter().map(|r| r.x).collect();
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    for x in xs {
        let idx: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].x == x).collect();
        let mut gs: Vec<usize> = idx.iter().map(|&i| rows[i].group).collect();
        gs.sort();
        gs.dedup();
        let n = gs.len();
        for &i in &idx {
            ns[i] = n;
            if n > 1 {
                let k = gs.iter().position(|&g| g == rows[i].group).unwrap() as f64 + 1.0;
                shift[i] = 0.6 * ((k - 0.5) / n as f64 - 0.5);
            }
        }
    }
    (shift, ns)
}

/// The layers of a BQR panel and its trained x positions: a line at 0 (inside the limits), the cohort's boxes
/// (width 0.5 split by the groups at a position), the sample's triangles (fill by sample type), and (SNV96) the
/// lines between the contexts on top.
fn bqr_layers(rows: Vec<BqrRow>, cohort: bool, zero_line: bool, snv96: bool) -> (Vec<f64>, Layers) {
    let (shift, ns) = dodge_total(&rows);
    let mut xc = Vec::new();
    for ((r, s), n) in rows.iter().zip(&shift).zip(&ns) {
        let x = r.x + s;
        // points: collide's x -+ 0.3, narrowed by pos_dodge to -+ 0.6 / n / 2; boxes -+ 0.5 / n / 2
        let hp = if *n > 1 { 0.6 / *n as f64 / 2.0 } else { 0.3 };
        let hb = if *n > 1 { 0.5 / *n as f64 / 2.0 } else { 0.25 };
        xc.push(x - hp);
        xc.push(x + hp);
        if cohort && r.pct.iter().any(|p| p.is_some()) {
            xc.push(x - hb);
            xc.push(x + hb);
        }
    }
    if snv96 {
        // the vlines train the x scale
        xc.extend((1..=15).map(|i| i as f64 + 0.5));
    }
    let layers: Layers = Rc::new(move |c: &mut Canvas, ctx: &Ctx, m: &Map| {
        let grey90 = col("grey90");
        if zero_line {
            let y = m.b.npc_y(m.v.npc(0.0));
            path(c, ctx, vec![(m.b.x0, y), (m.b.x1, y)], grey90, 0.21, "");
        }
        if cohort {
            for ((r, s), n) in rows.iter().zip(&shift).zip(&ns) {
                let [Some(mn), Some(lo), Some(md), Some(up), Some(mx)] = r.pct else { continue };
                let x = r.x + s;
                let hb = if *n > 1 { 0.5 / *n as f64 / 2.0 } else { 0.25 };
                let outline = col("grey70");
                let lw = BASE * 0.025;
                path(c, ctx, vec![m.at(x, up), m.at(x, mx)], outline, lw, "");
                path(c, ctx, vec![m.at(x, lo), m.at(x, mn)], outline, lw, "");
                let (a, b) = (m.at(x - hb, lo), m.at(x + hb, up));
                hmf_plot::ggp::rect(c, ctx, a.0, a.1, b.0, b.1, Some(type_colour(&r.st).with_alpha(0.3)), Some((outline, lw * PT)), Vec::new());
                path(c, ctx, vec![m.at(x - hb, md), m.at(x + hb, md)], outline, 2.0 * lw, "");
            }
        }
        for (r, s) in rows.iter().zip(&shift) {
            let Some(y) = r.y else { continue };
            let at = m.at(r.x + s, y);
            hmf_plot::ggp::geom::point(c, ctx, at.0, at.1, orientation_shape(&r.orientation), Rgba::BLACK, type_colour(&r.st), BASE * 0.15, LINE);
        }
        if snv96 {
            // & geom_vline: solid where the 5' base changes (every 4 contexts), dotted ("26", grey80) otherwise
            for i in [4, 8, 12] {
                let x = m.b.npc_x(m.h.npc(i as f64 + 0.5));
                path(c, ctx, vec![(x, m.b.y0), (x, m.b.y1)], grey90, 0.21, "");
            }
            for i in (1..=15).filter(|i| ![4, 8, 12].contains(i)) {
                let x = m.b.npc_x(m.h.npc(i as f64 + 0.5));
                path(c, ctx, vec![(x, m.b.y0), (x, m.b.y1)], col("grey80"), 0.21, "26");
            }
        }
    });
    (xc, layers)
}

/// The lower bound of a quality bin label, e.g. 30 of "HIGH (30+)" (sub(".*\\((\\d+).*", "\\1", x))
fn bin_lower_bound(s: &str) -> f64 {
    match s.find('(') {
        Some(i) => {
            let digits: String = s[i + 1..].chars().take_while(|c| c.is_ascii_digit()).collect();
            digits.parse().unwrap_or(f64::NAN)
        }
        None => f64::NAN,
    }
}

/// The BQR plots (upstream 883ab11bba's plot_bqr): BQR by original base quality (H, page 1: x the sample type,
/// rows the quality bins) and BQR by SNV96 context (page 2: x the context, rows the sample types); per read type a
/// facet grid with the mutations as columns (free x scales), the read types stacked with heights by their rows, the
/// y limits shared, the x axes and axis titles collected. Missing data: the missing-data plot.
fn bqr(th: &Theme, inp: &Inputs, ft: &str) -> Built {
    let snv96 = ft == "BQR_PER_SNV96_CONTEXT";
    let (title, xlab, ylab) = if snv96 { ("BQR by SNV96 context", "Mutation context", "Phred score adjustment") } else { ("BQR by original base quality", "Sample type", "Phred score adjustment") };
    let Some(d) = inp.plot_data(ft) else {
        return Built::Plot(missing(th, title, xlab, ylab));
    };
    let title = if snv96 { format!("{title}, original base quality: {}", d.rows.first().map(|r| r.field("OriginalQualBin")).unwrap_or("NA")) } else { title.to_string() };
    let x_field = if snv96 { "StandardTrinucContext" } else { "SampleType" };
    let row_field = if snv96 { "SampleType" } else { "OriginalQualBin" };
    let field = |r: &Row, f: &str| -> String { if f == "SampleType" { r.sample_type.clone() } else { r.field(f).to_string() } };
    let st_levels = d.levels("SampleType");
    let x_levels = if snv96 { d.levels(x_field) } else { st_levels.clone() };
    let mut row_levels = if snv96 { st_levels.clone() } else { d.levels(row_field) };
    if !snv96 {
        // levels ordered by the bins' lower bounds (a stable order)
        let mut keyed: Vec<(f64, String)> = row_levels.iter().map(|l| (bin_lower_bound(l), l.clone())).collect();
        keyed.sort_by(|a, b| a.0.total_cmp(&b.0));
        row_levels = keyed.into_iter().map(|(_, l)| l).collect();
    }
    let mut_levels = d.levels("StandardMutation");
    let orient_levels = d.levels("ReadOrientation");
    let rt_levels = d.levels("ReadType");
    // NaN values (too few counts) are NA, and so are their cohort percentiles
    let value = |r: &Row| r.value.filter(|v| !v.is_nan());
    let pct_of = |r: &Row| if value(r).is_some() { r.pct } else { [None; 5] };
    let mut all: Vec<f64> = Vec::new();
    for r in &d.rows {
        if let Some(v) = value(r) {
            all.push(v);
        }
        let p = pct_of(r);
        all.extend([p[0], p[4]].iter().filter_map(|x| *x));
    }
    let finite: Vec<f64> = all.iter().copied().filter(|v| v.is_finite()).collect();
    let y_limits = if finite.is_empty() { (0.0, 1.0) } else { (finite.iter().copied().fold(f64::INFINITY, f64::min), finite.iter().copied().fold(f64::NEG_INFINITY, f64::max)) };
    // scale_y_continuous(limits = y_limits, breaks = pretty breaks away from the edges, labels "+x"): the breaks
    // function is given the expanded range (view_scale's continuous range), and its breaks are kept inside it
    let range = scales::expand(y_limits, 0.05, 0.0);
    let margin = 0.1 * (range.1 - range.0);
    let breaks: Vec<f64> = if range.1 - range.0 == 0.0 {
        vec![range.0]
    } else {
        hmf_plot::ggp::scales::pretty(range.0, range.1, 5).into_iter().filter(|&b| b > range.0 + margin && b < range.1 - margin).collect()
    };
    let breaks: Vec<f64> = breaks.into_iter().filter(|&b| b >= range.0 && b <= range.1).collect();
    let v = Pos { range, log: false, breaks: breaks.clone(), labels: plus_labels(&breaks), logticks: Vec::new() };
    let zero_line = y_limits.0 <= 0.0 && 0.0 <= y_limits.1;
    let mut subplots: Vec<Patch> = Vec::new();
    let mut heights: Vec<Option<f64>> = Vec::new();
    for (ri, rt) in rt_levels.iter().enumerate() {
        let rt_rows: Vec<&Row> = d.rows.iter().filter(|r| r.field("ReadType") == rt).collect();
        if rt_rows.is_empty() {
            continue;
        }
        let rows_present: Vec<String> = row_levels.iter().filter(|l| rt_rows.iter().any(|r| field(r, row_field) == **l)).cloned().collect();
        let cols_present: Vec<String> = mut_levels.iter().filter(|l| rt_rows.iter().any(|r| r.field("StandardMutation") == *l)).cloned().collect();
        heights.push(Some(rows_present.len() as f64));
        // free x: each column's levels present
        let col_levels: Vec<Vec<String>> = cols_present.iter().map(|c| x_levels.iter().filter(|l| rt_rows.iter().any(|r| r.field("StandardMutation") == c && field(r, x_field) == **l)).cloned().collect()).collect();
        let mut col_x: Vec<Vec<f64>> = vec![Vec::new(); cols_present.len()];
        let mut grid: Vec<Vec<Layers>> = Vec::new();
        for rf in &rows_present {
            let mut row = Vec::new();
            for (ci, c) in cols_present.iter().enumerate() {
                let lv = &col_levels[ci];
                let rows: Vec<BqrRow> = rt_rows
                    .iter()
                    .filter(|r| field(r, row_field) == *rf && r.field("StandardMutation") == c)
                    .map(|r| {
                        let o = r.field("ReadOrientation").to_string();
                        BqrRow {
                            x: (lv.iter().position(|l| *l == field(r, x_field)).unwrap() + 1) as f64,
                            group: orient_levels.iter().position(|l| *l == o).unwrap_or(0),
                            st: r.sample_type.clone(),
                            orientation: o,
                            y: value(r),
                            pct: pct_of(r),
                        }
                    })
                    .collect();
                let (xc, layers) = bqr_layers(rows, inp.has_cohort, zero_line, snv96);
                col_x[ci].extend(xc);
                row.push(layers);
            }
            grid.push(row);
        }
        let hs: Vec<Pos> = col_levels.iter().zip(&col_x).map(|(lv, xs)| discrete(lv.len(), lv.clone(), xs)).collect();
        let mut spec = base_spec(th, if ri == 0 { Some(title.as_str()) } else { None }, if snv96 { Some(xlab) } else { None }, Some(ylab));
        spec.subtitle = Some((format!("Consensus type: {rt}"), BASE * 0.85));
        spec.strips_top = Some(cols_present.clone());
        spec.strips_right = Some(rows_present.iter().map(|l| if snv96 { if l == "NORMAL" { "Normal".to_string() } else if l == "TUMOR" { "Tumor".to_string() } else { l.clone() } } else { l.clone() }).collect());
        spec.spacing_x = th.pt(-0.5);
        spec.spacing_y = th.pt(-0.5);
        if snv96 {
            let mut h_text = th.axis_text_h();
            h_text.angle = 90.0;
            h_text.hjust = 1.0;
            h_text.vjust = 0.5;
            h_text.size = BASE * 0.5;
            spec.h_text = Some(h_text);
        } else {
            spec.h_text = None;
            spec.h_ticks = false;
        }
        spec.h = hs.clone();
        spec.v = vec![v.clone(); rows_present.len()];
        spec.layers = grid;
        let h_text = spec.h_text.clone();
        let h_ticks = spec.h_ticks;
        let mut gt = gg::table(th, spec);
        set_collect_keys(&mut gt, &hs, &h_text, h_ticks, if snv96 { Some(xlab) } else { None }, Some(ylab), None);
        subplots.push(patchwork::plot(gt));
    }
    let n = subplots.len();
    let layout = Layout { design: (1..=n).map(|i| patchwork::Area::cell(i, 1)).collect(), widths: vec![None], heights };
    Built::Patchwork(patchwork::build_collect(subplots, &layout, patchwork::Collect { axes_x: true, axes_y: false, titles_x: true, titles_y: true }))
}

/// Microsatellite indel error rates and bias (upstream 883ab11bba: the rates titled "Microsatellite indel errors"
/// with the subtitle "Rates", the bias with the subtitle "Bias", its y axis "Phred score diff." with the direction
/// of the bias written in the first column, and its y range expanded by 12%).
fn ms_indel(th: &Theme, inp: &Inputs, ft: &str) -> Gtable {
    let bias = ft == "MS_INDEL_ERROR_BIAS";
    let title = if bias { None } else { Some("Microsatellite indel errors") };
    let subtitle = if bias { "Bias" } else { "Rates" };
    let ylab = if bias { "Phred score diff." } else { "Phred score" };
    let Some(d) = inp.plot_data(ft) else {
        // plot_missing_data(plot_labels): the subtitle at the theme's size
        let mut spec = missing_spec(th, title, "Repeat units", ylab);
        spec.subtitle = Some((subtitle.to_string(), BASE));
        let mut gt = gg::table(th, spec);
        set_collect_keys(&mut gt, &[], &None, false, Some("Repeat units"), Some(ylab), None);
        return gt;
    };
    let cols = d.levels("RepeatUnitType");
    let rows_f = d.levels("ConsensusType");
    let st_levels = d.levels("SampleType");
    let mut spec = base_spec(th, title, Some("Repeat units"), Some(ylab));
    spec.subtitle = Some((subtitle.to_string(), BASE * 0.85));
    if bias {
        spec.grid_h = true;
    } else {
        spec.grid_h = true;
        spec.grid_v = true;
    }
    spec.ylab_r = Some("Consensus read type".into());
    spec.right_axis = true;
    spec.strips_top = Some(cols.clone());
    spec.strips_right = Some(rows_f.clone());
    spec.spacing_x = th.pt(-0.5);
    let mut all_x = Vec::new();
    let mut all_y = Vec::new();
    let mut grid: Vec<Vec<Layers>> = Vec::new();
    for rf in &rows_f {
        let mut row = Vec::new();
        for (ci, c) in cols.iter().enumerate() {
            let rows: Vec<PwRow> = d
                .rows
                .iter()
                .filter(|r| r.field("ConsensusType") == rf && r.field("RepeatUnitType") == c)
                .map(|r| PwRow { x: r.num("RefNumUnits"), st: r.sample_type.clone(), group: group_of(&st_levels, &r.sample_type), y: r.value, pct: r.pct })
                .collect();
            let hl = if bias { vec![0.0] } else { Vec::new() };
            let (xc, yc, layers) = pairwise_layers(rows, Pairwise::PointRange, inp.has_cohort, 1.0, hl, Vec::new(), 0, Vec::new());
            all_x.extend(xc);
            all_y.extend(yc);
            let layers: Layers = if bias && ci == 0 {
                // geom_text(x = -Inf, y = Inf / -Inf, vjust 1.5 / -0.5, hjust -0.05, size BASE * 0.22, grey40)
                Rc::new(move |cv: &mut Canvas, ctx: &Ctx, m: &Map| {
                    layers(cv, ctx, m);
                    let fs = BASE * 0.22 * PT;
                    let grey40 = col("grey40");
                    gg::text(cv, ctx, m.b.x0, m.b.y0, "More ins errors", fs, grey40, -0.05, 1.5);
                    gg::text(cv, ctx, m.b.x0, m.b.y1, "More del errors", fs, grey40, -0.05, -0.5);
                })
            } else {
                layers
            };
            row.push(layers);
        }
        grid.push(row);
    }
    // expand_limits(x = c(4, 15)); scale_x_continuous(breaks = breaks_width(3))
    all_x.extend([4.0, 15.0]);
    let h = continuous(&all_x, (None, None), (None, None), false, Breaks::Width(3.0), &default_labels, false);
    let v = if bias {
        // scale_y_continuous(expand = expansion(mult = 0.12))
        let finite: Vec<f64> = all_y.iter().copied().filter(|v| v.is_finite()).collect();
        let lim = if finite.is_empty() { (0.0, 1.0) } else { (finite.iter().copied().fold(f64::INFINITY, f64::min), finite.iter().copied().fold(f64::NEG_INFINITY, f64::max)) };
        let range = scales::expand(lim, 0.12, 0.0);
        let kept: Vec<f64> = rfmt::extended_breaks(range.0, range.1, 5).into_iter().filter(|&b| b >= range.0 && b <= range.1).collect();
        Pos { range, log: false, labels: plus_labels(&kept), breaks: kept, logticks: Vec::new() }
    } else {
        // scale_y_continuous(limits = c(0, NA)): values below 0 are dropped by the scale
        continuous(&all_y.iter().copied().filter(|v| *v >= 0.0).collect::<Vec<_>>(), (Some(0.0), None), (None, None), false, Breaks::Auto, &default_labels, false)
    };
    let hs = vec![h; cols.len()];
    spec.h = hs.clone();
    spec.v = vec![v; rows_f.len()];
    spec.layers = grid;
    let h_text = spec.h_text.clone();
    let mut gt = gg::table(th, spec);
    set_collect_keys(&mut gt, &hs, &h_text, true, Some("Repeat units"), Some(ylab), Some("Consensus read type"));
    gt
}

/// I (upstream 883ab11bba): wrap_plots(rates, bias, ncol = 1) + plot_layout(axes = "collect_x", axis_titles = "collect")
fn ms_indel_pair(th: &Theme, inp: &Inputs) -> Built {
    let p = vec![patchwork::plot(ms_indel(th, inp, "MS_INDEL_ERROR_RATES")), patchwork::plot(ms_indel(th, inp, "MS_INDEL_ERROR_BIAS"))];
    Built::Patchwork(patchwork::build_collect(p, &Layout::column(2, vec![None, None]), patchwork::Collect { axes_x: true, axes_y: false, titles_x: true, titles_y: true }))
}

// -------------------------------------------------------------------------------------------------------------
// the summary table (A)

struct SubTable {
    group: &'static str,
    format: &'static str,
    min_upper: f64,
    show_title: bool,
    show_sample_label: bool,
}

/// plot_sub_table: the values (labels coloured by sample type, boxed by QC status) beside the sample-vs-cohort
/// plot; a patchwork of the two.
fn sub_table(th: &Theme, d: &PlotData, s: &SubTable) -> Option<(Gtable, usize)> {
    let ctx = th.ctx;
    let rows: Vec<&Row> = d.rows.iter().filter(|r| r.field("FeatureGroup") == s.group && r.field("NumberFormat") == s.format).collect();
    // DisplayName and SampleType levels reversed (all levels of the summary's factors)
    let mut names: Vec<String> = d.levels("DisplayName");
    names.reverse();
    let mut st_levels: Vec<String> = d.levels("SampleType");
    st_levels.reverse();
    let present_names: Vec<String> = names.iter().filter(|n| rows.iter().any(|r| r.field("DisplayName") == **n)).cloned().collect();
    let n_rows = {
        let mut v: Vec<&str> = rows.iter().map(|r| r.field("DisplayName")).collect();
        v.sort();
        v.dedup();
        v.len()
    };
    let status = |r: &Row| -> &'static str {
        let q = r.field("QcStatus");
        if q.starts_with("FAIL") {
            "FAIL"
        } else if q.starts_with("WARN") {
            "WARN"
        } else {
            "PASS"
        }
    };
    let fmt_value = |v: f64| -> String {
        if s.format == "PERCENT" {
            scales::label_percent(&[v], Some(0.1), 100.0, false).remove(0)
        } else if v < 1.0 {
            rfmt::r_as_character(signif(v, 2))
        } else {
            rfmt::r_as_character(rfmt::r_round(v, 1))
        }
    };
    // the values plot: x = sample type (limits rev, drop = FALSE), y = display name (present levels)
    let x_limits: Vec<String> = st_levels.iter().rev().cloned().collect();
    let y_levels = present_names.clone();
    let div_lines: Vec<f64> = if n_rows >= 2 { (1..n_rows).map(|k| k as f64 + 0.5).collect() } else { vec![1.5, 0.5] };
    let div_col = if n_rows > 1 { col("grey90") } else { Rgba::WHITE };
    let labels: Vec<(f64, f64, String, String, String)> = rows
        .iter()
        .map(|r| {
            let x = (x_limits.iter().position(|l| *l == r.sample_type).unwrap() + 1) as f64;
            let y = (y_levels.iter().position(|l| l == r.field("DisplayName")).unwrap() + 1) as f64;
            (x, y, r.value.map(fmt_value).unwrap_or_else(|| "NA".into()), r.sample_type.clone(), status(r).to_string())
        })
        .collect();
    let mut yc: Vec<f64> = labels.iter().map(|l| l.1).collect();
    yc.extend(div_lines.iter().copied());
    let xc: Vec<f64> = labels.iter().map(|l| l.0).collect();
    // no rows (a group whose tool's output is missing): R's scales see only the div lines at seq(1.5, 0), so the
    // y scale is continuous (0.5 to 1.5, its labels drawn into the zero-height panel), the x scale untrained (no
    // breaks), and the dummy limit points of the comparison plot have no rows to train with
    let empty = rows.is_empty();
    let mut spec = base_spec(th, Some(s.group), None, None);
    spec.title_blank = !s.show_title;
    let names_map = |st: &str| if st == "NORMAL" { "Normal".to_string() } else { "Tumor".to_string() };
    if empty {
        spec.h = vec![discrete(0, Vec::new(), &[])];
        spec.v = vec![continuous(&yc, (None, None), (None, None), false, Breaks::Auto, &default_labels, false)];
    } else {
        spec.h = vec![discrete(x_limits.len(), x_limits.iter().map(|l| names_map(l)).collect(), &xc)];
        spec.v = vec![discrete(y_levels.len(), y_levels.clone(), &yc)];
    }
    spec.h_text = (s.show_sample_label && !empty).then(|| th.axis_text_h());
    spec.h_ticks = s.show_sample_label && !empty;
    let mut vt = th.axis_text_v();
    vt.size = BASE * 0.95;
    spec.v_text = Some(vt);
    spec.margin = [0.0, 0.0, th.pt(10.0), 0.0];
    let pad = th.pt(3.0);
    let dl = div_lines.clone();
    spec.layers = vec![vec![Rc::new(move |c: &mut Canvas, ctx: &Ctx, m: &Map| {
        for (x, y, t, st, q) in &labels {
            let border = if q == "PASS" { Rgba::TRANSPARENT } else { Rgba::BLACK };
            label(c, ctx, m.at(*x, *y), t, BASE * 0.3, type_colour(st), qc_colour(q), border, LINE * 0.5, pad);
        }
        for &h in &dl {
            let y = m.b.npc_y(m.v.npc(h));
            path(c, ctx, vec![(m.b.x0, y), (m.b.x1, y)], div_col, 0.21, "");
        }
    }) as Layers]];
    let values = gg::table(th, spec);
    let _ = ctx;

    // the sample-vs-cohort plot (coord_flip): x = display name, y = value; box width scaled by the sample types
    // present over the levels
    let present_types = {
        let mut v: Vec<&str> = rows.iter().map(|r| r.sample_type.as_str()).collect();
        v.sort();
        v.dedup();
        v.len()
    };
    let box_scale = present_types as f64 / st_levels.len() as f64;
    let pw: Vec<PwRow> = rows
        .iter()
        .map(|r| PwRow { x: (y_levels.iter().position(|l| l == r.field("DisplayName")).unwrap() + 1) as f64, st: r.sample_type.clone(), group: group_of(&st_levels, &r.sample_type), y: r.value, pct: r.pct })
        .collect();
    let log = s.format == "LOG10";
    let lower = if log { 1.0 } else { 0.0 };
    let dummies = if empty { Vec::new() } else { vec![(1.0, lower), (1.0, s.min_upper)] };
    let (xc, yc, layers) = pairwise_layers(pw, Pairwise::Box, th_has_cohort(d), box_scale, Vec::new(), div_lines.clone(), y_levels.len(), dummies);
    let mut spec = base_spec(th, None, None, None);
    spec.flip = true;
    if empty {
        // an untrained continuous scale: limits c(0, 1), no breaks, no axis (in R 4.0.3 the log scale's
        // guide_axis_logticks fails here under coord_flip: QS14)
        spec.v = vec![continuous(&xc, (None, None), (None, None), false, Breaks::Auto, &default_labels, false)];
        let mut h = continuous(&[], (None, None), (None, None), false, Breaks::Auto, &default_labels, false);
        h.breaks.clear();
        h.labels.clear();
        spec.h = vec![h];
        spec.h_text = None;
        spec.h_ticks = false;
    } else {
        spec.v = vec![discrete(y_levels.len(), y_levels.clone(), &xc)];
        spec.h = vec![match s.format {
            "PERCENT" => continuous(&yc, (None, None), (None, None), false, Breaks::Auto, &percent, false),
            "LOG10" => continuous(&yc, (None, None), (None, None), true, Breaks::Auto, &plain_labels, true),
            _ => continuous(&yc, (None, None), (None, None), false, Breaks::Auto, &default_labels, false),
        }];
    }
    spec.v_text = None;
    spec.v_ticks = false;
    spec.spacing_x = th.pt(-0.5);
    spec.margin = [0.0, th.pt(15.0), th.pt(10.0), th.pt(2.0)];
    spec.layers = vec![vec![layers]];
    let cmp = gg::table(th, spec);
    let both = patchwork::build(vec![patchwork::plot(values), patchwork::plot(cmp)], &Layout::row(2));
    Some((both, n_rows))
}

fn signif(x: f64, digits: i32) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let e = x.abs().log10().floor() as i32;
    rfmt::r_round(x, digits - 1 - e)
}

/// Whether R fails to draw the report (QS14): ggplot2 4.0.3's guide_axis_logticks cannot transform the untrained log
/// scale of an empty LOG10 sub-table (Copy number or Mutational burden, both from PURPLE) under coord_flip.
fn r_report_fails(inp: &Inputs) -> bool {
    let Some(d) = inp.plot_data("SUMMARY_TABLE") else { return false };
    [("Copy number", "LOG10"), ("Mutational burden", "LOG10")]
        .iter()
        .any(|(g, f)| !d.rows.iter().any(|r| r.field("FeatureGroup") == *g && r.field("NumberFormat") == *f))
}

/// The summary table: the seven sub-tables stacked, heights by their rows.
fn summary(th: &Theme, inp: &Inputs) -> Option<Gtable> {
    let d = inp.plot_data("SUMMARY_TABLE")?;
    let subs = [
        SubTable { group: "Mapping", format: "NUMBER", min_upper: 100.0, show_title: true, show_sample_label: false },
        SubTable { group: "Mapping", format: "PERCENT", min_upper: 1.0, show_title: false, show_sample_label: false },
        SubTable { group: "Copy number", format: "NUMBER", min_upper: 4.0, show_title: true, show_sample_label: false },
        SubTable { group: "Copy number", format: "PERCENT", min_upper: 1.0, show_title: false, show_sample_label: false },
        SubTable { group: "Copy number", format: "LOG10", min_upper: 1000.0, show_title: false, show_sample_label: false },
        SubTable { group: "Contamination", format: "PERCENT", min_upper: 1.0, show_title: true, show_sample_label: false },
        SubTable { group: "Mutational burden", format: "LOG10", min_upper: 1000.0, show_title: true, show_sample_label: true },
    ];
    let mut patches = Vec::new();
    let mut heights = Vec::new();
    for s in &subs {
        if let Some((gt, n)) = sub_table(th, &d, s) {
            patches.push(patchwork::nested(gt));
            heights.push(Some(n as f64));
        }
    }
    Some(patchwork::build(patches, &Layout::column(heights.len(), heights)))
}

// -------------------------------------------------------------------------------------------------------------
// the legend (Z) and the title

/// create_legend(show_qc_status) (upstream 883ab11bba): the sample types, the read strands (the BQR plots' triangles)
/// and, on page 1, the QC statuses, in one guide box at the legend's default position (right: the legends stacked,
/// left-aligned, the box centred), placed over the whole patch (wrap_elements(full = )).
fn legend_patch(th: &Theme, show_qc_status: bool) -> Patch {
    let ctx = th.ctx;
    let key = ctx.dev.bigpt(1.2 * 12.0 * 1.2);
    let st_entries: Vec<(String, KeyGlyph)> = [("Tumor sample", TUMOR_COL, 1.0), ("Tumor cohort", TUMOR_COL, 0.2), ("Normal sample", NORMAL_COL, 1.0), ("Normal cohort", NORMAL_COL, 0.2)]
        .iter()
        .map(|(l, c, a)| {
            let colour = Rgba { a: if *a < 1.0 { 0.2 } else { 1.0 }, ..col(c) };
            (l.to_string(), KeyGlyph::Point { shape: 15, colour, fill: colour, size: BASE * 0.5, stroke: LINE })
        })
        .collect();
    // guide_legend(override.aes = list(size = BASE * 0.3, fill = "grey50")): black outlines, the default stroke
    let rs_entries: Vec<(String, KeyGlyph)> =
        [("Forward", 24u8), ("Reverse", 25u8)].iter().map(|(l, sh)| (l.to_string(), KeyGlyph::Point { shape: *sh, colour: Rgba::BLACK, fill: col("grey50"), size: BASE * 0.3, stroke: 0.5 })).collect();
    let qc_entries: Vec<(String, KeyGlyph)> = ["FAIL", "WARN", "PASS"].iter().map(|q| (q.to_string(), KeyGlyph::Label { fill: qc_colour(q), border: Rgba::BLACK, size: BASE * 0.35, linewidth: LINE * 0.5 })).collect();
    let mut ttl = TextEl::new(BASE);
    ttl.hjust = 0.0;
    let mut txt = TextEl::new(0.8 * BASE);
    txt.hjust = 0.0;
    let lth = LegendTheme {
        title: ttl,
        text: txt,
        key_spacing: th.pt(HALF),
        margin: th.margin(),
        background: Some(RectEl { fill: Some(Rgba::WHITE), colour: None, linewidth: LINE }),
        key_background: Some(RectEl { fill: Some(Rgba::WHITE), colour: None, linewidth: LINE }),
    };
    let point_key = |size: f64| size / 10.0 * ctx.dev.cm(1.0);
    // plot_qc.R builds the legends (ggplot_gtable) at the top level, where pdf(NULL) is the open device: their sizes
    // are pdf()'s even in the PNG, where they are drawn
    let mctx = match ctx.dev {
        Device::Png { dpi } => Ctx { dev: Device::PdfAt { dpi }, v: ctx.v },
        _ => ctx,
    };
    let kst = key.max(point_key(BASE * 0.5));
    let krs = key.max(point_key(BASE * 0.3));
    let mut legends = vec![
        guides::legend_measured(ctx, mctx, &LegendSpec::Keys { title: "Sample type".into(), entries: st_entries, key_width: kst, key_height: kst }, &lth),
        guides::legend_measured(ctx, mctx, &LegendSpec::Keys { title: "Read strand".into(), entries: rs_entries, key_width: krs, key_height: krs }, &lth),
    ];
    if show_qc_status {
        legends.push(guides::legend_measured(ctx, mctx, &LegendSpec::Keys { title: "QC status".into(), entries: qc_entries, key_width: key, key_height: key }, &lth));
    }
    // legend.position "right": legend.box vertical, legend.box.just c("left", "top"), legend.justification centre
    let gb = guides::guide_box(ctx, GuideBox { legends, horizontal: false, just: (0.5, 0.5), box_just: (0.0, 1.0), spacing: th.pt(2.0 * HALF), box_margin: [0.0; 4] });
    let m = th.margin();
    patchwork::wrap_full(Grob::table(gb), m, Some(rect_grob(ctx, RectEl { fill: Some(Rgba::WHITE), colour: Some(Rgba::WHITE), linewidth: LINE })))
}

/// REPORT_TITLE: "<tumor>: <most severe status>" and the WARN and FAIL features per sample type.
fn title_texts(inp: &Inputs, tumor: &str) -> (String, Option<String>) {
    let Some(d) = inp.plot_data("SUMMARY_TABLE") else {
        return (format!("{tumor}: PASS"), None);
    };
    let status = |r: &Row| -> &'static str {
        let q = r.field("QcStatus");
        if q.starts_with("FAIL") {
            "FAIL"
        } else if q.starts_with("WARN") {
            "WARN"
        } else {
            "PASS"
        }
    };
    let mut strings = Vec::new();
    for (st, name) in [("TUMOR", "Tumor"), ("NORMAL", "Normal")] {
        let mut parts = Vec::new();
        for q in ["FAIL", "WARN"] {
            let items: Vec<String> = d.rows.iter().filter(|r| r.sample_type == st && status(r) == q).map(|r| format!("{} ({})", r.field("DisplayName"), r.field("QcThreshold"))).collect();
            if !items.is_empty() {
                parts.push(format!("{q}: {}", items.join(", ")));
            }
        }
        if !parts.is_empty() {
            strings.push(format!("{name} - {}", parts.join(" - ")));
        }
    }
    let worst = ["FAIL", "WARN", "PASS"].into_iter().find(|q| d.rows.iter().any(|r| status(r) == *q)).unwrap_or("PASS");
    (format!("{tumor}: {worst}"), if strings.is_empty() { None } else { Some(strings.join("\n")) })
}

// -------------------------------------------------------------------------------------------------------------
// the report

/// One page of the report (upstream 883ab11bba: two pages). Page 1: the summary table, coverage, fragment length,
/// GC bias, duplicate frequency, discordant fragments, missed variants, BQR by original base quality, the MS indel
/// errors and the legend in a design of 4 x 9 cells (columns 5 and 6 narrower); page 2: BQR by SNV96 context and the
/// legend without the QC statuses. Both with the report title.
fn report_table(ctx: Ctx, inp: &Inputs, tumor: &str, page: usize) -> Gtable {
    let th = Theme { ctx };
    let free_label = |gt: Gtable, nested: bool| -> Patch {
        let mut p = if nested { patchwork::nested(gt) } else { patchwork::plot(gt) };
        patchwork::free_label(&mut p.gt, [true; 4]);
        p
    };
    let free_panel_l = |gt: Gtable| -> Patch {
        let mut p = patchwork::plot(gt);
        patchwork::free_panel(&mut p.gt, [false, false, false, true]);
        p
    };
    let gt = if page == 1 {
        // free(<patchwork>) ends with merge(), which nests the free summary table once more (merge.patchwork); free() of
        // a ggplot (the missing-data plot) returns it as it is (merge.ggplot)
        let a = match summary(&th, inp) {
            Some(t) => patchwork::merge(free_label(t, true)),
            None => free_label(missing(&th, "", "", ""), false),
        };
        // in the design's order (A to I, then Z): wrap_plots matches the plots' names to the design's letters
        let patches = vec![
            a,
            free_label(distribution_plot(&th, inp, "COVERAGE_DISTRIBUTION", "Coverage", "Coverage", "Bases"), false),
            free_label(distribution_plot(&th, inp, "FRAG_LENGTH_DISTRIBUTION", "Fragment length", "Fragment length", "Fragments"), false),
            free_label(distribution_plot(&th, inp, "DUPLICATE_FREQ", "Duplicate frequency", "Duplicate read count", "Read groups"), false),
            free_label(distribution_plot(&th, inp, "GC_BIAS", "GC bias", "GC proportion", "Normalised median read depth"), false),
            free_panel_l(discordant(&th, inp)),
            free_panel_l(missed_variants(&th, inp)),
            // H and I are placed as they are (not freed)
            patch_of(bqr(&th, inp, "BQR_PER_ORIG_QUAL")),
            patch_of(ms_indel_pair(&th, inp)),
            {
                let mut p = legend_patch(&th, true);
                patchwork::free_label(&mut p.gt, [true; 4]);
                p
            },
        ];
        let mut layout = Layout::design("AABBEEDDZ\nAACCFFGGZ\nAAHHHIIIZ\nAAHHHIIIZ");
        layout.widths = [1.0, 1.0, 1.0, 1.0, 0.74, 0.74, 1.0, 1.0, 0.75].iter().map(|&w| Some(w)).collect();
        patchwork::build(patches, &layout)
    } else {
        let patches = vec![free_label_of(bqr(&th, inp, "BQR_PER_SNV96_CONTEXT")), legend_patch(&th, false)];
        let mut layout = Layout::row(2);
        layout.widths = vec![Some(1.0), Some(0.1)];
        patchwork::build(patches, &layout)
    };
    let (title, subtitle) = title_texts(inp, tumor);
    let mut tel = TextEl::new(BASE * 1.4);
    tel.face = Face::Bold;
    tel.hjust = 0.0;
    tel.vjust = 1.0;
    tel.margin = [0.0, 0.0, th.pt(HALF), 0.0];
    let mut sel = TextEl::new(BASE * 0.9);
    sel.hjust = 0.0;
    sel.vjust = 1.0;
    sel.margin = [0.0, 0.0, th.pt(HALF), 0.0];
    let t = TitleGrob::new(&title, &tel, true, true);
    let th_ = t.height(&ctx).abs();
    let sub = subtitle.map(|s| {
        let g = TitleGrob::new(&s, &sel, true, true);
        let h = g.height(&ctx).abs();
        (g.grob(ctx), h)
    });
    let m = th.margin();
    patchwork::annotate(gt, Annotation { margin: m, title: Some((t.grob(ctx), th_)), subtitle: sub.or(Some((Grob::Zero, 0.0))), background: Some(rect_grob(ctx, RectEl { fill: Some(Rgba::WHITE), colour: Some(Rgba::WHITE), linewidth: LINE })) })
}

/// The PNG of a page (upstream 883ab11bba): `<report>.page_<n>.png` beside the PDF (`sub("pdf$", "page_n.png", ...)`)
pub fn page_png_path(pdf_path: &str, page: usize) -> String {
    format!("{}page_{page}.png", pdf_path.strip_suffix("pdf").unwrap_or(pdf_path))
}

/// Writes the report: the PDF (two pages), and for a single patient a PNG per page (17 x 10 inches; 300 dpi); in the
/// default mode also page 1 as `<report>.png`, the name ORANGE reads (QS17). Returns the PDF's pages (for
/// -merge_plots).
pub fn report(tumor: &str, reference: Option<&str>, vis_path: &str, cohort: Option<&str>, pdf_path: &str, single_patient: bool) -> Result<Vec<Figure>> {
    let inp = Inputs::read(vis_path, cohort, tumor, reference)?;
    let (w_in, h_in) = (17.0, 10.0);
    hmf_common::info!("Writing PDF: {pdf_path}");
    if crate::compat::java("QS14") && r_report_fails(&inp) {
        // R opens the PDF, fails while building the plots ("attempt to set an attribute on NULL") and leaves one blank
        // page; no PNG; Java's QseePlotTask then throws
        let (w, h) = (w_in * 72.0, h_in * 72.0);
        let fig = Figure::new(w as u32, h as u32, Rgba::WHITE);
        hmf_plot::pdf_figure::write_r_pdf(pdf_path, &[fig], w, h).map_err(|e| anyhow::anyhow!("writing {pdf_path}: {e}"))?;
        anyhow::bail!("R: Error in names(object) <- nm : attempt to set an attribute on NULL (an empty LOG10 summary sub-table)");
    }
    let mut pages = Vec::new();
    for page in [1, 2] {
        let ctx = Ctx { dev: Device::Pdf, v: Version::V403 };
        let t = report_table(ctx, &inp, tumor, page);
        let (w, h) = (w_in * 72.0, h_in * 72.0);
        let mut fig = Figure::new(w as u32, h as u32, Rgba::WHITE);
        let mut c = Canvas::new(&mut fig);
        hmf_plot::ggp::record_drawn_boxes(boxes_wanted());
        t.draw(&mut c, Bx::new(0.0, 0.0, w, h));
        dump_boxes(&t, w, h, &format!("page_{page}.pdf"));
        pages.push(fig);
    }
    hmf_plot::pdf_figure::write_r_pdf(pdf_path, &pages, w_in * 72.0, h_in * 72.0).map_err(|e| anyhow::anyhow!("writing {pdf_path}: {e}"))?;
    if single_patient {
        for page in [1, 2] {
            let png = page_png_path(pdf_path, page);
            hmf_common::info!("Writing PNG: {png}");
            let ctx = Ctx { dev: Device::Png { dpi: 300.0 }, v: Version::V403 };
            let t = report_table(ctx, &inp, tumor, page);
            let (w, h) = ((w_in * 300.0) as u32, (h_in * 300.0) as u32);
            let mut fig = Figure::new(w, h, Rgba::WHITE);
            fig.meta = PngMeta { dpi: Some(300.0), transparent_key: None };
            let mut c = Canvas::new(&mut fig);
            hmf_plot::ggp::record_drawn_boxes(boxes_wanted());
            t.draw(&mut c, Bx::new(0.0, 0.0, w as f64, h as f64));
            dump_boxes(&t, w as f64, h as f64, &format!("page_{page}.png"));
            fig.write_png(&png).with_context(|| format!("writing {png}"))?;
            // QS17: the report's PNG under the name ORANGE (MiscToolFiles.QSEE_PLOT_FILE_ID) looks for, page 1
            if page == 1 && !crate::compat::java("QS17") {
                let legacy = format!("{}png", pdf_path.strip_suffix("pdf").unwrap_or(pdf_path));
                hmf_common::info!("Writing PNG: {legacy}");
                std::fs::copy(&png, &legacy).with_context(|| format!("writing {legacy}"))?;
            }
        }
    }
    let _ = (Item::Clip(None), Stroke::solid(Rgba::BLACK, 1.0), STROKE);
    Ok(pages)
}

/// QseePlot.mergePlots: the reports' pages in one PDF (`multisample.qsee.vis.report.pdf`), in sample order.
pub fn write_merged(path: &str, pages: &[Figure]) -> Result<()> {
    hmf_plot::pdf_figure::write_r_pdf(path, pages, 17.0 * 72.0, 10.0 * 72.0).map_err(|e| anyhow::anyhow!("writing {path}: {e}"))
}

fn boxes_wanted() -> bool {
    std::env::var_os("QSEE_PLOT_BOXES").is_some()
}

/// With QSEE_PLOT_BOXES=<prefix>, the layout's boxes (for comparing with R's viewports), and after them the boxes of
/// the label rectangles drawn (R's "rrvp" viewports; recorded while drawing).
fn dump_boxes(t: &Gtable, w: f64, h: f64, kind: &str) {
    let drawn = hmf_plot::ggp::take_drawn_boxes();
    if let Ok(p) = std::env::var("QSEE_PLOT_BOXES") {
        let mut out = Vec::new();
        t.boxes(Bx::new(0.0, 0.0, w, h), "", &mut out);
        out.extend(drawn);
        let text: String = out.iter().map(|(n, b)| format!("{:.4}\t{:.4}\t{:.4}\t{:.4}\t{n}\n", b.x0, b.y0, b.x1, b.y1)).collect();
        let _ = std::fs::write(format!("{p}.{kind}.tsv"), text);
        let f = |u: &[hmf_plot::gtable::Unit]| u.iter().map(|x| format!("{:.4}", x.abs())).collect::<Vec<_>>().join(" ");
        let _ = std::fs::write(format!("{p}.{kind}.sizes.txt"), format!("{}\n{}\n", f(&t.widths), f(&t.heights)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// upstream 883ab11bba: a PNG per page, `sub("pdf$", "page_<n>.png", <pdf>)`
    #[test]
    fn page_png_names() {
        assert_eq!(page_png_path("/out/T.qsee.vis.report.pdf", 1), "/out/T.qsee.vis.report.page_1.png");
        assert_eq!(page_png_path("T.qsee.vis.report.x.pdf", 2), "T.qsee.vis.report.x.page_2.png");
    }

    /// The BQR plot's quality bins in the order of their lower bounds, and position_dodge(0.6) of the read strands
    #[test]
    fn bqr_bins_and_dodge() {
        assert_eq!(bin_lower_bound("HIGH (30+)"), 30.0);
        assert_eq!(bin_lower_bound("LOW (0-29)"), 0.0);
        assert_eq!(bin_lower_bound("MEDIUM (19-29)"), 19.0);
        let row = |x: f64, group: usize| BqrRow { x, group, st: "TUMOR".into(), orientation: "FORWARD".into(), y: Some(1.0), pct: [None; 5] };
        let (shift, ns) = dodge_total(&[row(1.0, 1), row(1.0, 0), row(2.0, 1)]);
        assert_eq!(ns, vec![2, 2, 1]);
        assert!((shift[0] - 0.15).abs() < 1e-12 && (shift[1] + 0.15).abs() < 1e-12 && shift[2] == 0.0);
        assert_eq!((orientation_shape("FORWARD"), orientation_shape("REVERSE")), (24, 25));
    }

    /// geom_label in the summary table (ggplot2 4.0.3, png() at 300 dpi): the box is centred on its point (R's "rrvp"
    /// viewport has the centre of the label's anchor viewport), 51.9066 px high whatever the one-line text (R's boxes
    /// in the public HG008 report: e.g. (535.2227, 447.0872)-(679.1293, 498.9938) around the anchor (607.1760,
    /// 473.0405)), and the text's width plus the padding wide.
    #[test]
    fn label_box_centred_on_its_point_as_r() {
        let ctx = Ctx { dev: Device::Png { dpi: 300.0 }, v: Version::V403 };
        let th = Theme { ctx };
        let mut fig = Figure::new(8, 8, Rgba::WHITE);
        let mut c = Canvas::new(&mut fig);
        let pad = th.pt(3.0);
        hmf_plot::ggp::record_drawn_boxes(true);
        label(&mut c, &ctx, (607.1760, 473.0405), "0.0%", BASE * 0.3, Rgba::BLACK, Rgba::WHITE, Rgba::BLACK, LINE * 0.5, pad);
        let b = hmf_plot::ggp::take_drawn_boxes()[0].1;
        assert!(((b.y0 + b.y1) / 2.0 - 473.0405).abs() < 1e-9);
        assert!(((b.x0 + b.x1) / 2.0 - 607.1760).abs() < 1e-9);
        assert!((b.h() - (498.9938 - 447.0872)).abs() < 0.001, "height {}", b.h());
        let w = ctx.dev.str_width("0.0%", BASE * 0.3 * hmf_plot::ggp::PT, Face::Plain);
        assert!((b.w() - (w + 2.0 * pad)).abs() < 1e-9);
    }
}
