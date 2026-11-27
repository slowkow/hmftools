//! PURPLE's plots, drawn natively with hmf-plot-rs: the two Circos images (`circos.png`, `input.png`) from
//! the circos/ files, and the seven plots of the R scripts `copyNumberPlots.R` and `somaticVariantPlots.R`
//! from the same files the scripts read, with the same transformations (dplyr steps, ggplot2 stats, scales).
//!
//! Findings (corrected by default, Java's drawing with `-java_compat`):
//! - PL1: the Circos SNV track's `skip_run = yes` drops every point whose allele frequency text equals the
//!   previous point's.
//! - PL2: `size = 0.07` on geom_histogram and geom_bar is ignored by ggplot2 >= 3.4 (it wants `linewidth`), so
//!   the bars get 0.5 mm outlines.
//! - PL3: the rainfall plot's intermutation distance is taken across chromosome boundaries.
//! - PL4: the rainfall data stop after 100,000 PASS variants of any type, not 100,000 SNVs.

use std::collections::{BTreeMap, HashMap};

use anyhow::{Context, Result, bail};
use hmf_plot::color::{GradientN, Rgba, hue_pal, r_colour};
use hmf_plot::gg::{self, Key, Labels, Legend, Panel, Plot, PosScale, Trans, rfmt};
use hmf_plot::Figure;

use crate::compat::java;

/// One SNV of the rainfall plot (plot/<tumor>.somatic_data.tsv in Java).
#[derive(Clone, Debug)]
pub struct RainfallPoint {
    pub chromosome: String,
    pub position: i64,
    pub mutation: String,
    pub kataegis: String,
}

/// A tab-separated file with a header, as R's read.table sees it.
struct Table {
    cols: HashMap<String, usize>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn read(path: &str) -> Result<Table> {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
        let mut lines = text.lines();
        let header = lines.next().unwrap_or("");
        let cols = header.split('\t').enumerate().map(|(i, c)| (c.to_string(), i)).collect();
        let rows = lines.filter(|l| !l.is_empty()).map(|l| l.split('\t').map(|s| s.to_string()).collect()).collect();
        Ok(Table { cols, rows })
    }

    fn col(&self, name: &str) -> Result<usize> {
        self.cols.get(name).copied().with_context(|| format!("no column {name}"))
    }

    fn num(&self, row: &[String], c: usize) -> f64 {
        row.get(c).and_then(|v| v.parse::<f64>().ok()).unwrap_or(f64::NAN)
    }
}

/// R's `round()` to an integer, then ggplot2's paste0 of it.
fn r_round0(x: f64) -> f64 {
    rfmt::r_round(x, 0)
}

fn paste_int(x: f64) -> String {
    rfmt::r_as_character(if x == 0.0 { 0.0 } else { x })
}

fn hex(s: &str) -> Rgba {
    r_colour(s).unwrap_or(Rgba::BLACK)
}

/// The outline width of bars drawn with `size = 0.07` (PL2).
fn thin_bar_linewidth() -> f64 {
    if java("PL2") { 0.5 } else { 0.07 }
}

fn write(fig: &Figure, path: &str) -> Result<()> {
    fig.write_png(path).with_context(|| format!("writing {path}"))
}

// -------------------------------------------------------------------------------------------------------------
// stacked histograms and bars (stat_bin / stat_count with position = "stack")

/// One bar of a stacked bar chart, in data units, with its fill level.
#[derive(Clone, Debug)]
struct Bar {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    group: usize,
}

/// position_stack: within each x, the groups stacked from the last level at the bottom to the first at the top;
/// rows come out ordered by x, then by group descending (the order ggplot2 draws them).
fn stack(mut bars: Vec<(f64, f64, f64, usize)>) -> Vec<Bar> {
    // (xmin, xmax, y, group)
    bars.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.3.cmp(&a.3)));
    let mut out = Vec::with_capacity(bars.len());
    let mut i = 0;
    while i < bars.len() {
        let mut j = i;
        let mut h = 0.0f64;
        while j < bars.len() && bars[j].0 == bars[i].0 {
            let y = bars[j].2;
            out.push(Bar { xmin: bars[j].0, xmax: bars[j].1, ymin: h.min(h + y), ymax: h.max(h + y), group: bars[j].3 });
            h += y;
            j += 1;
        }
        i = j;
    }
    out
}

/// stat_bin (binwidth, weights) per fill level over the scale limits, stacked.
fn binned(values: &[(f64, f64, usize)], limits: (f64, f64), width: f64, levels: usize) -> Vec<Bar> {
    let breaks = gg::bin_breaks(limits, width);
    let nbins = breaks.len() - 1;
    let mut counts = vec![vec![0.0f64; nbins]; levels];
    let mut present = vec![false; levels];
    for &(x, w, g) in values {
        // continuous scale limits censor out-of-range values before the stat
        if !(x >= limits.0 && x <= limits.1) {
            continue;
        }
        if let Some(b) = gg::bin_index(&breaks, x) {
            counts[g][b] += w;
            present[g] = true;
        }
    }
    let mut rows = Vec::new();
    for (g, c) in counts.iter().enumerate() {
        if !present[g] {
            continue;
        }
        for (b, &n) in c.iter().enumerate() {
            rows.push((breaks[b], breaks[b + 1], n, g));
        }
    }
    stack(rows)
}

fn y_limits(bars: &[Bar]) -> (f64, f64) {
    let lo = bars.iter().map(|b| b.ymin).fold(f64::INFINITY, f64::min);
    let hi = bars.iter().map(|b| b.ymax).fold(f64::NEG_INFINITY, f64::max);
    if lo.is_finite() { (lo, hi) } else { (0.0, 1.0) }
}

/// dplyr: bucket = ceiling(v); per bucket sum the weights; cumulative proportion; the first bucket (in order
/// of proportion) whose proportion exceeds `q`.
fn quantile_bucket(rows: &[(f64, f64)], q: f64) -> Option<f64> {
    let mut sums: BTreeMap<i64, f64> = BTreeMap::new();
    for &(v, w) in rows {
        *sums.entry(v.ceil() as i64).or_insert(0.0) += w;
    }
    let mut cum = 0.0;
    let cums: Vec<(i64, f64)> = sums.iter().map(|(&b, &n)| {
        cum += n;
        (b, cum)
    }).collect();
    let max = cums.iter().map(|c| c.1).fold(f64::NEG_INFINITY, f64::max);
    let mut props: Vec<(i64, f64)> = cums.iter().map(|&(b, c)| (b, c / max)).collect();
    props.sort_by(|a, b| a.1.total_cmp(&b.1)); // arrange(proportion), stable
    props.into_iter().find(|p| p.1 > q).map(|p| p.0 as f64)
}

/// The bars, legend and scales of a stacked histogram plot (copynumber.png, map.png, somatic.png).
#[allow(clippy::too_many_arguments)]
fn stacked_plot(
    path: &str,
    title: &str,
    x_title: &str,
    y_title: &str,
    bars: &[Bar],
    levels: &[String],
    colours: &[(&str, &str)],
    x_limits: (f64, f64),
    x_breaks: Vec<f64>,
    outline_width: f64,
) -> Result<()> {
    let colour_of = |level: &str| colours.iter().find(|(l, _)| *l == level).map(|(_, c)| hex(c)).unwrap_or(hex("grey50"));
    // the y scale is trained on every bar; then the x scale's limits censor xmin and xmax too, removing the
    // bars that reach past them
    let ys = PosScale::new(y_limits(bars));
    let bars: Vec<Bar> = bars.iter().filter(|b| b.xmin >= x_limits.0 && b.xmax <= x_limits.1).cloned().collect();
    let bars = &bars;
    let mut xs = PosScale::new(x_limits);
    xs.breaks = Some(x_breaks);
    let present: Vec<usize> = (0..levels.len()).filter(|g| bars.iter().any(|b| b.group == *g)).collect();
    let legend = Legend::Keys {
        title: None,
        entries: present
            .iter()
            .map(|&g| (levels[g].clone(), Key::Polygon { fill: colour_of(&levels[g]), colour: Some(Rgba::BLACK), linewidth: outline_width }))
            .collect(),
    };
    let plot = Plot {
        title: Some(title.into()),
        x: xs.view(),
        y: ys.view(),
        x_title: Some(x_title.into()),
        y_title: Some(y_title.into()),
        legends: if present.is_empty() { vec![] } else { vec![legend] },
        ticks: false,
        minor_grid: false,
        panel_border: true,
    };
    let mut fig = gg::figure(4.8, 4.0);
    let (w, h) = (fig.width as f64, fig.height as f64);
    plot.draw(&mut fig, 0.0, 0.0, w, h, None, &mut |f, p| {
        for b in bars {
            gg::rect(f, p, b.xmin, b.xmax, b.ymin, b.ymax, Some(colour_of(&levels[b.group])), Some((Rgba::BLACK, outline_width)));
        }
    });
    write(&fig, path)
}

// -------------------------------------------------------------------------------------------------------------
// copyNumberPlots.R

struct CnvRow {
    copy_number: f64,
    baf_count: f64,
    minor: f64,
}

fn read_cnv(path: &str) -> Result<Vec<CnvRow>> {
    let t = Table::read(path)?;
    let (c_chr, c_cn, c_baf, c_minor) = (t.col("chromosome")?, t.col("copyNumber")?, t.col("bafCount")?, t.col("minorAlleleCopyNumber")?);
    Ok(t.rows
        .iter()
        .filter(|r| {
            let chr = r[c_chr].replace("chr", "");
            chr != "X" && chr != "Y"
        })
        .map(|r| CnvRow { copy_number: t.num(r, c_cn), baf_count: t.num(r, c_baf), minor: t.num(r, c_minor) })
        .filter(|r| r.baf_count > 0.0)
        .collect())
}

/// Levels of a character fill column: the sorted distinct values (R's factor levels).
fn levels_of(values: &[String]) -> Vec<String> {
    let mut l: Vec<String> = values.to_vec();
    l.sort();
    l.dedup();
    l
}

const CN_COLOURS: [(&str, &str); 7] = [("CN0", "#8dd3c7"), ("CN1", "#ffffb3"), ("CN2", "#bebada"), ("CN3", "#fb8072"), ("CN4", "#80b1d3"), ("CN5", "#fdb462"), ("CN6+", "#b3de69")];
const MAP_COLOURS: [(&str, &str); 6] = [("MACN0", "#8dd3c7"), ("MACN1", "#ffffb3"), ("MACN2", "#bebada"), ("MACN3", "#fb8072"), ("MACN4", "#80b1d3"), ("MACN5+", "#fdb462")];

/// copynumber_pdf: the copy number histogram weighted by BAF count, filled by minor allele copy number.
fn copynumber_png(rows: &[CnvRow], path: &str) -> Result<()> {
    let maps: Vec<String> = rows
        .iter()
        .map(|r| {
            let m = r_round0(r.minor);
            if m >= 5.0 { "MACN5+".to_string() } else { format!("MACN{}", paste_int(m)) }
        })
        .collect();
    let levels = levels_of(&maps);
    let max_cn = quantile_bucket(&rows.iter().map(|r| (r.copy_number, r.baf_count)).collect::<Vec<_>>(), 0.9).context("no copy number quantile")?;
    let min_cn = rows.iter().map(|r| r.copy_number).fold(f64::INFINITY, f64::min).floor();
    let limits = (min_cn - 0.1, max_cn + 0.1);
    let values: Vec<(f64, f64, usize)> = rows.iter().zip(&maps).map(|(r, m)| (r.copy_number, r.baf_count, levels.iter().position(|l| l == m).unwrap())).collect();
    let bars = binned(&values, limits, 0.1, levels.len());
    let breaks: Vec<f64> = ((min_cn - 1.0) as i64..=(max_cn + 1.0) as i64).map(|v| v as f64).collect();
    stacked_plot(path, "Copy Number PDF", "Copy Number", "Baf Count", &bars, &levels, &MAP_COLOURS, limits, breaks, thin_bar_linewidth())
}

/// minor_allele_ploidy_pdf: the minor allele copy number histogram, filled by copy number.
fn map_png(rows: &[CnvRow], path: &str) -> Result<()> {
    let cns: Vec<String> = rows
        .iter()
        .map(|r| {
            let cn = r_round0(r.copy_number).max(0.0);
            if cn >= 6.0 { "CN6+".to_string() } else { format!("CN{}", paste_int(cn)) }
        })
        .collect();
    let levels = levels_of(&cns);
    let max_map = quantile_bucket(&rows.iter().map(|r| (r.minor, r.baf_count)).collect::<Vec<_>>(), 0.9).context("no minor allele quantile")?;
    let limits = (-0.1, max_map + 0.1);
    let values: Vec<(f64, f64, usize)> = rows.iter().zip(&cns).map(|(r, c)| (r.minor, r.baf_count, levels.iter().position(|l| l == c).unwrap())).collect();
    let bars = binned(&values, limits, 0.1, levels.len());
    let breaks: Vec<f64> = (0..=10).map(|v| v as f64).collect();
    stacked_plot(path, "Minor Allele Copy Number PDF", "Minor Allele Copy Number", "Baf Count", &bars, &levels, &CN_COLOURS, limits, breaks, thin_bar_linewidth())
}

/// purity_ploidy_range_plot: the score of every fitted purity and ploidy as tiles, the best fit marked.
fn purity_range_png(purity_path: &str, range_path: &str, path: &str) -> Result<()> {
    let best = Table::read(purity_path)?;
    let row = best.rows.first().context("empty purity file")?;
    let best_purity = best.num(row, best.col("purity")?);
    let best_ploidy = best.num(row, best.col("ploidy")?);
    let t = Table::read(range_path)?;
    let (cp, cl, cs) = (t.col("purity")?, t.col("ploidy")?, t.col("score")?);
    let mut range: Vec<(f64, f64, f64)> = t.rows.iter().map(|r| (t.num(r, cp), t.num(r, cl), t.num(r, cs))).collect();
    let best_score = range.iter().map(|r| r.2).filter(|v| !v.is_nan()).fold(f64::INFINITY, f64::min);
    // arrange(purity, ploidy); group_by(purity); lag/lead of ploidy within the group
    range.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    struct Tile {
        purity: f64,
        score: f64,
        xmin: f64,
        xmax: f64,
    }
    let mut tiles: Vec<Tile> = Vec::with_capacity(range.len());
    let mut i = 0;
    while i < range.len() {
        let mut j = i;
        while j < range.len() && range[j].0 == range[i].0 {
            j += 1;
        }
        for k in i..j {
            let (purity, ploidy, score) = range[k];
            let left = if k > i { Some(range[k - 1].1) } else { None };
            let right = if k + 1 < j { Some(range[k + 1].1) } else { None };
            let xmin = left.map(|l| ploidy - (ploidy - l) / 2.0).unwrap_or(ploidy);
            let xmax = right.map(|r| ploidy + (r - ploidy) / 2.0).unwrap_or(ploidy);
            let rel = 1.0f64.min((score - best_score).abs() / score);
            tiles.push(Tile { purity, score: rel, xmin, xmax });
        }
        i = j;
    }
    // per purity: the xmax of the highest ploidy and the xmin of the lowest
    let mut per_purity: BTreeMap<u64, (f64, f64)> = BTreeMap::new();
    for tl in &tiles {
        let e = per_purity.entry(tl.purity.to_bits()).or_insert((f64::INFINITY, f64::NEG_INFINITY));
        e.0 = e.0.min(tl.xmin);
        e.1 = e.1.max(tl.xmax);
    }
    // the highest-ploidy row of a group has the largest xmax (sorted): min over groups; likewise the lowest
    let max_ploidy = per_purity.values().map(|v| v.1).fold(f64::INFINITY, f64::min).max(best_ploidy);
    let min_ploidy = per_purity.values().map(|v| v.0).fold(f64::NEG_INFINITY, f64::max).min(best_ploidy);
    let tiles: Vec<Tile> = tiles
        .into_iter()
        .filter(|t| t.xmin <= max_ploidy && t.xmax >= min_ploidy)
        .map(|t| Tile { xmax: t.xmax.min(max_ploidy), xmin: t.xmin.max(min_ploidy), ..t })
        .collect();
    let ramp = GradientN::new(
        &["blue", "blue", "green", "yellow", "orange", "red", "red2"].map(hex),
        Some(&[0.0, 0.0999, 0.1, 0.5, 0.8, 0.9, 1.0]),
    );
    // scales: x over the tiles, the segments and the label; y over the tiles and the segments
    let xlo = tiles.iter().map(|t| t.xmin).fold(f64::INFINITY, f64::min).min(min_ploidy).min(best_ploidy);
    let xhi = tiles.iter().map(|t| t.xmax).fold(f64::NEG_INFINITY, f64::max).max(max_ploidy + 0.4).max(best_ploidy);
    let ylo = tiles.iter().map(|t| t.purity - 0.005).fold(f64::INFINITY, f64::min).min(0.085).min(best_purity);
    let yhi = tiles.iter().map(|t| t.purity + 0.005).fold(f64::NEG_INFINITY, f64::max).max(1.05).max(best_purity);
    let xs = PosScale::new((xlo, xhi));
    let mut ys = PosScale::new((ylo, yhi));
    ys.breaks = Some(vec![0.25, 0.5, 0.75, 1.0]);
    ys.labels = Labels::Fixed(vec!["25%".into(), "50%".into(), "75%".into(), "100%".into()]);
    let legend = Legend::Bar {
        title: Some("Relative\nScore".into()),
        ramp: ramp.clone(),
        limits: (0.0, 1.0),
        breaks: vec![0.1, 0.5, 1.0],
        labels: vec!["10%".into(), "50%".into(), "100%".into()],
    };
    let plot = Plot {
        title: Some("Purity/Ploidy Scores".into()),
        x: xs.view(),
        y: ys.view(),
        x_title: Some("Ploidy".into()),
        y_title: Some("Purity".into()),
        legends: vec![legend],
        ticks: false,
        minor_grid: false,
        panel_border: true,
    };
    let mut fig = gg::figure(4.8, 4.0);
    let (w, h) = (fig.width as f64, fig.height as f64);
    plot.draw(&mut fig, 0.0, 0.0, w, h, None, &mut |f, p| {
        for t in &tiles {
            let c = if t.score.is_nan() { hex("grey50") } else { ramp.at(t.score).unwrap_or(hex("grey50")) };
            tile(f, p, t.xmin, t.xmax, t.purity - 0.005, t.purity + 0.005, c);
        }
        // geom_segment inherits the plot's data, so R draws each dashed line once per row of the range file;
        // their anti-aliased edges then add up to solid black. Drawn here up to 64 times, which looks the same.
        let copies = tiles.len().clamp(1, 64);
        let line = |f: &mut Figure, x1: f64, y1: f64, x2: f64, y2: f64| {
            for _ in 0..copies {
                gg::path(f, vec![(p.px(x1), p.py(y1)), (p.px(x2), p.py(y2))], Rgba::BLACK, 0.1, true);
            }
        };
        line(f, best_ploidy, 0.085, best_ploidy, 1.05);
        gg::label(f, p.px(best_ploidy), p.py(1.05), &rfmt::r_as_character(rfmt::r_round(best_ploidy, 2)), 2.5, 0.5, 0.5);
        line(f, min_ploidy, best_purity, max_ploidy + 0.4, best_purity);
        gg::label(f, p.px(max_ploidy + 0.4), p.py(best_purity), &format!("{}%", rfmt::r_as_character(best_purity * 100.0)), 2.5, 0.7, 0.5);
    });
    write(&fig, path)
}

/// A tile of a heat map (geom_rect without outline).
fn tile(f: &mut Figure, p: &Panel, xmin: f64, xmax: f64, ymin: f64, ymax: f64, c: Rgba) {
    gg::fill_rect_r(f, p.px(xmin), p.py(ymax), p.px(xmax), p.py(ymin), c);
}

/// fitted_segments_plot: diploid segments' major against minor allele copy number, sized by BAF count and
/// coloured by score.
fn segment_png(segment_path: &str, path: &str) -> Result<()> {
    let t = Table::read(segment_path)?;
    let (cg, cb, cbaf, ccn, cdp, cep) = (t.col("germlineStatus")?, t.col("bafCount")?, t.col("tumorBAF")?, t.col("tumorCopyNumber")?, t.col("deviationPenalty")?, t.col("eventPenalty")?);
    struct Seg {
        major: f64,
        minor: f64,
        score: f64,
        weight: f64,
    }
    let mut segs: Vec<Seg> = t
        .rows
        .iter()
        .filter(|r| r[cg] == "DIPLOID" && t.num(r, cb) > 0.0)
        .map(|r| {
            let major = t.num(r, cbaf) * t.num(r, ccn);
            Seg { major, minor: t.num(r, ccn) - major, score: t.num(r, cdp) * t.num(r, cep), weight: t.num(r, cb) }
        })
        .collect();
    segs.sort_by(|a, b| a.major.total_cmp(&b.major)); // arrange(major), stable
    let mut cum = 0.0;
    let cums: Vec<f64> = segs.iter().map(|s| {
        cum += s.weight * s.major;
        cum
    }).collect();
    let max_cum = cums.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let max_data: Vec<&Seg> = segs.iter().zip(&cums).filter(|(_, c)| **c / max_cum <= 0.9).map(|(s, _)| s).collect();
    if max_data.is_empty() {
        bail!("segment plot: no segments below the 90% cumulative proportion");
    }
    let max_score = max_data.iter().map(|s| s.score).fold(f64::NEG_INFINITY, f64::max).ceil();
    let min_score = max_data.iter().map(|s| s.score).fold(f64::INFINITY, f64::min).floor();
    let min_major = 0.0f64.min(max_data.iter().map(|s| s.major).fold(f64::INFINITY, f64::min).floor());
    let max_major = max_data.iter().map(|s| s.major).fold(f64::NEG_INFINITY, f64::max).ceil();
    let max_minor = max_major - 1.0;
    let ramp = GradientN::new(&["blue", "green", "yellow", "orange", "red"].map(hex), None);
    let (wmin, wmax) = segs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, s| (a.0.min(s.weight), a.1.max(s.weight)));
    let all_breaks: Vec<f64> = (-200..=200).map(|v| v as f64).collect();
    let mut xs = PosScale::new((min_major, max_major));
    xs.breaks = Some(all_breaks.clone());
    let mut ys = PosScale::new((0.0, max_minor));
    ys.breaks = Some(all_breaks);
    let score_breaks = rfmt::extended_breaks(min_score, max_score, 5).into_iter().filter(|b| *b >= min_score && *b <= max_score).collect::<Vec<_>>();
    let legend = Legend::Bar {
        title: Some("Score".into()),
        ramp: ramp.clone(),
        limits: (min_score, max_score),
        labels: rfmt::r_format(&score_breaks, true),
        breaks: score_breaks,
    };
    let plot = Plot {
        title: Some("Segment Scores".into()),
        x: xs.view(),
        y: ys.view(),
        x_title: Some("Major Allele".into()),
        y_title: Some("Minor Allele".into()),
        legends: vec![legend],
        ticks: false,
        minor_grid: false,
        panel_border: true,
    };
    let mut fig = gg::figure(4.8, 4.0);
    let (w, h) = (fig.width as f64, fig.height as f64);
    plot.draw(&mut fig, 0.0, 0.0, w, h, None, &mut |f, p| {
        for s in &segs {
            // the x and y limits remove what lies outside them
            if !(s.major >= min_major && s.major <= max_major && s.minor >= 0.0 && s.minor <= max_minor) {
                continue;
            }
            let x01 = if wmax > wmin { (s.weight - wmin) / (wmax - wmin) } else { 0.5 };
            let size = 1.0 + x01.sqrt() * 8.0; // scale_size(range = c(1, 9)): area palette
            let c01 = (s.score - min_score) / (max_score - min_score);
            let colour = ramp.at(c01).unwrap_or(hex("grey50"));
            gg::point(f, p.px(s.major), p.py(s.minor), size, colour.with_alpha(0.7), 0.5);
        }
    });
    write(&fig, path)
}

// -------------------------------------------------------------------------------------------------------------
// somaticVariantPlots.R

struct HistRow {
    vcn: f64,
    cn: f64,
    count: f64,
}

fn read_hist(path: &str) -> Result<Vec<HistRow>> {
    let t = Table::read(path)?;
    let (a, b, c) = (t.col("variantCopyNumberBucket")?, t.col("copyNumberBucket")?, t.col("count")?);
    Ok(t.rows.iter().map(|r| HistRow { vcn: t.num(r, a), cn: t.num(r, b), count: t.num(r, c) }).collect())
}

/// somatic_ploidy_pdf: the variant copy number histogram, filled by copy number.
fn somatic_png(hist: &[HistRow], path: &str) -> Result<()> {
    let max_ploidy = quantile_bucket(&hist.iter().map(|r| (r.vcn, r.count)).collect::<Vec<_>>(), 0.95);
    let cns: Vec<String> = hist
        .iter()
        .map(|r| {
            let cn = r_round0(r.cn).max(0.0);
            if cn >= 6.0 { "CN6+".to_string() } else { format!("CN{}", paste_int(cn)) }
        })
        .collect();
    let levels = levels_of(&cns);
    let hi = max_ploidy.map(|m| m + 1.1);
    // stat_count per fill level: weights summed at each x; bars 0.05 wide
    let mut sums: BTreeMap<(u64, usize), f64> = BTreeMap::new();
    for (r, c) in hist.iter().zip(&cns) {
        if r.vcn < -0.1 || hi.is_some_and(|h| r.vcn > h) {
            continue;
        }
        let g = levels.iter().position(|l| l == c).unwrap();
        *sums.entry((r.vcn.to_bits(), g)).or_insert(0.0) += r.count;
    }
    let rows: Vec<(f64, f64, f64, usize)> = sums.iter().map(|(&(x, g), &n)| (f64::from_bits(x) - 0.025, f64::from_bits(x) + 0.025, n, g)).collect();
    let bars = stack(rows);
    let x_hi = hi.unwrap_or_else(|| bars.iter().map(|b| b.xmax).fold(0.0, f64::max));
    let breaks: Vec<f64> = (0..=10).map(|v| v as f64).collect();
    stacked_plot(path, "Somatic Variant Copy Number", "Variant Copy Number", "Count", &bars, &levels, &CN_COLOURS, (-0.1, x_hi), breaks, 0.07)
}

struct PeakRow {
    peak: String,
    bucket: f64,
    weight: f64,
    subclonal: bool,
}

/// clonality_plot: the variant copy number histogram with the fitted peaks above, the subclonal likelihood
/// below (cowplot::plot_grid, rel_heights 5:1, aligned vertically).
fn clonality_png(hist: &[HistRow], clonality_path: &str, path: &str) -> Result<()> {
    let t = Table::read(clonality_path)?;
    let (cpk, cb, cw, cv, cs) = (t.col("peak")?, t.col("bucket")?, t.col("bucketWeight")?, t.col("isValid")?, t.col("isSubclonal")?);
    let model: Vec<PeakRow> = t
        .rows
        .iter()
        .filter(|r| r[cv] == "true")
        .map(|r| PeakRow { peak: rfmt::r_as_character(t.num(r, cpk)), bucket: t.num(r, cb), weight: t.num(r, cw), subclonal: r[cs] == "true" })
        .collect();
    // clonalityVariants: counts summed per variant copy number bucket
    let mut variants: BTreeMap<u64, f64> = BTreeMap::new();
    for r in hist {
        *variants.entry(r.vcn.to_bits()).or_insert(0.0) += r.count;
    }
    let mut variants: Vec<(f64, f64)> = variants.into_iter().map(|(k, v)| (f64::from_bits(k), v)).collect();
    variants.sort_by(|a, b| a.0.total_cmp(&b.0));
    // per bucket: total weight, and the subclonal weight over the total
    let subclonal_pct = |rows: &[&PeakRow]| -> Vec<(f64, f64, f64)> {
        let mut by: BTreeMap<u64, (f64, f64, bool)> = BTreeMap::new();
        for r in rows {
            let e = by.entry(r.bucket.to_bits()).or_insert((0.0, 0.0, false));
            e.0 += r.weight;
            if r.subclonal {
                e.1 += r.weight;
                e.2 = true;
            }
        }
        let mut out: Vec<(f64, f64, f64)> = by
            .into_iter()
            .filter(|(_, v)| v.2)
            .map(|(k, (total, sub, _))| (f64::from_bits(k), sub, if sub == 0.0 { 0.0 } else { sub / total }))
            .collect();
        out.sort_by(|a, b| a.0.total_cmp(&b.0));
        out
    };
    let all: Vec<&PeakRow> = model.iter().collect();
    let non_residual: Vec<&PeakRow> = model.iter().filter(|r| r.peak != "0").collect();
    let subclonal_all = subclonal_pct(&all);
    let subclonal_nr = subclonal_pct(&non_residual);
    let mut combined: BTreeMap<u64, f64> = BTreeMap::new();
    for r in &non_residual {
        *combined.entry(r.bucket.to_bits()).or_insert(0.0) += r.weight;
    }
    let mut combined: Vec<(f64, f64)> = combined.into_iter().map(|(k, v)| (f64::from_bits(k), v)).collect();
    combined.sort_by(|a, b| a.0.total_cmp(&b.0));
    let peaks = levels_of(&non_residual.iter().map(|r| r.peak.clone()).collect::<Vec<_>>());
    let peak_colours = hue_pal(peaks.len());

    let blue = hex("#6baed6");
    let red = hex("#d94701");
    let in_x = |x: f64| (0.0..=3.5).contains(&x);
    // top: y over the bars (from 0), the lines and the area
    let top_hi = variants.iter().filter(|v| in_x(v.0)).map(|v| v.1).chain(combined.iter().filter(|c| in_x(c.0)).map(|c| c.1)).chain(non_residual.iter().filter(|r| in_x(r.bucket)).map(|r| r.weight)).chain(subclonal_nr.iter().filter(|s| in_x(s.0)).map(|s| s.1)).fold(0.0, f64::max);
    let mut top_x = PosScale::new((0.0, 3.5));
    top_x.expand = (0.01, 0.01);
    let mut top_y = PosScale::new((0.0, top_hi));
    top_y.expand = (0.02, 0.02);
    let top = Plot {
        title: Some("Somatic Clonality".into()),
        x: top_x.view(),
        y: top_y.view(),
        x_title: Some("Variant Copy Number".into()),
        y_title: Some(String::new()),
        legends: vec![],
        ticks: false,
        minor_grid: false,
        panel_border: false,
    };
    let mut bottom_y = PosScale::new((0.0, 1.0));
    bottom_y.expand = (0.02, 0.02);
    bottom_y.breaks = Some(vec![0.0, 0.25, 0.5, 0.75, 1.0]);
    bottom_y.labels = Labels::Fixed(vec!["0%".into(), "25%".into(), "50%".into(), "75%".into(), "100%".into()]);
    let bottom = Plot { title: None, y: bottom_y.view(), x_title: Some(String::new()), ..top.clone() };
    let left = top.left_width().max(bottom.left_width());

    let mut fig = gg::figure(8.0, 6.0);
    let (w, h) = (fig.width as f64, fig.height as f64);
    let top_h = h * 5.0 / 6.0;
    let bar_lw = thin_bar_linewidth();
    top.draw(&mut fig, 0.0, 0.0, w, top_h, Some(left), &mut |f, p| {
        // bars whose edges (x -+ 0.025) leave the x limits are removed
        for &(x, n) in variants.iter().filter(|v| in_x(v.0 - 0.025) && in_x(v.0 + 0.025)) {
            gg::rect(f, p, x - 0.025, x + 0.025, 0.0, n, Some(blue.with_alpha(0.4)), Some((blue, bar_lw)));
        }
        let line = |f: &mut Figure, pts: &[(f64, f64)], c: Rgba| {
            let v: Vec<(f64, f64)> = pts.iter().filter(|q| in_x(q.0)).map(|q| (p.px(q.0), p.py(q.1))).collect();
            gg::path(f, v, c, 0.5, false);
        };
        line(f, &combined, Rgba::BLACK.with_alpha(0.8));
        for (peak, c) in peaks.iter().zip(&peak_colours) {
            let pts: Vec<(f64, f64)> = non_residual.iter().filter(|r| &r.peak == peak).map(|r| (r.bucket, r.weight)).collect();
            line(f, &pts, *c);
        }
        // geom_area: the polygon down to 0 filled at alpha 0.3, the upper line drawn
        let area: Vec<(f64, f64)> = subclonal_nr.iter().filter(|s| in_x(s.0)).map(|s| (s.0, s.1)).collect();
        if area.len() >= 2 {
            let mut poly: Vec<(f64, f64)> = area.iter().map(|q| (p.px(q.0), p.py(q.1))).collect();
            poly.push((p.px(area.last().unwrap().0), p.py(0.0)));
            poly.push((p.px(area[0].0), p.py(0.0)));
            f.polygon(poly, Some(red.with_alpha(0.3)), None);
            line(f, &area, red);
        }
    });
    bottom.draw(&mut fig, 0.0, top_h, w, h - top_h, Some(left), &mut |f, p| {
        for &(x, _, l) in subclonal_all.iter().filter(|s| in_x(s.0 - 0.025) && in_x(s.0 + 0.025)) {
            gg::rect(f, p, x - 0.025, x + 0.025, 0.0, l, Some(red.with_alpha(0.3)), Some((red, 0.5)));
        }
    });
    write(&fig, path)
}

/// standard_mutation: the pyrimidine context of a substitution.
fn standard_mutation(m: &str) -> String {
    match m {
        "G>T" => "C>A",
        "G>C" => "C>G",
        "G>A" => "C>T",
        "A>T" => "T>A",
        "A>G" => "T>C",
        "A>C" => "T>G",
        other => other,
    }
    .to_string()
}

const SBS_COLOURS: [(&str, &str); 6] = [("C>A", "#14B0EF"), ("C>G", "#060809"), ("C>T", "#E00714"), ("T>A", "#BFBEBF"), ("T>C", "#90CA4B"), ("T>G", "#E9BBB8")];

/// rainfall_plot: every SNV's distance to its nearest neighbour against its rank, kataegis clusters shaded.
fn rainfall_png(points: &[RainfallPoint], path: &str) -> Result<()> {
    let n = points.len();
    if n == 0 {
        // R draws the frame and the titles only: no breaks, no legends
        let empty = gg::View { range: (0.0, 1.0), trans: Trans::Identity, majors: vec![], labels: vec![], minors: vec![] };
        let plot = Plot {
            title: Some("Somatic Rainfall & Kataegis".into()),
            x: empty.clone(),
            y: empty,
            x_title: Some("Mutation number".into()),
            y_title: Some("Intermutation distance (bp)".into()),
            legends: vec![],
            ticks: true,
            minor_grid: true,
            panel_border: true,
        };
        let mut fig = gg::figure(8.0, 4.0);
        let (w, h) = (fig.width as f64, fig.height as f64);
        plot.draw(&mut fig, 0.0, 0.0, w, h, None, &mut |_, _| {});
        return write(&fig, path);
    }
    let pos: Vec<f64> = points.iter().map(|p| p.position as f64).collect();
    let same_chr = |a: usize, b: usize| points[a].chromosome == points[b].chromosome;
    // lag/lead of position, 0 where missing (PL3: R takes them across chromosome boundaries)
    let dist: Vec<f64> = (0..n)
        .map(|i| {
            let prev = if i > 0 && (java("PL3") || same_chr(i - 1, i)) { pos[i - 1] } else { 0.0 };
            let next = if i + 1 < n && (java("PL3") || same_chr(i, i + 1)) { pos[i + 1] } else { 0.0 };
            (pos[i] - prev).abs().min((next - pos[i]).abs())
        })
        .collect();
    let ymin = dist.iter().copied().fold(f64::INFINITY, f64::min);
    let ymax = dist.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    // kataegis clusters: rank range of each KT id, Forward for FWD ids
    let mut clusters: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for (i, p) in points.iter().enumerate() {
        if !p.kataegis.is_empty() {
            let e = clusters.entry(p.kataegis.clone()).or_insert((i + 1, i + 1));
            e.0 = e.0.min(i + 1);
            e.1 = e.1.max(i + 1);
        }
    }
    let strand_colours = [("Forward", hex("#6bd692")), ("Reverse", hex("#7e6bd6"))];
    let rects: Vec<(f64, f64, usize)> = clusters
        .iter()
        .map(|(k, &(a, b))| (a as f64, b as f64, if k.get(..3) == Some("FWD") { 0 } else { 1 }))
        .collect();
    let mutations: Vec<String> = points.iter().map(|p| standard_mutation(&p.mutation)).collect();
    let mut_levels = levels_of(&mutations);
    let colour_of = |m: &str| SBS_COLOURS.iter().find(|(k, _)| *k == m).map(|(_, c)| hex(c)).unwrap_or(hex("grey50"));

    let positive = |v: f64| v > 0.0;
    let ylo = dist.iter().copied().filter(|v| positive(*v)).fold(f64::INFINITY, f64::min).min(if positive(ymin) { ymin } else { f64::INFINITY });
    let yhi = ymax;
    let xhi = (n as f64).max(rects.iter().map(|r| r.1).fold(0.0, f64::max));
    let xs = PosScale::new((1.0f64.min(xhi), xhi));
    let mut ys = PosScale::new((ylo, yhi));
    ys.trans = Trans::Log10;
    ys.labels = Labels::FormatNoSci;
    let mut legends = Vec::new();
    let strands: Vec<usize> = [0usize, 1].into_iter().filter(|s| rects.iter().any(|r| r.2 == *s)).collect();
    if !strands.is_empty() {
        legends.push(Legend::Keys {
            title: Some("Kataegis Regions".into()),
            entries: strands.iter().map(|&s| (strand_colours[s].0.to_string(), Key::Polygon { fill: strand_colours[s].1.with_alpha(0.6), colour: None, linewidth: 0.0 })).collect(),
        });
    }
    if !mut_levels.is_empty() {
        legends.push(Legend::Keys {
            title: Some("Mutation".into()),
            entries: mut_levels.iter().map(|m| (m.clone(), Key::Point { colour: colour_of(m), size: 2.0 })).collect(),
        });
    }
    let plot = Plot {
        title: Some("Somatic Rainfall & Kataegis".into()),
        x: xs.view(),
        y: ys.view(),
        x_title: Some("Mutation number".into()),
        y_title: Some("Intermutation distance (bp)".into()),
        legends,
        ticks: true,
        minor_grid: true,
        panel_border: true,
    };
    let mut fig = gg::figure(8.0, 4.0);
    let (w, h) = (fig.width as f64, fig.height as f64);
    plot.draw(&mut fig, 0.0, 0.0, w, h, None, &mut |f, p| {
        for &(a, b, s) in &rects {
            // ymin may be 0: -Inf on the log scale, drawn to the panel's edge
            let y0 = if positive(ymin) { p.py(ymin) } else { p.y1 };
            gg::fill_rect_r(f, p.px(a), p.py(ymax), p.px(b), y0, strand_colours[s].1.with_alpha(0.6));
        }
        for (i, d) in dist.iter().enumerate() {
            if !positive(*d) {
                continue;
            }
            gg::point(f, p.px((i + 1) as f64), p.py(*d), 0.1, colour_of(&mutations[i]), 0.5);
        }
    });
    write(&fig, path)
}

// -------------------------------------------------------------------------------------------------------------

/// What the plots need beyond the files in the output folder.
pub struct PlotInputs<'a> {
    pub output_dir: &'a str,
    pub plot_dir: &'a str,
    pub tumor_id: &'a str,
    pub plot_somatics: bool,
    pub rainfall: &'a [RainfallPoint],
}

/// RCharts: copyNumberPlots.R (copynumber and map when there are autosomal segments with BAF support, purity
/// range, segment) and, with a somatic VCF, somaticVariantPlots.R (somatic, clonality, rainfall).
pub fn draw_r_plots(inp: &PlotInputs) -> Result<()> {
    let (out, plot, tumor) = (inp.output_dir, inp.plot_dir, inp.tumor_id);
    let cnv = read_cnv(&format!("{out}{tumor}.purple.cnv.somatic.tsv"))?;
    if !cnv.is_empty() {
        copynumber_png(&cnv, &format!("{plot}{tumor}.copynumber.png"))?;
        map_png(&cnv, &format!("{plot}{tumor}.map.png"))?;
    }
    purity_range_png(&format!("{out}{tumor}.purple.purity.tsv"), &format!("{out}{tumor}.purple.purity.range.tsv"), &format!("{plot}{tumor}.purity.range.png"))?;
    segment_png(&format!("{out}{tumor}.purple.segment.tsv"), &format!("{plot}{tumor}.segment.png"))?;
    if inp.plot_somatics {
        let hist = read_hist(&format!("{out}{tumor}.purple.somatic.hist.tsv"))?;
        somatic_png(&hist, &format!("{plot}{tumor}.somatic.png"))?;
        clonality_png(&hist, &format!("{out}{tumor}.purple.somatic.clonality.tsv"), &format!("{plot}{tumor}.somatic.clonality.png"))?;
        rainfall_png(inp.rainfall, &format!("{plot}{tumor}.somatic.rainfall.png"))?;
    }
    Ok(())
}

/// plot/<tumor>.somatic_data.tsv as Java writes it for R.
pub fn write_somatic_data(path: &str, points: &[RainfallPoint]) -> Result<()> {
    let mut text = String::from("chromosome\tposition\tmutation\tkataegis\n");
    for p in points {
        text.push_str(&format!("{}\t{}\t{}\t{}\n", p.chromosome, p.position, p.mutation, p.kataegis));
    }
    std::fs::write(path, text).with_context(|| format!("writing {path}"))
}

/// The Circos images, from the configuration and data files in circos/.
pub fn draw_circos(circos_dir: &str, plot_dir: &str, tumor: &str) -> Result<()> {
    let opts = hmf_plot::circos::Options { ignore_skip_run: !java("PL1") };
    for kind in ["input", "circos"] {
        let conf = format!("{circos_dir}{tumor}.{kind}.conf");
        let fig = hmf_plot::circos::draw_conf_file(std::path::Path::new(&conf), opts.clone()).map_err(|e| anyhow::anyhow!("circos {kind}: {e}"))?;
        write(&fig, &format!("{plot_dir}{tumor}.{kind}.png"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stacking_puts_the_first_level_on_top() {
        let bars = stack(vec![(0.0, 1.0, 2.0, 0), (0.0, 1.0, 3.0, 1), (1.0, 2.0, 1.0, 0)]);
        // at x 0: level 1 at the bottom (0..3), level 0 above (3..5); drawn bottom first
        assert_eq!((bars[0].group, bars[0].ymin, bars[0].ymax), (1, 0.0, 3.0));
        assert_eq!((bars[1].group, bars[1].ymin, bars[1].ymax), (0, 3.0, 5.0));
        assert_eq!(bars[2].xmin, 1.0);
    }

    #[test]
    fn quantile_bucket_as_dplyr() {
        // buckets 1: 10, 2: 80, 3: 10 -> proportions 0.1, 0.9, 1.0: the first above 0.9 is bucket 3
        let rows = [(0.5, 10.0), (1.2, 50.0), (1.9, 30.0), (2.5, 10.0)];
        assert_eq!(quantile_bucket(&rows, 0.9), Some(3.0));
        assert_eq!(quantile_bucket(&rows, 0.05), Some(1.0));
    }

    #[test]
    fn histogram_bins_weights_and_censors() {
        // x within limits (-0.1, 1.1), binwidth 0.1; 2.0 is censored
        let bars = binned(&[(0.0, 2.0, 0), (0.04, 1.0, 0), (0.06, 5.0, 0), (2.0, 7.0, 0)], (-0.1, 1.1), 0.1, 1);
        let total: f64 = bars.iter().map(|b| b.ymax - b.ymin).sum();
        assert_eq!(total, 8.0);
        let first = bars.iter().find(|b| b.ymax == 3.0).unwrap();
        assert!((first.xmin + 0.05).abs() < 1e-9 && (first.xmax - 0.05).abs() < 1e-9);
    }

    #[test]
    fn empty_rainfall_is_a_frame() {
        let dir = std::env::temp_dir().join(format!("purple_rs_plots_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("rain.png");
        rainfall_png(&[], path.to_str().unwrap()).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        // PNG signature and IHDR size 2400 x 1200
        assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n");
        assert_eq!(u32::from_be_bytes(bytes[16..20].try_into().unwrap()), 2400);
        assert_eq!(u32::from_be_bytes(bytes[20..24].try_into().unwrap()), 1200);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn mutations_standardised() {
        assert_eq!(standard_mutation("G>A"), "C>T");
        assert_eq!(standard_mutation("T>G"), "T>G");
    }
}
