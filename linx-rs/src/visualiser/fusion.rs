//! fusionPlot.R: under the Circos image, one panel per fusion (250 pixels high) and, when the fusions have
//! protein domains other than UTRs, a legend row (35 pixels), all 3000 pixels wide at 300 dpi:
//! - plot_fusion: the fused genes as thin bars, their exons as boxes (upstream blue, downstream red), the
//!   gene names above, UTRs as narrow bars on white, the other protein domains as translucent bars, the
//!   exons not in the fusion faded, a black line at the upstream gene's end, the exon ranks under their
//!   starts; ggplot2's theme_bw with the script's changes (no grid, border, ticks or y axis; 3-point side
//!   margins; an empty x title).
//! - plot_legend: the domains' colours with their names, as cowplot::get_legend extracts the legend.
//!
//! Finding LV3: R (ggplot2 4.0.3, cowplot 1.2.0) draws the legend's row black, so that only the colour keys
//! show and the domain names (black text) are hidden. Java mode draws the black row; by default the legend
//! is drawn on white.

use anyhow::Result;
use hmf_plot::color::Rgba;
use hmf_plot::compose::Canvas;
use hmf_plot::display::{Figure, LineCap, Stroke};
use hmf_plot::gg::{self, Panel, Trans};
use hmf_plot::rpdf;
use hmf_plot::text::Font;

use super::circos_data::{FusedExon, FusionData};
use super::files::{VisProteinDomain, PROTEIN_DOMAIN_UTR};
use super::VisOptions;
use crate::compat::java;

const SINGLE_BLUE: &str = "#6baed6";
const SINGLE_RED: &str = "#d6906b";

fn hex(s: &str) -> Rgba {
    Rgba::from_hex(s).unwrap_or(Rgba::BLACK)
}

/// A rectangle of a layer, in data units.
#[derive(Clone, Debug)]
struct RectL {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    fill: Rgba,
}

/// The layers of plot_fusion for one fusion, in drawing order (the script's dplyr steps).
struct FusionLayers {
    genes: Vec<RectL>,
    exons: Vec<RectL>,
    /// (x, label, hjust)
    labels: Vec<(f64, String, f64)>,
    utr_white: Vec<RectL>,
    utr: Vec<RectL>,
    other: Vec<RectL>,
    faded: Vec<RectL>,
    /// x of the black line at the upstream gene's end
    segments: Vec<f64>,
    /// (x, label) of the axis breaks
    breaks: Vec<(f64, String)>,
}

/// The domain names as fusionPlot.R reads them (`gsub(" domain", "", name)`).
fn domain_name(d: &VisProteinDomain) -> String {
    d.name().replace(" domain", "")
}

fn layers(exons: &[&FusedExon], domains: &[(&VisProteinDomain, String, Rgba)]) -> FusionLayers {
    let colour = |e: &FusedExon| hex(if e.is_upstream { SINGLE_BLUE } else { SINGLE_RED });
    // fusedGenes: group_by(gene, transcript, upGene, color), start = max(ifelse(rank == 1, start, geneStart),
    // geneStart), end = max(geneEnd); groups in sorted order
    let mut keys: Vec<(String, String, bool)> = Vec::new();
    for e in exons {
        let k = (e.gene.clone(), e.transcript.clone(), e.is_upstream);
        if !keys.contains(&k) {
            keys.push(k);
        }
    }
    keys.sort();
    let mut genes = Vec::new();
    let mut labels: Vec<(f64, String, f64)> = Vec::new();
    let mut segments = Vec::new();
    for (gene, transcript, up) in &keys {
        let g: Vec<&&FusedExon> = exons.iter().filter(|e| e.gene == *gene && e.transcript == *transcript && e.is_upstream == *up).collect();
        let mut start = f64::MIN;
        let mut end = f64::MIN;
        for e in &g {
            let s = if e.rank == 1 { e.start } else { e.gene_start } as f64;
            start = start.max(s).max(e.gene_start as f64);
            end = end.max(e.gene_end as f64);
        }
        let fill = hex(if *up { SINGLE_BLUE } else { SINGLE_RED });
        genes.push(RectL { xmin: start, xmax: end, ymin: 0.45, ymax: 0.55, fill });
        let label = format!("{gene} - {transcript}");
        if *up {
            labels.push((start, label, 0.0));
            segments.push(end);
        } else {
            labels.push((end, label, 1.0));
        }
    }
    // the upstream labels first (their layer comes before the downstream one's)
    labels.sort_by(|a, b| a.2.total_cmp(&b.2));
    let exon_rects: Vec<RectL> = exons.iter().map(|e| RectL { xmin: e.start as f64, xmax: e.end as f64, ymin: 0.0, ymax: 1.0, fill: colour(e) }).collect();
    let utr: Vec<&(&VisProteinDomain, String, Rgba)> = domains.iter().filter(|d| d.1 == PROTEIN_DOMAIN_UTR).collect();
    let mut other: Vec<&(&VisProteinDomain, String, Rgba)> = domains.iter().filter(|d| d.1 != PROTEIN_DOMAIN_UTR).collect();
    other.sort_by(|a, b| a.1.cmp(&b.1));
    // scale_fill_manual: a name's colour is the first in the data
    let fill_of = |name: &str| domains.iter().find(|d| d.1 == name).map(|d| d.2).unwrap_or(Rgba::BLACK);
    // fadedArea: the end of the upstream gene's included exons to the start of the downstream gene's
    let up_end = exons.iter().filter(|e| !e.skipped && e.is_upstream).map(|e| e.end).max();
    let down_start = exons.iter().filter(|e| !e.skipped && !e.is_upstream).map(|e| e.start).min();
    let faded = match (up_end, down_start) {
        (Some(a), Some(b)) => vec![RectL { xmin: a as f64, xmax: b as f64, ymin: 0.0, ymax: 1.0, fill: Rgba::WHITE.with_alpha(0.6) }],
        _ => Vec::new(),
    };
    FusionLayers {
        genes,
        exons: exon_rects,
        labels,
        utr_white: utr.iter().map(|d| RectL { xmin: d.0.start as f64, xmax: d.0.end as f64, ymin: 0.0, ymax: 1.0, fill: Rgba::WHITE }).collect(),
        utr: utr.iter().map(|d| RectL { xmin: d.0.start as f64, xmax: d.0.end as f64, ymin: 0.3, ymax: 0.7, fill: d.2 }).collect(),
        other: other.iter().map(|d| RectL { xmin: d.0.start as f64, xmax: d.0.end as f64, ymin: 0.3, ymax: 0.7, fill: fill_of(&d.1).with_alpha(0.8) }).collect(),
        faded,
        segments,
        breaks: exons.iter().map(|e| (e.start as f64, e.rank.to_string())).collect(),
    }
}

/// grid points to pixels at 300 dpi
fn pt(v: f64) -> f64 {
    gg::pt(v)
}

/// One fusion panel into `fig`, the cell at (0, top), `width` x `height` pixels.
fn draw_fusion(fig: &mut Figure, l: &FusionLayers, top: f64, width: f64, height: f64, font_size: f64) {
    let axis_size = font_size * 72.0 / 300.0; // points
    let text_size_mm = font_size * 25.4 / 300.0;
    let (asc, _) = gg::text_metrics(axis_size);
    // cowplot turns the plots into grobs on a pdf(NULL) device, so ggplot2 sizes the axis labels and the
    // empty x title with the PDF device's metrics (Helvetica at whole point sizes)
    let axis_h = pt(2.2) + gg::font_px(rpdf::text_height_bigpt(axis_size));
    let xlab_h = pt(2.75) + gg::font_px(rpdf::text_height_bigpt(5.0));
    // the x scale's limits from every layer, expanded by 5%
    let mut xs: Vec<f64> = Vec::new();
    for r in l.genes.iter().chain(&l.exons).chain(&l.utr_white).chain(&l.utr).chain(&l.other).chain(&l.faded) {
        xs.push(r.xmin);
        xs.push(r.xmax);
    }
    xs.extend(l.labels.iter().map(|t| t.0));
    xs.extend(l.segments.iter().copied());
    let (lo, hi) = (xs.iter().copied().fold(f64::MAX, f64::min), xs.iter().copied().fold(f64::MIN, f64::max));
    let d = hi - lo;
    let p = Panel {
        x0: pt(3.0),
        y0: top,
        x1: width - pt(3.0),
        y1: top + height - axis_h - xlab_h,
        xr: (lo - 0.05 * d, hi + 0.05 * d),
        yr: (-0.1, 2.1),
        xt: Trans::Identity,
        yt: Trans::Identity,
    };
    fig.clip(Some((p.x0, p.y0, p.x1 - p.x0, p.y1 - p.y0)));
    let rects = |fig: &mut Figure, rs: &[RectL]| {
        for r in rs {
            gg::rect(fig, &p, r.xmin, r.xmax, r.ymin, r.ymax, Some(r.fill), None);
        }
    };
    rects(fig, &l.genes);
    rects(fig, &l.exons);
    // geom_text, vjust 0: the baseline at y = 1.15
    let px = gg::font_px(text_size_mm * gg::PT);
    for (x, label, hjust) in &l.labels {
        let w = hmf_plot::text::advance(Font::Sans, px, label);
        fig.text((p.px(*x) - hjust * w).round(), p.py(1.15).round(), label, Font::Sans, px, Rgba::BLACK, 0.0);
    }
    rects(fig, &l.utr_white);
    rects(fig, &l.utr);
    rects(fig, &l.other);
    rects(fig, &l.faded);
    for x in &l.segments {
        let s = Stroke { color: Rgba::BLACK, width: gg::lw_px(0.5), dash: Vec::new(), cap: LineCap::Butt };
        fig.line(p.px(*x), p.py(-0.1), p.px(*x), p.py(1.1), s);
    }
    fig.clip(None);
    // the x axis labels (grey30), their tops 2.2 points under the panel
    let apx = gg::font_px(axis_size);
    let grey30 = hmf_plot::color::r_colour("grey30").unwrap_or(Rgba::BLACK);
    for (x, label) in &l.breaks {
        if *x < p.xr.0 || *x > p.xr.1 {
            continue;
        }
        let w = hmf_plot::text::advance(Font::Sans, apx, label);
        fig.text((p.px(*x) - w / 2.0).round(), (p.y1 + pt(2.2) + asc).round(), label, Font::Sans, apx, grey30, 0.0);
    }
}

/// The legend of plot_legend as cowplot::get_legend extracts it, centred in the row (x, top, width, height):
/// one row of keys (6 points, inset by the 0.5 mm line width) and names, 5.5 points apart.
fn draw_legend(fig: &mut Figure, entries: &[(String, Rgba)], top: f64, width: f64, height: f64, font_size: f64) {
    let size = font_size * 72.0 / 300.0;
    let px = gg::font_px(size);
    let to_pt = |v: f64| v * 72.27 / 300.0;
    let _ = to_pt;
    // label widths and the row's height as ggplot2 measures them on cowplot's pdf(NULL) device (points)
    let label_w: Vec<f64> = entries.iter().map(|e| 5.5 + rpdf::string_width_bigpt(&e.0, size) * 72.27 / 72.0).collect();
    let total_pt = 5.5 + 5.5 + label_w.iter().map(|w| 6.0 + w).sum::<f64>() + 5.5 * entries.len().saturating_sub(1) as f64 + 5.5;
    let row_pt = rpdf::text_height_bigpt(size) * 72.27 / 72.0;
    let legend_h = pt(5.5 + row_pt + 5.5);
    let left = (width - pt(total_pt)) / 2.0;
    let row_top = top + height / 2.0 - legend_h / 2.0 + pt(5.5);
    let row_h = pt(row_pt);
    let (asc, _) = gg::text_metrics(size);
    // draw_key_polygon: the key inset by half the 0.5 mm line width on each side
    let half = 0.25 * 300.0 / 25.4;
    let mut x = left + pt(5.5 + 5.5);
    for (i, (name, colour)) in entries.iter().enumerate() {
        let kw = pt(6.0);
        fig.rect(x + half, row_top + half, kw - 2.0 * half, row_h - 2.0 * half, *colour);
        x += kw;
        fig.text((x + pt(5.5)).round(), (row_top + row_h / 2.0 + asc / 2.0).round(), name, Font::Sans, px, Rgba::BLACK, 0.0);
        x += pt(label_w[i]);
        if i + 1 < entries.len() {
            x += pt(5.5);
        }
    }
}

/// The Circos image with the fusion panels (and legend) under it, as fusionPlot.R's plot_grid stacks them.
pub fn compose(circos: &Canvas, f: &FusionData, font_size: f64, opts: &VisOptions) -> Result<Canvas> {
    let mut fusions: Vec<String> = Vec::new();
    for e in &f.exons {
        if !fusions.contains(&e.fusion) {
            fusions.push(e.fusion.clone());
        }
    }
    let domains: Vec<(&VisProteinDomain, String, Rgba)> =
        f.domains.iter().map(|d| { let c = f.colours.colour(d); (d, domain_name(d), Rgba::rgb8(c.0, c.1, c.2)) }).collect();
    let non_utr: Vec<&(&VisProteinDomain, String, Rgba)> = domains.iter().filter(|d| d.1 != PROTEIN_DOMAIN_UTR).collect();
    let mut legend: Vec<(String, Rgba)> = Vec::new();
    for d in &non_utr {
        if !legend.iter().any(|l| l.0 == d.1) {
            legend.push((d.1.clone(), d.2));
        }
    }
    // the fill scale's levels: R's sort() of the names
    legend.sort_by(|a, b| r_collate(&a.0, &b.0));
    let width = circos.width as f64;
    let fusion_h = opts.fusion_height as f64 * fusions.len() as f64;
    let mut legend_h = if legend.is_empty() { 0.0 } else { opts.fusion_legend_height_per_row as f64 * opts.fusion_legend_rows as f64 };
    if !legend.is_empty() && !java("LV3") {
        // the legend's own height (its keys and names overflow the 35-pixel row in R)
        let size = font_size * 72.0 / 300.0;
        let full = gg::pt(5.5 + rpdf::text_height_bigpt(size) * 72.27 / 72.0 + 5.5).ceil();
        legend_h = legend_h.max(full);
    }
    let mut fig = Figure::new(circos.width, (fusion_h + legend_h).round() as u32, Rgba::WHITE);
    for (i, name) in fusions.iter().enumerate() {
        let exons: Vec<&FusedExon> = f.exons.iter().filter(|e| e.fusion == *name).collect();
        let ds: Vec<(&VisProteinDomain, String, Rgba)> = domains.iter().filter(|d| d.0.chromosome == *name).map(|d| (d.0, d.1.clone(), d.2)).collect();
        let l = layers(&exons, &ds);
        let h = opts.fusion_height as f64;
        draw_fusion(&mut fig, &l, i as f64 * h, width, h, font_size);
    }
    if !legend.is_empty() {
        if java("LV3") {
            fig.rect(0.0, fusion_h, width, legend_h, Rgba::BLACK);
        }
        draw_legend(&mut fig, &legend, fusion_h, width, legend_h, font_size);
    }
    let lower = Canvas::from_figure(&fig);
    let mut out = Canvas::new(circos.width, circos.height + lower.height, [255, 255, 255]);
    out.blit(circos, 0, 0);
    out.blit(&lower, 0, circos.height as i64);
    Ok(out)
}

/// R's sort() of strings with ICU's root collation: spaces and punctuation before digits before letters,
/// letters compared without case first (then lower case before upper case).
pub fn r_collate(a: &str, b: &str) -> std::cmp::Ordering {
    const PUNCT: &str = " _-,;:!?.'\"()[]{}@*/\\&#%`^+<=>|~$";
    let key = |c: char| -> (u32, u32) {
        if let Some(i) = PUNCT.find(c) {
            (0, i as u32)
        } else if c.is_ascii_digit() {
            (1, c as u32)
        } else if c.is_alphabetic() {
            (2, c.to_ascii_lowercase() as u32)
        } else {
            (3, c as u32)
        }
    };
    let primary = a.chars().map(key).cmp(b.chars().map(key));
    if primary != std::cmp::Ordering::Equal {
        return primary;
    }
    // tertiary: lower case first
    let case = |c: char| if c.is_uppercase() { 1 } else { 0 };
    a.chars().map(case).cmp(b.chars().map(case)).then(a.cmp(b))
}
