//! The visualiser's images, drawn natively with hmf-plot-rs as Java draws them with Circos and R:
//! - the Circos image (3000 x 3000) from the configuration and data files, by hmf-plot-rs's Circos engine;
//! - with fusions, `fusionPlot.R`'s panels below it (see `fusion.rs`);
//! - `chromosomeRangePlot.R`'s strip of Gviz ideograms of the plotted chromosomes below that, 0.85 times as
//!   wide as the image and 150 pixels high per row, centred.
//!
//! The composite is written at 300 dpi, as ggsave writes it.

use anyhow::{anyhow, Context, Result};
use hmf_plot::color::Rgba;
use hmf_plot::compose::Canvas;
use hmf_plot::display::Figure;
use hmf_plot::gviz;
use hmf_plot::png_out::PngMeta;

use super::circos_data::{CircosFiles, FusionData};
use super::VisOptions;
use crate::compat::java;

const CYTOBANDS_37: &str = include_str!("../../../hmf-common-rs/resources/refgenome/cytoBands.37.tsv");
const CYTOBANDS_38: &str = include_str!("../../../hmf-common-rs/resources/refgenome/cytoBands.38.tsv");

/// CytoBands.resourceAsString: the cytoband resource, written for the chromosome range plot
pub fn cytoband_resource(v37: bool) -> &'static str {
    if v37 { CYTOBANDS_37 } else { CYTOBANDS_38 }
}

pub struct PlotInputs<'a> {
    pub file_id: &'a str,
    pub conf: &'a str,
    pub files: &'a CircosFiles,
    pub fusion: Option<&'a FusionData>,
    pub chromosome_ranges: Option<&'a Vec<String>>,
    /// 1.2 times the Circos label size, as Java passes it to the R scripts
    pub label_size: f64,
    pub opts: &'a VisOptions,
    pub v37: bool,
    pub out: &'a str,
    /// -perf_debug: log the time of each step
    pub perf: bool,
    /// the last (or only) frame of the plot
    pub final_frame: bool,
}

/// The configuration as drawn. A link's frame is one more than the lowest frame of the segments it joins, so in
/// the last frame Java's rule (`var(frame) > conf(current_frame)` -> faint) greys out the links of the last
/// segment, and every link of a plot without segment tracks (LV4). Corrected: the last frame shows every
/// element in its colour (the configuration file written stays Java's).
fn drawn_conf(conf: &str, final_frame: bool) -> std::borrow::Cow<'_, str> {
    if !final_frame || java("LV4") {
        return std::borrow::Cow::Borrowed(conf);
    }
    let lines: Vec<String> = conf
        .lines()
        .map(|l| if l.trim_start().starts_with("current_frame") && l.contains('=') { "current_frame       = 2147483647".to_string() } else { l.to_string() })
        .collect();
    std::borrow::Cow::Owned(lines.join("\n") + "\n")
}

/// The Circos image from the configuration, its data files read from memory.
pub fn circos_figure(file_id: &str, conf: &str, files: &CircosFiles) -> Result<Figure> {
    let prefix = format!("{file_id}.");
    let load = |name: &str| -> std::result::Result<String, String> {
        let suffix = name.strip_prefix(&prefix).ok_or_else(|| format!("unknown data file {name}"))?;
        let lines = files.get(suffix).ok_or_else(|| format!("unknown data file {name}"))?;
        let mut s = String::new();
        for l in lines {
            s.push_str(l);
            s.push('\n');
        }
        Ok(s)
    };
    let c = hmf_plot::circos::Circos::new(conf, &load, hmf_plot::circos::Options::default()).map_err(|e| anyhow!("circos: {e}"))?.with_gd_fill_extent(true).with_gd_text_bounds(true);
    c.draw().map_err(|e| anyhow!("circos: {e}"))
}

pub fn draw(p: &PlotInputs) -> Result<()> {
    let t0 = std::time::Instant::now();
    let conf = drawn_conf(p.conf, p.final_frame);
    let fig = circos_figure(p.file_id, &conf, p.files)?;
    let t1 = std::time::Instant::now();
    let circos = Canvas::from_figure(&fig);
    let t2 = std::time::Instant::now();
    // fusionPlot.R: the fusion panels (and their legend) below the Circos image
    let upper = match p.fusion {
        Some(f) => super::fusion::compose(&circos, f, p.label_size, p.opts)?,
        None => circos,
    };
    let image = match p.chromosome_ranges {
        Some(ranges) => {
            let strip = chromosome_strip(ranges, p.label_size, p.opts, p.v37)?;
            stack_strip(&upper, &strip)
        }
        None => upper,
    };
    let t3 = std::time::Instant::now();
    image.write_png(p.out, PngMeta { dpi: Some(300.0), transparent_key: None }).with_context(|| format!("writing {}", p.out))?;
    if p.perf {
        hmf_common::info!(
            "perf: plot({}) circos display list {:.3} s ({} items), raster {:.3} s, fusion/strip/compose {:.3} s, png {:.3} s",
            p.file_id,
            (t1 - t0).as_secs_f64(),
            fig.items.len(),
            (t2 - t1).as_secs_f64(),
            (t3 - t2).as_secs_f64(),
            t3.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

/// chromosomeRangePlot.R's last steps: plot_grid(pCircos, pChr, rel_heights = c(circosWidth, chromosomeHeight)).
/// The upper image goes into a cell as tall as it is wide (Java: LV2) or as tall as it is; the strip, centred
/// between 3-point margins, below it.
fn stack_strip(upper: &Canvas, strip: &Canvas) -> Canvas {
    let w = upper.width;
    let cell_h = if java("LV2") { w } else { upper.height };
    let mut out = Canvas::new(w, cell_h + strip.height, [255, 255, 255]);
    if upper.height == cell_h {
        out.blit(upper, 0, 0);
    } else {
        // draw_image keeps the aspect ratio: the image fills the cell's height, centred
        let scale = cell_h as f64 / upper.height as f64;
        let dw = upper.width as f64 * scale;
        out.blit_scaled(upper, (w as f64 - dw) / 2.0, 0.0, dw, cell_h as f64);
    }
    // the strip's cell: 3-point margins left and right at 300 dpi; the strip at its size, centred
    let margin = 3.0 * 300.0 / 72.27;
    let panel_w = w as f64 - 2.0 * margin;
    let x = margin + (panel_w - strip.width as f64) / 2.0;
    if (x - x.round()).abs() < 1e-6 {
        out.blit(strip, x.round() as i64, cell_h as i64);
    } else {
        out.blit_scaled(strip, x, cell_h as f64, strip.width as f64, strip.height as f64);
    }
    out
}

struct ChrRange {
    chromosome: String,
    start: f64,
    end: f64,
    colour: Rgba,
}

/// The order of chromosomeRangePlot.R's factor levels: 1-22, X, Y
fn chr_level(c: &str) -> usize {
    match c {
        "X" => 23,
        "Y" => 24,
        n => n.parse().unwrap_or(99),
    }
}

/// chromosomeRangePlot.R's add_row_numbers: each chromosome's length relative to the first row's (at least
/// minRelativeLength), the row from the cumulative sum
fn add_row_numbers(columns: usize, lengths: &[f64], min_rel: f64) -> Vec<(f64, usize)> {
    let first_row: f64 = lengths.iter().take(columns).sum();
    let mut total = 0.0;
    lengths
        .iter()
        .map(|l| {
            let rel = f64::max(min_rel, l / first_row);
            total += rel;
            (rel, total.ceil() as usize)
        })
        .collect()
}

fn optimal_columns(max_columns: usize, max_rows: usize, min_rel: f64, lengths: &[f64]) -> usize {
    let mut rows = 0;
    for columns in 1..=max_columns {
        rows = add_row_numbers(columns, lengths, min_rel).iter().map(|r| r.1).max().unwrap_or(0);
        if rows <= max_rows {
            return columns;
        }
    }
    optimal_columns(max_columns, rows, min_rel, lengths)
}

/// The strip of chromosomeRangePlot.R: Gviz ideograms of the plotted chromosomes, laid out in rows.
fn chromosome_strip(ranges: &[String], label_size: f64, opts: &VisOptions, v37: bool) -> Result<Canvas> {
    let mut chr_ranges: Vec<ChrRange> = Vec::new();
    for l in ranges {
        let v: Vec<&str> = l.split('\t').collect();
        if v.len() < 4 {
            continue;
        }
        chr_ranges.push(ChrRange {
            chromosome: v[0].replace("chr", ""),
            start: v[1].parse().unwrap_or(0.0),
            end: v[2].parse().unwrap_or(0.0),
            colour: Rgba::from_hex(v[3]).unwrap_or(Rgba::BLACK),
        });
    }
    // the bands of every chromosome, and the lengths of those plotted, in level order
    let mut bands: Vec<(String, gviz::Band)> = Vec::new();
    for line in cytoband_resource(v37).lines().skip(1) {
        let v: Vec<&str> = line.split('\t').collect();
        if v.len() < 5 {
            continue;
        }
        bands.push((v[0].to_string(), gviz::Band { start: v[1].parse().unwrap_or(0), end: v[2].parse().unwrap_or(0), stain: v[4].to_string() }));
    }
    let mut chrs: Vec<(String, f64)> = Vec::new();
    for (chrom, b) in &bands {
        let c = chrom.replace("chr", "");
        if !chr_ranges.iter().any(|r| r.chromosome == c) {
            continue;
        }
        match chrs.iter_mut().find(|x| x.0 == c) {
            Some(x) => x.1 = x.1.max(b.end as f64),
            None => chrs.push((c, b.end as f64)),
        }
    }
    chrs.sort_by_key(|c| chr_level(&c.0));
    let lengths: Vec<f64> = chrs.iter().map(|c| c.1).collect();
    let max_columns = opts.chromosome_range_columns.max(1) as usize;
    let max_rows = chr_ranges.len().div_ceil(max_columns);
    let columns = optimal_columns(max_columns, max_rows, 0.1, &lengths);
    let rows = add_row_numbers(columns, &lengths, 0.1);
    let n_rows = rows.iter().map(|r| r.1).max().unwrap_or(1);
    // per row: columns, widths and x
    let mut row_len = vec![0.0f64; n_rows + 1];
    for (rel, row) in &rows {
        row_len[*row] += rel;
    }
    let max_row_len = row_len.iter().copied().fold(0.0, f64::max);
    let mut cells: Vec<(usize, f64, f64, usize)> = Vec::new(); // chromosome index, x, width, row
    let mut row_state: Vec<(usize, f64, f64, f64)> = vec![(0, 0.0, 0.0, 0.0); n_rows + 1]; // column count, cumsum prior width, cumsum spacing, last width
    for (i, (rel, row)) in rows.iter().enumerate() {
        let st = &mut row_state[*row];
        st.0 += 1;
        let spacing = if st.0 == 1 { (1.0 - row_len[*row]) / 2.0 } else { 0.0 };
        let width = rel / max_row_len;
        let prior = if st.0 == 1 { 0.0 } else { st.3 };
        st.1 += prior;
        st.2 += spacing;
        st.3 = width;
        cells.push((i, st.1 + st.2, width, *row));
    }
    let width_px = (0.85 * 3000.0f64).floor();
    let height_px = (opts.chromosome_range_height as usize * n_rows) as f64;
    let mut fig = Figure::new(width_px as u32, height_px as u32, Rgba::WHITE);
    for (i, x, w, row) in cells {
        let (c, _) = &chrs[i];
        let Some(r) = chr_ranges.iter().find(|r| r.chromosome == *c) else { continue };
        let chr_bands: Vec<gviz::Band> = bands.iter().filter(|(b, _)| b.replace("chr", "") == *c && *b == format!("chr{c}")).map(|(_, b)| b.clone()).collect();
        let h = 1.0 / n_rows as f64;
        let y = (n_rows - row) as f64 * h;
        // the viewport: (x, y) its bottom-left corner, in the device's y-up coordinates
        let cell = (x * width_px, height_px - (y + h) * height_px, w * width_px, h * height_px);
        let label = format!("CHR {c}");
        gviz::draw(&mut fig, cell, &gviz::Ideogram { label: &label, bands: &chr_bands, from: r.start, to: r.end, colour: r.colour, fontsize: label_size, lwd: 10.0 });
    }
    Ok(Canvas::from_figure(&fig))
}
