//! The parts of ggplot2 (4.0) the hmftools R scripts use, drawn into a display list: `theme_bw()` with base
//! size 11, continuous position scales (expansion, `extended` and log breaks, R's number formats), a panel
//! with grid lines and border, axes, titles, legends (keys and colour bars) on the right, and geoms (rect,
//! bar, point, segment, line, area, label). Sizes follow ggplot2's units: layout in points of 1/72.27 inch,
//! font sizes in big points (1/72 inch), line widths in mm (`* .pt` = 72.27 / 25.4, then R's lwd unit of
//! 1/96 inch).
//!
//! Layout is ggplot2's gtable: plot margin, title row, panel, axis rows and columns, axis titles and the guide
//! box, with every text height taken as the hinted ascent of "M" plus the descent of "gjpqyQ", as R's Cairo
//! device and ggplot2's `font_descent()` measure them.

use crate::fmath::FMath;

pub mod rfmt;

use crate::color::{GradientN, Rgba, r_colour};
use crate::display::{Figure, Item, LineCap, Stroke};
use crate::png_out::PngMeta;
use crate::text::{self, Font};

/// The device resolution of `ggsave()` (300 dpi).
pub const DPI: f64 = 300.0;

/// grid "points" (1/72.27 inch) to pixels
pub fn pt(v: f64) -> f64 {
    v * DPI / 72.27
}

/// font size in big points (1/72 inch) to pixels
pub fn font_px(size: f64) -> f64 {
    size * DPI / 72.0
}

/// a ggplot2 linewidth (mm) to pixels: lwd = linewidth * .pt, 1 lwd = 1/96 inch
pub fn lw_px(linewidth: f64) -> f64 {
    linewidth * 72.27 / 25.4 * DPI / 96.0
}

/// ggplot2's .pt
pub const PT: f64 = 72.27 / 25.4;
/// ggplot2's .stroke
pub const STROKE: f64 = 96.0 / 25.4;

/// theme_bw(base_size = 11)
pub mod theme {
    pub const HALF_LINE: f64 = 5.5;
    pub const TITLE_SIZE: f64 = 13.2;
    pub const AXIS_TITLE_SIZE: f64 = 11.0;
    pub const AXIS_TEXT_SIZE: f64 = 8.8;
    pub const LEGEND_TITLE_SIZE: f64 = 11.0;
    pub const LEGEND_TEXT_SIZE: f64 = 8.8;
    /// 1.2 lines of the default gpar (fontsize 12, lineheight 1.2), in big points
    pub const LEGEND_KEY_BIGPT: f64 = 1.2 * 12.0 * 1.2;
    pub const LINE_MM: f64 = 0.5;
}

/// The hinted ascent of "M" and descent of "gjpqyQ" at a font size (big points), in pixels.
pub fn text_metrics(size: f64) -> (f64, f64) {
    let px = font_px(size);
    (text::extents(Font::Sans, px, "M").ascent, text::extents(Font::Sans, px, "gjpqyQ").descent)
}

pub fn text_width(size: f64, s: &str) -> f64 {
    text::advance(Font::Sans, font_px(size), s)
}

/// The height ggplot2 gives a text row: ascent of M, plus line spacing for every extra line, plus the
/// descent.
pub fn text_height(size: f64, lines: usize) -> f64 {
    let (asc, desc) = text_metrics(size);
    asc + desc + (lines.saturating_sub(1)) as f64 * line_spacing(size)
}

/// grid's line spacing: lineheight 0.9 of 1.2 times the font size
pub fn line_spacing(size: f64) -> f64 {
    font_px(0.9 * 1.2 * size)
}

fn grey(n: u32) -> Rgba {
    r_colour(&format!("grey{n}")).unwrap()
}

/// Draws (possibly multi-line) left-aligned text whose block starts at `top`: the first baseline is the M
/// ascent below it, the next ones grid's line spacing apart.
pub fn draw_text_block(fig: &mut Figure, s: &str, x: f64, top: f64, size: f64, colour: Rgba) {
    let (asc, _) = text_metrics(size);
    let ls = line_spacing(size);
    for (i, line) in s.split('\n').enumerate() {
        fig.text(x.round(), (top + asc + i as f64 * ls).round(), line, Font::Sans, font_px(size), colour, 0.0);
    }
}

// -------------------------------------------------------------------------------------------------------------
// scales

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Trans {
    Identity,
    Log10,
}

impl Trans {
    pub fn apply(&self, v: f64) -> f64 {
        match self {
            Trans::Identity => v,
            Trans::Log10 => v.flog10(),
        }
    }
    pub fn inverse(&self, v: f64) -> f64 {
        match self {
            Trans::Identity => v,
            Trans::Log10 => 10f64.fpow(v),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Labels {
    /// `format(x, trim = TRUE)` of the in-range breaks (ggplot2's default)
    Auto,
    /// one label per break
    Fixed(Vec<String>),
    /// `format(x, scientific = FALSE)`: R's format without trimming, over the censored breaks
    FormatNoSci,
}

/// A continuous position scale.
#[derive(Clone, Debug)]
pub struct PosScale {
    pub trans: Trans,
    /// the limits in data space
    pub limits: (f64, f64),
    /// expansion: multiplicative and additive, on each side
    pub expand: (f64, f64),
    /// None: `extended` (or log) breaks of the expanded range
    pub breaks: Option<Vec<f64>>,
    pub labels: Labels,
    pub minor: bool,
}

impl PosScale {
    pub fn new(limits: (f64, f64)) -> PosScale {
        PosScale { trans: Trans::Identity, limits, expand: (0.05, 0.0), breaks: None, labels: Labels::Auto, minor: true }
    }

    /// The scale's view: the expanded range (transformed), the major breaks inside it with their labels, and
    /// the minor breaks.
    pub fn view(&self) -> View {
        let (a, b) = (self.trans.apply(self.limits.0), self.trans.apply(self.limits.1));
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let (m, add) = self.expand;
        let (r0, r1) = if hi > lo {
            (lo - (hi - lo) * m - add, hi + (hi - lo) * m + add)
        } else {
            // zero range: ggplot2 expands a single value by +-5% of its magnitude (or 0.5 for zero) -
            // expand_range(zero_width = 1)
            let w = if lo == 0.0 { 1.0 } else { lo.abs() * 0.1 };
            (lo - w / 2.0, hi + w / 2.0)
        };
        let breaks_data: Vec<f64> = match &self.breaks {
            Some(b) => b.clone(),
            None => match self.trans {
                Trans::Identity => rfmt::extended_breaks(r0, r1, 5),
                Trans::Log10 => rfmt::log_breaks(self.trans.inverse(r0), self.trans.inverse(r1), 5),
            },
        };
        let majors_t: Vec<f64> = breaks_data.iter().map(|&b| self.trans.apply(b)).collect();
        let inside = |v: f64| v >= r0 - 1e-10 * (r1 - r0).abs() && v <= r1 + 1e-10 * (r1 - r0).abs();
        // labels of the censored breaks
        let labels_all: Vec<Option<String>> = match &self.labels {
            Labels::Fixed(l) => majors_t.iter().zip(l).map(|(&v, s)| inside(v).then(|| s.clone())).collect(),
            Labels::Auto => {
                let vals: Vec<f64> = breaks_data.iter().zip(&majors_t).filter(|(_, t)| inside(**t)).map(|(b, _)| *b).collect();
                let f = rfmt::r_format(&vals, true);
                let mut it = f.into_iter();
                majors_t.iter().map(|&t| if inside(t) { it.next() } else { None }).collect()
            }
            Labels::FormatNoSci => {
                let vals: Vec<f64> = breaks_data.iter().zip(&majors_t).map(|(b, t)| if inside(*t) { *b } else { f64::NAN }).collect();
                let f = format_no_sci(&vals);
                f.into_iter().zip(&majors_t).map(|(s, &t)| inside(t).then_some(s)).collect()
            }
        };
        let mut majors = Vec::new();
        let mut labels = Vec::new();
        for (t, l) in majors_t.iter().zip(labels_all) {
            if inside(*t) {
                majors.push(*t);
                labels.push(l.unwrap_or_default());
            }
        }
        let minors = if self.minor { rfmt::regular_minor_breaks(&majors_t, (r0, r1)) } else { Vec::new() };
        View { range: (r0, r1), trans: self.trans, majors, labels, minors }
    }
}

/// R's `format(x, scientific = FALSE)` (NA for missing values), padded to a common width.
fn format_no_sci(xs: &[f64]) -> Vec<String> {
    let finite: Vec<f64> = xs.iter().copied().filter(|x| x.is_finite()).collect();
    let decimals = finite
        .iter()
        .map(|&x| {
            let s = rfmt::r_as_character(x);
            if s.contains('e') { 0 } else { s.split_once('.').map(|(_, d)| d.len()).unwrap_or(0) }
        })
        .max()
        .unwrap_or(0);
    let strs: Vec<String> = xs.iter().map(|&x| if x.is_finite() { format!("{:.*}", decimals, x) } else { "NA".into() }).collect();
    let w = strs.iter().map(|s| s.len()).max().unwrap_or(0);
    strs.into_iter().map(|s| format!("{s:>w$}")).collect()
}

/// A scale's view: what the panel and the axis draw.
#[derive(Clone, Debug)]
pub struct View {
    /// the expanded range, transformed
    pub range: (f64, f64),
    pub trans: Trans,
    /// major breaks (transformed) inside the range, and their labels
    pub majors: Vec<f64>,
    pub labels: Vec<String>,
    pub minors: Vec<f64>,
}

// -------------------------------------------------------------------------------------------------------------
// legends

#[derive(Clone, Debug)]
pub enum Key {
    /// draw_key_polygon: a filled square inset by the line width, outlined
    Polygon { fill: Rgba, colour: Option<Rgba>, linewidth: f64 },
    /// draw_key_point: a point of ggplot2 size `size`
    Point { colour: Rgba, size: f64 },
}

#[derive(Clone, Debug)]
pub enum Legend {
    Keys { title: Option<String>, entries: Vec<(String, Key)> },
    /// guide_colourbar: the ramp over `limits`, ticks and labels at `breaks`
    Bar { title: Option<String>, ramp: GradientN, limits: (f64, f64), breaks: Vec<f64>, labels: Vec<String> },
}

impl Legend {
    fn title(&self) -> Option<&str> {
        match self {
            Legend::Keys { title, .. } | Legend::Bar { title, .. } => title.as_deref(),
        }
    }

    fn title_size(&self) -> (f64, f64) {
        match self.title() {
            Some(t) => {
                let w = t.split('\n').map(|l| text_width(theme::LEGEND_TITLE_SIZE, l)).fold(0.0, f64::max);
                (w, text_height(theme::LEGEND_TITLE_SIZE, t.split('\n').count()) + pt(theme::HALF_LINE))
            }
            None => (0.0, 0.0),
        }
    }

    /// (width, height) of the legend, margins included.
    pub fn size(&self) -> (f64, f64) {
        let margin = pt(theme::HALF_LINE);
        let key = font_px(theme::LEGEND_KEY_BIGPT);
        let (tw, th) = self.title_size();
        let (body_w, body_h) = match self {
            Legend::Keys { entries, .. } => {
                let lw = entries.iter().map(|(l, _)| text_width(theme::LEGEND_TEXT_SIZE, l)).fold(0.0, f64::max);
                let kh = key.max(text_height(theme::LEGEND_TEXT_SIZE, 1));
                (key + margin + lw, kh * entries.len() as f64)
            }
            Legend::Bar { labels, .. } => {
                let lw = labels.iter().map(|l| text_width(theme::LEGEND_TEXT_SIZE, l)).fold(0.0, f64::max);
                (key + margin + lw, 5.0 * key)
            }
        };
        (2.0 * margin + body_w.max(tw), 2.0 * margin + th + body_h)
    }

    fn draw(&self, fig: &mut Figure, x0: f64, y0: f64) {
        let margin = pt(theme::HALF_LINE);
        let key = font_px(theme::LEGEND_KEY_BIGPT);
        let (_, th) = self.title_size();
        if let Some(t) = self.title() {
            let lines = t.split('\n').count();
            let h = text_height(theme::LEGEND_TITLE_SIZE, lines);
            let (_, desc) = text_metrics(theme::LEGEND_TITLE_SIZE);
            // vjust 0.5 in the text cell: the block (M ascent to last baseline) centred, the descent below
            let block = h - desc;
            draw_text_block(fig, t, x0 + margin, y0 + margin + (h - block) / 2.0, theme::LEGEND_TITLE_SIZE, Rgba::BLACK);
        }
        let top = y0 + margin + th;
        let (asc, _) = text_metrics(theme::LEGEND_TEXT_SIZE);
        match self {
            Legend::Keys { entries, .. } => {
                let kh = key.max(text_height(theme::LEGEND_TEXT_SIZE, 1));
                for (i, (label, k)) in entries.iter().enumerate() {
                    let (kx, ky) = (x0 + margin, top + i as f64 * kh + (kh - key) / 2.0);
                    match k {
                        Key::Polygon { fill, colour, linewidth } => {
                            let inset = linewidth / 2.0 * DPI / 25.4; // unit(lwd, "mm") off each dimension
                            let (x, y, s) = (kx + inset, ky + inset, key - 2.0 * inset);
                            fill_rect_r(fig, x, y, x + s, y + s, *fill);
                            if let Some(c) = colour {
                                if *linewidth > 0.0 {
                                    fig.rect_outline(x, y, s, s, Stroke { color: *c, width: lw_px(*linewidth), dash: Vec::new(), cap: LineCap::Round });
                                }
                            }
                        }
                        Key::Point { colour, size } => point(fig, kx + key / 2.0, ky + key / 2.0, *size, *colour, 0.5),
                    }
                    let cy = top + i as f64 * kh + kh / 2.0;
                    fig.text(x0 + margin + key + margin, (cy + asc / 2.0).round(), label, Font::Sans, font_px(theme::LEGEND_TEXT_SIZE), Rgba::BLACK, 0.0);
                }
            }
            Legend::Bar { ramp, limits, breaks, labels, .. } => {
                let (bx, by, bw, bh) = (x0 + margin, top, key, 5.0 * key);
                // the bar: nbin = 300 colours over the limits, bottom to top, drawn a pixel row at a time
                // (the raster ggplot2 draws, interpolated)
                let nbin = 300.0;
                let (ry0, ry1) = (by.round() as i64, (by + bh).round() as i64);
                for row in ry0..ry1 {
                    let f = 1.0 - ((row as f64 + 0.5) - by) / bh; // 0 at the bottom
                    let bin = (f * nbin - 0.5).clamp(0.0, nbin - 1.0);
                    let v = limits.0 + (limits.1 - limits.0) * bin / (nbin - 1.0);
                    let x = (v - limits.0) / (limits.1 - limits.0);
                    if let Some(c) = ramp.at(x.clamp(0.0, 1.0)) {
                        fig.rect(bx, row as f64, bw, 1.0, c);
                    }
                }
                let tick_len = 0.2 * bw;
                let tick = Stroke { color: Rgba::WHITE, width: 0.5 * DPI / 96.0, dash: Vec::new(), cap: LineCap::Butt };
                for (b, l) in breaks.iter().zip(labels) {
                    // the limits sit at the centres of the end bins
                    let f = (0.5 + (b - limits.0) / (limits.1 - limits.0) * (nbin - 1.0)) / nbin;
                    let y = by + bh * (1.0 - f);
                    fig.line(bx, y, bx + tick_len, y, tick.clone());
                    fig.line(bx + bw - tick_len, y, bx + bw, y, tick.clone());
                    fig.text(bx + bw + margin, (y + asc / 2.0).round(), l, Font::Sans, font_px(theme::LEGEND_TEXT_SIZE), Rgba::BLACK, 0.0);
                }
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------------------
// the plot

/// A plot: title, scales and axis titles, legends, theme switches.
#[derive(Clone, Debug)]
pub struct Plot {
    pub title: Option<String>,
    pub x: View,
    pub y: View,
    /// axis titles: None removes the row, Some("") keeps its space (`xlab("")`)
    pub x_title: Option<String>,
    pub y_title: Option<String>,
    pub legends: Vec<Legend>,
    pub ticks: bool,
    pub minor_grid: bool,
    pub panel_border: bool,
}

/// The panel's place and its data-to-pixel mapping.
#[derive(Clone, Copy, Debug)]
pub struct Panel {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    pub xt: Trans,
    pub yt: Trans,
}

impl Panel {
    pub fn px(&self, x: f64) -> f64 {
        let v = self.xt.apply(x);
        self.x0 + (v - self.xr.0) / (self.xr.1 - self.xr.0) * (self.x1 - self.x0)
    }
    pub fn py(&self, y: f64) -> f64 {
        let v = self.yt.apply(y);
        self.y1 - (v - self.yr.0) / (self.yr.1 - self.yr.0) * (self.y1 - self.y0)
    }
    fn clip(&self) -> (f64, f64, f64, f64) {
        (self.x0, self.y0, self.x1 - self.x0, self.y1 - self.y0)
    }
}

/// The fixed widths left of the panel: the y title column and the y axis column.
fn left_widths(p: &Plot) -> f64 {
    let mut w = 0.0;
    if p.y_title.is_some() {
        w += text_height(theme::AXIS_TITLE_SIZE, 1) + pt(theme::HALF_LINE / 2.0);
    }
    if !p.y.labels.is_empty() {
        let lw = p.y.labels.iter().map(|l| text_width(theme::AXIS_TEXT_SIZE, l)).fold(0.0, f64::max);
        w += lw + pt(0.8 * theme::HALF_LINE / 2.0) + if p.ticks { pt(theme::HALF_LINE / 2.0) } else { 0.0 };
    }
    w
}

impl Plot {
    /// The width of the plot's left columns (for aligning stacked plots, cowplot's align = "v").
    pub fn left_width(&self) -> f64 {
        left_widths(self)
    }

    /// Draws the plot into the box (x, y, w, h) of `fig`; `left` forces the width of the left columns.
    /// The layers are drawn by `layers` on the panel (clipped to it).
    pub fn draw(&self, fig: &mut Figure, bx: f64, by: f64, bw: f64, bh: f64, left: Option<f64>, layers: &mut dyn FnMut(&mut Figure, &Panel)) -> Panel {
        let margin = pt(theme::HALF_LINE);
        let left_w = left.unwrap_or_else(|| left_widths(self));
        let mut right = margin;
        let legend_sizes: Vec<(f64, f64)> = self.legends.iter().map(|l| l.size()).collect();
        let legend_w = legend_sizes.iter().map(|s| s.0).fold(0.0, f64::max);
        if !self.legends.is_empty() {
            right += pt(2.0 * theme::HALF_LINE) + legend_w;
        }
        let mut top = margin;
        if self.title.is_some() {
            top += text_height(theme::TITLE_SIZE, 1) + pt(theme::HALF_LINE);
        }
        let tick = if self.ticks { pt(theme::HALF_LINE / 2.0) } else { 0.0 };
        let axis_b = if self.x.labels.is_empty() { 0.0 } else { text_height(theme::AXIS_TEXT_SIZE, 1) + pt(0.8 * theme::HALF_LINE / 2.0) + tick };
        let xlab_h = if self.x_title.is_some() { text_height(theme::AXIS_TITLE_SIZE, 1) + pt(theme::HALF_LINE / 2.0) } else { 0.0 };
        let bottom = margin + axis_b + xlab_h;
        let panel = Panel {
            x0: bx + margin + left_w,
            y0: by + top,
            x1: bx + bw - right,
            y1: by + bh - bottom,
            xr: self.x.range,
            yr: self.y.range,
            xt: self.x.trans,
            yt: self.y.trans,
        };

        // panel: grid, layers, border
        fig.clip(Some(panel.clip()));
        let gl = |w: f64| Stroke { color: grey(92), width: lw_px(w), dash: Vec::new(), cap: LineCap::Butt };
        if self.minor_grid {
            for &v in &self.x.minors {
                let x = panel.x0 + (v - panel.xr.0) / (panel.xr.1 - panel.xr.0) * (panel.x1 - panel.x0);
                fig.line(x, panel.y0, x, panel.y1, gl(theme::LINE_MM / 2.0));
            }
            for &v in &self.y.minors {
                let y = panel.y1 - (v - panel.yr.0) / (panel.yr.1 - panel.yr.0) * (panel.y1 - panel.y0);
                fig.line(panel.x0, y, panel.x1, y, gl(theme::LINE_MM / 2.0));
            }
        }
        let xs: Vec<f64> = self.x.majors.iter().map(|v| panel.x0 + (v - panel.xr.0) / (panel.xr.1 - panel.xr.0) * (panel.x1 - panel.x0)).collect();
        let ys: Vec<f64> = self.y.majors.iter().map(|v| panel.y1 - (v - panel.yr.0) / (panel.yr.1 - panel.yr.0) * (panel.y1 - panel.y0)).collect();
        for &x in &xs {
            fig.line(x, panel.y0, x, panel.y1, gl(theme::LINE_MM));
        }
        for &y in &ys {
            fig.line(panel.x0, y, panel.x1, y, gl(theme::LINE_MM));
        }
        layers(fig, &panel);
        fig.clip(Some(panel.clip()));
        if self.panel_border {
            // a 0.5 mm grey20 line on the panel's edge, clipped to the panel: its inner half shows
            let h = lw_px(theme::LINE_MM) / 2.0;
            let c = grey(20);
            let (x0, y0, x1, y1) = (panel.x0, panel.y0, panel.x1, panel.y1);
            fig.rect(x0, y0, x1 - x0, h, c);
            fig.rect(x0, y1 - h, x1 - x0, h, c);
            fig.rect(x0, y0 + h, h, y1 - y0 - 2.0 * h, c);
            fig.rect(x1 - h, y0 + h, h, y1 - y0 - 2.0 * h, c);
        }
        fig.clip(None);

        // axes
        let (asc8, _) = text_metrics(theme::AXIS_TEXT_SIZE);
        let axis_colour = grey(30);
        let tick_stroke = Stroke { color: grey(20), width: lw_px(theme::LINE_MM), dash: Vec::new(), cap: LineCap::Butt };
        for (x, l) in xs.iter().zip(&self.x.labels) {
            if self.ticks {
                fig.line(*x, panel.y1, *x, panel.y1 + tick, tick_stroke.clone());
            }
            let w = text_width(theme::AXIS_TEXT_SIZE, l);
            let base = panel.y1 + tick + pt(0.8 * theme::HALF_LINE / 2.0) + asc8;
            fig.text((x - w / 2.0).round(), base.round(), l, Font::Sans, font_px(theme::AXIS_TEXT_SIZE), axis_colour, 0.0);
        }
        for (y, l) in ys.iter().zip(&self.y.labels) {
            if self.ticks {
                fig.line(panel.x0 - tick, *y, panel.x0, *y, tick_stroke.clone());
            }
            let w = text_width(theme::AXIS_TEXT_SIZE, l);
            let right_edge = panel.x0 - tick - pt(0.8 * theme::HALF_LINE / 2.0);
            fig.text((right_edge - w).round(), (y + asc8 / 2.0).round(), l, Font::Sans, font_px(theme::AXIS_TEXT_SIZE), axis_colour, 0.0);
        }
        // axis titles
        let (asc11, _) = text_metrics(theme::AXIS_TITLE_SIZE);
        if let Some(t) = self.x_title.as_deref().filter(|t| !t.is_empty()) {
            let w = text_width(theme::AXIS_TITLE_SIZE, t);
            let base = panel.y1 + axis_b + pt(theme::HALF_LINE / 2.0) + asc11;
            fig.text(((panel.x0 + panel.x1) / 2.0 - w / 2.0).round(), base.round(), t, Font::Sans, font_px(theme::AXIS_TITLE_SIZE), Rgba::BLACK, 0.0);
        }
        if let Some(t) = self.y_title.as_deref().filter(|t| !t.is_empty()) {
            let w = text_width(theme::AXIS_TITLE_SIZE, t);
            let base_x = bx + margin + asc11;
            let cy = (panel.y0 + panel.y1) / 2.0;
            fig.text(base_x.round(), (cy + w / 2.0).round(), t, Font::Sans, font_px(theme::AXIS_TITLE_SIZE), Rgba::BLACK, 90.0);
        }
        // title: left-aligned on the panel
        if let Some(t) = &self.title {
            let (asc, _) = text_metrics(theme::TITLE_SIZE);
            fig.text(panel.x0.round(), (by + margin + asc).round(), t, Font::Sans, font_px(theme::TITLE_SIZE), Rgba::BLACK, 0.0);
        }
        // legends: stacked, spaced by legend.spacing.y (11 pt), left-aligned, centred on the panel
        if !self.legends.is_empty() {
            let spacing = pt(2.0 * theme::HALF_LINE);
            let total_h: f64 = legend_sizes.iter().map(|s| s.1).sum::<f64>() + spacing * (self.legends.len() - 1) as f64;
            let x = panel.x1 + pt(2.0 * theme::HALF_LINE);
            let mut y = (panel.y0 + panel.y1) / 2.0 - total_h / 2.0;
            for (l, s) in self.legends.iter().zip(&legend_sizes) {
                l.draw(fig, x, y);
                y += s.1 + spacing;
            }
        }
        panel
    }
}

/// A figure of `w` x `h` inches at 300 dpi, white, as ggsave writes it.
pub fn figure(w_in: f64, h_in: f64) -> Figure {
    let mut f = Figure::new((w_in * DPI).round() as u32, (h_in * DPI).round() as u32, Rgba::WHITE);
    f.meta = PngMeta { dpi: Some(DPI), transparent_key: None };
    f
}

// -------------------------------------------------------------------------------------------------------------
// geoms

/// geom_point with shape 19: a disc of radius 0.375 * fontsize, outlined with lwd = stroke * .stroke / 2,
/// fill and outline in the same colour (R draws both, so a translucent point's rim is darker).
pub fn point(fig: &mut Figure, x: f64, y: f64, size: f64, colour: Rgba, stroke: f64) {
    let fs = size * PT + stroke * STROKE / 2.0;
    let r = 0.375 * font_px(fs);
    let lwd = stroke * STROKE / 2.0;
    let lw = lwd * DPI / 96.0;
    fig.push(Item::Circle { cx: x, cy: y, r, fill: Some(colour), stroke: Some(Stroke { color: colour, width: lw, dash: Vec::new(), cap: LineCap::Round }) });
}

/// R's Cairo device fills rectangles without anti-aliasing: the pixels whose centres lie inside.
pub fn fill_rect_r(fig: &mut Figure, x0: f64, y0: f64, x1: f64, y1: f64, colour: Rgba) {
    let snap = |v: f64| (v - 0.5).ceil();
    let (a, b) = (snap(x0.min(x1)), snap(x0.max(x1)));
    let (c, d) = (snap(y0.min(y1)), snap(y0.max(y1)));
    fig.rect(a, c, b - a, d - c, colour);
}

/// geom_rect / geom_bar: fill (pixel centres, as R's Cairo device), then an anti-aliased outline of
/// `linewidth` mm centred on the edges (mitre joins).
pub fn rect(fig: &mut Figure, p: &Panel, xmin: f64, xmax: f64, ymin: f64, ymax: f64, fill: Option<Rgba>, outline: Option<(Rgba, f64)>) {
    let (x0, x1) = (p.px(xmin).min(p.px(xmax)), p.px(xmin).max(p.px(xmax)));
    let (y0, y1) = (p.py(ymax).min(p.py(ymin)), p.py(ymax).max(p.py(ymin)));
    if let Some(f) = fill {
        fill_rect_r(fig, x0, y0, x1, y1, f);
    }
    if let Some((c, lw)) = outline {
        if lw > 0.0 && !c.is_transparent() {
            let s = Stroke { color: c, width: lw_px(lw), dash: Vec::new(), cap: LineCap::Round };
            if y1 - y0 <= 0.0 {
                fig.line(x0, y0, x1, y0, s);
            } else if x1 - x0 <= 0.0 {
                fig.line(x0, y0, x0, y1, s);
            } else {
                fig.rect_outline(x0, y0, x1 - x0, y1 - y0, s);
            }
        }
    }
}

/// geom_segment / geom_line: a polyline of `linewidth` mm; dashed is R's lty "44" (on and off 4 times the
/// line width, at least 1, in lwd units).
pub fn path(fig: &mut Figure, pts: Vec<(f64, f64)>, colour: Rgba, linewidth: f64, dashed: bool) {
    let lwd = linewidth * PT;
    let dash = if dashed {
        let d = 4.0 * lwd.max(1.0) * DPI / 96.0;
        vec![d, d]
    } else {
        Vec::new()
    };
    fig.polyline(pts, Stroke { color: colour, width: lw_px(linewidth), dash, cap: LineCap::Butt });
}

/// geom_label: a rounded box (padding 0.25 lines, radius 0.15 lines, border 0.25 mm) with the text inside.
#[allow(clippy::too_many_arguments)]
pub fn label(fig: &mut Figure, x: f64, y: f64, s: &str, size_mm: f64, hjust: f64, vjust: f64) {
    let size = size_mm * PT;
    let lines = font_px(size * 1.2);
    let pad = 0.25 * lines;
    let radius = 0.15 * lines;
    let (asc, desc) = text_metrics(size);
    let w = text_width(size, s) + 2.0 * pad;
    let h = asc + desc + 2.0 * pad;
    let (bx, by) = (x - hjust * w, y - (1.0 - vjust) * h);
    let border = Stroke { color: Rgba::BLACK, width: lw_px(0.25), dash: Vec::new(), cap: LineCap::Round };
    let pts = rounded_rect(bx, by, w, h, radius);
    fig.polygon(pts, Some(Rgba::WHITE), Some(border));
    fig.text((bx + pad).round(), (by + pad + asc).round(), s, Font::Sans, font_px(size), Rgba::BLACK, 0.0);
}

fn rounded_rect(x: f64, y: f64, w: f64, h: f64, r: f64) -> Vec<(f64, f64)> {
    let r = r.min(w / 2.0).min(h / 2.0);
    let mut pts = Vec::new();
    let corner = |pts: &mut Vec<(f64, f64)>, cx: f64, cy: f64, a0: f64| {
        for i in 0..=8 {
            let a = (a0 + 90.0 * i as f64 / 8.0).to_radians();
            pts.push((cx + r * a.fcos(), cy + r * a.fsin()));
        }
    };
    corner(&mut pts, x + w - r, y + r, -90.0);
    corner(&mut pts, x + w - r, y + h - r, 0.0);
    corner(&mut pts, x + r, y + h - r, 90.0);
    corner(&mut pts, x + r, y + r, 180.0);
    pts
}

/// ggplot2's stat_bin bins for a binwidth with neither centre nor boundary given (boundary = width / 2),
/// closed on the right, over `range`: the bin edges.
pub fn bin_breaks(range: (f64, f64), width: f64) -> Vec<f64> {
    let boundary = width / 2.0;
    let shift = ((range.0 - boundary) / width).floor();
    let origin = boundary + shift * width;
    let max_x = range.1 + (1.0 - 1e-8) * width;
    // seq(origin, max_x, width)
    let n = ((max_x - origin) / width + 1e-10).floor() as i64;
    let mut b: Vec<f64> = (0..=n).map(|i| origin + i as f64 * width).collect();
    if b.len() == 1 {
        b.push(b[0] + width);
    }
    b
}

/// The bin of x for stat_bin's right-closed bins (with the first bin closed on both sides, and the fuzz of
/// 1e-8 times the median width).
pub fn bin_index(breaks: &[f64], x: f64) -> Option<usize> {
    let n = breaks.len();
    if n < 2 {
        return None;
    }
    let fuzz = 1e-8 * (breaks[1] - breaks[0]);
    // fuzzy breaks: first moved down, the rest up (right-closed)
    let f = |i: usize| if i == 0 { breaks[0] - fuzz } else { breaks[i] + fuzz };
    if x < f(0) || x > f(n - 1) {
        return None;
    }
    (0..n - 1).find(|&i| x <= f(i + 1) && (i == 0 || x > f(i)))
}

#[cfg(test)]
mod tests;
