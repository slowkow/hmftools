//! ggplot2's plot tables, built as ggplot2 builds them (ggplot_gtable: the facet's panel table with axes and strips,
//! axis titles, the five guide-box cells, title, subtitle and caption rows, tag padding and the plot margin), so that
//! patchwork's alignment (`crate::patchwork`) works on the same table structure as in R. Sizes come from R's
//! devices (`crate::rdev`), so the same code lays out a PNG (Cairo, DejaVu Sans) and a PDF (Helvetica).
//!
//! This is a port of the parts of ggplot2 3.5.1 (with ggh4x 0.2.8's facet_grid2 and strips) and 4.0.3 that
//! CUPPA's `plot_predictions.R` and qsee's `plot_qc.R` use; where the two versions differ, [`Version`] says which.
//! The module [`gg`](crate::gg) (PURPLE's plots) is separate and unchanged.
//!
//! - this file: the context, theme elements, and text as R draws it (GEText's multi-line justification and rotation,
//!   titleGrob's margins and descent);
//! - [`guides`]: axes, colour bars and key legends, the guide box;
//! - [`facet`]: panel tables (facet_null, facet_grid with ggh4x's strips, facet_wrap);
//! - [`plot`]: the plot table;
//! - [`geom`]: geoms drawn in a panel.

pub mod facet;
pub mod geom;
pub mod guides;
pub mod plot;
pub mod plotmath;
pub mod scales;

use crate::color::Rgba;
use crate::display::{Item, LineCap, Stroke};
use crate::gtable::{Bx, Canvas, Grob, Unit};
use crate::rdev::{Device, Face};

/// The ggplot2 release whose layout rules apply.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Version {
    /// ggplot2 3.5.1 with patchwork 1.2.0 and ggh4x 0.2.8 (CUPPA 2.5.1's image)
    V351,
    /// ggplot2 4.0.3 with patchwork 1.3.2 (qsee's image)
    V403,
}

/// The device and the ggplot2 version a figure is laid out for.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ctx {
    pub dev: Device,
    pub v: Version,
}

impl Ctx {
    pub fn pt(&self, v: f64) -> f64 {
        self.dev.pt(v)
    }
    pub fn mm(&self, v: f64) -> f64 {
        self.dev.mm(v)
    }
    /// A line width of `linewidth` mm (ggplot2) in device units: lwd = linewidth * .pt, 1 lwd = 1/96 inch.
    pub fn lwd_mm(&self, linewidth: f64) -> f64 {
        linewidth * PT * self.dev.dpi() / 96.0
    }
    /// grid's lwd in device units
    pub fn lwd(&self, lwd: f64) -> f64 {
        lwd * self.dev.dpi() / 96.0
    }
    /// A font size (big points) in device units
    pub fn font_px(&self, size: f64) -> f64 {
        self.dev.bigpt(size)
    }
}

thread_local! {
    static DRAWN_BOXES: std::cell::RefCell<Option<Vec<(String, Bx)>>> = const { std::cell::RefCell::new(None) };
}

/// For layout tests against R's viewports: records, while on, the boxes of what is drawn without a table cell of its
/// own but in a viewport of its own in R (geom_label's and draw_key_label's rounded rectangles: R's "rrvp").
/// `record_drawn_boxes(true)` clears and starts; [`take_drawn_boxes`] stops and returns them (device units, y down).
pub fn record_drawn_boxes(on: bool) {
    DRAWN_BOXES.with(|d| *d.borrow_mut() = on.then(Vec::new));
}

/// Records a drawn box (when recording is on).
pub fn drawn_box(name: &str, b: Bx) {
    DRAWN_BOXES.with(|d| {
        if let Some(v) = d.borrow_mut().as_mut() {
            v.push((name.to_string(), b));
        }
    });
}

/// Stops recording and returns the boxes recorded.
pub fn take_drawn_boxes() -> Vec<(String, Bx)> {
    DRAWN_BOXES.with(|d| d.borrow_mut().take().unwrap_or_default())
}

/// ggplot2's .pt (points per mm)
pub const PT: f64 = 72.27 / 25.4;
/// ggplot2's .stroke
pub const STROKE: f64 = 96.0 / 25.4;

/// A text element (element_text), resolved. Sizes in big points, margins in device units (top, right, bottom,
/// left).
#[derive(Clone, Debug, PartialEq)]
pub struct TextEl {
    pub size: f64,
    pub face: Face,
    pub colour: Rgba,
    pub hjust: f64,
    pub vjust: f64,
    pub angle: f64,
    pub lineheight: f64,
    pub margin: [f64; 4],
}

impl TextEl {
    pub fn new(size: f64) -> TextEl {
        TextEl { size, face: Face::Plain, colour: Rgba::BLACK, hjust: 0.5, vjust: 0.5, angle: 0.0, lineheight: 0.9, margin: [0.0; 4] }
    }
}

/// A rectangle element (element_rect): fill, outline colour and width (mm).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectEl {
    pub fill: Option<Rgba>,
    pub colour: Option<Rgba>,
    pub linewidth: f64,
}

/// A line element (element_line): colour, width (mm), dash pattern (R's hex digits in lwd units).
#[derive(Clone, Debug, PartialEq)]
pub struct LineEl {
    pub colour: Rgba,
    pub linewidth: f64,
    pub dash: Vec<f64>,
}

impl LineEl {
    pub fn solid(colour: Rgba, linewidth: f64) -> LineEl {
        LineEl { colour, linewidth, dash: Vec::new() }
    }
}

/// R's dash pattern for a line type string of hex digits ("44", "11", "A1") at a line width (lwd) on Cairo's png():
/// each digit times max(lwd, 1), in 1/96 inch (CairoLineType). See [`dash_pattern_on`] for the device of a context.
pub fn dash_pattern(lty: &str, lwd: f64, dpi: f64) -> Vec<f64> {
    lty.chars().filter_map(|c| c.to_digit(16)).map(|d| d as f64 * lwd.max(1.0) * dpi / 96.0).collect()
}

/// R's dash pattern of a butt-ended line on the device of `ctx`: Cairo's png() scales the digits by max(lwd, 1)
/// ([`dash_pattern`]), pdf() by lwd itself (PDFSetLineTexture: digit * lwd * 0.75 pt; no cap correction for butt
/// ends), so a dotted line thinner than lwd 1 has shorter dashes in R's PDF than in its PNG.
pub fn dash_pattern_on(ctx: &Ctx, lty: &str, lwd: f64) -> Vec<f64> {
    match ctx.dev {
        Device::Pdf => lty.chars().filter_map(|c| c.to_digit(16)).map(|d| d as f64 * lwd * 0.75).collect(),
        _ => dash_pattern(lty, lwd, ctx.dev.dpi()),
    }
}

/// rotate_just: hjust and vjust turned into the justification of the rotated text's reference point.
pub fn rotate_just(angle: f64, hjust: f64, vjust: f64) -> (f64, f64) {
    let a = angle.rem_euclid(360.0);
    if a < 90.0 {
        (hjust, vjust)
    } else if a < 180.0 {
        (1.0 - vjust, hjust)
    } else if a < 270.0 {
        (1.0 - hjust, 1.0 - vjust)
    } else {
        (vjust, 1.0 - hjust)
    }
}

// -------------------------------------------------------------------------------------------------------------
// text

/// Draws a text at (x, y) (device units, y down) as R's graphics engine does (GEText): each line justified by its
/// own width and the ascent of "M", lines stacked by the line spacing with the block justified by vjust, the whole
/// rotated about (x, y) by `angle` degrees counter-clockwise.
#[allow(clippy::too_many_arguments)]
pub fn draw_text(c: &mut Canvas, ctx: &Ctx, x: f64, y: f64, label: &str, size: f64, face: Face, colour: Rgba, hjust: f64, vjust: f64, angle: f64, lineheight: f64) {
    if label.is_empty() || colour.is_transparent() {
        return;
    }
    let dev = ctx.dev;
    let lines: Vec<&str> = label.split('\n').collect();
    let n = lines.len() as f64;
    let rot = angle.to_radians();
    let (sn, cs) = (rot.sin(), rot.cos());
    let ls = dev.line_spacing(size, lineheight);
    let h = dev.m_ascent(size, face);
    for (i, line) in lines.iter().enumerate() {
        // offset of this line's reference point, upwards in the text's frame
        let yoff = if n > 1.0 { ((1.0 - vjust) * (n - 1.0) - i as f64) * ls } else { 0.0 };
        // in device space (y up): rotate (0, yoff)
        let (ox, oy) = (-yoff * sn, yoff * cs);
        let w = dev.str_width(line, size, face);
        // bottom left of the line (y up)
        let xl = ox - hjust * w * cs + vjust * h * sn;
        let yb = oy - hjust * w * sn - vjust * h * cs;
        let s = dev.shape(line);
        if s.trim().is_empty() {
            continue;
        }
        c.fig.push(Item::Text { x: x + xl, y: y - yb, text: s, font: dev.font(face), size: ctx.font_px(size), color: colour, angle });
    }
}

/// The bounding box of a text's (rotated) box relative to its anchor, as grid computes a textGrob's extent: the
/// box is [0, width] x [0, height] (height: GEStrHeight, the "M" ascent plus the extra lines) justified by hjust
/// and vjust and rotated. Returns (xmin, xmax, ymin, ymax), y up.
pub fn text_extent(ctx: &Ctx, label: &str, size: f64, face: Face, lineheight: f64, hjust: f64, vjust: f64, angle: f64) -> (f64, f64, f64, f64) {
    let w = ctx.dev.text_width(label, size, face);
    let h = ctx.dev.text_height(label, size, face, lineheight);
    let rot = angle.to_radians();
    let (sn, cs) = (rot.sin(), rot.cos());
    let (mut x0, mut x1, mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for (px, py) in [(0.0, 0.0), (w, 0.0), (w, h), (0.0, h)] {
        let (dx, dy) = (px - hjust * w, py - vjust * h);
        let (rx, ry) = (dx * cs - dy * sn, dx * sn + dy * cs);
        x0 = x0.min(rx);
        x1 = x1.max(rx);
        y0 = y0.min(ry);
        y1 = y1.max(ry);
    }
    (x0, x1, y0, y1)
}

/// ggplot2's titleGrob: labels with an element, at positions along x or y (npc of the cell) or at the rotated
/// justification, with margins on x and/or y. Its size counts the text's descent (font_descent) as ggplot2 does.
#[derive(Clone, Debug)]
pub struct TitleGrob {
    pub labels: Vec<String>,
    /// positions (npc) of the labels along x, or None for the rotated hjust
    pub x: Option<Vec<f64>>,
    /// positions (npc) along y, or None for the rotated vjust
    pub y: Option<Vec<f64>>,
    pub el: TextEl,
    pub margin_x: bool,
    pub margin_y: bool,
}

impl TitleGrob {
    pub fn new(label: &str, el: &TextEl, margin_x: bool, margin_y: bool) -> TitleGrob {
        TitleGrob { labels: vec![label.to_string()], x: None, y: None, el: el.clone(), margin_x, margin_y }
    }

    /// The text box of all labels together (grobWidth, grobHeight of the textGrob), with the descent.
    fn text_size(&self, ctx: &Ctx) -> (f64, f64) {
        if self.labels.is_empty() {
            return (0.0, 0.0);
        }
        let el = &self.el;
        let (mut x0, mut x1, mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for l in &self.labels {
            let (a, b, cc, d) = text_extent(ctx, l, el.size, el.face, el.lineheight, el.hjust, el.vjust, el.angle);
            x0 = x0.min(a);
            x1 = x1.max(b);
            y0 = y0.min(cc);
            y1 = y1.max(d);
        }
        // labels at different positions along an axis: the extent along it is not used by ggplot2 (null), only the
        // orthogonal one; the union of boxes at a common anchor is what grid gives for that dimension
        // font_descent(gp$fontfamily, gp$fontface, ...): grid's gpar keeps a face as `font`, so gp$fontface is NULL
        // and the descent is the plain face's, whatever the element's face
        let descent = ctx.dev.font_descent(el.size, Face::Plain);
        let rad = el.angle.rem_euclid(360.0).to_radians();
        let w = (x1 - x0) + rad.sin().abs() * descent;
        let h = (y1 - y0) + rad.cos().abs() * descent;
        (w, h)
    }

    /// The width of the grob: margins and text with margin_x; a null unit when only margin_y is set.
    pub fn width(&self, ctx: &Ctx) -> Unit {
        if self.labels.is_empty() {
            return Unit::ZERO;
        }
        if self.margin_x {
            let (w, _) = self.text_size(ctx);
            Unit::Abs(self.el.margin[3] + w + self.el.margin[1])
        } else if self.margin_y {
            Unit::Null(1.0)
        } else {
            Unit::Abs(self.text_size(ctx).0)
        }
    }

    pub fn height(&self, ctx: &Ctx) -> Unit {
        if self.labels.is_empty() {
            return Unit::ZERO;
        }
        if self.margin_y {
            let (_, h) = self.text_size(ctx);
            Unit::Abs(self.el.margin[0] + h + self.el.margin[2])
        } else if self.margin_x {
            Unit::Null(1.0)
        } else {
            Unit::Abs(self.text_size(ctx).1)
        }
    }

    /// Draws the labels in their cell `b`.
    pub fn draw(&self, c: &mut Canvas, ctx: &Ctx, b: Bx) {
        let el = &self.el;
        let (hj, vj) = rotate_just(el.angle, el.hjust, el.vjust);
        let [mt, mr, mb, ml] = el.margin;
        let n = self.labels.len();
        for (i, label) in self.labels.iter().enumerate() {
            // x: npc of the cell, then the margin shift (towards the inside of the margins)
            let xn = self.x.as_ref().map(|v| v[i.min(v.len() - 1)]).unwrap_or(hj);
            let yn = self.y.as_ref().map(|v| v[i.min(v.len() - 1)]).unwrap_or(vj);
            let mut x = b.npc_x(xn);
            let mut y_up = yn * b.h();
            if self.margin_x {
                x += -mr * hj + ml * (1.0 - hj);
            }
            if self.margin_y {
                y_up += -mt * vj + mb * (1.0 - vj);
            }
            let _ = n;
            draw_text(c, ctx, x, b.y1 - y_up, label, el.size, el.face, el.colour, el.hjust, el.vjust, el.angle, el.lineheight);
        }
    }

    pub fn grob(self, ctx: Ctx) -> Grob {
        if self.labels.is_empty() {
            return Grob::Zero;
        }
        Grob::draw_fn(move |c, b| self.draw(c, &ctx, b))
    }
}

// -------------------------------------------------------------------------------------------------------------
// shapes

/// A stroke of `lwd` (grid's line width) in a colour, with R's line ends.
pub fn stroke(ctx: &Ctx, colour: Rgba, lwd: f64, dash: Vec<f64>, cap: LineCap) -> Stroke {
    Stroke { color: colour, width: ctx.lwd(lwd), dash, cap }
}

/// Fills a rectangle as the device does: R's Cairo device fills the pixels whose centres are inside (no
/// anti-aliasing); a PDF takes the exact box.
pub fn fill_rect(c: &mut Canvas, ctx: &Ctx, x0: f64, y0: f64, x1: f64, y1: f64, fill: Rgba) {
    let (x0, x1) = (x0.min(x1), x0.max(x1));
    let (y0, y1) = (y0.min(y1), y0.max(y1));
    match ctx.dev {
        Device::Png { .. } | Device::PdfAt { .. } => {
            let snap = |v: f64| (v - 0.5).ceil();
            let (a, b) = (snap(x0), snap(x1));
            let (cc, d) = (snap(y0), snap(y1));
            c.fig.rect(a, cc, b - a, d - cc, fill);
        }
        Device::Pdf => c.fig.rect(x0, y0, x1 - x0, y1 - y0, fill),
    }
}

/// A rectangle (rectGrob): filled, then outlined with a stroke of `lwd` centred on its edges (mitre joins).
#[allow(clippy::too_many_arguments)]
pub fn rect(c: &mut Canvas, ctx: &Ctx, x0: f64, y0: f64, x1: f64, y1: f64, fill: Option<Rgba>, outline: Option<(Rgba, f64)>, dash: Vec<f64>) {
    if let Some(f) = fill {
        if !f.is_transparent() {
            fill_rect(c, ctx, x0, y0, x1, y1, f);
        }
    }
    if let Some((col, lwd)) = outline {
        if !col.is_transparent() && lwd > 0.0 {
            let (x0, x1) = (x0.min(x1), x0.max(x1));
            let (y0, y1) = (y0.min(y1), y0.max(y1));
            let s = stroke(ctx, col, lwd, dash, LineCap::Round);
            if x1 - x0 <= 0.0 || y1 - y0 <= 0.0 {
                c.fig.line(x0, y0, x1, y1, s);
            } else {
                c.fig.rect_outline(x0, y0, x1 - x0, y1 - y0, s);
            }
        }
    }
}

/// An element_rect grob filling its cell (the outline centred on the cell's edges).
pub fn rect_grob(ctx: Ctx, el: RectEl) -> Grob {
    Grob::draw_fn(move |c, b| {
        rect(c, &ctx, b.x0, b.y0, b.x1, b.y1, el.fill, el.colour.map(|col| (col, el.linewidth * PT)), Vec::new());
    })
}
