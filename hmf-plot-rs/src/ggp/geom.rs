//! Geoms as ggplot2 draws them in a panel, in device units: R's point symbols, rectangles (geom_tile, geom_rect,
//! geom_bar), text, segments and lines with R's dash patterns, ribbons, labels (geom_label), box plots.

use crate::color::Rgba;
use crate::display::{Item, LineCap, Stroke};
use crate::gtable::{Bx, Canvas};

use super::{Ctx, PT, STROKE, draw_text, rect, stroke};
use crate::rdev::Face;

/// A panel's mapping from data to device units (continuous ranges; discrete positions are 1..n).
#[derive(Clone, Copy, Debug)]
pub struct Panel {
    pub b: Bx,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    /// coord_flip: x data along the vertical axis
    pub flip: bool,
}

impl Panel {
    /// The device point of data (x, y).
    pub fn at(&self, x: f64, y: f64) -> (f64, f64) {
        if self.flip {
            (self.b.npc_x((y - self.yr.0) / (self.yr.1 - self.yr.0)), self.b.npc_y((x - self.xr.0) / (self.xr.1 - self.xr.0)))
        } else {
            (self.b.npc_x((x - self.xr.0) / (self.xr.1 - self.xr.0)), self.b.npc_y((y - self.yr.0) / (self.yr.1 - self.yr.0)))
        }
    }
}

/// R's point symbols (pch) as ggplot2 draws them: fontsize = size * .pt + stroke * .stroke / 2, the symbol's
/// radius 0.375 of it; lwd = stroke * .stroke / 2.
#[allow(clippy::too_many_arguments)]
pub fn point(c: &mut Canvas, ctx: &Ctx, x: f64, y: f64, shape: u8, colour: Rgba, fill: Rgba, size: f64, stroke_w: f64) {
    let fontsize = size * PT + stroke_w * STROKE / 2.0;
    let r = 0.375 * ctx.font_px(fontsize);
    let lwd = stroke_w * STROKE / 2.0;
    let lw = ctx.lwd(lwd);
    match shape {
        15 => {
            // filled square, no outline
            c.fig.rect(x - r, y - r, 2.0 * r, 2.0 * r, colour);
        }
        19 => {
            c.fig.push(Item::Circle { cx: x, cy: y, r, fill: Some(colour), stroke: Some(Stroke { color: colour, width: lw, dash: Vec::new(), cap: LineCap::Round }) });
        }
        21 => {
            c.fig.push(Item::Circle {
                cx: x,
                cy: y,
                r,
                fill: (!fill.is_transparent()).then_some(fill),
                stroke: (!colour.is_transparent()).then_some(Stroke { color: colour, width: lw, dash: Vec::new(), cap: LineCap::Round }),
            });
        }
        _ => {
            c.fig.push(Item::Circle { cx: x, cy: y, r, fill: Some(colour), stroke: None });
        }
    }
}

/// geom_tile / geom_rect: the rectangle between two data corners, filled and outlined (linewidth mm).
#[allow(clippy::too_many_arguments)]
pub fn data_rect(c: &mut Canvas, ctx: &Ctx, p: &Panel, x0: f64, x1: f64, y0: f64, y1: f64, fill: Option<Rgba>, outline: Option<(Rgba, f64)>) {
    let (ax, ay) = p.at(x0, y0);
    let (bx, by) = p.at(x1, y1);
    rect(c, ctx, ax, ay, bx, by, fill, outline.map(|(col, lw)| (col, lw * PT)), Vec::new());
}

/// geom_text: a label centred (by hjust, vjust) at a data point; size in mm (fontsize = size * .pt).
#[allow(clippy::too_many_arguments)]
pub fn text(c: &mut Canvas, ctx: &Ctx, p: &Panel, x: f64, y: f64, label: &str, size: f64, colour: Rgba, hjust: f64, vjust: f64) {
    let (dx, dy) = p.at(x, y);
    draw_text(c, ctx, dx, dy, label, size * PT, Face::Plain, colour, hjust, vjust, 0.0, 1.2);
}

/// A polyline of `linewidth` mm with a line type ("solid" or R's hex dash digits) and butt ends.
pub fn path(c: &mut Canvas, ctx: &Ctx, pts: Vec<(f64, f64)>, colour: Rgba, linewidth: f64, lty: &str) {
    let lwd = linewidth * PT;
    let dash = if lty == "solid" || lty.is_empty() { Vec::new() } else { super::dash_pattern_on(ctx, lty, lwd) };
    c.fig.polyline(pts, stroke(ctx, colour, lwd, dash, LineCap::Butt));
}

/// geom_hline: a line across the panel at y (data), as a segment from the panel's x range start to its end.
pub fn hline(c: &mut Canvas, ctx: &Ctx, p: &Panel, y: f64, colour: Rgba, linewidth: f64, lty: &str) {
    let a = p.at(p.xr.0, y);
    let b = p.at(p.xr.1, y);
    path(c, ctx, vec![a, b], colour, linewidth, lty);
}

/// geom_vline
pub fn vline(c: &mut Canvas, ctx: &Ctx, p: &Panel, x: f64, colour: Rgba, linewidth: f64, lty: &str) {
    let a = p.at(x, p.yr.0);
    let b = p.at(x, p.yr.1);
    path(c, ctx, vec![a, b], colour, linewidth, lty);
}

/// A filled polygon (no outline), e.g. geom_ribbon's area.
pub fn polygon(c: &mut Canvas, pts: Vec<(f64, f64)>, fill: Rgba) {
    if !fill.is_transparent() {
        c.fig.polygon(pts, Some(fill), None);
    }
}

/// The key of geom_label (draw_key_label, ggplot2 4.0.3): a label box around "a" (padding 0.25 lines, radius 0.15
/// lines, border `linewidth` mm) centred in the key, raised by font_descent(size) as measured when the legend was
/// built (on `mctx`; note that ggplot2 passes the size in mm as the font size there), its text drawn transparent.
#[allow(clippy::too_many_arguments)]
pub fn key_label(c: &mut Canvas, ctx: &Ctx, mctx: &Ctx, b: Bx, fill: Rgba, border: Rgba, size: f64, linewidth: f64) {
    let fs = size * PT;
    let dev = ctx.dev;
    let lines = dev.lines(1.0, fs, 1.2);
    let pad = 0.25 * lines;
    // labelGrob (ggplot2 4.0.3): top padding max(pad, descent), bottom max(pad - descent, 0), the box centred on its
    // point, which draw_key_label raises by font_descent(size) of the device the legend was built on (ggplot2 passes
    // the size in mm as the font size in points there)
    let desc = dev.font_descent(fs, Face::Plain);
    let w = dev.str_width("a", fs, Face::Plain) + 2.0 * pad;
    let h = dev.m_ascent(fs, Face::Plain) + desc + pad.max(desc) + (pad - desc).max(0.0);
    let raise = mctx.dev.font_descent(size, Face::Plain);
    let (cx, cy) = ((b.x0 + b.x1) / 2.0, (b.y0 + b.y1) / 2.0 - raise);
    let (x0, y0) = (cx - w / 2.0, cy - h / 2.0);
    let r = (0.15 * lines).min(w / 2.0).min(h / 2.0);
    super::drawn_box("rrvp", Bx::new(x0, y0, x0 + w, y0 + h));
    let mut pts = Vec::new();
    let corner = |pts: &mut Vec<(f64, f64)>, ccx: f64, ccy: f64, a0: f64| {
        for i in 0..=8 {
            let a = (a0 + 90.0 * i as f64 / 8.0).to_radians();
            pts.push((ccx + r * a.cos(), ccy + r * a.sin()));
        }
    };
    corner(&mut pts, x0 + w - r, y0 + r, -90.0);
    corner(&mut pts, x0 + w - r, y0 + h - r, 0.0);
    corner(&mut pts, x0 + r, y0 + h - r, 90.0);
    corner(&mut pts, x0 + r, y0 + r, 180.0);
    let st = (!border.is_transparent()).then(|| stroke(ctx, border, linewidth * PT, Vec::new(), LineCap::Round));
    c.fig.polygon(pts, (!fill.is_transparent()).then_some(fill), st);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ggp::{Version, record_drawn_boxes, take_drawn_boxes};
    use crate::rdev::Device;

    /// draw_key_label's box as R draws it (ggplot2 4.0.3, png() at 300 dpi, the legend built on pdf(NULL)): qsee's
    /// "QC status" key of size 10.5 * 0.35 in a 72 px key cell. R's viewports (the public HG008 report): the cell
    /// (4600.7401, 2495.0571)-(4672.7401, 2567.0571), the label box ("rrvp") (4610.1696, 2497.8199)-(4663.3105,
    /// 2556.9609): centred on the cell raised by one descent, not half a descent lower.
    #[test]
    fn key_label_box_as_r() {
        let ctx = Ctx { dev: Device::Png { dpi: 300.0 }, v: Version::V403 };
        let mctx = Ctx { dev: Device::PdfAt { dpi: 300.0 }, v: Version::V403 };
        let mut fig = crate::Figure::new(8, 8, Rgba::WHITE);
        let mut c = Canvas::new(&mut fig);
        record_drawn_boxes(true);
        key_label(&mut c, &ctx, &mctx, Bx::new(4600.7401, 2495.0571, 4672.7401, 2567.0571), Rgba::WHITE, Rgba::BLACK, 10.5 * 0.35, 10.5 / 22.0 * 0.5);
        let drawn = take_drawn_boxes();
        assert_eq!(drawn.len(), 1);
        let b = drawn[0].1;
        for (got, want) in [(b.x0, 4610.1696), (b.y0, 2497.8199), (b.x1, 4663.3105), (b.y1, 2556.9609)] {
            assert!((got - want).abs() < 0.001, "{got} vs R's {want}");
        }
        // dashes: Cairo scales by max(lwd, 1), pdf() by lwd (R's PDF of the qsee report: "0.64 w [ 0.64 0.64] 0 d")
        let lwd = 0.64 / 0.75;
        assert_eq!(crate::ggp::dash_pattern_on(&ctx, "11", lwd), vec![300.0 / 96.0; 2]);
        let pdf = Ctx { dev: Device::Pdf, v: Version::V403 };
        let d = crate::ggp::dash_pattern_on(&pdf, "11", lwd);
        assert!(d.iter().all(|&x| (x - 0.64).abs() < 1e-9), "{d:?}");
        // recording is off again
        key_label(&mut c, &ctx, &mctx, Bx::new(0.0, 0.0, 72.0, 72.0), Rgba::WHITE, Rgba::BLACK, 3.0, 0.25);
        assert!(take_drawn_boxes().is_empty());
    }
}
