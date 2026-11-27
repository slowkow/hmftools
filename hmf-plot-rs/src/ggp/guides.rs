//! Guides: position axes (GuideAxis), colour bars (GuideColourbar) and key legends (GuideLegend), and the guide box
//! that holds a plot's legends (Guides$package_box), with ggplot2's measurements.

use crate::color::Rgba;
use crate::display::LineCap;
use crate::gtable::{Bx, Canvas, Grob, Gtable, Pos, Unit, Vp, sum_abs};

use super::{Ctx, LineEl, PT, RectEl, TextEl, TitleGrob, Version, rect, stroke};

/// A side of a panel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Top,
    Right,
    Bottom,
    Left,
}

impl Side {
    pub fn vertical(&self) -> bool {
        matches!(self, Side::Left | Side::Right)
    }
}

/// An axis: breaks (npc along the axis) with labels, the label and tick elements, the tick length (device units).
#[derive(Clone, Debug)]
pub struct AxisSpec {
    pub side: Side,
    pub breaks: Vec<f64>,
    pub labels: Vec<String>,
    /// None: element_blank
    pub text: Option<TextEl>,
    pub ticks: Option<LineEl>,
    pub tick_len: f64,
    /// ticks at the breaks (guide_axis_logticks draws only its own)
    pub major_ticks: bool,
    /// more ticks (npc) with their lengths (guide_axis_logticks: long, mid and short)
    pub extra_ticks: Vec<(f64, f64)>,
}

/// An axis grob and its size across the axis (height of a horizontal axis, width of a vertical one).
pub struct Axis {
    pub grob: Grob,
    pub size: f64,
}

impl Axis {
    pub fn zero() -> Axis {
        Axis { grob: Grob::Zero, size: 0.0 }
    }
}

/// GuideAxis: [ticks, spacer, labels] outwards from the panel; the tick length counts even when the ticks are blank
/// (measure_grobs takes the theme's length), the labels only when they are drawn.
pub fn axis(ctx: Ctx, spec: AxisSpec) -> Axis {
    let vertical = spec.side.vertical();
    let tick_len = spec.extra_ticks.iter().map(|t| t.1).fold(spec.tick_len, f64::max).max(0.0);
    let labels = match &spec.text {
        Some(el) if !spec.labels.is_empty() => {
            let mut t = TitleGrob { labels: spec.labels.clone(), x: None, y: None, el: el.clone(), margin_x: vertical, margin_y: !vertical };
            if vertical {
                t.y = Some(spec.breaks.clone());
            } else {
                t.x = Some(spec.breaks.clone());
            }
            Some(t)
        }
        _ => None,
    };
    let label_size = labels.as_ref().map(|t| if vertical { t.width(&ctx).abs() } else { t.height(&ctx).abs() }).unwrap_or(0.0);
    let size = tick_len + label_size;
    let side = spec.side;
    let ticks = spec.ticks.clone();
    let breaks = spec.breaks.clone();
    let main_len = spec.tick_len;
    let major = spec.major_ticks;
    let extra = spec.extra_ticks.clone();
    let grob = Grob::draw_fn(move |c, b| {
        // the axis sits on the panel side of its cell
        let (tick_box, label_box) = match side {
            Side::Top => (Bx::new(b.x0, b.y1 - tick_len, b.x1, b.y1), Bx::new(b.x0, b.y1 - size, b.x1, b.y1 - tick_len)),
            Side::Bottom => (Bx::new(b.x0, b.y0, b.x1, b.y0 + tick_len), Bx::new(b.x0, b.y0 + tick_len, b.x1, b.y0 + size)),
            Side::Left => (Bx::new(b.x1 - tick_len, b.y0, b.x1, b.y1), Bx::new(b.x1 - size, b.y0, b.x1 - tick_len, b.y1)),
            Side::Right => (Bx::new(b.x0, b.y0, b.x0 + tick_len, b.y1), Bx::new(b.x0 + tick_len, b.y0, b.x0 + size, b.y1)),
        };
        if let Some(t) = &ticks {
            let s = stroke(&ctx, t.colour, t.linewidth * PT, t.dash.clone(), LineCap::Butt);
            let mut draw_tick = |p: f64, len: f64| match side {
                Side::Top => {
                    let x = tick_box.npc_x(p);
                    c.fig.line(x, tick_box.y1, x, tick_box.y1 - len, s.clone());
                }
                Side::Bottom => {
                    let x = tick_box.npc_x(p);
                    c.fig.line(x, tick_box.y0, x, tick_box.y0 + len, s.clone());
                }
                Side::Left => {
                    let y = tick_box.npc_y(p);
                    c.fig.line(tick_box.x1, y, tick_box.x1 - len, y, s.clone());
                }
                Side::Right => {
                    let y = tick_box.npc_y(p);
                    c.fig.line(tick_box.x0, y, tick_box.x0 + len, y, s.clone());
                }
            };
            if major {
                for &p in &breaks {
                    draw_tick(p, main_len);
                }
            }
            for &(p, len) in &extra {
                draw_tick(p, len);
            }
        }
        if let Some(t) = &labels {
            t.draw(c, &ctx, label_box);
        }
    });
    Axis { grob, size }
}

// -------------------------------------------------------------------------------------------------------------
// legends

/// A key of a key legend.
#[derive(Clone, Debug)]
pub enum KeyGlyph {
    /// draw_key_polygon: a rectangle inset by the outline width (mm), filled and outlined
    Polygon { fill: Rgba, colour: Option<Rgba>, linewidth: f64 },
    /// draw_key_point with a shape: 15 (square), 19 (disc)
    Point { shape: u8, colour: Rgba, fill: Rgba, size: f64, stroke: f64 },
    /// draw_key_label: a label box (geom_label's key); `linewidth` the border's (mm): the layer's linewidth, which in
    /// ggplot2 4.0.3 defaults to the theme's borderwidth / 2 (draw_key_label's own 0.25 applies only without one)
    Label { fill: Rgba, border: Rgba, size: f64, linewidth: f64 },
}

/// A legend to build.
#[derive(Clone, Debug)]
pub enum LegendSpec {
    /// guide_colourbar, vertical: the colours of the bar's bins from the low limit up (R's extract_decor: the
    /// scale's colours of nbin values from limit to limit), the breaks as fractions of the bar (key$.value) and
    /// their labels
    Bar { title: String, colours: Vec<Rgba>, breaks: Vec<f64>, labels: Vec<String>, bar_width: f64, bar_height: f64, frame: Option<(Rgba, f64)>, ticks: Option<(Rgba, f64)>, tick_len: f64 },
    /// guide_legend: keys with labels in one column (vertical) or as the caller arranges
    Keys { title: String, entries: Vec<(String, KeyGlyph)>, key_width: f64, key_height: f64 },
}

/// The theme parts a legend uses.
#[derive(Clone, Debug)]
pub struct LegendTheme {
    pub title: TextEl,
    pub text: TextEl,
    /// legend.key.spacing (the title gap and the key-label gap), device units
    pub key_spacing: f64,
    /// legend.margin (top, right, bottom, left)
    pub margin: [f64; 4],
    pub background: Option<RectEl>,
    /// legend.key: drawn in each key's cell before the key's glyph (R's "key-<i>-1-bg"); white without an outline in
    /// the themes CUPPA and qsee use. Invisible on a white legend, but it covers the edge of the key above it, as R's
    /// does (a hairline gap between outlined keys).
    pub key_background: Option<RectEl>,
}

/// A legend's gtable (GuideLegend / GuideColourbar assemble_drawing), title on top, text on the right.
pub fn legend(ctx: Ctx, spec: &LegendSpec, th: &LegendTheme) -> Gtable {
    legend_measured(ctx, ctx, spec, th)
}

/// A legend sized with the text metrics of `mctx` (the device that was open when R built it) and drawn on `ctx`.
pub fn legend_measured(ctx: Ctx, mctx: Ctx, spec: &LegendSpec, th: &LegendTheme) -> Gtable {
    let gap = th.key_spacing;
    let mut title_el = th.title.clone();
    // position_margin("top", margin, gap): the title's bottom margin grows by the gap
    title_el.margin[2] += gap;
    let title = TitleGrob::new(spec_title(spec), &title_el, true, true);
    let mut text_el = th.text.clone();
    text_el.margin[3] += gap;
    let mut gt = match spec {
        LegendSpec::Bar { colours, breaks, labels, bar_width, bar_height, frame, ticks, tick_len, .. } => {
            let values = breaks.clone();
            let lab = TitleGrob { labels: labels.clone(), x: None, y: Some(values.clone()), el: text_el.clone(), margin_x: true, margin_y: false };
            let label_w = lab.width(&mctx).abs();
            let mut gt = Gtable::new(vec![Unit::Abs(*bar_width), Unit::Abs(label_w)], vec![Unit::Abs(*bar_height)]);
            let (colours, frame, ticks, tick_len) = (colours.clone(), *frame, *ticks, *tick_len);
            let vals = values.clone();
            let bar = Grob::draw_fn(move |c, b| {
                draw_bar(c, &ctx, b, &colours);
                if let Some((col, lw)) = frame {
                    rect(c, &ctx, b.x0, b.y0, b.x1, b.y1, None, Some((col, lw * PT)), Vec::new());
                }
                if let Some((col, lw)) = ticks {
                    let s = stroke(&ctx, col, lw * PT, Vec::new(), LineCap::Butt);
                    for &v in &vals {
                        let y = b.npc_y(v);
                        c.fig.line(b.x1, y, b.x1 - tick_len, y, s.clone());
                        c.fig.line(b.x0, y, b.x0 + tick_len, y, s.clone());
                    }
                }
            });
            gt.add(bar, 1, 1, 1, 1, None, false, "bar");
            gt.add(lab.grob(ctx), 1, 2, 1, 2, None, false, "labels");
            gt
        }
        LegendSpec::Keys { entries, key_width, key_height, .. } => {
            let labs: Vec<TitleGrob> = entries.iter().map(|(l, _)| TitleGrob::new(l, &text_el, true, true)).collect();
            let label_w = labs.iter().map(|t| t.width(&mctx).abs()).fold(0.0, f64::max);
            let heights: Vec<Unit> = labs.iter().map(|t| Unit::Abs(t.height(&mctx).abs().max(*key_height))).collect();
            let mut gt = Gtable::new(vec![Unit::Abs(*key_width), Unit::Abs(label_w)], heights);
            for (i, ((_, key), lab)) in entries.iter().zip(labs).enumerate() {
                let key = key.clone();
                let kbg = th.key_background;
                gt.add(
                    Grob::draw_fn(move |c, b| {
                        if let Some(bg) = kbg {
                            rect(c, &ctx, b.x0, b.y0, b.x1, b.y1, bg.fill, bg.colour.map(|col| (col, bg.linewidth * PT)), Vec::new());
                        }
                        draw_key_m(c, &ctx, &mctx, b, &key)
                    }),
                    i + 1,
                    1,
                    i + 1,
                    1,
                    None,
                    false,
                    "key",
                );
                gt.add(lab.grob(ctx), i + 1, 2, i + 1, 2, None, false, "label");
            }
            gt
        }
    };
    // add_title: a row on top; a title wider than the table widens it on the right (hjust 0)
    let tw = title.width(&mctx).abs();
    let th_ = title.height(&mctx).abs();
    gt.add_rows(&[Unit::Abs(th_)], 0);
    let nc = gt.ncol();
    gt.add(title.grob(ctx), 1, 1, 1, nc, Some(f64::NEG_INFINITY), false, "title");
    let table_w = sum_abs(&gt.widths);
    let extra = (tw - table_w).max(0.0);
    if extra > 0.0 {
        gt.add_cols(&[Unit::Abs(0.0)], 0);
        gt.add_cols(&[Unit::Abs(extra)], usize::MAX);
    }
    let [t, r, b, l] = th.margin;
    gt.add_padding(t, r, b, l);
    if let Some(bg) = th.background {
        let (nr, nc) = (gt.nrow(), gt.ncol());
        gt.add(super::rect_grob(ctx, bg), 1, 1, nr, nc, Some(f64::NEG_INFINITY), false, "background");
    }
    gt
}

fn spec_title(spec: &LegendSpec) -> &str {
    match spec {
        LegendSpec::Bar { title, .. } | LegendSpec::Keys { title, .. } => title,
    }
}

/// The colour bar: R's raster of the bins' colours (bottom to top) stretched over the bar with interpolation. A PNG
/// row takes the colour at its centre, linearly between the two nearest bins (Cairo's bilinear filter); a PDF gets a
/// rectangle per bin.
fn draw_bar(c: &mut Canvas, ctx: &Ctx, b: Bx, colours: &[Rgba]) {
    let n = colours.len();
    if n == 0 {
        return;
    }
    match ctx.dev {
        crate::rdev::Device::Png { .. } | crate::rdev::Device::PdfAt { .. } => {
            let (r0, r1) = (b.y0.round() as i64, b.y1.round() as i64);
            for row in r0..r1 {
                let f = 1.0 - ((row as f64 + 0.5) - b.y0) / b.h();
                let pos = (f * n as f64 - 0.5).clamp(0.0, (n - 1) as f64);
                let i0 = pos.floor() as usize;
                let i1 = (i0 + 1).min(n - 1);
                let t = (pos - i0 as f64) as f32;
                let (a, z) = (colours[i0], colours[i1]);
                let mix = |p: f32, q: f32| ((p * 255.0) * (1.0 - t) + (q * 255.0) * t).round() / 255.0;
                let col = Rgba { r: mix(a.r, z.r), g: mix(a.g, z.g), b: mix(a.b, z.b), a: 1.0 };
                c.fig.rect(b.x0, row as f64, b.w(), 1.0, col);
            }
        }
        crate::rdev::Device::Pdf => {
            let h = b.h() / n as f64;
            for (k, col) in colours.iter().enumerate() {
                let y = b.y1 - (k + 1) as f64 * h;
                c.fig.rect(b.x0, y, b.w(), h + if k + 1 < n { h * 0.05 } else { 0.0 }, *col);
            }
        }
    }
}

/// Draws a legend key in its cell.
pub fn draw_key(c: &mut Canvas, ctx: &Ctx, b: Bx, key: &KeyGlyph) {
    draw_key_m(c, ctx, ctx, b, key)
}

/// Draws a legend key whose legend was built on another device (`mctx`).
pub fn draw_key_m(c: &mut Canvas, ctx: &Ctx, mctx: &Ctx, b: Bx, key: &KeyGlyph) {
    match key {
        KeyGlyph::Polygon { fill, colour, linewidth } => {
            // draw_key_polygon: lwd = min(linewidth, min(size) / 4) mm; the rectangle is the cell less lwd mm
            let size_mm = b.w().min(b.h()) / ctx.dev.mm(1.0);
            let lwd = linewidth.min(size_mm / 4.0);
            let inset = ctx.mm(lwd) / 2.0;
            let (x0, y0, x1, y1) = (b.x0 + inset, b.y0 + inset, b.x1 - inset, b.y1 - inset);
            rect(c, ctx, x0, y0, x1, y1, Some(*fill), colour.map(|col| (col, lwd * PT)), Vec::new());
        }
        KeyGlyph::Point { shape, colour, fill, size, stroke: st } => {
            super::geom::point(c, ctx, (b.x0 + b.x1) / 2.0, (b.y0 + b.y1) / 2.0, *shape, *colour, *fill, *size, *st);
        }
        KeyGlyph::Label { fill, border, size, linewidth } => {
            super::geom::key_label(c, ctx, mctx, b, *fill, *border, *size, *linewidth);
        }
    }
}

/// Guides$package_box for one or more legends at a side: a column (vertical box) or row (horizontal box) of
/// legends justified within the box, the box margin, and null margins that justify the box in the guide-box cell.
pub struct GuideBox {
    pub legends: Vec<Gtable>,
    /// legend.box: true for horizontal (legends side by side)
    pub horizontal: bool,
    /// legend.justification of the box in its cell (x, y)
    pub just: (f64, f64),
    /// legend.box.just of each legend in the box (x, y)
    pub box_just: (f64, f64),
    /// legend.spacing.x / y between legends
    pub spacing: f64,
    pub box_margin: [f64; 4],
}

pub fn guide_box(ctx: Ctx, gb: GuideBox) -> Gtable {
    let _ = ctx;
    let widths: Vec<f64> = gb.legends.iter().map(|g| sum_abs(&g.widths)).collect();
    let heights: Vec<f64> = gb.legends.iter().map(|g| sum_abs(&g.heights)).collect();
    let (bx, by) = gb.box_just;
    let mut guides;
    if gb.horizontal {
        let hmax = heights.iter().copied().fold(0.0, f64::max);
        let vp_w = widths.iter().sum::<f64>() + gb.spacing * (gb.legends.len().saturating_sub(1)) as f64;
        guides = Gtable::new(widths.iter().map(|&w| Unit::Abs(w)).collect(), vec![Unit::Abs(hmax)]);
        for (i, (g, h)) in gb.legends.into_iter().zip(&heights).enumerate() {
            let vp = Vp { x: Pos::npc(bx), y: Pos::npc(by), w: None, h: Some(*h), hjust: bx, vjust: by };
            guides.add_vp(Grob::table(g), vp, 1, i + 1, 1, i + 1, None, false, "guides");
        }
        guides.add_col_space(gb.spacing);
        guides.vp = Some(Vp { x: Pos::npc(gb.just.0), y: Pos::npc(gb.just.1), w: Some(vp_w), h: Some(hmax), hjust: gb.just.0, vjust: gb.just.1 });
    } else {
        let wmax = widths.iter().copied().fold(0.0, f64::max);
        let vp_h = heights.iter().sum::<f64>() + gb.spacing * (gb.legends.len().saturating_sub(1)) as f64;
        guides = Gtable::new(vec![Unit::Abs(wmax)], heights.iter().map(|&h| Unit::Abs(h)).collect());
        for (i, (g, w)) in gb.legends.into_iter().zip(&widths).enumerate() {
            let vp = Vp { x: Pos::npc(bx), y: Pos::npc(by), w: Some(*w), h: None, hjust: bx, vjust: by };
            guides.add_vp(Grob::table(g), vp, i + 1, 1, i + 1, 1, None, false, "guides");
        }
        guides.add_row_space(gb.spacing);
        guides.vp = Some(Vp { x: Pos::npc(gb.just.0), y: Pos::npc(gb.just.1), w: Some(wmax), h: Some(vp_h), hjust: gb.just.0, vjust: gb.just.1 });
    }
    // the box margin and the global margin are added to the same table, inside its justification viewport (as R
    // does: with a zero box margin the legends fill the viewport)
    let [t, r, b, l] = gb.box_margin;
    guides.add_padding(t, r, b, l);
    // global margin: null units that justify the box (t = 1 - yjust, b = yjust, r = 1 - xjust, l = xjust)
    let (xj, yj) = gb.just;
    guides.add_rows(&[Unit::Null(1.0 - yj)], 0);
    guides.add_cols(&[Unit::Null(1.0 - xj)], usize::MAX);
    guides.add_rows(&[Unit::Null(yj)], usize::MAX);
    guides.add_cols(&[Unit::Null(xj)], 0);
    guides
}

/// The version's legend.key.spacing default (5.5 pt in both for base size 11; 4.0.3 takes half_line)
pub fn key_spacing(ctx: &Ctx, half_line: f64) -> f64 {
    match ctx.v {
        Version::V351 => ctx.pt(5.5),
        Version::V403 => ctx.pt(half_line),
    }
}
