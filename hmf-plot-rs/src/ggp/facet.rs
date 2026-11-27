//! Panel tables: facet_null's 3 x 3 table, and a facet grid (ggplot2 4.0.3's FacetGrid, and ggh4x 0.2.8's
//! facet_grid2, whose extra axis rows and columns between panels are empty in the plots here and so take no space):
//! panels with null sizes and the panel spacing between them, axes seamed on the outside, strips outside the axes
//! (strip.placement "outside", with strip.switch.pad.grid between) or next to the panels.

use crate::color::Rgba;
use crate::gtable::{Grob, Gtable, Unit};

use super::guides::Axis;
use super::{Ctx, RectEl, TextEl, TitleGrob, rect};

/// A facet strip: its background and label, and its size across the strip.
pub struct Strip {
    pub grob: Grob,
    /// height of a horizontal strip, width of a vertical one
    pub size: f64,
}

/// A strip (render_strips / ggh4x's draw_labels): the background element, then the label as a titleGrob with
/// margins on both sides. `vertical`: a strip beside the panels (its size is a width).
pub fn strip(ctx: Ctx, label: &str, text: Option<&TextEl>, bg: Option<RectEl>, vertical: bool) -> Strip {
    let title = text.map(|el| TitleGrob::new(label, el, true, true));
    let size = match &title {
        Some(t) => {
            if vertical {
                t.width(&ctx).abs()
            } else {
                t.height(&ctx).abs()
            }
        }
        None => 0.0,
    };
    if title.is_none() && bg.is_none() {
        return Strip { grob: Grob::Zero, size: 0.0 };
    }
    let grob = Grob::draw_fn(move |c, b| {
        if let Some(bg) = bg {
            rect(c, &ctx, b.x0, b.y0, b.x1, b.y1, bg.fill, bg.colour.map(|col| (col, bg.linewidth * super::PT)), Vec::new());
        }
        if let Some(t) = &title {
            t.draw(c, &ctx, b);
        }
    });
    Strip { grob, size }
}

/// A transparent colour ("#00000000"): ggh4x's strip_themed hides strips with it.
pub const INVISIBLE: Rgba = Rgba::TRANSPARENT;

/// facet_null's table: spacer, axis-t, spacer / axis-l, panel, axis-r / spacer, axis-b, spacer, with ggplot2's z
/// values, the panel clipped.
pub fn facet_null(panel: Grob, top: Axis, right: Axis, bottom: Axis, left: Axis, clip: bool) -> Gtable {
    let mut gt = Gtable::new(vec![Unit::Abs(left.size), Unit::Null(1.0), Unit::Abs(right.size)], vec![Unit::Abs(top.size), Unit::Null(1.0), Unit::Abs(bottom.size)]);
    // gtable_matrix with z = (5 6 4 / 7 1 8 / 3 9 2), column by column
    let cells: [(usize, usize, Grob, f64, &str); 9] = [
        (1, 1, Grob::Zero, 5.0, "spacer"),
        (2, 1, left.grob, 7.0, "axis-l"),
        (3, 1, Grob::Zero, 3.0, "spacer"),
        (1, 2, top.grob, 6.0, "axis-t"),
        (2, 2, panel, 1.0, "panel"),
        (3, 2, bottom.grob, 9.0, "axis-b"),
        (1, 3, Grob::Zero, 4.0, "spacer"),
        (2, 3, right.grob, 8.0, "axis-r"),
        (3, 3, Grob::Zero, 2.0, "spacer"),
    ];
    for (r, c, g, z, name) in cells {
        let clip_this = name == "panel" && clip;
        gt.add(g, r, c, r, c, Some(z), clip_this, name);
    }
    gt
}

/// A facet grid.
pub struct Grid {
    /// null widths of the panel columns (space "free": the x ranges) and heights of the rows
    pub widths: Vec<f64>,
    pub heights: Vec<f64>,
    /// panels[row][col]
    pub panels: Vec<Vec<Grob>>,
    pub spacing_x: f64,
    pub spacing_y: f64,
    pub clip: bool,
    /// axes outside the grid: top and bottom per column, left and right per row
    pub axes_top: Vec<Axis>,
    pub axes_bottom: Vec<Axis>,
    pub axes_left: Vec<Axis>,
    pub axes_right: Vec<Axis>,
    /// strips: on top (or bottom, `switch_x`) per column, on the right (or left, `switch_y`) per row
    pub strips_x: Option<Vec<Strip>>,
    pub strips_y: Option<Vec<Strip>>,
    pub switch_x: bool,
    pub switch_y: bool,
    pub outside_x: bool,
    pub outside_y: bool,
    /// strip.switch.pad.grid
    pub padding: f64,
    /// ggh4x 0.2.8 adds the padding for outside strips whether or not there is an axis on that side; ggplot2
    /// 4.0.3 only when there is one
    pub pad_always: bool,
}

/// The rows and columns of the panels in a table (panel_rows / panel_cols).
pub fn panel_extent(gt: &Gtable) -> (usize, usize, usize, usize) {
    let ps: Vec<_> = gt.items.iter().filter(|p| p.name.starts_with("panel")).collect();
    let t = ps.iter().map(|p| p.t).min().unwrap_or(1);
    let b = ps.iter().map(|p| p.b).max().unwrap_or(1);
    let l = ps.iter().map(|p| p.l).min().unwrap_or(1);
    let r = ps.iter().map(|p| p.r).max().unwrap_or(1);
    (t, l, b, r)
}

fn distinct_rows(gt: &Gtable) -> Vec<usize> {
    let mut v: Vec<usize> = gt.items.iter().filter(|p| p.name.starts_with("panel")).map(|p| p.t).collect();
    v.sort();
    v.dedup();
    v
}

fn distinct_cols(gt: &Gtable) -> Vec<usize> {
    let mut v: Vec<usize> = gt.items.iter().filter(|p| p.name.starts_with("panel")).map(|p| p.l).collect();
    v.sort();
    v.dedup();
    v
}

fn max_size<'a>(xs: impl Iterator<Item = &'a f64>) -> f64 {
    xs.copied().fold(0.0, f64::max)
}

/// seam_table: a row or column of grobs added on one side of the panels, `shift` cells out, after an optional
/// spacing row or column.
fn seam(gt: &mut Gtable, grobs: Vec<Grob>, sizes: &[f64], side: char, shift: usize, name: &str, z: f64, clip: bool, spacing: Option<f64>) {
    let (t, l, b, r) = panel_extent(gt);
    let mut row = match side {
        'b' => b + shift - 1,
        't' => t - shift,
        _ => 0,
    };
    let mut col = match side {
        'r' => r + shift - 1,
        'l' => l - shift,
        _ => 0,
    };
    if let Some(sp) = spacing {
        match side {
            'b' | 't' => gt.add_rows(&[Unit::Abs(sp)], row),
            _ => gt.add_cols(&[Unit::Abs(sp)], col),
        }
        if side == 'b' {
            row += 1;
        }
        if side == 'r' {
            col += 1;
        }
    }
    let size = max_size(sizes.iter());
    match side {
        'b' | 't' => gt.add_rows(&[Unit::Abs(size)], row),
        _ => gt.add_cols(&[Unit::Abs(size)], col),
    }
    if side == 't' || side == 'b' {
        row += 1;
    } else {
        col += 1;
    }
    // the panels' rows and columns after the insertion
    let rows = distinct_rows(gt);
    let cols = distinct_cols(gt);
    for (i, g) in grobs.into_iter().enumerate() {
        let (tt, ll) = match side {
            't' | 'b' => (row, cols[i.min(cols.len() - 1)]),
            _ => (rows[i.min(rows.len() - 1)], col),
        };
        gt.add(g, tt, ll, tt, ll, Some(z), clip, &format!("{name}-{}", i + 1));
    }
}

impl Grid {
    pub fn build(self) -> Gtable {
        let nrow = self.heights.len();
        let ncol = self.widths.len();
        let mut gt = Gtable::new(self.widths.iter().map(|&w| Unit::Null(w)).collect(), self.heights.iter().map(|&h| Unit::Null(h)).collect());
        for (r, row) in self.panels.into_iter().enumerate() {
            for (c, p) in row.into_iter().enumerate() {
                gt.add(p, r + 1, c + 1, r + 1, c + 1, Some(1.0), self.clip, &format!("panel-{}-{}", r + 1, c + 1));
            }
        }
        let _ = (nrow, ncol);
        gt.add_col_space(self.spacing_x);
        gt.add_row_space(self.spacing_y);
        // axes
        let sizes = |a: &Vec<Axis>| a.iter().map(|x| x.size).collect::<Vec<f64>>();
        let (st, sb, sl, sr) = (sizes(&self.axes_top), sizes(&self.axes_bottom), sizes(&self.axes_left), sizes(&self.axes_right));
        let has = |a: &Vec<Axis>| a.iter().any(|x| !x.grob.is_zero());
        let (ht, hb, hl, hr) = (has(&self.axes_top), has(&self.axes_bottom), has(&self.axes_left), has(&self.axes_right));
        seam(&mut gt, self.axes_top.into_iter().map(|a| a.grob).collect(), &st, 't', 1, "axis-t", 3.0, false, None);
        seam(&mut gt, self.axes_bottom.into_iter().map(|a| a.grob).collect(), &sb, 'b', 1, "axis-b", 3.0, false, None);
        seam(&mut gt, self.axes_left.into_iter().map(|a| a.grob).collect(), &sl, 'l', 1, "axis-l", 3.0, false, None);
        seam(&mut gt, self.axes_right.into_iter().map(|a| a.grob).collect(), &sr, 'r', 1, "axis-r", 3.0, false, None);
        // strips
        if let Some(strips) = self.strips_x {
            let (side, has_axis) = if self.switch_x { ('b', hb) } else { ('t', ht) };
            let shift = if self.outside_x { 2 } else { 1 };
            let space = (self.outside_x && (self.pad_always || has_axis)).then_some(self.padding);
            let sizes: Vec<f64> = strips.iter().map(|s| s.size).collect();
            let name = if self.switch_x { "strip-b" } else { "strip-t" };
            seam(&mut gt, strips.into_iter().map(|s| s.grob).collect(), &sizes, side, shift, name, 2.0, true, space);
        }
        if let Some(strips) = self.strips_y {
            let (side, has_axis) = if self.switch_y { ('l', hl) } else { ('r', hr) };
            let shift = if self.outside_y { 2 } else { 1 };
            let space = (self.outside_y && (self.pad_always || has_axis)).then_some(self.padding);
            let sizes: Vec<f64> = strips.iter().map(|s| s.size).collect();
            let name = if self.switch_y { "strip-l" } else { "strip-r" };
            seam(&mut gt, strips.into_iter().map(|s| s.grob).collect(), &sizes, side, shift, name, 2.0, true, space);
        }
        gt
    }
}

/// Unused text element helper kept for strips drawn by callers.
pub fn strip_text(el: &TextEl) -> TextEl {
    el.clone()
}
