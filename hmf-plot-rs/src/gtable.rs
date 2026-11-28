//! gtable (0.3.6) and grid's layout, as ggplot2 and patchwork build figures with them: a table of column widths
//! and row heights (absolute lengths, or "null" units sharing what is left), grobs placed on cell ranges with a z
//! order, nested tables, viewports placing a grob inside its cell, and clipping. Lengths are in device units (pixels
//! for a PNG, points for a PDF); a [`Canvas`] draws the table into a display list.
//!
//! Only the operations the hmftools plots need are here: `gtable_add_rows/cols/grob/padding/row_space/col_space`,
//! subsetting (`gt[rows, cols]`), and drawing with grid's layout rules (absolute sizes first, null units sharing the
//! remainder in proportion, none when nothing is left).

use std::rc::Rc;

use crate::display::{Figure, Item};

/// A width or height: an absolute length in device units, or a share of the space left (grid's "null" unit).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Unit {
    Abs(f64),
    Null(f64),
}

impl Unit {
    pub const ZERO: Unit = Unit::Abs(0.0);

    /// The absolute length (a null unit converts to 0, as grid's convertUnit does outside a layout).
    pub fn abs(&self) -> f64 {
        match self {
            Unit::Abs(v) => *v,
            Unit::Null(_) => 0.0,
        }
    }

    pub fn is_null(&self) -> bool {
        matches!(self, Unit::Null(_))
    }

    /// patchwork's `as.numeric(unit) == 0`
    pub fn value(&self) -> f64 {
        match self {
            Unit::Abs(v) | Unit::Null(v) => *v,
        }
    }
}

/// The sum of the absolute parts of some units (gtable_width / gtable_height outside a layout).
pub fn sum_abs(units: &[Unit]) -> f64 {
    units.iter().map(|u| u.abs()).sum()
}

/// A box in device units, y downwards.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bx {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}

impl Bx {
    pub fn new(x0: f64, y0: f64, x1: f64, y1: f64) -> Bx {
        Bx { x0, y0, x1, y1 }
    }
    pub fn w(&self) -> f64 {
        self.x1 - self.x0
    }
    pub fn h(&self) -> f64 {
        self.y1 - self.y0
    }
    /// The x at a fraction of the width from the left (npc).
    pub fn npc_x(&self, f: f64) -> f64 {
        self.x0 + f * self.w()
    }
    /// The y at a fraction of the height from the bottom (grid's npc: 0 is the bottom).
    pub fn npc_y(&self, f: f64) -> f64 {
        self.y1 - f * self.h()
    }
}

/// Draws into a display list with a clip stack (a grob drawn with clip "on" restricts what it draws to its cell).
pub struct Canvas<'a> {
    pub fig: &'a mut Figure,
    clips: Vec<Option<(f64, f64, f64, f64)>>,
}

impl<'a> Canvas<'a> {
    pub fn new(fig: &'a mut Figure) -> Canvas<'a> {
        Canvas { fig, clips: vec![None] }
    }

    fn current_clip(&self) -> Option<(f64, f64, f64, f64)> {
        *self.clips.last().unwrap()
    }

    /// Restricts what follows to `b` (intersected with the current clip) until `pop_clip`.
    pub fn push_clip(&mut self, b: Bx) {
        let c = match self.current_clip() {
            None => (b.x0, b.y0, b.w(), b.h()),
            Some((x, y, w, h)) => {
                let x0 = x.max(b.x0);
                let y0 = y.max(b.y0);
                let x1 = (x + w).min(b.x1);
                let y1 = (y + h).min(b.y1);
                (x0, y0, (x1 - x0).max(0.0), (y1 - y0).max(0.0))
            }
        };
        self.clips.push(Some(c));
        self.fig.push(Item::Clip(Some(c)));
    }

    pub fn pop_clip(&mut self) {
        self.clips.pop();
        let c = self.current_clip();
        self.fig.push(Item::Clip(c));
    }
}

pub type DrawFn = dyn Fn(&mut Canvas, Bx);

/// What a table cell holds.
#[derive(Clone)]
pub enum Grob {
    Zero,
    Table(Rc<Gtable>),
    /// Draws itself in the box it is given.
    Draw(Rc<DrawFn>),
}

impl Grob {
    pub fn draw_fn(f: impl Fn(&mut Canvas, Bx) + 'static) -> Grob {
        Grob::Draw(Rc::new(f))
    }

    pub fn table(t: Gtable) -> Grob {
        Grob::Table(Rc::new(t))
    }

    pub fn is_zero(&self) -> bool {
        matches!(self, Grob::Zero)
    }

    pub fn draw(&self, c: &mut Canvas, b: Bx) {
        match self {
            Grob::Zero => {}
            Grob::Table(t) => t.draw(c, b),
            Grob::Draw(f) => f(c, b),
        }
    }
}

/// A position inside a cell: a fraction of its size (npc, 0 at the left / bottom) plus an absolute offset.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pos {
    pub npc: f64,
    pub abs: f64,
}

impl Pos {
    pub fn npc(npc: f64) -> Pos {
        Pos { npc, abs: 0.0 }
    }
}

/// A viewport placing a grob in its cell (grid's viewport(x, y, width, height, just)); None sizes are the cell's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vp {
    pub x: Pos,
    pub y: Pos,
    pub w: Option<f64>,
    pub h: Option<f64>,
    pub hjust: f64,
    pub vjust: f64,
}

impl Vp {
    pub fn full() -> Vp {
        Vp { x: Pos::npc(0.5), y: Pos::npc(0.5), w: None, h: None, hjust: 0.5, vjust: 0.5 }
    }

    /// The viewport's box inside `cell`.
    pub fn apply(&self, cell: Bx) -> Bx {
        let w = self.w.unwrap_or(cell.w());
        let h = self.h.unwrap_or(cell.h());
        let x = cell.x0 + self.x.npc * cell.w() + self.x.abs - self.hjust * w;
        // grid's y goes up: the reference point is npc from the bottom; vjust 0 puts the bottom there
        let y_ref = cell.y1 - (self.y.npc * cell.h() + self.y.abs);
        let y_top = y_ref - (1.0 - self.vjust) * h;
        Bx::new(x, y_top, x + w, y_top + h)
    }
}

/// A grob placed on the cells t..=b, l..=r (1-based, as gtable's layout).
#[derive(Clone)]
pub struct Placed {
    pub t: usize,
    pub l: usize,
    pub b: usize,
    pub r: usize,
    pub z: f64,
    pub clip: bool,
    pub name: String,
    pub grob: Grob,
    pub vp: Option<Vp>,
    /// What the grob is, for comparing grobs as patchwork's `hash(unname_grob(x))` does (collect_axes,
    /// collect_axis_titles): grobs with equal keys are the same; a grob without a key equals no other.
    pub key: Option<String>,
}

#[derive(Clone, Default)]
pub struct Gtable {
    pub widths: Vec<Unit>,
    pub heights: Vec<Unit>,
    pub items: Vec<Placed>,
    /// a viewport the whole table is drawn in, inside the box it is given (grob$vp)
    pub vp: Option<Vp>,
}

impl Gtable {
    pub fn new(widths: Vec<Unit>, heights: Vec<Unit>) -> Gtable {
        Gtable { widths, heights, items: Vec::new(), vp: None }
    }

    pub fn ncol(&self) -> usize {
        self.widths.len()
    }

    pub fn nrow(&self) -> usize {
        self.heights.len()
    }

    /// gtable_add_grob for one grob (t, l, b, r 1-based; z: None is Inf, Some(f64::NEG_INFINITY) is -Inf).
    #[allow(clippy::too_many_arguments)]
    pub fn add(&mut self, grob: Grob, t: usize, l: usize, b: usize, r: usize, z: Option<f64>, clip: bool, name: &str) {
        let z = self.resolve_z(z);
        self.items.push(Placed { t, l, b, r, z, clip, name: name.to_string(), grob, vp: None, key: None });
    }

    /// gtable_add_grob with a viewport for the grob.
    #[allow(clippy::too_many_arguments)]
    pub fn add_vp(&mut self, grob: Grob, vp: Vp, t: usize, l: usize, b: usize, r: usize, z: Option<f64>, clip: bool, name: &str) {
        let z = self.resolve_z(z);
        self.items.push(Placed { t, l, b, r, z, clip, name: name.to_string(), grob, vp: Some(vp), key: None });
    }

    /// gtable's z rule: Inf is one more than the largest finite z so far, -Inf one less than the smallest.
    fn resolve_z(&self, z: Option<f64>) -> f64 {
        let finite: Vec<f64> = self.items.iter().map(|p| p.z).collect();
        let (zmin, zmax) = if finite.is_empty() {
            (1.0, 0.0)
        } else {
            (finite.iter().copied().fold(f64::INFINITY, f64::min), finite.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        };
        match z {
            None => zmax + 1.0,
            Some(v) if v == f64::NEG_INFINITY => zmin - 1.0,
            Some(v) if v == f64::INFINITY => zmax + 1.0,
            Some(v) => v,
        }
    }

    /// gtable_add_rows: rows inserted after row `pos` (0: at the top; usize::MAX or nrow: at the bottom).
    pub fn add_rows(&mut self, heights: &[Unit], pos: usize) {
        let pos = pos.min(self.nrow());
        let n = heights.len();
        for (k, h) in heights.iter().enumerate() {
            self.heights.insert(pos + k, *h);
        }
        for p in &mut self.items {
            if p.t > pos {
                p.t += n;
            }
            if p.b > pos {
                p.b += n;
            }
        }
    }

    /// gtable_add_cols: columns inserted after column `pos` (0: at the left; usize::MAX: at the right).
    pub fn add_cols(&mut self, widths: &[Unit], pos: usize) {
        let pos = pos.min(self.ncol());
        let n = widths.len();
        for (k, w) in widths.iter().enumerate() {
            self.widths.insert(pos + k, *w);
        }
        for p in &mut self.items {
            if p.l > pos {
                p.l += n;
            }
            if p.r > pos {
                p.r += n;
            }
        }
    }

    /// gtable_add_padding (top, right, bottom, left)
    pub fn add_padding(&mut self, t: f64, r: f64, b: f64, l: f64) {
        self.add_rows(&[Unit::Abs(t)], 0);
        self.add_cols(&[Unit::Abs(r)], usize::MAX);
        self.add_rows(&[Unit::Abs(b)], usize::MAX);
        self.add_cols(&[Unit::Abs(l)], 0);
    }

    /// gtable_add_col_space
    pub fn add_col_space(&mut self, w: f64) {
        let n = self.ncol().saturating_sub(1);
        for i in (1..=n).rev() {
            self.add_cols(&[Unit::Abs(w)], i);
        }
    }

    /// gtable_add_row_space
    pub fn add_row_space(&mut self, h: f64) {
        let n = self.nrow().saturating_sub(1);
        for i in (1..=n).rev() {
            self.add_rows(&[Unit::Abs(h)], i);
        }
    }

    /// `gt[rows, cols]` for increasing 1-based index lists: the grobs whose cells lie inside both are kept.
    pub fn subset(&self, rows: &[usize], cols: &[usize]) -> Gtable {
        let heights = rows.iter().map(|&i| self.heights[i - 1]).collect();
        let widths = cols.iter().map(|&j| self.widths[j - 1]).collect();
        let row_pos = |i: usize| rows.iter().position(|&r| r == i);
        let col_pos = |j: usize| cols.iter().position(|&c| c == j);
        let mut items = Vec::new();
        for p in &self.items {
            if let (Some(t), Some(b), Some(l), Some(r)) = (row_pos(p.t), row_pos(p.b), col_pos(p.l), col_pos(p.r)) {
                // the cells between must be kept too (contiguous ranges, as patchwork uses it)
                let mut q = p.clone();
                q.t = t + 1;
                q.b = b + 1;
                q.l = l + 1;
                q.r = r + 1;
                items.push(q);
            }
        }
        Gtable { widths, heights, items, vp: None }
    }

    /// `gt[rows, cols]` for 1-based inclusive ranges.
    pub fn sub(&self, r0: usize, r1: usize, c0: usize, c1: usize) -> Gtable {
        let rows: Vec<usize> = (r0..=r1).collect();
        let cols: Vec<usize> = (c0..=c1).collect();
        self.subset(&rows, &cols)
    }

    /// Removes the grobs for which `remove` holds.
    pub fn retain(&mut self, keep: impl Fn(&Placed) -> bool) {
        self.items.retain(|p| keep(p));
    }

    /// Sets the identity key (see `Placed::key`) of every grob whose name starts with `prefix`.
    pub fn set_keys(&mut self, prefix: &str, key: &str) {
        for p in self.items.iter_mut().filter(|p| p.name.starts_with(prefix)) {
            p.key = Some(key.to_string());
        }
    }

    pub fn find(&self, name: &str) -> Option<&Placed> {
        self.items.iter().find(|p| p.name == name)
    }

    /// The column widths and row heights in a box of w x h (grid's layout without respect).
    pub fn resolve(units: &[Unit], total: f64) -> Vec<f64> {
        let abs: f64 = units.iter().map(|u| u.abs()).sum();
        let null: f64 = units.iter().map(|u| if let Unit::Null(v) = u { *v } else { 0.0 }).sum();
        let left = (total - abs).max(0.0);
        units
            .iter()
            .map(|u| match u {
                Unit::Abs(v) => *v,
                Unit::Null(v) => {
                    if null > 0.0 {
                        left * v / null
                    } else {
                        0.0
                    }
                }
            })
            .collect()
    }

    /// The box of every cell edge: x edges (ncol + 1) and y edges (nrow + 1) in the box.
    pub fn edges(&self, b: Bx) -> (Vec<f64>, Vec<f64>) {
        let ws = Gtable::resolve(&self.widths, b.w());
        let hs = Gtable::resolve(&self.heights, b.h());
        // grid centres the layout in its viewport when the sizes do not fill it (layout just "centre")
        let tw: f64 = ws.iter().sum();
        let th: f64 = hs.iter().sum();
        let mut xs = vec![b.x0 + (b.w() - tw) / 2.0];
        for w in ws {
            xs.push(xs.last().unwrap() + w);
        }
        let mut ys = vec![b.y0 + (b.h() - th) / 2.0];
        for h in hs {
            ys.push(ys.last().unwrap() + h);
        }
        (xs, ys)
    }

    /// Draws the table in box `b`: its grobs in z order (ties in the order added).
    pub fn draw(&self, c: &mut Canvas, b: Bx) {
        let b = match &self.vp {
            Some(vp) => vp.apply(b),
            None => b,
        };
        let (xs, ys) = self.edges(b);
        let mut order: Vec<usize> = (0..self.items.len()).collect();
        order.sort_by(|&i, &j| self.items[i].z.total_cmp(&self.items[j].z));
        for i in order {
            let p = &self.items[i];
            if p.grob.is_zero() {
                continue;
            }
            let cell = Bx::new(xs[p.l - 1], ys[p.t - 1], xs[p.r], ys[p.b]);
            let bx = match &p.vp {
                Some(vp) => vp.apply(cell),
                None => cell,
            };
            if p.clip {
                c.push_clip(bx);
                p.grob.draw(c, bx);
                c.pop_clip();
            } else {
                p.grob.draw(c, bx);
            }
        }
    }

    /// Every grob's box after layout (depth-first), for tests against R's viewports: (name path, box).
    pub fn boxes(&self, b: Bx, prefix: &str, out: &mut Vec<(String, Bx)>) {
        let b = match &self.vp {
            Some(vp) => vp.apply(b),
            None => b,
        };
        let (xs, ys) = self.edges(b);
        for p in &self.items {
            let cell = Bx::new(xs[p.l - 1], ys[p.t - 1], xs[p.r], ys[p.b]);
            let bx = match &p.vp {
                Some(vp) => vp.apply(cell),
                None => cell,
            };
            let name = if prefix.is_empty() { p.name.clone() } else { format!("{prefix}::{}", p.name) };
            out.push((name.clone(), bx));
            if let Grob::Table(t) = &p.grob {
                t.boxes(bx, &name, out);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout_and_z_as_gtable() {
        let mut g = Gtable::new(vec![Unit::Abs(10.0), Unit::Null(1.0), Unit::Null(3.0)], vec![Unit::Abs(5.0), Unit::Null(1.0)]);
        g.add(Grob::Zero, 1, 1, 1, 1, None, false, "a");
        g.add(Grob::Zero, 2, 2, 2, 3, Some(f64::NEG_INFINITY), false, "b");
        g.add(Grob::Zero, 2, 2, 2, 2, None, false, "c");
        assert_eq!(g.items.iter().map(|p| p.z).collect::<Vec<_>>(), vec![1.0, 0.0, 2.0]);
        let (xs, ys) = g.edges(Bx::new(0.0, 0.0, 110.0, 25.0));
        assert_eq!(xs, vec![0.0, 10.0, 35.0, 110.0]);
        assert_eq!(ys, vec![0.0, 5.0, 25.0]);
        g.add_rows(&[Unit::Abs(2.0)], 1);
        assert_eq!((g.items[1].t, g.items[1].b), (3, 3));
        let s = g.sub(1, 3, 2, 3);
        assert_eq!(s.items.len(), 2);
        assert_eq!((s.items[0].l, s.items[0].r), (1, 2));
    }

    #[test]
    fn viewport_justification() {
        let cell = Bx::new(0.0, 0.0, 100.0, 50.0);
        // viewport(y = 0, just = "bottom", height = 10)
        let vp = Vp { x: Pos::npc(0.5), y: Pos::npc(0.0), w: None, h: Some(10.0), hjust: 0.5, vjust: 0.0 };
        assert_eq!(vp.apply(cell), Bx::new(0.0, 40.0, 100.0, 50.0));
        // viewport(x = 1, just = "right", width = 30)
        let vp = Vp { x: Pos::npc(1.0), y: Pos::npc(0.5), w: Some(30.0), h: None, hjust: 1.0, vjust: 0.5 };
        assert_eq!(vp.apply(cell), Bx::new(70.0, 0.0, 100.0, 50.0));
    }
}
