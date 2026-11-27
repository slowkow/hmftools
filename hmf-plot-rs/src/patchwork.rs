//! patchwork (1.3.2, which agrees with 1.2.0 for what is used here) on gtables: a ggplot's table made compliant
//! (add_strips, simplify_gt: the facet panels nested into one cell), plots combined in a design grid with their
//! panels, axes and titles aligned (build_patchwork: table_dims, set_grob_sizes, set_panel_dimensions), nested
//! patchworks, `free()` (free_label, free_panel), `wrap_elements()` and `plot_annotation()` (annotate_table).
//!
//! Every plot's table has patchwork's 18 rows x 15 columns, the panel at row 10, column 8.

use std::rc::Rc;

use crate::gtable::{Grob, Gtable, Pos, Unit, Vp, sum_abs};

pub const TABLE_ROWS: usize = 18;
pub const TABLE_COLS: usize = 15;
pub const PANEL_ROW: usize = 10;
pub const PANEL_COL: usize = 8;

/// A patch area of a design (1-based, inclusive).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Area {
    pub t: usize,
    pub l: usize,
    pub b: usize,
    pub r: usize,
}

impl Area {
    pub fn cell(t: usize, l: usize) -> Area {
        Area { t, l, b: t, r: l }
    }
}

/// What kind of table a patch is (patchwork's classes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// a ggplot's table (simplified)
    Plot,
    /// a nested patchwork, wrapped by simplify_gt.gtable_patchwork
    PatchworkSimple,
    /// a patch (wrap_elements)
    PatchGrob,
}

/// A patch ready to be placed: its compliant table and kind.
#[derive(Clone)]
pub struct Patch {
    pub gt: Gtable,
    pub kind: Kind,
}

fn find_panel(gt: &Gtable) -> (usize, usize, usize, usize) {
    let ps: Vec<_> = gt.items.iter().filter(|p| p.name.starts_with("panel")).collect();
    (
        ps.iter().map(|p| p.t).min().unwrap_or(PANEL_ROW),
        ps.iter().map(|p| p.l).min().unwrap_or(PANEL_COL),
        ps.iter().map(|p| p.b).max().unwrap_or(PANEL_ROW),
        ps.iter().map(|p| p.r).max().unwrap_or(PANEL_COL),
    )
}

fn has_name(gt: &Gtable, prefix: &str) -> bool {
    gt.items.iter().any(|p| p.name.starts_with(prefix))
}

/// gt[-row, ] for one row (1-based)
fn drop_row(gt: &Gtable, row: usize) -> Gtable {
    let rows: Vec<usize> = (1..=gt.nrow()).filter(|&r| r != row).collect();
    let cols: Vec<usize> = (1..=gt.ncol()).collect();
    let mut s = gt.subset(&rows, &cols);
    s.vp = gt.vp;
    s
}

fn drop_col(gt: &Gtable, col: usize) -> Gtable {
    let rows: Vec<usize> = (1..=gt.nrow()).collect();
    let cols: Vec<usize> = (1..=gt.ncol()).filter(|&c| c != col).collect();
    let mut s = gt.subset(&rows, &cols);
    s.vp = gt.vp;
    s
}

/// add_strips: zero-sized rows and columns for missing strips; an outside strip's padding merged into its axis.
pub fn add_strips(mut gt: Gtable) -> Gtable {
    let (t, l, b, r) = find_panel(&gt);
    let names = |gt: &Gtable, p: &str| gt.items.iter().filter(|x| x.name.starts_with(p)).map(|x| (x.t, x.l, x.b, x.r)).collect::<Vec<_>>();
    // find_strip_pos
    let outside = {
        let st = names(&gt, "strip-t");
        let sr = names(&gt, "strip-r");
        let sb = names(&gt, "strip-b");
        let sl = names(&gt, "strip-l");
        (!st.is_empty() && t - st.iter().map(|x| x.0).min().unwrap() != 1)
            || (!sr.is_empty() && sr.iter().map(|x| x.3).max().unwrap() - r != 1)
            || (!sb.is_empty() && sb.iter().map(|x| x.2).max().unwrap() - b != 1)
            || (!sl.is_empty() && l - sl.iter().map(|x| x.1).min().unwrap() != 1)
    };
    let sp = if outside { 2 } else { 0 };
    if !has_name(&gt, "strip-b") {
        gt.add_rows(&[Unit::ZERO], b + sp);
    } else if sp == 2 && !gt.items.iter().any(|x| x.b == b + 2) {
        gt.heights[b] = Unit::Abs(gt.heights[b].abs() + gt.heights[b + 1].abs());
        gt = drop_row(&gt, b + 2);
    }
    if !has_name(&gt, "strip-t") {
        gt.add_rows(&[Unit::ZERO], t - 1 - sp);
    } else if sp == 2 && !gt.items.iter().any(|x| x.t == t - 2) {
        gt.heights[t - 2] = Unit::Abs(gt.heights[t - 2].abs() + gt.heights[t - 3].abs());
        gt = drop_row(&gt, t - 2);
    }
    if !has_name(&gt, "strip-r") {
        gt.add_cols(&[Unit::ZERO], r + sp);
    } else if sp == 2 && !gt.items.iter().any(|x| x.r == r + 2) {
        gt.widths[r] = Unit::Abs(gt.widths[r].abs() + gt.widths[r + 1].abs());
        gt = drop_col(&gt, r + 2);
    }
    if !has_name(&gt, "strip-l") {
        gt.add_cols(&[Unit::ZERO], l - 1 - sp);
    } else if sp == 2 && !gt.items.iter().any(|x| x.l == l - 2) {
        gt.widths[l - 2] = Unit::Abs(gt.widths[l - 2].abs() + gt.widths[l - 3].abs());
        gt = drop_col(&gt, l - 2);
    }
    gt
}

/// simplify_gt.gtable: a ggplot table with several panels (or strips) gets its panel area as one nested table.
pub fn simplify_plot(gt: Gtable) -> Gtable {
    let (t, l, b, r) = find_panel(&gt);
    if t == b && l == r && !has_name(&gt, "strip-") {
        return gt;
    }
    let p_rows: Vec<usize> = (t..=b).collect();
    let p_cols: Vec<usize> = (l..=r).collect();
    let panels = gt.subset(&p_rows, &p_cols);
    let o_rows: Vec<usize> = (1..=gt.nrow()).filter(|i| !p_rows.contains(i)).collect();
    let o_cols: Vec<usize> = (1..=gt.ncol()).filter(|j| !p_cols.contains(j)).collect();
    let mut gt_new = gt.subset(&o_rows, &o_cols);
    let w_abs = p_cols.iter().all(|&j| !gt.widths[j - 1].is_null());
    let h_abs = p_rows.iter().all(|&i| !gt.heights[i - 1].is_null());
    let new_w = if w_abs { Unit::Abs(p_cols.iter().map(|&j| gt.widths[j - 1].abs()).sum()) } else { Unit::Null(1.0) };
    let new_h = if h_abs { Unit::Abs(p_rows.iter().map(|&i| gt.heights[i - 1].abs()).sum()) } else { Unit::Null(1.0) };
    gt_new.add_rows(&[new_h], t - 1);
    gt_new.add_cols(&[new_w], l - 1);
    simplify_free(&gt, gt_new, panels, (t, b), (l, r))
}

fn max_z(t: &Gtable) -> f64 {
    t.items.iter().map(|p| p.z).fold(f64::NEG_INFINITY, f64::max)
}

fn joined_names(t: &Gtable) -> String {
    t.items.iter().map(|p| p.name.clone()).collect::<Vec<_>>().join(", ")
}

fn simplify_free(gt: &Gtable, mut gt_new: Gtable, panels: Gtable, rows: (usize, usize), cols: (usize, usize)) -> Gtable {
    let p_cols: Vec<usize> = (cols.0..=cols.1).collect();
    let p_rows: Vec<usize> = (rows.0..=rows.1).collect();
    if p_cols.len() == 1 {
        let c = p_cols[0];
        for p in gt.items.iter().filter(|p| p.l == c && p.r == c && p.b < rows.0) {
            gt_new.add(p.grob.clone(), p.t, c, p.b, c, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
        }
        let b_mod = rows.1 - rows.0;
        for p in gt.items.iter().filter(|p| p.l == c && p.r == c && p.t > rows.1) {
            gt_new.add(p.grob.clone(), p.t - b_mod, c, p.b - b_mod, c, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
        }
    } else {
        for i in 1..=gt.nrow() {
            let (ii, top) = if i >= rows.0 {
                if i <= rows.1 {
                    continue;
                }
                (i - (rows.1 - rows.0), false)
            } else {
                (i, true)
            };
            let mut table = gt.subset(&[i], &p_cols);
            if table.items.is_empty() {
                continue;
            }
            let name = joined_names(&table);
            let h = sum_abs(&table.heights);
            table.vp = Some(if top {
                Vp { x: Pos::npc(0.5), y: Pos::npc(0.0), w: None, h: Some(h), hjust: 0.5, vjust: 0.0 }
            } else {
                Vp { x: Pos::npc(0.5), y: Pos::npc(1.0), w: None, h: Some(h), hjust: 0.5, vjust: 1.0 }
            });
            let z = max_z(&table);
            gt_new.add(Grob::table(table), ii, cols.0, ii, cols.0, Some(z), false, &name);
        }
    }
    if p_rows.len() == 1 {
        let rr = p_rows[0];
        for p in gt.items.iter().filter(|p| p.t == rr && p.b == rr && p.r < cols.0) {
            gt_new.add(p.grob.clone(), rr, p.l, rr, p.r, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
        }
        let r_mod = cols.1 - cols.0;
        for p in gt.items.iter().filter(|p| p.t == rr && p.b == rr && p.l > cols.1) {
            gt_new.add(p.grob.clone(), rr, p.l - r_mod, rr, p.r - r_mod, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
        }
    } else {
        for i in 1..=gt.ncol() {
            let (ii, left) = if i >= cols.0 {
                if i <= cols.1 {
                    continue;
                }
                (i - (cols.1 - cols.0), false)
            } else {
                (i, true)
            };
            let mut table = gt.subset(&p_rows, &[i]);
            if table.items.is_empty() {
                continue;
            }
            let name = joined_names(&table);
            let w = sum_abs(&table.widths);
            table.vp = Some(if left {
                Vp { x: Pos::npc(1.0), y: Pos::npc(0.5), w: Some(w), h: None, hjust: 1.0, vjust: 0.5 }
            } else {
                Vp { x: Pos::npc(0.0), y: Pos::npc(0.5), w: Some(w), h: None, hjust: 0.0, vjust: 0.5 }
            });
            let z = max_z(&table);
            gt_new.add(Grob::table(table), rows.0, ii, rows.0, ii, Some(z), false, &name);
        }
    }
    let name = format!("panel; {}", joined_names(&panels));
    gt_new.add(Grob::table(panels), rows.0, cols.0, rows.0, cols.0, Some(1.0), false, &name);
    gt_new
}

/// plot_table.ggplot: add_strips, then (guides not collected, the five guide boxes kept) simplify_gt.
pub fn plot(gt: Gtable) -> Patch {
    Patch { gt: simplify_plot(add_strips(gt)), kind: Kind::Plot }
}

/// simplify_gt.gtable_patchwork: a nested patchwork as one table of patchwork's size, the whole patchwork drawn
/// over it.
pub fn nested(gt: Gtable) -> Patch {
    let (t, l, b, r) = find_panel(&gt);
    let new_w = if (l..=r).all(|j| !gt.widths[j - 1].is_null()) { Unit::Abs((l..=r).map(|j| gt.widths[j - 1].abs()).sum()) } else { Unit::Null(1.0) };
    let new_h = if (t..=b).all(|i| !gt.heights[i - 1].is_null()) { Unit::Abs((t..=b).map(|i| gt.widths[i - 1].abs()).sum()) } else { Unit::Null(1.0) };
    let mut widths: Vec<Unit> = gt.widths[..l - 1].to_vec();
    widths.push(new_w);
    widths.extend_from_slice(&gt.widths[r..]);
    let mut heights: Vec<Unit> = gt.heights[..t - 1].to_vec();
    heights.push(new_h);
    heights.extend_from_slice(&gt.heights[b..]);
    let mut gt_new = Gtable::new(widths, heights);
    gt_new.add(Grob::Zero, PANEL_ROW, PANEL_COL, PANEL_ROW, PANEL_COL, None, false, "panel-nested-patchwork");
    let (nr, nc) = (gt_new.nrow(), gt_new.ncol());
    gt_new.add(Grob::table(gt), 1, 1, nr, nc, None, false, "patchwork-table");
    Patch { gt: gt_new, kind: Kind::PatchworkSimple }
}

/// merge.patchwork, which `free()` ends with when given a patchwork: a new patchwork whose only patch is this one
/// (built alone in a 1 x 1 design), nested where it is used. It draws the same as the patch itself; the table tree
/// gets the extra level R's has (a `panel-area`, a `panel-nested-patchwork` and a `patchwork-table` around it).
/// `free()` on a ggplot calls merge.ggplot, which returns the plot unchanged: do not call this for plots.
pub fn merge(p: Patch) -> Patch {
    nested(build(vec![p], &Layout::column(1, vec![None])))
}

/// wrap_elements(panel = grob): patchwork's empty patch table with the grob in the panel cell (clipped), the
/// plot margins of the theme in the first and last rows and columns, and the plot background.
pub fn wrap_panel(grob: Grob, margin: [f64; 4], background: Option<Grob>) -> Patch {
    let mut widths = vec![Unit::ZERO; TABLE_COLS];
    widths[PANEL_COL - 1] = Unit::Null(1.0);
    let mut heights = vec![Unit::ZERO; TABLE_ROWS];
    heights[PANEL_ROW - 1] = Unit::Null(1.0);
    widths[0] = Unit::Abs(margin[3]);
    widths[TABLE_COLS - 1] = Unit::Abs(margin[1]);
    heights[0] = Unit::Abs(margin[0]);
    heights[TABLE_ROWS - 1] = Unit::Abs(margin[2]);
    let mut gt = Gtable::new(widths, heights);
    gt.add(Grob::Zero, PANEL_ROW, PANEL_COL, PANEL_ROW, PANEL_COL, Some(f64::NEG_INFINITY), false, "panel_patch");
    if let Some(bg) = background {
        gt.add(bg, 1, 1, TABLE_ROWS, TABLE_COLS, Some(-100.0), true, "background");
    }
    gt.add(grob, PANEL_ROW, PANEL_COL, PANEL_ROW, PANEL_COL, None, true, "panel");
    Patch { gt, kind: Kind::PatchGrob }
}

/// The layout of a patchwork: a design of areas (one per patch, in order), and relative widths and heights of
/// the design's columns and rows (None: 1).
pub struct Layout {
    pub design: Vec<Area>,
    pub widths: Vec<Option<f64>>,
    pub heights: Vec<Option<f64>>,
}

impl Layout {
    /// wrap_plots(ncol = 1) or (nrow = 1): one patch per row or column
    pub fn column(n: usize, heights: Vec<Option<f64>>) -> Layout {
        Layout { design: (1..=n).map(|i| Area::cell(i, 1)).collect(), widths: vec![None], heights }
    }
    pub fn row(n: usize) -> Layout {
        Layout { design: (1..=n).map(|i| Area::cell(1, i)).collect(), widths: vec![None; n], heights: vec![None] }
    }
    /// A design string as patchwork reads it: one letter per patch, rows of equal length ('#' empty).
    pub fn design(design: &str) -> Layout {
        let rows: Vec<Vec<char>> = design.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).map(|l| l.chars().collect()).collect();
        let mut letters: Vec<char> = rows.iter().flatten().copied().filter(|&c| c != '#').collect();
        letters.sort();
        letters.dedup();
        let areas = letters
            .iter()
            .map(|&ch| {
                let cells: Vec<(usize, usize)> = rows.iter().enumerate().flat_map(|(i, r)| r.iter().enumerate().filter(|(_, c)| **c == ch).map(move |(j, _)| (i + 1, j + 1))).collect();
                Area { t: cells.iter().map(|c| c.0).min().unwrap(), l: cells.iter().map(|c| c.1).min().unwrap(), b: cells.iter().map(|c| c.0).max().unwrap(), r: cells.iter().map(|c| c.1).max().unwrap() }
            })
            .collect();
        let (nr, nc) = (rows.len(), rows.first().map(|r| r.len()).unwrap_or(0));
        Layout { design: areas, widths: vec![None; nc], heights: vec![None; nr] }
    }
}

/// build_patchwork: the patches in a grid of 18 x 15 cells per design cell, widths and heights aligned across
/// patches sharing a design column or row, panel sizes from the layout.
pub fn build(patches: Vec<Patch>, layout: &Layout) -> Gtable {
    let nrow = layout.design.iter().map(|a| a.b).max().unwrap_or(1).max(layout.heights.len());
    let ncol = layout.design.iter().map(|a| a.r).max().unwrap_or(1).max(layout.widths.len());
    let mut gt_new = Gtable::new(vec![Unit::Null(0.0); TABLE_COLS * ncol], vec![Unit::Null(0.0); TABLE_ROWS * nrow]);
    let design = &layout.design[..patches.len()];
    let mut zoff = 0.0;
    for (i, p) in patches.iter().enumerate() {
        let loc = design[i];
        for it in &p.gt.items {
            let mut q = it.clone();
            if q.name != "background" {
                q.z += zoff;
            }
            q.t += if q.t <= PANEL_ROW { (loc.t - 1) * TABLE_ROWS } else { (loc.b - 1) * TABLE_ROWS };
            q.l += if q.l <= PANEL_COL { (loc.l - 1) * TABLE_COLS } else { (loc.r - 1) * TABLE_COLS };
            q.b += if q.b < PANEL_ROW { (loc.t - 1) * TABLE_ROWS } else { (loc.b - 1) * TABLE_ROWS };
            q.r += if q.r < PANEL_COL { (loc.l - 1) * TABLE_COLS } else { (loc.r - 1) * TABLE_COLS };
            q.name = format!("{}-{}", q.name, i + 1);
            gt_new.items.push(q);
        }
        zoff += max_z(&p.gt);
    }
    // table_dims
    let widths: Vec<f64> = (1..=ncol * TABLE_COLS)
        .map(|i| {
            let area = (i - 1) / TABLE_COLS + 1;
            let col_loc = match i % TABLE_COLS {
                0 => TABLE_COLS,
                c => c,
            };
            design
                .iter()
                .zip(&patches)
                .filter(|(a, _)| if col_loc <= PANEL_COL { a.l == area } else { a.r == area })
                .map(|(_, p)| p.gt.widths[col_loc - 1].abs())
                .fold(0.0, f64::max)
        })
        .collect();
    let heights: Vec<f64> = (1..=nrow * TABLE_ROWS)
        .map(|i| {
            let area = (i - 1) / TABLE_ROWS + 1;
            let row_loc = match i % TABLE_ROWS {
                0 => TABLE_ROWS,
                r => r,
            };
            design
                .iter()
                .zip(&patches)
                .filter(|(a, _)| if row_loc <= PANEL_ROW { a.t == area } else { a.b == area })
                .map(|(_, p)| p.gt.heights[row_loc - 1].abs())
                .fold(0.0, f64::max)
        })
        .collect();
    // set_grob_sizes: nested patchworks take the aligned border sizes
    let mut k = 0;
    for (i, p) in patches.iter().enumerate() {
        let n = p.gt.items.len();
        if p.kind == Kind::PatchworkSimple {
            let loc = design[i];
            let l = (loc.l - 1) * TABLE_COLS;
            let lw: Vec<Unit> = widths[l..l + PANEL_COL - 1].iter().map(|&w| Unit::Abs(w)).collect();
            let r = (loc.r - 1) * TABLE_COLS;
            let rw: Vec<Unit> = widths[r + PANEL_COL..r + TABLE_COLS].iter().map(|&w| Unit::Abs(w)).collect();
            let t = (loc.t - 1) * TABLE_ROWS;
            let th: Vec<Unit> = heights[t..t + PANEL_ROW - 1].iter().map(|&h| Unit::Abs(h)).collect();
            let b = (loc.b - 1) * TABLE_ROWS;
            let bh: Vec<Unit> = heights[b + PANEL_ROW..b + TABLE_ROWS].iter().map(|&h| Unit::Abs(h)).collect();
            for it in gt_new.items[k..k + n].iter_mut() {
                if it.name.starts_with("patchwork-table") {
                    if let Grob::Table(t) = &mut it.grob {
                        let tt = Rc::make_mut(t);
                        set_border_sizes(tt, Some(&lw), Some(&rw), Some(&th), Some(&bh));
                    }
                }
            }
        }
        k += n;
    }
    gt_new.widths = widths.iter().map(|&w| Unit::Abs(w)).collect();
    gt_new.heights = heights.iter().map(|&h| Unit::Abs(h)).collect();
    // set_panel_dimensions (no fixed aspect ratios here)
    for c in 0..ncol {
        let w = layout.widths.get(c).copied().flatten();
        if w.is_none() {
            // an absolute (non-zero) panel width of a patch alone in this column fixes it
            let fixed: Vec<f64> = design.iter().zip(&patches).filter(|(a, _)| a.l == c + 1 && a.r == c + 1).map(|(_, p)| p.gt.widths[PANEL_COL - 1]).filter(|u| !u.is_null() && u.abs() != 0.0).map(|u| u.abs()).collect();
            if let Some(m) = fixed.iter().copied().reduce(f64::max) {
                gt_new.widths[c * TABLE_COLS + PANEL_COL - 1] = Unit::Abs(m);
                continue;
            }
        }
        gt_new.widths[c * TABLE_COLS + PANEL_COL - 1] = Unit::Null(w.unwrap_or(1.0));
    }
    for r in 0..nrow {
        let h = layout.heights.get(r).copied().flatten();
        if h.is_none() {
            let fixed: Vec<f64> = design.iter().zip(&patches).filter(|(a, _)| a.t == r + 1 && a.b == r + 1).map(|(_, p)| p.gt.heights[PANEL_ROW - 1]).filter(|u| !u.is_null() && u.abs() != 0.0).map(|u| u.abs()).collect();
            if let Some(m) = fixed.iter().copied().reduce(f64::max) {
                gt_new.heights[r * TABLE_ROWS + PANEL_ROW - 1] = Unit::Abs(m);
                continue;
            }
        }
        gt_new.heights[r * TABLE_ROWS + PANEL_ROW - 1] = Unit::Null(h.unwrap_or(1.0));
    }
    let (b, r) = (PANEL_ROW + TABLE_ROWS * (nrow - 1), PANEL_COL + TABLE_COLS * (ncol - 1));
    gt_new.add(Grob::Zero, PANEL_ROW, PANEL_COL, b, r, Some(-1.0), false, "panel-area");
    gt_new
}

/// set_border_sizes: the outer widths and heights of a nested patchwork, recursively into patchworks touching the
/// border.
pub fn set_border_sizes(gt: &mut Gtable, l: Option<&[Unit]>, r: Option<&[Unit]>, t: Option<&[Unit]>, b: Option<&[Unit]>) {
    let (nr, nc) = (gt.nrow(), gt.ncol());
    if let Some(l) = l {
        gt.widths[..l.len()].copy_from_slice(l);
    }
    if let Some(r) = r {
        gt.widths[nc - r.len()..].copy_from_slice(r);
    }
    if let Some(t) = t {
        gt.heights[..t.len()].copy_from_slice(t);
    }
    if let Some(b) = b {
        gt.heights[nr - b.len()..].copy_from_slice(b);
    }
    for it in gt.items.iter_mut() {
        if !it.name.starts_with("patchwork-table") {
            continue;
        }
        let (il, ir, itt, ib) = (it.l == 1, it.r == nc, it.t == 1, it.b == nr);
        if let Grob::Table(tb) = &mut it.grob {
            let tt = Rc::make_mut(tb);
            set_border_sizes(tt, l.filter(|_| il), r.filter(|_| ir), t.filter(|_| itt), b.filter(|_| ib));
        }
    }
}

fn panel_positions(n: usize, by: usize, first: usize) -> Vec<usize> {
    (0..n / by).map(|k| first + k * by).collect()
}

/// liberate_area: the grobs inside rows top..=bottom, columns left..=right made one grob (a sub-table spanning the
/// area, drawn in `vp`).
fn liberate_area(gt: &mut Gtable, top: usize, right: usize, bottom: usize, left: usize, name: Option<&str>, vp: Option<Vp>) {
    let inside = |p: &crate::gtable::Placed| p.l >= left && p.t >= top && p.r <= right && p.b <= bottom;
    if !gt.items.iter().any(|p| inside(p)) {
        return;
    }
    let mut lib = gt.sub(top, bottom, left, right);
    lib.vp = vp;
    let nm = name.map(|s| s.to_string()).unwrap_or_else(|| lib.items.iter().map(|p| p.name.clone()).collect::<Vec<_>>().join("; "));
    gt.retain(|p| !inside(p));
    let z = max_z(&lib);
    gt.add(Grob::table(lib), top, left, bottom, right, Some(z), false, &nm);
}

/// free_label: the axis titles, strips and axes on each side made one grob drawn at the panel, so they keep their
/// own size and stay next to the panel when the patchwork aligns the rows and columns around it.
pub fn free_label(gt: &mut Gtable, has_side: [bool; 4]) {
    let (nr, nc) = (gt.nrow(), gt.ncol());
    for it in gt.items.iter_mut() {
        if !it.name.starts_with("patchwork-table") {
            continue;
        }
        let loc = [it.t == 1, it.r == nc, it.b == nr, it.l == 1];
        let sides = [loc[0] && has_side[0], loc[1] && has_side[1], loc[2] && has_side[2], loc[3] && has_side[3]];
        if !sides.iter().any(|&s| s) {
            continue;
        }
        if let Grob::Table(tb) = &mut it.grob {
            free_label(Rc::make_mut(tb), sides);
        }
    }
    let pcols = panel_positions(nc, TABLE_COLS, PANEL_COL);
    let prows = panel_positions(nr, TABLE_ROWS, PANEL_ROW);
    let pw: Vec<Unit> = pcols.iter().map(|&c| gt.widths[c - 1]).collect();
    let ph: Vec<Unit> = prows.iter().map(|&r| gt.heights[r - 1]).collect();
    for &c in &pcols {
        if gt.widths[c - 1].value() == 0.0 {
            gt.widths[c - 1] = Unit::Null(1.0);
        }
    }
    for &r in &prows {
        if gt.heights[r - 1].value() == 0.0 {
            gt.heights[r - 1] = Unit::Null(1.0);
        }
    }
    let top = PANEL_ROW;
    let right = nc - (TABLE_COLS - PANEL_COL);
    let bottom = nr - (TABLE_ROWS - PANEL_ROW);
    let left = PANEL_COL;
    let hsum = |gt: &Gtable, a: usize, b: usize| (a..=b).map(|i| gt.heights[i - 1].abs()).sum::<f64>();
    let wsum = |gt: &Gtable, a: usize, b: usize| (a..=b).map(|i| gt.widths[i - 1].abs()).sum::<f64>();
    if has_side[0] {
        let h = hsum(gt, top - 3, top - 1);
        liberate_area(gt, top - 3, right, top - 1, left, None, Some(Vp { x: Pos::npc(0.5), y: Pos::npc(0.0), w: None, h: Some(h), hjust: 0.5, vjust: 0.0 }));
    }
    if has_side[1] {
        let w = wsum(gt, right + 1, right + 3);
        liberate_area(gt, top, right + 3, bottom, right + 1, None, Some(Vp { x: Pos::npc(0.0), y: Pos::npc(0.5), w: Some(w), h: None, hjust: 0.0, vjust: 0.5 }));
    }
    if has_side[2] {
        let h = hsum(gt, bottom + 1, bottom + 3);
        liberate_area(gt, bottom + 1, right, bottom + 3, left, None, Some(Vp { x: Pos::npc(0.5), y: Pos::npc(1.0), w: None, h: Some(h), hjust: 0.5, vjust: 1.0 }));
    }
    if has_side[3] {
        let w = wsum(gt, left - 3, left - 1);
        liberate_area(gt, top, left - 1, bottom, left - 3, None, Some(Vp { x: Pos::npc(1.0), y: Pos::npc(0.5), w: Some(w), h: None, hjust: 1.0, vjust: 0.5 }));
    }
    for (k, &c) in pcols.iter().enumerate() {
        gt.widths[c - 1] = pw[k];
    }
    for (k, &r) in prows.iter().enumerate() {
        gt.heights[r - 1] = ph[k];
    }
}

/// free_panel: the panel area's sides not aligned with other patches (its axes and titles taken into the panel's
/// cell on those sides).
pub fn free_panel(gt: &mut Gtable, has_side: [bool; 4]) {
    let (nr, nc) = (gt.nrow(), gt.ncol());
    for it in gt.items.iter_mut() {
        if !it.name.starts_with("patchwork-table") {
            continue;
        }
        let loc = [it.t == 1, it.r == nc, it.b == nr, it.l == 1];
        let sides = [loc[0] && has_side[0], loc[1] && has_side[1], loc[2] && has_side[2], loc[3] && has_side[3]];
        if !sides.iter().any(|&s| s) {
            continue;
        }
        if let Grob::Table(tb) = &mut it.grob {
            free_panel(Rc::make_mut(tb), sides);
        }
    }
    let top = if has_side[0] { 3 } else { PANEL_ROW };
    let right = nc - if has_side[1] { 2 } else { TABLE_COLS - PANEL_COL };
    let bottom = nr - if has_side[2] { 2 } else { TABLE_ROWS - PANEL_ROW };
    let left = if has_side[3] { 3 } else { PANEL_COL };
    let pcols = panel_positions(nc, TABLE_COLS, PANEL_COL);
    let prows = panel_positions(nr, TABLE_ROWS, PANEL_ROW);
    let pw: Vec<Unit> = pcols.iter().map(|&c| gt.widths[c - 1]).collect();
    let ph: Vec<Unit> = prows.iter().map(|&r| gt.heights[r - 1]).collect();
    for &c in &pcols {
        if gt.widths[c - 1].value() == 0.0 {
            gt.widths[c - 1] = Unit::Null(1.0);
        }
    }
    for &r in &prows {
        if gt.heights[r - 1].value() == 0.0 {
            gt.heights[r - 1] = Unit::Null(1.0);
        }
    }
    liberate_area(gt, top, right, bottom, left, Some("free_panel"), None);
    let hsum = |gt: &Gtable, a: usize, b: usize| (a..=b).map(|i| gt.heights[i - 1].abs()).sum::<f64>();
    let wsum = |gt: &Gtable, a: usize, b: usize| (a..=b).map(|i| gt.widths[i - 1].abs()).sum::<f64>();
    let liberate_rows = |gt: &mut Gtable, t: usize, r: usize, b: usize, l: usize, align: f64| {
        let mut bands: Vec<(usize, usize)> = gt.items.iter().filter(|p| p.l >= l && p.t >= t && p.r <= r && p.b <= b).map(|p| (p.t, p.b)).collect();
        bands.dedup();
        let mut seen = Vec::new();
        for band in bands {
            if seen.contains(&band) {
                continue;
            }
            seen.push(band);
            let h = hsum(gt, band.0, band.1);
            liberate_area(gt, band.0, r, band.1, l, Some("free_row"), Some(Vp { x: Pos::npc(0.5), y: Pos::npc(align), w: None, h: Some(h), hjust: 0.5, vjust: align }));
        }
    };
    let liberate_cols = |gt: &mut Gtable, t: usize, r: usize, b: usize, l: usize, align: f64| {
        let bands: Vec<(usize, usize)> = gt.items.iter().filter(|p| p.l >= l && p.t >= t && p.r <= r && p.b <= b).map(|p| (p.l, p.r)).collect();
        let mut seen = Vec::new();
        for band in bands {
            if seen.contains(&band) {
                continue;
            }
            seen.push(band);
            let w = wsum(gt, band.0, band.1);
            liberate_area(gt, t, band.1, b, band.0, Some("free_col"), Some(Vp { x: Pos::npc(align), y: Pos::npc(0.5), w: Some(w), h: None, hjust: align, vjust: 0.5 }));
        }
    };
    if !has_side[0] && (has_side[1] || has_side[3]) {
        liberate_rows(gt, 3, right, top - 1, left, 0.0);
    }
    if !has_side[1] && (has_side[0] || has_side[2]) {
        let n = gt.ncol();
        liberate_cols(gt, top, n - 2, bottom, right + 1, 0.0);
    }
    if !has_side[2] && (has_side[1] || has_side[3]) {
        let n = gt.nrow();
        liberate_rows(gt, bottom + 1, right, n - 2, left, 1.0);
    }
    if !has_side[3] && (has_side[0] || has_side[2]) {
        liberate_cols(gt, top, left - 1, bottom, 3, 1.0);
    }
    let pcmin = *pcols.iter().min().unwrap();
    let pcmax = *pcols.iter().max().unwrap();
    for c in left..=right {
        if c < pcmin || c > pcmax {
            gt.widths[c - 1] = Unit::ZERO;
        }
    }
    for (k, &c) in pcols.iter().enumerate() {
        gt.widths[c - 1] = pw[k];
    }
    let prmin = *prows.iter().min().unwrap();
    let prmax = *prows.iter().max().unwrap();
    for r in top..=bottom {
        if r < prmin || r > prmax {
            gt.heights[r - 1] = Unit::ZERO;
        }
    }
    for (k, &r) in prows.iter().enumerate() {
        gt.heights[r - 1] = ph[k];
    }
}

/// The annotation of a patchwork (plot_annotation): the margins of the annotation theme around it, a title and
/// subtitle row (with their grobs) when given, and the background.
pub struct Annotation {
    /// plot.margin of the annotation theme (top, right, bottom, left)
    pub margin: [f64; 4],
    pub title: Option<(Grob, f64)>,
    pub subtitle: Option<(Grob, f64)>,
    pub background: Option<Grob>,
}

/// annotate_table
pub fn annotate(mut table: Gtable, a: Annotation) -> Gtable {
    let mz = max_z(&table);
    if a.title.is_some() || a.subtitle.is_some() {
        let th = a.title.as_ref().map(|t| t.1).unwrap_or(0.0);
        let sh = a.subtitle.as_ref().map(|t| t.1).unwrap_or(0.0);
        table.add_rows(&[Unit::Abs(a.margin[0]), Unit::Abs(th), Unit::Abs(sh)], 0);
        let nc = table.ncol();
        table.add(a.title.map(|t| t.0).unwrap_or(Grob::Zero), 2, 2, 2, nc - 1, Some(mz + 3.0), false, "title");
        table.add(a.subtitle.map(|t| t.0).unwrap_or(Grob::Zero), 3, 2, 3, nc - 1, Some(mz + 2.0), false, "subtitle");
    } else {
        table.add_rows(&[Unit::Abs(a.margin[0])], 0);
    }
    table.add_rows(&[Unit::Abs(a.margin[2])], usize::MAX);
    table.add_cols(&[Unit::Abs(a.margin[3])], 0);
    table.add_cols(&[Unit::Abs(a.margin[1])], usize::MAX);
    if let Some(bg) = a.background {
        let (nr, nc) = (table.nrow(), table.ncol());
        table.add(bg, 1, 1, nr, nc, Some(f64::NEG_INFINITY), true, "background");
    }
    table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtable::Bx;

    /// A plot table of patchwork's size: a null panel, a left axis column and a bottom axis row, a margin.
    fn plot_table(axis_w: f64, axis_h: f64) -> Patch {
        let mut widths = vec![Unit::ZERO; TABLE_COLS];
        let mut heights = vec![Unit::ZERO; TABLE_ROWS];
        widths[0] = Unit::Abs(5.0);
        widths[PANEL_COL - 2] = Unit::Abs(axis_w);
        widths[PANEL_COL - 1] = Unit::Null(1.0);
        heights[0] = Unit::Abs(5.0);
        heights[PANEL_ROW - 1] = Unit::Null(1.0);
        heights[PANEL_ROW] = Unit::Abs(axis_h);
        let mut gt = Gtable::new(widths, heights);
        gt.add(Grob::Zero, PANEL_ROW, PANEL_COL, PANEL_ROW, PANEL_COL, Some(1.0), false, "panel");
        gt.add(Grob::Zero, PANEL_ROW, PANEL_COL - 1, PANEL_ROW, PANEL_COL - 1, Some(3.0), false, "axis-l");
        gt.add(Grob::Zero, PANEL_ROW + 1, PANEL_COL, PANEL_ROW + 1, PANEL_COL, Some(3.0), false, "axis-b");
        Patch { gt, kind: Kind::Plot }
    }

    fn boxes(gt: &Gtable) -> Vec<(String, Bx)> {
        let mut out = Vec::new();
        gt.boxes(Bx::new(0.0, 0.0, 600.0, 400.0), "", &mut out);
        out
    }

    /// merge() (free() of a patchwork): one more level of tables named as R names them, every box of the patchwork
    /// where it was without it.
    #[test]
    fn merge_nests_once_more_and_moves_nothing() {
        let column = || nested(build(vec![plot_table(30.0, 10.0), plot_table(12.0, 20.0)], &Layout::column(2, vec![Some(1.0), Some(2.0)])));
        let free = |p: Patch| {
            let mut p = p;
            free_label(&mut p.gt, [true; 4]);
            p
        };
        let plain = build(vec![free(column()), plot_table(20.0, 10.0)], &Layout::row(2));
        let merged = build(vec![merge(free(column())), plot_table(20.0, 10.0)], &Layout::row(2));
        assert_eq!(plain.widths, merged.widths);
        assert_eq!(plain.heights, merged.heights);
        let (bp, bm) = (boxes(&plain), boxes(&merged));
        // the column's boxes one level deeper: patchwork-table-1::X becomes patchwork-table-1::patchwork-table-1::X
        for (name, b) in &bp {
            let want = match name.strip_prefix("patchwork-table-1::") {
                Some(rest) => format!("patchwork-table-1::patchwork-table-1::{rest}"),
                None => name.clone(),
            };
            assert!(bm.iter().any(|(n, x)| *n == want && x == b), "{name} {b:?} not found as {want}");
        }
        // the new level: its panel area, the nested patchwork's panel and table span the column
        let find = |n: &str| bm.iter().find(|(m, _)| m == n).map(|(_, b)| *b).unwrap();
        let col = find("patchwork-table-1::patchwork-table-1");
        assert_eq!(col, find("patchwork-table-1"));
        assert_eq!(find("patchwork-table-1::panel-area"), find("patchwork-table-1::panel-nested-patchwork-1"));
        assert_eq!(bm.len(), bp.len() + 3);
    }
}
