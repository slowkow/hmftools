//! patchwork (1.3.2, which agrees with 1.2.0 for what is used here) on gtables: a ggplot's table made compliant
//! (add_strips, simplify_gt: the facet panels nested into one cell), plots combined in a design grid with their
//! panels, axes and titles aligned (build_patchwork: table_dims, set_grob_sizes, set_panel_dimensions), nested
//! patchworks, `free()` (free_label, free_panel), `wrap_elements()` and `plot_annotation()` (annotate_table).
//!
//! Every plot's table has patchwork's 18 rows x 15 columns, the panel at row 10, column 8.

use std::rc::Rc;

use crate::gtable::{Grob, Gtable, Placed, Pos, Unit, Vp, sum_abs};

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

/// The identity of a table of grobs (simplify_gt's row and column tables): its grobs' keys (zero grobs equal),
/// positions and sizes; None when a grob has no key.
fn table_key(t: &Gtable) -> Option<String> {
    let mut parts = Vec::new();
    for p in &t.items {
        let k = if p.grob.is_zero() { "zero".to_string() } else { p.key.clone()? };
        parts.push(format!("{}:{}:{}:{}:{}", p.t, p.l, p.b, p.r, k));
    }
    let dims: Vec<String> = t.widths.iter().chain(&t.heights).map(|u| format!("{:?}", u)).collect();
    Some(format!("table[{}|{}]", parts.join(","), dims.join(",")))
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
            gt_new.items.last_mut().unwrap().key = p.key.clone();
        }
        let b_mod = rows.1 - rows.0;
        for p in gt.items.iter().filter(|p| p.l == c && p.r == c && p.t > rows.1) {
            gt_new.add(p.grob.clone(), p.t - b_mod, c, p.b - b_mod, c, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
            gt_new.items.last_mut().unwrap().key = p.key.clone();
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
            let key = table_key(&table);
            gt_new.add(Grob::table(table), ii, cols.0, ii, cols.0, Some(z), false, &name);
            gt_new.items.last_mut().unwrap().key = key;
        }
    }
    if p_rows.len() == 1 {
        let rr = p_rows[0];
        for p in gt.items.iter().filter(|p| p.t == rr && p.b == rr && p.r < cols.0) {
            gt_new.add(p.grob.clone(), rr, p.l, rr, p.r, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
            gt_new.items.last_mut().unwrap().key = p.key.clone();
        }
        let r_mod = cols.1 - cols.0;
        for p in gt.items.iter().filter(|p| p.t == rr && p.b == rr && p.l > cols.1) {
            gt_new.add(p.grob.clone(), rr, p.l - r_mod, rr, p.r - r_mod, Some(p.z), p.clip, &p.name);
            gt_new.items.last_mut().unwrap().vp = p.vp;
            gt_new.items.last_mut().unwrap().key = p.key.clone();
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
            let key = table_key(&table);
            gt_new.add(Grob::table(table), rows.0, ii, rows.0, ii, Some(z), false, &name);
            gt_new.items.last_mut().unwrap().key = key;
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

/// wrap_elements(full = grob): patchwork's empty patch table (as `wrap_panel`) with the grob over the whole table,
/// margins included (clipped).
pub fn wrap_full(grob: Grob, margin: [f64; 4], background: Option<Grob>) -> Patch {
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
    gt.add(grob, 1, 1, TABLE_ROWS, TABLE_COLS, None, true, "full");
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
    let (mut gt_new, nrow, ncol) = build_inner(patches, layout);
    add_panel_area(&mut gt_new, nrow, ncol);
    gt_new
}

fn add_panel_area(gt_new: &mut Gtable, nrow: usize, ncol: usize) {
    let (b, r) = (PANEL_ROW + TABLE_ROWS * (nrow - 1), PANEL_COL + TABLE_COLS * (ncol - 1));
    gt_new.add(Grob::Zero, PANEL_ROW, PANEL_COL, b, r, Some(-1.0), false, "panel-area");
}

/// plot_layout(axes = , axis_titles = ): which axes and axis titles build_patchwork collects.
#[derive(Clone, Copy, Debug, Default)]
pub struct Collect {
    pub axes_x: bool,
    pub axes_y: bool,
    pub titles_x: bool,
    pub titles_y: bool,
}

/// build_patchwork with plot_layout(axes, axis_titles): after the table is built, identical axes in a run of
/// patches are kept only at the run's outer side (collect_axes) and identical axis titles are kept once, a y title
/// stretched over the rows of its run (collect_axis_titles, merge = TRUE). Grobs are compared by `Placed::key`.
pub fn build_collect(patches: Vec<Patch>, layout: &Layout, collect: Collect) -> Gtable {
    let (mut gt, nrow, ncol) = build_inner(patches, layout);
    if collect.axes_x {
        gt = collect_axes(gt, 'x');
    }
    if collect.axes_y {
        gt = collect_axes(gt, 'y');
    }
    if collect.titles_x {
        gt = collect_axis_titles(gt, 'x');
    }
    if collect.titles_y {
        gt = collect_axis_titles(gt, 'y');
    }
    add_panel_area(&mut gt, nrow, ncol);
    gt
}

fn build_inner(patches: Vec<Patch>, layout: &Layout) -> (Gtable, usize, usize) {
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
    (gt_new, nrow, ncol)
}

// ---- collect_axes and collect_axis_titles (patchwork 1.3.2) ----

/// A matrix of item indices (None: NA), column-major as R's.
#[derive(Clone, Debug)]
struct Mat {
    nr: usize,
    nc: usize,
    v: Vec<Option<usize>>,
}

impl Mat {
    fn at(&self, r: usize, c: usize) -> Option<usize> {
        self.v[c * self.nr + r]
    }
    fn set(&mut self, r: usize, c: usize, x: Option<usize>) {
        self.v[c * self.nr + r] = x;
    }
}

/// grob_layout: the items' cells on the grid of their distinct edges (each cell holds the item covering it)
fn grob_layout(gt: &Gtable, idx: &[usize]) -> Mat {
    let mut tops: Vec<usize> = idx.iter().flat_map(|&i| [gt.items[i].t, gt.items[i].b]).collect();
    tops.sort();
    tops.dedup();
    let mut lefts: Vec<usize> = idx.iter().flat_map(|&i| [gt.items[i].l, gt.items[i].r]).collect();
    lefts.sort();
    lefts.dedup();
    let mut m = Mat { nr: tops.len(), nc: lefts.len(), v: vec![None; tops.len() * lefts.len()] };
    for &i in idx {
        let p = &gt.items[i];
        let (t, b) = (tops.iter().position(|&x| x == p.t).unwrap(), tops.iter().position(|&x| x == p.b).unwrap());
        let (l, r) = (lefts.iter().position(|&x| x == p.l).unwrap(), lefts.iter().position(|&x| x == p.r).unwrap());
        for rr in t..=b {
            for cc in l..=r {
                m.set(rr, cc, Some(i));
            }
        }
    }
    m
}

/// hash(unname_grob(x)) as a key: the item's key (with unpack, a one-grob table's grob's), or one of its own
fn grob_hash(gt: &Gtable, i: usize, unpack: bool) -> String {
    let p = &gt.items[i];
    if unpack {
        if let Grob::Table(t) = &p.grob {
            if t.items.len() == 1 {
                let inner = &t.items[0];
                if inner.grob.is_zero() {
                    return "zero".into();
                }
                // the inner grob's key, else the table's own (keys are given to the titles as placed)
                if let Some(k) = inner.key.clone().or_else(|| p.key.clone()) {
                    return k;
                }
                return format!("item-{i}");
            }
        }
    }
    if p.grob.is_zero() {
        return "zero".into();
    }
    p.key.clone().unwrap_or_else(|| format!("item-{i}"))
}

/// grob_id: each cell's grob as a number (first-seen order of the hashes); without merge, the hash includes the
/// span of the grob along the rows (byrow: its columns) so that only grobs in the same line are equal
fn grob_id(gt: &Gtable, layout: &Mat, byrow: bool, merge: bool, unpack: bool) -> Vec<Option<usize>> {
    let mut hashes: Vec<Option<String>> = vec![None; layout.v.len()];
    for c in 0..layout.nc {
        for r in 0..layout.nr {
            let Some(i) = layout.at(r, c) else { continue };
            let mut h = grob_hash(gt, i, unpack);
            if !merge {
                // ave(index, layout, FUN = min / max) over the cells of this grob
                let mut cells = Vec::new();
                for cc in 0..layout.nc {
                    for rr in 0..layout.nr {
                        if layout.at(rr, cc) == Some(i) {
                            cells.push(if byrow { cc + 1 } else { rr + 1 });
                        }
                    }
                }
                h = format!("{h}{};{}", cells.iter().min().unwrap(), cells.iter().max().unwrap());
            }
            hashes[c * layout.nr + r] = Some(h);
        }
    }
    let mut uniq: Vec<String> = Vec::new();
    hashes
        .into_iter()
        .map(|h| {
            h.map(|h| match uniq.iter().position(|u| *u == h) {
                Some(p) => p + 1,
                None => {
                    uniq.push(h);
                    uniq.len()
                }
            })
        })
        .collect()
}

/// A run of rle_2d (0-based, inclusive), with its value (None: NA)
#[derive(Clone, Debug)]
struct Run {
    col_start: usize,
    col_end: usize,
    row_start: usize,
    row_end: usize,
    value: Option<usize>,
}

/// rle_2d: runs down the columns (byrow: along the rows), merged with the next column's run of the same rows and
/// value
fn rle_2d(m: &[Option<usize>], nr: usize, nc: usize, byrow: bool, ignore_na: bool) -> Vec<Run> {
    let n = m.len();
    if n == 0 {
        return Vec::new();
    }
    // byrow: transpose (rows become columns), then swap the labels back
    let (m, nr, nc): (Vec<Option<usize>>, usize, usize) = if byrow {
        let mut t = vec![None; n];
        for c in 0..nc {
            for r in 0..nr {
                t[r * nc + c] = m[c * nr + r];
            }
        }
        (t, nc, nr)
    } else {
        (m.to_vec(), nr, nc)
    };
    let rename = |runs: Vec<Run>| -> Vec<Run> {
        if byrow {
            runs.into_iter().map(|r| Run { col_start: r.row_start, col_end: r.row_end, row_start: r.col_start, row_end: r.col_end, value: r.value }).collect()
        } else {
            runs
        }
    };
    let mut levels: Vec<Option<usize>> = Vec::new();
    for x in &m {
        if !levels.contains(x) {
            levels.push(*x);
        }
    }
    let non_na = levels.iter().filter(|x| x.is_some()).count();
    if (ignore_na && non_na == 1) || levels.len() == 1 {
        // sort(levels, na.last = TRUE)[1]
        let mut l: Vec<usize> = levels.iter().filter_map(|x| *x).collect();
        l.sort();
        let value = l.first().copied();
        return rename(vec![Run { col_start: 0, col_end: nc - 1, row_start: 0, row_end: nr - 1, value }]);
    }
    if levels.len() == n {
        let mut out = Vec::new();
        for c in 0..nc {
            for r in 0..nr {
                out.push(Run { col_start: c, col_end: c, row_start: r, row_end: r, value: m[c * nr + r] });
            }
        }
        return rename(out);
    }
    let code: Vec<usize> = m.iter().map(|x| levels.iter().position(|l| l == x).unwrap()).collect();
    // ends: where the value changes or a column ends (column-major), and the last cell
    let mut ends = Vec::new();
    for k in 0..n - 1 {
        if code[k + 1] != code[k] || k % nr == nr - 1 {
            ends.push(k);
        }
    }
    ends.push(n - 1);
    let mut starts = Vec::new();
    let mut prev = 0usize;
    for (j, &e) in ends.iter().enumerate() {
        starts.push(if j == 0 { 0 } else { prev + 1 });
        prev = e;
    }
    if nc == 1 {
        return rename(ends.iter().zip(&starts).map(|(&e, &s)| Run { col_start: 0, col_end: 0, row_start: s, row_end: e, value: levels[code[e]] }).collect());
    }
    let mut runs: Vec<Option<Run>> = ends
        .iter()
        .zip(&starts)
        .map(|(&e, &s)| Some(Run { col_start: e / nr, col_end: e / nr, row_start: s % nr, row_end: s % nr + (e - s), value: levels[code[e]] }))
        .collect();
    // runs of the same rows and value in adjacent columns merged into the first
    let mut seen: Vec<((usize, usize, usize), usize)> = Vec::new();
    for i in 0..runs.len() {
        let r = runs[i].clone().unwrap();
        let key = (r.row_start, r.row_end, code[ends[i]]);
        match seen.iter().find(|(k, _)| *k == key).map(|(_, h)| *h) {
            Some(h) if r.col_start == runs[h].as_ref().unwrap().col_end + 1 => {
                runs[h].as_mut().unwrap().col_end = r.col_start;
                runs[i] = None;
            }
            _ => {
                seen.retain(|(k, _)| *k != key);
                seen.push((key, i));
            }
        }
    }
    rename(runs.into_iter().flatten().collect())
}

fn items_named(gt: &Gtable, prefix: &str) -> Vec<usize> {
    gt.items.iter().enumerate().filter(|(_, p)| p.name.starts_with(prefix)).map(|(i, _)| i).collect()
}

fn nested_patch_items(gt: &Gtable) -> Vec<usize> {
    gt.items.iter().enumerate().filter(|(_, p)| p.name.contains("panel-nested-patchwork")).map(|(i, _)| i).collect()
}

/// delete_grobs(gt, idx, resize = TRUE): the grobs removed; their first rows and columns that no remaining
/// (non-zero) grob starts or ends in are given size 0
fn delete_grobs(mut gt: Gtable, idx: &[usize]) -> Gtable {
    if idx.is_empty() {
        return gt;
    }
    let rows: Vec<usize> = dedup(idx.iter().map(|&i| gt.items[i].t).collect());
    let cols: Vec<usize> = dedup(idx.iter().map(|&i| gt.items[i].l).collect());
    let keep: Vec<Placed> = gt.items.iter().enumerate().filter(|(i, _)| !idx.contains(i)).map(|(_, p)| p.clone()).collect();
    gt.items = keep;
    let used_rows: Vec<usize> = gt.items.iter().filter(|p| !p.grob.is_zero()).flat_map(|p| [p.t, p.b]).collect();
    let used_cols: Vec<usize> = gt.items.iter().filter(|p| !p.grob.is_zero()).flat_map(|p| [p.l, p.r]).collect();
    for r in rows.into_iter().filter(|r| !used_rows.contains(r)) {
        gt.heights[r - 1] = Unit::ZERO;
    }
    for c in cols.into_iter().filter(|c| !used_cols.contains(c)) {
        gt.widths[c - 1] = Unit::ZERO;
    }
    gt
}

/// a grob's size along rows (heights) or columns (widths) as grid's max_height/max_width give it for the grobs
/// here: a table's total; other grobs keep the current size of the line
fn grob_extent(p: &Placed, rows: bool, current: f64) -> f64 {
    match &p.grob {
        Grob::Table(t) => sum_abs(if rows { &t.heights } else { &t.widths }),
        _ => current,
    }
}

/// retrofit_rows / retrofit_cols: lines (of those given) whose remaining non-zero grobs all match `pattern`
/// (here: a prefix) take the largest size of those grobs
fn retrofit(mut gt: Gtable, lines: &[usize], prefix: &str, rows: bool) -> Gtable {
    if lines.is_empty() {
        return gt;
    }
    let first = |p: &Placed| if rows { p.t } else { p.l };
    let last = |p: &Placed| if rows { p.b } else { p.r };
    let mut groups: Vec<(usize, bool)> = Vec::new();
    for p in gt.items.iter().filter(|p| !p.grob.is_zero() && (lines.contains(&first(p)) || lines.contains(&last(p)))) {
        let all_pattern = p.name.starts_with(prefix);
        match groups.iter_mut().find(|(l, _)| *l == first(p)) {
            Some(g) => g.1 &= all_pattern,
            None => groups.push((first(p), all_pattern)),
        }
    }
    for (line, ok) in groups {
        if !ok {
            continue;
        }
        let current = if rows { gt.heights[line - 1].abs() } else { gt.widths[line - 1].abs() };
        let size = gt.items.iter().filter(|p| !p.grob.is_zero() && (first(p) == line || last(p) == line)).map(|p| grob_extent(p, rows, current)).fold(0.0, f64::max);
        if rows {
            gt.heights[line - 1] = Unit::Abs(size);
        } else {
            gt.widths[line - 1] = Unit::Abs(size);
        }
    }
    gt
}

/// collect_axes(gt, dir): of identical axes in a run of patches, only the outermost is kept (bottom for axis-b,
/// top for axis-t, left for axis-l, right for axis-r)
fn collect_axes(gt: Gtable, dir: char) -> Gtable {
    let names = if dir == 'x' { ["axis-b", "axis-t"] } else { ["axis-l", "axis-r"] };
    let mut delete: Vec<usize> = Vec::new();
    for name in names {
        let idx = items_named(&gt, name);
        if idx.len() < 2 || idx.iter().all(|&i| gt.items[i].grob.is_zero()) {
            continue;
        }
        let patch_index = nested_patch_items(&gt);
        let all: Vec<usize> = idx.iter().chain(&patch_index).copied().collect();
        let mut layout = grob_layout(&gt, &all);
        for x in layout.v.iter_mut() {
            if x.is_some_and(|i| patch_index.contains(&i)) {
                *x = None;
            }
        }
        let structure = grob_id(&gt, &layout, dir == 'x', false, false);
        let ids: Vec<usize> = structure.iter().filter_map(|x| *x).collect();
        if !has_duplicates(&ids) {
            continue;
        }
        let runs: Vec<Run> = rle_2d(&structure, layout.nr, layout.nc, dir == 'y', false).into_iter().filter(|r| r.value.is_some()).collect();
        let start: Vec<usize> = runs
            .iter()
            .filter_map(|r| {
                let row = if name == "axis-b" { r.row_end } else { r.row_start };
                let col = if name == "axis-r" { r.col_end } else { r.col_start };
                layout.at(row, col)
            })
            .collect();
        delete.extend(idx.iter().filter(|i| !start.contains(i)));
    }
    let deleted_rows: Vec<usize> = dedup(delete.iter().flat_map(|&i| [gt.items[i].t, gt.items[i].b]).collect());
    let deleted_cols: Vec<usize> = dedup(delete.iter().flat_map(|&i| [gt.items[i].l, gt.items[i].r]).collect());
    let new = delete_grobs(gt, &delete);
    let new = retrofit(new, &deleted_rows, "axis", true);
    retrofit(new, &deleted_cols, "axis", false)
}

fn has_duplicates(v: &[usize]) -> bool {
    let mut s = v.to_vec();
    s.sort();
    s.windows(2).any(|w| w[0] == w[1])
}

fn dedup(v: Vec<usize>) -> Vec<usize> {
    let mut out = Vec::new();
    for x in v {
        if !out.contains(&x) {
            out.push(x);
        }
    }
    out
}

/// collect_axis_titles(gt, dir, merge = TRUE): of identical axis titles in a run of patches one is kept (the
/// bottom x title, the top y title...); a y title is stretched over the rows of its run, an x title over its
/// columns
fn collect_axis_titles(mut gt: Gtable, dir: char) -> Gtable {
    let names = if dir == 'x' { ["xlab-t", "xlab-b"] } else { ["ylab-l", "ylab-r"] };
    let mut delete: Vec<usize> = Vec::new();
    for name in names {
        let idx = items_named(&gt, name);
        if idx.len() < 2 || idx.iter().all(|&i| gt.items[i].grob.is_zero()) {
            continue;
        }
        let patch_index = nested_patch_items(&gt);
        let all: Vec<usize> = idx.iter().chain(&patch_index).copied().collect();
        let mut layout = grob_layout(&gt, &all);
        let nested: Vec<bool> = layout.v.iter().map(|x| x.is_some_and(|i| patch_index.contains(&i))).collect();
        for (x, n) in layout.v.iter_mut().zip(&nested) {
            if *n {
                *x = None;
            }
        }
        let mut structure = grob_id(&gt, &layout, dir == 'x', true, true);
        let ids: Vec<usize> = structure.iter().filter_map(|x| *x).collect();
        if !has_duplicates(&ids) {
            continue;
        }
        for (s, n) in structure.iter_mut().zip(&nested) {
            if *n {
                *s = Some(0);
            }
        }
        let runs: Vec<Run> = rle_2d(&structure, layout.nr, layout.nc, dir == 'y', true).into_iter().filter(|r| r.value.is_some_and(|v| v != 0)).collect();
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for r in &runs {
            let first: Vec<Option<usize>> = match name {
                "xlab-t" => (r.col_start..=r.col_end).map(|c| layout.at(r.row_start, c)).collect(),
                "xlab-b" => (r.col_start..=r.col_end).map(|c| layout.at(r.row_end, c)).collect(),
                "ylab-l" => (r.row_start..=r.row_end).map(|rr| layout.at(rr, r.col_start)).collect(),
                _ => (r.row_start..=r.row_end).map(|rr| layout.at(rr, r.col_end)).collect(),
            };
            let first = first.into_iter().flatten().next();
            let mut members: Vec<usize> = first.into_iter().collect();
            for c in r.col_start..=r.col_end {
                for rr in r.row_start..=r.row_end {
                    if let Some(i) = layout.at(rr, c) {
                        if !members.contains(&i) {
                            members.push(i);
                        }
                    }
                }
            }
            groups.push(members);
        }
        let titles: Vec<usize> = groups.iter().filter_map(|g| g.first().copied()).collect();
        delete.extend(idx.iter().filter(|i| !titles.contains(i)));
        if (dir == 'x' && runs.iter().all(|r| r.col_start == r.col_end)) || (dir == 'y' && runs.iter().all(|r| r.row_start == r.row_end)) {
            continue;
        }
        let zmax = idx.iter().map(|&i| gt.items[i].z).fold(f64::NEG_INFINITY, f64::max);
        for g in &groups {
            let Some(&ti) = g.first() else { continue };
            if dir == 'y' {
                let t = g.iter().map(|&i| gt.items[i].t).min().unwrap();
                let b = g.iter().map(|&i| gt.items[i].b).max().unwrap();
                (gt.items[ti].t, gt.items[ti].b) = (t, b);
            } else {
                let l = g.iter().map(|&i| gt.items[i].l).min().unwrap();
                let r = g.iter().map(|&i| gt.items[i].r).max().unwrap();
                (gt.items[ti].l, gt.items[ti].r) = (l, r);
            }
            gt.items[ti].z = zmax;
        }
    }
    delete_grobs(gt, &dedup(delete))
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

    /// patchwork:::rle_2d on R's example (4 x 3, column-major 1 1 2 2 / 1 1 2 3 / 1 NA 2 3), by column and by row
    #[test]
    fn rle_2d_as_r() {
        let m: Vec<Option<usize>> = [1, 1, 2, 2, 1, 1, 2, 3, 1, 0, 2, 3].iter().map(|&v| if v == 0 { None } else { Some(v) }).collect();
        let runs = |byrow| rle_2d(&m, 4, 3, byrow, false).into_iter().map(|r| (r.col_start + 1, r.col_end + 1, r.row_start + 1, r.row_end + 1, r.value)).collect::<Vec<_>>();
        assert_eq!(runs(false), vec![(1, 2, 1, 2, Some(1)), (1, 1, 3, 4, Some(2)), (2, 3, 3, 3, Some(2)), (2, 3, 4, 4, Some(3)), (3, 3, 1, 1, Some(1)), (3, 3, 2, 2, None)]);
        assert_eq!(runs(true), vec![(1, 3, 1, 1, Some(1)), (1, 2, 2, 2, Some(1)), (3, 3, 2, 2, None), (1, 3, 3, 3, Some(2)), (1, 1, 4, 4, Some(2)), (2, 3, 4, 4, Some(3))]);
        let one = [Some(5), Some(5)];
        let r = rle_2d(&one, 2, 1, false, false);
        assert_eq!((r.len(), r[0].row_start, r[0].row_end), (1, 0, 1));
        let r = rle_2d(&one, 2, 1, true, false);
        assert_eq!((r.len(), r[0].row_start, r[0].row_end, r[0].col_start, r[0].col_end), (1, 0, 1, 0, 0));
    }

    /// plot_layout(axes = "collect_x", axis_titles = "collect") on two plots in a column: the top plot's x axis (equal
    /// to the bottom's) deleted and its row emptied; equal y titles kept once, the top one stretched over both
    /// plots' rows; different ones kept.
    #[test]
    fn collect_axes_and_titles_in_a_column() {
        let titled = |axis_key: &str, ylab: &str| {
            let mut p = plot_table(30.0, 10.0);
            p.gt.add(Grob::draw_fn(|_, _| {}), PANEL_ROW, 3, PANEL_ROW, 3, None, false, "ylab-l");
            p.gt.widths[2] = Unit::Abs(12.0);
            p.gt.items[2].grob = Grob::draw_fn(|_, _| {});
            p.gt.set_keys("axis-b", axis_key);
            p.gt.set_keys("ylab-l", ylab);
            p
        };
        let all = Collect { axes_x: true, axes_y: false, titles_x: true, titles_y: true };
        let gt = build_collect(vec![titled("x", "y"), titled("x", "y")], &Layout::column(2, vec![None, None]), all);
        let axes: Vec<&Placed> = gt.items.iter().filter(|p| p.name.starts_with("axis-b")).collect();
        assert_eq!(axes.len(), 1);
        assert_eq!(axes[0].name, "axis-b-2");
        assert_eq!(gt.heights[PANEL_ROW].abs(), 0.0);
        assert_eq!(gt.heights[TABLE_ROWS + PANEL_ROW].abs(), 10.0);
        let ylabs: Vec<&Placed> = gt.items.iter().filter(|p| p.name.starts_with("ylab-l")).collect();
        assert_eq!(ylabs.len(), 1);
        assert_eq!((ylabs[0].name.as_str(), ylabs[0].t, ylabs[0].b), ("ylab-l-1", PANEL_ROW, TABLE_ROWS + PANEL_ROW));
        // different axes and titles stay
        let gt = build_collect(vec![titled("x1", "y1"), titled("x2", "y2")], &Layout::column(2, vec![None, None]), all);
        assert_eq!(gt.items.iter().filter(|p| p.name.starts_with("axis-b")).count(), 2);
        assert_eq!(gt.items.iter().filter(|p| p.name.starts_with("ylab-l")).count(), 2);
        // without collecting, build_collect is build
        let a = build(vec![titled("x", "y"), titled("x", "y")], &Layout::column(2, vec![None, None]));
        let b = build_collect(vec![titled("x", "y"), titled("x", "y")], &Layout::column(2, vec![None, None]), Collect::default());
        assert_eq!((a.widths, a.heights, a.items.len()), (b.widths, b.heights, b.items.len()));
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
