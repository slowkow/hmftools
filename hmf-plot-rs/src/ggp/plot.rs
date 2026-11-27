//! The plot table (ggplot_gtable): the facet's panel table, the axis titles (Facet$draw_labels), the five guide-box
//! cells (table_add_legends), the subtitle, title and caption rows, the tag padding and the plot margin, and the
//! plot background. The result has ggplot2's 18 x 15 structure around a single panel, which patchwork relies on.

use crate::gtable::{Grob, Gtable, Unit, sum_abs};

use super::facet::panel_extent;
use super::{Ctx, RectEl, TitleGrob, rect_grob};

/// The parts of a plot, already built.
pub struct PlotParts {
    /// the panel table from `facet`
    pub panels: Gtable,
    /// axis titles: top, bottom (x) and left, right (y); None for no title (a zero-sized row or column)
    pub xlab_t: Option<TitleGrob>,
    pub xlab_b: Option<TitleGrob>,
    pub ylab_l: Option<TitleGrob>,
    /// an axis title on the left that is not a text (a plotmath expression): its grob and width
    pub ylab_l_custom: Option<(Grob, f64)>,
    pub ylab_r: Option<TitleGrob>,
    /// guide boxes (from `guides::guide_box`), by position
    pub legend_right: Option<Gtable>,
    pub legend_left: Option<Gtable>,
    pub legend_top: Option<Gtable>,
    pub legend_bottom: Option<Gtable>,
    /// legend.box.spacing
    pub legend_box_spacing: f64,
    pub title: Option<TitleGrob>,
    pub subtitle: Option<TitleGrob>,
    pub caption: Option<TitleGrob>,
    /// plot.title.position "plot" (title over the whole width) rather than "panel"
    pub title_on_plot: bool,
    /// plot.margin (top, right, bottom, left)
    pub margin: [f64; 4],
    pub background: Option<RectEl>,
}

fn opt_grob(ctx: Ctx, t: &Option<TitleGrob>) -> Grob {
    t.clone().map(|t| t.grob(ctx)).unwrap_or(Grob::Zero)
}

/// The plot table.
pub fn plot_table(ctx: Ctx, p: PlotParts) -> Gtable {
    let mut gt = p.panels;

    // Facet$draw_labels
    let (_, l, _, r) = panel_extent(&gt);
    let h = |t: &Option<TitleGrob>| t.as_ref().map(|t| t.height(&ctx).abs()).unwrap_or(0.0);
    let w = |t: &Option<TitleGrob>| t.as_ref().map(|t| t.width(&ctx).abs()).unwrap_or(0.0);
    gt.add_rows(&[Unit::Abs(h(&p.xlab_t))], 0);
    gt.add(opt_grob(ctx, &p.xlab_t), 1, l, 1, r, None, false, "xlab-t");
    gt.add_rows(&[Unit::Abs(h(&p.xlab_b))], usize::MAX);
    let n = gt.nrow();
    gt.add(opt_grob(ctx, &p.xlab_b), n, l, n, r, None, false, "xlab-b");
    let (t, _, b, _) = panel_extent(&gt);
    match p.ylab_l_custom {
        Some((g, width)) => {
            gt.add_cols(&[Unit::Abs(width)], 0);
            gt.add(g, t, 1, b, 1, None, false, "ylab-l");
        }
        None => {
            gt.add_cols(&[Unit::Abs(w(&p.ylab_l))], 0);
            gt.add(opt_grob(ctx, &p.ylab_l), t, 1, b, 1, None, false, "ylab-l");
        }
    }
    gt.add_cols(&[Unit::Abs(w(&p.ylab_r))], usize::MAX);
    let n = gt.ncol();
    gt.add(opt_grob(ctx, &p.ylab_r), t, n, b, n, None, false, "ylab-r");

    // table_add_legends
    let size = |g: &Option<Gtable>, horizontal: bool| g.as_ref().map(|g| if horizontal { sum_abs(&g.widths) } else { sum_abs(&g.heights) }).unwrap_or(0.0);
    let spacing = |g: &Option<Gtable>| if g.is_some() { p.legend_box_spacing } else { 0.0 };
    let grob = |g: Option<Gtable>| g.map(Grob::table).unwrap_or(Grob::Zero);
    let (t, l, b, r) = panel_extent(&gt);
    let (wr, sr) = (size(&p.legend_right, true), spacing(&p.legend_right));
    gt.add_cols(&[Unit::Abs(sr)], usize::MAX);
    gt.add_cols(&[Unit::Abs(wr)], usize::MAX);
    let n = gt.ncol();
    gt.add(grob(p.legend_right), t, n, b, n, None, false, "guide-box-right");
    let (wl, sl) = (size(&p.legend_left, true), spacing(&p.legend_left));
    gt.add_cols(&[Unit::Abs(sl)], 0);
    gt.add_cols(&[Unit::Abs(wl)], 0);
    gt.add(grob(p.legend_left), t, 1, b, 1, None, false, "guide-box-left");
    let _ = (l, r);
    let (_, l, _, r) = panel_extent(&gt);
    let (hb, sb) = (size(&p.legend_bottom, false), spacing(&p.legend_bottom));
    gt.add_rows(&[Unit::Abs(sb)], usize::MAX);
    gt.add_rows(&[Unit::Abs(hb)], usize::MAX);
    let n = gt.nrow();
    gt.add(grob(p.legend_bottom), n, l, n, r, None, false, "guide-box-bottom");
    let (ht, st) = (size(&p.legend_top, false), spacing(&p.legend_top));
    gt.add_rows(&[Unit::Abs(st)], 0);
    gt.add_rows(&[Unit::Abs(ht)], 0);
    gt.add(Grob::Zero, 1, l, 1, r, None, false, "guide-box-top");
    if p.legend_top.is_some() {
        let last = gt.items.len() - 1;
        gt.items[last].grob = grob(p.legend_top);
    }
    let (t, l, b, r) = panel_extent(&gt);
    gt.add(Grob::Zero, t, l, b, r, None, false, "guide-box-inside");

    // title, subtitle, caption
    let (tl, tr) = if p.title_on_plot { (1, gt.ncol()) } else { (l, r) };
    gt.add_rows(&[Unit::Abs(h(&p.subtitle))], 0);
    gt.add(opt_grob(ctx, &p.subtitle), 1, tl, 1, tr, None, false, "subtitle");
    gt.add_rows(&[Unit::Abs(h(&p.title))], 0);
    gt.add(opt_grob(ctx, &p.title), 1, tl, 1, tr, None, false, "title");
    gt.add_rows(&[Unit::Abs(h(&p.caption))], usize::MAX);
    let n = gt.nrow();
    gt.add(opt_grob(ctx, &p.caption), n, l, n, r, None, false, "caption");

    // table_add_tag: zero padding
    gt.add_padding(0.0, 0.0, 0.0, 0.0);
    // plot margin
    let [mt, mr, mb, ml] = p.margin;
    gt.add_padding(mt, mr, mb, ml);
    if let Some(bg) = p.background {
        let (nr, nc) = (gt.nrow(), gt.ncol());
        gt.add(rect_grob(ctx, bg), 1, 1, nr, nc, Some(f64::NEG_INFINITY), false, "background");
        // ggplot2 moves the background to the front of the layout (drawn first among equal z)
        let bgi = gt.items.pop().unwrap();
        gt.items.insert(0, bgi);
    }
    gt
}
