//! A ggplot2 4.0.3 plot as `plot_qc.R` builds them, on hmf-plot-rs's tables: one panel or a facet grid, position
//! scales already resolved (`Pos`), the theme of the script (theme_bw(10.5) with its changes), titles, axis titles
//! (also the secondary axis's) and the plot margin. The plots' layers are drawn by closures on each panel.

use std::rc::Rc;

use hmf_plot::color::{Rgba, r_colour};
use hmf_plot::ggp::facet::{self, Grid, Strip};
use hmf_plot::ggp::guides::{self, AxisSpec, Side};
use hmf_plot::ggp::plot::{PlotParts, plot_table};
use hmf_plot::ggp::{Ctx, LineEl, PT, RectEl, TextEl, TitleGrob, draw_text, rect};
use hmf_plot::gtable::{Bx, Canvas, Grob, Gtable, Unit};
use hmf_plot::rdev::Face;

pub const BASE: f64 = 10.5;
pub const HALF: f64 = BASE / 2.0;

pub fn col(s: &str) -> Rgba {
    r_colour(s).unwrap_or_else(|| panic!("colour {s}"))
}

/// R's col_mix(ink, paper, amount) for black ink and white paper: a grey of amount * 255
pub fn mix(amount: f64) -> Rgba {
    let v = (amount * 255.0).round() as u8;
    Rgba::rgb8(v, v, v)
}

/// base_line_size and base_rect_size of theme_bw(10.5): 10.5 / 22 mm
pub const LINE: f64 = BASE / 22.0;

/// A position scale as a panel draws it: the (transformed) range, the breaks (transformed) and labels, and the
/// ticks of guide_axis_logticks (transformed position, type 1 long, 2 mid, 3 short).
#[derive(Clone, Debug)]
pub struct Pos {
    pub range: (f64, f64),
    pub log: bool,
    pub breaks: Vec<f64>,
    pub labels: Vec<String>,
    pub logticks: Vec<(f64, u8)>,
}

impl Pos {
    pub fn t(&self, v: f64) -> f64 {
        if self.log { v.log10() } else { v }
    }
    /// npc of a data value; infinite values at the panel's edges (ggplot2 squishes them)
    pub fn npc(&self, v: f64) -> f64 {
        let t = self.t(v);
        if t == f64::INFINITY {
            return 1.0;
        }
        if t == f64::NEG_INFINITY {
            return 0.0;
        }
        (t - self.range.0) / (self.range.1 - self.range.0)
    }
    pub fn npc_t(&self, t: f64) -> f64 {
        (t - self.range.0) / (self.range.1 - self.range.0)
    }
}

/// The mapping of one panel: device box, horizontal and vertical scales; with `flip` the data's x is vertical.
#[derive(Clone)]
pub struct Map {
    pub b: Bx,
    pub h: Pos,
    pub v: Pos,
    pub flip: bool,
}

impl Map {
    /// The device point of data (x, y).
    pub fn at(&self, x: f64, y: f64) -> (f64, f64) {
        let (hx, vy) = if self.flip { (y, x) } else { (x, y) };
        (self.b.npc_x(self.h.npc(hx)), self.b.npc_y(self.v.npc(vy)))
    }
}

pub type Layers = Rc<dyn Fn(&mut Canvas, &Ctx, &Map)>;

/// An axis title: text, or plotmath's atop(a, atop(b)) (a over a smaller b) as the MS indel bias plot's.
#[derive(Clone, Debug)]
pub enum Label {
    Text(String),
    /// the MS indel bias plot's y title before upstream 883ab11bba
    #[allow(dead_code)]
    Atop(String, String),
}

/// What a plot shows besides its layers.
pub struct Spec {
    pub title: Option<String>,
    /// plot.subtitle (upstream 883ab11bba's BQR and MS indel plots): the text and its size in pt
    pub subtitle: Option<(String, f64)>,
    pub xlab: Option<Label>,
    pub ylab: Option<Label>,
    /// the secondary (duplicated) y axis's title, on the right
    pub ylab_r: Option<String>,
    /// the scales of the horizontal axis per panel column, of the vertical axis per panel row
    pub h: Vec<Pos>,
    pub v: Vec<Pos>,
    pub flip: bool,
    /// facet_grid strips: column labels (top) and row labels (right)
    pub strips_top: Option<Vec<String>>,
    pub strips_right: Option<Vec<String>>,
    pub layers: Vec<Vec<Layers>>,
    /// axis text elements (None: blank) and whether ticks are drawn
    pub h_text: Option<TextEl>,
    pub v_text: Option<TextEl>,
    pub h_ticks: bool,
    pub v_ticks: bool,
    /// a secondary y axis (its text and ticks blank, its tick length kept)
    pub right_axis: bool,
    /// major grid lines (grey90, 0.21 mm) across the horizontal axis's breaks (vertical lines) / the vertical's
    pub grid_h: bool,
    pub grid_v: bool,
    pub spacing_x: f64,
    pub spacing_y: f64,
    pub margin: [f64; 4],
    pub title_blank: bool,
}

pub struct Theme {
    pub ctx: Ctx,
}

impl Theme {
    pub fn pt(&self, v: f64) -> f64 {
        self.ctx.pt(v)
    }
    pub fn title_el(&self) -> TextEl {
        let mut t = TextEl::new(BASE);
        t.hjust = 0.0;
        t.vjust = 1.0;
        t.margin = [0.0, 0.0, self.pt(HALF), 0.0];
        t
    }
    pub fn axis_text_h(&self) -> TextEl {
        let mut t = TextEl::new(0.8 * BASE);
        t.colour = mix(0.302);
        t.vjust = 1.0;
        t.margin = [self.pt(0.8 * HALF / 2.0), 0.0, 0.0, 0.0];
        t
    }
    pub fn axis_text_v(&self) -> TextEl {
        let mut t = TextEl::new(0.8 * BASE);
        t.colour = mix(0.302);
        t.hjust = 1.0;
        t.margin = [0.0, self.pt(0.8 * HALF / 2.0), 0.0, 0.0];
        t
    }
    pub fn ticks(&self) -> LineEl {
        LineEl::solid(mix(0.2), LINE)
    }
    pub fn tick_len(&self) -> f64 {
        self.pt(HALF / 2.0)
    }
    pub fn margin(&self) -> [f64; 4] {
        let m = self.pt(HALF);
        [m, m, m, m]
    }
}

/// The axis title grobs as ggplot2 renders them (axis.title.x: margin top, vjust 1; .y: angle 90, margin right,
/// vjust 1; .y.right: angle -90, margin left, vjust 1).
fn axis_title(th: &Theme, l: &Label, side: Side) -> TitleGrob {
    let mut el = TextEl::new(BASE);
    let q = th.pt(HALF / 2.0);
    match side {
        Side::Bottom => {
            el.vjust = 1.0;
            el.margin = [q, 0.0, 0.0, 0.0];
        }
        Side::Left => {
            el.angle = 90.0;
            el.vjust = 1.0;
            el.margin = [0.0, q, 0.0, 0.0];
        }
        Side::Right => {
            el.angle = -90.0;
            el.vjust = 1.0;
            el.margin = [0.0, 0.0, 0.0, q];
        }
        Side::Top => {}
    }
    match l {
        Label::Text(s) => TitleGrob::new(s, &el, side.vertical(), !side.vertical()),
        Label::Atop(a, _) => TitleGrob::new(a, &el, side.vertical(), !side.vertical()),
    }
}

/// A plot's table (ggplot_gtable) from its spec.
pub fn table(th: &Theme, spec: Spec) -> Gtable {
    let ctx = th.ctx;
    let ncol = spec.h.len();
    let nrow = spec.v.len();
    let grey90 = col("grey90");
    let border = mix(0.2);
    // panels: background, grid lines, layers, border (clipped to the panel)
    let mut panels: Vec<Vec<Grob>> = Vec::new();
    for r in 0..nrow {
        let mut row = Vec::new();
        for c in 0..ncol {
            let (h, v) = (spec.h[c].clone(), spec.v[r].clone());
            let layers = spec.layers[r][c].clone();
            let (gh, gv, flip) = (spec.grid_h, spec.grid_v, spec.flip);
            row.push(Grob::draw_fn(move |cv, b| {
                rect(cv, &ctx, b.x0, b.y0, b.x1, b.y1, Some(Rgba::WHITE), None, Vec::new());
                let s = hmf_plot::ggp::stroke(&ctx, grey90, 0.21 * PT, Vec::new(), hmf_plot::LineCap::Butt);
                if gh {
                    for &t in &h.breaks {
                        let x = b.npc_x(h.npc_t(t));
                        cv.fig.line(x, b.y0, x, b.y1, s.clone());
                    }
                }
                if gv {
                    for &t in &v.breaks {
                        let y = b.npc_y(v.npc_t(t));
                        cv.fig.line(b.x0, y, b.x1, y, s.clone());
                    }
                }
                let m = Map { b, h: h.clone(), v: v.clone(), flip };
                layers(cv, &ctx, &m);
                rect(cv, &ctx, b.x0, b.y0, b.x1, b.y1, None, Some((border, LINE * PT)), Vec::new());
            }));
        }
        panels.push(row);
    }
    let axis = |pos: &Pos, side: Side, text: &Option<TextEl>, ticks: bool| -> guides::Axis {
        let breaks: Vec<f64> = pos.breaks.iter().map(|&t| pos.npc_t(t)).collect();
        let log = !pos.logticks.is_empty();
        // guide_axis_logticks: long 2.25, mid 1.5, short 0.75 times the tick length, none at the breaks; the
        // axis is as deep as its longest tick (drawn or not)
        let extra: Vec<(f64, f64)> = pos
            .logticks
            .iter()
            .map(|&(t, k)| (pos.npc_t(t), th.tick_len() * match k {
                1 => 2.25,
                2 => 1.5,
                _ => 0.75,
            }))
            .collect();
        guides::axis(
            ctx,
            AxisSpec {
                side,
                breaks,
                labels: pos.labels.clone(),
                text: text.clone(),
                ticks: ticks.then(|| th.ticks()),
                // ggplot2 4.0.3: blank ticks have no length (3.5.1 kept it)
                tick_len: if log || !ticks { 0.0 } else { th.tick_len() },
                major_ticks: !log,
                extra_ticks: if ticks { extra } else { Vec::new() },
            },
        )
    };
    let axes_bottom: Vec<guides::Axis> = spec.h.iter().map(|p| axis(p, Side::Bottom, &spec.h_text, spec.h_ticks)).collect();
    let axes_left: Vec<guides::Axis> = spec.v.iter().map(|p| axis(p, Side::Left, &spec.v_text, spec.v_ticks)).collect();
    let axes_right: Vec<guides::Axis> = spec.v.iter().map(|p| if spec.right_axis { axis(p, Side::Right, &None, false) } else { guides::Axis::zero() }).collect();
    let axes_top: Vec<guides::Axis> = spec.h.iter().map(|_| guides::Axis::zero()).collect();

    let mut strip_el = TextEl::new(0.8 * BASE);
    strip_el.colour = mix(0.1);
    let m = th.pt(0.8 * HALF);
    strip_el.margin = [m, m, m, m];
    let strip_bg = RectEl { fill: Some(mix(0.851)), colour: Some(mix(0.2)), linewidth: LINE };
    let strips_x: Option<Vec<Strip>> = spec.strips_top.as_ref().map(|labels| labels.iter().map(|l| facet::strip(ctx, l, Some(&strip_el), Some(strip_bg), false)).collect());
    let mut strip_y = strip_el.clone();
    strip_y.angle = -90.0;
    let strips_y: Option<Vec<Strip>> = spec.strips_right.as_ref().map(|labels| labels.iter().map(|l| facet::strip(ctx, l, Some(&strip_y), Some(strip_bg), true)).collect());

    let panel_table = if spec.strips_top.is_none() && spec.strips_right.is_none() && nrow == 1 && ncol == 1 {
        let mut a_b = axes_bottom;
        let mut a_l = axes_left;
        let mut a_r = axes_right;
        let mut a_t = axes_top;
        facet::facet_null(panels.remove(0).remove(0), a_t.remove(0), a_r.remove(0), a_b.remove(0), a_l.remove(0), true)
    } else {
        Grid {
            widths: vec![1.0; ncol],
            heights: vec![1.0; nrow],
            panels,
            spacing_x: spec.spacing_x,
            spacing_y: spec.spacing_y,
            clip: true,
            axes_top,
            axes_bottom,
            axes_left,
            axes_right,
            strips_x,
            strips_y,
            switch_x: false,
            switch_y: false,
            outside_x: false,
            outside_y: false,
            padding: th.pt(HALF / 2.0),
            pad_always: false,
        }
        .build()
    };
    let (ylab_l, ylab_l_custom) = match spec.ylab.as_ref() {
        Some(Label::Atop(a, b)) => (None, Some(atop_title(th, a, b))),
        Some(l) => (Some(axis_title(th, l, Side::Left)), None),
        None => (None, None),
    };
    let title = spec.title.as_ref().map(|t| TitleGrob::new(t, &th.title_el(), true, true));
    let subtitle = spec.subtitle.as_ref().map(|(t, size)| {
        let mut el = th.title_el();
        el.size = *size;
        TitleGrob::new(t, &el, true, true)
    });
    let parts = PlotParts {
        panels: panel_table,
        xlab_t: None,
        xlab_b: spec.xlab.as_ref().map(|l| axis_title(th, l, Side::Bottom)),
        ylab_l,
        ylab_l_custom,
        ylab_r: spec.ylab_r.as_ref().map(|t| axis_title(th, &Label::Text(t.clone()), Side::Right)),
        legend_right: None,
        legend_left: None,
        legend_top: None,
        legend_bottom: None,
        legend_box_spacing: th.pt(2.0 * HALF),
        title: if spec.title_blank { None } else { title },
        subtitle,
        caption: None,
        title_on_plot: false,
        margin: spec.margin,
        background: Some(RectEl { fill: Some(Rgba::WHITE), colour: Some(Rgba::WHITE), linewidth: LINE }),
    };
    plot_table(ctx, parts)
}

/// The MS indel bias plot's y title, plotmath's atop("...", atop("...")) as ggplot2's titleGrob places it at angle 90
/// (x at the cell's left, y at its middle, hjust 0.5 and vjust 1 for GEMathText) with the descent of the font and
/// the right margin added to its width.
fn atop_title(th: &Theme, a: &str, b: &str) -> (Grob, f64) {
    let ctx = th.ctx;
    let l = hmf_plot::ggp::plotmath::atop_atop(&ctx.dev, a, b, BASE);
    let descent = ctx.dev.font_descent(BASE, Face::Plain);
    let width = l.bbox.h + l.bbox.d + descent + th.pt(HALF / 2.0);
    let (a, b) = (a.to_string(), b.to_string());
    let grob = Grob::draw_fn(move |c, bx| {
        let (xref, yref) = (bx.x0, (bx.y0 + bx.y1) / 2.0);
        // the expression's frame (rotated 90 degrees): its left end at -W/2, its top at the reference line
        let (fx0, fy0) = (-0.5 * l.bbox.w, -l.bbox.h);
        for (s, off, size) in [(&a, l.a, l.size_a), (&b, l.b, l.size_b)] {
            let (fx, fy) = (fx0 + off.0, fy0 + off.1);
            draw_text(c, &ctx, xref - fy, yref - fx, s, size, Face::Plain, Rgba::BLACK, 0.0, 0.0, 90.0, 0.9);
        }
    });
    (grob, width)
}

/// Draws text with R's justification (for geoms).
#[allow(clippy::too_many_arguments)]
pub fn text(c: &mut Canvas, ctx: &Ctx, x: f64, y: f64, s: &str, size_pt: f64, colour: Rgba, hjust: f64, vjust: f64) {
    draw_text(c, ctx, x, y, s, size_pt, Face::Plain, colour, hjust, vjust, 0.0, 1.2);
}

/// facet_wrap(~ GeneGroup, scales = "free_y") with its panels in one row (the missed variants plot with more than
/// 10 genes), as ggplot2 4.0.3's FacetWrap lays it out: the panels with the panel spacing between them, axes woven
/// around every panel (a left axis per panel: the flipped genes' scale is free; a bottom axis per panel), and the
/// strips (blank here, so of no size) between the top axes and the panels.
pub fn table_wrap(th: &Theme, spec: Spec) -> Gtable {
    let ctx = th.ctx;
    let n = spec.h.len();
    // the panels as `table` draws them, through a one-panel spec each
    let mut cells: Vec<(Grob, guides::Axis, guides::Axis)> = Vec::new();
    for c in 0..n {
        let one = Spec {
            title: None,
            subtitle: None,
            xlab: None,
            ylab: None,
            ylab_r: None,
            h: vec![spec.h[c].clone()],
            v: vec![spec.v[c].clone()],
            flip: spec.flip,
            strips_top: None,
            strips_right: None,
            layers: vec![vec![spec.layers[0][c].clone()]],
            h_text: spec.h_text.clone(),
            v_text: spec.v_text.clone(),
            h_ticks: spec.h_ticks,
            v_ticks: spec.v_ticks,
            right_axis: false,
            grid_h: spec.grid_h,
            grid_v: spec.grid_v,
            spacing_x: spec.spacing_x,
            spacing_y: spec.spacing_y,
            margin: [0.0; 4],
            title_blank: true,
        };
        let (panel, bottom, left) = panel_parts(th, one);
        cells.push((panel, bottom, left));
    }
    let mut widths = Vec::new();
    let mut heights = vec![Unit::ZERO, Unit::ZERO, Unit::Null(1.0)];
    let bottom_h = cells.iter().map(|c| c.1.size).fold(0.0, f64::max);
    heights.push(Unit::Abs(bottom_h));
    for (i, c) in cells.iter().enumerate() {
        if i > 0 {
            widths.push(Unit::Abs(spec.spacing_x));
        }
        widths.push(Unit::Abs(c.2.size));
        widths.push(Unit::Null(1.0));
        widths.push(Unit::ZERO);
    }
    let mut gt = Gtable::new(widths, heights);
    for (i, (panel, bottom, left)) in cells.into_iter().enumerate() {
        let col = 1 + i * 4;
        gt.add(panel, 3, col + 1, 3, col + 1, Some(1.0), true, &format!("panel-{}-1", i + 1));
        gt.add(left.grob, 3, col, 3, col, Some(3.0), false, &format!("axis-l-1-{}", i + 1));
        gt.add(Grob::Zero, 3, col + 2, 3, col + 2, Some(3.0), false, &format!("axis-r-1-{}", i + 1));
        gt.add(bottom.grob, 4, col + 1, 4, col + 1, Some(3.0), false, &format!("axis-b-{}-1", i + 1));
        gt.add(Grob::Zero, 1, col + 1, 1, col + 1, Some(3.0), false, &format!("axis-t-{}-1", i + 1));
        gt.add(Grob::Zero, 2, col + 1, 2, col + 1, Some(2.0), false, &format!("strip-t-{}-1", i + 1));
    }
    let title = spec.title.as_ref().map(|t| TitleGrob::new(t, &th.title_el(), true, true));
    let parts = PlotParts {
        panels: gt,
        xlab_t: None,
        xlab_b: spec.xlab.as_ref().map(|l| axis_title(th, l, Side::Bottom)),
        ylab_l: spec.ylab.as_ref().map(|l| axis_title(th, l, Side::Left)),
        ylab_l_custom: None,
        ylab_r: None,
        legend_right: None,
        legend_left: None,
        legend_top: None,
        legend_bottom: None,
        legend_box_spacing: th.pt(2.0 * HALF),
        title,
        subtitle: None,
        caption: None,
        title_on_plot: false,
        margin: spec.margin,
        background: Some(RectEl { fill: Some(Rgba::WHITE), colour: Some(Rgba::WHITE), linewidth: LINE }),
    };
    plot_table(ctx, parts)
}

/// One panel's grob with its bottom and left axes (for table_wrap).
fn panel_parts(th: &Theme, spec: Spec) -> (Grob, guides::Axis, guides::Axis) {
    let t = table(th, spec);
    let find = |name: &str| t.items.iter().find(|p| p.name == name).map(|p| p.grob.clone()).unwrap_or(Grob::Zero);
    let (pt, pl, pb, pr) = hmf_plot::ggp::facet::panel_extent(&t);
    let _ = (pt, pr);
    let bottom_h = t.heights[pb].abs();
    let left_w = t.widths[pl - 2].abs();
    (find("panel"), guides::Axis { grob: find("axis-b"), size: bottom_h }, guides::Axis { grob: find("axis-l"), size: left_w })
}
