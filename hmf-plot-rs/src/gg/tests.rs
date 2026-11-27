//! ggplot2 layout and stat tests against what R 4.5 / ggplot2 4.0.3 computes (ggplotGrob widths in pt,
//! layer_data), on synthetic inputs.

use super::*;

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn scale_view_expands_and_breaks_like_ggplot() {
    // copy-number y scale: limits 0..575789, 5% each side, breaks of the expanded range
    let v = PosScale::new((0.0, 575789.0)).view();
    assert!(approx(v.range.0, -28789.45, 1e-6) && approx(v.range.1, 604578.45, 1e-6));
    assert_eq!(v.majors, [0.0, 2e5, 4e5, 6e5]);
    assert_eq!(v.labels, ["0e+00", "2e+05", "4e+05", "6e+05"]);
    // clonality x: limits 0..3.5, expand = c(0.01, 0.01)
    let mut s = PosScale::new((0.0, 3.5));
    s.expand = (0.01, 0.01);
    let v = s.view();
    assert!(approx(v.range.0, -0.045, 1e-12) && approx(v.range.1, 3.545, 1e-12));
    assert_eq!(v.majors, [0.0, 1.0, 2.0, 3.0]);
    // fixed breaks with labels, censored to the range
    let mut s = PosScale::new((0.0, 1.0));
    s.breaks = Some(vec![0.25, 0.5, 0.75, 1.0, 1.5]);
    s.labels = Labels::Fixed(vec!["25%".into(), "50%".into(), "75%".into(), "100%".into(), "150%".into()]);
    assert_eq!(s.view().labels, ["25%", "50%", "75%", "100%"]);
}

#[test]
fn log_scale_view() {
    let mut s = PosScale::new((20.0, 2.0e6));
    s.trans = Trans::Log10;
    s.labels = Labels::FormatNoSci;
    let v = s.view();
    // majors at even powers, labels padded to a common width, minors at the odd powers
    assert_eq!(v.majors.iter().map(|m| m.round() as i64).collect::<Vec<_>>(), [2, 3, 4, 5, 6]);
    assert_eq!(v.labels[0], "    100");
    assert_eq!(v.labels[4], "1000000");
}

#[test]
fn left_columns_match_ggplot_grob() {
    // R: ylab-l 13.3496 pt + axis-l 31.8307 pt for labels 0e+00 .. 6e+05 (no ticks)
    let v = PosScale::new((0.0, 575789.0)).view();
    let p = Plot { title: None, x: v.clone(), y: v, x_title: None, y_title: Some("Baf Count".into()), legends: vec![], ticks: false, minor_grid: false, panel_border: true };
    let expected = pt(13.3496 + 31.8307);
    assert!(approx(p.left_width(), expected, 2.5), "{} vs {}", p.left_width(), expected);
}

#[test]
fn legend_sizes_match_ggplot_grob() {
    // copy-number legend: widths 5.5 17.3448 37.2988 5.5 pt, heights 5.5 17.3448 0 17.3448 5.5 pt
    let key = |c: (u8, u8, u8)| Key::Polygon { fill: Rgba::rgb8(c.0, c.1, c.2), colour: Some(Rgba::BLACK), linewidth: 0.5 };
    let l = Legend::Keys { title: None, entries: vec![("MACN0".into(), key((141, 211, 199))), ("MACN1".into(), key((255, 255, 179)))] };
    let (w, h) = l.size();
    assert!(approx(w, pt(5.5 + 17.3448 + 37.2988 + 5.5), 3.0), "w {w}");
    assert!(approx(h, pt(5.5 + 2.0 * 17.3448 + 5.5), 0.5), "h {h}");
    // colour bar with a two-line title: widths 5.5 17.3448 30.5536 5.5, heights 5.5 28.0241 86.7240 5.5
    let bar = Legend::Bar {
        title: Some("Relative\nScore".into()),
        ramp: GradientN::new(&[Rgba::rgb8(0, 0, 255), Rgba::rgb8(255, 0, 0)], None),
        limits: (0.0, 1.0),
        breaks: vec![0.1, 0.5, 1.0],
        labels: vec!["10%".into(), "50%".into(), "100%".into()],
    };
    let (w, h) = bar.size();
    assert!(approx(w, pt(5.5 + 17.3448 + 30.5536 + 5.5), 3.0), "w {w}");
    assert!(approx(h, pt(5.5 + 28.0241 + 86.7240 + 5.5), 3.0), "h {h}");
}

#[test]
fn stat_bin_breaks() {
    // binwidth 0.1 over limits (-0.1, 3.1): bins centred on multiples of 0.1, closed on the right
    let b = bin_breaks((-0.1, 3.1), 0.1);
    assert!(approx(b[0], -0.15, 1e-12));
    assert!(approx(*b.last().unwrap(), 3.15, 1e-9));
    assert_eq!(b.len(), 34);
    // 0.05 belongs to the bin (-0.05, 0.05]; 0.0500001 to the next
    let i = bin_index(&b, 0.05).unwrap();
    assert!(approx(b[i], -0.05, 1e-9));
    assert_eq!(bin_index(&b, 0.0500001), Some(i + 1));
    assert_eq!(bin_index(&b, 5.0), None);
}

#[test]
fn point_radius() {
    // size 1.5 (ggplot2's default), stroke 0.5: fontsize 1.5 * .pt + 0.5 * .stroke / 2 = 5.212 pt
    let mut f = Figure::new(10, 10, Rgba::WHITE);
    point(&mut f, 5.0, 5.0, 1.5, Rgba::BLACK, 0.5);
    if let Item::Circle { r, stroke: Some(s), .. } = &f.items[0] {
        assert!(approx(*r, 0.375 * (1.5 * PT + 0.5 * STROKE / 2.0) * 300.0 / 72.0, 1e-9));
        assert!(approx(s.width, 0.5 * STROKE / 2.0 * 300.0 / 96.0, 1e-9));
    } else {
        panic!("no circle");
    }
}
