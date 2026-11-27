//! R's plotmath (src/main/plotmath.c) for the expressions the hmftools scripts use as axis titles:
//! `atop(a, atop(b))`, a fraction without a bar whose numerator `a` is in text style and whose denominator is
//! another atop with `b` (in script style, 0.7 times the size) over nothing. Boxes follow RenderFraction,
//! NumDenomHShift/VShift and the TeX parameters (sigma8, 9, 11, 12, 22, xi8) on the device's glyph metrics
//! (GEMetricInfo of 'x', 'X', '+', '0', 'g' and of each character), in device units.

use crate::rdev::{Device, Face};

/// A box: height above the baseline, depth below, width.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BBox {
    pub h: f64,
    pub d: f64,
    pub w: f64,
}

/// The layout of atop(a, atop(b)): its box, and each string's baseline-left offset from the box's left end on its
/// baseline (x right, y up), with their font sizes.
#[derive(Clone, Debug)]
pub struct AtopLayout {
    pub bbox: BBox,
    pub a: (f64, f64),
    pub b: (f64, f64),
    pub size_a: f64,
    pub size_b: f64,
}

/// RenderStr's box of a string: the tallest and deepest glyphs, the string's width.
fn str_box(dev: &Device, s: &str, size: f64) -> BBox {
    let (mut h, mut d) = (0.0f64, 0.0f64);
    for ch in s.chars() {
        let (a, de, _) = dev.char_metric(ch, size, Face::Plain);
        h = h.max(a);
        d = d.max(de);
    }
    BBox { h, d, w: dev.str_width(s, size, Face::Plain) }
}

struct Tex {
    axis: f64,
    rule: f64,
    x_upper: f64,
    fig: f64,
    desc: f64,
}

impl Tex {
    fn at(dev: &Device, size: f64) -> Tex {
        let m = |c| dev.char_metric(c, size, Face::Plain);
        Tex { axis: 0.5 * m('+').0, rule: 0.015 * dev.dpi(), x_upper: m('X').0, fig: m('0').0, desc: m('g').1 }
    }
    /// num1 and denom1 (display style) or num2 and denom2 (text style and smaller), and phi
    fn uv(&self, display: bool) -> (f64, f64, f64) {
        if display {
            (
                self.axis + 3.51 * self.rule + 0.15 * self.x_upper + 0.7 * self.desc,
                -self.axis + 3.51 * self.rule + 0.7 * self.fig + 0.344444 * self.x_upper,
                3.0 * self.rule,
            )
        } else {
            (self.axis + 1.51 * self.rule + 0.08333333 * self.x_upper, -self.axis + 1.51 * self.rule + 0.7 * self.fig + 0.08333333 * self.x_upper, self.rule)
        }
    }
}

/// NumDenomVShift
fn vshift(t: &Tex, display: bool, num: BBox, den: BBox) -> (f64, f64) {
    let (mut u, mut v, phi) = t.uv(display);
    let a = t.axis;
    let theta = t.rule;
    let delta = (u - num.d) - (a + 0.5 * theta);
    if delta < phi {
        u += phi - delta;
    }
    let delta = (a + 0.5 * theta) - (den.h - v);
    if delta < phi {
        v += phi - delta;
    }
    (u, v)
}

/// The layout of atop(a, atop(b)) at a font size, as GEMathText lays it out from display style.
pub fn atop_atop(dev: &Device, a: &str, b: &str, size: f64) -> AtopLayout {
    // inner atop(b), rendered in text' style (T1, cex 1): numerator b in S1 (0.7), denominator empty
    let size_b = size * 0.7;
    let tex = Tex::at(dev, size);
    let nb = str_box(dev, b, size_b);
    let empty = BBox { h: 0.0, d: 0.0, w: 0.0 };
    let (ui, vi) = vshift(&tex, false, nb, empty);
    // NumDenomHShift: the narrower centred
    let (nb_shift, _den_shift) = if nb.w > empty.w { (0.0, (nb.w - empty.w) / 2.0) } else { ((empty.w - nb.w) / 2.0, 0.0) };
    let inner = BBox { h: (nb.h + ui).max(empty.h - vi), d: (nb.d - ui).max(empty.d + vi), w: nb.w.max(empty.w) };
    // outer atop: numerator a in text style (T, cex 1), denominator the inner (T1, cex 1), display style shifts
    let na = str_box(dev, a, size);
    let (uo, vo) = vshift(&tex, true, na, inner);
    let (na_shift, inner_shift) = if na.w > inner.w { (0.0, (na.w - inner.w) / 2.0) } else { ((inner.w - na.w) / 2.0, 0.0) };
    let bbox = BBox { h: (na.h + uo).max(inner.h - vo), d: (na.d - uo).max(inner.d + vo), w: na.w.max(inner.w) };
    AtopLayout { bbox, a: (na_shift, uo), b: (inner_shift + nb_shift, -vo + ui), size_a: size, size_b }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atop_box_is_two_lines_deep() {
        let l = atop_atop(&Device::Pdf, "Phred score diff. (del - ins)", "more del errors <-> more ins errors", 10.5);
        assert!(l.bbox.h > 0.0 && l.bbox.d > 0.0);
        assert!(l.a.1 > 0.0 && l.b.1 < 0.0);
        assert!(l.size_b < l.size_a);
    }
}
