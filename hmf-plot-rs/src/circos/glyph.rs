//! Scatter glyphs other than circles (Circos.pm's scatter track): squares, rectangles and triangles as
//! polygons turned to the data point's angle, vertices rounded to whole pixels (rotate_xy). Added for LINX's
//! segment ends (squares at telomeres and centromeres, triangles at foldbacks).

use crate::fmath::FMath;

const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
const SQRT3_HALF: f64 = 0.866_025_403_784_438_6;

/// The polygon of a glyph of `size` pixels centred at (x, y) for a point at `angle` (degrees, with Circos's
/// angle offset `angle_offset`), or None for glyphs drawn otherwise (circle) or not known.
pub fn polygon(glyph: &str, x: f64, y: f64, size: f64, angle: f64, angle_offset: f64, angle_shift: f64) -> Option<Vec<(f64, f64)>> {
    let h = size / 2.0;
    let pts: Vec<(f64, f64)> = match glyph {
        "rectangle" | "square" => vec![(x - h, y - h), (x + h, y - h), (x + h, y + h), (x - h, y + h)],
        "triangle" => vec![(x, y - h * SQRT3_HALF), (x + h, y + h * SQRT3_HALF), (x - h, y + h * SQRT3_HALF)],
        _ => return None,
    };
    let a = (angle + angle_shift - angle_offset) * DEG2RAD;
    let (sa, ca) = (a.fsin(), a.fcos());
    Some(
        pts.into_iter()
            .map(|(px, py)| {
                let (xd, yd) = (px - x, py - y);
                (perl_round(xd * ca - yd * sa + x), perl_round(xd * sa + yd * ca + y))
            })
            .collect(),
    )
}

/// Math::Round's round: half away from zero
fn perl_round(v: f64) -> f64 {
    if v >= 0.0 { (v + 0.5).floor() } else { -((-v + 0.5).floor()) }
}

/// The radial and spiral pieces of a connector (Circos.pm's connector track) as (angle, radius) pairs of
/// line segments: radially out from r0 by dims[0..2], then from angle0 to angle1 (in steps of `step`
/// degrees) while moving out by dims[2], then radially by dims[3]. `rd` is |r1 - r0|.
pub fn connector_segments(angle0: f64, angle1: f64, r0: f64, rd: f64, dims: &[f64], step: f64) -> Vec<((f64, f64), (f64, f64))> {
    let d = |i: usize| dims.get(i).copied().unwrap_or(0.0);
    let mut out = vec![((angle0, r0 + d(0) * rd), (angle0, r0 + (d(0) + d(1)) * rd))];
    let rinit = r0 + (d(0) + d(1)) * rd;
    let rfinal = r0 + (d(0) + d(1) + d(2)) * rd;
    let rdiff = (rfinal - rinit).abs();
    if angle1 > angle0 {
        let adiff = angle1 - angle0;
        let mut a = angle0;
        while a + step <= angle1 {
            out.push(((a, rinit + rdiff * (a - angle0) / adiff), (a + step, rinit + rdiff * (a + step - angle0) / adiff)));
            a += step;
        }
        if a < angle1 {
            out.push(((a, rinit + rdiff * (a - angle0) / adiff), (angle1, rfinal)));
        }
    } else if angle1 < angle0 {
        let adiff = angle1 - angle0;
        let mut a = angle0;
        while a - step >= angle1 {
            out.push(((a, rinit + rdiff * (a - angle0) / adiff), (a - step, rinit + rdiff * (a - step - angle0) / adiff)));
            a -= step;
        }
        if a > angle1 {
            out.push(((a, rinit + rdiff * (a - angle0) / adiff), (angle1, rfinal)));
        }
    } else {
        out.push(((angle0, rinit), (angle1, rfinal)));
    }
    out.push(((angle1, rfinal), (angle1, r0 + (d(0) + d(1) + d(2) + d(3)) * rd)));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_turns_with_the_angle() {
        // at the top (angle -90 with offset -90) the square is not turned
        let p = polygon("square", 100.0, 50.0, 20.0, -90.0, -90.0, 0.0).unwrap();
        assert_eq!(p, vec![(90.0, 40.0), (110.0, 40.0), (110.0, 60.0), (90.0, 60.0)]);
        // a triangle at 3 o'clock points outwards (to the right)
        let t = polygon("triangle", 100.0, 100.0, 20.0, 0.0, -90.0, 0.0).unwrap();
        assert_eq!(t[0], (109.0, 100.0));
        assert!(polygon("circle", 0.0, 0.0, 1.0, 0.0, -90.0, 0.0).is_none());
    }

    #[test]
    fn straight_connector() {
        let s = connector_segments(10.0, 10.0, 100.0, 50.0, &[0.0, 0.3, 0.4, 0.3, 0.0], 0.5);
        assert_eq!(s.len(), 3);
        assert_eq!(s[0], ((10.0, 100.0), (10.0, 115.0)));
        assert_eq!(s[2].1, (10.0, 150.0));
    }
}
