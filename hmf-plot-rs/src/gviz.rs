//! Gviz's IdeogramTrack as `plotTracks(track, from, to, add = TRUE, panel.only = FALSE)` draws it into a
//! viewport (Gviz 1.54, drawGD for IdeogramTrack, .roundedCap, .getBioColorIdeo): the chromosome's bands as
//! rectangles between the bevel points, rounded caps at both ends cut into the bands they cover, centromere
//! triangles, black outlines, the chromosome's name on the left, and the plotted range as a filled box under
//! the bands and an outlined box over them. Sizes follow R's png() device (72 dpi: 1 point of text = 1 pixel,
//! line width 1 = 0.75 pixel). Added for LINX's chromosome range strip.

use crate::color::Rgba;
use crate::display::{Figure, LineCap, Stroke};
use crate::fmath::FMath;
use crate::text::{self, Font};

/// One band of the cytoband table (UCSC: chromStart, chromEnd, gieStain).
#[derive(Clone, Debug)]
pub struct Band {
    pub start: i64,
    pub end: i64,
    pub stain: String,
}

/// Gviz's colour of a band stain: gneg white, gposN on a 100-step ramp from white to black (truncated as
/// R's rgb() with maxColorValue 255 does), acen brown4, stalk brown3, gvar and gpos black.
pub fn band_colour(stain: &str) -> Rgba {
    match stain {
        "gneg" => Rgba::rgb8(255, 255, 255),
        "acen" => Rgba::rgb8(139, 35, 35),
        "stalk" => Rgba::rgb8(205, 51, 51),
        "gvar" | "gpos" => Rgba::rgb8(0, 0, 0),
        s => match s.strip_prefix("gpos").and_then(|n| n.parse::<i32>().ok()) {
            Some(n) if (1..=100).contains(&n) => {
                let v = (255.0 - 255.0 * (n - 1) as f64 / 99.0) as u8;
                Rgba::rgb8(v, v, v)
            }
            _ => Rgba::rgb8(0, 0, 0),
        },
    }
}

/// An IdeogramTrack to draw.
pub struct Ideogram<'a> {
    /// the track's name, shown on the left
    pub label: &'a str,
    /// the chromosome's bands, in position order
    pub bands: &'a [Band],
    /// the plotted range (plotTracks' from and to)
    pub from: f64,
    pub to: f64,
    /// fill and col of the range boxes
    pub colour: Rgba,
    /// the label's font size (points)
    pub fontsize: f64,
    /// lwd of the range box
    pub lwd: f64,
}

/// grid's 3 points of inner margin at 72 dpi
const INNER_MARGIN: f64 = 3.0 * 72.0 / 72.27;
/// plotTracks' margin (pixels)
const MARGIN: f64 = 6.0;
/// R's line width unit at 72 dpi: 1/96 inch
const LWD: f64 = 72.0 / 96.0;

/// .roundedCap's outline in [0, 1] x [0, 1] (n = 100, the bevel display parameter 0.45)
fn cap_coords(bevel: f64) -> Vec<(f64, f64)> {
    let n = 100.0f64;
    let bevel = f64::max(1.0 / n, f64::min(bevel, 0.5));
    let k = f64::max(2.0, (n * bevel).floor()) as usize;
    // R's seq(a, b, length.out = k): a + i * by, by = (b - a) / (k - 1)
    let seq = |a: f64, b: f64| -> Vec<f64> {
        let by = (b - a) / (k - 1) as f64;
        (0..k).map(|i| a + i as f64 * by).collect()
    };
    let pi = std::f64::consts::PI;
    let mut xs = vec![0.0];
    let mut ys = vec![1.0 + 0.5 - bevel];
    for t in seq(0.0, pi / 2.0) {
        xs.push(t.fsin() + bevel * 4.0);
        ys.push(t.fcos() + (0.5 - bevel));
    }
    for t in seq(pi / 2.0, pi) {
        xs.push(t.fsin() + bevel * 4.0);
        ys.push(t.fcos() - (0.5 - bevel));
    }
    xs.push(0.0);
    ys.push(pi.fcos() - (0.5 - bevel));
    xs.iter().zip(ys.iter()).map(|(x, y)| (x / (1.0 + 4.0 * bevel), (y + 1.0 + (0.5 - bevel)) / (2.0 + 2.0 * (0.5 - bevel)))).collect()
}

/// R's cut(x, breaks, right = TRUE, include.lowest = TRUE, labels = FALSE): the 1-based interval, or None
fn cut(x: f64, breaks: &[f64]) -> Option<usize> {
    for j in 0..breaks.len().saturating_sub(1) {
        let (lo, hi) = (breaks[j], breaks[j + 1]);
        if (x > lo || (j == 0 && x >= lo)) && x <= hi {
            return Some(j + 1);
        }
    }
    None
}

/// Draws the track into the plotTracks cell (x, y, w, h) (pixels, y down).
pub fn draw(fig: &mut Figure, cell: (f64, f64, f64, f64), t: &Ideogram) {
    if t.bands.is_empty() {
        return;
    }
    let (cx, cy, cw, ch) = cell;
    // vpBound (6 px margins), then the content viewport (no title, 3 points of spacing each side)
    let left = cx + MARGIN + INNER_MARGIN;
    let width = cw - 2.0 * MARGIN - 2.0 * INNER_MARGIN;
    let top = cy + MARGIN;
    let height = ch - 2.0 * MARGIN;
    let npc_y = |v: f64| top + (1.0 - v) * height;

    // the name: a viewport as wide as the text plus 0.2 inches, at most 0.75 of the track
    let px = t.fontsize; // points = pixels at 72 dpi
    let text_w_in = text::advance(Font::Sans, px, t.label) / 72.0;
    let nspace = f64::min((text_w_in + 0.2) / (width / 72.0), 0.75);
    let asc = text::extents(Font::Sans, px, "M").ascent;
    // grid.text(just = c("left", "centre")): the baseline half the ascent of "M" below the middle
    fig.text(left.round(), (npc_y(0.5) + asc / 2.0).round(), t.label, Font::Sans, px, Rgba::BLACK, 0.0);

    let il = left + nspace * width;
    let iw = (1.0 - nspace) * width;
    let nx = |v: f64| il + v * iw;
    let len = t.bands.iter().map(|b| b.end).max().unwrap_or(1) as f64;

    // the range, filled, under the bands
    let rx0 = t.from / len;
    let rw = f64::min(1.0, (t.to - t.from) / len);
    fig.rect(nx(rx0), npc_y(0.9), rw * iw, npc_y(0.1) - npc_y(0.9), t.colour);

    let cols: Vec<Rgba> = t.bands.iter().map(|b| band_colour(&b.stain)).collect();
    let st: Vec<f64> = t.bands.iter().map(|b| b.start as f64 / len).collect();
    let ed: Vec<f64> = t.bands.iter().map(|b| b.end as f64 / len).collect();
    let bevel = 0.02;
    // the bands overlapping the bevel points (IRanges at the truncated positions, width 1)
    let p1 = (bevel * len) as i64;
    let p2 = ((1.0 - bevel) * len) as i64;
    let ol: Vec<usize> = t.bands.iter().enumerate().filter(|(_, b)| (b.start <= p1 && b.end >= p1) || (b.start <= p2 && b.end >= p2)).map(|(i, _)| i).collect();
    let (stext, valsext): (Vec<f64>, Vec<usize>) = if t.bands.len() == 1 || ol.len() < 2 {
        (vec![0.0, bevel, 1.0 - bevel], vec![0, 0, 0])
    } else {
        let (o1, o2) = (ol[0], ol[1]);
        let mut s: Vec<f64> = st[..=o1].to_vec();
        let mut v: Vec<usize> = (0..=o1).collect();
        s.push(bevel);
        v.push(o1);
        s.extend_from_slice(&st[o1 + 1..=o2]);
        v.extend(o1 + 1..=o2);
        s.push(1.0 - bevel);
        v.push(o2);
        if o2 + 1 < st.len() {
            s.extend_from_slice(&st[o2 + 1..]);
            v.extend(o2 + 1..st.len());
        }
        (s, v)
    };
    let mut wd: Vec<f64> = Vec::with_capacity(stext.len());
    for i in 0..stext.len() {
        wd.push(if i + 1 < stext.len() { stext[i + 1] } else { 1.0 } - stext[i]);
    }
    let ls = if ol.len() >= 2 { ol[0] + 1 } else { 1 };
    let rs = if ol.len() >= 2 { ol[1] + 1 } else { 1 };
    let cent: Vec<usize> = valsext.iter().enumerate().filter(|(_, b)| t.bands[**b].stain == "acen").map(|(i, _)| i).collect();
    let (bef, aft): (Vec<usize>, Vec<usize>) = if let (Some(&cmin), Some(&cmax)) = (cent.first(), cent.last()) {
        ((ls..cmin.max(ls)).collect(), (cmax + 1..=rs).collect())
    } else {
        ((ls..=rs).collect(), Vec::new())
    };
    let margin = 0.3;
    let thin = |c: Rgba| Some(Stroke { color: c, width: LWD, dash: Vec::new(), cap: LineCap::Round });
    for &i in bef.iter().chain(aft.iter()) {
        let c = cols[valsext[i]];
        let (x0, x1) = (nx(stext[i]), nx(stext[i] + wd[i]));
        fig.polygon(vec![(x0, npc_y(margin)), (x1, npc_y(margin)), (x1, npc_y(1.0 - margin)), (x0, npc_y(1.0 - margin))], Some(c), thin(c));
    }
    let acen = band_colour("acen");
    if let (Some(&cmin), Some(&cmax)) = (cent.first(), cent.last()) {
        let (a, aw) = (stext[cmin], wd[cmin]);
        fig.polygon(vec![(nx(a), npc_y(margin)), (nx(a + aw), npc_y(0.5)), (nx(a), npc_y(1.0 - margin)), (nx(a), npc_y(margin))], Some(acen), thin(acen));
        let (b, bw) = (stext[cmax], wd[cmax]);
        fig.polygon(vec![(nx(b), npc_y(0.5)), (nx(b + bw), npc_y(margin)), (nx(b + bw), npc_y(1.0 - margin)), (nx(b), npc_y(0.5))], Some(acen), thin(acen));
    }

    // the caps: the outline split into the bands it covers, each part filled with its band's colour
    let unit = cap_coords(0.45);
    let breaks_l: Vec<f64> = if st.len() == 1 { vec![0.0, 1.0] } else { st.clone() };
    let breaks_r: Vec<f64> = if ed.len() == 1 { vec![1.0, 2.0] } else { ed.clone() };
    let cap = |side_left: bool, bl: (f64, f64), tr: (f64, f64), breaks: &[f64], fig: &mut Figure| -> Vec<(f64, f64)> {
        let (dx, dy) = ((tr.0 - bl.0).abs(), (tr.1 - bl.1).abs());
        let pts: Vec<(f64, f64)> =
            unit.iter().map(|&(x, y)| if side_left { ((1.0 - x) * dx + bl.0, (1.0 - y) * dy + bl.1) } else { (x * dx + bl.0, y * dy + bl.1) }).collect();
        // split by band interval, groups in interval order
        let mut groups: Vec<(usize, Vec<(f64, f64)>)> = Vec::new();
        for &p in &pts {
            if let Some(j) = cut(p.0, breaks) {
                match groups.iter_mut().find(|g| g.0 == j) {
                    Some(g) => g.1.push(p),
                    None => groups.push((j, vec![p])),
                }
            }
        }
        groups.sort_by_key(|g| g.0);
        let shift = if side_left { 0 } else { 1 };
        for (gi, (j, g)) in groups.iter().enumerate() {
            let mut poly: Vec<(f64, f64)> = Vec::new();
            // after the first, a piece is closed by the last point of interval j - 1 (none when that interval
            // has no points: R's lookup of the missing group gives NULL)
            let prev = if gi > 0 { groups.iter().find(|x| x.0 + 1 == *j).and_then(|x| x.1.last().copied()) } else { None };
            match prev {
                Some(last) => {
                    poly.push((last.0, 1.0 - last.1));
                    poly.extend(g.iter().copied());
                    poly.push(last);
                }
                None => poly.extend(g.iter().copied()),
            }
            let band = (j - 1 + shift).min(cols.len() - 1);
            let c = cols[band];
            fig.polygon(poly.iter().map(|&(x, y)| (nx(x), npc_y(y))).collect(), Some(c), thin(c));
        }
        pts
    };
    let lc = cap(true, (stext[0], margin), (stext[ls], 1.0 - margin), &breaks_l, fig);
    let rc = cap(false, (*stext.last().unwrap(), margin), (1.0, 1.0 - margin), &breaks_r, fig);
    let black = Stroke { color: Rgba::BLACK, width: LWD, dash: Vec::new(), cap: LineCap::Round };
    fig.polyline(lc.iter().map(|&(x, y)| (nx(x), npc_y(y))).collect(), black.clone());
    fig.polyline(rc.iter().map(|&(x, y)| (nx(x), npc_y(y))).collect(), black.clone());
    let lc_max = lc.iter().map(|p| p.0).fold(f64::MIN, f64::max);
    let mut segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    if let (Some(&cmin), Some(&cmax)) = (cent.first(), cent.last()) {
        let after = stext.get(cmax + 1).copied().unwrap_or(1.0);
        let last = *stext.last().unwrap();
        for y in [margin, 1.0 - margin] {
            segs.push(((lc_max, y), (stext[cmin], y)));
        }
        for y in [margin, 1.0 - margin] {
            segs.push(((after, y), (last, y)));
        }
        segs.push(((stext[cmin], margin), (stext[cmin] + wd[cmin], 0.5)));
        segs.push(((stext[cmin], 1.0 - margin), (stext[cmin] + wd[cmin], 0.5)));
        segs.push(((stext[cmax], 0.5), (stext[cmax] + wd[cmax], margin)));
        segs.push(((stext[cmax], 0.5), (stext[cmax] + wd[cmax], 1.0 - margin)));
    } else {
        let mx = stext.iter().copied().fold(f64::MIN, f64::max);
        for y in [margin, 1.0 - margin] {
            segs.push(((lc_max, y), (mx, y)));
        }
    }
    for (a, b) in segs {
        fig.line(nx(a.0), npc_y(a.1), nx(b.0), npc_y(b.1), black.clone());
    }

    // the range, outlined, over the bands
    let box_stroke = Stroke { color: t.colour, width: t.lwd * LWD, dash: Vec::new(), cap: LineCap::Round };
    let (bx0, bx1) = (nx(rx0), nx(rx0 + rw));
    fig.polygon(vec![(bx0, npc_y(0.9)), (bx1, npc_y(0.9)), (bx1, npc_y(0.1)), (bx0, npc_y(0.1))], None, Some(box_stroke));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stain_colours_as_gviz() {
        assert_eq!(band_colour("gpos15").to_hex(), "#DADADA");
        assert_eq!(band_colour("gpos20").to_hex(), "#CECECE");
        assert_eq!(band_colour("gpos25").to_hex(), "#C1C1C1");
        assert_eq!(band_colour("gpos30").to_hex(), "#B4B4B4");
        assert_eq!(band_colour("gpos31").to_hex(), "#B1B1B1");
        assert_eq!(band_colour("gpos100").to_hex(), "#000000");
        assert_eq!(band_colour("acen").to_hex(), "#8B2323");
    }

    #[test]
    fn cap_outline() {
        let c = cap_coords(0.45);
        assert_eq!(c.len(), 92);
        assert_eq!(c[0].0, 0.0);
        assert!((c[0].1 - 1.0).abs() < 1e-12);
        // the outermost point is in the middle
        let mx = c.iter().map(|p| p.0).fold(f64::MIN, f64::max);
        assert!((mx - 1.0).abs() < 1e-12);
        assert_eq!(cut(0.0, &[0.0, 0.5, 1.0]), Some(1));
        assert_eq!(cut(0.5, &[0.0, 0.5, 1.0]), Some(1));
        assert_eq!(cut(0.6, &[0.0, 0.5, 1.0]), Some(2));
        assert_eq!(cut(1.5, &[0.0, 0.5, 1.0]), None);
    }
}
