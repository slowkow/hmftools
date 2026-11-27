//! Track drawing added for LINX's configurations: connector tracks, rules on links, and axes placed by
//! `position` lists as well as `spacing` (Circos::Division::make_divisions).

use super::conf::{self, Block, Unit};
use super::{expr, glyph, seek, Circos, DatumResolver, Datum, Link, ANGLE_STEP};
use crate::display::{LineCap, Stroke};

impl Circos<'_> {
    /// Fills highlights as GD fills polygons, covering the pixels on both edges: half a pixel more on each
    /// side radially, which keeps bands a pixel or two thick (LINX's segment and gene tracks) as Circos shows
    /// them. Off by default.
    pub fn with_gd_fill_extent(mut self, on: bool) -> Self {
        self.gd_fill_extent = on;
        self
    }

    /// Sizes labels as Circos does from GD's string box (libgd's glyph metrics, see `gdtext`) for the faces
    /// whose metrics are known (TrueType: CMU Bright Roman), rather than by their rendered ink. Off by default.
    pub fn with_gd_text_bounds(mut self, on: bool) -> Self {
        self.gd_text_bounds = on;
        self
    }

    /// A label's (width, height) for placement: GD's string box with `with_gd_text_bounds`, else the ink.
    pub(super) fn text_size(&self, font: crate::text::Font, size_pt: f64, text: &str) -> (f64, f64) {
        if self.gd_text_bounds {
            if let Some(wh) = crate::gdtext::gd_label_size(font, size_pt, text) {
                return wh;
            }
        }
        Self::label_size(font, size_pt, text)
    }

    /// With `with_gd_fill_extent`, a filled slice's radii moved half a pixel apart.
    pub(super) fn gd_extent(&self, r0: f64, r1: f64, filled: bool) -> (f64, f64) {
        if !self.gd_fill_extent || !filled || r0 == r1 {
            return (r0, r1);
        }
        if r0 < r1 { (r0 - 0.5, r1 + 0.5) } else { (r0 + 0.5, r1 - 0.5) }
    }

    /// Applies a link track's `<rules>` (as for data tracks) and drops the links a rule hides.
    pub(super) fn apply_link_rules(&self, track: &Block, links: Vec<Link>) -> Result<Vec<Link>, String> {
        let Some(rules) = track.block("rules") else { return Ok(links) };
        let rules: Vec<&Block> = rules.blocks("rule").collect();
        let mut out = Vec::with_capacity(links.len());
        for mut l in links {
            let mut d = Datum { chr: l.a.0.clone(), start: l.a.1, end: l.a.2, value: String::new(), params: std::mem::take(&mut l.params), hidden: false };
            for rule in &rules {
                let Some(cond) = rule.get("condition") else { continue };
                let pass = expr::eval(cond, &DatumResolver { d: &d, conf: Some(&self.conf) }).map_err(|e| format!("rule '{cond}': {e}"))? != 0.0;
                if pass {
                    for (k, v) in &rule.params {
                        if k == "condition" || k == "flow" {
                            continue;
                        }
                        if k == "show" {
                            d.hidden = !conf::is_true(v);
                            continue;
                        }
                        d.params.push((k.clone(), v.clone()));
                    }
                }
                if pass && !rule.get("flow").is_some_and(|f| f.starts_with("continue")) {
                    break;
                }
            }
            if !d.hidden {
                l.params = d.params;
                out.push(l);
            }
        }
        Ok(out)
    }

    /// A connector track: per datum a line from r0 out to r1, bending from its start to its end angle in the
    /// middle part (`connector_dims`), drawn as GD draws lines (butt ends).
    pub(super) fn draw_connector_track(&mut self, track: &Block, data: Vec<Datum>) -> Result<(), String> {
        for d in &data {
            let Some(&i) = self.by_chr.get(&d.chr) else { continue };
            let id = &self.ideograms[i];
            if d.start < id.start || d.start > id.end || d.end < id.start || d.end > id.end {
                continue;
            }
            let p = Some(d.params.as_slice());
            let r0 = self.radius(seek(p, track, "connector", "r0").unwrap_or("0r"))?;
            let r1 = self.radius(seek(p, track, "connector", "r1").unwrap_or("0r"))?;
            let rd = (r0 - r1).abs();
            let a0 = self.angle(&d.chr, d.start as f64);
            let a1 = self.angle(&d.chr, d.end as f64);
            let dims: Vec<f64> = seek(p, track, "connector", "connector_dims").unwrap_or("0,0.3,0.4,0.3,0").split(',').map(|x| x.trim().parse().unwrap_or(0.0)).collect();
            let thickness = seek(p, track, "connector", "thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(1.0);
            let Some(colour) = self.colour(seek(p, track, "connector", "color")) else { continue };
            let segs = glyph::connector_segments(a0, a1, r0, rd, &dims, ANGLE_STEP);
            // consecutive pieces as one stroke, so that a translucent connector blends once
            let mut pts: Vec<(f64, f64)> = Vec::new();
            for ((sa, sr), (ea, er)) in segs {
                let s = self.xy(sa, sr);
                if pts.last() != Some(&s) {
                    pts.push(s);
                }
                pts.push(self.xy(ea, er));
            }
            pts.dedup();
            if pts.len() >= 2 {
                self.fig.polyline(pts, Stroke { color: colour, width: thickness.max(1.0), dash: Vec::new(), cap: LineCap::Butt });
            }
        }
        Ok(())
    }
}

/// Circos::Division::parse_value with a range: "<x>r" is min + x (max - min), with one range argument x times it.
fn parse_value(v: &str, min: f64, range: f64, two: bool) -> Option<f64> {
    let v = v.trim();
    if let Some(x) = v.strip_suffix('r') {
        let x: f64 = x.trim().parse().ok()?;
        Some(if two { min + x * range } else { x * range })
    } else {
        v.parse().ok()
    }
}

/// Perl's number stringification (%.15g), the key of Circos's division hashes
fn perl_key(x: f64) -> String {
    let s = format!("{:.14e}", x);
    let v: f64 = s.parse().unwrap_or(x);
    format!("{v}")
}

/// Circos::Division::make_divisions: the axis positions of a track's `<axes>` between min and max, each
/// with its block: from `position` lists and `spacing`, larger spacings first, a position drawn once.
pub(super) fn make_divisions<'b>(axes: &'b Block, min: f64, max: f64) -> Vec<(f64, &'b Block)> {
    // (spacing, position, block), first block of a (spacing, position) kept
    let mut divisions: Vec<(f64, f64, &Block)> = Vec::new();
    for axis in axes.blocks("axis") {
        if axis.flag("show") == Some(false) {
            continue;
        }
        let mut pos: Vec<f64> = Vec::new();
        if let Some(list) = axis.get("position") {
            for p in list.split(',') {
                if let Some(v) = parse_value(p, min, max - min, true) {
                    pos.push(v);
                }
            }
        }
        let spacing = axis.get("spacing").and_then(|s| parse_value(s, min, max - min, false)).unwrap_or(0.0);
        let ystart = axis.get("y0").and_then(|v| parse_value(v, min, max - min, true)).filter(|v| *v >= min).unwrap_or(min);
        let yend = axis.get("y1").and_then(|v| parse_value(v, min, max - min, true)).filter(|v| *v <= max).unwrap_or(max);
        if spacing != 0.0 {
            let tol = spacing / 1000.0;
            let mut y = min;
            while y <= max + tol {
                if y >= ystart - tol {
                    if y > yend + tol {
                        break;
                    }
                    pos.push(y);
                }
                y += spacing;
            }
        }
        let mut seen: Vec<String> = Vec::new();
        let mut uniq: Vec<f64> = Vec::new();
        for p in pos {
            let k = perl_key(p);
            if !seen.contains(&k) {
                seen.push(k);
                uniq.push(p);
            }
        }
        uniq.sort_by(|a, b| a.total_cmp(b));
        let tol = spacing / 1000.0;
        for p in uniq {
            if p < ystart - tol || p > yend + tol {
                continue;
            }
            if !divisions.iter().any(|(s, q, _)| *s == spacing && perl_key(*q) == perl_key(p)) {
                divisions.push((spacing, p, axis));
            }
        }
    }
    let mut spacings: Vec<f64> = Vec::new();
    for (s, _, _) in &divisions {
        if !spacings.contains(s) {
            spacings.push(*s);
        }
    }
    spacings.sort_by(|a, b| b.total_cmp(a));
    let mut seen: Vec<String> = Vec::new();
    let mut out = Vec::new();
    for s in spacings {
        let mut at: Vec<(f64, &Block)> = divisions.iter().filter(|(sp, _, _)| *sp == s).map(|(_, p, b)| (*p, *b)).collect();
        at.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (p, b) in at {
            let k = perl_key(p);
            if s != 0.0 && seen.contains(&k) {
                continue;
            }
            if s != 0.0 {
                seen.push(k);
            }
            out.push((p, b));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn divisions_from_positions_and_spacing() {
        let b = conf::parse("<axes>\n<axis>\nposition = 0.5r,1r\n</axis>\n<axis>\nspacing = 0.25r\n</axis>\n</axes>\n").unwrap();
        let axes = b.block("axes").unwrap();
        let d: Vec<f64> = make_divisions(axes, -2.0, 0.0).iter().map(|x| x.0).collect();
        // spacing 0.5 first (-2, -1.5, -1, -0.5, 0), then the listed positions (-1, 0), which spacing 0 does not
        // de-duplicate
        assert_eq!(d, vec![-2.0, -1.5, -1.0, -0.5, 0.0, -1.0, 0.0]);
        let single = conf::parse("<axes>\n<axis>\nspacing = 0.25r\n</axis>\n</axes>\n").unwrap();
        let d: Vec<f64> = make_divisions(single.block("axes").unwrap(), 0.0, 1.0).iter().map(|x| x.0).collect();
        assert_eq!(d, vec![0.0, 0.25, 0.5, 0.75, 1.0]);
    }
}
