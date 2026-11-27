//! Display lists written as vector PDF pages, as R's `pdf()` device writes its figures: the page in points, text in
//! the standard (not embedded) fonts Helvetica and Helvetica-Bold with R's ISO Latin-1 encoding (so that '-' is the
//! minus sign, as R's PDFs have it, and `pdftotext` reads the same characters), filled and stroked paths, circles as
//! four Bézier arcs, dash patterns, clipping rectangles, and transparency through ExtGState alpha. A figure laid out
//! for `rdev::Device::Pdf` has its device units in points with y downwards; pages are flipped to PDF's y upwards.
//!
//! Text in `text::Font::Sans` is written in Helvetica, `Font::SansBold` in Helvetica-Bold (R's "sans" on its pdf
//! device). No dates or random ids: the same figures give the same bytes.

use std::collections::BTreeMap;

use pdf_writer::types::LineCapStyle;
use pdf_writer::{Content, Filter, Finish, Name, Pdf, Rect, Ref, Str, TextStr};

use crate::color::Rgba;
use crate::display::{Figure, Item, LineCap, Stroke};
use crate::text::Font;

/// R's ISO Latin-1 names that differ from WinAnsiEncoding among the printable ASCII characters.
const DIFFERENCES: [(i32, &str); 3] = [(39, "quoteright"), (45, "minus"), (96, "quoteleft")];

/// A page's alpha states (fill alpha, stroke alpha) in 1/255 steps -> resource names.
struct Alphas {
    names: BTreeMap<(u8, u8), String>,
}

impl Alphas {
    fn name(&mut self, fill: f32, stroke: f32) -> Option<String> {
        let q = |a: f32| (a.clamp(0.0, 1.0) * 255.0).round() as u8;
        let key = (q(fill), q(stroke));
        if key == (255, 255) {
            return None;
        }
        let n = self.names.len();
        Some(self.names.entry(key).or_insert_with(|| format!("GS{}", n + 1)).clone())
    }
}

/// An axis-aligned rectangle given as its 4 corners in order (optionally with the first repeated last): (x0, y0, x1,
/// y1) with x0 < x1, y0 < y1 (device units, y down).
fn axis_rect(pts: &[(f64, f64)]) -> Option<(f64, f64, f64, f64)> {
    let pts = match pts {
        [a, b, c, d, e] if e == a => [*a, *b, *c, *d],
        [a, b, c, d] => [*a, *b, *c, *d],
        _ => return None,
    };
    let horizontal = |p: (f64, f64), q: (f64, f64)| p.1 == q.1 && p.0 != q.0;
    let vertical = |p: (f64, f64), q: (f64, f64)| p.0 == q.0 && p.1 != q.1;
    let edges = [(pts[0], pts[1]), (pts[1], pts[2]), (pts[2], pts[3]), (pts[3], pts[0])];
    let alternating = |first_h: bool| edges.iter().enumerate().all(|(i, &(p, q))| if (i % 2 == 0) == first_h { horizontal(p, q) } else { vertical(p, q) });
    if !(alternating(true) || alternating(false)) {
        return None;
    }
    let (xs, ys) = (pts.map(|p| p.0), pts.map(|p| p.1));
    let min = |v: [f64; 4]| v.iter().copied().fold(f64::INFINITY, f64::min);
    let max = |v: [f64; 4]| v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Some((min(xs), min(ys), max(xs), max(ys)))
}

fn rgb(c: Rgba) -> (f32, f32, f32) {
    let [r, g, b, _] = c.to_u8();
    (r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0)
}

/// Escapes a string for a PDF literal in the fonts' 8-bit encoding (ASCII; other characters as '?').
fn pdf_bytes(s: &str) -> Vec<u8> {
    s.chars().map(|c| if (c as u32) < 128 { c as u8 } else { b'?' }).collect()
}

struct Page<'a> {
    h: f32,
    content: Content,
    alphas: &'a mut Alphas,
    fonts_used: [bool; 2],
}

/// R's pdf() writes numbers with 2 decimals ("%.2f")
fn r2(v: f64) -> f32 {
    ((v * 100.0).round() / 100.0) as f32
}

impl Page<'_> {
    fn y(&self, y: f64) -> f32 {
        r2(self.h as f64 - y)
    }

    fn alpha(&mut self, fill: f32, stroke: f32) {
        if let Some(n) = self.alphas.name(fill, stroke) {
            self.content.set_parameters(Name(n.as_bytes()));
        }
    }

    fn stroke_state(&mut self, s: &Stroke) {
        let (r, g, b) = rgb(s.color);
        self.content.set_stroke_rgb(r, g, b);
        self.content.set_line_width(r2(s.width));
        self.content.set_line_cap(match s.cap {
            LineCap::Butt => LineCapStyle::ButtCap,
            LineCap::Round => LineCapStyle::RoundCap,
        });
        if s.dash.is_empty() {
            self.content.set_dash_pattern(std::iter::empty::<f32>(), 0.0);
        } else {
            self.content.set_dash_pattern(s.dash.iter().map(|&d| r2(d)), 0.0);
        }
    }

    fn path(&mut self, pts: &[(f64, f64)], close: bool) {
        let (x0, y0) = pts[0];
        let y0 = self.y(y0);
        self.content.move_to(r2(x0), y0);
        for &(x, y) in &pts[1..] {
            let yy = self.y(y);
            self.content.line_to(r2(x), yy);
        }
        if close {
            self.content.close_path();
        }
    }

    fn circle_path(&mut self, cx: f64, cy: f64, r: f64) {
        let k = 0.552_284_75 * r;
        let (cx, cy, r, k) = (r2(cx), self.y(cy), r2(r), r2(k));
        self.content.move_to(cx + r, cy);
        self.content.cubic_to(cx + r, cy + k, cx + k, cy + r, cx, cy + r);
        self.content.cubic_to(cx - k, cy + r, cx - r, cy + k, cx - r, cy);
        self.content.cubic_to(cx - r, cy - k, cx - k, cy - r, cx, cy - r);
        self.content.cubic_to(cx + k, cy - r, cx + r, cy - k, cx + r, cy);
        self.content.close_path();
    }

    fn item(&mut self, it: &Item) {
        match it {
            Item::Rect { x, y, w, h, fill } => {
                if fill.is_transparent() {
                    return;
                }
                self.content.save_state();
                self.alpha(fill.a, 1.0);
                let (r, g, b) = rgb(*fill);
                self.content.set_fill_rgb(r, g, b);
                let yy = self.y(y + h);
                self.content.rect(r2(*x), yy, r2(*w), r2(*h));
                self.content.fill_nonzero();
                self.content.restore_state();
            }
            Item::Polygon { points, fill, stroke } => {
                if points.len() < 2 {
                    return;
                }
                if let Some(f) = fill.filter(|f| !f.is_transparent()) {
                    self.content.save_state();
                    self.alpha(f.a, 1.0);
                    let (r, g, b) = rgb(f);
                    self.content.set_fill_rgb(r, g, b);
                    self.path(points, true);
                    self.content.fill_even_odd();
                    self.content.restore_state();
                }
                if let Some(s) = stroke.as_ref().filter(|s| !s.color.is_transparent() && s.width > 0.0) {
                    self.content.save_state();
                    self.alpha(1.0, s.color.a);
                    self.stroke_state(s);
                    self.content.set_line_join(pdf_writer::types::LineJoinStyle::MiterJoin);
                    // R's pdf() writes a rectangle (rectGrob: panel borders, tiles, keys) as "x y w h re", each number
                    // rounded to 2 decimals: an axis-aligned closed 4-point outline is written the same way, so that
                    // its far edges round as R's do (a path of corners rounds them separately, up to 0.01 pt off,
                    // which moves a border by a pixel where a renderer snaps strokes to pixels)
                    match axis_rect(points) {
                        Some((x0, y0, x1, y1)) => {
                            let yy = self.y(y1);
                            self.content.rect(r2(x0), yy, r2(x1 - x0), r2(y1 - y0));
                        }
                        None => self.path(points, true),
                    }
                    self.content.stroke();
                    self.content.restore_state();
                }
            }
            Item::Polyline { points, stroke } => {
                if points.len() < 2 || stroke.color.is_transparent() {
                    return;
                }
                self.content.save_state();
                self.alpha(1.0, stroke.color.a);
                self.stroke_state(stroke);
                self.content.set_line_join(pdf_writer::types::LineJoinStyle::RoundJoin);
                self.path(points, false);
                self.content.stroke();
                self.content.restore_state();
            }
            Item::Circle { cx, cy, r, fill, stroke } => {
                if let Some(f) = fill.filter(|f| !f.is_transparent()) {
                    self.content.save_state();
                    self.alpha(f.a, 1.0);
                    let (rr, g, b) = rgb(f);
                    self.content.set_fill_rgb(rr, g, b);
                    self.circle_path(*cx, *cy, *r);
                    self.content.fill_nonzero();
                    self.content.restore_state();
                }
                if let Some(s) = stroke.as_ref().filter(|s| !s.color.is_transparent() && s.width > 0.0) {
                    self.content.save_state();
                    self.alpha(1.0, s.color.a);
                    self.stroke_state(s);
                    self.circle_path(*cx, *cy, *r);
                    self.content.stroke();
                    self.content.restore_state();
                }
            }
            Item::Text { x, y, text, font, size, color, angle } => {
                if color.is_transparent() || text.trim().is_empty() {
                    return;
                }
                let fi = match font {
                    Font::SansBold => 1,
                    _ => 0,
                };
                self.fonts_used[fi] = true;
                self.content.save_state();
                self.alpha(color.a, 1.0);
                let (r, g, b) = rgb(*color);
                self.content.set_fill_rgb(r, g, b);
                self.content.begin_text();
                // R's pdf(): the font at 1 point scaled by the text matrix to the size rounded to whole points, numbers
                // with 2 decimals, and the AFM's kerning in a TJ array
                self.content.set_font(Name(if fi == 1 { b"F2" } else { b"F1" }), 1.0);
                let size_r = (size + 0.5).floor();
                let t = angle.to_radians();
                let (c, s) = (t.cos() * size_r, t.sin() * size_r);
                let yy = self.h as f64 - y;
                self.content.set_text_matrix([r2(c), r2(s), r2(-s), r2(c), r2(*x), r2(yy)]);
                let bytes = pdf_bytes(text);
                let face = if fi == 1 { crate::rdev::Face::Bold } else { crate::rdev::Face::Plain };
                let kerns: Vec<i16> = bytes.windows(2).map(|w| crate::rdev::afm_kern(w[0], w[1], face)).collect();
                if kerns.iter().all(|&k| k == 0) {
                    self.content.show(Str(&bytes));
                } else {
                    let mut tj = self.content.show_positioned();
                    let mut items = tj.items();
                    let mut start = 0;
                    for (i, &k) in kerns.iter().enumerate() {
                        if k != 0 {
                            items.show(Str(&bytes[start..=i]));
                            items.adjust(-(k as f32));
                            start = i + 1;
                        }
                    }
                    items.show(Str(&bytes[start..]));
                    items.finish();
                    tj.finish();
                }
                self.content.end_text();
                self.content.restore_state();
            }
            Item::Clip(c) => {
                // R's device: restore the page state, save it again, and clip
                self.content.restore_state();
                self.content.save_state();
                if let Some((x, y, w, h)) = c {
                    let yy = self.y(y + h);
                    self.content.rect(r2(*x), yy, r2(*w), r2(*h));
                    self.content.clip_nonzero();
                    self.content.end_path();
                }
            }
        }
    }
}

/// The PDF of a figure per page, each page `w` x `h` points.
pub fn to_r_pdf(figs: &[Figure], w: f64, h: f64) -> Result<Vec<u8>, String> {
    let mut pdf = Pdf::new();
    pdf.set_version(1, 4);
    let mut next = 1i32;
    let mut alloc = || {
        let r = Ref::new(next);
        next += 1;
        r
    };
    let catalog = alloc();
    let pages_id = alloc();
    let info = alloc();
    let font_ids = [alloc(), alloc()];
    let page_ids: Vec<Ref> = figs.iter().map(|_| alloc()).collect();
    let content_ids: Vec<Ref> = figs.iter().map(|_| alloc()).collect();
    let mut alphas = Alphas { names: BTreeMap::new() };
    let mut contents = Vec::new();
    for fig in figs {
        let mut page = Page { h: h as f32, content: Content::new(), alphas: &mut alphas, fonts_used: [false; 2] };
        page.content.save_state();
        // the background
        if !fig.background.is_transparent() {
            page.item(&Item::Rect { x: 0.0, y: 0.0, w, h, fill: fig.background });
        }
        for it in &fig.items {
            page.item(it);
        }
        page.content.restore_state();
        contents.push(page.content.finish());
    }
    let gs_ids: Vec<((u8, u8), String, Ref)> = alphas.names.iter().map(|(k, n)| (*k, n.clone(), alloc())).collect();

    pdf.catalog(catalog).pages(pages_id);
    pdf.pages(pages_id).kids(page_ids.iter().copied()).count(page_ids.len() as i32);
    pdf.document_info(info).title(TextStr("R Graphics Output")).producer(TextStr("hmf-plot-rs"));
    for (i, name) in [b"Helvetica".as_slice(), b"Helvetica-Bold".as_slice()].iter().enumerate() {
        let mut f = pdf.type1_font(font_ids[i]);
        f.base_font(Name(name));
        let mut enc = f.encoding_custom();
        enc.base_encoding(Name(b"WinAnsiEncoding"));
        let mut d = enc.differences();
        for (code, glyph) in DIFFERENCES {
            d.consecutive(code as u8, [Name(glyph.as_bytes())]);
        }
        d.finish();
        enc.finish();
        f.finish();
    }
    for (k, _, id) in &gs_ids {
        pdf.ext_graphics(*id).non_stroking_alpha(k.0 as f32 / 255.0).stroking_alpha(k.1 as f32 / 255.0);
    }
    for (i, data) in contents.iter().enumerate() {
        let mut p = pdf.page(page_ids[i]);
        p.parent(pages_id).media_box(Rect::new(0.0, 0.0, w as f32, h as f32)).contents(content_ids[i]);
        let mut res = p.resources();
        {
            let mut fonts = res.fonts();
            fonts.pair(Name(b"F1"), font_ids[0]);
            fonts.pair(Name(b"F2"), font_ids[1]);
        }
        if !gs_ids.is_empty() {
            let mut g = res.ext_g_states();
            for (_, n, id) in &gs_ids {
                g.pair(Name(n.as_bytes()), *id);
            }
        }
        res.finish();
        p.finish();
        let z = miniz_oxide::deflate::compress_to_vec_zlib(data, 6);
        pdf.stream(content_ids[i], &z).filter(Filter::FlateDecode);
    }
    Ok(pdf.finish())
}

/// Writes the PDF of the figures (a page each) to `path`.
pub fn write_r_pdf(path: &str, figs: &[Figure], w: f64, h: f64) -> Result<(), String> {
    let bytes = to_r_pdf(figs, w, h)?;
    std::fs::write(path, bytes).map_err(|e| format!("cannot write {path}: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangles_found_as_r_writes_them() {
        assert_eq!(axis_rect(&[(1.0, 2.0), (5.0, 2.0), (5.0, 7.0), (1.0, 7.0)]), Some((1.0, 2.0, 5.0, 7.0)));
        assert_eq!(axis_rect(&[(5.0, 7.0), (5.0, 2.0), (1.0, 2.0), (1.0, 7.0), (5.0, 7.0)]), Some((1.0, 2.0, 5.0, 7.0)));
        // not rectangles: a slanted edge, three corners, a degenerate edge
        assert_eq!(axis_rect(&[(1.0, 2.0), (5.0, 2.5), (5.0, 7.0), (1.0, 7.0)]), None);
        assert_eq!(axis_rect(&[(1.0, 2.0), (5.0, 2.0), (5.0, 7.0)]), None);
        assert_eq!(axis_rect(&[(1.0, 2.0), (1.0, 2.0), (5.0, 7.0), (1.0, 7.0)]), None);
    }

    #[test]
    fn deterministic_with_text_paths_and_alpha() {
        let mut f = Figure::new(100, 50, Rgba::WHITE);
        f.rect(1.0, 2.0, 10.0, 5.0, Rgba::rgb8(10, 20, 30));
        f.rect(1.0, 2.0, 10.0, 5.0, Rgba::rgba8(10, 20, 30, 50));
        f.line(0.0, 0.0, 50.0, 50.0, Stroke { color: Rgba::BLACK, width: 0.75, dash: vec![3.0, 3.0], cap: LineCap::Butt });
        f.circle(20.0, 20.0, 3.0, Some(Rgba::BLACK), None);
        f.clip(Some((0.0, 0.0, 40.0, 40.0)));
        f.text(5.0, 40.0, "Tumor - WARN", Font::Sans, 10.0, Rgba::BLACK, 90.0);
        f.clip(None);
        let a = to_r_pdf(&[f.clone(), f.clone()], 100.0, 50.0).unwrap();
        let b = to_r_pdf(&[f.clone(), f], 100.0, 50.0).unwrap();
        assert_eq!(a, b);
        let s = String::from_utf8_lossy(&a);
        assert!(s.contains("/BaseFont /Helvetica"));
        assert!(s.contains("/minus"));
        assert!(s.contains("/ca 0.19"));
    }
}
