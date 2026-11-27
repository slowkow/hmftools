//! Renders a display list with BUI's software raster (anti-aliased signed-distance primitives, even-odd
//! polygons with 4 vertical subsamples per row, FreeType text) into RGB pixels.
//!
//! What BUI gives directly: filled rectangles, polygons, circles (filled or stroked), capsule lines, text at
//! whole-pixel positions, textured triangle meshes. What this module adds on top:
//! - butt-capped and dashed strokes (a dash is a segment; horizontal and vertical ones are rectangles);
//! - translucent strokes of several segments drawn once (into an offscreen BUI raster, then composited as a
//!   textured quad), so that joints are not blended twice;
//! - rotated text: the run's coverage mask from BUI's text engine as a texture on a rotated quad.

use crate::fmath::FMath;
use std::ffi::c_int;

use crate::bui;
use crate::color::Rgba;
use crate::display::{Figure, Item, LineCap, Stroke};
use crate::text;

struct Raster {
    r: *mut bui::BuiRaster,
    w: i32,
    h: i32,
}

impl Raster {
    fn new(w: i32, h: i32) -> Raster {
        // SAFETY: plain allocation; freed in Drop.
        let r = unsafe { bui::BUI_raster_create(w, h) };
        assert!(!r.is_null(), "BUI_raster_create failed");
        unsafe { bui::BUI_raster_begin_frame(r, w, h) };
        Raster { r, w, h }
    }

    fn pixels(&self) -> &[u32] {
        // SAFETY: the raster holds w * h pixels while it lives.
        unsafe { std::slice::from_raw_parts(bui::BUI_raster_pixels(self.r), (self.w * self.h) as usize) }
    }

    fn set_clip(&self, clip: Option<(i32, i32, i32, i32)>) {
        let (x, y, w, h) = clip.unwrap_or((0, 0, 0, 0));
        unsafe { bui::BUI_raster_set_clip(self.r, x, y, w, h) };
    }

    fn fill_rect(&self, x: f64, y: f64, w: f64, h: f64, c: Rgba) {
        if w > 0.0 && h > 0.0 && !c.is_transparent() {
            let rgba = c.to_array();
            unsafe { bui::BUI_raster_fill_rect(self.r, x as f32, y as f32, w as f32, h as f32, rgba.as_ptr()) };
        }
    }

    fn polygon(&self, pts: &[(f64, f64)], c: Rgba) {
        if pts.len() < 3 || c.is_transparent() {
            return;
        }
        let p: Vec<[f32; 2]> = pts.iter().map(|&(x, y)| [x as f32, y as f32]).collect();
        let rgba = c.to_array();
        unsafe { bui::BUI_raster_polygon(self.r, p.as_ptr(), p.len() as c_int, rgba.as_ptr()) };
    }

    fn circle(&self, cx: f64, cy: f64, radius: f64, line_width: f64, c: Rgba) {
        if c.is_transparent() || radius <= 0.0 {
            return;
        }
        let rgba = c.to_array();
        unsafe { bui::BUI_raster_circle(self.r, cx as f32, cy as f32, radius as f32, line_width as f32, rgba.as_ptr()) };
    }

    /// One straight segment of a stroke in continuous coordinates.
    fn segment(&self, a: (f64, f64), b: (f64, f64), width: f64, cap: LineCap, c: Rgba) {
        let hw = width / 2.0;
        let ext = if cap == LineCap::Round { hw } else { 0.0 };
        if a.1 == b.1 {
            let (x0, x1) = (a.0.min(b.0) - ext, a.0.max(b.0) + ext);
            self.fill_rect(x0, a.1 - hw, (x1 - x0).max(0.0), width, c);
        } else if a.0 == b.0 {
            let (y0, y1) = (a.1.min(b.1) - ext, a.1.max(b.1) + ext);
            self.fill_rect(a.0 - hw, y0, width, (y1 - y0).max(0.0), c);
        } else {
            // BUI's capsule runs between pixel centres: integer coordinates address the pixel's centre
            let rgba = c.to_array();
            let (x1, y1, x2, y2) = if cap == LineCap::Round {
                (a.0, a.1, b.0, b.1)
            } else {
                // butt ends: pull the capsule's round ends in by half the width (exact for thin lines)
                let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                let len = (dx * dx + dy * dy).sqrt();
                let k = if len > width { hw / len } else { 0.0 };
                (a.0 + dx * k, a.1 + dy * k, b.0 - dx * k, b.1 - dy * k)
            };
            unsafe {
                bui::BUI_raster_line(self.r, (x1 - 0.5) as f32, (y1 - 0.5) as f32, (x2 - 0.5) as f32, (y2 - 0.5) as f32, width as f32, rgba.as_ptr())
            };
        }
    }

    /// The segments of a polyline, split by its dash pattern.
    fn stroke_segments(&self, pts: &[(f64, f64)], s: &Stroke, c: Rgba) {
        if s.dash.is_empty() {
            for w in pts.windows(2) {
                self.segment(w[0], w[1], s.width, s.cap, c);
            }
            return;
        }
        let total: f64 = s.dash.iter().sum();
        if total <= 0.0 {
            return;
        }
        let (mut k, mut left, mut on) = (0usize, s.dash[0], true);
        for w in pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            let len = ((b.0 - a.0).fpowi(2) + (b.1 - a.1).fpowi(2)).sqrt();
            let mut t = 0.0;
            while t < len {
                let step = left.min(len - t);
                if on && step > 0.0 {
                    let p = |u: f64| (a.0 + (b.0 - a.0) * u / len, a.1 + (b.1 - a.1) * u / len);
                    self.segment(p(t), p(t + step), s.width, s.cap, c);
                }
                t += step;
                left -= step;
                if left <= 1e-9 {
                    k = (k + 1) % s.dash.len();
                    left = s.dash[k];
                    on = !on;
                }
            }
        }
    }

    /// Composites `cov` (coverage 0..=255 of a w x h box at x0, y0) in colour `c`, as a textured quad.
    fn composite(&self, x0: i32, y0: i32, w: usize, h: usize, cov: &[u8], c: Rgba) {
        let rgba: Vec<u8> = cov.iter().flat_map(|&a| [255u8, 255, 255, a]).collect();
        let col = u32::from_le_bytes(c.to_u8());
        unsafe {
            let tex = bui::BUI_texture_create(w as c_int, h as c_int, rgba.as_ptr().cast(), 0);
            if tex == 0 {
                return;
            }
            let (x0, y0, x1, y1) = (x0 as f32, y0 as f32, x0 as f32 + w as f32, y0 as f32 + h as f32);
            let v = [
                bui::BuiMeshVertex { x: x0, y: y0, u: 0.0, v: 0.0, col },
                bui::BuiMeshVertex { x: x1, y: y0, u: 1.0, v: 0.0, col },
                bui::BuiMeshVertex { x: x1, y: y1, u: 1.0, v: 1.0, col },
                bui::BuiMeshVertex { x: x0, y: y1, u: 0.0, v: 1.0, col },
            ];
            let idx: [u16; 6] = [0, 1, 2, 0, 2, 3];
            bui::BUI_raster_mesh(self.r, v.as_ptr(), 4, idx.as_ptr().cast(), 6, bui::BUI_MESH_NEAREST, tex, std::ptr::null(), 0.0, 0.0);
            bui::BUI_texture_free(tex);
        }
    }

    /// A text run. Faces laid out glyph by glyph (R's text: see `text::advance`) are drawn a glyph at a time at
    /// their pen offsets along the baseline; the others as one BUI run.
    fn text(&self, x: f64, y: f64, s: &str, font: text::Font, size: f64, c: Rgba, angle: f64) {
        if let Some(glyphs) = text::glyph_offsets(font, size, s) {
            let t = angle.to_radians();
            for (ch, pen) in glyphs {
                if ch == ' ' {
                    continue;
                }
                let mut buf = [0u8; 4];
                let g = ch.encode_utf8(&mut buf);
                self.run(x + pen * t.fcos(), y - pen * t.fsin(), g, font, size, c, angle);
            }
            return;
        }
        self.run(x, y, s, font, size, c, angle);
    }

    /// One BUI text run: unrotated through BUI's own text path (whole-pixel origin, as Cairo and GD place
    /// glyphs), rotated as the run's mask on a rotated quad. The caller holds the text engine lock.
    fn run(&self, x: f64, y: f64, s: &str, font: text::Font, size: f64, c: Rgba, angle: f64) {
        if c.is_transparent() {
            return;
        }
        let Some(m) = text::mask_locked(font, size, s) else { return };
        let rgba = c.to_array();
        let a = angle.rem_euclid(360.0);
        if a.abs() < 1e-9 {
            // BUI_raster_text_run places the pen at floor(x + 0.5), floor(y + 0.5)
            let family = font.c_family();
            unsafe {
                let run = bui::BUI_text_prepare(family.as_ptr(), font.bui_flags(), size as f32, s.as_ptr().cast(), s.len() as c_int, std::ptr::null_mut());
                if !run.is_null() {
                    raster_text_run(self.r, x as f32, y as f32, run, rgba.as_ptr());
                }
            }
            return;
        }
        let right_angle = (a / 90.0 - (a / 90.0).round()).abs() < 1e-9;
        // the quad: mask pixel (i, j) has its pen origin at (ox, oy); rotate about the origin (x, y)
        let (ox, oy) = if right_angle { (x.round(), y.round()) } else { (x, y) };
        let t = a.to_radians();
        let (cs, sn) = (t.fcos(), t.fsin());
        // counter-clockwise on screen (y down): (dx, dy) -> (dx cos + dy sin, -dx sin + dy cos)
        let place = |mx: f64, my: f64| {
            let (dx, dy) = (mx - m.ox as f64, my - m.oy as f64);
            ((ox + dx * cs + dy * sn) as f32, (oy - dx * sn + dy * cs) as f32)
        };
        let (w, h) = (m.w as f64, m.h as f64);
        let col = u32::from_le_bytes(c.to_u8());
        let rgba_tex: Vec<u8> = m.coverage.iter().flat_map(|&a| [255u8, 255, 255, a]).collect();
        unsafe {
            let tex = bui::BUI_texture_create(m.w as c_int, m.h as c_int, rgba_tex.as_ptr().cast(), 0);
            if tex == 0 {
                return;
            }
            let corners = [(0.0, 0.0, 0.0, 0.0), (w, 0.0, 1.0, 0.0), (w, h, 1.0, 1.0), (0.0, h, 0.0, 1.0)];
            let v: Vec<bui::BuiMeshVertex> = corners
                .iter()
                .map(|&(mx, my, u, vv)| {
                    let (px, py) = place(mx, my);
                    bui::BuiMeshVertex { x: px, y: py, u, v: vv, col }
                })
                .collect();
            let idx: [u16; 6] = [0, 1, 2, 0, 2, 3];
            let flags = if right_angle { bui::BUI_MESH_NEAREST } else { 0 };
            bui::BUI_raster_mesh(self.r, v.as_ptr(), 4, idx.as_ptr().cast(), 6, flags, tex, std::ptr::null(), 0.0, 0.0);
            bui::BUI_texture_free(tex);
        }
    }
}

unsafe extern "C" {
    #[link_name = "BUI_raster_text_run"]
    fn raster_text_run(r: *mut bui::BuiRaster, x: f32, y: f32, run: *const bui::BuiTextRun, rgba: *const f32);
}

impl Drop for Raster {
    fn drop(&mut self) {
        unsafe { bui::BUI_raster_free(self.r) };
    }
}

fn clip_box(c: Option<(f64, f64, f64, f64)>) -> Option<(i32, i32, i32, i32)> {
    c.map(|(x, y, w, h)| {
        // an empty clip (a viewport of zero width or height, e.g. the panel of an empty ggplot2 sub-table) shows
        // nothing, as with Cairo: a one-pixel clip outside the canvas (BUI takes a 0 x 0 clip for none). Added for
        // CUPPA's and qsee's plots; PURPLE's panels are never empty.
        if w < 1e-6 || h < 1e-6 {
            return (-2, -2, 1, 1);
        }
        let (x0, y0) = (x.round() as i32, y.round() as i32);
        let (x1, y1) = ((x + w).round() as i32, (y + h).round() as i32);
        (x0, y0, (x1 - x0).max(1), (y1 - y0).max(1))
    })
}

/// A translucent stroke with more than one segment: drawn opaque into `scratch`, then composited once.
fn layered_stroke(main: &Raster, scratch: &Raster, pts: &[(f64, f64)], s: &Stroke, clip: Option<(i32, i32, i32, i32)>, closed: bool) {
    let pad = s.width + 2.0;
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for &(x, y) in pts {
        x0 = x0.min(x);
        y0 = y0.min(y);
        x1 = x1.max(x);
        y1 = y1.max(y);
    }
    let bx0 = ((x0 - pad).floor() as i32).max(0);
    let by0 = ((y0 - pad).floor() as i32).max(0);
    let bx1 = ((x1 + pad).ceil() as i32).min(main.w);
    let by1 = ((y1 + pad).ceil() as i32).min(main.h);
    if bx1 <= bx0 || by1 <= by0 {
        return;
    }
    scratch.set_clip(None);
    scratch.fill_rect(bx0 as f64, by0 as f64, (bx1 - bx0) as f64, (by1 - by0) as f64, Rgba::BLACK);
    scratch.set_clip(clip);
    let mut path = pts.to_vec();
    if closed {
        path.push(pts[0]);
    }
    scratch.stroke_segments(&path, s, Rgba::WHITE);
    scratch.set_clip(None);
    let (w, h) = ((bx1 - bx0) as usize, (by1 - by0) as usize);
    let px = scratch.pixels();
    let mut cov = Vec::with_capacity(w * h);
    for y in by0..by1 {
        let row = &px[(y * scratch.w) as usize..];
        for x in bx0..bx1 {
            cov.push(((row[x as usize] >> 8) & 0xff) as u8);
        }
    }
    main.composite(bx0, by0, w, h, &cov, s.color);
}

fn stroke(main: &Raster, scratch: &mut Option<Raster>, pts: &[(f64, f64)], s: &Stroke, clip: Option<(i32, i32, i32, i32)>, closed: bool) {
    if s.color.is_transparent() || s.width <= 0.0 || pts.len() < 2 {
        return;
    }
    let segments = pts.len() - 1 + usize::from(closed);
    if s.color.a < 1.0 && (segments > 1 || !s.dash.is_empty()) {
        let sc = scratch.get_or_insert_with(|| Raster::new(main.w, main.h));
        layered_stroke(main, sc, pts, s, clip, closed);
        return;
    }
    let mut path = pts.to_vec();
    if closed {
        path.push(pts[0]);
    }
    main.stroke_segments(&path, s, s.color);
}

/// Renders the figure; returns its pixels as RGB bytes, rows top to bottom.
pub fn render(fig: &Figure) -> Vec<u8> {
    let _lock = text::engine();
    let (w, h) = (fig.width as i32, fig.height as i32);
    let main = Raster::new(w, h);
    let mut scratch: Option<Raster> = None;
    main.set_clip(None);
    main.fill_rect(0.0, 0.0, w as f64, h as f64, Rgba { a: 1.0, ..fig.background });
    let mut clip = None;
    for item in &fig.items {
        match item {
            Item::Rect { x, y, w, h, fill } => main.fill_rect(*x, *y, *w, *h, *fill),
            Item::Polygon { points, fill, stroke: st } => {
                if let Some(f) = fill {
                    main.polygon(points, *f);
                }
                if let Some(s) = st {
                    stroke(&main, &mut scratch, points, s, clip, true);
                }
            }
            Item::Polyline { points, stroke: s } => stroke(&main, &mut scratch, points, s, clip, false),
            Item::Circle { cx, cy, r, fill, stroke: st } => {
                if let Some(f) = fill {
                    main.circle(*cx, *cy, *r, 0.0, *f);
                }
                if let Some(s) = st {
                    main.circle(*cx, *cy, *r, s.width, s.color);
                }
            }
            Item::Text { x, y, text: t, font, size, color, angle } => main.text(*x, *y, t, *font, *size, *color, *angle),
            Item::Clip(c) => {
                clip = clip_box(*c);
                main.set_clip(clip);
            }
        }
    }
    unsafe { bui::BUI_raster_end_frame(main.r) };
    let mut out = Vec::with_capacity((w * h * 3) as usize);
    for p in main.pixels() {
        out.extend_from_slice(&[(p >> 16) as u8, (p >> 8) as u8, *p as u8]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A clip of zero height (an empty panel) hides what is drawn in it; a clip of some size shows its area.
    #[test]
    fn empty_clip_shows_nothing() {
        let mut fig = Figure::new(20, 20, Rgba::WHITE);
        fig.clip(Some((2.0, 10.0, 16.0, 0.0)));
        fig.rect(0.0, 0.0, 20.0, 20.0, Rgba::BLACK);
        fig.clip(None);
        assert!(render(&fig).iter().all(|&v| v == 255), "an empty clip let a pixel through");
        let mut fig = Figure::new(20, 20, Rgba::WHITE);
        fig.clip(Some((2.0, 10.0, 16.0, 4.0)));
        fig.rect(0.0, 0.0, 20.0, 20.0, Rgba::BLACK);
        let px = render(&fig);
        let dark = px.chunks(3).filter(|p| p[0] == 0).count();
        assert_eq!(dark, 16 * 4);
    }
}
