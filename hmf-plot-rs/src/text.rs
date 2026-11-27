//! Fonts and text metrics. The fonts are embedded and registered with BUI's text engine (FreeType with
//! light hinting, as Cairo and GD use it) under names of their own, so the system's fonts are never used and
//! output is the same on every machine.
//!
//! BUI renders a font at a whole number of pixels (`round(size)`), as R's Cairo device does for its ppem.

use std::ffi::CString;
use std::sync::{Mutex, MutexGuard, Once};

use crate::bui;

/// The faces the plots use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Font {
    /// DejaVu Sans Book: R's "sans" in the pipeline's image (fontconfig), so ggplot2's text.
    Sans,
    /// R's bold "sans" in the pipeline's images, which have DejaVu Sans Book only: fontconfig asks for synthetic
    /// emboldening, so the glyphs are DejaVu Sans Book emboldened (BUI's synthetic bold here) at the advances of
    /// the regular face (Pango keeps them). Added for CUPPA's and qsee's titles.
    SansBold,
    /// CMU Bright Bold: Circos's `bold` (ideogram labels).
    CircosBold,
    /// CMU Sans Serif Demi Condensed: Circos's `condensed` (text tracks).
    CircosCondensed,
    /// CMU Bright Roman: Circos's `default`, `sans`, `normal` and `light` (LINX's labels).
    CircosRoman,
}

impl Font {
    fn family(&self) -> &'static str {
        match self {
            Font::Sans | Font::SansBold => "hmfplot-dejavu-sans",
            Font::CircosBold => "hmfplot-cmu-bright-bold",
            Font::CircosCondensed => "hmfplot-cmu-sans-demicondensed",
            Font::CircosRoman => "hmfplot-cmu-bright-roman",
        }
    }

    pub(crate) fn data(&self) -> &'static [u8] {
        match self {
            Font::Sans | Font::SansBold => include_bytes!("../fonts/DejaVuSans.ttf"),
            Font::CircosBold => include_bytes!("../fonts/cmunbbx.otf"),
            Font::CircosCondensed => include_bytes!("../fonts/cmunssdc.otf"),
            Font::CircosRoman => include_bytes!("../fonts/cmunbmr.ttf"),
        }
    }

    pub(crate) fn c_family(&self) -> CString {
        CString::new(self.family()).unwrap()
    }

    /// BUI's style flags for the face: bold is synthesized from the registered regular face.
    pub(crate) fn bui_flags(&self) -> i32 {
        match self {
            Font::SansBold => bui::BUI_FONT_BOLD,
            _ => 0,
        }
    }

    /// The faces laid out with R's Cairo advances (DejaVu Sans, regular and synthetic bold).
    pub(crate) fn is_sans(&self) -> bool {
        matches!(self, Font::Sans | Font::SansBold)
    }

    /// The Circos font names of etc/fonts.conf that map to an embedded face.
    pub fn circos(name: &str) -> Option<Font> {
        match name {
            "bold" => Some(Font::CircosBold),
            "condensed" | "condensedbold" | "boldcondensed" => Some(Font::CircosCondensed),
            "default" | "sans" | "sans_serif" | "sansserif" | "normal" | "light" => Some(Font::CircosRoman),
            _ => None,
        }
    }
}

static REGISTER: Once = Once::new();
/// BUI's text engine (run cache, faces) is global and not thread-safe: every call goes through this lock.
static TEXT_LOCK: Mutex<()> = Mutex::new(());

/// Registers the embedded fonts (once) and takes the text engine's lock.
pub(crate) fn engine() -> MutexGuard<'static, ()> {
    let guard = TEXT_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    REGISTER.call_once(|| {
        for font in [Font::Sans, Font::CircosBold, Font::CircosCondensed, Font::CircosRoman] {
            let data = font.data();
            let family = font.c_family();
            // SAFETY: BUI copies the font data; the family string outlives the call.
            let ok = unsafe { bui::BUI_text_add_font(family.as_ptr(), 0, data.as_ptr().cast(), data.len() as i32) };
            assert!(ok != 0, "BUI rejected the embedded font {font:?}");
        }
    });
    guard
}

/// The rendered (hinted) extents of a text run at a pixel size, relative to its baseline origin.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Extents {
    /// pen advance (hinted, whole pixels)
    #[allow(dead_code)]
    pub advance: f64,
    /// ink box: left edge (from the origin), right edge, top (above the baseline, positive) and bottom
    /// (below the baseline, positive); all zero for blank text
    pub left: f64,
    pub right: f64,
    pub ascent: f64,
    pub descent: f64,
}

/// A text run's coverage mask with the pen origin inside it.
pub(crate) struct Mask {
    pub w: usize,
    pub h: usize,
    pub ox: i32,
    pub oy: i32,
    pub coverage: Vec<u8>,
    #[allow(dead_code)]
    pub advance: f64,
}

/// The run's mask (None for blank text). The caller holds the engine lock.
pub(crate) fn mask_locked(font: Font, size: f64, text: &str) -> Option<Mask> {
    if text.trim().is_empty() {
        return None;
    }
    let family = font.c_family();
    let mut advance = 0.0f32;
    // SAFETY: BUI reads `text.len()` bytes; the run stays valid while no font is opened (we hold the lock).
    unsafe {
        let run = bui::BUI_text_prepare(family.as_ptr(), font.bui_flags(), size as f32, text.as_ptr().cast(), text.len() as i32, &mut advance);
        if run.is_null() {
            return None;
        }
        let (mut w, mut h, mut ox, mut oy) = (0, 0, 0, 0);
        let m = bui::BUI_text_run_mask(run, &mut w, &mut h, &mut ox, &mut oy);
        if m.is_null() || w <= 0 || h <= 0 {
            return None;
        }
        let coverage = std::slice::from_raw_parts(m, (w * h) as usize).to_vec();
        Some(Mask { w: w as usize, h: h as usize, ox, oy, coverage, advance: advance as f64 })
    }
}

/// The pen advance of `text` at `size` pixels. DejaVu Sans (R's text) follows R's Cairo device: each glyph's
/// advance at the fractional size, rounded to whole pixels; the Circos faces take BUI's run advance (GD, like
/// BUI, lays a run out with hinted advances).
pub fn advance(font: Font, size: f64, text: &str) -> f64 {
    if font.is_sans() {
        return advance_fractional(font, size, text);
    }
    bui_advance(font, size, text)
}

/// The pen offsets of a run's glyphs as `advance` lays them out (None: BUI lays the run out itself).
pub(crate) fn glyph_offsets(font: Font, size: f64, text: &str) -> Option<Vec<(char, f64)>> {
    if !font.is_sans() {
        return None;
    }
    let h = SANS_HMTX.get_or_init(|| Hmtx::parse(Font::Sans.data()).expect("DejaVu Sans tables"));
    let mut pen = 0.0;
    Some(
        h.pen_steps(size, text)
            .into_iter()
            .map(|(c, step)| {
                let at = pen;
                pen += step;
                (c, at)
            })
            .collect(),
    )
}

fn bui_advance(font: Font, size: f64, text: &str) -> f64 {
    if text.is_empty() {
        return 0.0;
    }
    let _g = engine();
    let family = font.c_family();
    let mut w = 0.0f32;
    // SAFETY: BUI reads `text.len()` bytes.
    unsafe { bui::BUI_text_measure(family.as_ptr(), font.bui_flags(), size as f32, text.as_ptr().cast(), text.len() as i32, &mut w) };
    w as f64
}

/// The advance and the ink box of `text` at `size` pixels, from the glyphs as BUI renders them.
pub fn extents(font: Font, size: f64, text: &str) -> Extents {
    // measured before taking the (non-reentrant) engine lock
    let advance_px = advance(font, size, text);
    let _g = engine();
    let Some(m) = mask_locked(font, size, text) else {
        return Extents { advance: advance_px, ..Extents::default() };
    };
    let (mut x0, mut x1, mut y0, mut y1) = (usize::MAX, 0, usize::MAX, 0);
    for y in 0..m.h {
        for x in 0..m.w {
            if m.coverage[y * m.w + x] != 0 {
                x0 = x0.min(x);
                x1 = x1.max(x + 1);
                y0 = y0.min(y);
                y1 = y1.max(y + 1);
            }
        }
    }
    if x0 == usize::MAX {
        return Extents { advance: advance_px, ..Extents::default() };
    }
    Extents {
        advance: advance_px,
        left: x0 as f64 - m.ox as f64,
        right: x1 as f64 - m.ox as f64,
        ascent: m.oy as f64 - y0 as f64,
        descent: y1 as f64 - m.oy as f64,
    }
}

/// Points to pixels at a resolution.
pub fn pt_to_px(pt: f64, dpi: f64) -> f64 {
    pt * dpi / 72.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metrics_are_whole_pixels_and_sane() {
        // R 4.5's Cairo png device at 300 dpi, DejaVu Sans (grobAscent, grobDescent, grobWidth in pt / 0.2409)
        let e = extents(Font::Sans, 55.0, "M");
        assert_eq!((e.ascent, e.descent), (40.0, 1.0), "'M' at 13.2 pt: 9.636 pt up, 0.2409 pt down");
        for (size, m, d) in [(45.83, 34.0, 10.0), (36.67, 27.0, 8.0), (29.63, 22.0, 6.0)] {
            assert_eq!(extents(Font::Sans, size, "M").ascent, m);
            assert_eq!(extents(Font::Sans, size, "gjpqyQ").descent, d);
        }
        assert_eq!(advance(Font::Sans, 55.0, "6e+05"), 185.0, "R: 44.5665 pt = 185 px");
        let d = extents(Font::Sans, 55.0, "gjpqyQ");
        assert_eq!(d.descent, 12.0, "font_descent at 13.2 pt: 2.8908 pt = 12 px");
        assert_eq!(advance(Font::Sans, 55.0, "Copy Number PDF"), 503.0, "R: 121.1727 pt = 503 px");
        assert_eq!(advance(Font::Sans, 8.8 * 300.0 / 72.0, "1000000"), 161.0, "R: 38.78 pt = 161 px");
    }
}


/// Advance widths from a TrueType/OpenType font's own tables (head, hhea, hmtx, cmap format 4).
pub(crate) struct Hmtx {
    units_per_em: f64,
    advances: Vec<u16>,
    cmap: std::collections::HashMap<u32, u16>,
    /// the 'kern' table's format 0 pairs (glyph, glyph) -> adjustment in font units
    kern: std::collections::HashMap<(u16, u16), i16>,
}

impl Hmtx {
    pub(crate) fn parse(data: &[u8]) -> Option<Hmtx> {
        let u16_at = |o: usize| -> Option<u16> { Some(u16::from_be_bytes([*data.get(o)?, *data.get(o + 1)?])) };
        let u32_at = |o: usize| -> Option<u32> { Some(u32::from_be_bytes([*data.get(o)?, *data.get(o + 1)?, *data.get(o + 2)?, *data.get(o + 3)?])) };
        let num_tables = u16_at(4)? as usize;
        let mut tables = std::collections::HashMap::new();
        for i in 0..num_tables {
            let rec = 12 + 16 * i;
            let tag = std::str::from_utf8(data.get(rec..rec + 4)?).ok()?.to_string();
            tables.insert(tag, u32_at(rec + 8)? as usize);
        }
        let head = *tables.get("head")?;
        let units_per_em = u16_at(head + 18)? as f64;
        let hhea = *tables.get("hhea")?;
        let n_metrics = u16_at(hhea + 34)? as usize;
        let hmtx = *tables.get("hmtx")?;
        let advances = (0..n_metrics).map(|i| u16_at(hmtx + 4 * i)).collect::<Option<Vec<u16>>>()?;
        let cmap = *tables.get("cmap")?;
        let mut map = std::collections::HashMap::new();
        let n_sub = u16_at(cmap + 2)? as usize;
        for i in 0..n_sub {
            let (pid, eid) = (u16_at(cmap + 4 + 8 * i)?, u16_at(cmap + 6 + 8 * i)?);
            let off = cmap + u32_at(cmap + 8 + 8 * i)? as usize;
            if !(pid == 3 && (eid == 1 || eid == 10) || pid == 0) || u16_at(off)? != 4 {
                continue;
            }
            let segs = u16_at(off + 6)? as usize / 2;
            let (ends, starts, deltas, ranges) = (off + 14, off + 16 + 2 * segs, off + 16 + 4 * segs, off + 16 + 6 * segs);
            for s in 0..segs {
                let (end, start) = (u16_at(ends + 2 * s)? as u32, u16_at(starts + 2 * s)? as u32);
                let delta = u16_at(deltas + 2 * s)?;
                let ro = u16_at(ranges + 2 * s)? as usize;
                // printable ASCII, and the ligatures fi and fl (U+FB01, U+FB02) that Pango's shaping substitutes
                for c in (start..=end).filter(|&c| c <= 0x7e || c == 0xfb01 || c == 0xfb02) {
                    let g = if ro == 0 {
                        (c as u16).wrapping_add(delta)
                    } else {
                        let gi = u16_at(ranges + 2 * s + ro + 2 * (c - start) as usize)?;
                        if gi == 0 { 0 } else { gi.wrapping_add(delta) }
                    };
                    map.insert(c, g);
                }
            }
            break;
        }
        let mut kern = std::collections::HashMap::new();
        if let Some(&k) = tables.get("kern") {
            let n_sub = u16_at(k + 2)? as usize;
            let mut off = k + 4;
            for _ in 0..n_sub {
                let len = u16_at(off + 2)? as usize;
                let coverage = u16_at(off + 4)?;
                // format 0, horizontal
                if coverage >> 8 == 0 && coverage & 1 == 1 {
                    let n_pairs = u16_at(off + 6)? as usize;
                    for p in 0..n_pairs {
                        let rec = off + 14 + 6 * p;
                        kern.insert((u16_at(rec)?, u16_at(rec + 2)?), u16_at(rec + 4)? as i16);
                    }
                }
                off += len;
            }
        }
        Some(Hmtx { units_per_em, advances, cmap: map, kern })
    }

    fn glyph(&self, c: char) -> u16 {
        self.cmap.get(&(c as u32)).copied().unwrap_or(0)
    }

    /// Each glyph's advance plus its kerning with the next, at `size` pixels, rounded to whole pixels.
    fn pen_steps(&self, size: f64, s: &str) -> Vec<(char, f64)> {
        let chars: Vec<char> = s.chars().collect();
        (0..chars.len())
            .map(|i| {
                let g = self.glyph(chars[i]);
                let k = chars.get(i + 1).and_then(|&n| self.kern.get(&(g, self.glyph(n)))).copied().unwrap_or(0) as f64;
                (chars[i], ((self.advance_units(chars[i]) + k) * size / self.units_per_em).round())
            })
            .collect()
    }

    fn advance_units(&self, c: char) -> f64 {
        let g = self.cmap.get(&(c as u32)).copied().unwrap_or(0) as usize;
        *self.advances.get(g.min(self.advances.len() - 1)).unwrap_or(&0) as f64
    }

    /// The glyphs' kerned advances at `size` pixels, each rounded to whole pixels, summed (hinted metrics).
    pub(crate) fn width(&self, size: f64, s: &str) -> f64 {
        self.pen_steps(size, s).iter().map(|p| p.1).sum()
    }
}

static SANS_HMTX: std::sync::OnceLock<Hmtx> = std::sync::OnceLock::new();

/// The advance of `text` as R's Cairo device measures it: the font's advance widths at the fractional
/// pixel size, each rounded to whole pixels (hinted metrics). BUI's own measure works at a whole-pixel size.
pub fn advance_fractional(font: Font, size: f64, text: &str) -> f64 {
    match font {
        Font::Sans | Font::SansBold => SANS_HMTX.get_or_init(|| Hmtx::parse(Font::Sans.data()).expect("DejaVu Sans tables")).width(size, text),
        _ => bui_advance(font, size, text),
    }
}

#[cfg(test)]
mod hmtx_tests {
    use super::*;

    #[test]
    fn widths_as_r_measures_them() {
        // grobWidth in pt / 0.2409 from R 4.5's Cairo png device at 300 dpi
        let cases = [
            (13.2, "Copy Number PDF", 503.0),
            (11.0, "Copy Number PDF", 419.0),
            (8.8, "Copy Number PDF", 335.0),
            (13.2, "6e+05", 185.0),
            (11.0, "6e+05", 153.0),
            (8.8, "6e+05", 123.0),
            (13.2, "Purity", 160.0),
            (8.8, "Purity", 105.0),
            (11.0, "Kataegis Regions", 393.0),
            (11.0, "Ka", 57.0),
            (11.0, "Regions", 181.0),
            (11.0, "Purity", 133.0),
            (11.0, "AV", 59.0),
            (11.0, "To", 48.0),
            (11.0, "Mutation", 203.0),
            (11.0, "Relative", 185.0),
            (8.8, "100%", 104.0),
            (8.8, "1000000", 161.0),
        ];
        for (pt, s, px) in cases {
            let w = advance_fractional(Font::Sans, pt * 300.0 / 72.0, s);
            eprintln!("{pt} {s}: {w} (R {px})");
        }
        for (pt, s, px) in cases {
            let w = advance_fractional(Font::Sans, pt * 300.0 / 72.0, s);
            assert!((w - px).abs() < 0.5, "{pt} {s}: {w} vs R {px}");
        }
    }
}
