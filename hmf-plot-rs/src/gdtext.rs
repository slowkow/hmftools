//! The string box GD's stringFT returns (libgd 2.3, gdft.c), from which Circos sizes and places its labels:
//! glyph metrics at 300 dpi (METRIC_RES) grid-fitted to whole pixels there, the box from the first glyph's
//! left bearing to the last glyph's advance and from the highest glyph top to the lowest glyph bottom,
//! scaled to 96 dpi and truncated to integers. Computed from a TrueType font's own tables (cmap, hmtx, kern,
//! loca, glyf). (Added for LINX's Circos labels in CMU Bright Roman.)

use std::collections::HashMap;

pub(crate) struct TtfMetrics {
    units_per_em: f64,
    advances: Vec<u16>,
    /// glyph bounding boxes (xMin, yMin, xMax, yMax) in font units
    bboxes: Vec<(i16, i16, i16, i16)>,
    cmap: HashMap<u32, u16>,
    kern: HashMap<(u16, u16), i16>,
}

impl TtfMetrics {
    pub(crate) fn parse(data: &[u8]) -> Option<TtfMetrics> {
        let u16_at = |o: usize| -> Option<u16> { Some(u16::from_be_bytes([*data.get(o)?, *data.get(o + 1)?])) };
        let u32_at = |o: usize| -> Option<u32> { Some(u32::from_be_bytes([*data.get(o)?, *data.get(o + 1)?, *data.get(o + 2)?, *data.get(o + 3)?])) };
        let mut tables = HashMap::new();
        for i in 0..u16_at(4)? as usize {
            let rec = 12 + 16 * i;
            tables.insert(std::str::from_utf8(data.get(rec..rec + 4)?).ok()?.to_string(), u32_at(rec + 8)? as usize);
        }
        let head = *tables.get("head")?;
        let units_per_em = u16_at(head + 18)? as f64;
        let long_loca = u16_at(head + 50)? == 1;
        let n_glyphs = u16_at(*tables.get("maxp")? + 4)? as usize;
        let hhea = *tables.get("hhea")?;
        let n_metrics = u16_at(hhea + 34)? as usize;
        let hmtx = *tables.get("hmtx")?;
        let advances = (0..n_metrics).map(|i| u16_at(hmtx + 4 * i)).collect::<Option<Vec<u16>>>()?;
        let (loca, glyf) = (*tables.get("loca")?, *tables.get("glyf")?);
        let mut bboxes = Vec::with_capacity(n_glyphs);
        for g in 0..n_glyphs {
            let (a, b) = if long_loca { (u32_at(loca + 4 * g)? as usize, u32_at(loca + 4 * g + 4)? as usize) } else { (2 * u16_at(loca + 2 * g)? as usize, 2 * u16_at(loca + 2 * g + 2)? as usize) };
            if b <= a {
                bboxes.push((0, 0, 0, 0));
                continue;
            }
            let o = glyf + a;
            bboxes.push((u16_at(o + 2)? as i16, u16_at(o + 4)? as i16, u16_at(o + 6)? as i16, u16_at(o + 8)? as i16));
        }
        let cmap_t = *tables.get("cmap")?;
        let mut cmap = HashMap::new();
        for i in 0..u16_at(cmap_t + 2)? as usize {
            let (pid, eid) = (u16_at(cmap_t + 4 + 8 * i)?, u16_at(cmap_t + 6 + 8 * i)?);
            let off = cmap_t + u32_at(cmap_t + 8 + 8 * i)? as usize;
            if !(pid == 3 && (eid == 1 || eid == 10) || pid == 0) || u16_at(off)? != 4 {
                continue;
            }
            let segs = u16_at(off + 6)? as usize / 2;
            let (ends, starts, deltas, ranges) = (off + 14, off + 16 + 2 * segs, off + 16 + 4 * segs, off + 16 + 6 * segs);
            for s in 0..segs {
                let (end, start) = (u16_at(ends + 2 * s)? as u32, u16_at(starts + 2 * s)? as u32);
                let delta = u16_at(deltas + 2 * s)?;
                let ro = u16_at(ranges + 2 * s)? as usize;
                for c in start..=end.min(0xff) {
                    let g = if ro == 0 {
                        (c as u16).wrapping_add(delta)
                    } else {
                        let gi = u16_at(ranges + 2 * s + ro + 2 * (c - start) as usize)?;
                        if gi == 0 { 0 } else { gi.wrapping_add(delta) }
                    };
                    cmap.insert(c, g);
                }
            }
            break;
        }
        let mut kern = HashMap::new();
        if let Some(&k) = tables.get("kern") {
            let mut off = k + 4;
            for _ in 0..u16_at(k + 2)? as usize {
                let len = u16_at(off + 2)? as usize;
                let coverage = u16_at(off + 4)?;
                if coverage >> 8 == 0 && coverage & 1 == 1 {
                    for p in 0..u16_at(off + 6)? as usize {
                        let rec = off + 14 + 6 * p;
                        kern.insert((u16_at(rec)?, u16_at(rec + 2)?), u16_at(rec + 4)? as i16);
                    }
                }
                off += len;
            }
        }
        Some(TtfMetrics { units_per_em, advances, bboxes, cmap, kern })
    }

    /// GD's string box at `size_pt` points: (width, height) as Circos's get_label_wh_from_bounds gives them
    /// (the box's extent minus one pixel).
    pub(crate) fn gd_label_size(&self, size_pt: f64, text: &str) -> (f64, f64) {
        // 26.6 units at 300 dpi
        let scale = size_pt * 300.0 / 72.0 / self.units_per_em * 64.0;
        let px = |v: f64| v * scale;
        let (mut pen, mut min_x, mut max_x, mut min_y, mut max_y) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut prev: Option<u16> = None;
        for (i, c) in text.chars().enumerate() {
            let g = self.cmap.get(&(c as u32)).copied().unwrap_or(0);
            if let Some(p) = prev {
                // FT_Get_Kerning, ft_kerning_default: rounded to whole pixels
                let k = self.kern.get(&(p, g)).copied().unwrap_or(0) as f64;
                pen += (px(k) / 64.0).round() * 64.0;
            }
            let (x0, y0, _x1, y1) = self.bboxes.get(g as usize).copied().unwrap_or((0, 0, 0, 0));
            let adv = *self.advances.get((g as usize).min(self.advances.len() - 1)).unwrap_or(&0) as f64;
            // grid-fitted metrics: bearings and extents to whole pixels, the advance rounded
            let bearing_x = (px(x0 as f64) / 64.0).floor() * 64.0;
            let top = (px(y1 as f64) / 64.0).ceil() * 64.0;
            let bottom = (px(y0 as f64) / 64.0).floor() * 64.0;
            let advance = (px(adv) / 64.0).round() * 64.0;
            let gmin_x = pen + bearing_x;
            let gmin_y = -top;
            let gmax_x = pen + advance;
            let gmax_y = gmin_y + (top - bottom);
            if i == 0 {
                (min_x, min_y, max_x, max_y) = (gmin_x, gmin_y, gmax_x, gmax_y);
            } else {
                min_x = min_x.min(gmin_x);
                min_y = min_y.min(gmin_y);
                max_x = max_x.max(gmax_x);
                max_y = max_y.max(gmax_y);
            }
            pen += advance;
            prev = Some(g);
        }
        // to 96 dpi pixels, truncated as the int brect is
        let s = 96.0 / (64.0 * 300.0);
        let t = |v: f64| (v * s).trunc();
        let w = (t(max_x) - t(min_x)).abs() - 1.0;
        let h = (t(min_y) - t(max_y)).abs() - 1.0;
        (w, h)
    }
}

static ROMAN: std::sync::OnceLock<Option<TtfMetrics>> = std::sync::OnceLock::new();

/// GD's label size (Circos's get_label_size) for the embedded TrueType faces; None for the others.
pub fn gd_label_size(font: crate::text::Font, size_pt: f64, text: &str) -> Option<(f64, f64)> {
    match font {
        crate::text::Font::CircosRoman => ROMAN.get_or_init(|| TtfMetrics::parse(font.data())).as_ref().map(|m| m.gd_label_size(size_pt, text)),
        _ => None,
    }
}
