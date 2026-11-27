//! R's graphics devices as the hmftools R scripts use them, for layout: the device's resolution (grid's units in
//! device units) and its text metrics for the family "sans", which decide every text-sized row and column of a
//! ggplot2 or patchwork layout.
//!
//! - [`Device::Png`]: R's `png()` with the Cairo type, as the pipeline's images draw (`ggsave()` without ragg, or
//!   `png(res = 300)`): Pango lays text out with DejaVu Sans (fontconfig's "sans"; the images have no other sans
//!   face), with hinted metrics: each glyph's advance at the fractional pixel size rounded to whole pixels, kerning,
//!   and the ligatures fi and fl. Bold is DejaVu Sans emboldened by FreeType (fontconfig's synthetic bold): the same
//!   advances, glyphs taller by ceil(ppem / 24) pixels. Heights are the hinted ink ascent of "M" (GEStrHeight) and
//!   descents the hinted ink descent of the string (GEStrMetric).
//! - [`Device::Pdf`]: R's `pdf()`: Helvetica and Helvetica-Bold from R's AFM files, measured at the font size
//!   rounded to whole points, with kerning; the ascent of "M" is its glyph box top (718 / 1000 em).
//!
//! Every value here was compared with R 4.5.3 and 4.2.3 in the pipeline's qsee and CUPPA images (`tests`).

use crate::rdev_afm as afm;
use crate::text::{self, Font};

/// A font face of the family "sans".
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Face {
    Plain,
    Bold,
}

/// The device a figure is laid out for.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Device {
    /// R's Cairo png() at `dpi` pixels per inch
    Png { dpi: f64 },
    /// R's pdf(): device units are big points (1/72 inch)
    Pdf,
    /// pdf()'s text metrics in the units of a device of `dpi` (R measures text on the current device: a gtable
    /// built while `pdf(NULL)` is open keeps those sizes when it is drawn on a PNG)
    PdfAt { dpi: f64 },
}

impl Device {
    pub fn dpi(&self) -> f64 {
        match self {
            Device::Png { dpi } => *dpi,
            Device::Pdf => 72.0,
            Device::PdfAt { dpi } => *dpi,
        }
    }

    pub fn is_pdf(&self) -> bool {
        matches!(self, Device::Pdf)
    }

    /// pdf()'s metrics in big points scaled to the device's units
    fn afm_scale(&self) -> f64 {
        self.dpi() / 72.0
    }

    /// grid's "pt" (1/72.27 inch) in device units
    pub fn pt(&self, v: f64) -> f64 {
        v * self.dpi() / 72.27
    }

    /// grid's "bigpts" (1/72 inch), the unit of font sizes
    pub fn bigpt(&self, v: f64) -> f64 {
        v * self.dpi() / 72.0
    }

    pub fn mm(&self, v: f64) -> f64 {
        v * self.dpi() / 25.4
    }

    pub fn cm(&self, v: f64) -> f64 {
        v * self.dpi() / 2.54
    }

    pub fn inch(&self, v: f64) -> f64 {
        v * self.dpi()
    }

    /// grid's "lines" at a font size and line height: fontsize * lineheight big points
    pub fn lines(&self, n: f64, fontsize: f64, lineheight: f64) -> f64 {
        self.bigpt(n * fontsize * lineheight)
    }

    /// The display font of a face (the PDF writer maps it to Helvetica or Helvetica-Bold).
    pub fn font(&self, face: Face) -> Font {
        match face {
            Face::Plain => Font::Sans,
            Face::Bold => Font::SansBold,
        }
    }

    /// The text as the device lays it out: Pango substitutes DejaVu Sans's ligatures fi and fl; pdf() draws the
    /// characters as they are.
    pub fn shape(&self, s: &str) -> String {
        match self {
            Device::Png { .. } => s.replace("fi", "\u{fb01}").replace("fl", "\u{fb02}"),
            Device::Pdf | Device::PdfAt { .. } => s.to_string(),
        }
    }

    /// The font size in pixels of the PNG device (big points at the resolution).
    fn px(&self, fontsize: f64) -> f64 {
        fontsize * self.dpi() / 72.0
    }

    /// GEStrWidth of one line of text, in device units.
    pub fn str_width(&self, s: &str, fontsize: f64, face: Face) -> f64 {
        match self {
            Device::Png { .. } => text::advance(Font::Sans, self.px(fontsize), &self.shape(s)),
            Device::Pdf | Device::PdfAt { .. } => {
                let size = r_round_size(fontsize);
                afm_width(s, face) as f64 * size / 1000.0 * self.afm_scale()
            }
        }
    }

    /// The width of a (possibly multi-line) text: its widest line.
    pub fn text_width(&self, s: &str, fontsize: f64, face: Face) -> f64 {
        s.split('\n').map(|l| self.str_width(l, fontsize, face)).fold(0.0, f64::max)
    }

    /// GEStrHeight of one line: the ascent of "M" (device metric info).
    pub fn m_ascent(&self, fontsize: f64, face: Face) -> f64 {
        match self {
            Device::Png { .. } => {
                let px = self.px(fontsize);
                let a = text::extents(Font::Sans, px, "M").ascent;
                match face {
                    Face::Plain => a,
                    // FreeType's synthetic bold grows the outline by ppem / 24 pixels
                    Face::Bold => a + (px.round() / 24.0).ceil(),
                }
            }
            Device::Pdf | Device::PdfAt { .. } => {
                let size = r_round_size(fontsize);
                let i = (b'M' - 32) as usize;
                let asc = if face == Face::Bold { afm::BOLD_ASCENT[i] } else { afm::REGULAR_ASCENT[i] };
                asc as f64 * size / 1000.0 * self.afm_scale()
            }
        }
    }

    /// The descent of a string (GEStrMetric: the deepest of its characters), in device units.
    pub fn descent(&self, s: &str, fontsize: f64, face: Face) -> f64 {
        match self {
            Device::Png { .. } => text::extents(Font::Sans, self.px(fontsize), &self.shape(s)).descent,
            Device::Pdf | Device::PdfAt { .. } => {
                let size = r_round_size(fontsize);
                let table = if face == Face::Bold { &afm::BOLD_DESCENT } else { &afm::REGULAR_DESCENT };
                let d = s.bytes().filter(|c| (32..=126).contains(c)).map(|c| table[(c - 32) as usize]).max().unwrap_or(0);
                d as f64 * size / 1000.0 * self.afm_scale()
            }
        }
    }

    /// GEMetricInfo of one character (plotmath's glyph boxes): ascent, descent and width in device units. Cairo:
    /// the hinted ink extents and the advance; pdf(): the AFM glyph box and width at the rounded size.
    pub fn char_metric(&self, ch: char, fontsize: f64, face: Face) -> (f64, f64, f64) {
        let s = ch.to_string();
        match self {
            Device::Png { .. } => {
                let e = text::extents(Font::Sans, self.px(fontsize), &s);
                let asc = if face == Face::Bold && e.ascent > 0.0 { e.ascent + (self.px(fontsize).round() / 24.0).ceil() } else { e.ascent };
                (asc, e.descent, self.str_width(&s, fontsize, face))
            }
            Device::Pdf | Device::PdfAt { .. } => {
                let size = r_round_size(fontsize) * self.afm_scale();
                let c = ch as u32;
                if !(32..=126).contains(&c) {
                    return (0.0, 0.0, 0.0);
                }
                let i = (c - 32) as usize;
                let (asc, desc) = if face == Face::Bold { (afm::BOLD_ASCENT[i], afm::BOLD_DESCENT[i]) } else { (afm::REGULAR_ASCENT[i], afm::REGULAR_DESCENT[i]) };
                (asc as f64 * size / 1000.0, desc as f64 * size / 1000.0, afm_width(&s, face) as f64 * size / 1000.0)
            }
        }
    }

    /// ggplot2's font_descent(): the descent of "gjpqyQ".
    pub fn font_descent(&self, fontsize: f64, face: Face) -> f64 {
        self.descent("gjpqyQ", fontsize, face)
    }

    /// The distance between the baselines of two lines: lineheight x the device's character height (cra[1], 1.2
    /// times the font size).
    pub fn line_spacing(&self, fontsize: f64, lineheight: f64) -> f64 {
        self.bigpt(lineheight * 1.2 * fontsize)
    }

    /// GEStrHeight of a multi-line text: the ascent of "M" plus a line spacing per extra line.
    pub fn text_height(&self, s: &str, fontsize: f64, face: Face, lineheight: f64) -> f64 {
        let n = s.split('\n').count();
        (n - 1) as f64 * self.line_spacing(fontsize, lineheight) + self.m_ascent(fontsize, face)
    }
}

/// R's PostScript and PDF devices measure text at floor(size + 0.5) points.
fn r_round_size(fontsize: f64) -> f64 {
    (fontsize + 0.5).floor()
}

/// The kerning R applies between two characters of Helvetica (or Helvetica-Bold), in 1/1000 em.
pub(crate) fn afm_kern(a: u8, b: u8, face: Face) -> i16 {
    let kerning: &[(u8, u8, i16)] = if face == Face::Bold { &afm::BOLD_KERNING } else { &afm::REGULAR_KERNING };
    kerning.iter().find(|p| p.0 == a && p.1 == b).map(|p| p.2).unwrap_or(0)
}

/// The width of a string in 1/1000 em of Helvetica (or Helvetica-Bold), kerning included.
fn afm_width(s: &str, face: Face) -> i64 {
    let (widths, kerning): (&[u16; 95], &[(u8, u8, i16)]) =
        if face == Face::Bold { (&afm::BOLD_WIDTHS, &afm::BOLD_KERNING) } else { (&afm::REGULAR_WIDTHS, &afm::REGULAR_KERNING) };
    let b: Vec<u8> = s.bytes().collect();
    let mut w = 0i64;
    for (i, &c) in b.iter().enumerate() {
        w += if (32..=126).contains(&c) { widths[(c - 32) as usize] as i64 } else { 556 };
        if let Some(&n) = b.get(i + 1) {
            if let Some(k) = kerning.iter().find(|p| p.0 == c && p.1 == n) {
                w += k.2 as i64;
            }
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    /// R's grobWidth, grobHeight (the ascent of "M") and grobDescent in device units, from the pipeline's images
    /// (R 4.5.3 with Pango 1.56 in qsee's, R 4.2.3 with Pango 1.54 in CUPPA's give the same values).
    const R_METRICS: &[(&str, u8, f64, &str, f64, f64, f64)] = &[
        ("png", 1, 13.2, "M", 47.0, 40.0, 1.0),
        ("png", 2, 13.2, "M", 47.0, 43.0, 1.0),
        ("png", 1, 13.2, "Probabilities by classifier", 682.0, 40.0, 12.0),
        ("png", 2, 13.2, "Probabilities by classifier", 682.0, 43.0, 12.0),
        ("png", 1, 13.2, "fi", 35.0, 40.0, 1.0),
        ("png", 1, 13.2, "classifier", 244.0, 40.0, 2.0),
        ("png", 1, 13.2, "AV", 72.0, 40.0, 1.0),
        ("png", 1, 7.68, "M", 28.0, 23.0, 0.0),
        ("png", 2, 7.68, "M", 28.0, 25.0, 0.0),
        ("png", 1, 7.68, "gjpqyQ", 113.0, 23.0, 7.0),
        ("png", 1, 8.8, "M", 32.0, 27.0, 0.0),
        ("png", 2, 8.8, "M", 32.0, 29.0, 0.0),
        ("png", 1, 8.8, "gjpqyQ", 130.0, 27.0, 8.0),
        ("png", 1, 9.45, "M", 34.0, 29.0, 0.0),
        ("png", 2, 9.45, "M", 34.0, 31.0, 0.0),
        ("png", 1, 9.45, "gjpqyQ", 140.0, 29.0, 8.0),
        ("png", 1, 10.5, "M", 38.0, 33.0, 0.0),
        ("png", 2, 10.5, "M", 38.0, 35.0, 0.0),
        ("png", 1, 10.5, "gjpqyQ", 156.0, 33.0, 9.0),
        ("png", 1, 11.0, "M", 40.0, 34.0, 0.0),
        ("png", 2, 11.0, "M", 40.0, 36.0, 0.0),
        ("png", 1, 11.0, "gjpqyQ", 163.0, 34.0, 10.0),
        ("png", 1, 14.7, "M", 53.0, 45.0, 0.0),
        ("png", 2, 14.7, "M", 53.0, 48.0, 0.0),
        ("png", 1, 14.7, "gjpqyQ", 218.0, 45.0, 13.0),
        ("pdf", 1, 13.2, "M", 10.829, 9.334, 0.0),
        ("pdf", 1, 13.2, "gjpqyQ", 41.184, 9.334, 2.86),
        ("pdf", 1, 13.2, "Probabilities by classifier", 142.051, 9.334, 2.782),
        ("pdf", 1, 13.2, "fi", 6.5, 9.334, 0.0),
        ("pdf", 1, 13.2, "AV", 16.432, 9.334, 0.0),
        ("pdf", 1, 14.7, "Probabilities by classifier", 163.905, 10.77, 3.21),
        ("pdf", 1, 10.5, "gjpqyQ", 34.848, 7.898, 2.42),
        ("pdf", 2, 13.2, "gjpqyQ", 44.785, 9.334, 2.821),
        ("pdf", 2, 13.2, "Probabilities by classifier", 156.013, 9.334, 2.782),
        ("pdf", 2, 14.7, "fi", 9.165, 10.77, 0.0),
        ("pdf", 2, 10.5, "AV", 14.399, 7.898, 0.0),
    ];

    /// R's widths (1/1000 em) on pdf(): kerning pairs with a space and the last pair of each first character
    /// (y period, P period) are not applied, the others are; '-' is the minus sign.
    #[test]
    fn pdf_kerning_as_r_applies_it() {
        let r = [
            ("y,", 678), ("y.", 778), ("v.", 778), ("v,", 698), ("w.", 1000), ("r.", 561), ("P.", 945), ("F.", 739), ("V.", 820),
            ("T.", 769), ("f.", 526), ("Y.", 805), ("T:", 869), ("V:", 905), ("LY", 1083), ("Lyy", 1526), ("L'", 618), ("r-", 917),
            ("Ty", 991), ("rv", 863), (": ", 556), (". ", 556), (" T", 889), ("AV", 1264), ("Fe", 1137), ("P,", 765),
            ("Non-small cell: LUSC", 9752),
        ];
        for (s, w) in r {
            assert_eq!(afm_width(s, Face::Plain), w, "{s}");
        }
        assert_eq!(afm_width("AV", Face::Bold), 1309);
        assert_eq!(afm_width("Ty", Face::Bold), 1107);
    }

    #[test]
    fn metrics_as_r_measures_them() {
        for &(dev, face, size, s, w, h, d) in R_METRICS {
            let device = if dev == "png" { Device::Png { dpi: 300.0 } } else { Device::Pdf };
            let face = if face == 2 { Face::Bold } else { Face::Plain };
            let tol = if dev == "png" { 0.5 } else { 1e-6 };
            assert!((device.str_width(s, size, face) - w).abs() < tol, "{dev} {face:?} {size} {s}: width {} vs R {w}", device.str_width(s, size, face));
            assert!((device.m_ascent(size, face) - h).abs() < tol, "{dev} {face:?} {size}: M ascent {} vs R {h}", device.m_ascent(size, face));
            assert!((device.descent(s, size, face) - d).abs() < tol, "{dev} {face:?} {size} {s}: descent {} vs R {d}", device.descent(s, size, face));
        }
    }
}
