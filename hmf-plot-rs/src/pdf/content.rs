//! Content streams with PDFBox's operators and number format: a [`Stream`] is one PDPageContentStream opened on a
//! page in append mode with a reset context (its operators run from the default graphics state); numbers are
//! written as PDFBox's NumberFormatUtil.formatFloatFast writes them (at most 5 decimals, half up), so that a
//! viewer computes the same positions from the same operators as from Java's file.

use std::collections::BTreeSet;

use super::{Document, FontId, ImageId, MissingGlyph, Rgb};

/// NumberFormatUtil.formatFloatFast(value, 5): the integer part, and up to 5 decimals rounded half up from the
/// float's exact value, trailing zeros left out
pub fn format_number(value: f32, out: &mut Vec<u8>) {
    const POW: i64 = 100_000;
    if !value.is_finite() {
        out.extend_from_slice(b"0");
        return;
    }
    let mut integer = value as i64;
    if value < 0.0 {
        out.push(b'-');
        integer = -integer;
    }
    let mut fraction = (((value as f64).abs() - integer as f64) * POW as f64 + 0.5) as i64;
    if fraction >= POW {
        integer += 1;
        fraction -= POW;
    }
    out.extend_from_slice(integer.to_string().as_bytes());
    if fraction > 0 {
        out.push(b'.');
        let digits = format!("{fraction:05}");
        out.extend_from_slice(digits.trim_end_matches('0').as_bytes());
    }
}

/// A content stream being written on one page
pub struct Stream {
    page: usize,
    buf: Vec<u8>,
    pub(crate) fonts: BTreeSet<FontId>,
    pub(crate) images: BTreeSet<ImageId>,
}

impl Stream {
    pub(crate) fn new(page: usize) -> Stream {
        Stream { page, buf: Vec::with_capacity(4096), fonts: BTreeSet::new(), images: BTreeSet::new() }
    }

    pub fn page(&self) -> usize {
        self.page
    }

    pub(crate) fn into_parts(self) -> (usize, Vec<u8>, BTreeSet<FontId>, BTreeSet<ImageId>) {
        (self.page, self.buf, self.fonts, self.images)
    }

    fn num(&mut self, v: f32) {
        format_number(v, &mut self.buf);
        self.buf.push(b' ');
    }

    fn op(&mut self, op: &str) {
        self.buf.extend_from_slice(op.as_bytes());
        self.buf.push(b'\n');
    }

    /// BT
    pub fn begin_text(&mut self) {
        self.op("BT");
    }

    /// ET
    pub fn end_text(&mut self) {
        self.op("ET");
    }

    /// Tf
    pub fn set_font(&mut self, font: FontId, size: f32) {
        self.fonts.insert(font);
        self.buf.extend_from_slice(format!("/F{} ", font + 1).as_bytes());
        self.num(size);
        self.op("Tf");
    }

    /// Td (newLineAtOffset)
    pub fn new_line_at_offset(&mut self, tx: f32, ty: f32) {
        self.num(tx);
        self.num(ty);
        self.op("Td");
    }

    /// Tj with the text encoded in the font's subset (2-byte glyph codes)
    pub fn show_text(&mut self, doc: &mut Document, font: FontId, text: &str) -> Result<(), MissingGlyph> {
        let codes = doc.fonts[font].encode(text)?;
        self.fonts.insert(font);
        self.buf.push(b'<');
        for b in codes {
            self.buf.extend_from_slice(format!("{b:02X}").as_bytes());
        }
        self.buf.extend_from_slice(b"> ");
        self.op("Tj");
        Ok(())
    }

    fn rgb(&mut self, c: Rgb) {
        // Color.getRGBColorComponents: component / 255f
        self.num(c.0 as f32 / 255f32);
        self.num(c.1 as f32 / 255f32);
        self.num(c.2 as f32 / 255f32);
    }

    /// setNonStrokingColor
    pub fn set_non_stroking_color(&mut self, c: Rgb) {
        self.rgb(c);
        self.op("rg");
    }

    /// setStrokingColor
    pub fn set_stroking_color(&mut self, c: Rgb) {
        self.rgb(c);
        self.op("RG");
    }

    /// addRect (re)
    pub fn add_rect(&mut self, x: f32, y: f32, w: f32, h: f32) {
        self.num(x);
        self.num(y);
        self.num(w);
        self.num(h);
        self.op("re");
    }

    /// f
    pub fn fill(&mut self) {
        self.op("f");
    }

    /// S
    pub fn stroke(&mut self) {
        self.op("S");
    }

    /// m
    pub fn move_to(&mut self, x: f32, y: f32) {
        self.num(x);
        self.num(y);
        self.op("m");
    }

    /// l
    pub fn line_to(&mut self, x: f32, y: f32) {
        self.num(x);
        self.num(y);
        self.op("l");
    }

    /// w
    pub fn set_line_width(&mut self, w: f32) {
        self.num(w);
        self.op("w");
    }

    /// J
    pub fn set_line_cap(&mut self, cap: i32) {
        self.buf.extend_from_slice(format!("{cap} ").as_bytes());
        self.op("J");
    }

    /// d with an empty dash array (solid lines)
    pub fn set_no_dash(&mut self) {
        self.op("[] 0 d");
    }

    /// drawImage(image, x, y, w, h): q, w 0 0 h x y cm, Do, Q
    pub fn draw_image(&mut self, image: ImageId, x: f32, y: f32, w: f32, h: f32) {
        self.images.insert(image);
        self.op("q");
        self.num(w);
        self.num(0.0);
        self.num(0.0);
        self.num(h);
        self.num(x);
        self.num(y);
        self.op("cm");
        self.buf.extend_from_slice(format!("/Im{} ", image + 1).as_bytes());
        self.op("Do");
        self.op("Q");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f(v: f32) -> String {
        let mut b = Vec::new();
        format_number(v, &mut b);
        String::from_utf8(b).unwrap()
    }

    #[test]
    fn numbers_as_pdfbox() {
        assert_eq!(f(0.0), "0");
        assert_eq!(f(11.0), "11");
        assert_eq!(f(61.7), "61.7");
        assert_eq!(f(771.88977), "771.88977");
        assert_eq!(f(242.0 / 255.0), "0.94902");
        assert_eq!(f(255.0 / 255.0), "1");
        assert_eq!(f(-1.5), "-1.5");
        assert_eq!(f(0.25), "0.25");
        assert_eq!(f(-0.000001), "-0");
        assert_eq!(f(99.999999), "100");
    }
}
