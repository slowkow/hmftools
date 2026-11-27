//! TrueType fonts for the PDF writer: metrics as PDFBox computes them for a `PDType0Font` (so that text widths,
//! and with them every table layout, equal the Java tools'), the glyphs a document uses, and the embedded
//! subset (a Type0 font with Identity-H encoding over a CIDFontType2, its widths and a ToUnicode map).

use std::collections::BTreeMap;

use subsetter::GlyphRemapper;

/// The error PDFBox raises for a character the font has no glyph for (IllegalArgumentException)
#[derive(Debug, Clone, PartialEq)]
pub struct MissingGlyph {
    pub ch: char,
    pub font: String,
}

impl std::fmt::Display for MissingGlyph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "No glyph for U+{:04X} in font {}", self.ch as u32, self.font)
    }
}

impl std::error::Error for MissingGlyph {}

/// A TrueType font, its PDFBox metrics, and the glyphs used so far
pub struct Font {
    data: Vec<u8>,
    /// the PostScript name (name ID 6), PDFont.getName() before subsetting
    name: String,
    units_per_em: u16,
    /// char -> glyph id, from the font's Unicode cmap
    cmap: BTreeMap<char, u16>,
    /// advance widths in 1000 units of text space, rounded as PDFBox writes the CIDFont's /W array
    widths: Vec<f32>,
    ascent: f32,
    descent: f32,
    x_height: f32,
    cap_height: f32,
    italic_angle: f32,
    bbox: [f32; 4],
    average_width: f32,
    /// glyphs used by the document: old glyph id -> new (subset) glyph id
    remapper: GlyphRemapper,
    /// new glyph id -> the character it was used for (for ToUnicode)
    used: BTreeMap<u16, char>,
}

impl Font {
    /// Reads a TrueType font (.ttf). Metrics follow PDFBox 3 (TrueTypeEmbedder / PDCIDFontType2Embedder):
    /// scaling 1000 / unitsPerEm in float; ascent and descent from 'hhea'; x-height and cap height from 'OS/2'
    /// (version 2 and later; else the bounding boxes of 'x' and 'H'); widths rounded to integers.
    pub fn from_bytes(data: Vec<u8>) -> Result<Font, String> {
        let face = ttf_parser::Face::parse(&data, 0).map_err(|e| format!("cannot parse font: {e}"))?;
        let units_per_em = face.units_per_em();
        let scaling = 1000f32 / units_per_em as f32;
        let name = face
            .names()
            .into_iter()
            .find(|n| n.name_id == ttf_parser::name_id::POST_SCRIPT_NAME && n.is_unicode())
            .and_then(|n| n.to_string())
            .or_else(|| {
                face.names()
                    .into_iter()
                    .find(|n| n.name_id == ttf_parser::name_id::POST_SCRIPT_NAME)
                    .map(|n| String::from_utf8_lossy(n.name).to_string())
            })
            .unwrap_or_else(|| "Font".to_string());
        let mut cmap = BTreeMap::new();
        if let Some(table) = face.tables().cmap {
            for sub in table.subtables.into_iter().filter(|s| s.is_unicode()) {
                sub.codepoints(|cp| {
                    if let (Some(ch), Some(g)) = (char::from_u32(cp), sub.glyph_index(cp)) {
                        cmap.entry(ch).or_insert(g.0);
                    }
                });
            }
        }
        let n = face.number_of_glyphs();
        let widths: Vec<f32> = (0..n)
            .map(|g| {
                let adv = face.glyph_hor_advance(ttf_parser::GlyphId(g)).unwrap_or(0) as f32;
                // Math.round(float) on the scaled advance
                java_round_f32(adv * scaling) as f32
            })
            .collect();
        let ascent = face.ascender() as f32 * scaling;
        let descent = face.descender() as f32 * scaling;
        let (x_height, cap_height) = match (face.x_height(), face.capital_height()) {
            (Some(x), Some(c)) => (x as f32 * scaling, c as f32 * scaling),
            _ => {
                let gh = |ch: char| {
                    face.glyph_index(ch).and_then(|g| face.glyph_bounding_box(g)).map(|b| b.y_max as f32 * scaling).unwrap_or(0.0)
                };
                (gh('x'), gh('H'))
            }
        };
        let gb = face.global_bounding_box();
        let bbox = [gb.x_min as f32 * scaling, gb.y_min as f32 * scaling, gb.x_max as f32 * scaling, gb.y_max as f32 * scaling];
        let italic_angle = face.italic_angle();
        // PDCIDFont.getAverageFontWidth: the mean of the positive /W widths
        let (sum, count) = widths.iter().filter(|w| **w > 0.0).fold((0f32, 0u32), |(s, c), w| (s + w, c + 1));
        let average_width = if count > 0 { sum / count as f32 } else { 1000.0 };
        drop(face);
        Ok(Font {
            data,
            name,
            units_per_em,
            cmap,
            widths,
            ascent,
            descent,
            x_height,
            cap_height,
            italic_angle,
            bbox,
            average_width,
            remapper: GlyphRemapper::new(),
            used: BTreeMap::new(),
        })
    }

    /// PDFont.getName()
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn units_per_em(&self) -> u16 {
        self.units_per_em
    }
    /// PDFontDescriptor.getAscent() (1000 units)
    pub fn ascent(&self) -> f32 {
        self.ascent
    }
    /// PDFontDescriptor.getDescent() (1000 units, negative)
    pub fn descent(&self) -> f32 {
        self.descent
    }
    /// PDFontDescriptor.getXHeight() (1000 units)
    pub fn x_height(&self) -> f32 {
        self.x_height
    }
    pub fn cap_height(&self) -> f32 {
        self.cap_height
    }
    /// PDFont.getAverageFontWidth() (1000 units)
    pub fn average_width(&self) -> f32 {
        self.average_width
    }

    /// The glyph of a character, or PDFBox's error
    pub fn glyph(&self, ch: char) -> Result<u16, MissingGlyph> {
        self.cmap.get(&ch).copied().ok_or_else(|| MissingGlyph { ch, font: self.name.clone() })
    }

    /// Whether every character of `text` has a glyph
    pub fn can_encode(&self, text: &str) -> bool {
        text.chars().all(|c| self.cmap.contains_key(&c))
    }

    /// The width of one character in 1000 units (PDFont.getWidth of its code)
    pub fn char_width(&self, ch: char) -> Result<f32, MissingGlyph> {
        Ok(self.widths[self.glyph(ch)? as usize])
    }

    /// PDFont.getStringWidth: the float sum of the characters' widths, in 1000 units of text space
    pub fn string_width(&self, text: &str) -> Result<f32, MissingGlyph> {
        let mut w = 0f32;
        for ch in text.chars() {
            w += self.char_width(ch)?;
        }
        Ok(w)
    }

    /// Marks the glyphs of `text` as used and returns the encoded string (2-byte codes of the subset's
    /// glyph ids, Identity-H)
    pub(crate) fn encode(&mut self, text: &str) -> Result<Vec<u8>, MissingGlyph> {
        let mut out = Vec::with_capacity(text.len() * 2);
        for ch in text.chars() {
            let old = self.glyph(ch)?;
            let new = self.remapper.remap(old);
            self.used.entry(new).or_insert(ch);
            out.extend_from_slice(&new.to_be_bytes());
        }
        Ok(out)
    }

    pub(crate) fn is_used(&self) -> bool {
        !self.used.is_empty()
    }

    /// The subset font program and its glyphs: (subset bytes, [(new gid, width, char)])
    pub(crate) fn subset(&self) -> Result<(Vec<u8>, Vec<(u16, f32, char)>), String> {
        let bytes = subsetter::subset(&self.data, 0, &self.remapper).map_err(|e| format!("font subset of {}: {e}", self.name))?;
        let mut glyphs = Vec::new();
        for old in self.remapper.remapped_gids() {
            let new = self.remapper.get(old).unwrap_or(0);
            if let Some(ch) = self.used.get(&new) {
                glyphs.push((new, self.widths[old as usize], *ch));
            }
        }
        glyphs.sort_by_key(|g| g.0);
        Ok((bytes, glyphs))
    }

    /// A six-letter subset tag derived from the glyphs used (as PDFBox's, it only has to differ between subsets)
    pub(crate) fn subset_tag(&self) -> String {
        let mut h: u64 = 0xcbf29ce484222325;
        for (g, c) in &self.used {
            for b in g.to_le_bytes().iter().chain((*c as u32).to_le_bytes().iter()) {
                h ^= *b as u64;
                h = h.wrapping_mul(0x100000001b3);
            }
        }
        (0..6)
            .map(|i| {
                let v = (h >> (i * 5)) % 26;
                (b'A' + v as u8) as char
            })
            .collect()
    }

    pub(crate) fn italic_angle(&self) -> f32 {
        self.italic_angle
    }
    pub(crate) fn bbox(&self) -> [f32; 4] {
        self.bbox
    }
}

/// java.lang.Math.round(float): floor(x + 0.5) as an int
fn java_round_f32(x: f32) -> i32 {
    if x.is_nan() {
        return 0;
    }
    (x + 0.5).floor() as i32
}

/// The font metrics boxable derives from a font descriptor (FontUtils.createFontMetrics), per point of font size
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoxMetrics {
    pub ascent: f32,
    pub descent: f32,
    pub height: f32,
}

impl BoxMetrics {
    /// base = xHeight / 1000; ascent = ascent / 1000 - base; descent = descent / 1000; height = base + ascent - descent
    pub fn of(font: &Font) -> BoxMetrics {
        let base = font.x_height() / 1000f32;
        let ascent = font.ascent() / 1000f32 - base;
        let descent = font.descent() / 1000f32;
        BoxMetrics { ascent, descent, height: base + ascent - descent }
    }
    /// FontUtils.getHeight(font, fontSize)
    pub fn height(&self, size: f32) -> f32 {
        self.height * size
    }
    /// FontUtils.getDescent(font, fontSize) (negative)
    pub fn descent(&self, size: f32) -> f32 {
        self.descent * size
    }
    /// FontUtils.getAscent(font, fontSize)
    pub fn ascent(&self, size: f32) -> f32 {
        self.ascent * size
    }
}
