//! A PDF writer for reports drawn directly in PDF coordinates (points, origin bottom left): pages of any size,
//! filled and stroked rectangles and lines, text in embedded TrueType font subsets (Type0, Identity-H, with a
//! ToUnicode map so the text can be extracted), PNG images (with soft masks), document outline (bookmarks) and
//! document information. [`table`] adds a table layout engine that reproduces boxable's, the library ORANGE's
//! Java report draws its tables with.
//!
//! The file structure is written by `pdf-writer`, font subsets by `subsetter`; both are pure Rust.
//!
//! Drawing goes through content [`Stream`]s, each like a PDFBox PDPageContentStream opened in append mode with a
//! reset context, with PDFBox's operators and number format (so a viewer computes the same positions as from
//! the Java tools' files). The output is deterministic: no dates, no random ids.

pub mod content;
pub mod font;
pub mod image;
pub mod table;

#[cfg(test)]
mod tests;

use std::collections::{BTreeSet, HashMap};

use pdf_writer::types::{CidFontType, FontFlags, SystemInfo, UnicodeCmap};
use pdf_writer::{Filter, Finish, Name, Pdf, Rect, Ref, Str, TextStr};

pub use content::Stream;
pub use font::{BoxMetrics, Font, MissingGlyph};

/// An sRGB colour (java.awt.Color without alpha)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Rgb(pub u8, pub u8, pub u8);

impl Rgb {
    pub const BLACK: Rgb = Rgb(0, 0, 0);
    pub const WHITE: Rgb = Rgb(255, 255, 255);
}

/// Page sizes as PDFBox computes them (PDRectangle.A4: 210 x 297 mm in float points)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PageSize {
    pub width: f32,
    pub height: f32,
}

impl PageSize {
    pub fn a4() -> PageSize {
        // POINTS_PER_MM = 1 / (10 * 2.54f) * POINTS_PER_INCH
        let per_mm = 1f32 / (10f32 * 2.54f32) * 72f32;
        PageSize { width: 210f32 * per_mm, height: 297f32 * per_mm }
    }
    pub fn a4_landscape() -> PageSize {
        let a4 = PageSize::a4();
        PageSize { width: a4.height, height: a4.width }
    }
}

pub type FontId = usize;
pub type ImageId = usize;

struct Page {
    size: PageSize,
    /// the page's content streams, each wrapped in q ... Q (PDFBox's reset context)
    content: Vec<u8>,
    fonts: BTreeSet<FontId>,
    images: BTreeSet<ImageId>,
}

/// A PDF document being built
pub struct Document {
    fonts: Vec<Font>,
    images: Vec<image::Image>,
    image_keys: HashMap<String, ImageId>,
    pages: Vec<Page>,
    outline: Vec<(String, usize)>,
    title: Option<String>,
    author: Option<String>,
    version: (u8, u8),
}

impl Default for Document {
    fn default() -> Self {
        Document::new()
    }
}

impl Document {
    pub fn new() -> Document {
        Document {
            fonts: Vec::new(),
            images: Vec::new(),
            image_keys: HashMap::new(),
            pages: Vec::new(),
            outline: Vec::new(),
            title: None,
            author: None,
            version: (1, 7),
        }
    }

    /// Adds a TrueType font (embedded as a subset of the glyphs used)
    pub fn add_font(&mut self, data: Vec<u8>) -> Result<FontId, String> {
        self.fonts.push(Font::from_bytes(data)?);
        Ok(self.fonts.len() - 1)
    }

    pub fn font(&self, id: FontId) -> &Font {
        &self.fonts[id]
    }

    /// Adds a PNG image; images added again under the same key are embedded once
    pub fn add_png(&mut self, key: &str, bytes: &[u8]) -> Result<ImageId, String> {
        if let Some(id) = self.image_keys.get(key) {
            return Ok(*id);
        }
        let img = image::Image::from_png(bytes).map_err(|e| format!("{key}: {e}"))?;
        self.images.push(img);
        let id = self.images.len() - 1;
        self.image_keys.insert(key.to_string(), id);
        Ok(id)
    }

    /// Adds the PNG file at `path` (once per path)
    pub fn add_png_file(&mut self, path: &str) -> Result<ImageId, String> {
        if let Some(id) = self.image_keys.get(path) {
            return Ok(*id);
        }
        let bytes = std::fs::read(path).map_err(|e| format!("Could not read image from {path}: {e}"))?;
        self.add_png(path, &bytes)
    }

    /// The pixel size of an image
    pub fn image_size(&self, id: ImageId) -> (u32, u32) {
        (self.images[id].width, self.images[id].height)
    }

    /// Appends a page and returns its index
    pub fn add_page(&mut self, size: PageSize) -> usize {
        self.pages.push(Page { size, content: Vec::new(), fonts: BTreeSet::new(), images: BTreeSet::new() });
        self.pages.len() - 1
    }

    pub fn page_count(&self) -> usize {
        self.pages.len()
    }

    pub fn page_size(&self, page: usize) -> PageSize {
        self.pages[page].size
    }

    /// A new content stream on page `page` (PDPageContentStream in append mode with a reset context)
    pub fn open_stream(&self, page: usize) -> Stream {
        Stream::new(page)
    }

    /// Appends a content stream to its page
    pub fn close_stream(&mut self, stream: Stream) {
        let (page, buf, fonts, images) = stream.into_parts();
        let p = &mut self.pages[page];
        p.content.extend_from_slice(b"q\n");
        p.content.extend_from_slice(&buf);
        p.content.extend_from_slice(b"Q\n");
        p.fonts.extend(fonts);
        p.images.extend(images);
    }

    /// Adds an outline entry (bookmark) at the top level pointing to `page`
    pub fn add_bookmark(&mut self, title: &str, page: usize) {
        self.outline.push((title.to_string(), page));
    }

    pub fn set_title(&mut self, title: &str) {
        self.title = Some(title.to_string());
    }

    pub fn set_author(&mut self, author: &str) {
        self.author = Some(author.to_string());
    }

    /// The PDF version in the header (1.7 by default)
    pub fn set_version(&mut self, major: u8, minor: u8) {
        self.version = (major, minor);
    }

    /// The PDF file
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        let mut pdf = Pdf::new();
        pdf.set_version(self.version.0, self.version.1);
        let mut next = 1i32;
        let mut alloc = || {
            let r = Ref::new(next);
            next += 1;
            r
        };
        let catalog_id = alloc();
        let pages_id = alloc();
        let info_id = alloc();
        let outline_id = if self.outline.is_empty() { None } else { Some(alloc()) };
        let page_ids: Vec<Ref> = self.pages.iter().map(|_| alloc()).collect();
        let content_ids: Vec<Ref> = self.pages.iter().map(|_| alloc()).collect();
        let font_ids: Vec<Ref> = self.fonts.iter().map(|_| alloc()).collect();
        let image_ids: Vec<(Ref, Option<Ref>)> = self.images.iter().map(|i| (alloc(), i.smask.as_ref().map(|_| alloc()))).collect();
        let outline_item_ids: Vec<Ref> = self.outline.iter().map(|_| alloc()).collect();

        let mut catalog = pdf.catalog(catalog_id);
        catalog.pages(pages_id);
        if let Some(o) = outline_id {
            catalog.outlines(o);
        }
        catalog.finish();
        pdf.pages(pages_id).kids(page_ids.iter().copied()).count(page_ids.len() as i32);

        let mut info = pdf.document_info(info_id);
        if let Some(t) = &self.title {
            info.title(TextStr(t));
        }
        if let Some(a) = &self.author {
            info.author(TextStr(a));
        }
        info.finish();

        for (i, page) in self.pages.iter().enumerate() {
            let mut p = pdf.page(page_ids[i]);
            p.parent(pages_id).media_box(Rect::new(0.0, 0.0, page.size.width, page.size.height)).contents(content_ids[i]);
            let mut res = p.resources();
            if !page.fonts.is_empty() {
                let mut fonts = res.fonts();
                for f in &page.fonts {
                    fonts.pair(Name(format!("F{}", f + 1).as_bytes()), font_ids[*f]);
                }
                fonts.finish();
            }
            if !page.images.is_empty() {
                let mut xo = res.x_objects();
                for im in &page.images {
                    xo.pair(Name(format!("Im{}", im + 1).as_bytes()), image_ids[*im].0);
                }
                xo.finish();
            }
            res.finish();
            p.finish();
            let z = image::zlib(&page.content, 6);
            pdf.stream(content_ids[i], &z).filter(Filter::FlateDecode);
        }

        for (fi, font) in self.fonts.iter().enumerate() {
            if !font.is_used() {
                // never drawn: an empty Type0 placeholder is not needed, but the id is referenced by no page
                pdf.indirect(font_ids[fi]).primitive(pdf_writer::Null);
                continue;
            }
            write_font(&mut pdf, font_ids[fi], font, &mut alloc)?;
        }

        for (ii, img) in self.images.iter().enumerate() {
            let (id, smask_id) = image_ids[ii];
            let (data, colors, predicted) = match &img.pixels {
                image::Pixels::PngPredicted { zlib, colors } => (zlib, *colors, true),
                image::Pixels::Flate { zlib, colors } => (zlib, *colors, false),
            };
            let mut x = pdf.image_xobject(id, data);
            x.filter(Filter::FlateDecode);
            x.width(img.width as i32).height(img.height as i32).bits_per_component(8);
            if colors == 1 {
                x.color_space().device_gray();
            } else {
                x.color_space().device_rgb();
            }
            if predicted {
                x.decode_parms().predictor(pdf_writer::types::Predictor::PngOptimum).colors(colors as i32).bits_per_component(8).columns(img.width as i32);
            }
            if let Some(s) = smask_id {
                x.s_mask(s);
            }
            x.finish();
            if let (Some(s), Some(mask)) = (smask_id, &img.smask) {
                let mut m = pdf.image_xobject(s, mask);
                m.filter(Filter::FlateDecode);
                m.width(img.width as i32).height(img.height as i32).bits_per_component(8);
                m.color_space().device_gray();
                m.finish();
            }
        }

        if let Some(oid) = outline_id {
            let n = outline_item_ids.len();
            let mut o = pdf.outline(oid);
            o.first(outline_item_ids[0]).last(outline_item_ids[n - 1]).count(n as i32);
            o.finish();
            for (k, (title, page)) in self.outline.iter().enumerate() {
                let mut item = pdf.outline_item(outline_item_ids[k]);
                item.title(TextStr(title)).parent(oid);
                if k > 0 {
                    item.prev(outline_item_ids[k - 1]);
                }
                if k + 1 < n {
                    item.next(outline_item_ids[k + 1]);
                }
                // PDPageFitWidthDestination without a top: [page /FitH null]
                item.dest().page(page_ids[*page]).fit_horizontal(self.pages[*page].size.height);
                item.finish();
            }
        }

        Ok(pdf.finish())
    }

    /// Writes the PDF file
    pub fn save(&self, path: &str) -> Result<(), String> {
        let bytes = self.to_bytes()?;
        std::fs::write(path, bytes).map_err(|e| format!("cannot write {path}: {e}"))
    }
}

fn write_font(pdf: &mut Pdf, type0_id: Ref, font: &Font, alloc: &mut impl FnMut() -> Ref) -> Result<(), String> {
    let cid_id = alloc();
    let descriptor_id = alloc();
    let file_id = alloc();
    let cmap_id = alloc();
    let (subset, glyphs) = font.subset()?;
    let base_name = format!("{}+{}", font.subset_tag(), font.name());
    let system_info = SystemInfo { registry: Str(b"Adobe"), ordering: Str(b"Identity"), supplement: 0 };

    pdf.type0_font(type0_id)
        .base_font(Name(base_name.as_bytes()))
        .encoding_predefined(Name(b"Identity-H"))
        .descendant_font(cid_id)
        .to_unicode(cmap_id);

    let mut cid = pdf.cid_font(cid_id);
    cid.subtype(CidFontType::Type2).base_font(Name(base_name.as_bytes())).system_info(system_info).font_descriptor(descriptor_id);
    cid.default_width(1000.0);
    {
        let mut w = cid.widths();
        for (gid, width, _) in &glyphs {
            w.consecutive(*gid, [*width]);
        }
        w.finish();
    }
    cid.cid_to_gid_map_predefined(Name(b"Identity"));
    cid.finish();

    let b = font.bbox();
    pdf.font_descriptor(descriptor_id)
        .name(Name(base_name.as_bytes()))
        .flags(FontFlags::SYMBOLIC)
        .bbox(Rect::new(b[0], b[1], b[2], b[3]))
        .italic_angle(font.italic_angle())
        .ascent(font.ascent())
        .descent(font.descent())
        .cap_height(font.cap_height())
        .x_height(font.x_height())
        .stem_v(80.0)
        .font_file2(file_id);

    let z = image::zlib(&subset, 6);
    pdf.stream(file_id, &z).filter(Filter::FlateDecode).pair(Name(b"Length1"), subset.len() as i32);

    let mut cmap = UnicodeCmap::new(Name(b"Custom"), system_info);
    for (gid, _, ch) in &glyphs {
        cmap.pair(*gid, *ch);
    }
    pdf.cmap(cmap_id, &cmap.finish());
    Ok(())
}
