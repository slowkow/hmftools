//! ORANGE's PDF report (`orange/report`: ReportWriter, DocumentContext, PageEventHandler, the components, chapters
//! and tables), drawn with hmf-plot-rs's PDF backend from the same record and plots as Java. The page layout follows
//! Java's PDFBox and boxable code step by step (cursor, spacing, page breaks, table layout in float arithmetic), with
//! the same fonts (Nimbus Sans L, ORANGE's resources), so that pages break at the same rows and text sits at the
//! same places.

mod chapters;
mod common;
mod tables;

use anyhow::{anyhow, Result};
use hmf_plot::pdf::table::Table;
use hmf_plot::pdf::{Document, FontId, ImageId, MissingGlyph, PageSize, Rgb};

use crate::compat::java;
use crate::config::OrangeConfig;
use crate::model::OrangeRecord;

pub use common::{Resources, TextStyle};

const FONT_REGULAR: &[u8] = include_bytes!("../../resources/fonts/NimbusSansL-Regular.ttf");
const FONT_BOLD: &[u8] = include_bytes!("../../resources/fonts/NimbusSansL-Bold.ttf");
const ORANGE_CIRCOS: &[u8] = include_bytes!("../../resources/orange_circos.png");

pub const METADATA_TITLE: &str = "HMF ORANGE Report v5.1";
pub const METADATA_AUTHOR: &str = "Hartwig Pipeline";

pub const PAGE_MARGIN_TOP: f32 = 100.0;
pub const PAGE_MARGIN_LEFT: f32 = 30.0;
pub const PAGE_MARGIN_RIGHT: f32 = 30.0;
pub const PAGE_MARGIN_BOTTOM: f32 = 40.0;

const HEADER_ORANGE_HEIGHT: f32 = 60.0;

/// A Java exception message for a character without a glyph (PDFBox's IllegalArgumentException)
pub(crate) fn glyph_error(e: MissingGlyph) -> anyhow::Error {
    anyhow!("java.lang.IllegalArgumentException: {e}")
}

/// PlotPathResolver: a plot path as it is when it exists (or without an output directory), else relative to the
/// output directory
pub(crate) fn resolve_plot(output_dir: Option<&str>, path: &str) -> String {
    match output_dir {
        Some(dir) if !std::path::Path::new(path).exists() => format!("{dir}/{path}"),
        _ => path.to_string(),
    }
}

/// Header, side panel and footer (PageEventHandler with Header, SidePanel and Footer)
struct PageEvents {
    sample_id: String,
    add_disclaimer: bool,
    logo: ImageId,
    chapter_title: String,
    first_page_of_chapter: bool,
}

impl PageEvents {
    /// onPageStart: header, the chapter's bookmark on its first page, side panel (the footer only records the page)
    fn on_page_start(&mut self, doc: &mut Document, res: &Resources, page: usize) -> Result<()> {
        self.render_header(doc, res, page)?;
        if self.first_page_of_chapter {
            self.first_page_of_chapter = false;
            doc.add_bookmark(&self.chapter_title, page);
        }
        self.render_side_panel(doc, res, page)?;
        Ok(())
    }

    /// Header.renderHeader
    fn render_header(&self, doc: &mut Document, res: &Resources, page: usize) -> Result<()> {
        let page_height = doc.page_size(page).height;
        let (iw, ih) = doc.image_size(self.logo);
        let img_scale = HEADER_ORANGE_HEIGHT / ih as f32;
        let img_width = iw as f32 * img_scale;
        let y = page_height - HEADER_ORANGE_HEIGHT - 10.0;
        let mut s = doc.open_stream(page);
        s.draw_image(self.logo, 50.0, y, img_width, HEADER_ORANGE_HEIGHT);
        let left = if self.add_disclaimer { 150.0 } else { 180.0 };
        let text_y = page_height - 40.0;
        let colors = [common::PALETTE_ORANGE_1, common::PALETTE_ORANGE_2, common::PALETTE_ORANGE_3, common::PALETTE_ORANGE_4, common::PALETTE_ORANGE_5, common::PALETTE_ORANGE_6];
        let letters = ["O", "R", "A", "N", "G", "E"];
        s.begin_text();
        s.set_font(res.bold, 11.0);
        s.new_line_at_offset(left, text_y);
        for (l, c) in letters.iter().zip(colors.iter()) {
            s.set_non_stroking_color(*c);
            s.show_text(doc, res.bold, l).map_err(glyph_error)?;
        }
        s.set_non_stroking_color(Rgb::BLACK);
        s.show_text(doc, res.bold, " Report").map_err(glyph_error)?;
        if self.add_disclaimer {
            s.set_font(res.bold, 9.0);
            s.show_text(doc, res.bold, " (Research Use Only)").map_err(glyph_error)?;
        }
        s.end_text();
        doc.close_stream(s);
        Ok(())
    }

    /// SidePanel.renderSidePanel
    fn render_side_panel(&self, doc: &mut Document, res: &Resources, page: usize) -> Result<()> {
        const ROW_SPACING: f32 = 35.0;
        const VALUE_TEXT_Y_OFFSET: f32 = 18.0;
        const MAX_WIDTH: f32 = 120.0;
        const RECTANGLE_WIDTH: f32 = 170.0;
        const RECTANGLE_HEIGHT: f32 = 60.0;
        let size = doc.page_size(page);
        let mut s = doc.open_stream(page);
        s.set_non_stroking_color(common::PALETTE_ORANGE);
        s.add_rect(size.width - RECTANGLE_WIDTH, size.height - RECTANGLE_HEIGHT, RECTANGLE_WIDTH, RECTANGLE_HEIGHT);
        s.fill();
        let x = size.width - RECTANGLE_WIDTH + 15.0;
        let mut y = (size.height + 15.0) - ROW_SPACING;
        let label = res.side_panel_label_style();
        s.begin_text();
        s.set_font(label.font, label.size);
        s.set_non_stroking_color(label.color);
        s.new_line_at_offset(x, y);
        s.show_text(doc, label.font, "SAMPLE").map_err(glyph_error)?;
        s.end_text();
        let value_size = max_point_size_for_width(doc, res.bold, 11.0, 6.0, &self.sample_id, MAX_WIDTH);
        y -= VALUE_TEXT_Y_OFFSET;
        s.begin_text();
        s.set_font(res.bold, value_size);
        s.set_non_stroking_color(Rgb::WHITE);
        s.new_line_at_offset(x, y);
        s.show_text(doc, res.bold, &self.sample_id).map_err(glyph_error)?;
        s.end_text();
        doc.close_stream(s);
        Ok(())
    }

    /// Footer.writeFooters: page numbers (and the disclaimer) on every page
    fn write_footers(&self, doc: &mut Document, res: &Resources) -> Result<()> {
        let total = doc.page_count();
        for i in 0..total {
            let mut s = doc.open_stream(i);
            let st = res.page_number_style();
            s.begin_text();
            s.set_font(st.font, st.size);
            s.set_non_stroking_color(st.color);
            s.new_line_at_offset(58.0, 20.0);
            s.show_text(doc, st.font, &format!("{}/{}", i + 1, total)).map_err(glyph_error)?;
            s.end_text();
            if self.add_disclaimer {
                let d = res.disclaimer_style();
                let text = "This report is for 'Research Use Only (RUO)' and is not suitable for diagnostic or clinical applications. No rights can be derived from the contents of this report.";
                s.begin_text();
                s.set_font(d.font, d.size);
                s.set_non_stroking_color(d.color);
                s.new_line_at_offset(98.0, 20.0);
                s.show_text(doc, d.font, text).map_err(glyph_error)?;
                s.end_text();
            }
            doc.close_stream(s);
        }
        Ok(())
    }
}

/// SidePanel.maxPointSizeForWidth: the font size shrunk by 0.1 until the text fits
fn max_point_size_for_width(doc: &Document, font: FontId, initial: f32, min: f32, text: &str, max_width: f32) -> f32 {
    let increment = 0.1f32;
    let mut size = initial;
    let Ok(w) = doc.font(font).string_width(text) else { return size };
    let mut width = w / 1000f32 * size;
    while width > max_width && size > min {
        size -= increment;
        width = w / 1000f32 * size;
    }
    size
}

/// DocumentContext: the document, its current page and the vertical cursor
pub struct DocCtx {
    pub doc: Document,
    pub res: Resources,
    current_page: Option<usize>,
    page_size: PageSize,
    cursor_y: f32,
    events: PageEvents,
    output_dir: Option<String>,
}

impl DocCtx {
    pub fn page_size(&self) -> PageSize {
        self.page_size
    }
    pub fn cursor_y(&self) -> f32 {
        self.cursor_y
    }
    pub fn set_cursor_y(&mut self, y: f32) {
        self.cursor_y = y;
    }
    pub fn content_width(&self) -> f32 {
        self.page_size.width - PAGE_MARGIN_LEFT - PAGE_MARGIN_RIGHT
    }
    pub fn content_start_y(&self) -> f32 {
        self.page_size.height - PAGE_MARGIN_TOP
    }
    pub fn content_end_y(&self) -> f32 {
        PAGE_MARGIN_BOTTOM
    }
    pub fn margin_left(&self) -> f32 {
        PAGE_MARGIN_LEFT
    }
    pub fn current_page(&self) -> usize {
        self.current_page.expect("a page")
    }
    pub fn resolve_plot(&self, path: &str) -> String {
        resolve_plot(self.output_dir.as_deref(), path)
    }

    /// newPage(pageSize)
    pub fn new_page_sized(&mut self, size: PageSize) -> Result<()> {
        self.page_size = size;
        let page = self.doc.add_page(size);
        self.current_page = Some(page);
        self.cursor_y = self.content_start_y();
        self.events.on_page_start(&mut self.doc, &self.res, page)
    }

    pub fn new_page(&mut self) -> Result<()> {
        self.new_page_sized(self.page_size)
    }

    fn ensure_page(&mut self) -> Result<()> {
        if self.current_page.is_none() {
            self.new_page()?;
        }
        Ok(())
    }

    fn has_space_for(&self, height: f32) -> bool {
        self.current_page.is_some() && (self.cursor_y - height) >= self.content_end_y()
    }

    /// createTable(width, columnWidths): a table at the left margin from the cursor
    pub fn create_table(&mut self, width: f32) -> Result<Table> {
        self.create_table_at_x(width, PAGE_MARGIN_LEFT)
    }

    /// createTableAtX
    pub fn create_table_at_x(&mut self, width: f32, x: f32) -> Result<Table> {
        self.ensure_page()?;
        Ok(Table::new(&self.doc, self.cursor_y, self.content_start_y(), self.content_end_y(), width, x, self.current_page()))
    }

    /// addTable: draws the table; when it ended on another page, that page gets the header and side panel. OR16:
    /// Java decorates only the last page, so the pages in between of a table over three or more pages have no
    /// header, side panel (sample id); corrected: every page the table adds
    pub fn add_table(&mut self, mut table: Table) -> Result<()> {
        self.ensure_page()?;
        let pages_before = self.doc.page_count();
        let final_y = table.draw(&mut self.doc).map_err(glyph_error)?;
        if java("OR16") {
            if table.current_page() != self.current_page() {
                self.current_page = Some(table.current_page());
                let page = self.current_page();
                self.events.on_page_start(&mut self.doc, &self.res, page)?;
            }
        } else {
            self.decorate_new_pages(pages_before)?;
            self.current_page = Some(table.current_page());
        }
        self.cursor_y = final_y - 5.0;
        Ok(())
    }

    /// The page events of the pages added since `pages_before` (OR16 corrected)
    fn decorate_new_pages(&mut self, pages_before: usize) -> Result<()> {
        for page in pages_before..self.doc.page_count() {
            self.events.on_page_start(&mut self.doc, &self.res, page)?;
        }
        Ok(())
    }

    /// addTableNoAdvance: draws the table and returns the y below it, the cursor unchanged (OR16: Java decorates
    /// none of the pages the table adds)
    pub fn add_table_no_advance(&mut self, mut table: Table) -> Result<f32> {
        self.ensure_page()?;
        let pages_before = self.doc.page_count();
        let final_y = table.draw(&mut self.doc).map_err(glyph_error)?;
        if !java("OR16") {
            self.decorate_new_pages(pages_before)?;
        }
        self.current_page = Some(table.current_page());
        Ok(final_y)
    }

    /// addParagraph: one line of text at the margin, on a new page when the line does not fit
    pub fn add_paragraph(&mut self, text: &str, style: TextStyle) -> Result<()> {
        self.ensure_page()?;
        let line_height = style.size * 1.4f32;
        if !self.has_space_for(line_height) {
            self.new_page()?;
        }
        let (page, y) = (self.current_page(), self.cursor_y - style.size);
        let mut s = self.doc.open_stream(page);
        s.begin_text();
        s.set_font(style.font, style.size);
        s.set_non_stroking_color(style.color);
        s.new_line_at_offset(PAGE_MARGIN_LEFT, y);
        s.show_text(&mut self.doc, style.font, text).map_err(glyph_error)?;
        s.end_text();
        self.doc.close_stream(s);
        self.cursor_y -= line_height;
        Ok(())
    }

    fn load_image(&mut self, path: &str) -> Result<ImageId> {
        self.doc.add_png_file(path).map_err(|e| anyhow!("java.io.IOException: {e}"))
    }

    /// The scale of an image to fit (max_w, max_h), never above 1
    fn fit(&self, image: ImageId, max_w: f32, max_h: f32) -> (f32, f32) {
        let (w, h) = self.doc.image_size(image);
        let (w, h) = (w as f32, h as f32);
        let mut scale = (max_w / w).min(max_h / h);
        if scale > 1.0f32 {
            scale = 1.0;
        }
        (w * scale, h * scale)
    }

    /// addImage(path, maxWidth, maxHeight): centred, on a new page when it does not fit
    pub fn add_image(&mut self, path: &str, max_w: f32, max_h: f32) -> Result<()> {
        self.ensure_page()?;
        let image = self.load_image(path)?;
        let (dw, dh) = self.fit(image, max_w, max_h);
        if !self.has_space_for(dh) {
            self.new_page()?;
        }
        let x = PAGE_MARGIN_LEFT + (self.content_width() - dw) / 2f32;
        let y = self.cursor_y - dh;
        self.draw_image(image, x, y, dw, dh);
        self.cursor_y = y - 5.0;
        Ok(())
    }

    /// addImageAt: centred in a cell of max_w from (x, y) down, the cursor unchanged
    pub fn add_image_at(&mut self, path: &str, x: f32, y: f32, max_w: f32, max_h: f32) -> Result<f32> {
        self.ensure_page()?;
        let image = self.load_image(path)?;
        let (dw, dh) = self.fit(image, max_w, max_h);
        self.draw_image(image, x + (max_w - dw) / 2f32, y - dh, dw, dh);
        Ok(dh)
    }

    /// addTitledImage: a title line and the image below it
    pub fn add_titled_image(&mut self, text: &str, style: TextStyle, path: &str, x: f32, y: f32, max_w: f32, max_h: f32) -> Result<f32> {
        self.ensure_page()?;
        let line_height = style.size * 1.4f32;
        let image = self.load_image(path)?;
        let (dw, dh) = self.fit(image, max_w, max_h);
        let xpos = x + (max_w - dw) / 2f32;
        let ypos = y - dh - line_height * 2f32;
        let mut s = self.doc.open_stream(self.current_page());
        s.begin_text();
        s.set_font(style.font, style.size);
        s.new_line_at_offset(x + 10.0, y - line_height);
        s.set_non_stroking_color(style.color);
        s.show_text(&mut self.doc, style.font, text).map_err(glyph_error)?;
        s.end_text();
        s.draw_image(image, xpos, ypos, dw, dh);
        self.doc.close_stream(s);
        Ok(dh)
    }

    /// addSpacing
    pub fn add_spacing(&mut self, height: f32) -> Result<()> {
        self.ensure_page()?;
        if !self.has_space_for(height) {
            self.new_page()?;
        } else {
            self.cursor_y -= height;
        }
        Ok(())
    }

    /// An image in its own content stream (drawImage)
    fn draw_image(&mut self, image: ImageId, x: f32, y: f32, w: f32, h: f32) {
        let mut s = self.doc.open_stream(self.current_page());
        s.draw_image(image, x, y, w, h);
        self.doc.close_stream(s);
    }

    /// A stroked rectangle on the current page (FrontPageChapter.drawBoxBorder)
    pub fn draw_box_border(&mut self, x: f32, y: f32, w: f32, h: f32) {
        let mut s = self.doc.open_stream(self.current_page());
        s.set_stroking_color(Rgb::BLACK);
        s.set_line_width(0.5);
        s.add_rect(x, y, w, h);
        s.stroke();
        self.doc.close_stream(s);
    }
}

impl DocCtx {
    /// A document with ORANGE's fonts, header image and metadata (ReportResources.create, PageEventHandler.create,
    /// the DocumentContext of ReportWriter.writePdfChapters)
    pub fn new(sample_id: &str, add_disclaimer: bool, output_dir: Option<&str>) -> Result<DocCtx> {
        let mut doc = Document::new();
        let regular = doc.add_font(FONT_REGULAR.to_vec()).map_err(|e| anyhow!(e))?;
        let bold = doc.add_font(FONT_BOLD.to_vec()).map_err(|e| anyhow!(e))?;
        let logo = doc.add_png("orange_circos.png", ORANGE_CIRCOS).map_err(|e| anyhow!("Failed to load header image: {e}"))?;
        doc.set_title(METADATA_TITLE);
        doc.set_author(METADATA_AUTHOR);
        // PDFBox 3's default
        doc.set_version(1, 6);
        let events = PageEvents {
            sample_id: sample_id.to_string(),
            add_disclaimer,
            logo,
            chapter_title: "Undefined".into(),
            first_page_of_chapter: true,
        };
        Ok(DocCtx {
            doc,
            res: Resources { regular, bold },
            current_page: None,
            page_size: PageSize::a4(),
            cursor_y: 0.0,
            events,
            output_dir: output_dir.map(|d| d.to_string()),
        })
    }

    /// The start of a chapter: its page size, title for the bookmark, and a new page
    pub fn start_chapter(&mut self, name: &str, size: PageSize) -> Result<()> {
        self.page_size = size;
        self.events.chapter_title = name.to_string();
        self.events.first_page_of_chapter = true;
        self.new_page_sized(size)
    }

    /// The footers on every page (second pass), and the file
    pub fn finish(self) -> Result<Document> {
        let DocCtx { mut doc, res, events, .. } = self;
        events.write_footers(&mut doc, &res)?;
        Ok(doc)
    }
}

/// ReportWriter.writePdf: the chapters in order and the footers
pub fn render_pdf(record: &OrangeRecord, config: &OrangeConfig) -> Result<Document> {
    let mut ctx = DocCtx::new(&config.display_sample_id, config.add_disclaimer, Some(&config.output_dir))?;
    for chapter in chapters::chapters(record) {
        ctx.start_chapter(chapter.name(), chapter.page_size())?;
        chapter.render(&mut ctx, record, config)?;
    }
    ctx.finish()
}

/// The PDF report written to `path` (`<output_dir><tumor>.orange[.<id>].pdf`)
pub fn write_pdf(record: &OrangeRecord, config: &OrangeConfig, path: &str) -> Result<()> {
    let t0 = std::time::Instant::now();
    let doc = render_pdf(record, config)?;
    hmf_common::debug!("PDF: {} pages laid out and drawn (images read) in {:.3} s", doc.page_count(), t0.elapsed().as_secs_f64());
    hmf_common::info!("writing PDF report to {path}");
    let bytes = doc.to_bytes().map_err(|e| anyhow!(e))?;
    hmf_common::debug!("PDF: file of {} bytes made at {:.3} s", bytes.len(), t0.elapsed().as_secs_f64());
    std::fs::write(path, bytes).map_err(|e| anyhow!("cannot write {path}: {e}"))?;
    hmf_common::debug!("PDF: written at {:.3} s", t0.elapsed().as_secs_f64());
    Ok(())
}

#[cfg(test)]
mod tests;
