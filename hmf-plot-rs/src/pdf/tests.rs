//! Tests of the PDF backend on synthetic content (the crate's DejaVu Sans): font metrics as PDFBox computes them,
//! boxable's line wrapping and table layout, PNG embedding, and the file structure.

use super::table::{Paragraph, Table, VAlign};
use super::*;

const DEJAVU: &[u8] = include_bytes!("../../fonts/DejaVuSans.ttf");

fn doc_with_font() -> (Document, FontId) {
    let mut doc = Document::new();
    let f = doc.add_font(DEJAVU.to_vec()).unwrap();
    (doc, f)
}

#[test]
fn font_metrics_as_pdfbox() {
    let (doc, f) = doc_with_font();
    let font = doc.font(f);
    assert_eq!(font.name(), "DejaVuSans");
    assert_eq!(font.units_per_em(), 2048);
    // widths: Math.round(advance * 1000 / 2048) in float (space 651 units, 'A' 1401 units)
    assert_eq!(font.char_width(' ').unwrap(), 318.0);
    assert_eq!(font.char_width('A').unwrap(), 684.0);
    // getStringWidth: the float sum of the character widths
    let w: f32 = "AA A".chars().map(|c| font.char_width(c).unwrap()).sum();
    assert_eq!(font.string_width("AA A").unwrap(), w);
    // ascent and descent from 'hhea' scaled in float
    assert_eq!(font.ascent(), 1901f32 * (1000f32 / 2048f32));
    assert_eq!(font.descent(), -483f32 * (1000f32 / 2048f32));
    // a character without a glyph: PDFBox's message
    let e = font.string_width("a\u{E000}").unwrap_err();
    assert_eq!(e.to_string(), "No glyph for U+E000 in font DejaVuSans");
}

fn line_texts(p: &Paragraph) -> Vec<String> {
    p.lines.iter().map(|l| l.tokens.concat()).collect()
}

#[test]
fn paragraphs_wrap_as_boxable() {
    let (doc, f) = doc_with_font();
    let size = 10f32;
    let w = |s: &str| doc.font(f).string_width(s).unwrap() / 1000f32 * size;
    // one line when it fits; tokens split after wrap characters, spaces kept at the token's end
    let p = Paragraph::new(&doc, f, size, 500.0, "c.35G>T [p.Gly12Val]").unwrap();
    assert_eq!(p.lines.len(), 1);
    assert_eq!(p.lines[0].tokens, vec!["c.", "35G>T ", "[p.", "Gly12Val]"]);
    // (the width adds up the tokens in float as PipelineLayer does)
    assert!((p.lines[0].width - w("c.35G>T [p.Gly12Val]")).abs() < 1e-3);
    // wrapped at the last wrap point that fits; the line's width leaves out trailing spaces
    let width = w("alpha beta") + 0.01;
    let p = Paragraph::new(&doc, f, size, width, "alpha beta alpha beta-alpha").unwrap();
    // (a hyphen is a wrap point; "alpha beta-" is wider than the cell)
    assert_eq!(line_texts(&p), vec!["alpha beta ", "alpha ", "beta-", "alpha"]);
    assert!((p.lines[0].width - w("alpha beta")).abs() < 1e-3);
    // a word wider than the cell is broken where it reaches the width
    let p = Paragraph::new(&doc, f, size, w("ABCD") + 0.01, "ABCDEFGHIJ").unwrap();
    assert_eq!(line_texts(&p), vec!["ABCD", "EFGH", "IJ"]);
    // empty text has no lines (height 0)
    let p = Paragraph::new(&doc, f, size, 100.0, "").unwrap();
    assert!(p.lines.is_empty());
    assert_eq!(p.height_from(0, 1.0), 0.0);
    // line heights are the font's boxable height
    let p = Paragraph::new(&doc, f, size, 100.0, "x").unwrap();
    assert_eq!(p.lines[0].height, BoxMetrics::of(doc.font(f)).height(size));
}

fn table_with_rows(doc: &mut Document, f: FontId, n: usize, header: bool) -> Table {
    let page = doc.add_page(PageSize::a4());
    let mut t = Table::new(doc, 700.0, 741.89, 40.0, 500.0, 30.0, page);
    if header {
        let r = t.create_row(16.0);
        for name in ["A", "B"] {
            let c = r.create_cell(50.0, name, f).unwrap();
            c.valign = VAlign::Bottom;
        }
        let i = t.last_row_index();
        t.add_header_row(i);
    }
    for k in 0..n {
        let r = t.create_row(15.0);
        r.create_cell(50.0, &format!("row {k}"), f).unwrap();
        r.create_cell(50.0, "x", f).unwrap();
    }
    t
}

#[test]
fn tables_break_pages_and_repeat_headers() {
    let (mut doc, f) = doc_with_font();
    // boxable's default cells: padding 5, borders 1, font size 8: a row is max(15, text + 12) high
    let h = 8f32 * BoxMetrics::of(doc.font(f)).height + 12f32;
    let mut t = table_with_rows(&mut doc, f, 100, true);
    let header_h = 16f32.max(h);
    let row_h = 15f32.max(h);
    let y = t.draw(&mut doc).unwrap();
    // rows per page: from 700 on the first page and from 741.89 on the next, down to 40, each page starting with
    // the header
    let mut pages = 1;
    let mut cursor = 700f32 - header_h;
    for _ in 0..100 {
        if cursor - row_h < 40.0 {
            pages += 1;
            cursor = 741.89 - header_h;
        }
        cursor -= row_h;
    }
    assert_eq!(doc.page_count(), pages);
    assert!((y - cursor).abs() < 0.01, "{y} {cursor}");
    assert_eq!(t.current_page(), pages - 1);
}

#[test]
fn keep_together_moves_the_first_rows() {
    let (mut doc, f) = doc_with_font();
    let page = doc.add_page(PageSize::a4());
    // two rows of 15 at y 70: only one fits above 40
    let mut t = Table::new(&doc, 70.0, 741.89, 40.0, 500.0, 30.0, page);
    for k in 0..2 {
        t.create_row(15.0).create_cell(100.0, &format!("r{k}"), f).unwrap();
    }
    let mut t2 = t.clone();
    t.draw(&mut doc).unwrap();
    assert_eq!(t.current_page(), 1, "boxable: the second row moves to a new page");
    t2.set_keep_together(2);
    let pages = doc.page_count();
    t2.draw(&mut doc).unwrap();
    assert_eq!(doc.page_count(), pages + 1);
    assert_eq!(t2.current_page(), pages, "both rows on the new page");
}

fn png(color: png::ColorType, w: u32, h: u32, data: &[u8], trns: Option<&[u8]>) -> Vec<u8> {
    let mut out = Vec::new();
    {
        let mut e = png::Encoder::new(&mut out, w, h);
        e.set_color(color);
        e.set_depth(png::BitDepth::Eight);
        if let Some(t) = trns {
            e.set_trns(t.to_vec());
        }
        let mut wr = e.write_header().unwrap();
        wr.write_image_data(data).unwrap();
    }
    out
}

#[test]
fn png_images_as_pdfbox_embeds_them() {
    // RGB: the IDAT data as it is, no soft mask
    let rgb = png(png::ColorType::Rgb, 2, 1, &[255, 0, 0, 0, 0, 255], None);
    let img = image::Image::from_png(&rgb).unwrap();
    assert!(matches!(img.pixels, image::Pixels::PngPredicted { colors: 3, .. }));
    assert!(img.smask.is_none());
    // RGBA: colour and alpha separated, a soft mask even when opaque (as LosslessFactory)
    let rgba = png(png::ColorType::Rgba, 2, 1, &[1, 2, 3, 255, 4, 5, 6, 0], None);
    let img = image::Image::from_png(&rgba).unwrap();
    let image::Pixels::Flate { zlib, colors } = &img.pixels else { panic!() };
    assert_eq!(*colors, 3);
    assert_eq!(miniz_oxide::inflate::decompress_to_vec_zlib(zlib).unwrap(), vec![1, 2, 3, 4, 5, 6]);
    assert_eq!(miniz_oxide::inflate::decompress_to_vec_zlib(img.smask.as_ref().unwrap()).unwrap(), vec![255, 0]);
    // RGB with a tRNS colour key: the colour data as it is, the key colour transparent
    let keyed = png(png::ColorType::Rgb, 2, 1, &[255, 255, 255, 9, 9, 9], Some(&[0, 255, 0, 255, 0, 255]));
    let img = image::Image::from_png(&keyed).unwrap();
    assert!(matches!(img.pixels, image::Pixels::PngPredicted { colors: 3, .. }));
    assert_eq!(miniz_oxide::inflate::decompress_to_vec_zlib(img.smask.as_ref().unwrap()).unwrap(), vec![0, 255]);
}

fn count(hay: &[u8], needle: &[u8]) -> usize {
    hay.windows(needle.len()).filter(|w| *w == needle).count()
}

#[test]
fn document_structure() {
    let (mut doc, f) = doc_with_font();
    let p0 = doc.add_page(PageSize::a4());
    let p1 = doc.add_page(PageSize::a4_landscape());
    let img = doc.add_png("x", &png(png::ColorType::Rgb, 1, 1, &[0, 0, 0], None)).unwrap();
    assert_eq!(doc.add_png("x", &[]).unwrap(), img, "the same key is embedded once");
    for (page, text) in [(p0, "Hello"), (p1, "World")] {
        let mut s = doc.open_stream(page);
        s.begin_text();
        s.set_font(f, 10.0);
        s.new_line_at_offset(30.0, 700.0);
        s.show_text(&mut doc, f, text).unwrap();
        s.end_text();
        s.draw_image(img, 10.0, 10.0, 20.0, 20.0);
        doc.close_stream(s);
    }
    doc.add_bookmark("First", p0);
    doc.set_title("T");
    let bytes = doc.to_bytes().unwrap();
    assert!(bytes.starts_with(b"%PDF-"));
    assert_eq!(count(&bytes, b"/Type /Page\n"), 2);
    assert_eq!(count(&bytes, b"/MediaBox [0 0 595.27"), 1);
    assert_eq!(count(&bytes, b"/MediaBox [0 0 841.8"), 1);
    assert_eq!(count(&bytes, b"/Subtype /Type0"), 1);
    assert_eq!(count(&bytes, b"/ToUnicode"), 1);
    assert_eq!(count(&bytes, b"/Subtype /Image"), 1);
    assert_eq!(count(&bytes, b"/Type /Outlines"), 1);
    // deterministic
    assert_eq!(bytes, doc.to_bytes().unwrap());
    // A4 as PDFBox computes it
    assert_eq!(PageSize::a4().width, 210f32 * (1f32 / (10f32 * 2.54f32) * 72f32));
}
