//! Tables laid out and drawn as boxable 1.8.2 does (be.quodlibet.boxable: BaseTable, Row, Cell, Paragraph,
//! Tokenizer, PipelineLayer, FontUtils), for reports whose Java version draws its tables with boxable: the same
//! cell widths, line wrapping, row heights, page breaks (with header rows repeated), fills, borders and text
//! positions, computed in the same float arithmetic.
//!
//! What is ported: left-aligned rows of text cells, horizontal and vertical alignment within a cell, padding,
//! borders, fill colours, header rows, rows split across pages, words wider than their cell. Not ported (no
//! report uses them): boxable's HTML-like markup in cell text (`<b>`, `<i>`, `<p>`, `<br>`, lists, headings:
//! drawn as plain text here; boxable would switch to Helvetica or break lines), rotated text, image and table
//! cells, fixed row heights, outer borders, bookmarks and the "start a new page if the table does not fit"
//! option.

use super::{BoxMetrics, Document, FontId, MissingGlyph, PageSize, Rgb, Stream};

/// boxable's LineStyle (solid lines only)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LineStyle {
    pub color: Rgb,
    pub width: f32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HAlign {
    Left,
    Center,
    Right,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VAlign {
    Top,
    Middle,
    Bottom,
}

// ------------------------------------------------------------------------------------------------
// Tokenizer and Paragraph

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Text(String),
    PossibleWrapPoint,
}

/// Tokenizer.isWrapPointChar
fn is_wrap_point_char(c: char) -> bool {
    matches!(c, ' ' | ',' | '.' | '-' | '@' | ':' | ';' | '\n' | '\t' | '\r' | '\u{c}' | '\u{b}')
}

/// Tokenizer.tokenize without a wrapping function, for text without markup: a possible wrap point after every
/// wrap character (but the last character), text between them
fn tokenize(text: &str) -> Vec<Token> {
    let chars: Vec<char> = text.chars().collect();
    let n = chars.len();
    // findWrapPoints: a stack whose top is the first wrap point
    let mut points: Vec<usize> = vec![n];
    if n >= 2 {
        for i in (0..=n - 2).rev() {
            if is_wrap_point_char(chars[i]) {
                points.push(i + 1);
            }
        }
    }
    let mut tokens = Vec::new();
    let mut sb = String::new();
    let mut current = points.pop().unwrap_or(n);
    for (i, c) in chars.iter().enumerate() {
        if i == current {
            if !sb.is_empty() {
                tokens.push(Token::Text(std::mem::take(&mut sb)));
            }
            tokens.push(Token::PossibleWrapPoint);
            current = points.pop().unwrap_or(n);
        }
        sb.push(*c);
    }
    if !sb.is_empty() {
        tokens.push(Token::Text(sb));
    }
    tokens.push(Token::PossibleWrapPoint);
    tokens
}

/// Whether `text` holds something boxable would read as markup (`<b>`, `<i>`, `<p>`, `<br>`, `<u>`, `<s>`,
/// `<sup>`, `<sub>`, `<h1>`..`<h6>`, `<ol>`, `<ul>`, `<li>` or their closing tags), which this port draws as text
pub fn has_markup(text: &str) -> bool {
    const TAGS: [&str; 19] =
        ["b>", "i>", "p>", "u>", "s>", "br>", "br/>", "br />", "sup>", "sub>", "ol>", "ul>", "li>", "h1>", "h2>", "h3>", "h4>", "h5>", "h6>"];
    text.match_indices('<').any(|(i, _)| {
        let rest = &text[i + 1..];
        let rest = rest.strip_prefix('/').unwrap_or(rest);
        TAGS.iter().any(|t| rest.starts_with(t))
    })
}

fn rtrim(s: &str) -> &str {
    s.trim_end_matches(|c: char| c <= ' ')
}

/// PipelineLayer: text gathered since a wrap point, or on the current line
#[derive(Default)]
struct Layer {
    text: String,
    last: String,
    trimmed_last: String,
    width: f32,
    w_last: f32,
    w_trimmed_last: f32,
    w_current_text: f32,
    tokens: Vec<String>,
}

impl Layer {
    fn is_empty(&self) -> bool {
        self.tokens.is_empty()
    }
    fn push_text(&mut self, m: &Measure, word: &str) -> Result<(), MissingGlyph> {
        self.text.push_str(&self.last);
        self.width += self.w_last;
        self.last = word.to_string();
        self.trimmed_last = rtrim(word).to_string();
        self.w_last = m.width(word)?;
        self.w_trimmed_last = if self.trimmed_last.len() == self.last.len() { self.w_last } else { m.width(&self.trimmed_last)? };
        self.w_current_text = if self.text.is_empty() { 0.0 } else { m.width(&self.text)? };
        self.tokens.push(word.to_string());
        Ok(())
    }
    fn push_layer(&mut self, other: &mut Layer) {
        self.text.push_str(&self.last);
        self.width += self.w_last;
        self.text.push_str(&other.text);
        if !other.text.is_empty() {
            self.width += other.w_current_text;
        }
        self.last = std::mem::take(&mut other.last);
        self.trimmed_last = std::mem::take(&mut other.trimmed_last);
        self.w_last = other.w_last;
        self.w_trimmed_last = other.w_trimmed_last;
        self.tokens.append(&mut other.tokens);
        other.reset();
    }
    fn reset(&mut self) {
        self.text.clear();
        self.width = 0.0;
        self.last.clear();
        self.trimmed_last.clear();
        self.w_last = 0.0;
        self.w_trimmed_last = 0.0;
        self.tokens.clear();
    }
    fn width(&self) -> f32 {
        self.width + self.w_last
    }
    fn trimmed_width(&self) -> f32 {
        self.width + self.w_trimmed_last
    }
}

/// A font at a size: string widths in points as boxable computes them (getStringWidth / 1000f * fontSize)
struct Measure<'a> {
    font: &'a super::Font,
    size: f32,
}

impl Measure<'_> {
    fn width(&self, s: &str) -> Result<f32, MissingGlyph> {
        Ok(self.font.string_width(s)? / 1000f32 * self.size)
    }
}

/// One line of a paragraph: the text tokens drawn on it, its trimmed width and its height
#[derive(Clone, Debug, PartialEq)]
pub struct Line {
    pub tokens: Vec<String>,
    pub width: f32,
    pub height: f32,
}

/// Paragraph.getLines for plain text
#[derive(Clone, Debug, PartialEq)]
pub struct Paragraph {
    pub lines: Vec<Line>,
    pub max_line_width: f32,
    pub font_height: f32,
}

impl Paragraph {
    pub fn new(doc: &Document, font: FontId, size: f32, width: f32, text: &str) -> Result<Paragraph, MissingGlyph> {
        let f = doc.font(font);
        let m = Measure { font: f, size };
        let metrics = BoxMetrics::of(f);
        let default_height = metrics.height(size);
        let mut lines = Vec::new();
        let mut max_line_width = i32::MIN as f32;
        let mut in_line = Layer::default();
        let mut since = Layer::default();
        let mut in_line_h = 0f32;
        let mut since_h = 0f32;
        let emit = |layer: &mut Layer, h: &mut f32, lines: &mut Vec<Line>, max_w: &mut f32| {
            let w = layer.trimmed_width();
            lines.push(Line { tokens: layer.tokens.clone(), width: w, height: if *h == 0.0 { default_height } else { *h } });
            *max_w = max_w.max(w);
            layer.reset();
            *h = 0.0;
        };
        for token in tokenize(text) {
            match token {
                Token::Text(word) => {
                    let mut word = word;
                    let mut word_width = f.string_width(&word)?;
                    if word_width / 1000f32 * size > width && width > f.average_width() / 1000f32 * size {
                        // a word wider than the cell: broken where it reaches the width
                        let mut already_text_in_line = in_line.trimmed_width() > 0.0;
                        while word_width / 1000f32 * size > width {
                            let mut w = 0f32;
                            let mut first = String::new();
                            let mut rest = String::new();
                            let chars: Vec<char> = word.chars().collect();
                            for (i, c) in chars.iter().enumerate() {
                                w += f.char_width(*c)? / 1000f32 * size;
                                if already_text_in_line {
                                    if w < width - in_line.trimmed_width() {
                                        first.push(*c);
                                    } else {
                                        rest.push(*c);
                                    }
                                } else if w < width {
                                    first.push(*c);
                                } else if i == 0 {
                                    first.push(*c);
                                    rest.extend(chars[1..].iter());
                                    break;
                                } else {
                                    rest.push(*c);
                                }
                            }
                            already_text_in_line = false;
                            since.push_text(&m, &first)?;
                            since_h = since_h.max(default_height);
                            in_line.push_layer(&mut since);
                            in_line_h = in_line_h.max(since_h);
                            since_h = 0.0;
                            emit(&mut in_line, &mut in_line_h, &mut lines, &mut max_line_width);
                            word = rest;
                            word_width = f.string_width(&word)?;
                        }
                        since.push_text(&m, &word)?;
                        since_h = since_h.max(default_height);
                    } else {
                        since.push_text(&m, &word)?;
                        since_h = since_h.max(default_height);
                    }
                }
                Token::PossibleWrapPoint => {
                    if in_line.width() + since.trimmed_width() > width && !in_line.is_empty() {
                        emit(&mut in_line, &mut in_line_h, &mut lines, &mut max_line_width);
                    }
                    in_line.push_layer(&mut since);
                    in_line_h = in_line_h.max(since_h);
                    since_h = 0.0;
                }
            }
        }
        if since.trimmed_width() + in_line.trimmed_width() > 0.0 {
            in_line.push_layer(&mut since);
            in_line_h = in_line_h.max(since_h);
            emit(&mut in_line, &mut in_line_h, &mut lines, &mut max_line_width);
        }
        Ok(Paragraph { lines, max_line_width, font_height: default_height })
    }

    /// Paragraph.getHeight(lineStart) (line spacing 1: every line's height)
    pub fn height_from(&self, line_start: usize, line_spacing: f32) -> f32 {
        if self.lines.len() <= line_start {
            return 0.0;
        }
        let mut h = 0f32;
        for (i, l) in self.lines.iter().enumerate().skip(line_start) {
            h += if i == line_start { l.height } else { l.height * line_spacing };
        }
        h
    }
}

// ------------------------------------------------------------------------------------------------
// Cells, rows, tables

/// boxable's Cell (text cells)
#[derive(Clone, Debug)]
pub struct Cell {
    pub width: f32,
    pub text: String,
    pub font: FontId,
    pub font_size: f32,
    pub fill: Option<Rgb>,
    pub text_color: Rgb,
    pub left_padding: f32,
    pub right_padding: f32,
    pub top_padding: f32,
    pub bottom_padding: f32,
    pub left_border: Option<LineStyle>,
    pub right_border: Option<LineStyle>,
    pub top_border: Option<LineStyle>,
    pub bottom_border: Option<LineStyle>,
    pub align: HAlign,
    pub valign: VAlign,
    pub line_spacing: f32,
    line_start: usize,
    paragraph: Option<Paragraph>,
}

fn bw(b: &Option<LineStyle>) -> f32 {
    b.map(|l| l.width).unwrap_or(0.0)
}

impl Cell {
    /// Cell(row, width%, text, isCalculated = true): boxable's defaults (padding 5, black borders of width 1,
    /// font size 8, black text, top left alignment)
    fn new(row_width: f32, width_pct: f32, text: &str, font: FontId) -> Result<Cell, String> {
        let width = row_width * (width_pct / 100f32);
        if width > row_width {
            return Err(format!("Cell Width={} can't be bigger than row width={}", width, row_width));
        }
        let black = Some(LineStyle { color: Rgb::BLACK, width: 1.0 });
        Ok(Cell {
            width,
            text: text.to_string(),
            font,
            font_size: 8.0,
            fill: None,
            text_color: Rgb::BLACK,
            left_padding: 5.0,
            right_padding: 5.0,
            top_padding: 5.0,
            bottom_padding: 5.0,
            left_border: black,
            right_border: black,
            top_border: black,
            bottom_border: black,
            align: HAlign::Left,
            valign: VAlign::Top,
            line_spacing: 1.0,
            line_start: 0,
            paragraph: None,
        })
    }

    pub fn set_padding(&mut self, top: f32, bottom: f32, left: f32, right: f32) {
        self.top_padding = top;
        self.bottom_padding = bottom;
        self.left_padding = left;
        self.right_padding = right;
        self.paragraph = None;
    }

    pub fn set_borders(&mut self, top: Option<LineStyle>, bottom: Option<LineStyle>, left: Option<LineStyle>, right: Option<LineStyle>) {
        self.top_border = top;
        self.bottom_border = bottom;
        self.left_border = left;
        self.right_border = right;
        self.paragraph = None;
    }

    pub fn set_font(&mut self, font: FontId, size: f32, color: Rgb) {
        self.font = font;
        self.font_size = size;
        self.text_color = color;
        self.paragraph = None;
    }

    /// getInnerWidth
    pub fn inner_width(&self) -> f32 {
        self.width - self.left_padding - self.right_padding - bw(&self.left_border) - bw(&self.right_border)
    }

    fn paragraph(&mut self, doc: &Document) -> Result<&Paragraph, MissingGlyph> {
        if self.paragraph.is_none() {
            self.paragraph = Some(Paragraph::new(doc, self.font, self.font_size, self.inner_width(), &self.text)?);
        }
        Ok(self.paragraph.as_ref().unwrap())
    }

    /// getTextHeight
    fn text_height(&mut self, doc: &Document) -> Result<f32, MissingGlyph> {
        let (ls, sp) = (self.line_start, self.line_spacing);
        Ok(self.paragraph(doc)?.height_from(ls, sp))
    }

    /// getCellHeight
    fn cell_height(&mut self, doc: &Document) -> Result<f32, MissingGlyph> {
        Ok(self.text_height(doc)? + self.top_padding + self.bottom_padding + bw(&self.top_border) + bw(&self.bottom_border))
    }

    /// getInnerHeight, for the row's height
    fn inner_height(&self, row_height: f32) -> f32 {
        row_height - self.bottom_padding - self.top_padding - bw(&self.top_border) - bw(&self.bottom_border)
    }
}

/// boxable's Row (left-aligned, not of fixed height)
#[derive(Clone, Debug)]
pub struct Row {
    pub height: f32,
    pub cells: Vec<Cell>,
    width: f32,
}

impl Row {
    /// createCell(width%, text): the cell's left border is removed unless it is the row's first
    pub fn create_cell(&mut self, width_pct: f32, text: &str, font: FontId) -> Result<&mut Cell, String> {
        let mut cell = Cell::new(self.width, width_pct, text, font)?;
        if !self.cells.is_empty() {
            cell.left_border = None;
        }
        self.cells.push(cell);
        Ok(self.cells.last_mut().unwrap())
    }

    /// getHeight: the row's height or its tallest cell
    fn get_height(&mut self, doc: &Document) -> Result<f32, MissingGlyph> {
        let mut max = 0f32;
        for c in &mut self.cells {
            let h = c.cell_height(doc)?;
            if h > max {
                max = h;
            }
        }
        Ok(self.height.max(max))
    }

    fn total_cells_width(&self) -> f32 {
        self.cells.iter().fold(0f32, |s, c| s + c.width)
    }

    /// getLastCellLayoutWidth for a left-aligned row: the last cell stretched to the table's width
    fn last_cell_layout_width(&self, table_width: f32) -> f32 {
        if self.cells.is_empty() {
            return 0.0;
        }
        let mut used = 0f32;
        for c in &self.cells[..self.cells.len() - 1] {
            used += c.width;
        }
        table_width - used
    }

    fn has_bottom_border(&self) -> bool {
        !self.cells.is_empty() && self.cells.iter().all(|c| c.bottom_border.is_some())
    }
}

/// boxable's BaseTable: drawn from `y_start` on `page`, continuing on new pages (of the first page's size) from
/// `y_start_new_page` down to `bottom_margin`
#[derive(Clone, Debug)]
pub struct Table {
    y_start: f32,
    y_start_new_page: f32,
    page_top_margin: f32,
    page_bottom_margin: f32,
    width: f32,
    margin: f32,
    rows: Vec<Row>,
    header: Vec<usize>,
    page: usize,
    page_size: PageSize,
    provider_index: Option<usize>,
    remove_top_borders: bool,
    previous_drawn_row: Option<usize>,
    keep_together: usize,
}

impl Table {
    /// new BaseTable(yStart, yStartNewPage, bottomMargin, width, margin, document, currentPage, true, true)
    pub fn new(doc: &Document, y_start: f32, y_start_new_page: f32, bottom_margin: f32, width: f32, margin: f32, page: usize) -> Table {
        Table {
            y_start,
            y_start_new_page,
            page_top_margin: 0.0,
            page_bottom_margin: bottom_margin,
            width,
            margin,
            rows: Vec::new(),
            header: Vec::new(),
            page,
            page_size: doc.page_size(page),
            provider_index: None,
            remove_top_borders: false,
            previous_drawn_row: None,
            keep_together: 0,
        }
    }

    /// Not boxable: the first `rows` rows start on a new page when they do not fit together on this one (and the
    /// table is not at the top of a page already), so that a title row is not left alone at the end of a page
    pub fn set_keep_together(&mut self, rows: usize) {
        self.keep_together = rows;
    }

    pub fn width(&self) -> f32 {
        self.width
    }

    /// The page the table ended on
    pub fn current_page(&self) -> usize {
        self.page
    }

    /// createRow(height)
    pub fn create_row(&mut self, height: f32) -> &mut Row {
        self.rows.push(Row { height, cells: Vec::new(), width: self.width });
        self.rows.last_mut().unwrap()
    }

    /// The index of the last row created
    pub fn last_row_index(&self) -> usize {
        self.rows.len() - 1
    }

    pub fn row_mut(&mut self, index: usize) -> &mut Row {
        &mut self.rows[index]
    }

    /// addHeaderRow: a row repeated at the top of every page the table continues on
    pub fn add_header_row(&mut self, index: usize) {
        self.header.push(index);
    }

    /// getHeaderAndDataHeight
    pub fn height(&mut self, doc: &Document) -> Result<f32, MissingGlyph> {
        let mut h = 0f32;
        for r in &mut self.rows {
            h += r.get_height(doc)?;
        }
        Ok(h)
    }

    /// getMinimumHeight: the header rows and the row at the header rows' count (the first data row when the
    /// header rows come first)
    fn minimum_height(&mut self, doc: &Document) -> Result<f32, MissingGlyph> {
        let mut h = 0f32;
        let header = self.header.clone();
        for r in &header {
            h += self.rows[*r].get_height(doc)?;
        }
        if self.rows.len() > header.len() {
            h += self.rows[header.len()].get_height(doc)?;
        }
        Ok(h)
    }

    fn is_end_of_page(&self, free_space: f32) -> bool {
        self.y_start - free_space <= self.page_bottom_margin
    }

    /// The page part of pageBreak (DefaultPageProvider.nextPage): the next page, appended when it does not exist
    fn next_page(&mut self, doc: &mut Document) {
        self.y_start = self.y_start_new_page - self.page_top_margin;
        let index = match self.provider_index {
            None => doc.page_count(),
            Some(i) => i + 1,
        };
        self.provider_index = Some(index);
        if index >= doc.page_count() {
            doc.add_page(self.page_size);
        }
        self.page = index;
    }

    /// draw(): returns the y below the table
    pub fn draw(&mut self, doc: &mut Document) -> Result<f32, MissingGlyph> {
        let mut out = Optimized::new(doc.open_stream(self.page));
        self.remove_top_borders = false;
        self.previous_drawn_row = None;
        let keep = self.keep_together.min(self.rows.len());
        if keep > 0 {
            let mut h = 0f32;
            for r in 0..keep {
                h += self.rows[r].get_height(doc)?;
            }
            let full_page = self.y_start_new_page - self.page_top_margin - self.page_bottom_margin;
            let at_top = self.y_start >= self.y_start_new_page - self.page_top_margin;
            if self.is_end_of_page(h) && h <= full_page && !at_top {
                self.page_break(doc, &mut out);
            }
        }
        for r in 0..self.rows.len() {
            if !self.header.is_empty() && r == self.header[0] {
                let min_height = self.minimum_height(doc)?;
                if self.is_end_of_page(min_height) {
                    let total_page_height = self.y_start_new_page - self.page_top_margin - self.page_bottom_margin;
                    if min_height > total_page_height {
                        let available = self.y_start - self.page_bottom_margin;
                        if available < total_page_height - 10f32 {
                            self.page_break(doc, &mut out);
                        }
                    } else {
                        self.page_break(doc, &mut out);
                    }
                }
            }
            self.draw_row(r, doc, &mut out)?;
        }
        out.close(doc);
        Ok(self.y_start)
    }

    /// pageBreak: the content stream closed, the next page, a new stream
    fn page_break(&mut self, doc: &mut Document, out: &mut Optimized) {
        let old = std::mem::replace(out, Optimized::new(doc.open_stream(self.page)));
        old.close(doc);
        self.next_page(doc);
        *out = Optimized::new(doc.open_stream(self.page));
    }

    fn draw_row(&mut self, r: usize, doc: &mut Document, out: &mut Optimized) -> Result<(), MissingGlyph> {
        let is_header = self.header.contains(&r);
        if is_header {
            for c in &mut self.rows[r].cells {
                c.line_start = 0;
            }
        }
        loop {
            let row_height = self.rows[r].get_height(doc)?;
            let available = self.y_start - self.page_bottom_margin;
            let mut split = false;
            let mut height_to_draw = row_height;
            if row_height > available && !is_header {
                let full_page_height = self.y_start_new_page - self.page_top_margin - self.page_bottom_margin;
                let mut header_height = 0f32;
                for h in self.header.clone() {
                    header_height += self.rows[h].get_height(doc)?;
                }
                if row_height > full_page_height - header_height {
                    split = true;
                    height_to_draw = available;
                }
            }
            if !split && row_height > available && !is_header {
                self.page_break(doc, out);
                self.remove_top_borders = false;
                for h in self.header.clone() {
                    if h != r {
                        self.draw_row(h, doc, out)?;
                    }
                }
                continue;
            }
            if self.remove_top_borders {
                for c in &mut self.rows[r].cells {
                    c.top_border = None;
                }
            }
            self.fill_row_cell_colors(r, height_to_draw, out);
            self.draw_vertical_lines(r, height_to_draw, out);
            self.draw_cell_content(r, height_to_draw, doc, out)?;
            self.previous_drawn_row = Some(r);
            self.y_start -= height_to_draw;
            if split {
                self.page_break(doc, out);
                self.remove_top_borders = true;
                self.previous_drawn_row = None;
                for h in self.header.clone() {
                    self.draw_row(h, doc, out)?;
                }
            } else {
                self.remove_top_borders = self.rows[r].has_bottom_border();
                break;
            }
        }
        Ok(())
    }

    fn fill_row_cell_colors(&self, r: usize, row_height: f32, out: &mut Optimized) {
        let row = &self.rows[r];
        let mut x = self.margin;
        let n = row.cells.len();
        for (i, cell) in row.cells.iter().enumerate() {
            let layout_width = if i + 1 < n { cell.width } else { row.last_cell_layout_width(self.width) };
            let fill_width = cell.width.min(layout_width);
            if let Some(fill) = cell.fill {
                let y = self.y_start - row_height;
                let h = row_height - bw(&cell.top_border);
                out.set_non_stroking_color(fill);
                out.add_rect(x, y, fill_width, h);
                out.fill();
            }
            x += layout_width;
        }
    }

    /// drawLine with PDStreamUtils.setLineStyles
    fn line(&self, out: &mut Optimized, x1: f32, y1: f32, x2: f32, y2: f32, style: LineStyle) {
        if style.width > 0.0 {
            out.set_stroking_color(style.color);
            out.set_line_width(style.width);
            out.set_line_cap(0);
            out.set_no_dash();
            out.move_to(x1, y1);
            out.line_to(x2, y2);
            out.stroke();
        }
    }

    fn draw_vertical_lines(&self, r: usize, row_height: f32, out: &mut Optimized) {
        let x_start = self.margin;
        // the bottom borders of the previous row continued where its cells do not reach the table's width
        if let Some(p) = self.previous_drawn_row {
            let prev = &self.rows[p];
            if !prev.cells.is_empty() {
                let prev_end = self.margin + prev.total_cells_width();
                let table_end = self.margin + self.width;
                if prev_end < table_end {
                    if let Some(b) = prev.cells[prev.cells.len() - 1].bottom_border {
                        self.line(out, prev_end, self.y_start + b.width / 2f32, table_end, self.y_start + b.width / 2f32, b);
                    }
                }
            }
        }
        let row = &self.rows[r];
        let n = row.cells.len();
        let mut x = x_start;
        for (i, cell) in row.cells.iter().enumerate() {
            let cell_width = if i + 1 < n { cell.width } else { row.last_cell_layout_width(self.width) };
            let draw_left = i == 0 || row.cells[i - 1].right_border.is_none();
            self.draw_cell_borders(out, row_height, cell, x, draw_left);
            x += cell_width;
        }
    }

    fn draw_cell_borders(&self, out: &mut Optimized, row_height: f32, cell: &Cell, x: f32, draw_left: bool) {
        let y_end = self.y_start - row_height;
        if let Some(t) = cell.top_border {
            let y = self.y_start - t.width / 2f32;
            self.line(out, x, y, x + cell.width, y, t);
        }
        if let Some(rb) = cell.right_border {
            let xx = x + cell.width - rb.width / 2f32;
            self.line(out, xx, self.y_start - bw(&cell.top_border), xx, y_end, rb);
        }
        if let Some(b) = cell.bottom_border {
            let y = y_end + b.width / 2f32;
            self.line(out, x, y, x + cell.width - bw(&cell.right_border), y, b);
        }
        if draw_left {
            if let Some(l) = cell.left_border {
                let xx = x + l.width / 2f32;
                self.line(out, xx, self.y_start, xx, y_end + bw(&cell.bottom_border), l);
            }
        }
    }

    fn draw_cell_content(&mut self, r: usize, row_height: f32, doc: &mut Document, out: &mut Optimized) -> Result<(), MissingGlyph> {
        let mut cursor_x = self.margin;
        let y_start = self.y_start;
        for ci in 0..self.rows[r].cells.len() {
            let cell_start_x = cursor_x;
            // the row's height as boxable reads it while drawing this cell: the cells drawn before it are finished
            // (their text no longer counts), which matters for bottom or middle alignment
            let full_row_height = self.rows[r].get_height(doc)?;
            let (font, size, color, width) = {
                let c = &self.rows[r].cells[ci];
                (c.font, c.font_size, c.text_color, c.width)
            };
            let metrics = BoxMetrics::of(doc.font(font));
            let text_height = self.rows[r].cells[ci].text_height(doc)?;
            let cell = &mut self.rows[r].cells[ci];
            let mut cursor_y = y_start - cell.top_padding - metrics.height(size) - metrics.descent(size) - bw(&cell.top_border);
            cursor_x += cell.left_padding + bw(&cell.left_border);
            let inner_width = cell.inner_width();
            let vertical_free = cell.inner_height(full_row_height) - text_height;
            let (align, valign, line_spacing, line_start) = (cell.align, cell.valign, cell.line_spacing, cell.line_start);
            let paragraph = cell.paragraph(doc)?.clone();
            let horizontal_free = inner_width - paragraph.max_line_width;
            match align {
                HAlign::Center => cursor_x += horizontal_free / 2f32,
                HAlign::Right => cursor_x += horizontal_free,
                HAlign::Left => {}
            }
            match valign {
                VAlign::Middle => cursor_y -= vertical_free / 2f32,
                VAlign::Bottom => cursor_y -= vertical_free,
                VAlign::Top => {}
            }
            let line_start_x = cursor_x;
            out.set_non_stroking_color(color);
            let mut finished = true;
            let mut new_line_start = line_start;
            for (li, line) in paragraph.lines.iter().enumerate() {
                if li < line_start {
                    continue;
                }
                if cursor_y <= y_start - row_height {
                    new_line_start = li;
                    finished = false;
                    break;
                }
                let free = paragraph.max_line_width - line.width;
                cursor_x = line_start_x;
                match align {
                    HAlign::Center => cursor_x += free / 2f32,
                    HAlign::Right => cursor_x += free,
                    HAlign::Left => {}
                }
                for token in &line.tokens {
                    out.set_font(font, size);
                    out.new_line_at(cursor_x, cursor_y + 0f32);
                    out.show_text(doc, font, token)?;
                    cursor_x += doc.font(font).string_width(token)? / 1000f32 * size;
                }
                cursor_y -= line.height * line_spacing;
            }
            let cell = &mut self.rows[r].cells[ci];
            cell.line_start = if finished { usize::MAX } else { new_line_start };
            cursor_x = cell_start_x + width;
        }
        Ok(())
    }
}

/// boxable's PageContentStreamOptimized: text mode stays open until a path operator or a change of line style,
/// text positions are offsets from the previous one, and repeated settings are left out
struct Optimized {
    s: Stream,
    text_mode: bool,
    cx: f32,
    cy: f32,
    font: Option<(FontId, f32)>,
    fill: Option<Rgb>,
    stroke: Option<Rgb>,
    line_width: f32,
    line_cap: i32,
    dash_set: bool,
}

impl Optimized {
    fn new(s: Stream) -> Optimized {
        Optimized { s, text_mode: false, cx: 0.0, cy: 0.0, font: None, fill: None, stroke: None, line_width: -1.0, line_cap: -1, dash_set: false }
    }
    fn close(mut self, doc: &mut Document) {
        self.end_text();
        doc.close_stream(self.s);
    }
    fn begin_text(&mut self) {
        if !self.text_mode {
            self.s.begin_text();
            self.text_mode = true;
            self.cx = 0.0;
            self.cy = 0.0;
        }
    }
    fn end_text(&mut self) {
        if self.text_mode {
            self.s.end_text();
            self.text_mode = false;
        }
    }
    fn set_font(&mut self, font: FontId, size: f32) {
        if self.font != Some((font, size)) {
            self.s.set_font(font, size);
            self.font = Some((font, size));
        }
    }
    fn show_text(&mut self, doc: &mut Document, font: FontId, text: &str) -> Result<(), MissingGlyph> {
        self.begin_text();
        self.s.show_text(doc, font, text)
    }
    fn new_line_at(&mut self, tx: f32, ty: f32) {
        self.begin_text();
        let (dx, dy) = (tx - self.cx, ty - self.cy);
        self.s.new_line_at_offset(dx, dy);
        self.cx = tx;
        self.cy = ty;
    }
    fn set_stroking_color(&mut self, c: Rgb) {
        if self.stroke != Some(c) {
            self.s.set_stroking_color(c);
            self.stroke = Some(c);
        }
    }
    fn set_non_stroking_color(&mut self, c: Rgb) {
        if self.fill != Some(c) {
            self.s.set_non_stroking_color(c);
            self.fill = Some(c);
        }
    }
    fn add_rect(&mut self, x: f32, y: f32, w: f32, h: f32) {
        self.end_text();
        self.s.add_rect(x, y, w, h);
    }
    fn move_to(&mut self, x: f32, y: f32) {
        self.end_text();
        self.s.move_to(x, y);
    }
    fn line_to(&mut self, x: f32, y: f32) {
        self.end_text();
        self.s.line_to(x, y);
    }
    fn stroke(&mut self) {
        self.end_text();
        self.s.stroke();
    }
    fn fill(&mut self) {
        self.end_text();
        self.s.fill();
    }
    fn set_line_width(&mut self, w: f32) {
        if w != self.line_width {
            self.end_text();
            self.s.set_line_width(w);
            self.line_width = w;
        }
    }
    fn set_line_cap(&mut self, cap: i32) {
        if cap != self.line_cap {
            self.end_text();
            self.s.set_line_cap(cap);
            self.line_cap = cap;
        }
    }
    fn set_no_dash(&mut self) {
        if !self.dash_set {
            self.end_text();
            self.s.set_no_dash();
            self.dash_set = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens() {
        let t = tokenize("6 (ARID1B, KRAS");
        let texts: Vec<String> = t
            .iter()
            .map(|t| match t {
                Token::Text(s) => s.clone(),
                Token::PossibleWrapPoint => "|".into(),
            })
            .collect();
        assert_eq!(texts, vec!["6 ", "|", "(ARID1B,", "|", " ", "|", "KRAS", "|"]);
        assert_eq!(tokenize(""), vec![Token::PossibleWrapPoint]);
        // a wrap character at the end is no wrap point
        let t = tokenize("a.");
        assert_eq!(t, vec![Token::Text("a.".into()), Token::PossibleWrapPoint]);
        assert!(has_markup("x <b>y</b>"));
        assert!(has_markup("</li>"));
        assert!(!has_markup("c.35G>T <1%"));
    }
}
