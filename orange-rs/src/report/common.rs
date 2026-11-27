//! ReportResources (fonts, colours, text styles), Cells and Tables (styled boxable cells and the standard tables),
//! and TableCommon (column names, number formats, column widths).

use anyhow::{anyhow, Result};
use hmf_common::fmt::{decimal_format, decimal_format_fixed, java_fixed, java_sci};
use hmf_plot::pdf::table::{LineStyle, Table, VAlign};
use hmf_plot::pdf::{FontId, Rgb};

use super::DocCtx;
use crate::compat::java;

/// OR15: a table's title row is an ordinary first row, so boxable's rule that keeps the header rows with the first
/// data row leaves the title at the end of a page when the rest moves to the next; corrected: the title, the header
/// row and the first data row (or the title and its one row) start on a new page together
fn keep_title_with_rows(table: &mut Table, rows: usize) {
    if !java("OR15") {
        table.set_keep_together(rows);
    }
}

pub const NOT_AVAILABLE: &str = "NA";
pub const NONE: &str = "NONE";

pub const FRONT_CIRCOS_IMAGE_HEIGHT: f32 = 400.0;
pub const FULL_PAGE_IMAGE_WIDTH: f32 = 750.0;
pub const FULL_PAGE_IMAGE_HEIGHT: f32 = 430.0;
pub const MAX_TABLE_WIDTH: f64 = 450.0;
pub const MAX_TABLE_COLUMN_BUFFER: f64 = 5.0;

pub const PALETTE_WHITE: Rgb = Rgb(255, 255, 255);
pub const PALETTE_DARK_GREY: Rgb = Rgb(39, 47, 50);
pub const PALETTE_MID_GREY: Rgb = Rgb(101, 106, 108);
pub const PALETTE_SMOKE_GREY: Rgb = Rgb(240, 240, 240);
pub const PALETTE_ORANGE: Rgb = Rgb(242, 139, 31);
pub const PALETTE_ORANGE_1: Rgb = Rgb(255, 165, 0);
pub const PALETTE_ORANGE_2: Rgb = Rgb(235, 155, 0);
pub const PALETTE_ORANGE_3: Rgb = Rgb(215, 145, 0);
pub const PALETTE_ORANGE_4: Rgb = Rgb(195, 135, 0);
pub const PALETTE_ORANGE_5: Rgb = Rgb(175, 125, 0);
pub const PALETTE_ORANGE_6: Rgb = Rgb(155, 115, 0);

/// ReportResources.TextStyle
#[derive(Clone, Copy, Debug)]
pub struct TextStyle {
    pub font: FontId,
    pub size: f32,
    pub color: Rgb,
}

/// ReportResources: the two fonts and the styles
#[derive(Clone, Copy, Debug)]
pub struct Resources {
    pub regular: FontId,
    pub bold: FontId,
}

impl Resources {
    pub fn chapter_title_style(&self) -> TextStyle {
        TextStyle { font: self.bold, size: 10.0, color: PALETTE_ORANGE }
    }
    pub fn table_title_style(&self) -> TextStyle {
        TextStyle { font: self.bold, size: 8.0, color: PALETTE_ORANGE }
    }
    pub fn table_header_style(&self) -> TextStyle {
        TextStyle { font: self.regular, size: 7.0, color: PALETTE_MID_GREY }
    }
    pub fn table_content_style(&self) -> TextStyle {
        TextStyle { font: self.regular, size: 7.0, color: PALETTE_DARK_GREY }
    }
    pub fn key_style(&self) -> TextStyle {
        TextStyle { font: self.regular, size: 7.0, color: PALETTE_MID_GREY }
    }
    pub fn value_style(&self) -> TextStyle {
        TextStyle { font: self.regular, size: 7.0, color: PALETTE_MID_GREY }
    }
    pub fn page_number_style(&self) -> TextStyle {
        TextStyle { font: self.bold, size: 7.0, color: PALETTE_ORANGE }
    }
    pub fn disclaimer_style(&self) -> TextStyle {
        TextStyle { font: self.regular, size: 6.0, color: PALETTE_MID_GREY }
    }
    pub fn qc_warning_style(&self) -> TextStyle {
        TextStyle { font: self.bold, size: 7.0, color: PALETTE_DARK_GREY }
    }
    pub fn side_panel_label_style(&self) -> TextStyle {
        TextStyle { font: self.bold, size: 7.0, color: PALETTE_WHITE }
    }
}

// ------------------------------------------------------------------------------------------------
// Cells

pub const COMPACT_ROW_HEIGHT: f32 = 15.0;
pub const HEADER_ROW_HEIGHT: f32 = 16.0;

/// The cell styles of Cells
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Style {
    Header,
    Content,
    Key,
    Value,
    Warning,
    Title,
}

/// createCell on the table's row `row` with a style applied (Cells.apply*Style: compact padding of 2, no borders
/// but a thin grey line under content cells, header text bottom-aligned)
pub fn add_cell(res: &Resources, table: &mut Table, row: usize, pct: f32, text: &str, style: Style) -> Result<()> {
    let st = match style {
        Style::Header => res.table_header_style(),
        Style::Content => res.table_content_style(),
        Style::Key => res.key_style(),
        Style::Value => res.value_style(),
        Style::Warning => res.qc_warning_style(),
        Style::Title => res.table_title_style(),
    };
    let text = if style == Style::Header { text.to_uppercase() } else { text.to_string() };
    if hmf_plot::pdf::table::has_markup(&text) {
        // boxable reads <b>, <i>, <p>, <br> ... in cell text as markup; the port draws them as text
        hmf_common::warn!("PDF cell text with markup drawn as text (Java's boxable would format it): {text}");
    }
    let cell = table.row_mut(row).create_cell(pct, &text, st.font).map_err(|e| anyhow!("java.lang.IllegalArgumentException: {e}"))?;
    cell.set_font(st.font, st.size, st.color);
    let bottom = if style == Style::Content { Some(LineStyle { color: PALETTE_MID_GREY, width: 0.25 }) } else { None };
    cell.set_borders(None, bottom, None, None);
    if style == Style::Header {
        cell.valign = VAlign::Bottom;
    }
    cell.set_padding(2.0, 2.0, 2.0, 2.0);
    Ok(())
}

/// Cells.addRow: a row of content cells; returns the row's index
pub fn add_row(res: &Resources, table: &mut Table, pcts: &[f32], values: &[String]) -> Result<usize> {
    table.create_row(COMPACT_ROW_HEIGHT);
    let r = table.last_row_index();
    for (i, v) in values.iter().enumerate() {
        let p = if i < pcts.len() { pcts[i] } else { pcts[pcts.len() - 1] };
        add_cell(res, table, r, p, v, Style::Content)?;
    }
    Ok(r)
}

/// ReportResources.shadeCandidateCells
pub fn shade_row(table: &mut Table, row: usize) {
    for c in &mut table.row_mut(row).cells {
        c.fill = Some(PALETTE_SMOKE_GREY);
    }
}

// ------------------------------------------------------------------------------------------------
// Tables

/// Tables.createNonContent: a title row and one content row
pub(crate) fn keep_title(table: &mut Table, rows: usize) {
    keep_title_with_rows(table, rows);
}

pub fn create_non_content(ctx: &mut DocCtx, title: &str, width: f32, value: &str) -> Result<Table> {
    let res = ctx.res;
    let mut table = ctx.create_table(width)?;
    table.create_row(HEADER_ROW_HEIGHT);
    let r = table.last_row_index();
    add_cell(&res, &mut table, r, 100.0, title, Style::Title)?;
    table.create_row(COMPACT_ROW_HEIGHT);
    let r = table.last_row_index();
    add_cell(&res, &mut table, r, 100.0, value, Style::Content)?;
    keep_title_with_rows(&mut table, 2);
    Ok(table)
}

/// TableCommon.createEmptyTable: the title and "NONE"
pub fn create_empty_table(ctx: &mut DocCtx, title: &str, width: f32) -> Result<Table> {
    create_non_content(ctx, title, width, NONE)
}

fn header_row(res: &Resources, table: &mut Table, widths: &[f32], headers: &[String]) -> Result<()> {
    let mut total = 0f32;
    for w in widths {
        total += *w;
    }
    table.create_row(HEADER_ROW_HEIGHT);
    let r = table.last_row_index();
    for (i, h) in headers.iter().enumerate() {
        let pct = (widths[i] / total) * 100f32;
        add_cell(res, table, r, pct, h, Style::Header)?;
    }
    table.add_header_row(r);
    Ok(())
}

/// Tables.createContent: a header row (repeated on new pages)
pub fn create_content(ctx: &mut DocCtx, width: f32, widths: &[f32], headers: &[String]) -> Result<Table> {
    let res = ctx.res;
    let mut table = ctx.create_table(width)?;
    header_row(&res, &mut table, widths, headers)?;
    Ok(table)
}

/// TableCommon.createStandardTable (Tables.createWithTitle): a title row, then a header row
pub fn create_standard_table(ctx: &mut DocCtx, title: &str, width: f32, widths: &[f32], headers: &[String]) -> Result<Table> {
    let res = ctx.res;
    let mut table = ctx.create_table(width)?;
    table.create_row(HEADER_ROW_HEIGHT);
    let r = table.last_row_index();
    add_cell(&res, &mut table, r, 100.0, title, Style::Title)?;
    header_row(&res, &mut table, widths, headers)?;
    keep_title_with_rows(&mut table, 3);
    Ok(table)
}

/// TableCommon.toPercentages (float)
pub fn to_percentages(widths: &[f32]) -> Vec<f32> {
    let mut total = 0f32;
    for w in widths {
        total += *w;
    }
    widths.iter().map(|w| (w / total) * 100f32).collect()
}

// ------------------------------------------------------------------------------------------------
// TableCommon: number formats

/// formatSingleDigitDecimal: DecimalFormat("0.0")
pub fn single_digit(v: f64) -> String {
    decimal_format_fixed(v, 1)
}

/// formatTwoDigitDecimal: DecimalFormat("0.00")
pub fn two_digit(v: f64) -> String {
    decimal_format_fixed(v, 2)
}

/// formatPercentage: DecimalFormat("0'%'") of v * 100
pub fn percentage(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    format!("{}%", decimal_format(v * 100.0, 0))
}

/// formatPercentageField: DecimalFormat("#'%'") of v * 100 (the same digits as "0'%'")
pub fn percentage_field(v: f64) -> String {
    percentage(v)
}

/// formatTpmField
pub fn tpm_field(v: f64) -> String {
    single_digit(v)
}

/// formatFoldChangeField
pub fn fold_change_field(v: f64) -> String {
    if v > 1000.0 { ">1000".into() } else { single_digit(v) }
}

/// formatSupportField: "%d/%d"
pub fn support_field(fragments: i32, depth: i32) -> String {
    format!("{fragments}/{depth}")
}

/// String.format("%.<n>f")
pub fn fixed(v: f64, n: usize) -> String {
    java_fixed(v, n)
}

/// String.format("%4.2e", value) for a Double: "null" truncated to the precision ("nu") when null, padded to 4
pub fn sci_4_2(v: Option<f64>) -> String {
    let s = match v {
        // Formatter prints the sign of -0.0, which hmf-common's java_sci leaves out
        Some(v) if v == 0.0 && v.is_sign_negative() => format!("-{}", java_sci(0.0, 2)),
        Some(v) => java_sci(v, 2),
        None => "nu".to_string(),
    };
    format!("{s:>4}")
}

/// TableCommon.zeroPrefixed: "0" before a chromosome number below 10, for sorting
pub fn zero_prefixed(location: &str) -> String {
    let mut arm_start = location.find('q');
    if arm_start.is_none() {
        arm_start = location.find('p');
    }
    let chromosome = match arm_start {
        Some(i) if i > 0 => &location[..i],
        _ => location,
    };
    match java_parse_int(chromosome) {
        Some(i) if i < 10 => format!("0{location}"),
        _ => location.to_string(),
    }
}

/// Integer.parseInt: an optional sign and decimal digits
fn java_parse_int(s: &str) -> Option<i32> {
    let digits = s.strip_prefix('-').or_else(|| s.strip_prefix('+')).unwrap_or(s);
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    s.parse::<i32>().ok()
}

// ------------------------------------------------------------------------------------------------
// TableCommon: column widths scaled to the text (scaleTableWidths, convertWidths)

/// TableHeaderData
#[derive(Clone, Debug)]
pub struct HeaderData {
    pub text: String,
    pub wrap_at_percent: f64,
    pub can_wrap_data: bool,
    pub header_width: f64,
    pub data_width: f64,
    pub excess_width: f64,
}

impl HeaderData {
    pub fn new(text: &str, wrap_at_percent: f64, can_wrap_data: bool) -> HeaderData {
        HeaderData { text: text.to_string(), wrap_at_percent, can_wrap_data, header_width: 0.0, data_width: 0.0, excess_width: 0.0 }
    }
    pub fn adjusted_width(&self) -> f64 {
        self.header_width.max(self.data_width) + self.excess_width
    }
}

/// calculateTextWidth: the text's width in the header (upper case) or content style, in float
fn text_width(ctx: &DocCtx, text: &str, header: bool) -> f32 {
    let st = if header { ctx.res.table_header_style() } else { ctx.res.table_content_style() };
    let t = if header { text.to_uppercase() } else { text.to_string() };
    match ctx.doc.font(st.font).string_width(&t) {
        Ok(w) => (w / 1000f32) * st.size,
        Err(_) => 0.0,
    }
}

/// scaleTableWidths
pub fn scale_table_widths(ctx: &DocCtx, headers: &mut [HeaderData], rows: &[Vec<String>], max_table_width: f64) {
    scale_table_widths_by(headers, rows, max_table_width, |t, header| text_width(ctx, t, header));
}

/// scaleTableWidths with the text widths of `width(text, is_header)` (Java without ReportResources: the length)
pub fn scale_table_widths_by(headers: &mut [HeaderData], rows: &[Vec<String>], max_table_width: f64, width: impl Fn(&str, bool) -> f32) {
    let n = headers.len();
    let mut max_header_data_length = 0f64;
    for c in 0..n {
        for row in rows {
            let w = width(&row[c], false);
            headers[c].data_width = headers[c].data_width.max(w as f64);
        }
        headers[c].header_width = width(&headers[c].text, true) as f64;
        max_header_data_length += headers[c].adjusted_width();
    }
    if max_header_data_length > max_table_width {
        let mut required = max_header_data_length - max_table_width;
        for h in headers.iter_mut() {
            if h.wrap_at_percent > 0.0 {
                let header_length = h.header_width;
                let min_header_length = (header_length * h.wrap_at_percent) as f32 as f64;
                if min_header_length >= h.data_width {
                    let reduction = header_length - min_header_length;
                    h.header_width = min_header_length;
                    required -= reduction;
                }
                if required <= 0.0 {
                    break;
                }
            }
        }
        if required > 0.0 {
            for h in headers.iter_mut() {
                let max_data_length = h.data_width;
                if h.can_wrap_data {
                    let new_length = h.header_width.max(max_data_length * 0.5);
                    let reduction = h.data_width - new_length;
                    h.data_width = new_length;
                    required -= reduction;
                    if required <= 0.0 {
                        break;
                    }
                }
            }
        }
    }
    let total: f64 = headers.iter().map(|h| h.adjusted_width()).sum();
    if total < max_table_width {
        let mut excess = max_table_width - total;
        excess = excess.min(MAX_TABLE_COLUMN_BUFFER * n as f64);
        let per_column = excess / n as f64;
        for h in headers.iter_mut() {
            h.excess_width = per_column;
        }
    }
}

/// convertWidths: (header texts, column widths, column percentages)
pub fn convert_widths(headers: &[HeaderData]) -> (Vec<String>, Vec<f32>, Vec<f32>) {
    let total: f64 = headers.iter().map(|h| h.adjusted_width()).sum();
    let mut texts = Vec::new();
    let mut widths = Vec::new();
    let mut pcts = Vec::new();
    for h in headers {
        let w = h.adjusted_width() as f32;
        widths.push(w);
        texts.push(h.text.clone());
        pcts.push((w as f64 / total * 100.0) as f32);
    }
    (texts, widths, pcts)
}

/// The cell padding of the table cells (2 left and 2 right), which scaleTableWidths leaves out
const CELL_PADDING: f64 = 4.0;

/// Whether a column cannot hold its widest text on one line where it must not wrap: a header without
/// WrapAtPercent, or data that cannot wrap (boxable then breaks the text inside, "1." over "00")
fn breaks_values(headers: &[HeaderData], pcts: &[f32], table_width: f32) -> bool {
    headers.iter().zip(pcts).any(|(h, p)| {
        let inner = (table_width * (p / 100f32)) as f64 - CELL_PADDING;
        (h.wrap_at_percent == 0.0 && inner < h.header_width) || (!h.can_wrap_data && inner < h.data_width)
    })
}

/// OR14 corrected: column widths that hold each column's text and padding: the text widths of scaleTableWidths
/// plus the padding, the space left shared in proportion to the text widths; when the table is too narrow even
/// so, the columns whose data or header may wrap give up width first (down to their header, or to 60% of it), then
/// all columns shrink in proportion as Java's do
fn padded_widths(headers: &[HeaderData], table_width: f32) -> Vec<f32> {
    let tw = table_width as f64;
    let need: Vec<f64> = headers.iter().map(|h| h.adjusted_width() + CELL_PADDING).collect();
    let total_need: f64 = need.iter().sum();
    let mut w: Vec<f64> = if total_need <= tw {
        let sum_adj: f64 = headers.iter().map(|h| h.adjusted_width()).sum();
        let slack = tw - total_need;
        need.iter().zip(headers).map(|(n, h)| n + if sum_adj > 0.0 { slack * h.adjusted_width() / sum_adj } else { slack / headers.len() as f64 }).collect()
    } else {
        let floor: Vec<f64> = headers
            .iter()
            .zip(&need)
            .map(|(h, n)| {
                let mut f = *n;
                if h.can_wrap_data {
                    f = f.min(h.header_width + CELL_PADDING);
                }
                if h.wrap_at_percent > 0.0 {
                    f = f.min(h.header_width.max(if h.can_wrap_data { 0.0 } else { h.data_width }) + CELL_PADDING);
                }
                f
            })
            .collect();
        let reducible: f64 = need.iter().zip(&floor).map(|(n, f)| n - f).sum();
        let overflow = total_need - tw;
        let take = overflow.min(reducible);
        need.iter().zip(&floor).map(|(n, f)| if reducible > 0.0 { n - take * (n - f) / reducible } else { *n }).collect()
    };
    let sum: f64 = w.iter().sum();
    if sum > tw {
        for x in &mut w {
            *x *= tw / sum;
        }
    }
    w.iter().map(|x| *x as f32).collect()
}

/// A table whose columns are sized to its text (the scaled tables): title, header row, the rows, shaded where asked.
/// OR14: Java sizes the columns to their text without the cells' padding, so in a table whose text is wider than
/// about 450 pt the narrow columns cannot hold their values: boxable breaks them inside ("1." over "00" for an allele
/// frequency of 1.00, the header JCN as "JC" over "N"); corrected (only where Java breaks a value): widths with the
/// padding
pub fn scaled_table(ctx: &mut DocCtx, title: &str, width: f32, mut headers: Vec<HeaderData>, rows: Vec<Vec<String>>, shade: &[bool]) -> Result<Table> {
    scale_table_widths(ctx, &mut headers, &rows, MAX_TABLE_WIDTH);
    let (texts, mut widths, mut pcts) = convert_widths(&headers);
    if !java("OR14") && breaks_values(&headers, &pcts, width) {
        widths = padded_widths(&headers, width);
        pcts = to_percentages(&widths);
    }
    let mut table = create_standard_table(ctx, title, width, &widths, &texts)?;
    let res = ctx.res;
    for (i, row) in rows.iter().enumerate() {
        let r = add_row(&res, &mut table, &pcts, row)?;
        if shade.get(i).copied().unwrap_or(false) {
            shade_row(&mut table, r);
        }
    }
    Ok(table)
}

/// A table with fixed relative column widths (addEntry widths): title, header row, the rows
pub fn fixed_table(ctx: &mut DocCtx, title: &str, width: f32, widths: &[f32], headers: &[&str], rows: Vec<Vec<String>>, shade: &[bool]) -> Result<Table> {
    let headers: Vec<String> = headers.iter().map(|s| s.to_string()).collect();
    let mut table = create_standard_table(ctx, title, width, widths, &headers)?;
    let pcts = to_percentages(widths);
    let res = ctx.res;
    for (i, row) in rows.iter().enumerate() {
        let r = add_row(&res, &mut table, &pcts, row)?;
        if shade.get(i).copied().unwrap_or(false) {
            shade_row(&mut table, r);
        }
    }
    Ok(table)
}
