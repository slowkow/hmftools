//! The visualiser's drawing layer: a writer that serialises shapes as JFreeSVG 5.0.5's
//! `SVGGraphics2D.getSVGElement()` does, the font metrics Java's AWT reports for the "SansSerif"
//! logical font, and the ports of `hmf-common` `vis/SvgRender` (base sequences, coordinates and
//! chromosome labels), `SvgUtil`, `BaseSeqViewModel`, `BaseViewModel` and `ColorUtil`.

use super::{VisError, VisResult};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color(pub u8, pub u8, pub u8);

pub const BLACK: Color = Color(0, 0, 0);
pub const WHITE: Color = Color(255, 255, 255);
pub const LIGHT_GRAY: Color = Color(192, 192, 192);
pub const PURPLE: Color = Color(128, 0, 128);
pub const FORWARD_STRAND_COLOR: Color = Color(81, 144, 207);
pub const REVERSE_STRAND_COLOR: Color = Color(242, 167, 121);
const INSERT_COLOR: Color = Color(153, 50, 204);

pub const BOX_PADDING: i32 = 1;
pub const BASE_BOX_SIZE: f64 = 30.0;
const DEL_CONNECTOR_BOX_PROPORTION: f64 = 0.2;
const BOUNDARY_BOX_PROPORTION: f64 = 0.2;
const INDEL_FONT_SIZE: i32 = 18;
const BASE_FONT_SIZE: i32 = 26;
pub const COORD_FONT_SIZE: i32 = 30;
const MAX_BASEQ_SHADING_CUTTOFF: i32 = 37;
const INSERT_LABEL: &str = "INS";

/// Math.round(double)
fn jround(v: f64) -> i64 {
    if v.is_nan() {
        return 0;
    }
    (v + 0.5).floor() as i64
}

/// ColorUtil.interpolateColors
pub fn interpolate_colors(start: Color, end: Color, proportion: f64) -> Color {
    let p = proportion.clamp(0.0, 1.0);
    let mix = |a: u8, b: u8| -> u8 { jround((1.0 - p) * a as f64 + p * b as f64).clamp(0, 255) as u8 };
    Color(mix(start.0, end.0), mix(start.1, end.1), mix(start.2, end.2))
}

/// ColorUtil.lighten
pub fn lighten(color: Color, proportion: f64) -> Color {
    interpolate_colors(color, WHITE, proportion)
}

/// SvgRender.OVERLAPPING_FRAGMENT_BORDER_COLOR
pub fn overlapping_fragment_border_color() -> Color {
    lighten(PURPLE, 0.5f32 as f64)
}

/// hmf-common BaseRegion: closed interval
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BaseRegion {
    pub start: i32,
    pub end: i32,
}

impl BaseRegion {
    pub fn new(start: i32, end: i32) -> BaseRegion {
        BaseRegion { start, end }
    }
    pub fn base_length(&self) -> i32 {
        self.end - self.start + 1
    }
    pub fn overlaps(&self, other: &BaseRegion) -> bool {
        !(self.start > other.end || self.end < other.start)
    }
}

/// Double.toString
pub fn java_double(d: f64) -> String {
    if d.is_nan() {
        return "NaN".into();
    }
    if d.is_infinite() {
        return if d > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    if d == 0.0 {
        return if d.is_sign_negative() { "-0.0".into() } else { "0.0".into() };
    }
    let a = d.abs();
    let mut out = String::new();
    if d < 0.0 {
        out.push('-');
    }
    if (1e-3..1e7).contains(&a) {
        // Rust's shortest round-trip form; Display never uses an exponent
        let s = format!("{a}");
        out.push_str(&s);
        if !s.contains('.') {
            out.push_str(".0");
        }
    } else {
        let s = format!("{a:e}");
        let (mant, exp) = s.split_once('e').unwrap();
        out.push_str(mant);
        if !mant.contains('.') {
            out.push_str(".0");
        }
        out.push('E');
        out.push_str(exp);
    }
    out
}

/// Advance widths of the printable ASCII characters (32 to 126) in 1/2048 em, as Java's AWT reports
/// them for the logical font "SansSerif" on macOS (which maps it to Lucida Grande). Dumped from Java at
/// sizes 18, 26 and 30; the advances scale exactly with the size.
const ADVANCE_UNITS: [u16; 95] = [
    648, 648, 765, 1295, 1295, 1369, 1428, 469, 666, 666, 987, 1628, 648, 1185, 648, 1074, 1295, 1295, 1295, 1295, 1295,
    1295, 1295, 1295, 1295, 1295, 648, 648, 1628, 1628, 1628, 864, 1758, 1413, 1178, 1417, 1534, 1110, 1098, 1480, 1505,
    590, 637, 1337, 1092, 1764, 1513, 1591, 1132, 1591, 1295, 1103, 1295, 1419, 1339, 1752, 1282, 1276, 1238, 666, 1074,
    666, 1295, 1024, 1258, 1131, 1289, 1049, 1289, 1141, 753, 1277, 1271, 592, 623, 1197, 592, 1912, 1271, 1258, 1289,
    1289, 838, 1044, 766, 1271, 1060, 1578, 1256, 1070, 1174, 666, 765, 666, 1295,
];
/// ascent + descent + leading of the same font, in 1/2048 em
const LINE_HEIGHT_UNITS: f32 = 2412.0;

/// FontMetrics.getDescent for the three sizes the visualiser uses
pub fn font_descent(size: i32) -> i32 {
    match size {
        18 => 4,
        26 => 6,
        30 => 7,
        // FontDesignMetrics rounds the descent (432/2048 em) up from 0.05
        _ => (432.0 * size as f32 / 2048.0 + 0.95) as i32,
    }
}

/// SvgUtil.getStringBounds: (width, height); the height has the integer descent added, as in Java.
/// The x of the bounds is always zero. Characters outside printable ASCII are given the width of a space.
pub fn string_bounds(size: i32, text: &str) -> (f64, f64) {
    let mut width = 0f32;
    for c in text.chars() {
        let units = if (' '..='~').contains(&c) { ADVANCE_UNITS[c as usize - 32] } else { ADVANCE_UNITS[0] };
        width += units as f32 * size as f32 / 2048.0;
    }
    let height = (LINE_HEIGHT_UNITS * size as f32 / 2048.0) as f64 + font_descent(size) as f64;
    (width as f64, height)
}

/// SVGUtils.escapeForXML
fn escape_xml(text: &str, out: &mut String) {
    for c in text.chars() {
        match c {
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '&' => out.push_str("&amp;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
}

/// The part of `SVGGraphics2D` the visualiser uses. The transform is an `AffineTransform` restricted to
/// the operations SvgRender performs, each done with the arithmetic `AffineTransform` uses for it.
pub struct Canvas {
    pub width: f64,
    pub height: f64,
    body: String,
    m00: f64,
    m10: f64,
    m01: f64,
    m11: f64,
    m02: f64,
    m12: f64,
    /// set by a quarter turn: AffineTransform is then in its "shear" state
    sheared: bool,
    color: Color,
    font_size: i32,
}

#[derive(Clone, Copy)]
pub struct Transform(f64, f64, f64, f64, f64, f64, bool);

impl Canvas {
    pub fn new(width: f64, height: f64) -> Canvas {
        Canvas { width, height, body: String::new(), m00: 1.0, m10: 0.0, m01: 0.0, m11: 1.0, m02: 0.0, m12: 0.0, sheared: false, color: BLACK, font_size: 12 }
    }
    pub fn set_color(&mut self, color: Color) {
        self.color = color;
    }
    pub fn set_font_size(&mut self, size: i32) {
        self.font_size = size;
    }
    pub fn font_size(&self) -> i32 {
        self.font_size
    }
    pub fn transform(&self) -> Transform {
        Transform(self.m00, self.m10, self.m01, self.m11, self.m02, self.m12, self.sheared)
    }
    pub fn set_transform(&mut self, t: Transform) {
        (self.m00, self.m10, self.m01, self.m11, self.m02, self.m12, self.sheared) = (t.0, t.1, t.2, t.3, t.4, t.5, t.6);
    }
    /// AffineTransform.scale in the identity or scale-only state
    pub fn scale(&mut self, sx: f64, sy: f64) {
        assert!(!self.sheared && self.m02 == 0.0 && self.m12 == 0.0);
        self.m00 *= sx;
        self.m11 *= sy;
    }
    /// AffineTransform.rotate(Math.PI / 2), which Java performs as an exact quarter turn
    pub fn rotate_quarter(&mut self) {
        assert!(!self.sheared && self.m02 == 0.0 && self.m12 == 0.0);
        let m0 = self.m00;
        self.m00 = self.m01;
        self.m01 = -m0;
        let m0 = self.m10;
        self.m10 = self.m11;
        self.m11 = -m0;
        self.sheared = true;
    }
    /// AffineTransform.translate in the shear-only state that follows a quarter turn of a scale
    pub fn translate(&mut self, tx: f64, ty: f64) {
        assert!(self.sheared && self.m02 == 0.0 && self.m12 == 0.0);
        self.m02 = ty * self.m01;
        self.m12 = tx * self.m10;
    }
    fn is_identity(&self) -> bool {
        self.m00 == 1.0 && self.m11 == 1.0 && self.m01 == 0.0 && self.m10 == 0.0 && self.m02 == 0.0 && self.m12 == 0.0
    }
    fn push_transform(&mut self) {
        if !self.is_identity() {
            self.body.push_str(" transform='");
            self.push_matrix();
            self.body.push('\'');
        }
    }
    fn push_matrix(&mut self) {
        self.body.push_str("matrix(");
        for (i, v) in [self.m00, self.m10, self.m01, self.m11, self.m02, self.m12].into_iter().enumerate() {
            if i > 0 {
                self.body.push(',');
            }
            self.body.push_str(&java_double(v));
        }
        self.body.push(')');
    }
    fn push_rgb(&mut self) {
        use std::fmt::Write;
        let Color(r, g, b) = self.color;
        let _ = write!(self.body, "rgb({r},{g},{b})");
    }
    /// fill(Rectangle2D)
    pub fn fill_rect(&mut self, x: f64, y: f64, w: f64, h: f64) {
        if w <= 0.0 || h <= 0.0 {
            return;
        }
        self.body.push_str("<rect x='");
        self.body.push_str(&java_double(x));
        self.body.push_str("' y='");
        self.body.push_str(&java_double(y));
        self.body.push_str("' width='");
        self.body.push_str(&java_double(w));
        self.body.push_str("' height='");
        self.body.push_str(&java_double(h));
        self.body.push_str("' style='fill:");
        self.push_rgb();
        self.body.push('\'');
        self.push_transform();
        self.body.push_str("/>");
    }
    /// fill(Path2D) of a closed triangle
    pub fn fill_triangle(&mut self, points: [(f64, f64); 3]) {
        self.body.push_str("<g style='fill:");
        self.push_rgb();
        self.body.push_str(";stroke:none'");
        self.push_transform();
        self.body.push_str("><path d='");
        for (i, (x, y)) in points.into_iter().enumerate() {
            self.body.push(if i == 0 { 'M' } else { 'L' });
            self.body.push_str(&java_double(x));
            self.body.push(',');
            self.body.push_str(&java_double(y));
        }
        self.body.push_str("Z'/></g>");
    }
    /// drawString(String, int, int)
    pub fn draw_string(&mut self, text: &str, x: i32, y: i32) {
        use std::fmt::Write;
        if text.is_empty() {
            return;
        }
        let wrapped = !self.is_identity();
        if wrapped {
            self.body.push_str("<g transform='");
            self.push_matrix();
            self.body.push_str("'>");
        }
        let _ = write!(self.body, "<text x='{x}.0' y='{y}.0' style='fill: ");
        self.push_rgb();
        let _ = write!(self.body, "; fill-opacity: 1.0; font-family: sans-serif; font-size: {}px;'>", self.font_size);
        escape_xml(text, &mut self.body);
        self.body.push_str("</text>");
        if wrapped {
            self.body.push_str("</g>");
        }
    }
    /// getSVGElement
    pub fn svg_element(&self) -> String {
        let mut s = String::with_capacity(self.body.len() + 200);
        s.push_str(
            "<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' \
             xmlns:jfreesvg='http://www.jfree.org/jfreesvg/svg' width='",
        );
        s.push_str(&java_double(self.width));
        s.push_str("' height='");
        s.push_str(&java_double(self.height));
        s.push_str("'>");
        s.push_str(&self.body);
        s.push_str("</svg>");
        s
    }
}

/// SvgUtil.Alignment
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Alignment {
    Left,
    Right,
    Center,
}

/// SvgUtil.drawStringFromCenter
fn draw_string_from_center(canvas: &mut Canvas, text: &str, center_x: f64, center_y: f64) {
    let size = canvas.font_size();
    let (width, height) = string_bounds(size, text);
    let x = (center_x - 0.5 * width) as i32;
    let y = (center_y + 0.5 * height - font_descent(size) as f64) as i32;
    canvas.draw_string(text, x, y);
}

/// SvgRender.renderText
#[allow(clippy::too_many_arguments)]
fn render_text(
    canvas: &mut Canvas, offset: (f64, f64), box_x: f64, box_y: f64, font_size: i32, fg: Color, text: &str, bg: Option<Color>, expand: f64,
    alignment: Alignment,
) {
    let current_font = canvas.font_size();
    canvas.set_font_size(font_size);
    let current_transform = canvas.transform();
    canvas.scale(1.0 / BASE_BOX_SIZE, 1.0 / BASE_BOX_SIZE);
    let mut bounds_width = 0.0;
    if alignment == Alignment::Center {
        let center_x = (box_x + offset.0) * BASE_BOX_SIZE;
        let center_y = (box_y + offset.1) * BASE_BOX_SIZE;
        if let Some(bg) = bg {
            // SvgUtil.scaleRectangleFromCenter(SvgUtil.getStringBoundsFromCenter(...), factor)
            let (width, height) = string_bounds(font_size, text);
            let (rx, ry) = (center_x - 0.5 * width, center_y - 0.5 * height);
            let (scaled_w, scaled_h) = (expand * width, expand * height);
            let (rcx, rcy) = (rx + width / 2.0, ry + height / 2.0);
            canvas.set_color(bg);
            canvas.fill_rect(rcx - 0.5 * scaled_w, rcy - 0.5 * scaled_h, scaled_w, scaled_h);
        }
        canvas.set_color(fg);
        draw_string_from_center(canvas, text, center_x, center_y);
    } else {
        bounds_width = string_bounds(font_size, text).0;
    }
    canvas.set_transform(current_transform);
    canvas.set_font_size(current_font);
    // Rectangle2D.getCenterX of bounds that start at x = 0
    let bounds_center_x = 0.0 + bounds_width / 2.0;
    match alignment {
        Alignment::Left => {
            render_text(canvas, offset, box_x + bounds_center_x / BASE_BOX_SIZE, box_y, font_size, fg, text, bg, expand, Alignment::Center)
        }
        Alignment::Right => {
            render_text(canvas, offset, box_x - bounds_center_x / BASE_BOX_SIZE, box_y, font_size, fg, text, bg, expand, Alignment::Center)
        }
        Alignment::Center => {}
    }
}

fn draw_top_box_border(canvas: &mut Canvas, o: (f64, f64), box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + o.0, o.1, 1.0, prop);
}
fn draw_right_box_border(canvas: &mut Canvas, o: (f64, f64), box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + 1.0 - prop + o.0, o.1, prop, 1.0);
}
fn draw_bottom_box_border(canvas: &mut Canvas, o: (f64, f64), box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + o.0, 1.0 - prop + o.1, 1.0, prop);
}
fn draw_left_box_border(canvas: &mut Canvas, o: (f64, f64), box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + o.0, o.1, prop, 1.0);
}
fn draw_right_insert_indicator(canvas: &mut Canvas, o: (f64, f64), box_idx: f64) {
    draw_right_box_border(canvas, o, box_idx, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + 1.0 - 1.0 / 3.0 + o.0, o.1, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + 1.0 - 1.0 / 3.0 + o.0, 1.0 - BOUNDARY_BOX_PROPORTION + o.1, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
}
fn draw_left_insert_indicator(canvas: &mut Canvas, o: (f64, f64), box_idx: f64) {
    draw_left_box_border(canvas, o, box_idx, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + o.0, o.1, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + o.0, 1.0 - BOUNDARY_BOX_PROPORTION + o.1, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
}
fn draw_del(canvas: &mut Canvas, o: (f64, f64), label_del_len: Option<i32>, box_idx_start: i32, box_idx_end: i32) {
    let del_len = box_idx_end - box_idx_start + 1;
    canvas.set_color(BLACK);
    canvas.fill_rect(box_idx_start as f64 + o.0, 0.5 - 0.5 * DEL_CONNECTOR_BOX_PROPORTION + o.1, del_len as f64, DEL_CONNECTOR_BOX_PROPORTION);
    let Some(label) = label_del_len else { return };
    render_text(
        canvas,
        o,
        (box_idx_end + 1) as f64 - 0.5 * del_len as f64,
        0.5,
        INDEL_FONT_SIZE,
        BLACK,
        &label.to_string(),
        Some(WHITE),
        1.5,
        Alignment::Center,
    );
}
fn draw_forward_arrow(canvas: &mut Canvas, o: (f64, f64), left: f64, top: f64, width: f64, height: f64) {
    canvas.fill_triangle([(left + o.0, top + o.1), (left + width + o.0, top + 0.5 * height + o.1), (left + o.0, top + height + o.1)]);
}
fn draw_reverse_arrow(canvas: &mut Canvas, o: (f64, f64), left: f64, top: f64, width: f64, height: f64) {
    canvas.fill_triangle([(left + o.0, top + 0.5 * height + o.1), (left + width + o.0, top + o.1), (left + width + o.0, top + height + o.1)]);
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BaseKind {
    Char(u8),
    Missing,
    Del,
}

/// BaseViewModel. The inserted bases and their qualities, which Java keeps beside the count, are not
/// drawn and are not kept.
#[derive(Clone, Copy, Debug)]
pub struct BaseViewModel {
    pub kind: BaseKind,
    pub base_q: i32,
    pub soft_clip: bool,
    pub overlapped: bool,
    pub right_insert_count: i32,
}

impl BaseViewModel {
    pub fn base(base: u8) -> Self {
        BaseViewModel { kind: BaseKind::Char(base), base_q: 0, soft_clip: false, overlapped: false, right_insert_count: 0 }
    }
    pub fn with_qual(base: u8, base_q: i32, soft_clip: bool, overlapped: bool) -> Self {
        BaseViewModel { kind: BaseKind::Char(base), base_q, soft_clip, overlapped, right_insert_count: 0 }
    }
    pub fn missing() -> Self {
        BaseViewModel { kind: BaseKind::Missing, base_q: 0, soft_clip: false, overlapped: false, right_insert_count: 0 }
    }
    pub fn del(overlapped: bool) -> Self {
        BaseViewModel { kind: BaseKind::Del, base_q: 0, soft_clip: false, overlapped, right_insert_count: 0 }
    }
    pub fn is_missing(&self) -> bool {
        self.kind == BaseKind::Missing
    }
    pub fn is_del(&self) -> bool {
        self.kind == BaseKind::Del
    }
    pub fn char_base(&self) -> Option<u8> {
        match self.kind {
            BaseKind::Char(c) => Some(c),
            _ => None,
        }
    }
}

/// htsjdk CigarOperator, as far as the view models distinguish them
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VisCigarOp {
    /// M, EQ, X
    Match,
    S,
    I,
    D,
    /// H, N
    Skip,
    /// P: no case in Java's switch
    Pad,
}

/// BaseSeqViewModel
#[derive(Clone, Debug)]
pub struct BaseSeqViewModel {
    bases: Vec<BaseViewModel>,
    pos_start: i32,
    pub left_is_forward_strand: Option<bool>,
    pub right_is_forward_strand: Option<bool>,
    pub first_base_pos: i32,
    pub last_base_pos: i32,
}

impl BaseSeqViewModel {
    pub fn new(bases: Vec<BaseViewModel>, pos_start: i32, left: Option<bool>, right: Option<bool>) -> Self {
        let mut first_base_pos = -1;
        let mut last_base_pos = -1;
        for (i, b) in bases.iter().enumerate() {
            if b.is_missing() {
                continue;
            }
            if first_base_pos == -1 {
                first_base_pos = i as i32 + pos_start;
            }
            last_base_pos = i as i32 + pos_start;
        }
        BaseSeqViewModel { bases, pos_start, left_is_forward_strand: left, right_is_forward_strand: right, first_base_pos, last_base_pos }
    }

    /// fromStr
    pub fn from_bases(bases: &[u8], pos_start: i32) -> Self {
        Self::new(bases.iter().map(|b| BaseViewModel::base(*b)).collect(), pos_start, None, None)
    }

    /// create(unclippedStart, cigarElements, bases, baseQuals, readNegativeStrandFlag). The qualities are
    /// Java bytes (signed). An index beyond the bases or the qualities is Java's
    /// ArrayIndexOutOfBoundsException.
    pub fn create(unclipped_start: i32, cigar: &[(i32, VisCigarOp)], bases: &[u8], base_quals: &[u8], read_negative_strand: bool) -> VisResult<Self> {
        let mut indexed: Vec<BaseViewModel> = Vec::new();
        let mut base_idx: usize = 0;
        let at = |idx: usize| -> VisResult<(u8, i32)> {
            match (bases.get(idx), base_quals.get(idx)) {
                (Some(b), Some(q)) => Ok((*b, *q as i8 as i32)),
                _ => Err(VisError::new(format!("java.lang.ArrayIndexOutOfBoundsException: Index {idx} out of bounds for length {}", bases.len().min(base_quals.len())))),
            }
        };
        for (len, op) in cigar {
            for _ in 0..*len {
                match op {
                    VisCigarOp::Match => {
                        let (b, q) = at(base_idx)?;
                        indexed.push(BaseViewModel::with_qual(b, q, false, false));
                        base_idx += 1;
                    }
                    VisCigarOp::S => {
                        let (b, q) = at(base_idx)?;
                        indexed.push(BaseViewModel::with_qual(b, q, true, false));
                        base_idx += 1;
                    }
                    VisCigarOp::I => {
                        if !indexed.is_empty() {
                            at(base_idx)?;
                            indexed.last_mut().unwrap().right_insert_count += 1;
                        }
                        base_idx += 1;
                    }
                    VisCigarOp::D => indexed.push(BaseViewModel::del(false)),
                    VisCigarOp::Skip => indexed.push(BaseViewModel::missing()),
                    VisCigarOp::Pad => {}
                }
            }
        }
        Ok(Self::new(indexed, unclipped_start, Some(!read_negative_strand), Some(!read_negative_strand)))
    }

    /// fromStringWithCigar; an index beyond the string is Java's StringIndexOutOfBoundsException
    pub fn from_string_with_cigar(base_str: &[u8], cigar: &[(i32, VisCigarOp)], pos_start: i32) -> VisResult<Self> {
        let mut indexed: Vec<BaseViewModel> = Vec::new();
        let mut base_idx: usize = 0;
        let at = |idx: usize| -> VisResult<u8> {
            base_str.get(idx).copied().ok_or_else(|| {
                VisError::new(format!("java.lang.StringIndexOutOfBoundsException: Index {idx} out of bounds for length {}", base_str.len()))
            })
        };
        for (len, op) in cigar {
            for _ in 0..*len {
                match op {
                    VisCigarOp::Match => {
                        indexed.push(BaseViewModel::base(at(base_idx)?));
                        base_idx += 1;
                    }
                    VisCigarOp::S => {
                        indexed.push(BaseViewModel::with_qual(at(base_idx)?, 0, true, false));
                        base_idx += 1;
                    }
                    VisCigarOp::I => {
                        if !indexed.is_empty() {
                            at(base_idx)?;
                            indexed.last_mut().unwrap().right_insert_count += 1;
                        }
                        base_idx += 1;
                    }
                    VisCigarOp::D => indexed.push(BaseViewModel::del(false)),
                    VisCigarOp::Skip => indexed.push(BaseViewModel::missing()),
                    VisCigarOp::Pad => {}
                }
            }
        }
        Ok(Self::new(indexed, pos_start, None, None))
    }

    pub fn get_base(&self, pos: i32) -> BaseViewModel {
        let idx = pos - self.pos_start;
        if idx < 0 || idx >= self.bases.len() as i32 {
            return BaseViewModel::missing();
        }
        self.bases[idx as usize]
    }

    pub fn has_orientation(&self) -> bool {
        self.left_is_forward_strand.is_some() && self.right_is_forward_strand.is_some()
    }

    pub fn mismatch_count(&self, ref_bases: &BaseSeqViewModel) -> i32 {
        let mut mismatches = 0;
        for i in ref_bases.first_base_pos..=ref_bases.last_base_pos {
            let ref_base = ref_bases.get_base(i);
            let base = self.get_base(i);
            if let (false, Some(b), Some(r)) = (base.soft_clip, base.char_base(), ref_base.char_base()) {
                if r != b {
                    mismatches += 1;
                }
            }
        }
        mismatches
    }

    pub fn clear_soft_clips(&mut self) {
        for b in &mut self.bases {
            b.soft_clip = false;
        }
    }
}

fn base_colours(base: u8) -> Option<(Color, Color)> {
    match base {
        b'G' => Some((Color(255, 232, 145), Color(151, 83, 9))),
        b'C' => Some((Color(212, 221, 240), Color(11, 105, 186))),
        b'A' => Some((Color(191, 236, 199), Color(15, 95, 16))),
        b'T' => Some((Color(255, 192, 199), Color(145, 0, 7))),
        _ => None,
    }
}

/// SvgRender.renderBaseSeq(svgCanvas, boxOffset, baseBoxSizePx, renderRegion, bases, shadeQuals,
/// posToBordersMap, refBases, renderLeftOrientationMarker, renderRightOrientationMarker,
/// reverseOrientationMarkers), for an empty border map (the assembly visualiser never passes borders).
#[allow(clippy::too_many_arguments)]
pub fn render_base_seq(
    canvas: &mut Canvas, o: (f64, f64), base_box_size_px: f64, render_region: BaseRegion, bases: &BaseSeqViewModel, shade_quals: bool,
    ref_bases: Option<&BaseSeqViewModel>, render_left_marker: bool, render_right_marker: bool, reverse_markers: bool,
) {
    canvas.scale(base_box_size_px, base_box_size_px);
    canvas.set_font_size(BASE_FONT_SIZE);

    let mut first_base_idx: Option<i32> = None;
    let mut last_base_idx: Option<i32> = None;
    let mut first_overlapping: Option<i32> = None;
    let mut last_overlapping: Option<i32> = None;
    let mut prev_base = bases.get_base(render_region.start - 1);
    let mut del_len = 0;
    for i in render_region.start..=render_region.end {
        let box_idx = i - render_region.start + BOX_PADDING;
        let box_f = box_idx as f64;
        let ref_base = ref_bases.map_or_else(BaseViewModel::missing, |r| r.get_base(i));
        let base = bases.get_base(i);
        if !base.is_missing() {
            if first_base_idx.is_none() {
                first_base_idx = Some(box_idx);
            }
            last_base_idx = Some(box_idx);
        }
        if base.is_del() {
            del_len += 1;
        } else if let BaseKind::Char(ch) = base.kind {
            let matches_ref = !base.soft_clip && ref_base.kind == BaseKind::Char(ch);
            let mut bg = WHITE;
            let mut fg = BLACK;
            if matches_ref {
                bg = LIGHT_GRAY;
            } else if let Some((b, f)) = base_colours(ch) {
                bg = b;
                fg = f;
            }
            if shade_quals {
                let lighten_factor = 1.0 - 1.0 * base.base_q as f64 / MAX_BASEQ_SHADING_CUTTOFF as f64;
                fg = lighten(fg, lighten_factor);
                bg = lighten(bg, lighten_factor);
            }
            canvas.set_color(bg);
            canvas.fill_rect(box_f + o.0, o.1, 1.0, 1.0);
            if !matches_ref {
                let mut buf = [0u8; 4];
                let text: &str = (ch as char).encode_utf8(&mut buf);
                let size = canvas.font_size();
                render_text(canvas, o, box_f + 0.5, 0.5, size, fg, text, None, 1.0, Alignment::Center);
            }
        }
        if base.overlapped {
            if first_overlapping.is_none() {
                first_overlapping = Some(box_idx);
            }
            last_overlapping = Some(box_idx);
            canvas.set_color(overlapping_fragment_border_color());
            draw_top_box_border(canvas, o, box_f, BOUNDARY_BOX_PROPORTION);
            draw_bottom_box_border(canvas, o, box_f, BOUNDARY_BOX_PROPORTION);
        }
        if !base.is_del() && del_len > 0 {
            draw_del(canvas, o, Some(del_len), box_idx - del_len, box_idx - 1);
            del_len = 0;
        }
        canvas.set_color(INSERT_COLOR);
        if base.right_insert_count > 0 {
            draw_right_insert_indicator(canvas, o, box_f);
        }
        if prev_base.right_insert_count > 0 {
            draw_left_insert_indicator(canvas, o, box_f);
            render_text(canvas, o, box_f, 0.5, INDEL_FONT_SIZE, WHITE, &prev_base.right_insert_count.to_string(), Some(INSERT_COLOR), 1.0, Alignment::Center);
        }
        prev_base = base;
    }
    if del_len > 0 {
        let box_end_idx = render_region.end - render_region.start + BOX_PADDING;
        draw_del(canvas, o, Some(del_len), box_end_idx + 1 - del_len, box_end_idx);
    }
    if let (Some(first), Some(last)) = (first_overlapping, last_overlapping) {
        canvas.set_color(overlapping_fragment_border_color());
        draw_left_box_border(canvas, o, first as f64, BOUNDARY_BOX_PROPORTION);
        draw_right_box_border(canvas, o, last as f64, BOUNDARY_BOX_PROPORTION);
    }
    let (Some(left_forward), Some(right_forward), Some(first), Some(last)) =
        (bases.left_is_forward_strand, bases.right_is_forward_strand, first_base_idx, last_base_idx)
    else {
        return;
    };
    if render_left_marker {
        if left_forward != reverse_markers {
            canvas.set_color(FORWARD_STRAND_COLOR);
            draw_right_box_border(canvas, o, (first - 1) as f64, 0.5);
        } else {
            canvas.set_color(REVERSE_STRAND_COLOR);
            draw_reverse_arrow(canvas, o, first as f64 - 0.5, 0.0, 0.5, 1.0);
        }
    }
    if render_right_marker {
        if right_forward != reverse_markers {
            canvas.set_color(FORWARD_STRAND_COLOR);
            draw_forward_arrow(canvas, o, (last + 1) as f64, 0.0, 0.5, 1.0);
        } else {
            canvas.set_color(REVERSE_STRAND_COLOR);
            draw_left_box_border(canvas, o, (last + 1) as f64, 0.5);
        }
    }
}

/// String.format(Locale.US, "%,d", n)
pub fn grouped(n: i32) -> String {
    let digits = (n as i64).abs().to_string();
    let mut out = String::new();
    if n < 0 {
        out.push('-');
    }
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// SvgRender.renderCoords(svgCanvas, boxOffset, canvasSize, baseBoxSizePx, renderRegion, centerPosition,
/// displayEveryNthCoord, maxStringWidth, reverse)
#[allow(clippy::too_many_arguments)]
pub fn render_coords(
    canvas: &mut Canvas, o: (f64, f64), canvas_size: (f64, f64), base_box_size_px: f64, render_region: BaseRegion, center_position: i32,
    display_every_nth: i32, max_string_width: f64, reverse: bool,
) {
    let scaling_factor = base_box_size_px / BASE_BOX_SIZE;
    let char_width = string_bounds(COORD_FONT_SIZE, "9").0;
    let widths: Vec<f64> = (render_region.start..=render_region.end).map(|i| string_bounds(COORD_FONT_SIZE, &grouped(i)).0).collect();
    canvas.set_font_size(COORD_FONT_SIZE);
    canvas.set_color(LIGHT_GRAY);
    canvas.fill_rect(base_box_size_px * o.0, base_box_size_px * o.1, canvas_size.0 + base_box_size_px * o.0, canvas_size.1 + base_box_size_px * o.1);
    canvas.scale(scaling_factor, scaling_factor);
    let scaling_transform = canvas.transform();
    canvas.set_color(BLACK);
    for b in 0..render_region.base_length() {
        canvas.set_transform(scaling_transform);
        let box_idx = (b + BOX_PADDING) as f64;
        let base_pos = if reverse { render_region.end - b } else { render_region.start + b };
        let base_pos_offset = base_pos - center_position;
        let x = BASE_BOX_SIZE * (box_idx + 0.5 - 0.25 * BOUNDARY_BOX_PROPORTION + o.0);
        if base_pos_offset % display_every_nth != 0 {
            canvas.fill_rect(
                x,
                canvas_size.1 / scaling_factor - 0.75 * char_width + BASE_BOX_SIZE * o.1,
                BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
                0.75 * char_width,
            );
            continue;
        }
        canvas.fill_rect(x, canvas_size.1 / scaling_factor - 1.5 * char_width + BASE_BOX_SIZE * o.1, BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION, 1.5 * char_width);
        let string_idx = if reverse { widths.len() as i32 - 1 - b } else { b };
        canvas.rotate_quarter();
        canvas.translate(max_string_width - 0.5 * widths[string_idx as usize] + BASE_BOX_SIZE * o.1, -BASE_BOX_SIZE * (box_idx + 0.5 + o.0));
        draw_string_from_center(canvas, &grouped(base_pos), 0.0, 0.0);
    }
}

/// SvgRender.ChrLabel; a label without a chromosome marks inserted bases
#[derive(Clone, Debug)]
pub struct ChrLabel {
    pub chromosome: Option<String>,
    pub position: i32,
    pub view_region: BaseRegion,
    pub ref_region: Option<BaseRegion>,
    pub is_reversed: bool,
    pub alignments: Vec<Alignment>,
}

/// SvgRender.renderChrLabels
pub fn render_chr_labels(base_box_size_px: f64, labels: &[ChrLabel]) -> VisResult<Canvas> {
    let (Some(first), Some(last)) = (labels.first(), labels.last()) else {
        return Err(VisError::new("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"));
    };
    let start_idx = first.view_region.start;
    let end_idx = last.view_region.end;
    let box_length = end_idx - start_idx + 1;
    let mut canvas = Canvas::new(base_box_size_px * box_length as f64, base_box_size_px);
    canvas.scale(base_box_size_px, base_box_size_px);
    canvas.set_font_size(BASE_FONT_SIZE);
    for label in labels {
        let view_region = label.view_region;
        let Some(chromosome) = &label.chromosome else {
            let box_x = view_region.start as f64 + 0.5 * view_region.base_length() as f64;
            let size = canvas.font_size();
            render_text(&mut canvas, (0.0, 0.0), box_x, 0.5, size, BLACK, INSERT_LABEL, None, 1.0, Alignment::Center);
            continue;
        };
        // RefGenomeFunctions.enforceChrPrefix
        let chromosome = if chromosome.starts_with("chr") { chromosome.clone() } else { format!("chr{chromosome}") };
        for alignment in &label.alignments {
            let ref_region = label.ref_region.ok_or_else(|| VisError::new("java.lang.NullPointerException"))?;
            let strand_label = if label.is_reversed { "<<" } else { ">>" };
            let left_position = if label.is_reversed { ref_region.end } else { ref_region.start };
            let right_position = if label.is_reversed { ref_region.start } else { ref_region.end };
            let position = if *alignment == Alignment::Left { left_position } else { right_position };
            let chr_position_label = format!("{chromosome}:{position}");
            let (full_label, box_x) = if *alignment == Alignment::Left {
                (format!("{chr_position_label} {strand_label}"), (view_region.start + BOX_PADDING) as f64)
            } else {
                (format!("{strand_label} {chr_position_label}"), (view_region.start + view_region.base_length()) as f64)
            };
            let size = canvas.font_size();
            render_text(&mut canvas, (0.0, 0.0), box_x, 0.5, size, BLACK, &full_label, None, 1.0, *alignment);
        }
    }
    Ok(canvas)
}
