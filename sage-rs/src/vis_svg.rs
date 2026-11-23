//! The visualiser's drawing layer: a writer that serialises shapes as JFreeSVG 5.0.5's
//! `SVGGraphics2D.getSVGElement()` does, the font metrics Java's AWT reported for the "SansSerif"
//! logical font, and the ports of `hmf-common` `vis/SvgRender`, `SvgUtil`, `BaseSeqViewModel`,
//! `BaseViewModel`, `GeneRegionViewModel` and `ColorUtil`.

use crate::constants::jround;
use crate::types::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color(pub u8, pub u8, pub u8);

pub const BLACK: Color = Color(0, 0, 0);
pub const WHITE: Color = Color(255, 255, 255);
pub const RED: Color = Color(255, 0, 0);
pub const GREEN: Color = Color(0, 255, 0);
pub const BLUE: Color = Color(0, 0, 255);
pub const LIGHT_GRAY: Color = Color(192, 192, 192);
pub const DARK_BLUE: Color = Color(0, 0, 139);
pub const DARK_GREEN: Color = Color(0, 139, 0);
pub const PURPLE: Color = Color(128, 0, 128);
pub const FORWARD_STRAND_COLOR: Color = Color(81, 144, 207);
pub const REVERSE_STRAND_COLOR: Color = Color(242, 167, 121);
const INSERT_COLOR: Color = Color(153, 50, 204);
const AMINO_ACID_MATCH_BG: [Color; 2] = [Color(92, 92, 164), Color(12, 12, 120)];
const AMINO_ACID_MISMATCH_BG: [Color; 2] = [Color(179, 130, 77), Color(182, 103, 18)];

pub const BOX_PADDING: i32 = 1;
pub const BASE_BOX_SIZE: f64 = 30.0;
const DEL_CONNECTOR_BOX_PROPORTION: f64 = 0.2;
const BOUNDARY_BOX_PROPORTION: f64 = 0.2;
const INDEL_FONT_SIZE: i32 = 18;
const BASE_FONT_SIZE: i32 = 26;
const COORD_FONT_SIZE: i32 = 30;
const MAX_BASEQ_SHADING_CUTTOFF: i32 = 37;

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
/// Characters outside printable ASCII are given the width of a space.
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
    stroke_width: f32,
}

#[derive(Clone, Copy)]
pub struct Transform(f64, f64, f64, f64, f64, f64, bool);

impl Canvas {
    pub fn new(width: f64, height: f64) -> Canvas {
        Canvas {
            width,
            height,
            body: String::new(),
            m00: 1.0,
            m10: 0.0,
            m01: 0.0,
            m11: 1.0,
            m02: 0.0,
            m12: 0.0,
            sheared: false,
            color: BLACK,
            font_size: 12,
            stroke_width: 1.0,
        }
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
        let Color(r, g, b) = self.color;
        self.body.push_str(&format!("rgb({r},{g},{b})"));
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
    pub fn set_stroke_width(&mut self, width: f32) {
        self.stroke_width = width;
    }
    pub fn stroke_width(&self) -> f32 {
        self.stroke_width
    }
    /// draw(Line2D) with a BasicStroke of the current width (square cap, mitre join, limit 10)
    pub fn draw_line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64) {
        self.body.push_str("<line x1='");
        self.body.push_str(&java_double(x1));
        self.body.push_str("' y1='");
        self.body.push_str(&java_double(y1));
        self.body.push_str("' x2='");
        self.body.push_str(&java_double(x2));
        self.body.push_str("' y2='");
        self.body.push_str(&java_double(y2));
        self.body.push_str("' style='stroke-width:");
        self.body.push_str(&java_double(self.stroke_width as f64));
        self.body.push_str(";stroke:");
        self.push_rgb();
        self.body.push_str(";stroke-opacity:1.0;stroke-linecap:square;stroke-miterlimit:10.0'");
        self.push_transform();
        self.body.push_str("/>");
    }
    /// drawString(String, int, int)
    pub fn draw_string(&mut self, text: &str, x: i32, y: i32) {
        if text.is_empty() {
            return;
        }
        let wrapped = !self.is_identity();
        if wrapped {
            self.body.push_str("<g transform='");
            self.push_matrix();
            self.body.push_str("'>");
        }
        self.body.push_str(&format!("<text x='{x}.0' y='{y}.0' style='fill: "));
        self.push_rgb();
        self.body.push_str(&format!("; fill-opacity: 1.0; font-family: sans-serif; font-size: {}px;'>", self.font_size));
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

/// SvgUtil.drawStringFromCenter
fn draw_string_from_center(canvas: &mut Canvas, text: &str, center_x: f64, center_y: f64) {
    let size = canvas.font_size();
    let (width, height) = string_bounds(size, text);
    let x = (center_x - 0.5 * width) as i32;
    let y = (center_y + 0.5 * height - font_descent(size) as f64) as i32;
    canvas.draw_string(text, x, y);
}

/// SvgRender.renderText with CENTER alignment (the only one SAGE uses); the box offset is always zero
fn render_text(
    canvas: &mut Canvas, box_x: f64, box_y: f64, font_size: i32, fg: Color, text: &str, bg: Option<Color>, expand: f64,
) {
    let current_font = canvas.font_size();
    canvas.set_font_size(font_size);
    let current_transform = canvas.transform();
    canvas.scale(1.0 / BASE_BOX_SIZE, 1.0 / BASE_BOX_SIZE);
    let center_x = (box_x + 0.0) * BASE_BOX_SIZE;
    let center_y = (box_y + 0.0) * BASE_BOX_SIZE;
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
    canvas.set_transform(current_transform);
    canvas.set_font_size(current_font);
}

fn draw_top_box_border(canvas: &mut Canvas, box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + 0.0, 0.0, 1.0, prop);
}
fn draw_right_box_border(canvas: &mut Canvas, box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + 1.0 - prop + 0.0, 0.0, prop, 1.0);
}
fn draw_bottom_box_border(canvas: &mut Canvas, box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + 0.0, 1.0 - prop + 0.0, 1.0, prop);
}
fn draw_left_box_border(canvas: &mut Canvas, box_idx: f64, prop: f64) {
    canvas.fill_rect(box_idx + 0.0, 0.0, prop, 1.0);
}
fn draw_right_insert_indicator(canvas: &mut Canvas, box_idx: f64) {
    draw_right_box_border(canvas, box_idx, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + 1.0 - 1.0 / 3.0 + 0.0, 0.0, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + 1.0 - 1.0 / 3.0 + 0.0, 1.0 - BOUNDARY_BOX_PROPORTION + 0.0, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
}
fn draw_left_insert_indicator(canvas: &mut Canvas, box_idx: f64) {
    draw_left_box_border(canvas, box_idx, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + 0.0, 0.0, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
    canvas.fill_rect(box_idx + 0.0, 1.0 - BOUNDARY_BOX_PROPORTION + 0.0, 1.0 / 3.0, BOUNDARY_BOX_PROPORTION);
}
fn draw_del(canvas: &mut Canvas, label_del_len: Option<i32>, box_idx_start: i32, box_idx_end: i32) {
    let del_len = box_idx_end - box_idx_start + 1;
    canvas.set_color(BLACK);
    canvas.fill_rect(
        box_idx_start as f64 + 0.0,
        0.5 - 0.5 * DEL_CONNECTOR_BOX_PROPORTION + 0.0,
        del_len as f64,
        DEL_CONNECTOR_BOX_PROPORTION,
    );
    let Some(label) = label_del_len else { return };
    render_text(
        canvas,
        (box_idx_end + 1) as f64 - 0.5 * del_len as f64,
        0.5,
        INDEL_FONT_SIZE,
        BLACK,
        &label.to_string(),
        Some(WHITE),
        1.5,
    );
}
fn draw_forward_arrow(canvas: &mut Canvas, left: f64, top: f64, width: f64, height: f64) {
    canvas.fill_triangle([
        (left + 0.0, top + 0.0),
        (left + width + 0.0, top + 0.5 * height + 0.0),
        (left + 0.0, top + height + 0.0),
    ]);
}
fn draw_reverse_arrow(canvas: &mut Canvas, left: f64, top: f64, width: f64, height: f64) {
    canvas.fill_triangle([
        (left + 0.0, top + 0.5 * height + 0.0),
        (left + width + 0.0, top + 0.0),
        (left + width + 0.0, top + height + 0.0),
    ]);
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BaseKind {
    Char(u8),
    Missing,
    Del,
}

/// BaseViewModel
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
}

/// BaseSeqViewModel
#[derive(Clone, Debug)]
pub struct BaseSeqViewModel {
    pub bases: Vec<BaseViewModel>,
    pub pos_start: i32,
    pub left_is_forward_strand: Option<bool>,
    pub right_is_forward_strand: Option<bool>,
    pub first_base_pos: i32,
    pub last_base_pos: i32,
}

/// SAMRecord.getUnclippedStart: soft and hard clips both count
pub fn unclipped_start(read: &Read) -> i32 {
    let mut start = read.start;
    for el in &read.cigar {
        if matches!(el.op, Op::S | Op::H) {
            start -= el.len;
        } else {
            break;
        }
    }
    start
}

/// SAMRecord.getUnclippedEnd
pub fn unclipped_end(read: &Read) -> i32 {
    let mut end = read.end;
    for el in read.cigar.iter().rev() {
        if matches!(el.op, Op::S | Op::H) {
            end += el.len;
        } else {
            break;
        }
    }
    end
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

    /// fromRead
    pub fn from_read(read: &Read) -> Self {
        Self::from_consensus_fragment(read, None, None)
    }

    /// fromConsensusFragment / create.
    pub fn from_consensus_fragment(read: &Read, first: Option<&BaseSeqViewModel>, second: Option<&BaseSeqViewModel>) -> Self {
        // S44, Java: the mates' bases are looked up at `baseIdx + unclippedStart`, where baseIdx counts
        // read bases and not reference positions, so after an insert or a delete in the merged read
        // the overlap border is drawn against the wrong columns. Corrected: the column of the base.
        let java_overlap_index = crate::compat::java("S44");
        let unclipped_start = unclipped_start(read);
        let mut indexed: Vec<BaseViewModel> = Vec::new();
        let mut base_idx: usize = 0;
        for el in &read.cigar {
            for _ in 0..el.len {
                // a record whose cigar reads more bases than it holds is drawn with N at quality 0
                let base = read.bases.get(base_idx).copied().unwrap_or(b'N');
                let qual = read.quals.get(base_idx).map_or(0, |q| *q as i8 as i32);
                let pos = if java_overlap_index { base_idx as i32 } else { indexed.len() as i32 } + unclipped_start;
                let overlapped = match (first, second) {
                    (Some(f), Some(s)) => !f.get_base(pos).is_missing() && !s.get_base(pos).is_missing(),
                    _ => false,
                };
                match el.op {
                    Op::M | Op::Eq | Op::X => {
                        indexed.push(BaseViewModel::with_qual(base, qual, false, overlapped));
                        base_idx += 1;
                    }
                    Op::S => {
                        indexed.push(BaseViewModel::with_qual(base, qual, true, overlapped));
                        base_idx += 1;
                    }
                    Op::I => {
                        if let Some(last) = indexed.last_mut() {
                            last.right_insert_count += 1;
                        }
                        base_idx += 1;
                    }
                    Op::D => indexed.push(BaseViewModel::del(overlapped)),
                    Op::H | Op::N => indexed.push(BaseViewModel::missing()),
                    Op::P => {}
                }
            }
        }
        let forward = !read.negative_strand();
        match (first, second) {
            (Some(f), Some(s)) => {
                let left = if f.first_base_pos <= s.first_base_pos { f.left_is_forward_strand } else { s.left_is_forward_strand };
                let right = if f.last_base_pos >= s.last_base_pos { f.right_is_forward_strand } else { s.right_is_forward_strand };
                Self::new(indexed, unclipped_start, left, right)
            }
            _ => Self::new(indexed, unclipped_start, Some(forward), Some(forward)),
        }
    }

    pub fn get_base(&self, pos: i32) -> BaseViewModel {
        let idx = pos - self.pos_start;
        if idx < 0 || idx >= self.bases.len() as i32 {
            return BaseViewModel::missing();
        }
        self.bases[idx as usize]
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BorderLocation {
    Right,
    Left,
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

/// SvgRender.renderBaseSeq(baseBoxSizePx, renderRegion, bases, shadeQuals, posToBordersMap, refBases)
pub fn render_base_seq(
    base_box_size_px: f64, render_region: BaseRegion, bases: &BaseSeqViewModel, shade_quals: bool,
    borders: &[(i32, Vec<(BorderLocation, Color)>)], ref_bases: Option<&BaseSeqViewModel>,
) -> Canvas {
    let mut canvas =
        Canvas::new(base_box_size_px * (render_region.base_length() + 2 * BOX_PADDING) as f64, base_box_size_px);
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
            canvas.fill_rect(box_f + 0.0, 0.0, 1.0, 1.0);
            if !matches_ref {
                let mut buf = [0u8; 4];
                let text: &str = (ch as char).encode_utf8(&mut buf);
                let size = canvas.font_size();
                render_text(&mut canvas, box_f + 0.5, 0.5, size, fg, text, None, 1.0);
            }
        }
        if base.overlapped {
            if first_overlapping.is_none() {
                first_overlapping = Some(box_idx);
            }
            last_overlapping = Some(box_idx);
            canvas.set_color(overlapping_fragment_border_color());
            draw_top_box_border(&mut canvas, box_f, BOUNDARY_BOX_PROPORTION);
            draw_bottom_box_border(&mut canvas, box_f, BOUNDARY_BOX_PROPORTION);
        }
        if !base.is_del() && del_len > 0 {
            draw_del(&mut canvas, Some(del_len), box_idx - del_len, box_idx - 1);
            del_len = 0;
        }
        if let Some((_, list)) = borders.iter().find(|(pos, _)| *pos == i) {
            for (location, colour) in list {
                canvas.set_color(*colour);
                match location {
                    BorderLocation::Right => draw_right_box_border(&mut canvas, box_f, 0.5 * BOUNDARY_BOX_PROPORTION),
                    BorderLocation::Left => draw_left_box_border(&mut canvas, box_f, 0.5 * BOUNDARY_BOX_PROPORTION),
                }
            }
        }
        canvas.set_color(INSERT_COLOR);
        if base.right_insert_count > 0 {
            draw_right_insert_indicator(&mut canvas, box_f);
        }
        if prev_base.right_insert_count > 0 {
            draw_left_insert_indicator(&mut canvas, box_f);
            render_text(
                &mut canvas,
                box_f,
                0.5,
                INDEL_FONT_SIZE,
                WHITE,
                &prev_base.right_insert_count.to_string(),
                Some(INSERT_COLOR),
                1.0,
            );
        }
        prev_base = base;
    }
    if del_len > 0 {
        let box_end_idx = render_region.end - render_region.start + BOX_PADDING;
        draw_del(&mut canvas, Some(del_len), box_end_idx + 1 - del_len, box_end_idx);
    }
    if let (Some(first), Some(last)) = (first_overlapping, last_overlapping) {
        canvas.set_color(overlapping_fragment_border_color());
        draw_left_box_border(&mut canvas, first as f64, BOUNDARY_BOX_PROPORTION);
        draw_right_box_border(&mut canvas, last as f64, BOUNDARY_BOX_PROPORTION);
    }
    let (Some(left_forward), Some(right_forward), Some(first), Some(last)) =
        (bases.left_is_forward_strand, bases.right_is_forward_strand, first_base_idx, last_base_idx)
    else {
        return canvas;
    };
    if left_forward {
        canvas.set_color(FORWARD_STRAND_COLOR);
        draw_right_box_border(&mut canvas, (first - 1) as f64, 0.5);
    } else {
        canvas.set_color(REVERSE_STRAND_COLOR);
        draw_reverse_arrow(&mut canvas, first as f64 - 0.5, 0.0, 0.5, 1.0);
    }
    if right_forward {
        canvas.set_color(FORWARD_STRAND_COLOR);
        draw_forward_arrow(&mut canvas, (last + 1) as f64, 0.0, 0.5, 1.0);
    } else {
        canvas.set_color(REVERSE_STRAND_COLOR);
        draw_left_box_border(&mut canvas, (last + 1) as f64, 0.5);
    }
    canvas
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

/// SvgRender.renderCoords(baseBoxSizePx, renderRegion, centerPosition, displayEveryNthCoord)
pub fn render_coords(base_box_size_px: f64, render_region: BaseRegion, center_position: i32, display_every_nth: i32) -> Canvas {
    let scaling_factor = base_box_size_px / BASE_BOX_SIZE;
    let char_width = string_bounds(COORD_FONT_SIZE, "9").0;
    let widths: Vec<f64> =
        (render_region.start..=render_region.end).map(|i| string_bounds(COORD_FONT_SIZE, &grouped(i)).0).collect();
    let max_string_width = widths.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut canvas = Canvas::new(
        scaling_factor * BASE_BOX_SIZE * (render_region.base_length() + 2 * BOX_PADDING) as f64,
        scaling_factor * (max_string_width + 2.0 * char_width),
    );
    let (canvas_x, canvas_y) = (canvas.width, canvas.height);
    canvas.set_font_size(COORD_FONT_SIZE);
    canvas.set_color(LIGHT_GRAY);
    canvas.fill_rect(base_box_size_px * 0.0, base_box_size_px * 0.0, canvas_x + base_box_size_px * 0.0, canvas_y + base_box_size_px * 0.0);
    canvas.scale(scaling_factor, scaling_factor);
    let scaling_transform = canvas.transform();
    canvas.set_color(BLACK);
    for b in 0..render_region.base_length() {
        canvas.set_transform(scaling_transform);
        let box_idx = (b + BOX_PADDING) as f64;
        let base_pos = render_region.start + b;
        let base_pos_offset = base_pos - center_position;
        let x = BASE_BOX_SIZE * (box_idx + 0.5 - 0.25 * BOUNDARY_BOX_PROPORTION + 0.0);
        if base_pos_offset % display_every_nth != 0 {
            canvas.fill_rect(
                x,
                canvas_y / scaling_factor - 0.75 * char_width + BASE_BOX_SIZE * 0.0,
                BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
                0.75 * char_width,
            );
            continue;
        }
        canvas.fill_rect(
            x,
            canvas_y / scaling_factor - 1.5 * char_width + BASE_BOX_SIZE * 0.0,
            BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
            1.5 * char_width,
        );
        canvas.rotate_quarter();
        canvas.translate(max_string_width - 0.5 * widths[b as usize] + BASE_BOX_SIZE * 0.0, -BASE_BOX_SIZE * (box_idx + 0.5 + 0.0));
        draw_string_from_center(&mut canvas, &grouped(base_pos), 0.0, 0.0);
    }
    canvas
}

/// SvgRender.renderColoredBox
pub fn render_colored_box(size_px: f64, color: Color) -> Canvas {
    let mut canvas = Canvas::new(size_px, size_px);
    canvas.set_color(color);
    canvas.fill_rect(0.0, 0.0, size_px, size_px);
    canvas
}

/// GeneRegionViewModel
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeneRegionKind {
    AminoAcid { amino_acid_pos: i32, ref_: u8, alt: u8, is_start: bool },
    NonCodingExonic,
    Intronic,
    Del,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GeneRegionViewModel {
    pub region: BaseRegion,
    pub kind: GeneRegionKind,
    /// not defined for a Del (Java throws NotImplementedException)
    pub bases_before_right_insert: i32,
}

impl GeneRegionViewModel {
    pub fn amino_acid(region: BaseRegion, amino_acid_pos: i32, ref_: u8, alt: u8, is_start: bool) -> Self {
        GeneRegionViewModel { region, kind: GeneRegionKind::AminoAcid { amino_acid_pos, ref_, alt, is_start }, bases_before_right_insert: 0 }
    }
    pub fn non_coding(region: BaseRegion) -> Self {
        GeneRegionViewModel { region, kind: GeneRegionKind::NonCodingExonic, bases_before_right_insert: 0 }
    }
    pub fn intronic(region: BaseRegion) -> Self {
        GeneRegionViewModel { region, kind: GeneRegionKind::Intronic, bases_before_right_insert: 0 }
    }
    pub fn del(region: BaseRegion) -> Self {
        GeneRegionViewModel { region, kind: GeneRegionKind::Del, bases_before_right_insert: 0 }
    }
    pub fn is_amino_acid(&self) -> bool {
        matches!(self.kind, GeneRegionKind::AminoAcid { .. })
    }
    pub fn is_del(&self) -> bool {
        self.kind == GeneRegionKind::Del
    }
    pub fn with_region(mut self, region: BaseRegion) -> Self {
        self.region = region;
        self
    }
}

/// SvgRender.renderGeneData: the reference and the predicted amino acid tracks
pub fn render_gene_data(
    base_box_size_px: f64, render_region: BaseRegion, pos_strand: bool, ref_models: &[GeneRegionViewModel],
    alt_models: &[GeneRegionViewModel],
) -> (Canvas, Canvas) {
    let width = base_box_size_px * (render_region.base_length() + 2 * BOX_PADDING) as f64;
    let mut ref_canvas = Canvas::new(width, base_box_size_px);
    render_gene_regions(&mut ref_canvas, base_box_size_px, render_region, pos_strand, ref_models, true);
    let mut alt_canvas = Canvas::new(width, base_box_size_px);
    render_gene_regions(&mut alt_canvas, base_box_size_px, render_region, pos_strand, alt_models, false);
    (ref_canvas, alt_canvas)
}

fn render_gene_regions(
    canvas: &mut Canvas, base_box_size_px: f64, render_region: BaseRegion, pos_strand: bool, models: &[GeneRegionViewModel],
    render_ref: bool,
) {
    canvas.scale(base_box_size_px, base_box_size_px);
    canvas.set_font_size(BASE_FONT_SIZE);
    for (i, model) in models.iter().enumerate() {
        if model.region.start > render_region.end {
            break;
        }
        if model.region.end < render_region.start {
            continue;
        }
        match model.kind {
            GeneRegionKind::Intronic => {
                canvas.set_color(AMINO_ACID_MATCH_BG[0]);
                let current = canvas.stroke_width();
                canvas.set_stroke_width((2.0 / BASE_BOX_SIZE) as f32);
                let left = (model.region.start.max(render_region.start) - render_region.start + BOX_PADDING) as f64;
                let right = (model.region.end.min(render_region.end) + 1 - render_region.start + BOX_PADDING) as f64;
                canvas.draw_line(left, 0.5, right, 0.5);
                canvas.set_stroke_width(current);
                continue;
            }
            GeneRegionKind::NonCodingExonic => {
                canvas.set_color(AMINO_ACID_MATCH_BG[0]);
                let left = (model.region.start.max(render_region.start) - render_region.start + BOX_PADDING) as f64;
                let right = (model.region.end.min(render_region.end) - render_region.start + BOX_PADDING) as f64;
                canvas.fill_rect(left, 0.25, right - left + 1.0, 0.5);
                continue;
            }
            _ => {}
        }
        let start_pos = model.region.start.max(render_region.start);
        let end_pos = model.region.end.min(render_region.end);
        let box_width = end_pos - start_pos + 1;
        let box_idx = start_pos - render_region.start + BOX_PADDING;
        let GeneRegionKind::AminoAcid { ref_, alt, is_start, .. } = model.kind else {
            draw_del(canvas, None, box_idx, box_idx + box_width - 1);
            continue;
        };
        let aa = if render_ref { ref_ } else { alt };
        let matches_ref = render_ref || ref_ == alt;
        let bg = if aa == b'*' {
            RED
        } else if is_start {
            GREEN
        } else if matches_ref {
            AMINO_ACID_MATCH_BG[i % 2]
        } else {
            AMINO_ACID_MISMATCH_BG[i % 2]
        };
        canvas.set_color(bg);
        canvas.fill_rect(box_idx as f64, 0.0, box_width as f64, 1.0);
        let size = canvas.font_size();
        render_text(canvas, box_idx as f64 + 0.5 * box_width as f64, 0.5, size, WHITE, &(aa as char).to_string(), None, 1.0);
    }

    canvas.set_color(BLUE);
    let right_edge = (render_region.base_length() + BOX_PADDING) as f64;
    if pos_strand {
        draw_forward_arrow(canvas, 0.5, 0.0, 0.5, 1.0);
        draw_forward_arrow(canvas, right_edge, 0.0, 0.5, 1.0);
    } else {
        draw_reverse_arrow(canvas, 0.5, 0.0, 0.5, 1.0);
        draw_reverse_arrow(canvas, right_edge, 0.0, 0.5, 1.0);
    }
    if render_ref {
        return;
    }
    for model in models {
        if model.region.start > render_region.end {
            break;
        }
        if model.region.end < render_region.start {
            continue;
        }
        if model.is_del() || model.bases_before_right_insert <= 0 {
            continue;
        }
        let start_pos = model.region.start.max(render_region.start);
        let end_pos = model.region.end.min(render_region.end);
        let box_width = end_pos - start_pos + 1;
        let box_idx = start_pos - render_region.start + BOX_PADDING;
        let indel_pos = if pos_strand {
            model.region.start + model.bases_before_right_insert - 1
        } else {
            model.region.end - model.bases_before_right_insert + 1
        };
        if !render_region.contains(indel_pos) {
            continue;
        }
        let indel_box_idx = (indel_pos - render_region.start + BOX_PADDING) as f64;
        canvas.set_color(INSERT_COLOR);
        if pos_strand {
            draw_right_insert_indicator(canvas, indel_box_idx);
            draw_left_insert_indicator(canvas, indel_box_idx + 1.0);
        } else {
            draw_left_insert_indicator(canvas, indel_box_idx);
            draw_right_insert_indicator(canvas, indel_box_idx - 1.0);
        }
        let GeneRegionKind::AminoAcid { alt, .. } = model.kind else { continue };
        let size = canvas.font_size();
        render_text(canvas, box_idx as f64 + 0.5 * box_width as f64, 0.5, size, WHITE, &(alt as char).to_string(), None, 1.0);
    }
}
