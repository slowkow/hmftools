//! The visualiser's HTML layer: ports of `hmf-common` `vis/CssBuilder`, `CssSize` and `HtmlUtil`, and element
//! writers that produce what j2html 1.6.0 `render()` does (no whitespace, attributes in the order they
//! were set, text and attribute values escaped).

use super::svg::Color;
use crate::common::fmt::java_fixed;
use std::collections::BTreeMap;

/// CssSize
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CssSize {
    pub magnitude: f64,
    pub unit: &'static str,
}

impl CssSize {
    pub const ZERO: CssSize = CssSize { magnitude: 0.0, unit: "px" };
    pub fn px(magnitude: f64) -> CssSize {
        CssSize { magnitude, unit: "px" }
    }
    pub fn em(magnitude: f64) -> CssSize {
        CssSize { magnitude, unit: "em" }
    }
    pub fn is_zero(&self) -> bool {
        (10000.0 * self.magnitude + 0.5).floor() as i64 == 0
    }
    pub fn scale(&self, factor: f64) -> CssSize {
        CssSize { magnitude: factor * self.magnitude, unit: self.unit }
    }
}

impl std::fmt::Display for CssSize {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        write!(f, "{}{}", java_fixed(self.magnitude, 4), self.unit)
    }
}

fn hex(color: Color) -> String {
    format!("#{:02x}{:02x}{:02x}", color.0, color.1, color.2)
}

/// CssBuilder: an immutable, sorted property map. Box shadows are keyed by direction
/// (top, right, bottom, left).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CssBuilder {
    properties: BTreeMap<&'static str, String>,
    box_shadows: BTreeMap<usize, (CssSize, Color)>,
}

const SHADOW_SCALES: [(f64, f64); 4] = [(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)];

impl CssBuilder {
    pub fn empty() -> CssBuilder {
        CssBuilder::default()
    }
    fn plus(&self, property: &'static str, value: String) -> CssBuilder {
        let mut next = self.clone();
        next.properties.insert(property, value);
        next
    }
    /// the other builder's values win
    pub fn merge(&self, other: &CssBuilder) -> CssBuilder {
        let mut next = self.clone();
        for (k, v) in &other.properties {
            next.properties.insert(k, v.clone());
        }
        for (k, v) in &other.box_shadows {
            next.box_shadows.insert(*k, *v);
        }
        next
    }
    pub fn text_align(&self, value: &str) -> CssBuilder {
        self.plus("text-align", value.to_string())
    }
    pub fn border_collapse(&self, value: &str) -> CssBuilder {
        self.plus("border-collapse", value.to_string())
    }
    pub fn margin_right(&self, size: CssSize) -> CssBuilder {
        self.plus("margin-right", size.to_string())
    }
    pub fn writing_mode(&self, value: &str) -> CssBuilder {
        self.plus("writing-mode", value.to_string())
    }
    pub fn font_family(&self, value: &str) -> CssBuilder {
        self.plus("font-family", value.to_string())
    }
    pub fn display(&self, value: &str) -> CssBuilder {
        self.plus("display", value.to_string())
    }
    pub fn font_weight(&self, value: &str) -> CssBuilder {
        self.plus("font-weight", value.to_string())
    }
    pub fn height(&self, size: CssSize) -> CssBuilder {
        self.plus("height", size.to_string())
    }
    pub fn padding(&self, size: CssSize) -> CssBuilder {
        self.plus("padding", size.to_string())
    }
    pub fn padding_right(&self, size: CssSize) -> CssBuilder {
        self.plus("padding-right", size.to_string())
    }
    pub fn padding_left(&self, size: CssSize) -> CssBuilder {
        self.plus("padding-left", size.to_string())
    }
    pub fn margin(&self, size: CssSize) -> CssBuilder {
        self.plus("margin", size.to_string())
    }
    pub fn border_spacing(&self, size: CssSize) -> CssBuilder {
        self.plus("border-spacing", size.to_string())
    }
    pub fn background_color(&self, color: Color) -> CssBuilder {
        self.plus("background-color", hex(color))
    }
    pub fn color(&self, color: Color) -> CssBuilder {
        self.plus("color", hex(color))
    }
    pub fn border_bottom(&self, size: CssSize, style: &str, color: Color) -> CssBuilder {
        self.plus("border-bottom", format!("{} {} {}", size, style, hex(color)))
    }
    pub fn border_top(&self, size: CssSize, style: &str, color: Color) -> CssBuilder {
        self.plus("border-top", format!("{} {} {}", size, style, hex(color)))
    }
    pub fn border_left(&self, size: CssSize, style: &str, color: Color) -> CssBuilder {
        self.plus("border-left", format!("{} {} {}", size, style, hex(color)))
    }
    pub fn border_right(&self, size: CssSize, style: &str, color: Color) -> CssBuilder {
        self.plus("border-right", format!("{} {} {}", size, style, hex(color)))
    }
    pub fn transform(&self, value: &str) -> CssBuilder {
        self.plus("transform", value.to_string())
    }
    pub fn border(&self, size: CssSize, style: &str, color: Color) -> CssBuilder {
        self.plus("border", format!("{} {} {}", size, style, hex(color)))
    }
    pub fn no_border(&self) -> CssBuilder {
        self.plus("border", "none".to_string())
    }
    pub fn vertical_align(&self, value: &str) -> CssBuilder {
        self.plus("vertical-align", value.to_string())
    }
    pub fn font_size_pt(&self, pt: i32) -> CssBuilder {
        self.plus("font-size", format!("{pt}pt"))
    }
    pub fn font_size_pct(&self, pct: f64) -> CssBuilder {
        self.plus("font-size", format!("{}%", java_fixed(100.0 * pct, 2)))
    }
    pub fn left_box_shadow(&self, size: CssSize, color: Color) -> CssBuilder {
        let mut next = self.clone();
        next.box_shadows.insert(3, (size, color));
        next
    }
    pub fn right_box_shadow(&self, size: CssSize, color: Color) -> CssBuilder {
        let mut next = self.clone();
        next.box_shadows.insert(1, (size, color));
        next
    }
}

impl std::fmt::Display for CssBuilder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (property, value) in &self.properties {
            write!(f, "{property}:{value};")?;
        }
        if !self.box_shadows.is_empty() {
            let shadows: Vec<String> = self
                .box_shadows
                .iter()
                .map(|(direction, (size, color))| {
                    let (h, v) = SHADOW_SCALES[*direction];
                    format!("{} {} 0 0 {}", size.scale(h), size.scale(v), hex(*color))
                })
                .collect();
            write!(f, "box-shadow:{};", shadows.join(", "))?;
        }
        Ok(())
    }
}

/// j2html's text and attribute escaping
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#x27;"),
            _ => out.push(c),
        }
    }
    out
}

/// `<name attr="value"...>inner</name>`; `inner` is already rendered
pub fn tag(name: &str, attrs: &[(&str, &str)], inner: &str) -> String {
    let mut out = String::with_capacity(inner.len() + 64);
    push_tag(&mut out, name, attrs, inner);
    out
}

pub fn push_tag(out: &mut String, name: &str, attrs: &[(&str, &str)], inner: &str) {
    out.push('<');
    out.push_str(name);
    for (key, value) in attrs {
        out.push(' ');
        out.push_str(key);
        out.push_str("=\"");
        out.push_str(&escape(value));
        out.push('"');
    }
    out.push('>');
    out.push_str(inner);
    out.push_str("</");
    out.push_str(name);
    out.push('>');
}

/// `td(text)`
pub fn td(text: &str) -> String {
    tag("td", &[], &escape(text))
}

/// `td(...).withStyle(style)`; `inner` is already rendered
pub fn td_styled(inner: &str, style: &CssBuilder) -> String {
    tag("td", &[("style", &style.to_string())], inner)
}

/// HtmlUtil.BASE_FONT_STYLE
pub fn base_font_style() -> CssBuilder {
    CssBuilder::empty().font_size_pt(10).font_family("sans-serif")
}

/// HtmlUtil.styledTable
pub fn styled_table(rows: &str, style: &CssBuilder) -> String {
    tag("table", &[("style", &base_font_style().merge(style).to_string())], rows)
}

/// HtmlUtil.getJavascript: the script element with the lines of `vis/sagevis.js` joined by line feeds
pub fn javascript() -> String {
    let content = include_str!("sagevis.js").lines().collect::<Vec<_>>().join("\n");
    tag("script", &[("type", "text/javascript")], &content)
}

/// The arguments of HtmlUtil.renderReadInfoTable; `None` is Java's null.
#[derive(Clone, Debug, Default)]
pub struct ReadInfo<'a> {
    pub read_name: &'a str,
    pub chromosome: &'a str,
    pub alignment_start: i32,
    pub alignment_end: i32,
    /// chromosome, start and end of the mate; an end of zero (NO_POSITION) is shown as "?"
    pub mate_alignment: Option<(&'a str, i32, i32)>,
    pub cigar: &'a str,
    pub mate_cigar: Option<&'a str>,
    pub insert_size: Option<i32>,
    pub orientation: Option<&'a str>,
    pub map_q: Option<i32>,
    pub read_nm: Option<i32>,
    pub consensus_type: Option<&'a str>,
    pub consensus_read: Option<&'a str>,
    pub second_map_q: Option<i32>,
    pub second_read_nm: Option<i32>,
    pub extra_info: &'a [(&'a str, String)],
}

/// HtmlUtil.renderReadInfoTable
pub fn render_read_info_table(info: &ReadInfo) -> String {
    let base_div_style = CssBuilder::empty().padding(CssSize::ZERO).margin(CssSize::ZERO);
    let read_info_style = base_div_style.display("none");
    let row = |key: &str, value: &str| tag("tr", &[], &(td(key) + &td(value)));
    let mut rows = String::new();
    rows.push_str(&row("Read name:", info.read_name));
    let alignment_str = format!("{}:{}-{}", info.chromosome, info.alignment_start, info.alignment_end);
    let mate_alignment_str = match info.mate_alignment {
        Some((chromosome, start, end)) => format!("{}:{}-{}", chromosome, start, if end == 0 { "?".to_string() } else { end.to_string() }),
        None => "unmapped".to_string(),
    };
    rows.push_str(&row("Alignment:", &format!("{alignment_str}, {mate_alignment_str}")));
    rows.push_str(&row("Cigar:", &format!("{}, {}", info.cigar, info.mate_cigar.unwrap_or("missing"))));
    rows.push_str(&row("Insert size:", &info.insert_size.map_or("unknown".to_string(), |v| v.to_string())));
    rows.push_str(&row("Orientation:", info.orientation.unwrap_or("unknown")));
    let first_map_q = info.map_q.map_or("unknown".to_string(), |v| v.to_string());
    rows.push_str(&row("MapQ:", &match info.second_map_q {
        Some(second) => format!("{first_map_q}, {second}"),
        None => first_map_q,
    }));
    let first_nm = info.read_nm.map_or("unknown".to_string(), |v| v.to_string());
    rows.push_str(&row("NM:", &match info.second_read_nm {
        Some(second) => format!("{first_nm}, {second}"),
        None => first_nm,
    }));
    if let Some(consensus_type) = info.consensus_type {
        rows.push_str(&row("Dup type:", consensus_type));
    }
    // String.split(CONSENSUS_INFO_DELIM, 2)[0]
    let dup_count = info.consensus_read.map_or("0", |v| v.split(';').next().unwrap_or(""));
    rows.push_str(&row("Dup count:", dup_count));
    for (key, value) in info.extra_info {
        rows.push_str(&row(key, value));
    }
    let table = styled_table(&rows, &CssBuilder::empty());
    tag("div", &[("class", "read-info"), ("style", &read_info_style.to_string())], &table)
}
