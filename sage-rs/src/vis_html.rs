//! The visualiser's HTML layer: ports of `hmf-common` `vis/CssBuilder` and `CssSize`, and element
//! writers that produce what j2html 1.6.0 `render()` does (no whitespace, attributes in the order they
//! were set, text and attribute values escaped).

use crate::fmt::java_fixed;
use crate::vis_svg::Color;
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
        crate::constants::jround(10000.0 * self.magnitude) == 0
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
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn font_size_pct(&self, pct: f64) -> CssBuilder {
        self.plus("font-size", format!("{}%", java_fixed(100.0 * pct, 2)))
    }
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn left_box_shadow(&self, size: CssSize, color: Color) -> CssBuilder {
        let mut next = self.clone();
        next.box_shadows.insert(3, (size, color));
        next
    }
    #[cfg_attr(not(test), allow(dead_code))]
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
