//! Ports of the `hmf-common` `vis` package that the assembly visualiser draws with: `svg` writes what
//! JFreeSVG's `SVGGraphics2D` and `SvgRender` produce, `html` what j2html, `CssBuilder` and `HtmlUtil` do.

pub mod html;
pub mod svg;

/// A Java exception thrown while a visualisation is built; the message is the exception's `toString()`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VisError(pub String);

impl VisError {
    pub fn new(message: impl Into<String>) -> VisError {
        VisError(message.into())
    }
}

impl std::fmt::Display for VisError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

pub type VisResult<T> = Result<T, VisError>;

#[cfg(test)]
mod tests;
