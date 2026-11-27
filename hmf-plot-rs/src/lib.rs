//! hmf-plot-rs: native drawing for the hmftools Rust ports.
//!
//! - [`display`]: a figure is a pixel size and a display list (rectangles, polygons, polylines, circles,
//!   text, clips) in device pixels. Nothing in it depends on the renderer.
//! - [`raster`]: BUI's software raster turns a display list into pixels; [`png_out`] writes them.
//! - [`text`]: embedded fonts registered with BUI's text engine, and hinted text metrics.
//! - [`color`]: R's and Circos's colours, R's `hcl`, and scales' Lab colour ramps.
//! - [`gg`]: the parts of ggplot2 (theme_bw, scales, breaks, axes, legends, geoms, layout) the hmftools R
//!   scripts use.
//! - [`circos`]: a Circos ring engine driven by Circos configuration files: karyotype, ideograms, highlights,
//!   scatter, histogram, heatmap and text tracks, axes, backgrounds, links.
//!
//! - [`pdf`] (feature "pdf"): a PDF writer (pages, vector shapes, text in embedded font subsets, PNG images,
//!   outlines) and a table layout engine that reproduces boxable's (the Java tables library ORANGE uses).
//!
//! BUI is used only through `bui` and `raster`: the plotting code builds display lists. Everything that needs
//! BUI is behind the default feature "raster"; without it (`default-features = false, features = ["pdf"]`)
//! the crate is pure Rust and builds without CMake or a C compiler.

#[cfg(feature = "raster")]
pub mod bui;
pub mod color;
#[cfg(feature = "raster")]
pub mod display;
pub mod fmath;
pub mod png_out;
#[cfg(feature = "raster")]
pub mod raster;
#[cfg(feature = "raster")]
pub mod text;
#[cfg(feature = "raster")]
pub mod circos;
#[cfg(feature = "raster")]
pub mod gg;
#[cfg(feature = "pdf")]
pub mod pdf;

pub use color::Rgba;
#[cfg(feature = "raster")]
pub use display::{Figure, Item, LineCap, Stroke};
#[cfg(feature = "raster")]
pub use text::Font;
