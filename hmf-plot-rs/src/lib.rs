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
//! BUI is used only through `bui` and `raster`: the plotting code builds display lists.

pub mod bui;
pub mod color;
pub mod display;
pub mod fmath;
pub mod png_out;
pub mod raster;
pub mod text;
pub mod circos;
pub mod gg;

pub use color::Rgba;
pub use display::{Figure, Item, LineCap, Stroke};
pub use text::Font;
