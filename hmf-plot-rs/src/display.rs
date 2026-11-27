//! The display list: what a figure draws, in device pixels, independent of the renderer. The BUI raster
//! (`raster.rs`) turns it into pixels; a vector backend (PDF, SVG) can take the same list later.
//!
//! Coordinates are continuous: pixel (i, j) covers [i, i + 1) x [j, j + 1), y grows downwards.

use crate::color::Rgba;
use crate::text::Font;

/// How the ends of an open stroke look.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LineCap {
    Butt,
    Round,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Stroke {
    pub color: Rgba,
    /// line width in pixels
    pub width: f64,
    /// on/off lengths in pixels (empty: solid)
    pub dash: Vec<f64>,
    pub cap: LineCap,
}

impl Stroke {
    pub fn solid(color: Rgba, width: f64) -> Stroke {
        Stroke { color, width, dash: Vec::new(), cap: LineCap::Butt }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Item {
    /// An axis-aligned filled rectangle.
    Rect { x: f64, y: f64, w: f64, h: f64, fill: Rgba },
    /// A filled polygon (even-odd rule), optionally outlined.
    Polygon { points: Vec<(f64, f64)>, fill: Option<Rgba>, stroke: Option<Stroke> },
    /// An open polyline.
    Polyline { points: Vec<(f64, f64)>, stroke: Stroke },
    /// A circle, filled and/or outlined (the outline centred on the radius).
    Circle { cx: f64, cy: f64, r: f64, fill: Option<Rgba>, stroke: Option<Stroke> },
    /// A text run. (x, y) is the origin of the baseline (its left end before rotation); `angle` turns the
    /// run counter-clockwise about it, in degrees.
    Text { x: f64, y: f64, text: String, font: Font, size: f64, color: Rgba, angle: f64 },
    /// Restricts what follows to a rectangle (None: no clip).
    Clip(Option<(f64, f64, f64, f64)>),
}

/// A figure: its pixel size, background and display list.
#[derive(Clone, Debug)]
pub struct Figure {
    pub width: u32,
    pub height: u32,
    pub background: Rgba,
    pub items: Vec<Item>,
    /// what the PNG records about the image (resolution, transparent colour)
    pub meta: crate::png_out::PngMeta,
}

impl Figure {
    pub fn new(width: u32, height: u32, background: Rgba) -> Figure {
        Figure { width, height, background, items: Vec::new(), meta: crate::png_out::PngMeta { dpi: None, transparent_key: None } }
    }

    pub fn push(&mut self, item: Item) {
        self.items.push(item);
    }

    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: Rgba) {
        if !fill.is_transparent() && w > 0.0 && h > 0.0 {
            self.items.push(Item::Rect { x, y, w, h, fill });
        }
    }

    /// A rectangle outline, centred on the rectangle's edges as R draws it.
    pub fn rect_outline(&mut self, x: f64, y: f64, w: f64, h: f64, stroke: Stroke) {
        self.items.push(Item::Polygon {
            points: vec![(x, y), (x + w, y), (x + w, y + h), (x, y + h)],
            fill: None,
            stroke: Some(stroke),
        });
    }

    pub fn polygon(&mut self, points: Vec<(f64, f64)>, fill: Option<Rgba>, stroke: Option<Stroke>) {
        if points.len() >= 2 {
            self.items.push(Item::Polygon { points, fill, stroke });
        }
    }

    pub fn polyline(&mut self, points: Vec<(f64, f64)>, stroke: Stroke) {
        if points.len() >= 2 && !stroke.color.is_transparent() && stroke.width > 0.0 {
            self.items.push(Item::Polyline { points, stroke });
        }
    }

    pub fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: Stroke) {
        self.polyline(vec![(x1, y1), (x2, y2)], stroke);
    }

    pub fn circle(&mut self, cx: f64, cy: f64, r: f64, fill: Option<Rgba>, stroke: Option<Stroke>) {
        if r > 0.0 {
            self.items.push(Item::Circle { cx, cy, r, fill, stroke });
        }
    }

    pub fn text(&mut self, x: f64, y: f64, text: &str, font: Font, size: f64, color: Rgba, angle: f64) {
        if !text.is_empty() {
            self.items.push(Item::Text { x, y, text: text.to_string(), font, size, color, angle });
        }
    }

    pub fn clip(&mut self, clip: Option<(f64, f64, f64, f64)>) {
        self.items.push(Item::Clip(clip));
    }

    /// Renders with BUI and writes an RGB PNG.
    pub fn write_png(&self, path: &str) -> std::io::Result<()> {
        let pixels = crate::raster::render(self);
        crate::png_out::write_rgb(path, self.width, self.height, &pixels, self.meta)
    }
}
