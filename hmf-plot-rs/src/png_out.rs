//! PNG output: 8-bit RGB, no interlacing, fixed compression and filter settings and no time or text
//! chunks, so that the same pixels always give the same bytes. The resolution (pHYs) is written as R's png
//! device and Circos's GD write it: 300 dpi for the R plots, 96 dpi for Circos (which also marks the
//! colour 0,0,1 transparent with a tRNS chunk).

use std::io::{BufWriter, Write};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PngMeta {
    /// pixels per inch written in pHYs (None: no pHYs)
    pub dpi: Option<f64>,
    /// an RGB colour marked transparent (tRNS)
    pub transparent_key: Option<[u8; 3]>,
}

pub fn encode_rgb<W: Write>(out: W, width: u32, height: u32, rgb: &[u8], meta: PngMeta) -> std::io::Result<()> {
    let mut enc = png::Encoder::new(out, width, height);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.set_compression(png::Compression::Default);
    enc.set_filter(png::FilterType::Paeth);
    if let Some(dpi) = meta.dpi {
        let ppm = (dpi / 0.0254).round() as u32;
        enc.set_pixel_dims(Some(png::PixelDimensions { xppu: ppm, yppu: ppm, unit: png::Unit::Meter }));
    }
    if let Some([r, g, b]) = meta.transparent_key {
        enc.set_trns(vec![0, r, 0, g, 0, b]);
    }
    let mut w = enc.write_header().map_err(std::io::Error::other)?;
    w.write_image_data(rgb).map_err(std::io::Error::other)?;
    w.finish().map_err(std::io::Error::other)
}

pub fn write_rgb(path: &str, width: u32, height: u32, rgb: &[u8], meta: PngMeta) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut out = BufWriter::new(file);
    encode_rgb(&mut out, width, height, rgb, meta)?;
    out.flush()
}
