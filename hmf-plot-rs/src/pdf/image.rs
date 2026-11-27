//! PNG images as PDF image XObjects. The pixels are those of the PNG: the colour data of an 8-bit,
//! non-interlaced RGB or grey PNG is embedded as it is (its zlib data is a PDF Flate stream with PNG predictors,
//! no recompression), with a soft mask from its tRNS colour key if it has one; any other PNG (alpha channel,
//! palette, 16 bits, interlaced) is decoded and embedded as 8-bit RGB or grey with its alpha as a soft mask. As
//! PDFBox's LosslessFactory, every PNG with transparency (alpha or tRNS) gets a soft mask.

use std::io::Cursor;

/// The pixel data of one image XObject
pub(crate) enum Pixels {
    /// IDAT data of the PNG: Flate with /DecodeParms << /Predictor 15 /Colors n /BitsPerComponent 8 /Columns w >>
    PngPredicted { zlib: Vec<u8>, colors: u8 },
    /// zlib-compressed samples, 8 bits per component, rows without filter bytes
    Flate { zlib: Vec<u8>, colors: u8 },
}

pub(crate) struct Image {
    pub width: u32,
    pub height: u32,
    pub pixels: Pixels,
    /// the zlib-compressed 8-bit alpha samples, when the PNG has an alpha channel or a tRNS key
    pub smask: Option<Vec<u8>>,
}

/// Flate level for decoded images: their streams are written once, and big plots (6000 x 3300) are common, so a
/// fast level (the PDF is a few percent larger than at the default level 6)
const FLATE_LEVEL: u8 = 3;

pub(crate) fn zlib(data: &[u8], level: u8) -> Vec<u8> {
    miniz_oxide::deflate::compress_to_vec_zlib(data, level)
}

impl Image {
    pub(crate) fn from_png(bytes: &[u8]) -> Result<Image, String> {
        // an 8-bit RGB or grey PNG: its colour data as it is (and with a tRNS key, the key's soft mask)
        if let Some(mut img) = pass_through(bytes) {
            if img.smask.is_some() {
                let (alpha, _, _) = decode(bytes)?;
                img.smask = Some(zlib(&alpha.unwrap_or_default(), FLATE_LEVEL));
            }
            return Ok(img);
        }
        let (alpha, color, (w, h, colors)) = decode(bytes)?;
        // as PDFBox's LosslessFactory: a soft mask for every image with an alpha channel (or a tRNS key)
        let smask = alpha.map(|a| zlib(&a, FLATE_LEVEL));
        Ok(Image { width: w, height: h, pixels: Pixels::Flate { zlib: zlib(&color, FLATE_LEVEL), colors }, smask })
    }
}

/// The PNG decoded to 8 bits (palettes and tRNS keys expanded): (alpha samples when it has transparency, colour
/// samples, (width, height, colour components))
#[allow(clippy::type_complexity)]
fn decode(bytes: &[u8]) -> Result<(Option<Vec<u8>>, Vec<u8>, (u32, u32, u8)), String> {
    let mut decoder = png::Decoder::new(Cursor::new(bytes));
    decoder.set_transformations(png::Transformations::EXPAND | png::Transformations::STRIP_16);
    let mut reader = decoder.read_info().map_err(|e| format!("cannot read PNG: {e}"))?;
    let mut buf = vec![0; reader.output_buffer_size()];
    let info = reader.next_frame(&mut buf).map_err(|e| format!("cannot decode PNG: {e}"))?;
    let (w, h) = (info.width, info.height);
    let n = (w as usize) * (h as usize);
    let (colors, has_alpha) = match info.color_type {
        png::ColorType::Grayscale => (1u8, false),
        png::ColorType::GrayscaleAlpha => (1, true),
        png::ColorType::Rgb => (3, false),
        png::ColorType::Rgba => (3, true),
        png::ColorType::Indexed => return Err("unexpected indexed PNG after expansion".into()),
    };
    let stride = colors as usize + has_alpha as usize;
    buf.truncate(n * stride);
    if !has_alpha {
        return Ok((None, buf, (w, h, colors)));
    }
    let mut color = Vec::with_capacity(n * colors as usize);
    let mut alpha = Vec::with_capacity(n);
    for px in buf.chunks_exact(stride) {
        color.extend_from_slice(&px[..colors as usize]);
        alpha.push(px[colors as usize]);
    }
    Ok((Some(alpha), color, (w, h, colors)))
}

/// An 8-bit, non-interlaced RGB or grey PNG: its IDAT data as it is; with a tRNS chunk the soft mask is marked
/// as needed (an empty one, filled by the caller)
fn pass_through(bytes: &[u8]) -> Option<Image> {
    if bytes.len() < 33 || &bytes[..8] != b"\x89PNG\r\n\x1a\n" {
        return None;
    }
    let mut pos = 8;
    let mut idat = Vec::new();
    let (mut width, mut height, mut colors, mut trns) = (0u32, 0u32, 0u8, false);
    while pos + 8 <= bytes.len() {
        let len = u32::from_be_bytes(bytes[pos..pos + 4].try_into().ok()?) as usize;
        let kind = &bytes[pos + 4..pos + 8];
        let data = bytes.get(pos + 8..pos + 8 + len)?;
        match kind {
            b"IHDR" => {
                width = u32::from_be_bytes(data[0..4].try_into().ok()?);
                height = u32::from_be_bytes(data[4..8].try_into().ok()?);
                let (depth, color_type, interlace) = (data[8], data[9], data[12]);
                if depth != 8 || interlace != 0 {
                    return None;
                }
                colors = match color_type {
                    0 => 1,
                    2 => 3,
                    _ => return None,
                };
            }
            b"tRNS" => trns = true,
            b"IDAT" => idat.extend_from_slice(data),
            b"IEND" => break,
            _ => {}
        }
        pos += 12 + len;
    }
    if colors == 0 || idat.is_empty() {
        return None;
    }
    Some(Image { width, height, pixels: Pixels::PngPredicted { zlib: idat, colors }, smask: if trns { Some(Vec::new()) } else { None } })
}
