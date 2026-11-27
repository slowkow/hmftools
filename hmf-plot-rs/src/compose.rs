//! Composing rendered images, as R scripts do with cowplot's `draw_image` and `plot_grid` (LINX's visualiser
//! stacks the Circos image, fusion panels and the chromosome strip). A `Canvas` is an RGB pixel buffer;
//! images are placed at whole pixels unscaled, or scaled with bilinear filtering (Cairo's filter for R's
//! interpolated rasters). Added for LINX.

use crate::display::Figure;
use crate::png_out::{self, PngMeta};

/// An RGB image.
#[derive(Clone, Debug)]
pub struct Canvas {
    pub width: u32,
    pub height: u32,
    pub rgb: Vec<u8>,
}

impl Canvas {
    /// A canvas filled with one colour.
    pub fn new(width: u32, height: u32, fill: [u8; 3]) -> Canvas {
        let mut rgb = Vec::with_capacity(width as usize * height as usize * 3);
        for _ in 0..(width as usize * height as usize) {
            rgb.extend_from_slice(&fill);
        }
        Canvas { width, height, rgb }
    }

    /// A figure rendered by BUI's raster.
    pub fn from_figure(fig: &Figure) -> Canvas {
        Canvas { width: fig.width, height: fig.height, rgb: crate::raster::render(fig) }
    }

    /// Copies `src` with its top-left pixel at (x, y), clipped to the canvas.
    pub fn blit(&mut self, src: &Canvas, x: i64, y: i64) {
        for sy in 0..src.height as i64 {
            let dy = y + sy;
            if dy < 0 || dy >= self.height as i64 {
                continue;
            }
            let x0 = x.max(0);
            let x1 = (x + src.width as i64).min(self.width as i64);
            if x1 <= x0 {
                continue;
            }
            let s = ((sy * src.width as i64 + (x0 - x)) * 3) as usize;
            let d = ((dy * self.width as i64 + x0) * 3) as usize;
            let n = ((x1 - x0) * 3) as usize;
            self.rgb[d..d + n].copy_from_slice(&src.rgb[s..s + n]);
        }
    }

    /// Draws `src` scaled into the rectangle (x, y, w, h) (device pixels, fractional) with bilinear filtering
    /// and padded edges: each covered pixel samples the source at its centre.
    pub fn blit_scaled(&mut self, src: &Canvas, x: f64, y: f64, w: f64, h: f64) {
        if w <= 0.0 || h <= 0.0 || src.width == 0 || src.height == 0 {
            return;
        }
        let (sx, sy) = (src.width as f64 / w, src.height as f64 / h);
        let px0 = x.floor().max(0.0) as i64;
        let px1 = ((x + w).ceil() as i64).min(self.width as i64);
        let py0 = y.floor().max(0.0) as i64;
        let py1 = ((y + h).ceil() as i64).min(self.height as i64);
        let get = |ix: i64, iy: i64| -> [f64; 3] {
            let ix = ix.clamp(0, src.width as i64 - 1) as usize;
            let iy = iy.clamp(0, src.height as i64 - 1) as usize;
            let i = (iy * src.width as usize + ix) * 3;
            [src.rgb[i] as f64, src.rgb[i + 1] as f64, src.rgb[i + 2] as f64]
        };
        for py in py0..py1 {
            // coverage of the destination pixel by the rectangle (edges partly covered blend with what is there)
            let cy = ((py as f64 + 1.0).min(y + h) - (py as f64).max(y)).clamp(0.0, 1.0);
            let v = (py as f64 + 0.5 - y) * sy - 0.5;
            let (v0, fv) = (v.floor(), v - v.floor());
            for px in px0..px1 {
                let cx = ((px as f64 + 1.0).min(x + w) - (px as f64).max(x)).clamp(0.0, 1.0);
                let cov = cx * cy;
                if cov <= 0.0 {
                    continue;
                }
                let u = (px as f64 + 0.5 - x) * sx - 0.5;
                let (u0, fu) = (u.floor(), u - u.floor());
                let (a, b, c, d) = (get(u0 as i64, v0 as i64), get(u0 as i64 + 1, v0 as i64), get(u0 as i64, v0 as i64 + 1), get(u0 as i64 + 1, v0 as i64 + 1));
                let di = ((py as usize) * self.width as usize + px as usize) * 3;
                for k in 0..3 {
                    let top = a[k] * (1.0 - fu) + b[k] * fu;
                    let bottom = c[k] * (1.0 - fu) + d[k] * fu;
                    let val = top * (1.0 - fv) + bottom * fv;
                    let old = self.rgb[di + k] as f64;
                    self.rgb[di + k] = (old * (1.0 - cov) + val * cov).round().clamp(0.0, 255.0) as u8;
                }
            }
        }
    }

    /// Writes the canvas as an RGB PNG.
    pub fn write_png(&self, path: &str, meta: PngMeta) -> std::io::Result<()> {
        png_out::write_rgb(path, self.width, self.height, &self.rgb, meta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blit_and_scale() {
        let mut c = Canvas::new(4, 4, [255, 255, 255]);
        let s = Canvas::new(2, 2, [0, 0, 0]);
        c.blit(&s, 3, 3);
        assert_eq!(&c.rgb[(3 * 4 + 3) * 3..(3 * 4 + 3) * 3 + 3], &[0, 0, 0]);
        assert_eq!(&c.rgb[0..3], &[255, 255, 255]);
        let mut c = Canvas::new(4, 4, [255, 255, 255]);
        c.blit_scaled(&s, 0.0, 0.0, 4.0, 4.0);
        assert!(c.rgb.iter().all(|v| *v == 0));
        // an unscaled copy at whole pixels equals blit
        let mut a = Canvas::new(5, 5, [9, 9, 9]);
        let mut b = a.clone();
        let mut src = Canvas::new(3, 2, [0, 0, 0]);
        src.rgb[0] = 200;
        src.rgb[10] = 100;
        a.blit(&src, 1, 2);
        b.blit_scaled(&src, 1.0, 2.0, 3.0, 2.0);
        assert_eq!(a.rgb, b.rgb);
    }
}
