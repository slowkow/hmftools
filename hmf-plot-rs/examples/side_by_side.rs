//! Puts two PNGs side by side (left: reference, right: ours) under a caption band:
//! `side_by_side <left.png> <right.png> <out.png> [x y w h zoom]` (the optional crop, magnified `zoom` times
//! with nearest pixels, is taken from both).
use hmf_plot::{Figure, Font, Rgba};

fn read(path: &str) -> (usize, usize, Vec<u8>) {
    let dec = png::Decoder::new(std::io::BufReader::new(std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"))));
    let mut dec = dec;
    dec.set_transformations(png::Transformations::EXPAND | png::Transformations::STRIP_16);
    let mut r = dec.read_info().unwrap();
    let mut buf = vec![0; r.output_buffer_size()];
    let info = r.next_frame(&mut buf).unwrap();
    let (w, h) = (info.width as usize, info.height as usize);
    let ch = info.color_type.samples();
    let mut rgb = Vec::with_capacity(w * h * 3);
    for p in buf[..w * h * ch].chunks(ch) {
        match ch {
            1 | 2 => rgb.extend_from_slice(&[p[0], p[0], p[0]]),
            _ => rgb.extend_from_slice(&p[..3]),
        }
    }
    (w, h, rgb)
}

fn crop(img: &(usize, usize, Vec<u8>), x: usize, y: usize, w: usize, h: usize, zoom: usize) -> (usize, usize, Vec<u8>) {
    let (iw, ih, px) = img;
    let mut out = Vec::with_capacity(w * h * zoom * zoom * 3);
    for yy in 0..h * zoom {
        for xx in 0..w * zoom {
            let (sx, sy) = (x + xx / zoom, y + yy / zoom);
            if sx < *iw && sy < *ih {
                let i = (sy * iw + sx) * 3;
                out.extend_from_slice(&px[i..i + 3]);
            } else {
                out.extend_from_slice(&[255, 255, 255]);
            }
        }
    }
    (w * zoom, h * zoom, out)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (mut l, mut r) = (read(&a[1]), read(&a[2]));
    let mut note = String::new();
    if a.len() >= 9 {
        let n = |i: usize| a[i].parse::<usize>().unwrap();
        let (x, y, w, h, z) = (n(4), n(5), n(6), n(7), n(8));
        l = crop(&l, x, y, w, h, z);
        r = crop(&r, x, y, w, h, z);
        note = format!("  crop {x},{y} {w}x{h} x{z}");
    }
    let gap = 20;
    let band = 60;
    let width = l.0 + gap + r.0;
    let height = band + l.1.max(r.1);
    let mut cap = Figure::new(width as u32, band as u32, Rgba::rgb8(235, 235, 235));
    cap.text(10.0, 42.0, &format!("Java (R / Circos){note}"), Font::Sans, 32.0, Rgba::BLACK, 0.0);
    cap.text((l.0 + gap + 10) as f64, 42.0, &format!("Rust (hmf-plot-rs){}", if l.0 != r.0 || l.1 != r.1 { format!("  size differs: {}x{} vs {}x{}", l.0, l.1, r.0, r.1) } else { String::new() }), Font::Sans, 32.0, Rgba::BLACK, 0.0);
    let capx = hmf_plot::raster::render(&cap);
    let mut out = vec![255u8; width * height * 3];
    out[..width * band * 3].copy_from_slice(&capx);
    for (img, x0) in [(&l, 0usize), (&r, l.0 + gap)] {
        for y in 0..img.1 {
            let src = &img.2[y * img.0 * 3..(y + 1) * img.0 * 3];
            let o = ((band + y) * width + x0) * 3;
            out[o..o + img.0 * 3].copy_from_slice(src);
        }
    }
    hmf_plot::png_out::write_rgb(&a[3], width as u32, height as u32, &out, hmf_plot::png_out::PngMeta { dpi: None, transparent_key: None }).unwrap();
}
