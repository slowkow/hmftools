//! Colours: straight-alpha RGBA, R's colour names and functions (grey levels, `hcl`, the Lab interpolation
//! of `scales::gradient_n_pal`), Circos's named colours with `_aN` transparency, and HSV/hex parsing.

use crate::fmath::FMath;

/// A colour with straight (not premultiplied) alpha, every channel in 0..=1.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rgba {
    pub r: f32,
    pub g: f32,
    pub b: f32,
    pub a: f32,
}

impl Rgba {
    pub const BLACK: Rgba = Rgba { r: 0.0, g: 0.0, b: 0.0, a: 1.0 };
    pub const WHITE: Rgba = Rgba { r: 1.0, g: 1.0, b: 1.0, a: 1.0 };
    pub const TRANSPARENT: Rgba = Rgba { r: 0.0, g: 0.0, b: 0.0, a: 0.0 };

    pub fn rgb8(r: u8, g: u8, b: u8) -> Rgba {
        Rgba { r: r as f32 / 255.0, g: g as f32 / 255.0, b: b as f32 / 255.0, a: 1.0 }
    }

    pub fn rgba8(r: u8, g: u8, b: u8, a: u8) -> Rgba {
        Rgba { r: r as f32 / 255.0, g: g as f32 / 255.0, b: b as f32 / 255.0, a: a as f32 / 255.0 }
    }

    /// The same colour at `alpha` times its opacity (ggplot2's `alpha` aesthetic replaces the alpha of an
    /// opaque colour: `scales::alpha(colour, alpha)`).
    pub fn with_alpha(self, alpha: f32) -> Rgba {
        Rgba { a: alpha, ..self }
    }

    pub fn is_transparent(&self) -> bool {
        self.a <= 0.0
    }

    /// 8-bit channels, rounded as R does when it writes a colour as "#RRGGBB".
    pub fn to_u8(&self) -> [u8; 4] {
        let c = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
        [c(self.r), c(self.g), c(self.b), c(self.a)]
    }

    pub fn to_array(&self) -> [f32; 4] {
        [self.r, self.g, self.b, self.a]
    }

    /// "#RRGGBB" or "#RRGGBBAA" (6 or 8 hex digits), as R reads them.
    pub fn from_hex(s: &str) -> Option<Rgba> {
        let h = s.strip_prefix('#')?;
        let byte = |i: usize| u8::from_str_radix(h.get(i..i + 2)?, 16).ok();
        match h.len() {
            6 => Some(Rgba::rgb8(byte(0)?, byte(2)?, byte(4)?)),
            8 => Some(Rgba::rgba8(byte(0)?, byte(2)?, byte(4)?, byte(6)?)),
            _ => None,
        }
    }

    pub fn to_hex(&self) -> String {
        let [r, g, b, _] = self.to_u8();
        format!("#{r:02X}{g:02X}{b:02X}")
    }
}

/// R's grey levels: `greyN` / `grayN` is round(N * 255 / 100) in every channel (grDevices colors.c).
pub fn r_grey(level: u32) -> Rgba {
    let v = ((level as f64 * 255.0 / 100.0) + 0.5) as u8;
    Rgba::rgb8(v, v, v)
}

/// The R colour names the PURPLE scripts and ggplot2's themes use, plus hex strings and grey levels.
pub fn r_colour(name: &str) -> Option<Rgba> {
    if name.starts_with('#') {
        return Rgba::from_hex(name);
    }
    let lower = name.to_ascii_lowercase();
    for prefix in ["grey", "gray"] {
        if let Some(n) = lower.strip_prefix(prefix) {
            if n.is_empty() {
                return Some(Rgba::rgb8(190, 190, 190));
            }
            if let Ok(level) = n.parse::<u32>() {
                if level <= 100 {
                    return Some(r_grey(level));
                }
            }
        }
    }
    let rgb = match lower.as_str() {
        "black" => (0, 0, 0),
        "white" => (255, 255, 255),
        "red" | "red1" => (255, 0, 0),
        "red2" => (238, 0, 0),
        "red3" => (205, 0, 0),
        "green" | "green1" => (0, 255, 0),
        "blue" | "blue1" => (0, 0, 255),
        "yellow" | "yellow1" => (255, 255, 0),
        "orange" | "orange1" => (255, 165, 0),
        "purple" => (160, 32, 240),
        "cyan" => (0, 255, 255),
        "magenta" => (255, 0, 255),
        _ => return None,
    };
    Some(Rgba::rgb8(rgb.0, rgb.1, rgb.2))
}

// ---------------------------------------------------------------------------------------------------------
// sRGB, CIE XYZ, Lab and polar Luv (R's hcl)

/// sRGB gamma decoding of one channel in 0..=1.
fn srgb_to_linear(c: f64) -> f64 {
    if c <= 0.04045 { c / 12.92 } else { ((c + 0.055) / 1.055).fpow(2.4) }
}

fn linear_to_srgb(c: f64) -> f64 {
    if c <= 0.0031308 { 12.92 * c } else { 1.055 * c.fpow(1.0 / 2.4) - 0.055 }
}

/// D65 white, X and Z scaled to Y = 100 (farver's and R's reference white).
const WHITE_D65: [f64; 3] = [95.047, 100.0, 108.883];

fn rgb_to_xyz(rgb: [f64; 3]) -> [f64; 3] {
    let [r, g, b] = rgb.map(srgb_to_linear);
    [
        (r * 0.4124 + g * 0.3576 + b * 0.1805) * 100.0,
        (r * 0.2126 + g * 0.7152 + b * 0.0722) * 100.0,
        (r * 0.0193 + g * 0.1192 + b * 0.9505) * 100.0,
    ]
}

fn xyz_to_rgb(xyz: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = xyz.map(|v| v / 100.0);
    [x * 3.2406 + y * -1.5372 + z * -0.4986, x * -0.9689 + y * 1.8758 + z * 0.0415, x * 0.0557 + y * -0.2040 + z * 1.0570]
        .map(linear_to_srgb)
}

fn lab_f(t: f64) -> f64 {
    if t > 0.008856 { t.fcbrt() } else { 7.787 * t + 16.0 / 116.0 }
}

fn lab_f_inv(t: f64) -> f64 {
    let t3 = t * t * t;
    if t3 > 0.008856 { t3 } else { (t - 16.0 / 116.0) / 7.787 }
}

/// sRGB (0..=1) to CIE Lab with D65 white, as the ColorSpace library inside farver computes it.
pub fn rgb_to_lab(rgb: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = rgb_to_xyz(rgb);
    let fx = lab_f(x / WHITE_D65[0]);
    let fy = lab_f(y / WHITE_D65[1]);
    let fz = lab_f(z / WHITE_D65[2]);
    [116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

pub fn lab_to_rgb(lab: [f64; 3]) -> [f64; 3] {
    let fy = (lab[0] + 16.0) / 116.0;
    let fx = lab[1] / 500.0 + fy;
    let fz = fy - lab[2] / 200.0;
    xyz_to_rgb([lab_f_inv(fx) * WHITE_D65[0], lab_f_inv(fy) * WHITE_D65[1], lab_f_inv(fz) * WHITE_D65[2]])
}

/// farver's encode_colour: channels rounded to the nearest of 0..=255, clamped.
fn encode_rgb(rgb: [f64; 3]) -> Rgba {
    let c = |v: f64| (v * 255.0).round().clamp(0.0, 255.0) as u8;
    Rgba::rgb8(c(rgb[0]), c(rgb[1]), c(rgb[2]))
}

/// R's `hcl(h, c, l)`: polar CIE-Luv with D65 white, gamma-corrected sRGB, channels rounded to 8 bits
/// (grDevices colors.c hcl2rgb).
pub fn hcl(h: f64, c: f64, l: f64) -> Rgba {
    if l <= 0.0 {
        return Rgba::BLACK;
    }
    let h = h.to_radians();
    let (u, v) = (c * h.fcos(), c * h.fsin());
    // CIE-Luv to CIE-XYZ
    let y = WHITE_D65[1] * if l > 7.999592 { ((l + 16.0) / 116.0).fpowi(3) } else { l / 903.3 };
    let t = WHITE_D65[0] + WHITE_D65[1] + WHITE_D65[2];
    let x = WHITE_D65[0] / t;
    let yw = WHITE_D65[1] / t;
    let uw = 2.0 * x / (6.0 * yw - x + 1.5);
    let vw = 4.5 * yw / (6.0 * yw - x + 1.5);
    let uu = u / (13.0 * l) + uw;
    let vv = v / (13.0 * l) + vw;
    let xx = 9.0 * y * uu / (4.0 * vv);
    let zz = -xx / 3.0 - 5.0 * y + 3.0 * y / vv;
    // grDevices gtrans: sRGB gamma
    let g = |u: f64| if u > 0.00304 { 1.055 * u.fpow(1.0 / 2.4) - 0.055 } else { 12.92 * u };
    let r = g((3.240479 * xx - 1.537150 * y - 0.498535 * zz) / WHITE_D65[1]);
    let gg = g((-0.969256 * xx + 1.875992 * y + 0.041556 * zz) / WHITE_D65[1]);
    let b = g((0.055648 * xx - 0.204043 * y + 1.057311 * zz) / WHITE_D65[1]);
    // R: ScaleColor = (unsigned int)(255 * x + 0.5) after fixing to [0, 1]
    let s = |v: f64| (255.0 * v.clamp(0.0, 1.0) + 0.5) as u8;
    Rgba::rgb8(s(r), s(gg), s(b))
}

/// ggplot2's default discrete colour scale (`scales::hue_pal()(n)`): hues from 15 to 375 degrees in n + 1
/// steps (the last dropped when the range is a full circle), chroma 100, luminance 65.
pub fn hue_pal(n: usize) -> Vec<Rgba> {
    if n == 0 {
        return Vec::new();
    }
    let (h0, h1) = (15.0, 375.0);
    // scales::hue_pal: if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360 / n
    let h1 = if ((h1 - h0) % 360.0) < 1.0 { h1 - 360.0 / n as f64 } else { h1 };
    (0..n)
        .map(|i| {
            let h = if n == 1 { h0 } else { h0 + (h1 - h0) * i as f64 / (n - 1) as f64 };
            hcl(h % 360.0, 100.0, 65.0)
        })
        .collect()
}

/// `scales::gradient_n_pal(colours, values)`: a colour ramp interpolated linearly in CIE Lab between the
/// colours at positions `values` (0..=1; evenly spaced when absent).
#[derive(Clone, Debug)]
pub struct GradientN {
    stops: Vec<f64>,
    labs: Vec<[f64; 3]>,
}

impl GradientN {
    pub fn new(colours: &[Rgba], values: Option<&[f64]>) -> GradientN {
        let n = colours.len();
        let stops = match values {
            Some(v) => v.to_vec(),
            None => (0..n).map(|i| if n == 1 { 0.0 } else { i as f64 / (n - 1) as f64 }).collect(),
        };
        let labs = colours.iter().map(|c| rgb_to_lab([c.r as f64, c.g as f64, c.b as f64])).collect();
        GradientN { stops, labs }
    }

    /// The colour at x (0..=1); outside, NA (None), as scales' ramp returns na.value.
    pub fn at(&self, x: f64) -> Option<Rgba> {
        if !x.is_finite() || x < 0.0 || x > 1.0 {
            return None;
        }
        // stats::approxfun(values, channel): linear, ties averaged (scales passes sorted values)
        let s = &self.stops;
        let mut i = 0;
        while i + 1 < s.len() && x > s[i + 1] {
            i += 1;
        }
        if i + 1 >= s.len() {
            return Some(encode_rgb(lab_to_rgb(self.labs[s.len() - 1])));
        }
        let (x0, x1) = (s[i], s[i + 1]);
        let t = if x1 > x0 { (x - x0) / (x1 - x0) } else { 0.0 };
        let l = |k: usize| self.labs[i][k] + (self.labs[i + 1][k] - self.labs[i][k]) * t;
        Some(encode_rgb(lab_to_rgb([l(0), l(1), l(2)])))
    }
}

// ---------------------------------------------------------------------------------------------------------
// Circos colours

/// The Circos 0.69 colours of etc/colors.conf (with its Brewer palettes) that PURPLE's and LINX's
/// configurations name, as r, g, b.
const CIRCOS_COLOURS: &[(&str, (u8, u8, u8))] = &[
    ("white", (255, 255, 255)),
    ("black", (0, 0, 0)),
    ("vvvvlgrey", (255, 255, 255)),
    ("vvvlgrey", (255, 255, 255)),
    ("vvlgrey", (240, 240, 240)),
    ("vlgrey", (217, 217, 217)),
    ("lgrey", (189, 189, 189)),
    ("grey", (150, 150, 150)),
    ("dgrey", (115, 115, 115)),
    ("vdgrey", (82, 82, 82)),
    ("vvdgrey", (37, 37, 37)),
    ("vvvdgrey", (0, 0, 0)),
    ("vvvvdgrey", (0, 0, 0)),
    ("vlred", (252, 187, 161)),
    ("lred", (252, 146, 114)),
    ("red", (251, 106, 74)),
    ("dred", (239, 59, 44)),
    ("vdred", (203, 24, 29)),
    ("vlgreen", (199, 233, 192)),
    ("lgreen", (161, 217, 155)),
    ("green", (116, 196, 118)),
    ("dgreen", (65, 171, 93)),
    ("vdgreen", (35, 139, 69)),
    ("vlblue", (198, 219, 239)),
    ("lblue", (158, 202, 225)),
    ("blue", (107, 174, 214)),
    ("dblue", (66, 146, 198)),
    ("vdblue", (33, 113, 181)),
    ("vlpurple", (218, 218, 235)),
    ("lpurple", (188, 189, 220)),
    ("purple", (158, 154, 200)),
    ("dpurple", (128, 125, 186)),
    ("vdpurple", (106, 81, 163)),
    ("vlorange", (253, 208, 162)),
    ("lorange", (253, 174, 107)),
    ("orange", (253, 141, 60)),
    ("dorange", (241, 105, 19)),
    ("vdorange", (217, 72, 1)),
    ("yellow", (255, 244, 80)),
    ("vdyellow", (227, 210, 0)),
];

/// The Brewer palettes Circos names by "<palette>-<n>-<type>-<i>" that the colour lists here use.
const BREWER: &[(&str, &[(u8, u8, u8)])] = &[(
    "spectral-11-div",
    &[
        (158, 1, 66),
        (213, 62, 79),
        (244, 109, 67),
        (253, 174, 97),
        (254, 224, 139),
        (255, 255, 191),
        (230, 245, 152),
        (171, 221, 164),
        (102, 194, 165),
        (50, 136, 189),
        (94, 79, 162),
    ],
)];

/// A Circos colour table: the defaults, plus the `<colors>` a configuration adds (`name = r,g,b`).
#[derive(Clone, Debug, Default)]
pub struct CircosColours {
    extra: Vec<(String, Rgba)>,
    auto_alpha_steps: u32,
}

impl CircosColours {
    pub fn new(auto_alpha_steps: u32) -> CircosColours {
        CircosColours { extra: Vec::new(), auto_alpha_steps }
    }

    /// Adds or replaces a named colour, given as "r,g,b", "r,g,b,a" or another colour's name.
    pub fn define(&mut self, name: &str, value: &str) {
        if let Some(c) = self.resolve(value) {
            let name = name.to_ascii_lowercase();
            self.extra.retain(|(n, _)| *n != name);
            self.extra.push((name, c));
        }
    }

    /// The colour a Circos colour expression names: "r,g,b" (0..255), "r,g,b,a" (a: opacity 0..1), a name,
    /// or a name with "_aN" (opacity 1 - N / (auto_alpha_steps + 1)); None when unknown.
    pub fn resolve(&self, spec: &str) -> Option<Rgba> {
        // Circos's colour names are case-insensitive (the karyotypes name "chrx", configurations "chrX")
        let lower = spec.trim().to_ascii_lowercase();
        let spec = lower.as_str();
        if spec.contains(',') {
            let parts: Vec<&str> = spec.split(',').map(|s| s.trim()).collect();
            if parts.len() == 3 || parts.len() == 4 {
                let ch = |s: &str| s.parse::<f64>().ok().map(|v| v.clamp(0.0, 255.0).round() as u8);
                let mut c = Rgba::rgb8(ch(parts[0])?, ch(parts[1])?, ch(parts[2])?);
                if parts.len() == 4 {
                    c.a = parts[3].parse::<f32>().ok()?.clamp(0.0, 1.0);
                }
                return Some(c);
            }
            return None;
        }
        if let Some((base, n)) = spec.rsplit_once("_a") {
            if let Ok(n) = n.parse::<u32>() {
                let steps = self.auto_alpha_steps.max(1);
                let base = self.resolve(base)?;
                return Some(base.with_alpha(1.0 - n as f32 / (steps + 1) as f32));
            }
        }
        if let Some((_, c)) = self.extra.iter().rev().find(|(n, _)| n == spec) {
            return Some(*c);
        }
        if let Some((_, (r, g, b))) = CIRCOS_COLOURS.iter().find(|(n, _)| *n == spec) {
            return Some(Rgba::rgb8(*r, *g, *b));
        }
        for (name, list) in BREWER {
            if let Some(i) = spec.strip_prefix(name).and_then(|s| s.strip_prefix('-')).and_then(|s| s.parse::<usize>().ok()) {
                if (1..=list.len()).contains(&i) {
                    let (r, g, b) = list[i - 1];
                    return Some(Rgba::rgb8(r, g, b));
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_greys() {
        assert_eq!(r_colour("grey92").unwrap().to_hex(), "#EBEBEB");
        assert_eq!(r_colour("grey20").unwrap().to_hex(), "#333333");
        assert_eq!(r_colour("grey30").unwrap().to_hex(), "#4D4D4D");
        assert_eq!(r_colour("red2").unwrap().to_hex(), "#EE0000");
    }

    #[test]
    fn hue_palette_matches_scales() {
        // scales::hue_pal()(3) and (4) in R 4.5 / scales 1.4
        let p: Vec<String> = hue_pal(3).iter().map(|c| c.to_hex()).collect();
        assert_eq!(p, ["#F8766D", "#00BA38", "#619CFF"]);
        let p: Vec<String> = hue_pal(4).iter().map(|c| c.to_hex()).collect();
        assert_eq!(p, ["#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"]);
    }

    #[test]
    fn lab_round_trip() {
        for c in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.2, 0.4, 0.6], [1.0, 1.0, 1.0]] {
            let back = lab_to_rgb(rgb_to_lab(c));
            for k in 0..3 {
                assert!((back[k] - c[k]).abs() < 1e-3, "{c:?} {back:?}");
            }
        }
    }

    #[test]
    fn gradient_ends_and_middle() {
        let g = GradientN::new(&[Rgba::rgb8(0, 0, 255), Rgba::rgb8(255, 0, 0)], None);
        assert_eq!(g.at(0.0).unwrap().to_hex(), "#0000FF");
        assert_eq!(g.at(1.0).unwrap().to_hex(), "#FF0000");
        assert!(g.at(1.5).is_none());
        // scales::gradient_n_pal(c("blue", "red"))(0.5) == "#C1009D" (Lab interpolation; R 4.5, scales 1.4)
        assert_eq!(g.at(0.5).unwrap().to_hex(), "#CA0088");
    }

    #[test]
    fn circos_colours() {
        let c = CircosColours::new(5);
        assert_eq!(c.resolve("dgrey").unwrap().to_hex(), "#737373");
        let a = c.resolve("black_a3").unwrap();
        assert!((a.a - 0.5).abs() < 1e-6);
        assert_eq!(c.resolve("224,7,20").unwrap().to_hex(), "#E00714");
        assert!((c.resolve("vlpurple_a5").unwrap().a - 1.0 / 6.0).abs() < 1e-6);
    }
}
