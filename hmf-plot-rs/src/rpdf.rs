//! Text sizes as R's pdf() device reports them for the family "sans": Helvetica's AFM metrics at the font
//! size rounded to whole points (R's PDF device measures at integer sizes), with kerning. cowplot turns
//! plots into grobs on a pdf(NULL) device, so the ggplot2 layouts it composes (LINX's fusion plots) are
//! sized with these metrics, while the image is then drawn with Cairo and DejaVu Sans. (Added for LINX.)

use crate::rpdf_helvetica::{KERNING, WIDTHS};

fn width_units(s: &str) -> i64 {
    let b: Vec<u8> = s.bytes().collect();
    let mut w: i64 = 0;
    for (i, &c) in b.iter().enumerate() {
        w += if (32..=126).contains(&c) { WIDTHS[(c - 32) as usize] as i64 } else { 556 };
        if let Some(&n) = b.get(i + 1) {
            if let Some(k) = KERNING.iter().find(|p| p.0 == c && p.1 == n) {
                w += k.2 as i64;
            }
        }
    }
    w
}

/// The width of `s` in big points (1/72 inch) at `fontsize` points.
pub fn string_width_bigpt(s: &str, fontsize: f64) -> f64 {
    width_units(s) as f64 * fontsize.round() / 1000.0
}

/// The height of a text row (a label's grob height) in big points at `fontsize` points: 0.938 em of the
/// rounded size, as R 4.5's pdf device gives it.
pub fn text_height_bigpt(fontsize: f64) -> f64 {
    0.938 * fontsize.round()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn as_r_measures() {
        // convertWidth(grobWidth(textGrob(s, gp = gpar(fontsize = 11.52))), "pt") on pdf(NULL) in R 4.5
        let pt = |bigpt: f64| bigpt * 72.27 / 72.0;
        assert!((pt(string_width_bigpt("3-Prime UTR", 11.52)) - 73.462455).abs() < 1e-5);
        assert!((pt(string_width_bigpt("ANK_REPEAT", 11.52)) - 78.208185).abs() < 1e-5);
        assert!((pt(string_width_bigpt("ANK_REP_REGION", 11.52)) - 110.440605).abs() < 1e-5);
        assert!((pt(string_width_bigpt("ZINC_FINGER_C2H2_2", 11.52)) - 131.182095).abs() < 1e-5);
        assert!((pt(text_height_bigpt(11.52)) - 11.29821).abs() < 1e-4);
        assert!((pt(text_height_bigpt(5.0)) - 4.7075875).abs() < 1e-6);
    }
}
