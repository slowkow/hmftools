//! Circos's unit_parse for arithmetic of values with units ("0.526r-0.5p", "1r-50p"): each number with a
//! unit (b, u, p, r) becomes pixels or bases before the expression is evaluated. A relative value below 1 is
//! a fraction of the ideogram's inner radius, from 1 up of its outer radius. (Added for LINX's configurations.)

/// Replaces every `<number><unit>` term of `expr` by its value in pixels (r, p) or bases (u, b), as
/// Circos::Unit::unit_parse does (each match replaced once, left to right).
pub fn substitute(expr: &str, radius_inner: f64, radius_outer: f64, chromosomes_units: f64) -> String {
    let b = expr.as_bytes();
    let mut out = String::with_capacity(expr.len() + 16);
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        let starts_number = (c.is_ascii_digit() || c == b'.') && (i == 0 || !(b[i - 1].is_ascii_alphanumeric() || b[i - 1] == b'_'));
        if starts_number {
            let mut j = i;
            while j < b.len() && (b[j].is_ascii_digit() || b[j] == b'.') {
                j += 1;
            }
            if j < b.len() && matches!(b[j], b'b' | b'u' | b'p' | b'r') && (j + 1 >= b.len() || !(b[j + 1].is_ascii_alphanumeric() || b[j + 1] == b'_')) {
                if let Ok(v) = expr[i..j].parse::<f64>() {
                    let converted = match b[j] {
                        b'r' => v * if v < 1.0 { radius_inner } else { radius_outer },
                        b'u' => v * chromosomes_units,
                        _ => v,
                    };
                    out.push_str(&format!("{converted:?}"));
                    i = j + 1;
                    continue;
                }
            }
            out.push_str(&expr[i..j]);
            i = j;
            continue;
        }
        out.push(c as char);
        i += 1;
    }
    out
}

/// Whether `expr` has a term with a unit (so that it needs `substitute`).
pub fn has_unit_terms(expr: &str) -> bool {
    substitute(expr, 1.0, 1.0, 1.0) != expr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn terms_become_pixels() {
        assert_eq!(substitute("0.5r-2p", 100.0, 200.0, 1.0), "50.0-2.0");
        assert_eq!(substitute("1r-50p", 100.0, 200.0, 1.0), "200.0-50.0");
        assert_eq!(substitute("(dims(ideogram,radius_outer)+dims(ideogram,radius_inner))/2-15", 1.0, 2.0, 1.0), "(dims(ideogram,radius_outer)+dims(ideogram,radius_inner))/2-15");
        assert!(!has_unit_terms("var(value) < 0.001"));
        assert!(has_unit_terms("0.526r+1.5p"));
    }
}
