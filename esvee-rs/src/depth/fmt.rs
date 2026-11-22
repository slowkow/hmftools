//! Number formatting that reproduces java.util.Formatter and htsjdk's VCF encoder.

/// Shortest round-trip decimal digits of a positive finite double: (digits, exponent) where the
/// value is 0.d1d2d3... x 10^exponent.
fn digits(d: f64) -> (Vec<u8>, i32) {
    let s = format!("{:e}", d);
    let (mant, exp) = s.split_once('e').unwrap();
    let digits: Vec<u8> = mant.bytes().filter(|b| b.is_ascii_digit()).map(|b| b - b'0').collect();
    (digits, exp.parse::<i32>().unwrap() + 1)
}

/// Round the digit string half-up so that `keep` digits remain; may carry into a new leading digit.
fn round_half_up(digits: &mut Vec<u8>, exp: &mut i32, keep: i32) {
    if keep < 0 {
        digits.clear();
        return;
    }
    let keep = keep as usize;
    if digits.len() <= keep {
        return;
    }
    let round_up = digits[keep] >= 5;
    digits.truncate(keep);
    if round_up {
        let mut i = keep;
        loop {
            if i == 0 {
                digits.insert(0, 1);
                *exp += 1;
                break;
            }
            i -= 1;
            if digits[i] == 9 {
                digits[i] = 0;
            } else {
                digits[i] += 1;
                break;
            }
        }
    }
}

/// String.format("%.<prec>f", d)
pub fn java_fixed(d: f64, prec: usize) -> String {
    if d.is_nan() {
        return "NaN".into();
    }
    if d.is_infinite() {
        return if d > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    let neg = d < 0.0 || (d == 0.0 && d.is_sign_negative());
    let a = d.abs();
    let mut int_part: Vec<u8> = vec![0];
    let mut frac: Vec<u8> = vec![0; prec];
    if a != 0.0 {
        let (mut ds, mut exp) = digits(a);
        let keep = exp + prec as i32;
        round_half_up(&mut ds, &mut exp, keep);
        // value = 0.ds x 10^exp
        let mut all: Vec<u8> = Vec::new();
        if exp <= 0 {
            all.push(0);
            all.extend(std::iter::repeat(0).take((-exp) as usize));
            all.extend(&ds);
            int_part = vec![0];
            frac = all[1..].to_vec();
        } else {
            let e = exp as usize;
            let mut ds2 = ds.clone();
            while ds2.len() < e {
                ds2.push(0);
            }
            int_part = ds2[..e].to_vec();
            frac = ds2[e..].to_vec();
        }
        frac.resize(prec, 0);
    }
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    s.extend(int_part.iter().map(|d| (d + b'0') as char));
    if prec > 0 {
        s.push('.');
        s.extend(frac.iter().map(|d| (d + b'0') as char));
    }
    s
}

/// String.format("%.<prec>e", d)
pub fn java_sci(d: f64, prec: usize) -> String {
    if d.is_nan() {
        return "NaN".into();
    }
    if d.is_infinite() {
        return if d > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    let neg = d < 0.0;
    let a = d.abs();
    let (mut ds, mut exp) = if a == 0.0 { (vec![0], 1) } else { digits(a) };
    round_half_up(&mut ds, &mut exp, prec as i32 + 1);
    ds.resize(prec + 1, 0);
    let e = exp - 1;
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    s.push((ds[0] + b'0') as char);
    if prec > 0 {
        s.push('.');
        s.extend(ds[1..].iter().map(|d| (d + b'0') as char));
    }
    s.push('e');
    s.push(if e < 0 { '-' } else { '+' });
    s.push_str(&format!("{:02}", e.abs()));
    s
}

/// htsjdk VCFEncoder.formatVCFDouble
pub fn vcf_double(d: f64) -> String {
    if d < 1.0 {
        if d < 0.01 {
            if d.abs() >= 1e-20 {
                java_sci(d, 3)
            } else {
                "0.00".into()
            }
        } else {
            java_fixed(d, 3)
        }
    } else {
        java_fixed(d, 2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn formats() {
        assert_eq!(java_fixed(0.125, 2), "0.13");
        assert_eq!(java_fixed(0.9996, 3), "1.000");
        assert_eq!(java_fixed(0.0, 3), "0.000");
        assert_eq!(java_sci(0.009668, 3), "9.668e-03");
        assert_eq!(vcf_double(0.0), "0.00");
        assert_eq!(vcf_double(0.192), "0.192");
        assert_eq!(vcf_double(1.0), "1.00");
        assert_eq!(vcf_double(f64::NAN), "NaN");
    }
}
