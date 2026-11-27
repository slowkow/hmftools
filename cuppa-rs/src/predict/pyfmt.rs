//! Numbers and text as Python and pandas write them: `repr(float)` (which is also numpy's `str` of a
//! float64 and what `DataFrame.to_csv` writes for float columns without a float_format), and the
//! `%(...)s` log format of Python's logging module.

/// The shortest round-trip digits of a positive finite double and the position of the decimal point
/// (value = 0.d1d2... x 10^decpt, as _Py_dg_dtoa's decpt)
fn shortest(a: f64) -> (String, i32) {
    let s = format!("{:e}", a);
    let (mant, exp) = s.split_once('e').unwrap();
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let digits = digits.trim_end_matches('0').to_string();
    let digits = if digits.is_empty() { "0".to_string() } else { digits };
    (digits, exp.parse::<i32>().unwrap() + 1)
}

/// repr(x) for a Python float (float_repr_style 'short'): exponential if decpt <= -4 or decpt > 16
pub fn repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let sign = if x.is_sign_negative() { "-" } else { "" };
    if x == 0.0 {
        return format!("{sign}0.0");
    }
    let (digits, decpt) = shortest(x.abs());
    let nd = digits.len() as i32;
    if decpt <= -4 || decpt > 16 {
        let exp = decpt - 1;
        let mant = if nd == 1 { digits.clone() } else { format!("{}.{}", &digits[..1], &digits[1..]) };
        let es = if exp < 0 { format!("-{:02}", -exp) } else { format!("+{:02}", exp) };
        format!("{sign}{mant}e{es}")
    } else if decpt <= 0 {
        format!("{sign}0.{}{}", "0".repeat((-decpt) as usize), digits)
    } else if decpt >= nd {
        format!("{sign}{}{}.0", digits, "0".repeat((decpt - nd) as usize))
    } else {
        format!("{sign}{}.{}", &digits[..decpt as usize], &digits[decpt as usize..])
    }
}

/// a float cell of DataFrame.to_csv: repr, NaN as the empty na_rep
pub fn csv_float(x: f64) -> String {
    if x.is_nan() { String::new() } else { repr(x) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_repr() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (0.1, "0.1"),
            (128.0, "128.0"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (1.5e-5, "1.5e-05"),
            (1e16, "1e+16"),
            (1234567890123456.0, "1234567890123456.0"),
            (1e22, "1e+22"),
            (1e-300, "1e-300"),
            (2.5e300, "2.5e+300"),
            (0.39593293399613794, "0.39593293399613794"),
            (9227.0, "9227.0"),
            (0.4444444444444444, "0.4444444444444444"),
            (f64::NAN, "nan"),
            (f64::INFINITY, "inf"),
        ];
        for (x, s) in cases {
            assert_eq!(repr(*x), *s, "{x:e}");
        }
    }
}
