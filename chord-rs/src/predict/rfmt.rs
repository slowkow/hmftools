//! R's ways with numbers and strings that reach CHORD's prediction file:
//! - `write.table` writes each double with R's formatReal and EncodeReal0 at 15 significant digits
//!   (R_print.digits = DBL_DIG), fixed or scientific notation, whichever is narrower (scipen 0);
//! - `read.delim(..., row.names = 1)` converts the sample id column as type.convert does (logical,
//!   integer, double or else character), and the row names are then the values as.character gives
//!   (finding CH3);
//! - `set.seed` and `unif_rand` of R's default Mersenne-Twister, which `max.col` draws from to break
//!   ties (finding CH1).
//!
//! formatReal follows R's double-precision code path (R on arm64, where long double is double); on
//! x86-64 R scales in long double. The two agree on the values CHORD writes (vote fractions k/ntree and
//! their sums, which are within a few units in the last place of a short decimal).

/// tbl[k + 1] = 10^k as in format.c (tbl[0] = 1e-1), up to KP_MAX
const KP_MAX: i32 = 27;

fn tbl(i: i32) -> f64 {
    // format.c's table of exact powers of ten, indexed from 1e-1
    const T: [f64; 29] = [
        1e-1, 1e00, 1e01, 1e02, 1e03, 1e04, 1e05, 1e06, 1e07, 1e08, 1e09, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20,
        1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27,
    ];
    T[i as usize]
}

const R_DEC_MIN_EXPONENT: i32 = -308;

/// format.c scientific(): sign, decimal exponent, significant digits (at most `digits`) and whether
/// rounding to `digits` widens the number in fixed notation
fn scientific(x: f64, digits: i32) -> (bool, i32, i32, bool) {
    if x == 0.0 {
        return (false, 0, 1, false);
    }
    let (neg, r) = if x < 0.0 { (true, -x) } else { (false, x) };
    let mut kp = r.log10().floor() as i32 - digits + 1;
    let mut r_prec = r;
    if kp.abs() < 10 {
        if kp > 0 {
            r_prec /= tbl(kp + 1);
        } else if kp < 0 {
            r_prec *= tbl(-kp + 1);
        }
    } else if kp <= R_DEC_MIN_EXPONENT {
        r_prec = (r * 1e303) / 10f64.powf((kp + 303) as f64);
    } else {
        r_prec /= 10f64.powf(kp as f64);
    }
    if r_prec < tbl(digits) {
        r_prec *= 10.0;
        kp -= 1;
    }
    // nearbyint: round half to even
    let mut alpha = r_prec.round_ties_even();
    let mut nsig = digits;
    for _ in 1..=digits {
        alpha /= 10.0;
        if alpha == alpha.floor() {
            nsig -= 1;
        } else {
            break;
        }
    }
    if nsig == 0 && digits > 0 {
        nsig = 1;
        kp += 1;
    }
    let kpower = kp + digits - 1;
    let rgt = (digits - kpower).clamp(0, KP_MAX);
    let fuzz = 0.5 / tbl(1 + rgt);
    let rounding_widens = kpower > 0 && kpower <= KP_MAX && r < tbl(kpower + 1) - fuzz;
    (neg, kpower, nsig, rounding_widens)
}

/// formatReal for one value: width, digits after the point (or of the mantissa) and exponent mode
/// (0 fixed, 1 or 2 scientific)
fn format_real(x: f64, digits: i32) -> (i32, i32, i32) {
    if !x.is_finite() {
        return (if x.is_nan() { 3 } else if x > 0.0 { 3 } else { 4 }, 0, 0);
    }
    let (neg_i, kpower, nsig, rounding_widens) = scientific(x, digits);
    let mut left = kpower + 1;
    if rounding_widens {
        left -= 1;
    }
    let neg = neg_i as i32;
    let sleft = neg + if left <= 0 { 1 } else { left };
    let right = nsig - left;
    let mut rgt = right;
    let mxl = left;
    let mut mxsl = sleft;
    let mxns = nsig;
    if mxl < 0 {
        mxsl = 1 + neg;
    }
    if rgt < 0 {
        rgt = 0;
    }
    // the sign is in mxsl already
    let wf = mxsl + rgt + (rgt != 0) as i32;
    // a three-digit exponent from 1e100 and 1e-100 on (C's %e writes at least two digits)
    let mut e = if kpower >= 100 || kpower <= -100 { 2 } else { 1 };
    let mut d = mxns - 1;
    let mut w = neg + (d > 0) as i32 + d + 4 + e;
    // R_print.scipen is 0
    if wf <= w {
        e = 0;
        d = rgt;
        w = wf;
    }
    (w, d, e)
}

/// C's "%.<d>e": mantissa, 'e', sign and at least two exponent digits
fn c_exp(x: f64, d: usize, alt: bool) -> String {
    let s = format!("{:.*e}", d, x);
    let (mant, exp) = s.split_once('e').unwrap();
    let mut mant = mant.to_string();
    if alt && !mant.contains('.') {
        mant.push('.');
    }
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp),
    };
    format!("{mant}e{sign}{digits:0>2}")
}

/// EncodeReal0 with the width, digits and exponent mode of formatReal
fn encode_real(x: f64, w: i32, d: i32, e: i32) -> String {
    let x = if x == 0.0 { 0.0 } else { x };
    let s = if !x.is_finite() {
        if x.is_nan() {
            "NaN".to_string()
        } else if x > 0.0 {
            "Inf".to_string()
        } else {
            "-Inf".to_string()
        }
    } else if e != 0 {
        c_exp(x, d as usize, d != 0)
    } else {
        format!("{:.*}", d as usize, x)
    };
    format!("{:>width$}", s, width = w.max(0) as usize)
}

/// A double as write.table writes it (15 significant digits)
pub fn write_table_double(x: f64) -> String {
    let (w, d, e) = format_real(x, 15);
    encode_real(x, w, d, e)
}

/// as.character of a double, as far as the sample ids need it: R's 15 significant digits (R 4.3 and
/// later print integer-valued doubles of more than 15 digits in full; such ids are not reproduced)
pub fn as_character_double(x: f64) -> String {
    let (_, d, e) = format_real(x, 15);
    encode_real(x, 0, d, e)
}

/// The type type.convert gives a column of strings
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ColumnType {
    Logical,
    Integer,
    Double,
    Character,
}

/// R's Strtoi(s, 10): strtol in base 10 (leading white space, a sign), the whole string, within int
/// range and not NA_INTEGER
fn r_strtoi(s: &str) -> Option<i32> {
    let t = s.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let (neg, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let v: i64 = digits.parse::<i64>().ok()?;
    let v = if neg { -v } else { v };
    if v > i32::MAX as i64 || v <= i32::MIN as i64 {
        return None;
    }
    Some(v as i32)
}

/// R_strtod's decimal forms with only white space after them: [+-]digits[.digits][e[+-]digits],
/// "NaN", "Inf", "infinity" (any case); hexadecimal forms are not converted
fn r_strtod(s: &str) -> Option<f64> {
    let t = s.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']).trim_end_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let (neg, body) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let lower = body.to_ascii_lowercase();
    let v = if lower == "nan" {
        f64::NAN
    } else if lower == "inf" || lower == "infinity" {
        f64::INFINITY
    } else {
        let b = body.as_bytes();
        let mut i = 0;
        let mut mant_digits = 0;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            mant_digits += 1;
        }
        if i < b.len() && b[i] == b'.' {
            i += 1;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
                mant_digits += 1;
            }
        }
        if mant_digits == 0 {
            return None;
        }
        if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
            let mut j = i + 1;
            if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
                j += 1;
            }
            let start = j;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            if j > start {
                i = j;
            }
        }
        if i != b.len() {
            return None;
        }
        body[..i].parse::<f64>().ok()?
    };
    Some(if neg { -v } else { v })
}

/// read.delim's conversion of the sample id column: the row names as R holds them, or Err where R stops
/// (an id "NA", or two ids that convert to the same value)
pub fn r_row_names(ids: &[String]) -> Result<Vec<String>, String> {
    if ids.iter().any(|s| s == "NA") {
        return Err("missing values in 'row.names' are not allowed".to_string());
    }
    let is_logical = |s: &str| matches!(s, "F" | "T" | "FALSE" | "TRUE");
    let ty = if ids.iter().all(|s| is_logical(s)) {
        ColumnType::Logical
    } else if ids.iter().all(|s| r_strtoi(s).is_some()) {
        ColumnType::Integer
    } else if ids.iter().all(|s| r_strtod(s).is_some()) {
        ColumnType::Double
    } else {
        ColumnType::Character
    };
    let names: Vec<String> = ids
        .iter()
        .map(|s| match ty {
            ColumnType::Logical => if s.starts_with('T') { "TRUE" } else { "FALSE" }.to_string(),
            ColumnType::Integer => r_strtoi(s).unwrap().to_string(),
            ColumnType::Double => as_character_double(r_strtod(s).unwrap()),
            ColumnType::Character => s.clone(),
        })
        .collect();
    let mut seen = std::collections::HashSet::new();
    if names.iter().any(|n| !seen.insert(n.clone())) {
        return Err("duplicate 'row.names' are not allowed".to_string());
    }
    Ok(names)
}

/// R's Mersenne-Twister (RNG.c), as set.seed(seed) leaves it
pub struct RRng {
    mt: [u32; 624],
    mti: usize,
}

impl RRng {
    /// set.seed(seed): RNG_Init's scrambling, then FixupSeeds (mti = 624)
    pub fn set_seed(seed: u32) -> RRng {
        let mut s = seed;
        for _ in 0..50 {
            s = s.wrapping_mul(69069).wrapping_add(1);
        }
        let mut dummy = [0u32; 625];
        for d in dummy.iter_mut() {
            s = s.wrapping_mul(69069).wrapping_add(1);
            *d = s;
        }
        let mut mt = [0u32; 624];
        mt.copy_from_slice(&dummy[1..]);
        RRng { mt, mti: 624 }
    }

    fn genrand(&mut self) -> f64 {
        const N: usize = 624;
        const M: usize = 397;
        const MATRIX_A: u32 = 0x9908b0df;
        const UPPER_MASK: u32 = 0x80000000;
        const LOWER_MASK: u32 = 0x7fffffff;
        let mag01 = [0u32, MATRIX_A];
        if self.mti >= N {
            let mt = &mut self.mt;
            let mut kk = 0;
            while kk < N - M {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[(y & 1) as usize];
                kk += 1;
            }
            while kk < N - 1 {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ mag01[(y & 1) as usize];
                kk += 1;
            }
            let y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y as f64 * 2.3283064365386963e-10
    }

    /// unif_rand: in (0, 1)
    pub fn unif_rand(&mut self) -> f64 {
        const I2_32M1: f64 = 2.328306437080797e-10;
        let v = self.genrand();
        if v <= 0.0 {
            0.5 * I2_32M1
        } else if 1.0 - v <= 0.0 {
            1.0 - 0.5 * I2_32M1
        } else {
            v
        }
    }
}

/// Parses a C99 hexadecimal float as R's sprintf("%a") writes it (0x1.881723c1f2245p-4, 0x0p+0)
pub fn parse_hex_float(s: &str) -> Option<f64> {
    let (neg, t) = match s.strip_prefix('-') {
        Some(t) => (true, t),
        None => (false, s),
    };
    let t = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X"))?;
    let (mant, exp) = t.split_once(['p', 'P'])?;
    let exp: i32 = exp.parse().ok()?;
    let (int_part, frac_part) = mant.split_once('.').unwrap_or((mant, ""));
    if int_part.is_empty() || int_part.len() + frac_part.len() > 15 {
        return None;
    }
    let mut m: u64 = 0;
    for c in int_part.chars().chain(frac_part.chars()) {
        m = m * 16 + c.to_digit(16)? as u64;
    }
    let shift = exp - 4 * frac_part.len() as i32;
    if m >= (1u64 << 53) {
        return None;
    }
    // m * 2^shift, exact for normal doubles
    let v = (m as f64) * 2f64.powi(shift);
    Some(if neg { -v } else { v })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a hex float of any exponent, subnormals included (for the fixture)
    fn parse_hex_float_any(s: &str) -> Option<f64> {
        let (neg, t) = match s.strip_prefix('-') {
            Some(t) => (true, t),
            None => (false, s),
        };
        let t = t.strip_prefix("0x")?;
        let (mant, exp) = t.split_once('p')?;
        let exp: i32 = exp.parse().ok()?;
        let (int_part, frac_part) = mant.split_once('.').unwrap_or((mant, ""));
        let mut m: u64 = 0;
        for c in int_part.chars().chain(frac_part.chars()) {
            m = m * 16 + c.to_digit(16)? as u64;
        }
        // m * 2^shift in two steps, so that neither factor overflows or underflows on its own
        let shift = exp - 4 * frac_part.len() as i32;
        let half = shift / 2;
        let v = (m as f64) * 2f64.powi(half) * 2f64.powi(shift - half);
        Some(if neg { -v } else { v })
    }

    #[test]
    fn write_table_doubles() {
        let cases = [
            (0.0, "0"),
            (1.0, "1"),
            (0.004, "0.004"),
            (0.002, "0.002"),
            (0.006, "0.006"),
            (0.586, "0.586"),
            (0.5, "0.5"),
            (0.1 + 0.2, "0.3"),
            (1.0 / 3.0, "0.333333333333333"),
            (2.0 / 3.0, "0.666666666666667"),
            (1e5, "1e+05"),
            (123456.0, "123456"),
            (0.0001, "1e-04"),
            (0.00012, "0.00012"),
            (1e-300, "1e-300"),
            (123456789012.0, "123456789012"),
            (-0.25, "-0.25"),
            (1e15, "1e+15"),
            (9.9999999999999995, "10"),
        ];
        for (x, want) in cases {
            assert_eq!(write_table_double(x), want, "{x:e}");
        }
    }

    /// write.table's text for 1,742 doubles, written by R 4.5.2 (scripts/r_format_fixture.R): special
    /// values, every vote fraction k/500, sums of two, random magnitudes from 1e-320 to 1e308 and random
    /// decimals of up to 15 digits
    #[test]
    fn write_table_doubles_as_r() {
        let text = std::fs::read_to_string(format!("{}/testdata/r/write_table_doubles.tsv", env!("CARGO_MANIFEST_DIR"))).unwrap();
        let mut n = 0;
        for line in text.lines() {
            let (hex, want) = line.split_once('\t').unwrap();
            let x = parse_hex_float_any(hex).unwrap_or_else(|| panic!("{hex}"));
            assert_eq!(write_table_double(x), want, "{hex}");
            n += 1;
        }
        assert!(n > 1700);
    }

    /// Every vote pair of a 500-tree forest as R writes it (scripts/r_vote_pairs.R makes the file, given
    /// in CHORD_R_VOTE_PAIRS; skipped without it)
    #[test]
    fn vote_pairs_as_r() {
        let Ok(path) = std::env::var("CHORD_R_VOTE_PAIRS") else {
            eprintln!("CHORD_R_VOTE_PAIRS not set; skipped");
            return;
        };
        let text = std::fs::read_to_string(path).unwrap();
        let mut n = 0;
        for line in text.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            let (k1, k2): (u32, u32) = (f[0].parse().unwrap(), f[1].parse().unwrap());
            let (a, b) = (k1 as f64 / 500.0, k2 as f64 / 500.0);
            assert_eq!([write_table_double(a), write_table_double(b), write_table_double(a + b)], [f[2], f[3], f[4]], "{k1} {k2}");
            n += 1;
        }
        assert_eq!(n, 125751);
    }

    #[test]
    fn row_names_as_r_converts_them() {
        let one = |s: &str| r_row_names(&[s.to_string()]).unwrap()[0].clone();
        assert_eq!(one("HG008-T"), "HG008-T");
        assert_eq!(one("0123"), "123");
        assert_eq!(one("T"), "TRUE");
        assert_eq!(one("true"), "true");
        assert_eq!(one("1e5"), "1e+05");
        assert_eq!(one("1e3"), "1000");
        assert_eq!(one("100000"), "100000");
        assert_eq!(one("12.50"), "12.5");
        assert_eq!(one(" 12"), "12");
        assert_eq!(one("12 "), "12");
        assert_eq!(one("+5"), "5");
        assert_eq!(one("-2147483648"), "-2147483648");
        assert_eq!(one(".5"), "0.5");
        assert_eq!(one("5."), "5");
        assert_eq!(one("inf"), "Inf");
        assert_eq!(one("nan"), "NaN");
        assert_eq!(one("1d5"), "1d5");
        assert_eq!(one("TRUE "), "TRUE ");
        let many = |v: &[&str]| r_row_names(&v.iter().map(|s| s.to_string()).collect::<Vec<_>>());
        assert_eq!(many(&["007", "8.5"]).unwrap(), ["7", "8.5"]);
        assert_eq!(many(&["007", "A"]).unwrap(), ["007", "A"]);
        assert_eq!(many(&["T", "1"]).unwrap(), ["T", "1"]);
        assert!(many(&["NA"]).is_err());
        assert!(many(&["1", "01"]).is_err());
    }

    #[test]
    fn r_mersenne_twister() {
        // set.seed(1); sprintf("%.17g", runif(5)) in R
        let mut rng = RRng::set_seed(1);
        let want = [0.26550866314209998, 0.37212389963679016, 0.57285336335189641, 0.90820778999477625, 0.2016819310374558];
        for w in want {
            assert_eq!(rng.unif_rand(), w);
        }
    }

    #[test]
    fn hex_floats() {
        assert_eq!(parse_hex_float("0x0p+0"), Some(0.0));
        assert_eq!(parse_hex_float("0x1p+0"), Some(1.0));
        assert_eq!(parse_hex_float("0x1.8p+1"), Some(3.0));
        assert_eq!(parse_hex_float("0x1.3f8ed9737145cp-3"), Some(0x13f8ed9737145c_u64 as f64 * 2f64.powi(-55)));
        assert_eq!(parse_hex_float("-0x1p-2"), Some(-0.25));
        assert_eq!(parse_hex_float("0.5"), None);
    }
}
