//! R's number formatting (`format`, `as.character`, `round`) and the break algorithms of the scales and
//! labeling packages (`extended_breaks`, `breaks_log`, regular minor breaks).

use crate::fmath::FMath;

/// R's `scientific()` for one number at 7 significant digits: (negative, exponent, significant digits).
fn sci_parts(x: f64, digits: usize) -> (bool, i32, usize) {
    if x == 0.0 {
        return (false, 0, 1);
    }
    let s = format!("{:.*e}", digits - 1, x.abs());
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let digits_str: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let nsig = digits_str.trim_end_matches('0').len().max(1);
    (x < 0.0, exp, nsig)
}

/// R's `format(x)` of a numeric vector with digits = 7: one common notation (fixed unless scientific is
/// narrower), one number of decimals. `trim`: no padding to a common width.
pub fn r_format(xs: &[f64], trim: bool) -> Vec<String> {
    let finite: Vec<f64> = xs.iter().copied().filter(|x| x.is_finite()).collect();
    if finite.is_empty() {
        return xs.iter().map(|_| "NA".to_string()).collect();
    }
    let digits = 7;
    let (mut neg, mut mxsl, mut rgt, mut mxe, mut mf) = (false, 1i32, 0i32, 1i32, 1usize);
    let mut min_e = i32::MAX;
    for &x in &finite {
        let (n, kp, nsig) = sci_parts(x, digits);
        neg |= n;
        let left = if kp >= 0 { kp + 1 } else { 1 };
        mxsl = mxsl.max(left);
        rgt = rgt.max((nsig as i32 - kp - 1).max(0));
        mxe = mxe.max(kp);
        min_e = min_e.min(kp);
        mf = mf.max(nsig);
    }
    let neg_w = i32::from(neg);
    let fixed_w = neg_w + mxsl + if rgt > 0 { rgt + 1 } else { 0 };
    let e_digits = if mxe.abs() >= 100 || min_e.abs() >= 100 { 3 } else { 2 };
    let sci_w = neg_w + if mf > 1 { mf as i32 + 1 } else { mf as i32 } + 2 + e_digits;
    let fixed = fixed_w <= sci_w;
    let mut out: Vec<String> = xs
        .iter()
        .map(|&x| {
            if !x.is_finite() {
                return "NA".to_string();
            }
            if fixed {
                format!("{:.*}", rgt as usize, x)
            } else {
                let s = format!("{:.*e}", mf - 1, x);
                let (m, e) = s.split_once('e').unwrap();
                let e: i32 = e.parse().unwrap();
                format!("{m}e{}{:0w$}", if e < 0 { '-' } else { '+' }, e.abs(), w = e_digits as usize)
            }
        })
        .collect();
    if !trim {
        let w = out.iter().map(|s| s.len()).max().unwrap_or(0);
        for s in &mut out {
            *s = format!("{s:>w$}");
        }
    }
    out
}

/// R's `as.character(x)` for a double: 15 significant digits, the shortest form ("1.46", "100", "1e+05").
pub fn r_as_character(x: f64) -> String {
    if x.is_nan() {
        return "NaN".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    if x == 0.0 {
        return "0".into();
    }
    let s = format!("{:.14e}", x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    let neg = m.starts_with('-');
    let digits: String = m.chars().filter(|c| c.is_ascii_digit()).collect();
    let digits = digits.trim_end_matches('0');
    let digits = if digits.is_empty() { "0" } else { digits };
    let nsig = digits.len() as i32;
    // fixed when it is no wider than scientific (R's rule, scipen 0)
    let fixed_w = if e >= 0 { (e + 1).max(nsig) + if nsig > e + 1 { 1 } else { 0 } } else { 1 + 1 + (-e - 1) + nsig };
    let sci_w = nsig + if nsig > 1 { 1 } else { 0 } + 4;
    let body = if fixed_w <= sci_w {
        if e >= 0 {
            let int_len = (e + 1) as usize;
            if digits.len() <= int_len {
                format!("{digits}{}", "0".repeat(int_len - digits.len()))
            } else {
                format!("{}.{}", &digits[..int_len], &digits[int_len..])
            }
        } else {
            format!("0.{}{digits}", "0".repeat((-e - 1) as usize))
        }
    } else {
        let mant = if digits.len() > 1 { format!("{}.{}", &digits[..1], &digits[1..]) } else { digits.to_string() };
        format!("{mant}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
    };
    if neg { format!("-{body}") } else { body }
}

/// R's `round(x, digits)`: IEC 60559 rounding of the decimal value (R >= 4.0 picks the closest of the two
/// candidates, half to even on exact ties).
pub fn r_round(x: f64, digits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let p = 10f64.fpowi(digits);
    let y = x * p;
    let lo = y.floor();
    let hi = y.ceil();
    let (a, b) = (lo / p, hi / p);
    let (da, db) = ((x - a).abs(), (b - x).abs());
    if da < db {
        a
    } else if db < da {
        b
    } else if (lo as i64) % 2 == 0 {
        a
    } else {
        b
    }
}

/// R's `%%` for doubles.
fn r_mod(a: f64, b: f64) -> f64 {
    a - (a / b).floor() * b
}

/// labeling::extended (Talbot, Lin and Hanrahan), as scales::extended_breaks calls it: Q = 1, 5, 2, 2.5, 4, 3,
/// weights 0.25, 0.2, 0.5, 0.05, not only loose.
pub fn extended_breaks(dmin: f64, dmax: f64, m: usize) -> Vec<f64> {
    if !(dmin.is_finite() && dmax.is_finite()) {
        return Vec::new();
    }
    if (dmax - dmin).abs() < 1e-10 * dmax.abs().max(1.0) {
        return vec![dmin];
    }
    let q_list = [1.0, 5.0, 2.0, 2.5, 4.0, 3.0];
    let w = [0.25, 0.2, 0.5, 0.05];
    let n = q_list.len() as f64;
    let m = m as f64;
    let eps = f64::EPSILON * 100.0;
    let mut best: Option<(f64, f64, f64)> = None;
    let mut best_score = -2.0;
    let simplicity = |i: usize, j: f64, lmin: f64, lmax: f64, lstep: f64| {
        let v = if (r_mod(lmin, lstep) < eps || lstep - r_mod(lmin, lstep) < eps) && lmin <= 0.0 && lmax >= 0.0 { 1.0 } else { 0.0 };
        1.0 - i as f64 / (n - 1.0) - j + v
    };
    let simplicity_max = |i: usize, j: f64| 1.0 - i as f64 / (n - 1.0) - j + 1.0;
    let coverage = |lmin: f64, lmax: f64| {
        let range = dmax - dmin;
        1.0 - 0.5 * ((dmax - lmax).fpowi(2) + (dmin - lmin).fpowi(2)) / (0.1 * range).fpowi(2)
    };
    let coverage_max = |span: f64| {
        let range = dmax - dmin;
        if span > range {
            let half = (span - range) / 2.0;
            1.0 - 0.5 * (2.0 * half * half) / (0.1 * range).fpowi(2)
        } else {
            1.0
        }
    };
    let density = |k: f64, lmin: f64, lmax: f64| {
        let r = (k - 1.0) / (lmax - lmin);
        let rt = (m - 1.0) / (lmax.max(dmax) - dmin.min(lmin));
        2.0 - (r / rt).max(rt / r)
    };
    let density_max = |k: f64| if k >= m { 2.0 - (k - 1.0) / (m - 1.0) } else { 1.0 };
    let mut j = 1.0f64;
    'outer: loop {
        for (i, &q) in q_list.iter().enumerate() {
            let sm = simplicity_max(i, j);
            if w[0] * sm + w[1] + w[2] + w[3] < best_score {
                break 'outer;
            }
            let mut k = 2.0f64;
            loop {
                let dm = density_max(k);
                if w[0] * sm + w[1] + w[2] * dm + w[3] < best_score {
                    break;
                }
                let delta = (dmax - dmin) / (k + 1.0) / j / q;
                let mut z = delta.flog10().ceil();
                loop {
                    let step = j * q * 10f64.fpow(z);
                    let cm = coverage_max(step * (k - 1.0));
                    if w[0] * sm + w[1] * cm + w[2] * dm + w[3] < best_score {
                        break;
                    }
                    let min_start = (dmax / step).floor() * j - (k - 1.0) * j;
                    let max_start = (dmin / step).ceil() * j;
                    if min_start > max_start {
                        z += 1.0;
                        continue;
                    }
                    let mut start = min_start;
                    while start <= max_start {
                        let lmin = start * (step / j);
                        let lmax = lmin + step * (k - 1.0);
                        let lstep = step;
                        let s = simplicity(i, j, lmin, lmax, lstep);
                        let c = coverage(lmin, lmax);
                        let g = density(k, lmin, lmax);
                        let score = w[0] * s + w[1] * c + w[2] * g + w[3];
                        if score > best_score {
                            best_score = score;
                            best = Some((lmin, lmax, lstep));
                        }
                        start += 1.0;
                    }
                    z += 1.0;
                }
                k += 1.0;
            }
        }
        j += 1.0;
    }
    let Some((lmin, lmax, lstep)) = best else { return Vec::new() };
    // seq(from, to, by): from + (0:n) * by, n = floor((to - from) / by + 1e-10)
    let count = ((lmax - lmin) / lstep + 1e-10).floor() as i64;
    (0..=count).map(|i| lmin + i as f64 * lstep).collect()
}

/// scales::breaks_log(n = 5, base = 10) for a positive range.
pub fn log_breaks(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let rng = (lo.flog10(), hi.flog10());
    let min = rng.0.floor();
    let max = rng.1.ceil();
    if max == min {
        return vec![10f64.fpow(min)];
    }
    let relevant = |b: &[f64]| b.iter().filter(|&&x| lo <= x && x <= hi).count();
    let seq = |by: f64| {
        let count = ((max - min) / by + 1e-10).floor() as i64;
        (0..=count).map(|i| 10f64.fpowi((min + i as f64 * by) as i32)).collect::<Vec<f64>>()
    };
    let mut by = ((max - min) / n as f64).floor() + 1.0;
    let breaks = seq(by);
    if relevant(&breaks) + 2 >= n {
        return breaks;
    }
    while by > 1.0 {
        by -= 1.0;
        let breaks = seq(by);
        if relevant(&breaks) + 2 >= n {
            return breaks;
        }
    }
    log_sub_breaks(rng, lo, hi, n)
}

/// scales:::log_sub_breaks (base 10): adds the multiples 2..9 one at a time, each the one that best evens
/// out the spacing on the log scale, until n - 2 breaks fall in the range.
fn log_sub_breaks(rng: (f64, f64), lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let (min, max) = (rng.0.floor(), rng.1.ceil());
    let mut steps: Vec<f64> = vec![1.0];
    let mut candidates: Vec<f64> = (2..10).map(|x| x as f64).collect();
    let relevant = |b: &[f64]| b.iter().filter(|&&x| lo <= x && x <= hi).count();
    let mut breaks: Vec<f64> = Vec::new();
    while !candidates.is_empty() {
        let delta = |x: f64| {
            let mut v: Vec<f64> = steps.iter().copied().chain([x, 10.0]).collect();
            v.sort_by(f64::total_cmp);
            v.windows(2).map(|w| w[1].flog10() - w[0].flog10()).fold(f64::INFINITY, f64::min)
        };
        // which.max: the first of the largest
        let mut best = 0;
        for i in 1..candidates.len() {
            if delta(candidates[i]) > delta(candidates[best]) {
                best = i;
            }
        }
        steps.push(candidates.remove(best));
        // as.vector(outer(base^seq(min, max), steps))
        breaks.clear();
        for s in &steps {
            let mut e = min;
            while e <= max {
                breaks.push(10f64.fpowi(e as i32) * s);
                e += 1.0;
            }
        }
        if relevant(&breaks) + 2 >= n {
            break;
        }
    }
    if relevant(&breaks) + 2 >= n {
        breaks.sort_by(f64::total_cmp);
        let first = breaks.iter().position(|&x| lo <= x).unwrap_or(0);
        let last = breaks.iter().rposition(|&x| x <= hi).unwrap_or(breaks.len() - 1);
        let (a, b) = (first.saturating_sub(1), (last + 1).min(breaks.len() - 1));
        breaks[a..=b].to_vec()
    } else {
        extended_breaks(lo, hi, n)
    }
}

/// scales::regular_minor_breaks()(b, limits, n = 2): the major breaks (extended by one step towards the
/// limits) and the midpoints between them; ggplot2 then keeps those in the range that are not major.
pub fn regular_minor_breaks(major: &[f64], limits: (f64, f64)) -> Vec<f64> {
    let mut b: Vec<f64> = major.iter().copied().filter(|x| x.is_finite()).collect();
    if b.len() < 2 {
        return Vec::new();
    }
    let bd = b[1] - b[0];
    if limits.0 < b[0] {
        b.insert(0, b[0] - bd);
    }
    if limits.1 > *b.last().unwrap() {
        b.push(b.last().unwrap() + bd);
    }
    let mut out = Vec::new();
    for w in b.windows(2) {
        out.push(w[0]);
        out.push(w[0] + (w[1] - w[0]) / 2.0);
    }
    out.push(*b.last().unwrap());
    out.into_iter().filter(|x| *x >= limits.0 && *x <= limits.1 && !major.contains(x)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_format_vectors() {
        // as R 4.5 prints them
        assert_eq!(r_format(&[0.0, 2e5, 4e5, 6e5], true), ["0e+00", "2e+05", "4e+05", "6e+05"]);
        assert_eq!(r_format(&[0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0], true), ["0", "500", "1000", "1500", "2000", "2500"]);
        assert_eq!(r_format(&[0.25, 0.5, 1.0], true), ["0.25", "0.50", "1.00"]);
        assert_eq!(r_format(&[100.0, 10000.0, 1e6], false), ["1e+02", "1e+04", "1e+06"]);
        assert_eq!(r_format(&[0.0, 1.0, 2.0, 3.0], true), ["0", "1", "2", "3"]);
    }

    #[test]
    fn as_character() {
        assert_eq!(r_as_character(1.46), "1.46");
        assert_eq!(r_as_character(100.0), "100");
        assert_eq!(r_as_character(0.57 * 100.0), "57");
        assert_eq!(r_as_character(0.98), "0.98");
        assert_eq!(r_as_character(1e5), "1e+05");
        assert_eq!(r_as_character(123456.0), "123456");
    }

    #[test]
    fn rounding() {
        assert_eq!(r_round(1.455, 2), 1.46);
        assert_eq!(r_round(2.5, 0), 2.0);
        assert_eq!(r_round(1.4600, 2), 1.46);
    }

    #[test]
    fn extended_breaks_as_r() {
        // scales::extended_breaks()(c(...)) in R 4.5
        assert_eq!(extended_breaks(0.0, 575789.0, 5), [0.0, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5]);
        assert_eq!(extended_breaks(-28789.45, 604578.4, 5), [0.0, 2e5, 4e5, 6e5]);
        assert_eq!(extended_breaks(0.0, 2478.0, 5), [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0]);
        assert_eq!(extended_breaks(1.0, 8.4, 5), [0.0, 2.0, 4.0, 6.0, 8.0]);
        assert_eq!(extended_breaks(0.0, 3.0, 5), [0.0, 1.0, 2.0, 3.0]);
        assert_eq!(extended_breaks(0.0, 9200.0, 5), [0.0, 2500.0, 5000.0, 7500.0, 10000.0]);
        assert_eq!(extended_breaks(-0.1, 3.1, 5), [0.0, 1.0, 2.0, 3.0]);
        let b = extended_breaks(0.0001, 0.37, 5);
        assert_eq!(b.len(), 5);
        assert!((b[4] - 0.4).abs() < 1e-12);
    }

    #[test]
    fn log_breaks_as_scales() {
        // scales::breaks_log()(...) in R 4.5 / scales 1.4
        assert_eq!(log_breaks(10.0, 3.4e6, 5), [10.0, 1e3, 1e5, 1e7]);
        assert_eq!(log_breaks(1.0, 3.4e6, 5), [1.0, 100.0, 1e4, 1e6]);
        assert_eq!(log_breaks(20.0, 2e6, 5), [10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]);
        assert_eq!(log_breaks(2.0, 50.0, 5), [1.0, 3.0, 10.0, 30.0, 100.0]);
        // minor breaks of majors 0, 2, 4, 6 (log10) over 0.5..6.5: the odd powers
        assert_eq!(regular_minor_breaks(&[0.0, 2.0, 4.0, 6.0], (0.5, 6.5)), [1.0, 3.0, 5.0]);
    }
}
