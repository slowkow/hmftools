//! The scales package (1.4.0) as ggplot2's position scales use it: label functions (`number()`, `label_percent()`
//! with `precision()`'s accuracy, R's `format()` with `nsmall`, big marks and `drop0trailing`), break functions
//! (`breaks_width`, `breaks_log`, `extended_breaks` from [`rfmt`](crate::gg::rfmt)), the log ticks of
//! `guide_axis_logticks` (`minor_breaks_log`), and continuous range expansion.

use crate::gg::rfmt;

/// scales::precision: the accuracy for labelling a set of numbers (10^(floor(log10(smallest gap)) - 1), one
/// order coarser when every value's last digit would be 0, at most 1).
pub fn precision(xs: &[f64]) -> f64 {
    let mut x: Vec<f64> = xs.iter().copied().filter(|v| v.is_finite()).collect();
    x.sort_by(f64::total_cmp);
    x.dedup();
    if x.len() <= 1 {
        return 1.0;
    }
    let smallest = x.windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min);
    if smallest < f64::EPSILON.sqrt() {
        return 1.0;
    }
    let mut p = 10f64.powf((smallest.log10()).floor() - 1.0);
    if x.iter().all(|v| r_mod((v / p).round_ties_even_r(), 10.0) == 0.0) {
        p *= 10.0;
    }
    p.min(1.0)
}

/// R's round() for x / precision as `round()` (half to even on representable ties, as R >= 4.0)
trait RRound {
    fn round_ties_even_r(self) -> f64;
}

impl RRound for f64 {
    fn round_ties_even_r(self) -> f64 {
        rfmt::r_round(self, 0)
    }
}

fn r_mod(a: f64, b: f64) -> f64 {
    a - (a / b).floor() * b
}

/// R's `format(x, nsmall = n, scientific = FALSE, trim = TRUE, big.mark = mark)` of a numeric vector (7 significant
/// digits, one common number of decimals, at least nsmall), with `drop0trailing`.
pub fn format_fixed(xs: &[f64], nsmall: usize, big_mark: &str, drop0trailing: bool) -> Vec<String> {
    let mut rgt = 0usize;
    for &x in xs.iter().filter(|x| x.is_finite()) {
        if x == 0.0 {
            continue;
        }
        let s = format!("{:.6e}", x.abs());
        let (m, e) = s.split_once('e').unwrap();
        let e: i32 = e.parse().unwrap();
        let digits: String = m.chars().filter(|c| c.is_ascii_digit()).collect();
        let nsig = digits.trim_end_matches('0').len().max(1) as i32;
        rgt = rgt.max((nsig - e - 1).max(0) as usize);
    }
    let rgt = rgt.max(nsmall);
    xs.iter()
        .map(|&x| {
            if x.is_nan() {
                return "NaN".into();
            }
            if x.is_infinite() {
                return if x > 0.0 { "Inf".into() } else { "-Inf".into() };
            }
            let mut s = format!("{:.*}", rgt, x);
            if !big_mark.is_empty() {
                let (neg, body) = if let Some(b) = s.strip_prefix('-') { (true, b.to_string()) } else { (false, s.clone()) };
                let (int, frac) = match body.split_once('.') {
                    Some((i, f)) => (i.to_string(), Some(f.to_string())),
                    None => (body.clone(), None),
                };
                let mut grouped = String::new();
                for (k, ch) in int.chars().enumerate() {
                    if k > 0 && (int.len() - k) % 3 == 0 {
                        grouped.push_str(big_mark);
                    }
                    grouped.push(ch);
                }
                s = format!("{}{}{}", if neg { "-" } else { "" }, grouped, frac.map(|f| format!(".{f}")).unwrap_or_default());
            }
            if drop0trailing && s.contains('.') {
                // R's dropTrailing0: the zeros after the decimals, and the point itself
                let t = s.trim_end_matches('0');
                s = t.trim_end_matches('.').to_string();
            }
            s
        })
        .collect()
}

/// scales::number (and label_number / label_percent): rounded to the accuracy (precision of the scaled values
/// when not given), formatted with the decimals the accuracy needs, big mark " ", prefix and suffix, a hyphen for
/// negatives.
pub fn number(xs: &[f64], accuracy: Option<f64>, scale: f64, suffix: &str, drop0trailing: bool) -> Vec<String> {
    if xs.is_empty() {
        return Vec::new();
    }
    let scaled: Vec<f64> = xs.iter().map(|x| x * scale).collect();
    let acc = accuracy.unwrap_or_else(|| precision(&scaled));
    // round_any(x, accuracy / scale)
    let step = acc / scale;
    let rounded: Vec<f64> = xs.iter().map(|&x| rfmt::r_round(x / step, 0) * step).collect();
    let nsmall = (-(acc.log10().floor())).clamp(0.0, 20.0) as usize;
    let abs_scaled: Vec<f64> = rounded.iter().map(|x| scale * x.abs()).collect();
    let body = format_fixed(&abs_scaled, nsmall, " ", drop0trailing);
    rounded
        .iter()
        .zip(body)
        .map(|(&x, b)| {
            if x.is_nan() {
                return "NA".into();
            }
            let s = format!("{b}{suffix}");
            if x < 0.0 { format!("-{s}") } else { s }
        })
        .collect()
}

/// label_percent(accuracy, scale = 100, drop0trailing)
pub fn label_percent(xs: &[f64], accuracy: Option<f64>, scale: f64, drop0trailing: bool) -> Vec<String> {
    number(xs, accuracy, scale, "%", drop0trailing)
}

/// scales::breaks_width(width): fullseq over the range, from round_any(lo, width, floor) to round_any(hi, width,
/// ceiling).
pub fn breaks_width(range: (f64, f64), width: f64) -> Vec<f64> {
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let start = (lo / width).floor() * width;
    let end = (hi / width).ceil() * width;
    let n = ((end - start) / width + 1e-10).floor() as i64;
    (0..=n).map(|k| start + k as f64 * width).collect()
}

/// scales::breaks_log(n = 5, base = 10) (as `rfmt::log_breaks`)
pub fn breaks_log(range: (f64, f64)) -> Vec<f64> {
    rfmt::log_breaks(range.0, range.1, 5)
}

/// minor_breaks_log(smallest = 0.1)(limits) as guide_axis_logticks uses it: the ticks (data values) and their
/// type, 1 for the powers of ten (long), 2 for the fives (mid), 3 for the other digits (short).
pub fn log_ticks(limits: (f64, f64)) -> Vec<(f64, u8)> {
    let (mut lo, mut hi) = (limits.0.min(limits.1), limits.0.max(limits.1));
    let has_negatives = lo <= 0.0;
    let small: f64 = 0.1;
    if has_negatives {
        let large = lo.abs().max(hi.abs());
        let a = small * 10.0;
        lo = a.min(large);
        hi = a.max(large);
    }
    let start = lo.log10().floor() as i32 - 1;
    let end = hi.log10().ceil() as i32 + 1;
    let span = (end - start).abs();
    // findInterval(abs(end - start), c(8, 15), left.open = TRUE): detail 1, 5 or 10
    let detail = if span <= 8 { 1 } else if span <= 15 { 5 } else { 10 };
    let ladder: Vec<f64> = (start..=end).map(|e| 10f64.powi(e)).collect();
    let mut tens = ladder.clone();
    let mut fives: Vec<f64> = if detail <= 5 { ladder.iter().map(|v| 5.0 * v).collect() } else { Vec::new() };
    let mut ones: Vec<f64> = Vec::new();
    if detail == 1 {
        for v in &ladder {
            for d in 1..=9 {
                let o = d as f64 * v;
                if !tens.contains(&o) && !fives.contains(&o) {
                    ones.push(o);
                }
            }
        }
    }
    if has_negatives {
        tens.retain(|&t| t >= small);
        let neg: Vec<f64> = tens.iter().map(|t| -t).collect();
        tens.extend(neg);
        tens.push(0.0);
        fives.retain(|&t| t >= small);
        let neg: Vec<f64> = fives.iter().map(|t| -t).collect();
        fives.extend(neg);
        ones.retain(|&t| t >= small);
        let neg: Vec<f64> = ones.iter().map(|t| -t).collect();
        ones.extend(neg);
    }
    let mut out: Vec<(f64, u8)> = tens.into_iter().map(|t| (t, 1)).collect();
    out.extend(fives.into_iter().map(|t| (t, 2)));
    out.extend(ones.into_iter().map(|t| (t, 3)));
    out
}

/// A continuous range expanded by a multiplicative and an additive part on each side (expand_range4).
pub fn expand(range: (f64, f64), mult: f64, add: f64) -> (f64, f64) {
    let (lo, hi) = range;
    if hi == lo {
        // zero range: ggplot2 widens a single value by +-5% of its magnitude (0.5 for zero)
        let w = if lo == 0.0 { 1.0 } else { lo.abs() * 0.1 };
        return (lo - w / 2.0 - add, hi + w / 2.0 + add);
    }
    let d = hi - lo;
    (lo - d * mult - add, hi + d * mult + add)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scales_labels() {
        assert_eq!(label_percent(&[0.0, 0.25, 0.5, 0.75, 1.0], None, 100.0, false), vec!["0%", "25%", "50%", "75%", "100%"]);
        assert_eq!(label_percent(&[0.0, 0.025, 0.05], None, 100.0, false), vec!["0.0%", "2.5%", "5.0%"]);
        assert_eq!(label_percent(&[0.996717], Some(0.1), 100.0, false), vec!["99.7%"]);
        assert_eq!(label_percent(&[0.00001, 0.0001, 0.001, 0.01], None, 100.0, true), vec!["0.001%", "0.01%", "0.1%", "1%"]);
        assert_eq!(label_percent(&[30.0, 40.0, 50.0, 60.0], None, 1.0, false), vec!["30%", "40%", "50%", "60%"]);
        assert_eq!(format_fixed(&[1234567.0], 0, " ", false), vec!["1 234 567"]);
        assert_eq!(breaks_width((4.0, 15.0), 3.0), vec![3.0, 6.0, 9.0, 12.0, 15.0]);
        let t = log_ticks((1.0, 1000.0));
        assert!(t.contains(&(100.0, 1)) && t.contains(&(50.0, 2)) && t.contains(&(300.0, 3)));
    }
}
