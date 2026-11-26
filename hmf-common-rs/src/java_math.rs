//! Java numeric library behaviour that differs from the obvious Rust equivalent.

/// java.util.stream.DoubleStream.sum() (JDK 17 and later): Kahan compensated summation, falling back to the
/// simple sum when the compensated result is NaN because of infinities.
pub fn double_stream_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0f64;
    let mut compensation = 0.0f64;
    let mut simple_sum = 0.0f64;
    for value in values {
        let tmp = value - compensation;
        let velvel = sum + tmp;
        compensation = (velvel - sum) - tmp;
        sum = velvel;
        simple_sum += value;
    }
    let tmp = sum - compensation;
    if tmp.is_nan() && simple_sum.is_infinite() { simple_sum } else { tmp }
}

/// Math.max(double, double): NaN if either is NaN, and 0.0 is larger than -0.0.
pub fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && a.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a >= b { a } else { b }
}

/// Double.doubleToLongBits: the bits, with every NaN folded to the canonical NaN.
pub fn double_to_long_bits(d: f64) -> i64 {
    if d.is_nan() { 0x7ff8_0000_0000_0000 } else { d.to_bits() as i64 }
}

/// Double.compare: numeric order, with -0.0 before 0.0 and NaN (any NaN) after everything.
pub fn double_compare(a: f64, b: f64) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    if a < b {
        Less
    } else if a > b {
        Greater
    } else {
        double_to_long_bits(a).cmp(&double_to_long_bits(b))
    }
}

/// Arrays.sort(double[]): ascending in the order of Double.compare.
pub fn sort_doubles(values: &mut [f64]) {
    values.sort_by(|a, b| double_compare(*a, *b));
}

/// Math.round(double): the closest long, ties towards positive infinity; NaN gives 0 and values beyond the
/// long range are clamped.
pub fn java_round(d: f64) -> i64 {
    if d.is_nan() {
        return 0;
    }
    let f = d.floor();
    // exact: the fraction of a double is representable
    let r = if d - f >= 0.5 { f + 1.0 } else { f };
    r as i64
}

#[cfg(test)]
mod tests {
    use super::double_stream_sum;

    #[test]
    fn compensated_sum_differs_from_naive_sum() {
        // ten 0.1s: the naive sum is 0.9999999999999999, Java's DoubleStream.sum() gives 1.0
        let values = [0.1f64; 10];
        assert_eq!(1.0, double_stream_sum(values));
        assert_ne!(1.0, values.iter().sum::<f64>());
        assert_eq!(0.0, double_stream_sum(std::iter::empty()));
        assert_eq!(f64::INFINITY, double_stream_sum([f64::INFINITY, 1.0]));
    }

    #[test]
    fn max_compare_round() {
        use super::*;
        assert!(java_max(f64::NAN, 1.0).is_nan() && java_max(1.0, f64::NAN).is_nan());
        assert!(java_max(-0.0, 0.0).is_sign_positive() && java_max(0.0, -0.0).is_sign_positive());
        assert_eq!(java_max(2.0, 3.0), 3.0);
        let mut v = [f64::NAN, 1.0, 0.0, -0.0, f64::NEG_INFINITY, f64::INFINITY, -2.5];
        sort_doubles(&mut v);
        assert_eq!(&format!("{v:?}"), "[-inf, -2.5, -0.0, 0.0, 1.0, inf, NaN]");
        assert_eq!(java_round(2.5), 3);
        assert_eq!(java_round(-2.5), -2);
        assert_eq!(java_round(0.49999999999999994), 0);
        assert_eq!(java_round(f64::NAN), 0);
        assert_eq!(java_round(1e300), i64::MAX);
    }
}
