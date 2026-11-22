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
}
