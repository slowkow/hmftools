//! hmf-common `segmentation.copynumber`: the least-cost piecewise constant fit of a series of values (the
//! dynamic programme of the R copynumber package's pcf, Nilsen et al. 2012), its segment penalty
//! (`Gamma`), R's running median (`Runmed`, `WindowedMedian`) and the median absolute deviation (`Stats`).
//!
//! The floating point operations are Java's, in Java's order: the cumulative sums, the cost of every
//! candidate segment, the first least cost winning ties, `Arrays.sort` and `Math.max` semantics for
//! -0.0 and NaN, and the two heaps of `WindowedMedian` with `PriorityQueue.remove`'s equality.

use crate::java_math::{double_compare, double_to_long_bits, java_max, java_round, sort_doubles};
use anyhow::{bail, Result};

/// Doubles.mean: the plain sum divided by the count; NaN when empty.
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sum = 0.0;
    for v in values {
        sum += v;
    }
    sum / values.len() as f64
}

/// Doubles.median(double[]): of a sorted copy; the mean of the two middle values when the count is even.
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sort_doubles(&mut sorted);
    let n = sorted.len();
    if n.is_multiple_of(2) { (sorted[n / 2] + sorted[n / 2 - 1]) / 2.0 } else { sorted[n / 2] }
}

/// Stats.medianAbsoluteDeviation: 1.4826 times the median absolute deviation from the median.
pub fn median_absolute_deviation(values: &[f64]) -> f64 {
    let m = median(values);
    let deviations: Vec<f64> = values.iter().map(|v| (v - m).abs()).collect();
    median(&deviations) * 1.4826
}

/// The segments of a fit: their lengths, first indices and means rounded to three decimals.
#[derive(Clone, Debug, PartialEq)]
pub struct PiecewiseConstantFit {
    pub lengths: Vec<usize>,
    pub start_positions: Vec<usize>,
    pub means: Vec<f64>,
}

/// Gamma: the segment penalty for a series. Normalised, it is `gamma` times the square of the median
/// absolute deviation of the values from their running median (window up to 51), or `gamma / 100` when
/// that deviation is 0. Java throws on an empty series.
pub fn segment_penalty(y: &[f64], gamma: f64, normalise: bool) -> Result<f64> {
    if y.is_empty() {
        bail!("Input array must not be empty");
    }
    if !normalise {
        return Ok(gamma);
    }
    let running_medians = runmed(y, filter_width(y.len() as i32), true)?;
    let diffs: Vec<f64> = y.iter().zip(&running_medians).map(|(a, b)| a - b).collect();
    let sd = median_absolute_deviation(&diffs);
    Ok(if sd == 0.0 { 0.01 * gamma } else { sd * sd * gamma })
}

/// Gamma.filterWidth
pub fn filter_width(n: i32) -> i32 {
    if n >= 51 { 51 } else if n % 2 == 0 { n - 1 } else { n }
}

/// Runmed.runmed: R's running median of odd window `k`, the ends smoothed as R's smoothEnds.
pub fn runmed(data: &[f64], k: i32, smooth: bool) -> Result<Vec<f64>> {
    if k < 1 {
        bail!("Window size k must be >= 1");
    }
    if k % 2 != 1 {
        bail!("Window size k must be odd");
    }
    if data.len() < 2 {
        return Ok(data.to_vec());
    }
    let result = windowed_median(data, k as usize)?;
    if smooth { smooth_ends(&result, k) } else { Ok(result) }
}

/// Runmed.smoothEnds
pub fn smooth_ends(y: &[f64], k: i32) -> Result<Vec<f64>> {
    if k < 0 || k % 2 == 0 {
        bail!("bandwidth 'k' must be >= 1 and odd!");
    }
    let half_k = (k / 2) as usize;
    if half_k < 1 {
        return Ok(y.to_vec());
    }
    let n = y.len();
    let mut sm = y.to_vec();
    if half_k >= 2 {
        sm[1] = med3(y[0], y[1], y[2]);
        sm[n - 2] = med3(y[n - 1], y[n - 2], y[n - 3]);
        if half_k >= 3 {
            for i in 3..=half_k {
                let j = 2 * i - 1;
                sm[i - 1] = median_odd(&y[0..j])?;
                sm[n - i] = median_odd(&y[n - j..n])?;
            }
        }
    }
    sm[0] = med3(y[0], sm[1], 3.0 * sm[1] - 2.0 * sm[2]);
    sm[n - 1] = med3(y[n - 1], sm[n - 2], 3.0 * sm[n - 2] - 2.0 * sm[n - 3]);
    Ok(sm)
}

/// Runmed.med3
pub fn med3(a: f64, b: f64, c: f64) -> f64 {
    if a < b {
        if c < b { java_max(a, c) } else { b }
    } else if c < a {
        java_max(c, b)
    } else {
        a
    }
}

/// Runmed.medianOdd
pub fn median_odd(x: &[f64]) -> Result<f64> {
    let n = x.len();
    if n % 2 != 1 {
        bail!("n must be odd for medianOdd");
    }
    let mut sorted = x.to_vec();
    sort_doubles(&mut sorted);
    Ok(sorted[n.div_ceil(2) - 1])
}

/// One of WindowedMedian's PriorityQueue<Double>s, kept as a vector sorted in Double.compare order.
/// Only the multiset of values matters: peek and poll give the largest (max heap) or smallest (min heap)
/// value in that order, and remove(Object) takes out one value that Double.equals the argument.
#[derive(Default)]
struct SortedHeap {
    values: Vec<f64>,
}

impl SortedHeap {
    fn len(&self) -> usize {
        self.values.len()
    }
    fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
    fn add(&mut self, v: f64) {
        let at = self.values.partition_point(|x| double_compare(*x, v) != std::cmp::Ordering::Greater);
        self.values.insert(at, v);
    }
    fn last(&self) -> f64 {
        self.values[self.values.len() - 1]
    }
    fn remove_first(&mut self) -> f64 {
        self.values.remove(0)
    }
    fn remove_last(&mut self) -> f64 {
        self.values.pop().unwrap()
    }
    /// PriorityQueue.remove(Object): one element equal by Double.equals, if there is one
    fn remove(&mut self, v: f64) {
        let key = double_to_long_bits(v);
        let at = self.values.partition_point(|x| double_compare(*x, v) == std::cmp::Ordering::Less);
        if at < self.values.len() && double_to_long_bits(self.values[at]) == key {
            self.values.remove(at);
        }
    }
}

/// WindowedMedian: running medians of an odd window by a max heap of the lower half and a min heap of the
/// upper half; the first and last `window / 2` values are copied from the data.
pub fn windowed_median(data: &[f64], window_size: usize) -> Result<Vec<f64>> {
    if window_size % 2 != 1 {
        bail!("Window size must be odd");
    }
    let n = data.len();
    let half = window_size / 2;
    // SG3: Java's heaps lose a 0.0 or -0.0 that compares equal to its other zero but is not equal to it;
    // corrected, both are taken as 0.0 (a zero median may then lose its sign)
    let fold_zeros = !crate::compat::java("SG3") && data.iter().any(|v| *v == 0.0 && v.is_sign_negative());
    let folded: Vec<f64>;
    let values: &[f64] = if fold_zeros {
        folded = data.iter().map(|v| if *v == 0.0 { 0.0 } else { *v }).collect();
        &folded
    } else {
        data
    };
    let mut lower = SortedHeap::default(); // max heap: its peek is the last value
    let mut upper = SortedHeap::default(); // min heap: its peek is the first value
    let mut result = vec![0.0; n];

    fn balance(lower: &mut SortedHeap, upper: &mut SortedHeap) {
        if lower.len() > upper.len() + 1 {
            let v = lower.remove_last();
            upper.add(v);
        } else if lower.len() < upper.len() {
            let v = upper.remove_first();
            lower.add(v);
        }
    }
    fn add(lower: &mut SortedHeap, upper: &mut SortedHeap, v: f64) {
        if lower.is_empty() || v <= lower.last() {
            lower.add(v);
        } else {
            upper.add(v);
        }
        balance(lower, upper);
    }

    for &v in values.iter().take(window_size) {
        add(&mut lower, &mut upper, v);
    }
    if half < n {
        result[half] = lower.last();
    }
    for i in window_size..n {
        let old = values[i - window_size];
        // Java's maxHeap.peek() cannot be null here: the window is never empty
        if old <= lower.last() {
            lower.remove(old);
        } else {
            upper.remove(old);
        }
        balance(&mut lower, &mut upper);
        add(&mut lower, &mut upper, values[i]);
        result[i - half] = lower.last();
    }
    let m = half.min(n);
    result[..m].copy_from_slice(&data[..m]);
    if half > n {
        bail!("arraycopy: last index {} out of bounds for length {n}", n as i64 - half as i64);
    }
    result[n - half..].copy_from_slice(&data[n - half..]);
    Ok(result)
}

/// The least-cost segmentation of `y` for a segment penalty, as segment lengths (Segmenter).
///
/// For every end index the cost of each start is `leastCost[start] + penalty + segmentCost(start, end)`
/// with the segment cost `sum of squares - square of sum / length` taken from cumulative sums; the first
/// start with the smallest cost wins. When no cost is below Double.MAX_VALUE (only with infinite or NaN
/// values) Java records 0 as the previous segment end; if that happens for the first value its
/// backtracking never ends, which is an error here.
pub fn least_cost_segment_lengths(y: &[f64], penalty: f64) -> Result<Vec<usize>> {
    let n = y.len();
    if n == 0 {
        // a single empty segment, as Java's singletonList(y)
        return Ok(vec![0]);
    }
    // Java's cumulative sums, with a leading 0 so that start 0 needs no branch: x - 0.0 is x, bit for bit
    let mut sums = vec![0.0f64; n + 1];
    let mut squares = vec![0.0f64; n + 1];
    sums[1] = y[0];
    squares[1] = y[0] * y[0];
    for i in 1..n {
        sums[i + 1] = sums[i] + y[i];
        squares[i + 1] = squares[i] + (y[i] * y[i]);
    }
    // leastCostEndingJustBefore[start] + segmentPenalty, which Java adds first
    let mut least_plus_penalty = vec![0.0f64; n + 1];
    least_plus_penalty[0] = 0.0 + penalty;
    let mut previous_end = vec![0i64; n];
    for end in 0..n {
        let (min_cost, start) = cheapest_start(&sums[..=end], &squares[..=end], &least_plus_penalty[..=end], sums[end + 1], squares[end + 1]);
        previous_end[end] = match start {
            Some(s) => s as i64 - 1,
            None => 0,
        };
        least_plus_penalty[end + 1] = min_cost + penalty;
    }
    let mut ends = Vec::new();
    let mut last = n as i64 - 1;
    while last >= 0 {
        ends.push(last as usize);
        let next = previous_end[last as usize];
        if next == last {
            bail!("segmentation does not terminate (Java loops for ever here): no finite cost for a first value");
        }
        last = next;
    }
    ends.reverse();
    let mut lengths = Vec::with_capacity(ends.len());
    let mut start = 0;
    for e in ends {
        lengths.push(e + 1 - start);
        start = e + 1;
    }
    Ok(lengths)
}

/// The least cost over the starts `0..=end` and the first start that has it (None when no cost is
/// below Double.MAX_VALUE). `prior_sums[s]` and `prior_squares[s]` are the cumulative sums before `s`.
/// Eight lanes keep their own first minimum; merging them by (cost, start) gives the sequential answer.
#[inline]
fn cheapest_start(prior_sums: &[f64], prior_squares: &[f64], least_plus_penalty: &[f64], sum_end: f64, square_end: f64) -> (f64, Option<usize>) {
    const L: usize = 8;
    let m = prior_sums.len(); // starts 0..m, segment length m - start
    let chunks = m / L;
    let mut best = [f64::MAX; L];
    let mut best_at = [usize::MAX; L];
    let mut lens = [0.0f64; L];
    for (l, len) in lens.iter_mut().enumerate() {
        *len = (m - l) as f64;
    }
    for c in 0..chunks {
        let base = c * L;
        let ps = &prior_sums[base..base + L];
        let pq = &prior_squares[base..base + L];
        let lp = &least_plus_penalty[base..base + L];
        for l in 0..L {
            let d = sum_end - ps[l];
            let sum_squared = d * d;
            let sum_of_squares = square_end - pq[l];
            let segment_cost = sum_of_squares - (sum_squared / lens[l]);
            let cost = lp[l] + segment_cost;
            let better = cost < best[l];
            best[l] = if better { cost } else { best[l] };
            best_at[l] = if better { base + l } else { best_at[l] };
        }
        for len in lens.iter_mut() {
            *len -= L as f64;
        }
    }
    let mut min_cost = f64::MAX;
    let mut min_at = usize::MAX;
    for l in 0..L {
        if best_at[l] != usize::MAX && (best[l] < min_cost || (best[l] == min_cost && best_at[l] < min_at)) {
            min_cost = best[l];
            min_at = best_at[l];
        }
    }
    for s in chunks * L..m {
        let d = sum_end - prior_sums[s];
        let sum_squared = d * d;
        let sum_of_squares = square_end - prior_squares[s];
        let segment_cost = sum_of_squares - (sum_squared / (m - s) as f64);
        let cost = least_plus_penalty[s] + segment_cost;
        if cost < min_cost {
            min_cost = cost;
            min_at = s;
        }
    }
    (min_cost, if min_at == usize::MAX { None } else { Some(min_at) })
}

/// Segmenter(y, penalty).pcf()
pub fn pcf(y: &[f64], penalty: f64) -> Result<PiecewiseConstantFit> {
    let lengths = least_cost_segment_lengths(y, penalty)?;
    Ok(Segmentation::from_lengths(y, &lengths).pcf())
}

/// Segmentation: a series cut into consecutive segments.
#[derive(Clone, Debug, PartialEq)]
pub struct Segmentation {
    pub segments: Vec<Vec<f64>>,
}

impl Segmentation {
    pub fn new(segments: Vec<Vec<f64>>) -> Self {
        Segmentation { segments }
    }

    pub fn from_lengths(y: &[f64], lengths: &[usize]) -> Self {
        let mut segments = Vec::with_capacity(lengths.len());
        let mut start = 0;
        for &len in lengths {
            segments.push(y[start..start + len].to_vec());
            start += len;
        }
        Segmentation { segments }
    }

    /// Segmenter.segmentBy: cut after each of the given (increasing) indices, the last of which must be
    /// the last index of `y`.
    pub fn segment_by(y: &[f64], segment_end_indices: &[usize]) -> Result<Self> {
        if segment_end_indices.windows(2).any(|w| w[1] <= w[0]) {
            bail!("Segment endpoint indices must be increasing");
        }
        if segment_end_indices.last().copied() != Some(y.len().wrapping_sub(1)) {
            bail!("Last segment endpoint index must be the last index of the input array");
        }
        let mut segments = Vec::new();
        let mut start = 0;
        for &e in segment_end_indices {
            segments.push(y[start..=e].to_vec());
            start = e + 1;
        }
        Ok(Segmentation { segments })
    }

    pub fn complete_series(&self) -> Vec<f64> {
        self.segments.iter().flatten().copied().collect()
    }

    pub fn pcf(&self) -> PiecewiseConstantFit {
        let mut fit = PiecewiseConstantFit { lengths: Vec::new(), start_positions: Vec::new(), means: Vec::new() };
        let mut cursor = 0;
        for s in &self.segments {
            fit.lengths.push(s.len());
            fit.start_positions.push(cursor);
            cursor += s.len();
            fit.means.push(java_round(mean(s) * 1000.0) as f64 / 1000.0);
        }
        fit
    }

    pub fn cost(&self, gamma: f64) -> f64 {
        let mut sum = 0.0;
        for s in &self.segments {
            sum += segment_cost(s);
        }
        sum + gamma * self.segments.len() as f64
    }
}

fn segment_cost(segment: &[f64]) -> f64 {
    if segment.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0;
    let mut sum_squared = 0.0;
    for v in segment {
        sum += v;
        sum_squared += v * v;
    }
    let mean = sum / segment.len() as f64;
    sum_squared - segment.len() as f64 * mean * mean
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    const D: f64 = 1e-10;

    fn assert_close(expected: &[f64], actual: &[f64], delta: f64) {
        assert_eq!(expected.len(), actual.len(), "{expected:?} {actual:?}");
        for (e, a) in expected.iter().zip(actual) {
            assert!((e - a).abs() <= delta, "{expected:?} {actual:?}");
        }
    }

    fn fit(lengths: &[usize], starts: &[usize], means: &[f64]) -> PiecewiseConstantFit {
        PiecewiseConstantFit { lengths: lengths.to_vec(), start_positions: starts.to_vec(), means: means.to_vec() }
    }

    fn seg(parts: &[&[f64]]) -> Segmentation {
        Segmentation::new(parts.iter().map(|p| p.to_vec()).collect())
    }

    fn round4(v: f64) -> f64 {
        java_round(v * 10000.0) as f64 / 10000.0
    }

    // ---- StatsTest
    #[test]
    fn median_absolute_deviation_test() {
        assert!((median_absolute_deviation(&[1.0]) - 0.0).abs() < 0.001);
        assert!((median_absolute_deviation(&[1.0, 3.0, 5.0]) - 2.9652).abs() < 0.001);
        assert!((median_absolute_deviation(&[1.0, 3.0, 5.0, 2.0, 8.0, 7.0]) - 3.7064).abs() < 0.001);
        // DoubleArrayTest (its other tests cover helpers of the Java test code only)
        assert!((median_absolute_deviation(&[1.0, 2.0, 1.0]) - 0.0).abs() < 1e-10);
        assert!((median_absolute_deviation(&[10.0, 12.0, 11.0, 14.0]) - 1.4826).abs() < 1e-10);
    }

    // ---- WindowedMedianTest
    #[test]
    fn windowed_median_compare_with_r() {
        let check = |data: &[f64], expected: &[f64], w: usize| assert_close(expected, &windowed_median(data, w).unwrap(), D);
        check(&[1., 2., 3., 4., 5.], &[1., 2., 3., 4., 5.], 1);
        check(&[1., 2., 3., 4., 5.], &[1., 2., 3., 4., 5.], 3);
        check(&[1., 2., 3., 4., 5.], &[1., 2., 3., 4., 5.], 5);
        check(&[1., 3., 5., 2., 8., 7.], &[1., 3., 5., 2., 8., 7.], 1);
        check(&[1., 3., 5., 2., 8., 7.], &[1., 3., 3., 5., 7., 7.], 3);
        check(&[1., 3., 5., 2., 8., 7.], &[1., 3., 3., 5., 8., 7.], 5);
        check(&[1., 0., 0., 4., 5., 6., 7., 8., 9.], &[1., 0., 0., 4., 5., 6., 7., 8., 9.], 1);
        check(&[1., 0., 0., 4., 5., 6., 7., 8., 9.], &[1., 0., 0., 4., 5., 6., 7., 8., 9.], 3);
        check(&[1., 0., 0., 4., 5., 6., 7., 8., 9.], &[1., 0., 1., 4., 5., 6., 7., 8., 9.], 5);
        check(&[1., 0., -1., 2., 5., 6., 1., 8., 9.], &[1., 0., 0., 2., 5., 5., 6., 8., 9.], 3);
        check(&[1., 0., -1., 2., 5., 6., 1., 8., 9.], &[1., 0., 1., 2., 2., 5., 6., 8., 9.], 5);
    }

    /// The heaps against a plain sorted window, on series with many ties
    #[test]
    fn windowed_median_matches_sorting() {
        let mut state = 12345u64;
        for case in 0..200 {
            let n = 1 + case % 97;
            let w = [1usize, 3, 5, 7, 11, 51][case % 6];
            if w > n {
                continue;
            }
            let data: Vec<f64> = (0..n)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    ((state >> 33) % 7) as f64 * 0.5 - 1.0
                })
                .collect();
            let medians = windowed_median(&data, w).unwrap();
            let half = w / 2;
            for i in half..n - half {
                let mut window = data[i - half..=i + half].to_vec();
                sort_doubles(&mut window);
                assert_eq!(medians[i], window[half]);
            }
        }
    }

    /// Java's heaps compare with `<=` (-0.0 equals 0.0) but remove by Double.equals (-0.0 is not 0.0): with
    /// both zeros in a window a value can stay behind, and the median goes wrong. Not reachable from AMBER or
    /// COBALT, whose segmented values are never -0.0; reproduced anyway (Java's Runmed gives the same, checked
    /// with the harness).
    #[test]
    fn signed_zeros_as_java() {
        // the window [-0.0, 1, 1] at index 3 has median 1; Java's heaps still hold the first 0.0 and give 0
        let data = [0.0, -0.0, -0.0, 1.0, 1.0];
        let medians = crate::compat::with_java(true, || windowed_median(&data, 3).unwrap());
        assert_eq!(format!("{medians:?}"), "[0.0, -0.0, -0.0, 0.0, 1.0]");
        // SG3 corrected: zeros folded in the heaps; the copied ends keep the data's own values
        let medians = crate::compat::with_java(false, || windowed_median(&data, 3).unwrap());
        assert_eq!(format!("{medians:?}"), "[0.0, 0.0, 0.0, 1.0, 1.0]");
        let medians = crate::compat::with_java(false, || windowed_median(&[-0.0, -0.0, 0.0, 1.0, 1.0], 3).unwrap());
        assert_eq!(format!("{medians:?}"), "[-0.0, 0.0, 0.0, 1.0, 1.0]");
    }

    // ---- RunmedTest
    #[test]
    fn runmed_compare_with_r() {
        let check = |data: &[f64], expected: &[f64], k: i32| assert_close(expected, &runmed(data, k, true).unwrap(), D);
        check(&[1., 2., 3., 4., 5.], &[1., 2., 3., 4., 5.], 1);
        check(&[1., 2., 3., 4., 5.], &[1., 2., 3., 4., 5.], 3);
        check(&[1., 2., 3., 4., 5.], &[1., 2., 3., 4., 5.], 5);
        check(&[1., 1., 2., 2., 3., 3.], &[1., 1., 2., 2., 3., 3.], 3);
        check(&[1., 3., 5., 2., 8., 7., 8., 3., 12.], &[3., 3., 3., 5., 7., 8., 7., 8., 10.], 3);
        check(&[1., 3., 5., 2., 8., 7., 8., 3., 12.], &[3., 3., 3., 5., 7., 7., 8., 8., 8.], 5);
        let x = [3., 1., 4., 1., 5., 9., 2., 6., 5., 3., 5., 9.];
        check(&x, &[3., 3., 1., 4., 5., 5., 6., 5., 5., 5., 5., 5.], 3);
        check(&x, &[3., 3., 3., 4., 4., 5., 5., 5., 5., 5., 5., 5.], 5);
        check(&x, &[3., 3., 3., 3., 4., 5., 5., 5., 5., 5., 5., 5.], 7);
        check(&[1., 1., -5., -6., 1., 2.], &[1., 1., -5., -6., 1., 2.], 1);
        check(&[1., 1., -5., -6., 1., 2.], &[1., 1., -5., -5., 1., 2.], 3);
        check(&[1., 1., -5., -6., 1., 2.], &[1., 1., 1., 1., 1., 1.], 5);
    }

    #[test]
    fn runmed_empty_single_and_invalid() {
        assert!(runmed(&[], 1, true).unwrap().is_empty());
        assert_eq!(runmed(&[42.0], 1, true).unwrap(), vec![42.0]);
        for k in [2, 0, -1] {
            assert!(runmed(&[1.0, 2.0, 3.0], k, true).is_err());
        }
    }

    // ---- GammaTest
    #[test]
    fn gamma_minimum_cost() {
        assert!((segment_penalty(&[1., 2., 3., 4.], 28.0, true).unwrap() - 0.01 * 28.0).abs() < 0.0001);
        assert!(segment_penalty(&[], 28.0, true).is_err());
        assert_eq!(segment_penalty(&[1.0, 5.0], 28.0, false).unwrap(), 28.0);
    }

    #[test]
    fn gamma_filter_width() {
        let cases = [(0, -1), (1, 1), (2, 1), (3, 3), (4, 3), (5, 5), (6, 5), (49, 49), (50, 49), (51, 51), (52, 51), (53, 51), (54, 51), (10000, 51)];
        for (n, w) in cases {
            assert_eq!(filter_width(n), w, "{n}");
        }
    }

    #[test]
    fn gamma_compare_with_r() {
        let Some(dir) = super::super::tests::resources_dir() else { return };
        let gamma_for = |f: &str| round4(segment_penalty(&read_column(&dir, f, "ratio", false), 50.0, true).unwrap());
        assert!((gamma_for("ratios1.tsv") - 0.2361).abs() < 0.0001);
        assert!((gamma_for("ratios4.tsv") - 0.3325).abs() < 0.0001);
        assert!((gamma_for("ratios5.csv") - 0.5).abs() < 0.0001);
        assert!((gamma_for("ratios_chr13_1.tsv") - 0.3099).abs() < 0.0001);
    }

    /// SegmentationTestBase.readResourcesFile and Double.parseDouble of one column (rounded to ten decimals
    /// as ExamplesTest.readDoubles does when `round10`)
    pub(crate) fn read_column(dir: &std::path::Path, file: &str, column: &str, round10: bool) -> Vec<f64> {
        let text = std::fs::read_to_string(dir.join("segmentation").join(file)).unwrap();
        let sep = if file.ends_with("csv") { ',' } else { '\t' };
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split(sep).collect();
        let col = header.iter().position(|h| *h == column).unwrap();
        lines
            .map(|l| {
                let v: f64 = l.split(sep).nth(col).unwrap().trim().parse().unwrap();
                if round10 { java_round(v * 10_000_000_000.0) as f64 / 10_000_000_000.0 } else { v }
            })
            .collect()
    }

    // ---- SegmentationTest
    #[test]
    fn segmentation_complete_series() {
        assert_close(&[1., 1.], &seg(&[&[1., 1.]]).complete_series(), D);
        assert_close(&[1., 1., 2., 2.1], &seg(&[&[1., 1.], &[2., 2.1]]).complete_series(), D);
        assert_close(&[10., 1., 2., 2.1], &seg(&[&[10.], &[1.], &[2., 2.1]]).complete_series(), D);
    }

    #[test]
    fn segmentation_cost() {
        let c = |parts: &[&[f64]], expected: f64| assert!((seg(parts).cost(50.0) - expected).abs() < D, "{parts:?}");
        c(&[&[2.], &[3.], &[1.], &[12.], &[13.], &[11.]], 300.0);
        c(&[&[2., 3., 1.]], 52.0);
        c(&[&[2., 3., 1.], &[12., 13., 11.]], 104.0);
        c(&[&[2., 3.], &[1.], &[12.], &[13.], &[11.]], 250.5);
        c(&[&[2., 3., 1.], &[12.], &[13.], &[11.]], 202.0);
        c(&[&[2., 3., 1.], &[12., 13.], &[11.]], 152.5);
        c(&[&[2., 3., 1., 12.], &[13., 11.]], 179.0);
        c(&[&[2., 3.], &[1., 12.], &[13., 11.]], 213.0);
        c(&[&[2., 3.], &[1., 12., 13., 11.]], 193.25);
        c(&[&[2., 3., 1., 12.], &[13.], &[11.]], 227.0);
        c(&[&[2., 3., 1., 12., 13., 11.]], 204.0);
        c(&[&[2.], &[3.], &[1.], &[12.], &[13., 11.]], 252.0);
    }

    // ---- SegmenterTest
    fn segmenter_cost(y: &[f64], penalty: f64) -> f64 {
        Segmentation::from_lengths(y, &least_cost_segment_lengths(y, penalty).unwrap()).cost(penalty)
    }

    #[test]
    fn compare_algorithm_with_exhaustive_search() {
        let cases: [(&[f64], f64); 10] = [
            (&[0.5, 0.51, 0.49, 10.09, 9.92, 10.02], 50.0),
            (&[0.5, 0.51, 0.49, 0.52], 50.0),
            (&[1., 1., 2., 2.], 50.0),
            (&[1., 1., 2., 2.], 10.0),
            (&[-1., 1., -2., 2.], 10.0),
            (&[2., 3., 1., 12., 13., 11.], 50.0),
            (&[5., 4., 1., 1.], 0.0),
            (&[5., 4., 1., 1.], 50.0),
            (&[1., 1., -5., -6., 1., 2.], 50.0),
            (&[55., 55., -276., -275., 55., 56.], 50.0),
        ];
        for (y, gamma) in cases {
            let penalty = segment_penalty(y, gamma, true).unwrap();
            let fast = segmenter_cost(y, penalty);
            let exhaustive = cheapest_by_exhaustive_search(y, penalty).cost(penalty);
            assert!((fast - exhaustive).abs() < 0.001, "{y:?}");
        }
    }

    #[test]
    fn segment_by() {
        let y = [2., 3., 1., 12., 13., 11.];
        assert_eq!(seg(&[&[2.], &[3.], &[1.], &[12.], &[13.], &[11.]]), Segmentation::segment_by(&y, &[0, 1, 2, 3, 4, 5]).unwrap());
        assert_eq!(seg(&[&[2.], &[3.], &[1.], &[12.], &[13., 11.]]), Segmentation::segment_by(&y, &[0, 1, 2, 3, 5]).unwrap());
        assert_eq!(seg(&[&[2.], &[3., 1., 12.], &[13., 11.]]), Segmentation::segment_by(&y, &[0, 3, 5]).unwrap());
        assert_eq!(seg(&[&[2.], &[3., 1., 12., 13., 11.]]), Segmentation::segment_by(&y, &[0, 5]).unwrap());
        assert_eq!(seg(&[&[2., 3., 1., 12., 13., 11.]]), Segmentation::segment_by(&y, &[5]).unwrap());
        assert!(Segmentation::segment_by(&y, &[3, 2, 5]).is_err());
        assert!(Segmentation::segment_by(&y, &[2, 4]).is_err());
    }

    #[test]
    fn set_minimum_cost() {
        let y = [0.3, 1.5, 1.8, 2.4];
        let penalty = segment_penalty(&y, 50.0, true).unwrap();
        assert_eq!(fit(&[1, 3], &[0, 1], &[0.3, 1.9]), pcf(&y, penalty).unwrap());
        assert!((penalty - 0.5).abs() < 0.0001);
    }

    #[test]
    fn pcf_test() {
        assert_eq!(fit(&[1], &[0], &[1.0]), pcf(&[1.0], 50.0).unwrap());
        assert_eq!(fit(&[2], &[0], &[1.0]), pcf(&[1.0, 1.0], 50.0).unwrap());
        assert_eq!(fit(&[3], &[0], &[1.0]), pcf(&[1.0, 1.0, 1.0], 50.0).unwrap());
        assert_eq!(fit(&[3, 3], &[0, 3], &[2.0, 12.0]), pcf(&[2., 3., 1., 12., 13., 11.], 50.0).unwrap());
        assert_eq!(fit(&[3, 3, 3], &[0, 3, 6], &[2.0, 12.0, 2.0]), pcf(&[2., 3., 1., 12., 13., 11., 2., 3., 1.], 50.0).unwrap());
        assert_eq!(fit(&[3, 3, 3], &[0, 3, 6], &[2.0, 12.0, 22.0]), pcf(&[2., 3., 1., 12., 13., 11., 22., 23., 21.], 50.0).unwrap());
    }

    /// The first value infinite: no finite cost for it, where Java's backtracking never ends
    #[test]
    fn non_terminating_case_is_an_error() {
        assert!(least_cost_segment_lengths(&[f64::INFINITY, 1.0], 1.0).is_err());
        // later infinite values give Java's odd but finite answer: the first value alone, then the rest
        assert_eq!(least_cost_segment_lengths(&[1.0, 2.0, f64::INFINITY, 3.0], 1.0).unwrap(), vec![1, 3]);
    }

    /// The eight-lane minimum against the plain sequential loop, with many equal costs
    #[test]
    fn lanes_match_sequential_loop() {
        fn sequential(y: &[f64], penalty: f64) -> Vec<usize> {
            let n = y.len();
            let mut cs = vec![0.0; n];
            let mut css = vec![0.0; n];
            cs[0] = y[0];
            css[0] = y[0] * y[0];
            for i in 1..n {
                cs[i] = cs[i - 1] + y[i];
                css[i] = css[i - 1] + (y[i] * y[i]);
            }
            let mut least = vec![0.0; n + 1];
            let mut prev = vec![0i64; n];
            for end in 0..n {
                let mut min_cost = f64::MAX;
                let mut best = 0i64;
                for start in 0..=end {
                    let sum_squared = if start == 0 { cs[end] * cs[end] } else { (cs[end] - cs[start - 1]) * (cs[end] - cs[start - 1]) };
                    let sum_of_squares = if start == 0 { css[end] } else { css[end] - css[start - 1] };
                    let len = (end - start + 1) as f64;
                    let cost = least[start] + penalty + (sum_of_squares - (sum_squared / len));
                    if cost < min_cost {
                        min_cost = cost;
                        best = start as i64 - 1;
                    }
                }
                least[end + 1] = min_cost;
                prev[end] = best;
            }
            let mut ends = Vec::new();
            let mut last = n as i64 - 1;
            while last >= 0 {
                ends.push(last as usize);
                last = prev[last as usize];
            }
            ends.reverse();
            let mut start = 0;
            ends.iter().map(|e| { let l = e + 1 - start; start = e + 1; l }).collect()
        }
        let mut state = 99u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            state >> 11
        };
        for case in 0..300 {
            let n = 1 + (next() % 70) as usize;
            let levels = 1 + next() % 4;
            let y: Vec<f64> = (0..n).map(|_| ((next() % levels) as f64) * if case % 3 == 0 { 1.0 } else { 0.25 } + if case % 2 == 0 { 0.0 } else { (next() % 1000) as f64 * 1e-3 }).collect();
            let penalty = [0.0, 0.01, 0.5, 2.0][case % 4];
            assert_eq!(least_cost_segment_lengths(&y, penalty).unwrap(), sequential(&y, penalty), "{y:?} {penalty}");
        }
    }

    // ---- ExhaustiveSearchSegmenter and its tests
    fn all_possible_segmentations(y: &[f64]) -> Vec<Segmentation> {
        let n = y.len();
        let mut result = Vec::new();
        for mask in 0..(1u32 << (n - 1)) {
            let mut ends: Vec<usize> = (0..n - 1).filter(|i| mask & (1 << i) != 0).collect();
            ends.push(n - 1);
            result.push(Segmentation::segment_by(y, &ends).unwrap());
        }
        result
    }

    fn cheapest_by_exhaustive_search(y: &[f64], penalty: f64) -> Segmentation {
        let mut cheapest = None;
        let mut min_cost = f64::MAX;
        for s in all_possible_segmentations(y) {
            let c = s.cost(penalty);
            if c < min_cost {
                min_cost = c;
                cheapest = Some(s);
            }
        }
        cheapest.unwrap()
    }

    #[test]
    fn exhaustive_search_segmenter_test() {
        let set = |y: &[f64]| -> HashSet<String> { all_possible_segmentations(y).iter().map(|s| format!("{:?}", s.segments)).collect() };
        let has = |s: &HashSet<String>, parts: &[&[f64]]| s.contains(&format!("{:?}", seg(parts).segments));
        let s1 = set(&[1.]);
        assert_eq!(1, s1.len());
        assert!(has(&s1, &[&[1.]]));
        let s2 = set(&[1., 1.]);
        assert_eq!(2, s2.len());
        assert!(has(&s2, &[&[1.], &[1.]]) && has(&s2, &[&[1., 1.]]));
        let s3 = set(&[1., 2., 3.]);
        assert_eq!(4, s3.len());
        assert!(has(&s3, &[&[1.], &[2.], &[3.]]) && has(&s3, &[&[1., 2.], &[3.]]) && has(&s3, &[&[1.], &[2., 3.]]) && has(&s3, &[&[1., 2., 3.]]));
        let s4 = set(&[1., 2., 3., 1.]);
        assert_eq!(8, s4.len());
        for parts in [
            vec![vec![1.], vec![2.], vec![3.], vec![1.]],
            vec![vec![1., 2.], vec![3.], vec![1.]],
            vec![vec![1.], vec![2., 3.], vec![1.]],
            vec![vec![1.], vec![2.], vec![3., 1.]],
            vec![vec![1., 2.], vec![3., 1.]],
            vec![vec![1.], vec![2., 3., 1.]],
            vec![vec![1., 2., 3.], vec![1.]],
            vec![vec![1., 2., 3., 1.]],
        ] {
            assert!(s4.contains(&format!("{parts:?}")));
        }
    }

    #[test]
    fn exhaustive_search_cheapest_test() {
        let c = |y: &[f64]| cheapest_by_exhaustive_search(y, 50.0).pcf();
        assert_eq!(fit(&[1], &[0], &[1.0]), c(&[1.0]));
        assert_eq!(fit(&[2], &[0], &[1.0]), c(&[1.0, 1.0]));
        assert_eq!(fit(&[3], &[0], &[1.0]), c(&[1.0, 1.0, 1.0]));
        assert_eq!(fit(&[3, 3], &[0, 3], &[2.0, 12.0]), c(&[2., 3., 1., 12., 13., 11.]));
        assert_eq!(fit(&[3, 3, 3], &[0, 3, 6], &[2.0, 12.0, 2.0]), c(&[2., 3., 1., 12., 13., 11., 2., 3., 1.]));
        assert_eq!(fit(&[3, 3, 3], &[0, 3, 6], &[2.0, 12.0, 22.0]), c(&[2., 3., 1., 12., 13., 11., 22., 23., 21.]));
    }

    // ---- ExamplesTest (ignored in Java too: the large arms take a while)
    fn pcf_for_file(file: &str, gamma: f64, column: &str) -> Option<PiecewiseConstantFit> {
        let dir = super::super::tests::resources_dir()?;
        let y = read_column(&dir, file, column, true);
        let penalty = segment_penalty(&y, gamma, true).unwrap();
        Some(pcf(&y, penalty).unwrap())
    }

    #[test]
    fn examples_amber() {
        let Some(p) = pcf_for_file("amber1.tsv", 100.0, "tumorModifiedBAF") else { return };
        assert_eq!(fit(&[5, 5, 5], &[0, 5, 10], &[0.515, 0.962, 0.511]), p);
    }

    #[test]
    fn examples_amber2() {
        let Some(dir) = super::super::tests::resources_dir() else { return };
        let y = read_column(&dir, "amber2.tsv", "tumorModifiedBAF", true);
        let penalty = 0.036533517458912736;
        let ours = segmenter_cost(&y, penalty);
        let r = Segmentation::from_lengths(&y, &[1, 1, 108, 1, 19]).cost(penalty);
        assert!(ours < r);
    }

    #[test]
    fn examples_minima() {
        let Some(p) = pcf_for_file("ratios1.tsv", 50.0, "ratio") else { return };
        assert_eq!(fit(&[16], &[0], &[0.883]), p);
        assert_eq!(fit(&[36], &[0], &[0.902]), pcf_for_file("ratios2.tsv", 50.0, "ratio").unwrap());
        assert_eq!(fit(&[75], &[0], &[0.899]), pcf_for_file("ratios3.tsv", 50.0, "ratio").unwrap());
        assert_eq!(fit(&[109], &[0], &[0.889]), pcf_for_file("ratios4.tsv", 50.0, "ratio").unwrap());
    }

    #[test]
    fn examples_panel_9p() {
        let Some(p) = pcf_for_file("ratios_tii_9p.tsv", 35.0, "logRatio") else { return };
        assert_eq!(fit(&[86, 4, 19], &[0, 86, 90], &[-0.153, -0.764, -0.193]), p);
        assert_eq!(
            fit(&[7, 2, 77, 4, 3, 16], &[0, 7, 9, 86, 90, 93], &[-0.159, -0.594, -0.141, -0.764, 0.123, -0.252]),
            pcf_for_file("ratios_tii_9p.tsv", 10.0, "logRatio").unwrap()
        );
        assert_eq!(fit(&[109], &[0], &[-0.227]), pcf_for_file("ratios_27t_9p.tsv", 50.0, "logRatio").unwrap());
        assert_eq!(fit(&[86, 7, 16], &[0, 86, 93], &[-0.195, -0.604, -0.236]), pcf_for_file("ratios_27t_9p.tsv", 18.0, "logRatio").unwrap());
    }

    #[test]
    fn examples_chr13() {
        let Some(p) = pcf_for_file("ratios_chr13_1.tsv", 50.0, "ratio") else { return };
        assert_eq!(fit(&[108, 25, 48], &[0, 108, 133], &[0.546, 0.185, 0.519]), p);
    }

    #[test]
    #[ignore = "about 10^10 candidate segments; run with --ignored"]
    fn examples_large_arms() {
        let Some(p) = pcf_for_file("chr1p.ratios.tsv", 50.0, "logRatio") else { return };
        assert_eq!(146, p.lengths.len());
        assert_eq!(&p.lengths[..4], &[7, 4, 3, 54]);
        assert_eq!((p.lengths[141], p.lengths[143], p.lengths[145]), (2767, 2, 1));

        let p = pcf_for_file("chr2q.ratios.tsv", 100.0, "logRatio").unwrap();
        assert_eq!(59, p.lengths.len());
        assert_eq!(&p.lengths[..5], &[5976, 1, 1032, 1, 1559]);
        assert_eq!(&p.start_positions[..5], &[0, 5976, 5977, 7009, 7010]);

        assert_eq!(
            fit(&[35543, 89, 51], &[0, 35543, 35632], &[-0.007, 1.394, 2.459]),
            pcf_for_file("chr10p.ratios.tsv", 2000.0, "logRatio").unwrap()
        );
        assert_eq!(
            fit(
                &[24, 104, 5, 39, 2849, 174, 72458, 15, 6684],
                &[0, 24, 128, 133, 172, 3021, 3195, 75653, 75668],
                &[5.135, 0.264, 6.383, 1.076, -0.021, 0.961, -0.016, 2.88, -0.016]
            ),
            pcf_for_file("chr10q.ratios.tsv", 5000.0, "logRatio").unwrap()
        );
    }
}
