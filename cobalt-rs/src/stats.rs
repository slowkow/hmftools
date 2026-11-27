//! The commons-math3 3.6.1 statistics COBALT uses (DescriptiveStatistics.getMean and getPercentile(50)), and
//! COBALT's RollingMedian.

use std::collections::BTreeMap;

/// Mean.evaluate: the sum divided by the count, plus the mean of the residuals (a second pass); NaN for
/// no values. The values are summed in their order.
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let n = values.len() as f64;
    let mut sum = 0.0;
    for &v in values {
        sum += v;
    }
    let xbar = sum / n;
    let mut correction = 0.0;
    for &v in values {
        correction += v - xbar;
    }
    xbar + correction / n
}

/// Percentile(50) with the LEGACY estimation and NaN values removed, as DescriptiveStatistics.getPercentile(50):
/// position p (n + 1); below 1 the minimum, from n on the maximum, otherwise
/// `lower + (pos - floor(pos)) * (upper - lower)` of the order statistics around it. NaN for no values.
/// The slice is reordered.
pub fn median(values: &mut [f64]) -> f64 {
    percentile(values, 50.0)
}

pub fn percentile(values: &mut [f64], p: f64) -> f64 {
    let length = values.len();
    if length == 0 {
        return f64::NAN;
    }
    if length == 1 {
        return values[0];
    }
    // NaNStrategy.REMOVED
    let mut n = length;
    if values.iter().any(|v| v.is_nan()) {
        let mut k = 0;
        for i in 0..length {
            if !values[i].is_nan() {
                values[k] = values[i];
                k += 1;
            }
        }
        n = k;
        if n == 0 {
            return f64::NAN;
        }
    }
    let work = &mut values[..n];
    let q = p / 100.0;
    let pos = if q == 0.0 {
        0.0
    } else if q == 1.0 {
        n as f64
    } else {
        q * (n as f64 + 1.0)
    };
    let fpos = pos.floor();
    let int_pos = fpos as usize;
    let dif = pos - fpos;
    if pos < 1.0 {
        return select(work, 0);
    }
    if pos >= n as f64 {
        return select(work, n - 1);
    }
    let lower = select(work, int_pos - 1);
    // the smallest value above the lower order statistic's place
    let upper = work[int_pos..].iter().copied().fold(f64::INFINITY, f64::min);
    lower + dif * (upper - lower)
}

/// The k-th smallest value (KthSelector gives the order statistic, whatever its pivoting). The values hold
/// no NaN here; -0.0 and 0.0 are taken as equal, as Java's `<` does.
fn select(work: &mut [f64], k: usize) -> f64 {
    let (_, v, _) = work.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
    *v
}

/// RollingMedian: a max-heap below the median and a min-heap above it, with Java's PriorityQueue semantics
/// kept as sorted multisets (only the multisets and their extremes reach the median).
///
/// CO2: Java removes a value from the max-heap when it is at most the median, else from the min-heap, and
/// PriorityQueue.remove does nothing when the value is not in that heap. When the min-heap holds the extra
/// value and the value to remove is the median (the min-heap's top), the max-heap is searched and the value
/// stays in the window for good. Corrected: the value is taken from the heap that holds it.
#[derive(Default)]
pub struct RollingMedian {
    /// values at most the median: count by bit pattern (positive doubles order as their bits)
    max_heap: Heap,
    min_heap: Heap,
}

#[derive(Default)]
struct Heap {
    counts: BTreeMap<u64, u32>,
    len: usize,
}

impl Heap {
    fn add(&mut self, v: f64) {
        *self.counts.entry(key(v)).or_insert(0) += 1;
        self.len += 1;
    }
    /// PriorityQueue.remove(Object): one value equal (Double.equals) to `v`, if there is one
    fn remove(&mut self, v: f64) -> bool {
        let k = key(v);
        match self.counts.get_mut(&k) {
            Some(c) => {
                *c -= 1;
                if *c == 0 {
                    self.counts.remove(&k);
                }
                self.len -= 1;
                true
            }
            None => false,
        }
    }
    fn min(&self) -> f64 {
        unkey(*self.counts.keys().next().unwrap())
    }
    fn max(&self) -> f64 {
        unkey(*self.counts.keys().next_back().unwrap())
    }
    fn pop_min(&mut self) -> f64 {
        let v = self.min();
        self.remove(v);
        v
    }
    fn pop_max(&mut self) -> f64 {
        let v = self.max();
        self.remove(v);
        v
    }
}

/// The values are ratios above 1e-10: their bit patterns order as the values do.
fn key(v: f64) -> u64 {
    debug_assert!(v > 0.0 && v.is_finite());
    v.to_bits()
}
fn unkey(k: u64) -> f64 {
    f64::from_bits(k)
}

impl RollingMedian {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, n: f64) {
        if self.is_empty() {
            self.min_heap.add(n);
        } else if hmf_common::java_math::double_compare(n, self.median()).is_le() {
            self.max_heap.add(n);
        } else {
            self.min_heap.add(n);
        }
        self.fix_chaos();
    }

    pub fn remove(&mut self, n: f64) {
        if !self.is_empty() {
            let below = hmf_common::java_math::double_compare(n, self.median()).is_le();
            let removed = if below { self.max_heap.remove(n) } else { self.min_heap.remove(n) };
            if !removed && !crate::compat::java("CO2") {
                if below {
                    self.min_heap.remove(n);
                } else {
                    self.max_heap.remove(n);
                }
            }
        }
        self.fix_chaos();
    }

    pub fn size(&self) -> usize {
        self.max_heap.len + self.min_heap.len
    }

    fn is_empty(&self) -> bool {
        self.size() == 0
    }

    fn fix_chaos(&mut self) {
        let (a, b) = (self.max_heap.len, self.min_heap.len);
        if a.abs_diff(b) > 1 {
            if a > b {
                let v = self.max_heap.pop_max();
                self.min_heap.add(v);
            } else {
                let v = self.min_heap.pop_min();
                self.max_heap.add(v);
            }
        }
    }

    pub fn median(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let (a, b) = (self.max_heap.len, self.min_heap.len);
        if a == b {
            (self.max_heap.max() + self.min_heap.min()) / 2.0
        } else if a > b {
            self.max_heap.max()
        } else {
            self.min_heap.min()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_percentile() {
        assert!(median(&mut []).is_nan());
        assert_eq!(median(&mut [3.0]), 3.0);
        assert_eq!(median(&mut [3.0, 1.0]), 2.0);
        assert_eq!(median(&mut [3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&mut [4.0, 1.0, 3.0, 2.0]), 2.5);
        // lower + 0.5 * (upper - lower), not (lower + upper) / 2
        let (a, b) = (0.1f64, 0.7f64);
        assert_eq!(median(&mut [b, a]), a + 0.5 * (b - a));
        assert_eq!(median(&mut [f64::NAN, 2.0, 4.0]), 3.0);
        assert_eq!(percentile(&mut [1.0, 2.0, 3.0, 4.0], 25.0), 1.25);
    }

    #[test]
    fn two_pass_mean() {
        assert!(mean(&[]).is_nan());
        assert_eq!(mean(&[1.0, 2.0, 3.0, 4.0]), 2.5);
    }

    /// RollingMedianTest
    #[test]
    fn rolling_median_java_tests() {
        let mut m = RollingMedian::new();
        assert_eq!(m.median(), 0.0);
        m.add(1.0);
        assert_eq!(m.median(), 1.0);
        m.add(2.0);
        assert_eq!(m.median(), 1.5);
        m.add(3.0);
        assert_eq!(m.median(), 2.0);
        m.add(4.0);
        assert_eq!(m.median(), 2.5);
        m.remove(1.0);
        assert_eq!(m.median(), 3.0);
        m.remove(4.0);
        assert_eq!(m.median(), 2.5);
        m.remove(2.0);
        assert_eq!(m.median(), 3.0);
        m.remove(3.0);
        assert_eq!(m.median(), 0.0);
    }

    #[test]
    fn rolling_median_lost_removal() {
        // {1, 2, 3}: max-heap {1}, min-heap {2, 3}, median 2; removing 2 looks in the max-heap
        let run = |java: bool| {
            crate::compat::with_java(java, || {
                let mut m = RollingMedian::new();
                for v in [1.0, 2.0, 3.0] {
                    m.add(v);
                }
                m.remove(2.0);
                (m.size(), m.median())
            })
        };
        assert_eq!(run(true), (3, 2.0)); // Java: 2 stays
        assert_eq!(run(false), (2, 2.0)); // corrected: {1, 3}
    }
}
