//! CobaltConstants and the constants COBALT takes from hmf-common.

pub const DEFAULT_GC_RATIO_MIN: f64 = 0.24;
pub const DEFAULT_GC_RATIO_MAX: f64 = 0.68;

/// The window of a read depth (CobaltConstants.WINDOW_SIZE, GCProfileFactory.WINDOW_SIZE)
pub const WINDOW_SIZE: i32 = 1000;
/// Java reads the BAM in partitions of 100 Mb (CobaltConstants.PARTITION_SIZE). The port reads 10 Mb
/// partitions, for an even load over the threads: the counts do not depend on the partition length, as each
/// partition counts exactly the bases inside it and every partition starts at a window start.
pub const PARTITION_SIZE: i32 = 10_000_000;

pub const ROLLING_MEDIAN_MAX_DISTANCE: i64 = 5000;
pub const ROLLING_MEDIAN_MIN_COVERAGE: usize = 1000;

/// Sparse windows are not consolidated over a longer stretch (the centromeres)
pub const MAX_SPARSE_CONSOLIDATE_DISTANCE: i32 = 3_000_000;

/// GCProfile.MIN_MAPPABLE_PERCENTAGE
pub const MIN_MAPPABLE_PERCENTAGE: f64 = 0.85;

/// CobaltChromosomes cut-offs
pub const MIN_Y_COUNT: i32 = 1000;
pub const TWO_X_CUTOFF: f64 = 0.65;
pub const TWO_Y_CUTOFF: f64 = 0.75;
pub const Y_CUTOFF: f64 = 0.05;
pub const MOSAIC_X_CUTOFF: f64 = 0.8;
pub const TRISOMY_CUTOFF: f64 = 1.35;

/// GCBucket.calcGcBucket: (int) Math.round(gcContent * 100)
pub fn calc_gc_bucket(gc_content: f64) -> i32 {
    hmf_common::java_math::java_round(gc_content * 100.0) as i32
}

/// Doubles (hmf-common): comparisons with a tolerance of 1e-10
pub mod doubles {
    const EPSILON: f64 = 1e-10;
    pub fn equal(a: f64, b: f64) -> bool {
        (a - b).abs() < EPSILON
    }
    pub fn is_zero(v: f64) -> bool {
        equal(v, 0.0)
    }
    pub fn greater_than(v: f64, r: f64) -> bool {
        v - r > EPSILON
    }
    pub fn greater_or_equal(v: f64, r: f64) -> bool {
        v - r > -EPSILON
    }
    pub fn less_than(v: f64, r: f64) -> bool {
        v - r < -EPSILON
    }
    pub fn positive(v: f64) -> bool {
        greater_than(v, 0.0)
    }
    /// Doubles.median(List<Double>): sorted, the mean of the middle two for an even count, 0 for none
    pub fn median(values: &mut [f64]) -> f64 {
        if values.is_empty() {
            return 0.0;
        }
        values.sort_by(|a, b| hmf_common::java_math::double_compare(*a, *b));
        let n = values.len();
        if n.is_multiple_of(2) { (values[n / 2] + values[n / 2 - 1]) / 2.0 } else { values[n / 2] }
    }
}
