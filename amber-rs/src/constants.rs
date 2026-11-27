//! AmberConstants (the values the core uses; the purity and LOH constants live with those parts).

pub const DEFAULT_TUMOR_ONLY_MIN_DEPTH: i32 = 25;
pub const DEFAULT_TUMOR_MIN_DEPTH: i32 = 8;

pub const CONTAMINATON_MIN_NORMAL_READ_DEPTH: i32 = 7;
pub const CONTAMINATON_THREE_PLUS_READS_MIN: i64 = 10;
pub const CONTAMINATON_THREE_PLUS_READS_SITE_PERC: f64 = 0.03;
pub const CONTAMINATON_THREE_PLUS_READS_SITE_LOW_VAF_PERC: f64 = 0.002;
pub const CONTAMINATON_THREE_PLUS_READS_VAF_MIN: f64 = 0.05;

pub const QUAL_FILTERED_THRESHOLD: f64 = 0.15;

pub const HOMOZYGOUS_REGION_MIN_SIZE: i32 = 500_000;
pub const HOMOZYGOUS_REGION_MIN_SNP_LOCI_COUNT: i32 = 50;
pub const HOMOZYGOUS_REGION_WINDOW_SIZE: i32 = 50;
pub const HOMOZYGOUS_REGION_MAX_HET_IN_WINDOW: i32 = 5;
pub const HOMOZYGOUS_REGION_MAX_HET_RATIO: f64 = 0.05;
pub const HOMOZYGOUS_REGION_LONG_SIZE: i32 = 3_000_000;
pub const UNIPARENTAL_DISOMY_MIN_LENGTH: i64 = 10_000_000;

pub const CRAM_MIN_GAP_START: i32 = 10000;
pub const BAM_MIN_GAP_START: i32 = 4000;

pub const TARGET_REGION_SITE_BUFFER: i32 = 300;

/// hmf-common Doubles: comparisons with a tolerance of 1e-10
pub mod doubles {
    const EPSILON: f64 = 1e-10;
    pub fn less_than(value: f64, reference: f64) -> bool {
        value - reference < -EPSILON
    }
    pub fn less_or_equal(value: f64, reference: f64) -> bool {
        value - reference < EPSILON
    }
    pub fn greater_than(value: f64, reference: f64) -> bool {
        value - reference > EPSILON
    }
    pub fn greater_or_equal(value: f64, reference: f64) -> bool {
        value - reference > -EPSILON
    }
}

/// Integers.medianPositiveValue: the median of the values above zero, the mean of the middle two
/// (rounded down) for an even count, 0 for none
pub fn median_positive_value(values: impl IntoIterator<Item = i32>) -> i32 {
    let mut v: Vec<i32> = values.into_iter().filter(|&x| x > 0).collect();
    if v.is_empty() {
        return 0;
    }
    v.sort_unstable();
    let n = v.len();
    if n % 2 == 0 { (v[n / 2] + v[n / 2 - 1]) / 2 } else { v[n / 2] }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median() {
        assert_eq!(median_positive_value([0, 0, 0]), 0);
        assert_eq!(median_positive_value([3, 0, 1, 2]), 2);
        assert_eq!(median_positive_value([4, 1, 0, 2, 3]), 2); // (2+3)/2
        assert_eq!(median_positive_value([5]), 5);
    }
}
