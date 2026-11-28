//! BiallelicModel (upstream bc9b187516, "Purple biallelic model changes (AUS-504)"): the probability that
//! no wildtype allele remains at a somatic variant, from the sample purity, the variant's read counts and
//! variant copy number, and the copy number and minor-allele copy number of its segment.
//!
//! The probability is P(LOH) * P(no wildtype | LOH) + (1 - P(LOH)) * P(no wildtype | no LOH):
//! - P(LOH): a logistic step in the minor-allele copy number at 0.5, steeper with purity
//!   (growth rate purity * 70, at least 20);
//! - P(no wildtype | LOH): the chance that a Poisson count, with mean the variant's alt read count moved by
//!   at most 2.5% of it towards the read count it would have at the no-wildtype copy number threshold
//!   (and dampened above 300 reads), exceeds that threshold count;
//! - P(no wildtype | no LOH): from the odds of the above and the floored P(LOH), then moved on a
//!   log-odds-like scale by ln(distance of the variant copy number to the nearest allele copy number /
//!   its distance to the segment copy number).

use hmf_common::commons_math::Poisson;
use hmf_common::java_math::double_compare;

use crate::constants::*;
use crate::types::*;

const PROB_UPPER_LIMIT: f64 = 1e10;
const PROB_LOWER_LIMIT: f64 = 1.0 / PROB_UPPER_LIMIT;

/// BiallelicModel.classifyIsBiallelic
pub fn classify_is_biallelic(probability: f64) -> bool {
    probability >= 0.5
}

/// BiallelicModel.calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN).
/// The caller returns 0 when the alt read count is 0, as Java's variant overload does.
pub fn calculate_probability(sample_purity: f64, allele_read_count: i32, depth: i32, segment_cn: f64, variant_cn: f64, minor_allele_cn: f64) -> f64 {
    if depth == 0 {
        return 0.0;
    }
    let ad = allele_read_count as f64;
    let vaf = ad / depth as f64;

    let loh_exponent = jmax(sample_purity * BIALLELIC_LOH_GROWTH_RATE as f64, BIALLELIC_LOH_MIN_PURITY_EXP);
    let probability_loh = 1.0 - 1.0 / (1.0 + jexp(-loh_exponent * (minor_allele_cn - 0.5)));

    let vcn_threshold = jmin(segment_cn - 0.5, jmax(BIALLELIC_VCN_THRESHOLD_PI, segment_cn - BIALLELIC_VCN_THRESHOLD_P2));
    let mut read_count_at_threshold = vcn_threshold / variant_cn * ad;

    let af_uncertainty = vaf * BIALLELIC_AF_UNCERTAINTY_FACTOR;
    let mut conservative_ad = jround(median3(
        read_count_at_threshold,
        ad - depth as f64 * af_uncertainty,
        ad + depth as f64 * af_uncertainty,
    )) as f64;

    if conservative_ad > BIALLELIC_DAMPENING_AD_LIMIT as f64 {
        let factor = jpow(conservative_ad - BIALLELIC_DAMPENING_AD_LIMIT as f64, BIALLELIC_DAMPENING_AD_EXP);
        let dampening_factor = (BIALLELIC_DAMPENING_AD_LIMIT as f64 + factor) / conservative_ad;
        conservative_ad = jround(conservative_ad * dampening_factor) as f64;
        read_count_at_threshold *= dampening_factor;
    }

    let threshold_count = threshold_read_count(jround(read_count_at_threshold));
    let prob_no_wildtype_given_loh = 1.0 - Poisson::new(conservative_ad).cumulative_probability(threshold_count);

    let mut prob_wildtype_given_loh_relative = if prob_no_wildtype_given_loh == 0.0 {
        PROB_UPPER_LIMIT
    } else {
        0.5 * (1.0 - prob_no_wildtype_given_loh) / prob_no_wildtype_given_loh
    };
    prob_wildtype_given_loh_relative = cap_probability(prob_wildtype_given_loh_relative);

    let probability_loh_floored = jmax(probability_loh, BIALLELIC_LOH_BASE_ERROR_RATE);
    let mut prob_no_wildtype_given_no_loh_raw = probability_loh_floored / (prob_wildtype_given_loh_relative + probability_loh_floored);
    prob_no_wildtype_given_no_loh_raw = cap_probability(prob_no_wildtype_given_no_loh_raw);

    let min_allele_diff = jmax(
        jmin((variant_cn - minor_allele_cn).abs(), (variant_cn - (segment_cn - minor_allele_cn)).abs()),
        PROB_LOWER_LIMIT,
    );
    let min_cn_diff = jmax((variant_cn - segment_cn).abs(), 1.0 / PROB_UPPER_LIMIT);

    let adjustment = jlog(min_allele_diff / min_cn_diff);

    let mut prob_no_wildtype_given_no_loh = calc_scaled_adjustment(prob_no_wildtype_given_no_loh_raw, adjustment);
    prob_no_wildtype_given_no_loh = jmin(jmax(prob_no_wildtype_given_no_loh, 0.0), 1.0);

    probability_loh * prob_no_wildtype_given_loh + (1.0 - probability_loh) * prob_no_wildtype_given_no_loh
}

/// `(int)Math.round(readCountAtThreshold)`. Java narrows the long to an int by keeping its low 32 bits, so a
/// count beyond the int range wraps (-infinity, from a variant copy number of 0, becomes 0). PU10: corrected,
/// the count saturates at the int range, as Java's double-to-int cast of the old model did.
fn threshold_read_count(rounded: i64) -> i32 {
    if crate::compat::java("PU10") {
        rounded as i32
    } else {
        rounded.clamp(i32::MIN as i64, i32::MAX as i64) as i32
    }
}

/// Doubles.median(List.of(a, b, c)): the middle value in the order of Double.compareTo
fn median3(a: f64, b: f64, c: f64) -> f64 {
    let mut v = [a, b, c];
    v.sort_by(|x, y| double_compare(*x, *y));
    v[1]
}

fn cap_probability(probability: f64) -> f64 {
    jmin(jmax(probability, PROB_LOWER_LIMIT), PROB_UPPER_LIMIT)
}

fn scaled(orig_value: f64, adjustment: f64, upper: bool) -> f64 {
    if upper {
        1.0 - (1.0 - orig_value) / jpow(std::f64::consts::E, adjustment)
    } else {
        orig_value * jpow(std::f64::consts::E, adjustment)
    }
}

/// BiallelicModel.calcScaledAdjustment: moves a probability by `adjustment` on the scale ln(2p) below 0.5 and
/// -ln(2(1 - p)) above it, so values near 0 and 1 move least.
fn calc_scaled_adjustment(orig_value: f64, orig_adjustment: f64) -> f64 {
    if orig_value > 0.5 {
        if orig_adjustment > 0.0 {
            scaled(orig_value, orig_adjustment, true)
        } else {
            let area_to_midpoint = -jlog(0.5) + jlog(1.0 - orig_value);
            if orig_adjustment - area_to_midpoint > 0.0 {
                scaled(orig_value, orig_adjustment, true)
            } else {
                scaled(0.5, orig_adjustment - area_to_midpoint, false)
            }
        }
    } else if orig_adjustment < 0.0 {
        scaled(orig_value, orig_adjustment, false)
    } else {
        let area_to_midpoint = jlog(0.5) - jlog(orig_value);
        if orig_adjustment - area_to_midpoint < 0.0 {
            scaled(orig_value, orig_adjustment, false)
        } else {
            scaled(0.5, orig_adjustment - area_to_midpoint, true)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// BiallelicProbabilityTest.testBiallelicCalculations (upstream's cases and tolerance)
    #[test]
    fn biallelic_calculations() {
        let cases: [(f64, i32, i32, f64, f64, f64, f64); 5] = [
            (0.34, 869, 1932, 4.04, 3.56, 2.0, 0.7136),
            // extreme subclonality: the relative calculation takes its fallback, giving 0
            (0.32, 8, 755, 1.58, 0.062, 0.552, 0.0),
            // depth dampening and AF uncertainty lower P(no wildtype | LOH); the adjustment lowers P(no wildtype | no LOH)
            (0.73, 1337, 2144, 2.74, 2.17, 0.889, 0.2989),
            // the variant copy number is closer to the minor allele's than to the segment's: moved down
            (0.26, 11, 107, 1.61, 0.751, 0.489, 0.0399),
            // the variant copy number is very close to the segment's and far from both alleles': moved far up
            (0.45, 20, 58, 1.39, 1.32, 0.63, 0.9138),
        ];
        for (purity, ad, depth, cn, vcn, macn, expected) in cases {
            let p = calculate_probability(purity, ad, depth, cn, vcn, macn);
            assert!((p - expected).abs() < 0.001, "purity {purity} ad {ad} depth {depth}: {p} != {expected}");
        }
    }

    #[test]
    fn zero_depth_gives_zero() {
        assert_eq!(0.0, calculate_probability(0.5, 0, 0, 2.0, 1.0, 1.0));
    }

    #[test]
    fn scaled_adjustment_is_continuous_and_monotone() {
        // the two halves meet at 0.5 whichever side the value starts on
        for p in [0.01, 0.2, 0.49, 0.5, 0.51, 0.8, 0.99] {
            let to_mid = if p > 0.5 { -jlog(0.5) + jlog(1.0 - p) } else { jlog(0.5) - jlog(p) };
            assert!((calc_scaled_adjustment(p, to_mid) - 0.5).abs() < 1e-12, "{p}");
            assert_eq!(p, calc_scaled_adjustment(p, 0.0));
            let mut last = 0.0;
            for i in -40..=40 {
                let v = calc_scaled_adjustment(p, i as f64 * 0.25);
                assert!(v >= last && (0.0..=1.0).contains(&v), "{p} {i}");
                last = v;
            }
        }
    }

    #[test]
    fn het_variant_in_diploid_region_is_not_biallelic() {
        // CN 2, minor allele CN 1, the variant on one copy
        let p = calculate_probability(0.6, 30, 60, 2.0, 1.0, 1.0);
        assert!(p < 0.01, "{p}");
        // the same in an LOH region (minor allele CN 0) with the variant on both copies
        let p = calculate_probability(0.6, 60, 60, 2.0, 2.0, 0.0);
        assert!(p > 0.95, "{p}");
    }

    /// PU10: a variant copy number of 0 in a deleted segment (CN <= 0) makes the threshold count -infinity.
    /// Java's long-to-int narrowing turns it into 0, so P(no wildtype | LOH) is 1 - exp(-mean); saturated, it is 1.
    #[test]
    fn threshold_count_saturates() {
        let java = crate::compat::with_java(true, || calculate_probability(0.5, 1, 10, -0.1, 0.0, 0.0));
        let ours = crate::compat::with_java(false, || calculate_probability(0.5, 1, 10, -0.1, 0.0, 0.0));
        assert!((java - (1.0 - (-1.0f64).exp())).abs() < 1e-6, "{java}");
        assert!((ours - 1.0).abs() < 1e-6, "{ours}");
        crate::compat::with_java(true, || assert_eq!(-1, threshold_read_count(i64::MAX)));
        crate::compat::with_java(true, || assert_eq!(0, threshold_read_count(i64::MIN)));
        crate::compat::with_java(false, || assert_eq!(i32::MAX, threshold_read_count(i64::MAX)));
        crate::compat::with_java(false, || assert_eq!(i32::MIN, threshold_read_count(-3_000_000_000)));
    }

    /// BiallelicModel.calculateProbability bit for bit on the synthetic inputs of the harness
    /// (testdata/biallelic_java.tsv, made by Java at upstream 883ab11bba): upstream's cases, edge cases and
    /// 5,000 random ones. In Java mode every result is Java's; by default only the PU10 cases differ.
    #[test]
    fn matches_java_vectors() {
        let text = include_str!("testdata/biallelic_java.tsv");
        let f = |h: &str| f64::from_bits(u64::from_str_radix(h, 16).unwrap());
        let mut n = 0;
        let mut corrected = 0;
        for line in text.lines() {
            let c: Vec<&str> = line.split('\t').collect();
            let (purity, ad, depth, cn, vcn, macn, expected) =
                (f(c[0]), c[1].parse::<i32>().unwrap(), c[2].parse::<i32>().unwrap(), f(c[3]), f(c[4]), f(c[5]), f(c[6]));
            let java = crate::compat::with_java(true, || calculate_probability(purity, ad, depth, cn, vcn, macn));
            assert_eq!(expected.to_bits(), java.to_bits(), "{line}: {java} != {expected}");
            let ours = crate::compat::with_java(false, || calculate_probability(purity, ad, depth, cn, vcn, macn));
            if ours.to_bits() != java.to_bits() {
                corrected += 1;
                // only where the threshold count is beyond the int range
                let threshold = jmin(cn - 0.5, jmax(BIALLELIC_VCN_THRESHOLD_PI, cn - BIALLELIC_VCN_THRESHOLD_P2));
                assert!((threshold / vcn * ad as f64).abs() > i32::MAX as f64 / 4.0, "{line}");
            }
            n += 1;
        }
        assert!(n > 5000, "{n}");
        assert!(corrected >= 2, "{corrected}");
    }
}
