//! commons-math3's PoissonDistribution (probability, cumulativeProbability) with Gamma.logGamma and the
//! regularized gamma functions: now in `commons_math` (with PURPLE's additions), re-exported here for
//! AMBER and LINX; and hmf-common's PoissonCalcs (NEO).

pub use crate::commons_math::{Poisson, log_gamma, regularized_gamma_p, regularized_gamma_q};

/// PoissonCalcs.calcPoissonNoiseValue(observedCount, requiredProb): the mean of a Poisson distribution at
/// which the observed count falls at the required cumulative probability, by bisection (NEO, ISOFOX)
pub fn calc_poisson_noise_value(observed_count: i32, required_prob: f64) -> f64 {
    let threshold = if required_prob < 0.01 { required_prob * 0.1 } else { MIN_PROB_DIFF_PERC_HIGH };
    calc_poisson_noise_value_with(observed_count, required_prob, threshold, DEFAULT_MAX_ITERATIONS)
}

const DEFAULT_MAX_ITERATIONS: i32 = 20;
const MIN_PROB_DIFF_PERC_HIGH: f64 = 0.005;

/// PoissonCalcs.calcPoissonNoiseValue(observedCount, requiredProb, probDiffThreshold, maxIterations)
pub fn calc_poisson_noise_value_with(observed_count: i32, required_prob: f64, prob_diff_threshold: f64, max_iterations: i32) -> f64 {
    if observed_count < 0 {
        return 0.0;
    }
    if observed_count == 0 && required_prob > 0.5 {
        return 0.0;
    }
    let mut iterations = 0;
    let mut current_value;
    let mut test_value_upper;
    let mut test_value_lower;
    if required_prob > 0.5 {
        current_value = observed_count as f64 * 0.5;
        test_value_upper = observed_count as f64;
        test_value_lower = current_value * 0.01;
    } else {
        test_value_upper = find_upper_start_value(observed_count, required_prob, max_iterations);
        test_value_lower = observed_count as f64;
        current_value = (test_value_upper + test_value_lower) * 0.5;
    }
    let mut current_prob = Poisson::new(current_value).cumulative_probability(observed_count);
    let mut prob_diff = 0.0;
    while iterations < max_iterations {
        prob_diff = (required_prob - current_prob).abs();
        if prob_diff <= prob_diff_threshold {
            return current_value;
        }
        if current_prob < required_prob {
            test_value_upper = current_value;
            current_value = (current_value + test_value_lower) * 0.5;
        } else {
            test_value_lower = current_value;
            current_value = (current_value + test_value_upper) * 0.5;
        }
        current_prob = Poisson::new(current_value).cumulative_probability(observed_count);
        iterations += 1;
    }
    if iterations >= max_iterations {
        crate::warn!(
            "max iterations reached: value({}) test({}) prob({} diff={})",
            observed_count,
            crate::fmt::java_fixed(current_value, 1),
            crate::fmt::java_fixed(current_prob, 4),
            crate::fmt::java_fixed(prob_diff, 4)
        );
    }
    current_value
}

/// PoissonCalcs.findUpperStartValue
pub fn find_upper_start_value(observed_count: i32, required_prob: f64, max_iterations: i32) -> f64 {
    let mut current_value = (observed_count as f64 * 2.0).max(1.0);
    for _ in 0..max_iterations {
        if Poisson::new(current_value).cumulative_probability(observed_count) < required_prob {
            break;
        }
        current_value *= 2.0;
    }
    current_value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poisson_noise_values() {
        // PoissonCalcs.calcPoissonNoiseValue as NEO's TpmCalcsTest checks it (the means at which 10, 1 and 2
        // fragments fall at the 25% and 75% probabilities)
        let near = |a: f64, b: f64| assert!((a - b).abs() < 0.1, "{a} vs {b}");
        near(calc_poisson_noise_value(10, 0.25), 13.0);
        near(calc_poisson_noise_value(10, 0.75), 8.6);
        near(calc_poisson_noise_value(1, 0.25), 2.7);
        near(calc_poisson_noise_value(1, 0.75), 0.876);
        near(calc_poisson_noise_value(2, 0.25), 3.9);
        near(calc_poisson_noise_value(2, 0.75), 1.7);
        assert_eq!(calc_poisson_noise_value(-1, 0.25), 0.0);
        assert_eq!(calc_poisson_noise_value(0, 0.75), 0.0);
        // the upper start value doubles until the probability falls below the required one
        assert_eq!(find_upper_start_value(0, 0.25, 20), 2.0);
        assert_eq!(find_upper_start_value(10, 0.25, 20), 20.0);
    }

    #[test]
    fn probabilities_sum_to_one() {
        for mean in [0.05, 0.5, 3.0, 17.5, 60.0] {
            let p = Poisson::new(mean);
            let total: f64 = (0..400).map(|k| p.probability(k)).sum();
            assert!((total - 1.0).abs() < 1e-9, "{mean} {total}");
            let mut cdf = 0.0;
            for k in 0..60 {
                cdf += p.probability(k);
                assert!((cdf - p.cumulative_probability(k)).abs() < 1e-9, "{mean} {k}");
            }
        }
        assert_eq!(Poisson::new(2.0).probability(-1), 0.0);
        assert_eq!(Poisson::new(2.0).probability(0), (-2.0f64).exp());
    }
}
