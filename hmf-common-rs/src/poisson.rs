//! commons-math3's PoissonDistribution (probability, cumulativeProbability) with Gamma.logGamma and the
//! regularized gamma functions: now in `commons_math` (with PURPLE's additions), re-exported here for
//! AMBER and LINX.

pub use crate::commons_math::{Poisson, log_gamma, regularized_gamma_p, regularized_gamma_q};

#[cfg(test)]
mod tests {
    use super::*;

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
