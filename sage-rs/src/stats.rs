//! Binomial distribution, following commons-math3 3.6.1 (BinomialDistribution, Beta.regularizedBeta,
//! ContinuedFraction). logBeta uses lgamma rather than the commons-math logGamma series, so results
//! can differ from Java in the last few bits.

fn log_beta(a: f64, b: f64) -> f64 {
    libm::lgamma(a) + libm::lgamma(b) - libm::lgamma(a + b)
}

fn precision_equals(x: f64, y: f64, eps: f64) -> bool {
    (y - x).abs() <= eps
}

pub fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    const EPSILON: f64 = 1e-14;
    if x.is_nan() || a.is_nan() || b.is_nan() || x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }
    if x > (a + 1.0) / (2.0 + b + a) && 1.0 - x <= (b + 1.0) / (2.0 + b + a) {
        return 1.0 - regularized_beta(1.0 - x, b, a);
    }
    let get_b = |n: i32| -> f64 {
        if n % 2 == 0 {
            let m = n as f64 / 2.0;
            (m * (b - m) * x) / ((a + (2.0 * m) - 1.0) * (a + (2.0 * m)))
        } else {
            let m = (n as f64 - 1.0) / 2.0;
            -((a + m) * (a + b + m) * x) / ((a + (2.0 * m)) * (a + (2.0 * m) + 1.0))
        }
    };
    // ContinuedFraction.evaluate, with getA(n, x) = 1
    let small = 1e-50;
    let mut h_prev = 1.0;
    let mut n = 1;
    let mut d_prev = 0.0;
    let mut c_prev = h_prev;
    let mut h_n = h_prev;
    while n < i32::MAX {
        let a_n = 1.0;
        let b_n = get_b(n);
        let mut d_n = a_n + b_n * d_prev;
        if precision_equals(d_n, 0.0, small) {
            d_n = small;
        }
        let mut c_n = a_n + b_n / c_prev;
        if precision_equals(c_n, 0.0, small) {
            c_n = small;
        }
        d_n = 1.0 / d_n;
        let delta = c_n * d_n;
        h_n = h_prev * delta;
        if h_n.is_infinite() || h_n.is_nan() {
            return f64::NAN;
        }
        if (delta - 1.0).abs() < EPSILON {
            break;
        }
        d_prev = d_n;
        c_prev = c_n;
        h_prev = h_n;
        n += 1;
    }
    (a * x.ln() + b * (-x).ln_1p() - a.ln() - log_beta(a, b)).exp() * 1.0 / h_n
}

/// BinomialDistribution(n, p).cumulativeProbability(x)
pub fn binomial_cdf(n: i32, p: f64, x: i32) -> f64 {
    if x < 0 {
        0.0
    } else if x >= n {
        1.0
    } else {
        1.0 - regularized_beta(p, x as f64 + 1.0, (n - x) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cdf_matches_direct_sum() {
        // P(X <= 3), X ~ Binom(10, 0.3) = 0.6496107184
        assert!((binomial_cdf(10, 0.3, 3) - 0.6496107184).abs() < 1e-9);
        assert!((binomial_cdf(100, 0.5, 40) - 0.02844396682).abs() < 1e-9);
        assert_eq!(binomial_cdf(5, 0.2, 5), 1.0);
        assert_eq!(binomial_cdf(5, 0.2, -1), 0.0);
    }
}
