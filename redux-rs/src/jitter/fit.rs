//! The jitter model fit (Java JitterModelFitter, JitterModelLoss, JitterModelCalc) and the Apache Commons Math
//! least squares regression it uses.

use super::table::{JitterCountsTable, JitterModelParams, JitterTableRow};
use anyhow::{bail, Result};

/// Doubles.EPSILON
const EPSILON: f64 = 1e-10;

/// Math.max / Math.min: NaN if either argument is NaN (f64::max returns the other argument).
#[inline]
pub fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() { f64::NAN } else if a > b { a } else { b }
}
#[inline]
pub fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() { f64::NAN } else if a < b { a } else { b }
}

/// Math.log10 as the JVM computes it where there is no intrinsic (fdlibm e_log10.c): the result differs from a
/// correctly rounded log10 in the last bit for some arguments.
pub fn java_log10(x: f64) -> f64 {
    const TWO54: f64 = 1.80143985094819840000e+16;
    const IVLN10: f64 = 4.34294481903251816668e-01;
    const LOG10_2HI: f64 = 3.01029995663611771306e-01;
    const LOG10_2LO: f64 = 3.69423907715893078616e-13;
    let mut x = x;
    let mut hx = (x.to_bits() >> 32) as i32;
    let lx = x.to_bits() as u32;
    let mut k: i32 = 0;
    if hx < 0x0010_0000 {
        if ((hx & 0x7fff_ffff) as u32 | lx) == 0 {
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            return f64::NAN;
        }
        k -= 54;
        x *= TWO54;
        hx = (x.to_bits() >> 32) as i32;
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    let i = ((k as u32) & 0x8000_0000) >> 31;
    hx = (hx & 0x000f_ffff) | ((0x3ff - i as i32) << 20);
    let y = (k + i as i32) as f64;
    x = f64::from_bits(((hx as u32 as u64) << 32) | (x.to_bits() & 0xffff_ffff));
    let z = y * LOG10_2LO + IVLN10 * fdlibm_log(x);
    z + y * LOG10_2HI
}

/// fdlibm e_log.c (__ieee754_log), which Math.log10 calls. The libm crate's log sums the same terms in another
/// order and differs in the last bit for about 2% of arguments, so the original is written out here.
fn fdlibm_log(x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const TWO54: f64 = 1.80143985094819840000e+16;
    const LG1: f64 = 6.666666666666735130e-01;
    const LG2: f64 = 3.999999999940941908e-01;
    const LG3: f64 = 2.857142874366239149e-01;
    const LG4: f64 = 2.222219843214978396e-01;
    const LG5: f64 = 1.818357216161805012e-01;
    const LG6: f64 = 1.531383769920937332e-01;
    const LG7: f64 = 1.479819860511658591e-01;
    let mut x = x;
    let mut hx = (x.to_bits() >> 32) as i32;
    let lx = x.to_bits() as u32;
    let mut k: i32 = 0;
    if hx < 0x0010_0000 {
        if ((hx & 0x7fff_ffff) as u32 | lx) == 0 {
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            return f64::NAN;
        }
        k -= 54;
        x *= TWO54;
        hx = (x.to_bits() >> 32) as i32;
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000f_ffff;
    let mut i = (hx + 0x95f64) & 0x10_0000;
    // normalize x or x/2
    x = f64::from_bits((((hx | (i ^ 0x3ff0_0000)) as u32 as u64) << 32) | (x.to_bits() & 0xffff_ffff));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000f_ffff & (2 + hx)) < 3 {
        // |f| < 2**-20
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            let dk = k as f64;
            return dk * LN2_HI + dk * LN2_LO;
        }
        let r = f * f * (0.5 - 0.33333333333333333 * f);
        if k == 0 {
            return f - r;
        }
        let dk = k as f64;
        return dk * LN2_HI - ((r - dk * LN2_LO) - f);
    }
    let s = f / (2.0 + f);
    let dk = k as f64;
    let z = s * s;
    i = hx - 0x6147a;
    let w = z * z;
    let j = 0x6b851 - hx;
    let t1 = w * (LG2 + w * (LG4 + w * LG6));
    let t2 = z * (LG1 + w * (LG3 + w * (LG5 + w * LG7)));
    i |= j;
    let r = t2 + t1;
    if i > 0 {
        let hfsq = 0.5 * f * f;
        if k == 0 { f - (hfsq - s * (hfsq + r)) } else { dk * LN2_HI - ((hfsq - (s * (hfsq + r) + dk * LN2_LO)) - f) }
    } else if k == 0 {
        f - s * (f - r)
    } else {
        dk * LN2_HI - ((s * (f - r) - dk * LN2_LO) - f)
    }
}

/// Math.exp (fdlibm e_exp.c; the libm crate is a port of the same code).
#[inline]
pub fn java_exp(x: f64) -> f64 {
    libm::exp(x)
}

/// Math.pow (fdlibm e_pow.c; the libm crate is a port of the same code).
#[inline]
pub fn java_pow(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

#[inline]
fn sgn(x: f64) -> f64 {
    if x >= 0.0 { 1.0 } else { -1.0 }
}

/// JitterModelCalc.rawAsymmetricLaplace
pub fn raw_asymmetric_laplace(x: f64, scale: f64, skew: f64) -> f64 {
    1.0 / (scale * (skew + 1.0 / skew)) * java_exp((-x * sgn(x) * java_pow(skew, sgn(x))) / scale)
}

/// JitterModelCalc.modifiedAsymmetricLaplace for the jitters -5..5: index 0 is jitter -5.
pub fn modified_asymmetric_laplace_all(scale: f64, raw_skew: f64) -> [f64; 11] {
    let skew = 1.0 + (raw_skew - 1.0) * java_min(scale, 1.0);
    let mut raw = [0.0; 11];
    let mut pdf_sum = 0.0;
    let mut x1 = -5.0;
    let mut i = 0;
    while x1 <= 5.0 {
        raw[i] = raw_asymmetric_laplace(x1, scale, skew);
        pdf_sum += raw[i];
        x1 += 1.0;
        i += 1;
    }
    let mut out = [0.0; 11];
    for i in 0..11 {
        out[i] = raw[i] / pdf_sum;
    }
    out
}

/// JitterModelCalc.modifiedAsymmetricLaplace
#[cfg(test)]
pub fn modified_asymmetric_laplace(x: i32, scale: f64, raw_skew: f64) -> f64 {
    modified_asymmetric_laplace_all(scale, raw_skew)[(x + 5) as usize]
}

/// JitterModelFitter.searchIncrement
fn search_increment(last: f64) -> f64 {
    if last - 0.25 < EPSILON {
        return last + 0.01;
    }
    if last - 0.5 < -EPSILON {
        return last + 0.02;
    }
    if last - 2.0 < EPSILON {
        return last + 0.05;
    }
    last + 0.2
}

/// The scale and skew values of JitterModelFitter.gridSearch with the model frequencies at every pair. Java
/// evaluates the model eleven times for each jitter of each pair of each table row; the values depend only on the
/// pair, so they are computed once for the whole run.
pub struct SearchGrid {
    /// The values the scale loop and the skew loop both step through, with their accumulated rounding.
    pub values: Vec<f64>,
    /// max(model frequency, 1e-4) for jitters -5..5, indexed [scale * values.len() + skew]
    model_freq: Vec<[f64; 11]>,
}

impl SearchGrid {
    pub fn new() -> Self {
        let mut values = Vec::new();
        let mut v = 0.05;
        while v <= 5.001 {
            values.push(v);
            v = search_increment(v);
        }
        let mut model_freq = Vec::with_capacity(values.len() * values.len());
        for &scale in &values {
            for &skew in &values {
                let mut f = modified_asymmetric_laplace_all(scale, skew);
                for x in f.iter_mut() {
                    *x = java_max(*x, 1e-4);
                }
                model_freq.push(f);
            }
        }
        SearchGrid { values, model_freq }
    }
}

/// JitterModelFitter.ScaleSkew
#[derive(Clone, Copy, Debug)]
pub struct ScaleSkew {
    pub scale: f64,
    pub skew: f64,
    #[allow(dead_code)]
    pub loss: f64,
}

/// JitterModelFitter.gridSearch with JitterModelLoss.totalLoss. `length_minus_one_scale` is the fitted scale of the
/// previous repeat count, or 0.
pub fn grid_search(grid: &SearchGrid, num_repeats: i32, row: &JitterTableRow, fallback_scale: f64, length_minus_one_scale: f64) -> ScaleSkew {
    let mut min_loss = f64::MAX;
    let mut best_scale = f64::NAN;
    let mut best_skew = f64::NAN;

    // JitterModelLoss.freq: NaN for every jitter when the row has no reads, which leaves the best values NaN
    let mut freq = [0.0; 11];
    let mut floored = [0.0; 11];
    for x in -5..=5 {
        let i = (x + 5) as usize;
        freq[i] = row.jitter_read_count(x) as f64 / row.total_read_count as f64;
        floored[i] = java_max(freq[i], 1e-4);
    }

    let n = grid.values.len();
    for (si, &scale) in grid.values.iter().enumerate() {
        for (ki, &skew) in grid.values.iter().enumerate() {
            let model = &grid.model_freq[si * n + ki];
            let mut loss = 0.0;
            for i in 0..11 {
                if i != 5 {
                    loss += -10.0 * java_log10(java_min(floored[i], model[i]) / java_max(floored[i], model[i])) * freq[i];
                }
            }
            loss += (java_max(skew, 1.0 / skew) - 1.0) * 0.01;
            if num_repeats >= 7 {
                // Math.pow(x, 2) is x * x
                let d = (scale - length_minus_one_scale).abs();
                loss += d * d;
            }
            if loss < min_loss {
                min_loss = loss;
                best_scale = scale;
                best_skew = skew;
            }
        }
    }
    let scale = if row.total_read_count < 20000 && !fallback_scale.is_nan() { fallback_scale } else { best_scale };
    ScaleSkew { scale, skew: best_skew, loss: min_loss }
}

/// JitterModelFitter.FittingDefaults: scale at 4, 5 and 6 repeats, and the slope.
fn defaults(repeat_unit_length: usize) -> (f64, f64, f64, f64) {
    match repeat_unit_length {
        1 => (0.1, 0.15, 0.2, 0.06),
        2 => (0.15, 0.2, 0.25, 0.06),
        _ => (0.2, 0.25, 0.3, 0.06),
    }
}

/// JitterModelFitter.performFit
pub fn perform_fit(grid: &SearchGrid, table: &JitterCountsTable) -> Result<JitterModelParams> {
    let (scale_four, scale_five, scale_six, slope) = defaults(table.repeat_unit_length());
    // fitted values for repeat counts 4..15, index 0 is 4
    let mut fitted: [Option<ScaleSkew>; 12] = [None; 12];
    let mut total_shortened: i32 = 0;
    let mut total_lengthened: i32 = 0;

    for num_repeats in 4..=15 {
        let default_scale = match num_repeats {
            4 => scale_four,
            5 => scale_five,
            6 => scale_six,
            _ => f64::NAN,
        };
        if let Some(row) = table.row(num_repeats) {
            for (&jitter, &count) in &row.jitter_counts {
                if jitter < 0 {
                    total_shortened = total_shortened.wrapping_add(count);
                }
                if jitter > 0 {
                    total_lengthened = total_lengthened.wrapping_add(count);
                }
            }
            let previous = if num_repeats > 4 { fitted[(num_repeats - 5) as usize].map(|s| s.scale).unwrap_or(0.0) } else { 0.0 };
            fitted[(num_repeats - 4) as usize] = Some(grid_search(grid, num_repeats, row, default_scale, previous));
        }
    }

    let mut read_count_sum_7_15 = 0.0;
    let mut num_repeat_counts_over_min = 0;
    for num_repeats in 7..=15 {
        let read_count = table.read_count(num_repeats);
        read_count_sum_7_15 += read_count as f64;
        if read_count >= 20000 {
            num_repeat_counts_over_min += 1;
        }
    }

    let mut scale_fit_gradient = f64::NAN;
    let mut scale_fit_intercept = f64::NAN;
    let mut microsatellite_skew = f64::NAN;

    if num_repeat_counts_over_min > 1 {
        // weighted least squares through the scales of repeat counts 7 to 15
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut weights = Vec::new();
        let mut skews = Vec::new();
        for num_repeats in 7..=15 {
            let Some(scale_skew) = fitted[(num_repeats - 4) as usize] else { continue };
            let w = java_min(table.read_count(num_repeats) as f64 / read_count_sum_7_15, 0.2);
            xs.push(num_repeats as f64);
            ys.push(scale_skew.scale);
            weights.push(w);
            skews.push(scale_skew.skew);
        }
        // two usable repeat counts pass this test and then fail in the regression (see ols_no_intercept)
        if xs.len() >= 2 {
            let y: Vec<f64> = (0..ys.len()).map(|i| weights[i] * ys[i]).collect();
            let x: Vec<[f64; 2]> = (0..xs.len()).map(|i| [weights[i] * xs[i], weights[i]]).collect();
            let params = ols_no_intercept(&y, &x)?;
            // the gradient must be positive
            if params[0] - 0.0 > EPSILON {
                scale_fit_gradient = params[0];
                scale_fit_intercept = params[1];
                microsatellite_skew = 0.0;
                let mut weight_sum = 0.0;
                for i in 0..skews.len() {
                    microsatellite_skew += weights[i] * skews[i];
                    weight_sum += weights[i];
                }
                microsatellite_skew /= weight_sum;
            }
        }
    }

    if scale_fit_gradient.is_nan() {
        scale_fit_gradient = slope;
        scale_fit_intercept = scale_six - scale_fit_gradient * 6.0;
    }
    if microsatellite_skew.is_nan() {
        if total_shortened.min(total_lengthened) < 50 {
            microsatellite_skew = 1.0;
        } else {
            microsatellite_skew = java_max(java_log10(total_shortened as f64 / total_lengthened as f64) + 1.0, 0.01);
        }
    }

    Ok(JitterModelParams {
        repeat_unit: table.repeat_unit.clone(),
        consensus_type: table.consensus_type,
        optimal_scale_repeat4: fitted[0].map(|s| s.scale).unwrap_or(scale_four),
        optimal_scale_repeat5: fitted[1].map(|s| s.scale).unwrap_or(scale_five),
        optimal_scale_repeat6: fitted[2].map(|s| s.scale).unwrap_or(scale_six),
        scale_fit_gradient,
        scale_fit_intercept,
        microsatellite_skew,
    })
}

/// OLSMultipleLinearRegression with setNoIntercept(true): newSampleData(y, x) then estimateRegressionParameters,
/// for two regressors. This is Commons Math 3.6.1 QRDecomposition (Householder reflections on the transposed
/// matrix, threshold 0) and its solver, operation for operation.
pub fn ols_no_intercept(y: &[f64], x: &[[f64; 2]]) -> Result<[f64; 2]> {
    let m = x.len();
    if m != y.len() {
        bail!("regression: {} observations but {} rows", y.len(), m);
    }
    // AbstractMultipleLinearRegression.validateSampleData requires more rows than predictors plus one, also
    // without an intercept: two rows are an exception in Java, which ends the run
    if m < 3 {
        bail!("org.apache.commons.math3.exception.MathIllegalArgumentException: not enough data ({m} rows) for this many predictors (2 predictors)");
    }
    // qrt[column][row]
    let mut qrt: [Vec<f64>; 2] = [x.iter().map(|r| r[0]).collect(), x.iter().map(|r| r[1]).collect()];
    let mut r_diag = [0.0; 2];

    // QRDecomposition.performHouseholderReflection
    for minor in 0..2 {
        let mut x_norm_sqr = 0.0;
        for row in minor..m {
            let c = qrt[minor][row];
            x_norm_sqr += c * c;
        }
        let a = if qrt[minor][minor] > 0.0 { -x_norm_sqr.sqrt() } else { x_norm_sqr.sqrt() };
        r_diag[minor] = a;
        if a != 0.0 {
            qrt[minor][minor] -= a;
            for col in minor + 1..2 {
                let (head, tail) = qrt.split_at_mut(col);
                let qrt_minor = &head[minor];
                let qrt_col = &mut tail[0];
                let mut alpha = 0.0;
                for row in minor..m {
                    alpha -= qrt_col[row] * qrt_minor[row];
                }
                alpha /= a * qrt_minor[minor];
                for row in minor..m {
                    qrt_col[row] -= alpha * qrt_minor[row];
                }
            }
        }
    }

    // QRDecomposition.Solver.isNonSingular; Java throws SingularMatrixException
    if r_diag.iter().any(|d| d.abs() <= 0.0) {
        bail!("org.apache.commons.math3.linear.SingularMatrixException: matrix is singular");
    }

    // QRDecomposition.Solver.solve(RealVector): apply the reflections to y, then back-substitute
    let mut yv = y.to_vec();
    for minor in 0..2 {
        let qrt_minor = &qrt[minor];
        let mut dot_product = 0.0;
        for row in minor..m {
            dot_product += yv[row] * qrt_minor[row];
        }
        dot_product /= r_diag[minor] * qrt_minor[minor];
        for row in minor..m {
            yv[row] += dot_product * qrt_minor[row];
        }
    }
    let mut out = [0.0; 2];
    for row in (0..2).rev() {
        yv[row] /= r_diag[row];
        let y_row = yv[row];
        out[row] = y_row;
        for i in 0..row {
            yv[i] -= y_row * qrt[row][i];
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jitter::read::ConsensusType;

    #[test]
    fn grid_values() {
        let grid = SearchGrid::new();
        // 0.05..0.25 by 0.01, ..0.5 by 0.02, ..2.0 by 0.05, then by 0.2 up to 5.001
        assert_eq!(grid.values.len(), 79);
        assert_eq!(grid.values[0], 0.05);
        assert_eq!(grid.values[1], 0.060000000000000005);
        assert!((grid.values[78] - 4.85).abs() < 1e-9);
        assert!(grid.values.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn laplace_sums_to_one() {
        for (scale, skew) in [(0.1, 1.0), (0.5, 0.7), (2.0, 3.0)] {
            let all = modified_asymmetric_laplace_all(scale, skew);
            let sum: f64 = all.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12);
            assert_eq!(modified_asymmetric_laplace(-2, scale, skew), all[3]);
        }
        // symmetric when the skew is 1
        let all = modified_asymmetric_laplace_all(0.4, 1.0);
        assert_eq!(all[4], all[6]);
    }

    #[test]
    fn java_math() {
        assert!(java_max(f64::NAN, 1.0).is_nan());
        assert!(java_min(1.0, f64::NAN).is_nan());
        assert_eq!(java_log10(1000.0), 3.0);
        assert_eq!(java_log10(1.0), 0.0);
        assert_eq!(java_log10(0.0), f64::NEG_INFINITY);
        assert!(java_log10(-1.0).is_nan());
        // values from the Java runtime (Double.doubleToLongBits of Math.log10, Math.exp, Math.pow)
        assert_eq!(java_log10(0.3).to_bits(), JAVA_LOG10_0_3);
        assert_eq!(java_log10(7.25e-5).to_bits(), JAVA_LOG10_7_25E_5);
        assert_eq!(java_exp(-3.7).to_bits(), JAVA_EXP_M3_7);
        assert_eq!(java_pow(0.0123, 3.0).to_bits(), JAVA_POW_0_0123_3);
    }
    const JAVA_LOG10_0_3: u64 = 0xbfe0bb6c34d81502;
    const JAVA_LOG10_7_25E_5: u64 = 0xc0108f038db91645;
    const JAVA_EXP_M3_7: u64 = 0x3f99511fc6871044;
    const JAVA_POW_0_0123_3: u64 = 0x3ebf385ce77c93ac;

    #[test]
    fn regression_recovers_line() {
        // y = 0.05 x - 0.1 with weights
        let w = [0.2, 0.2, 0.15, 0.1, 0.05];
        let xs = [7.0, 8.0, 9.0, 10.0, 11.0];
        let y: Vec<f64> = (0..5).map(|i| w[i] * (0.05 * xs[i] - 0.1)).collect();
        let x: Vec<[f64; 2]> = (0..5).map(|i| [w[i] * xs[i], w[i]]).collect();
        let p = ols_no_intercept(&y, &x).unwrap();
        assert!((p[0] - 0.05).abs() < 1e-12);
        assert!((p[1] + 0.1).abs() < 1e-12);
        // singular, and too few rows
        assert!(ols_no_intercept(&[0.0, 0.0, 0.0], &[[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]).is_err());
        assert!(ols_no_intercept(&[0.1, 0.2], &[[1.0, 0.1], [2.0, 0.2]]).is_err());
    }

    fn table_with(rows: &[(i32, &[(i32, i32)])]) -> JitterCountsTable {
        let mut t = JitterCountsTable::new("A/T", ConsensusType::None);
        for (units, counts) in rows {
            let row = t.get_or_create_row(*units);
            for (j, c) in counts.iter() {
                row.add_reads(*j, *c);
            }
        }
        t
    }

    #[test]
    fn fit_defaults_for_low_counts() {
        let grid = SearchGrid::new();
        // too few reads everywhere: default scales, default slope, skew 1
        let t = table_with(&[(4, &[(0, 1000), (-1, 3), (1, 2)]), (8, &[(0, 500), (-1, 20)])]);
        let p = perform_fit(&grid, &t).unwrap();
        assert_eq!(p.optimal_scale_repeat4, 0.1);
        assert_eq!(p.optimal_scale_repeat5, 0.15);
        assert_eq!(p.optimal_scale_repeat6, 0.2);
        assert_eq!(p.scale_fit_gradient, 0.06);
        assert_eq!(p.scale_fit_intercept, 0.2 - 0.06 * 6.0);
        assert_eq!(p.microsatellite_skew, 1.0);
        // an empty table gives the same
        let p = perform_fit(&grid, &JitterCountsTable::new("AT/TA", ConsensusType::Dual)).unwrap();
        assert_eq!((p.optimal_scale_repeat4, p.optimal_scale_repeat5, p.optimal_scale_repeat6), (0.15, 0.2, 0.25));
        assert_eq!(p.scale_fit_intercept, 0.25 - 0.06 * 6.0);
    }

    #[test]
    fn fit_skew_from_counts() {
        let grid = SearchGrid::new();
        // enough shortened and lengthened reads for the log ratio, not enough rows for the regression
        let t = table_with(&[(4, &[(0, 30000), (-1, 400), (1, 100)])]);
        let p = perform_fit(&grid, &t).unwrap();
        assert_eq!(p.microsatellite_skew, java_log10(4.0) + 1.0);
        // 30000 reads: the grid search scale is used instead of the default
        assert_ne!(p.optimal_scale_repeat4, 0.1);
        assert!(grid.values.contains(&p.optimal_scale_repeat4));
    }

    #[test]
    fn empty_row_gives_nan_then_fallback() {
        let grid = SearchGrid::new();
        let row = JitterTableRow::new(8);
        let s = grid_search(&grid, 8, &row, f64::NAN, 0.0);
        assert!(s.scale.is_nan() && s.skew.is_nan());
        assert_eq!(s.loss, f64::MAX);
        let s = grid_search(&grid, 4, &row, 0.1, 0.0);
        assert_eq!(s.scale, 0.1);
    }

    #[test]
    fn grid_search_finds_generating_parameters() {
        let grid = SearchGrid::new();
        // reads distributed by the model at scale 0.3, skew 1.5
        let (scale, skew) = (grid.values[23], grid.values[55]);
        let freq = modified_asymmetric_laplace_all(scale, skew);
        let mut row = JitterTableRow::new(5);
        for x in -5..=5 {
            row.add_reads(x, (freq[(x + 5) as usize] * 1e7).round() as i32);
        }
        let s = grid_search(&grid, 5, &row, f64::NAN, 0.0);
        assert_eq!(s.scale, scale);
        assert_eq!(s.skew, skew);
    }
}
