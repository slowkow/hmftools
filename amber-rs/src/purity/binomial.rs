//! BinomialDistribution(n, p).cumulativeProbability(x) as commons-math3 3.6.1 computes it (Beta.regularizedBeta,
//! Beta.logBeta, Gamma.logGamma, ContinuedFraction; Apache License 2.0), so that the capture tests of the
//! candidate peaks, which compare it against 0.16 and 0.84, come out as in Java. Taken from sage-rs (stats.rs),
//! which checked it against commons-math3. exp and log are FastMath's (hmf_common::java_fastmath), log1p
//! the correctly rounded one: on the 5,111,983 values of the peak search's grid (`java_grid`) the result is
//! bit-identical to Java's.

#![allow(clippy::excessive_precision)]


/// FastMath's exp and log (ported), and a correctly rounded log1p (FastMath.log1p is accurate to well under
/// one ulp, so correct rounding reproduces it far more often than the platform libm does). This matters
/// where a result is subtracted from 1: BinomialDistribution(54, 0.5).cumulativeProbability(0) is exactly
/// 0 in Java.
trait CorrectlyRounded {
    fn cr_ln(self) -> f64;
    fn cr_exp(self) -> f64;
    fn cr_ln_1p(self) -> f64;
}

impl CorrectlyRounded for f64 {
    fn cr_ln(self) -> f64 {
        hmf_common::java_fastmath::log(self)
    }
    fn cr_exp(self) -> f64 {
        hmf_common::java_fastmath::exp(self)
    }
    fn cr_ln_1p(self) -> f64 {
        core_math::log1p(self)
    }
}

const HALF_LOG_TWO_PI: f64 = 0.918_938_533_204_672_7;

// --- Gamma ---

const LANCZOS_G: f64 = 607.0 / 128.0;
const LANCZOS: [f64; 15] = [
    0.999_999_999_999_997_091_82,
    57.156_235_665_862_923_517,
    -59.597_960_355_475_491_248,
    14.136_097_974_741_747_174,
    -0.491_913_816_097_620_199_78,
    0.339_946_499_848_118_886_99e-4,
    0.465_236_289_270_485_756_65e-4,
    -0.983_744_753_048_795_646_77e-4,
    0.158_088_703_224_912_488_84e-3,
    -0.210_264_441_724_104_883_19e-3,
    0.217_439_618_115_212_643_20e-3,
    -0.164_318_106_536_763_890_22e-3,
    0.844_182_239_838_527_432_93e-4,
    -0.261_908_384_015_814_086_70e-4,
    0.368_991_826_595_316_227_04e-5,
];

fn lanczos(x: f64) -> f64 {
    let mut sum = 0.0;
    for i in (1..LANCZOS.len()).rev() {
        sum += LANCZOS[i] / (x + i as f64);
    }
    sum + LANCZOS[0]
}

/// Gamma.invGamma1pm1, for -0.5 <= x <= 1.5
fn inv_gamma_1p_m1(x: f64) -> f64 {
    const A0: f64 = 0.611_609_510_448_141_581_788E-08;
    const A1: f64 = 0.624_730_830_116_465_516_210E-08;
    const B: [f64; 8] = [
        0.203_610_414_066_806_987_300E+00,
        0.266_205_348_428_949_217_746E-01,
        0.493_944_979_382_446_875_238E-03,
        -0.851_419_432_440_314_906_588E-05,
        -0.643_045_481_779_353_022_248E-05,
        0.992_641_840_672_773_722_196E-06,
        -0.607_761_895_722_825_260_739E-07,
        0.195_755_836_614_639_731_882E-09,
    ];
    const P: [f64; 7] = [
        0.611_609_510_448_141_581_786_1E-08,
        0.687_167_411_306_719_873_615_2E-08,
        0.682_016_166_849_617_065_791_8E-09,
        0.468_684_332_294_884_803_108_0E-10,
        0.157_283_302_771_044_628_699_5E-11,
        -0.124_944_157_227_636_621_322_2E-12,
        0.434_352_993_740_859_425_517_8E-14,
    ];
    const Q: [f64; 4] = [
        0.305_696_107_836_522_102_500_9E+00,
        0.546_421_308_604_229_653_601_6E-01,
        0.495_683_009_382_588_731_202_0E-02,
        0.269_236_946_618_636_119_287_6E-03,
    ];
    const C: f64 = -0.422_784_335_098_467_139_393_487_909_917_598E+00;
    const C0: f64 = 0.577_215_664_901_532_860_606_512_090_082_402E+00;
    // C1..C13
    const CS: [f64; 13] = [
        -0.655_878_071_520_253_881_077_019_515_145_390E+00,
        -0.420_026_350_340_952_355_290_039_348_754_298E-01,
        0.166_538_611_382_291_489_501_700_795_102_105E+00,
        -0.421_977_345_555_443_367_482_083_012_891_874E-01,
        -0.962_197_152_787_697_356_211_492_167_234_820E-02,
        0.721_894_324_666_309_954_239_501_034_044_657E-02,
        -0.116_516_759_185_906_511_211_397_108_401_839E-02,
        -0.215_241_674_114_950_972_815_729_963_053_648E-03,
        0.128_050_282_388_116_186_153_198_626_328_164E-03,
        -0.201_348_547_807_882_386_556_893_914_210_218E-04,
        -0.125_049_348_214_267_065_734_535_947_383_309E-05,
        0.113_302_723_198_169_588_237_412_962_033_074E-05,
        -0.205_633_841_697_760_710_345_015_413_002_057E-06,
    ];
    let t = if x <= 0.5 { x } else { (x - 0.5) - 0.5 };
    // c = C13 + first; then C12..C1 by Horner
    let horner = |first: f64| {
        let mut c = CS[12] + first;
        for i in (0..12).rev() {
            c = CS[i] + t * c;
        }
        c
    };
    if t < 0.0 {
        let a = A0 + t * A1;
        let mut b = B[7];
        for i in (0..7).rev() {
            b = B[i] + t * b;
        }
        b = 1.0 + t * b;
        let c = C + t * horner(t * (a / b));
        if x > 0.5 {
            t * c / x
        } else {
            x * ((c + 0.5) + 0.5)
        }
    } else {
        let mut p = P[6];
        for i in (0..6).rev() {
            p = P[i] + t * p;
        }
        let mut q = Q[3];
        for i in (0..3).rev() {
            q = Q[i] + t * q;
        }
        q = 1.0 + t * q;
        let c = C0 + t * horner((p / q) * t);
        if x > 0.5 {
            (t / x) * ((c - 0.5) - 0.5)
        } else {
            x * c
        }
    }
}

fn log_gamma_1p(x: f64) -> f64 {
    -inv_gamma_1p_m1(x).cr_ln_1p()
}

/// Gamma.logGamma
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        f64::NAN
    } else if x < 0.5 {
        log_gamma_1p(x) - x.cr_ln()
    } else if x <= 2.5 {
        log_gamma_1p((x - 0.5) - 0.5)
    } else if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut prod = 1.0;
        for i in 1..=n {
            prod *= x - i as f64;
        }
        log_gamma_1p(x - (n + 1) as f64) + prod.cr_ln()
    } else {
        let sum = lanczos(x);
        let tmp = x + LANCZOS_G + 0.5;
        ((x + 0.5) * tmp.cr_ln()) - tmp + HALF_LOG_TWO_PI + (sum / x).cr_ln()
    }
}

// --- Beta.logBeta ---

const DELTA: [f64; 15] = [
    0.833_333_333_333_333_333_333_333_333_333E-01,
    -0.277_777_777_777_777_777_777_777_752_282E-04,
    0.793_650_793_650_793_650_791_732_130_419E-07,
    -0.595_238_095_238_095_232_389_839_236_182E-09,
    0.841_750_841_750_832_853_294_451_671_990E-11,
    -0.191_752_691_751_854_612_334_149_171_243E-12,
    0.641_025_640_510_325_475_730_918_472_625E-14,
    -0.295_506_514_125_338_232_839_867_823_991E-15,
    0.179_643_716_359_402_238_723_287_696_452E-16,
    -0.139_228_964_661_627_791_231_203_060_395E-17,
    0.133_802_855_014_020_915_603_275_339_093E-18,
    -0.154_246_009_867_966_094_273_710_216_533E-19,
    0.197_701_992_980_957_427_278_370_133_333E-20,
    -0.234_065_664_793_997_056_856_992_426_667E-21,
    0.171_348_014_966_398_575_409_015_466_667E-22,
];

/// for 1 <= a, b <= 2
fn log_gamma_sum(a: f64, b: f64) -> f64 {
    let x = (a - 1.0) + (b - 1.0);
    if x <= 0.5 {
        log_gamma_1p(1.0 + x)
    } else if x <= 1.5 {
        log_gamma_1p(x) + x.cr_ln_1p()
    } else {
        log_gamma_1p(x - 1.0) + (x * (1.0 + x)).cr_ln()
    }
}

/// for 0 <= a <= b, b >= 10
fn delta_minus_delta_sum(a: f64, b: f64) -> f64 {
    let h = a / b;
    let p = h / (1.0 + h);
    let q = 1.0 / (1.0 + h);
    let q2 = q * q;
    let mut s = [0.0; 15];
    s[0] = 1.0;
    for i in 1..s.len() {
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }
    let sqrt_t = 10.0 / b;
    let t = sqrt_t * sqrt_t;
    let mut w = DELTA[14] * s[14];
    for i in (0..14).rev() {
        w = t * w + DELTA[i] * s[i];
    }
    w * p / b
}

/// for a >= 0, b >= 10
fn log_gamma_minus_log_gamma_sum(a: f64, b: f64) -> f64 {
    let (d, w) = if a <= b { (b + (a - 0.5), delta_minus_delta_sum(a, b)) } else { (a + (b - 0.5), delta_minus_delta_sum(b, a)) };
    let u = d * (a / b).cr_ln_1p();
    let v = a * (b.cr_ln() - 1.0);
    if u <= v {
        (w - u) - v
    } else {
        (w - v) - u
    }
}

/// for p, q >= 10
fn sum_delta_minus_delta_sum(p: f64, q: f64) -> f64 {
    let a = p.min(q);
    let b = p.max(q);
    let sqrt_t = 10.0 / a;
    let t = sqrt_t * sqrt_t;
    let mut z = DELTA[14];
    for i in (0..14).rev() {
        z = t * z + DELTA[i];
    }
    z / a + delta_minus_delta_sum(a, b)
}

/// Beta.logBeta
pub fn log_beta(p: f64, q: f64) -> f64 {
    if p.is_nan() || q.is_nan() || p <= 0.0 || q <= 0.0 {
        return f64::NAN;
    }
    let a = p.min(q);
    let b = p.max(q);
    if a >= 10.0 {
        let w = sum_delta_minus_delta_sum(a, b);
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -(a - 0.5) * c.cr_ln();
        let v = b * h.cr_ln_1p();
        if u <= v {
            (((-0.5 * b.cr_ln() + HALF_LOG_TWO_PI) + w) - u) - v
        } else {
            (((-0.5 * b.cr_ln() + HALF_LOG_TWO_PI) + w) - v) - u
        }
    } else if a > 2.0 {
        if b > 1000.0 {
            let n = (a - 1.0).floor() as i32;
            let mut prod = 1.0;
            let mut ared = a;
            for _ in 0..n {
                ared -= 1.0;
                prod *= ared / (1.0 + ared / b);
            }
            (prod.cr_ln() - n as f64 * b.cr_ln()) + (log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b))
        } else {
            let mut prod1 = 1.0;
            let mut ared = a;
            while ared > 2.0 {
                ared -= 1.0;
                let h = ared / b;
                prod1 *= h / (1.0 + h);
            }
            if b < 10.0 {
                let mut prod2 = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod2 *= bred / (ared + bred);
                }
                prod1.cr_ln() + prod2.cr_ln() + (log_gamma(ared) + (log_gamma(bred) - log_gamma_sum(ared, bred)))
            } else {
                prod1.cr_ln() + log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b)
            }
        }
    } else if a >= 1.0 {
        if b > 2.0 {
            if b < 10.0 {
                let mut prod = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod *= bred / (a + bred);
                }
                prod.cr_ln() + (log_gamma(a) + (log_gamma(bred) - log_gamma_sum(a, bred)))
            } else {
                log_gamma(a) + log_gamma_minus_log_gamma_sum(a, b)
            }
        } else {
            log_gamma(a) + log_gamma(b) - log_gamma_sum(a, b)
        }
    } else if b >= 10.0 {
        log_gamma(a) + log_gamma_minus_log_gamma_sum(a, b)
    } else {
        // Java uses Gamma.gamma here; not reachable from the binomial CDF, where both arguments are >= 1
        libm::lgamma(a) + libm::lgamma(b) - libm::lgamma(a + b)
    }
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
    (a * x.cr_ln() + b * (-x).cr_ln_1p() - a.cr_ln() - log_beta(a, b)).cr_exp() * 1.0 / h_n
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
        // n = 0: every x >= 0 is >= n
        assert_eq!(binomial_cdf(0, 0.2, 2), 1.0);
    }

    /// (n, p, x, cumulative probability bits) from commons-math3 3.6.1 (taken from sage-rs)
    const JAVA_CDF: [(i32, f64, i32, u64); 21] = [
        (30, 0.003, 0, 0x3fed3de95d237a21),
        (30, 0.05, 0, 0x3fcb794873cb5d1c),
        (30, 0.3, 0, 0x3ef7a25babb68000),
        (30, 0.75, 0, 0x0),
        (30, 0.003, 3, 0x3feffffba056ee01),
        (30, 0.05, 3, 0x3fee0e28cc74dde7),
        (30, 0.3, 3, 0x3f8314907f4cc6c0),
        (30, 0.75, 3, 0x3d3bc00000000000),
        (100, 0.003, 7, 0x3fefffffff7cb289),
        (100, 0.05, 7, 0x3febe7bf6d24a38c),
        (100, 0.3, 7, 0x3e51c1f320000000),
        (100, 0.75, 7, 0x0),
        (45, 0.05, 10, 0x3fefffeab52cb521),
        (45, 0.3, 10, 0x3fc515df03322004),
        (45, 0.75, 10, 0x3d47940000000000),
        (12, 0.3, 12, 0x3ff0000000000000),
        (200, 0.003, 1, 0x3fec1b026601795e),
        (200, 0.05, 1, 0x3f3a7a777a63c800),
        (200, 0.3, 1, 0x0),
        (200, 0.75, 1, 0x0),
        (200, 1.0, 1, 0x0),
    ];

    /// The Java grid of BinomialGrid (work/amber/tumoronly/harness), in AMBER_BINOMIAL_GRID with the grid's sizes
    /// in AMBER_BINOMIAL_GRID_N ("maxN,maxDepthN"): counts the values that differ from Java and those on the other
    /// side of the 0.16 / 0.84 bounds the peak search compares them with
    #[test]
    #[ignore]
    fn java_grid() {
        let path = std::env::var("AMBER_BINOMIAL_GRID").expect("AMBER_BINOMIAL_GRID");
        let sizes = std::env::var("AMBER_BINOMIAL_GRID_N").expect("AMBER_BINOMIAL_GRID_N");
        let (max_n, max_depth_n) = sizes.split_once(',').map(|(a, b)| (a.parse::<i32>().unwrap(), b.parse::<i32>().unwrap())).unwrap();
        let bytes = std::fs::read(path).unwrap();
        let mut values = bytes.chunks_exact(8).map(|c| f64::from_bits(u64::from_be_bytes(c.try_into().unwrap())));
        let (mut total, mut differ, mut max_ulps, mut crossings) = (0u64, 0u64, 0u64, 0u64);
        let mut check = |java: f64, rust: f64, bounds: &[f64]| {
            total += 1;
            if java.to_bits() != rust.to_bits() {
                differ += 1;
                max_ulps = max_ulps.max(java.to_bits().abs_diff(rust.to_bits()));
                if bounds.iter().any(|&b| (java < b) != (rust < b) || (java > b) != (rust > b)) {
                    crossings += 1;
                }
            }
        };
        for (level, _) in crate::purity::peak::search_values_and_steps() {
            for p in [level, level / 2.0] {
                for n in 0..=max_n {
                    for k in 0..=n / 2 {
                        check(values.next().unwrap(), binomial_cdf(n, p, k), &[0.16, 0.84]);
                    }
                }
            }
            for n in 0..=max_depth_n {
                check(values.next().unwrap(), binomial_cdf(n, level / 2.0, 2), &[0.16]);
            }
        }
        assert!(values.next().is_none());
        println!("{total} values: {differ} differ from Java (at most {max_ulps} ulps), {crossings} on the other side of a bound");
        assert_eq!(crossings, 0);
    }

    #[test]
    fn cdf_matches_commons_math() {
        for (n, p, x, cdf) in JAVA_CDF {
            let got = binomial_cdf(n, p, x);
            assert_eq!(got.to_bits(), cdf, "cdf({n}, {p}, {x}) = {got:e}, Java {:e}", f64::from_bits(cdf));
        }
    }
}
