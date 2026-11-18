//! Binomial distribution, following commons-math3 3.6.1 (BinomialDistribution, SaddlePointExpansion,
//! Beta.regularizedBeta / logBeta, Gamma.logGamma, ContinuedFraction; Apache License 2.0), so that
//! probabilities compared against thresholds come out as they do in Java. The remaining difference is
//! Java's FastMath exp/log/log1p against the platform libm, which is at the level of the last bit.

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
    -inv_gamma_1p_m1(x).ln_1p()
}

/// Gamma.logGamma
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        f64::NAN
    } else if x < 0.5 {
        log_gamma_1p(x) - x.ln()
    } else if x <= 2.5 {
        log_gamma_1p((x - 0.5) - 0.5)
    } else if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut prod = 1.0;
        for i in 1..=n {
            prod *= x - i as f64;
        }
        log_gamma_1p(x - (n + 1) as f64) + prod.ln()
    } else {
        let sum = lanczos(x);
        let tmp = x + LANCZOS_G + 0.5;
        ((x + 0.5) * tmp.ln()) - tmp + HALF_LOG_TWO_PI + (sum / x).ln()
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
        log_gamma_1p(x) + x.ln_1p()
    } else {
        log_gamma_1p(x - 1.0) + (x * (1.0 + x)).ln()
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
    let u = d * (a / b).ln_1p();
    let v = a * (b.ln() - 1.0);
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
        let u = -(a - 0.5) * c.ln();
        let v = b * h.ln_1p();
        if u <= v {
            (((-0.5 * b.ln() + HALF_LOG_TWO_PI) + w) - u) - v
        } else {
            (((-0.5 * b.ln() + HALF_LOG_TWO_PI) + w) - v) - u
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
            (prod.ln() - n as f64 * b.ln()) + (log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b))
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
                prod1.ln() + prod2.ln() + (log_gamma(ared) + (log_gamma(bred) - log_gamma_sum(ared, bred)))
            } else {
                prod1.ln() + log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b)
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
                prod.ln() + (log_gamma(a) + (log_gamma(bred) - log_gamma_sum(a, bred)))
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

// --- SaddlePointExpansion (commons-math3), used by BinomialDistribution.probability ---

const EXACT_STIRLING_ERRORS: [f64; 31] = [
    0.0,
    0.153_426_409_720_027_345_291_384_8,
    0.081_061_466_795_327_258_219_670_2,
    0.054_814_121_051_917_653_896_139_0,
    0.041_340_695_955_409_294_093_822_1,
    0.033_162_873_519_936_287_485_110_48,
    0.027_677_925_684_998_339_148_789_29,
    0.023_746_163_656_297_495_971_329_20,
    0.020_790_672_103_765_093_111_522_77,
    0.018_488_450_532_673_185_230_779_34,
    0.016_644_691_189_821_192_163_194_87,
    0.015_134_973_221_917_378_873_512_55,
    0.013_876_128_823_070_747_998_745_73,
    0.012_810_465_242_920_226_924_249_86,
    0.011_896_709_945_891_770_095_055_72,
    0.011_104_559_758_206_917_326_629_91,
    0.010_411_265_261_972_096_497_478_567,
    0.009_799_416_126_158_803_298_389_475,
    0.009_255_462_182_712_732_917_728_637,
    0.008_768_700_134_139_385_462_952_823,
    0.008_330_563_433_362_871_256_469_318,
    0.007_934_114_564_314_020_547_248_100,
    0.007_573_675_487_951_840_794_972_024,
    0.007_244_554_301_320_383_179_543_912,
    0.006_942_840_107_209_529_865_664_152,
    0.006_665_247_032_707_682_442_354_394,
    0.006_408_994_188_004_207_068_439_631,
    0.006_171_712_263_039_457_647_532_867,
    0.005_951_370_112_758_847_735_624_416,
    0.005_746_216_513_010_115_682_023_589,
    0.005_554_733_551_962_801_371_038_690,
];

fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            EXACT_STIRLING_ERRORS[z2 as usize]
        } else {
            log_gamma(z + 1.0) - (z + 0.5) * z.ln() + z - HALF_LOG_TWO_PI
        }
    } else {
        let z2 = z * z;
        (0.083_333_333_333_333_333_333
            - (0.002_777_777_777_777_777_777_78
                - (0.000_793_650_793_650_793_650_793_65
                    - (0.000_595_238_095_238_095_238_095_238 - 0.000_841_750_841_750_841_750_841_750_8 / z2) / z2)
                    / z2)
                / z2)
            / z
    }
}

fn deviance_part(x: f64, mu: f64) -> f64 {
    if (x - mu).abs() < 0.1 * (x + mu) {
        let d = x - mu;
        let mut v = d / (x + mu);
        let mut s1 = v * d;
        let mut s = f64::NAN;
        let mut ej = 2.0 * x * v;
        v *= v;
        let mut j = 1;
        while s1 != s {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1) as f64;
            j += 1;
        }
        s1
    } else {
        x * (x / mu).ln() + mu - x
    }
}

fn log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    let (xf, nf) = (x as f64, n as f64);
    if x == 0 {
        if p < 0.1 {
            -deviance_part(nf, nf * q) - nf * p
        } else {
            nf * q.ln()
        }
    } else if x == n {
        if q < 0.1 {
            -deviance_part(nf, nf * p) - nf * q
        } else {
            nf * p.ln()
        }
    } else {
        let ret = stirling_error(nf) - stirling_error(xf) - stirling_error(nf - xf)
            - deviance_part(xf, nf * p)
            - deviance_part(nf - xf, nf * q);
        let f = (2.0 * std::f64::consts::PI * xf * (nf - xf)) / nf;
        -0.5 * f.ln() + ret
    }
}

/// BinomialDistribution(n, p).probability(x)
pub fn binomial_probability(n: i32, p: f64, x: i32) -> f64 {
    let log_p = if n == 0 {
        if x == 0 { 0.0 } else { f64::NEG_INFINITY }
    } else if x < 0 || x > n {
        f64::NEG_INFINITY
    } else {
        log_binomial_probability(x, n, p, 1.0 - p)
    };
    if log_p == f64::NEG_INFINITY { 0.0 } else { log_p.exp() }
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

    /// (n, p, x, probability bits, cumulative probability bits) from commons-math3 3.6.1
    const JAVA_BINOMIAL: [(i32, f64, i32, u64, u64); 42] = [
        (30, 0.0, 0, 0x3ff0000000000000, 0x3ff0000000000000),
        (30, 0.003, 0, 0x3fed3de95d237a21, 0x3fed3de95d237a21),
        (30, 0.05, 0, 0x3fcb794873cb5d2c, 0x3fcb794873cb5d1c),
        (30, 0.3, 0, 0x3ef7a25babb65698, 0x3ef7a25babb68000),
        (30, 0.75, 0, 0x3c30000000000008, 0x0),
        (30, 1.0, 0, 0x0, 0x0),
        (30, 0.0, 3, 0x0, 0x3ff0000000000000),
        (30, 0.003, 3, 0x3f1a7f448abc7f11, 0x3feffffba056ee01),
        (30, 0.05, 3, 0x3fc0432981574525, 0x3fee0e28cc74dde7),
        (30, 0.3, 3, 0x3f7d814d0df558aa, 0x3f8314907f4cc6c0),
        (30, 0.75, 3, 0x3d3ac3400000001d, 0x3d3bc00000000000),
        (30, 1.0, 3, 0x0, 0x0),
        (30, 0.0, 30, 0x0, 0x3ff0000000000000),
        (30, 0.003, 30, 0x3037d72f39816293, 0x3ff0000000000000),
        (30, 0.05, 30, 0x37d4484bfeebc2b2, 0x3ff0000000000000),
        (30, 0.3, 30, 0x3cadac0ac33ce14c, 0x3ff0000000000000),
        (30, 0.75, 30, 0x3f276838794f1724, 0x3ff0000000000000),
        (30, 1.0, 30, 0x3ff0000000000000, 0x3ff0000000000000),
        (100, 0.0, 7, 0x0, 0x3ff0000000000000),
        (100, 0.003, 7, 0x3e5c6d2fba32099c, 0x3fefffffff7cb289),
        (100, 0.05, 7, 0x3fbb247d5785c8d3, 0x3febe7bf6d24a38c),
        (100, 0.3, 7, 0x3e4d86f34d0d36d2, 0x3e51c1f320000000),
        (100, 0.75, 7, 0x363fd70f521e1f96, 0x0),
        (100, 1.0, 7, 0x0, 0x0),
        (45, 0.0, 10, 0x0, 0x3ff0000000000000),
        (45, 0.003, 10, 0x3ca8702d2905caa3, 0x3ff0000000000000),
        (45, 0.05, 10, 0x3f0b20ae4b5a4ef6, 0x3fefffeab52cb521),
        (45, 0.3, 10, 0x3fb244b4b4d5e245, 0x3fc515df03322004),
        (45, 0.75, 10, 0x3d456a80c54ff0ba, 0x3d47940000000000),
        (45, 1.0, 10, 0x0, 0x0),
        (12, 0.0, 12, 0x0, 0x3ff0000000000000),
        (12, 0.003, 12, 0x39a58ecc82b7b8b3, 0x3ff0000000000000),
        (12, 0.05, 12, 0x3cb19799812dea11, 0x3ff0000000000000),
        (12, 0.3, 12, 0x3ea1d50b1e323882, 0x3ff0000000000000),
        (12, 0.75, 12, 0x3fa037e200000001, 0x3ff0000000000000),
        (12, 1.0, 12, 0x3ff0000000000000, 0x3ff0000000000000),
        (200, 0.0, 1, 0x0, 0x3ff0000000000000),
        (200, 0.003, 1, 0x3fd51e64de337129, 0x3fec1b026601795e),
        (200, 0.05, 1, 0x3f382e616c43be13, 0x3f3a7a777a63c800),
        (200, 0.3, 1, 0x39e6bc1a9bb9f13e, 0x0),
        (200, 0.75, 1, 0x2782bfffffffff75, 0x0),
        (200, 1.0, 1, 0x0, 0x0),
    ];

    fn ulps(a: f64, b: u64) -> u64 {
        a.to_bits().abs_diff(b)
    }

    #[test]
    fn binomial_matches_commons_math() {
        // exp/log/lgamma are not the same implementations as Java's FastMath, so allow a few ulps
        for (n, p, x, prob, cdf) in JAVA_BINOMIAL {
            let got = binomial_probability(n, p, x);
            assert!(ulps(got, prob) <= 8, "probability({n}, {p}, {x}) = {got:e}, Java {:e}", f64::from_bits(prob));
            let got = binomial_cdf(n, p, x);
            assert!(ulps(got, cdf) <= 8, "cdf({n}, {p}, {x}) = {got:e}, Java {:e}", f64::from_bits(cdf));
        }
    }
}
