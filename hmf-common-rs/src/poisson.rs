//! (Moved here from amber-rs/src/poisson.rs for LINX; the same code.)
//! PoissonDistribution of commons-math3 3.6.1 (Apache License 2.0): probability (SaddlePointExpansion)
//! and cumulativeProbability (Gamma.regularizedGammaQ / P with the ContinuedFraction), with Gamma.logGamma.
//! The Gamma and saddle-point code follows sage-rs/src/stats.rs. exp and log are FastMath's (ported in
//! crate::java_fastmath); log1p is core-math's correctly rounded one
//! (FastMath.log1p needs FastMath's high-precision log, not ported), used only by logGamma below 8.
//! Tests compare with values from Java.

#[inline]
fn ln(x: f64) -> f64 {
    crate::java_fastmath::log(x)
}
#[inline]
fn exp(x: f64) -> f64 {
    crate::java_fastmath::exp(x)
}
#[inline]
fn ln_1p(x: f64) -> f64 {
    core_math::log1p(x)
}

const HALF_LOG_TWO_PI: f64 = 0.918_938_533_204_672_7;

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
        if x > 0.5 { t * c / x } else { x * ((c + 0.5) + 0.5) }
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
        if x > 0.5 { (t / x) * ((c - 0.5) - 0.5) } else { x * c }
    }
}

fn log_gamma_1p(x: f64) -> f64 {
    -ln_1p(inv_gamma_1p_m1(x))
}

/// Gamma.logGamma
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        f64::NAN
    } else if x < 0.5 {
        log_gamma_1p(x) - ln(x)
    } else if x <= 2.5 {
        log_gamma_1p((x - 0.5) - 0.5)
    } else if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut prod = 1.0;
        for i in 1..=n {
            prod *= x - i as f64;
        }
        log_gamma_1p(x - (n + 1) as f64) + ln(prod)
    } else {
        let sum = lanczos(x);
        let tmp = x + LANCZOS_G + 0.5;
        ((x + 0.5) * ln(tmp)) - tmp + HALF_LOG_TWO_PI + ln(sum / x)
    }
}

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

/// SaddlePointExpansion.getStirlingError
fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            EXACT_STIRLING_ERRORS[z2 as usize]
        } else {
            log_gamma(z + 1.0) - (z + 0.5) * ln(z) + z - HALF_LOG_TWO_PI
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

/// SaddlePointExpansion.getDeviancePart
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
        x * ln(x / mu) + mu - x
    }
}

const EPSILON: f64 = 1e-12;
const MAX_ITERATIONS: i32 = 10_000_000;

/// PoissonDistribution(mean) with the default accuracy (1e-12, 10 million iterations). Java throws
/// NotStrictlyPositiveException for a mean that is not above 0; callers make sure of it.
#[derive(Clone, Copy, Debug)]
pub struct Poisson {
    mean: f64,
}

impl Poisson {
    pub fn new(mean: f64) -> Poisson {
        debug_assert!(mean > 0.0);
        Poisson { mean }
    }
    /// logProbability
    pub fn log_probability(&self, x: i32) -> f64 {
        if x < 0 || x == i32::MAX {
            f64::NEG_INFINITY
        } else if x == 0 {
            -self.mean
        } else {
            let xf = x as f64;
            -stirling_error(xf) - deviance_part(xf, self.mean) - 0.5 * ln(std::f64::consts::TAU) - 0.5 * ln(xf)
        }
    }
    /// probability
    pub fn probability(&self, x: i32) -> f64 {
        let lp = self.log_probability(x);
        if lp == f64::NEG_INFINITY { 0.0 } else { exp(lp) }
    }
    /// cumulativeProbability
    pub fn cumulative_probability(&self, x: i32) -> f64 {
        if x < 0 {
            return 0.0;
        }
        if x == i32::MAX {
            return 1.0;
        }
        regularized_gamma_q(x as f64 + 1.0, self.mean)
    }
}

fn precision_equals(x: f64, y: f64, eps: f64) -> bool {
    (y - x).abs() <= eps
}

/// Gamma.regularizedGammaP(a, x, 1e-12, 10 million)
pub fn regularized_gamma_p(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        f64::NAN
    } else if x == 0.0 {
        0.0
    } else if x >= a + 1.0 {
        1.0 - regularized_gamma_q(a, x)
    } else {
        let mut n = 0.0;
        let mut an = 1.0 / a;
        let mut sum = an;
        while (an / sum).abs() > EPSILON && n < MAX_ITERATIONS as f64 && sum < f64::INFINITY {
            n += 1.0;
            an *= x / (a + n);
            sum += an;
        }
        if sum.is_infinite() { 1.0 } else { exp(-x + (a * ln(x)) - log_gamma(a)) * sum }
    }
}

/// Gamma.regularizedGammaQ(a, x, 1e-12, 10 million), with ContinuedFraction.evaluate
pub fn regularized_gamma_q(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        f64::NAN
    } else if x == 0.0 {
        1.0
    } else if x < a + 1.0 {
        1.0 - regularized_gamma_p(a, x)
    } else {
        let get_a = |n: i32| ((2.0 * n as f64) + 1.0) - a + x;
        let get_b = |n: i32| n as f64 * (a - n as f64);
        let small = 1e-50;
        let mut h_prev = get_a(0);
        if precision_equals(h_prev, 0.0, small) {
            h_prev = small;
        }
        let mut n = 1;
        let mut d_prev = 0.0;
        let mut c_prev = h_prev;
        let mut h_n = h_prev;
        while n < MAX_ITERATIONS {
            let a_n = get_a(n);
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
            let delta_n = c_n * d_n;
            h_n = h_prev * delta_n;
            if h_n.is_infinite() || h_n.is_nan() {
                return f64::NAN;
            }
            if (delta_n - 1.0).abs() < EPSILON {
                break;
            }
            d_prev = d_n;
            c_prev = c_n;
            h_prev = h_n;
            n += 1;
        }
        let ret = 1.0 / h_n;
        exp(-x + (a * ln(x)) - log_gamma(a)) * ret
    }
}

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
