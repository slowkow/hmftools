//! Java numeric behaviour LILAC depends on: Math.min/max on doubles, Math.round, and StrictMath.log10 (the
//! fdlibm algorithm, which the JVM uses for Math.log10 where it has no intrinsic; the allele frequency
//! penalty is a sum of these logarithms).

/// Math.min(double, double)
pub fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    if a == 0.0 && b == 0.0 {
        return if a.is_sign_negative() || b.is_sign_negative() { -0.0 } else { 0.0 };
    }
    if a <= b { a } else { b }
}

/// Math.max(double, double)
pub fn java_max(a: f64, b: f64) -> f64 {
    hmf_common::java_math::java_max(a, b)
}

/// Math.round(double) cast to int, as `(int) round(x)`
pub fn round_int(x: f64) -> i32 {
    hmf_common::java_math::java_round(x) as i32
}

#[inline]
fn hi(x: f64) -> i32 {
    (x.to_bits() >> 32) as i32
}
#[inline]
fn lo(x: f64) -> u32 {
    x.to_bits() as u32
}
#[inline]
fn with_hi(x: f64, h: i32) -> f64 {
    f64::from_bits(((h as u32 as u64) << 32) | lo(x) as u64)
}

const TWO54: f64 = 1.80143985094819840000e+16;

/// fdlibm __ieee754_log
pub fn fdlibm_log(mut x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const LG1: f64 = 6.666666666666735130e-01;
    const LG2: f64 = 3.999999999940941908e-01;
    const LG3: f64 = 2.857142874366239149e-01;
    const LG4: f64 = 2.222219843214978396e-01;
    const LG5: f64 = 1.818357216161805012e-01;
    const LG6: f64 = 1.531383769920937332e-01;
    const LG7: f64 = 1.479819860511658591e-01;
    let mut hx = hi(x);
    let lx = lo(x);
    let mut k: i32 = 0;
    if hx < 0x00100000 {
        if ((hx & 0x7fffffff) as u32 | lx) == 0 {
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            return f64::NAN;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= 0x7ff00000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000fffff;
    let i = (hx + 0x95f64) & 0x100000;
    x = with_hi(x, hx | (i ^ 0x3ff00000));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000fffff & (2 + hx)) < 3 {
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
    let mut i = hx - 0x6147a;
    let w = z * z;
    let j = 0x6b851 - hx;
    let t1 = w * (LG2 + w * (LG4 + w * LG6));
    let t2 = z * (LG1 + w * (LG3 + w * (LG5 + w * LG7)));
    i |= j;
    let r = t2 + t1;
    if i > 0 {
        let hfsq = 0.5 * f * f;
        if k == 0 {
            f - (hfsq - s * (hfsq + r))
        } else {
            dk * LN2_HI - ((hfsq - (s * (hfsq + r) + dk * LN2_LO)) - f)
        }
    } else if k == 0 {
        f - s * (f - r)
    } else {
        dk * LN2_HI - ((s * (f - r) - dk * LN2_LO) - f)
    }
}

/// StrictMath.log10 (fdlibm __ieee754_log10)
pub fn log10(mut x: f64) -> f64 {
    const IVLN10: f64 = 4.34294481903251816668e-01;
    const LOG10_2HI: f64 = 3.01029995663611771306e-01;
    const LOG10_2LO: f64 = 3.69423907715893078616e-13;
    let mut hx = hi(x);
    let lx = lo(x);
    let mut k: i32 = 0;
    if hx < 0x00100000 {
        if ((hx & 0x7fffffff) as u32 | lx) == 0 {
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            return f64::NAN;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= 0x7ff00000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    let i = ((k as u32) & 0x80000000) >> 31;
    hx = (hx & 0x000fffff) | ((0x3ff - i as i32) << 20);
    let y = (k + i as i32) as f64;
    x = with_hi(x, hx);
    let z = y * LOG10_2LO + IVLN10 * fdlibm_log(x);
    z + y * LOG10_2HI
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn logs() {
        assert_eq!(log10(1000.0), 3.0);
        assert_eq!(log10(0.0001), -4.0);
        assert_eq!(log10(1.0), 0.0);
        assert!(log10(0.0).is_infinite());
        for v in [0.06167487625085214, 0.0017203552857809527, 0.5, 2.5, 1e-300] {
            assert!((log10(v) - v.log10()).abs() < 1e-15);
            assert!((fdlibm_log(v) - v.ln()).abs() < 1e-12);
        }
        assert_eq!(round_int(2.5), 3);
        assert_eq!(round_int(-2.5), -2);
        assert_eq!(java_min(0.7, 0.5), 0.5);
    }
}

#[cfg(test)]
mod java_log10 {
    /// Math.log10 of the JVM on the frequencies of the cohort file, from /tmp/lilac_hash/log10_java.txt
    /// (written by a Java harness: input bits and result bits, hex); skipped without the file
    #[test]
    fn matches_java_math_log10() {
        let Ok(text) = std::fs::read_to_string("/tmp/lilac_hash/log10_java.txt") else { return };
        let (mut n, mut port_diff, mut std_diff) = (0, 0, 0);
        for l in text.lines() {
            let mut it = l.split(' ');
            let x = f64::from_bits(u64::from_str_radix(it.next().unwrap(), 16).unwrap());
            let y = f64::from_bits(u64::from_str_radix(it.next().unwrap(), 16).unwrap());
            n += 1;
            if super::log10(x) != y {
                port_diff += 1;
            }
            if x.log10() != y {
                std_diff += 1;
            }
        }
        eprintln!("log10 values {n}: port differs {port_diff}, std differs {std_diff}");
        assert_eq!(port_diff, 0);
    }
}
