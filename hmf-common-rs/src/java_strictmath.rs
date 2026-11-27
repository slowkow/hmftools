//! java.lang.StrictMath (the JDK's FdLibm.java, a port of fdlibm 5.3): `exp`, `log`, `log10` and `pow`.
//!
//! On the arm64 JVM of the reference machine `Math.exp`, `Math.log`, `Math.log10` and `Math.pow` give
//! exactly these results (checked on two million random arguments each, after JIT compilation), while
//! the platform's libm and Rust's `f64` methods differ from them in the last bit for some arguments. A
//! port whose decisions depend on such values (PURPLE's purity fit compares sums of logs and powers)
//! uses these. (On x86-64, HotSpot's intrinsics for these functions can differ from StrictMath.)

#[inline]
fn hi(x: f64) -> i32 {
    (x.to_bits() >> 32) as i32
}
#[inline]
fn lo(x: f64) -> i32 {
    x.to_bits() as i32
}
#[inline]
fn with_hi(x: f64, high: i32) -> f64 {
    f64::from_bits((x.to_bits() & 0x0000_0000_FFFF_FFFF) | ((high as u32 as u64) << 32))
}
#[inline]
fn with_lo(x: f64, low: i32) -> f64 {
    f64::from_bits((x.to_bits() & 0xFFFF_FFFF_0000_0000) | (low as u32 as u64))
}

const TWO54: f64 = 18014398509481984.0; // 0x1.0p54
const EXP_BITS: i32 = 0x7ff0_0000;
const EXP_SIGNIF_BITS: i32 = 0x7fff_ffff;

/// StrictMath.exp
pub fn exp(x: f64) -> f64 {
    const HALF: [f64; 2] = [0.5, -0.5];
    const HUGE: f64 = 1.0e+300;
    let twom1000: f64 = f64::from_bits(0x0170000000000000); // 2^-1000
    const O_THRESHOLD: f64 = 7.09782712893383973096e+02;
    const U_THRESHOLD: f64 = -7.45133219101941108420e+02;
    const LN2_HI: [f64; 2] = [6.93147180369123816490e-01, -6.93147180369123816490e-01];
    const LN2_LO: [f64; 2] = [1.90821492927058770002e-10, -1.90821492927058770002e-10];
    const INVLN2: f64 = 1.44269504088896338700e+00;
    const P1: f64 = 1.66666666666666019037e-01;
    const P2: f64 = -2.77777777770155933842e-03;
    const P3: f64 = 6.61375632143793436117e-05;
    const P4: f64 = -1.65339022054652515390e-06;
    const P5: f64 = 4.13813679705723846039e-08;

    let mut x = x;
    let mut hi_v = 0.0;
    let mut lo_v = 0.0;
    let mut k: i32 = 0;
    let mut hx = hi(x);
    let xsb = ((hx >> 31) & 1) as usize;
    hx &= EXP_SIGNIF_BITS;
    if hx >= 0x40862E42 {
        if hx >= 0x7ff00000 {
            if ((hx & 0xfffff) | lo(x)) != 0 {
                return x + x;
            } else {
                return if xsb == 0 { x } else { 0.0 };
            }
        }
        if x > O_THRESHOLD {
            return HUGE * HUGE;
        }
        if x < U_THRESHOLD {
            return twom1000 * twom1000;
        }
    }
    if hx > 0x3fd62e42 {
        if hx < 0x3FF0A2B2 {
            hi_v = x - LN2_HI[xsb];
            lo_v = LN2_LO[xsb];
            k = 1 - xsb as i32 - xsb as i32;
        } else {
            k = (INVLN2 * x + HALF[xsb]) as i32;
            let t = k as f64;
            hi_v = x - t * LN2_HI[0];
            lo_v = t * LN2_LO[0];
        }
        x = hi_v - lo_v;
    } else if hx < 0x3e300000 {
        if HUGE + x > 1.0 {
            return 1.0 + x;
        }
    } else {
        k = 0;
    }
    let t = x * x;
    let c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    if k == 0 {
        return 1.0 - ((x * c) / (c - 2.0) - x);
    }
    let y = 1.0 - ((lo_v - (x * c) / (2.0 - c)) - hi_v);
    if k >= -1021 {
        with_hi(y, hi(y).wrapping_add(k << 20))
    } else {
        with_hi(y, hi(y).wrapping_add((k + 1000) << 20)) * twom1000
    }
}

/// StrictMath.log
pub fn log(x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const LG1: f64 = 6.666666666666735130e-01;
    const LG2: f64 = 3.999999999940941908e-01;
    const LG3: f64 = 2.857142874366239149e-01;
    const LG4: f64 = 2.222219843214978396e-01;
    const LG5: f64 = 1.818357216161805012e-01;
    const LG6: f64 = 1.531383769920937332e-01;
    const LG7: f64 = 1.479819860511658591e-01;

    let mut x = x;
    let mut hx = hi(x);
    let lx = lo(x);
    let mut k: i32 = 0;
    if hx < 0x0010_0000 {
        if ((hx & EXP_SIGNIF_BITS) | lx) == 0 {
            return -TWO54 / 0.0;
        }
        if hx < 0 {
            return (x - x) / 0.0;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= EXP_BITS {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000f_ffff;
    let i = (hx + 0x9_5f64) & 0x10_0000;
    x = with_hi(x, hx | (i ^ 0x3ff0_0000));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000f_ffff & (2 + hx)) < 3 {
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
    let mut i = hx - 0x6_147a;
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

/// StrictMath.log10
pub fn log10(x: f64) -> f64 {
    const IVLN10: f64 = 4.34294481903251816668e-01;
    const LOG10_2HI: f64 = 3.01029995663611771306e-01;
    const LOG10_2LO: f64 = 3.69423907715893078616e-13;
    let mut x = x;
    let mut hx = hi(x);
    let lx = lo(x);
    let mut k: i32 = 0;
    if hx < 0x0010_0000 {
        if ((hx & EXP_SIGNIF_BITS) | lx) == 0 {
            return -TWO54 / 0.0;
        }
        if hx < 0 {
            return (x - x) / 0.0;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= EXP_BITS {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    let i = ((k as u32) & 0x8000_0000) >> 31;
    let i = i as i32;
    hx = (hx & 0x000f_ffff) | ((0x3ff - i) << 20);
    let y = (k + i) as f64;
    x = with_hi(x, hx);
    let z = y * LOG10_2LO + IVLN10 * log(x);
    z + y * LOG10_2HI
}

/// Math.scalb(double, int)
fn scalb(d: f64, scale_factor: i32) -> f64 {
    const MAX_SCALE: i32 = 1023 + 1022 + 53 + 1;
    let (mut scale_factor, scale_increment, exp_delta) = if scale_factor < 0 {
        (scale_factor.max(-MAX_SCALE), -512, f64::from_bits(((1023 - 512) as u64) << 52))
    } else {
        (scale_factor.min(MAX_SCALE), 512, f64::from_bits(((1023 + 512) as u64) << 52))
    };
    let t = ((scale_factor >> 8) as u32 >> 23) as i32;
    let exp_adjust = ((scale_factor + t) & 511) - t;
    let mut d = d * f64::from_bits(((exp_adjust + 1023) as u64) << 52);
    scale_factor -= exp_adjust;
    while scale_factor != 0 {
        d *= exp_delta;
        scale_factor -= scale_increment;
    }
    d
}

/// StrictMath.pow
pub fn pow(x: f64, y: f64) -> f64 {
    if y == 0.0 {
        return 1.0;
    }
    if x.is_nan() || y.is_nan() {
        return x + y;
    }
    let y_abs = y.abs();
    let mut x_abs = x.abs();
    if y == 2.0 {
        return x * x;
    } else if y == 0.5 {
        if x >= -f64::MAX {
            return (x + 0.0).sqrt();
        }
    } else if y_abs == 1.0 {
        return if y == 1.0 { x } else { 1.0 / x };
    } else if y_abs == f64::INFINITY {
        if x_abs == 1.0 {
            return y - y;
        } else if x_abs > 1.0 {
            return if y >= 0.0 { y } else { 0.0 };
        } else {
            return if y < 0.0 { -y } else { 0.0 };
        }
    }
    let hx = hi(x);
    let mut ix = hx & EXP_SIGNIF_BITS;
    let mut y_is_int: i32 = 0;
    if hx < 0 {
        if y_abs >= 9007199254740992.0 {
            y_is_int = 2;
        } else if y_abs >= 1.0 {
            let y_abs_as_long = y_abs as i64;
            if (y_abs_as_long as f64) == y_abs {
                y_is_int = 2 - (y_abs_as_long & 0x1) as i32;
            }
        }
    }
    if x_abs == 0.0 || x_abs == f64::INFINITY || x_abs == 1.0 {
        let mut z = x_abs;
        if y < 0.0 {
            z = 1.0 / z;
        }
        if hx < 0 {
            if ((ix - 0x3ff00000) | y_is_int) == 0 {
                z = (z - z) / (z - z);
            } else if y_is_int == 1 {
                z = -1.0 * z;
            }
        }
        return z;
    }
    let mut n: i32 = (hx >> 31) + 1;
    if (n | y_is_int) == 0 {
        return (x - x) / (x - x);
    }
    let mut s = 1.0;
    if (n | (y_is_int - 1)) == 0 {
        s = -1.0;
    }
    let t1: f64;
    let t2: f64;
    if y_abs > f64::from_bits(0x41e00000ffffffff) {
        const INV_LN2: f64 = 1.44269504088896338700e+00;
        const INV_LN2_H: f64 = 1.44269502162933349609e+00;
        const INV_LN2_L: f64 = 1.92596299112661746887e-08;
        if x_abs < f64::from_bits(0x3fefffff00000000) {
            return if y < 0.0 { s * f64::INFINITY } else { s * 0.0 };
        }
        if x_abs > f64::from_bits(0x3ff00000ffffffff) {
            return if y > 0.0 { s * f64::INFINITY } else { s * 0.0 };
        }
        let t = x_abs - 1.0;
        let w = (t * t) * (0.5 - t * (0.3333333333333333333333 - t * 0.25));
        let u = INV_LN2_H * t;
        let v = t * INV_LN2_L - w * INV_LN2;
        let tt1 = with_lo(u + v, 0);
        t1 = tt1;
        t2 = v - (tt1 - u);
    } else {
        const CP: f64 = 9.61796693925975554329e-01;
        const CP_H: f64 = 9.61796700954437255859e-01;
        const CP_L: f64 = -7.02846165095275826516e-09;
        n = 0;
        if ix < 0x00100000 {
            x_abs *= 9007199254740992.0;
            n -= 53;
            ix = hi(x_abs);
        }
        n += (ix >> 20) - 0x3ff;
        let j = ix & 0x000fffff;
        ix = j | 0x3ff00000;
        let k: usize;
        if j <= 0x3988E {
            k = 0;
        } else if j < 0xBB67A {
            k = 1;
        } else {
            k = 0;
            n += 1;
            ix -= 0x00100000;
        }
        x_abs = with_hi(x_abs, ix);
        const BP: [f64; 2] = [1.0, 1.5];
        const DP_H: [f64; 2] = [0.0, 5.84962487220764160156e-01];
        const DP_L: [f64; 2] = [0.0, 1.35003920212974897128e-08];
        const L1: f64 = 5.99999999999994648725e-01;
        const L2: f64 = 4.28571428578550184252e-01;
        const L3: f64 = 3.33333329818377432918e-01;
        const L4: f64 = 2.72728123808534006489e-01;
        const L5: f64 = 2.30660745775561754067e-01;
        const L6: f64 = 2.06975017800338417784e-01;
        let mut u = x_abs - BP[k];
        let mut v = 1.0 / (x_abs + BP[k]);
        let ss = u * v;
        let s_h = with_lo(ss, 0);
        let mut t_h = with_hi(0.0, ((ix >> 1) | 0x20000000) + 0x00080000 + ((k as i32) << 18));
        let mut t_l = x_abs - (t_h - BP[k]);
        let s_l = v * ((u - s_h * t_h) - s_h * t_l);
        let mut s2 = ss * ss;
        let mut r = s2 * s2 * (L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
        r += s_l * (s_h + ss);
        s2 = s_h * s_h;
        t_h = with_lo(3.0 + s2 + r, 0);
        t_l = r - ((t_h - 3.0) - s2);
        u = s_h * t_h;
        v = s_l * t_h + t_l * ss;
        let p_h = with_lo(u + v, 0);
        let p_l = v - (p_h - u);
        let z_h = CP_H * p_h;
        let z_l = CP_L * p_h + p_l * CP + DP_L[k];
        let t = n as f64;
        let tt1 = with_lo(((z_h + z_l) + DP_H[k]) + t, 0);
        t1 = tt1;
        t2 = z_l - (((tt1 - t) - DP_H[k]) - z_h);
    }
    let y1 = with_lo(y, 0);
    let p_l = (y - y1) * t1 + y * t2;
    let mut p_h = y1 * t1;
    let mut z = p_l + p_h;
    let mut j = hi(z);
    let i = lo(z);
    if j >= 0x40900000 {
        if ((j - 0x40900000) | i) != 0 {
            return s * f64::INFINITY;
        } else {
            const OVT: f64 = 8.0085662595372944372e-0017;
            if p_l + OVT > z - p_h {
                return s * f64::INFINITY;
            }
        }
    } else if (j & EXP_SIGNIF_BITS) >= 0x4090cc00 {
        if ((j.wrapping_sub(0xc090cc00u32 as i32)) | i) != 0 {
            return s * 0.0;
        } else if p_l <= z - p_h {
            return s * 0.0;
        }
    }
    const P1: f64 = 1.66666666666666019037e-01;
    const P2: f64 = -2.77777777770155933842e-03;
    const P3: f64 = 6.61375632143793436117e-05;
    const P4: f64 = -1.65339022054652515390e-06;
    const P5: f64 = 4.13813679705723846039e-08;
    const LG2: f64 = 6.93147180559945286227e-01;
    const LG2_H: f64 = 6.93147182464599609375e-01;
    const LG2_L: f64 = -1.90465429995776804525e-09;
    let i = j & EXP_SIGNIF_BITS;
    let mut k = (i >> 20) - 0x3ff;
    let mut n: i32 = 0;
    if i > 0x3fe00000 {
        n = j + (0x00100000 >> (k + 1));
        k = ((n & EXP_SIGNIF_BITS) >> 20) - 0x3ff;
        let t = with_hi(0.0, n & !(0x000fffff >> k));
        n = ((n & 0x000fffff) | 0x00100000) >> (20 - k);
        if j < 0 {
            n = -n;
        }
        p_h -= t;
    }
    let t = with_lo(p_l + p_h, 0);
    let u = t * LG2_H;
    let v = (p_l - (t - p_h)) * LG2 + t * LG2_L;
    z = u + v;
    let w = v - (z - u);
    let t = z * z;
    let t1 = z - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    let r = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    j = hi(z);
    j = j.wrapping_add(n << 20);
    if (j >> 20) <= 0 {
        z = scalb(z, n);
    } else {
        let z_hi = hi(z).wrapping_add(n << 20);
        z = with_hi(z, z_hi);
    }
    s * z
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_values() {
        assert_eq!(exp(0.0), 1.0);
        assert_eq!(log(1.0), 0.0);
        assert_eq!(pow(2.0, 10.0), 1024.0);
        assert_eq!(pow(2.0, -1.0), 0.5);
        assert!(pow(-8.0, 1.0 / 3.0).is_nan());
        assert_eq!(pow(-2.0, 3.0), -8.0);
    }

    /// 20,000 argument sets with Java's Math.exp, log, log10 and pow results (as raw bits), written by
    /// a Java program on the arm64 reference JVM
    #[test]
    fn java_vectors() {
        let text = include_str!("testdata/strictmath_vectors.txt");
        let b = |s: &str| f64::from_bits(u64::from_str_radix(s, 16).unwrap());
        let mut bad = Vec::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split(' ').collect();
            let (x, y, p) = (b(f[0]), b(f[1]), b(f[2]));
            let checks = [("exp", exp(x), b(f[3])), ("log", log(y), b(f[4])), ("log10", log10(y), b(f[5])), ("pow", pow(y, p), b(f[6])), ("powx", pow(x, p), b(f[7]))];
            for (name, got, want) in checks {
                if got.to_bits() != want.to_bits() && !(got.is_nan() && want.is_nan()) {
                    bad.push(format!("{name} x={x} y={y} p={p} got={got:e} want={want:e}"));
                }
            }
        }
        assert!(bad.is_empty(), "{} differences, first: {:?}", bad.len(), &bad[..bad.len().min(5)]);
    }
}
