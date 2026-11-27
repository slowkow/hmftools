//! Number formatting that reproduces java.util.Formatter and htsjdk's VCF encoder.

/// Shortest round-trip decimal digits of a positive finite double: (digits, exponent) where the
/// value is 0.d1d2d3... x 10^exponent.
fn digits(d: f64) -> (Vec<u8>, i32) {
    let s = format!("{:e}", d);
    let (mant, exp) = s.split_once('e').unwrap();
    let digits: Vec<u8> = mant.bytes().filter(|b| b.is_ascii_digit()).map(|b| b - b'0').collect();
    (digits, exp.parse::<i32>().unwrap() + 1)
}

/// Round the digit string half-up so that `keep` digits remain; may carry into a new leading digit.
fn round_half_up(digits: &mut Vec<u8>, exp: &mut i32, keep: i32) {
    if keep < 0 {
        digits.clear();
        return;
    }
    let keep = keep as usize;
    if digits.len() <= keep {
        return;
    }
    let round_up = digits[keep] >= 5;
    digits.truncate(keep);
    if round_up {
        let mut i = keep;
        loop {
            if i == 0 {
                digits.insert(0, 1);
                *exp += 1;
                break;
            }
            i -= 1;
            if digits[i] == 9 {
                digits[i] = 0;
            } else {
                digits[i] += 1;
                break;
            }
        }
    }
}

/// String.format("%.<prec>f", d)
pub fn java_fixed(d: f64, prec: usize) -> String {
    if d.is_nan() {
        return "NaN".into();
    }
    if d.is_infinite() {
        return if d > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    let neg = d < 0.0 || (d == 0.0 && d.is_sign_negative());
    let a = d.abs();
    let mut int_part: Vec<u8> = vec![0];
    let mut frac: Vec<u8> = vec![0; prec];
    if a != 0.0 {
        let (mut ds, mut exp) = digits(a);
        let keep = exp + prec as i32;
        round_half_up(&mut ds, &mut exp, keep);
        // value = 0.ds x 10^exp
        let mut all: Vec<u8> = Vec::new();
        if exp <= 0 {
            all.push(0);
            all.extend(std::iter::repeat(0).take((-exp) as usize));
            all.extend(&ds);
            int_part = vec![0];
            frac = all[1..].to_vec();
        } else {
            let e = exp as usize;
            let mut ds2 = ds.clone();
            while ds2.len() < e {
                ds2.push(0);
            }
            int_part = ds2[..e].to_vec();
            frac = ds2[e..].to_vec();
        }
        frac.resize(prec, 0);
    }
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    s.extend(int_part.iter().map(|d| (d + b'0') as char));
    if prec > 0 {
        s.push('.');
        s.extend(frac.iter().map(|d| (d + b'0') as char));
    }
    s
}

/// String.format("%.<prec>e", d)
pub fn java_sci(d: f64, prec: usize) -> String {
    if d.is_nan() {
        return "NaN".into();
    }
    if d.is_infinite() {
        return if d > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    let neg = d < 0.0;
    let a = d.abs();
    let (mut ds, mut exp) = if a == 0.0 { (vec![0], 1) } else { digits(a) };
    round_half_up(&mut ds, &mut exp, prec as i32 + 1);
    ds.resize(prec + 1, 0);
    let e = exp - 1;
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    s.push((ds[0] + b'0') as char);
    if prec > 0 {
        s.push('.');
        s.extend(ds[1..].iter().map(|d| (d + b'0') as char));
    }
    s.push('e');
    s.push(if e < 0 { '-' } else { '+' });
    s.push_str(&format!("{:02}", e.abs()));
    s
}

/// htsjdk VCFEncoder.formatVCFDouble
pub fn vcf_double(d: f64) -> String {
    if d < 1.0 {
        if d < 0.01 {
            if d.abs() >= 1e-20 {
                java_sci(d, 3)
            } else {
                "0.00".into()
            }
        } else {
            java_fixed(d, 3)
        }
    } else {
        java_fixed(d, 2)
    }
}

/// java.text.DecimalFormat("#.###...") with English symbols and `max_fraction` digits after the point:
/// HALF_EVEN rounding, no trailing zeros, no grouping, at least one integer digit, "-0" for negative values
/// that round to zero, "NaN" and "∞".
///
/// Java does not round the exact binary value: it rounds the decimal digits of the old
/// `FloatingDecimal` conversion (as `Double.toString` printed them before JDK 19, an E-form with one
/// digit given a second, "5.0E-5"), helped by that conversion's flags for whether its last digit was
/// rounded up or is exact. So 0.00005 (just above the tie) gives "0" where exact rounding gives
/// "0.0001". This reproduces the conversion (`java_double_digits`) and DigitList's rounding.
pub fn decimal_format(v: f64, max_fraction: usize) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    let negative = v < 0.0 || (v == 0.0 && v.is_sign_negative());
    if v.is_infinite() {
        return if negative { "-∞".into() } else { "∞".into() };
    }
    let a = v.abs();
    // DigitList.set(isNegative, double, maximumDigits = max_fraction, fixedPoint = true)
    let (mut digits, mut decimal_at) = if a == 0.0 { (Vec::new(), 0) } else { digit_list(a, max_fraction as i32) };
    if digits.is_empty() {
        decimal_at = 0; // subformat: digitList.isZero() normalises
    }
    let mut out = String::new();
    if negative {
        out.push('-');
    }
    // integer part: at least one digit
    let count = decimal_at.max(1);
    let mut index = 0usize;
    for i in (0..count).rev() {
        if i < decimal_at && index < digits.len() {
            out.push(digits[index] as char);
            index += 1;
        } else {
            out.push('0');
        }
    }
    if index < digits.len() {
        out.push('.');
        for i in 0..max_fraction as i32 {
            if index >= digits.len() {
                break;
            }
            if -1 - i > decimal_at - 1 {
                out.push('0');
                continue;
            }
            out.push(digits[index] as char);
            index += 1;
        }
    }
    digits.clear();
    out
}

/// DigitList.set(boolean, double, int, true) for a positive finite value: the kept digits (ASCII) and
/// decimalAt, the number of digits before the point.
fn digit_list(a: f64, maximum_digits: i32) -> (Vec<u8>, i32) {
    let conv = java_double_digits(a);
    // parse the Java format string as DigitList does
    let s = conv.java_format_string();
    let mut digits: Vec<u8> = Vec::new();
    let mut decimal_at: i32 = -1;
    let mut exponent = 0i32;
    let mut leading_zeros_after_decimal = 0i32;
    let mut non_zero_seen = false;
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        i += 1;
        if c == b'.' {
            decimal_at = digits.len() as i32;
        } else if c == b'e' || c == b'E' {
            exponent = s[i..].parse::<i32>().unwrap();
            break;
        } else {
            if !non_zero_seen {
                non_zero_seen = c != b'0';
                if !non_zero_seen && decimal_at != -1 {
                    leading_zeros_after_decimal += 1;
                }
            }
            if non_zero_seen {
                digits.push(c);
            }
        }
    }
    if decimal_at == -1 {
        decimal_at = digits.len() as i32;
    }
    if non_zero_seen {
        decimal_at += exponent - leading_zeros_after_decimal;
    }
    if -decimal_at > maximum_digits {
        return (Vec::new(), decimal_at);
    } else if -decimal_at == maximum_digits {
        if should_round_up(&digits, 0, conv.rounded_up, conv.exact) {
            return (vec![b'1'], decimal_at + 1);
        }
        return (Vec::new(), decimal_at);
    }
    while digits.len() > 1 && digits[digits.len() - 1] == b'0' {
        digits.pop();
    }
    // round(maximumDigits + decimalAt)
    let max = maximum_digits + decimal_at;
    if max >= 0 && (max as usize) < digits.len() {
        let mut max = max as usize;
        if should_round_up(&digits, max, conv.rounded_up, conv.exact) {
            // roundUp
            loop {
                if max == 0 {
                    digits[0] = b'1';
                    decimal_at += 1;
                    max = 0;
                    break;
                }
                max -= 1;
                digits[max] += 1;
                if digits[max] <= b'9' {
                    break;
                }
            }
            max += 1;
        }
        digits.truncate(max);
        while digits.len() > 1 && digits[digits.len() - 1] == b'0' {
            digits.pop();
        }
    }
    (digits, decimal_at)
}

/// DigitList.shouldRoundUp for HALF_EVEN
fn should_round_up(digits: &[u8], maximum_digits: usize, already_rounded: bool, exact: bool) -> bool {
    if maximum_digits >= digits.len() {
        return false;
    }
    let d = digits[maximum_digits];
    if d > b'5' {
        return true;
    }
    if d == b'5' {
        if maximum_digits == digits.len() - 1 {
            if exact {
                return maximum_digits > 0 && !(digits[maximum_digits - 1] - b'0').is_multiple_of(2);
            }
            return !already_rounded;
        }
        return digits[maximum_digits + 1..].iter().any(|&c| c != b'0');
    }
    false
}

/// The result of jdk.internal.math.FloatingDecimal's binary to decimal conversion of a positive finite
/// double: digits (ASCII), decExponent (value = 0.d1d2... x 10^decExponent), and whether the last digit was
/// rounded up and whether the digits are exact.
#[derive(Debug, Clone)]
pub struct JavaDoubleDigits {
    pub digits: Vec<u8>,
    pub dec_exponent: i32,
    pub rounded_up: bool,
    pub exact: bool,
}

impl JavaDoubleDigits {
    /// BinaryToASCIIBuffer.getChars without the sign: the string the pre-JDK 19 Double.toString printed.
    pub fn java_format_string(&self) -> String {
        let d = &self.digits;
        let n = d.len() as i32;
        let e = self.dec_exponent;
        let mut s = String::new();
        let push = |s: &mut String, b: &[u8]| s.push_str(std::str::from_utf8(b).unwrap());
        if e > 0 && e < 8 {
            let char_length = n.min(e) as usize;
            push(&mut s, &d[..char_length]);
            if (char_length as i32) < e {
                for _ in 0..(e - char_length as i32) {
                    s.push('0');
                }
                s.push_str(".0");
            } else {
                s.push('.');
                if (char_length as i32) < n {
                    push(&mut s, &d[char_length..]);
                } else {
                    s.push('0');
                }
            }
        } else if e <= 0 && e > -3 {
            s.push_str("0.");
            for _ in 0..(-e) {
                s.push('0');
            }
            push(&mut s, d);
        } else {
            s.push(d[0] as char);
            s.push('.');
            if n > 1 {
                push(&mut s, &d[1..]);
            } else {
                s.push('0');
            }
            s.push('E');
            let x = if e <= 0 {
                s.push('-');
                -e + 1
            } else {
                e - 1
            };
            s.push_str(&x.to_string());
        }
        s
    }
}

const N_5_BITS: [i32; 27] = [0, 3, 5, 7, 10, 12, 14, 17, 19, 21, 24, 26, 28, 31, 33, 35, 38, 40, 42, 45, 47, 49, 52, 54, 56, 59, 61];
const INSIGNIFICANT_DIGITS_NUMBER: [i32; 64] = [
    0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12,
    12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19,
];

fn pow5_i64(n: i32) -> i64 {
    5i64.pow(n as u32)
}

/// FloatingDecimal.getBinaryToASCIIConverter(d) and BinaryToASCIIBuffer.dtoa for a positive finite d.
pub fn java_double_digits(d: f64) -> JavaDoubleDigits {
    const EXP_SHIFT: i32 = 52;
    const FRACT_HOB: i64 = 1 << 52;
    let bits = d.to_bits() as i64;
    let mut fract_bits = bits & 0x000f_ffff_ffff_ffff;
    let mut bin_exp = ((bits >> 52) & 0x7ff) as i32;
    let n_significant_bits;
    if bin_exp == 0 {
        let leading_zeros = fract_bits.leading_zeros() as i32;
        let shift = leading_zeros - (63 - EXP_SHIFT);
        fract_bits <<= shift;
        bin_exp = 1 - shift;
        n_significant_bits = 64 - leading_zeros;
    } else {
        fract_bits |= FRACT_HOB;
        n_significant_bits = EXP_SHIFT + 1;
    }
    bin_exp -= 1023;

    let tail_zeros = fract_bits.trailing_zeros() as i32;
    let n_fract_bits = EXP_SHIFT + 1 - tail_zeros;
    let n_tiny_bits = 0.max(n_fract_bits - bin_exp - 1);
    if (-21..=62).contains(&bin_exp) && n_tiny_bits < 27 && n_fract_bits + N_5_BITS[n_tiny_bits as usize] < 64 && n_tiny_bits == 0 {
        let insignificant = if bin_exp > n_significant_bits {
            let p2 = bin_exp - n_significant_bits - 1;
            if p2 > 1 && p2 < INSIGNIFICANT_DIGITS_NUMBER.len() as i32 { INSIGNIFICANT_DIGITS_NUMBER[p2 as usize] } else { 0 }
        } else {
            0
        };
        let v = if bin_exp >= EXP_SHIFT { fract_bits << (bin_exp - EXP_SHIFT) } else { ((fract_bits as u64) >> (EXP_SHIFT - bin_exp)) as i64 };
        return develop_long_digits(0, v, insignificant);
    }
    let mut dec_exp = estimate_dec_exp(fract_bits, bin_exp);
    let b5 = 0.max(-dec_exp);
    let mut b2 = b5 + n_tiny_bits + bin_exp;
    let s5 = 0.max(dec_exp);
    let mut s2 = s5 + n_tiny_bits;
    let m5 = b5;
    let mut m2 = b2 - n_significant_bits;
    fract_bits = ((fract_bits as u64) >> tail_zeros) as i64;
    b2 -= n_fract_bits - 1;
    let common = b2.min(s2);
    b2 -= common;
    s2 -= common;
    m2 -= common;
    if n_fract_bits == 1 {
        m2 -= 1;
    }
    if m2 < 0 {
        b2 -= m2;
        s2 -= m2;
        m2 = 0;
    }
    let b_bits = n_fract_bits + b2 + if (b5 as usize) < N_5_BITS.len() { N_5_BITS[b5 as usize] } else { b5 * 3 };
    let ten_s_bits = s2 + 1 + if ((s5 + 1) as usize) < N_5_BITS.len() { N_5_BITS[(s5 + 1) as usize] } else { (s5 + 1) * 3 };
    let mut digits: Vec<u8> = Vec::with_capacity(20);
    let low;
    let high;
    let low_digit_difference: i64;
    let exact;
    if b_bits < 64 && ten_s_bits < 64 {
        // the int and the long paths do the same in their own width (with the same overflow hack)
        if b_bits < 32 && ten_s_bits < 32 {
            let mut b: i32 = (fract_bits as i32).wrapping_mul(pow5_i64(b5) as i32) << b2;
            let s: i32 = (pow5_i64(s5) as i32) << s2;
            let mut m: i32 = (pow5_i64(m5) as i32) << m2;
            let tens: i32 = s.wrapping_mul(10);
            let mut q = b / s;
            b = 10i32.wrapping_mul(b % s);
            m = m.wrapping_mul(10);
            let mut lo = b < m;
            let mut hi = b.wrapping_add(m) > tens;
            if q == 0 && !hi {
                dec_exp -= 1;
            } else {
                digits.push(b'0' + q as u8);
            }
            if !(-3..8).contains(&dec_exp) {
                hi = false;
                lo = false;
            }
            while !lo && !hi {
                q = b / s;
                b = 10i32.wrapping_mul(b % s);
                m = m.wrapping_mul(10);
                if m > 0 {
                    lo = b < m;
                    hi = b.wrapping_add(m) > tens;
                } else {
                    lo = true;
                    hi = true;
                }
                digits.push(b'0' + q as u8);
            }
            low_digit_difference = ((b << 1).wrapping_sub(tens)) as i64;
            exact = b == 0;
            low = lo;
            high = hi;
        } else {
            let mut b: i64 = fract_bits.wrapping_mul(pow5_i64(b5)) << b2;
            let s: i64 = pow5_i64(s5) << s2;
            let mut m: i64 = pow5_i64(m5) << m2;
            let tens: i64 = s.wrapping_mul(10);
            let mut q = (b / s) as i32;
            b = 10i64.wrapping_mul(b % s);
            m = m.wrapping_mul(10);
            let mut lo = b < m;
            let mut hi = b.wrapping_add(m) > tens;
            if q == 0 && !hi {
                dec_exp -= 1;
            } else {
                digits.push(b'0' + q as u8);
            }
            if !(-3..8).contains(&dec_exp) {
                hi = false;
                lo = false;
            }
            while !lo && !hi {
                q = (b / s) as i32;
                b = 10i64.wrapping_mul(b % s);
                m = m.wrapping_mul(10);
                if m > 0 {
                    lo = b < m;
                    hi = b.wrapping_add(m) > tens;
                } else {
                    lo = true;
                    hi = true;
                }
                digits.push(b'0' + q as u8);
            }
            low_digit_difference = (b << 1).wrapping_sub(tens);
            exact = b == 0;
            low = lo;
            high = hi;
        }
    } else {
        let s_val = Big::pow52(s5, s2);
        let mut b_val = Big::from_u64(fract_bits as u64).mul_pow5(b5).shl(b2 as usize);
        let mut m_val = Big::pow52(m5 + 1, m2 + 1);
        let ten_s_val = Big::pow52(s5 + 1, s2 + 1);
        let mut q = b_val.quo_rem_iteration(&s_val);
        let mut lo = b_val.cmp(&m_val).is_lt();
        let mut hi = ten_s_val.cmp(&b_val.add(&m_val)).is_le();
        if q == 0 && !hi {
            dec_exp -= 1;
        } else {
            digits.push(b'0' + q as u8);
        }
        if !(-3..8).contains(&dec_exp) {
            hi = false;
            lo = false;
        }
        while !lo && !hi {
            q = b_val.quo_rem_iteration(&s_val);
            m_val = m_val.mul_small(10);
            lo = b_val.cmp(&m_val).is_lt();
            hi = ten_s_val.cmp(&b_val.add(&m_val)).is_le();
            digits.push(b'0' + q as u8);
        }
        if hi && lo {
            b_val = b_val.shl(1);
            low_digit_difference = match b_val.cmp(&ten_s_val) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            };
        } else {
            low_digit_difference = 0;
        }
        exact = b_val.is_zero();
        low = lo;
        high = hi;
    }
    let mut r = JavaDoubleDigits { digits, dec_exponent: dec_exp + 1, rounded_up: false, exact };
    if high {
        if low {
            if low_digit_difference == 0 {
                if (r.digits[r.digits.len() - 1] & 1) != 0 {
                    round_up(&mut r);
                }
            } else if low_digit_difference > 0 {
                round_up(&mut r);
            }
        } else {
            round_up(&mut r);
        }
    }
    r
}

/// BinaryToASCIIBuffer.roundup
fn round_up(r: &mut JavaDoubleDigits) {
    let mut i = r.digits.len() - 1;
    let mut q = r.digits[i];
    if q == b'9' {
        while q == b'9' && i > 0 {
            r.digits[i] = b'0';
            i -= 1;
            q = r.digits[i];
        }
        if q == b'9' {
            r.dec_exponent += 1;
            r.digits[0] = b'1';
            return;
        }
    }
    r.digits[i] = q + 1;
    r.rounded_up = true;
}

/// BinaryToASCIIBuffer.developLongDigits
fn develop_long_digits(mut dec_exponent: i32, mut lvalue: i64, insignificant: i32) -> JavaDoubleDigits {
    if insignificant != 0 {
        let pow10 = pow5_i64(insignificant) << insignificant;
        let residue = lvalue % pow10;
        lvalue /= pow10;
        dec_exponent += insignificant;
        if residue >= (pow10 >> 1) {
            lvalue += 1;
        }
    }
    // digits from the right; trailing zeros dropped
    let mut rev: Vec<u8> = Vec::new();
    let mut c = lvalue % 10;
    lvalue /= 10;
    while c == 0 {
        dec_exponent += 1;
        c = lvalue % 10;
        lvalue /= 10;
    }
    while lvalue != 0 {
        rev.push(b'0' + c as u8);
        dec_exponent += 1;
        c = lvalue % 10;
        lvalue /= 10;
    }
    rev.push(b'0' + c as u8);
    rev.reverse();
    JavaDoubleDigits { digits: rev, dec_exponent: dec_exponent + 1, rounded_up: false, exact: false }
}

/// BinaryToASCIIBuffer.estimateDecExp (its constants are Java's literals)
#[allow(clippy::approx_constant)]
fn estimate_dec_exp(fract_bits: i64, bin_exp: i32) -> i32 {
    let d2 = f64::from_bits((0x3ff0_0000_0000_0000u64) | (fract_bits as u64 & 0x000f_ffff_ffff_ffff));
    let d = (d2 - 1.5) * 0.289529654 + 0.176091259 + bin_exp as f64 * 0.301029995663981;
    let d_bits = d.to_bits() as i64;
    let exponent = (((d_bits as u64) & 0x7ff0_0000_0000_0000) >> 52) as i32 - 1023;
    let is_negative = (d_bits as u64 & 0x8000_0000_0000_0000) != 0;
    if (0..52).contains(&exponent) {
        let mask = 0x000f_ffff_ffff_ffffi64 >> exponent;
        let r = (((d_bits & 0x000f_ffff_ffff_ffff) | (1i64 << 52)) >> (52 - exponent)) as i32;
        if is_negative {
            if (mask & d_bits) == 0 { -r } else { -r - 1 }
        } else {
            r
        }
    } else if exponent < 0 {
        if (d_bits & !(0x8000_0000_0000_0000u64 as i64)) == 0 {
            0
        } else if is_negative {
            -1
        } else {
            0
        }
    } else {
        d as i32
    }
}

/// An unsigned big integer, little-endian 32-bit limbs: enough of FDBigInteger for dtoa's hard case.
#[derive(Clone, Debug)]
struct Big(Vec<u32>);

impl Big {
    fn from_u64(v: u64) -> Big {
        let mut b = Big(vec![v as u32, (v >> 32) as u32]);
        b.trim();
        b
    }
    fn trim(&mut self) {
        while self.0.len() > 1 && *self.0.last().unwrap() == 0 {
            self.0.pop();
        }
    }
    fn is_zero(&self) -> bool {
        self.0.iter().all(|&x| x == 0)
    }
    fn mul_small(&self, m: u32) -> Big {
        let mut carry = 0u64;
        let mut out = Vec::with_capacity(self.0.len() + 1);
        for &x in &self.0 {
            let p = x as u64 * m as u64 + carry;
            out.push(p as u32);
            carry = p >> 32;
        }
        if carry != 0 {
            out.push(carry as u32);
        }
        let mut b = Big(out);
        b.trim();
        b
    }
    fn mul_pow5(&self, mut n: i32) -> Big {
        let mut b = self.clone();
        while n >= 13 {
            b = b.mul_small(1_220_703_125); // 5^13
            n -= 13;
        }
        if n > 0 {
            b = b.mul_small(5u32.pow(n as u32));
        }
        b
    }
    fn pow52(p5: i32, p2: i32) -> Big {
        Big::from_u64(1).mul_pow5(p5).shl(p2 as usize)
    }
    fn shl(&self, n: usize) -> Big {
        let words = n / 32;
        let bits = n % 32;
        let mut out = vec![0u32; words];
        if bits == 0 {
            out.extend_from_slice(&self.0);
        } else {
            let mut carry = 0u32;
            for &x in &self.0 {
                out.push((x << bits) | carry);
                carry = x >> (32 - bits);
            }
            out.push(carry);
        }
        let mut b = Big(out);
        b.trim();
        b
    }
    fn add(&self, o: &Big) -> Big {
        let n = self.0.len().max(o.0.len());
        let mut out = Vec::with_capacity(n + 1);
        let mut carry = 0u64;
        for i in 0..n {
            let s = *self.0.get(i).unwrap_or(&0) as u64 + *o.0.get(i).unwrap_or(&0) as u64 + carry;
            out.push(s as u32);
            carry = s >> 32;
        }
        if carry != 0 {
            out.push(carry as u32);
        }
        let mut b = Big(out);
        b.trim();
        b
    }
    fn sub_assign(&mut self, o: &Big) {
        let mut borrow = 0i64;
        for i in 0..self.0.len() {
            let mut d = self.0[i] as i64 - *o.0.get(i).unwrap_or(&0) as i64 - borrow;
            if d < 0 {
                d += 1 << 32;
                borrow = 1;
            } else {
                borrow = 0;
            }
            self.0[i] = d as u32;
        }
        self.trim();
    }
    fn cmp(&self, o: &Big) -> std::cmp::Ordering {
        let (a, b) = (self.0.len(), o.0.len());
        let n = a.max(b);
        for i in (0..n).rev() {
            let x = *self.0.get(i).unwrap_or(&0);
            let y = *o.0.get(i).unwrap_or(&0);
            if x != y {
                return x.cmp(&y);
            }
        }
        std::cmp::Ordering::Equal
    }
    /// q = this / s (below 10), this = 10 * (this mod s)
    fn quo_rem_iteration(&mut self, s: &Big) -> u32 {
        let mut q = 0;
        while self.cmp(s).is_ge() {
            self.sub_assign(s);
            q += 1;
        }
        *self = self.mul_small(10);
        q
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn formats() {
        assert_eq!(java_fixed(0.125, 2), "0.13");
        assert_eq!(java_fixed(0.5555, 3), "0.556");
        assert_eq!(java_fixed(0.9996, 3), "1.000");
        assert_eq!(java_fixed(33.925, 2), "33.93");
        assert_eq!(java_fixed(0.0, 3), "0.000");
        assert_eq!(java_fixed(1234.5, 0), "1235");
        assert_eq!(java_sci(0.009668, 3), "9.668e-03");
        assert_eq!(java_sci(2.9394e-5, 3), "2.939e-05");
        assert_eq!(java_sci(-12.345, 3), "-1.235e+01");
        assert_eq!(vcf_double(0.0), "0.00");
        assert_eq!(vcf_double(0.192), "0.192");
        assert_eq!(vcf_double(33.93), "33.93");
    }

    /// Values and results printed by Java's DecimalFormat("#.####", English symbols)
    #[test]
    fn decimal_format_as_java() {
        assert_eq!(java_double_digits(5.0E-5).java_format_string(), "5.0E-5");
        assert_eq!(java_double_digits(100.0).java_format_string(), "100.0");
        assert_eq!(java_double_digits(0.001).java_format_string(), "0.001");
        assert_eq!(java_double_digits(1.0E7).java_format_string(), "1.0E7");
        assert_eq!(java_double_digits(1234567.5).java_format_string(), "1234567.5");
        assert_eq!(java_double_digits(0.1 + 0.2).java_format_string(), "0.30000000000000004");
        assert_eq!(java_double_digits(f64::MIN_POSITIVE / 8.0).java_format_string(), "2.781342323134002E-309");
        assert_eq!(java_double_digits(1e300).java_format_string(), "1.0E300");
        // the old FloatingDecimal's well-known extra digit (Double.toString since JDK 19 gives 2.0E-3)
        assert_eq!(java_double_digits(2.0E-3).java_format_string(), "0.002");
        let cases: [(f64, &str); 23] = [
            (0.5, "0.5"), (0.03125, "0.0312"), (0.09375, "0.0938"), (0.15625, "0.1562"), (-0.00001, "-0"),
            (-0.0, "-0"), (0.0, "0"), (1.0, "1"), (0.99995, "1"), (0.999951, "1"), (12345.678912, "12345.6789"),
            (1e20, "100000000000000000000"), (1e-20, "0"), (f64::NAN, "NaN"), (f64::INFINITY, "∞"),
            (f64::NEG_INFINITY, "-∞"), (0.00005, "0"), (0.00015, "0.0001"), (2.5, "2.5"), (0.12345, "0.1235"),
            (1234567.5, "1234567.5"), (3.0E-5, "0"), (123456789012.12345, "123456789012.1234"),
        ];
        for (v, expected) in cases {
            assert_eq!(decimal_format(v, 4), expected, "{v}");
        }
    }
}

/// Double.toString as Java 19 and later print it: the shortest digits that round-trip; plain notation
/// with at least one fraction digit for magnitudes from 1e-3 up to 1e7, else d.dddE<exponent>.
/// (Added 2026-11-26 for the driver gene panel writer; `java_double_digits` gives the older form.)
pub fn java_double_to_string(d: f64) -> String {
    if d.is_nan() {
        return "NaN".into();
    }
    if d.is_infinite() {
        return if d > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    if d == 0.0 {
        return if d.is_sign_negative() { "-0.0".into() } else { "0.0".into() };
    }
    let neg = d < 0.0;
    let (ds, exp) = digits(d.abs()); // value = 0.ds x 10^exp
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    let a = d.abs();
    if (1e-3..1e7).contains(&a) {
        if exp <= 0 {
            s.push_str("0.");
            s.extend(std::iter::repeat_n('0', (-exp) as usize));
            s.extend(ds.iter().map(|x| (x + b'0') as char));
        } else {
            let e = exp as usize;
            for i in 0..e {
                s.push(if i < ds.len() { (ds[i] + b'0') as char } else { '0' });
            }
            s.push('.');
            if ds.len() > e {
                s.extend(ds[e..].iter().map(|x| (x + b'0') as char));
            } else {
                s.push('0');
            }
        }
    } else {
        s.push((ds[0] + b'0') as char);
        s.push('.');
        if ds.len() > 1 {
            s.extend(ds[1..].iter().map(|x| (x + b'0') as char));
        } else {
            s.push('0');
        }
        s.push('E');
        s.push_str(&(exp - 1).to_string());
    }
    s
}

#[cfg(test)]
mod double_to_string_tests {
    use super::java_double_to_string as t;

    #[test]
    fn java_19_double_to_string() {
        assert_eq!(t(3.0), "3.0");
        assert_eq!(t(0.6), "0.6");
        assert_eq!(t(-2.5), "-2.5");
        assert_eq!(t(0.001), "0.001");
        assert_eq!(t(0.0001), "1.0E-4");
        assert_eq!(t(1e7), "1.0E7");
        assert_eq!(t(1234567.0), "1234567.0");
        assert_eq!(t(2e-3), "0.002");
        assert_eq!(t(1.0E23), "1.0E23");
        assert_eq!(t(123.456), "123.456");
        assert_eq!(t(-0.0), "-0.0");
    }
}
