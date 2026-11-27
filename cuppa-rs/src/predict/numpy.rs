//! The arithmetic of numpy 1.26 and its OpenBLAS (0.3.23, the armv8 kernels of the numpy wheel) as
//! pycuppa's prediction meets it, operation for operation, so that the results equal numpy's to the last
//! bit (checked against numpy in the tests and on whole predictions):
//!
//! - `pairwise_sum`: numpy's float64 add reduction (np.sum, DataFrame.sum along a contiguous row);
//! - `dot`: OpenBLAS dgemv_t (kernel/arm64/gemv_t.S), which numpy uses for a row vector times a matrix
//!   whose columns are contiguous (`x @ coef_.T`, `np.dot(x, profiles)`);
//! - `interp`: np.interp (no fused multiply-add);
//! - `argsort`: np.argsort's default introsort (quicksort with insertion sort below 17 items and heapsort
//!   past the depth limit), whose order of ties reaches the outputs;
//! - the elementwise functions come from the C library (exp, log1p, pow, log10), as numpy's float64 loops
//!   on arm64 call them.

/// numpy's DOUBLE_pairwise_sum for a contiguous array: sequential below 8 items, 8 accumulators up to 128,
/// halves above
pub fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for v in a {
            res += *v;
        }
        res
    } else if n <= 128 {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

/// The sum of one row of a DataFrame of `n_rows` rows, as `DataFrame.sum(axis=1)` gives it: pandas sums
/// the transposed blocks, so for one row numpy reduces a contiguous row (the pairwise sum), and for several
/// rows it adds column after column (a sequential sum from the first value). The cohort order is taken in
/// Java mode only (CU7); corrected, each sample is summed as if alone.
pub fn row_sum(a: &[f64], n_rows: usize) -> f64 {
    if n_rows > 1 && crate::compat::java("CU7") {
        let Some((first, rest)) = a.split_first() else { return 0.0 };
        let mut s = *first;
        for v in rest {
            s += *v;
        }
        s
    } else {
        pairwise_sum(a)
    }
}

/// `DataFrame.sum(axis=1, min_count=1)` of one row: NaN taken as 0, NaN if every value is NaN
pub fn row_sum_skipna(a: &[f64], n_rows: usize) -> f64 {
    if a.iter().all(|v| v.is_nan()) {
        return f64::NAN;
    }
    if a.iter().any(|v| v.is_nan()) {
        let filled: Vec<f64> = a.iter().map(|v| if v.is_nan() { 0.0 } else { *v }).collect();
        return row_sum(&filled, n_rows);
    }
    row_sum(a, n_rows)
}

/// OpenBLAS 0.3.23 dgemm (driver/level3/level3.c with the armv8 dgemm_kernel_8x4) for one element of
/// `X @ W.T`: K split in blocks of GEMM_Q = 128 (the last two balanced, rounded up to GEMM_UNROLL_M = 8),
/// each block a fused chain over k starting with a product, added to C (zeroed by beta 0) as
/// C = alpha * acc + C with alpha 1
pub fn gemm_dot(x: &[f64], w: &[f64]) -> f64 {
    const Q: usize = 128;
    const UNROLL_M: usize = 8;
    let k = x.len();
    let mut c = 0.0f64;
    let mut ls = 0;
    while ls < k {
        let mut min_l = k - ls;
        if min_l >= Q * 2 {
            min_l = Q;
        } else if min_l > Q {
            min_l = (min_l / 2 + UNROLL_M - 1) / UNROLL_M * UNROLL_M;
        }
        let mut acc = x[ls] * w[ls];
        for i in ls + 1..ls + min_l {
            acc = x[i].mul_add(w[i], acc);
        }
        c = acc.mul_add(1.0, c);
        ls += min_l;
    }
    c
}

/// `X @ W.T` for several rows of X: as numpy computes it for a cohort (OpenBLAS dgemm) in Java mode (CU7);
/// corrected, each sample as if predicted alone (dgemv), so that its numbers do not depend on the others
pub fn gemm_rows(xs: &[&[f64]], w_rows: &[&[f64]]) -> Vec<Vec<f64>> {
    if crate::compat::java("CU7") {
        xs.iter().map(|x| w_rows.iter().map(|w| gemm_dot(x, w)).collect()).collect()
    } else {
        xs.iter().map(|x| w_rows.iter().map(|w| dot(w, x)).collect()).collect()
    }
}

/// OpenBLAS 0.3.23 dgemv_t for one output (kernel/arm64/gemv_t.S, unit strides): eight fused partial sums
/// over blocks of 32 (lanes of v1-v4) folded as ((v1 + v2) + v3) + v4, then pairs of 4 into v1, the two
/// lanes added, the last 1-3 products fused one by one; y = alpha * t + y with alpha 1 and y 0 (dscal by
/// beta 0 first).
pub fn dot(a: &[f64], x: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), x.len());
    let m = a.len();
    let mut v = [[0.0f64; 2]; 4];
    let mut i = 0;
    let blocks = m >> 5;
    if blocks > 0 {
        for _ in 0..blocks {
            for g in 0..4 {
                for (r, acc) in v.iter_mut().enumerate() {
                    for (l, lane) in acc.iter_mut().enumerate() {
                        let idx = i + g * 8 + r * 2 + l;
                        *lane = a[idx].mul_add(x[idx], *lane);
                    }
                }
            }
            i += 32;
        }
        for l in 0..2 {
            v[0][l] = ((v[0][l] + v[1][l]) + v[2][l]) + v[3][l];
        }
    }
    for _ in 0..((m & 31) >> 2) {
        for l in 0..2 {
            v[0][l] = a[i + l].mul_add(x[i + l], v[0][l]);
        }
        for l in 0..2 {
            v[0][l] = a[i + 2 + l].mul_add(x[i + 2 + l], v[0][l]);
        }
        i += 4;
    }
    let mut t = v[0][0] + v[0][1];
    for _ in 0..(m & 3) {
        t = a[i].mul_add(x[i], t);
        i += 1;
    }
    1.0f64.mul_add(t, 0.0)
}

/// OpenBLAS 0.3.23 dgemv_n for one output (kernel/arm64/gemv_n.S): y += (alpha * x[j]) * A[i, j] fused,
/// column after column, from y = 0
pub fn dot_n(a: &[f64], x: &[f64]) -> f64 {
    let mut y = 0.0f64;
    for (aj, xj) in a.iter().zip(x) {
        y = (1.0 * xj).mul_add(*aj, y);
    }
    y
}

/// np.interp(x, xp, fp) for one value, with the default left and right values fp[0] and fp[-1]
pub fn interp(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if n == 1 {
        // left, right and only value are all fp[0], NaN included
        return fp[0];
    }
    if x.is_nan() {
        return x;
    }
    // binary_search_with_guess: the last index j with xp[j] <= x, -1 below, n above
    if x > xp[n - 1] {
        return fp[n - 1];
    }
    if x < xp[0] {
        return fp[0];
    }
    let j = xp.partition_point(|v| *v <= x) - 1;
    if j == n - 1 {
        return fp[j];
    }
    if xp[j] == x {
        return fp[j];
    }
    let slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    let mut res = slope * (x - xp[j]) + fp[j];
    if res.is_nan() {
        res = slope * (x - xp[j + 1]) + fp[j + 1];
        if res.is_nan() && fp[j] == fp[j + 1] {
            res = fp[j];
        }
    }
    res
}

/// np.clip(x, lo, hi) for floats: maximum then minimum, NaN kept
pub fn clip(x: f64, lo: f64, hi: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    let y = if x < lo || lo.is_nan() { lo } else { x };
    if y > hi || hi.is_nan() { hi } else { y }
}

/// npysort's Tag::less for doubles: NaNs sort last
#[inline]
fn less(a: f64, b: f64) -> bool {
    a < b || (b != b && a == a)
}

/// np.argsort(v) with the default kind (npysort aquicksort for doubles)
pub fn argsort(v: &[f64]) -> Vec<usize> {
    let num = v.len();
    let mut tosort: Vec<usize> = (0..num).collect();
    if num < 2 {
        return tosort;
    }
    const SMALL_QUICKSORT: isize = 15;
    let t = &mut tosort;
    let mut pl: isize = 0;
    let mut pr: isize = num as isize - 1;
    let mut stack: Vec<(isize, isize, i32)> = Vec::new();
    // npy_get_msb(num) * 2
    let mut msb = 0;
    let mut u = num >> 1;
    while u > 0 {
        msb += 1;
        u >>= 1;
    }
    let mut cdepth: i32 = msb * 2;
    loop {
        let mut do_insertion = true;
        if cdepth < 0 {
            aheapsort(v, &mut t[pl as usize..=pr as usize]);
            do_insertion = false;
        } else {
            while (pr - pl) > SMALL_QUICKSORT {
                let pm = pl + ((pr - pl) >> 1);
                if less(v[t[pm as usize]], v[t[pl as usize]]) {
                    t.swap(pm as usize, pl as usize);
                }
                if less(v[t[pr as usize]], v[t[pm as usize]]) {
                    t.swap(pr as usize, pm as usize);
                }
                if less(v[t[pm as usize]], v[t[pl as usize]]) {
                    t.swap(pm as usize, pl as usize);
                }
                let vp = v[t[pm as usize]];
                let mut pi = pl;
                let mut pj = pr - 1;
                t.swap(pm as usize, pj as usize);
                loop {
                    loop {
                        pi += 1;
                        if !less(v[t[pi as usize]], vp) {
                            break;
                        }
                    }
                    loop {
                        pj -= 1;
                        if !less(vp, v[t[pj as usize]]) {
                            break;
                        }
                    }
                    if pi >= pj {
                        break;
                    }
                    t.swap(pi as usize, pj as usize);
                }
                let pk = pr - 1;
                t.swap(pi as usize, pk as usize);
                cdepth -= 1;
                // push the larger partition on the stack
                if pi - pl < pr - pi {
                    stack.push((pi + 1, pr, cdepth));
                    pr = pi - 1;
                } else {
                    stack.push((pl, pi - 1, cdepth));
                    pl = pi + 1;
                }
            }
        }
        if do_insertion {
            let mut pi = pl + 1;
            while pi <= pr {
                let vi = t[pi as usize];
                let vp = v[vi];
                let mut pj = pi;
                let mut pk = pi - 1;
                while pj > pl && less(vp, v[t[pk as usize]]) {
                    t[pj as usize] = t[pk as usize];
                    pj -= 1;
                    pk -= 1;
                }
                t[pj as usize] = vi;
                pi += 1;
            }
        }
        match stack.pop() {
            None => break,
            Some((l, r, d)) => {
                pl = l;
                pr = r;
                cdepth = d;
            }
        }
    }
    tosort
}

/// npysort aheapsort for doubles on a slice of indices
fn aheapsort(v: &[f64], tosort: &mut [usize]) {
    let mut n = tosort.len();
    // 1-based indexing as in the C code
    let a = |k: usize| k - 1;
    let mut l = n >> 1;
    while l > 0 {
        let tmp = tosort[a(l)];
        let mut i = l;
        let mut j = l << 1;
        while j <= n {
            if j < n && less(v[tosort[a(j)]], v[tosort[a(j + 1)]]) {
                j += 1;
            }
            if less(v[tmp], v[tosort[a(j)]]) {
                tosort[a(i)] = tosort[a(j)];
                i = j;
                j += j;
            } else {
                break;
            }
        }
        tosort[a(i)] = tmp;
        l -= 1;
    }
    while n > 1 {
        let tmp = tosort[a(n)];
        tosort[a(n)] = tosort[a(1)];
        n -= 1;
        let mut i = 1;
        let mut j = 2;
        while j <= n {
            if j < n && less(v[tosort[a(j)]], v[tosort[a(j + 1)]]) {
                j += 1;
            }
            if less(v[tmp], v[tosort[a(j)]]) {
                tosort[a(i)] = tosort[a(j)];
                i = j;
                j += j;
            } else {
                break;
            }
        }
        tosort[a(i)] = tmp;
    }
}

/// np.round(x, decimals) for float64: rint(x * 10^d) / 10^d (decimals >= 0)
pub fn round(x: f64, decimals: i32) -> f64 {
    if decimals == 0 {
        return x.round_ties_even();
    }
    let f = 10f64.powi(decimals);
    (x * f).round_ties_even() / f
}

/// softmax of one row as sklearn.utils.extmath.softmax: subtract the maximum, exp, divide by the pairwise
/// sum
pub fn softmax(row: &mut [f64]) {
    let mut max = f64::NEG_INFINITY;
    let mut nan = false;
    for v in row.iter() {
        if v.is_nan() {
            nan = true;
        } else if *v > max {
            max = *v;
        }
    }
    // np.max propagates NaN
    if nan {
        max = f64::NAN;
    }
    for v in row.iter_mut() {
        *v = (*v - max).exp();
    }
    let s = pairwise_sum(row);
    for v in row.iter_mut() {
        *v /= s;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_sums() {
        assert_eq!(pairwise_sum(&[]), 0.0);
        assert_eq!(pairwise_sum(&[1.0, 2.0, 3.0]), 6.0);
        let v: Vec<f64> = (0..300).map(|i| i as f64 * 0.1).collect();
        // halves of 144 and 156, each with 8 accumulators
        let expect = pairwise_sum(&v[..144]) + pairwise_sum(&v[144..]);
        assert_eq!(pairwise_sum(&v), expect);
    }

    #[test]
    fn argsort_ties() {
        // numpy: np.argsort([0.0, 1.0, 0.0, 2.0, 0.0]) == [0, 2, 4, 1, 3] (insertion sort, stable for small n)
        assert_eq!(argsort(&[0.0, 1.0, 0.0, 2.0, 0.0]), vec![0, 2, 4, 1, 3]);
        assert_eq!(argsort(&[f64::NAN, 1.0, 0.5]), vec![2, 1, 0]);
    }

    #[test]
    fn interp_edges() {
        let xp = [0.0, 1.0, 2.0];
        let fp = [0.0, 10.0, 30.0];
        assert_eq!(interp(-1.0, &xp, &fp), 0.0);
        assert_eq!(interp(3.0, &xp, &fp), 30.0);
        assert_eq!(interp(2.0, &xp, &fp), 30.0);
        assert_eq!(interp(1.5, &xp, &fp), 20.0);
        assert!(interp(f64::NAN, &xp, &fp).is_nan());
    }
}
