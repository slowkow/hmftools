//! java.util.TimSort, as used by Collections.sort / List.sort / Stream.sorted on objects.
//!
//! A stable sort gives one result for a consistent comparator, but some Java comparators are not consistent
//! (for example two different contigs of the same rank each compare as greater than the other), and then the
//! result depends on the exact sequence of comparisons. This is a direct port so the outcome is the same.

use std::cmp::Ordering;

const MIN_MERGE: usize = 32;
const MIN_GALLOP: i32 = 7;

/// Sorts like `list.sort(comparator)` in Java.
pub fn java_sort<T: Clone, F: FnMut(&T, &T) -> Ordering>(a: &mut [T], mut c: F) {
    let hi = a.len();
    let mut lo = 0usize;
    let mut n_remaining = hi;
    if n_remaining < 2 {
        return;
    }
    if n_remaining < MIN_MERGE {
        let init_run_len = count_run_and_make_ascending(a, lo, hi, &mut c);
        binary_sort(a, lo, hi, lo + init_run_len, &mut c);
        return;
    }
    let mut ts = TimSort { min_gallop: MIN_GALLOP, tmp: Vec::new(), run_base: Vec::new(), run_len: Vec::new() };
    let min_run = min_run_length(n_remaining);
    loop {
        let mut run_len = count_run_and_make_ascending(a, lo, hi, &mut c);
        if run_len < min_run {
            let force = if n_remaining <= min_run { n_remaining } else { min_run };
            binary_sort(a, lo, lo + force, lo + run_len, &mut c);
            run_len = force;
        }
        ts.run_base.push(lo);
        ts.run_len.push(run_len);
        ts.merge_collapse(a, &mut c);
        lo += run_len;
        n_remaining -= run_len;
        if n_remaining == 0 {
            break;
        }
    }
    ts.merge_force_collapse(a, &mut c);
}

/// Sorts by a key with natural ordering (Comparator.comparing / comparingInt).
pub fn java_sort_by_key<T: Clone, K: Ord, F: FnMut(&T) -> K>(a: &mut [T], mut key: F) {
    java_sort(a, |x, y| key(x).cmp(&key(y)));
}

fn lt<T, F: FnMut(&T, &T) -> Ordering>(c: &mut F, x: &T, y: &T) -> bool {
    c(x, y) == Ordering::Less
}

fn binary_sort<T: Clone, F: FnMut(&T, &T) -> Ordering>(a: &mut [T], lo: usize, hi: usize, mut start: usize, c: &mut F) {
    if start == lo {
        start += 1;
    }
    while start < hi {
        let pivot = a[start].clone();
        let mut left = lo;
        let mut right = start;
        while left < right {
            let mid = (left + right) >> 1;
            if lt(c, &pivot, &a[mid]) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        a[left..=start].rotate_right(1);
        start += 1;
    }
}

fn count_run_and_make_ascending<T, F: FnMut(&T, &T) -> Ordering>(a: &mut [T], lo: usize, hi: usize, c: &mut F) -> usize {
    let mut run_hi = lo + 1;
    if run_hi == hi {
        return 1;
    }
    let descending = lt(c, &a[run_hi], &a[lo]);
    run_hi += 1;
    if descending {
        while run_hi < hi && lt(c, &a[run_hi], &a[run_hi - 1]) {
            run_hi += 1;
        }
        a[lo..run_hi].reverse();
    } else {
        while run_hi < hi && !lt(c, &a[run_hi], &a[run_hi - 1]) {
            run_hi += 1;
        }
    }
    run_hi - lo
}

fn min_run_length(mut n: usize) -> usize {
    let mut r = 0;
    while n >= MIN_MERGE {
        r |= n & 1;
        n >>= 1;
    }
    n + r
}

struct TimSort<T> {
    min_gallop: i32,
    tmp: Vec<T>,
    run_base: Vec<usize>,
    run_len: Vec<usize>,
}

impl<T: Clone> TimSort<T> {
    fn merge_collapse<F: FnMut(&T, &T) -> Ordering>(&mut self, a: &mut [T], c: &mut F) {
        while self.run_len.len() > 1 {
            let mut n = self.run_len.len() as isize - 2;
            let rl = |i: isize| self.run_len[i as usize];
            if (n > 0 && rl(n - 1) <= rl(n) + rl(n + 1)) || (n > 1 && rl(n - 2) <= rl(n) + rl(n - 1)) {
                if rl(n - 1) < rl(n + 1) {
                    n -= 1;
                }
            } else if n < 0 || rl(n) > rl(n + 1) {
                break;
            }
            self.merge_at(a, n as usize, c);
        }
    }

    fn merge_force_collapse<F: FnMut(&T, &T) -> Ordering>(&mut self, a: &mut [T], c: &mut F) {
        while self.run_len.len() > 1 {
            let mut n = self.run_len.len() - 2;
            if n > 0 && self.run_len[n - 1] < self.run_len[n + 1] {
                n -= 1;
            }
            self.merge_at(a, n, c);
        }
    }

    fn merge_at<F: FnMut(&T, &T) -> Ordering>(&mut self, a: &mut [T], i: usize, c: &mut F) {
        let mut base1 = self.run_base[i];
        let mut len1 = self.run_len[i];
        let base2 = self.run_base[i + 1];
        let mut len2 = self.run_len[i + 1];
        self.run_len[i] = len1 + len2;
        self.run_base.remove(i + 1);
        self.run_len.remove(i + 1);
        let k = gallop_right(&a[base2].clone(), a, base1, len1, 0, c);
        base1 += k;
        len1 -= k;
        if len1 == 0 {
            return;
        }
        len2 = gallop_left(&a[base1 + len1 - 1].clone(), a, base2, len2, len2 - 1, c);
        if len2 == 0 {
            return;
        }
        if len1 <= len2 {
            self.merge_lo(a, base1, len1, base2, len2, c);
        } else {
            self.merge_hi(a, base1, len1, base2, len2, c);
        }
    }

    fn merge_lo<F: FnMut(&T, &T) -> Ordering>(&mut self, a: &mut [T], base1: usize, mut len1: usize, base2: usize, mut len2: usize, c: &mut F) {
        self.tmp.clear();
        self.tmp.extend_from_slice(&a[base1..base1 + len1]);
        let tmp = std::mem::take(&mut self.tmp);
        let mut cursor1 = 0usize;
        let mut cursor2 = base2;
        let mut dest = base1;
        a[dest] = a[cursor2].clone();
        dest += 1;
        cursor2 += 1;
        len2 -= 1;
        if len2 == 0 {
            a[dest..dest + len1].clone_from_slice(&tmp[cursor1..cursor1 + len1]);
            self.tmp = tmp;
            return;
        }
        if len1 == 1 {
            copy_within(a, cursor2, dest, len2);
            a[dest + len2] = tmp[cursor1].clone();
            self.tmp = tmp;
            return;
        }
        let mut min_gallop = self.min_gallop;
        'outer: loop {
            let mut count1: i32 = 0;
            let mut count2: i32 = 0;
            loop {
                if lt(c, &a[cursor2], &tmp[cursor1]) {
                    a[dest] = a[cursor2].clone();
                    dest += 1;
                    cursor2 += 1;
                    count2 += 1;
                    count1 = 0;
                    len2 -= 1;
                    if len2 == 0 {
                        break 'outer;
                    }
                } else {
                    a[dest] = tmp[cursor1].clone();
                    dest += 1;
                    cursor1 += 1;
                    count1 += 1;
                    count2 = 0;
                    len1 -= 1;
                    if len1 == 1 {
                        break 'outer;
                    }
                }
                if (count1 | count2) >= min_gallop {
                    break;
                }
            }
            loop {
                count1 = gallop_right(&a[cursor2].clone(), &tmp, cursor1, len1, 0, c) as i32;
                if count1 != 0 {
                    let n = count1 as usize;
                    a[dest..dest + n].clone_from_slice(&tmp[cursor1..cursor1 + n]);
                    dest += n;
                    cursor1 += n;
                    len1 -= n;
                    if len1 <= 1 {
                        break 'outer;
                    }
                }
                a[dest] = a[cursor2].clone();
                dest += 1;
                cursor2 += 1;
                len2 -= 1;
                if len2 == 0 {
                    break 'outer;
                }
                count2 = gallop_left(&tmp[cursor1], a, cursor2, len2, 0, c) as i32;
                if count2 != 0 {
                    let n = count2 as usize;
                    copy_within(a, cursor2, dest, n);
                    dest += n;
                    cursor2 += n;
                    len2 -= n;
                    if len2 == 0 {
                        break 'outer;
                    }
                }
                a[dest] = tmp[cursor1].clone();
                dest += 1;
                cursor1 += 1;
                len1 -= 1;
                if len1 == 1 {
                    break 'outer;
                }
                min_gallop -= 1;
                if !(count1 >= MIN_GALLOP || count2 >= MIN_GALLOP) {
                    break;
                }
            }
            if min_gallop < 0 {
                min_gallop = 0;
            }
            min_gallop += 2;
        }
        self.min_gallop = if min_gallop < 1 { 1 } else { min_gallop };
        if len1 == 1 {
            copy_within(a, cursor2, dest, len2);
            a[dest + len2] = tmp[cursor1].clone();
        } else if len1 == 0 {
            panic!("Comparison method violates its general contract!");
        } else {
            a[dest..dest + len1].clone_from_slice(&tmp[cursor1..cursor1 + len1]);
        }
        self.tmp = tmp;
    }

    fn merge_hi<F: FnMut(&T, &T) -> Ordering>(&mut self, a: &mut [T], base1: usize, mut len1: usize, base2: usize, mut len2: usize, c: &mut F) {
        self.tmp.clear();
        self.tmp.extend_from_slice(&a[base2..base2 + len2]);
        let tmp = std::mem::take(&mut self.tmp);
        // cursors are kept as isize since they can step to one before the start
        let mut cursor1 = (base1 + len1) as isize - 1;
        let mut cursor2 = len2 as isize - 1;
        let mut dest = (base2 + len2) as isize - 1;
        a[dest as usize] = a[cursor1 as usize].clone();
        dest -= 1;
        cursor1 -= 1;
        len1 -= 1;
        if len1 == 0 {
            let d = (dest - (len2 as isize - 1)) as usize;
            a[d..d + len2].clone_from_slice(&tmp[0..len2]);
            self.tmp = tmp;
            return;
        }
        if len2 == 1 {
            dest -= len1 as isize;
            cursor1 -= len1 as isize;
            copy_within(a, (cursor1 + 1) as usize, (dest + 1) as usize, len1);
            a[dest as usize] = tmp[cursor2 as usize].clone();
            self.tmp = tmp;
            return;
        }
        let mut min_gallop = self.min_gallop;
        'outer: loop {
            let mut count1: i32 = 0;
            let mut count2: i32 = 0;
            loop {
                if lt(c, &tmp[cursor2 as usize], &a[cursor1 as usize]) {
                    a[dest as usize] = a[cursor1 as usize].clone();
                    dest -= 1;
                    cursor1 -= 1;
                    count1 += 1;
                    count2 = 0;
                    len1 -= 1;
                    if len1 == 0 {
                        break 'outer;
                    }
                } else {
                    a[dest as usize] = tmp[cursor2 as usize].clone();
                    dest -= 1;
                    cursor2 -= 1;
                    count2 += 1;
                    count1 = 0;
                    len2 -= 1;
                    if len2 == 1 {
                        break 'outer;
                    }
                }
                if (count1 | count2) >= min_gallop {
                    break;
                }
            }
            loop {
                count1 = (len1 - gallop_right(&tmp[cursor2 as usize], a, base1, len1, len1 - 1, c)) as i32;
                if count1 != 0 {
                    let n = count1 as usize;
                    dest -= n as isize;
                    cursor1 -= n as isize;
                    len1 -= n;
                    copy_within(a, (cursor1 + 1) as usize, (dest + 1) as usize, n);
                    if len1 == 0 {
                        break 'outer;
                    }
                }
                a[dest as usize] = tmp[cursor2 as usize].clone();
                dest -= 1;
                cursor2 -= 1;
                len2 -= 1;
                if len2 == 1 {
                    break 'outer;
                }
                count2 = (len2 - gallop_left(&a[cursor1 as usize].clone(), &tmp, 0, len2, len2 - 1, c)) as i32;
                if count2 != 0 {
                    let n = count2 as usize;
                    dest -= n as isize;
                    cursor2 -= n as isize;
                    len2 -= n;
                    let d = (dest + 1) as usize;
                    let s = (cursor2 + 1) as usize;
                    a[d..d + n].clone_from_slice(&tmp[s..s + n]);
                    if len2 <= 1 {
                        break 'outer;
                    }
                }
                a[dest as usize] = a[cursor1 as usize].clone();
                dest -= 1;
                cursor1 -= 1;
                len1 -= 1;
                if len1 == 0 {
                    break 'outer;
                }
                min_gallop -= 1;
                if !(count1 >= MIN_GALLOP || count2 >= MIN_GALLOP) {
                    break;
                }
            }
            if min_gallop < 0 {
                min_gallop = 0;
            }
            min_gallop += 2;
        }
        self.min_gallop = if min_gallop < 1 { 1 } else { min_gallop };
        if len2 == 1 {
            dest -= len1 as isize;
            cursor1 -= len1 as isize;
            copy_within(a, (cursor1 + 1) as usize, (dest + 1) as usize, len1);
            a[dest as usize] = tmp[cursor2 as usize].clone();
        } else if len2 == 0 {
            panic!("Comparison method violates its general contract!");
        } else {
            let d = (dest - (len2 as isize - 1)) as usize;
            a[d..d + len2].clone_from_slice(&tmp[0..len2]);
        }
        self.tmp = tmp;
    }
}

/// System.arraycopy within one array (ranges may overlap)
fn copy_within<T: Clone>(a: &mut [T], src: usize, dest: usize, len: usize) {
    if src == dest || len == 0 {
        return;
    }
    if dest < src {
        for i in 0..len {
            a[dest + i] = a[src + i].clone();
        }
    } else {
        for i in (0..len).rev() {
            a[dest + i] = a[src + i].clone();
        }
    }
}

fn gallop_left<T, F: FnMut(&T, &T) -> Ordering>(key: &T, a: &[T], base: usize, len: usize, hint: usize, c: &mut F) -> usize {
    let (len, hint) = (len as isize, hint as isize);
    let at = |i: isize| &a[(base as isize + i) as usize];
    let mut last_ofs: isize = 0;
    let mut ofs: isize = 1;
    if c(key, at(hint)) == Ordering::Greater {
        let max_ofs = len - hint;
        while ofs < max_ofs && c(key, at(hint + ofs)) == Ordering::Greater {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        if ofs > max_ofs {
            ofs = max_ofs;
        }
        last_ofs += hint;
        ofs += hint;
    } else {
        let max_ofs = hint + 1;
        while ofs < max_ofs && c(key, at(hint - ofs)) != Ordering::Greater {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        if ofs > max_ofs {
            ofs = max_ofs;
        }
        let tmp = last_ofs;
        last_ofs = hint - ofs;
        ofs = hint - tmp;
    }
    last_ofs += 1;
    while last_ofs < ofs {
        let m = last_ofs + ((ofs - last_ofs) >> 1);
        if c(key, at(m)) == Ordering::Greater {
            last_ofs = m + 1;
        } else {
            ofs = m;
        }
    }
    ofs as usize
}

fn gallop_right<T, F: FnMut(&T, &T) -> Ordering>(key: &T, a: &[T], base: usize, len: usize, hint: usize, c: &mut F) -> usize {
    let (len, hint) = (len as isize, hint as isize);
    let at = |i: isize| &a[(base as isize + i) as usize];
    let mut ofs: isize = 1;
    let mut last_ofs: isize = 0;
    if c(key, at(hint)) == Ordering::Less {
        let max_ofs = hint + 1;
        while ofs < max_ofs && c(key, at(hint - ofs)) == Ordering::Less {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        if ofs > max_ofs {
            ofs = max_ofs;
        }
        let tmp = last_ofs;
        last_ofs = hint - ofs;
        ofs = hint - tmp;
    } else {
        let max_ofs = len - hint;
        while ofs < max_ofs && c(key, at(hint + ofs)) != Ordering::Less {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        if ofs > max_ofs {
            ofs = max_ofs;
        }
        last_ofs += hint;
        ofs += hint;
    }
    last_ofs += 1;
    while last_ofs < ofs {
        let m = last_ofs + ((ofs - last_ofs) >> 1);
        if c(key, at(m)) == Ordering::Less {
            ofs = m;
        } else {
            last_ofs = m + 1;
        }
    }
    ofs as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_stable_sort() {
        let mut seed = 12345u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        for n in [0usize, 1, 2, 5, 31, 32, 33, 64, 100, 257, 1000, 5000] {
            for modulus in [3u32, 50, 100000] {
                let v: Vec<(u32, usize)> = (0..n).map(|i| (next() % modulus, i)).collect();
                let mut a = v.clone();
                let mut b = v.clone();
                java_sort(&mut a, |x, y| x.0.cmp(&y.0));
                b.sort_by(|x, y| x.0.cmp(&y.0));
                assert_eq!(a, b, "n={n} mod={modulus}");
                // runs and partially sorted input
                let mut a2 = a.clone();
                a2.reverse();
                let mut b2 = a2.clone();
                java_sort(&mut a2, |x, y| x.0.cmp(&y.0));
                b2.sort_by(|x, y| x.0.cmp(&y.0));
                assert_eq!(a2, b2);
            }
        }
    }
}
