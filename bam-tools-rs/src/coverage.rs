//! BaseCoverage (the depth of each base of a partition, and the bases each filter removed) and
//! CoverageMetrics (the combined depth distribution, its statistics and the filtered percentages).

use hmf_common::regions::BaseRegion;

use crate::compat;
use crate::record::{Rec, OP_D, OP_EQ, OP_I, OP_M, OP_N, OP_S, OP_X};

/// FilterType, in its order
pub const UNFILTERED: usize = 0;
pub const LOW_MAP_QUAL: usize = 1;
pub const DUPLICATE: usize = 2;
pub const MATE_UNMAPPED: usize = 3;
pub const LOW_BASE_QUAL: usize = 4;
pub const OVERLAPPED: usize = 5;
pub const MAX_COVERAGE: usize = 6;
pub const FILTER_TYPES: usize = 7;
pub const FILTER_NAMES: [&str; FILTER_TYPES] = ["UNFILTERED", "LOW_MAP_QUAL", "DUPLICATE", "MATE_UNMAPPED", "LOW_BASE_QUAL", "OVERLAPPED", "MAX_COVERAGE"];

/// The thresholds BaseCoverage applies
#[derive(Clone, Copy)]
pub struct CoverageParams {
    pub map_qual_threshold: i32,
    pub base_qual_threshold: i32,
    pub max_coverage: i32,
    pub expect_duplicates: bool,
}

/// BaseCoverage for one partition at a time (the buffer is reused)
pub struct BaseCoverage {
    params: CoverageParams,
    region_start: i32,
    region_end: i32,
    depth: Vec<u32>,
    pub filter_counts: [i64; FILTER_TYPES],
    /// the unmappable regions overlapping the partition, sorted by start
    unmappable: Vec<BaseRegion>,
    /// BT11: = and X taken as M
    eq_x_as_match: bool,
}

impl BaseCoverage {
    pub fn new(params: CoverageParams) -> Self {
        BaseCoverage {
            params,
            region_start: 1,
            region_end: 0,
            depth: Vec::new(),
            filter_counts: [0; FILTER_TYPES],
            unmappable: Vec::new(),
            eq_x_as_match: !compat::java("BT11"),
        }
    }

    /// starts a partition: depths and counts cleared
    pub fn reset(&mut self, region_start: i32, region_end: i32, unmappable: Vec<BaseRegion>) {
        self.region_start = region_start;
        self.region_end = region_end;
        let size = (region_end - region_start + 1).max(0) as usize;
        self.depth.clear();
        self.depth.resize(size, 0);
        self.filter_counts = [0; FILTER_TYPES];
        self.unmappable = unmappable;
    }

    #[cfg(test)]
    pub fn depth(&self) -> &[u32] {
        &self.depth
    }

    #[inline]
    fn is_match_op(&self, op: u32) -> bool {
        op == OP_M || (self.eq_x_as_match && (op == OP_EQ || op == OP_X))
    }

    /// The aligned bases the read-level filters count: the M operators (and = and X, BT11)
    #[inline]
    pub fn aligned_bases(&self, rec: &Rec) -> i64 {
        let mut n = 0i64;
        for &c in rec.cigar {
            if self.is_match_op(c & 0xf) {
                n += (c >> 4) as i64;
            }
        }
        n
    }

    /// BaseCoverage.processRead. `count_read` is false for a read that started in an earlier partition
    /// (BT1), whose read-level filter counts were made there; only its bases in this partition count.
    /// Returns an error where Java's array access would throw (more CIGAR bases than qualities).
    pub fn process_read(&mut self, rec: &Rec, mate_coords: Option<&[(i32, i32)]>, is_consensus: bool, count_read: bool) -> Result<(), String> {
        let p = self.params;
        if (rec.mapq as i32) < p.map_qual_threshold {
            if !is_consensus && count_read {
                self.filter_counts[LOW_MAP_QUAL] += self.aligned_bases(rec);
            }
            return Ok(());
        }
        if rec.duplicate() {
            if count_read {
                self.filter_counts[DUPLICATE] += self.aligned_bases(rec);
            }
            return Ok(());
        } else if is_consensus && p.expect_duplicates && count_read {
            // lower duplicate data to account for the additionally marked duplicate
            self.filter_counts[DUPLICATE] -= self.aligned_bases(rec);
        }
        if rec.paired() && rec.mate_unmapped() {
            if count_read {
                // BT13: Java counts all the read's bases here (soft clips included) where the other
                // read-level filters count aligned bases. Corrected: aligned bases.
                self.filter_counts[MATE_UNMAPPED] += if compat::java("BT13") { rec.seq_len as i64 } else { self.aligned_bases(rec) };
            }
            return Ok(());
        }
        let mut position = rec.start;
        let mut read_index = 0usize;
        for &c in rec.cigar {
            let op = c & 0xf;
            let len = (c >> 4) as i32;
            if self.is_match_op(op) {
                self.process_matched_bases(rec, position, read_index, len, mate_coords)?;
                position += len;
                read_index += len as usize;
            } else if op == OP_S || op == OP_I {
                read_index += len as usize;
            } else if op == OP_D || op == OP_N {
                position += len;
            }
            // H and P consume nothing; = and X nothing in Java (BT11)
        }
        Ok(())
    }

    fn process_matched_bases(&mut self, rec: &Rec, pos_start: i32, ri_start: usize, len: i32, mate_coords: Option<&[(i32, i32)]>) -> Result<(), String> {
        let lo = pos_start.max(self.region_start);
        let hi = (pos_start + len - 1).min(self.region_end);
        if lo > hi {
            return Ok(());
        }
        if ri_start + (hi - pos_start) as usize >= rec.quals.len() {
            return Err(format!(
                "ArrayIndexOutOfBoundsException: read {} has fewer base qualities than its CIGAR",
                String::from_utf8_lossy(rec.name)
            ));
        }
        let check_unmappable = !self.unmappable.is_empty() && self.unmappable.iter().any(|u| u.start <= rec.end && u.end >= pos_start);
        if !check_unmappable {
            self.count_range(rec, pos_start, ri_start, lo, hi, mate_coords);
            return Ok(());
        }
        // the parts of [lo, hi] outside the unmappable regions
        let mut from = lo;
        for i in 0..self.unmappable.len() {
            let u = self.unmappable[i];
            if u.end < from {
                continue;
            }
            if u.start > hi {
                break;
            }
            if u.start > from {
                self.count_range(rec, pos_start, ri_start, from, u.start - 1, mate_coords);
            }
            from = from.max(u.end + 1);
            if from > hi {
                return Ok(());
            }
        }
        if from <= hi {
            self.count_range(rec, pos_start, ri_start, from, hi, mate_coords);
        }
        Ok(())
    }

    /// counts the bases of [lo, hi] (within the partition, outside unmappable regions); bases inside the
    /// mate's (or the fragment's other reads') aligned coordinates are OVERLAPPED
    #[inline]
    fn count_range(&mut self, rec: &Rec, pos_start: i32, ri_start: usize, lo: i32, hi: i32, mate_coords: Option<&[(i32, i32)]>) {
        match mate_coords {
            Some(coords) if !coords.is_empty() => {
                // split [lo, hi] at the union of the coordinates
                let mut ivs: [(i32, i32); 16] = [(0, 0); 16];
                let mut n = 0;
                let mut heap: Vec<(i32, i32)> = Vec::new();
                let list: &mut [(i32, i32)] = if coords.len() <= 16 {
                    for &(a, b) in coords {
                        ivs[n] = (a, b);
                        n += 1;
                    }
                    &mut ivs[..n]
                } else {
                    heap.extend_from_slice(coords);
                    &mut heap[..]
                };
                list.sort_unstable();
                let mut from = lo;
                for &(a, b) in list.iter() {
                    if b < from || a > b {
                        continue;
                    }
                    if a > hi {
                        break;
                    }
                    if a > from {
                        self.count_plain(rec, pos_start, ri_start, from, a - 1);
                    }
                    let ov_end = b.min(hi);
                    let ov_start = a.max(from);
                    if ov_end >= ov_start {
                        self.filter_counts[OVERLAPPED] += (ov_end - ov_start + 1) as i64;
                    }
                    from = from.max(b + 1);
                    if from > hi {
                        return;
                    }
                }
                if from <= hi {
                    self.count_plain(rec, pos_start, ri_start, from, hi);
                }
            }
            _ => self.count_plain(rec, pos_start, ri_start, lo, hi),
        }
    }

    /// the inner loop: a base below the base quality threshold is LOW_BASE_QUAL, else one at the depth
    /// cap is MAX_COVERAGE, else it adds to the depth
    #[inline]
    fn count_plain(&mut self, rec: &Rec, pos_start: i32, ri_start: usize, lo: i32, hi: i32) {
        let n = (hi - lo + 1) as usize;
        let q0 = ri_start + (lo - pos_start) as usize;
        let quals = &rec.quals[q0..q0 + n];
        let d0 = (lo - self.region_start) as usize;
        let depth = &mut self.depth[d0..d0 + n];
        let thr = self.params.base_qual_threshold;
        let max = self.params.max_coverage.max(0) as u32;
        let (mut added, mut low) = (0u32, 0u32);
        if thr <= 0 {
            for d in depth.iter_mut() {
                let ok = (*d < max) as u32;
                *d += ok;
                added += ok;
            }
        } else if thr > 255 {
            low = n as u32;
        } else {
            let t = thr as u8;
            for (d, &q) in depth.iter_mut().zip(quals) {
                let is_low = (q < t) as u32;
                let ok = (1 - is_low) & ((*d < max) as u32);
                *d += ok;
                added += ok;
                low += is_low;
            }
        }
        self.filter_counts[UNFILTERED] += added as i64;
        self.filter_counts[LOW_BASE_QUAL] += low as i64;
        self.filter_counts[MAX_COVERAGE] += (n as u32 - added - low) as i64;
    }

    /// BaseCoverage.createMetrics: the depth distribution of the partition's bases outside unmappable
    /// regions (and inside the target regions, when the run has any), added to `metrics`
    pub fn add_metrics(&self, target_regions: Option<&[BaseRegion]>, exclude_zero_coverage: bool, metrics: &mut CoverageMetrics) {
        let (s, e) = (self.region_start, self.region_end);
        let counted = |lo: i32, hi: i32, m: &mut CoverageMetrics| {
            let depths = &self.depth[(lo - s) as usize..=(hi - s) as usize];
            let zero = histogram(depths, &mut m.coverage_frequency);
            let n = (hi - lo + 1) as i64;
            if exclude_zero_coverage {
                m.coverage_frequency[0] -= zero;
            }
            m.zero_coverage_bases += zero;
            m.coverage_bases += n - zero;
        };
        // the positions to count: [s, e], within the target regions (if any), outside unmappable regions
        let mut ranges: Vec<(i32, i32)> = Vec::new();
        match target_regions {
            Some(t) => {
                let mut ts: Vec<(i32, i32)> = t.iter().map(|r| (r.start.max(s), r.end.min(e))).filter(|r| r.0 <= r.1).collect();
                ts.sort_unstable();
                for (a, b) in ts {
                    match ranges.last_mut() {
                        Some(last) if a <= last.1 + 1 => last.1 = last.1.max(b),
                        _ => ranges.push((a, b)),
                    }
                }
            }
            None => {
                if s <= e {
                    ranges.push((s, e));
                }
            }
        }
        for (a, b) in ranges {
            let mut from = a;
            for u in &self.unmappable {
                if u.end < from {
                    continue;
                }
                if u.start > b {
                    break;
                }
                if u.start > from {
                    counted(from, u.start - 1, metrics);
                }
                from = from.max(u.end + 1);
                if from > b {
                    break;
                }
            }
            if from <= b {
                counted(from, b, metrics);
            }
        }
        for i in 0..FILTER_TYPES {
            metrics.filter_counts[i] += self.filter_counts[i];
        }
    }
}

/// Adds the depths to the frequencies (all depths are within them) and returns the number of zeros.
/// Blocks of zeros are counted at once; other blocks go through four sub-histograms, so that runs of
/// equal depths do not wait on the same counter.
fn histogram(depths: &[u32], freq: &mut [i64]) -> i64 {
    const BLOCK: usize = 64;
    let bins = freq.len();
    let mut lanes = vec![0u32; bins * 4];
    let (l0, rest) = lanes.split_at_mut(bins);
    let (l1, rest) = rest.split_at_mut(bins);
    let (l2, l3) = rest.split_at_mut(bins);
    let mut zero_blocks = 0i64;
    let mut blocks = depths.chunks_exact(BLOCK);
    for b in &mut blocks {
        if b.iter().fold(0u32, |a, &d| a | d) == 0 {
            zero_blocks += 1;
            continue;
        }
        for q in b.chunks_exact(4) {
            l0[q[0] as usize] += 1;
            l1[q[1] as usize] += 1;
            l2[q[2] as usize] += 1;
            l3[q[3] as usize] += 1;
        }
    }
    for &d in blocks.remainder() {
        l0[d as usize] += 1;
    }
    let zeros = zero_blocks * BLOCK as i64 + (l0[0] + l1[0] + l2[0] + l3[0]) as i64;
    freq[0] += zero_blocks * BLOCK as i64;
    for i in 0..bins {
        freq[i] += (l0[i] + l1[i] + l2[i] + l3[i]) as i64;
    }
    zeros
}

/// Statistics of the depth distribution
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Statistics {
    pub mean: f64,
    pub median: f64,
    pub mad: i32,
    pub std_deviation: f64,
}

/// CoverageMetrics
#[derive(Clone, Debug)]
pub struct CoverageMetrics {
    pub filter_counts: [i64; FILTER_TYPES],
    pub coverage_frequency: Vec<i64>,
    pub coverage_bases: i64,
    pub zero_coverage_bases: i64,
    total_bases: i64,
    statistics: Option<Statistics>,
}

impl CoverageMetrics {
    pub fn new(max_coverage: i32) -> Self {
        CoverageMetrics {
            filter_counts: [0; FILTER_TYPES],
            coverage_frequency: vec![0; max_coverage.max(0) as usize + 1],
            coverage_bases: 0,
            zero_coverage_bases: 0,
            total_bases: -1,
            statistics: None,
        }
    }

    /// genomeTerritory: bases with and without coverage (unmappable regions not included)
    pub fn genome_territory(&self) -> i64 {
        self.zero_coverage_bases + self.coverage_bases
    }

    pub fn merge(&mut self, other: &CoverageMetrics) {
        for i in 0..FILTER_TYPES {
            self.filter_counts[i] += other.filter_counts[i];
        }
        for (a, b) in self.coverage_frequency.iter_mut().zip(&other.coverage_frequency) {
            *a += *b;
        }
        self.coverage_bases += other.coverage_bases;
        self.zero_coverage_bases += other.zero_coverage_bases;
    }

    pub fn finalise(&mut self, exclude_zero_coverage: bool) {
        self.calc_total_filtered();
        self.calc_statistics(exclude_zero_coverage);
    }

    pub fn statistics(&self) -> Statistics {
        self.statistics.expect("finalised")
    }

    fn calc_total_filtered(&mut self) {
        if self.total_bases >= 0 {
            return;
        }
        self.total_bases = self.filter_counts.iter().sum();
    }

    /// calcFilteredPercentage: the type's bases over all bases, filtered and unfiltered
    pub fn calc_filtered_percentage(&mut self, filter: usize) -> f64 {
        self.calc_total_filtered();
        if self.total_bases == 0 {
            return 0.0;
        }
        self.filter_counts[filter] as f64 / self.total_bases as f64
    }

    /// calcCoverageFrequency: the share of the territory with at least this coverage
    pub fn calc_coverage_frequency(&self, level: i32, genome_territory: f64) -> f64 {
        if genome_territory <= 0.0 {
            return 0.0;
        }
        let total: i64 = self.coverage_frequency.iter().enumerate().filter(|(i, _)| *i as i64 >= level as i64).map(|(_, f)| *f).sum();
        total as f64 / genome_territory
    }

    fn calc_statistics(&mut self, exclude_zero: bool) {
        if self.statistics.is_some() {
            return;
        }
        let cf = &self.coverage_frequency;
        let mut total: i64 = 0;
        let mut total_frequency: i64 = 0;
        for (i, &f) in cf.iter().enumerate() {
            if exclude_zero && i == 0 {
                continue;
            }
            total += f;
            total_frequency = total_frequency.wrapping_add(f.wrapping_mul(i as i64));
        }
        if total == 0 {
            self.statistics = Some(Statistics { mean: 0.0, median: 0.0, mad: 0, std_deviation: 0.0 });
            return;
        }
        let mean = total_frequency as f64 / total as f64;
        // BT14: Java's rank total / 2 gives the value one below the median for an odd count (0 for a
        // single base). Corrected: rank (total + 1) / 2, the median for an odd count and the lower one of
        // the middle two for an even count.
        let median_total = if compat::java("BT14") { total / 2 } else { (total + 1) / 2 };
        let mut median: f64 = -1.0;
        let mut variance_total = 0.0f64;
        let mut cumulative: i64 = 0;
        for (i, &f) in cf.iter().enumerate() {
            if exclude_zero && i == 0 {
                continue;
            }
            if median < 0.0 {
                if cumulative + f >= median_total {
                    median = i as f64;
                } else {
                    cumulative += f;
                }
            }
            let diff = i as f64 - mean;
            // Math.pow(x, 2) is x * x
            variance_total += diff * diff * f as f64;
        }
        let std_deviation = (variance_total / total as f64).sqrt();
        // the median absolute deviation: frequencies by |i - median| (every i), ascending
        let max_diff = (cf.len() as f64 - 1.0 - median).max(median).max(0.0) as usize;
        let mut by_diff = vec![0i64; max_diff + 1];
        for (i, &f) in cf.iter().enumerate() {
            // (int)Math.round(abs(i - median)), exact for these integer values
            let d = (i as f64 - median).abs() as usize;
            by_diff[d] += f;
        }
        cumulative = 0;
        let mut mad = 0i32;
        for (d, &f) in by_diff.iter().enumerate() {
            if cumulative + f >= median_total {
                mad = d as i32;
                break;
            }
            cumulative += f;
        }
        crate::debug!(
            "mean({}) median({}) mad({mad}) totalFrequency={total_frequency} total={total}) stdDeviation({})",
            hmf_common::fmt::java_fixed(mean, 3),
            hmf_common::fmt::java_fixed(median, 1),
            hmf_common::fmt::java_fixed(std_deviation, 3)
        );
        self.statistics = Some(Statistics { mean, median, mad, std_deviation });
    }

    /// toString: counts(UNFILTERED=n LOW_MAP_QUAL=n ...)
    pub fn describe(&self) -> String {
        let parts: Vec<String> = (0..FILTER_TYPES).map(|i| format!("{}={}", FILTER_NAMES[i], self.filter_counts[i])).collect();
        format!("counts({})", parts.join(" "))
    }
}

/// CoverageMetrics.createCoverageBuckets: levels growing by 1% (at least 1) from 1, the first one being 2
pub fn create_coverage_buckets(distribution_length: usize) -> Vec<i32> {
    let mut current: i32 = 1;
    let mut buckets = Vec::new();
    while (current as usize) < distribution_length {
        current += (current as f64 * 0.01).ceil() as i32;
        buckets.push(current);
    }
    buckets
}

/// CoverageMetrics.getNextCoverageBucket: the bucket a coverage level goes to, the highest one not above
/// it. BT15: a level below the first bucket (coverage 1) goes to the first bucket, 2, in Java, so the
/// collapsed distribution has no row for 1x and counts those bases under 2x. Corrected: level 1 has its
/// own row.
pub fn next_coverage_bucket(coverage: i32, buckets: &[i32]) -> i32 {
    if coverage < buckets[0] {
        return if compat::java("BT15") { buckets[0] } else { coverage };
    }
    for i in 0..buckets.len() {
        if buckets[i] > coverage {
            return buckets[i - 1];
        }
    }
    buckets[buckets.len() - 1]
}

/// MetricsWriter.writeCoverageFrequency's rows: every level, or for a distribution longer than 3000
/// levels (a panel's max coverage) the levels collapsed into ~1% buckets
pub fn coverage_rows(coverage_frequency: &[i64]) -> Vec<hmf_common::metrics::ValueFrequency> {
    use hmf_common::metrics::ValueFrequency;
    const LONG_DISTRIBUTION_THRESHOLD: usize = 3000;
    let mut rows = Vec::new();
    if coverage_frequency.len() > LONG_DISTRIBUTION_THRESHOLD {
        let buckets = create_coverage_buckets(coverage_frequency.len());
        let mut current = 0i32;
        let mut total = coverage_frequency[0];
        for (i, &f) in coverage_frequency.iter().enumerate().skip(1) {
            let next = next_coverage_bucket(i as i32, &buckets);
            if current < next {
                rows.push(ValueFrequency::new(current, total));
                current = next;
                total = 0;
            }
            total += f;
        }
        rows.push(ValueFrequency::new(current, total));
    } else {
        for (i, &f) in coverage_frequency.iter().enumerate() {
            rows.push(ValueFrequency::new(i as i32, f));
        }
    }
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metrics_of(freq: &[i64]) -> CoverageMetrics {
        let mut m = CoverageMetrics::new(freq.len() as i32 - 1);
        m.coverage_frequency.copy_from_slice(freq);
        m
    }

    #[test]
    fn statistics_java_and_corrected() {
        // depths 1, 2, 3: Java's median rank 3 / 2 = 1 gives 1; the median is 2
        let mut m = metrics_of(&[0, 1, 1, 1]);
        compat::with_java(true, || m.finalise(false));
        assert_eq!(m.statistics().median, 1.0);
        let mut m = metrics_of(&[0, 1, 1, 1]);
        compat::with_java(false, || m.finalise(false));
        let s = m.statistics();
        assert_eq!((s.mean, s.median, s.mad), (2.0, 2.0, 1));
        assert!((s.std_deviation - (2.0f64 / 3.0).sqrt()).abs() < 1e-15);
        // a single base at depth 5: Java gives median 0
        let mut m = metrics_of(&[0, 0, 0, 0, 0, 1]);
        compat::with_java(true, || m.finalise(false));
        assert_eq!(m.statistics().median, 0.0);
        let mut m = metrics_of(&[0, 0, 0, 0, 0, 1]);
        compat::with_java(false, || m.finalise(false));
        assert_eq!(m.statistics().median, 5.0);
        let mut empty = metrics_of(&[5, 0]);
        empty.finalise(true);
        assert_eq!(empty.statistics(), Statistics { mean: 0.0, median: 0.0, mad: 0, std_deviation: 0.0 });
    }

    #[test]
    fn histograms() {
        let mut depths: Vec<u32> = vec![0; 300];
        for (i, d) in depths.iter_mut().enumerate().skip(130) {
            *d = (i % 7) as u32;
        }
        let mut freq = vec![0i64; 8];
        let zeros = histogram(&depths, &mut freq);
        let mut expect = vec![0i64; 8];
        for &d in &depths {
            expect[d as usize] += 1;
        }
        assert_eq!(freq, expect);
        assert_eq!(zeros, expect[0]);
    }

    #[test]
    fn coverage_buckets() {
        let b = create_coverage_buckets(10001);
        assert_eq!(&b[..3], &[2, 3, 4]);
        assert!(*b.last().unwrap() >= 10001);
        // 1% steps above 100
        assert_eq!(b[b.iter().position(|x| *x == 100).unwrap() + 1], 101);
        assert_eq!(b[b.iter().position(|x| *x == 101).unwrap() + 1], 103);
        compat::with_java(true, || {
            assert_eq!(next_coverage_bucket(1, &b), 2);
            assert_eq!(next_coverage_bucket(2, &b), 2);
            assert_eq!(next_coverage_bucket(102, &b), 101);
            let mut freq = vec![0i64; 10001];
            freq[0] = 7;
            freq[1] = 3;
            freq[2] = 4;
            let rows = coverage_rows(&freq);
            assert_eq!((rows[0].value, rows[0].count), (0, 7));
            assert_eq!((rows[1].value, rows[1].count), (2, 7));
        });
        compat::with_java(false, || {
            let mut freq = vec![0i64; 10001];
            freq[1] = 3;
            freq[2] = 4;
            let rows = coverage_rows(&freq);
            assert_eq!((rows[1].value, rows[1].count), (1, 3));
            assert_eq!((rows[2].value, rows[2].count), (2, 4));
        });
        let short = coverage_rows(&[1, 2, 3]);
        assert_eq!(short.len(), 3);
    }
}

#[cfg(test)]
mod java_harness {
    //! Replays the cases of harness/BtHarness.java (testdata/java_harness.tsv): Java's CoverageMetrics
    //! statistics, filtered percentages and coverage levels, and its coverage buckets, bit for bit.
    use super::*;

    fn nums<T: std::str::FromStr>(s: &str) -> Vec<T>
    where
        T::Err: std::fmt::Debug,
    {
        s.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect()
    }

    #[test]
    fn statistics_and_buckets_match_java() {
        let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/java_harness.tsv")).unwrap();
        let (mut stats, mut buckets) = (0, 0);
        compat::with_java(true, || {
            for line in text.lines() {
                let f: Vec<&str> = line.split('\t').collect();
                match f[0] {
                    "STATS" => {
                        let max_cov: i32 = f[1].parse().unwrap();
                        let exclude_zero = f[2] == "true";
                        let mut m = CoverageMetrics::new(max_cov);
                        m.zero_coverage_bases = f[3].parse().unwrap();
                        m.coverage_bases = f[4].parse().unwrap();
                        m.coverage_frequency = nums(f[5]);
                        let fc: Vec<i64> = nums(f[6]);
                        m.filter_counts.copy_from_slice(&fc);
                        m.finalise(exclude_zero);
                        let s = m.statistics();
                        assert_eq!(s.mean.to_bits() as i64, f[7].parse::<i64>().unwrap(), "mean: {line:.80}");
                        assert_eq!(s.median.to_bits() as i64, f[8].parse::<i64>().unwrap(), "median");
                        assert_eq!(s.mad, f[9].parse::<i32>().unwrap(), "mad");
                        assert_eq!(s.std_deviation.to_bits() as i64, f[10].parse::<i64>().unwrap(), "sd");
                        let pcts: Vec<i64> = nums(f[11]);
                        for (t, p) in pcts.iter().enumerate() {
                            assert_eq!(m.calc_filtered_percentage(t).to_bits() as i64, *p, "filter {t}");
                        }
                        let levels: Vec<i64> = nums(f[12]);
                        let territory = m.genome_territory() as f64;
                        for (k, l) in crate::metrics::COVERAGE_LEVELS.iter().enumerate() {
                            assert_eq!(m.calc_coverage_frequency(*l, territory).to_bits() as i64, levels[k], "level {l}");
                        }
                        stats += 1;
                    }
                    "BUCKETS" => {
                        let n: usize = f[1].parse().unwrap();
                        let b = create_coverage_buckets(n);
                        assert_eq!(b, nums::<i32>(f[2]));
                        let next: Vec<i32> = nums(f[3]);
                        for i in 1..n {
                            assert_eq!(next_coverage_bucket(i as i32, &b), next[i - 1], "bucket of {i}");
                        }
                        buckets += 1;
                    }
                    _ => {}
                }
            }
        });
        assert_eq!((stats, buckets), (400, 4));
    }
}
