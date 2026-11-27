//! The reseg package (new in PURPLE 4.5): the primary tumor-ratio peak (RatioPeakAnalyser), GC
//! normalisation of tumor ratios (GcNormaliser), the segmentation penalty from ratio-change troughs
//! (PenaltyCalculator), supersegments of adjacent NONE-supported diploid segments
//! (SupersegmentIdentifier) and their least-penalty re-partition (SupersegmentResegmenter,
//! SkippedSegmentReinjector).

use crate::constants::*;
use crate::region::ObservedRegion;
use crate::types::*;

#[derive(Clone, Copy, Debug)]
pub struct Bucket {
    pub level: f64,
    pub value: f64,
}

/// RatioBucketSeries
pub struct RatioBucketSeries {
    min: f64,
    step: f64,
    counts: Vec<f64>,
    max: f64,
}

impl RatioBucketSeries {
    pub fn new() -> RatioBucketSeries {
        Self::with_range(RESEG_RATIO_BUCKET_MIN, RESEG_RATIO_BUCKET_MAX, RESEG_RATIO_BUCKET_STEP)
    }
    pub fn with_range(min: f64, max: f64, step: f64) -> RatioBucketSeries {
        let n = jround((max - min) / step) as usize + 1;
        RatioBucketSeries { min, step, counts: vec![0.0; n], max }
    }
    pub fn add_count(&mut self, value: f64, weight: f64) {
        if value < self.min || value > self.max {
            return;
        }
        let n = self.counts.len() as i64;
        let index = jround((value - self.min) / self.step).clamp(0, n - 1) as usize;
        self.counts[index] += weight;
    }
    pub fn build_smoothed_series(&self) -> Vec<Bucket> {
        let window = RESEG_ROLLING_AVG_WINDOW as usize;
        let half = window / 2;
        let mut smoothed: Vec<Bucket> = Vec::new();
        let n = self.counts.len();
        if n > half {
            for i in half..n - half {
                let mut sum = 0.0;
                for k in i - half..=i + half {
                    sum += self.counts[k];
                }
                let avg = sum / window as f64;
                let level = jround((self.min + i as f64 * self.step) / self.step) as f64 * self.step;
                smoothed.push(Bucket { level, value: avg });
            }
        }
        let mut collapsed: Vec<Bucket> = Vec::new();
        for b in smoothed {
            if collapsed.last().is_some_and(|l| l.value == b.value) {
                continue;
            }
            collapsed.push(b);
        }
        collapsed
    }
}

/// PeakTroughData
#[derive(Clone, Copy, Debug)]
pub struct PeakTroughData {
    pub level: f64,
    pub value: f64,
    pub support_below_level: f64,
    pub support_above_level: f64,
}

fn find_local_indices(series: &[Bucket], peaks: bool) -> Vec<usize> {
    let mut out = Vec::new();
    for i in 1..series.len().saturating_sub(1) {
        let (prev, curr, next) = (series[i - 1].value, series[i].value, series[i + 1].value);
        if peaks {
            if curr > prev && curr > next {
                out.push(i);
            }
        } else if curr < prev && curr < next {
            out.push(i);
        }
    }
    out
}

fn nearest_below(sorted: &[usize], center: usize, ok: impl Fn(usize) -> bool) -> Option<usize> {
    sorted.iter().rev().copied().filter(|&c| c < center).find(|&c| ok(c))
}
fn nearest_above(sorted: &[usize], center: usize, ok: impl Fn(usize) -> bool) -> Option<usize> {
    sorted.iter().copied().filter(|&c| c > center).find(|&c| ok(c))
}

/// PeakTroughFinder.consolidateResults
fn consolidate(sorted: &[PeakTroughData], keep_both: impl Fn(&PeakTroughData, &PeakTroughData) -> bool, keep_higher: bool) -> Vec<PeakTroughData> {
    let mut current: Vec<PeakTroughData> = sorted.to_vec();
    loop {
        let mut any = false;
        let mut drop = vec![false; current.len()];
        for i in 0..current.len().saturating_sub(1) {
            let (l, r) = (&current[i], &current[i + 1]);
            if keep_both(l, r) {
                continue;
            }
            any = true;
            let drop_left = if keep_higher { r.value > l.value } else { r.value < l.value };
            drop[if drop_left { i } else { i + 1 }] = true;
        }
        if !any {
            return current;
        }
        current = current.iter().enumerate().filter(|(i, _)| !drop[*i]).map(|(_, p)| *p).collect();
    }
}

pub struct RatioPeakResult {
    pub left_bound: f64,
    pub right_bound: f64,
}

/// RatioPeakAnalyser.findRatioPeak
pub fn find_ratio_peak(segments: &[ObservedRegion]) -> Option<RatioPeakResult> {
    let mut raw = RatioBucketSeries::new();
    for s in segments {
        if s.germline_status != GermlineStatus::Diploid || s.baf_count <= 0 {
            continue;
        }
        raw.add_count(round(s.observed_tumor_ratio, 2), s.baf_count as f64);
    }
    let ratios = raw.build_smoothed_series();
    if ratios.len() < 2 {
        return None;
    }
    let peaks = find_local_indices(&ratios, true);
    let total: f64 = hmf_common::java_math::double_stream_sum(ratios.iter().map(|b| b.value));
    let all: Vec<usize> = (0..ratios.len()).collect();
    let valid_support = |peak: f64, other: f64| other == 0.0 || peak / other >= RESEG_PEAK_MIN_RATIO;
    let mut valid: Vec<PeakTroughData> = Vec::new();
    for &pi in &peaks {
        let pv = ratios[pi].value;
        if total <= 0.0 || pv / total < RESEG_PEAK_MIN_PROPORTION {
            continue;
        }
        let left = nearest_below(&all, pi, |o| valid_support(pv, ratios[o].value));
        let right = nearest_above(&all, pi, |o| valid_support(pv, ratios[o].value));
        let (Some(l), Some(r)) = (left, right) else { continue };
        valid.push(PeakTroughData { level: ratios[pi].level, value: pv, support_below_level: ratios[l].level, support_above_level: ratios[r].level });
    }
    let consolidated = consolidate(
        &valid,
        |l, r| r.level - l.support_above_level >= RESEG_PEAK_MIN_GAP && r.support_below_level - l.level >= RESEG_PEAK_MIN_GAP,
        true,
    );
    let mut top = consolidated.clone();
    // sorted by value descending (stable)
    top.sort_by(|a, b| dcmp(b.value, a.value));
    top.truncate(RESEG_MAX_PEAKS as usize);
    if top.len() < 2 {
        return None;
    }
    let primary = top[0];
    let others = &top[1..];
    // Stream.max / min with comparingDouble (BinaryOperator.maxBy / minBy): the first maximum, the first minimum
    let mut left: Option<PeakTroughData> = None;
    for p in others.iter().filter(|p| p.level < primary.level) {
        if left.is_none_or(|l| dcmp(p.level, l.level) == std::cmp::Ordering::Greater) {
            left = Some(*p);
        }
    }
    let mut right: Option<PeakTroughData> = None;
    for p in others.iter().filter(|p| p.level > primary.level) {
        if right.is_none_or(|r| dcmp(p.level, r.level) == std::cmp::Ordering::Less) {
            right = Some(*p);
        }
    }
    let left_bound = left.map(|p| round((p.support_above_level + 2.0 * primary.support_below_level) / 3.0, 4)).unwrap_or(RESEG_RATIO_BUCKET_MIN);
    let right_bound = right.map(|p| round((p.support_below_level + 2.0 * primary.support_above_level) / 3.0, 4)).unwrap_or(RESEG_RATIO_BUCKET_MAX);
    Some(RatioPeakResult { left_bound, right_bound })
}

// ---------------------------------------------------------------------------------------------
// GC normalisation

#[derive(Clone, Copy, Debug, PartialEq)]
struct GcBin {
    lower: f64,
    upper: f64,
}
impl GcBin {
    fn contains(&self, gc: f64) -> bool {
        gc > self.lower && gc <= self.upper
    }
}

fn gc_bins() -> Vec<GcBin> {
    let n = jround((RESEG_GC_BIN_MAX - RESEG_GC_BIN_MIN) / RESEG_GC_BIN_WIDTH) as i32;
    (0..n)
        .map(|i| {
            let lower = round(RESEG_GC_BIN_MIN + i as f64 * RESEG_GC_BIN_WIDTH, 2);
            GcBin { lower, upper: round(lower + RESEG_GC_BIN_WIDTH, 2) }
        })
        .collect()
}

fn find_bin(gc: f64, bins: &[GcBin]) -> Option<usize> {
    bins.iter().position(|b| b.contains(gc))
}

/// (ratio, baf count) sorted by ratio, equal ratios merged
fn ratio_baf_counts(segments: &[&ObservedRegion]) -> Vec<(f64, i32)> {
    let mut out: Vec<(f64, i32)> = Vec::new();
    for s in segments {
        match out.iter_mut().find(|(r, _)| *r == s.observed_tumor_ratio) {
            Some(e) => e.1 += s.baf_count,
            None => out.push((s.observed_tumor_ratio, s.baf_count)),
        }
    }
    out.sort_by(|a, b| dcmp(a.0, b.0));
    out
}

fn weight_ratio_across_percentile(counts: &[(f64, i32)], percentile: f64) -> Option<f64> {
    let total: f64 = hmf_common::java_math::double_stream_sum(counts.iter().map(|c| c.1 as f64));
    if total <= 0.0 {
        return None;
    }
    let mut cumulative = 0.0;
    let mut lower: Option<f64> = None;
    let mut upper: Option<f64> = None;
    for (ratio, count) in counts {
        cumulative += *count as f64;
        let prop = cumulative / total;
        if lower.is_none() && prop >= percentile {
            lower = Some(*ratio);
        }
        if upper.is_none() && prop > percentile {
            upper = Some(*ratio);
        }
        if lower.is_some() && upper.is_some() {
            break;
        }
    }
    Some((lower? + upper?) / 2.0)
}

fn suppress_sparse_values(grid: &mut [Vec<f64>], bin_proportion: &[f64]) {
    let bins = grid.len();
    let percentiles = grid[0].len();
    let mut any_nan = true;
    let mut pass = 0;
    while any_nan && pass < 20 {
        pass += 1;
        for p in 0..percentiles {
            for i in 1..bins {
                let (curr, prev) = (grid[i][p], grid[i - 1][p]);
                if !curr.is_nan() && !prev.is_nan() && (curr - prev).abs() > RESEG_GC_SPARSE_MIN_JUMP && bin_proportion[i] < RESEG_GC_SPARSE_MAX_PROPORTION {
                    grid[i][p] = f64::NAN;
                }
            }
        }
        any_nan = grid.iter().any(|row| row.iter().any(|v| v.is_nan()));
        for p in 0..percentiles {
            for i in 1..bins {
                if grid[i][p].is_nan() {
                    grid[i][p] = grid[i - 1][p];
                }
            }
        }
    }
}

fn gc_adjustments(regions: &[&ObservedRegion], bins: &[GcBin]) -> Vec<f64> {
    let total_weight: f64 = hmf_common::java_math::double_stream_sum(regions.iter().map(|r| r.baf_count as f64));
    let mut bin_proportion = vec![0.0; bins.len()];
    let mut by_bin: Vec<Vec<&ObservedRegion>> = Vec::with_capacity(bins.len());
    for (i, bin) in bins.iter().enumerate() {
        let mut list = Vec::new();
        let mut weight = 0.0;
        for r in regions {
            if bin.contains(r.gc_content) {
                weight += r.baf_count as f64;
                list.push(*r);
            }
        }
        bin_proportion[i] = if total_weight > 0.0 { weight / total_weight } else { 0.0 };
        by_bin.push(list);
    }
    let all_counts = ratio_baf_counts(regions);
    let mut grid = vec![vec![0.0; RESEG_GC_PERCENTILES.len()]; bins.len()];
    for (p, &percentile) in RESEG_GC_PERCENTILES.iter().enumerate() {
        for i in 0..bins.len() {
            let gc_counts = ratio_baf_counts(&by_bin[i]);
            let all = weight_ratio_across_percentile(&all_counts, percentile);
            let gc = weight_ratio_across_percentile(&gc_counts, percentile);
            grid[i][p] = match (all, gc) {
                (Some(a), Some(g)) => g / a,
                _ => f64::NAN,
            };
        }
    }
    suppress_sparse_values(&mut grid, &bin_proportion);
    grid.iter()
        .map(|row| {
            let valid: Vec<f64> = row.iter().copied().filter(|v| !v.is_nan()).collect();
            if valid.is_empty() { f64::NAN } else { median(&valid) }
        })
        .collect()
}

/// GcNormaliser.normaliseRegions
pub fn normalise_regions(all: &[ObservedRegion], peak: Option<&RatioPeakResult>) -> Vec<ObservedRegion> {
    let bins = gc_bins();
    let mut lower: Option<usize> = None;
    let mut high: Option<usize> = None;
    let mut observed: Vec<usize> = Vec::new();
    for s in all {
        if let Some(b) = find_bin(s.gc_content, &bins) {
            if !observed.contains(&b) {
                observed.push(b);
            }
            if lower.is_none_or(|l| bins[b].lower < bins[l].lower) {
                lower = Some(b);
            }
            if high.is_none_or(|h| bins[b].upper > bins[h].upper) {
                high = Some(b);
            }
            if observed.len() == bins.len() {
                break;
            }
        }
    }
    let lower = lower.unwrap_or(0);
    let high = high.unwrap_or(bins.len() - 1);
    let valid: Vec<&ObservedRegion> = all
        .iter()
        .filter(|r| {
            r.baf_count > 0
                && r.germline_status == GermlineStatus::Diploid
                && peak.is_none_or(|p| !(r.observed_tumor_ratio < p.left_bound || r.observed_tumor_ratio > p.right_bound))
        })
        .collect();
    let adjustments = gc_adjustments(&valid, &bins);
    all.iter()
        .map(|s| {
            let bin = find_bin(s.gc_content, &bins).unwrap_or(if s.gc_content <= bins[lower].lower { lower } else { high });
            let mut adjustment = adjustments[bin];
            if adjustment.is_nan() || adjustment <= 0.0 {
                adjustment = 1.0;
            }
            let mut r = s.clone();
            r.observed_tumor_ratio = s.observed_tumor_ratio / adjustment;
            r
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// penalty

fn calc_penalty(ratio: f64) -> f64 {
    round(jpow(ratio, 2.0) * 2.0, 4)
}

/// PenaltyCalculator.calculatePenalty
pub fn calculate_penalty(segments: &[ObservedRegion]) -> f64 {
    let mut raw = RatioBucketSeries::new();
    for i in 1..segments.len() {
        let (prev, curr) = (&segments[i - 1], &segments[i]);
        if curr.observed_tumor_ratio <= 0.0 || prev.observed_tumor_ratio <= 0.0 {
            continue;
        }
        let change = jmin(round((curr.observed_tumor_ratio - prev.observed_tumor_ratio).abs(), 2), RESEG_RATIO_BUCKET_MAX);
        raw.add_count(change, 1.0);
    }
    let diffs = raw.build_smoothed_series();
    if diffs.len() < 2 {
        return calc_penalty(RESEG_MAX_SEGMENTATION_PENALTY_RATIO);
    }
    let mut peaks = find_local_indices(&diffs, true);
    let mut troughs = find_local_indices(&diffs, false);
    if diffs[0].value < diffs[1].value {
        troughs.insert(0, 0);
    } else if diffs[0].value > diffs[1].value {
        peaks.insert(0, 0);
    }
    let valid_support = |trough: f64, peak: f64| {
        if peak - trough < RESEG_TROUGH_MIN_DIFF {
            return false;
        }
        trough == 0.0 || peak / trough >= RESEG_TROUGH_MIN_RATIO
    };
    let mut valid: Vec<PeakTroughData> = Vec::new();
    for &ti in &troughs {
        let tv = diffs[ti].value;
        let prev = nearest_below(&peaks, ti, |pi| valid_support(tv, diffs[pi].value));
        let next = nearest_above(&peaks, ti, |pi| valid_support(tv, diffs[pi].value));
        let (Some(p), Some(n)) = (prev, next) else { continue };
        valid.push(PeakTroughData { level: diffs[ti].level, value: tv, support_below_level: diffs[p].level, support_above_level: diffs[n].level });
    }
    let consolidated = consolidate(&valid, |l, r| r.level - l.level >= RESEG_TROUGH_MIN_GAP, false);
    let lowest = consolidated.iter().map(|t| t.level).fold(None, |m: Option<f64>, v| Some(m.map_or(v, |m| jmin(m, v))));
    let lowest = lowest.unwrap_or(RESEG_MAX_SEGMENTATION_PENALTY_RATIO);
    calc_penalty(jmin(lowest, RESEG_MAX_SEGMENTATION_PENALTY_RATIO))
}

// ---------------------------------------------------------------------------------------------
// supersegments

pub struct Supersegment {
    pub both_none: Vec<ObservedRegion>,
    pub skippable: Vec<ObservedRegion>,
}

/// SupersegmentIdentifier.identify: (supersegments, excluded segments)
pub fn identify_supersegments(segments: &[ObservedRegion]) -> (Vec<Supersegment>, Vec<ObservedRegion>) {
    let n = segments.len();
    let mut both_none = vec![false; n];
    let mut can_skip = vec![false; n];
    for i in 0..n {
        let s = &segments[i];
        let other = if i + 1 < n { Some(segments[i + 1].support) } else { None };
        both_none[i] = s.support == SegmentSupport::None && other == Some(SegmentSupport::None) && s.germline_status == GermlineStatus::Diploid;
        let is_deletion = matches!(s.germline_status, GermlineStatus::HomDeletion | GermlineStatus::HetDeletion);
        let unknown_bridge = s.germline_status == GermlineStatus::Unknown
            && other.is_some()
            && Some(s.support) == other
            && matches!(s.support, SegmentSupport::Del | SegmentSupport::Dup);
        can_skip[i] = is_deletion || unknown_bridge;
    }
    let raw_new: Vec<bool> = (0..n).map(|i| both_none[i] && !(i > 0 && both_none[i - 1])).collect();
    let skip_over: Vec<bool> = (0..n).map(|i| can_skip[i] && i > 0 && both_none[i - 1] && i + 1 < n && raw_new[i + 1]).collect();
    let mut new_super = raw_new.clone();
    for i in 0..n {
        if skip_over[i] && i + 1 < n {
            new_super[i + 1] = false;
        }
    }
    let mut supers: Vec<Supersegment> = Vec::new();
    let mut excluded: Vec<ObservedRegion> = Vec::new();
    let mut current: Option<Supersegment> = None;
    for i in 0..n {
        let s = &segments[i];
        if new_super[i] {
            if let Some(c) = current.take() {
                supers.push(c);
            }
            current = Some(Supersegment { both_none: Vec::new(), skippable: Vec::new() });
        }
        if both_none[i] {
            current.as_mut().expect("both-none segment before any supersegment").both_none.push(s.clone());
        } else {
            excluded.push(s.clone());
            if skip_over[i] {
                current.as_mut().expect("skippable segment before any supersegment").skippable.push(s.clone());
            }
        }
    }
    if let Some(c) = current {
        supers.push(c);
    }
    (supers, excluded)
}

/// SupersegmentResegmenter.deviationPenalty
pub fn deviation_penalty(members: &[ObservedRegion]) -> f64 {
    let m = members.len();
    let mean = if m == 0 { 0.0 } else { hmf_common::java_math::double_stream_sum(members.iter().map(|x| x.observed_tumor_ratio)) / m as f64 };
    let mut total = 0.0;
    for i in 0..m {
        let ratio = members[i].observed_tumor_ratio;
        let prev = if i > 0 { (ratio - members[i - 1].observed_tumor_ratio).abs() } else { 0.0 };
        let next = if i + 1 < m { (ratio - members[i + 1].observed_tumor_ratio).abs() } else { 0.0 };
        let mean_diff = (ratio - mean).abs();
        let dev = jmax(prev, jmax(next, mean_diff));
        total += dev * dev;
    }
    total
}

fn aggregate(regions: &[ObservedRegion]) -> ObservedRegion {
    let use_baf = regions.iter().any(|r| r.baf_count >= RESEG_BAF_WEIGHT_THRESHOLD);
    let (mut baf_total, mut depth_total) = (0i32, 0i32);
    let (mut w_baf, mut w_ratio) = (0.0, 0.0);
    for r in regions {
        let count = if use_baf { r.baf_count } else { r.depth_window_count };
        depth_total += r.depth_window_count;
        baf_total += r.baf_count;
        w_baf += r.observed_baf * count as f64;
        w_ratio += r.observed_tumor_ratio * count as f64;
    }
    let counts = if use_baf { baf_total } else { depth_total } as f64;
    let mut out = regions[0].clone();
    out.end = regions[regions.len() - 1].end;
    out.baf_count = baf_total;
    out.observed_baf = w_baf / counts;
    out.depth_window_count = depth_total;
    out.observed_tumor_ratio = w_ratio / counts;
    out
}

fn combinations(range: usize, k: usize) -> Vec<Vec<usize>> {
    fn helper(out: &mut Vec<Vec<usize>>, combo: &mut Vec<usize>, start: usize, range: usize, k: usize) {
        if combo.len() == k {
            out.push(combo.clone());
            return;
        }
        let remaining = k - combo.len() - 1;
        let mut v = start;
        while v + remaining <= range {
            combo.push(v);
            helper(out, combo, v + 1, range, k);
            combo.pop();
            v += 1;
        }
    }
    let mut out = Vec::new();
    helper(&mut out, &mut Vec::with_capacity(k), 1, range, k);
    out
}

fn find_best_partition(members: &[ObservedRegion], penalty: f64) -> Option<Vec<usize>> {
    let n = members.len();
    let mut lowest = deviation_penalty(members);
    let mut best: Option<Vec<usize>> = None;
    for k in 1..n {
        let partition_cost = k as f64 * penalty;
        if partition_cost >= lowest {
            break;
        }
        for cuts in combinations(n - 1, k) {
            let mut padded = Vec::with_capacity(k + 2);
            padded.push(0);
            padded.extend(cuts);
            padded.push(n);
            let mut deviation = 0.0;
            for j in 0..padded.len() - 1 {
                deviation += deviation_penalty(&members[padded[j]..padded[j + 1]]);
                if partition_cost + deviation >= lowest {
                    break;
                }
            }
            if partition_cost + deviation < lowest {
                lowest = partition_cost + deviation;
                best = Some(padded);
            }
        }
    }
    best
}

/// SupersegmentResegmenter.resegment
pub fn resegment(members: &[ObservedRegion], penalty: f64) -> Vec<ObservedRegion> {
    if members.len() as i32 > RESEG_MAX_SUPERSEGMENT_SIZE_FOR_BRUTE_FORCE {
        let max = members.iter().map(|m| m.observed_tumor_ratio).fold(f64::NEG_INFINITY, jmax);
        let min = members.iter().map(|m| m.observed_tumor_ratio).fold(f64::INFINITY, jmin);
        if max - min > (penalty / 2.0).sqrt() {
            let mut best_index = 1;
            let mut best_diff = -1.0;
            for i in 1..members.len() {
                let diff = (members[i].observed_tumor_ratio - members[i - 1].observed_tumor_ratio).abs();
                if diff >= best_diff {
                    best_diff = diff;
                    best_index = i;
                }
            }
            let mut out = resegment(&members[..best_index], penalty);
            out.extend(resegment(&members[best_index..], penalty));
            return out;
        }
        return vec![aggregate(members)];
    }
    match find_best_partition(members, penalty) {
        None => vec![aggregate(members)],
        Some(points) => (0..points.len() - 1).map(|j| aggregate(&members[points[j]..points[j + 1]])).collect(),
    }
}

/// SkippedSegmentReinjector.formBlendedRegion
fn form_blended_region(regions: &[&ObservedRegion], sub: &ObservedRegion) -> ObservedRegion {
    let (mut baf_total, mut depth_total) = (0i32, 0i32);
    let (mut w_norm, mut w_unnorm, mut gc_total) = (0.0, 0.0, 0.0);
    let mut total_length = 0i32;
    let mut ratio_support = false;
    let mut min_start = -1;
    let mut max_start = 0;
    for r in regions {
        baf_total += r.baf_count;
        depth_total += r.depth_window_count;
        w_norm += r.depth_window_count as f64 * r.observed_normal_ratio;
        w_unnorm += r.depth_window_count as f64 * r.unnormalised_observed_normal_ratio;
        ratio_support |= r.ratio_support;
        let len = r.end - r.start + 1;
        total_length += len;
        gc_total += r.gc_content * len as f64;
        min_start = if min_start < 0 { r.min_start } else { min_start.min(r.min_start) };
        max_start = max_start.max(r.max_start);
    }
    let mut out = regions[0].clone();
    out.end = regions[regions.len() - 1].end;
    out.min_start = min_start;
    out.max_start = max_start;
    out.observed_tumor_ratio = sub.observed_tumor_ratio;
    out.observed_baf = sub.observed_baf;
    out.baf_count = baf_total;
    out.depth_window_count = depth_total;
    out.observed_normal_ratio = w_norm / depth_total as f64;
    out.unnormalised_observed_normal_ratio = w_unnorm / depth_total as f64;
    out.ratio_support = ratio_support;
    out.gc_content = gc_total / total_length as f64;
    out
}

fn build_reinjected(sub: &ObservedRegion, ss: &Supersegment, start: i32, end: i32) -> ObservedRegion {
    let relevant: Vec<&ObservedRegion> = ss.both_none.iter().filter(|s| s.start >= start && s.end <= end).collect();
    if relevant.is_empty() {
        panic!("no bothNone members found for reinjected range({start}-{end})");
    }
    if relevant.len() == 1 {
        let mut r = relevant[0].clone();
        r.observed_tumor_ratio = sub.observed_tumor_ratio;
        r.observed_baf = sub.observed_baf;
        return r;
    }
    form_blended_region(&relevant, sub)
}

/// SkippedSegmentReinjector.reinjectSkippedRegions
pub fn reinject_skipped_regions(resegmented: Vec<ObservedRegion>, ss: &Supersegment) -> Vec<ObservedRegion> {
    if ss.skippable.is_empty() {
        if resegmented.len() == ss.both_none.len() {
            return resegmented;
        }
        if resegmented.len() == 1 {
            let all: Vec<&ObservedRegion> = ss.both_none.iter().collect();
            return vec![form_blended_region(&all, &resegmented[0])];
        }
    }
    let mut result = Vec::new();
    for sub in &resegmented {
        let mut skips: Vec<&ObservedRegion> = ss.skippable.iter().filter(|s| s.start > sub.start && s.end < sub.end).collect();
        skips.sort_by_key(|s| s.start);
        if skips.is_empty() {
            result.push(build_reinjected(sub, ss, sub.start, sub.end));
            continue;
        }
        let mut ranges = vec![(sub.start, skips[0].start - 1)];
        for i in 0..skips.len() - 1 {
            ranges.push((skips[i].end + 1, skips[i + 1].start - 1));
        }
        ranges.push((skips[skips.len() - 1].end + 1, sub.end));
        for (s, e) in ranges {
            result.push(build_reinjected(sub, ss, s, e));
        }
    }
    result
}

/// Resegmentation.run
pub fn run(input: &[ObservedRegion], threads: usize) -> Vec<ObservedRegion> {
    let peak = find_ratio_peak(input);
    let normalised = normalise_regions(input, peak.as_ref());
    let penalty = calculate_penalty(&normalised);
    let (supers, excluded) = identify_supersegments(&normalised);
    let results: Vec<Vec<ObservedRegion>> = crate::parallel::map(&supers, threads, |ss| reinject_skipped_regions(resegment(&ss.both_none, penalty), ss));
    let excluded_count = excluded.len();
    let mut final_regions = excluded;
    for r in results {
        final_regions.extend(r);
    }
    final_regions.sort_by(|a, b| region_compare(&a.chromosome, a.start, &b.chromosome, b.start));
    crate::debug!(
        "resegmentation segments({} -> {}) penalty({}) superSegments({} excluded={})",
        input.len(),
        final_regions.len(),
        hmf_common::fmt::java_fixed(penalty, 4),
        supers.len(),
        excluded_count
    );
    final_regions
}
