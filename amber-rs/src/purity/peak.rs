//! CandidatePeak, CandidatePeakEvaluation, CandidatePeakResult, SearchGrid, LocalMaximaFinder and PeakSearch.
//!
//! A candidate peak at VAF level L tests each point of sufficient depth: with n = ref + alt reads and k the
//! smaller of the two counts, the point is in the peak's homozygous band where the binomial CDF of k for
//! (n, L) lies between 0.16 and 0.84 or k/n is within the grid step of L, in the heterozygous band likewise
//! for L/2 (a point in both goes to the nearer band). The score of a level is the share of the tested points
//! in either band; the peaks are the local maxima of the scores over the grid of levels.
//!
//! Everything a test depends on is (L, step, n, k), so the tests are memoised per level; the bands of a level
//! are only kept for the levels that turn out to be maxima.

use super::binomial::binomial_cdf;
use crate::types::PositionEvidence;
use std::collections::HashMap;

pub const LOWER_CDF_BOUND_FOR_CAPTURE: f64 = 0.16;
pub const UPPER_CDF_BOUND_FOR_CAPTURE: f64 = 0.84;
pub const MINIMUM_CAPTURED_POINTS: usize = 15;
pub const PEAK_SEARCH_START: f64 = 0.005;
pub const PEAK_SEARCH_END: f64 = 0.37;
pub const PEAK_SEARCH_STEP_RATIO: f64 = 1.05;
pub const PEAK_SEARCH_INITIAL_STEP: f64 = 0.001;
pub const PEAK_SEARCH_OVERSHOOT: f64 = 0.0001;

/// Which points of a list are equal as Java's PositionEvidence.equals sees them (chromosome, position, ref,
/// alt), so that the sets Java builds from a peak's bands can be counted: `canonical[i]` is the first index
/// equal to i, or None when the list has no two equal points.
pub struct EvidenceKeys {
    canonical: Option<Vec<u32>>,
}

impl EvidenceKeys {
    pub fn new<'e>(points: &[&'e PositionEvidence]) -> EvidenceKeys {
        // points in genome order have equal points next to each other (within one position); check the
        // order, and fall back to a map when it does not hold
        let mut sorted = true;
        let mut seen_chromosomes: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for i in 1..points.len() {
            let (a, b) = (points[i - 1], points[i]);
            if a.chromosome == b.chromosome {
                if b.position < a.position {
                    sorted = false;
                    break;
                }
            } else {
                seen_chromosomes.insert(&a.chromosome);
                if seen_chromosomes.contains(b.chromosome.as_str()) {
                    sorted = false;
                    break;
                }
            }
        }
        let same = |a: &PositionEvidence, b: &PositionEvidence| a.position == b.position && a.ref_ == b.ref_ && a.alt == b.alt && a.chromosome == b.chromosome;
        let mut canonical: Option<Vec<u32>> = None;
        if sorted {
            let mut start = 0;
            while start < points.len() {
                let mut end = start + 1;
                while end < points.len() && points[end].position == points[start].position && points[end].chromosome == points[start].chromosome {
                    end += 1;
                }
                for j in start + 1..end {
                    if let Some(first) = (start..j).find(|&k| same(points[k], points[j])) {
                        let c = canonical.get_or_insert_with(|| (0..points.len() as u32).collect());
                        c[j] = c[first];
                    }
                }
                start = end;
            }
        } else {
            let mut first: HashMap<(&str, i32, u8, u8), u32> = HashMap::new();
            for (i, p) in points.iter().enumerate() {
                let key = (p.chromosome.as_str(), p.position, p.ref_ as u8, p.alt as u8);
                match first.get(&key) {
                    Some(&f) => {
                        let c = canonical.get_or_insert_with(|| (0..points.len() as u32).collect());
                        c[i] = f;
                    }
                    None => {
                        first.insert(key, i as u32);
                    }
                }
            }
        }
        EvidenceKeys { canonical }
    }

    pub fn has_duplicates(&self) -> bool {
        self.canonical.is_some()
    }

    /// the points of `bands` (concatenated in order) as a set: the first of each group of equal points
    pub fn distinct(&self, bands: &[&[u32]]) -> Vec<u32> {
        match &self.canonical {
            None => bands.iter().flat_map(|b| b.iter().copied()).collect(),
            Some(c) => {
                let mut seen = std::collections::HashSet::new();
                bands.iter().flat_map(|b| b.iter().copied()).filter(|&i| seen.insert(c[i as usize])).collect()
            }
        }
    }
}

/// The band a point falls in at one level
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Band {
    None,
    Homozygous,
    Heterozygous,
}

/// CandidatePeak.RangeStepPeak.isCaptured
pub fn is_captured(depth: i32, count: i32, vaf: f64, step: f64, symmetric_vaf: f64) -> bool {
    let cdf = binomial_cdf(depth, vaf, count);
    let by_cdf = cdf > LOWER_CDF_BOUND_FOR_CAPTURE && cdf < UPPER_CDF_BOUND_FOR_CAPTURE;
    let by_step = (symmetric_vaf - vaf).abs() < step;
    by_cdf || by_step
}

/// n = ref + alt and k, the alt count or the ref count if that is smaller (CandidatePeak.test)
fn depth_and_count(e: &PositionEvidence) -> (i32, i32) {
    let n = e.ref_support.wrapping_add(e.alt_support);
    let mut k = e.alt_support;
    if k > n / 2 {
        k = n - k;
    }
    (n, k)
}

/// CandidatePeak: a VAF level and the points in its homozygous and heterozygous bands, as indices into the
/// list of points the peak was tested on. After the purity analysis the bands are dropped; the counts and
/// the statistics that decided the peak's class remain.
#[derive(Clone, Debug)]
pub struct CandidatePeak {
    level: f64,
    step_to_next_level: f64,
    tested: usize,
    homozygous: Vec<u32>,
    heterozygous: Vec<u32>,
    /// the statistics of the purity analysis (None for a peak that has not been through it)
    pub(crate) stats: Option<PeakStats>,
}

/// What TumorOnlyPurityAnalysis worked out for a peak
#[derive(Clone, Copy, Debug, Default)]
pub struct PeakStats {
    pub score: f64,
    pub homozygous_count: usize,
    pub heterozygous_count: usize,
    pub captured: usize,
    pub gnomad: super::gnomad::GnomadCheck,
    pub homozygous_proportion: f64,
    pub chr_arm_auc: f64,
    pub mutation_auc: f64,
    pub chr_arm_lower_bound: f64,
    /// whether Java's int product overflowed in either AUC (AM53)
    pub auc_overflow: bool,
}

impl CandidatePeak {
    /// CandidatePeak(level): step 0.00001
    pub fn new(level: f64) -> CandidatePeak {
        CandidatePeak::with_step(level, 0.00001)
    }
    pub fn with_step(level: f64, step_to_next_level: f64) -> CandidatePeak {
        CandidatePeak { level, step_to_next_level, tested: 0, homozygous: Vec::new(), heterozygous: Vec::new(), stats: None }
    }
    pub fn vaf(&self) -> f64 {
        self.level
    }
    pub fn step_to_next_level(&self) -> f64 {
        self.step_to_next_level
    }
    pub fn stats(&self) -> Option<&PeakStats> {
        self.stats.as_ref()
    }

    /// hasSufficientDepthForEventDetection: het points above the depth-2 noise floor
    pub fn has_sufficient_depth_for_event_detection(&self, e: &PositionEvidence) -> bool {
        sufficient_depth(self.level, e.ref_support.wrapping_add(e.alt_support))
    }

    /// test: puts the point `index` in a band, or in none
    pub fn test(&mut self, index: u32, e: &PositionEvidence) {
        self.tested += 1;
        let (n, k) = depth_and_count(e);
        match band(self.level, self.step_to_next_level, n, k, e.symmetric_vaf()) {
            Band::Homozygous => self.homozygous.push(index),
            Band::Heterozygous => self.heterozygous.push(index),
            Band::None => {}
        }
    }

    pub fn tested(&self) -> usize {
        self.tested
    }
    /// the homozygous band, a list (Java's homozygousEvidencePoints() makes a set of it)
    pub fn homozygous_points(&self) -> &[u32] {
        &self.homozygous
    }
    pub fn heterozygous_points(&self) -> &[u32] {
        &self.heterozygous
    }
    /// allCapturedPoints, as a list of the distinct points: the homozygous band, then the heterozygous
    pub fn all_captured_points(&self, keys: &EvidenceKeys) -> Vec<u32> {
        keys.distinct(&[&self.homozygous, &self.heterozygous])
    }
    pub fn number_of_captured_evidence_points(&self, keys: &EvidenceKeys) -> usize {
        if keys.has_duplicates() { self.all_captured_points(keys).len() } else { self.homozygous.len() + self.heterozygous.len() }
    }
    /// homozygousProportion: the homozygous band's size (with any repeats) over the number of distinct
    /// captured points; NaN with none
    pub fn homozygous_proportion(&self, keys: &EvidenceKeys) -> f64 {
        self.homozygous.len() as f64 / self.number_of_captured_evidence_points(keys) as f64
    }

    pub(crate) fn drop_bands(&mut self) {
        self.homozygous = Vec::new();
        self.heterozygous = Vec::new();
    }
}

impl std::fmt::Display for CandidatePeak {
    /// toString, in Locale.ENGLISH
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        use hmf_common::fmt::java_fixed;
        write!(
            f,
            "VafLevel{{vaf={}, step={}, tested: {}, homozygous: {}, heterozygous: {}}}",
            java_fixed(self.level, 2),
            java_fixed(self.step_to_next_level, 2),
            self.tested,
            self.homozygous.len(),
            self.heterozygous.len()
        )
    }
}

fn sufficient_depth(level: f64, n: i32) -> bool {
    binomial_cdf(n, level / 2.0, 2) < LOWER_CDF_BOUND_FOR_CAPTURE
}

/// CandidatePeak.test's choice of band for a point with n reads, k of the minor allele
fn band(level: f64, step: f64, n: i32, k: i32, symmetric_vaf: f64) -> Band {
    let hom = is_captured(n, k, level, step, symmetric_vaf);
    let het = is_captured(n, k, level / 2.0, step / 2.0, symmetric_vaf);
    if hom {
        if het {
            // add it to the band to which it is closest
            if symmetric_vaf > level * 0.75 { Band::Homozygous } else { Band::Heterozygous }
        } else {
            Band::Homozygous
        }
    } else if het {
        Band::Heterozygous
    } else {
        Band::None
    }
}

/// SearchGrid.searchValuesAndSteps: (level, step) pairs; the step is the one after the level's increment
pub fn search_values_and_steps() -> Vec<(f64, f64)> {
    let mut step = PEAK_SEARCH_INITIAL_STEP;
    let mut result = Vec::with_capacity(60);
    let mut current = PEAK_SEARCH_START;
    while current <= PEAK_SEARCH_END + PEAK_SEARCH_OVERSHOOT {
        // Math.round of a positive value is floor(x + 0.5), as f64::round for positive values
        let current_value = (current * 1000.0).round() / 1000.0;
        current += step;
        step *= PEAK_SEARCH_STEP_RATIO;
        result.push((current_value, step));
    }
    result
}

/// SearchGrid.ValueScore.compareTo: by score, then value
pub fn compare_value_score(a: (f64, f64), b: (f64, f64)) -> std::cmp::Ordering {
    a.1.total_cmp(&b.1).then(a.0.total_cmp(&b.0))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Direction {
    Up,
    Down,
    Flat,
}

const FIRST_PEAK_ITEM_COUNT: usize = 3;
const FIRST_PEAK_MIN_RATIO: f64 = 1.25;
const PEAK_VS_TROUGH_MIN_RATIO: f64 = 0.03;

/// LocalMaximaFinder.maxima over the scores: the indices of the maxima
pub fn local_maxima(scores: &[f64]) -> Vec<usize> {
    let mut maxima: Vec<usize> = Vec::new();
    let Some(first_non_zero) = scores.iter().position(|&s| s != 0.0) else { return maxima };
    let mut previous_direction: Option<Direction> = None;
    let mut current_direction: Option<Direction> = None;
    let mut previous: Option<usize> = None;
    let mut first_non_zero_score = 0.0;
    let mut last_trough_score = -1.0;
    let mut non_zero_score_count = 0usize;
    for (i, &score) in scores.iter().enumerate() {
        if let Some(p) = previous {
            let previous_score = scores[p];
            current_direction = Some(if score > previous_score {
                Direction::Up
            } else if score < previous_score {
                Direction::Down
            } else {
                Direction::Flat
            });
            if current_direction == Some(Direction::Down) {
                // down after peak
                let mut add_peak = false;
                if previous_direction == Some(Direction::Up) {
                    // avoid creating a peak close to another at the start
                    if non_zero_score_count <= FIRST_PEAK_ITEM_COUNT {
                        if previous_score >= first_non_zero_score * FIRST_PEAK_MIN_RATIO {
                            add_peak = true;
                        }
                    } else {
                        add_peak = true;
                    }
                }
                if add_peak {
                    if last_trough_score > 0.0 && !maxima.is_empty() {
                        let last_peak = scores[*maxima.last().unwrap()];
                        let peak_to_trough_ratio = (last_peak - last_trough_score) / last_trough_score;
                        if peak_to_trough_ratio < PEAK_VS_TROUGH_MIN_RATIO {
                            maxima.pop();
                        }
                    }
                    maxima.push(p);
                    last_trough_score = -1.0;
                }
            } else if previous_direction == Some(Direction::Down) && matches!(current_direction, Some(Direction::Up) | Some(Direction::Flat)) {
                last_trough_score = previous_score;
            }
        }
        if score > 0.0 {
            non_zero_score_count += 1;
            if first_non_zero_score == 0.0 {
                first_non_zero_score = score;
            }
        }
        if current_direction != Some(Direction::Flat) {
            previous_direction = current_direction;
            previous = Some(i);
        }
    }
    // maxima.remove(firstNonZero)
    if let Some(pos) = maxima.iter().position(|&m| m == first_non_zero) {
        maxima.remove(pos);
    }
    maxima
}

/// CandidatePeakResult without the peak: the score and message of one level
#[derive(Clone, Debug)]
pub struct LevelResult {
    pub level: f64,
    pub step: f64,
    pub testable: usize,
    pub homozygous: usize,
    pub heterozygous: usize,
    pub captured: usize,
    pub score: f64,
    pub message: String,
}

/// the memoised tests of one level
struct LevelMemo {
    level: f64,
    step: f64,
    depth: Vec<u8>,
    depth_large: HashMap<i32, bool>,
    bands: Vec<u8>,
    bands_large: HashMap<(i32, i32), Band>,
}

const DENSE_DEPTH: i32 = 1024;
const DENSE_ROW: usize = (DENSE_DEPTH as usize) / 2 + 1;

impl LevelMemo {
    fn new() -> LevelMemo {
        LevelMemo {
            level: 0.0,
            step: 0.0,
            depth: vec![0; DENSE_DEPTH as usize],
            depth_large: HashMap::new(),
            bands: vec![0; DENSE_DEPTH as usize * DENSE_ROW],
            bands_large: HashMap::new(),
        }
    }
    fn reset(&mut self, level: f64, step: f64) {
        self.level = level;
        self.step = step;
        self.depth.fill(0);
        self.depth_large.clear();
        self.bands.fill(0);
        self.bands_large.clear();
    }
    fn sufficient_depth(&mut self, n: i32) -> bool {
        if (0..DENSE_DEPTH).contains(&n) {
            let slot = &mut self.depth[n as usize];
            if *slot == 0 {
                *slot = if sufficient_depth(self.level, n) { 2 } else { 1 };
            }
            *slot == 2
        } else {
            let level = self.level;
            *self.depth_large.entry(n).or_insert_with(|| sufficient_depth(level, n))
        }
    }
    fn band(&mut self, e: &PositionEvidence) -> Band {
        let (n, k) = depth_and_count(e);
        // k <= n / 2 always holds for counts >= 0; anything else is not memoised
        if (0..DENSE_DEPTH).contains(&n) && (0..=n / 2).contains(&k) {
            let slot = &mut self.bands[n as usize * DENSE_ROW + k as usize];
            if *slot == 0 {
                *slot = match band(self.level, self.step, n, k, e.symmetric_vaf()) {
                    Band::None => 1,
                    Band::Homozygous => 2,
                    Band::Heterozygous => 3,
                };
            }
            match *slot {
                2 => Band::Homozygous,
                3 => Band::Heterozygous,
                _ => Band::None,
            }
        } else if n >= DENSE_DEPTH && (0..=n / 2).contains(&k) {
            let (level, step) = (self.level, self.step);
            *self.bands_large.entry((n, k)).or_insert_with(|| band(level, step, n, k, e.symmetric_vaf()))
        } else {
            band(self.level, self.step, n, k, e.symmetric_vaf())
        }
    }
}

/// CandidatePeakEvaluation.call for one level, counting only (or keeping the bands, for `peak`)
fn evaluate(points: &[&PositionEvidence], keys: &EvidenceKeys, memo: &mut LevelMemo, mut peak: Option<&mut CandidatePeak>) -> LevelResult {
    let (level, step) = (memo.level, memo.step);
    let testable: Vec<u32> = (0..points.len() as u32).filter(|&i| memo.sufficient_depth(points[i as usize].ref_support.wrapping_add(points[i as usize].alt_support))).collect();
    let mut result = LevelResult { level, step, testable: testable.len(), homozygous: 0, heterozygous: 0, captured: 0, score: 0.0, message: String::new() };
    if testable.len() < MINIMUM_CAPTURED_POINTS {
        // CandidatePeakResult.tooFewPointsOfSufficientDepth(peak, peak.numberOfCapturedEvidencePoints()), with no
        // point tested yet
        result.message = "Too few points of sufficient depth to test contamination level: 0".to_string();
        return result;
    }
    let keep = peak.is_some() || keys.has_duplicates();
    let mut local = CandidatePeak::with_step(level, step);
    for &i in &testable {
        let e = points[i as usize];
        match memo.band(e) {
            Band::Homozygous => {
                result.homozygous += 1;
                if keep {
                    local.homozygous.push(i);
                }
            }
            Band::Heterozygous => {
                result.heterozygous += 1;
                if keep {
                    local.heterozygous.push(i);
                }
            }
            Band::None => {}
        }
    }
    local.tested = testable.len();
    result.captured = if keep { local.number_of_captured_evidence_points(keys) } else { result.homozygous + result.heterozygous };
    if result.captured < MINIMUM_CAPTURED_POINTS {
        result.message = format!("Too few points in bands: {}", result.captured);
    } else {
        result.score = result.captured as f64 / testable.len() as f64;
    }
    if let Some(p) = peak.as_deref_mut() {
        *p = local;
    }
    result
}

/// PeakSearch: every level of the grid evaluated (on `threads` threads), then the local maxima of the
/// scores, with their bands
pub struct PeakSearch {
    pub levels: Vec<LevelResult>,
    pub peaks: Vec<CandidatePeak>,
}

impl PeakSearch {
    pub fn new(points: &[&PositionEvidence], keys: &EvidenceKeys, threads: usize) -> PeakSearch {
        let grid = search_values_and_steps();
        let threads = threads.max(1).min(grid.len());
        let next = std::sync::atomic::AtomicUsize::new(0);
        let mut levels: Vec<Option<LevelResult>> = vec![None; grid.len()];
        std::thread::scope(|scope| {
            let workers: Vec<_> = (0..threads)
                .map(|_| {
                    let (grid, next) = (&grid, &next);
                    scope.spawn(move || {
                        let mut memo = LevelMemo::new();
                        let mut out = Vec::new();
                        loop {
                            let g = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if g >= grid.len() {
                                break;
                            }
                            memo.reset(grid[g].0, grid[g].1);
                            out.push((g, evaluate(points, keys, &mut memo, None)));
                        }
                        out
                    })
                })
                .collect();
            for w in workers {
                for (g, r) in w.join().unwrap() {
                    levels[g] = Some(r);
                }
            }
        });
        let levels: Vec<LevelResult> = levels.into_iter().map(|r| r.unwrap()).collect();
        let scores: Vec<f64> = levels.iter().map(|r| r.score).collect();
        let maxima = local_maxima(&scores);
        let mut memo = LevelMemo::new();
        let mut peaks = Vec::new();
        for g in maxima {
            memo.reset(grid[g].0, grid[g].1);
            let mut peak = CandidatePeak::with_step(grid[g].0, grid[g].1);
            let r = evaluate(points, keys, &mut memo, Some(&mut peak));
            crate::debug!("actual peak({}) with score({})", hmf_common::fmt::java_fixed(peak.vaf(), 3), hmf_common::fmt::java_fixed(r.score, 3));
            peak.stats = Some(PeakStats { score: r.score, ..Default::default() });
            peaks.push(peak);
        }
        PeakSearch { levels, peaks }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::AmberBase;

    /// PurityTestBase.evidenceWithDepthAndAltCount (A>T)
    pub(crate) fn evidence(chromosome: &str, position: i32, depth: i32, alt: i32) -> PositionEvidence {
        let mut e = PositionEvidence::new(chromosome, position, AmberBase::A, AmberBase::T);
        e.read_depth = depth;
        e.ref_support = depth - alt;
        e.alt_support = alt;
        e
    }

    fn evidence_with_depth(depth: i32) -> PositionEvidence {
        evidence("1", 1234, depth, 0)
    }

    fn keys_of(points: &[PositionEvidence]) -> EvidenceKeys {
        EvidenceKeys::new(&points.iter().collect::<Vec<_>>())
    }

    // CandidatePeakTest
    #[test]
    fn should_go_in_bucket_to_which_it_is_closest() {
        let mut peak = CandidatePeak::with_step(0.006, 0.01);
        peak.test(0, &evidence("chr1", 1000, 1746, 7));
        assert_eq!((peak.heterozygous_points().len(), peak.homozygous_points().len()), (1, 0));
        let mut peak = CandidatePeak::with_step(0.006, 0.01);
        peak.test(0, &evidence("chr1", 1000, 1746, 8));
        assert_eq!((peak.heterozygous_points().len(), peak.homozygous_points().len()), (0, 1));
    }

    #[test]
    fn has_sufficient_depth_for_event_detection() {
        let level35 = CandidatePeak::new(0.35);
        for (depth, expected) in [(0, false), (1, false), (20, false), (25, false), (26, true), (30, true)] {
            assert_eq!(level35.has_sufficient_depth_for_event_detection(&evidence_with_depth(depth)), expected, "depth {depth}");
        }
        let level10 = CandidatePeak::new(0.1);
        for (depth, expected) in [(0, false), (1, false), (91, false), (92, true), (1000, true)] {
            assert_eq!(level10.has_sufficient_depth_for_event_detection(&evidence_with_depth(depth)), expected, "depth {depth}");
        }
    }

    fn captured(level: f64, step: Option<f64>, depth: i32, alt: i32) -> Band {
        let mut peak = match step {
            Some(s) => CandidatePeak::with_step(level, s),
            None => CandidatePeak::new(level),
        };
        peak.test(0, &evidence("chr1", 1234, depth, alt));
        match (peak.homozygous_points().len(), peak.heterozygous_points().len()) {
            (1, 0) => Band::Homozygous,
            (0, 1) => Band::Heterozygous,
            (0, 0) => Band::None,
            _ => unreachable!(),
        }
    }

    #[test]
    fn capture_by_binomial() {
        let (none, hom, het) = (Band::None, Band::Homozygous, Band::Heterozygous);
        let cases = [
            (0.1, 1000, 20, none),
            (0.1, 1000, 42, none),
            (0.1, 1000, 50, het),
            (0.1, 1000, 56, het),
            (0.1, 1000, 57, none),
            (0.1, 1000, 79, none),
            (0.1, 1000, 90, none),
            (0.1, 1000, 91, hom),
            (0.1, 1000, 100, hom),
            (0.1, 1000, 108, hom),
            (0.1, 1000, 110, none),
            (0.30, 100, 5, none),
            (0.30, 100, 10, none),
            (0.30, 100, 11, het),
            (0.30, 100, 15, het),
            (0.30, 100, 18, het),
            (0.30, 100, 19, none),
            (0.30, 100, 20, none),
            (0.30, 100, 24, none),
            (0.30, 100, 25, hom),
            (0.30, 100, 30, hom),
            (0.30, 100, 34, hom),
            (0.30, 100, 35, none),
        ];
        for (level, depth, alt, expected) in cases {
            assert_eq!(captured(level, None, depth, alt), expected, "level {level} depth {depth} alt {alt}");
        }
    }

    #[test]
    fn capture_by_next_level_step() {
        let (none, hom, het) = (Band::None, Band::Homozygous, Band::Heterozygous);
        let cases = [
            (0.018, 40, none),
            (0.019, 40, none),
            (0.020, 40, none),
            (0.020001, 40, het),
            (0.018, 60, none),
            (0.019, 60, none),
            (0.019999, 60, none),
            (0.02, 60, het),
            (0.02, 79, none),
            (0.019, 80, none),
            (0.02, 80, none),
            (0.020001, 80, hom),
            (0.02, 120, hom),
            (0.02, 121, none),
        ];
        for (gap, alt, expected) in cases {
            assert_eq!(captured(0.1, Some(gap), 1000, alt), expected, "gap {gap} alt {alt}");
        }
        // when the VAF is close to 1
        let cases = [
            (0.018, 960, none),
            (0.019, 960, none),
            (0.020, 960, none),
            (0.020001, 960, het),
            (0.018, 940, none),
            (0.019, 940, none),
            (0.019999, 940, none),
            (0.02, 940, het),
            (0.02, 921, none),
            (0.019, 920, none),
            (0.02, 920, none),
            (0.020001, 920, hom),
            (0.02, 880, hom),
            (0.02, 879, none),
        ];
        for (gap, alt, expected) in cases {
            assert_eq!(captured(0.1, Some(gap), 1000, alt), expected, "gap {gap} alt {alt}");
        }
    }

    #[test]
    fn alt_and_ref_depth_are_used() {
        let mut pe = PositionEvidence::new("1", 1000, AmberBase::A, AmberBase::T);
        pe.alt_support = 10;
        pe.ref_support = 40;
        pe.read_depth = 10000;
        let mut peak = CandidatePeak::new(0.2);
        peak.test(0, &pe);
        assert_eq!(peak.all_captured_points(&keys_of(&[pe.clone()])), vec![0]);
        assert!(!CandidatePeak::new(0.02).has_sufficient_depth_for_event_detection(&pe));
    }

    #[test]
    fn number_of_events_captured() {
        let points = [evidence("chr3", 1_000_000, 1000, 100), evidence("chr3", 2_000_000, 1000, 50), evidence("chr3", 3_000_000, 1000, 10)];
        let keys = keys_of(&points);
        let mut level = CandidatePeak::new(0.1);
        assert_eq!(level.number_of_captured_evidence_points(&keys), 0);
        level.test(0, &points[0]);
        assert_eq!(level.number_of_captured_evidence_points(&keys), 1);
        level.test(1, &points[1]);
        assert_eq!(level.number_of_captured_evidence_points(&keys), 2);
        level.test(2, &points[2]);
        assert_eq!(level.number_of_captured_evidence_points(&keys), 2);
        assert_eq!(level.vaf(), 0.1);
    }

    #[test]
    fn homozygous_proportion() {
        let points = [
            evidence("chr3", 3_000_000, 1000, 10),
            evidence("chr3", 4_000_000, 1000, 50),
            evidence("chr3", 5_000_000, 1000, 50),
            evidence("chr3", 6_000_000, 1000, 100),
            evidence("chr3", 7_000_000, 1000, 100),
        ];
        let keys = keys_of(&points);
        let mut level = CandidatePeak::new(0.1);
        assert!(level.homozygous_proportion(&keys).is_nan());
        let expected = [f64::NAN, 0.0, 0.0, 1.0 / 3.0, 0.5];
        for (i, e) in points.iter().enumerate() {
            level.test(i as u32, e);
            let got = level.homozygous_proportion(&keys);
            assert!(got == expected[i] || (got.is_nan() && expected[i].is_nan()), "{i}: {got}");
        }
    }

    #[test]
    fn duplicates_count_once() {
        // two equal points (same site) are one point of the set; the band lists keep both
        let points = [evidence("chr3", 1_000_000, 1000, 100), evidence("chr3", 1_000_000, 1000, 100), evidence("chr3", 2_000_000, 1000, 50)];
        let keys = keys_of(&points);
        assert!(keys.has_duplicates());
        let mut level = CandidatePeak::new(0.1);
        for (i, e) in points.iter().enumerate() {
            level.test(i as u32, e);
        }
        assert_eq!(level.homozygous_points().len(), 2);
        assert_eq!(level.number_of_captured_evidence_points(&keys), 2);
        assert_eq!(level.homozygous_proportion(&keys), 1.0);
        // out of order, found through the map
        let shuffled = [points[2].clone(), points[0].clone(), points[1].clone()];
        assert!(keys_of(&shuffled).has_duplicates());
        assert!(!keys_of(&[points[0].clone(), points[2].clone()]).has_duplicates());
    }

    #[test]
    fn to_string() {
        let mut level = CandidatePeak::with_step(0.1, 0.02);
        level.test(0, &evidence("chr3", 1_000_000, 1000, 100));
        level.test(1, &evidence("chr3", 2_000_000, 1000, 50));
        level.test(2, &evidence("chr3", 2_001_000, 1000, 50));
        level.test(3, &evidence("chr3", 3_000_000, 1000, 10));
        assert_eq!(level.to_string(), "VafLevel{vaf=0.10, step=0.02, tested: 4, homozygous: 1, heterozygous: 2}");
    }

    // RangeStepPeakTest
    #[test]
    fn range_step_peak() {
        let (ref_support, alt_support) = (1739, 7);
        let depth = ref_support + alt_support;
        let symmetric_vaf = alt_support as f64 / depth as f64;
        assert!(is_captured(depth, alt_support, 0.006, 0.001103, symmetric_vaf));
        assert!(!is_captured(depth, alt_support, 0.003, 0.001103 / 2.0, symmetric_vaf));
    }

    // CandidatePeakEvaluationTest
    struct Points {
        position: i32,
    }
    impl Points {
        fn make(&mut self, n: usize, depth: i32, ref_depth: i32, alt_depth: i32) -> Vec<PositionEvidence> {
            (0..n)
                .map(|_| {
                    let mut e = PositionEvidence::new("1", self.position, AmberBase::A, AmberBase::T);
                    e.read_depth = depth;
                    e.ref_support = ref_depth;
                    e.alt_support = alt_depth;
                    self.position += 100;
                    e
                })
                .collect()
        }
        fn hom(&mut self, n: usize, depth: i32) -> Vec<PositionEvidence> {
            let alt = (0.1 * depth as f64) as i32;
            self.make(n, depth, depth - alt, alt)
        }
        fn het(&mut self, n: usize, depth: i32) -> Vec<PositionEvidence> {
            let alt = (0.1 / 2.0 * depth as f64) as i32;
            self.make(n, depth, depth - alt, alt)
        }
        fn uncaptured(&mut self, n: usize, depth: i32) -> Vec<PositionEvidence> {
            self.make(n, depth, depth, 0)
        }
    }

    fn result_of(points: &[PositionEvidence]) -> LevelResult {
        let refs: Vec<&PositionEvidence> = points.iter().collect();
        let keys = EvidenceKeys::new(&refs);
        let mut memo = LevelMemo::new();
        memo.reset(0.1, 0.00001);
        evaluate(&refs, &keys, &mut memo, None)
    }

    #[test]
    fn need_to_capture_minimum_number_of_deep_points() {
        let mut p = Points { position: 1000 };
        let deep = p.hom(10, 1000);
        let shallow = p.hom(5, 91);
        let all: Vec<PositionEvidence> = deep.iter().chain(shallow.iter()).cloned().collect();
        let r = result_of(&all);
        assert_eq!(r.score, 0.0);
        assert_eq!(r.message, "Too few points of sufficient depth to test contamination level: 0");
        let more = p.hom(5, 920);
        let all: Vec<PositionEvidence> = deep.iter().chain(more.iter()).cloned().collect();
        assert!(result_of(&all).score > 0.0);
    }

    #[test]
    fn need_to_capture_a_minimum_number_of_points() {
        let mut p = Points { position: 1000 };
        let captured = p.hom(14, 1000);
        let missed = p.uncaptured(14, 1000);
        let all: Vec<PositionEvidence> = captured.iter().chain(missed.iter()).cloned().collect();
        let r = result_of(&all);
        assert_eq!(r.score, 0.0);
        assert_eq!(r.message, "Too few points in bands: 14");
        let more = p.hom(1, 1000);
        let all: Vec<PositionEvidence> = captured.iter().chain(more.iter()).cloned().collect();
        assert!(result_of(&all).score > 0.0);
    }

    #[test]
    fn captured_points_can_be_in_het_peak() {
        let mut p = Points { position: 1000 };
        assert!(result_of(&p.het(15, 1000)).score > 0.0);
    }

    #[test]
    fn score_is_ratio_of_captured_points() {
        let mut p = Points { position: 1000 };
        let mut all = p.hom(15, 1000);
        all.extend(p.het(16, 1000));
        all.extend(p.uncaptured(69, 1000));
        assert!((result_of(&all).score - 0.31).abs() < 1e-4);
    }

    // SearchGridTest
    #[test]
    fn search_values() {
        let values = search_values_and_steps();
        let d = 0.00001;
        assert!((values[0].0 - 0.005).abs() < d && (values[0].1 - 0.00105).abs() < d);
        assert!((values[1].0 - 0.006).abs() < d && (values[1].1 - 0.0011).abs() < d);
        assert!((values[2].0 - 0.007).abs() < d && (values[2].1 - 0.00116).abs() < d);
        let n = values.len();
        assert!((values[n - 1].0 - 0.359).abs() < d);
        assert!((values[n - 2].0 - 0.341).abs() < d);
        assert!((values[n - 1].1 / values[n - 2].1 - 1.05).abs() < d);
        assert!(n >= 60);
        for w in values.windows(2) {
            assert!(w[1].0 - w[0].0 >= 0.001 - 1e-12);
        }
    }

    #[test]
    fn value_score_comparison() {
        use std::cmp::Ordering::*;
        let vs1 = (0.05, 100.0);
        assert_eq!(compare_value_score(vs1, (0.05, 100.0)), Equal);
        assert_eq!(compare_value_score(vs1, (0.05, 101.0)), Less);
        assert_eq!(compare_value_score(vs1, (0.04, 100.0)), Greater);
        assert_eq!(compare_value_score(vs1, (0.04, 101.0)), Less);
        assert_eq!(compare_value_score(vs1, (0.04, 99.0)), Greater);
        assert_eq!(compare_value_score(vs1, (0.06, 99.0)), Greater);
    }

    // LocalMaximaFinderTest
    #[test]
    fn maxima() {
        assert!(local_maxima(&[1.0, 2.0, 3.0]).is_empty());
        assert!(local_maxima(&[10.0, 5.0, 3.0]).is_empty());
        assert!(local_maxima(&[10.0, 10.0, 10.0]).is_empty());
        assert!(local_maxima(&[]).is_empty());
        assert!(local_maxima(&[10.0]).is_empty());
        assert_eq!(local_maxima(&[1.0, 10.0, 10.0, 5.0]), vec![1]);
        assert!(local_maxima(&[100.0, 10.0, 10.0, 50.0]).is_empty());
        assert!(local_maxima(&[0.0, 0.0, 0.001, 0.0]).is_empty());
        assert!(local_maxima(&[0.0, 0.0, 0.002, 0.001, 0.0]).is_empty());
        assert_eq!(local_maxima(&[0.0, 0.0, 0.001, 0.002, 0.0]), vec![3]);
        let typical = [100.0, 50.0, 20.0, 18.0, 19.0, 20.0, 25.0, 30.0, 29.0, 25.0, 21.0, 22.1, 22.0, 25.0, 30.0, 28.0, 25.0, 30.0];
        assert_eq!(local_maxima(&typical), vec![7, 14]);
    }

    // PeaksSearchTest
    #[test]
    fn peak_search() {
        let mut evidence_list = Vec::new();
        let start = 2_000_000;
        let iterations = 4000;
        for i in 0..iterations {
            let position = start + i * 1000;
            let multiplier = 1.0 + (i as f64 - 0.5 * iterations as f64) / (20 * iterations) as f64;
            let alt = |base: f64| (base * multiplier) as i32;
            let e = match i % 16 {
                0 => evidence("chr3", position, 1000, 0),
                1 | 3 => evidence("chr3", position, 1000, alt(100.0)),
                2 => evidence("chr3", position, 1000, alt(200.0)),
                4 | 8 => evidence("chr3", position, 1000, alt(400.0)),
                5 | 7 | 9 | 11 => evidence("chr3", position, 1000, alt(500.0)),
                6 | 10 => evidence("chr3", position, 1000, alt(600.0)),
                12 => evidence("chr3", position, 1000, alt(800.0)),
                13 | 14 => evidence("chr3", position, 1000, alt(900.0)),
                _ => evidence("chr3", position, 1000, 1000),
            };
            evidence_list.push(e);
        }
        let run = |list: &[PositionEvidence]| {
            let refs: Vec<&PositionEvidence> = list.iter().collect();
            let keys = EvidenceKeys::new(&refs);
            PeakSearch::new(&refs, &keys, 2).peaks.iter().map(|p| p.vaf()).collect::<Vec<f64>>()
        };
        let peaks = run(&evidence_list);
        assert_eq!(peaks.len(), 2);
        assert!((peaks[0] - 0.1).abs() < 0.01 && (peaks[1] - 0.2).abs() < 0.01);
        for i in 0..100 {
            evidence_list.push(evidence("chr3", 1_000_000 + i * 1000, 1000, 300));
        }
        let peaks = run(&evidence_list);
        assert_eq!(peaks.len(), 3);
        assert!((peaks[0] - 0.1).abs() < 0.01 && (peaks[1] - 0.2).abs() < 0.01 && (peaks[2] - 0.3).abs() < 0.01);
    }
}
