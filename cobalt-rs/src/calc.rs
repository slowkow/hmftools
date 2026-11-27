//! CobaltCalculator and its parts: BamCalculation (TumorCalculation, ReferenceCalculation), CobaltWindow,
//! BamRatio, GCPail and GcBucketStatistics, the normalisers (ReadDepthStatisticsNormaliser,
//! DiploidNormaliser with DiploidRatioNormaliser, UnityNormaliser), the consolidators (LowCoverageConsolidator)
//! and ResultsCollator.

use anyhow::{bail, Result};

use crate::chromosome::{Chr, PerChr};
use crate::compat::java;
use crate::constants::*;
use crate::count::SampleCounts;
use crate::gc::WindowStatuses;
use crate::stats;
use crate::targeted::Scope;

/// BamRatio: a window's read depth, GC content, ratio and diploid-adjusted ratio
#[derive(Clone, Debug, PartialEq)]
pub struct BamRatio {
    pub position: i32,
    pub read_depth: f64,
    pub gc_content: f64,
    ratio: f64,
    diploid: f64,
    included: bool,
}

impl BamRatio {
    /// BamRatio(chromosome, position, readDepth, gcContent, included)
    pub fn new(position: i32, read_depth: f64, gc_content: f64, included: bool) -> BamRatio {
        let read_depth = if read_depth.is_finite() { read_depth } else { -1.0 };
        let gc_content = if gc_content.is_finite() { gc_content } else { -1.0 };
        BamRatio { position, read_depth, gc_content, ratio: if included { read_depth } else { -1.0 }, diploid: -1.0, included }
    }

    /// BamRatio(chromosome, position, readDepth, gcContent): the ratio is the read depth
    pub fn simple(position: i32, read_depth: f64, gc_content: f64) -> BamRatio {
        BamRatio { position, read_depth, gc_content, ratio: read_depth, diploid: -1.0, included: true }
    }

    pub fn is_valid(&self) -> bool {
        self.ratio >= 0.0
    }
    pub fn ratio(&self) -> f64 {
        if self.included { self.ratio } else { -1.0 }
    }
    pub fn diploid_adjusted_ratio(&self) -> f64 {
        self.diploid
    }

    pub fn normalise(&mut self, factor: f64) {
        if factor.is_nan() {
            self.included = false;
            self.ratio = -1.0;
            return;
        }
        // CO4: Java returns here for a ratio of 0 whatever the factor, so a window without reads keeps ratio 0
        // where the factor masks every other window (a GC bucket outside the limits or without depth).
        // Corrected: a ratio of 0 is kept only when the factor applies.
        if doubles::is_zero(self.ratio) && (java("CO4") || (factor > 0.0 && self.included)) {
            return;
        }
        if factor <= 0.0 || !self.included || self.ratio < 0.0 {
            self.included = false;
            self.ratio = -1.0;
        } else {
            self.ratio /= factor;
        }
    }

    pub fn normalise_diploid_adjusted_ratio(&mut self, factor: f64) {
        if self.diploid == 0.0 {
            return;
        }
        if factor <= 0.0 || factor.is_nan() || self.diploid < 0.0 {
            self.diploid = -1.0;
        } else {
            self.diploid /= factor;
        }
    }

    pub fn set_diploid_adjusted_ratio(&mut self, ratio: f64) {
        self.diploid = if self.ratio == 0.0 { 0.0 } else { ratio };
    }

    pub fn override_ratio(&mut self, ratio: f64) {
        self.ratio = ratio;
        if ratio > 0.0 {
            self.included = true;
        }
    }
}

/// The ratios of one chromosome, in position order
#[derive(Clone, Debug)]
pub struct ChrRatios {
    pub chr: Chr,
    pub ratios: Vec<BamRatio>,
}

/// MedianRatio
#[derive(Clone, Debug, PartialEq)]
pub struct MedianRatio {
    pub chromosome: String,
    pub median_ratio: f64,
    pub count: i32,
}

/// GcMedianReadDepth
#[derive(Clone, Debug)]
pub struct GcMedianReadDepth {
    pub mean: f64,
    pub median: f64,
    /// the bucket medians above 0 (bucketToMedianReadDepth)
    pub bucket_medians: Vec<(usize, f64)>,
}

impl GcMedianReadDepth {
    pub fn no_results() -> Self {
        GcMedianReadDepth { mean: 0.0, median: 0.0, bucket_medians: Vec::new() }
    }
}

/// CobaltWindow: a window's read depth and GC content, whether it is masked (excluded) or on target, and
/// its GC bucket once bucketed (none for a masked window)
#[derive(Clone, Debug, PartialEq)]
pub struct CobaltWindow {
    pub position: i32,
    pub read_depth: f64,
    pub gc_content: f64,
    pub bucket: Option<usize>,
    pub excluded: bool,
    pub in_target: bool,
}

impl CobaltWindow {
    pub fn new(position: i32, read_depth: f64, gc_content: f64, excluded: bool, in_target: bool) -> CobaltWindow {
        CobaltWindow { position, read_depth, gc_content, bucket: None, excluded, in_target }
    }
    pub fn include(&self) -> bool {
        !self.excluded && self.in_target
    }
    /// correctedByReferenceValue: an included window of depth below 1 takes the GC profile's GC content
    /// (unchanged where there is no profile window, where Java catches the exception)
    pub fn corrected_by_reference_value(self, reference_gc: impl FnOnce() -> Option<f64>) -> CobaltWindow {
        if self.read_depth >= 1.0 || self.read_depth < 0.0 || !self.include() {
            return self;
        }
        match reference_gc() {
            Some(gc_content) => CobaltWindow { gc_content, excluded: false, ..self },
            None => self,
        }
    }
    /// GCPailsList.getGCPail, then bucketed: the GC bucket of a window that is not masked, and its depth in
    /// the bucket for an autosome window on target with reads. Java's list access throws outside 0-100.
    pub fn bucketed(self, chr: Chr, pails: &mut [Vec<f64>; 101]) -> Result<CobaltWindow> {
        let index = gc_bucket_index(self.gc_content);
        if !(0..=100).contains(&index) {
            bail!("java.lang.IndexOutOfBoundsException: GC bucket {index}");
        }
        if self.excluded {
            return Ok(self);
        }
        if chr.is_autosome() && self.in_target && self.read_depth > 0.0 {
            pails[index as usize].push(self.read_depth);
        }
        Ok(CobaltWindow { bucket: Some(index as usize), ..self })
    }
    pub fn to_bam_ratio(&self) -> BamRatio {
        BamRatio::new(self.position, self.read_depth, self.gc_content, self.include())
    }
}

/// GCPail.bucketIndex: (int) Math.round(gcContent * 100)
pub fn gc_bucket_index(gc_content: f64) -> i64 {
    hmf_common::java_math::java_round(gc_content * 100.0)
}

/// GcBucketStatistics: for each allowed bucket the mean of the medians of it and its two neighbours
pub struct GcBucketStatistics {
    min_allowed: i32,
    max_allowed: i32,
    mean_depths: [f64; 101],
}

impl GcBucketStatistics {
    pub fn new(pails: &mut [Vec<f64>; 101], min_allowed: i32, max_allowed: i32) -> Result<GcBucketStatistics> {
        if min_allowed <= 0 || max_allowed >= 100 {
            bail!("java.lang.IllegalArgumentException (GcBucketStatistics: GC buckets {min_allowed} to {max_allowed})");
        }
        let medians: Vec<f64> = pails.iter_mut().map(|p| if p.is_empty() { 0.0 } else { stats::median(p) }).collect();
        let mut s = GcBucketStatistics { min_allowed, max_allowed, mean_depths: [0.0; 101] };
        let mut window = [-1.0f64; 3];
        for i in 0..100 {
            window[0] = window[1];
            window[1] = window[2];
            window[2] = medians[i + 1];
            // DoubleStream.average: the compensated sum over the count
            s.mean_depths[i] = if s.is_allowed(Some(i + 1)) { hmf_common::java_math::double_stream_sum(window) / 3.0 } else { -1.0 };
        }
        s.mean_depths[100] = -1.0;
        Ok(s)
    }

    fn is_allowed(&self, bucket: Option<usize>) -> bool {
        bucket.is_some_and(|b| b as i32 > self.min_allowed && b as i32 <= self.max_allowed)
    }

    pub fn median_read_depth(&self, bucket: Option<usize>) -> f64 {
        if !self.is_allowed(bucket) { -1.0 } else { self.mean_depths[bucket.unwrap()] }
    }

    pub fn bucket_medians(&self) -> Vec<(usize, f64)> {
        (0..101).filter(|&i| self.mean_depths[i] > 0.0).map(|i| (i, self.mean_depths[i])).collect()
    }
}

/// ResultsConsolidator: none, or LowCoverageConsolidator with its buckets, made from the first ratios it
/// consolidates (the tumor's) and kept for the reference.
pub enum Consolidator {
    NoOp,
    LowCoverage { count: i32, boundaries: Vec<(Chr, Vec<LowCovBucket>)> },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LowCovBucket {
    pub start: i32,
    pub end: i32,
    pub position: i32,
}

/// ResultsConsolidator.calcConsolidationCount: 80 / median read depth, 1 below 10, at most 1000, else
/// rounded to one significant digit
pub fn calc_consolidation_count(median_read_depth: f64) -> i32 {
    let c = 80.0 / median_read_depth;
    if c < 10.0 {
        return 1;
    }
    if c >= 1000.0 {
        return 1000;
    }
    let round_by = 10f64.powf(c.log10().floor());
    (hmf_common::java_math::java_round(c / round_by) as f64 * round_by) as i32
}

fn round_down_to_window_boundary(p: f64) -> i32 {
    ((p / WINDOW_SIZE as f64).floor() * WINDOW_SIZE as f64) as i32 + 1
}

/// ResultsConsolidator.consolidateIntoBuckets
pub fn consolidate_into_buckets(positions: &[i32], count: i32) -> Vec<LowCovBucket> {
    let mut buckets = Vec::new();
    if positions.is_empty() {
        return buckets;
    }
    let mut window_count = 0;
    let mut bucket_start = positions[0];
    for i in 0..positions.len() {
        let position = positions[i];
        if position - bucket_start >= MAX_SPARSE_CONSOLIDATE_DISTANCE {
            if i > 0 {
                let bucket_end = positions[i - 1] + WINDOW_SIZE;
                let bucket_pos = round_down_to_window_boundary((bucket_start + bucket_end) as f64 * 0.5);
                buckets.push(LowCovBucket { start: bucket_start, end: bucket_end, position: bucket_pos });
            }
            bucket_start = position;
            window_count = 0;
        }
        if window_count == count {
            let bucket_end = round_down_to_window_boundary((positions[i - 1] + position) as f64 * 0.5);
            let bucket_pos = round_down_to_window_boundary((bucket_start + bucket_end) as f64 * 0.5);
            buckets.push(LowCovBucket { start: bucket_start, end: bucket_end, position: bucket_pos });
            bucket_start = bucket_end + WINDOW_SIZE;
            window_count = 0;
        }
        window_count += 1;
    }
    if window_count > 0 {
        let bucket_end = positions[positions.len() - 1] + WINDOW_SIZE;
        let bucket_pos = round_down_to_window_boundary((bucket_start + bucket_end) as f64 * 0.5);
        buckets.push(LowCovBucket { start: bucket_start, end: bucket_end, position: bucket_pos });
    }
    buckets
}

impl Consolidator {
    #[cfg(test)]
    pub fn consolidated_for_test(&mut self, ratios: &[ChrRatios]) -> Vec<(Chr, Vec<BamRatio>)> {
        self.consolidate(ratios).unwrap_or_default()
    }

    /// WholeGenome.resultsConsolidator
    pub fn for_whole_genome(median_read_depth: f64) -> Consolidator {
        if median_read_depth.is_nan() {
            return Consolidator::NoOp;
        }
        let count = calc_consolidation_count(median_read_depth);
        if count == 1 {
            crate::info!("median read depth: {}, not using sparse consolidation", java_double(median_read_depth));
            return Consolidator::NoOp;
        }
        crate::info!("median read depth: {}, sparse consolidation count: {count}", java_double(median_read_depth));
        Consolidator::LowCoverage { count, boundaries: Vec::new() }
    }

    /// consolidate: the consolidated ratios of each chromosome with buckets (LowCoverageConsolidator), or
    /// None for no change
    fn consolidate(&mut self, ratios: &[ChrRatios]) -> Option<Vec<(Chr, Vec<BamRatio>)>> {
        let Consolidator::LowCoverage { count, boundaries } = self else { return None };
        if boundaries.is_empty() {
            for c in ratios {
                let positions: Vec<i32> = c.ratios.iter().filter(|r| r.ratio() >= 0.0).map(|r| r.position).collect();
                let buckets = consolidate_into_buckets(&positions, *count);
                crate::info!("chromosome: {}, low cov buckets count: {}", c.chr.short_name(), buckets.len());
                if !buckets.is_empty() {
                    boundaries.push((c.chr, buckets));
                }
            }
        }
        let empty = Vec::new();
        let mut result = Vec::new();
        for (chr, buckets) in boundaries.iter() {
            let list = ratios.iter().find(|c| c.chr == *chr).map(|c| &c.ratios).unwrap_or(&empty);
            let mut out = Vec::new();
            let mut it = buckets.iter();
            let mut bucket = it.next();
            let mut values: Vec<f64> = Vec::new();
            let mut gcs: Vec<f64> = Vec::new();
            for r in list {
                let Some(b) = bucket else { continue };
                if r.ratio() >= 0.0 {
                    if r.position > b.end {
                        out.push(BamRatio::simple(b.position, stats::mean(&values), stats::mean(&gcs)));
                        bucket = it.next();
                        values.clear();
                        gcs.clear();
                        if bucket.is_some() {
                            values.push(r.ratio());
                            gcs.push(r.gc_content);
                        }
                    } else {
                        values.push(r.ratio());
                        gcs.push(r.gc_content);
                    }
                }
            }
            if let Some(b) = bucket {
                out.push(BamRatio::simple(b.position, stats::mean(&values), stats::mean(&gcs)));
            }
            result.push((*chr, out));
        }
        Some(result)
    }
}

/// BamRatios.consolidate: each original window at a consolidated position takes its ratio, the others -1
pub(crate) fn apply_consolidation(ratios: &mut [ChrRatios], consolidator: &mut Consolidator) -> Result<()> {
    if ratios.is_empty() {
        return Ok(());
    }
    let Some(consolidated) = consolidator.consolidate(ratios) else { return Ok(()) };
    apply_consolidated(ratios, &consolidated)
}

/// BamRatios.consolidate with the consolidated ratios of each chromosome
pub fn apply_consolidated(ratios: &mut [ChrRatios], consolidated: &[(Chr, Vec<BamRatio>)]) -> Result<()> {
    if ratios.iter().all(|c| c.ratios.is_empty()) {
        return Ok(());
    }
    let total: usize = ratios.iter().map(|c| c.ratios.len()).sum();
    let consolidated_total: usize = consolidated.iter().map(|c| c.1.len()).sum();
    if total < consolidated_total {
        bail!("java.lang.IllegalStateException (BamRatios.consolidate: more consolidated ratios than ratios)");
    }
    let mut keys: Vec<Chr> = ratios.iter().map(|c| c.chr).collect();
    let mut ckeys: Vec<Chr> = consolidated.iter().map(|c| c.0).collect();
    keys.sort();
    ckeys.sort();
    if keys != ckeys {
        bail!("java.lang.IllegalStateException (BamRatios.consolidate: a chromosome without consolidated ratios)");
    }
    if consolidated_total == total {
        return Ok(());
    }
    for c in ratios.iter_mut() {
        let cons = &consolidated.iter().find(|x| x.0 == c.chr).unwrap().1;
        if c.ratios.len() <= cons.len() {
            bail!("java.lang.IllegalStateException (BamRatios.consolidate: chromosome {})", c.chr.short_name());
        }
        let mut originals = c.ratios.iter_mut();
        for cr in cons {
            for original in originals.by_ref() {
                if original.position == cr.position {
                    original.override_ratio(cr.ratio());
                    break;
                }
                original.override_ratio(-1.0);
            }
        }
        for original in originals {
            original.override_ratio(-1.0);
        }
    }
    Ok(())
}

/// The order in which Java sums the values of a sample's chromosomes (CO3): the chromosomes' order in its
/// ratio map. Java's map is a HashMap keyed by enum, whose order follows identity hashes and so differs
/// from run to run; `java_order` is the order under -XX:hashCode=2. Corrected: genome order.
fn summation_order(ratios: &[ChrRatios], java_order: &[Chr]) -> Vec<usize> {
    if java("CO3") {
        java_order.iter().filter_map(|c| ratios.iter().position(|r| r.chr == *c)).collect()
    } else {
        let mut idx: Vec<usize> = (0..ratios.len()).collect();
        idx.sort_by_key(|&i| ratios[i].chr);
        idx
    }
}

/// ReadDepthStatisticsNormaliser: the mean and median read depth of autosome windows of valid ratio within
/// the GC range
struct ReadDepthStatistics {
    mean: f64,
    median: f64,
}

fn read_depth_statistics(ratios: &[ChrRatios], order: &[usize], gc_min: f64, gc_max: f64) -> ReadDepthStatistics {
    let mut depths: Vec<f64> = Vec::new();
    for &i in order {
        let c = &ratios[i];
        if !c.chr.is_autosome() {
            continue;
        }
        for r in &c.ratios {
            if r.gc_content < gc_min || r.gc_content > gc_max {
                continue;
            }
            if r.ratio() > 0.0 {
                depths.push(r.read_depth);
            }
        }
    }
    let mean = stats::mean(&depths);
    let median = stats::median(&mut depths);
    ReadDepthStatistics { mean, median }
}

#[cfg(test)]
pub fn read_depth_statistics_pub(ratios: &[ChrRatios], gc_min: f64, gc_max: f64) -> (f64, f64) {
    let order: Vec<usize> = (0..ratios.len()).collect();
    let s = read_depth_statistics(ratios, &order, gc_min, gc_max);
    (s.mean, s.median)
}

#[cfg(test)]
pub fn unity_normaliser_pub(ratios: &mut [ChrRatios]) {
    let order: Vec<usize> = (0..ratios.len()).collect();
    unity_normaliser(ratios, &order)
}

/// What differs between the tumor and the reference calculation
pub enum Sample<'a> {
    Tumor,
    /// with the tumor's consolidator; `diploid_normalise` (whole genome) gives the median ratios
    Reference { consolidator: &'a mut Consolidator, diploid_normalise: bool, v38: bool },
}

/// The result of BamCalculation.calculateRatios
pub struct Calculation {
    pub ratios: Vec<ChrRatios>,
    pub gc_median: GcMedianReadDepth,
    /// the tumor's consolidator (TumorCalculation.consolidator)
    pub consolidator: Option<Consolidator>,
    /// DiploidNormaliser.medianRatios (reference, whole genome)
    pub median_ratios: Option<Vec<MedianRatio>>,
}

pub struct CalcConfig {
    pub gc_bucket_min: i32,
    pub gc_bucket_max: i32,
    pub gc_ratio_min: f64,
    pub gc_ratio_max: f64,
}

/// BamCalculation: addReading for every window, then calculateRatios
pub fn calculate(counts: &SampleCounts, statuses: &WindowStatuses, scope: &Scope, cfg: &CalcConfig, sample: Sample) -> Result<Calculation> {
    // the windows (addReading), and the GC pails; a chromosome without a whole window is not in Java's maps
    let chromosomes: Vec<&crate::count::ChrCounts> = counts.chromosomes.iter().filter(|c| !c.bases.is_empty()).collect();
    let mut pails: [Vec<f64>; 101] = std::array::from_fn(|_| Vec::new());
    let mut windows: Vec<Vec<CobaltWindow>> = Vec::with_capacity(chromosomes.len());
    for c in &chromosomes {
        let mut list = Vec::with_capacity(c.bases.len());
        for i in 0..c.bases.len() {
            let position = i as i32 * WINDOW_SIZE + 1;
            let bases = c.bases[i] as f64;
            let read_depth = bases / WINDOW_SIZE as f64;
            // DepthReading: GC content 0 for no depth
            let gc_content = if read_depth == 0.0 { 0.0 } else { c.gc[i] as f64 / bases };
            let excluded = statuses.exclude(c.chr, position)?;
            let in_target = scope.on_target(c.chr, position)?;
            let window = CobaltWindow::new(position, read_depth, gc_content, excluded, in_target)
                .corrected_by_reference_value(|| statuses.reference_gc(c.chr, position))
                .bucketed(c.chr, &mut pails)
                .map_err(|e| anyhow::anyhow!("{e} (window {}:{position})", c.contig))?;
            list.push(window);
        }
        windows.push(list);
    }
    let bucket_stats = GcBucketStatistics::new(&mut pails, cfg.gc_bucket_min, cfg.gc_bucket_max)?;
    drop(pails);
    let mut ratios: Vec<ChrRatios> = Vec::with_capacity(windows.len());
    for (c, list) in chromosomes.iter().zip(windows) {
        let mut out = Vec::with_capacity(list.len());
        for w in list {
            let mut r = w.to_bam_ratio();
            r.normalise(bucket_stats.median_read_depth(w.bucket));
            r.normalise(scope.region_enrichment(c.chr, w.position)?);
            out.push(r);
        }
        ratios.push(ChrRatios { chr: c.chr, ratios: out });
    }
    let header_order: Vec<Chr> = chromosomes.iter().map(|c| c.chr).collect();
    let java_order = crate::chromosome::java_multimap_order(&header_order, 3);
    let order = summation_order(&ratios, &java_order);

    // the mean normaliser (no-op for a targeted panel, whose statistics are still taken)
    let depth_stats = read_depth_statistics(&ratios, &order, cfg.gc_ratio_min, cfg.gc_ratio_max);
    if matches!(scope, Scope::WholeGenome) && !ratios.iter().all(|c| c.ratios.is_empty()) {
        let factor = depth_stats.mean / depth_stats.median;
        if factor.is_nan() || factor <= 0.0 {
            bail!("java.lang.IllegalStateException (ReadDepthStatisticsNormaliser: mean/median read depth {factor}, no autosome window with a valid ratio)");
        }
        for c in ratios.iter_mut() {
            for r in c.ratios.iter_mut() {
                r.normalise(factor);
            }
        }
    }

    let mut median_ratios = None;
    let mut own_consolidator = None;
    match sample {
        Sample::Tumor => {
            let mut consolidator = match scope {
                Scope::WholeGenome => Consolidator::for_whole_genome(depth_stats.median),
                Scope::Targeted(_) => Consolidator::NoOp,
            };
            apply_consolidation(&mut ratios, &mut consolidator)?;
            own_consolidator = Some(consolidator);
        }
        Sample::Reference { consolidator, diploid_normalise, v38 } => {
            apply_consolidation(&mut ratios, consolidator)?;
            if diploid_normalise {
                median_ratios = Some(diploid_normaliser(&mut ratios, v38));
            }
        }
    }
    if let Scope::Targeted(_) = scope {
        unity_normaliser(&mut ratios, &order);
    }
    Ok(Calculation {
        ratios,
        gc_median: GcMedianReadDepth { mean: depth_stats.mean, median: depth_stats.median, bucket_medians: bucket_stats.bucket_medians() },
        consolidator: own_consolidator,
        median_ratios,
    })
}

/// UnityNormaliser: the ratios over the mean autosome ratio, the diploid-adjusted ratios over theirs
fn unity_normaliser(ratios: &mut [ChrRatios], order: &[usize]) {
    let mut values = Vec::new();
    let mut diploid = Vec::new();
    for &i in order {
        let c = &ratios[i];
        if !c.chr.is_autosome() {
            continue;
        }
        for r in &c.ratios {
            if r.ratio() > 0.0 {
                values.push(r.ratio());
            }
            if r.diploid_adjusted_ratio() > 0.0 {
                diploid.push(r.diploid_adjusted_ratio());
            }
        }
    }
    let mean = stats::mean(&values);
    let diploid_mean = stats::mean(&diploid);
    for c in ratios.iter_mut() {
        for r in c.ratios.iter_mut() {
            r.normalise(mean);
            if !diploid.is_empty() {
                r.normalise_diploid_adjusted_ratio(diploid_mean);
            }
        }
    }
}

/// Java's Math.min for doubles: NaN if either is NaN, -0.0 below 0.0
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a == 0.0 && b == 0.0 {
        if a.is_sign_negative() { a } else { b }
    } else {
        a.min(b)
    }
}

/// CobaltChromosomes: the actual (expected) ratio of each chromosome it keeps
pub fn expected_ratios(median_ratios: &[MedianRatio]) -> PerChr<Option<f64>> {
    let chr_of = |m: &MedianRatio| Chr::from_name(&m.chromosome).expect("human chromosome");
    let (mut min_autosome, mut y_median, mut x_median) = (0.0f64, 0.0f64, 0.0f64);
    for m in median_ratios {
        let c = chr_of(m);
        if c.is_autosome() {
            min_autosome = if min_autosome == 0.0 { m.median_ratio } else { java_min(min_autosome, m.median_ratio) };
        } else if c == Chr::Y && m.count >= MIN_Y_COUNT {
            y_median = m.median_ratio;
        } else if c == Chr::X {
            x_median = m.median_ratio;
        }
    }
    let female = doubles::greater_or_equal(x_median, TWO_X_CUTOFF) && doubles::less_than(y_median, Y_CUTOFF);
    let mut result: PerChr<Option<f64>> = PerChr::default();
    for m in median_ratios {
        let c = chr_of(m);
        if c == Chr::Y && m.count < MIN_Y_COUNT {
            continue;
        }
        let median = m.median_ratio;
        let trisomy = doubles::greater_or_equal(median, TRISOMY_CUTOFF) && matches!(c.short_name(), "13" | "15" | "18" | "21" | "X");
        let xyy = !female && c == Chr::Y && doubles::greater_or_equal(median, TWO_Y_CUTOFF);
        let mosaic_x = female && doubles::less_than(median, java_min(min_autosome, MOSAIC_X_CUTOFF)) && c == Chr::X;
        let klinefelter = !female && c == Chr::X && doubles::greater_or_equal(median, TWO_X_CUTOFF);
        let typical = if female {
            if c == Chr::Y { 0.0 } else { 1.0 }
        } else if c.is_autosome() {
            1.0
        } else {
            0.5
        };
        let actual = if trisomy {
            1.5
        } else if xyy || klinefelter {
            1.0
        } else if mosaic_x {
            median
        } else {
            typical
        };
        if doubles::positive(typical) {
            *result.get_mut(c) = Some(actual);
        }
    }
    result
}

/// DiploidNormaliser over the reference ratios: per chromosome the median of the valid ratios, the
/// expected ratio from CobaltChromosomes, then each valid ratio times the expected ratio over the rolling
/// median of the valid ratios within 5,000 of it (by index) when at least 1,000 are; Y keeps its ratio.
/// Returns the median ratios (MedianRatio list in chromosome order).
fn diploid_normaliser(ratios: &mut [ChrRatios], v38: bool) -> Vec<MedianRatio> {
    diploid_normaliser_with(ratios, v38, ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE)
}

pub fn diploid_normaliser_with(ratios: &mut [ChrRatios], v38: bool, max_distance: i64, min_coverage: usize) -> Vec<MedianRatio> {
    let valid = |r: &BamRatio| doubles::greater_than(r.ratio(), 0.0);
    let mut median_ratios: Vec<MedianRatio> = ratios
        .iter()
        .map(|c| {
            let mut values: Vec<f64> = c.ratios.iter().filter(|r| valid(r)).map(|r| r.ratio()).collect();
            let count = values.len() as i32;
            MedianRatio { chromosome: c.chr.versioned(v38), median_ratio: stats::median(&mut values), count }
        })
        .collect();
    median_ratios.sort_by_key(|m| Chr::from_name(&m.chromosome).unwrap());
    let expected = expected_ratios(&median_ratios);
    for c in ratios.iter_mut() {
        if c.chr == Chr::Y {
            for r in c.ratios.iter_mut() {
                let v = r.ratio();
                r.set_diploid_adjusted_ratio(v);
            }
            continue;
        }
        let expected_ratio = expected.get(c.chr).unwrap_or(0.0);
        let values: Vec<f64> = c.ratios.iter().map(|r| r.ratio()).collect();
        let normalised = diploid_ratio_normalise(&values, expected_ratio, max_distance, min_coverage);
        for (r, v) in c.ratios.iter_mut().zip(normalised) {
            r.set_diploid_adjusted_ratio(v);
        }
    }
    median_ratios
}

/// DiploidRatioNormaliser.normalise over one chromosome's ratios in order: a valid ratio (above 1e-10) times
/// the expected ratio over the rolling median of the valid ratios at most `max_distance` valid ratios away,
/// when the window holds at least `min_coverage`; any other ratio as it is.
pub fn diploid_ratio_normalise(ratios: &[f64], expected_ratio: f64, max_distance: i64, min_coverage: usize) -> Vec<f64> {
    let values: Vec<f64> = ratios.iter().copied().filter(|r| doubles::greater_than(*r, 0.0)).collect();
    let mut rolling = stats::RollingMedian::new();
    let (mut start, mut end): (usize, i64) = (0, -1);
    let mut current = 0usize;
    let mut out = Vec::with_capacity(ratios.len());
    for &r in ratios {
        if !doubles::greater_than(r, 0.0) {
            out.push(r);
            continue;
        }
        let value = values[current];
        // removeExpiredRatios
        while start < current && (current - start) as i64 > max_distance {
            rolling.remove(values[start]);
            start += 1;
        }
        // addNewRatios
        let mut later = (end + 1) as usize;
        while later < values.len() && (later as i64 - current as i64).abs() <= max_distance {
            end += 1;
            rolling.add(values[later]);
            later += 1;
        }
        current += 1;
        let median = rolling.median();
        out.push(if rolling.size() >= min_coverage { expected_ratio * value / median } else { value });
    }
    out
}

/// ReferenceCalculation.medianRatios for a targeted panel: every human chromosome, the median of the
/// positive ratios (which become the diploid-adjusted ratios) and the number of windows
pub fn panel_median_ratios(ratios: &mut [ChrRatios], v38: bool) -> Vec<MedianRatio> {
    let mut out = Vec::new();
    for chr in Chr::all() {
        let mut values = Vec::new();
        let mut count = 0;
        if let Some(c) = ratios.iter_mut().find(|c| c.chr == chr) {
            count = c.ratios.len() as i32;
            for r in c.ratios.iter_mut() {
                if r.ratio() > 0.0 {
                    let v = r.ratio();
                    r.set_diploid_adjusted_ratio(v);
                    values.push(v);
                }
            }
        }
        out.push(MedianRatio { chromosome: chr.versioned(v38), median_ratio: doubles::median(&mut values), count });
    }
    out
}

/// TargetRegions.mergeTumorRatios: the first panel's ratios, then for each window not valid yet the next
/// panel's valid ratio
pub fn merge_panel_ratios(all: &mut Vec<ChrRatios>, panel: Vec<ChrRatios>) -> Result<()> {
    if all.iter().all(|c| c.ratios.is_empty()) {
        *all = panel;
        return Ok(());
    }
    for p in panel {
        let existing = all.iter_mut().find(|c| c.chr == p.chr);
        let existing_len = existing.as_ref().map_or(0, |e| e.ratios.len());
        if existing_len != p.ratios.len() {
            crate::error!("chromosome({}) inconsistent tumor ratio array size(existing={} new={})", p.chr.short_name(), existing_len, p.ratios.len());
            std::process::exit(1);
        }
        let Some(existing) = existing else { continue };
        for (e, n) in existing.ratios.iter_mut().zip(&p.ratios) {
            if !e.is_valid() && n.is_valid() {
                e.override_ratio(n.ratio());
            }
        }
    }
    Ok(())
}

/// One row of the ratio file (CobaltRatio)
#[derive(Clone, Debug, PartialEq)]
pub struct CobaltRatio {
    pub chr: Chr,
    pub position: i32,
    pub reference_read_depth: f64,
    pub reference_gc_ratio: f64,
    pub reference_gc_content: f64,
    pub reference_gc_diploid_ratio: f64,
    pub tumor_read_depth: f64,
    pub tumor_gc_ratio: f64,
    pub tumor_gc_content: f64,
}

/// The tumor and reference ratios of one chromosome side by side (ResultsCollator), either may be empty
pub struct CollatedChr<'a> {
    pub chr: Chr,
    pub tumor: &'a [BamRatio],
    pub reference: &'a [BamRatio],
}

impl CollatedChr<'_> {
    pub fn len(&self) -> usize {
        if self.tumor.is_empty() { self.reference.len() } else { self.tumor.len() }
    }
    /// ResultsCollator.merge
    pub fn row(&self, i: usize) -> CobaltRatio {
        let (tr, rr) = (self.tumor.get(i), self.reference.get(i));
        CobaltRatio {
            chr: self.chr,
            position: tr.or(rr).unwrap().position,
            reference_read_depth: rr.map_or(-1.0, |x| x.read_depth),
            reference_gc_ratio: rr.map_or(-1.0, |x| x.ratio()),
            reference_gc_content: rr.map_or(-1.0, |x| x.gc_content),
            reference_gc_diploid_ratio: rr.map_or(-1.0, |x| x.diploid_adjusted_ratio()),
            tumor_read_depth: tr.map_or(-1.0, |x| x.read_depth),
            tumor_gc_ratio: tr.map_or(-1.0, |x| x.ratio()),
            tumor_gc_content: tr.map_or(-1.0, |x| x.gc_content),
        }
    }
}

/// ResultsCollator.collateResults, in the order CobaltRatioFile.write sorts the rows (chromosome rank, then
/// position)
pub fn collate<'a>(tumor: &'a [ChrRatios], reference: &'a [ChrRatios]) -> Result<Vec<CollatedChr<'a>>> {
    let empty = |s: &[ChrRatios]| s.iter().all(|c| c.ratios.is_empty());
    if empty(tumor) && empty(reference) {
        bail!("java.lang.IllegalArgumentException (ResultsCollator: no ratios)");
    }
    let keys = |s: &[ChrRatios]| {
        let mut k: Vec<Chr> = s.iter().filter(|c| !c.ratios.is_empty()).map(|c| c.chr).collect();
        k.sort();
        k
    };
    if !empty(tumor) && !empty(reference) && keys(tumor) != keys(reference) {
        bail!("java.lang.IllegalArgumentException (ResultsCollator: the tumor and reference BAMs have different chromosomes)");
    }
    let primary = if empty(tumor) { reference } else { tumor };
    let mut out = Vec::new();
    for c in primary.iter().filter(|c| !c.ratios.is_empty()) {
        let t: &[BamRatio] = tumor.iter().find(|x| x.chr == c.chr).map(|x| &x.ratios[..]).unwrap_or(&[]);
        let r: &[BamRatio] = reference.iter().find(|x| x.chr == c.chr).map(|x| &x.ratios[..]).unwrap_or(&[]);
        if !t.is_empty() && !r.is_empty() {
            if t.len() != r.len() {
                bail!("java.lang.IllegalStateException (ResultsCollator: chromosome {} has {} tumor and {} reference windows)", c.chr.short_name(), t.len(), r.len());
            }
            if t.iter().zip(r).any(|(a, b)| a.position != b.position) {
                bail!("java.lang.IllegalArgumentException (ResultsCollator: positions differ)");
            }
        }
        out.push(CollatedChr { chr: c.chr, tumor: t, reference: r });
    }
    out.sort_by_key(|c| c.chr);
    Ok(out)
}

/// A double as Java's logs print it ({} placeholders: Double.toString), close enough for log lines
pub fn java_double(d: f64) -> String {
    format!("{d:?}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consolidation_counts() {
        // ResultsConsolidatorTest
        assert_eq!(calc_consolidation_count(100.0), 1);
        assert_eq!(calc_consolidation_count(8.0), 10);
        assert_eq!(calc_consolidation_count(4.0), 20);
        assert_eq!(calc_consolidation_count(0.3), 300);
        assert_eq!(calc_consolidation_count(0.0281), 1000);
        assert_eq!(calc_consolidation_count(0.08), 1000);
        assert_eq!(calc_consolidation_count(0.21), 400);
    }
}
