//! The fitting package: PurityAdjuster, PloidyDeviation, RegionFitCalculator, AmbiguousBaf, the purity and
//! ploidy grid (FittedPurityFactory), FittedPurityScoreFactory, BestFit, AneuploidyDetector and
//! PurityPloidyFitter.

use std::collections::BTreeMap;

use hmf_common::commons_math::{normal_cdf, Binomial};

use crate::config::{ambiguous_baf_threshold, FittingConfig, PurpleConfig};
use crate::constants::*;
use crate::inputs::{AmberBaf, CobaltChromosomes};
use crate::region::ObservedRegion;
use crate::reseg::RatioBucketSeries;
use crate::types::*;
use crate::variant::{AllelicDepth, Variant};

#[derive(Clone, Debug)]
pub struct PurityAdjuster {
    pub purity: f64,
    pub norm_factor: f64,
    /// CobaltChromosome.actualRatio by HumanChromosome ordinal (0 where the chromosome is absent)
    pub ratios: [f64; 24],
}

impl PurityAdjuster {
    pub fn new(purity: f64, norm_factor: f64, cobalt: &CobaltChromosomes) -> PurityAdjuster {
        let mut ratios = [0.0; 24];
        for c in cobalt.chromosomes.values() {
            ratios[c.chromosome.ordinal()] = c.actual_ratio;
        }
        PurityAdjuster { purity, norm_factor, ratios }
    }
    pub fn with_ratios(ratios: [f64; 24], purity: f64, norm_factor: f64) -> PurityAdjuster {
        PurityAdjuster { purity, norm_factor, ratios }
    }
    fn germline_ratio(&self, chromosome: &str) -> f64 {
        self.ratios[chr(chromosome).ordinal()]
    }
    pub fn germline_copy_number(&self, chromosome: &str) -> f64 {
        self.germline_ratio(chromosome) * 2.0
    }
    pub fn purity_adjusted_copy_number_chr(&self, chromosome: &str, ratio: f64) -> f64 {
        self.purity_adjusted_copy_number(ratio, self.germline_ratio(chromosome))
    }
    pub fn purity_adjusted_copy_number(&self, tumor_ratio: f64, normal_ratio: f64) -> f64 {
        if is_zero(tumor_ratio) {
            0.0
        } else {
            2.0 * normal_ratio + 2.0 * (tumor_ratio - normal_ratio * self.norm_factor) / self.purity / self.norm_factor
        }
    }
    pub fn purity_adjusted_vaf(&self, chromosome: &str, copy_number: f64, observed: f64) -> f64 {
        self.purity_adjusted_vaf_het(chromosome, copy_number, observed, false)
    }
    pub fn purity_adjusted_vaf_het(&self, chromosome: &str, copy_number: f64, observed: f64, germline_het_deletion: bool) -> f64 {
        let typical = if !germline_het_deletion { self.germline_copy_number(chromosome) } else { 1.0 };
        self.purity_adjusted_frequency(typical, 0.0, copy_number, observed)
    }
    pub fn purity_adjusted_baf_simple(&self, chromosome: &str, copy_number: f64, observed: f64) -> f64 {
        let typical = self.germline_copy_number(chromosome);
        if typical < 2.0 || less_or_equal(copy_number, 1.0) {
            return 1.0;
        }
        self.purity_adjusted_frequency(2.0, 1.0, copy_number, observed)
    }
    pub fn purity_adjusted_frequency(&self, normal_cn: f64, normal_ploidy: f64, tumor_cn: f64, observed: f64) -> f64 {
        self.purity_adjusted_ploidy(normal_cn, normal_ploidy, tumor_cn, observed) / tumor_cn
    }
    pub fn purity_adjusted_ploidy(&self, normal_cn: f64, normal_ploidy: f64, tumor_cn: f64, observed: f64) -> f64 {
        let total = self.purity * tumor_cn + normal_cn * (1.0 - self.purity);
        let normal_obs = normal_ploidy * (1.0 - self.purity);
        (observed * total - normal_obs) / self.purity
    }
    pub fn expected_frequency(&self, normal_cn: f64, normal_ploidy: i32, tumor_cn: f64, tumor_ploidy: f64) -> f64 {
        if less_or_equal(tumor_cn, 0.0) {
            return 0.0;
        }
        let total = self.purity * tumor_cn + normal_cn * (1.0 - self.purity);
        let normal_obs = normal_ploidy as f64 * (1.0 - self.purity);
        let tumor_obs = tumor_ploidy * self.purity;
        (normal_obs + tumor_obs) / total
    }
    pub fn purity_adjusted_vaf_with_heterozygous_normal(&self, chromosome: &str, copy_number: f64, observed: f64) -> f64 {
        let typical = self.germline_copy_number(chromosome);
        self.purity_adjusted_frequency(typical, typical / 2.0, copy_number, observed)
    }
    pub fn purity_adjusted_vaf_with_homozygous_normal(&self, chromosome: &str, copy_number: f64, observed: f64) -> f64 {
        let typical = self.germline_copy_number(chromosome);
        self.purity_adjusted_frequency(typical, typical, copy_number, observed)
    }
    pub fn implied_ploidy(average_ratio: f64, purity: f64, norm_factor: f64) -> f64 {
        (average_ratio - norm_factor) / purity / norm_factor * 2.0 + 2.0
    }
    pub fn implied_norm_factor(average_ratio: f64, purity: f64, tumor_ploidy: f64) -> f64 {
        2.0 * average_ratio / (2.0 - 2.0 * purity + tumor_ploidy * purity)
    }
}

/// PloidyDeviation
pub struct PloidyDeviation {
    sd: f64,
    min_sd_per_ploidy: f64,
    major_sub_one_multiplier: f64,
    major_sub_min_additional: f64,
    min_deviation: f64,
}

impl PloidyDeviation {
    pub fn new(c: &FittingConfig) -> PloidyDeviation {
        PloidyDeviation {
            sd: c.ploidy_penalty_standard_deviation,
            min_sd_per_ploidy: c.ploidy_penalty_min_standard_deviation_per_ploidy,
            major_sub_one_multiplier: c.ploidy_penalty_major_allele_sub_one_multiplier.abs(),
            major_sub_min_additional: c.ploidy_penalty_major_allele_sub_one_additional.abs(),
            min_deviation: c.ploidy_penalty_min_deviation,
        }
    }
    pub fn major_allele_deviation(&self, purity: f64, norm_factor: f64, ploidy: f64) -> f64 {
        let multiplier = if greater_than(ploidy, 0.0) && less_than(ploidy, 1.0) { jmax(1.0, self.major_sub_one_multiplier * (1.0 - ploidy)) } else { 1.0 };
        let deviation = multiplier * self.allele_deviation(purity, norm_factor, ploidy) + self.sub_min_additional(1.0, ploidy);
        jmax(deviation, self.min_deviation)
    }
    pub fn minor_allele_deviation(&self, purity: f64, norm_factor: f64, ploidy: f64) -> f64 {
        let deviation = self.allele_deviation(purity, norm_factor, ploidy) + self.sub_min_additional(0.0, ploidy);
        jmax(deviation, self.min_deviation)
    }
    fn allele_deviation(&self, purity: f64, norm_factor: f64, ploidy: f64) -> f64 {
        let distance = if less_than(ploidy, -0.5) { 0.5 } else { abs_distance_from_integer(ploidy) };
        let sd_per_ploidy = jmax(self.min_sd_per_ploidy, purity * norm_factor / 2.0 / self.sd);
        let value = distance * sd_per_ploidy;
        2.0 * normal_cdf(0.0, 1.0, value) - 1.0 + jmax(-0.5 - ploidy, 0.0)
    }
    fn sub_min_additional(&self, min_ploidy: f64, ploidy: f64) -> f64 {
        jmin(self.major_sub_min_additional, jmax(0.0, -self.major_sub_min_additional * (ploidy - min_ploidy)))
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RegionFitCalcs {
    pub tumor_copy_number: f64,
    pub tumor_baf: f64,
    pub ref_normalised_copy_number: f64,
    pub minor_allele_copy_number_deviation: f64,
    pub major_allele_copy_number_deviation: f64,
    pub event_penalty: f64,
    pub deviation_penalty: f64,
}

impl RegionFitCalcs {
    pub fn major_allele_copy_number(&self) -> f64 {
        self.tumor_baf * self.tumor_copy_number
    }
    pub fn minor_allele_copy_number(&self) -> f64 {
        self.tumor_copy_number - self.major_allele_copy_number()
    }
    pub fn is_diploid(&self) -> bool {
        ObservedRegion::is_diploid_cn(self.major_allele_copy_number()) && ObservedRegion::is_diploid_cn(self.minor_allele_copy_number())
    }
}

pub struct RegionFitCalculator<'a> {
    pub ambiguous_baf: f64,
    pub config: &'a FittingConfig,
    pub deviation: PloidyDeviation,
    pub cobalt: &'a CobaltChromosomes,
}

/// RegionFitCalculator.calculateEventPenalty
pub fn event_penalty(factor: f64, major: f64, minor: f64) -> f64 {
    let wgd = 1.0 + (major - 2.0).abs() + (minor - 2.0).abs();
    let single = (major - 1.0).abs() + (minor - 1.0).abs();
    1.0 + factor * jmin(single, wgd)
}

const NO_CALC_BAF: f64 = -1.0;

/// RegionFitCalculator.estimateMinMaxBaf
pub fn estimate_min_max_baf(cn: f64, min_baf: f64, max_baf: f64) -> f64 {
    let major_min = min_baf * cn;
    let major_max = max_baf * cn;
    let minor_min = cn - major_min;
    let minor_max = cn - major_max;
    let minor_max_ceil = minor_max.ceil();
    let signum = |v: f64| if v > 0.0 { 1.0 } else if v < 0.0 { -1.0 } else { v };
    let minor_diff = !equal(signum(minor_max_ceil - minor_min), signum(minor_max_ceil - minor_max));
    let major_min_ceil = major_min.ceil();
    let major_diff = !equal(signum(major_min_ceil - major_min), signum(major_min_ceil - major_max));
    if !minor_diff && !major_diff {
        return NO_CALC_BAF;
    }
    if !minor_diff {
        return major_min_ceil / cn;
    } else if !major_diff {
        return 1.0 - minor_max_ceil / cn;
    }
    let half = cn * 0.5;
    let major_max_floor = major_max.floor();
    let est_high = if major_max_floor >= half { major_max_floor } else { major_max };
    let mut est_low = est_high - 1.0;
    let lower_bound = jmax(half, half);
    let mut minor_int = 0;
    while (minor_int as f64) < cn {
        let major = cn - minor_int as f64;
        if major > minor_int as f64 && major - minor_int as f64 <= 0.5 {
            est_low = major;
            break;
        }
        if major >= lower_bound {
            est_low = major;
        } else if major < minor_int as f64 {
            let previous_floor = (major + 1.0).floor();
            if previous_floor >= lower_bound {
                est_low = previous_floor;
            }
            break;
        }
        minor_int += 1;
    }
    est_low = jmin(est_low, est_high);
    let minor_high = cn - est_high;
    let minor_low = cn - est_low;
    let diff_high = est_high - minor_high;
    let diff_low = est_low - minor_low;
    let est = if diff_high <= diff_low { est_high } else { est_low };
    jmax(est / cn, BAF_PNT_5)
}

impl RegionFitCalculator<'_> {
    pub fn new<'a>(cobalt: &'a CobaltChromosomes, config: &'a FittingConfig, ambiguous_baf: f64) -> RegionFitCalculator<'a> {
        RegionFitCalculator { ambiguous_baf, config, deviation: PloidyDeviation::new(config), cobalt }
    }

    /// RegionFitCalculator.fitRegion over regions on chromosomes COBALT has
    pub fn fit_regions(&self, purity: f64, norm_factor: f64, regions: &[ObservedRegion]) -> Vec<ObservedRegion> {
        regions.iter().filter(|r| self.cobalt.has_chromosome(&r.chromosome)).map(|r| self.fit_region(purity, norm_factor, r)).collect()
    }

    pub fn fit_region(&self, purity: f64, norm_factor: f64, region: &ObservedRegion) -> ObservedRegion {
        let c = self.calculate_region_fit(purity, norm_factor, region);
        let mut r = region.clone();
        r.tumor_copy_number = c.tumor_copy_number;
        r.tumor_baf = c.tumor_baf;
        r.ref_normalised_copy_number = c.ref_normalised_copy_number;
        r.minor_allele_copy_number_deviation = c.minor_allele_copy_number_deviation;
        r.major_allele_copy_number_deviation = c.major_allele_copy_number_deviation;
        r.deviation_penalty = c.deviation_penalty;
        r.event_penalty = c.event_penalty;
        r
    }

    pub fn calculate_region_fit(&self, purity: f64, norm_factor: f64, region: &ObservedRegion) -> RegionFitCalcs {
        let adjuster = PurityAdjuster::new(purity, norm_factor, self.cobalt);
        self.calculate_region_fit_with(&adjuster, region)
    }

    pub fn calculate_region_fit_with(&self, adjuster: &PurityAdjuster, region: &ObservedRegion) -> RegionFitCalcs {
        let (purity, norm_factor) = (adjuster.purity, adjuster.norm_factor);
        let observed_tumor_ratio = region.observed_tumor_ratio;
        let implied_cn = adjuster.purity_adjusted_copy_number_chr(&region.chromosome, observed_tumor_ratio);
        let observed_baf = region.observed_baf;
        let implied_baf = self.implied_baf(adjuster, &region.chromosome, implied_cn, observed_baf);
        let ref_normalised = adjuster.purity_adjusted_copy_number(observed_tumor_ratio, region.observed_normal_ratio);
        let major = implied_baf * implied_cn;
        let minor = implied_cn - major;
        let major_dev = self.deviation.major_allele_deviation(purity, norm_factor, major);
        let minor_dev = self.deviation.minor_allele_deviation(purity, norm_factor, minor);
        let event = event_penalty(self.config.ploidy_penalty_factor, major, minor);
        let mut deviation = (minor_dev + major_dev) * observed_baf;
        if self.config.gc_ratio_exponent > 0.0 || self.config.deviation_penalty_gc_min_adjust > 0.0 {
            let adj = if self.config.gc_ratio_exponent > 0.0 { jpow(observed_tumor_ratio, self.config.gc_ratio_exponent) } else { 1.0 };
            deviation /= jmax(self.config.deviation_penalty_gc_min_adjust, adj);
        }
        RegionFitCalcs {
            tumor_copy_number: implied_cn,
            tumor_baf: implied_baf,
            ref_normalised_copy_number: if ref_normalised.is_nan() { 0.0 } else { ref_normalised },
            minor_allele_copy_number_deviation: minor_dev,
            major_allele_copy_number_deviation: major_dev,
            event_penalty: event,
            deviation_penalty: deviation,
        }
    }

    fn implied_baf(&self, adjuster: &PurityAdjuster, chromosome: &str, cn: f64, observed_baf: f64) -> f64 {
        let Some(c) = self.cobalt.get(chromosome) else { return 1.0 };
        if !c.is_normal() || !c.is_diploid() || less_or_equal(cn, 0.1) {
            return 1.0;
        }
        if less_or_equal(observed_baf, self.ambiguous_baf) {
            self.baf_to_minimise_deviation(adjuster, chromosome, cn, observed_baf)
        } else {
            adjuster.purity_adjusted_baf_simple(chromosome, cn, observed_baf)
        }
    }

    pub fn baf_to_minimise_deviation(&self, adjuster: &PurityAdjuster, chromosome: &str, cn: f64, observed_baf: f64) -> f64 {
        let min_baf = jmax(0.0, jmin(1.0, adjuster.purity_adjusted_baf_simple(chromosome, cn, BAF_PNT_5)));
        let max_baf = jmax(0.0, jmin(1.0, adjuster.purity_adjusted_baf_simple(chromosome, cn, observed_baf)));
        let est = estimate_min_max_baf(cn, min_baf, max_baf);
        if est != NO_CALC_BAF {
            return est;
        }
        let (purity, norm) = (adjuster.purity, adjuster.norm_factor);
        let (major_min, major_max) = (min_baf * cn, max_baf * cn);
        let (minor_min, minor_max) = (cn - major_min, cn - major_max);
        let min_dev = self.deviation.major_allele_deviation(purity, norm, major_min) + self.deviation.minor_allele_deviation(purity, norm, minor_min);
        let max_dev = self.deviation.major_allele_deviation(purity, norm, major_max) + self.deviation.minor_allele_deviation(purity, norm, minor_max);
        if less_than(min_dev, max_dev) { BAF_PNT_5 } else { observed_baf }
    }
}

// ---------------------------------------------------------------------------------------------
// AmbiguousBaf

/// AmbiguousBaf.depthBasedBAF
pub fn depth_based_baf(average_depth: i32) -> f64 {
    let min_depth = average_depth / 2;
    let max_depth = average_depth * 3 / 2;
    let mut total = 0.0;
    for i in min_depth..max_depth {
        let d = Binomial::new(i, BAF_PNT_5);
        total += d.inverse_cumulative_probability(STANDARD_BAF_TARGET_PERC) as f64 / i as f64;
    }
    let expected = total / (max_depth - min_depth) as f64;
    jmax(expected, ambiguous_baf_threshold())
}

/// AmbiguousBaf.findAmbiguousBaf
pub fn find_ambiguous_baf(segments: &[ObservedRegion], average_tumor_depth: i32) -> f64 {
    let depth_based = depth_based_baf(average_tumor_depth);
    let mut raw = RatioBucketSeries::with_range(0.5, 0.6, 0.005);
    for s in segments {
        if s.germline_status != GermlineStatus::Diploid || s.baf_count <= 0 {
            continue;
        }
        if s.observed_tumor_ratio <= 0.0 || s.observed_tumor_ratio >= AMBIGUOUS_BAF_OBSERVED_RATIO_MAX {
            continue;
        }
        let ratio = jround(s.observed_baf / 0.005) as f64 * 0.005;
        raw.add_count(ratio, s.baf_count as f64);
    }
    let bafs = raw.build_smoothed_series();
    if bafs.len() < 2 {
        return -1.0;
    }
    let mut peaks: Vec<usize> = (1..bafs.len() - 1).filter(|&i| bafs[i].value > bafs[i - 1].value && bafs[i].value > bafs[i + 1].value).collect();
    if bafs[0].value > bafs[1].value {
        peaks.insert(0, 0);
    }
    if let Some((index, peak)) = closest_strongly_supported_peak(&bafs, &peaks) {
        for b in bafs.iter().skip(index + 1) {
            if b.value < AMBIGUOUS_BAF_EDGE_RATIO * peak.value {
                let ambiguous = jmax(jmin((b.level + peak.level) * 0.5, AMBIGUOUS_BAF_MAX), AMBIGUOUS_BAF_MIN);
                return jmax(ambiguous, depth_based);
            }
        }
    }
    depth_based
}

#[derive(Clone, Copy)]
struct PeakCandidate {
    index: usize,
    baf: crate::reseg::Bucket,
    lower_edge: f64,
    upper_edge: f64,
}

fn closest_strongly_supported_peak(bafs: &[crate::reseg::Bucket], peaks: &[usize]) -> Option<(usize, crate::reseg::Bucket)> {
    let mut candidates: Vec<PeakCandidate> = Vec::new();
    for &pi in peaks {
        let b = bafs[pi];
        if b.value < AMBIGUOUS_BAF_MIN_VALUE as f64 {
            continue;
        }
        let (mut lower, mut upper): (Option<f64>, Option<f64>) = (None, None);
        for up in [true, false] {
            let mut i = pi as i64 + if up { 1 } else { -1 };
            while i >= 0 && (i as usize) < bafs.len() {
                let n = bafs[i as usize];
                if n.value < b.value * AMBIGUOUS_BAF_EDGE_RATIO {
                    if up {
                        upper = Some(n.level);
                    } else {
                        lower = Some(n.level);
                    }
                    break;
                }
                i += if up { 1 } else { -1 };
            }
        }
        if pi == 0 {
            lower = Some(b.level);
        }
        if let (Some(l), Some(u)) = (lower, upper) {
            candidates.push(PeakCandidate { index: pi, baf: b, lower_edge: l, upper_edge: u });
        }
    }
    if candidates.is_empty() {
        return None;
    }
    // consolidatePeaks
    while candidates.len() > 1 {
        let mut found = false;
        let mut index = 0;
        while index + 1 < candidates.len() {
            let (c, n) = (candidates[index], candidates[index + 1]);
            if c.baf.level < n.lower_edge && c.upper_edge < n.baf.level {
                index += 1;
                continue;
            }
            let (lv, uv) = (c.baf.value, n.baf.value);
            let mut lowest = jmax(lv, uv);
            for b in &bafs[c.index + 1..n.index] {
                lowest = jmin(lowest, b.value);
            }
            let keep = lv > lowest * AMBIGUOUS_BAF_CONSOLIDATE_RATIO && uv > lowest * AMBIGUOUS_BAF_CONSOLIDATE_RATIO;
            if keep {
                index += 1;
            } else {
                found = true;
                if c.baf.value < n.baf.value {
                    candidates.remove(index);
                } else {
                    candidates.remove(index + 1);
                }
            }
        }
        if !found {
            break;
        }
    }
    let mut closest: Option<PeakCandidate> = None;
    let mut closest_diff = 0.0;
    for c in &candidates {
        let diff = (c.baf.level - BAF_PNT_5).abs();
        if closest.is_none() || diff < closest_diff {
            closest_diff = diff;
            closest = Some(*c);
        }
    }
    closest.map(|c| (c.index, c.baf))
}

// ---------------------------------------------------------------------------------------------
// FittedPurity

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FittedPurity {
    pub purity: f64,
    pub norm_factor: f64,
    pub ploidy: f64,
    pub score: f64,
    pub diploid_proportion: f64,
    pub somatic_penalty: f64,
}

/// Collections.sort of FittedPurity: by score (Double.compare), stable
pub fn sort_by_score(v: &mut [FittedPurity]) {
    v.sort_by(|a, b| dcmp(a.score, b.score));
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FittedPurityScore {
    pub min_purity: f64,
    pub max_purity: f64,
    pub min_ploidy: f64,
    pub max_ploidy: f64,
    pub min_diploid_proportion: f64,
    pub max_diploid_proportion: f64,
}

impl FittedPurityScore {
    pub fn purity_spread(&self) -> f64 {
        self.max_purity - self.min_purity
    }
}

/// FittedPurityScoreFactory.score
pub fn score(purities: &[FittedPurity]) -> FittedPurityScore {
    let (mut max_ploidy, mut min_ploidy) = (0.0, f64::MAX);
    let (mut max_purity, mut min_purity) = (0.0, f64::MAX);
    let (mut max_dp, mut min_dp) = (0.0, f64::MAX);
    for f in purities {
        max_ploidy = jmax(max_ploidy, f.ploidy);
        min_ploidy = jmin(min_ploidy, f.ploidy);
        max_purity = jmax(max_purity, f.purity);
        min_purity = jmin(min_purity, f.purity);
        max_dp = jmax(max_dp, f.diploid_proportion);
        min_dp = jmin(min_dp, f.diploid_proportion);
    }
    FittedPurityScore { min_purity, max_purity, min_ploidy, max_ploidy, min_diploid_proportion: min_dp, max_diploid_proportion: max_dp }
}

/// FittedPurityScoreFactory.isPolyclonal / polyclonalProportion
pub fn polyclonal_proportion(copy_numbers: &[crate::copynumber::PurpleCopyNumber]) -> f64 {
    let (mut poly, mut total) = (0, 0);
    for c in copy_numbers {
        total += c.baf_count;
        let remainder = (c.average_tumor_copy_number - jround(c.average_tumor_copy_number) as f64).abs();
        if greater_than(remainder, 0.25) {
            poly += c.baf_count;
        }
    }
    if total == 0 { 0.0 } else { poly as f64 / total as f64 }
}

/// BestFit.mostDiploidPerPurity
pub fn most_diploid_per_purity(all: &[FittedPurity]) -> Vec<FittedPurity> {
    let mut sortable = all.to_vec();
    sortable.sort_by(|a, b| dcmp((DIPLOID_PLOIDY - a.ploidy).abs(), (DIPLOID_PLOIDY - b.ploidy).abs()));
    // TreeSet<Double> with Doubles.comparator(): purities within 1e-10 are one
    let mut purities: Vec<f64> = Vec::new();
    let mut result = Vec::new();
    for f in sortable {
        if !purities.iter().any(|p| equal(*p, f.purity)) {
            purities.push(f.purity);
            result.push(f);
        }
    }
    sort_by_score(&mut result);
    result
}

/// PurityPloidyFitter.findMatchedFittedPurity
pub fn find_matched_fitted_purity(purity: f64, candidates: &[FittedPurity]) -> Option<FittedPurity> {
    let mut closest: Option<FittedPurity> = None;
    let mut closest_diff = 0.0;
    let epsilon = PURITY_INCREMENT_DEFAULT * 0.25;
    for f in candidates {
        if (f.ploidy - 2.0).abs() > 0.005 {
            continue;
        }
        let diff = (f.purity - purity).abs();
        if closest.is_none() || diff < closest_diff {
            if diff < epsilon {
                return Some(*f);
            }
            closest_diff = diff;
            closest = Some(*f);
        }
    }
    closest
}

/// FittedPurityFactory.ploidyRange
pub fn ploidy_range(min_ploidy: f64, max_ploidy: f64) -> Vec<f64> {
    let sequence = |min: f64, max: f64, inc: f64| {
        let mut out = Vec::new();
        let mut p = min;
        while less_than(p, max) {
            out.push(p);
            p = round(p + inc, 2);
        }
        out
    };
    let mut out = sequence(jmax(0.0, min_ploidy), jmin(3.0, max_ploidy), 0.02);
    out.extend(sequence(jmax(3.0, min_ploidy), jmin(5.0, max_ploidy), 0.05));
    out.extend(sequence(jmax(5.0, min_ploidy), max_ploidy, 0.1));
    out.push(max_ploidy);
    out
}

/// FittedPurityFactory.useRegionToFitPurity
pub fn use_region_to_fit_purity(tumor_only: bool, cobalt: &CobaltChromosomes, region: &ObservedRegion) -> bool {
    if region.baf_count <= 0 || !positive_or_zero(region.observed_tumor_ratio) || region.germline_status != GermlineStatus::Diploid {
        return false;
    }
    if greater_than(region.observed_tumor_ratio, 3.0) {
        return false;
    }
    let Some(c) = cobalt.get(&region.chromosome) else { return false };
    if tumor_only && cobalt.gender == Gender::Male && !c.is_autosome() {
        return false;
    }
    if tumor_only && cobalt.gender == Gender::Female && is_y(c.chromosome) {
        return false;
    }
    c.is_normal() && c.is_diploid()
}

/// Downsample.downsample
pub fn downsample<T: Clone>(max_entries: usize, input: &[T]) -> Vec<T> {
    if input.len() <= max_entries {
        return input.to_vec();
    }
    let scale = jround((input.len() as f64 / max_entries as f64).ceil()) as usize;
    input.iter().enumerate().filter(|(i, _)| i % scale == 0).map(|(_, x)| x.clone()).collect()
}

/// a fitting region and its somatic variants (allelic depth only)
struct RegionData<'a> {
    region: &'a ObservedRegion,
    variants: Vec<AllelicDepth>,
}

pub struct FittedPurityFactory<'a> {
    config: &'a PurpleConfig,
    cobalt: &'a CobaltChromosomes,
    calculator: &'a RegionFitCalculator<'a>,
    regions: Vec<RegionData<'a>>,
    total_baf_count: i32,
    pub average_fitting_ratio: f64,
    ploidy_range: Vec<f64>,
}

impl<'a> FittedPurityFactory<'a> {
    pub fn new(
        config: &'a PurpleConfig,
        cobalt: &'a CobaltChromosomes,
        calculator: &'a RegionFitCalculator<'a>,
        observed: &'a [ObservedRegion],
        variants: &[&Variant],
    ) -> FittedPurityFactory<'a> {
        let tumor_only = config.tumor_only_mode();
        let mut regions: Vec<RegionData> = Vec::new();
        let mut filtered_variants: Vec<&Variant> = Vec::new();
        let (mut baf_count, mut weighted_ratio) = (0i32, 0.0);
        // GenomePositionSelectorFactory.create(List): the variants within each fitting region
        let in_region = |v: &Variant, r: &ObservedRegion| v.chromosome() == r.chromosome && v.position() >= r.start && v.position() <= r.end;
        for r in observed {
            if use_region_to_fit_purity(tumor_only, cobalt, r) {
                for v in variants.iter().filter(|v| in_region(v, r)) {
                    filtered_variants.push(v);
                }
                regions.push(RegionData { region: r, variants: Vec::new() });
                baf_count += r.baf_count;
                weighted_ratio += r.baf_count as f64 * r.observed_tumor_ratio;
            }
        }
        let average = if baf_count > 0 { weighted_ratio / baf_count as f64 } else { 0.0 };
        let downsampled = downsample(1000, &filtered_variants);
        for rd in regions.iter_mut() {
            for v in downsampled.iter().filter(|v| in_region(v, rd.region)) {
                rd.variants.push(v.tumor_depth.unwrap_or_default());
            }
        }
        FittedPurityFactory {
            config,
            cobalt,
            calculator,
            regions,
            total_baf_count: baf_count,
            average_fitting_ratio: average,
            ploidy_range: ploidy_range(config.fitting.min_ploidy, config.fitting.max_ploidy),
        }
    }

    pub fn valid_data_for_fit(&self) -> bool {
        self.average_fitting_ratio > 0.0
    }

    /// FittedPurityFactory.fitPurity: every purity and ploidy, sorted by score
    pub fn fit_purity(&self) -> Vec<FittedPurity> {
        let mut purities = Vec::new();
        let mut p = self.config.fitting.min_purity;
        while p <= self.config.fitting.max_purity {
            purities.push(p);
            p = round(p + self.config.fitting.purity_increment, 2);
        }
        let threads = self.config.threads.max(1) as usize;
        let per_purity: Vec<Vec<FittedPurity>> = crate::parallel::map(&purities, threads, |&purity| {
            self.ploidy_range
                .iter()
                .map(|&ploidy| {
                    let norm = PurityAdjuster::implied_norm_factor(self.average_fitting_ratio, purity, ploidy);
                    self.fit(purity, norm)
                })
                .collect()
        });
        let mut all: Vec<FittedPurity> = per_purity.into_iter().flatten().collect();
        sort_by_score(&mut all);
        all
    }

    fn weight(&self, value: f64, baf_count: i32) -> f64 {
        value * baf_count as f64 / self.total_baf_count as f64
    }

    fn fit(&self, purity: f64, norm_factor: f64) -> FittedPurity {
        let (mut event, mut deviation, mut diploid, mut ploidy, mut somatic_total) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let mut somatic_count = 0;
        let weight = self.config.somatic_fitting.somatic_penalty_weight;
        let adjuster = PurityAdjuster::new(purity, norm_factor, self.cobalt);
        for rd in &self.regions {
            let r = rd.region;
            let c = self.calculator.calculate_region_fit_with(&adjuster, r);
            event += self.weight(c.event_penalty, r.baf_count);
            deviation += self.weight(c.deviation_penalty, r.baf_count);
            ploidy += self.weight(c.tumor_copy_number, r.baf_count);
            if c.is_diploid() {
                diploid += self.weight(1.0, r.baf_count);
            }
            for depth in &rd.variants {
                somatic_count += 1;
                let penalty = crate::fittingsnv::deviation_from_max(&adjuster, &r.chromosome, c.major_allele_copy_number(), c.tumor_copy_number, depth);
                somatic_total += weight * penalty;
            }
        }
        let somatic_penalty = if weight > 0.0 && somatic_count > 0 { somatic_total / somatic_count as f64 } else { 0.0 };
        FittedPurity { purity, norm_factor, ploidy, score: event * deviation + somatic_penalty, diploid_proportion: diploid, somatic_penalty }
    }
}

// ---------------------------------------------------------------------------------------------
// AneuploidyDetector

pub struct AneuploidyDetector<'a> {
    regions: Vec<&'a ObservedRegion>,
    bafs: &'a BTreeMap<Chromosome, Vec<AmberBaf>>,
    pub has_high_baf_region: bool,
    total_baf_count: i32,
}

impl<'a> AneuploidyDetector<'a> {
    pub fn new(regions: &'a [ObservedRegion], bafs: &'a BTreeMap<Chromosome, Vec<AmberBaf>>, cobalt: &CobaltChromosomes) -> AneuploidyDetector<'a> {
        let regions = regions
            .iter()
            .filter(|r| r.baf_count > 0 && positive_or_zero(r.observed_tumor_ratio) && r.germline_status == GermlineStatus::Diploid && cobalt.has_chromosome(&r.chromosome))
            .collect();
        AneuploidyDetector { regions, bafs, has_high_baf_region: false, total_baf_count: 0 }
    }

    pub fn has_aneuploidy(&mut self) -> bool {
        let mut high = 0;
        for r in &self.regions {
            self.total_baf_count += r.baf_count;
            if r.observed_baf > SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF && r.baf_count > SOMATIC_FIT_ANEUPLOIDIC_REGION_MIN_BAF_COUNT {
                high += r.baf_count;
                if r.observed_baf > HIGHLY_ANEUPLOIDIC_REFIT_BAF_CUTOFF {
                    self.has_high_baf_region = true;
                }
            }
        }
        let ratio = high as f64 / self.total_baf_count as f64;
        crate::debug!("aneuploidy ratio({})", hmf_common::fmt::java_fixed(ratio, 3));
        if greater_or_equal(ratio, SOMATIC_FIT_ANEUPLOIDIC_RATIO_CUTOFF) {
            return true;
        }
        let required = jmax(HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT, self.total_baf_count as f64 * HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF);
        // AmberPointsProvider: an iterator selector that only moves forward
        let mut cursors: BTreeMap<Chromosome, usize> = BTreeMap::new();
        let mut found = false;
        for r in &self.regions {
            if r.observed_baf < HIGHLY_ANEUPLOIDIC_REGION_CUTOFF || (r.baf_count as f64) < required {
                continue;
            }
            let c = chr(&r.chromosome);
            let Some(list) = self.bafs.get(&c) else { continue };
            let cur = cursors.entry(c).or_insert(0);
            while *cur < list.len() && list[*cur].position < r.start {
                *cur += 1;
            }
            let mut points = Vec::new();
            while *cur < list.len() && list[*cur].position <= r.end {
                points.push(&list[*cur]);
                *cur += 1;
            }
            let (mut has_high, mut has_low) = (false, false);
            let mut aneuploidic = false;
            for b in points {
                if b.tumor_baf > HIGHLY_ANEUPLOIDIC_REGION_CUTOFF {
                    has_high = true;
                } else if b.tumor_baf < 1.0 - HIGHLY_ANEUPLOIDIC_REGION_CUTOFF {
                    has_low = true;
                }
                if has_high && has_low {
                    aneuploidic = true;
                    break;
                }
            }
            if aneuploidic {
                found = true;
            }
        }
        found
    }

    pub fn calculate_baf_purity(&mut self) -> Option<f64> {
        if self.regions.is_empty() {
            return None;
        }
        self.regions.sort_by(|a, b| dcmp(a.observed_baf, b.observed_baf));
        let mut cumulative = 0;
        let max_threshold = self.total_baf_count as f64 * (1.0 - HIGHLY_ANEUPLOIDIC_RATIO_CUTOFF);
        let mut last_valid: Option<&ObservedRegion> = None;
        for r in &self.regions {
            cumulative += r.baf_count;
            let valid = r.baf_count as f64 >= HIGHLY_ANEUPLOIDIC_REGION_MIN_BAF_COUNT && r.observed_baf >= SOMATIC_FIT_ANEUPLOIDIC_REGION_CUTOFF;
            if (cumulative as f64) < max_threshold {
                if valid {
                    last_valid = Some(r);
                }
            } else if valid {
                return Some(r.observed_baf * 2.0 - 1.0);
            } else {
                return last_valid.map(|l| l.observed_baf * 2.0 - 1.0);
            }
        }
        None
    }
}

#[cfg(test)]
mod purity_adjuster_tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    /// buildPurityAdjuster: CobaltChromosomes of a typical sample of the gender
    fn adjuster(gender: Gender, purity: f64, norm_factor: f64) -> PurityAdjuster {
        let mut ratios = [1.0; 24];
        let (x, y) = if gender == Gender::Female { (1.0, 0.0) } else { (0.5, 0.5) };
        ratios[22] = x;
        ratios[23] = y;
        PurityAdjuster::with_ratios(ratios, purity, norm_factor)
    }

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= EPSILON, "{a} != {b}");
    }

    /// PurityAdjuster.impliedAverageRatio
    fn implied_average_ratio(purity: f64, norm_factor: f64, ploidy: f64) -> f64 {
        norm_factor * (1.0 - purity + ploidy / 2.0 * purity)
    }

    #[test]
    fn germline_copy_number() {
        let a = adjuster(Gender::Female, 1.0, 1.0);
        close(a.germline_copy_number("1"), 2.0);
        close(a.germline_copy_number("X"), 2.0);
        close(a.germline_copy_number("Y"), 0.0);
    }

    #[test]
    fn purity_adjusted_copy_number() {
        let a = adjuster(Gender::Male, 1.0, 1.0);
        close(a.purity_adjusted_copy_number_chr("1", 1.0), 2.0);
        close(a.purity_adjusted_copy_number_chr("1", 1.2), 2.4);
        close(a.purity_adjusted_copy_number_chr("1", 2.0), 4.0);
        let a = adjuster(Gender::Male, 0.8, 1.0);
        close(a.purity_adjusted_copy_number_chr("1", 1.0), 2.0);
        close(a.purity_adjusted_copy_number_chr("1", 1.5), 3.25);
        close(a.purity_adjusted_copy_number_chr("1", 2.0), 4.5);
        close(adjuster(Gender::Male, 0.85, 1.04).purity_adjusted_copy_number_chr("1", 0.0), 0.0);
        close(adjuster(Gender::Male, 0.85, 1.0).purity_adjusted_copy_number_chr("1", 0.575), 1.0);
    }

    #[test]
    fn purity_adjusted_copy_number_given_ratio() {
        let a = adjuster(Gender::Male, 1.0, 1.0);
        close(a.purity_adjusted_copy_number(1.0, 1.0), 2.0);
        close(a.purity_adjusted_copy_number(2.0, 1.0), 4.0);
        close(a.purity_adjusted_copy_number(1.0, 2.0), 2.0);
        let a = adjuster(Gender::Male, 1.0, 0.5);
        close(a.purity_adjusted_copy_number(1.0, 1.0), 4.0);
        close(a.purity_adjusted_copy_number(2.0, 1.0), 8.0);
        close(a.purity_adjusted_copy_number(1.0, 2.0), 4.0);
        let a = adjuster(Gender::Male, 0.5, 1.0);
        close(a.purity_adjusted_copy_number(1.0, 1.0), 2.0);
        close(a.purity_adjusted_copy_number(2.0, 1.0), 6.0);
        close(a.purity_adjusted_copy_number(1.0, 2.0), 0.0);
        close(a.purity_adjusted_copy_number(1.0, 4.0), -4.0);
    }

    fn frequency_matches_ploidy(a: &PurityAdjuster, normal_cn: i32, normal_ploidy: i32, tumor_cn: f64, tumor_ploidy: f64) {
        let f = a.expected_frequency(normal_cn as f64, normal_ploidy, tumor_cn, tumor_ploidy);
        close(a.purity_adjusted_ploidy(normal_cn as f64, normal_ploidy as f64, tumor_cn, f), tumor_ploidy);
    }

    #[test]
    fn purity_adjusted_frequency() {
        let a = adjuster(Gender::Female, 0.9, 1.0);
        frequency_matches_ploidy(&a, 2, 1, 3.0, 2.0);
        frequency_matches_ploidy(&a, 2, 0, 3.0, 2.0);
        let a = adjuster(Gender::Female, 1.0, 1.0);
        frequency_matches_ploidy(&a, 2, 1, 3.0, 2.0);
        frequency_matches_ploidy(&a, 2, 0, 3.0, 2.0);
    }

    #[test]
    fn expected_frequency_with_negative_tumor_copy_number() {
        close(adjuster(Gender::Female, 1.0, 1.0).expected_frequency(2.0, 1, -0.2, -0.1), 0.0);
    }

    #[test]
    fn norm_factor_ploidy_average_ratio_consistency() {
        let (average_ratio, norm_factor, purity) = (0.990374738586366, 0.8782, 0.08);
        let ploidy = PurityAdjuster::implied_ploidy(average_ratio, purity, norm_factor);
        close(PurityAdjuster::implied_norm_factor(average_ratio, purity, ploidy), norm_factor);
        close(implied_average_ratio(purity, norm_factor, ploidy), average_ratio);
    }
}

#[cfg(test)]
mod ambiguous_baf_tests {
    use super::*;
    use crate::test_support::observed_region;

    const EPSILON: f64 = 0.001;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= EPSILON, "{a} != {b}");
    }

    fn region(baf: f64, ratio: f64, count: i32) -> ObservedRegion {
        observed_region("1", 1001, 200, baf, ratio, GermlineStatus::Diploid, count)
    }

    fn around_peak(regions: &mut Vec<ObservedRegion>, peak: &ObservedRegion, count_diff: i32) {
        for i in 0..=10 {
            let baf = peak.observed_baf - (5 - i) as f64 * 0.001;
            let count = peak.baf_count - (5i32 - i).abs() * count_diff;
            regions.push(region(baf, 2.0, count));
        }
    }

    fn sorted(mut regions: Vec<ObservedRegion>) -> Vec<ObservedRegion> {
        regions.sort_by(|a, b| a.observed_baf.total_cmp(&b.observed_baf));
        regions
    }

    // AmbiguousBAFTest
    #[test]
    fn depth_based() {
        close(depth_based_baf(30), 0.581);
        close(depth_based_baf(60), 0.5564);
        close(depth_based_baf(90), 0.546);
        close(depth_based_baf(120), 0.5398);
    }

    #[test]
    fn ambiguous() {
        let (valid, invalid, n) = (1.0, 3.1, 100);
        let regions: Vec<ObservedRegion> =
            [(0.49, valid, n), (0.50, invalid, n), (0.51, valid, n * 2), (0.52, valid, n * 3), (0.53, valid, n * 4), (0.54, valid, n * 6), (0.55, valid, n * 4), (0.56, valid, n * 2), (0.57, valid, n), (0.58, valid, n)]
                .iter()
                .map(|&(b, r, c)| region(b, r, c))
                .collect();
        close(find_ambiguous_baf(&regions, 150), 0.5475);
        let regions: Vec<ObservedRegion> =
            [(0.52, valid, n - 1), (0.53, invalid, n), (0.54, valid, n), (0.55, valid, n - 1), (0.56, invalid, n), (0.57, valid, n - 1), (0.58, valid, n - 1)].iter().map(|&(b, r, c)| region(b, r, c)).collect();
        close(find_ambiguous_baf(&regions, 150), depth_based_baf(150));
        let mut regions: Vec<ObservedRegion> = [0.51, 0.515, 0.575, 0.58].iter().map(|&b| region(b, valid, n)).collect();
        around_peak(&mut regions, &region(0.53, valid, n * 5), 20);
        around_peak(&mut regions, &region(0.56, valid, n * 5), 20);
        let regions = sorted(regions);
        close(find_ambiguous_baf(&regions, 150), 0.535);
        close(find_ambiguous_baf(&regions, 40), depth_based_baf(40));
        let first = region(0.515, valid, n * 5);
        let second = region(0.56, valid, n * 5);
        let mut regions = vec![first.clone(), region(0.54, valid, n), region(0.545, valid, n), second.clone(), region(0.57, valid, n)];
        around_peak(&mut regions, &first, 20);
        around_peak(&mut regions, &second, 20);
        close(find_ambiguous_baf(&sorted(regions), 400), 0.525);
    }
}

#[cfg(test)]
mod region_fit_tests {
    use super::*;
    use crate::inputs::MedianRatio;
    use crate::test_support::{cobalt_chromosomes, default_fitting_config};

    // RegionFitCalculatorTest
    #[test]
    fn allowed_regions() {
        // RegionFitCalculator.isAllowedRegion: the chromosome is one of COBALT's
        assert!(cobalt_chromosomes(false).has_chromosome("Y"));
        assert!(!cobalt_chromosomes(true).has_chromosome("Y"));
        assert!(cobalt_chromosomes(false).has_chromosome("X") && cobalt_chromosomes(true).has_chromosome("X"));
    }

    #[test]
    fn baf_deviation_calcs() {
        let config = default_fitting_config();
        let cobalt = CobaltChromosomes::new(&[MedianRatio { chromosome: "1".into(), median_ratio: 0.5, count: 1 }], true);
        let calculator = RegionFitCalculator::new(&cobalt, &config, depth_based_baf(100));
        let adjuster = PurityAdjuster::new(0.08, 0.986, &cobalt);
        let close = |a: f64, b: f64| assert!((a - b).abs() <= 0.001, "{a} != {b}");
        close(calculator.baf_to_minimise_deviation(&adjuster, "1", 2.018, 0.53), 0.504);
        close(calculator.baf_to_minimise_deviation(&adjuster, "1", 1.5, 0.53), 0.667);
        close(calculator.baf_to_minimise_deviation(&adjuster, "1", 2.02, 0.51), 0.504);
        close(calculator.baf_to_minimise_deviation(&adjuster, "1", 1.9, 0.53), 0.526);
    }

    #[test]
    fn min_max_baf_estimates() {
        let expected_major = [1.0, 1.0, 1.0, 1.0, 1.1, 1.25, 1.5, 1.75, 2.0, 2.0, 2.0, 2.0, 2.0, 2.25, 2.5, 2.75];
        let cns = [1.25, 1.5, 1.75, 2.0, 2.1, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75];
        for max_baf in [1.0, 0.7] {
            for (i, (&cn, &major)) in cns.iter().zip(&expected_major).enumerate() {
                let major = if max_baf == 0.7 && i == 0 { -1.0 } else { major };
                let calc = estimate_min_max_baf(cn, BAF_PNT_5, max_baf);
                if calc == -1.0 {
                    assert_eq!(major, -1.0, "cn {cn} max {max_baf}");
                } else {
                    assert!((calc - major / cn).abs() <= 0.001, "cn {cn} max {max_baf}: {calc}");
                }
            }
        }
    }
}
