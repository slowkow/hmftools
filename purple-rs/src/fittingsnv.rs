//! The fittingsnv package: SomaticPurityFitter (fitting variants, kernel-density peak purity, hotspot
//! purity, VAF-based purity of tumor-only and targeted samples), SomaticReadjustmentFit, SomaticDeviation,
//! WeightedPloidyHistogram and PeakModelFactory (the clonality model written to somatic.clonality.tsv).

use std::cell::RefCell;
use std::collections::HashMap;

use hmf_common::commons_math::{Binomial, Poisson};
use hmf_common::kernel_estimator::KernelEstimator;

use crate::constants::*;
use crate::copynumber::{chromosome_map, PurpleCopyNumber};
use crate::fitting::{find_matched_fitted_purity, FittedPurity, PurityAdjuster};
use crate::reference::DriverGeneResource;
use crate::region::ObservedRegion;
use crate::sv::StructuralVariant;
use crate::types::*;
use crate::variant::{AllelicDepth, Variant};

// ---------------------------------------------------------------------------------------------
// SomaticDeviation

thread_local! {
    /// SomaticDeviation's cache of maxConceivableReads by the binomial p (a pure function of p)
    static MAX_CONCEIVABLE: RefCell<HashMap<u64, i32>> = RefCell::new(HashMap::new());
}

fn max_conceivable_for_p(p: f64) -> i32 {
    MAX_CONCEIVABLE.with(|m| {
        *m.borrow_mut().entry(p.to_bits()).or_insert_with(|| Binomial::new(10_000, jmin(1.0, p)).inverse_cumulative_probability(0.999))
    })
}

/// SomaticDeviation.deviationFromMax(purityAdjuster, chromosome, major, tumorCopyNumber, variant)
pub fn deviation_from_max(adjuster: &PurityAdjuster, chromosome: &str, major_cn: f64, tumor_cn: f64, depth: &AllelicDepth) -> f64 {
    let normal_cn = adjuster.germline_copy_number(chromosome);
    deviation_from_max_with(adjuster, normal_cn, depth, jmax(0.0, tumor_cn), jmax(0.0, major_cn))
}

pub fn deviation_from_max_with(adjuster: &PurityAdjuster, normal_cn: f64, depth: &AllelicDepth, tumor_cn: f64, major: f64) -> f64 {
    let max_ploidy = max_conceivable_ploidy(adjuster, normal_cn, depth, tumor_cn, major);
    let somatic_ploidy = adjuster.purity_adjusted_ploidy(normal_cn, 0.0, tumor_cn, depth.allele_frequency());
    jmax(0.0, somatic_ploidy - max_ploidy)
}

pub fn max_conceivable_ploidy(adjuster: &PurityAdjuster, normal_cn: f64, depth: &AllelicDepth, tumor_cn: f64, major: f64) -> f64 {
    let reads = max_conceivable_reads(adjuster, normal_cn, depth, tumor_cn, major);
    let vaf = reads as f64 / depth.total_read_count as f64;
    adjuster.purity_adjusted_ploidy(normal_cn, 0.0, tumor_cn, vaf)
}

pub fn max_conceivable_reads(adjuster: &PurityAdjuster, normal_cn: f64, depth: &AllelicDepth, tumor_cn: f64, major: f64) -> i32 {
    let expected = adjuster.expected_frequency(normal_cn, 0, tumor_cn, major);
    let p = jround(expected * depth.total_read_count as f64 * 100.0) as f64 / 100.0 / 10_000.0;
    max_conceivable_for_p(p)
}

// ---------------------------------------------------------------------------------------------
// fitting variants

/// SomaticPurityFitter.findFittingVariants
pub fn find_fitting_variants<'a>(tumor_only: bool, variants: &'a [Variant], observed: &[ObservedRegion], immune: &[ChrRegion]) -> Vec<&'a Variant> {
    // GenomeRegionSelectorFactory.createImproved(fromRegions(observed)): the region containing the variant
    let mut by_chr: HashMap<String, Vec<&ObservedRegion>> = HashMap::new();
    for r in observed {
        if Chromosome::contains(&r.chromosome) {
            by_chr.entry(chr(&r.chromosome).short_name().to_string()).or_default().push(r);
        }
    }
    for list in by_chr.values_mut() {
        list.sort_by_key(|r| r.start);
    }
    let region_of = |v: &Variant| -> Option<&ObservedRegion> {
        let list = by_chr.get(chr(v.chromosome()).short_name())?;
        let i = list.partition_point(|r| r.start <= v.position());
        if i == 0 {
            return None;
        }
        let r = list[i - 1];
        if v.position() <= r.end { Some(r) } else { None }
    };
    let mut out = Vec::new();
    for v in variants {
        // CompoundFilter: SGT, human chromosome, NT
        let sgt_ok = match v.rec.info.get("SGT") {
            Some(val) => {
                let s = val.java_string();
                let parts: Vec<&str> = s.split("->").collect();
                !(parts.len() == 2 && parts[0] == parts[1])
            }
            None => true,
        };
        let nt_ok = match v.rec.info.get("NT") {
            None => true,
            Some(val) => val.java_string() == "ref",
        };
        if !v.is_pass() || !sgt_ok || !Chromosome::contains(v.chromosome()) || !nt_ok {
            continue;
        }
        if v.vtype == VariantType::Mnp {
            continue;
        } else if tumor_only && v.vtype == VariantType::Indel && v.repeat_count() > 0 {
            continue;
        }
        if !is_fitting_candidate(v, immune) {
            continue;
        }
        let status = region_of(v).map(|r| r.germline_status).unwrap_or(GermlineStatus::Unknown);
        if !v.is_hotspot_type() && status != GermlineStatus::Diploid {
            continue;
        }
        out.push(v);
    }
    out
}

fn is_fitting_candidate(v: &Variant, immune: &[ChrRegion]) -> bool {
    if immune.iter().any(|r| r.contains_position(v.chromosome(), v.position())) {
        return false;
    }
    if v.tumor_depth.is_none_or(|d| d.total_read_count == 0) {
        return false;
    }
    let gnomad = v.rec.info.as_double("GND_FREQ", -1.0);
    if gnomad > 0.0 {
        if v.is_hotspot() && gnomad < HOTSPOT_GNOMAD_FREQ_THRESHOLD {
            return true;
        }
        if v.tier == VariantTier::Panel && gnomad < PANEL_GNOMAD_FREQ_THRESHOLD {
            return true;
        }
    }
    if v.is_hotspot_type() {
        return true;
    }
    if gnomad > 0.0 {
        return false;
    }
    if matches!(v.tier, VariantTier::LowConfidence | VariantTier::Unknown) {
        return false;
    }
    if v.repeat_count() > SNV_FITTING_MAX_REPEATS {
        return false;
    }
    if v.rec.info.has("MAPPABILITY") && v.rec.info.as_double("MAPPABILITY", 0.0) < SNV_FITTING_MAPPABILITY {
        return false;
    }
    if v.reference_allele_read_count() > 0 {
        return false;
    }
    true
}

// ---------------------------------------------------------------------------------------------
// kernel density peaks

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SomaticPeak {
    pub allele_frequency: f64,
    pub count: i32,
}

/// SomaticKernelDensityPeaks.findPeaks
pub fn find_peaks(sample: &[f64]) -> Vec<SomaticPeak> {
    let mut est = KernelEstimator::new(0.001, 0.03);
    for x in sample {
        est.add_value(*x, 1.0);
    }
    let vafs: Vec<f64> = (0..=51).map(|x| x as f64 / 100.0).collect();
    let densities: Vec<f64> = vafs.iter().map(|v| est.probability(*v)).collect();
    let mut out = Vec::new();
    for i in 1..densities.len() - 1 {
        let d = densities[i];
        if greater_than(d, densities[i - 1]) && greater_than(d, densities[i + 1]) {
            let af = vafs[i];
            let count = sample.iter().filter(|vaf| greater_or_equal(af, *vaf - 0.015) && less_or_equal(af, *vaf + 0.015)).count() as i32;
            crate::debug!("discovered somatic peak: count({count}) alleleFrequency({})", hmf_common::fmt::java_fixed(af, 3));
            out.push(SomaticPeak { allele_frequency: af, count });
        }
    }
    out
}

/// SomaticKernelDensityPeaks.fitPurity
pub fn kd_fit_purity(candidates: &[FittedPurity], variants: &[&Variant], min_variants: i32, min_peak: f64, min_purity: f64, max_purity: f64) -> Option<FittedPurity> {
    if (variants.len() as i32) < min_variants {
        crate::info!("somatic variants count({}) too low for somatic fit", variants.len());
        return None;
    }
    let peaks = find_peaks(&variants.iter().map(|v| v.allele_frequency()).collect::<Vec<_>>());
    let in_range = |p: f64| greater_or_equal(p, min_purity) && less_or_equal(p, max_purity);
    let mut max_peak = 0;
    for peak in peaks.iter().rev() {
        let implied = peak.allele_frequency * 2.0;
        if in_range(implied) {
            if peak.count as f64 >= min_peak {
                match find_matched_fitted_purity(implied, candidates) {
                    Some(f) => {
                        crate::debug!("somatic implied purity({})", hmf_common::fmt::java_double_to_string(implied));
                        return Some(f);
                    }
                    None => crate::warn!("unable to find diploid solution for implied purity: {}", hmf_common::fmt::java_double_to_string(implied)),
                }
            }
            max_peak = max_peak.max(peak.count);
        }
    }
    if max_peak > 0 {
        for peak in peaks.iter().rev() {
            let implied = peak.allele_frequency * 2.0;
            if peak.count == max_peak {
                match find_matched_fitted_purity(implied, candidates) {
                    Some(f) => return Some(f),
                    None => crate::warn!("unable to find diploid solution for implied purity: {}", hmf_common::fmt::java_double_to_string(implied)),
                }
            }
        }
    }
    crate::debug!("unable to determine somatic implied purity.");
    None
}

// ---------------------------------------------------------------------------------------------
// SomaticPurityFitter

pub struct SomaticPurityFitter {
    pub min_read_count: i32,
    pub max_read_count: i32,
    pub min_purity: f64,
    pub max_purity: f64,
    pub kd_min_peak: i32,
    pub kd_min_somatics: i32,
}

/// SomaticPurityFitter.calc75thPercentileValue
pub fn calc_75th_percentile(vafs: &mut [f64]) -> f64 {
    if vafs.is_empty() {
        return 0.0;
    }
    if vafs.len() == 1 {
        return vafs[0];
    }
    vafs.sort_by(|a, b| dcmp(*a, *b));
    if vafs.len() <= 3 {
        let top = vafs.len() - 1;
        return if vafs[top] < SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF { vafs[top] } else { vafs[top.saturating_sub(1)] };
    }
    let idx = 0.75 * (vafs.len() - 1) as f64;
    let lower = (idx.floor() as usize).min(vafs.len() - 2);
    let fraction = idx - lower as f64;
    vafs[lower] + fraction * (vafs[lower + 1] - vafs[lower])
}

impl SomaticPurityFitter {
    pub fn fit_from_somatics(
        &self,
        somatics: &[&Variant],
        hotspot_svs: &[&StructuralVariant],
        diploid_candidates: &[FittedPurity],
        copy_numbers: &[PurpleCopyNumber],
        gender: Gender,
    ) -> Option<FittedPurity> {
        let filtered: Vec<&Variant> = somatics
            .iter()
            .copied()
            .filter(|x| x.is_hotspot_type() || (x.total_read_count() >= self.min_read_count && x.total_read_count() <= self.max_read_count))
            .collect();
        let mut peak_purity = 0.0;
        let mut fit: Option<FittedPurity> = None;
        if filtered.len() as i32 >= self.kd_min_somatics {
            crate::info!("looking for peak somatic allelic frequencies");
            if let Some(kd) = kd_fit_purity(diploid_candidates, &filtered, self.kd_min_somatics, self.kd_min_peak as f64, self.min_purity, self.max_purity) {
                peak_purity = kd.purity;
                fit = Some(kd);
                crate::info!("peak somatic purity({})", hmf_common::fmt::java_fixed(peak_purity, 3));
            }
        } else {
            crate::info!("somatic variants count({}) too low for somatic fit", filtered.len());
        }
        let hotspot_purity = self.find_hotspot_purity(&filtered, hotspot_svs, peak_purity, gender);
        if hotspot_purity > peak_purity {
            if let Some(m) = find_matched_fitted_purity(hotspot_purity, diploid_candidates) {
                return Some(m);
            }
        }
        if let Some(f) = fit {
            let reassessment = calc_readjustment_purity(&filtered, &f, copy_numbers);
            if reassessment > peak_purity {
                if let Some(m) = find_matched_fitted_purity(reassessment, diploid_candidates) {
                    return Some(m);
                }
            }
        }
        fit
    }

    pub fn fit_from_somatics_only(&self, driver_genes: &DriverGeneResource, variants: &[&Variant], candidates: &[FittedPurity], tumor_only: bool) -> Option<FittedPurity> {
        let mut vafs = Vec::new();
        for v in variants {
            if !v.is_hotspot_type() {
                let effect = v.impact.canonical_coding_effect;
                if !matches!(effect, CodingEffect::NonsenseOrFrameshift | CodingEffect::Missense | CodingEffect::Splice) {
                    continue;
                }
                let Some(g) = driver_genes.gene(&v.impact.gene_name) else { continue };
                if effect == CodingEffect::NonsenseOrFrameshift && !g.report_nonsense_and_frameshift {
                    continue;
                } else if effect == CodingEffect::Missense && !g.report_missense_and_inframe {
                    continue;
                } else if effect == CodingEffect::Splice && !g.report_splice {
                    continue;
                }
            }
            let mut vaf = v.allele_frequency();
            if vaf < SOMATIC_FIT_TUMOR_ONLY_VAF_MIN {
                continue;
            }
            if v.is_hotspot_type() {
                vaf = jmin(vaf, SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF);
            } else if tumor_only && vaf > SOMATIC_FIT_TUMOR_ONLY_VAF_MAX {
                continue;
            }
            vafs.push(vaf);
        }
        if vafs.is_empty() {
            return None;
        }
        let n = vafs.len();
        let purity = calc_75th_percentile(&mut vafs) * 2.0;
        crate::info!("somatic VAF-based purity({}) from {} variants", hmf_common::fmt::java_fixed(purity, 3), n);
        find_matched_fitted_purity(purity, candidates)
    }

    fn find_hotspot_purity(&self, somatics: &[&Variant], hotspot_svs: &[&StructuralVariant], peak_purity: f64, gender: Gender) -> f64 {
        if somatics.len() as i32 > SNV_HOTSPOT_MAX_SNV_COUNT {
            return 0.0;
        }
        let mut max_vaf = 0.0;
        for v in somatics {
            if !v.is_hotspot_type() {
                continue;
            }
            let c = chr(v.chromosome());
            if !is_autosome(c) && !(gender == Gender::Female && is_x(c)) {
                continue;
            }
            let af = v.allele_frequency();
            if af * 2.0 <= peak_purity || af > 0.5 || af < max_vaf {
                continue;
            }
            if !below_required_probability(peak_purity, v.total_read_count(), v.allele_read_count()) {
                continue;
            }
            crate::info!(
                "hotspot({}:{}) vaf({} {}/{})",
                v.chromosome(),
                v.position(),
                hmf_common::fmt::java_fixed(af, 3),
                v.allele_read_count(),
                v.total_read_count()
            );
            max_vaf = jmax(af, max_vaf);
        }
        for sv in hotspot_svs {
            if !sv.hotspot || sv.is_filtered() || sv.end.is_none() {
                continue;
            }
            let af = jmin(
                sv.start.allele_frequency.expect("hotspot SV without allele frequency"),
                sv.end.as_ref().unwrap().allele_frequency.expect("hotspot SV without allele frequency"),
            );
            if af < max_vaf {
                continue;
            }
            crate::info!(
                "hotspotSV({} {}:{}-{}:{}) vaf({})",
                sv.sv_type,
                sv.start.chromosome,
                sv.start.position,
                sv.end.as_ref().unwrap().chromosome,
                sv.end.as_ref().unwrap().position,
                hmf_common::fmt::java_fixed(af, 3)
            );
            max_vaf = af;
        }
        max_vaf * 2.0
    }
}

fn below_required_probability(peak_purity: f64, total: i32, allele: i32) -> bool {
    if peak_purity <= 0.0 {
        return true;
    }
    let expected = peak_purity * 0.5 * total as f64;
    let p = 1.0 - Poisson::new(expected).cumulative_probability(allele - 1);
    p < SNV_HOTSPOT_VAF_PROBABILITY
}

// ---------------------------------------------------------------------------------------------
// SomaticReadjustmentFit

fn vaf_probabilities(variants: &[&Variant], test_vaf: f64) -> Vec<f64> {
    let mut out: Vec<f64> = variants
        .iter()
        .map(|v| {
            let expected = v.total_read_count() as f64 * test_vaf;
            if v.allele_read_count() as f64 > expected {
                1.0 - Binomial::new(v.total_read_count(), test_vaf).cumulative_probability(v.allele_read_count() - 1)
            } else {
                1.0
            }
        })
        .collect();
    out.sort_by(|a, b| dcmp(*a, *b));
    out
}

fn probability_threshold(probabilities: &[f64], min: f64) -> f64 {
    if (probabilities.len() as i32) < SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS {
        return min;
    }
    jmax(probabilities[(SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS - 1) as usize], min)
}

/// SomaticReadjustmentFit.calcReadjustmentPurity
pub fn calc_readjustment_purity(variants: &[&Variant], standard: &FittedPurity, copy_numbers: &[PurpleCopyNumber]) -> f64 {
    let map = chromosome_map(copy_numbers);
    let mut candidates: Vec<&Variant> = Vec::new();
    let mut max_vaf = 0.0;
    for v in variants {
        if !is_autosome(chr(v.chromosome())) {
            continue;
        }
        let Some(list) = map.get(v.chromosome()) else { continue };
        let Some(cn) = list.iter().find(|x| position_within(v.position(), x.start, x.end)) else { continue };
        if !(cn.average_tumor_copy_number >= SNV_READJUST_CN_MIN
            && cn.average_tumor_copy_number <= SNV_READJUST_CN_MAX
            && cn.minor_allele_copy_number() >= SNV_READJUST_MINOR_ALLELE_MIN)
        {
            continue;
        }
        candidates.push(v);
        max_vaf = jmax(max_vaf, v.allele_frequency());
    }
    if (candidates.len() as i32) < SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS {
        return INVALID_PURITY;
    }
    let standard_purity = standard.purity;
    let standard_vaf = standard_purity * 0.5;
    let probs = vaf_probabilities(&candidates, standard_vaf);
    let threshold = probability_threshold(&probs, SNV_READJUST_PROB_THRESHOLD);
    let outlier_vafs = candidates.iter().filter(|x| x.allele_frequency() >= standard_vaf).count() as i32;
    let expected = outlier_vafs as f64 * threshold / 0.5;
    let outliers = probs.iter().filter(|x| **x >= threshold).count() as i32;
    let min_outliers = if standard_purity < SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY_LEVEL {
        SNV_READJUST_EXPECTED_VARIANT_COUNT
    } else {
        SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY
    };
    let required = jmax(expected * SNV_READJUST_EXPECTED_VARIANT_RATIO, min_outliers as f64);
    if (outliers as f64) < required {
        return INVALID_PURITY;
    }
    // findClosestExpectedPurity
    let max_implied = jmin(max_vaf * 2.0, 1.0);
    let mut min_valid = INVALID_PURITY;
    let mut test = standard_purity;
    while test < max_implied {
        let test_vaf = test * 0.5;
        let probs = vaf_probabilities(&candidates, test_vaf);
        let threshold = probability_threshold(&probs, SNV_READJUST_PROB_THRESHOLD);
        let outlier_vafs = candidates.iter().filter(|x| x.allele_frequency() >= test_vaf).count() as i32;
        let expected = outlier_vafs as f64 * threshold / 0.5;
        let outliers = probs.iter().filter(|x| **x <= threshold).count() as i32;
        if outliers as f64 <= expected && (min_valid == INVALID_PURITY || test < min_valid) {
            min_valid = test;
        }
        test += SNV_READJUST_PURITY_INCREMENT;
    }
    crate::debug!("reassessment purity({}) from candidate({})", hmf_common::fmt::java_fixed(min_valid, 4), candidates.len());
    min_valid
}

// ---------------------------------------------------------------------------------------------
// WeightedPloidyHistogram and PeakModelFactory

#[derive(Clone, Debug)]
pub struct WeightedPloidy {
    pub depth: AllelicDepth,
    pub ploidy: f64,
    pub weight: f64,
}

pub struct WeightedPloidyHistogram {
    pub max_ploidy: f64,
    pub bin_width: f64,
    pub offset: f64,
}

impl WeightedPloidyHistogram {
    pub fn new(max_ploidy: f64, bin_width: f64, offset: f64) -> WeightedPloidyHistogram {
        WeightedPloidyHistogram { max_ploidy, bin_width, offset }
    }
    pub fn bucket(&self, ploidy: f64) -> i32 {
        jround((ploidy - self.offset) / self.bin_width) as i32
    }
    pub fn ploidy(&self, bucket: i32) -> f64 {
        bucket as f64 * self.bin_width + self.offset
    }
    pub fn histogram(&self, ploidies: &[WeightedPloidy]) -> Vec<f64> {
        let max_bucket = self.bucket(self.max_ploidy);
        let mut result = vec![0.0; (max_bucket + 1).max(0) as usize];
        for p in ploidies {
            let b = self.bucket(p.ploidy);
            if b <= max_bucket {
                result[b as usize] += p.weight;
            }
        }
        result
    }
    pub fn model_histogram(&self, peaks: &[&PeakModelData]) -> Vec<f64> {
        let max_bucket = self.bucket(self.max_ploidy);
        let mut result = vec![0.0; (max_bucket + 1).max(0) as usize];
        for p in peaks {
            let b = self.bucket(p.bucket);
            if b <= max_bucket {
                result[b as usize] += p.bucket_weight;
            }
        }
        result
    }
    pub fn peak_ploidy(&self, peak_bin_count: i32, ploidies: &[WeightedPloidy]) -> f64 {
        self.ploidy(peak_bucket(peak_bin_count, &self.histogram(ploidies)))
    }
}

/// WeightedPloidyHistogram.peakBucket
pub fn peak_bucket(buffer: i32, histogram: &[f64]) -> i32 {
    let mut max_bucket = -1;
    let mut max_weight = 0.0;
    for i in 0..histogram.len() as i32 {
        if greater_than(histogram[i as usize], 0.0) {
            let mut w = 0.0;
            let lo = (i - buffer).max(0);
            let hi = (i + buffer).min(histogram.len() as i32 - 1);
            for j in lo..=hi {
                w += histogram[j as usize];
            }
            if greater_than(w, max_weight) {
                max_weight = w;
                max_bucket = i;
            }
        }
    }
    max_bucket
}

#[derive(Clone, Debug, PartialEq)]
pub struct PeakModelData {
    pub peak: f64,
    pub peak_avg_weight: f64,
    pub bucket: f64,
    pub bucket_weight: f64,
    pub is_valid: bool,
    pub is_subclonal: bool,
}

pub struct PeakModelFactory {
    max_ploidy: f64,
    model_width: f64,
    precise: WeightedPloidyHistogram,
    binomials: RefCell<HashMap<(i32, i32), Binomial>>,
}

impl PeakModelFactory {
    pub fn new(max_ploidy: f64, model_width: f64) -> PeakModelFactory {
        PeakModelFactory { max_ploidy, model_width, precise: WeightedPloidyHistogram::new(max_ploidy, PEAK_BIN_WIDTH, 0.0), binomials: RefCell::new(HashMap::new()) }
    }

    pub fn model(&self, ploidies: &mut [WeightedPloidy]) -> Vec<PeakModelData> {
        let mut has_valid_subclonal = false;
        let residual = WeightedPloidyHistogram::new(0.85, self.model_width, 0.0);
        let residual_actual = residual.histogram(ploidies);
        let mut model: Vec<PeakModelData> = Vec::new();
        let initial = positive_weight(ploidies);
        for _ in 0..10 {
            let peak = self.precise.peak_ploidy(PEAK_BIN_COUNT, ploidies);
            let offset = self.offset(peak);
            let hist = WeightedPloidyHistogram::new(self.max_ploidy, self.model_width, offset);
            let peak_ploidies: Vec<WeightedPloidy> = ploidies
                .iter()
                .filter(|x| greater_than(x.ploidy, peak - self.model_width / 2.0) && less_than(x.ploidy, peak + self.model_width / 2.0))
                .cloned()
                .collect();
            let avg = if peak_ploidies.is_empty() {
                0.0
            } else {
                hmf_common::java_math::double_stream_sum(peak_ploidies.iter().map(|x| x.weight)) / peak_ploidies.len() as f64
            };
            let peak_hist = self.model_peak_histogram(peak, &peak_ploidies);
            let current = hist.histogram(ploidies);
            for p in ploidies.iter_mut() {
                let b = hist.bucket(p.ploidy) as usize;
                let bucket_weight = current[b];
                let peak_weight = peak_hist[b];
                p.weight = if is_zero(bucket_weight) { 0.0 } else { p.weight - (peak_weight / bucket_weight).abs() };
            }
            let valid = greater_or_equal(avg, PEAK_BIN_MIN_AVERAGE_WEIGHT) && greater_than(peak, 0.0);
            let subclonal = less_than(peak, PEAK_BIN_CLONAL_PLOIDY);
            has_valid_subclonal |= subclonal && valid;
            for (b, w) in peak_hist.iter().enumerate() {
                model.push(PeakModelData { peak, peak_avg_weight: avg, bucket: b as f64 * self.model_width, bucket_weight: *w, is_valid: valid, is_subclonal: subclonal });
            }
            let remaining = positive_weight(ploidies);
            if less_than(remaining / initial, MAX_UNEXPLAINED_WEIGHT_PERCENT) {
                break;
            }
        }
        if has_valid_subclonal {
            return model;
        }
        let valid: Vec<&PeakModelData> = model.iter().filter(|x| x.is_valid).collect();
        let residual_model = residual.model_histogram(&valid);
        let mut extra = Vec::new();
        for (i, actual) in residual_actual.iter().enumerate() {
            let pct = if is_zero(*actual) { 1.0 } else { (actual - residual_model[i]) / actual };
            if greater_than(pct, 0.0) {
                extra.push(PeakModelData { peak: 0.0, peak_avg_weight: 1.0, bucket: i as f64 * self.model_width, bucket_weight: pct, is_valid: true, is_subclonal: true });
            }
        }
        model.extend(extra);
        model
    }

    fn offset(&self, peak: f64) -> f64 {
        peak - jround(peak / self.model_width) as f64 * self.model_width
    }

    fn bucket(&self, ploidy: f64) -> i32 {
        jround(ploidy / self.model_width) as i32
    }

    pub fn model_peak_histogram(&self, peak: f64, peak_ploidies: &[WeightedPloidy]) -> Vec<f64> {
        let offset = self.offset(peak);
        let max_bucket = self.bucket(self.max_ploidy);
        let mut result = vec![0.0; (max_bucket + 1) as usize];
        let scaling: Vec<f64> = peak_ploidies.iter().map(|p| p.weight / self.ploidy_likelihood(peak, p)).collect();
        let start = self.bucket(peak - offset);
        let likelihood = |ploidy: f64| -> f64 {
            let mut r = 0.0;
            for (i, p) in peak_ploidies.iter().enumerate() {
                r += scaling[i] * self.ploidy_likelihood(ploidy, p);
            }
            r
        };
        let mut i = start;
        while i <= max_bucket {
            let l = likelihood(i as f64 * self.model_width + offset);
            result[i as usize] = l;
            if is_zero(l) {
                break;
            }
            i += 1;
        }
        let mut i = start - 1;
        while i >= 0 {
            let l = likelihood(i as f64 * self.model_width + offset);
            result[i as usize] = l;
            if is_zero(l) {
                break;
            }
            i -= 1;
        }
        result
    }

    pub fn ploidy_likelihood(&self, ploidy: f64, w: &WeightedPloidy) -> f64 {
        let key = (w.depth.allele_read_count, w.depth.total_read_count);
        let d = *self.binomials.borrow_mut().entry(key).or_insert_with(|| Binomial::new(w.depth.total_read_count, w.depth.allele_frequency()));
        let lower = jmax(0.0, ploidy - self.model_width / 2.0) / w.ploidy * w.depth.allele_read_count as f64;
        let lower_r = jround(lower) as i32;
        let lower_add = lower_r as f64 + 0.5 - lower;
        let upper = jmax(0.0, ploidy + self.model_width / 2.0) / w.ploidy * w.depth.allele_read_count as f64;
        let upper_r = jround(upper) as i32;
        let upper_sub = upper_r as f64 + 0.5 - upper;
        let raw = d.cumulative_probability(upper_r) - d.cumulative_probability(lower_r) + lower_add * d.probability(lower_r) - upper_sub * d.probability(upper_r);
        jround(raw * 100.0) as f64 / 100.0
    }
}

fn positive_weight(ploidies: &[WeightedPloidy]) -> f64 {
    hmf_common::java_math::double_stream_sum(ploidies.iter().map(|x| jmax(0.0, x.weight)))
}

/// SomaticPeakStream.generateModelPeaks
pub fn generate_model_peaks(variants: &[Variant]) -> Vec<PeakModelData> {
    let mut ploidies: Vec<WeightedPloidy> = Vec::new();
    for v in variants {
        if !v.is_pass() || v.copy_number() >= CLONALITY_MAX_PLOIDY {
            continue;
        }
        if !Chromosome::contains(v.chromosome()) || !is_autosome(chr(v.chromosome())) {
            continue;
        }
        if let Some(d) = v.tumor_depth {
            ploidies.push(WeightedPloidy { depth: d, ploidy: v.copy_number(), weight: 1.0 });
        }
    }
    crate::debug!("somatic peak uses {} variants", ploidies.len());
    PeakModelFactory::new(CLONALITY_MAX_PLOIDY, CLONALITY_BIN_WIDTH).model(&mut ploidies)
}

/// PeakModelFile.write
pub fn write_peak_model(path: &str, model: &[PeakModelData]) -> anyhow::Result<()> {
    let f = |v: f64| hmf_common::fmt::java_fixed(v, 4);
    let mut out = String::from("peak\tbucket\tbucketWeight\tpeakAvgWeight\tisValid\tisSubclonal\n");
    for m in model {
        out.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", f(m.peak), f(m.bucket), f(m.bucket_weight), f(m.peak_avg_weight), m.is_valid, m.is_subclonal));
    }
    std::fs::write(path, out).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}

/// SubclonalLikelihood
pub struct SubclonalLikelihood {
    histogram: WeightedPloidyHistogram,
    likelihood: Vec<f64>,
}

impl SubclonalLikelihood {
    pub fn new(bin_width: f64, model: &[PeakModelData]) -> SubclonalLikelihood {
        let histogram = WeightedPloidyHistogram::new(2.0, bin_width, 0.0);
        let valid: Vec<&PeakModelData> = model.iter().filter(|x| x.is_valid).collect();
        let clonal = histogram.model_histogram(&valid.iter().copied().filter(|x| !x.is_subclonal).collect::<Vec<_>>());
        let subclonal = histogram.model_histogram(&valid.iter().copied().filter(|x| x.is_subclonal).collect::<Vec<_>>());
        let mut likelihood = vec![0.0; clonal.len()];
        for i in 0..clonal.len() {
            let total = clonal[i] + subclonal[i];
            if greater_than(total, 0.0) {
                likelihood[i] = subclonal[i] / (clonal[i] + subclonal[i]);
            }
        }
        SubclonalLikelihood { histogram, likelihood }
    }
    pub fn subclonal_likelihood(&self, ploidy: f64) -> f64 {
        let b = self.histogram.bucket(ploidy);
        if b >= 0 && (b as usize) < self.likelihood.len() {
            self.likelihood[b as usize]
        } else if b < 0 {
            panic!("ArrayIndexOutOfBoundsException: Index {b} out of bounds for length {}", self.likelihood.len())
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod weighted_ploidy_tests {
    use super::*;

    fn wp(ploidy: f64, alt: i32, total: i32) -> WeightedPloidy {
        WeightedPloidy { depth: AllelicDepth { total_read_count: total, allele_read_count: alt }, ploidy, weight: 1.0 }
    }

    // WeightedPloidyHistogramTest
    #[test]
    fn max_bucket() {
        WeightedPloidyHistogram::new(5.0, 0.01, 0.0).histogram(&[wp(5.0, 20, 40)]);
    }

    #[test]
    fn bucket_with_offset() {
        let neg = WeightedPloidyHistogram::new(10.0, 0.05, -0.02);
        assert_eq!([0.05, 0.06, 0.08, 0.10, 0.11].map(|p| neg.bucket(p)), [1, 2, 2, 2, 3]);
        let pos = WeightedPloidyHistogram::new(10.0, 0.05, 0.02);
        assert_eq!([0.04, 0.05, 0.07, 0.09, 0.10].map(|p| pos.bucket(p)), [0, 1, 1, 1, 2]);
    }

    #[test]
    fn bucket_no_offset() {
        let h = WeightedPloidyHistogram::new(10.0, 0.01, 0.0);
        assert_eq!([0.00499, 0.005, 0.01499, 0.015].map(|p| h.bucket(p)), [0, 1, 1, 2]);
        for (b, p) in [(0, 0.0), (1, 0.01), (2, 0.02)] {
            assert!((h.ploidy(b) - p).abs() <= 0.01);
        }
    }

    #[test]
    fn histogram_construction() {
        // ploidies.tsv of the Java test resources: ploidy, allele reads, total reads
        let ploidies: Vec<WeightedPloidy> = include_str!("testdata/ploidies.tsv")
            .lines()
            .map(|l| {
                let v: Vec<&str> = l.split('\t').collect();
                wp(v[0].parse().unwrap(), v[1].parse().unwrap(), v[2].parse().unwrap())
            })
            .filter(|w| less_than(w.ploidy, 10.0))
            .collect();
        let close = |a: f64, b: f64| assert!((a - b).abs() <= 0.1, "{a} != {b}");
        let h = WeightedPloidyHistogram::new(4.0, 0.02, 0.0);
        let hist = h.histogram(&ploidies);
        for (p, e) in [(0.5, 3.0), (1.0, 24.0), (2.0, 8.0), (3.0, 5.0)] {
            close(hist[h.bucket(p) as usize], e);
        }
        let h = WeightedPloidyHistogram::new(4.0, 0.05, 0.0);
        let hist = h.histogram(&ploidies);
        for (p, e) in [(1.0, 65.0), (2.0, 20.0), (3.0, 8.0), (4.0, 0.0)] {
            close(hist[h.bucket(p) as usize], e);
        }
    }

    #[test]
    fn histogram_peak() {
        assert_eq!(peak_bucket(1, &[10.0, 11.0, 1.0, 18.0, 8.0]), 2);
        assert_eq!(peak_bucket(1, &[10.0, 11.0, 0.0, 18.0, 8.0]), 3);
        assert_eq!(peak_bucket(2, &[10.0, 11.0, 0.0, 18.0, 8.0]), 1);
    }
}
