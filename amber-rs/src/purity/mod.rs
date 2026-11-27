//! Tumor-only purity analysis (com.hartwig.hmftools.amber.purity): the noise floor below which tumor-only
//! sites are dropped, and the contamination peaks.
//!
//! The evidence outside the HLA, IG and TR regions is searched for peaks of symmetric VAF over a grid of
//! levels (peak.rs). A peak whose captured points have the gnomAD frequencies of the heterozygous baseline
//! (gnomad.rs) is a contamination peak when 25-75% of its points are in its homozygous band and the points
//! spread over chromosome arms and SNV types as the baseline does (auc.rs); it is a copy number peak when
//! its points are concentrated in a few arms. The cutoff is a third of the lowest copy number peak, at most
//! 0.04.
//!
//! Java faults found here are AM53 to AM55 (REPORT.md); the LOH ones, AM51 and AM52, are in crate::loh.

#![allow(dead_code)]

pub mod auc;
pub mod binomial;
#[cfg(test)]
mod driver;
pub mod genome;
pub mod gnomad;
pub mod hash_order;
pub mod peak;
pub mod regions_filter;

use crate::types::{AmberSite, PositionEvidence};
use auc::{AucCalculator, ChrArm, CanonicalSnvType};
use hmf_common::fmt::java_fixed;
pub use peak::CandidatePeak;
use peak::{EvidenceKeys, LevelResult, PeakSearch, PeakStats};

pub const HET_VAF_LOWER_BOUND: f64 = 0.35;
pub const HET_VAF_UPPER_BOUND: f64 = 0.65;
pub const MIN_CUTOFF: f64 = 0.04;
pub const HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION: f64 = 0.25;
pub const HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION: f64 = 0.75;
pub const MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION: f64 = 0.825;
pub const MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION_PEAK_FACTOR: f64 = 100.0;
pub const CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION: f64 = 0.4;
pub const CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS: f64 = 0.1;

/// PurityAnalysisConfig
pub struct PurityAnalysisConfig {
    pub tumor_id: String,
    pub v38: bool,
    pub output_dir: String,
    pub threads: usize,
}

/// What a peak was found to be
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeakClass {
    /// failed the gnomAD frequency check
    FailedGnomad,
    Contamination,
    CopyNumber,
    /// neither
    Other,
}

/// The steps of the analysis, for the comparison with Java (purity/driver.rs)
#[derive(Clone, Debug, Default)]
pub struct PurityReport {
    pub evidence: usize,
    pub filtered: usize,
    pub het_variants: usize,
    pub baseline_frequency: f64,
    pub duplicates: bool,
    pub levels: Vec<LevelResult>,
    /// every maximum of the peak search, in order, with its class
    pub peaks: Vec<(CandidatePeak, PeakClass)>,
}

/// TumorOnlyPurityAnalysis
pub struct TumorOnlyPurityAnalysis {
    cutoff: f64,
    contamination_peaks: Vec<CandidatePeak>,
    copy_number_peaks: Vec<CandidatePeak>,
    report: PurityReport,
}

impl TumorOnlyPurityAnalysis {
    /// `evidence`: the tumor evidence that passed the depth, quality and support filters, without the
    /// LOH-recovered sites, in genome order. `sites`: every loci-file site, in file order (Java: the sites
    /// AMBER loaded, after any -specific_regions / target-region filtering).
    pub fn new(evidence: &[PositionEvidence], sites: &[AmberSite], config: &PurityAnalysisConfig) -> TumorOnlyPurityAnalysis {
        let refs: Vec<&PositionEvidence> = evidence.iter().collect();
        TumorOnlyPurityAnalysis::from_refs(&refs, sites, config)
    }

    /// as `new`, for evidence the caller holds elsewhere (no copy of the points is needed)
    pub fn from_refs(evidence: &[&PositionEvidence], sites: &[AmberSite], config: &PurityAnalysisConfig) -> TumorOnlyPurityAnalysis {
        let v38 = config.v38;
        crate::debug!("evaluating noise floor with {} evidence points", evidence.len());

        // filterOutExcludedRegions: HLA, IG and TR regions
        let mut immune = genome::hla_regions(v38);
        immune.extend(genome::ig_regions(v38));
        immune.extend(genome::tr_regions(v38));
        let kept = regions_filter::RegionsFilter::new(&immune).filter(evidence.iter().map(|e| (e.chromosome.as_str(), e.position)));
        let points: Vec<&PositionEvidence> = kept.iter().map(|&i| evidence[i as usize]).collect();
        drop(kept);
        crate::debug!("immune region data removed leaving {} evidence points", points.len());

        // getBaselineHetVariants
        let het: Vec<&PositionEvidence> = points
            .iter()
            .copied()
            .filter(|e| {
                let vaf = e.symmetric_vaf();
                vaf > HET_VAF_LOWER_BOUND && vaf < HET_VAF_UPPER_BOUND
            })
            .collect();
        let supplier = gnomad::DefaultGnomadFrequencySupplier::new(sites, v38);
        let baseline_frequency = gnomad::mean_frequency(&supplier, het.iter().copied());
        crate::debug!("baseline het gnomad frequency: {}", java_fixed(baseline_frequency, 3));

        let chr_arm_auc = AucCalculator::<ChrArm>::new(het.iter().map(|e| auc::chr_arm(v38, &e.chromosome, e.position)));
        let snv_type_auc = AucCalculator::<CanonicalSnvType>::new(het.iter().map(|e| auc::snv_type(e)));
        let het_count = het.len();
        drop(het);

        let keys = EvidenceKeys::new(&points);
        let search = PeakSearch::new(&points, &keys, config.threads);

        let mut report = PurityReport {
            evidence: evidence.len(),
            filtered: points.len(),
            het_variants: het_count,
            baseline_frequency,
            duplicates: keys.has_duplicates(),
            levels: search.levels,
            peaks: Vec::new(),
        };
        let mut contamination_peaks = Vec::new();
        let mut copy_number_peaks = Vec::new();
        let java_order = crate::compat::java("AM54");

        for mut peak in search.peaks {
            let captured = peak.all_captured_points(&keys);
            // the order Java iterates allCapturedPoints(), a HashSet<PositionEvidence> (AM54)
            let ordered: Vec<u32> = if java_order { java_hash_set_order(&captured, &points) } else {
                let mut sorted = captured.clone();
                sorted.sort_unstable();
                sorted
            };
            let gnomad = gnomad::check_gnomad_frequencies(&supplier, ordered.iter().map(|&i| points[i as usize]), baseline_frequency);
            let mut stats: PeakStats = peak.stats.unwrap_or_default();
            stats.homozygous_count = peak.homozygous_points().len();
            stats.heterozygous_count = peak.heterozygous_points().len();
            stats.captured = captured.len();
            stats.gnomad = gnomad;
            if !gnomad.pass {
                crate::debug!("peak {} failed Gnomad frequency check", java_fixed(peak.vaf(), 3));
                peak.stats = Some(stats);
                peak.drop_bands();
                report.peaks.push((peak, PeakClass::FailedGnomad));
                continue;
            }

            let homozygous_proportion = peak.homozygous_proportion(&keys);
            let (chr_arm, arm_overflow) = chr_arm_auc.integral(captured.iter().map(|&i| auc::chr_arm(v38, &points[i as usize].chromosome, points[i as usize].position)));
            let (mutation, mutation_overflow) = snv_type_auc.integral(captured.iter().map(|&i| auc::snv_type(points[i as usize])));
            let peak_points = captured.len();
            let chr_arm_lower_bound =
                CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION * (peak_points as f64 / MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION_PEAK_FACTOR).min(1.0).sqrt();
            crate::debug!(
                "peak {}, homProportion: {}, chrArmAUC: {}, mutationAUC: {}",
                java_fixed(peak.vaf(), 3),
                java_fixed(homozygous_proportion, 3),
                java_fixed(chr_arm, 3),
                java_fixed(mutation, 3)
            );
            stats.homozygous_proportion = homozygous_proportion;
            stats.chr_arm_auc = chr_arm;
            stats.mutation_auc = mutation;
            stats.chr_arm_lower_bound = chr_arm_lower_bound;
            stats.auc_overflow = arm_overflow || mutation_overflow;
            peak.stats = Some(stats);
            peak.drop_bands();

            // HOM_PROP_RANGE.contains: the open range (0.25, 0.75)
            let in_hom_range = homozygous_proportion > HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION
                && homozygous_proportion < HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION;
            if in_hom_range && mutation > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION && chr_arm > chr_arm_lower_bound {
                crate::debug!("peak {} is a contamination peak", java_fixed(peak.vaf(), 3));
                report.peaks.push((peak.clone(), PeakClass::Contamination));
                contamination_peaks.push(peak);
                continue;
            }
            if chr_arm < CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS || mutation > MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION && chr_arm < chr_arm_lower_bound {
                crate::debug!("peak {} is a copy number peak", java_fixed(peak.vaf(), 3));
                report.peaks.push((peak.clone(), PeakClass::CopyNumber));
                copy_number_peaks.push(peak);
            } else {
                report.peaks.push((peak, PeakClass::Other));
            }
        }

        // cutoff(): a third of the lowest copy number peak, at most MIN_CUTOFF
        let min_copy_number_peak = copy_number_peaks.iter().map(|p| p.vaf()).min_by(|a, b| a.total_cmp(b)).unwrap_or(0.0);
        // Doubles.isZero: within 1e-10
        let cutoff = if min_copy_number_peak.abs() < 1e-10 { MIN_CUTOFF } else { MIN_CUTOFF.min(min_copy_number_peak / 3.0) };

        TumorOnlyPurityAnalysis { cutoff, contamination_peaks, copy_number_peaks, report }
    }
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }
    pub fn contamination_peaks(&self) -> &[CandidatePeak] {
        &self.contamination_peaks
    }
    pub fn copy_number_peaks(&self) -> &[CandidatePeak] {
        &self.copy_number_peaks
    }
    /// AmberApplication.runTumorOnly's contamination: the highest contamination peak, 0 with none
    pub fn contamination(&self) -> f64 {
        self.contamination_peaks.iter().map(|p| p.vaf()).max_by(|a, b| a.total_cmp(b)).unwrap_or(0.0)
    }
    pub fn report(&self) -> &PurityReport {
        &self.report
    }
}

/// PositionEvidence.hashCode on a JVM with constant identity hashes (-XX:hashCode=2), where the AmberBase
/// enum constants hash to 1: Objects.hash(Chromosome, Position, Ref, Alt)
pub fn position_evidence_hash(e: &PositionEvidence) -> i32 {
    let mut h: i32 = 1;
    h = h.wrapping_mul(31).wrapping_add(hmf_common::java_hash::java_string_hash(&e.chromosome));
    h = h.wrapping_mul(31).wrapping_add(e.position);
    h = h.wrapping_mul(31).wrapping_add(1);
    h.wrapping_mul(31).wrapping_add(1)
}

/// The iteration order of Java's `new HashSet<>()` filled with `captured` (distinct points, in insertion order)
fn java_hash_set_order(captured: &[u32], points: &[&PositionEvidence]) -> Vec<u32> {
    let mut set = hash_order::JavaHashSetOrder::new();
    let mut last_chromosome: Option<(&str, i32)> = None;
    for &i in captured {
        let e = points[i as usize];
        // the chromosome's String hash, cached across a run of one chromosome
        let chr_hash = match last_chromosome {
            Some((c, h)) if c == e.chromosome => h,
            _ => {
                let h = hmf_common::java_hash::java_string_hash(&e.chromosome);
                last_chromosome = Some((&e.chromosome, h));
                h
            }
        };
        let mut h: i32 = 1;
        h = h.wrapping_mul(31).wrapping_add(chr_hash);
        h = h.wrapping_mul(31).wrapping_add(e.position);
        h = h.wrapping_mul(31).wrapping_add(1);
        h = h.wrapping_mul(31).wrapping_add(1);
        set.insert(i, h);
    }
    set.iteration_order()
}

#[cfg(test)]
mod tests;
