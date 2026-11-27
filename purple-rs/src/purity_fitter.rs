//! PurityPloidyFitter and VariantPurityFitter: the copy-number fit, the somatic fit where the sample is
//! highly diploid (or tumor-only / targeted), and the copy numbers of the final fit.

use crate::config::PurpleConfig;
use crate::constants::*;
use crate::copynumber::{validate_copy_numbers, CopyNumberFactory, PurpleCopyNumber};
use crate::fitting::*;
use crate::fittingsnv::{find_fitting_variants, SomaticPurityFitter};
use crate::inputs::{AmberData, CobaltData};
use crate::reference::ReferenceData;
use crate::region::ObservedRegion;
use crate::sv::StructuralVariant;
use crate::types::*;
use crate::variant::Variant;

#[derive(Clone, Debug)]
pub struct BestFit {
    pub fit: FittedPurity,
    pub score: FittedPurityScore,
    pub method: FittedPurityMethod,
    pub all_fits: Vec<FittedPurity>,
}

impl BestFit {
    /// PurityPloidyFitter.buildGermlineBestFit
    pub fn germline() -> BestFit {
        let fit = FittedPurity { purity: 1.0, norm_factor: 1.0, ploidy: 2.0, score: 0.0, diploid_proportion: 1.0, somatic_penalty: 0.0 };
        let score = FittedPurityScore {
            min_ploidy: 2.0,
            max_ploidy: 2.0,
            min_purity: 1.0,
            max_purity: 1.0,
            min_diploid_proportion: 1.0,
            max_diploid_proportion: 1.0,
        };
        BestFit { fit, score, method: FittedPurityMethod::Normal, all_fits: vec![fit] }
    }
}

pub struct FitResult {
    pub best_fit: BestFit,
    pub purity_adjuster: PurityAdjuster,
    pub copy_numbers: Vec<PurpleCopyNumber>,
    pub fitted_regions: Vec<ObservedRegion>,
}

/// VariantPurityFitter state
struct VariantPurityFitter<'a> {
    fitting_somatics: Vec<&'a Variant>,
    hotspot_svs: Vec<&'a StructuralVariant>,
    has_tumor: bool,
    fitter: SomaticPurityFitter,
}

impl<'a> VariantPurityFitter<'a> {
    fn new(config: &PurpleConfig, amber: &AmberData) -> VariantPurityFitter<'a> {
        VariantPurityFitter {
            fitting_somatics: Vec::new(),
            hotspot_svs: Vec::new(),
            has_tumor: false,
            fitter: SomaticPurityFitter {
                min_read_count: amber.min_somatic_total_read_count(),
                max_read_count: amber.max_somatic_total_read_count(),
                min_purity: config.fitting.min_purity,
                max_purity: config.fitting.max_purity,
                kd_min_peak: config.somatic_fitting.min_peak_variants,
                kd_min_somatics: config.somatic_fitting.min_total_variants,
            },
        }
    }

    fn set_state(&mut self, config: &PurpleConfig, somatics: &'a [Variant], svs: &'a [StructuralVariant], observed: &[ObservedRegion], immune: &[ChrRegion]) {
        self.fitting_somatics = find_fitting_variants(config.tumor_only_mode(), somatics, observed, immune);
        if !self.fitting_somatics.is_empty() {
            crate::debug!("somatic fitting variants({})", self.fitting_somatics.len());
        }
        let (mut sv_hotspots, mut sv_fragments) = (0, 0);
        for v in svs {
            if v.is_filtered() {
                continue;
            }
            if v.hotspot && v.end.is_some() {
                sv_hotspots += 1;
                self.hotspot_svs.push(v);
            }
            if v.end.is_some() {
                if let Some(f) = v.start.tumor_variant_fragment_count {
                    sv_fragments += f;
                }
            }
        }
        let (mut hotspots, mut allele_reads) = (0, 0);
        for v in &self.fitting_somatics {
            if !v.is_pass() {
                continue;
            }
            if v.is_hotspot_type() {
                hotspots += 1;
            }
            if v.vtype == VariantType::Snp {
                allele_reads += v.allele_read_count();
            }
        }
        if hotspots > 0 || allele_reads >= MIN_TOTAL_SOMATIC_VAR_ALLELE_READ_COUNT {
            crate::debug!("tumor evidence: somaticHotspotCount({hotspots}) alleleReadCountTotal({allele_reads})");
            self.has_tumor = true;
            return;
        }
        if sv_hotspots > 0 || sv_fragments >= MIN_TOTAL_SV_FRAGMENT_COUNT {
            crate::debug!("tumor evidence: svHotspotCount({sv_hotspots}) svFragmentReadCount({sv_fragments})");
            self.has_tumor = true;
            return;
        }
        let baf_total: i32 = observed
            .iter()
            .filter(|x| x.germline_status == GermlineStatus::Diploid)
            .filter(|x| x.observed_tumor_ratio < NO_TUMOR_DEPTH_RATIO_MIN || x.observed_tumor_ratio > NO_TUMOR_DEPTH_RATIO_MAX)
            .map(|x| x.baf_count)
            .sum();
        if baf_total >= NO_TUMOR_BAF_TOTAL {
            crate::debug!("tumor evidence: tumorEvidenceBafCountTotal({baf_total})");
            self.has_tumor = true;
            return;
        }
        self.has_tumor = false;
    }
}

/// VariantPurityFitter.somaticFitIsWorse
fn somatic_fit_is_worse(lowest: &FittedPurity, somatic: &FittedPurity, min_somatic_purity: f64) -> bool {
    less_than(lowest.purity, min_somatic_purity) && less_than(somatic.purity, min_somatic_purity) && greater_than(somatic.purity, lowest.purity)
}

fn in_range_of_lowest(lowest: f64, purities: &[FittedPurity]) -> Vec<FittedPurity> {
    purities
        .iter()
        .filter(|f| {
            let abs = (f.score - lowest).abs();
            let rel = (abs / lowest).abs();
            less_or_equal(abs, 0.0005) || less_or_equal(rel, 0.1)
        })
        .copied()
        .collect()
}

pub struct PurityPloidyFitter<'a> {
    pub config: &'a PurpleConfig,
    pub reference: &'a ReferenceData,
    pub amber: &'a AmberData,
    pub cobalt: &'a CobaltData,
    pub somatics: &'a [Variant],
    pub svs: &'a [StructuralVariant],
    pub calculator: &'a RegionFitCalculator<'a>,
    pub observed: &'a [ObservedRegion],
    pub gender: Gender,
}

impl PurityPloidyFitter<'_> {
    fn build_copy_numbers(&self, fit: &FittedPurity) -> (PurityAdjuster, Vec<ObservedRegion>, Vec<PurpleCopyNumber>, bool) {
        let adjuster = PurityAdjuster::new(fit.purity, fit.norm_factor, &self.cobalt.cobalt_chromosomes);
        let factory = CopyNumberFactory {
            min_tumor_ratio_count: self.config.fitting.min_diploid_tumor_ratio_count,
            min_tumor_ratio_count_at_centromere: self.config.fitting.min_diploid_tumor_ratio_count_at_centromere,
            targeted: self.config.target_regions_mode,
            average_read_depth: self.amber.average_tumor_depth,
            ploidy: fit.ploidy,
            purity_adjuster: &adjuster,
            cobalt: &self.cobalt.cobalt_chromosomes,
        };
        crate::debug!("building copy numbers");
        let fitted = self.calculator.fit_regions(fit.purity, fit.norm_factor, self.observed);
        let copy_numbers = factory.build(&fitted, self.svs);
        let valid = validate_copy_numbers(&copy_numbers);
        if !valid {
            crate::warn!("invalid copy numbers");
        }
        (adjuster, fitted, copy_numbers, valid)
    }

    /// PurityPloidyFitter.run; None where Java marks the fitter invalid
    pub fn run(&self) -> Option<FitResult> {
        let mut vpf = VariantPurityFitter::new(self.config, self.amber);
        vpf.set_state(self.config, self.somatics, self.svs, self.observed, &self.reference.special_regions.excluded_immune);

        // performCopyNumberFit
        let fitting_variants: Vec<&Variant> = if !self.config.tumor_only_mode() { vpf.fitting_somatics.clone() } else { Vec::new() };
        let factory = FittedPurityFactory::new(self.config, &self.cobalt.cobalt_chromosomes, self.calculator, self.observed, &fitting_variants);
        if !factory.valid_data_for_fit() {
            return None;
        }
        let mut candidates = factory.fit_purity();
        sort_by_score(&mut candidates);
        let copy_number_fit = *candidates.first()?;
        let fit_score = score(&in_range_of_lowest(copy_number_fit.score, &candidates));

        let (mut adjuster, mut fitted, mut copy_numbers, mut valid) = self.build_copy_numbers(&copy_number_fit);

        // performSomaticFit
        let diploid_candidates = most_diploid_per_purity(&candidates);
        let lowest_purity_fit = if !diploid_candidates.is_empty() {
            let mut m = diploid_candidates[0];
            for f in &diploid_candidates {
                if dcmp(f.purity, m.purity) == std::cmp::Ordering::Less {
                    m = *f;
                }
            }
            m
        } else {
            copy_number_fit
        };
        let highly_diploid = greater_or_equal(fit_score.max_diploid_proportion, self.config.somatic_fitting.highly_diploid_percentage);
        let targeted = self.reference.target_regions.has_target_regions();
        let (final_fit, method) = if self.config.tumor_only_mode() || targeted {
            let mut detector = AneuploidyDetector::new(self.observed, &self.amber.chromosome_bafs, &self.cobalt.cobalt_chromosomes);
            let has_aneuploidy = detector.has_aneuploidy();
            let diploid_high_purity = copy_number_fit.purity > SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN
                && copy_number_fit.ploidy > SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN
                && copy_number_fit.ploidy < SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX;
            let significant_contamination = self.amber.contamination > SOMATIC_FIT_CONTAMINATION_CUTOFF;
            if significant_contamination || (!has_aneuploidy && (diploid_high_purity || highly_diploid)) {
                match vpf.fitter.fit_from_somatics_only(&self.reference.driver_genes, &vpf.fitting_somatics, &candidates, self.config.tumor_only_mode()) {
                    Some(f) => {
                        crate::debug!("somatic fit found");
                        (f, FittedPurityMethod::Somatic)
                    }
                    None => {
                        crate::debug!("somatic fit not found, reverting to lowest purity fit");
                        (
                            FittedPurity {
                                purity: MIN_PURITY_DEFAULT,
                                norm_factor: lowest_purity_fit.norm_factor,
                                ploidy: DIPLOID_PLOIDY,
                                score: lowest_purity_fit.score,
                                diploid_proportion: lowest_purity_fit.diploid_proportion,
                                somatic_penalty: 0.0,
                            },
                            FittedPurityMethod::NoTumor,
                        )
                    }
                }
            } else {
                let mut baf_fit = None;
                if has_aneuploidy && diploid_high_purity && !detector.has_high_baf_region {
                    if let Some(p) = detector.calculate_baf_purity() {
                        baf_fit = find_matched_fitted_purity(p, &candidates);
                    }
                }
                (baf_fit.unwrap_or(copy_number_fit), FittedPurityMethod::Normal)
            }
        } else {
            let has_tumor = !highly_diploid || vpf.has_tumor;
            crate::info!(
                "maxDiploidProportion({}) diploidCandidates({}) purityRange({} - {}) hasTumor({})",
                hmf_common::fmt::java_fixed(fit_score.max_diploid_proportion, 3),
                diploid_candidates.len(),
                hmf_common::fmt::java_fixed(fit_score.min_purity, 3),
                hmf_common::fmt::java_fixed(fit_score.max_purity, 3),
                has_tumor
            );
            if !has_tumor {
                (lowest_purity_fit, FittedPurityMethod::NoTumor)
            } else if diploid_candidates.is_empty() {
                crate::warn!("unable to use somatic fit as there are no diploid candidates");
                (copy_number_fit, FittedPurityMethod::Normal)
            } else {
                let exceeds_spread = greater_or_equal(fit_score.purity_spread(), self.config.somatic_fitting.min_somatic_purity_spread);
                if !(exceeds_spread && highly_diploid) {
                    (copy_number_fit, FittedPurityMethod::Normal)
                } else {
                    match vpf.fitter.fit_from_somatics(&vpf.fitting_somatics, &vpf.hotspot_svs, &diploid_candidates, &copy_numbers, self.gender) {
                        None => (lowest_purity_fit, FittedPurityMethod::NoTumor),
                        Some(s) => {
                            if !somatic_fit_is_worse(&copy_number_fit, &s, self.config.somatic_fitting.min_somatic_purity) {
                                (s, FittedPurityMethod::Somatic)
                            } else {
                                (copy_number_fit, FittedPurityMethod::Normal)
                            }
                        }
                    }
                }
            }
        };
        let mut final_fit = final_fit;
        if method == FittedPurityMethod::NoTumor {
            // convertNoTumorFit
            let diploid100 = candidates.iter().find(|x| (x.ploidy - DIPLOID_PLOIDY).abs() < 0.001 && (x.purity - 1.0).abs() < 0.001).copied();
            final_fit = match diploid100 {
                Some(d) => FittedPurity { purity: 1.0, norm_factor: final_fit.norm_factor, ploidy: DIPLOID_PLOIDY, score: d.score, diploid_proportion: d.diploid_proportion, somatic_penalty: d.somatic_penalty },
                None => FittedPurity { purity: 1.0, ploidy: DIPLOID_PLOIDY, ..final_fit },
            };
        }
        if final_fit != copy_number_fit {
            let r = self.build_copy_numbers(&final_fit);
            adjuster = r.0;
            fitted = r.1;
            copy_numbers = r.2;
            valid = r.3;
        }
        if !valid {
            return None;
        }
        Some(FitResult { best_fit: BestFit { fit: final_fit, score: fit_score, method, all_fits: candidates }, purity_adjuster: adjuster, copy_numbers, fitted_regions: fitted })
    }
}
