//! The copynumber package: PurpleCopyNumber, CombinedRegion and the extensions that turn fitted regions
//! into copy-number segments (ExtendDiploid, ExtendNonDiploid, SV-implied copy numbers, ExtendLongArm,
//! populateUnknown, ExtendDiploidBAF), the SV leg ploidies and copy-number changes
//! (copynumber.sv, CopyNumberEnrichedStructuralVariantFactory), ChromosomeCopyNumbers and LohCalcs.

use std::collections::{BTreeMap, HashMap};

use crate::constants::*;
use crate::fitting::PurityAdjuster;
use crate::inputs::CobaltChromosomes;
use crate::region::ObservedRegion;
use crate::sv::{StructuralVariant, SvLeg};
use crate::types::*;

#[derive(Clone, Debug, PartialEq)]
pub struct PurpleCopyNumber {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub baf_count: i32,
    pub average_actual_baf: f64,
    pub average_observed_baf: f64,
    pub average_tumor_copy_number: f64,
    pub depth_window_count: i32,
    pub segment_start_support: SegmentSupport,
    pub segment_end_support: SegmentSupport,
    pub method: CopyNumberMethod,
    pub gc_content: f64,
    pub min_start: i32,
    pub max_start: i32,
}

impl PurpleCopyNumber {
    pub fn minor_allele_copy_number(&self) -> f64 {
        if less_than(self.average_actual_baf, 0.50) { 0.0 } else { jmax(0.0, (1.0 - self.average_actual_baf) * self.average_tumor_copy_number) }
    }
    pub fn major_allele_copy_number(&self) -> f64 {
        self.average_tumor_copy_number - self.minor_allele_copy_number()
    }
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
    pub fn sv_support(&self) -> bool {
        self.segment_start_support.is_sv() || self.segment_end_support.is_sv()
    }
    pub fn chr(&self) -> Chromosome {
        chr(&self.chromosome)
    }
}

pub const COPY_NUMBER_HEADER: &str = "chromosome\tstart\tend\tcopyNumber\tbafCount\tobservedBAF\tbaf\tsegmentStartSupport\tsegmentEndSupport\tmethod\tdepthWindowCount\tgcContent\tminStart\tmaxStart\tminorAlleleCopyNumber\tmajorAlleleCopyNumber";

fn f4(v: f64) -> String {
    hmf_common::fmt::decimal_format_min(v, 4, 4)
}

/// PurpleCopyNumberFile.write
pub fn write_copy_numbers(path: &str, copy_numbers: &[PurpleCopyNumber]) -> anyhow::Result<()> {
    let mut out = String::with_capacity(copy_numbers.len() * 160 + 200);
    out.push_str(COPY_NUMBER_HEADER);
    out.push('\n');
    for c in copy_numbers {
        out.push_str(
            &[
                c.chromosome.clone(),
                c.start.to_string(),
                c.end.to_string(),
                f4(c.average_tumor_copy_number),
                c.baf_count.to_string(),
                f4(c.average_observed_baf),
                f4(c.average_actual_baf),
                c.segment_start_support.name().into(),
                c.segment_end_support.name().into(),
                c.method.name().into(),
                c.depth_window_count.to_string(),
                f4(c.gc_content),
                c.min_start.to_string(),
                c.max_start.to_string(),
                f4(c.minor_allele_copy_number()),
                f4(c.major_allele_copy_number()),
            ]
            .join("\t"),
        );
        out.push('\n');
    }
    std::fs::write(path, out).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}

/// PurpleCopyNumber.buildChromosomeMap: consecutive copy numbers of a chromosome
pub fn chromosome_map(copy_numbers: &[PurpleCopyNumber]) -> HashMap<String, Vec<PurpleCopyNumber>> {
    let mut map: HashMap<String, Vec<PurpleCopyNumber>> = HashMap::new();
    let mut current = String::new();
    for c in copy_numbers {
        if c.chromosome != current {
            current = c.chromosome.clone();
            map.insert(current.clone(), Vec::new());
        }
        map.get_mut(&current).unwrap().push(c.clone());
    }
    map
}

/// Multimaps.fromRegions for copy numbers: by chromosome, in list order
pub fn by_chromosome(copy_numbers: &[PurpleCopyNumber]) -> BTreeMap<Chromosome, Vec<PurpleCopyNumber>> {
    let mut map: BTreeMap<Chromosome, Vec<PurpleCopyNumber>> = BTreeMap::new();
    for c in copy_numbers {
        if Chromosome::contains(&c.chromosome) {
            map.entry(c.chr()).or_default().push(c.clone());
        }
    }
    map
}

// ---------------------------------------------------------------------------------------------
// CombinedRegion

#[derive(Clone, Debug)]
pub struct CombinedRegion {
    pub combined: ObservedRegion,
    pub method: CopyNumberMethod,
    pub inferred_baf: bool,
    pub regions: Vec<ObservedRegion>,
    unweighted_count: i64,
    is_baf_weighted: bool,
}

fn weighted_average(current_weight: i64, current_value: f64, new_weight: i64, new_value: f64) -> f64 {
    if is_zero(current_value) {
        return new_value;
    }
    let total = current_weight + new_weight;
    (current_weight as f64 * current_value + new_weight as f64 * new_value) / total as f64
}

impl CombinedRegion {
    pub fn new(baf_weighted: bool, region: &ObservedRegion) -> CombinedRegion {
        let mut combined = region.clone();
        if region.germline_status != GermlineStatus::Diploid {
            combined.baf_count = 0;
        }
        CombinedRegion {
            combined,
            method: CopyNumberMethod::Unknown,
            inferred_baf: false,
            regions: vec![region.clone()],
            unweighted_count: 1,
            is_baf_weighted: baf_weighted,
        }
    }
    pub fn start(&self) -> i32 {
        self.combined.start
    }
    pub fn end(&self) -> i32 {
        self.combined.end
    }
    pub fn chromosome(&self) -> &str {
        &self.combined.chromosome
    }
    pub fn tumor_copy_number(&self) -> f64 {
        self.combined.tumor_copy_number
    }
    pub fn tumor_baf(&self) -> f64 {
        self.combined.tumor_baf
    }
    pub fn baf_count(&self) -> i32 {
        self.combined.baf_count
    }
    pub fn is_processed(&self) -> bool {
        self.method != CopyNumberMethod::Unknown
    }
    pub fn support(&self) -> SegmentSupport {
        self.combined.support
    }

    pub fn extend(&mut self, region: &ObservedRegion) {
        let c = &mut self.combined;
        c.start = c.start.min(region.start);
        c.end = c.end.max(region.end);
        c.min_start = c.min_start.min(region.min_start);
        c.max_start = c.max_start.min(region.max_start);
        if region.start <= c.start {
            self.regions.insert(0, region.clone());
            c.support = region.support;
            c.ratio_support = region.ratio_support;
        } else {
            self.regions.push(region.clone());
        }
    }

    pub fn extend_with_unweighted_average(&mut self, region: &ObservedRegion) {
        self.extend(region);
        let w = self.unweighted_count;
        self.apply_depth_weights(region, w, 1);
        self.apply_baf_weights(region, w, 1);
        self.unweighted_count += 1;
    }

    pub fn extend_with_weighted_average(&mut self, region: &ObservedRegion) {
        self.combined.germline_status = GermlineStatus::Diploid;
        let mut current = self.combined.depth_window_count as i64;
        let mut new = region.depth_window_count as i64;
        self.apply_depth_weights(region, current, new);
        if self.is_baf_weighted && (self.combined.baf_count > 0 || region.baf_count > 0) {
            current = self.combined.baf_count as i64;
            new = region.baf_count as i64;
        }
        self.apply_baf_weights(region, current, new);
        self.extend(region);
    }

    fn apply_baf_weights(&mut self, region: &ObservedRegion, cw: i64, nw: i64) {
        let c = &mut self.combined;
        if !is_zero(region.observed_baf) {
            c.observed_baf = weighted_average(cw, c.observed_baf, nw, region.observed_baf);
        }
        if !is_zero(region.tumor_baf) {
            c.tumor_baf = weighted_average(cw, c.tumor_baf, nw, region.tumor_baf);
        }
        c.baf_count += region.baf_count;
    }

    fn apply_depth_weights(&mut self, region: &ObservedRegion, cw: i64, nw: i64) {
        let c = &mut self.combined;
        if !is_zero(region.tumor_copy_number) {
            c.tumor_copy_number = weighted_average(cw, c.tumor_copy_number, nw, region.tumor_copy_number);
        }
        if !is_zero(region.ref_normalised_copy_number) {
            c.ref_normalised_copy_number = weighted_average(cw, c.ref_normalised_copy_number, nw, region.ref_normalised_copy_number);
        }
        if !is_zero(region.gc_content) {
            c.gc_content = weighted_average(cw, c.gc_content, nw, region.gc_content);
        }
        c.depth_window_count += region.depth_window_count;
    }

    pub fn set_tumor_copy_number(&mut self, method: CopyNumberMethod, cn: f64) {
        self.method = method;
        self.combined.tumor_copy_number = cn;
    }

    pub fn set_inferred_tumor_baf(&mut self, baf: f64) {
        self.inferred_baf = true;
        self.combined.tumor_baf = baf;
        self.combined.baf_count = 0;
        self.combined.observed_baf = 0.0;
    }
}

// ---------------------------------------------------------------------------------------------
// AlleleTolerance and ExtendDiploid

pub struct AlleleTolerance {
    pub purity: f64,
    pub targeted: bool,
}

/// AlleleTolerance.relativeCopyNumberChange
pub fn relative_copy_number_change(a: f64, b: f64) -> f64 {
    let diff = (jmax(a, b) - jmin(a, b)).abs();
    if is_zero(diff) {
        return 0.0;
    }
    if is_zero(a) || is_zero(b) {
        return 1.0;
    }
    diff / jmin(a, b).abs()
}

impl AlleleTolerance {
    pub fn in_tolerance(&self, first: &ObservedRegion, second: &ObservedRegion) -> bool {
        let purity_adjustment = jmax(1.0, MAX_DEVIATION_ADJUSTMENT / self.purity);
        let tol = |min: f64, additional: f64, samples: i32| min + additional / (samples as f64).sqrt();
        let min_baf_count = first.baf_count.min(second.baf_count);
        if min_baf_count > 0 {
            let baf_dev = (first.observed_baf - second.observed_baf).abs();
            let max_cn = jmax(first.tumor_copy_number, second.tumor_copy_number);
            let max_minor_dev = purity_adjustment * tol(MIN_ABS_CN_TOLERANCE, 0.5 * max_cn, min_baf_count);
            let minor_dev = (first.minor_allele_copy_number() - second.minor_allele_copy_number()).abs();
            if greater_than(minor_dev, max_minor_dev) && greater_than(baf_dev, MIN_OBSERVED_BAF_CHANGE) {
                return false;
            }
        }
        let min_depth = first.depth_window_count.min(second.depth_window_count);
        let abs_add = if self.targeted { MIN_ABS_CN_ADDITION_TARGETED } else { MIN_ABS_CN_ADDITION };
        let rel_add = if self.targeted { MIN_REL_CN_ADDITION_TARGETED } else { MIN_REL_CN_ADDITION };
        let abs_tol = purity_adjustment * tol(MIN_ABS_CN_TOLERANCE, abs_add, min_depth);
        let rel_tol = purity_adjustment * tol(MIN_REL_CN_TOLERANCE, rel_add, min_depth);
        let in_abs = |a: f64, b: f64| less_or_equal((a - b).abs(), abs_tol);
        let in_rel = |a: f64, b: f64| less_or_equal(relative_copy_number_change(a, b), rel_tol);
        let cn = in_abs(first.tumor_copy_number, second.tumor_copy_number) || in_rel(first.tumor_copy_number, second.tumor_copy_number);
        let ref_cn = in_abs(first.ref_normalised_copy_number, second.ref_normalised_copy_number)
            || in_rel(first.ref_normalised_copy_number, second.ref_normalised_copy_number);
        cn || ref_cn
    }
}

/// ExtendUtils.doNotExtend
fn do_not_extend(target: &CombinedRegion, neighbour: &ObservedRegion) -> bool {
    let centromere = if target.start() < neighbour.start {
        neighbour.support == SegmentSupport::Centromere
    } else {
        target.combined.support == SegmentSupport::Centromere
    };
    let sv = if target.start() < neighbour.start { neighbour.support.is_sv() } else { target.combined.support.is_sv() };
    centromere || sv
}

pub struct ExtendDiploid<'a> {
    pub tolerance: &'a AlleleTolerance,
    pub min_tumor_count: i32,
    pub centromere_min_tumor_count: i32,
}

fn mv(index: i64, right: bool) -> i64 {
    if right { index + 1 } else { index - 1 }
}

impl ExtendDiploid<'_> {
    pub fn extend_diploid(&self, fitted: &[ObservedRegion]) -> Vec<CombinedRegion> {
        let baf_weighted = fitted.iter().any(|x| x.baf_count >= 50);
        let mut regions: Vec<CombinedRegion> = fitted.iter().map(|r| CombinedRegion::new(baf_weighted, r)).collect();
        let mut index = next_index(&regions);
        while index > -1 {
            let i = index as usize;
            regions[i].method = CopyNumberMethod::BafWeighted;
            self.extend_right(&mut regions, i);
            self.extend_left(&mut regions, i);
            index = next_index(&regions);
        }
        regions
    }

    fn extend_right(&self, regions: &mut Vec<CombinedRegion>, target: usize) {
        let neighbour = target + 1;
        while neighbour < regions.len() {
            if !self.merge(regions, true, target) {
                return;
            }
            regions.remove(neighbour);
        }
    }

    fn extend_left(&self, regions: &mut Vec<CombinedRegion>, target: usize) {
        let mut neighbour = target as i64 - 1;
        while neighbour >= 0 {
            if !self.merge(regions, false, (neighbour + 1) as usize) {
                return;
            }
            regions.remove(neighbour as usize);
            neighbour -= 1;
        }
    }

    fn is_valid(min_tumor_count: i32, r: &ObservedRegion) -> bool {
        r.germline_status == GermlineStatus::Diploid && (r.support.is_sv() || r.depth_window_count >= min_tumor_count)
    }
    fn is_dubious(min_tumor_count: i32, r: &ObservedRegion) -> bool {
        r.germline_status == GermlineStatus::Diploid && !r.support.is_sv() && r.depth_window_count < min_tumor_count
    }

    fn next_big_break_is_centromere_or_telomere(regions: &[CombinedRegion], right: bool, target: usize) -> bool {
        let mut i = mv(target as i64, right);
        while i >= 0 && (i as usize) < regions.len() {
            let n = &regions[i as usize].combined;
            if n.support == SegmentSupport::Centromere {
                return true;
            }
            if n.support.is_sv() {
                return false;
            }
            i = mv(i, right);
        }
        true
    }

    fn push_through_dubious(&self, min_tumor_count: i32, regions: &[CombinedRegion], right: bool, target_index: usize) -> bool {
        let mut dubious = 0;
        let target = &regions[target_index];
        let mut i = mv(target_index as i64, right);
        while i >= 0 && (i as usize) < regions.len() {
            let n = &regions[i as usize].combined;
            if n.start > target.start() {
                if n.support == SegmentSupport::Centromere {
                    return dubious < min_tumor_count;
                }
                if n.support.is_sv() {
                    return false;
                }
            }
            let is_dubious = Self::is_dubious(min_tumor_count, n);
            let in_tol = self.tolerance.in_tolerance(&target.combined, n);
            if is_dubious && !in_tol {
                dubious += n.depth_window_count;
                if dubious >= min_tumor_count {
                    return false;
                }
            }
            if Self::is_valid(min_tumor_count, n) {
                return in_tol;
            }
            if n.start < target.start() {
                if n.support == SegmentSupport::Centromere {
                    return dubious < min_tumor_count;
                }
                if n.support.is_sv() {
                    return false;
                }
            }
            i = mv(i, right);
        }
        dubious < min_tumor_count
    }

    fn merge(&self, regions: &mut [CombinedRegion], right: bool, target_index: usize) -> bool {
        let neighbour = regions[mv(target_index as i64, right) as usize].combined.clone();
        if do_not_extend(&regions[target_index], &neighbour) {
            return false;
        }
        let min_tumor_count = if Self::next_big_break_is_centromere_or_telomere(regions, right, target_index) {
            self.centromere_min_tumor_count
        } else {
            self.min_tumor_count
        };
        if Self::is_dubious(min_tumor_count, &neighbour) {
            if self.tolerance.in_tolerance(&regions[target_index].combined, &neighbour) {
                regions[target_index].extend_with_weighted_average(&neighbour);
                return true;
            } else if self.push_through_dubious(min_tumor_count, regions, right, target_index) {
                regions[target_index].extend(&neighbour);
                return true;
            }
            return false;
        }
        if !Self::is_valid(min_tumor_count, &neighbour) {
            regions[target_index].extend(&neighbour);
            return true;
        } else if self.tolerance.in_tolerance(&regions[target_index].combined, &neighbour) {
            regions[target_index].extend_with_weighted_average(&neighbour);
            return true;
        }
        false
    }
}

fn next_index(regions: &[CombinedRegion]) -> i64 {
    let (mut idx_baf, mut idx_depth) = (-1i64, -1i64);
    let mut largest_baf = 0;
    let mut largest_depth: i64 = 0;
    for (i, c) in regions.iter().enumerate() {
        let r = &c.combined;
        if !c.is_processed() && r.germline_status == GermlineStatus::Diploid {
            if r.baf_count > largest_baf {
                largest_baf = r.baf_count;
                idx_baf = i as i64;
            }
            if r.depth_window_count as i64 > largest_depth {
                largest_depth = r.depth_window_count as i64;
                idx_depth = i as i64;
            }
        }
    }
    if idx_baf > -1 { idx_baf } else if idx_depth > -1 { idx_depth } else { -1 }
}

// ---------------------------------------------------------------------------------------------
// ExtendRegion (non-diploid, structural variant)

#[derive(Clone, Copy, PartialEq)]
enum ExtendKind {
    NonDiploid,
    StructuralVariant,
}

fn extend_region(kind: ExtendKind, regions: &mut Vec<CombinedRegion>) {
    let method = match kind {
        ExtendKind::NonDiploid => CopyNumberMethod::NonDiploid,
        ExtendKind::StructuralVariant => CopyNumberMethod::StructuralVariant,
    };
    let extend_with = |target: &mut CombinedRegion, n: &ObservedRegion| match kind {
        ExtendKind::NonDiploid => target.extend_with_weighted_average(n),
        ExtendKind::StructuralVariant => target.extend_with_unweighted_average(n),
    };
    let mut i: i64 = 0;
    while (i as usize) < regions.len() {
        let iu = i as usize;
        if regions[iu].method == method {
            // extendRight
            let ti = iu + 1;
            while ti < regions.len() {
                let n = regions[ti].combined.clone();
                if do_not_extend(&regions[iu], &n) {
                    break;
                }
                if !regions[ti].is_processed() {
                    regions[iu].extend(&n);
                } else if regions[ti].method == method {
                    extend_with(&mut regions[iu], &n);
                } else {
                    break;
                }
                regions.remove(ti);
            }
            // extendLeft
            let mut start = iu;
            let mut ti = iu as i64 - 1;
            while ti >= 0 {
                let n = regions[ti as usize].combined.clone();
                if do_not_extend(&regions[start], &n) {
                    break;
                }
                if !regions[ti as usize].is_processed() {
                    regions[start].extend(&n);
                } else if regions[ti as usize].method == method {
                    extend_with(&mut regions[start], &n);
                } else {
                    break;
                }
                regions.remove(ti as usize);
                start -= 1;
                ti -= 1;
            }
            i -= (iu - start) as i64;
        }
        i += 1;
    }
}

/// ExtendNonDiploid.isEligible
fn non_diploid_eligible(region: &ObservedRegion, neighbour: Option<&ObservedRegion>) -> bool {
    let bound_by_centromere = region.support == SegmentSupport::Centromere || neighbour.is_some_and(|n| n.support == SegmentSupport::Centromere);
    let bound_by_sv = region.support.is_sv() || neighbour.is_some_and(|n| n.support.is_sv());
    greater_than(region.observed_tumor_ratio, region.observed_normal_ratio)
        && !matches!(region.germline_status, GermlineStatus::Excluded | GermlineStatus::Diploid | GermlineStatus::Noise)
        && region.depth_window_count > 0
        && !bound_by_centromere
        && bound_by_sv
}

/// ExtendNonDiploid.nonDiploid
pub fn extend_non_diploid(mut regions: Vec<CombinedRegion>) -> Vec<CombinedRegion> {
    for i in 0..regions.len() {
        let next = regions.get(i + 1).map(|r| r.combined.clone());
        if regions[i].method == CopyNumberMethod::Unknown && non_diploid_eligible(&regions[i].combined, next.as_ref()) {
            let cn = regions[i].combined.ref_normalised_copy_number;
            regions[i].set_tumor_copy_number(CopyNumberMethod::NonDiploid, cn);
        }
    }
    extend_region(ExtendKind::NonDiploid, &mut regions);
    regions
}

// ---------------------------------------------------------------------------------------------
// ExtendLongArm and populateUnknown

fn short_arm(c: Chromosome) -> bool {
    matches!(c.short_name(), "13" | "14" | "15" | "21" | "22")
}

/// ExtendLongArm.extendLongArm
pub fn extend_long_arm(regions: &mut Vec<CombinedRegion>) {
    let mut centromere_index: i64 = -1;
    for (i, r) in regions.iter().enumerate() {
        let contig = r.chromosome();
        if !Chromosome::contains(contig) || !(short_arm(chr(contig)) || is_y(chr(contig))) {
            break;
        }
        if r.combined.support == SegmentSupport::Centromere {
            centromere_index = i as i64;
            break;
        }
    }
    if centromere_index > 0 {
        let ci = centromere_index as usize;
        if regions[ci].method != CopyNumberMethod::Unknown {
            let cn = regions[ci].tumor_copy_number();
            long_arm_extend_left(cn, ci as i64 - 1, regions);
        }
    }
}

fn long_arm_is_processed(r: &CombinedRegion) -> bool {
    if short_arm(chr(r.chromosome())) {
        return false;
    }
    r.is_processed()
}

fn long_arm_extend_left(cn: f64, target_index: i64, regions: &mut Vec<CombinedRegion>) {
    if target_index < 0 {
        return;
    }
    let ti = target_index as usize;
    if long_arm_is_processed(&regions[ti]) {
        return;
    }
    regions[ti].set_tumor_copy_number(CopyNumberMethod::LongArm, cn);
    if regions[ti].combined.support != SegmentSupport::None {
        long_arm_extend_left(cn, target_index - 1, regions);
        return;
    }
    let mut target = ti;
    let mut idx = target_index - 1;
    while idx >= 0 {
        if long_arm_is_processed(&regions[idx as usize]) {
            return;
        }
        let n = regions[idx as usize].combined.clone();
        regions[target].extend(&n);
        regions.remove(idx as usize);
        target -= 1;
        idx -= 1;
        if regions[target].combined.support != SegmentSupport::None {
            long_arm_extend_left(cn, idx, regions);
            return;
        }
    }
}

/// ExtendUtils.populateUnknown
pub fn populate_unknown(regions: &mut Vec<CombinedRegion>, cobalt: &CobaltChromosomes) {
    let mut i: i64 = 0;
    while (i as usize) < regions.len() {
        let iu = i as usize;
        if regions[iu].method == CopyNumberMethod::Unknown {
            let ratio = cobalt.get(regions[iu].chromosome()).map(|c| c.actual_ratio).unwrap_or_else(|| panic!("NullPointerException: no COBALT chromosome {}", regions[iu].chromosome()));
            regions[iu].set_tumor_copy_number(CopyNumberMethod::Unknown, 2.0 * ratio);
            if regions[iu].support() == SegmentSupport::None && iu > 0 && regions[iu - 1].method == CopyNumberMethod::Unknown {
                let r = regions[iu].combined.clone();
                regions[iu - 1].extend(&r);
                regions.remove(iu);
                i -= 1;
            }
        }
        i += 1;
    }
}

// ---------------------------------------------------------------------------------------------
// ExtendDiploidBAF

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InferRegion {
    pub left_source: i64,
    pub left_target: i64,
    pub right_target: i64,
    pub right_source: i64,
}

impl InferRegion {
    fn valid(&self) -> bool {
        self.left_valid() || self.right_valid()
    }
    fn left_valid(&self) -> bool {
        self.left_source != -1
    }
    fn right_valid(&self) -> bool {
        self.right_source != -1
    }
}

const INVALID_PAIR: InferRegion = InferRegion { left_source: -1, left_target: -1, right_target: -1, right_source: -1 };

pub struct ExtendDiploidBaf {
    /// start cnaPosition -> end cnaPosition of simple DUPs (HashMap.put: a later DUP at the same position wins)
    pub simple_dups: HashMap<i32, i32>,
}

fn ignore_support(s: SegmentSupport) -> bool {
    matches!(s, SegmentSupport::Centromere | SegmentSupport::Telomere | SegmentSupport::Unknown)
}
fn major_different(a: &ObservedRegion, b: &ObservedRegion) -> bool {
    greater_than((a.major_allele_copy_number() - b.major_allele_copy_number()).abs(), 0.5)
}
fn minor_different(a: &ObservedRegion, b: &ObservedRegion) -> bool {
    greater_than((a.minor_allele_copy_number() - b.minor_allele_copy_number()).abs(), 0.5)
}
fn either_different(a: &ObservedRegion, b: &ObservedRegion) -> bool {
    major_different(a, b) || minor_different(a, b)
}

/// ExtendDiploidBAF.bafForTargetAllele
pub fn baf_for_target_allele(target_allele: f64, copy_number: f64) -> f64 {
    if less_or_equal(copy_number, 1.0) {
        return 1.0;
    }
    let major = jmin(target_allele, copy_number);
    let minor = copy_number - major;
    jmax(major, minor) / (major + minor)
}

/// ExtendDiploidBAF.minorOrMajorMovedTargetPloidy
pub fn minor_or_major_moved_target_ploidy(source: &ObservedRegion, primary: &ObservedRegion, secondary: &ObservedRegion) -> f64 {
    let minor_diff = minor_different(primary, secondary);
    let major_diff = major_different(primary, secondary);
    let closest = |target: f64| {
        let dmin = (source.minor_allele_copy_number() - target).abs();
        let dmaj = (source.major_allele_copy_number() - target).abs();
        if less_than(dmin, dmaj) { source.minor_allele_copy_number() } else { source.major_allele_copy_number() }
    };
    if minor_diff && major_diff {
        if less_than((primary.minor_allele_copy_number() - secondary.major_allele_copy_number()).abs(), 0.5) {
            return closest((primary.minor_allele_copy_number() + secondary.major_allele_copy_number()) / 2.0);
        }
        if less_than((primary.major_allele_copy_number() - secondary.minor_allele_copy_number()).abs(), 0.5) {
            return closest((primary.major_allele_copy_number() + secondary.minor_allele_copy_number()) / 2.0);
        }
    }
    if major_diff {
        return source.minor_allele_copy_number();
    }
    source.major_allele_copy_number()
}

/// ExtendDiploidBAF.nextRegion
pub fn next_region(cross_centromere: bool, regions: &[CombinedRegion]) -> InferRegion {
    let mut left_source: i64 = -1;
    let mut should_infer = false;
    let mut left_target: i64 = i32::MAX as i64;
    for (i, s) in regions.iter().enumerate() {
        let i = i as i64;
        if s.support() == SegmentSupport::Centromere && !cross_centromere {
            if should_infer && left_source > -1 {
                return InferRegion { left_source, left_target, right_target: i - 1, right_source: -1 };
            } else {
                left_source = -1;
                should_infer = false;
                left_target = i32::MAX as i64;
            }
        }
        if s.baf_count() == 0 && !s.inferred_baf {
            should_infer = true;
            left_target = left_target.min(i);
        } else if should_infer {
            return InferRegion { left_source, left_target, right_target: i - 1, right_source: i };
        } else {
            left_source = i;
        }
    }
    if should_infer { InferRegion { left_source, left_target, right_target: regions.len() as i64 - 1, right_source: -1 } } else { INVALID_PAIR }
}

impl ExtendDiploidBaf {
    pub fn new(simple_variants: &[&StructuralVariant]) -> ExtendDiploidBaf {
        let mut simple_dups = HashMap::new();
        for v in simple_variants {
            if v.sv_type == SvType::Dup {
                if let Some(end) = &v.end {
                    simple_dups.insert(v.start.cna_position(), end.cna_position());
                }
            }
        }
        ExtendDiploidBaf { simple_dups }
    }

    pub fn extend_baf(&self, regions: &mut [CombinedRegion]) {
        for cross in [false, true] {
            let mut ir = next_region(cross, regions);
            while ir.valid() {
                self.infer_between(ir, regions);
                ir = next_region(cross, regions);
            }
        }
    }

    fn r(regions: &[CombinedRegion], i: i64) -> &ObservedRegion {
        &regions[i as usize].combined
    }

    fn look_right(ir: InferRegion, regions: &[CombinedRegion]) -> Option<ObservedRegion> {
        let src = Self::r(regions, ir.right_source);
        for i in (ir.right_source + 1) as usize..regions.len() {
            let right = &regions[i].combined;
            if right.support == SegmentSupport::Centromere {
                return None;
            }
            if right.baf_count > 0 && !ignore_support(right.support) && either_different(right, src) {
                return Some(right.clone());
            }
        }
        None
    }

    fn look_left(ir: InferRegion, regions: &[CombinedRegion]) -> Option<ObservedRegion> {
        let src = Self::r(regions, ir.left_source);
        let mut i = ir.left_source - 1;
        while i >= 0 {
            if regions[(i + 1) as usize].combined.support == SegmentSupport::Centromere {
                return None;
            }
            let left = &regions[i as usize].combined;
            if left.baf_count > 0 && !ignore_support(left.support) && either_different(left, src) {
                return Some(left.clone());
            }
            i -= 1;
        }
        None
    }

    fn is_small_region(max: i32, ir: InferRegion, regions: &[CombinedRegion]) -> bool {
        regions[ir.right_target as usize].end() - regions[ir.left_target as usize].start() + 1 <= max
    }

    fn single_tiny_with_greater_cn(ir: InferRegion, regions: &[CombinedRegion], source: &ObservedRegion) -> bool {
        ir.left_target == ir.right_target
            && Self::is_small_region(30, ir, regions)
            && greater_than(regions[ir.left_target as usize].tumor_copy_number(), source.tumor_copy_number)
    }

    fn single_source_target_ploidy(ir: InferRegion, regions: &[CombinedRegion]) -> f64 {
        let (source, different) = if ir.left_valid() {
            (Self::r(regions, ir.left_source).clone(), Self::look_left(ir, regions))
        } else {
            (Self::r(regions, ir.right_source).clone(), Self::look_right(ir, regions))
        };
        if Self::single_tiny_with_greater_cn(ir, regions, &source) {
            return source.minor_allele_copy_number();
        }
        if let Some(d) = different {
            if minor_different(&source, &d) {
                return source.major_allele_copy_number();
            }
        }
        source.minor_allele_copy_number()
    }

    fn is_simple_dup_surrounded_by_loh(&self, ir: InferRegion, regions: &[CombinedRegion]) -> bool {
        if !(ir.left_target == ir.right_target && ir.left_source != -1 && ir.right_source != -1) {
            return false;
        }
        let target = Self::r(regions, ir.left_target);
        let right = Self::r(regions, ir.right_source);
        if !(target.support == SegmentSupport::Dup && right.support == SegmentSupport::Dup) {
            return false;
        }
        if self.simple_dups.get(&target.start) != Some(&right.start) {
            return false;
        }
        let left = Self::r(regions, ir.left_source);
        less_than(left.minor_allele_copy_number(), 0.5) && less_than(right.minor_allele_copy_number(), 0.5)
    }

    fn multi_source_target_ploidy(&self, ir: InferRegion, regions: &[CombinedRegion]) -> f64 {
        let (primary, secondary) = if regions[ir.left_source as usize].baf_count() > regions[ir.right_source as usize].baf_count() {
            (Self::r(regions, ir.left_source).clone(), Self::r(regions, ir.right_source).clone())
        } else {
            (Self::r(regions, ir.right_source).clone(), Self::r(regions, ir.left_source).clone())
        };
        // minCopyNumberLessThanSourceMajorAllelePloidies
        let mut min_cn = f64::MAX;
        for i in ir.left_target..=ir.right_target {
            min_cn = jmin(min_cn, regions[i as usize].tumor_copy_number());
        }
        if greater_than(Self::r(regions, ir.left_source).major_allele_copy_number() - min_cn, 0.5)
            && greater_than(Self::r(regions, ir.right_source).major_allele_copy_number() - min_cn, 0.5)
        {
            return primary.minor_allele_copy_number();
        }
        if minor_different(&primary, &secondary) || major_different(&primary, &secondary) {
            return minor_or_major_moved_target_ploidy(&primary, &primary, &secondary);
        }
        if Self::single_tiny_with_greater_cn(ir, regions, &primary) {
            return primary.minor_allele_copy_number();
        }
        if self.is_simple_dup_surrounded_by_loh(ir, regions) {
            return primary.minor_allele_copy_number();
        }
        // nearestDifferentRegion
        let left_source = Self::r(regions, ir.left_source).clone();
        let right_source = Self::r(regions, ir.right_source).clone();
        let opt_right = Self::look_right(ir, regions);
        let opt_left = Self::look_left(ir, regions);
        let nearest = match (opt_left, opt_right) {
            (l, None) => l,
            (None, r) => r,
            (Some(l), Some(r)) => {
                let ld = left_source.end - l.start;
                let rd = r.start - right_source.start;
                if ld < rd { Some(l) } else { Some(r) }
            }
        };
        if let Some(nearest) = nearest {
            // isSingleSmallRegionFlankedByLargeLOH
            let lt = &regions[ir.left_target as usize];
            let distance = if nearest.start > lt.end() { nearest.start - lt.end() } else { lt.start() - nearest.end };
            let is_loh = |r: &CombinedRegion| less_or_equal(r.combined.minor_allele_copy_number(), 0.5);
            let flanked = is_loh(&regions[ir.left_source as usize]) && is_loh(&regions[ir.right_source as usize]);
            if Self::is_small_region(1000, ir, regions) && flanked && distance > 1_000_000 {
                return primary.minor_allele_copy_number();
            }
            let other = if nearest.start > right_source.start { &right_source } else { &left_source };
            return minor_or_major_moved_target_ploidy(&primary, &nearest, other);
        }
        primary.minor_allele_copy_number()
    }

    fn infer_between(&self, ir: InferRegion, regions: &mut [CombinedRegion]) {
        let target = if ir.left_valid() ^ ir.right_valid() {
            Self::single_source_target_ploidy(ir, regions)
        } else {
            self.multi_source_target_ploidy(ir, regions)
        };
        for i in ir.left_target..=ir.right_target {
            let r = &mut regions[i as usize];
            let cn = r.tumor_copy_number();
            r.set_inferred_tumor_baf(baf_for_target_allele(target, cn));
        }
    }
}

// ---------------------------------------------------------------------------------------------
// SV leg copy numbers and ploidies

#[derive(Clone, Debug)]
pub struct LegPloidy {
    pub leg: SvLeg,
    pub left_copy_number: Option<f64>,
    pub right_copy_number: Option<f64>,
    pub observed_vaf: f64,
    pub adjusted_vaf: f64,
    pub unweighted_implied_ploidy: f64,
    pub average_implied_ploidy: f64,
    pub weight: f64,
}

impl LegPloidy {
    fn decreasing_from_zero(&self, orientation: i8, cn: Option<f64>) -> bool {
        self.leg.orientation == orientation && cn.is_some_and(is_zero)
    }
    pub fn implied_right_copy_number_weight(&self) -> f64 {
        if self.decreasing_from_zero(1, self.left_copy_number) || self.left_copy_number.is_some() { self.weight } else { 0.0 }
    }
    pub fn implied_left_copy_number_weight(&self) -> f64 {
        if self.decreasing_from_zero(-1, self.right_copy_number) || self.right_copy_number.is_some() { self.weight } else { 0.0 }
    }
    pub fn implied_right_copy_number(&self) -> f64 {
        if self.decreasing_from_zero(1, self.left_copy_number) {
            return 0.0;
        }
        match self.left_copy_number {
            Some(l) => l - self.leg.orientation as f64 * self.average_implied_ploidy,
            None => 0.0,
        }
    }
    pub fn implied_left_copy_number(&self) -> f64 {
        if self.decreasing_from_zero(-1, self.right_copy_number) {
            return 0.0;
        }
        match self.right_copy_number {
            Some(r) => r + self.leg.orientation as f64 * self.average_implied_ploidy,
            None => 0.0,
        }
    }
    /// StructuralVariantLegCopyNumber.adjustedCopyNumber
    pub fn adjusted_copy_number(&self) -> f64 {
        if self.leg.orientation == 1 { self.left_copy_number.unwrap_or(0.0) } else { self.right_copy_number.unwrap_or(0.0) }
    }
}

/// copy numbers by chromosome as a region selector: the region containing a position
pub trait RegionCopyNumbers {
    /// the copy number (after the extractor; None where the extractor gives null) of the region containing the position
    fn copy_number_at(&self, chromosome: &str, position: i32) -> Option<Option<f64>>;
}

pub struct CopyNumberSelector<'a> {
    /// per chromosome: (start, end, copy number) sorted by start
    pub regions: BTreeMap<Chromosome, Vec<(i32, i32, f64)>>,
    pub _p: std::marker::PhantomData<&'a ()>,
}

impl<'a> CopyNumberSelector<'a> {
    pub fn from_copy_numbers(copy_numbers: &BTreeMap<Chromosome, Vec<PurpleCopyNumber>>) -> CopyNumberSelector<'a> {
        let mut regions = BTreeMap::new();
        for (c, list) in copy_numbers {
            let mut v: Vec<(i32, i32, f64)> = list.iter().map(|x| (x.start, x.end, x.average_tumor_copy_number)).collect();
            v.sort_by_key(|x| x.0);
            regions.insert(*c, v);
        }
        CopyNumberSelector { regions, _p: std::marker::PhantomData }
    }
    pub fn from_combined(combined: &BTreeMap<Chromosome, Vec<CombinedRegion>>) -> CopyNumberSelector<'a> {
        let mut regions = BTreeMap::new();
        for (c, list) in combined {
            let mut v: Vec<(i32, i32, f64)> =
                list.iter().filter(|x| x.is_processed()).map(|x| (x.start(), x.end(), x.tumor_copy_number())).collect();
            v.sort_by_key(|x| x.0);
            regions.insert(*c, v);
        }
        CopyNumberSelector { regions, _p: std::marker::PhantomData }
    }
    pub fn at(&self, chromosome: &str, position: i32) -> Option<f64> {
        let list = self.regions.get(&chr(chromosome))?;
        let i = list.partition_point(|r| r.0 <= position);
        if i == 0 {
            return None;
        }
        let r = list[i - 1];
        if position >= r.0 && position <= r.1 { Some(r.2) } else { None }
    }
}

pub struct LegPloidyFactory<'a> {
    pub average_copy_number: f64,
    pub average_read_depth: i32,
    pub purity_adjuster: &'a PurityAdjuster,
}

/// StructuralVariantLegs
#[derive(Clone, Debug)]
pub struct Legs {
    pub start: Option<SvLeg>,
    pub end: Option<SvLeg>,
}

fn reciprocal_vaf(vaf: f64) -> f64 {
    vaf / (1.0 - vaf)
}

impl LegPloidyFactory<'_> {
    /// SvLegPloidyFactory.create(leg, left, right)
    fn create_leg(&self, leg: &SvLeg, left: Option<f64>, right: Option<f64>) -> Option<LegPloidy> {
        let (larger, smaller) = if leg.orientation == 1 { (left, right) } else { (right, left) };
        if larger.is_none() && smaller.is_none() {
            return None;
        }
        let observed_vaf = leg.allele_frequency?;
        let (adjusted_vaf, ploidy, weight);
        if let Some(cn) = larger {
            adjusted_vaf = self.purity_adjuster.purity_adjusted_vaf(&leg.chromosome, jmax(0.001, cn), observed_vaf);
            ploidy = adjusted_vaf * cn;
            weight = 1.0;
        } else {
            let cn = smaller.unwrap();
            let rvaf = reciprocal_vaf(observed_vaf);
            adjusted_vaf = self.purity_adjuster.purity_adjusted_vaf(&leg.chromosome, jmax(0.001, cn), rvaf);
            if self.average_read_depth > 0 && (adjusted_vaf.is_infinite() || greater_than(adjusted_vaf, reciprocal_vaf(0.75))) {
                match leg.tumor_variant_fragment_count {
                    Some(f) if f > 0 => ploidy = self.average_copy_number * f as f64 / self.average_read_depth as f64,
                    _ => return None,
                }
            } else if !rvaf.is_finite() {
                return None;
            } else {
                ploidy = adjusted_vaf * cn;
            }
            weight = 1.0 / (1.0 + jpow(jmax(cn, 2.0) / jmin(jmax(cn, 0.01), 2.0), 2.0));
        }
        Some(LegPloidy {
            leg: leg.clone(),
            left_copy_number: left,
            right_copy_number: right,
            observed_vaf,
            adjusted_vaf,
            unweighted_implied_ploidy: ploidy,
            average_implied_ploidy: 0.0,
            weight,
        })
    }

    fn leg_copy_numbers(leg: &SvLeg, selector: &CopyNumberSelector) -> (Option<f64>, Option<f64>) {
        let left = selector.at(&leg.chromosome, leg.cna_position() - 1).map(|x| jmax(0.0, x));
        let right = selector.at(&leg.chromosome, leg.cna_position()).map(|x| jmax(0.0, x));
        (left, right)
    }

    /// SvLegPloidyFactory.create(StructuralVariantLegs, copyNumbers)
    pub fn create_legs(&self, legs: &Legs, selector: &CopyNumberSelector) -> Vec<LegPloidy> {
        let start = legs.start.as_ref().and_then(|l| {
            let (a, b) = Self::leg_copy_numbers(l, selector);
            self.create_leg(l, a, b)
        });
        let end = legs.end.as_ref().and_then(|l| {
            let (a, b) = Self::leg_copy_numbers(l, selector);
            self.create_leg(l, a, b)
        });
        if start.is_none() && end.is_none() {
            return Vec::new();
        }
        let sw = start.as_ref().map(|x| x.weight).unwrap_or(0.0);
        let sp = start.as_ref().map(|x| x.unweighted_implied_ploidy).unwrap_or(0.0);
        let ew = end.as_ref().map(|x| x.weight).unwrap_or(0.0);
        let ep = end.as_ref().map(|x| x.unweighted_implied_ploidy).unwrap_or(0.0);
        let total = sw + ew;
        let average = (sw * sp + ew * ep) / total;
        let mut out = Vec::new();
        for mut p in [start, end].into_iter().flatten() {
            p.weight = total;
            p.average_implied_ploidy = average;
            out.push(p);
        }
        sort_legs(&mut out);
        out
    }

    /// SvLegPloidyFactory.create(variant, copyNumbers): one variant, inserts allowed
    pub fn create_variant(&self, variant: &StructuralVariant, selector: &CopyNumberSelector) -> Vec<LegPloidy> {
        let mut out = Vec::new();
        for legs in create_legs_list(true, &[variant]) {
            out.extend(self.create_legs(&legs, selector));
        }
        sort_legs(&mut out);
        out
    }

    /// SvLegPloidyFactory.create(variants, copyNumbers)
    pub fn create_all(&self, variants: &[StructuralVariant], selector: &CopyNumberSelector) -> Vec<LegPloidy> {
        let refs: Vec<&StructuralVariant> = variants.iter().collect();
        let mut out = Vec::new();
        for legs in sv_legs(&refs) {
            out.extend(self.create_legs(&legs, selector));
        }
        sort_legs(&mut out);
        out
    }
}

/// Collections.sort of leg ploidies: GenomePosition order (chromosome rank, position), stable
fn sort_legs(v: &mut [LegPloidy]) {
    v.sort_by(|a, b| position_compare(&a.leg.chromosome, a.leg.position, &b.leg.chromosome, b.leg.position));
}

fn create_legs_list(allow_inserts: bool, variants: &[&StructuralVariant]) -> Vec<Legs> {
    let mut out = Vec::new();
    for v in variants {
        if allow_inserts || v.sv_type != SvType::Ins {
            let start = if v.start.allele_frequency.is_some() { Some(v.start.clone()) } else { None };
            let end = v.end.as_ref().filter(|e| e.allele_frequency.is_some()).cloned();
            if start.is_some() || end.is_some() {
                out.push(Legs { start, end });
            }
        }
    }
    out
}

type CnaPos = (String, i32);
fn cna_pos(l: &SvLeg) -> CnaPos {
    (l.chromosome.clone(), l.cna_position())
}

/// SvLegsFactory.reduce
pub fn reduce_legs(position: &CnaPos, legs: &[SvLeg]) -> SvLeg {
    let (mut max_pos, mut max_neg) = (0.0f64, 0.0f64);
    for l in legs {
        let af = l.allele_frequency.unwrap_or_else(|| panic!("NullPointerException: leg without allele frequency"));
        if l.orientation == 1 {
            max_pos = jmax(max_pos, af);
        } else {
            max_neg = jmax(max_neg, af);
        }
    }
    if is_zero(max_neg) {
        if let Some(l) = legs.iter().find(|l| equal(max_pos, l.allele_frequency.unwrap())) {
            return l.clone();
        }
    }
    if is_zero(max_pos) {
        if let Some(l) = legs.iter().find(|l| equal(max_neg, l.allele_frequency.unwrap())) {
            return l.clone();
        }
    }
    let mut max_ref = 0;
    let mut cumulative = 0;
    for l in legs {
        if let Some(t) = l.tumor_variant_fragment_count {
            cumulative += l.orientation as i32 * t;
        }
        if let Some(r) = l.tumor_reference_fragment_count {
            max_ref = max_ref.max(r);
        }
    }
    let orientation: i8 = if greater_than(max_pos, max_neg) { 1 } else { -1 };
    let vaf = (max_pos - max_neg).abs();
    SvLeg {
        chromosome: position.0.clone(),
        position: if orientation == -1 { position.1 } else { position.1 - 1 },
        orientation,
        start_offset: 0,
        end_offset: 0,
        homology: String::new(),
        allele_frequency: Some(vaf.abs()),
        inexact_homology_offset_start: None,
        inexact_homology_offset_end: None,
        tumor_variant_fragment_count: if cumulative == 0 { None } else { Some(cumulative.abs()) },
        tumor_reference_fragment_count: if max_ref == 0 { None } else { Some(max_ref.abs()) },
        normal_variant_fragment_count: None,
        normal_reference_fragment_count: None,
        anchoring_support_distance: 0,
    }
}

/// SvLegsFactory.create(variants): legs of every non-insert variant, legs at a shared copy-number position
/// reduced to one
pub fn sv_legs(variants: &[&StructuralVariant]) -> Vec<Legs> {
    let mut result = create_legs_list(false, variants);
    // findDuplicatePositions (a HashMap key set; each position is handled independently)
    let mut counts: Vec<(CnaPos, usize)> = Vec::new();
    for l in &result {
        for leg in [&l.start, &l.end].into_iter().flatten() {
            let p = cna_pos(leg);
            match counts.iter_mut().find(|(k, _)| *k == p) {
                Some(e) => e.1 += 1,
                None => counts.push((p, 1)),
            }
        }
    }
    let duplicates: Vec<CnaPos> = counts.into_iter().filter(|(_, n)| *n > 1).map(|(k, _)| k).collect();
    for dp in &duplicates {
        let matching: Vec<usize> = (0..result.len())
            .filter(|&i| {
                let l = &result[i];
                l.start.as_ref().is_some_and(|s| cna_pos(s) == *dp) || l.end.as_ref().is_some_and(|e| cna_pos(e) == *dp)
            })
            .collect();
        let mut dup_legs: Vec<SvLeg> = Vec::new();
        for &i in &matching {
            if let Some(s) = &result[i].start {
                if cna_pos(s) == *dp {
                    dup_legs.push(s.clone());
                }
            }
            if let Some(e) = &result[i].end {
                if cna_pos(e) == *dp {
                    dup_legs.push(e.clone());
                }
            }
        }
        let approved = reduce_legs(dp, &dup_legs);
        let mut matched = false;
        for &i in &matching {
            if let Some(s) = result[i].start.clone() {
                if cna_pos(&s) == *dp {
                    if !matched && s == approved {
                        matched = true;
                    } else {
                        result[i].start = None;
                    }
                }
            }
            if let Some(e) = result[i].end.clone() {
                if cna_pos(&e) == *dp {
                    if !matched && e == approved {
                        matched = true;
                    } else {
                        result[i].end = None;
                    }
                }
            }
        }
        if !matched {
            result.push(Legs { start: Some(approved), end: None });
        }
    }
    result.into_iter().filter(|l| l.start.is_some() || l.end.is_some()).collect()
}

/// SvImpliedCopyNumber.inferCopyNumberFromStructuralVariants
pub fn infer_copy_number_from_svs(start: Option<&LegPloidy>, end: Option<&LegPloidy>) -> f64 {
    let sw = start.map(|x| x.implied_right_copy_number_weight()).unwrap_or(0.0);
    let scn = start.map(|x| x.implied_right_copy_number()).unwrap_or(0.0);
    let ew = end.map(|x| x.implied_left_copy_number_weight()).unwrap_or(0.0);
    let ecn = end.map(|x| x.implied_left_copy_number()).unwrap_or(0.0);
    let weighted = (scn * sw + ecn * ew) / (sw + ew);
    jmin(jmax(0.0, weighted), SV_MAX_INFERRED_COPY_NUMBER as f64)
}

/// SvImpliedCopyNumber.svImpliedCopyNumber
pub fn sv_implied_copy_number(
    factory: &LegPloidyFactory,
    variants: &[StructuralVariant],
    copy_numbers: &mut BTreeMap<Chromosome, Vec<CombinedRegion>>,
) {
    let missing = |cn: &BTreeMap<Chromosome, Vec<CombinedRegion>>| cn.values().flatten().filter(|x| !x.is_processed()).count() as i64;
    let mut previous: i64 = copy_numbers.values().map(|v| v.len() as i64).sum();
    let mut current = missing(copy_numbers);
    while current < previous && current > 0 {
        let selector = CopyNumberSelector::from_combined(copy_numbers);
        let ploidies = factory.create_all(variants, &selector);
        let mut map: HashMap<CnaPos, LegPloidy> = HashMap::new();
        for p in ploidies {
            let k = (p.leg.chromosome.clone(), p.leg.cna_position());
            if map.contains_key(&k) {
                panic!("IllegalStateException: Duplicate key {}:{}", k.0, k.1);
            }
            map.insert(k, p);
        }
        for c in Chromosome::ALL {
            let Some(list) = copy_numbers.get_mut(&c) else { continue };
            let mut inferred = false;
            for cn in list.iter_mut() {
                if !cn.is_processed() {
                    let start = map.get(&(cn.chromosome().to_string(), cn.start())).filter(|x| x.implied_right_copy_number_weight() > 0.0);
                    let end = map.get(&(cn.chromosome().to_string(), cn.end() + 1)).filter(|x| x.implied_left_copy_number_weight() > 0.0);
                    if start.is_some() || end.is_some() {
                        inferred = true;
                        cn.combined.depth_window_count = 0;
                        cn.set_tumor_copy_number(CopyNumberMethod::StructuralVariant, infer_copy_number_from_svs(start, end));
                    }
                }
            }
            if inferred {
                extend_region(ExtendKind::StructuralVariant, list);
            }
        }
        previous = current;
        current = missing(copy_numbers);
    }
}

// ---------------------------------------------------------------------------------------------
// PurpleCopyNumberFactory

pub struct CopyNumberFactory<'a> {
    pub min_tumor_ratio_count: i32,
    pub min_tumor_ratio_count_at_centromere: i32,
    pub targeted: bool,
    pub average_read_depth: i32,
    pub ploidy: f64,
    pub purity_adjuster: &'a PurityAdjuster,
    pub cobalt: &'a CobaltChromosomes,
}

impl CopyNumberFactory<'_> {
    pub fn build(&self, fitted: &[ObservedRegion], svs: &[StructuralVariant]) -> Vec<PurpleCopyNumber> {
        let tolerance = AlleleTolerance { purity: self.purity_adjuster.purity, targeted: self.targeted };
        let ed = ExtendDiploid { tolerance: &tolerance, min_tumor_count: self.min_tumor_ratio_count, centromere_min_tumor_count: self.min_tumor_ratio_count_at_centromere };
        let mut diploid: BTreeMap<Chromosome, Vec<CombinedRegion>> = BTreeMap::new();
        for c in Chromosome::ALL {
            let regions: Vec<ObservedRegion> = fitted.iter().filter(|x| Chromosome::from_name(&x.chromosome) == Some(c)).cloned().collect();
            let extended = extend_non_diploid(ed.extend_diploid(&regions));
            if !extended.is_empty() {
                diploid.insert(c, extended);
            }
        }
        let factory = LegPloidyFactory { average_copy_number: self.ploidy, average_read_depth: self.average_read_depth, purity_adjuster: self.purity_adjuster };
        sv_implied_copy_number(&factory, svs, &mut diploid);
        let mut out = Vec::new();
        for c in Chromosome::ALL {
            let simple: Vec<&StructuralVariant> = svs
                .iter()
                .filter(|x| {
                    x.end.as_ref().is_some_and(|e| {
                        Chromosome::contains(&x.start.chromosome)
                            && Chromosome::contains(&e.chromosome)
                            && x.start.chromosome == e.chromosome
                            && chr(&x.start.chromosome) == c
                    })
                })
                .collect();
            let baf = ExtendDiploidBaf::new(&simple);
            let mut regions = diploid.remove(&c).unwrap_or_default();
            extend_long_arm(&mut regions);
            populate_unknown(&mut regions, self.cobalt);
            baf.extend_baf(&mut regions);
            out.extend(to_copy_numbers(&regions));
        }
        out
    }
}

fn to_copy_numbers(regions: &[CombinedRegion]) -> Vec<PurpleCopyNumber> {
    let mut out: Vec<PurpleCopyNumber> = Vec::with_capacity(regions.len());
    for (i, region) in regions.iter().enumerate() {
        let mut start = region.start();
        if i > 0 {
            let first = &region.regions[0];
            let prev_end = out.last().unwrap().end;
            if first.germline_status == GermlineStatus::Diploid && first.min_start < first.max_start && first.min_start == prev_end + 1 {
                start = (first.min_start + first.max_start) / 2;
                out.last_mut().unwrap().end = start - 1;
            }
        }
        let trailing = if i + 1 < regions.len() { regions[i + 1].combined.support } else { SegmentSupport::Telomere };
        let r = &region.combined;
        out.push(PurpleCopyNumber {
            chromosome: r.chromosome.clone(),
            start,
            end: region.end(),
            baf_count: region.baf_count(),
            average_actual_baf: region.tumor_baf(),
            average_observed_baf: r.observed_baf,
            average_tumor_copy_number: region.tumor_copy_number(),
            depth_window_count: r.depth_window_count,
            segment_start_support: r.support,
            segment_end_support: trailing,
            method: region.method,
            gc_content: r.gc_content,
            min_start: r.min_start,
            max_start: r.max_start,
        });
    }
    out
}

/// PurpleCopyNumberFactory.validateCopyNumbers
pub fn validate_copy_numbers(copy_numbers: &[PurpleCopyNumber]) -> bool {
    let mut valid = true;
    for i in 1..copy_numbers.len() {
        let (c, p) = (&copy_numbers[i], &copy_numbers[i - 1]);
        if c.chromosome != p.chromosome {
            continue;
        }
        if c.start <= p.end {
            crate::error!("purple copy-number({}:{}-{}) overlaps previous({}-{})", c.chromosome, c.start, c.end, p.start, p.end);
            valid = false;
        }
    }
    valid
}


// ---------------------------------------------------------------------------------------------
// SV enrichment (CopyNumberEnrichedStructuralVariantFactory, SvLegCopyNumberChangeFactory)

/// SvCopyNumberChange
struct CopyNumberChange {
    down_scale: f64,
    up_scale: f64,
    down_offset: f64,
    up_offset: f64,
}

impl CopyNumberChange {
    fn new(legs: &[LegPloidy]) -> CopyNumberChange {
        let template = &legs[0];
        let difference = template.right_copy_number.unwrap_or(0.0) - template.left_copy_number.unwrap_or(0.0);
        let (mut up_count, mut down_count) = (0, 0);
        let (mut up, mut down) = (0.0, 0.0);
        let mut min_ploidy = f64::MAX;
        let mut have_min = false;
        for l in legs {
            let p = jmax(0.0, l.average_implied_ploidy);
            if l.leg.orientation == 1 {
                down_count += 1;
                down += p;
            } else {
                up_count += 1;
                up += p;
            }
            if less_than(p, min_ploidy) {
                min_ploidy = p;
                have_min = true;
            }
        }
        let down_scale = (up - difference) / (up + down);
        let up_scale = 1.0 - down_scale;
        let (down_offset, up_offset) = if up_count > 0 && down_count > 0 && have_min {
            let target_down = jmin(template.left_copy_number.unwrap_or(0.0), down);
            let target_up = jmin(template.right_copy_number.unwrap_or(0.0), up);
            let target = target_up + target_down;
            let explained = down_scale * down + up_scale * up;
            let offset = target - explained;
            let d = offset / (2 * down_count) as f64;
            (d, d * down_count as f64 / up_count as f64)
        } else {
            (0.0, 0.0)
        };
        CopyNumberChange { down_scale, up_scale, down_offset, up_offset }
    }
    fn change(&self, leg: &LegPloidy) -> f64 {
        let p = jmax(0.0, leg.average_implied_ploidy);
        if leg.leg.orientation == 1 { self.down_offset + self.down_scale * p } else { self.up_offset + self.up_scale * p }
    }
}

fn change_simple(orientation: i8, left: Option<f64>, right: Option<f64>) -> f64 {
    let (l, r) = (left.unwrap_or(0.0), right.unwrap_or(0.0));
    if orientation == 1 { l - r } else { r - l }
}

#[derive(Clone, Debug, Default)]
pub struct EnrichedLeg {
    pub adjusted_allele_frequency: Option<f64>,
    pub adjusted_copy_number: Option<f64>,
    pub adjusted_copy_number_change: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct EnrichedSv {
    pub index: usize,
    pub start: EnrichedLeg,
    pub end: Option<EnrichedLeg>,
    pub junction_copy_number: Option<f64>,
}

/// CopyNumberEnrichedStructuralVariantFactory.enrich
pub fn enrich_svs(purity_adjuster: &PurityAdjuster, copy_numbers: &BTreeMap<Chromosome, Vec<PurpleCopyNumber>>, variants: &[StructuralVariant]) -> Vec<EnrichedSv> {
    let selector = CopyNumberSelector::from_copy_numbers(copy_numbers);
    let factory = LegPloidyFactory { average_copy_number: 0.0, average_read_depth: 0, purity_adjuster };
    // SvLegCopyNumberChangeFactory: variants with a leg at a position shared with another
    let mut at: Vec<(CnaPos, usize)> = Vec::new();
    for (i, v) in variants.iter().enumerate() {
        at.push((cna_pos(&v.start), i));
        if let Some(e) = &v.end {
            at.push((cna_pos(e), i));
        }
    }
    let mut shared: Vec<usize> = Vec::new();
    for (k, _) in &at {
        let here: Vec<usize> = at.iter().filter(|(k2, _)| k2 == k).map(|(_, i)| *i).collect();
        if here.len() > 1 {
            for i in here {
                if !shared.contains(&i) {
                    shared.push(i);
                }
            }
        }
    }
    shared.sort();
    let mut ploidy_map: Vec<(CnaPos, Vec<LegPloidy>)> = Vec::new();
    for &i in &shared {
        for p in factory.create_variant(&variants[i], &selector) {
            let k = cna_pos(&p.leg);
            match ploidy_map.iter_mut().find(|(k2, _)| *k2 == k) {
                Some(e) => e.1.push(p),
                None => ploidy_map.push((k, vec![p])),
            }
        }
    }
    let changes: HashMap<CnaPos, CopyNumberChange> =
        ploidy_map.iter().filter(|(_, v)| v.len() > 1).map(|(k, v)| (k.clone(), CopyNumberChange::new(v))).collect();
    let change_of_ploidy = |p: &LegPloidy| match changes.get(&cna_pos(&p.leg)) {
        Some(c) => c.change(p),
        None => change_simple(p.leg.orientation, p.left_copy_number, p.right_copy_number),
    };
    let leg_cn = |leg: &SvLeg| {
        let left = selector.at(&leg.chromosome, leg.cna_position() - 1).map(|x| jmax(0.0, x));
        let right = selector.at(&leg.chromosome, leg.cna_position()).map(|x| jmax(0.0, x));
        (left, right)
    };
    let mut out = Vec::with_capacity(variants.len());
    for (i, v) in variants.iter().enumerate() {
        let mut start = EnrichedLeg::default();
        let mut end = v.end.as_ref().map(|_| EnrichedLeg::default());
        let mut junction = None;
        let ploidies = factory.create_variant(v, &selector);
        if !ploidies.is_empty() {
            junction = Some(ploidies[0].average_implied_ploidy);
            let sp = &ploidies[0];
            start.adjusted_allele_frequency = Some(sp.adjusted_vaf);
            start.adjusted_copy_number = Some(sp.adjusted_copy_number());
            start.adjusted_copy_number_change = Some(change_of_ploidy(sp));
            if let Some(ep) = ploidies.get(1) {
                // Java: the end builder may be null (no end leg) where a second ploidy exists; it is not
                let e = end.get_or_insert_with(EnrichedLeg::default);
                e.adjusted_allele_frequency = Some(ep.adjusted_vaf);
                e.adjusted_copy_number = Some(ep.adjusted_copy_number());
                e.adjusted_copy_number_change = Some(change_of_ploidy(ep));
            }
        } else {
            let (l, r) = leg_cn(&v.start);
            start.adjusted_copy_number = Some(if v.start.orientation == 1 { l.unwrap_or(0.0) } else { r.unwrap_or(0.0) });
            let simple_change = if changes.contains_key(&cna_pos(&v.start)) { None } else { Some(change_simple(v.start.orientation, l, r)) };
            start.adjusted_copy_number_change = simple_change;
            if v.sv_type == SvType::Inf {
                junction = simple_change;
            }
            if let Some(el) = &v.end {
                let (l, r) = leg_cn(el);
                let e = end.as_mut().unwrap();
                e.adjusted_copy_number = Some(if el.orientation == 1 { l.unwrap_or(0.0) } else { r.unwrap_or(0.0) });
                e.adjusted_copy_number_change = if changes.contains_key(&cna_pos(el)) { None } else { Some(change_simple(el.orientation, l, r)) };
            }
        }
        out.push(EnrichedSv { index: i, start, end, junction_copy_number: junction });
    }
    out
}

// ---------------------------------------------------------------------------------------------
// ChromosomeCopyNumbers and LohCalcs

#[derive(Clone, Debug, PartialEq)]
pub struct ChrArmCopyNumber {
    pub chromosome: Chromosome,
    pub arm: hmf_common::segmentation::Arm,
    pub mean: f64,
    pub median: f64,
    pub min: f64,
    pub max: f64,
}

/// ChromosomeCopyNumbers.data
pub fn chromosome_arm_copy_numbers(copy_numbers: &[PurpleCopyNumber], v38: bool) -> Vec<ChrArmCopyNumber> {
    let locator = hmf_common::segmentation::default_locator(v38);
    arm_copy_numbers_with(copy_numbers, |chromosome, position| locator(chromosome, position).map(|a| a.arm).unwrap_or(hmf_common::segmentation::Arm::Q))
}

/// ChromosomeCopyNumbers with a ChrArmLocator
pub fn arm_copy_numbers_with(copy_numbers: &[PurpleCopyNumber], arm_of: impl Fn(&str, i32) -> hmf_common::segmentation::Arm) -> Vec<ChrArmCopyNumber> {
    use hmf_common::segmentation::Arm;
    let stats = |list: &[&PurpleCopyNumber]| -> (f64, f64, f64, f64) {
        let (mut total_len, mut total_weight) = (0.0, 0.0);
        let mut min = f64::MAX;
        // PU6: Java starts the maximum at Double.MIN_VALUE, the smallest positive double
        let mut max = if crate::compat::java("PU6") { f64::from_bits(1) } else { f64::NEG_INFINITY };
        let mut sorted: Vec<&PurpleCopyNumber> = Vec::new();
        for c in list {
            let len = c.length() as f64;
            total_len += len;
            total_weight += len * c.average_tumor_copy_number;
            let pos = sorted.partition_point(|x| {
                dcmp(x.average_tumor_copy_number, c.average_tumor_copy_number).then((x.start as f64).total_cmp(&(c.start as f64))) == std::cmp::Ordering::Less
            });
            let dup = sorted.get(pos).is_some_and(|x| {
                dcmp(x.average_tumor_copy_number, c.average_tumor_copy_number) == std::cmp::Ordering::Equal && x.start == c.start
            });
            if !dup {
                sorted.insert(pos, c);
            }
            min = jmin(min, c.average_tumor_copy_number);
            max = jmax(max, c.average_tumor_copy_number);
        }
        let mean = total_weight / total_len;
        let half = (total_len / 2.0) as i32;
        let mut median = 0.0;
        let mut running = 0;
        for c in &sorted {
            running += c.length();
            if running > half {
                median = c.average_tumor_copy_number;
                break;
            }
        }
        (mean, median, min, max)
    };
    let mut result = Vec::new();
    let mut current: Vec<&PurpleCopyNumber> = Vec::new();
    let mut current_chr: Option<Chromosome> = None;
    let mut current_arm: Option<Arm> = None;
    for c in copy_numbers {
        let ch = Chromosome::from_name(&c.chromosome);
        let arm = arm_of(&c.chromosome, c.start);
        if ch != current_chr || Some(arm) != current_arm {
            if !current.is_empty() {
                let (mean, median, min, max) = stats(&current);
                result.push(ChrArmCopyNumber { chromosome: current_chr.unwrap(), arm: current_arm.unwrap(), mean, median, min, max });
            }
            current_chr = ch;
            current_arm = Some(arm);
            current = Vec::new();
        }
        current.push(c);
    }
    if !current.is_empty() {
        let (mean, median, min, max) = stats(&current);
        result.push(ChrArmCopyNumber { chromosome: current_chr.unwrap(), arm: current_arm.unwrap(), mean, median, min, max });
    }
    result.into_iter().filter(|a| !short_arm(a.chromosome) || a.arm == Arm::Q).collect()
}

/// ChrArmCopyNumbersFile.write
pub fn write_chromosome_arms(path: &str, arms: &[ChrArmCopyNumber]) -> anyhow::Result<()> {
    let mut out = String::from("chromosome\tarm\tmeanCopyNumber\tmedianCopyNumber\tminCopyNumber\tmaxCopyNumber\n");
    for a in arms {
        let arm = match a.arm {
            hmf_common::segmentation::Arm::P => "P",
            hmf_common::segmentation::Arm::Q => "Q",
        };
        out.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", a.chromosome.short_name(), arm, f4(a.mean), f4(a.median), f4(a.min), f4(a.max)));
    }
    std::fs::write(path, out).map_err(|e| anyhow::anyhow!("{path}: {e}"))
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LohCalcData {
    pub loh_base_count: i64,
    pub total_base_count: i64,
    pub segments: i32,
}

impl LohCalcData {
    pub fn loh_percent(&self) -> f64 {
        if self.total_base_count > 0 { self.loh_base_count as f64 / self.total_base_count as f64 } else { 0.0 }
    }
}

/// LohCalcs.calcLohSegments
pub fn calc_loh_segments(chromosome: &str, copy_numbers: &[PurpleCopyNumber], min_loh_length: f64, allow_whole_arm: bool, count_by_arm: bool) -> LohCalcData {
    let c = chr(chromosome);
    if !is_autosome(c) {
        return LohCalcData::default();
    }
    let (mut loh_count, mut loh_bases, mut total_bases) = (0i32, 0i32, 0i32);
    let mut in_loh = false;
    let mut loh_start = 0;
    let mut before_centromere = false;
    let mut has_non_loh = false;
    for i in 0..copy_numbers.len() {
        let cn = &copy_numbers[i];
        if cn.segment_start_support == SegmentSupport::Telomere {
            in_loh = false;
            before_centromere = true;
        } else if cn.segment_start_support == SegmentSupport::Centromere {
            before_centromere = false;
        }
        if before_centromere && short_arm(c) {
            continue;
        }
        total_bases += cn.length();
        let next = copy_numbers.get(i + 1);
        let is_loh = cn.minor_allele_copy_number() < 0.5;
        has_non_loh |= !is_loh;
        let next_loh = next.is_some_and(|n| n.minor_allele_copy_number() < 0.5);
        if !in_loh && is_loh {
            in_loh = true;
            loh_start = cn.start;
        }
        let end_of_arm = cn.segment_end_support == SegmentSupport::Telomere || (count_by_arm && cn.segment_end_support == SegmentSupport::Centromere);
        let mut check_end = false;
        if in_loh {
            if end_of_arm {
                check_end = true;
            } else if !next_loh {
                check_end = !next.is_some_and(|n| n.length() <= 1000);
            }
        }
        if check_end {
            let length = cn.end - loh_start + 1;
            if length as f64 > min_loh_length {
                loh_count += 1;
                loh_bases += length;
            }
            in_loh = false;
        }
    }
    if allow_whole_arm {
        if loh_count == 0 {
            return LohCalcData::default();
        }
        if !has_non_loh {
            loh_count = 1;
        }
        return LohCalcData { loh_base_count: loh_bases as i64, total_base_count: total_bases as i64, segments: loh_count };
    }
    if has_non_loh {
        LohCalcData { loh_base_count: loh_bases as i64, total_base_count: total_bases as i64, segments: loh_count }
    } else {
        LohCalcData::default()
    }
}

/// LohCalcs.calcLohData
pub fn calc_loh_data(copy_numbers: &[PurpleCopyNumber]) -> LohCalcData {
    let mut total = LohCalcData::default();
    for (chromosome, list) in chromosome_map(copy_numbers) {
        let d = calc_loh_segments(&chromosome, &list, 1000.0, false, true);
        total.loh_base_count += d.loh_base_count;
        total.total_base_count += d.total_base_count;
        total.segments += d.segments;
    }
    total
}

#[cfg(test)]
mod chromosome_arm_tests {
    use super::*;
    use hmf_common::segmentation::Arm;

    const EPSILON: f64 = 1e-9;

    fn pcn(chromosome: &str, start: i32, end: i32, cn: f64, support: SegmentSupport) -> PurpleCopyNumber {
        crate::test_support::copy_number(chromosome, start, end, cn, support, SegmentSupport::None)
    }

    /// DummyChrArmLocator: P up to the given centre (default 1e9)
    fn arms(cns: &[PurpleCopyNumber], centres: &[(&str, i32)]) -> Vec<ChrArmCopyNumber> {
        let centres: Vec<(String, i32)> = centres.iter().map(|(c, p)| (c.to_string(), *p)).collect();
        arm_copy_numbers_with(cns, move |chromosome, position| {
            let centre = centres.iter().find(|(c, _)| c == chromosome).map(|(_, p)| *p).unwrap_or(1_000_000_000);
            if position <= centre { Arm::P } else { Arm::Q }
        })
    }

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= EPSILON, "{a} != {b}");
    }

    #[test]
    fn one_purple_copy_number() {
        let r = arms(&[pcn("1", 10_000_000, 100_000_000, 0.55, SegmentSupport::Telomere)], &[]);
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].chromosome, r[0].arm), (chr("1"), Arm::P));
        for v in [r[0].mean, r[0].median, r[0].max, r[0].min] {
            close(v, 0.55);
        }
    }

    #[test]
    fn one_arm() {
        let n = SegmentSupport::None;
        let cns = [
            pcn("1", 1, 10_000_000, 0.5, SegmentSupport::Telomere),
            pcn("1", 10_000_001, 20_000_000, 0.5, n),
            pcn("1", 20_000_001, 30_000_000, 0.5, n),
            pcn("1", 30_000_001, 40_000_000, 0.7, n),
            pcn("1", 40_000_001, 50_000_000, 0.8, n),
        ];
        let r = arms(&cns, &[]);
        assert_eq!(r.len(), 1);
        close(r[0].mean, 0.6);
        close(r[0].median, 0.5);
        close(r[0].max, 0.8);
        close(r[0].min, 0.5);
    }

    #[test]
    fn mean_and_median_weighted_by_length() {
        let n = SegmentSupport::None;
        let cns = [pcn("1", 1, 100, 0.4, n), pcn("1", 101, 400, 0.5, n), pcn("1", 401, 500, 0.6, n), pcn("1", 501, 1000, 0.7, n), pcn("1", 1001, 2000, 0.8, n)];
        let r = arms(&cns, &[]);
        close(r[0].mean, 0.7);
        close(r[0].median, 0.8);
        let cns = [
            pcn("1", 1, 100, 0.4, n),
            pcn("1", 101, 400, 0.5, n),
            pcn("1", 401, 600, 0.4, n),
            pcn("1", 601, 1000, 0.7, n),
            pcn("1", 1001, 1200, 0.3, n),
            pcn("1", 1201, 1600, 0.7, n),
            pcn("1", 1601, 1700, 0.6, n),
            pcn("1", 1701, 1900, 0.8, n),
        ];
        close(arms(&cns, &[])[0].median, 0.7);
    }

    #[test]
    fn multiple_arms() {
        let (t, c, n) = (SegmentSupport::Telomere, SegmentSupport::Centromere, SegmentSupport::None);
        let cns = [
            pcn("1", 1, 50_000_000, 0.5, t),
            pcn("1", 50_000_001, 100_000_000, 0.7, n),
            pcn("1", 160_000_001, 200_000_000, 0.4, c),
            pcn("1", 200_000_001, 240_000_000, 0.6, n),
            pcn("2", 1, 40_000_000, 0.55, t),
            pcn("2", 40_000_001, 80_000_000, 0.75, n),
            pcn("2", 100_000_001, 150_000_000, 0.45, c),
            pcn("2", 150_000_001, 200_000_000, 0.65, n),
        ];
        let r = arms(&cns, &[("1", 150_000_000), ("2", 90_000_000)]);
        assert_eq!(r.len(), 4);
        let expect = [("1", Arm::P, 0.6, 0.7, 0.7, 0.5), ("1", Arm::Q, 0.5, 0.6, 0.6, 0.4), ("2", Arm::P, 0.65, 0.75, 0.75, 0.55), ("2", Arm::Q, 0.55, 0.65, 0.65, 0.45)];
        for (a, (c, arm, mean, median, max, min)) in r.iter().zip(expect) {
            assert_eq!((a.chromosome, a.arm), (chr(c), arm));
            close(a.mean, mean);
            close(a.median, median);
            close(a.max, max);
            close(a.min, min);
        }
    }

    #[test]
    fn non_consecutive_arms() {
        let (t, c, n) = (SegmentSupport::Telomere, SegmentSupport::Centromere, SegmentSupport::None);
        let cns = [
            pcn("1", 1, 50_000_000, 0.5, t),
            pcn("1", 50_000_001, 100_000_000, 0.7, n),
            pcn("1", 100_000_001, 120_000_000, 0.4, n),
            pcn("2", 80_000_001, 90_000_000, 0.55, c),
            pcn("2", 90_000_001, 120_000_000, 0.55, n),
            pcn("3", 60_000_001, 80_000_000, 0.75, c),
            pcn("3", 80_000_001, 100_000_000, 0.75, n),
        ];
        let r = arms(&cns, &[("1", 150_000_000), ("2", 80_000_000), ("3", 60_000_000)]);
        let got: Vec<(Chromosome, Arm)> = r.iter().map(|a| (a.chromosome, a.arm)).collect();
        assert_eq!(got, vec![(chr("1"), Arm::P), (chr("2"), Arm::Q), (chr("3"), Arm::Q)]);
    }

    #[test]
    fn centromeric_copy_numbers() {
        let (t, c, n) = (SegmentSupport::Telomere, SegmentSupport::Centromere, SegmentSupport::None);
        let cns = [
            pcn("10", 1, 50_000_000, 0.5, t),
            pcn("10", 50_000_001, 60_000_000, 0.7, n),
            pcn("10", 70_000_001, 71_000_000, 0.4, c),
            pcn("10", 71_000_001, 72_000_000, 0.3, c),
            pcn("10", 72_000_001, 100_000_000, 0.4, n),
            pcn("11", 1, 40_000_000, 0.55, t),
            pcn("11", 50_000_001, 60_000_000, 0.55, c),
        ];
        let r = arms(&cns, &[("10", 65_000_000), ("11", 45_000_000)]);
        assert_eq!(r.len(), 4);
        assert_eq!((r[1].chromosome, r[1].arm), (chr("10"), Arm::Q));
        close(r[1].min, 0.3);
        let cns = [pcn("10", 1, 50_000_000, 0.5, t), pcn("10", 50_000_001, 60_000_000, 0.7, n), pcn("11", 40_000_001, 60_000_000, 0.55, t), pcn("12", 50_000_001, 60_000_000, 0.55, c)];
        let r = arms(&cns, &[("10", 70_000_000), ("11", 70_000_000), ("12", 50_000_000)]);
        let got: Vec<(Chromosome, Arm)> = r.iter().map(|a| (a.chromosome, a.arm)).collect();
        assert_eq!(got, vec![(chr("10"), Arm::P), (chr("11"), Arm::P), (chr("12"), Arm::Q)]);
    }

    #[test]
    fn exclude_non_reportable_arms() {
        let mut cns = Vec::new();
        for c in Chromosome::ALL {
            cns.push(pcn(c.short_name(), 1, 100, 0.5, SegmentSupport::None));
            cns.push(pcn(c.short_name(), 1, 100, 0.5, SegmentSupport::Centromere));
        }
        // the short arms of the acrocentric chromosomes are left out
        for a in arms(&cns, &[]) {
            assert!(!short_arm(a.chromosome) || a.arm == Arm::Q);
        }
    }
}

#[cfg(test)]
mod extend_diploid_tests {
    use super::*;
    use crate::test_support::default_fitted_region;

    const MIN_TUMOR_COUNT: i32 = 30;
    const MIN_TUMOR_COUNT_AT_CENTROMERE: i32 = 50;

    fn extend(regions: &[ObservedRegion]) -> Vec<CombinedRegion> {
        // AlleleTolerance of a pure tumor (PurityAdjuster of purity 1)
        let tolerance = AlleleTolerance { purity: 1.0, targeted: false };
        ExtendDiploid { tolerance: &tolerance, min_tumor_count: MIN_TUMOR_COUNT, centromere_min_tumor_count: MIN_TUMOR_COUNT_AT_CENTROMERE }.extend_diploid(regions)
    }

    fn region(start: i32, end: i32, cn: f64, ratio_count: i32, baf_count: i32, support: SegmentSupport) -> ObservedRegion {
        let mut r = default_fitted_region("1", start, end);
        r.germline_status = GermlineStatus::Diploid;
        r.tumor_copy_number = cn;
        r.ref_normalised_copy_number = cn;
        r.observed_normal_ratio = 0.5;
        r.baf_count = baf_count;
        r.support = support;
        r.depth_window_count = ratio_count;
        r
    }

    fn somatic(start: i32, end: i32, cn: f64, baf_count: i32, support: SegmentSupport) -> ObservedRegion {
        region(start, end, cn, MIN_TUMOR_COUNT_AT_CENTROMERE, baf_count, support)
    }

    fn dubious(start: i32, end: i32, cn: f64, ratio_count: i32) -> ObservedRegion {
        region(start, end, cn, ratio_count, 1, SegmentSupport::None)
    }

    fn fitted(start: i32, end: i32, cn: f64, status: GermlineStatus, support: SegmentSupport) -> ObservedRegion {
        let mut r = default_fitted_region("1", start, end);
        r.germline_status = status;
        r.tumor_copy_number = cn;
        r.ref_normalised_copy_number = cn;
        r.observed_normal_ratio = 1.1;
        r.baf_count = 50;
        r.support = support;
        r
    }

    fn check(r: &CombinedRegion, start: i32, end: i32, cn: f64) {
        assert_eq!((r.start(), r.end()), (start, end));
        assert!((r.tumor_copy_number() - cn).abs() <= 1e-10, "{} != {cn}", r.tumor_copy_number());
    }

    // ExtendDiploidTest
    #[test]
    fn favour_tumor_ratio_count_over_length() {
        let r = extend(&[region(1, 50000, 2.0, 10, 0, SegmentSupport::None), region(50001, 60000, 3.0, 20, 0, SegmentSupport::None)]);
        assert_eq!(r.len(), 1);
        check(&r[0], 1, 60000, 3.0);
    }

    #[test]
    fn centromere_is_dubious() {
        let r = extend(&[
            somatic(1, 10000, 1.0, 1, SegmentSupport::Telomere),
            region(10001, 20000, 4.0, 8, 0, SegmentSupport::Centromere),
            region(20001, 30000, 4.0, 9, 0, SegmentSupport::None),
            region(30001, 40000, 4.0, 10, 0, SegmentSupport::None),
        ]);
        assert_eq!(r.len(), 2);
        check(&r[0], 1, 10000, 1.0);
        check(&r[1], 10001, 40000, 4.0);
    }

    #[test]
    fn dubious_region_included_or_excluded() {
        use SegmentSupport::*;
        for left in [50, 30] {
            let r = extend(&[somatic(1, 10000, 3.0, left, Bnd), dubious(10001, 20000, 2.0, 10), somatic(20001, 30000, 3.0, 40, None)]);
            assert_eq!(r.len(), 1);
            check(&r[0], 1, 30000, 3.0);
        }
        let r = extend(&[somatic(1, 10000, 3.0, 50, Bnd), dubious(10001, 20000, 2.0, 10), somatic(20001, 30000, 3.0, 40, Bnd)]);
        assert_eq!(r.len(), 3);
    }

    #[test]
    fn multiple_dubious_regions() {
        use SegmentSupport::*;
        let r = extend(&[somatic(1, 10000, 3.0, 50, Bnd), dubious(10001, 11000, 2.0, 20), dubious(11001, 12000, 2.0, 9), dubious(12001, 20000, 2.0, 20), somatic(20001, 30000, 3.0, 40, None)]);
        assert_eq!(r.len(), 1);
        check(&r[0], 1, 30000, 3.0);
        let r = extend(&[somatic(1, 10000, 3.0, 50, Bnd), dubious(10001, 11000, 2.0, 20), dubious(11001, 12000, 2.0, 10), dubious(12001, 20000, 2.0, 20), somatic(20001, 30000, 3.0, 40, None)]);
        assert_eq!(r.len(), 3);
        check(&r[0], 1, 10000, 3.0);
        check(&r[1], 10001, 20000, 2.0);
        check(&r[2], 20001, 30000, 3.0);
    }

    #[test]
    fn towards_centromere() {
        use SegmentSupport::*;
        let r = extend(&[somatic(1, 10000, 3.0, 50, Bnd), dubious(10001, 11000, 2.0, 11), dubious(11001, 12000, 2.0, 11), dubious(12001, 20000, 2.0, 11), somatic(20001, 30000, 3.0, 40, Centromere)]);
        assert_eq!(r.len(), 2);
        check(&r[0], 1, 20000, 3.0);
        check(&r[1], 20001, 30000, 3.0);
        let r = extend(&[somatic(1, 10000, 3.0, 50, Bnd), dubious(10001, 11000, 2.0, 21), dubious(11001, 12000, 2.0, 21), dubious(12001, 20000, 2.0, 21), somatic(20001, 30000, 3.0, 40, Centromere)]);
        assert_eq!(r.len(), 3);
        check(&r[0], 1, 10000, 3.0);
        check(&r[1], 10001, 20000, 2.0);
        check(&r[2], 20001, 30000, 3.0);
    }

    #[test]
    fn invalid_germline() {
        use SegmentSupport::*;
        let r = extend(&[fitted(1, 10000, 3.0, GermlineStatus::Diploid, Bnd), fitted(10001, 20000, 4.0, GermlineStatus::Amplification, None)]);
        assert_eq!(r.len(), 1);
        check(&r[0], 1, 20000, 3.0);
        let r = extend(&[fitted(1, 10000, 3.0, GermlineStatus::Diploid, Bnd), fitted(10001, 20000, 4.0, GermlineStatus::Amplification, Bnd)]);
        assert_eq!(r.len(), 2);
        check(&r[0], 1, 10000, 3.0);
        check(&r[1], 10001, 20000, 4.0);
    }
}

#[cfg(test)]
mod combined_region_tests {
    use super::*;
    use crate::test_support::default_fitted_region;

    fn create(start: i32, end: i32, baf_count: i32, baf: f64, cn: f64) -> ObservedRegion {
        let mut r = default_fitted_region("1", start, end);
        r.baf_count = baf_count;
        r.observed_baf = baf;
        r.tumor_copy_number = cn;
        r.ref_normalised_copy_number = cn;
        r
    }

    fn averages(r: &CombinedRegion, baf: f64, cn: f64) {
        assert!((r.combined.observed_baf - baf).abs() <= 1e-10, "baf {} != {baf}", r.combined.observed_baf);
        assert!((r.combined.tumor_copy_number - cn).abs() <= 1e-10, "cn {} != {cn}", r.combined.tumor_copy_number);
    }

    // CombinedRegionTest (new CombinedRegion(region) is BAF weighted)
    #[test]
    fn depth_window_count_summation_only_applies_to_somatic() {
        let mut somatic = default_fitted_region("1", 2001, 3000);
        somatic.depth_window_count = 2;
        somatic.germline_status = GermlineStatus::Diploid;
        let mut region = CombinedRegion::new(true, &somatic);
        assert_eq!(region.combined.depth_window_count, 2);
        for (s, e) in [(1, 1000), (1001, 2000)] {
            let mut amp = default_fitted_region("1", s, e);
            amp.depth_window_count = 2;
            amp.germline_status = GermlineStatus::Amplification;
            region.extend(&amp);
            assert_eq!(region.combined.depth_window_count, 2);
        }
    }

    #[test]
    fn average_on_length_for_non_zero_ratio() {
        let mut r = CombinedRegion::new(true, &create(1, 100, 0, 0.0, 3.0));
        averages(&r, 0.0, 3.0);
        r.extend_with_weighted_average(&create(101, 200, 0, 0.0, 0.0));
        averages(&r, 0.0, 3.0);
    }

    #[test]
    fn do_not_include_zero_copy_number() {
        let mut r = CombinedRegion::new(true, &create(1, 100, 200, 0.5, 0.0));
        averages(&r, 0.5, 0.0);
        r.extend_with_weighted_average(&create(201, 300, 200, 1.0, 2.0));
        averages(&r, 0.75, 2.0);
    }
}
