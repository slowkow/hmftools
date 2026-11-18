//! Soft filters (VariantFilters, StrandBiasCalcs, FilterUtils) and the SageVariant record.

use crate::candidate::Candidate;
use crate::config::Config;
use crate::constants::*;
use crate::evidence::{ReadContextCounter, StrandBias};
use crate::stats::binomial_cdf;
use crate::types::*;
use std::collections::HashSet;
use std::sync::OnceLock;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub enum SoftFilter {
    MinTumorQual,
    MinMapQualFactor,
    MinTumorVaf,
    MinTumorSupport,
    MinGermlineDepth,
    MaxGermlineVaf,
    MaxGermlineRelativeQual,
    MaxGermlineAltSupport,
    MinAvgBaseQuality,
    MinAvgHpQual,
    MinAvgT0Qual,
    FragmentStrandBias,
    ReadStrandBias,
    FragmentCoords,
    MaxEdgeDistance,
    RealignedFreq,
    FragmentLength,
    MapQualRefAltDifference,
    Jitter,
    DedupMnv,
    DedupMixedGermlineSomatic,
    DedupSnvMnv,
    DedupIndel,
    DedupMatch,
}

pub const ALL_FILTERS: [(SoftFilter, &str, &str); 24] = [
    (SoftFilter::MinTumorQual, "minTumorQual", "Insufficient tumor quality"),
    (SoftFilter::MinMapQualFactor, "minMQF", "Insufficient map qual factor"),
    (SoftFilter::MinTumorVaf, "minTumorVAF", "Insufficient tumor VAF"),
    (SoftFilter::MinTumorSupport, "minTumorSupport", "Insufficient tumor strong support"),
    (SoftFilter::MinGermlineDepth, "minGermlineDepth", "Insufficient germline depth"),
    (SoftFilter::MaxGermlineVaf, "maxGermlineVAF", "Excess germline VAF"),
    (SoftFilter::MaxGermlineRelativeQual, "maxGermlineRelQual", "Excess germline relative qual"),
    (SoftFilter::MaxGermlineAltSupport, "maxGermlineAltSupport", "Excess germline alt support"),
    (SoftFilter::MinAvgBaseQuality, "minAvgRawBaseQual", "Variant average raw base quality below limit"),
    (SoftFilter::MinAvgHpQual, "minAvgHpQual", "Average base quality within core homopolymer below limit"),
    (SoftFilter::MinAvgT0Qual, "minAvgT0Qual", "Average T0 quality within core homopolymer below limit"),
    (SoftFilter::FragmentStrandBias, "fragmentStrandBias", "Variant exceeds fragment strand bias limit"),
    (SoftFilter::ReadStrandBias, "readStrandBias", "Variant exceeds read strand bias limit"),
    (SoftFilter::FragmentCoords, "minFragmentCoords", "Insufficient fragment coordinate variation"),
    (SoftFilter::MaxEdgeDistance, "maxEdgeDistance", "Variant close to read edge"),
    (SoftFilter::RealignedFreq, "realignFrequency", "Support is predominantly realigned"),
    (SoftFilter::FragmentLength, "fragmentLength", "Alt max fragment length well below average non-alt"),
    (SoftFilter::MapQualRefAltDifference, "mapQualRefAltDiff", "Alt support map qual well below ref support"),
    (SoftFilter::Jitter, "jitter", "Jitter filter"),
    (SoftFilter::DedupMnv, "dedupMnv", "Filter MNV duplicate"),
    (SoftFilter::DedupMixedGermlineSomatic, "dedupMixedGermlineSomatic", "Variant duplicate mixed somatic/germline"),
    (SoftFilter::DedupSnvMnv, "dedupSnvMnv", "Variant duplicate MNV vs SNV"),
    (SoftFilter::DedupIndel, "dedupIndel", "Variant duplicate SNV/MNV vs INDEL"),
    (SoftFilter::DedupMatch, "dedupMatch", "Variant duplicate with different read contexts"),
];

impl SoftFilter {
    pub fn name(self) -> &'static str {
        ALL_FILTERS.iter().find(|(f, _, _)| *f == self).unwrap().1
    }
    pub fn is_germline(self) -> bool {
        matches!(
            self,
            SoftFilter::MinGermlineDepth
                | SoftFilter::MaxGermlineVaf
                | SoftFilter::MaxGermlineRelativeQual
                | SoftFilter::MaxGermlineAltSupport
        )
    }
    pub fn is_tumor(self) -> bool {
        !self.is_germline()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SoftFilterConfig {
    pub qual_p_score: f64,
    pub map_qual_factor: f64,
    pub min_tumor_vaf: f64,
    pub min_germline_coverage: i32,
    pub min_germline_coverage_allosome: i32,
    pub max_germline_vaf: f64,
}

pub const HOTSPOT_FILTER: SoftFilterConfig = SoftFilterConfig {
    qual_p_score: 1e-2,
    map_qual_factor: -6.0,
    min_tumor_vaf: 0.01,
    min_germline_coverage: 0,
    min_germline_coverage_allosome: 0,
    max_germline_vaf: 0.05,
};
pub const PANEL_FILTER: SoftFilterConfig = SoftFilterConfig { qual_p_score: 1e-5, min_tumor_vaf: 0.02, ..HOTSPOT_FILTER };
pub const HIGH_CONFIDENCE_FILTER: SoftFilterConfig = SoftFilterConfig {
    qual_p_score: 1e-8,
    map_qual_factor: 0.0,
    min_tumor_vaf: 0.025,
    min_germline_coverage: 10,
    min_germline_coverage_allosome: 6,
    max_germline_vaf: 0.05,
};
pub const LOW_CONFIDENCE_FILTER: SoftFilterConfig = SoftFilterConfig { qual_p_score: 1e-14, ..HIGH_CONFIDENCE_FILTER };

pub fn tier_config(tier: Tier) -> &'static SoftFilterConfig {
    match tier {
        Tier::Hotspot => &HOTSPOT_FILTER,
        Tier::Panel => &PANEL_FILTER,
        Tier::HighConfidence => &HIGH_CONFIDENCE_FILTER,
        Tier::LowConfidence => &LOW_CONFIDENCE_FILTER,
    }
}

pub struct SageVariant {
    pub candidate: Candidate,
    pub ref_counters: Vec<ReadContextCounter>,
    pub tumor_counters: Vec<ReadContextCounter>,
    pub filters: HashSet<SoftFilter>,
    pub mixed_impact: i32,
    pub near_indel: bool,
    pub near_hotspot: bool,
    pub near_multi_base_indel: bool,
}

impl SageVariant {
    pub fn new(candidate: Candidate, ref_counters: Vec<ReadContextCounter>, tumor_counters: Vec<ReadContextCounter>) -> Self {
        SageVariant {
            candidate,
            ref_counters,
            tumor_counters,
            filters: HashSet::new(),
            mixed_impact: 0,
            near_indel: false,
            near_hotspot: false,
            near_multi_base_indel: false,
        }
    }
    pub fn variant(&self) -> &SimpleVariant {
        self.candidate.variant()
    }
    pub fn position(&self) -> i32 {
        self.candidate.position()
    }
    pub fn tier(&self) -> Tier {
        self.candidate.tier
    }
    pub fn is_passing(&self) -> bool {
        self.filters.is_empty()
    }
    pub fn is_indel(&self) -> bool {
        self.variant().is_indel()
    }
    pub fn local_phase_sets(&self) -> &[i32] {
        &self.tumor_counters[0].local_phase_sets
    }
    pub fn has_local_phase_sets(&self) -> bool {
        !self.tumor_counters[0].local_phase_sets.is_empty()
    }
    pub fn has_matching_lps(&self, other: &[i32]) -> bool {
        self.local_phase_sets().iter().any(|x| other.contains(x))
    }
    pub fn has_lps(&self, lps: i32) -> bool {
        self.local_phase_sets().contains(&lps)
    }
    pub fn lps_read_count(&self, lps: i32) -> i32 {
        let t = &self.tumor_counters[0];
        t.local_phase_sets.iter().position(|&x| x == lps).map_or(0, |i| t.lps_counts[i])
    }
    pub fn remove_lps(&mut self, lps: i32) {
        let t = &mut self.tumor_counters[0];
        if let Some(i) = t.local_phase_sets.iter().position(|&x| x == lps) {
            t.local_phase_sets.remove(i);
            t.lps_counts.remove(i);
        }
    }
    pub fn total_quality(&self) -> i32 {
        self.tumor_counters.iter().map(|x| jround(x.tumor_quality()) as i32).sum()
    }
}

// Doubles, with EPSILON = 1e-10
fn d_less_than(v: f64, r: f64) -> bool {
    v - r < -1e-10
}
fn d_less_or_equal(v: f64, r: f64) -> bool {
    v - r < 1e-10
}
fn d_greater_than(v: f64, r: f64) -> bool {
    v - r > 1e-10
}
fn d_greater_or_equal(v: f64, r: f64) -> bool {
    v - r > -1e-10
}

const STRAND_BIAS_CHECK_THRESHOLD: f64 = 0.15;
const STRAND_BIAS_NON_ALT_MIN_DEPTH: i32 = 5;
const STRAND_BIAS_NON_ALT_MIN_BIAS: f64 = 0.25;

/// StrandBiasCalcs probability cache, built once.
fn strand_bias_values() -> &'static [f64] {
    static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
    TABLE.get_or_init(|| {
        const MAX_AD: usize = 1000;
        const PROB: f64 = 0.0005;
        let mut values = vec![0.0; MAX_AD + 1];
        let mut current_observed = 0;
        for depth in 1..=MAX_AD as i32 {
            if depth < 100 {
                let max_observed = depth / 2;
                let mut observed = max_observed;
                let mut o = current_observed;
                while o <= max_observed {
                    if binomial_cdf(depth, 0.5, o) > PROB {
                        observed = if o == 0 { -1 } else { o };
                        break;
                    }
                    o += 1;
                }
                if observed == -1 {
                    values[depth as usize] = -1.0;
                } else {
                    values[depth as usize] = (observed as f64 - 1e-6) / depth as f64;
                    current_observed = observed;
                }
            } else {
                // calcDepthVsStrandBias
                let (mut lower, mut upper) = (0.0f64, 0.5f64);
                let mut current = (lower + upper) * 0.5;
                let mut last_observed = 0;
                let mut result = None;
                for _ in 0..30 {
                    let observed = jround(depth as f64 * current) as i32;
                    if observed == 0 {
                        result = Some(-1.0);
                        break;
                    }
                    if last_observed > 0 && observed == last_observed {
                        result = Some(current);
                        break;
                    }
                    let prob = binomial_cdf(depth, 0.5, observed);
                    if (prob - PROB).abs() < 0.0001 {
                        result = Some(current);
                        break;
                    }
                    if prob < PROB {
                        lower = current;
                    } else {
                        upper = current;
                    }
                    if (upper - lower).abs() < 0.001 {
                        result = Some(current);
                        break;
                    }
                    current = (lower + upper) * 0.5;
                    last_observed = observed;
                }
                values[depth as usize] = result.unwrap_or(current);
            }
        }
        values
    })
}

fn is_depth_below_probability(alt: &StrandBias, non_alt: &StrandBias) -> bool {
    let bias = alt.bias();
    let min_bias = jmin(bias, 1.0 - bias);
    if min_bias > STRAND_BIAS_CHECK_THRESHOLD {
        return false;
    }
    if non_alt.forward < STRAND_BIAS_NON_ALT_MIN_DEPTH || non_alt.reverse < STRAND_BIAS_NON_ALT_MIN_DEPTH {
        return false;
    }
    let non_alt_bias = non_alt.bias();
    let min_non_alt = jmin(non_alt_bias, 1.0 - non_alt_bias);
    if min_non_alt < STRAND_BIAS_NON_ALT_MIN_BIAS && (non_alt_bias < 0.5) == (bias < 0.5) {
        return false;
    }
    let table = strand_bias_values();
    let depth = alt.depth() as usize;
    let required = if depth < table.len() { table[depth] } else { table[table.len() - 1] };
    if required == -1.0 {
        return false;
    }
    d_less_or_equal(min_bias, required)
}

fn all_one_side(alt: &StrandBias) -> bool {
    if alt.depth() < STRAND_BIAS_NON_ALT_MIN_DEPTH {
        return false;
    }
    alt.forward == 0 || alt.reverse == 0
}

fn read_edge_distance_threshold(tier: Tier) -> f64 {
    if matches!(tier, Tier::Hotspot | Tier::Panel) {
        0.2
    } else {
        0.33
    }
}

fn phred_to_prob(q: i8) -> f64 {
    10f64.powf(-(q as f64) / 10.0)
}

/// ReadContextCounter.phredQual
fn phred_qual(rc: &ReadContextCounter) -> f64 {
    -10.0 * jmin(jmax(rc.tumor_qual_probability, 1e-20), 1.0).log10()
}

pub struct VariantFilters<'a> {
    pub cfg: &'a Config,
    pub chr: &'a str,
    pub hla: &'a dyn Fn(i32) -> bool,
}

fn boost_novel_indel(t: &ReadContextCounter) -> bool {
    (matches!(t.tier, Tier::Hotspot | Tier::Panel) || is_illumina())
        && t.is_indel
        && !t.use_msi_error_rate()
        && t.jitter.lengthened + t.jitter.shortened == 0
        && (t.non_alt_nm_count_total as f64) / ((t.depth() - t.alt_support()) as f64) < 0.75
}

impl<'a> VariantFilters<'a> {
    /// applySoftFilters, for a single tumor sample
    pub fn apply_soft_filters(&self, v: &mut SageVariant) {
        let tier = v.tier();
        let config = tier_config(tier);
        let max_ref = v.ref_counters.len().min(self.cfg.ref_sample_count);
        let near_multi = v.near_multi_base_indel;
        for ti in 0..v.tumor_counters.len() {
            let mut filters = HashSet::new();
            self.apply_tumor_filters(near_multi, tier, config, &mut v.tumor_counters[ti], &mut filters);
            for ri in 0..if self.cfg.germline { 0 } else { max_ref } {
                self.apply_tumor_germline_filters(tier, config, &v.ref_counters[ri], &v.tumor_counters[ti], &mut filters);
            }
            if filters.is_empty() {
                v.filters.clear();
                break;
            } else {
                v.filters.extend(filters);
            }
        }
    }

    fn apply_tumor_filters(
        &self, near_multi: bool, tier: Tier, config: &SoftFilterConfig, t: &mut ReadContextCounter,
        filters: &mut HashSet<SoftFilter>,
    ) {
        let skip_qual = tier == Tier::Hotspot
            && t.alt_support() >= 8
            && d_greater_or_equal(t.vaf(), 0.08)
            && t.qual_counters.alt_recalibrated_base_quality_total >= 150;
        if !skip_qual {
            if self.below_min_tumor_qual(config, t) {
                filters.insert(SoftFilter::MinTumorQual);
            }
            if self.below_min_map_qual_factor(near_multi, config, tier, t) {
                filters.insert(SoftFilter::MinMapQualFactor);
            }
            if d_less_than(t.vaf(), config.min_tumor_vaf) && below_min_tumor_vaf_probability(t) {
                filters.insert(SoftFilter::MinTumorVaf);
            }
        }
        let skip_edge = t.variant().is_long_insert() && !near_multi;
        // belowMaxEdgeDistance
        if !skip_edge {
            let alt_med = t.read_edge_distance.max_alt_distance_from_edge;
            if !(alt_med >= read_edge_distance_threshold(tier)) {
                let max_med = t.read_edge_distance.max_distance_from_edge;
                let med_prob = (2.0 * alt_med / (2.0 * max_med)).powf(t.alt_support() as f64);
                if med_prob < 0.001 {
                    filters.insert(SoftFilter::MaxEdgeDistance);
                }
            }
        }
        // exceedsAltVsRefMapQual
        {
            let depth = t.depth() as f64;
            let alt = t.alt_support() as f64;
            if depth != 0.0 && alt != 0.0 {
                let avg = t.qual_counters.map_quality_total as f64 / depth;
                let avg_alt = t.qual_counters.alt_map_quality_total as f64 / alt;
                if avg - avg_alt > 15.0 {
                    filters.insert(SoftFilter::MapQualRefAltDifference);
                }
            }
        }
        // belowMinFragmentCoords
        {
            let required = if t.alt_support() >= 5 {
                3
            } else if t.alt_support() >= 3 {
                2
            } else {
                1
            };
            if t.fragment_coords.min_count() < required {
                filters.insert(SoftFilter::FragmentCoords);
            }
        }
        // belowMinStrongSupport
        {
            let threshold = if tier == Tier::Hotspot { 2 } else { 3 };
            let strong = if t.use_msi_error_rate() {
                t.jitter.valid_qual_full_support
            } else if is_ultima() {
                // depending on HP lengths, may not have any high qual reads
                t.strong_alt_support()
            } else {
                t.high_qual_strong_support
            };
            if strong < threshold {
                filters.insert(SoftFilter::MinTumorSupport);
            }
        }
        // exceedsRealignedPercentage
        if !t.is_long_indel() && t.alt_support() != 0 {
            let perc = t.counts.realigned as f64 / t.alt_support() as f64;
            if perc > 0.7 {
                filters.insert(SoftFilter::RealignedFreq);
            }
        }
        if tier != Tier::Hotspot {
            if is_depth_below_probability(&t.alt_fragment_strand_bias, &t.non_alt_fragment_strand_bias) {
                filters.insert(SoftFilter::FragmentStrandBias);
            }
            if is_depth_below_probability(&t.alt_read_strand_bias, &t.non_alt_read_strand_bias)
                || (t.is_indel && all_one_side(&t.alt_read_strand_bias))
            {
                filters.insert(SoftFilter::ReadStrandBias);
            }
        }
        if t.ctx.max_repeat.is_some() && t.jitter.filter_on_noise {
            filters.insert(SoftFilter::Jitter);
        }
        if is_ultima() {
            let panel_indel_repeat =
                is_panel_indel_repeat_variant(t.tier, phred_qual(t), t.ctx.variant.is_indel(), t.ctx.max_repeat_count());
            if crate::ultima::below_expected_hp_quals(t, self.cfg.germline, panel_indel_repeat) {
                filters.insert(SoftFilter::MinAvgHpQual);
            }
            if crate::ultima::below_expected_t0_quals(t, all_one_side(&t.alt_read_strand_bias)) {
                filters.insert(SoftFilter::MinAvgT0Qual);
            }
        } else {
            // belowMinAverageBaseQuality
            let threshold = if tier == Tier::Hotspot { self.cfg.min_avg_base_qual_hotspot } else { self.cfg.min_avg_base_qual };
            if !t.use_msi_error_rate() && d_less_than(t.average_alt_seq_tech_base_quality(), threshold as f64) {
                filters.insert(SoftFilter::MinAvgBaseQuality);
            }
        }
    }

    fn below_min_tumor_qual(&self, config: &SoftFilterConfig, t: &mut ReadContextCounter) -> bool {
        let depth = t.depth();
        let strong = t.strong_alt_support();
        if strong == 0 {
            t.tumor_qual_probability = 1.0;
            return true;
        }
        let q = &t.qual_counters;
        let strong_non_medium = strong - t.medium_qual_strong_support;
        let medium_total = q.strong_alt_medium_recalibrated_base_quality_total;
        let non_medium_total = q.strong_alt_recalibrated_base_quality_total - medium_total;
        let medium_contribution = (strong_non_medium as f64 * medium_total / non_medium_total) as i32;
        let adjusted_strong = strong_non_medium + medium_contribution;
        let alt_non_medium_final = q.alt_final_base_quality_total - q.alt_medium_final_base_quality_total;
        let mut qual_per_read = jround(alt_non_medium_final / strong_non_medium as f64) as i32;
        if boost_novel_indel(t) {
            qual_per_read += DEFAULT_BASE_QUAL_FIXED_PENALTY;
        }
        let min_qual = if t.qual_cache.is_msi_sample_and_variant { 20 } else { 15 };
        let floored = qual_per_read.max(min_qual) as i8;
        let mut prob = 1.0 - binomial_cdf(depth, phred_to_prob(floored), adjusted_strong - 1);
        if self.cfg.germline {
            // GERMLINE_HET_MIN_EXPECTED_VAF, GERMLINE_HET_MIN_SAMPLING_PROB
            let het_prob = binomial_cdf(depth, 0.4, t.alt_support());
            if het_prob < prob {
                prob /= jmax(het_prob, 1e-3);
                prob = jmin(prob, 1.0);
            }
        }
        t.tumor_qual_probability = prob;
        prob >= config.qual_p_score
    }

    fn below_min_map_qual_factor(&self, near_multi: bool, config: &SoftFilterConfig, tier: Tier, t: &mut ReadContextCounter) -> bool {
        let depth = t.depth();
        let alt_support = t.alt_support();
        if depth == 0 || alt_support == 0 {
            return false;
        }
        let q = &t.qual_counters;
        let mut avg_alt_final_map_quality = q.alt_final_map_quality_total / t.strong_alt_support() as f64;
        let avg_map_qual = q.map_quality_total as f64 / depth as f64;
        let avg_alt_map_qual = q.alt_map_quality_total as f64 / alt_support as f64;
        let mut map_qual_diff_penalty = jmax(2.0 * (avg_map_qual - avg_alt_map_qual), 0.0);
        let polymorphic = (self.hla)(t.position());
        if polymorphic {
            avg_alt_final_map_quality = jmin((25 + HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY) as f64, avg_alt_map_qual);
            if avg_alt_map_qual >= 40.0 {
                map_qual_diff_penalty = 0.0;
            }
        }
        let alt = &t.alt_read_strand_bias;
        let non_alt = &t.non_alt_read_strand_bias;
        let read_strand_bias_penalty = if alt.min_bias() > STRAND_BIAS_CHECK_THRESHOLD
            || non_alt.forward < STRAND_BIAS_NON_ALT_MIN_DEPTH
            || non_alt.reverse < STRAND_BIAS_NON_ALT_MIN_DEPTH
            || (non_alt.min_bias() < STRAND_BIAS_NON_ALT_MIN_BIAS && (non_alt.bias() < 0.5) == (alt.bias() < 0.5))
        {
            0.0
        } else {
            let p = 2.0 * binomial_cdf(alt.depth(), 0.5, jround(alt.depth() as f64 * alt.min_bias()) as i32);
            if p > 0.0 {
                -10.0 * p.log10()
            } else {
                50.0
            }
        };
        let avg_edge = t.read_edge_distance.avg_distance_from_edge();
        let avg_alt_edge = t.read_edge_distance.avg_alt_distance_from_edge();
        let mut edge_penalty = 0.0;
        let ratio = avg_alt_edge / avg_edge;
        let skip_edge = t.variant().is_long_insert() && !near_multi;
        if !polymorphic && ratio < 2.0 * read_edge_distance_threshold(tier) && !skip_edge {
            edge_penalty = 10.0 * alt_support as f64 * (avg_edge / jmax(avg_alt_edge, 0.001)).log10();
            if !is_sbx() && !near_multi && !self.cfg.high_depth_mode {
                edge_penalty = jmin(edge_penalty, 10.0);
            }
        }
        let mut repeat_penalty = 0.0;
        if let Some(r) = &t.ctx.max_repeat {
            if r.repeat_length() > 1 && r.repeat_length() * r.count >= 15 {
                let max_penalty = if t.is_indel { 18 } else { 24 };
                repeat_penalty = (3 * r.count).min(max_penalty) as f64;
            }
        }
        let mqf = avg_alt_final_map_quality - 25.0 - map_qual_diff_penalty - read_strand_bias_penalty - edge_penalty - repeat_penalty;
        t.map_qual_factor = mqf;
        let mut threshold = config.map_qual_factor;
        if is_sbx() && t.min_number_of_events == 1 {
            threshold -= 10.0; // SbxUtils.MQF_NM_1_THRESHOLD_DEDUCTION
        }
        mqf < threshold
    }

    fn apply_tumor_germline_filters(
        &self, tier: Tier, config: &SoftFilterConfig, r: &ReadContextCounter, t: &ReadContextCounter,
        filters: &mut HashSet<SoftFilter>,
    ) {
        // belowMinGermlineCoverage
        let allosome = matches!(chromosome_rank(self.chr), Some(23 | 24));
        let min_cov = if allosome { config.min_germline_coverage_allosome } else { config.min_germline_coverage };
        if r.depth() < min_cov {
            filters.insert(SoftFilter::MinGermlineDepth);
        }
        let valid_indel_repeat =
            is_panel_or_long_indel_repeat_variant(t.tier, phred_qual(t), t.is_indel, t.ctx.max_repeat_count());
        // aboveMaxGermlineVaf
        let tumor_vaf = t.vaf();
        if tumor_vaf != 0.0 {
            let mut adjusted_ref_alt = r.counts.alt_support() + r.simple_alt_matches;
            let max_repeat_bases: &[u8] = t.ctx.max_repeat.as_ref().map_or(&[], |r| &r.bases);
            if t.is_long_indel()
                || include_ref_jitter_in_msi_indel(
                    t.ctx.max_repeat_count(), max_repeat_bases.len() as i32, t.variant().indel_length_abs(),
                )
            {
                adjusted_ref_alt += r.jitter.shortened + r.jitter.lengthened;
            }
            if above_max_germline_vaf(tier, valid_indel_repeat, tumor_vaf, adjusted_ref_alt, r.counts.total, config.max_germline_vaf) {
                filters.insert(SoftFilter::MaxGermlineVaf);
            }
        }
        // aboveMaxGermlineRelativeQual
        {
            let tumor_avg = jround(t.average_alt_recalibrated_base_quality()) as i32 as f64;
            let ref_avg = jround(r.average_alt_recalibrated_base_quality()) as i32 as f64;
            let tumor_qual = (t.qualities.full + t.qualities.partial_core + t.qualities.realigned) as f64;
            let ref_qual = (r.qualities.full + r.qualities.partial_core + r.qualities.realigned) as f64;
            if above_max_germline_relative_qual(
                tier, tumor_qual, t.vaf(), t.depth(), tumor_avg, ref_qual, r.depth(), r.alt_support(), ref_avg, valid_indel_repeat,
            ) {
                filters.insert(SoftFilter::MaxGermlineRelativeQual);
            }
        }
        // aboveMaxMnvIndelGermlineAltSupport
        if tier != Tier::Hotspot && (r.variant().is_mnv() || r.variant().is_long_insert()) {
            let perc = if r.depth() > 0 { r.alt_support() as f64 / r.depth() as f64 } else { 0.0 };
            if perc >= 0.01 {
                filters.insert(SoftFilter::MaxGermlineAltSupport);
            }
        }
    }
}

fn below_min_tumor_vaf_probability(t: &ReadContextCounter) -> bool {
    let depth = t.depth();
    let base_qual_avg = t.qual_counters.recalibrated_base_quality_total as f64 / depth as f64;
    let support = t.strong_alt_support();
    let mut alt_base_qual_avg = t.average_alt_recalibrated_base_quality();
    if boost_novel_indel(t) {
        alt_base_qual_avg += DEFAULT_BASE_QUAL_FIXED_PENALTY as f64;
    }
    let adjusted = support * jround(jmin(alt_base_qual_avg / base_qual_avg, 1.0)) as i32;
    let probability = 10f64.powf(-base_qual_avg / 10.0);
    let p_value = 1.0 - binomial_cdf(depth, probability, adjusted - 1);
    if t.tier == Tier::Hotspot {
        p_value > 1e-9
    } else {
        p_value > 1e-14
    }
}

/// VariantFilters.includeRefJitterInMsiIndel: Ultima only. Takes the length of the read context's max
/// repeat unit (0 if there is no repeat) and the absolute indel length.
pub fn include_ref_jitter_in_msi_indel(rc_repeat_count: i32, rc_repeat_bases_length: i32, indel_length_abs: i32) -> bool {
    // a repeat count without repeat bases (possible for TINC, which reads both from the VCF) divides by
    // zero in Java
    if !is_ultima() || indel_length_abs == 0 || rc_repeat_count == 0 || rc_repeat_bases_length == 0 {
        return false;
    }
    // DOUBLE_JITTER_MIN_REPEAT_CHANGE, JITTER_MIN_REPEAT_CHANGE
    let min_repeat_change = if rc_repeat_count >= DEFAULT_DOUBLE_JITTER_REPEAT_COUNT { 4 } else { 3 };
    indel_length_abs / rc_repeat_bases_length >= min_repeat_change
}

/// VariantFilters.isPanelIndelRepeatVariant
pub fn is_panel_indel_repeat_variant(tier: Tier, qual: f64, is_indel: bool, repeat_length: i32) -> bool {
    // PANEL_REPEAT_MIN_QUAL
    matches!(tier, Tier::Panel | Tier::Hotspot) && qual >= 100.0 && is_indel && repeat_length > 0
}

/// VariantFilters.isPanelOrLongIndelRepeatVariant: only affects the SBX and Ultima germline thresholds
pub fn is_panel_or_long_indel_repeat_variant(tier: Tier, qual: f64, is_indel: bool, repeat_length: i32) -> bool {
    // PANEL_REPEAT_MIN_QUAL, LONG_INDEL_REPEAT_MIN, LONG_INDEL_REPEAT_MAX
    let is_indel_repeat = qual >= 100.0 && is_indel && repeat_length > 0;
    (matches!(tier, Tier::Panel | Tier::Hotspot) && is_indel_repeat)
        || ((7..=18).contains(&repeat_length) && is_indel_repeat)
}

/// VariantFilters.aboveMaxGermlineVaf (static form, also used by TINC)

pub fn above_max_germline_vaf(
    tier: Tier, is_valid_indel_repeat: bool, tumor_vaf: f64, adjusted_ref_alt: i32, ref_total: i32,
    mut max_germline_vaf: f64,
) -> bool {
    if tumor_vaf == 0.0 {
        return false;
    }
    let is_valid_ultima_indel_repeat = is_ultima() && is_valid_indel_repeat;
    if matches!(tier, Tier::Panel | Tier::Hotspot) {
        let mut threshold = tumor_vaf;
        if tier == Tier::Panel {
            // MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_VAF_FACTOR, MAX_GERMLINE_VAF_PANEL_VAF_FACTOR
            threshold /= if is_valid_ultima_indel_repeat { 1.5 } else { 2.0 };
        }
        // MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_THRESHOLD_FACTOR, MAX_GERMLINE_VAF_THRESHOLD_FACTOR
        let min_threshold_factor = if is_valid_ultima_indel_repeat { 5.0 } else { 2.0 };
        max_germline_vaf = jmax(jmin(threshold, max_germline_vaf * min_threshold_factor), max_germline_vaf);
    } else if is_valid_ultima_indel_repeat {
        let threshold = tumor_vaf / crate::ultima::GERMLINE_VAF_INDEL_REPEAT_THRESHOLD;
        let min_threshold_factor = crate::ultima::GERMLINE_VAF_INDEL_REPEAT_MIN_THRESHOLD_FACTOR;
        max_germline_vaf = jmax(jmin(threshold, max_germline_vaf * min_threshold_factor), max_germline_vaf);
    } else if is_sbx() && is_valid_indel_repeat {
        // SbxUtils.GERMLINE_VAF_INDEL_REPEAT_THRESHOLD, GERMLINE_VAF_INDEL_REPEAT_MIN_THRESHOLD_FACTOR
        let threshold = tumor_vaf / 3.5;
        max_germline_vaf = jmax(jmin(threshold, max_germline_vaf * 2.0), max_germline_vaf);
    }
    let adjusted_ref_vaf = adjusted_ref_alt as f64 / ref_total as f64;
    d_greater_than(adjusted_ref_vaf, max_germline_vaf)
}

/// VariantFilters.aboveMaxMnvIndelGermlineAltSupport (static form)
pub fn above_max_mnv_indel_germline_alt_support(
    tier: Tier, is_mnv: bool, is_long_insert: bool, ref_depth: i32, ref_alt_support: i32,
) -> bool {
    if tier == Tier::Hotspot {
        return false;
    }
    if is_mnv || is_long_insert {
        let perc = if ref_depth > 0 { ref_alt_support as f64 / ref_depth as f64 } else { 0.0 };
        return perc >= 0.01;
    }
    false
}

#[allow(clippy::too_many_arguments)]
pub fn above_max_germline_relative_qual(
    tier: Tier, tumor_qual: f64, tumor_vaf: f64, tumor_depth: i32, tumor_avg_base_qual: f64, ref_qual: f64,
    ref_depth: i32, ref_alt_support: i32, ref_avg_base_qual: f64, is_valid_indel_repeat: bool,
) -> bool {
    if tumor_qual == 0.0 {
        return false;
    }
    if ref_depth == 0 || tumor_depth == 0 {
        return false;
    }
    let adjusted_tumor_vaf = jmin(jmax(tumor_vaf - 10f64.powf(-tumor_avg_base_qual / 10.0), 0.0), 0.5);
    let adjusted_ref_ad =
        (jround(ref_alt_support as f64 - ref_depth as f64 * 10f64.powf(-ref_avg_base_qual / 10.0)) as i32).max(0);
    let depth_ratio = ref_depth as f64 / tumor_depth as f64;
    let adjusted_ref_qual_ratio = ref_qual / tumor_qual / depth_ratio;
    let prob = binomial_cdf(ref_depth, adjusted_tumor_vaf, adjusted_ref_ad);
    let prob_threshold = match tier {
        Tier::Hotspot => 0.01,
        Tier::Panel => 0.00025,
        _ => 1e-7,
    };
    let ratio_threshold = if tier == Tier::Hotspot {
        0.5
    } else if is_sbx() && is_valid_indel_repeat {
        0.2 // SbxUtils.GERMLINE_VAF_REL_QUAL_RATIO_THRESHOLD
    } else if is_ultima() && is_valid_indel_repeat {
        crate::ultima::GERMLINE_VAF_REL_QUAL_RATIO_THRESHOLD
    } else {
        0.12
    };
    prob > prob_threshold && adjusted_ref_qual_ratio > ratio_threshold
}

const NEAR_INDEL_MIN_VAF: f64 = 0.15;
const NEAR_INDEL_PROXIMITY: i32 = 30;

/// FilterUtils.setNearByIndelStatusPreFilter
pub fn set_near_by_indel_status_pre_filter(variants: &mut [SageVariant]) {
    for index in 0..variants.len() {
        let v = &variants[index];
        if !v.is_indel() || v.variant().indel_length_abs() < 2 {
            continue;
        }
        if v.tumor_counters[0].vaf() < NEAR_INDEL_MIN_VAF {
            continue;
        }
        let indel_start = v.position() - NEAR_INDEL_PROXIMITY;
        let indel_end = v.variant().position_end() + NEAR_INDEL_PROXIMITY;
        for search_up in [true, false] {
            let mut other = if search_up { index as i64 + 1 } else { index as i64 - 1 };
            while other >= 0 && (other as usize) < variants.len() {
                let o = &mut variants[other as usize];
                let other_end = o.variant().position_end();
                if positions_overlap(indel_start, indel_end, o.position(), other_end) {
                    o.near_multi_base_indel = true;
                } else if (search_up && o.position() > indel_end) || (!search_up && other_end < indel_start) {
                    break;
                }
                other += if search_up { 1 } else { -1 };
            }
        }
    }
}

/// FilterUtils.setNearByIndelStatusPostFilter
pub fn set_near_by_indel_status_post_filter(variants: &mut [SageVariant]) {
    for index in 0..variants.len() {
        let v = &variants[index];
        if !v.is_indel() {
            continue;
        }
        if !v.is_passing() && v.filters.iter().any(|f| f.is_tumor()) {
            continue;
        }
        let (start, end) = {
            let ctx = &v.tumor_counters[0].ctx;
            (ctx.alignment_start, ctx.alignment_end)
        };
        for search_up in [true, false] {
            let mut other = if search_up { index as i64 + 1 } else { index as i64 - 1 };
            while other >= 0 && (other as usize) < variants.len() {
                let o = &mut variants[other as usize];
                if position_within(o.position(), start, end) {
                    o.near_indel = true;
                } else if (search_up && o.position() > end) || (!search_up && o.position() < start) {
                    break;
                }
                other += if search_up { 1 } else { -1 };
            }
        }
    }
}

/// VariantFilters.checkFinalFilters (somatic, not panel-only)
pub fn check_final_filters(v: &SageVariant, passing_phase_sets: &[i32], cfg: &Config, chr: &str) -> bool {
    if cfg.panel_only && !matches!(v.tier(), Tier::Hotspot | Tier::Panel) {
        return false;
    }
    if v.is_passing() {
        return true;
    }
    if cfg.disable_hard_filter {
        return true;
    }
    if matches!(v.tier(), Tier::Hotspot | Tier::Panel) {
        return true;
    }
    if v.mixed_impact > 0 {
        return true;
    }
    if !cfg.germline && !v.ref_counters.is_empty() && !is_mito(chr) && !v.has_matching_lps(passing_phase_sets) {
        if v.ref_counters[0].alt_support() > cfg.filtered_max_germline_alt_support {
            return false;
        }
    }
    true
}
