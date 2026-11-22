//! Variant filters: FilterConstants, VariantFilters, SeqTechUtils, LineChecker and Deduplication.

use super::types::*;
use super::variant::*;
use super::vcf::Genotype;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

impl SequencingType {
    pub fn value_of(s: &str) -> Option<SequencingType> {
        match s {
            "ILLUMINA" => Some(SequencingType::Illumina),
            "ULTIMA" => Some(SequencingType::Ultima),
            "SBX" => Some(SequencingType::Sbx),
            _ => None,
        }
    }
}

#[cfg(test)]
pub const DEFAULT_MIN_QUAL: i32 = 30;
pub const DEFAULT_MIN_QUAL_PANEL: i32 = 60;
#[cfg(test)]
pub const DEFAULT_MIN_SUPPORT_JUNCTION: i32 = 4;
#[cfg(test)]
pub const DEFAULT_MIN_SUPPORT_SGL: i32 = 6;
#[cfg(test)]
pub const DEFAULT_MIN_SUPPORT_HOTSPOT: i32 = 2;
pub const DEFAULT_MIN_AF_SGL: f64 = 0.05;
pub const DEFAULT_MIN_AF_HOTSPOT: f64 = 0.001;
pub const DEFAULT_MIN_AF_JUNCTION: f64 = 0.002;
#[cfg(test)]
pub const DEFAULT_PON_DISTANCE: i32 = 3;
#[cfg(test)]
pub const DEFAULT_SGL_PON_DISTANCE: i32 = 10;

pub const PON_SHORT_INDEL_LENGTH: i32 = 200;
pub const PON_MAX_INS_SEQ_LENGTH: i32 = 100;
pub const PON_SHORT_INDEL_PON_DISTANCE: i32 = 10;
pub const PON_SHORT_INDEL_MAX_REPEAT_PON_DISTANCE: i32 = 20;
pub const PON_SHORT_INDEL_MAX_REPEATS: i32 = 11;

const MIN_VARIANT_LENGTH: i32 = 32;
const MIN_ANCHOR_LENGTH: i32 = 50;
const MIN_TRIMMED_ANCHOR_LENGTH: i32 = MIN_ANCHOR_LENGTH;
const MIN_AVG_FRAG_FACTOR: f64 = 3.0;
const MIN_AVG_FRAG_STD_DEV_FACTOR: f64 = 0.6;

const DEL_ARTEFACT_SHORT_LENGTH: i32 = 3000;
const DEL_ARTEFACT_MAX_HOMOLOGY: i32 = 10;
const DEL_ARTEFACT_LENGTH_FACTOR: f64 = 1.5;
const DEL_ARTEFACT_MIN_AF: f64 = 0.05;

const INV_SHORT_LENGTH: i32 = 3000;
const INV_SHORT_MIN_AF_LOWER: f64 = 0.1;
const INV_SHORT_MIN_AF_HIGHER: f64 = 0.2;
const INV_SHORT_MAX_HOMOLOGY_LOWER: i32 = 3;
const INV_SHORT_MAX_HOMOLOGY_HIGHER: i32 = 6;
const INV_SHORT_RATE_LOWER: f64 = 200.0;
const INV_SHORT_RATE_HIGHER: f64 = 400.0;

const INV_SHORT_FRAGMENT_LENGTH: i32 = 300;
const INV_SHORT_FRAGMENT_MIN_AF: f64 = 0.05;
const INV_SHORT_FRAGMENT_AF_RATIO: f64 = 50.0;

const PRIME_MAX_PERMITTED_RANGE: i32 = 50;
const PRIME_MAX_BASE_FACTOR: i32 = 10;
const PRIME_MAX_SGL_FACTOR: i32 = 5;

const INV_ADJACENT_LENGTH: i32 = 100;
const INV_ADJACENT_MIN_UPS: i32 = 4;
const INV_ADJACENT_EXCLUDED_FILTERS: [FilterType; 4] =
    [FilterType::Pon, FilterType::InvShortFragLowVaf, FilterType::InvShortLowVafHom, FilterType::InvShortIsolated];

const MIN_AF_WEAK_JUNCTION_DISC_RATE_MIN: f64 = 0.005;
const MIN_AF_WEAK_JUNCTION_STRAND_BIAS_MIN: f64 = 0.1;
const MIN_AF_WEAK_JUNCTION_SUPP_REMOTE_THRESHOLD: f64 = 0.5;
const MIN_AF_WEAK_JUNCTION_PROX_JUNC_READ_THRESHOLD: f64 = 0.1;
const MIN_AF_WEAK_JUNCTION_DISC_RATE_FACTOR_LONG: f64 = 1.25;
const MIN_AF_WEAK_JUNCTION_DISC_RATE_FACTOR_SGL: f64 = 1.0;
const MIN_AF_WEAK_JUNCTION_DISC_RATE_FACTOR_OTHER: f64 = 0.75;
const WEAK_ASSEMBLY_LONG_LENGTH: i32 = 100_000;

const PON_INS_SEQ_FWD_STRAND_1: &str = "GCCGTATCATTAAAAA";
const PON_INS_SEQ_FWD_STRAND_2: &str = "GTAGATCTCGGTGGTC";

pub const PANEL_INCLUSION_BUFFER: i32 = 1000;

// LINE elements
const LINE_POLY_AT_TEST_LEN: i32 = 18;
const LINE_POLY_AT_REQ: i32 = 16;
const LINE_INDEL_MAX_OVERLAP: i32 = 40;
const LINE_INDEL_MAX_GAP: i32 = 30;

const LOCAL_JUNCTION_MAX_DISTANCE: i32 = 50;

#[derive(Clone, Debug)]
pub struct FilterConstants {
    pub min_qual: i32,
    pub min_qual_hotspot: i32,
    pub min_length: i32,
    pub min_support_junction: i32,
    pub min_support_sgl: i32,
    pub min_support_hotspot: i32,
    pub min_af_junction: f64,
    pub min_af_sgl: f64,
    pub min_af_hotspot: f64,
    pub filter_sgls: bool,
    pub low_qual_region: ChrBaseRegion,
}

/// the PMS2 region, where the quality threshold is halved
fn low_qual_region(version: RefGenomeVersion) -> ChrBaseRegion {
    if version.is37() { ChrBaseRegion::new("7", 6002870, 6058756) } else { ChrBaseRegion::new("chr7", 5960925, 6019106) }
}

impl FilterConstants {
    pub fn new(
        min_qual: i32, min_qual_hotspot: i32, min_support: i32, min_support_sgl: i32, min_support_hotspot: i32, filter_sgls: bool,
        version: RefGenomeVersion,
    ) -> FilterConstants {
        FilterConstants {
            min_qual,
            min_qual_hotspot,
            min_length: MIN_VARIANT_LENGTH,
            min_support_junction: min_support,
            min_support_sgl,
            min_support_hotspot,
            min_af_junction: DEFAULT_MIN_AF_JUNCTION,
            min_af_sgl: DEFAULT_MIN_AF_SGL,
            min_af_hotspot: DEFAULT_MIN_AF_HOTSPOT,
            filter_sgls,
            low_qual_region: low_qual_region(version),
        }
    }

    /// FilterConstants.from(filterSgls, refGenomeVersion, ponDistance), as used by the tests
    #[cfg(test)]
    pub fn defaults(filter_sgls: bool, version: RefGenomeVersion) -> FilterConstants {
        Self::new(
            DEFAULT_MIN_QUAL,
            DEFAULT_MIN_QUAL,
            DEFAULT_MIN_SUPPORT_JUNCTION,
            DEFAULT_MIN_SUPPORT_SGL,
            DEFAULT_MIN_SUPPORT_HOTSPOT,
            filter_sgls,
            version,
        )
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FragmentLengthBounds {
    pub lower_bound: i32,
    pub upper_bound: i32,
    pub median: i32,
    pub std_deviation: f64,
}

pub const DISCORDANT_FRAG_TYPES: usize = 11;
pub const DISC_INV_LT_1K: usize = 1;
pub const DISC_INV_1_TO_5K: usize = 2;

#[derive(Clone, Debug, Default)]
pub struct DiscordantStats {
    pub total_reads: i64,
    pub type_counts: [i64; DISCORDANT_FRAG_TYPES],
}

impl DiscordantStats {
    pub fn short_inversion_rate(&self) -> f64 {
        let count = self.type_counts[DISC_INV_LT_1K] + self.type_counts[DISC_INV_1_TO_5K];
        if self.total_reads > 0 { count as f64 / self.total_reads as f64 } else { 0.0 }
    }
    pub fn short_inversion_fragment_rate(&self) -> f64 {
        let count = self.type_counts[DISC_INV_LT_1K];
        if self.total_reads > 0 { count as f64 / self.total_reads as f64 } else { 0.0 }
    }
    pub fn discordant_rate(&self) -> f64 {
        // Java sums the counts as doubles
        let total: f64 = self.type_counts.iter().map(|&c| c as f64).sum();
        if self.total_reads > 0 { total / self.total_reads as f64 } else { 0.0 }
    }
}

pub struct VariantFilters {
    pub constants: FilterConstants,
    pub fragment_length_bounds: FragmentLengthBounds,
    pub discordant_stats: DiscordantStats,
    pub sequencing_type: SequencingType,
    short_inversion_rate: f64,
    short_fragment_inversion_rate: f64,
}

fn adjusted_strand_bias(genotype: &Genotype) -> f64 {
    let strand_bias = genotype.attr_double(STRAND_BIAS, 0.5);
    if strand_bias > 0.5 { 1.0 - strand_bias } else { strand_bias }
}

impl VariantFilters {
    pub fn new(
        constants: FilterConstants, fragment_length_bounds: FragmentLengthBounds, discordant_stats: DiscordantStats,
        sequencing_type: SequencingType,
    ) -> VariantFilters {
        VariantFilters {
            short_inversion_rate: discordant_stats.short_inversion_rate(),
            short_fragment_inversion_rate: discordant_stats.short_inversion_fragment_rate(),
            constants,
            fragment_length_bounds,
            discordant_stats,
            sequencing_type,
        }
    }

    fn has_paired_reads(&self) -> bool {
        self.sequencing_type == SequencingType::Illumina
    }
    fn is_sbx(&self) -> bool {
        self.sequencing_type == SequencingType::Sbx
    }

    /// VariantFilters.applyFilters for the variant at `v`
    pub fn apply_filters(&self, variants: &mut [Variant], v: usize) {
        let mut add = FilterSet::default();
        {
            let var = &variants[v];
            // keep existing filters eg from assembly
            for breakend in var.breakends.iter().flatten() {
                for filter in &breakend.context.filters {
                    if let Some(filter_type) = FilterType::from_vcf_tag(filter) {
                        if filter_type != FilterType::Duplicate {
                            add.add(filter_type);
                        }
                    }
                }
            }

            if self.constants.filter_sgls && var.is_sgl() && !var.is_hotspot {
                add.add(FilterType::Sgl);
            }
            if self.below_min_af(var) {
                add.add(FilterType::MinAf);
            }
            if self.below_weak_junction_min_af(var) {
                add.add(FilterType::WeakJuncMinAf);
            }
            if self.below_min_support(variants, var) {
                add.add(FilterType::MinSupport);
            }
            if self.below_min_quality(variants, var) {
                add.add(FilterType::MinQuality);
            }
            if self.below_min_length(var) {
                add.add(FilterType::MinLength);
            }
            if self.below_min_anchor_length(var) {
                add.add(FilterType::MinAnchorLength);
            }
            if self.below_min_fragment_length(var) {
                add.add(FilterType::ShortFragLength);
            }
            if self.is_inversion_short_low_vaf_homology(var) {
                add.add(FilterType::InvShortLowVafHom);
            }
            if self.is_inversion_short_fragment_low_vaf(var) {
                add.add(FilterType::InvShortFragLowVaf);
            }
            if self.is_deletion_short_low_vaf(var) {
                add.add(FilterType::DelShortLowVaf);
            }
            if fails_strand_bias(var) {
                add.add(FilterType::StrandBias);
            }
            if self.low_three_prime_position_range(var) {
                add.add(FilterType::UnpairedThreePrimeRange);
            }
            if self.is_sbx_strand_bias(var) {
                add.add(FilterType::SbxStrandBias);
            }
            if self.is_sbx_artefact(var) {
                add.add(FilterType::SbxArtefact);
            }
            if self.is_sbx_zero_length_inversion(var) {
                add.add(FilterType::SbxInvZeroMotif);
            }
            if var.breakends.iter().flatten().any(|b| b.line_remote_source_site) {
                add.add(FilterType::LineSource);
            }
        }
        for filter in add.iter() {
            variants[v].add_filter(filter);
        }
    }

    fn below_min_support(&self, variants: &[Variant], var: &Variant) -> bool {
        if var.is_saga_matched {
            return false;
        }
        let support_threshold = if var.is_hotspot {
            self.constants.min_support_hotspot
        } else if var.is_sgl() && !var.is_line_site() {
            self.constants.min_support_sgl
        } else {
            self.constants.min_support_junction
        } as f64;

        let breakend = var.start();
        let java_c2 = crate::common::compat::java("C2");
        for genotype in &breakend.context.genotypes {
            let mut fragment_count = Breakend::fragment_count_of(genotype);
            if let (true, Some(line)) = (var.is_sgl(), breakend.line_site_breakend) {
                if java_c2 {
                    // Java (C2): asks the paired LINE breakend for the count of this same genotype object,
                    // so the breakend's own count is added a second time
                    fragment_count += Breakend::fragment_count_of(genotype);
                } else {
                    // the support of the paired breakend in the same sample
                    let pair = variants[line.var as usize].breakend(line.se as usize);
                    fragment_count += pair.context.genotype(&genotype.sample).map(Breakend::fragment_count_of).unwrap_or(0);
                }
            }
            if fragment_count as f64 >= support_threshold {
                return false;
            }
        }
        true
    }

    fn below_min_af(&self, var: &Variant) -> bool {
        if var.is_saga_matched {
            return false;
        }
        let af_threshold = if var.is_hotspot {
            self.constants.min_af_hotspot
        } else if var.is_sgl() && !var.is_line_site() {
            self.constants.min_af_sgl
        } else {
            self.constants.min_af_junction
        };
        !any_samples_above_af_threshold(var, af_threshold)
    }

    fn below_weak_junction_min_af(&self, var: &Variant) -> bool {
        let discordant_rate = self.discordant_stats.discordant_rate();
        if discordant_rate < MIN_AF_WEAK_JUNCTION_DISC_RATE_MIN {
            return false;
        }
        if !apply_weak_junction_min_af_filter(var) {
            return false;
        }
        let mut af_threshold = if var.is_sgl() && !var.is_line_site() { self.constants.min_af_sgl } else { self.constants.min_af_junction };
        let disc_rate_threshold = MIN_AF_WEAK_JUNCTION_DISC_RATE_MIN;
        if discordant_rate > disc_rate_threshold {
            if var.sv_type == SvType::Sgl {
                af_threshold += (discordant_rate - disc_rate_threshold) * MIN_AF_WEAK_JUNCTION_DISC_RATE_FACTOR_SGL;
            } else if (var.sv_type == SvType::Bnd || var.sv_length() >= WEAK_ASSEMBLY_LONG_LENGTH) && var.original_junctions.len() == 1 {
                af_threshold += (discordant_rate - disc_rate_threshold) * MIN_AF_WEAK_JUNCTION_DISC_RATE_FACTOR_LONG;
            } else {
                af_threshold += (discordant_rate - disc_rate_threshold) * MIN_AF_WEAK_JUNCTION_DISC_RATE_FACTOR_OTHER;
            }
        }
        !any_samples_above_af_threshold(var, af_threshold)
    }

    fn below_min_quality(&self, variants: &[Variant], var: &Variant) -> bool {
        if var.is_saga_matched {
            return false;
        }
        let mut qual_threshold = if var.is_hotspot { self.constants.min_qual_hotspot } else { self.constants.min_qual } as f64;
        let breakend = var.start();
        let region = &self.constants.low_qual_region;
        if region.contains_position(&breakend.chromosome, breakend.position) {
            qual_threshold *= 0.5;
        } else if let Some(other) = var.end() {
            if region.contains_position(&other.chromosome, other.position) {
                qual_threshold *= 0.5;
            }
        }
        let mut variant_qual = var.qual;
        if let Some(line) = breakend.line_site_breakend {
            variant_qual += variants[line.var as usize].qual;
        }
        variant_qual < qual_threshold
    }

    fn below_min_length(&self, var: &Variant) -> bool {
        if matches!(var.sv_type, SvType::Del | SvType::Dup | SvType::Ins) {
            return var.adjusted_length() < self.constants.min_length;
        }
        false
    }

    fn below_min_anchor_length(&self, var: &Variant) -> bool {
        if var.is_saga_matched || var.is_line_site() || var.in_chained_assembly() {
            return false;
        }
        if var.start().anchor_length() < MIN_TRIMMED_ANCHOR_LENGTH {
            return true;
        }
        if var.is_sgl() {
            let insert_sequence = var.insert_sequence.as_bytes();
            let repeats = find_repeats(insert_sequence);
            let trimmed = calc_trimmed_base_length(0, insert_sequence.len() as i32 - 1, &repeats);
            trimmed < MIN_TRIMMED_ANCHOR_LENGTH
        } else {
            var.breakend(SE_END).anchor_length() < MIN_TRIMMED_ANCHOR_LENGTH
        }
    }

    fn below_min_fragment_length(&self, var: &Variant) -> bool {
        if var.is_saga_matched || !self.has_paired_reads() || var.is_sgl() || var.is_line_site() {
            return false;
        }
        let sv_avg_length = var.average_fragment_length();
        let split_frag_count = var.split_fragment_count();
        if sv_avg_length == 0 {
            return false;
        }
        let median_length = self.fragment_length_bounds.median;
        let std_deviation = self.fragment_length_bounds.std_deviation;
        let lower_std_dev_threshold = MIN_AVG_FRAG_STD_DEV_FACTOR * std_deviation;
        let lower_length_limit =
            median_length as f64 - f64::max(MIN_AVG_FRAG_FACTOR * std_deviation / (split_frag_count as f64).sqrt(), lower_std_dev_threshold);
        (sv_avg_length as f64) < lower_length_limit
    }

    fn is_inversion_short_low_vaf_homology(&self, var: &Variant) -> bool {
        if var.is_saga_matched || var.sv_type != SvType::Inv || var.adjusted_length() > INV_SHORT_LENGTH {
            return false;
        }
        let inexact_homology = var.start().inexact_homology.length();
        if inexact_homology < INV_SHORT_MAX_HOMOLOGY_LOWER {
            return false;
        }
        let vaf_threshold = if inexact_homology >= INV_SHORT_MAX_HOMOLOGY_HIGHER {
            f64::min(INV_SHORT_MIN_AF_HIGHER, self.short_inversion_rate * INV_SHORT_RATE_HIGHER)
        } else {
            f64::min(INV_SHORT_MIN_AF_LOWER, self.short_inversion_rate * INV_SHORT_RATE_LOWER)
        };
        !any_samples_above_af_threshold(var, vaf_threshold)
    }

    fn is_inversion_short_fragment_low_vaf(&self, var: &Variant) -> bool {
        if var.is_saga_matched || var.sv_type != SvType::Inv || var.adjusted_length() > INV_SHORT_FRAGMENT_LENGTH {
            return false;
        }
        let vaf_threshold = f64::min(INV_SHORT_FRAGMENT_MIN_AF, INV_SHORT_FRAGMENT_AF_RATIO * self.short_fragment_inversion_rate);
        !any_samples_above_af_threshold(var, vaf_threshold)
    }

    fn is_deletion_short_low_vaf(&self, var: &Variant) -> bool {
        if var.is_saga_matched || var.sv_type != SvType::Del || var.adjusted_length() > DEL_ARTEFACT_SHORT_LENGTH {
            return false;
        }
        let inexact_homology = var.start().inexact_homology.length();
        if inexact_homology < DEL_ARTEFACT_MAX_HOMOLOGY {
            return false;
        }
        let inferred_fragment_length = var.average_fragment_length() + var.adjusted_length();
        let length_threshold = self.fragment_length_bounds.upper_bound as f64 * DEL_ARTEFACT_LENGTH_FACTOR;
        if inferred_fragment_length as f64 >= length_threshold {
            return false;
        }
        !any_samples_above_af_threshold(var, DEL_ARTEFACT_MIN_AF)
    }

    fn low_three_prime_position_range(&self, var: &Variant) -> bool {
        if self.has_paired_reads() {
            return false;
        }
        let split_fragment_factor = if var.original_junctions.len() == 1 { PRIME_MAX_SGL_FACTOR } else { 1 };
        for breakend in var.breakends.iter().flatten() {
            let min_position_range = breakend.min_orientation_position_range();
            if min_position_range < 0 {
                continue;
            }
            let mut split_fragments = 0;
            let mut max_strand_bias = 0.0f64;
            for genotype in &breakend.context.genotypes {
                split_fragments += genotype.attr_int(SPLIT_FRAGS, 0);
                max_strand_bias = f64::max(adjusted_strand_bias(genotype), max_strand_bias);
            }
            if max_strand_bias > 0.0 {
                continue;
            }
            let max_permitted_range = i32::min(PRIME_MAX_BASE_FACTOR + split_fragments * split_fragment_factor, PRIME_MAX_PERMITTED_RANGE) as f64;
            if (min_position_range as f64) < max_permitted_range {
                return true;
            }
        }
        false
    }

    // SeqTechUtils

    fn is_sbx_strand_bias(&self, var: &Variant) -> bool {
        if !self.is_sbx() {
            return false;
        }
        if !matches!(var.sv_type, SvType::Bnd | SvType::Inv | SvType::Sgl) {
            return false;
        }
        for breakend in var.breakends.iter().flatten() {
            let mut split_fragments = 0;
            let mut max_strand_bias = 0.0f64;
            for genotype in &breakend.context.genotypes {
                split_fragments += genotype.attr_int(SPLIT_FRAGS, 0);
                max_strand_bias = f64::max(adjusted_strand_bias(genotype), max_strand_bias);
            }
            if max_strand_bias > 0.0 {
                return false;
            }
            if var.sv_type != SvType::Bnd && split_fragments < SBX_STRAND_BIAS_NON_BND_MIN_FRAGS {
                return false;
            }
        }
        true
    }

    fn is_sbx_artefact(&self, var: &Variant) -> bool {
        if !self.is_sbx() || var.is_hotspot {
            return false;
        }
        let mut max_strand_bias = 0.0f64;
        let full_assembly_length = var.start().context.attr_int(ASM_LENGTH, 0);
        let mut inexact_hom_length = 0;
        let mut local_junction_count = 0;

        for se in SE_START..=SE_END {
            let Some(breakend) = &var.breakends[se] else { continue };
            if var.close_to_original_junction(se, LOCAL_JUNCTION_MAX_DISTANCE) {
                local_junction_count += 1;
            }
            inexact_hom_length = inexact_hom_length.max(breakend.inexact_homology.length());
            for genotype in &breakend.context.genotypes {
                max_strand_bias = f64::max(adjusted_strand_bias(genotype), max_strand_bias);
            }
        }

        let non_line_stranded = !var.is_line_site() && max_strand_bias == 0.0;
        let mut heuristic_value = var.qual * f64::min(full_assembly_length as f64 / SBX_HEURISTIC_ASM_LENGTH_FACTOR, 1.0);
        let mut heuristic_threshold = SBX_HEURISTIC_THRESHOLD as f64;
        let inexact_hom_penalty = i32::min(inexact_hom_length * SBX_HEURISTIC_INEXACT_HOM_FACTOR, SBX_HEURISTIC_INEXACT_HOM_MAX) as f64;

        let mut var_length = 0;
        if var.sv_type.has_length() {
            var_length = var.sv_length();
            if var_length > SBX_HEURISTIC_SHORT_LENGTH {
                return false;
            }
        }

        match var.sv_type {
            SvType::Del => {
                if non_line_stranded {
                    return true;
                }
                heuristic_value -= inexact_hom_penalty;
                if local_junction_count < 2 {
                    heuristic_value -= SBX_HEURISTIC_LOCATION_JUNC_PENALTY as f64;
                }
                if var.is_line_site() {
                    heuristic_value += SBX_HEURISTIC_SGL_LINE_BONUS as f64;
                }
            }
            SvType::Dup | SvType::Ins => {
                if non_line_stranded {
                    return true;
                }
                heuristic_value -= inexact_hom_penalty;
                if local_junction_count < 2 {
                    heuristic_value -= SBX_HEURISTIC_LOCATION_JUNC_PENALTY as f64;
                }
                if var.is_line_site() {
                    heuristic_value += SBX_HEURISTIC_SGL_LINE_BONUS as f64;
                } else if var_length < SBX_HEURISTIC_DUP_INS_LENGTH {
                    heuristic_value -= SBX_HEURISTIC_SHORT_DUP_INS_PENALTY as f64;
                }
            }
            SvType::Inv => {
                if var_length == 0 {
                    return true;
                }
                if var_length < SBX_HEURISTIC_INV_SHORT_LENGTH && var.insert_sequence.len() as i32 > SBX_HEURISTIC_INV_INS_LENGTH {
                    return true;
                }
            }
            SvType::Bnd => {
                if local_junction_count < 2 {
                    heuristic_value -= SBX_HEURISTIC_BND_LOCATION_JUNC_PENALTY as f64;
                }
                let is_line_insert =
                    has_line_sequence(&var.insert_sequence, Orient::Forward) || has_line_sequence(&var.insert_sequence, Orient::Reverse);
                if is_line_insert {
                    heuristic_value += SBX_HEURISTIC_BND_LINE_BONUS as f64;
                } else if var.insert_sequence.len() as i32 >= LINE_POLY_AT_REQ {
                    heuristic_value -= SBX_HEURISTIC_BND_INS_PENALTY as f64;
                }
            }
            SvType::Sgl => {
                if non_line_stranded {
                    return true;
                }
                heuristic_threshold = SBX_HEURISTIC_SGL_THRESHOLD as f64;
                if local_junction_count == 0 {
                    heuristic_value -= SBX_HEURISTIC_LOCATION_JUNC_PENALTY as f64;
                }
                if var.is_line_site() {
                    heuristic_value += SBX_HEURISTIC_SGL_LINE_BONUS as f64;
                }
            }
            SvType::Inf => {}
        }
        heuristic_value < heuristic_threshold
    }

    fn is_sbx_zero_length_inversion(&self, var: &Variant) -> bool {
        if !self.is_sbx() || var.sv_type != SvType::Inv {
            return false;
        }
        if var.sv_length() != 0 {
            return false;
        }
        let insert_sequence = &var.insert_sequence;
        let length = insert_sequence.len() as i32;
        if length == 0 {
            return true;
        }
        if !(SBX_INV_INSERT_MOTIF_MIN_LENGTH..=SBX_INV_INSERT_MOTIF_MAX_LENGTH).contains(&length) {
            return false;
        }
        let reversed = reverse_complement(insert_sequence);
        let mut min_sequence_diffs = -1;
        for motif in SBX_INV_INSERT_MOTIFS {
            let length_diff = length - motif.len() as i32;
            for i in 0..=length_diff.abs() {
                let s1_offset = if length_diff > 0 { i } else { 0 };
                let s2_offset = if length_diff < 0 { i } else { 0 };
                let mut diffs = sequence_mismatches(insert_sequence, motif, s1_offset, s2_offset);
                diffs = diffs.min(sequence_mismatches(&reversed, motif, s1_offset, s2_offset));
                if min_sequence_diffs < 0 || diffs < min_sequence_diffs {
                    min_sequence_diffs = diffs;
                }
            }
        }
        min_sequence_diffs <= SBX_INV_INSERT_MOTIF_MAX_DIFFS
    }

    /// VariantFilters.applyAdjacentFilters: short INVs that are isolated from other variants
    pub fn apply_adjacent_filters(&self, variants: &mut [Variant], map: &BreakendMap) {
        for (_, breakends) in map.iter() {
            for i in 0..breakends.len() {
                let b = breakends[i];
                let var = &variants[b.var as usize];
                if var.filters.contains(FilterType::InvShortIsolated) {
                    continue;
                }
                if !is_candidate_isolated_inv(var) {
                    continue;
                }
                let position = var.breakend(b.se as usize).position;
                let mut has_valid_adjacent = false;
                for search_up in [true, false] {
                    let mut k = if search_up { i as i64 + 1 } else { i as i64 - 1 };
                    while k >= 0 && (k as usize) < breakends.len() {
                        let other = breakends[k as usize];
                        let other_position = variants[other.var as usize].breakend(other.se as usize).position;
                        if (other_position - position).abs() > INV_ADJACENT_LENGTH * 2 {
                            break;
                        }
                        if has_adjacent_non_artefact(variants, b, other) {
                            has_valid_adjacent = true;
                            break;
                        }
                        k += if search_up { 1 } else { -1 };
                    }
                    if has_valid_adjacent {
                        break;
                    }
                }
                if has_valid_adjacent {
                    continue;
                }
                variants[b.var as usize].add_filter(FilterType::InvShortIsolated);
            }
        }
    }
}

/// VariantFilters.applyWeakJunctionMinAfFilter
pub fn apply_weak_junction_min_af_filter(var: &Variant) -> bool {
    if var.is_saga_matched || var.is_hotspot {
        return false;
    }
    if !var.original_junctions.iter().any(|j| j.soft_clip_based()) {
        return false;
    }
    let breakend = var.start();
    if var.sv_type == SvType::Sgl || var.sv_type == SvType::Inv {
        let strand_bias = breakend.strand_bias();
        if strand_bias < MIN_AF_WEAK_JUNCTION_STRAND_BIAS_MIN || strand_bias > (1.0 - MIN_AF_WEAK_JUNCTION_STRAND_BIAS_MIN) {
            return true;
        }
    }
    if breakend.inexact_homology.length() >= INV_SHORT_MAX_HOMOLOGY_HIGHER {
        return true;
    }
    let proximate_junc_ratio = breakend.context.attr_double(PROX_JUNC_READ_RATIO, -1.0);
    if proximate_junc_ratio >= 0.0 && proximate_junc_ratio <= MIN_AF_WEAK_JUNCTION_PROX_JUNC_READ_THRESHOLD {
        return true;
    }
    let supp_remote_region_ratio = breakend.context.attr_double(SUPP_REMOTE_REGION_RATIO, -1.0);
    if supp_remote_region_ratio >= 0.0 && supp_remote_region_ratio <= MIN_AF_WEAK_JUNCTION_SUPP_REMOTE_THRESHOLD {
        return true;
    }
    false
}

/// both breakends of some sample must be at or above the threshold
fn any_samples_above_af_threshold(var: &Variant, af_threshold: f64) -> bool {
    let genotype_count = var.start().context.genotypes.len();
    for g in 0..genotype_count {
        let mut above_threshold = true;
        for breakend in var.breakends.iter().flatten() {
            let Some(genotype) = breakend.context.genotypes.get(g) else {
                panic!("IndexOutOfBoundsException: genotype {g} of {}", breakend.context.genotypes.len());
            };
            let af = var.calc_allelic_frequency(genotype);
            above_threshold &= af >= af_threshold;
        }
        if above_threshold {
            return true;
        }
    }
    false
}

fn fails_strand_bias(var: &Variant) -> bool {
    let breakend = var.start();
    let mut has_passing = false;
    let mut has_failing = false;
    for genotype in &breakend.context.genotypes {
        let fragment_count = Breakend::fragment_count_of(genotype);
        if fragment_count <= 1 {
            continue;
        }
        let strand_bias = genotype.attr_double(STRAND_BIAS, 0.5);
        let ins_sequence = &var.insert_sequence;
        if strand_bias == 1.0 && (ins_sequence.contains(PON_INS_SEQ_FWD_STRAND_1) || ins_sequence.contains(PON_INS_SEQ_FWD_STRAND_2)) {
            has_failing = true;
        } else if strand_bias == 0.0
            && (ins_sequence.contains(&reverse_complement(PON_INS_SEQ_FWD_STRAND_1)) || ins_sequence.contains(&reverse_complement(PON_INS_SEQ_FWD_STRAND_2)))
        {
            has_failing = true;
        } else {
            has_passing = true;
        }
    }
    has_failing && !has_passing
}

fn is_candidate_isolated_inv(var: &Variant) -> bool {
    var.sv_type == SvType::Inv
        && var.sv_length() <= INV_ADJACENT_LENGTH
        && (!var.in_chained_assembly() || var.max_unique_fragment_positions() < INV_ADJACENT_MIN_UPS)
}

fn has_adjacent_non_artefact(variants: &[Variant], b: BreakendRef, other: BreakendRef) -> bool {
    let var = &variants[b.var as usize];
    if !var.is_sgl() && other == b.other() {
        return false;
    }
    let breakend = var.breakend(b.se as usize);
    let other_var = &variants[other.var as usize];
    let other_breakend = other_var.breakend(other.se as usize);
    if breakend.orient == other_breakend.orient {
        return false;
    }
    if is_candidate_isolated_inv(other_var) {
        return false;
    }
    if INV_ADJACENT_EXCLUDED_FILTERS.iter().any(|f| other_var.filters.contains(*f)) {
        return false;
    }
    let (max_lower_position, min_upper_position) = if breakend.position < other_breakend.position {
        (breakend.position + breakend.inexact_homology.end, other_breakend.position + other_breakend.inexact_homology.start)
    } else {
        (other_breakend.position + other_breakend.inexact_homology.end, breakend.position + breakend.inexact_homology.start)
    };
    max_lower_position >= min_upper_position - INV_ADJACENT_LENGTH
}

// ---------------------------------------------------------------------------------------------
// SBX constants

const SBX_STRAND_BIAS_NON_BND_MIN_FRAGS: i32 = 10;
const SBX_HEURISTIC_THRESHOLD: i32 = 50;
const SBX_HEURISTIC_SGL_THRESHOLD: i32 = 30;
const SBX_HEURISTIC_SHORT_LENGTH: i32 = 300;
const SBX_HEURISTIC_ASM_LENGTH_FACTOR: f64 = 600.0;
const SBX_HEURISTIC_INEXACT_HOM_FACTOR: i32 = 2;
const SBX_HEURISTIC_INEXACT_HOM_MAX: i32 = 50;
const SBX_HEURISTIC_DUP_INS_LENGTH: i32 = 40;
const SBX_HEURISTIC_SHORT_DUP_INS_PENALTY: i32 = 20;
const SBX_HEURISTIC_INV_SHORT_LENGTH: i32 = 100;
const SBX_HEURISTIC_INV_INS_LENGTH: i32 = 20;
const SBX_HEURISTIC_LOCATION_JUNC_PENALTY: i32 = 30;
const SBX_HEURISTIC_BND_LOCATION_JUNC_PENALTY: i32 = 20;
const SBX_HEURISTIC_BND_INS_PENALTY: i32 = 20;
const SBX_HEURISTIC_BND_LINE_BONUS: i32 = 20;
const SBX_HEURISTIC_SGL_LINE_BONUS: i32 = 30;
const SBX_HEURISTIC_BND_LINE_PERC: f64 = 0.9;

pub const SBX_INV_INSERT_MOTIFS: [&str; 16] = [
    "ACCATCAATCCGGATGTATGCCGGATTGATGGT",
    "AGGAGTAACATCCATGTATGGGATGTTACTCCT",
    "ACCCGTAGTATTAATGTATGTAATACTACGGGT",
    "ACTATTGAAGGCTATGTATGAGCCTTCAATAGT",
    "AGGATTATGGCAAATGTATGTTGCCATAATCCT",
    "AGGCGAGTATTCCATACATGGAATACTCGCCT",
    "AGGTTATCAGCTTATGTATGAAGCTGATAACCT",
    "AGGGACATTACCAATGTATGTGGTAATGTCCCT",
    "ACAACAGCCGAAGATGTATGCTTCGGCTGTTGT",
    "AGGTATTGCCTGGATGTATGCCAGGCAATACCT",
    "AGAATCGCAGTTAATGTATGTAACTGCGATTCT",
    "ACCTGCGGTTGAAATGTATGTTCAACCGCAGGT",
    "ACCACTTCGCCAAATGTATGTTGGCGAAGTGGT",
    "AGGCTTACTACGGATGTATGCCGTAGTAAGCCT",
    "TCTTCTGTATTTCTAGTGTTCAATAGAAATACAGA",
    "ACATAGAAGGTAGATGTATGCTACCTTCTATGT",
];
const SBX_INV_INSERT_MOTIF_MIN_LENGTH: i32 = 29;
const SBX_INV_INSERT_MOTIF_MAX_LENGTH: i32 = 34;
const SBX_INV_INSERT_MOTIF_MAX_DIFFS: i32 = 10;

/// SeqTechUtils.sequenceMismatches
pub fn sequence_mismatches(s1: &str, s2: &str, s1_offset: i32, s2_offset: i32) -> i32 {
    compare_sequences_mismatch_no_quals(
        s1.as_bytes(),
        s1_offset,
        s1.len() as i32 - 1,
        s2.as_bytes(),
        s2_offset,
        s2.len() as i32 - 1,
        SBX_INV_INSERT_MOTIF_MAX_DIFFS,
    )
}

const MISMATCH_RECOVERY_BASE_COUNT: i32 = 3;

/// SequenceCompare.canRecoverMatch, forwards
fn can_recover_match(first: &[u8], first_current: i32, second: &[u8], second_current: i32) -> bool {
    let first_end = first_current + MISMATCH_RECOVERY_BASE_COUNT;
    let second_end = second_current + MISMATCH_RECOVERY_BASE_COUNT;
    if first_end >= first.len() as i32 || second_end >= second.len() as i32 {
        return false;
    }
    let (mut i, mut j) = (first_current, second_current);
    loop {
        if first[i as usize] != second[j as usize] {
            return false;
        }
        if i == first_end && j == second_end {
            break;
        }
        i += 1;
        j += 1;
    }
    true
}

/// SequenceCompare.checkSkipShortIndel, forwards
fn check_skip_short_indel(first: &[u8], first_index: i32, second: &[u8], second_index: i32) -> i32 {
    for diff in 1..=2 {
        if can_recover_match(first, first_index + diff, second, second_index) {
            return diff;
        }
    }
    for diff in 1..=2 {
        if can_recover_match(first, first_index, second, second_index + diff) {
            return -diff;
        }
    }
    0
}

/// SequenceCompare.compareSequencesMismatchNoQuals
pub fn compare_sequences_mismatch_no_quals(
    first: &[u8], first_start: i32, first_end: i32, second: &[u8], second_start: i32, second_end: i32, max_mismatch_penalty: i32,
) -> i32 {
    if first_start < 0 || first_end >= first.len() as i32 || second_start < 0 || second_end >= second.len() as i32 {
        return if max_mismatch_penalty < 0 { max_mismatch_penalty } else { max_mismatch_penalty + 1 };
    }
    let mut mismatch_total = 0;
    let mut first_index = first_start;
    let mut second_index = second_start;
    // never set in the Java either, so a skip only cancels out when it is zero
    let last_repeat_skip_bases = 0;

    while first_index >= 0 && first_index <= first_end && second_index >= 0 && second_index <= second_end {
        if first[first_index as usize] == second[second_index as usize] {
            first_index += 1;
            second_index += 1;
            continue;
        }
        let recovery_skip_bases = check_skip_short_indel(first, first_index, second, second_index);
        if recovery_skip_bases != 0 {
            if recovery_skip_bases > 0 {
                first_index += recovery_skip_bases;
            } else {
                second_index += -recovery_skip_bases;
            }
            if last_repeat_skip_bases + recovery_skip_bases == 0 {
                continue;
            }
            mismatch_total += 1;
            continue;
        }
        mismatch_total += 1;
        if max_mismatch_penalty >= 0 && mismatch_total > max_mismatch_penalty {
            return mismatch_total;
        }
        first_index += 1;
        second_index += 1;
    }
    mismatch_total
}

// ---------------------------------------------------------------------------------------------
// repeats (RepeatInfo)

const REPEAT_MIN_COUNT: i32 = 4;

/// (index, repeat unit length, count)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RepeatInfo {
    pub index: i32,
    pub unit_length: i32,
    pub count: i32,
}

impl RepeatInfo {
    fn total_length(&self) -> i32 {
        self.count * self.unit_length
    }
}

fn find_multi_base_repeat(bases: &[u8], index: usize, repeat_length: usize, min_count: i32) -> Option<RepeatInfo> {
    if index + min_count as usize * repeat_length - 1 >= bases.len() {
        return None;
    }
    let mut i = index + repeat_length;
    while i + (repeat_length - 1) < bases.len() {
        if bases[i..i + repeat_length] != bases[index..index + repeat_length] {
            break;
        }
        i += repeat_length;
    }
    let count = ((i - index) / repeat_length) as i32;
    if count < min_count {
        return None;
    }
    Some(RepeatInfo { index: index as i32, unit_length: repeat_length as i32, count })
}

fn find_dual_dual_repeat(bases: &[u8], index: usize, min_count: i32) -> Option<RepeatInfo> {
    if index + min_count as usize * 4 - 1 >= bases.len() {
        return None;
    }
    if bases[index] != bases[index + 1] || bases[index + 2] != bases[index + 3] {
        return None;
    }
    find_multi_base_repeat(bases, index, 4, min_count)
}

/// RepeatInfo.findRepeats: single, dual, 3, 4 and 5 base repeats, then pairs of doubled bases
pub fn find_repeats(bases: &[u8]) -> Vec<RepeatInfo> {
    let mut repeats = Vec::new();
    let mut index: i64 = 0;
    while index <= bases.len() as i64 - REPEAT_MIN_COUNT as i64 {
        let i = index as usize;
        // a single-base repeat counts every following equal base, the others whole units
        let single = {
            if i + REPEAT_MIN_COUNT as usize - 1 >= bases.len() {
                None
            } else {
                let mut j = i + 1;
                while j < bases.len() && bases[j] == bases[i] {
                    j += 1;
                }
                let count = (j - i) as i32;
                if count >= REPEAT_MIN_COUNT { Some(RepeatInfo { index: i as i32, unit_length: 1, count }) } else { None }
            }
        };
        let repeat = single
            .or_else(|| find_multi_base_repeat(bases, i, 2, REPEAT_MIN_COUNT))
            .or_else(|| find_multi_base_repeat(bases, i, 3, REPEAT_MIN_COUNT))
            .or_else(|| find_multi_base_repeat(bases, i, 4, REPEAT_MIN_COUNT))
            .or_else(|| find_multi_base_repeat(bases, i, 5, REPEAT_MIN_COUNT))
            .or_else(|| find_dual_dual_repeat(bases, i, REPEAT_MIN_COUNT));
        match repeat {
            Some(r) => {
                index += r.total_length() as i64;
                repeats.push(r);
            }
            None => index += 1,
        }
    }
    repeats
}

/// RepeatInfo.calcTrimmedBaseLength: each repeat counts as two of its units
pub fn calc_trimmed_base_length(seq_start: i32, seq_end: i32, repeats: &[RepeatInfo]) -> i32 {
    if repeats.is_empty() {
        return seq_end - seq_start + 1;
    }
    let mut current_index = seq_start;
    let mut current_repeat_index = 0;
    let mut trimmed = 0;
    let mut current_repeat = Some(repeats[0]);
    while current_index <= seq_end {
        while let Some(r) = current_repeat {
            if r.index >= current_index {
                break;
            }
            current_repeat_index += 1;
            current_repeat = repeats.get(current_repeat_index).copied();
        }
        match current_repeat {
            Some(r) if r.index == current_index => {
                current_index += r.total_length();
                trimmed += r.unit_length * 2;
            }
            _ => {
                current_index += 1;
                trimmed += 1;
            }
        }
    }
    trimmed
}

// ---------------------------------------------------------------------------------------------
// LINE sites (LineChecker)

/// LineUtils.findLineSequenceCount
fn find_line_sequence_count(bases: &[u8], index_start: i32, index_end: i32, line_base: u8) -> i32 {
    if index_start < 0 || index_end >= bases.len() as i32 {
        return 0;
    }
    if index_end - index_start + 1 < LINE_POLY_AT_REQ {
        return 0;
    }
    let mut line_base_count = 0;
    let mut other_count = 0;
    for i in index_start..=index_end {
        if bases[i as usize] == line_base {
            line_base_count += 1;
        } else {
            other_count += 1;
            if other_count > LINE_POLY_AT_TEST_LEN - LINE_POLY_AT_REQ {
                break;
            }
        }
    }
    if line_base_count >= LINE_POLY_AT_REQ { line_base_count } else { 0 }
}

/// LineChecker.hasLineSequence
pub fn has_line_sequence(insert_sequence: &str, orientation: Orient) -> bool {
    let line_base = if orientation.is_forward() { b'T' } else { b'A' };
    let count = find_line_sequence_count(insert_sequence.as_bytes(), 0, insert_sequence.len() as i32 - 1, line_base);
    let line_perc = count as f64 / insert_sequence.len() as f64;
    count >= LINE_POLY_AT_REQ && line_perc >= SBX_HEURISTIC_BND_LINE_PERC
}

/// CommonUtils.withinLineProximity
fn within_line_proximity(pos1: i32, pos2: i32, orient1: Orient, orient2: Orient) -> bool {
    if orient1 == orient2 {
        return false;
    }
    let pos_diff = (pos2 - pos1).abs();
    if (pos1 < pos2) == orient1.is_forward() { pos_diff <= LINE_INDEL_MAX_GAP } else { pos_diff <= LINE_INDEL_MAX_OVERLAP }
}

/// LineChecker.markLineSites
pub fn mark_line_sites(variants: &mut [Variant], map: &BreakendMap) {
    let max_proximate_line_distance = LINE_INDEL_MAX_OVERLAP.max(LINE_INDEL_MAX_GAP);

    for (_, breakends) in map.iter() {
        if breakends.is_empty() {
            continue;
        }
        for i in 0..breakends.len() - 1 {
            let b = breakends[i];
            let (bv, bse) = (b.var as usize, b.se as usize);

            // first search for likely remote line source sites
            for &next in &breakends[i + 1..] {
                let (nv, nse) = (next.var as usize, next.se as usize);
                let breakend = variants[bv].breakend(bse);
                let next_breakend = variants[nv].breakend(nse);
                if (breakend.position - next_breakend.position).abs() > max_proximate_line_distance {
                    break;
                }
                if !variants[bv].is_sgl() && b.other() == next {
                    continue;
                }
                if !breakend.is_line() && !next_breakend.is_line() {
                    continue;
                }
                if !variants[bv].is_line_site() && !variants[bv].close_to_original_junction(bse, max_proximate_line_distance) {
                    variants[bv].breakend_mut(bse).line_remote_source_site = true;
                } else if !variants[nv].is_line_site() && !variants[nv].close_to_original_junction(nse, max_proximate_line_distance) {
                    variants[nv].breakend_mut(nse).line_remote_source_site = true;
                }
            }

            // otherwise check for pairs of line and the other insertion breakend
            let next = breakends[i + 1];
            let (nv, nse) = (next.var as usize, next.se as usize);
            if !variants[bv].is_sgl() && b.other() == next {
                continue;
            }
            let breakend = variants[bv].breakend(bse);
            let next_breakend = variants[nv].breakend(nse);
            if !breakend.is_line() && !next_breakend.is_line() {
                continue;
            }
            if breakend.line_remote_source_site || next_breakend.line_remote_source_site {
                continue;
            }
            if !within_line_proximity(breakend.position, next_breakend.position, breakend.orient, next_breakend.orient) {
                continue;
            }
            variants[bv].breakend_mut(bse).line_site_breakend = Some(next);
            variants[nv].breakend_mut(nse).line_site_breakend = Some(b);
        }
    }
}

/// LineChecker.adjustLineSites
pub fn adjust_line_sites(variants: &mut [Variant], v: usize) {
    let var = &variants[v];
    if !var.is_line_site() {
        return;
    }
    let is_germline = var.germline;
    let is_filtered = var.is_filtered();
    let mut has_line_germline = false;
    let mut has_line_pass = false;
    for breakend in var.breakends.iter().flatten() {
        if let Some(line) = breakend.line_site_breakend {
            has_line_germline |= variants[line.var as usize].germline;
            has_line_pass |= variants[line.var as usize].is_pass();
        }
    }
    if !is_germline && has_line_germline {
        variants[v].germline = true;
    }
    if is_filtered && has_line_pass {
        variants[v].filters.clear();
    }
}

// ---------------------------------------------------------------------------------------------
// Deduplication

fn is_duplicate(variants: &[Variant], a: BreakendRef, b: BreakendRef) -> bool {
    let (first_var, second_var) = (&variants[a.var as usize], &variants[b.var as usize]);
    let first = first_var.breakend(a.se as usize);
    let second = second_var.breakend(b.se as usize);
    if first.position != second.position || first.orient != second.orient || first.insert_sequence != second.insert_sequence {
        return false;
    }
    if first_var.is_sgl() != second_var.is_sgl() {
        return false;
    }
    if first_var.is_sgl() {
        return true;
    }
    let first_other = first_var.breakend(1 - a.se as usize);
    let second_other = second_var.breakend(1 - b.se as usize);
    first_other.position == second_other.position && first_other.orient == second_other.orient && first_other.chromosome == second_other.chromosome
}

fn keep_first(variants: &[Variant], a: BreakendRef, b: BreakendRef) -> bool {
    let (first_var, second_var) = (&variants[a.var as usize], &variants[b.var as usize]);
    let first_support = first_var.breakend(a.se as usize).fragment_count();
    let second_support = second_var.breakend(b.se as usize).fragment_count();
    if first_support != second_support {
        return first_support > second_support;
    }
    first_var.qual >= second_var.qual
}

/// Deduplication.deduplicateVariants
pub fn deduplicate_variants(variants: &mut [Variant], map: &BreakendMap) {
    for (_, breakends) in map.iter() {
        if breakends.is_empty() {
            continue;
        }
        for i in 0..breakends.len() - 1 {
            let b = breakends[i];
            if variants[b.var as usize].is_filtered() {
                continue;
            }
            for &next in &breakends[i + 1..] {
                if variants[next.var as usize].is_filtered() {
                    continue;
                }
                let breakend = variants[b.var as usize].breakend(b.se as usize);
                let next_breakend = variants[next.var as usize].breakend(next.se as usize);
                if breakend.chromosome != next_breakend.chromosome || next_breakend.position > breakend.position {
                    break;
                }
                if is_duplicate(variants, b, next) {
                    if keep_first(variants, b, next) {
                        variants[next.var as usize].add_filter(FilterType::Duplicate);
                    } else {
                        variants[b.var as usize].add_filter(FilterType::Duplicate);
                    }
                }
            }
        }
    }
}
