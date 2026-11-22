//! AlignData: one BWA alignment of an assembly sequence, with the adjusted values derived from it.

use crate::assembly::alignment::alternative_alignment::AlternativeAlignment;
use crate::assembly::alignment::bwa::BwaMemAlignment;
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::P;
use crate::assembly::sam::*;
use crate::assembly::types::assembly_stats::java_round;
use crate::assembly::types::remote_region::compare_chr_regions;
use crate::assembly::types::repeat_info::{calc_trimmed_base_length, find_repeats};

pub type AlignP = P<AlignData>;

pub const HUMAN_CHROMOSOME_NAMES: [&str; 24] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y",
];

/// (chromosome, start, end) of the DUX4 regions where the default alignment is listed with the alternatives
const MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37: [(&str, i32, i32); 2] = [("4", 190930000, 191030000), ("10", 135420000, 135520000)];
const MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38: [(&str, i32, i32); 2] = [("chr4", 190060000, 190190000), ("chr10", 133660000, 133770000)];

pub struct AlignData {
    chromosome: String,
    position_start: i32,
    position_end: i32,
    map_qual: i32,
    n_matches: i32,
    score: i32,
    flags: i32,
    cigar: String,
    xa_tag: String,
    md_tag: String,
    cigar_elements: Vec<CigarElement>,
    soft_clip_left: i32,
    soft_clip_right: i32,
    orientation: Orientation,
    aligned_bases: i32,
    raw_sequence_start: i32,
    raw_sequence_end: i32,
    sequence_start: i32,
    sequence_end: i32,
    dropped_on_requery: bool,
    has_sequence_coords_set: bool,
    is_requery: bool,
    adjusted_alignment: i32,
    modified_map_qual: i32,
    raw_alt_alignments: Option<Vec<AlternativeAlignment>>,
    selected_alt_alignment: Option<AlternativeAlignment>,
    unselected_alt_alignments: Option<Vec<AlternativeAlignment>>,
    has_low_map_qual_short_sv_link: bool,
}

impl AlignData {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        chromosome: &str, position_start: i32, position_end: i32, sequence_start: i32, sequence_end: i32, map_qual: i32, score: i32, flags: i32, cigar: &str,
        n_matches: i32, xa_tag: Option<&str>, md_tag: Option<&str>,
    ) -> AlignData {
        let cigar_elements = cigar_elements_from_str(cigar);
        let (soft_clip_left, soft_clip_right) = if !cigar_elements.is_empty() {
            (left_soft_clip_length(&cigar_elements), right_soft_clip_length(&cigar_elements))
        } else {
            (0, 0)
        };
        let aligned_bases: i32 = cigar_elements.iter().filter(|x| x.op == CigarOp::M).map(|x| x.len).sum();
        AlignData {
            chromosome: chromosome.to_string(),
            position_start,
            position_end,
            raw_sequence_start: sequence_start,
            raw_sequence_end: sequence_end,
            map_qual,
            n_matches,
            score,
            flags,
            cigar: cigar.to_string(),
            xa_tag: xa_tag.unwrap_or("").to_string(),
            md_tag: md_tag.unwrap_or("").to_string(),
            cigar_elements,
            soft_clip_left,
            soft_clip_right,
            orientation: if flags & FLAG_READ_REVERSE_STRAND != 0 { Orientation::Reverse } else { Orientation::Forward },
            aligned_bases,
            sequence_start,
            sequence_end: (sequence_end - 1).max(0),
            has_sequence_coords_set: false,
            is_requery: false,
            dropped_on_requery: false,
            adjusted_alignment: aligned_bases,
            modified_map_qual: 0,
            raw_alt_alignments: None,
            selected_alt_alignment: None,
            unselected_alt_alignments: None,
            has_low_map_qual_short_sv_link: false,
        }
    }

    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }
    pub fn position_start(&self) -> i32 {
        self.position_start
    }
    pub fn position_end(&self) -> i32 {
        self.position_end
    }
    pub fn ref_location_str(&self) -> String {
        format!("{}:{}-{}", self.chromosome, self.position_start, self.position_end)
    }
    pub fn map_qual(&self) -> i32 {
        self.map_qual
    }
    pub fn n_matches(&self) -> i32 {
        self.n_matches
    }
    pub fn score(&self) -> i32 {
        self.score
    }
    pub fn flags(&self) -> i32 {
        self.flags
    }
    pub fn cigar(&self) -> &str {
        &self.cigar
    }
    pub fn xa_tag(&self) -> &str {
        &self.xa_tag
    }
    pub fn md_tag(&self) -> &str {
        &self.md_tag
    }
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }
    pub fn is_forward(&self) -> bool {
        self.orientation.is_forward()
    }
    pub fn is_reverse(&self) -> bool {
        self.orientation.is_reverse()
    }
    pub fn is_supplementary(&self) -> bool {
        self.flags & FLAG_SUPPLEMENTARY != 0
    }
    pub fn left_soft_clip_length(&self) -> i32 {
        self.soft_clip_left
    }
    pub fn right_soft_clip_length(&self) -> i32 {
        self.soft_clip_right
    }
    pub fn aligned_bases(&self) -> i32 {
        self.aligned_bases
    }
    pub fn segment_length(&self) -> i32 {
        self.sequence_end - self.sequence_start + 1
    }
    pub fn set_soft_clip_lengths(&mut self, left: i32, right: i32) {
        self.soft_clip_left = left;
        self.soft_clip_right = right;
    }

    pub fn set_full_sequence_data(&mut self, full_sequence: &str, full_sequence_length: i32) {
        if self.has_sequence_coords_set {
            return;
        }
        if self.orientation.is_reverse() {
            let new_sequence_start = (full_sequence_length - 1) - (self.raw_sequence_end - 1);
            let new_sequence_end = (full_sequence_length - 1) - self.sequence_start;
            self.sequence_start = new_sequence_start.max(0);
            self.sequence_end = new_sequence_end;
        }
        if self.sequence_start < 0 || self.sequence_start > self.sequence_end || self.sequence_end > full_sequence.len() as i32 {
            crate::error!("alignment({}) invalid subsequence request({}-{}) vs fullSequenceLength({})", self, self.sequence_start, self.sequence_end, full_sequence_length);
        }
        self.has_sequence_coords_set = true;
    }

    pub fn mark_dropped_on_requery(&mut self) {
        self.dropped_on_requery = true;
    }
    pub fn dropped_on_requery(&self) -> bool {
        self.dropped_on_requery
    }
    pub fn is_requery(&self) -> bool {
        self.is_requery
    }
    pub fn set_requeried_sequence_coords(&mut self, sequence_start: i32, sequence_end: i32) {
        self.sequence_start = sequence_start;
        self.sequence_end = sequence_end;
        self.has_sequence_coords_set = true;
        self.is_requery = true;
    }
    pub fn sequence_start(&self) -> i32 {
        self.sequence_start
    }
    pub fn sequence_end(&self) -> i32 {
        self.sequence_end
    }
    pub fn raw_sequence_start(&self) -> i32 {
        self.raw_sequence_start
    }
    pub fn raw_sequence_end(&self) -> i32 {
        self.raw_sequence_end
    }
    pub fn cigar_elements(&self) -> &Vec<CigarElement> {
        &self.cigar_elements
    }
    pub fn adjusted_alignment(&self) -> i32 {
        self.adjusted_alignment
    }
    pub fn modified_map_qual(&self) -> f64 {
        self.modified_map_qual as f64
    }

    fn calc_adjusted_indel_score(&self) -> i32 {
        if self.map_qual < ALIGNMENT_INDEL_RECOVERY_MIN_MAP_QUAL {
            return 0;
        }
        let mut indel_length = 0;
        for element in &self.cigar_elements {
            if element.op == CigarOp::S {
                return 0;
            }
            if element.op.is_indel() && element.len >= MIN_VARIANT_LENGTH {
                if indel_length != 0 {
                    return 0;
                }
                indel_length = element.len;
                if element.op == CigarOp::D {
                    indel_length = -indel_length;
                }
            }
        }
        indel_length
    }

    pub fn set_adjusted_alignment(&mut self, full_sequence: &str, inexact_homology_start: i32, inexact_homology_end: i32) {
        let mut sequence_length = self.sequence_end - self.sequence_start + 1;
        let mut adjusted_score = self.score;
        let adjusted_indel_score = self.calc_adjusted_indel_score();
        if adjusted_indel_score != 0 {
            if adjusted_indel_score < 0 {
                sequence_length += adjusted_indel_score;
            } else {
                sequence_length -= adjusted_indel_score;
            }
            adjusted_score += adjusted_indel_score.abs() + BWA_GAP_OPEN_PENALTY;
        }
        let score_adjustment = sequence_length - adjusted_score;
        let seq_start = if inexact_homology_start > 0 { self.sequence_start + inexact_homology_start } else { self.sequence_start };
        let seq_end = if inexact_homology_end > 0 { self.sequence_end - inexact_homology_end } else { self.sequence_end };
        if seq_start < 0 || seq_start > seq_end || seq_end >= full_sequence.len() as i32 {
            self.adjusted_alignment = 0;
            return;
        }
        let aligned_bases = &full_sequence.as_bytes()[seq_start as usize..(seq_end + 1) as usize];
        let repeats = find_repeats(aligned_bases);
        let trimmed_alignment_length = calc_trimmed_base_length(0, aligned_bases.len() as i32 - 1, &repeats);
        self.adjusted_alignment = (trimmed_alignment_length - score_adjustment).max(0);
        let inexact_homology_length = inexact_homology_start + inexact_homology_end;
        let adjusted_sequence_length = (sequence_length - inexact_homology_length).max(0) as f64;
        if (adjusted_score + ALIGNMENT_CALC_SCORE_FACTOR) as f64 / adjusted_sequence_length < ALIGNMENT_CALC_SCORE_THRESHOLD {
            self.modified_map_qual = 0;
        } else {
            let length_factor = (self.adjusted_alignment as f64 / 100.max(self.aligned_bases - inexact_homology_length) as f64).powi(2);
            self.modified_map_qual = java_round(self.map_qual as f64 * length_factor.min(1.0)) as i32;
        }
    }

    pub fn has_alt_alignments(&self) -> bool {
        !self.xa_tag.is_empty()
    }

    pub fn raw_alt_alignments(&mut self) -> Vec<AlternativeAlignment> {
        if self.xa_tag.is_empty() {
            return Vec::new();
        }
        if let Some(a) = &self.raw_alt_alignments {
            return a.clone();
        }
        let a = AlternativeAlignment::from_location_tag(&self.xa_tag);
        self.raw_alt_alignments = Some(a.clone());
        a
    }

    pub fn all_alignments(&mut self) -> Vec<AlternativeAlignment> {
        let alt_alignments = self.raw_alt_alignments();
        let default_alignment = AlternativeAlignment::new(&self.chromosome, self.position_start, self.orientation, &self.cigar, self.map_qual);
        if alt_alignments.is_empty() {
            let regions: &[(&str, i32, i32)] = if self.chromosome.starts_with("chr") { &MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38 } else { &MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37 };
            if !regions.iter().any(|(c, s, e)| *c == self.chromosome && *s <= self.position_end && *e >= self.position_start) {
                return alt_alignments;
            }
        }
        let mut all_alignments = Vec::with_capacity(alt_alignments.len() + 1);
        all_alignments.push(default_alignment);
        all_alignments.extend(alt_alignments);
        all_alignments
    }

    pub fn has_low_map_qual_alignment(&self) -> bool {
        self.selected_alt_alignment.is_some() || self.unselected_alt_alignments.is_some()
    }
    pub fn has_selected_alt_alignment(&self) -> bool {
        self.selected_alt_alignment.is_some()
    }
    pub fn has_low_map_qual_short_sv_link(&self) -> bool {
        self.has_low_map_qual_short_sv_link
    }
    pub fn selected_alt_alignment(&self) -> Option<&AlternativeAlignment> {
        self.selected_alt_alignment.as_ref()
    }
    pub fn set_selected_alt_alignments(&mut self, selected_alignment: Option<AlternativeAlignment>, other_alignments: Option<Vec<AlternativeAlignment>>, has_short_sv_link: bool) {
        self.selected_alt_alignment = selected_alignment;
        self.unselected_alt_alignments = other_alignments;
        self.has_low_map_qual_short_sv_link = has_short_sv_link;
    }
    pub fn set_selected_low_map_qual_alt_alignment(&mut self, selected_alignment: Option<AlternativeAlignment>) {
        self.selected_alt_alignment = selected_alignment;
        self.modified_map_qual = ALIGNMENT_MIN_MOD_MAP_QUAL;
        self.raw_alt_alignments = None;
        self.soft_clip_right = 0;
        self.soft_clip_left = 0;
    }
    pub fn set_modified_map_qual(&mut self, map_qual: i32) {
        self.modified_map_qual = map_qual;
    }
    pub fn unselected_alt_alignments(&self) -> &[AlternativeAlignment] {
        self.unselected_alt_alignments.as_deref().unwrap_or(&[])
    }

    pub fn is_lower_alignment(&self, other: &AlignData) -> bool {
        compare_chr_regions(&self.chromosome, self.position_start, &other.chromosome, other.position_start) == std::cmp::Ordering::Less
    }

    /// AlignData.fromRef: alignments to anything but the first 24 contigs of the index are dropped
    pub fn from_ref(alignment: &BwaMemAlignment, v38: bool) -> Option<AlignData> {
        let chr_index = alignment.ref_id;
        if chr_index < 0 || chr_index as usize >= HUMAN_CHROMOSOME_NAMES.len() {
            return None;
        }
        let name = HUMAN_CHROMOSOME_NAMES[chr_index as usize];
        let chromosome = if v38 { format!("chr{name}") } else { name.to_string() };
        Some(AlignData::new(
            &chromosome,
            alignment.ref_start + 1,
            alignment.ref_end,
            alignment.seq_start,
            alignment.seq_end,
            alignment.map_qual,
            alignment.aligner_score,
            alignment.sam_flag,
            &alignment.cigar,
            alignment.n_mismatches,
            alignment.xa_tag.as_deref(),
            alignment.md_tag.as_deref(),
        ))
    }
}

impl AlignData {
    /// AlignData.fromSaga: the alignment of an assembly to a SAGA variant assembly, located on the variant id
    pub fn from_saga(alignment: &crate::assembly::saga::SagaAlignment) -> AlignData {
        let raw = &alignment.raw_alignment;
        AlignData::new(
            &alignment.saga_variant().id,
            alignment.saga_start() + 1,
            alignment.saga_end(),
            alignment.query_start(),
            alignment.query_end(),
            raw.map_qual,
            raw.aligner_score,
            raw.sam_flag,
            &raw.cigar,
            raw.n_mismatches,
            raw.xa_tag.as_deref(),
            raw.md_tag.as_deref(),
        )
    }
}

impl std::fmt::Display for AlignData {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}:{} {} seq({}-{} adj={}-{} {}) score({}) flags({}) mapQual({} adj={}) aligned({} adj={}) md({})",
            self.ref_location_str(),
            self.orientation.as_byte(),
            self.cigar,
            self.raw_sequence_start,
            self.raw_sequence_end,
            self.sequence_start,
            self.sequence_end,
            self.orientation.as_byte(),
            self.score,
            self.flags,
            self.map_qual,
            self.modified_map_qual,
            self.aligned_bases,
            self.adjusted_alignment,
            self.md_tag
        )
    }
}
