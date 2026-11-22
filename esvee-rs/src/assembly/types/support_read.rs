//! SupportRead: the properties of a read kept by an assembly it supports.

use crate::assembly::common_utils::IndelCoords;
use crate::assembly::constants::MIN_INDEL_LENGTH;
use crate::assembly::globals::is_sbx;
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::A;
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::is_discordant_fragment;
use crate::assembly::read_parse_state::ReadParseState;
use crate::assembly::sam::*;
use crate::assembly::seq_tech_utils::set_sbx_medium_qual_count;
use crate::assembly::types::SupportType;

pub type SupP = A<SupportRead>;

/// value of the reference mismatch count before it is set (Java null)
const NO_REF_BASE_MISMATCHES: i32 = i32::MIN;

pub struct SupportRead {
    type_: SupportType,
    id: Box<str>,
    chromosome: &'static str,
    alignment_start: i32,
    alignment_end: i32,
    cigar: Box<str>,
    unclipped_start: i32,
    unclipped_end: i32,
    mate_chromosome: Option<&'static str>,
    mate_alignment_start: i32,
    mate_alignment_end: i32,
    mate_fragment_end: i32,
    flags: i32,
    base_length: i32,
    is_discordant: bool,
    is_reference: bool,
    sample_index: i32,
    supplementary_data: Option<Box<SupplementaryReadData>>,
    map_qual: i32,
    insert_size: i32,
    trim_count_start: i32,
    trim_count_end: i32,
    has_indel: bool,
    indel_coords: Option<A<IndelCoords>>,
    has_line_tail: bool,
    junction_read_start_distance: i32,
    full_assembly_index: i32,
    full_assembly_orientation: Orientation,
    inferred_fragment_length: i32,
    breakend_type: Option<SupportType>,
    ext_base_matches: i32,
    ext_base_mismatches: i32,
    ext_base_overlap: i32,
    medium_qual_count: i32,
    /// NO_REF_BASE_MISMATCHES until set
    ref_base_mismatches: i32,
    mismatch_info: Box<str>,
    read: Option<ReadP>,
}

impl Drop for SupportRead {
    fn drop(&mut self) {
        crate::assembly::mem_stats::SUPPORT_READS_LIVE.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
    }
}

impl SupportRead {
    pub fn new(read: &ReadP, type_: SupportType, junct_read_start_distance: i32, matches: i32, mismatches: i32) -> SupP {
        let (mate_chromosome, mate_alignment_start, mate_alignment_end, mate_fragment_end) = if read.is_paired_read() {
            (read.mate_chromosome(), read.mate_alignment_start(), read.mate_alignment_end(), read.mate_fragment_end())
        } else {
            (Some(NO_CHROMOSOME_NAME), NO_POSITION, NO_POSITION, NO_POSITION)
        };
        let indel = read.indel_coords();
        crate::assembly::mem_stats::SUPPORT_READS_LIVE.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        A::new(SupportRead {
            type_,
            id: read.id().into(),
            chromosome: read.chromosome(),
            alignment_start: read.alignment_start(),
            alignment_end: read.alignment_end(),
            flags: read.get_flags(),
            cigar: read.cigar_string().into(),
            unclipped_start: read.unclipped_start(),
            unclipped_end: read.unclipped_end(),
            mate_chromosome,
            mate_alignment_start,
            mate_alignment_end,
            mate_fragment_end,
            is_reference: read.is_reference(),
            sample_index: read.sample_index(),
            is_discordant: is_discordant_fragment(read),
            supplementary_data: read.supplementary_data().cloned().map(Box::new),
            base_length: read.bases_length(),
            insert_size: inferred_insert_size_abs(read.bam_record()),
            trim_count_start: read.trim_count_start(),
            trim_count_end: read.trim_count_end(),
            map_qual: read.mapping_quality(),
            has_indel: indel.is_some(),
            indel_coords: indel.filter(|c| c.length >= MIN_INDEL_LENGTH),
            has_line_tail: read.has_line_tail_start() || read.has_line_tail_end(),
            ext_base_matches: matches,
            ext_base_mismatches: mismatches,
            ext_base_overlap: 0,
            ref_base_mismatches: NO_REF_BASE_MISMATCHES,
            mismatch_info: Box::default(),
            medium_qual_count: if is_sbx() && type_ == SupportType::Junction { set_sbx_medium_qual_count(read) } else { 0 },
            junction_read_start_distance: junct_read_start_distance,
            full_assembly_index: -1,
            full_assembly_orientation: Orientation::Forward,
            inferred_fragment_length: -1,
            breakend_type: None,
            read: Some(read.clone()),
        })
    }

    /// estimate of the bytes held, for the memory diagnostics
    pub fn approx_bytes(&self) -> usize {
        std::mem::size_of::<SupportRead>() + 16 + self.id.len() + self.cigar.len() + self.mismatch_info.len()
            + self.supplementary_data.as_ref().map_or(0, |s| std::mem::size_of::<SupplementaryReadData>() + s.cigar.capacity())
    }
    pub fn type_(&self) -> SupportType {
        self.type_
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn chromosome(&self) -> &'static str {
        self.chromosome
    }
    pub fn alignment_start(&self) -> i32 {
        self.alignment_start
    }
    pub fn alignment_end(&self) -> i32 {
        self.alignment_end
    }
    pub fn unclipped_start(&self) -> i32 {
        self.unclipped_start
    }
    pub fn unclipped_end(&self) -> i32 {
        self.unclipped_end
    }
    pub fn untrimmed_start(&self) -> i32 {
        self.unclipped_start - self.trim_count_start
    }
    pub fn untrimmed_end(&self) -> i32 {
        self.unclipped_end + self.trim_count_end
    }
    pub fn is_left_clipped(&self) -> bool {
        self.unclipped_start < self.alignment_start
    }
    pub fn is_right_clipped(&self) -> bool {
        self.unclipped_end > self.alignment_end
    }
    pub fn left_clip_length(&self) -> i32 {
        (self.alignment_start - self.unclipped_start).max(0)
    }
    pub fn right_clip_length(&self) -> i32 {
        (self.unclipped_end - self.alignment_end).max(0)
    }
    /// null in Java when the mate is unmapped
    pub fn mate_chromosome(&self) -> Option<&'static str> {
        self.mate_chromosome
    }
    pub fn mate_alignment_start(&self) -> i32 {
        self.mate_alignment_start
    }
    pub fn mate_alignment_end(&self) -> i32 {
        self.mate_alignment_end
    }
    pub fn mate_fragment_end(&self) -> i32 {
        self.mate_fragment_end
    }
    pub fn base_length(&self) -> i32 {
        self.base_length
    }
    pub fn insert_size(&self) -> i32 {
        self.insert_size
    }
    pub fn trim_count(&self) -> i32 {
        self.trim_count_start + self.trim_count_end
    }
    pub fn trim_count_start(&self) -> i32 {
        self.trim_count_start
    }
    pub fn trim_count_end(&self) -> i32 {
        self.trim_count_end
    }
    pub fn has_indel(&self) -> bool {
        self.has_indel
    }
    pub fn indel_coords(&self) -> Option<&A<IndelCoords>> {
        self.indel_coords.as_ref()
    }
    pub fn cigar(&self) -> &str {
        &self.cigar
    }
    pub fn orientation(&self) -> Orientation {
        if self.is_flag_set(FLAG_READ_REVERSE_STRAND) { Orientation::Reverse } else { Orientation::Forward }
    }
    pub fn mate_orientation(&self) -> Orientation {
        if self.is_flag_set(FLAG_MATE_REVERSE_STRAND) { Orientation::Reverse } else { Orientation::Forward }
    }
    pub fn supplementary_data(&self) -> Option<&SupplementaryReadData> {
        self.supplementary_data.as_deref()
    }
    pub fn map_qual(&self) -> i32 {
        self.map_qual
    }
    pub fn has_line_tail(&self) -> bool {
        self.has_line_tail
    }
    pub fn is_reference(&self) -> bool {
        self.is_reference
    }
    pub fn sample_index(&self) -> i32 {
        self.sample_index
    }
    pub fn flags(&self) -> i32 {
        self.flags
    }
    pub fn is_supplementary(&self) -> bool {
        self.is_flag_set(FLAG_SUPPLEMENTARY)
    }
    pub fn is_paired_read(&self) -> bool {
        self.is_flag_set(FLAG_READ_PAIRED)
    }
    pub fn first_in_pair(&self) -> bool {
        self.is_flag_set(FLAG_FIRST_OF_PAIR)
    }
    pub fn is_unmapped(&self) -> bool {
        self.is_flag_set(FLAG_READ_UNMAPPED)
    }
    pub fn is_mate_unmapped(&self) -> bool {
        self.is_flag_set(FLAG_MATE_UNMAPPED)
    }
    pub fn is_mate_mapped(&self) -> bool {
        self.is_flag_set(FLAG_READ_PAIRED) && !self.is_flag_set(FLAG_MATE_UNMAPPED)
    }
    pub fn five_prime_fragment_position(&self) -> i32 {
        if self.orientation().is_forward() { self.unclipped_start } else { self.unclipped_end }
    }
    pub fn is_discordant(&self) -> bool {
        self.is_discordant
    }
    #[inline]
    pub fn is_flag_set(&self, flag: i32) -> bool {
        self.flags & flag != 0
    }
    pub fn extension_base_mismatches(&self) -> i32 {
        self.ext_base_mismatches
    }
    pub fn extension_base_matches(&self) -> i32 {
        self.ext_base_matches
    }
    pub fn extension_base_overlap(&self) -> i32 {
        self.ext_base_overlap
    }
    pub fn medium_qual_count(&self) -> i32 {
        self.medium_qual_count
    }
    pub fn reference_mismatches(&self) -> i32 {
        if self.ref_base_mismatches == NO_REF_BASE_MISMATCHES { -1 } else { self.ref_base_mismatches }
    }
    pub fn has_reference_mismatches(&self) -> bool {
        self.ref_base_mismatches != NO_REF_BASE_MISMATCHES
    }
    pub fn mismatch_info(&self) -> &str {
        &self.mismatch_info
    }
    pub fn set_extension_match_info(&mut self, read_info: &ReadParseState) {
        self.mismatch_info = read_info.mismatch_info().into();
        self.ext_base_overlap = read_info.evaluated_base_count() - 1;
    }
    pub fn set_reference_mismatches(&mut self, mismatches: i32) {
        self.ref_base_mismatches = mismatches;
    }
    pub fn cached_read(&self) -> Option<&ReadP> {
        self.read.as_ref()
    }
    pub fn clear_cached_read(&mut self) {
        self.read = None;
    }
    pub fn junction_read_start_distance(&self) -> i32 {
        self.junction_read_start_distance
    }
    pub fn junction_extension_length(&self, junction_orientation: Orientation) -> i32 {
        if self.type_ == SupportType::Junction || self.type_ == SupportType::Indel {
            return if junction_orientation.is_forward() { self.base_length - self.junction_read_start_distance - 1 } else { self.junction_read_start_distance };
        }
        0
    }
    pub fn set_full_assembly_info(&mut self, assembly_index: i32, orientation: Orientation) {
        self.full_assembly_index = assembly_index;
        self.full_assembly_orientation = orientation;
    }
    pub fn full_assembly_index_start(&self) -> i32 {
        self.full_assembly_index
    }
    pub fn full_assembly_index_end(&self) -> i32 {
        self.full_assembly_index + self.base_length - 1
    }
    pub fn full_assembly_orientation(&self) -> Orientation {
        self.full_assembly_orientation
    }
    pub fn inferred_fragment_length(&self) -> i32 {
        self.inferred_fragment_length
    }
    pub fn set_inferred_fragment_length(&mut self, length: i32) {
        self.inferred_fragment_length = length;
    }
    pub fn breakend_support_type(&self) -> Option<SupportType> {
        self.breakend_type
    }
    pub fn set_breakend_support_type(&mut self, breakend_type: Option<SupportType>) {
        self.breakend_type = breakend_type;
    }
    pub fn matches_fragment(&self, other: &SupportRead, allow_read_match: bool) -> bool {
        if self.id != other.id {
            return false;
        }
        allow_read_match || self.flags != other.flags
    }
    pub fn matches_fragment_read(&self, other: &Read, allow_read_match: bool) -> bool {
        if &*self.id != other.id() {
            return false;
        }
        allow_read_match || self.flags != other.get_flags()
    }
    pub fn has_fragment_other_read(support: &[SupP], read: &Read) -> bool {
        support.iter().any(|x| x.matches_fragment_read(read, false))
    }
}

impl std::fmt::Display for SupportRead {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "type({}) read({} {}:{}-{} {} {}) index(juncDist={} asm={}:{}) hqMatch({}) mismatch(junc={} ref={})",
            self.type_.name(),
            self.id,
            self.chromosome,
            self.alignment_start,
            self.alignment_end,
            self.cigar,
            self.orientation().as_byte(),
            self.junction_read_start_distance,
            self.full_assembly_index,
            self.full_assembly_orientation.as_byte(),
            self.ext_base_matches,
            self.ext_base_mismatches,
            self.reference_mismatches()
        )
    }
}
