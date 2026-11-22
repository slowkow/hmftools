//! esvee.common: IndelCoords and CommonUtils.

use std::cmp::Ordering;

use crate::assembly::constants::*;
use crate::assembly::globals::{is_illumina, is_sbx};
use crate::assembly::orientation::Orientation;
use crate::assembly::sam::*;
use crate::assembly::types::junction_assembly::JunctionAssembly;
use crate::common::regions::chromosome_rank;

pub struct IndelCoords {
    pub pos_start: i32,
    pub pos_end: i32,
    pub length: i32,
    inserted_bases: Option<String>,
}

impl IndelCoords {
    pub fn new(pos_start: i32, pos_end: i32, length: i32) -> Self {
        IndelCoords { pos_start, pos_end, length, inserted_bases: None }
    }
    pub fn is_insert(&self) -> bool {
        self.pos_end == self.pos_start + 1
    }
    pub fn is_delete(&self) -> bool {
        !self.is_insert()
    }
    pub fn inserted_bases(&self) -> &str {
        self.inserted_bases.as_deref().unwrap_or("")
    }
    pub fn set_inserted_bases(&mut self, bases: String) {
        self.inserted_bases = Some(bases);
    }
    pub fn matches_junction(&self, position: i32, orientation: Orientation) -> bool {
        if orientation.is_forward() { self.pos_start == position } else { self.pos_end == position }
    }
    pub fn matches(&self, other: &IndelCoords) -> bool {
        self.pos_start == other.pos_start && self.pos_end == other.pos_end && self.length == other.length
    }

    pub fn find_indel_coords(read_start: i32, cigar_elements: &[CigarElement], min_indel_length: i32) -> Option<IndelCoords> {
        let mut max_indel_length = 0;
        let mut indel_start_pos = read_start - 1;
        let mut indel_end_pos = 0;
        let mut ref_position = read_start;
        for element in cigar_elements {
            if element.op.is_indel() && element.len > max_indel_length {
                max_indel_length = element.len;
                indel_start_pos = ref_position - 1;
                indel_end_pos = if element.op == CigarOp::D { indel_start_pos + element.len + 1 } else { indel_start_pos + 1 };
            }
            if element.op.consumes_reference_bases() {
                ref_position += element.len;
            }
        }
        if indel_end_pos <= indel_start_pos || max_indel_length < min_indel_length {
            return None;
        }
        Some(IndelCoords::new(indel_start_pos, indel_end_pos, max_indel_length))
    }

    pub fn find_matching_indel_coords(read_start: i32, cigar_elements: &[CigarElement], specific_indel: &CigarElement) -> Option<IndelCoords> {
        let mut indel_start_pos = -1;
        let mut indel_end_pos = 0;
        let mut matched_indel_length = 0;
        let mut ref_position = read_start;
        for element in cigar_elements {
            if element.op == specific_indel.op && (specific_indel.len - element.len).abs() <= 1 {
                indel_start_pos = ref_position - 1;
                indel_end_pos = if element.op == CigarOp::D { indel_start_pos + element.len + 1 } else { indel_start_pos + 1 };
                matched_indel_length = element.len;
                break;
            }
            if element.op.consumes_reference_bases() {
                ref_position += element.len;
            }
        }
        if matched_indel_length > 0 { Some(IndelCoords::new(indel_start_pos, indel_end_pos, matched_indel_length)) } else { None }
    }
}

impl std::fmt::Display for IndelCoords {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({} - {}) len({})", if self.is_delete() { "delete" } else { "insert" }, self.pos_start, self.pos_end, self.length)
    }
}

// BaseQualAdjustment and the per-technology high / medium quality tests
#[inline]
pub fn above_min_qual(qual: u8) -> bool {
    qual > LOW_BASE_QUAL_THRESHOLD
}
#[inline]
pub fn below_min_qual(qual: u8) -> bool {
    qual <= LOW_BASE_QUAL_THRESHOLD
}
#[inline]
pub fn is_uncertain_base_qual(qual: u8) -> bool {
    qual <= BASE_QUAL_MINIMUM
}
pub fn is_high_base_qual(qual: u8) -> bool {
    if is_illumina() {
        crate::assembly::seq_tech_utils::illumina_is_high_base_qual(qual)
    } else if is_sbx() {
        crate::assembly::seq_tech_utils::sbx_is_high_base_qual(qual)
    } else {
        crate::assembly::seq_tech_utils::ultima_is_high_base_qual(qual)
    }
}
pub fn is_medium_base_qual(qual: u8) -> bool {
    if is_sbx() { crate::assembly::seq_tech_utils::sbx_is_medium_base_qual(qual) } else { false }
}

pub fn is_higher_base_qual_category(qual1: u8, qual2: u8) -> bool {
    if qual1 == qual2 {
        return false;
    }
    let low1 = below_min_qual(qual1);
    let low2 = below_min_qual(qual2);
    if low1 != low2 {
        return low2;
    }
    let high1 = is_high_base_qual(qual1);
    let high2 = is_high_base_qual(qual2);
    if high1 != high2 {
        return high1;
    }
    qual1 > qual2
}

pub fn is_discordant_fragment(read: &SamRecord, max_concordant_fragment_length: i32, supp_data: Option<&SupplementaryReadData>) -> bool {
    if read.read_unmapped() || !read.read_paired() || read.mate_unmapped() {
        return false;
    }
    if read.supplementary() && supp_data.is_some() {
        if supp_data.unwrap().chromosome != read.mate_reference_name {
            return true;
        }
    } else if read.reference_name != read.mate_reference_name {
        return true;
    }
    if read.read_negative_strand() == read.mate_negative_strand() {
        return true;
    }
    let fragment_size = read.inferred_insert_size.abs();
    if fragment_size == 0 || (max_concordant_fragment_length > 0 && fragment_size >= max_concordant_fragment_length) {
        return true;
    }
    is_duplication_fragment(read, fragment_size)
}

pub fn is_duplication_fragment(read: &SamRecord, fragment_size: i32) -> bool {
    if fragment_size > 2 * read.read_bases.len() as i32 {
        if !read.read_negative_strand() {
            if read.mate_alignment_start < read.alignment_start {
                return true;
            }
        } else if read.alignment_start < read.mate_alignment_start {
            return true;
        }
    }
    false
}

pub fn is_line_insert_pair(assembly1: &JunctionAssembly, assembly2: &JunctionAssembly) -> bool {
    let (j1, j2) = (assembly1.junction(), assembly2.junction());
    if j1.chromosome != j2.chromosome {
        return false;
    }
    within_line_proximity(j1.position, j2.position, j1.orient, j2.orient)
}

pub fn within_line_proximity(pos1: i32, pos2: i32, orient1: Orientation, orient2: Orientation) -> bool {
    if orient1 == orient2 {
        return false;
    }
    let pos_diff = (pos2 - pos1).abs();
    if (pos1 < pos2) == orient1.is_forward() { pos_diff <= LINE_INDEL_MAX_GAP } else { pos_diff <= LINE_INDEL_MAX_OVERLAP }
}

/// CommonUtils.compareJunctions. Two different contigs of equal rank each compare as greater than the other,
/// as in Java.
pub fn compare_junctions(chr1: &str, chr2: &str, pos1: i32, pos2: i32, orient1: Orientation, orient2: Orientation) -> Ordering {
    if chr1 != chr2 {
        return if chromosome_rank(chr1) < chromosome_rank(chr2) { Ordering::Less } else { Ordering::Greater };
    }
    if pos1 == pos2 {
        if orient1 == orient2 {
            return Ordering::Equal;
        }
        return if orient1.is_forward() { Ordering::Less } else { Ordering::Greater };
    }
    if pos1 < pos2 { Ordering::Less } else { Ordering::Greater }
}

pub fn is_valid_dna_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

/// CommonUtils.checkStandardNucleotides
pub fn check_standard_nucleotides(read: &mut SamRecord) {
    for i in 0..read.read_bases.len() {
        if !is_valid_dna_base(read.read_bases[i]) {
            read.read_bases[i] = b'N';
            if i < read.base_qualities.len() {
                read.base_qualities[i] = BASE_QUAL_MINIMUM;
            }
        }
    }
}

pub fn create_byte_array(length: usize, value: u8) -> Vec<u8> {
    vec![value; length]
}

/// BaseRegion.positionsOverlap
#[inline]
pub fn positions_overlap(start1: i32, end1: i32, start2: i32, end2: i32) -> bool {
    !(start1 > end2 || end1 < start2)
}
/// BaseRegion.positionWithin
#[inline]
pub fn position_within(position: i32, start: i32, end: i32) -> bool {
    position >= start && position <= end
}
/// BaseRegion.positionsWithin
#[inline]
pub fn positions_within(inner_start: i32, inner_end: i32, outer_start: i32, outer_end: i32) -> bool {
    inner_start <= inner_end && inner_start >= outer_start && inner_end <= outer_end
}

/// HumanChromosome.contains: 1-22, X or Y, with or without the chr prefix
pub fn human_chromosome_contains(chromosome: &str) -> bool {
    let trimmed = chromosome.strip_prefix("chr").unwrap_or(chromosome);
    if !trimmed.is_empty() && trimmed.bytes().all(|b| b.is_ascii_digit()) {
        return matches!(trimmed.parse::<i64>(), Ok(v) if (1..=22).contains(&v));
    }
    trimmed == "X" || trimmed == "Y"
}
