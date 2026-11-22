//! Breakend, BreakendSegment and BreakendSupport: a breakend called from the alignments of an assembly sequence.

use std::cmp::Ordering;
use std::collections::HashSet;

use crate::common::compat;
use crate::assembly::alignment::align_data::AlignP;
use crate::assembly::alignment::alternative_alignment::AlternativeAlignment;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::homology_data::HomologyData;
use crate::assembly::common_utils::{compare_junctions, position_within};
use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::P;
use crate::assembly::types::assembly_stats::java_round;
use crate::assembly::types::support_read::SupportRead;
use crate::assembly::types::{form_sv_type, StructuralVariantType, SupportType};

pub type BreakendP = P<Breakend>;

pub struct BreakendSegment {
    pub assembly_alignment_id: i32,
    pub sequence_index: i32,
    pub index: i32,
    pub alignment: AlignP,
    indel_sequence_indices: Option<[i32; 2]>,
}

impl BreakendSegment {
    pub fn new(assembly_alignment_id: i32, sequence_index: i32, index: i32, alignment: AlignP) -> Self {
        BreakendSegment { assembly_alignment_id, sequence_index, index, alignment, indel_sequence_indices: None }
    }
    pub fn new_indel(assembly_alignment_id: i32, sequence_index: i32, index: i32, alignment: AlignP, indel_indices: Option<[i32; 2]>) -> Self {
        BreakendSegment { assembly_alignment_id, sequence_index, index, alignment, indel_sequence_indices: indel_indices }
    }
    pub fn indel_sequence_indices(&self) -> Option<[i32; 2]> {
        self.indel_sequence_indices
    }
    pub fn unique_id(&self) -> String {
        format!("{}_{}", self.assembly_alignment_id, self.index)
    }
}

#[derive(Default, Clone, Debug)]
pub struct BreakendSupport {
    pub split_fragments: i32,
    pub discordant_fragments: i32,
    pub forward_reads: i32,
    pub reverse_reads: i32,
}

impl BreakendSupport {
    pub fn total_support(&self) -> i32 {
        self.split_fragments + self.discordant_fragments
    }
    pub fn strand_bias(&self) -> f64 {
        let total_reads = self.forward_reads + self.reverse_reads;
        if total_reads > 0 { self.forward_reads as f64 / total_reads as f64 } else { 0.0 }
    }
}

pub struct Breakend {
    pub chromosome: String,
    pub position: i32,
    pub orient: Orientation,
    pub inserted_bases: String,
    pub homology: HomologyData,
    id: i32,
    assembly: AsmAlignP,
    other_breakend: Option<BreakendP>,
    facing_breakends: Vec<BreakendP>,
    alternative_alignments: Option<Vec<AlternativeAlignment>>,
    segments: Vec<BreakendSegment>,
    breakend_support: Vec<BreakendSupport>,
    fragment_length_total: i32,
    fragment_length_count: i32,
    incomplete_fragment_count: i32,
    non_primary_assembly_fragment_count: i32,
    fragment_positions: Option<HashSet<i32>>,
    five_prime_range_start: i32,
    five_prime_range_end: i32,
    three_prime_range_start: i32,
    three_prime_range_end: i32,
    max_local_repeat: i32,
    saga_inferred: bool,
}

impl Breakend {
    pub fn new(assembly: AsmAlignP, chromosome: &str, position: i32, orientation: Orientation, inserted_bases: &str, homology: Option<HomologyData>) -> BreakendP {
        P::new(Breakend {
            chromosome: chromosome.to_string(),
            position,
            orient: orientation,
            inserted_bases: inserted_bases.to_string(),
            homology: homology.unwrap_or_else(HomologyData::no_homology),
            id: -1,
            assembly,
            other_breakend: None,
            facing_breakends: Vec::new(),
            alternative_alignments: None,
            segments: Vec::new(),
            breakend_support: Vec::new(),
            fragment_length_total: 0,
            fragment_length_count: 0,
            incomplete_fragment_count: 0,
            non_primary_assembly_fragment_count: 0,
            fragment_positions: None,
            five_prime_range_start: 0,
            five_prime_range_end: 0,
            three_prime_range_start: 0,
            three_prime_range_end: 0,
            max_local_repeat: 0,
            saga_inferred: false,
        })
    }
    pub fn id(&self) -> i32 {
        self.id
    }
    pub fn set_id(&mut self, id: i32) {
        self.id = id;
    }
    pub fn assembly(&self) -> AsmAlignP {
        self.assembly
    }
    pub fn other_breakend(&self) -> Option<BreakendP> {
        self.other_breakend
    }
    pub fn set_other_breakend(&mut self, breakend: Option<BreakendP>) {
        self.other_breakend = breakend;
    }
    pub fn facing_breakends(&self) -> &Vec<BreakendP> {
        &self.facing_breakends
    }
    pub fn add_facing_breakend(&mut self, breakend: BreakendP) {
        self.facing_breakends.push(breakend);
    }
    pub fn segments(&self) -> &Vec<BreakendSegment> {
        &self.segments
    }
    pub fn add_segment(&mut self, segment: BreakendSegment) {
        self.segments.push(segment);
    }
    pub fn add_fragment_position(&mut self, position: i32) {
        self.fragment_positions.get_or_insert_with(HashSet::new).insert(position);
    }
    pub fn unique_fragment_position_count(&self) -> i32 {
        self.fragment_positions.as_ref().map(|s| s.len() as i32).unwrap_or(0)
    }
    pub fn add_orientation_positions(&mut self, read_orientation: Orientation, read_start: i32, read_end: i32) {
        let five_prime_position = if read_orientation.is_forward() { read_start } else { read_end };
        let three_prime_position = if read_orientation.is_forward() { read_end } else { read_start };
        if self.five_prime_range_start == 0 {
            self.five_prime_range_start = five_prime_position;
        } else {
            self.five_prime_range_start = self.five_prime_range_start.min(five_prime_position);
        }
        self.five_prime_range_end = self.five_prime_range_end.max(five_prime_position);
        if self.three_prime_range_start == 0 {
            self.three_prime_range_start = three_prime_position;
        } else {
            self.three_prime_range_start = self.three_prime_range_start.min(three_prime_position);
        }
        self.three_prime_range_end = self.three_prime_range_end.max(three_prime_position);
    }
    pub fn read_orientation_range(&self) -> [i32; 2] {
        [
            if self.five_prime_range_start > 0 { self.five_prime_range_end - self.five_prime_range_start } else { -1 },
            if self.three_prime_range_start > 0 { self.three_prime_range_end - self.three_prime_range_start } else { -1 },
        ]
    }
    pub fn set_max_local_repeat(&mut self, max_local_repeat: i32) {
        self.max_local_repeat = max_local_repeat;
    }
    pub fn max_local_repeat(&self) -> i32 {
        self.max_local_repeat
    }
    pub fn alternative_alignments(&self) -> &[AlternativeAlignment] {
        self.alternative_alignments.as_deref().unwrap_or(&[])
    }
    pub fn has_alternative_alignments_set(&self) -> bool {
        self.alternative_alignments.is_some()
    }
    pub fn set_alternative_alignments(&mut self, alt_alignments: Option<Vec<AlternativeAlignment>>) {
        self.alternative_alignments = alt_alignments;
    }
    pub fn low_qual_alt_alignments(&self) -> Vec<AlternativeAlignment> {
        if !self.segments.is_empty() && self.segments[0].alignment.has_low_map_qual_alignment() {
            return self.segments[0].alignment.unselected_alt_alignments().to_vec();
        }
        Vec::new()
    }
    pub fn sample_support(&self) -> &Vec<BreakendSupport> {
        &self.breakend_support
    }
    pub fn sample_support_mut(&mut self) -> &mut Vec<BreakendSupport> {
        &mut self.breakend_support
    }
    pub fn update_breakend_support(&mut self, sample_index: i32, is_split_fragment: bool) {
        if (sample_index as usize) < self.breakend_support.len() {
            let s = &mut self.breakend_support[sample_index as usize];
            if is_split_fragment {
                s.split_fragments += 1;
            } else {
                s.discordant_fragments += 1;
            }
        }
    }
    pub fn update_breakend_read_strand(&mut self, sample_index: i32, read: &SupportRead) {
        if read.type_() != SupportType::Junction || !crate::assembly::alignment::alignment_fragments::read_alignment_spans_junction(self, read) {
            return;
        }
        if (sample_index as usize) < self.breakend_support.len() {
            let s = &mut self.breakend_support[sample_index as usize];
            if read.orientation().is_forward() {
                s.forward_reads += 1;
            } else {
                s.reverse_reads += 1;
            }
        }
    }
    pub fn average_fragment_length(&self) -> i32 {
        if self.fragment_length_count > 0 { java_round(self.fragment_length_total as f64 / self.fragment_length_count as f64) as i32 } else { 0 }
    }
    pub fn valid_fragment_length_percent(&self) -> f64 {
        if self.fragment_length_count > 0 { self.fragment_length_count as f64 / (self.fragment_length_count + self.incomplete_fragment_count) as f64 } else { 0.0 }
    }
    pub fn incomplete_fragment_count(&self) -> i32 {
        self.incomplete_fragment_count
    }
    pub fn add_inferred_fragment_length(&mut self, length: i32, is_complete: bool) {
        if is_complete && length > 0 {
            self.fragment_length_count += 1;
            self.fragment_length_total += length;
        } else {
            self.incomplete_fragment_count += 1;
        }
    }
    pub fn add_non_primary_assembly_fragment_count(&mut self) {
        self.non_primary_assembly_fragment_count += 1;
    }
    pub fn non_primary_assembly_fragment_count(&self) -> i32 {
        self.non_primary_assembly_fragment_count
    }
    pub fn is_single(&self) -> bool {
        self.other_breakend.is_none()
    }
    pub fn sv_type(&self) -> StructuralVariantType {
        match self.other_breakend {
            None => StructuralVariantType::Sgl,
            Some(o) => form_sv_type(&self.chromosome, &o.chromosome, self.position, o.position, self.orient, o.orient, !self.inserted_bases.is_empty()),
        }
    }
    pub fn sv_length(&self) -> i32 {
        match self.other_breakend {
            Some(o) if o.chromosome == self.chromosome => {
                let pos_length = (self.position - o.position).abs();
                if self.sv_type() == StructuralVariantType::Dup { pos_length + 1 } else { pos_length }
            }
            _ => 0,
        }
    }
    pub fn is_short_local_del_dup_ins(&self) -> bool {
        matches!(self.sv_type(), StructuralVariantType::Del | StructuralVariantType::Dup | StructuralVariantType::Ins) && self.sv_length() <= 1000
    }
    pub fn min_position(&self) -> i32 {
        self.position + self.homology.exact_start
    }
    pub fn max_position(&self) -> i32 {
        self.position + self.homology.exact_end
    }
    pub fn calc_sv_qual(&self) -> i32 {
        let mut final_qual = self.calc_qual();
        if let Some(o) = self.other_breakend {
            final_qual += o.calc_qual();
        }
        let total_support: i32 = self.breakend_support.iter().map(|x| x.total_support()).sum();
        let support_factor = total_support as f64 / (total_support as f64 + QUAL_CALC_FRAG_SUPPORT_FACTOR);
        java_round(final_qual as f64 * support_factor) as i32
    }
    pub fn calc_qual(&self) -> i32 {
        let mut max_segment_qual = 0;
        for segment in &self.segments {
            let mut segment_qual = java_round(segment.alignment.modified_map_qual()) as i32;
            if segment.alignment.has_low_map_qual_short_sv_link() {
                segment_qual += ALIGNMENT_LOW_MOD_MQ_QUAL_BOOST;
            }
            max_segment_qual = segment_qual.max(max_segment_qual);
        }
        max_segment_qual
    }
    pub fn matches_coordinates(&self, chromosome: &str, position: i32, orientation: Orientation) -> bool {
        if self.chromosome != chromosome || self.orient != orientation {
            return false;
        }
        if !self.homology.is_no_homology() {
            position_within(position, self.position + self.homology.exact_start, self.position + self.homology.exact_end)
        } else {
            self.position == position
        }
    }
    pub fn is_saga_inferred(&self) -> bool {
        self.saga_inferred
    }
    pub fn set_saga_inferred(&mut self) {
        self.saga_inferred = true;
    }

    /// Breakend.compareTo. As in Java, the last comparison passes the same chromosome twice, so mates on
    /// different chromosomes are compared by position alone.
    pub fn compare_to(&self, other: &Breakend) -> Ordering {
        let breakend_compare = compare_junctions(&self.chromosome, &other.chromosome, self.position, other.position, self.orient, other.orient);
        if breakend_compare != Ordering::Equal {
            return breakend_compare;
        }
        let (t1, t2) = (self.sv_type(), other.sv_type());
        if t1 == StructuralVariantType::Sgl || t2 == StructuralVariantType::Sgl {
            if t1 == t2 {
                return self.inserted_bases.len().cmp(&other.inserted_bases.len());
            }
            return if t1 == StructuralVariantType::Sgl { Ordering::Less } else { Ordering::Greater };
        }
        let (o1, o2) = (self.other_breakend.unwrap(), other.other_breakend.unwrap());
        // Java: passes this breakend's mate chromosome twice, so the mates' chromosomes are never compared
        let o2_chromosome = if compat::java("A18") { &o1.chromosome } else { &o2.chromosome };
        compare_junctions(&o1.chromosome, o2_chromosome, o1.position, o2.position, o1.orient, o2.orient)
    }
}

impl std::fmt::Display for Breakend {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}: {}:{}:{} {} hom({}) otherId({}) segs({}) alts({})",
            self.id,
            self.chromosome,
            self.position,
            self.orient.as_byte(),
            self.sv_type().name(),
            if !self.homology.is_no_homology() { self.homology.to_string() } else { String::new() },
            self.other_breakend.map(|o| o.id()).unwrap_or(-1),
            self.segments.len(),
            self.alternative_alignments.as_ref().map(|a| a.len()).unwrap_or(0)
        )
    }
}
