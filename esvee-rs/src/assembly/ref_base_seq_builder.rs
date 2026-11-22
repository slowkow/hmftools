//! RefBaseSeqBuilder: the consensus of the reference-side bases of a junction assembly, built from the CIGARs
//! of its junction reads walking away from the junction.

use crate::assembly::assembly_utils::{self, DNA_BASE_BYTES, DNA_BASE_COUNT, NO_BASE};
use crate::assembly::common_utils::{above_min_qual, below_min_qual};
use crate::assembly::constants::*;
use crate::assembly::read::read::Read;
use crate::assembly::read::read_utils::get_read_index_at_junction;
use crate::assembly::read_parse_state::{ReadParseState, RpsP};
use crate::assembly::sam::{cigar_elements_to_str, CigarElement, CigarOp};
use crate::assembly::sequence_builder::exceeds_read_mismatches;
use crate::assembly::sequence_diff::SequenceDiffType;
use crate::assembly::types::junction_assembly::AsmP;

pub struct RefBaseSeqBuilder {
    assembly: AsmP,
    reads: Vec<RpsP>,
    is_forward_junction: bool,
    junction_position: i32,
    bases: Vec<u8>,
    base_quals: Vec<u8>,
    cigar_elements: Vec<CigarElement>,
    ref_base_position: i32,
}

impl RefBaseSeqBuilder {
    pub fn new(assembly: AsmP) -> RefBaseSeqBuilder {
        let junction = assembly.junction();
        let is_forward_junction = junction.is_forward();
        let junction_position = junction.position;
        let mut max_read_base_length = 0;
        let mut reads: Vec<RpsP> = Vec::with_capacity(assembly.support_count());
        for support in assembly.support() {
            let read = support.cached_read().expect("support read without cached read");
            let read_junction_index = get_read_index_at_junction(read, &junction, true);
            let mut has_valid_junction_overlap;
            let indel_coords = read.indel_coords();
            if indel_coords.is_some_and(|c| c.matches_junction(junction_position, junction.orient)) {
                has_valid_junction_overlap = true;
            } else if is_forward_junction {
                has_valid_junction_overlap = assembly.discordant_only() || read.is_right_clipped();
                has_valid_junction_overlap &= read.unclipped_end() - junction_position >= ASSEMBLY_REF_READ_MIN_SOFT_CLIP;
            } else {
                has_valid_junction_overlap = assembly.discordant_only() || read.is_left_clipped();
                has_valid_junction_overlap &= junction_position - read.unclipped_start() >= ASSEMBLY_REF_READ_MIN_SOFT_CLIP;
            }
            let read_ref_base_length = Self::read_ref_base_length(read, read_junction_index, is_forward_junction);
            let read_state = ReadParseState::new_cigar(!is_forward_junction, read, read_junction_index, true);
            if has_valid_junction_overlap {
                max_read_base_length = max_read_base_length.max(read_ref_base_length);
            } else {
                read_state.m().mark_invalid();
            }
            reads.push(read_state);
        }
        let base_length = (max_read_base_length + 1) as usize;
        let mut b = RefBaseSeqBuilder {
            assembly,
            reads,
            is_forward_junction,
            junction_position,
            bases: vec![0u8; base_length],
            base_quals: vec![0u8; base_length],
            cigar_elements: Vec::new(),
            ref_base_position: junction_position,
        };
        b.build_sequence();
        b.trim_final_sequence();
        b
    }

    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
    pub fn base_qualities(&self) -> &[u8] {
        &self.base_quals
    }
    pub fn ref_base_length(&self) -> i32 {
        self.bases.len() as i32
    }
    pub fn ref_base_position(&self) -> i32 {
        self.ref_base_position
    }
    pub fn cigar_elements(&self) -> &[CigarElement] {
        &self.cigar_elements
    }
    pub fn cigar_str(&self) -> String {
        cigar_elements_to_str(&self.cigar_elements)
    }
    pub fn reads(&self) -> &Vec<RpsP> {
        &self.reads
    }
    pub fn ref_base_sequence(&self) -> String {
        self.bases.iter().map(|b| *b as char).collect()
    }

    pub fn read_ref_base_length(read: &Read, read_junction_index: i32, is_forward_junction: bool) -> i32 {
        let insert_length: i32 = read.cigar_elements().iter().filter(|x| x.op == CigarOp::I).map(|x| x.len).sum();
        let aligned_length = if is_forward_junction {
            (read_junction_index - read.left_clip_length()).max(0)
        } else {
            (read.bases_length() - read_junction_index - 1 - read.right_clip_length()).max(0)
        };
        aligned_length + insert_length
    }

    fn build_sequence(&mut self) {
        let len = self.bases.len() as i32;
        let mut current_index: i32 = if self.is_forward_junction { len - 1 } else { 0 };
        let mut ref_position = self.junction_position;
        self.bases[current_index as usize] = self.assembly.bases()[self.assembly.junction_index() as usize];
        self.base_quals[current_index as usize] = self.assembly.base_quals()[self.assembly.junction_index() as usize];
        let mut current_element_type = CigarOp::M;
        let mut current_element_length = 0;
        let mut active_reads: Vec<RpsP> = self.reads.iter().filter(|x| x.is_valid()).cloned().collect();

        while !active_reads.is_empty() {
            if current_element_type == CigarOp::M || current_element_type == CigarOp::I {
                current_index += if self.is_forward_junction { -1 } else { 1 };
            }
            current_element_length += 1;
            if current_index < 0 || current_index >= len {
                break;
            }
            let mut consensus_base: u8 = 0;
            let mut consensus_max_qual: i32 = 0;
            let mut consensus_qual_total: i32 = 0;
            let mut consensus_read_count = 0;
            let mut has_arrays = false;
            let mut read_counts = [0i32; DNA_BASE_COUNT];
            let mut total_quals = [0i32; DNA_BASE_COUNT];
            let mut max_quals = [0i32; DNA_BASE_COUNT];

            let prior_active_reads = active_reads.len();
            progress_read_state(&mut active_reads, current_element_type);
            if active_reads.is_empty() {
                break;
            }
            let soft_clipped_reads = self.reads.iter().filter(|x| x.exhausted() && x.operator() == Some(CigarOp::S)).count();
            if (active_reads.len() as f64 / prior_active_reads as f64) < 0.1 && soft_clipped_reads as f64 / prior_active_reads as f64 > 0.9 {
                break;
            }
            let next_element_type = find_next_operator(&active_reads);
            if next_element_type != current_element_type {
                self.cigar_elements.push(CigarElement::new(current_element_length, current_element_type));
                current_element_length = 0;
                current_element_type = next_element_type;
            }
            if current_element_type == CigarOp::D {
                ref_position += if self.is_forward_junction { -1 } else { 1 };
                mark_read_base_matches(&active_reads, current_element_type, current_index, NO_BASE);
                continue;
            }
            for read in active_reads.iter() {
                if read.operator() != Some(current_element_type) {
                    continue;
                }
                let base = read.current_base();
                let qual = read.current_qual();
                if !has_arrays {
                    if consensus_base == NO_BASE || (base != consensus_base && below_min_qual(consensus_max_qual as u8) && above_min_qual(qual)) {
                        consensus_base = base;
                        consensus_max_qual = qual as i32;
                        consensus_qual_total = qual as i32;
                        consensus_read_count = 1;
                        continue;
                    } else if base == consensus_base {
                        consensus_max_qual = consensus_max_qual.max(qual as i32);
                        consensus_qual_total += qual as i32;
                        consensus_read_count += 1;
                        continue;
                    }
                    has_arrays = true;
                    let base_index = assembly_utils::base_index(consensus_base);
                    total_quals[base_index] = consensus_qual_total;
                    max_quals[base_index] = consensus_max_qual;
                    read_counts[base_index] = consensus_read_count;
                }
                let base_index = assembly_utils::base_index(base);
                total_quals[base_index] += qual as i32;
                max_quals[base_index] = max_quals[base_index].max(qual as i32);
                read_counts[base_index] += 1;
            }
            if has_arrays {
                let mut max_qual = 0;
                let mut max_base_index = 0;
                for b in 0..DNA_BASE_COUNT {
                    if total_quals[b] > max_qual {
                        max_qual = total_quals[b];
                        max_base_index = b;
                    }
                }
                consensus_base = DNA_BASE_BYTES[max_base_index];
                consensus_max_qual = max_quals[max_base_index];
            }
            self.bases[current_index as usize] = consensus_base;
            self.base_quals[current_index as usize] = consensus_max_qual as u8;
            mark_read_base_matches(&active_reads, current_element_type, current_index, consensus_base);
            if current_element_type == CigarOp::M || current_element_type == CigarOp::D {
                ref_position += if self.is_forward_junction { -1 } else { 1 };
            }
        }
        self.cigar_elements.push(CigarElement::new(current_element_length, current_element_type));
        self.ref_base_position = ref_position;
        if self.cigar_elements.iter().any(|x| x.op == CigarOp::I && x.len >= MIN_INDEL_SUPPORT_LENGTH) {
            let mut outer_ref_position = -1;
            for read in &self.reads {
                if !read.is_valid() {
                    continue;
                }
                if self.is_forward_junction {
                    if outer_ref_position < 0 || read.ref_position() < outer_ref_position {
                        outer_ref_position = read.ref_position();
                    }
                } else {
                    outer_ref_position = outer_ref_position.max(read.ref_position());
                }
            }
            self.ref_base_position = outer_ref_position;
        }
        if self.is_forward_junction {
            self.cigar_elements.reverse();
        }
    }

    fn trim_final_sequence(&mut self) {
        let len = self.bases.len() as i32;
        let mut current_index = if self.is_forward_junction { len - 1 } else { 0 };
        let mut valid_length = 0;
        while current_index >= 0 && current_index < len {
            if self.bases[current_index as usize] == 0 {
                break;
            }
            valid_length += 1;
            current_index += if self.is_forward_junction { -1 } else { 1 };
        }
        if valid_length == len {
            return;
        }
        let reduction = len - valid_length;
        let start_index = if self.is_forward_junction { reduction } else { 0 };
        let end_index = start_index + valid_length;
        self.bases = self.bases[start_index as usize..end_index as usize].to_vec();
        self.base_quals = self.base_quals[start_index as usize..end_index as usize].to_vec();
    }
}

impl std::fmt::Display for RefBaseSeqBuilder {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "refPosition({} len={}) cigar({}) reads({})", self.ref_base_position, self.bases.len(), self.cigar_str(), self.reads.len())
    }
}

fn find_next_operator(reads: &[RpsP]) -> CigarOp {
    let mut inserts = 0;
    let mut deletes = 0;
    let mut aligned = 0;
    for read in reads {
        match read.operator() {
            Some(CigarOp::M) => aligned += 1,
            Some(CigarOp::D) => deletes += 1,
            Some(CigarOp::I) => inserts += 1,
            _ => {}
        }
    }
    if aligned >= inserts && aligned >= deletes {
        return CigarOp::M;
    }
    if inserts >= deletes { CigarOp::I } else { CigarOp::D }
}

fn progress_read_state(reads: &mut Vec<RpsP>, current_element_type: CigarOp) {
    let mut allow_non_indel_match = false;
    if current_element_type == CigarOp::I {
        let mut indel_count = 0;
        for read in reads.iter() {
            if read.operator() == Some(CigarOp::I) && read.element_length() >= MIN_INDEL_SUPPORT_LENGTH {
                indel_count += 1;
                if indel_count > reads.len() / 2 {
                    allow_non_indel_match = true;
                    break;
                }
            }
        }
    }
    let mut index = 0;
    while index < reads.len() {
        let read = reads[index].clone();
        if current_element_type == CigarOp::M || current_element_type == CigarOp::D {
            if read.operator() == Some(CigarOp::I) && read.element_length() < MIN_INDEL_SUPPORT_LENGTH {
                read.m().skip_insert();
            } else {
                read.m().move_next();
            }
        } else if read.operator() == Some(CigarOp::I) || (read.operator() == Some(CigarOp::M) && allow_non_indel_match) {
            read.m().move_next();
        }
        if exceeds_read_mismatches(&read) {
            read.m().mark_mismatched();
        }
        if read.exhausted() || read.mismatched() {
            reads.remove(index);
        } else {
            index += 1;
        }
    }
}

fn mark_read_base_matches(reads: &[RpsP], current_element_type: CigarOp, current_index: i32, consensus_base: u8) {
    for read in reads {
        if current_element_type == CigarOp::D {
            if read.operator() != Some(current_element_type) {
                let type_ = if read.operator() == Some(CigarOp::M) { SequenceDiffType::Base } else { SequenceDiffType::Insert };
                read.m().add_mismatch(current_index, type_);
            }
            continue;
        }
        let hq = above_min_qual(read.current_qual());
        let operator_mismatch = read.operator() != Some(current_element_type);
        if operator_mismatch && read.operator() == Some(CigarOp::D) {
            read.m().add_mismatch(current_index, SequenceDiffType::Delete);
            continue;
        }
        if read.current_base() != consensus_base {
            read.m().add_mismatch(current_index, SequenceDiffType::Base);
        } else {
            read.m().add_base_match(hq);
        }
    }
}
