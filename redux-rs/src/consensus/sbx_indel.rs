//! SBX consensus over reads whose cigars differ (Java consensus/SbxIndelConsensus): the consensus cigar is
//! chosen element by element from the reads, then the bases are built over it.

use super::*;
use crate::sbx::{DEFAULT_SBX_MAX_DUPLICATE_DISTANCE, RAW_SIMPLEX_QUAL, SBX_DUPLEX_MISMATCH_QUAL};
use crate::sbx_dup::check_supplementary_consistency;

/// ConsensusState.setFromRead(read, true): the read's cigar, boundaries, bases and qualities
fn set_from_read_with_bases(state: &mut ConsensusState, read: &ReadView) {
    state.set_from_read(read);
    state.bases.copy_from_slice(&read.bases);
    for (q, r) in state.base_qualities.iter_mut().zip(read.quals.iter()) {
        *q = *r as i8 as i32;
    }
}

/// ReadValidReason.hasValidCigar
pub fn has_valid_cigar(cigar: &[u32]) -> bool {
    let count = cigar.len();
    if count > 1 {
        for (i, &e) in cigar.iter().enumerate() {
            let o = op(e);
            if i == 0 || i == count - 1 {
                if !aligned_or_clipped(o) && o != OP_I {
                    return false;
                }
            } else if o == OP_S || o == OP_H {
                return false;
            }
        }
        true
    } else {
        op(cigar[0]) == OP_M
    }
}

fn soft_clip_overlaps_indels(soft_clip_length: i32, other_cigar: &[u32], check_left: bool) -> bool {
    let mut read_bases = 0;
    let count = other_cigar.len() as i32;
    let mut cigar_index = if check_left { 0 } else { count - 1 };
    while cigar_index >= 0 && cigar_index < count {
        let e = other_cigar[cigar_index as usize];
        if consumes_read(op(e)) {
            if op(e) == OP_I && read_bases <= soft_clip_length {
                return true;
            }
            read_bases += oplen(e);
            if read_bases > soft_clip_length {
                return false;
            }
        }
        cigar_index += if check_left { 1 } else { -1 };
    }
    false
}

/// SbxIndelConsensus.checkSoftClipIndelPair: of two reads, the one whose soft clip covers an insert of the other
pub fn check_soft_clip_indel_pair(cigars: &[&[u32]]) -> Option<usize> {
    if cigars.len() != 2 {
        return None;
    }
    for i in 0..=1 {
        let soft_clip_read = cigars[i];
        let indel_read = cigars[1 - i];
        for check_left in [true, false] {
            let soft_clip_length = if check_left { left_soft_clip(soft_clip_read) } else { right_soft_clip(soft_clip_read) };
            if soft_clip_length == 0 {
                continue;
            }
            if soft_clip_overlaps_indels(soft_clip_length, indel_read, check_left) {
                return Some(i);
            }
        }
    }
    None
}

/// SbxIndelConsensus.checkNonOverlappingSupplementaries: when a read's ends are out of range of the template's,
/// the first read that is in range
/// The reads that Java's List.contains compares by value on the way, and so decodes, go to `eager_decoded`.
fn check_non_overlapping_supplementaries(reads: &[ReadView], names: &[String], eager_decoded: &mut Vec<usize>) -> Option<usize> {
    let records: Vec<&bam::Record> = reads.iter().map(|r| r.rec).collect();
    let non_consensus = check_supplementary_consistency(&records, DEFAULT_SBX_MAX_DUPLICATE_DISTANCE, names)?;
    for i in 0..reads.len() {
        let mut contained = false;
        for &n in &non_consensus {
            if n == i {
                contained = true;
                break;
            }
            if crate::sbx::fixed_fields_equal(records[i], records[n]) {
                eager_decoded.push(records[i] as *const bam::Record as usize);
                eager_decoded.push(records[n] as *const bam::Record as usize);
            }
        }
        if !contained {
            return Some(i);
        }
    }
    None
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SbxQualType {
    Simplex,
    DuplexLowQual,
    DuplexHighQual,
}

fn is_valid_operator_transition(current: u32, proposed: u32) -> bool {
    current == proposed || current == OP_M || proposed == OP_M
}

/// SbxIndelConsensus.getMajorityOperatorByType over the counts by (operator, quality type). Among the entries
/// with the highest count the operator that comes first in the CigarOperator enum is taken, so the order of the
/// Java map does not matter.
fn get_majority_operator_by_type(counts: &[(u32, SbxQualType, i32)], required: Option<SbxQualType>, min_count: f64) -> Option<u32> {
    let mut max_operator = 0;
    let mut max_count = 0;
    for &(operator, qual_type, count) in counts {
        if required.is_some_and(|r| r != qual_type) {
            continue;
        }
        if max_count == 0 || count > max_count || (count == max_count && operator <= max_operator) {
            max_count = count;
            max_operator = operator;
        }
    }
    if max_count as f64 > min_count { Some(max_operator) } else { None }
}

/// SbxIndelConsensus.determineConsensusCigarOperator: the operator of the reads at this position when they
/// agree; else the operator of more than half of the duplex reads with a high quality there (of the simplex
/// reads when no read is duplex there), else aligned when any read is aligned, else the most common
fn determine_consensus_cigar_operator(
    read_states: &[ReadParseState], active: &[usize], in_aligned_bases: bool, current_operator: Option<u32>, ref_position: i32,
    is_forward_consensus: bool,
) -> Option<u32> {
    let mut operator: Option<u32> = None;
    let mut matched = true;
    for &r in active {
        let read = &read_states[r];
        if read.before_unclipped_position(ref_position) {
            continue;
        }
        let read_operator = read.element_type;
        if current_operator.is_some_and(|c| !is_valid_operator_transition(c, read_operator)) {
            continue;
        }
        match operator {
            None => operator = Some(read_operator),
            Some(o) if o != read_operator => {
                matched = false;
                break;
            }
            _ => {}
        }
    }
    if matched {
        return operator;
    }

    let mut counts: Vec<(u32, SbxQualType, i32)> = Vec::with_capacity(4);
    let mut duplex_count = 0;
    let mut valid_read_count = 0;
    let mut has_aligned = false;
    for &r in active {
        let read = &read_states[r];
        if read.before_unclipped_position(ref_position) {
            continue;
        }
        let read_operator = read.element_type;
        if current_operator.is_some_and(|c| !is_valid_operator_transition(c, read_operator)) {
            continue;
        }
        if in_aligned_bases && read_operator == OP_S {
            // the first soft clip of a read is ignored once the aligned section has started
            if is_forward_consensus && read.cigar_index == 0 {
                continue;
            } else if !is_forward_consensus && read.cigar_index > 0 {
                continue;
            }
        }
        valid_read_count += 1;
        has_aligned |= read_operator == OP_M;
        let qual_type = if read.base_qual() == RAW_SIMPLEX_QUAL as i32 {
            SbxQualType::Simplex
        } else {
            duplex_count += 1;
            // an element without read bases is in both strands of the duplex read
            if consumes_read(read_operator) && read.base_qual() <= SBX_DUPLEX_MISMATCH_QUAL as i32 {
                SbxQualType::DuplexLowQual
            } else {
                SbxQualType::DuplexHighQual
            }
        };
        match counts.iter_mut().find(|c| c.0 == read_operator && c.1 == qual_type) {
            Some(c) => c.2 += 1,
            None => counts.push((read_operator, qual_type, 1)),
        }
    }
    let consensus_operator = if duplex_count == 0 {
        get_majority_operator_by_type(&counts, Some(SbxQualType::Simplex), valid_read_count as f64 * 0.5)
    } else {
        get_majority_operator_by_type(&counts, Some(SbxQualType::DuplexHighQual), duplex_count as f64 * 0.5)
    };
    if consensus_operator.is_some() {
        return consensus_operator;
    }
    if has_aligned {
        return Some(OP_M);
    }
    get_majority_operator_by_type(&counts, None, 0.0)
}

/// SbxIndelConsensus.determineConsensusCigar: walks all reads from the outermost unclipped position, choosing
/// the operator at each step. None when no operator can be chosen.
fn determine_consensus_cigar(read_states: &mut [ReadParseState], is_forward_consensus: bool, initial_unclipped_position: i32) -> Option<Vec<u32>> {
    let mut cigar: Vec<u32> = Vec::new();
    let mut active: Vec<usize> = (0..read_states.len()).collect();
    let mut element_length: i32 = 0;
    let mut element_operator: Option<u32> = None;
    // no way back from aligned to the first soft clip, nor from the last soft clip to aligned
    let mut in_aligned_bases = false;
    let mut in_final_soft_clip = false;
    let mut ref_position = initial_unclipped_position;
    let add = |cigar: &mut Vec<u32>, len: i32, o: u32| {
        let e = (len as u32) << 4 | o;
        if is_forward_consensus {
            cigar.push(e);
        } else {
            cigar.insert(0, e);
        }
    };
    loop {
        let next_operator = if in_final_soft_clip {
            OP_S
        } else {
            determine_consensus_cigar_operator(read_states, &active, in_aligned_bases, element_operator, ref_position, is_forward_consensus)?
        };
        if element_operator != Some(next_operator) {
            if let Some(o) = element_operator {
                add(&mut cigar, element_length, o);
            }
            element_length = 0;
            element_operator = Some(next_operator);
            if !in_aligned_bases && next_operator == OP_M {
                in_aligned_bases = true;
            } else if in_aligned_bases && next_operator == OP_S {
                in_final_soft_clip = true;
            }
        }
        element_length += 1;

        // each read moves on unless an insert was chosen that it does not have
        let mut has_exhausted = false;
        for &r in &active {
            let read = &mut read_states[r];
            if read.before_unclipped_position(ref_position) {
                continue;
            }
            if consumes_ref(next_operator) && read.element_type == OP_I {
                read.skip_insert_element();
            }
            let move_next = !(read.element_type != next_operator && next_operator == OP_I);
            if move_next {
                read.move_next_base();
            }
            has_exhausted |= read.exhausted;
        }
        if has_exhausted {
            active.retain(|&r| !read_states[r].exhausted);
            if active.is_empty() {
                break;
            }
        }
        if consumes_ref_or_unclipped(next_operator) {
            ref_position += if is_forward_consensus { 1 } else { -1 };
        }
    }
    add(&mut cigar, element_length, element_operator.unwrap());
    Some(cigar)
}

impl ConsensusReads {
    /// SbxIndelConsensus.determineIndelConsensus
    pub(super) fn sbx_determine_indel_consensus(&mut self, reads: &[ReadView], state: &mut ConsensusState, names: &[String]) {
        let cigars: Vec<&[u32]> = reads.iter().map(|r| r.cigar).collect();
        if let Some(i) = check_soft_clip_indel_pair(&cigars) {
            set_from_read_with_bases(state, &reads[i]);
            state.outcome = ConsensusOutcome::IndelSoftclip;
            return;
        }
        if let Some(i) = check_non_overlapping_supplementaries(reads, names, &mut self.sbx_eager_decoded) {
            set_from_read_with_bases(state, &reads[i]);
            state.outcome = ConsensusOutcome::IndelMismatch;
            return;
        }
        let is_forward = state.is_forward;
        let mut read_states: Vec<ReadParseState> = reads.iter().map(|r| ReadParseState::new(r, is_forward)).collect();
        let mut outer_unclipped_position = -1;
        for r in reads {
            if is_forward {
                let p = r.unclipped_start();
                if outer_unclipped_position < 0 || p < outer_unclipped_position {
                    outer_unclipped_position = p;
                }
            } else {
                outer_unclipped_position = outer_unclipped_position.max(r.unclipped_end());
            }
        }
        let initial_ref_position = outer_unclipped_position;
        for rs in read_states.iter_mut() {
            rs.move_to_ref_position(initial_ref_position);
        }
        let Some(cigar) = determine_consensus_cigar(&mut read_states, is_forward, initial_ref_position) else {
            // Java logs this, then fails on the empty cigar and exits when the group has no consensus read
            crate::error!(
                "indel consensus, invalid next operator, sample read: {}", String::from_utf8_lossy(reads[0].rec.qname())
            );
            std::process::exit(1);
        };
        if !has_valid_cigar(&cigar) {
            state.outcome = ConsensusOutcome::IndelFail;
            return;
        }
        let base_length: i32 = cigar.iter().filter(|e| consumes_read(op(**e))).map(|e| oplen(*e)).sum();
        state.set_base_length(base_length as usize);
        let mut base_index: i32 = if is_forward { 0 } else { base_length - 1 };
        let mut ref_position = initial_ref_position;
        let cigar_count = cigar.len() as i32;
        let mut cigar_index: i32 = if is_forward { 0 } else { cigar_count - 1 };
        for rs in read_states.iter_mut() {
            *rs = ReadParseState::new(rs.read, is_forward);
            rs.move_to_ref_position(initial_ref_position);
        }
        let chromosome_length = if self.chromosome_length == 0 {
            self.refs.chromosome_length(&names[reads[0].rec.tid() as usize])
        } else {
            self.chromosome_length
        };
        let dir = if is_forward { 1 } else { -1 };
        while cigar_index >= 0 && cigar_index < cigar_count {
            let element = cigar[cigar_index as usize];
            self.add_element_bases(state, &mut read_states, element, base_index, ref_position, false, &[], chromosome_length);
            if state.outcome == ConsensusOutcome::IndelFail {
                break;
            }
            let (o, n) = (op(element), oplen(element));
            if consumes_read(o) {
                base_index += dir * n;
            }
            if consumes_ref_or_unclipped(o) {
                ref_position += dir * n;
            }
            cigar_index += dir;
        }
        state.finalise_cigar();

        let ref_base_length = ref_length(&cigar);
        let (first, last) = (cigar[0], cigar[cigar.len() - 1]);
        let (alignment_start, alignment_end, unclipped_start, unclipped_end);
        if is_forward {
            unclipped_start = outer_unclipped_position;
            alignment_start = outer_unclipped_position + if op(first) == OP_S { oplen(first) } else { 0 };
            alignment_end = alignment_start + ref_base_length - 1;
            unclipped_end = alignment_end + if op(last) == OP_S { oplen(last) } else { 0 };
        } else {
            unclipped_end = outer_unclipped_position;
            alignment_end = unclipped_end - if op(last) == OP_S { oplen(last) } else { 0 };
            alignment_start = alignment_end - ref_base_length + 1;
            unclipped_start = alignment_start - if op(first) == OP_S { oplen(first) } else { 0 };
        }
        state.unclipped_pos_start = unclipped_start;
        state.unclipped_pos_end = unclipped_end;
        state.alignment_start = alignment_start;
        state.alignment_end = alignment_end;
        if state.outcome != ConsensusOutcome::IndelFail {
            state.outcome = ConsensusOutcome::IndelMismatch;
        }
    }
}
