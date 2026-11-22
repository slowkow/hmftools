//! JunctionReadTypes: sorts the reads at a soft-clip junction into junction, extension and other reads.

use crate::assembly::constants::*;
use crate::assembly::extension_seq_builder::{JunctionThresholdState, ReadJunctionInfo};
use crate::assembly::read::read::{Read, ReadP};
use crate::assembly::read::read_utils::{read_index_from_position, read_junction_extension_length, record_soft_clips_at_junction};
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::sam::NUM_MUTATONS_ATTRIBUTE;
use crate::assembly::types::junction::{JuncP, Junction};

const MIN_SNV_TYPE_EXTENSION: i32 = 3;

pub struct JunctionReadTypes<'a> {
    junction: JuncP,
    candidate_junction_reads: Vec<ReadJunctionInfo>,
    non_junction_reads: Vec<ReadP>,
    ref_genome: Option<&'a dyn RefGenomeInterface>,
}

impl<'a> JunctionReadTypes<'a> {
    pub fn new(junction: JuncP, raw_reads: &[ReadP], ref_genome: Option<&'a dyn RefGenomeInterface>) -> JunctionReadTypes<'a> {
        let mut j = JunctionReadTypes { junction, candidate_junction_reads: Vec::new(), non_junction_reads: Vec::new(), ref_genome };
        for read in raw_reads {
            match j.establish_junction_info(read) {
                Some(info) => j.candidate_junction_reads.push(info),
                None => j.non_junction_reads.push(read.clone()),
            }
        }
        j
    }
    pub fn non_junction_reads(&self) -> &Vec<ReadP> {
        &self.non_junction_reads
    }
    pub fn junction_reads(&self) -> Vec<ReadP> {
        self.candidate_junction_reads.iter().map(|x| x.read.clone()).collect()
    }
    pub fn candidate_junction_reads(&self) -> &Vec<ReadJunctionInfo> {
        &self.candidate_junction_reads
    }

    pub fn establish_extension_reads(&self, extension_reads: &mut Vec<ReadP>, junc_state: &mut JunctionThresholdState) {
        let candidates: Vec<&ReadJunctionInfo> = self.candidate_junction_reads.iter().filter(|x| x.matches_junction).collect();
        let candidate_count = candidates.len();
        let line_count = candidates.iter().filter(|x| x.read.has_line_tail(self.junction.is_forward())).count();
        if line_count as f64 >= 0.5 * candidate_count as f64 {
            junc_state.is_line = true;
            junc_state.min_extension_length = LINE_MIN_EXTENSION_LENGTH;
            junc_state.min_second_extension_length = LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH;
        } else {
            junc_state.min_extension_length = ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            junc_state.min_second_extension_length = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
        }
        let mut apply = |extension_reads: &mut Vec<ReadP>, junc_state: &mut JunctionThresholdState| {
            for read_info in &candidates {
                if read_info.extension_length < junc_state.min_second_extension_length {
                    continue;
                }
                if !junc_state.extension_length_valid && read_info.extension_length >= junc_state.min_extension_length {
                    junc_state.extension_length_valid = true;
                } else if !junc_state.second_extension_length_valid {
                    junc_state.second_extension_length_valid = true;
                }
                extension_reads.push(read_info.read.clone());
            }
        };
        apply(extension_reads, junc_state);
        if junc_state.is_line {
            return;
        }
        if junc_state.extension_length_valid {
            return;
        }
        junc_state.min_extension_length = ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER;
        junc_state.min_second_extension_length = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER;
        junc_state.uses_lower_saga_limits = true;
        extension_reads.clear();
        apply(extension_reads, junc_state);
    }

    fn establish_junction_info(&self, read: &ReadP) -> Option<ReadJunctionInfo> {
        let junction = &*self.junction;
        if record_soft_clips_at_junction(read, junction) {
            let extension_length = read_junction_extension_length(read, junction);
            return Some(ReadJunctionInfo { read: read.clone(), extension_length, matches_junction: true });
        }
        let mut extension_length = 0;
        let pos = junction.position;
        if junction.is_forward() {
            if read.is_right_clipped() && read.unclipped_end() > pos && (read.alignment_end() - pos).abs() <= ASSEMBLY_MAX_JUNC_POS_DIFF {
                extension_length = read.unclipped_end() - pos;
            } else if read.has_indel_implied_unclipped_end() && read.max_unclipped_end() > pos {
                extension_length = read.max_unclipped_end() - pos;
            } else if !read.is_right_clipped() && read.alignment_end() - pos >= MIN_SNV_TYPE_EXTENSION && self.read_mismatches_past_junction(read, junction) {
                extension_length = read.alignment_end() - pos;
            }
        } else if read.is_left_clipped() && read.unclipped_start() < pos && (read.alignment_start() - pos).abs() <= ASSEMBLY_MAX_JUNC_POS_DIFF {
            extension_length = pos - read.unclipped_start();
        } else if read.has_indel_implied_unclipped_start() && read.min_unclipped_start() < pos {
            extension_length = pos - read.min_unclipped_start();
        } else if !read.is_left_clipped() && pos - read.alignment_start() >= MIN_SNV_TYPE_EXTENSION && self.read_mismatches_past_junction(read, junction) {
            extension_length = pos - read.alignment_start();
        }
        if extension_length > 0 { Some(ReadJunctionInfo { read: read.clone(), extension_length, matches_junction: false }) } else { None }
    }

    pub fn read_mismatches_past_junction(&self, read: &Read, junction: &Junction) -> bool {
        let outer_alignment_pos = if junction.is_forward() { read.alignment_end() } else { read.alignment_start() };
        let extension_length = (outer_alignment_pos - junction.position).abs();
        if extension_length > ASSEMBLY_MAX_JUNC_POS_DIFF {
            return false;
        }
        match read.bam_record().get_integer_attribute(NUM_MUTATONS_ATTRIBUTE) {
            Some(n) if n >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH => {}
            _ => return false,
        }
        let pos_start;
        let pos_end;
        if junction.is_forward() {
            if read.alignment_end() <= junction.position {
                return false;
            }
            pos_start = junction.position + 1;
            pos_end = read.alignment_end();
        } else {
            if read.alignment_start() >= junction.position {
                return false;
            }
            pos_start = read.alignment_start();
            pos_end = junction.position - 1;
        }
        let ref_length = pos_end - pos_start + 1;
        if ref_length < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH {
            return false;
        }
        if ref_length > MIN_VARIANT_LENGTH / 2 {
            return false;
        }
        let Some(ref_bases) = self.ref_genome.expect("ref genome").get_bases(read.chromosome(), pos_start, pos_end) else { return false };
        if (ref_bases.len() as i32) < ref_length {
            return false;
        }
        let read_index_start = read_index_from_position(read, pos_start);
        let mut mismatches = 0;
        let mut i = 0;
        while i < ref_length && (i as usize) < ref_bases.len() {
            if read_index_start + i >= read.bases_length() {
                return false;
            }
            if read.get_bases()[(read_index_start + i) as usize] != ref_bases[i as usize] {
                mismatches += 1;
            }
            if mismatches >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH {
                return true;
            }
            i += 1;
        }
        false
    }
}

pub fn calc_proximate_junction_read_ratio(junction: &Junction, raw_reads: &[ReadP]) -> f64 {
    if !junction.soft_clip_based() {
        return 0.0;
    }
    let junc_position = junction.position;
    let mut proximate_count = 0;
    let mut exact_extension_count = 0;
    for read in raw_reads {
        if junction.is_forward() {
            if read.is_right_clipped() && read.unclipped_end() > junc_position && (read.alignment_end() - junc_position).abs() <= JUNCTION_PROXIMATE_READ_DISTANCE {
                proximate_count += 1;
                if read.alignment_end() == junc_position && read.unclipped_end() - junc_position >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH {
                    exact_extension_count += 1;
                }
            }
        } else if read.is_left_clipped() && read.unclipped_start() < junc_position && (read.alignment_start() - junc_position).abs() <= JUNCTION_PROXIMATE_READ_DISTANCE {
            proximate_count += 1;
            if read.alignment_start() == junc_position && junc_position - read.unclipped_start() >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH {
                exact_extension_count += 1;
            }
        }
    }
    if proximate_count > 0 { exact_extension_count as f64 / proximate_count as f64 } else { 0.0 }
}
