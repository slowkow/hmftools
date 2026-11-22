//! DiscordantJunctionBuilder: places a junction for a group of discordant reads without soft clips.

use std::collections::HashSet;

use crate::assembly::constants::*;
use crate::assembly::globals::max_observed_concordant_frag_length;
use crate::assembly::ptr::P;
use crate::assembly::read::read::ReadP;
use crate::assembly::remote_region_finder::add_or_create_mate_remote_region;
use crate::assembly::types::junction::{JuncP, Junction};
use crate::assembly::types::remote_region::{RemoteRegion, RemoteRegionP};
use crate::common::java_sort::java_sort;

pub struct DiscordantJunctionBuilder {
    junction: JuncP,
    new_junction: Option<JuncP>,
}

impl DiscordantJunctionBuilder {
    pub fn new(junction: JuncP) -> Self {
        DiscordantJunctionBuilder { junction, new_junction: None }
    }
    pub fn has_new_junction(&self) -> bool {
        self.new_junction.is_some()
    }
    pub fn new_junction(&self) -> Option<JuncP> {
        self.new_junction
    }

    /// `non_junction_reads` is the list the Java constructor is given
    pub fn assess_discordant_junction(&mut self, raw_reads: &[ReadP], extension_reads: &mut Vec<ReadP>, junction_reads: &mut Vec<ReadP>, non_junction_reads: &mut Vec<ReadP>) {
        if (raw_reads.len() as i32) < ASSEMBLY_MIN_READ_SUPPORT {
            return;
        }
        let junction = self.junction;
        let original_junc_position = junction.position;
        let candidate_junction_reads: Vec<ReadP> = raw_reads
            .iter()
            .filter(|x| (junction.is_forward() && x.alignment_end() == original_junc_position) || (junction.is_reverse() && x.alignment_start() == original_junc_position))
            .cloned()
            .collect();
        if candidate_junction_reads.is_empty() {
            return;
        }
        let mut remote_regions: Vec<RemoteRegionP> = Vec::new();
        for x in &candidate_junction_reads {
            add_or_create_mate_remote_region(&mut remote_regions, x, true);
        }
        if remote_regions.is_empty() {
            return;
        }
        if remote_regions.len() > 1 {
            RemoteRegion::merge_regions(&mut remote_regions);
            java_sort(&mut remote_regions, |a, b| (-a.read_count()).cmp(&(-b.read_count())));
        }
        let main_remote_region = remote_regions[0];
        let mut extension_junc_positions: Vec<i32> = Vec::with_capacity(raw_reads.len());
        let mut candidate_read_ids: HashSet<String> = HashSet::new();
        let min_read_pos_start;
        let max_read_pos_end;
        if junction.is_forward() {
            min_read_pos_start = 1.max(original_junc_position - max_observed_concordant_frag_length());
            max_read_pos_end = original_junc_position;
        } else {
            min_read_pos_start = original_junc_position;
            max_read_pos_end = original_junc_position + max_observed_concordant_frag_length();
        }
        for read in raw_reads {
            if read.mapping_quality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY {
                continue;
            }
            match read.mate_chromosome() {
                Some(mc) if main_remote_region.overlaps(mc, read.mate_alignment_start(), read.mate_alignment_end()) => {}
                _ => continue,
            }
            if read.alignment_start() < min_read_pos_start || read.alignment_end() > max_read_pos_end {
                continue;
            }
            candidate_read_ids.insert(read.id().to_string());
            if junction.is_forward() {
                if read.alignment_end() > original_junc_position {
                    continue;
                }
                extension_junc_positions.push(read.unclipped_end());
            } else {
                if read.alignment_start() < original_junc_position {
                    continue;
                }
                extension_junc_positions.push(read.unclipped_start());
            }
        }
        if (extension_junc_positions.len() as i32) < ASSEMBLY_MIN_READ_SUPPORT {
            return;
        }
        if junction.is_forward() {
            extension_junc_positions.sort_by(|a, b| b.cmp(a));
        } else {
            extension_junc_positions.sort();
        }
        let mut outer_junc_position = extension_junc_positions[0];
        let mut second_junc_position = extension_junc_positions[(ASSEMBLY_MIN_READ_SUPPORT - 1) as usize];
        let adjusted_junc_position;
        if junction.is_forward() {
            outer_junc_position -= ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            second_junc_position -= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
            adjusted_junc_position = outer_junc_position.min(second_junc_position);
        } else {
            outer_junc_position += ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            second_junc_position += ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
            adjusted_junc_position = outer_junc_position.max(second_junc_position);
        }
        let mut new_junction = Junction::new(junction.chromosome, adjusted_junc_position, junction.orient, true, false, false, 0, None);
        new_junction.set_raw_discordant_position(original_junc_position);
        self.new_junction = Some(P::new(new_junction));

        for read in raw_reads {
            if read.alignment_start() < min_read_pos_start || read.alignment_end() > max_read_pos_end {
                continue;
            }
            if junction.is_forward() {
                if read.has_indel_implied_unclipped_end() {
                    continue;
                }
            } else if read.has_indel_implied_unclipped_start() {
                continue;
            }
            if !candidate_read_ids.contains(read.id()) || read.mapping_quality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY {
                non_junction_reads.push(read.clone());
                continue;
            }
            let extension_length;
            if junction.is_forward() {
                if read.alignment_end() > original_junc_position {
                    continue;
                }
                extension_length = if read.alignment_end() >= adjusted_junc_position { read.unclipped_end() - adjusted_junc_position } else { -1 };
            } else {
                if read.alignment_start() < original_junc_position {
                    continue;
                }
                extension_length = if read.alignment_start() <= adjusted_junc_position { adjusted_junc_position - read.unclipped_start() } else { -1 };
            }
            if extension_length >= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH {
                extension_reads.push(read.clone());
            } else if extension_length > 0 {
                junction_reads.push(read.clone());
            } else {
                non_junction_reads.push(read.clone());
            }
        }
    }
}
