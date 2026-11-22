//! RemoteRegion: a region away from the junction where mates or supplementaries of an assembly's reads align.

use std::cmp::Ordering;

use crate::assembly::common_utils::positions_overlap;
use crate::assembly::constants::*;
use crate::assembly::ptr::P;
use crate::assembly::types::junction_assembly::JunctionAssembly;
use crate::assembly::types::RemoteReadType;
use crate::common::compat;
use crate::common::java_hash::JavaHashSet;
use crate::common::java_sort::java_sort;
use crate::common::regions::compare_chromosomes;

pub type RemoteRegionP = P<RemoteRegion>;

pub struct RemoteRegion {
    pub chromosome: &'static str,
    start: i32,
    end: i32,
    read_ids: JavaHashSet<String>,
    read_type_count: [i32; 3],
    soft_clip_map_qual_total: i32,
    has_high_qual_discordant_read: bool,
}

/// ChrBaseRegion.compareTo
pub fn compare_chr_regions(chr1: &str, start1: i32, chr2: &str, start2: i32) -> Ordering {
    if chr1 == chr2 {
        return start1.cmp(&start2);
    }
    compare_chromosomes(chr1, chr2)
}

impl RemoteRegion {
    pub fn new(chromosome: &'static str, start: i32, end: i32, read_id: &str, read_type: RemoteReadType) -> RemoteRegion {
        let mut read_type_count = [0; 3];
        read_type_count[read_type.ordinal()] += 1;
        RemoteRegion {
            chromosome,
            start,
            end,
            read_ids: JavaHashSet::of_one(read_id.to_string()),
            read_type_count,
            soft_clip_map_qual_total: 0,
            has_high_qual_discordant_read: false,
        }
    }
    pub fn start(&self) -> i32 {
        self.start
    }
    pub fn end(&self) -> i32 {
        self.end
    }
    pub fn set_start(&mut self, v: i32) {
        self.start = v;
    }
    pub fn set_end(&mut self, v: i32) {
        self.end = v;
    }
    pub fn add_read_details(&mut self, read_id: &str, pos_start: i32, pos_end: i32, read_type: RemoteReadType) {
        self.start = self.start.min(pos_start);
        self.end = self.end.max(pos_end);
        self.read_ids.insert(read_id.to_string());
        self.read_type_count[read_type.ordinal()] += 1;
    }
    pub fn read_ids(&self) -> &JavaHashSet<String> {
        &self.read_ids
    }
    pub fn read_ids_mut(&mut self) -> &mut JavaHashSet<String> {
        &mut self.read_ids
    }
    pub fn read_count(&self) -> i32 {
        self.read_ids.len() as i32
    }
    pub fn read_type_counts(&self) -> &[i32; 3] {
        &self.read_type_count
    }
    pub fn read_type_counts_mut(&mut self) -> &mut [i32; 3] {
        &mut self.read_type_count
    }
    pub fn non_supp_read_count(&self) -> i32 {
        self.read_type_count[RemoteReadType::JunctionMate.ordinal()] + self.read_type_count[RemoteReadType::Discordant.ordinal()]
    }
    pub fn has_discordant_reads(&self) -> bool {
        self.read_type_count[RemoteReadType::Discordant.ordinal()] > 0
    }
    pub fn set_has_high_qual_discordant_read(&mut self) {
        self.has_high_qual_discordant_read = true;
    }
    pub fn has_high_qual_discordant_read(&self) -> bool {
        self.has_high_qual_discordant_read
    }
    pub fn is_supp_only_region(&self) -> bool {
        self.non_supp_read_count() == 0
    }
    pub fn soft_clip_map_qual_total(&self) -> i32 {
        self.soft_clip_map_qual_total
    }
    pub fn add_soft_clip_map_qual(&mut self, soft_clip_length: i32, map_qual: i32) {
        self.soft_clip_map_qual_total += soft_clip_length * map_qual;
    }
    pub fn matches(&self, other: &RemoteRegion) -> bool {
        self.overlaps(other.chromosome, other.start, other.end)
    }
    pub fn overlaps(&self, other_chr: &str, other_pos_start: i32, other_pos_end: i32) -> bool {
        self.chromosome == other_chr && positions_overlap(self.start, self.end, other_pos_start, other_pos_end)
    }
    pub fn overlaps_assembly(&self, assembly: &JunctionAssembly) -> bool {
        self.overlaps(assembly.junction().chromosome, assembly.min_aligned_position(), assembly.max_aligned_position())
    }
    pub fn has_read_id(&self, full_read_id: &str) -> bool {
        self.read_ids.iter_unordered().any(|x| x == full_read_id)
    }

    pub fn merge_regions(regions: &mut Vec<RemoteRegionP>) {
        Self::merge_regions_buffer(regions, REMOTE_REGION_MERGE_MARGIN)
    }

    pub fn merge_regions_buffer(regions: &mut Vec<RemoteRegionP>, distance_buffer: i32) {
        java_sort(regions, |a, b| compare_chr_regions(a.chromosome, a.start, b.chromosome, b.start));
        let java_a9 = compat::java("A9");
        let mut index = 0usize;
        while index + 1 < regions.len() {
            let region = regions[index];
            let mut next_index = index + 1;
            while next_index < regions.len() {
                let next_region = regions[next_index];
                let merge = region.chromosome == next_region.chromosome
                    && (positions_overlap(region.start, region.end, next_region.start, next_region.end) || region.end >= next_region.start - distance_buffer);
                if merge {
                    regions.remove(next_index);
                    let r = region.m();
                    r.end = r.end.max(next_region.end);
                    for id in next_region.read_ids.iter() {
                        r.read_ids.insert(id.clone());
                    }
                    for i in 0..r.read_type_count.len() {
                        r.read_type_count[i] += next_region.read_type_count[i];
                        if java_a9 {
                            // Java: inside the loop, so the other region's total is added once per read type
                            r.add_soft_clip_map_qual(1, next_region.soft_clip_map_qual_total);
                        }
                    }
                    if !java_a9 {
                        r.add_soft_clip_map_qual(1, next_region.soft_clip_map_qual_total);
                    }
                } else {
                    next_index += 1;
                }
            }
            index += 1;
        }
    }

    pub fn purge_weak_supplementary_regions(regions: &mut Vec<RemoteRegionP>) {
        let mut non_supp_support = 0;
        let mut max_soft_clip_map_qual = 0;
        let mut has_supp_only_regions = false;
        for region in regions.iter() {
            has_supp_only_regions |= region.is_supp_only_region();
            non_supp_support += region.non_supp_read_count();
            max_soft_clip_map_qual = max_soft_clip_map_qual.max(region.soft_clip_map_qual_total());
        }
        if !has_supp_only_regions {
            return;
        }
        let mut index = 0;
        while index < regions.len() {
            let region = regions[index];
            if region.is_supp_only_region()
                && ((region.soft_clip_map_qual_total() as f64) < REMOTE_REGION_WEAK_SUPP_PERCENT * max_soft_clip_map_qual as f64
                    || (region.read_count() as f64) < REMOTE_REGION_WEAK_SUPP_PERCENT * non_supp_support as f64)
            {
                regions.remove(index);
                continue;
            }
            index += 1;
        }
    }

    pub fn purge_low_qual_discordant_only_regions(regions: &mut Vec<RemoteRegionP>) {
        let mut index = 0;
        while index < regions.len() {
            let region = regions[index];
            let is_discordant_only =
                region.read_type_count[RemoteReadType::JunctionMate.ordinal()] == 0 && region.read_type_count[RemoteReadType::Supplementary.ordinal()] == 0;
            if is_discordant_only && !region.has_high_qual_discordant_read() {
                regions.remove(index);
                continue;
            }
            index += 1;
        }
    }
}

impl std::fmt::Display for RemoteRegion {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} reads({}) counts(mate={} supp={} disc={}) softClipMapQual({})",
            self.chromosome,
            self.start,
            self.end,
            self.read_ids.len(),
            self.read_type_count[RemoteReadType::JunctionMate.ordinal()],
            self.read_type_count[RemoteReadType::Supplementary.ordinal()],
            self.read_type_count[RemoteReadType::Discordant.ordinal()],
            self.soft_clip_map_qual_total
        )
    }
}
