//! Discordant read groups (DiscordantGroups, DiscordantGroup, DiscordantRemoteRegion) and the counts of
//! discordant fragment types (DiscordantStats).

use super::depth::DepthTracker;
use super::filters::*;
use super::hotspot::KnownHotspot;
use super::rec::*;
use super::types::*;
use std::sync::Arc;

pub const DISCORDANT_GROUP_MIN_FRAGMENTS: usize = 3;
pub const DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT: usize = 5;
pub const DISCORDANT_GROUP_MAX_LOCAL_LENGTH: i32 = 1_000_000;
pub const DISCORDANT_GROUP_MIN_MAP_QUAL: i32 = 40;
pub const DISCORDANT_GROUP_MIN_ALIGN_SCORE: i32 = 75;
pub const REMOTE_REGION_MERGE_MARGIN: i32 = 500;
/// PrepConstants: the depth from which a site's support is tested against it, and the share of the depth a
/// discordant group must have
pub const DEPTH_MIN_CHECK: i32 = 100;
pub const DEPTH_MIN_SUPPORT_RATIO_DISCORDANT: f64 = 0.01;

pub const DISCORDANT_FRAG_TYPES: [&str; 11] =
    ["Translocation", "InvLt1K", "Inv1To5K", "Inv5To100K", "InvGt100K", "Del1To5K", "Del5To100K", "DelGt100K", "Dup1To5K", "Dup5To100K", "DupGt100K"];

const TRANSLOCATION: usize = 0;
const INV_LT_1K: usize = 1;
const INV_1_TO_5K: usize = 2;
const INV_5_TO_100K: usize = 3;
const INV_GT_100K: usize = 4;
const DEL_1_TO_5K: usize = 5;
const DEL_5_TO_100K: usize = 6;
const DEL_GT_100K: usize = 7;
const DUP_1_TO_5K: usize = 8;
const DUP_5_TO_100K: usize = 9;
const DUP_GT_100K: usize = 10;

#[derive(Default, Clone, Debug, PartialEq, Eq)]
pub struct DiscordantStats {
    pub type_counts: [i64; 11],
}

impl DiscordantStats {
    pub fn add(&mut self, other: &DiscordantStats) {
        for i in 0..self.type_counts.len() {
            self.type_counts[i] += other.type_counts[i];
        }
    }

    pub fn has_primary_with_supplementary(g: &ReadGroup, reads: &[PrepRead]) -> bool {
        g.reads.iter().any(|&i| !reads[i as usize].is_supplementary() && reads[i as usize].has_supp_alignment())
    }

    /// (read position, supplementary position) used for the span of a primary and its supplementary
    fn form_position_pairing(read: &PrepRead, contigs: &Contigs) -> (i32, i32) {
        let supp = read.supp.as_ref().unwrap();
        if contigs.name(read.chr) == supp.chr_name {
            if read.sc_left > read.sc_right {
                (read.end, supp.position)
            } else {
                (read.start, read_boundary_position(supp.position, supp.cigar.as_bytes(), false, false))
            }
        } else {
            (read.start, supp.position)
        }
    }

    fn primary_supp_inferred_insert_size(read: &PrepRead, pairing: (i32, i32), contigs: &Contigs) -> i32 {
        if contigs.name(read.chr) == read.supp.as_ref().unwrap().chr_name { (pairing.0 - pairing.1).wrapping_abs() } else { 0 }
    }

    pub fn process_read_group(&mut self, g: &ReadGroup, reads: &[PrepRead], is_discordant_pair: bool, contigs: &Contigs) {
        let Some(read) = g.reads.iter().map(|&i| &reads[i as usize]).find(|r| !r.is_supplementary()) else { return };
        if read.is_paired() {
            if is_discordant_pair {
                self.add_to_counts(
                    read.chr == read.mate_chr,
                    read.start,
                    read.mate_pos_start,
                    read.orientation(),
                    read.mate_orientation(),
                    inferred_insert_size_abs(&read.rec),
                );
            } else {
                // the longest primary to supplementary pairing
                let mut max_insert_size = -1;
                let mut max_read: Option<(&PrepRead, (i32, i32))> = None;
                for &i in &g.reads {
                    let primary = &reads[i as usize];
                    if !primary.has_supp_alignment() || primary.is_unmapped() {
                        continue;
                    }
                    let pairing = Self::form_position_pairing(primary, contigs);
                    let size = Self::primary_supp_inferred_insert_size(primary, pairing, contigs);
                    if size > max_insert_size {
                        max_read = Some((primary, pairing));
                        max_insert_size = size;
                    }
                }
                if let Some((primary, pairing)) = max_read {
                    let supp = primary.supp.as_ref().unwrap();
                    self.add_to_counts(
                        contigs.name(primary.chr) == supp.chr_name,
                        pairing.0,
                        pairing.1,
                        primary.orientation(),
                        supp.orient().opposite(),
                        max_insert_size,
                    );
                }
            }
        } else if let Some(supp) = &read.supp {
            let pairing = Self::form_position_pairing(read, contigs);
            let size = Self::primary_supp_inferred_insert_size(read, pairing, contigs);
            self.add_to_counts(contigs.name(read.chr) == supp.chr_name, pairing.0, pairing.1, read.orientation(), supp.orient().opposite(), size);
        }
    }

    fn add_to_counts(&mut self, same_chromosome: bool, pos_start1: i32, pos_start2: i32, orientation1: Orientation, orientation2: Orientation, insert_size: i32) {
        const LENGTH_SHORT: i32 = 1000;
        const LENGTH_5K: i32 = 5000;
        const LENGTH_100K: i32 = 100_000;
        let c = &mut self.type_counts;
        if !same_chromosome {
            c[TRANSLOCATION] += 1;
            return;
        }
        let distance = insert_size;
        if orientation1 == orientation2 {
            if distance < LENGTH_SHORT {
                c[INV_LT_1K] += 2; // the mate is in the same group and is not counted again
            } else if distance < LENGTH_5K {
                c[INV_1_TO_5K] += 1;
            } else if distance < LENGTH_100K {
                c[INV_5_TO_100K] += 1;
            } else {
                c[INV_GT_100K] += 1;
            }
        } else if (pos_start1 < pos_start2) == orientation1.is_forward() {
            if distance < LENGTH_5K {
                c[DEL_1_TO_5K] += 1;
            } else if distance < LENGTH_100K {
                c[DEL_5_TO_100K] += 1;
            } else {
                c[DEL_GT_100K] += 1;
            }
        } else if distance < LENGTH_5K {
            c[DUP_1_TO_5K] += 1;
        } else if distance < LENGTH_100K {
            c[DUP_5_TO_100K] += 1;
        } else {
            c[DUP_GT_100K] += 1;
        }
    }

    /// The disc_stats file (EsveeDiscordantStats.write)
    pub fn file_text(&self, total_reads: i64, prep_reads: i64) -> String {
        let mut header = vec!["TotalReads".to_string(), "PrepReads".to_string()];
        header.extend(DISCORDANT_FRAG_TYPES.iter().map(|s| s.to_string()));
        let mut values = vec![total_reads.to_string(), prep_reads.to_string()];
        values.extend(self.type_counts.iter().map(|c| c.to_string()));
        format!("{}\n{}\n", header.join("\t"), values.join("\t"))
    }
}

struct RemoteRegion {
    chr: i32,
    start: i32,
    end: i32,
    orient: Orientation,
    read_groups: Vec<GroupId>,
}

struct DiscordantGroup {
    region: Region,
    orient: Orientation,
    innermost_read: Option<ReadId>,
    read_groups: Vec<GroupId>,
    remote_regions: Vec<RemoteRegion>,
}

fn first_primary_read(g: &ReadGroup, reads: &[PrepRead]) -> Option<ReadId> {
    g.reads.iter().copied().find(|&i| !reads[i as usize].is_supplementary())
}

impl DiscordantGroup {
    fn from_read_group(gid: GroupId, groups: &[ReadGroup], reads: &[PrepRead]) -> DiscordantGroup {
        let first = &reads[first_primary_read(&groups[gid as usize], reads).expect("primary read") as usize];
        let mut dg = DiscordantGroup {
            region: Region::new(first.chr, first.start, first.end),
            orient: first.orientation(),
            innermost_read: None,
            read_groups: vec![gid],
            remote_regions: Vec::new(),
        };
        dg.add_remote_region_read(gid, first);
        dg
    }

    fn inner_position(&self, reads: &[PrepRead]) -> i32 {
        match self.innermost_read {
            None => {
                if self.orient.is_forward() {
                    self.region.end
                } else {
                    self.region.start
                }
            }
            Some(r) => {
                if self.orient.is_forward() {
                    reads[r as usize].end
                } else {
                    reads[r as usize].start
                }
            }
        }
    }

    fn try_add_read_group(&mut self, gid: GroupId, read: &PrepRead) -> bool {
        if read.orientation() != self.orient {
            return false;
        }
        if !self.region.overlaps(read.chr, read.start, read.end) {
            return false;
        }
        self.read_groups.push(gid);
        if read.end > self.region.end {
            self.region.end = read.end;
        }
        if read.start < self.region.start {
            self.region.start = read.start;
        }
        self.add_remote_region_read(gid, read);
        true
    }

    fn add_remote_region_read(&mut self, gid: GroupId, read: &PrepRead) {
        let remote_chr = read.mate_chr;
        let remote_start = read.mate_pos_start;
        let remote_end = remote_start + read.bases.len() as i32;
        let remote_orient = read.mate_orientation();
        let found = self.remote_regions.iter_mut().find(|x| x.chr == remote_chr && positions_overlap(x.start, x.end, remote_start, remote_end) && x.orient == remote_orient);
        match found {
            Some(r) => {
                r.start = r.start.min(remote_start);
                r.end = r.end.max(remote_end);
                r.read_groups.push(gid);
            }
            None => self.remote_regions.push(RemoteRegion { chr: remote_chr, start: remote_start, end: remote_end, orient: remote_orient, read_groups: vec![gid] }),
        }
    }

    /// purgeReads: drops the given groups from the group and from the first remote region that holds each
    fn purge_reads(&mut self, excluded: &mut Vec<GroupId>) {
        if excluded.is_empty() {
            return;
        }
        let mut index = 0;
        while index < self.read_groups.len() {
            let gid = self.read_groups[index];
            if let Some(p) = excluded.iter().position(|&x| x == gid) {
                self.read_groups.remove(index);
                excluded.remove(p);
                if let Some(r) = self.remote_regions.iter_mut().find(|x| x.read_groups.contains(&gid)) {
                    let p = r.read_groups.iter().position(|&x| x == gid).unwrap();
                    r.read_groups.remove(p);
                }
                if excluded.is_empty() {
                    return;
                }
            } else {
                index += 1;
            }
        }
    }
}

/// DiscordantRemoteRegion.mergeRemoteRegions: sorted by chromosome and start, regions of one orientation that
/// overlap or lie within the margin become one
fn merge_remote_regions(regions: &mut Vec<RemoteRegion>, contigs: &Contigs) {
    regions.sort_by(|a, b| {
        if a.chr == b.chr {
            a.start.cmp(&b.start)
        } else {
            let (na, nb) = (contigs.name(a.chr), contigs.name(b.chr));
            let (ra, rb) = (chromosome_rank(na), chromosome_rank(nb));
            if ra == rb { na.cmp(nb) } else { ra.cmp(&rb) }
        }
    });
    let mut index = 0;
    while index + 1 < regions.len() {
        let mut next_index = index + 1;
        while next_index < regions.len() {
            let (region, next) = (&regions[index], &regions[next_index]);
            let merge = region.chr == next.chr
                && region.orient == next.orient
                && (positions_overlap(region.start, region.end, next.start, next.end) || region.end >= next.start - REMOTE_REGION_MERGE_MARGIN);
            if merge {
                let next = regions.remove(next_index);
                let region = &mut regions[index];
                region.end = region.end.max(next.end);
                region.read_groups.extend(next.read_groups);
            } else {
                next_index += 1;
            }
        }
        index += 1;
    }
}

pub struct DiscordantGroups {
    region: Region,
    max_concordant_fragment_length: i32,
    min_discordant_fragment_length: i32,
    known_hotspots: Vec<Arc<KnownHotspot>>,
    contigs: Arc<Contigs>,
    /// junctions made at a remote region that lies outside the partition, and those of them on another
    /// chromosome; Java writes both with the partition's chromosome (finding P8)
    pub remote_outside_partition: std::cell::Cell<u32>,
    pub remote_other_chromosome: std::cell::Cell<u32>,
}

impl DiscordantGroups {
    pub fn new(region: Region, observed_max_fragment_length: i32, known_hotspots: Vec<Arc<KnownHotspot>>, contigs: Arc<Contigs>) -> Self {
        DiscordantGroups {
            region,
            max_concordant_fragment_length: observed_max_fragment_length,
            min_discordant_fragment_length: max_concordant_fragment_length(observed_max_fragment_length) * 2,
            known_hotspots,
            contigs,
            remote_outside_partition: std::cell::Cell::new(0),
            remote_other_chromosome: std::cell::Cell::new(0),
        }
    }

    /// isDiscordantGroup: by the first read that is not a supplementary
    pub fn is_discordant_group(&self, g: &ReadGroup, reads: &[PrepRead]) -> bool {
        match first_primary_read(g, reads) {
            Some(r) => is_discordant_fragment(&reads[r as usize].rec, self.min_discordant_fragment_length),
            None => false,
        }
    }

    /// KnownHotspot.readGroupMatchesHotspot
    fn read_group_matches_hotspot(&self, g: &ReadGroup, reads: &[PrepRead]) -> bool {
        for h in &self.known_hotspots {
            let (mut matches_start, mut matches_end) = (false, false);
            for &i in &g.reads {
                let r = &reads[i as usize];
                let mate_pos = mate_start(&r.rec);
                matches_start |= h.region_start.overlaps(r.chr, r.start, r.end);
                matches_start |= h.region_start.contains(r.mate_chr, mate_pos);
                matches_end |= h.region_end.overlaps(r.chr, r.start, r.end);
                matches_end |= h.region_end.contains(r.mate_chr, mate_pos);
                if !matches_start || !matches_end {
                    continue;
                }
                if (r.orientation() == r.mate_orientation()) == h.same_orientation() {
                    return true;
                }
            }
        }
        false
    }

    /// isRelevantDiscordantGroup: in a known hotspot pair, or a local deletion or duplication within range
    pub fn is_relevant_discordant_group(&self, g: &ReadGroup, reads: &[PrepRead]) -> bool {
        let (mut has_non_supp, mut has_passing_map_qual) = (false, false);
        for &i in &g.reads {
            let r = &reads[i as usize];
            has_non_supp |= !r.is_supplementary();
            has_passing_map_qual |= r.map_quality() >= MIN_MAP_QUALITY;
        }
        if !has_non_supp || !has_passing_map_qual {
            return false;
        }
        let Some(first) = first_primary_read(g, reads) else { return false };
        let first = &reads[first as usize];
        if self.read_group_matches_hotspot(g, reads) {
            return true;
        }
        if first.chr != first.mate_chr || first.orientation() == first.mate_orientation() {
            return false;
        }
        (first.start - mate_start(&first.rec)).wrapping_abs() <= DISCORDANT_GROUP_MAX_LOCAL_LENGTH
    }

    fn has_remote_required_fragments(&self, dg: &DiscordantGroup, remote: &RemoteRegion) -> bool {
        let mut is_short_local = false;
        if dg.region.chr == remote.chr {
            let min_distance = if remote.start > dg.region.end { remote.start - dg.region.end } else { dg.region.start - remote.end };
            is_short_local = min_distance <= self.min_discordant_fragment_length;
        }
        remote.read_groups.len() >= if is_short_local { DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT } else { DISCORDANT_GROUP_MIN_FRAGMENTS }
    }

    fn is_valid_discordant_group(&self, dg: &mut DiscordantGroup, groups: &[ReadGroup], reads: &[PrepRead], contigs: &Contigs) -> bool {
        if dg.read_groups.len() < DISCORDANT_GROUP_MIN_FRAGMENTS {
            return false;
        }
        if !self.region.overlaps_region(&dg.region) {
            return false;
        }
        merge_remote_regions(&mut dg.remote_regions, contigs);
        dg.remote_regions.sort_by_key(|x| std::cmp::Reverse(x.read_groups.len()));
        if dg.remote_regions[0].read_groups.len() < DISCORDANT_GROUP_MIN_FRAGMENTS {
            return false;
        }
        let (mut above_min_map_qual, mut above_min_align_score) = (false, false);
        for &gid in &dg.read_groups {
            let read = &reads[first_primary_read(&groups[gid as usize], reads).unwrap() as usize];
            above_min_map_qual |= read.map_quality() >= DISCORDANT_GROUP_MIN_MAP_QUAL;
            above_min_align_score |= above_repeat_trimmed_alignment_threshold(read, DISCORDANT_GROUP_MIN_ALIGN_SCORE, false);
            if above_min_align_score && above_min_map_qual {
                break;
            }
        }
        if !above_min_align_score || !above_min_map_qual {
            return false;
        }
        dg.remote_regions.iter().any(|r| self.has_remote_required_fragments(dg, r))
    }

    /// Corrected (finding P13): whether the group's fragments reach the share of the depth at its junction
    /// that the README and `junctionAboveMinAF` set for discordant sites, 1% where the depth is 100 or more.
    /// The support counted is that of the junction `add_junctions` would make: its supporting groups.
    fn above_min_af(&self, dg: &DiscordantGroup, groups: &[ReadGroup], reads: &[PrepRead], depth: &DepthTracker) -> bool {
        let (_, inner_position) = Self::inner_read(dg, groups, reads);
        let site_depth = depth.calc_depth(inner_position) as i32;
        if site_depth < DEPTH_MIN_CHECK {
            return true;
        }
        let range = self.supporting_read_range(dg.orient, inner_position);
        let support = dg.read_groups.iter().filter(|&&gid| Self::in_supporting_range(dg.orient, range, &groups[gid as usize], reads)).count();
        support as f64 >= site_depth as f64 * DEPTH_MIN_SUPPORT_RATIO_DISCORDANT
    }

    /// The innermost read of the remote region with the most support and its position, the junction's.
    fn inner_read(dg: &DiscordantGroup, groups: &[ReadGroup], reads: &[PrepRead]) -> (Option<ReadId>, i32) {
        let mut inner_read: Option<ReadId> = None;
        let mut inner_position = 0;
        for &gid in &dg.remote_regions[0].read_groups {
            let rid = first_primary_read(&groups[gid as usize], reads).unwrap();
            let read = &reads[rid as usize];
            if dg.orient.is_forward() {
                if inner_read.is_none_or(|r| read.end > reads[r as usize].end) {
                    inner_read = Some(rid);
                    inner_position = read.end;
                }
            } else if inner_read.is_none_or(|r| read.start < reads[r as usize].start) {
                inner_read = Some(rid);
                inner_position = read.start;
            }
        }
        (inner_read, inner_position)
    }

    /// The positions within the concordant fragment length of the junction, on its side.
    fn supporting_read_range(&self, orient: Orientation, inner_position: i32) -> (i32, i32) {
        if orient.is_forward() {
            (1.max(inner_position - self.max_concordant_fragment_length), inner_position)
        } else {
            (inner_position, inner_position + self.max_concordant_fragment_length)
        }
    }

    /// Whether a read of the group starts (forward junction) or ends (reverse) within the range.
    fn in_supporting_range(orient: Orientation, (start, end): (i32, i32), g: &ReadGroup, reads: &[PrepRead]) -> bool {
        g.reads.iter().any(|&i| {
            let r = &reads[i as usize];
            position_within(if orient.is_forward() { r.start } else { r.end }, start, end)
        })
    }

    /// formDiscordantJunctions: groups of discordant fragments that start close together and share a remote
    /// region make a junction at their innermost read, and one at the remote region. New junctions are
    /// appended to `jstore`. With `depth` (finding P13 corrected) a group whose fragments are too few for the
    /// depth at its junction makes none.
    pub fn form_discordant_junctions(
        &self,
        candidates: &[GroupId],
        groups: &mut [ReadGroup],
        reads: &mut [PrepRead],
        jstore: &mut Vec<JunctionData>,
        depth: Option<&DepthTracker>,
    ) {
        let contigs: &Contigs = &self.contigs;
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let mut read_groups: Vec<GroupId> = candidates
                .iter()
                .copied()
                .filter(|&g| first_primary_read(&groups[g as usize], reads).is_some_and(|r| reads[r as usize].orientation() == orientation))
                .collect();
            let mut assigned = vec![false; groups.len()];
            read_groups.sort_by_key(|&g| reads[first_primary_read(&groups[g as usize], reads).unwrap() as usize].start);
            let mut i = 0usize;
            while i + DISCORDANT_GROUP_MIN_FRAGMENTS <= read_groups.len() {
                let first_group = read_groups[i];
                if assigned[first_group as usize] {
                    i += 1;
                    continue;
                }
                let mut dg = DiscordantGroup::from_read_group(first_group, groups, reads);
                let mut j = i + 1;
                while j < read_groups.len() {
                    let next_group = read_groups[j];
                    if assigned[next_group as usize] {
                        j += 1;
                        continue;
                    }
                    let next_read = &reads[first_primary_read(&groups[next_group as usize], reads).unwrap() as usize];
                    if next_read.start > dg.region.end {
                        break;
                    }
                    dg.try_add_read_group(next_group, next_read);
                    j += 1;
                }
                if self.is_valid_discordant_group(&mut dg, groups, reads, contigs) && depth.is_none_or(|d| self.above_min_af(&dg, groups, reads, d)) {
                    self.add_junctions(&mut dg, groups, reads, jstore);
                    for &g in &dg.read_groups {
                        assigned[g as usize] = true;
                    }
                }
                i = j;
            }
        }
    }

    fn add_junctions(&self, dg: &mut DiscordantGroup, groups: &mut [ReadGroup], reads: &mut [PrepRead], jstore: &mut Vec<JunctionData>) {
        // the junction is at the innermost read of the remote region with the most support
        let (inner_read, inner_position) = Self::inner_read(dg, groups, reads);
        let mut junction = JunctionData::new(inner_position, dg.orient, inner_read);
        dg.innermost_read = inner_read;
        junction.discordant_group = true;

        let range = self.supporting_read_range(dg.orient, inner_position);
        let mut excluded: Vec<GroupId> = Vec::new();
        for &gid in &dg.read_groups {
            // only reads within the concordant fragment length of the junction are of its group
            if !Self::in_supporting_range(dg.orient, range, &groups[gid as usize], reads) {
                if !excluded.contains(&gid) {
                    excluded.push(gid);
                }
                continue;
            }
            junction.supporting_groups.push(gid);
            groups[gid as usize].add_junction_position(junction.position, junction.orient);
            for &i in &groups[gid as usize].reads {
                reads[i as usize].set_read_type_ranked(ReadType::Support);
            }
        }
        dg.purge_reads(&mut excluded);
        // every remaining group is kept, even if its remote reads may not be
        for &gid in &dg.read_groups {
            groups[gid as usize].has_remote_junction_reads = true;
        }

        // the main remote region is the first of the list, by position in memory in Java
        let mut remote_junctions: Vec<JunctionData> = Vec::new();
        for (k, remote) in dg.remote_regions.iter().enumerate() {
            let is_main = k == 0;
            if is_main {
                junction.add_remote_junction(RemoteJunction { chr: remote.chr, position: remote.start, orient: Orientation::Forward, fragments: remote.read_groups.len() as i32 });
            }
            if !self.has_remote_required_fragments(dg, remote) {
                continue;
            }
            let remote_position = if remote.orient.is_forward() { remote.end } else { remote.start };
            if !self.region.contains(remote.chr, remote_position) {
                self.remote_outside_partition.set(self.remote_outside_partition.get() + 1);
                if remote.chr != self.region.chr {
                    self.remote_other_chromosome.set(self.remote_other_chromosome.get() + 1);
                }
            }
            let first_read = groups[remote.read_groups[0] as usize].reads[0];
            let mut remote_junction = JunctionData::new(remote_position, remote.orient, Some(first_read));
            remote_junction.discordant_group = true;
            if !self.region.contains(remote.chr, remote_position) {
                remote_junction.outside_chr = Some(remote.chr);
            }
            for &gid in &remote.read_groups {
                remote_junction.supporting_groups.push(gid);
                groups[gid as usize].add_junction_position(remote_junction.position, remote_junction.orient);
                for &i in &groups[gid as usize].reads {
                    reads[i as usize].set_read_type_ranked(ReadType::Support);
                }
            }
            if is_main {
                remote_junction.add_remote_junction(RemoteJunction {
                    chr: dg.region.chr,
                    position: dg.inner_position(reads),
                    orient: dg.orient,
                    fragments: dg.read_groups.len() as i32,
                });
            }
            remote_junctions.push(remote_junction);
        }
        jstore.push(junction);
        jstore.extend(remote_junctions);
    }
}
