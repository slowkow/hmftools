//! SBX duplicate collapsing (Java duplicate/SbxDuplicateCollapser): fragments whose two unclipped ends are within
//! a few bases of each other are duplicates of each other, since the ends of SBX reads jitter.

use crate::dup::{DuplicateGroup, FragmentCoordReads, FragmentCoords, ReadInfo};
use crate::rec::*;
use rust_htslib::bam;

const MAX_CHAINED_MULTIPLE: i32 = 3;

fn within_range(first_pos_lower: i32, first_pos_upper: i32, second_pos_lower: i32, second_pos_upper: i32, max_distance: i32) -> bool {
    (first_pos_lower - second_pos_lower).abs() + (first_pos_upper - second_pos_upper).abs() <= max_distance
}

#[derive(Clone, Copy)]
enum Source {
    Single(usize),
    Group(usize),
}

/// SbxDuplicateCollapser.GroupInfo: a single read or a duplicate group, with the range of positions of all that
/// has been collapsed into it
struct GroupInfo {
    source: Source,
    position_lower: i32,
    position_upper: i32,
    position_lower_min: i32,
    position_lower_max: i32,
    position_upper_min: i32,
    position_upper_max: i32,
    size: i32,
    collapsed: bool,
    collapsed_groups: Vec<usize>,
}

fn min_distance(first_min: i32, first_max: i32, second_min: i32, second_max: i32) -> i32 {
    if first_max < second_min {
        second_min - first_max
    } else if second_max < first_min {
        first_min - second_max
    } else {
        0
    }
}

pub struct SbxDuplicateCollapser {
    max_duplicate_distance: i32,
}

impl SbxDuplicateCollapser {
    pub fn new(max_duplicate_distance: i32) -> Self {
        SbxDuplicateCollapser { max_duplicate_distance }
    }

    fn distance(a: &GroupInfo, b: &GroupInfo) -> i32 {
        min_distance(a.position_lower_min, a.position_lower_max, b.position_lower_min, b.position_lower_max)
            + min_distance(a.position_upper_min, a.position_upper_max, b.position_upper_min, b.position_upper_max)
    }

    /// GroupInfo.addGroup: the position ranges widen, by no more than three times the distance from the original
    fn add_group(&self, infos: &mut [GroupInfo], main: usize, other: usize) {
        let (o_size, o_lmin, o_lmax, o_umin, o_umax) = {
            let o = &infos[other];
            (o.size, o.position_lower_min, o.position_lower_max, o.position_upper_min, o.position_upper_max)
        };
        let g = &mut infos[main];
        g.collapsed_groups.push(other);
        g.size += o_size;
        g.position_lower_min = g.position_lower_min.min(o_lmin);
        g.position_lower_max = g.position_lower_max.max(o_lmax);
        g.position_upper_min = g.position_upper_min.min(o_umin);
        g.position_upper_max = g.position_upper_max.max(o_umax);
        let max_distance_from_original = MAX_CHAINED_MULTIPLE * self.max_duplicate_distance;
        g.position_lower_min = g.position_lower_min.max(g.position_lower - max_distance_from_original);
        g.position_lower_max = g.position_lower_max.min(g.position_lower + max_distance_from_original);
        g.position_upper_min = g.position_upper_min.max(g.position_upper - max_distance_from_original);
        g.position_upper_max = g.position_upper_max.min(g.position_upper + max_distance_from_original);
    }

    /// SbxDuplicateCollapser.findProximateGroups: the groups in range of this one that are not collapsed and not
    /// larger. `groups` is sorted by position and `index` is this group's place in it.
    fn find_proximate_groups(&self, infos: &[GroupInfo], groups: &[usize], index: usize, out: &mut Vec<usize>) -> bool {
        out.clear();
        let mut found = false;
        let g = &infos[groups[index]];
        // a position range reaches at most this far from the group's own lower position
        let max_reach = (MAX_CHAINED_MULTIPLE + 1) * self.max_duplicate_distance;
        for search_up in [true, false] {
            let mut j = if search_up { index as i64 + 1 } else { index as i64 - 1 };
            while j >= 0 && (j as usize) < groups.len() {
                let next = &infos[groups[j as usize]];
                if Self::distance(g, next) > self.max_duplicate_distance {
                    if next.position_lower > g.position_lower + self.max_duplicate_distance {
                        break;
                    }
                    // Java walks down to the start of the list; no group below this point can be in range
                    if !search_up && next.position_lower + max_reach < g.position_lower_min {
                        break;
                    }
                } else if !next.collapsed && next.size <= g.size {
                    out.push(groups[j as usize]);
                    found = true;
                }
                j += if search_up { 1 } else { -1 };
            }
        }
        found
    }

    /// SbxDuplicateCollapser.findGroupsToCollapse
    fn find_groups_to_collapse(&self, infos: &[GroupInfo], groups: &[usize], scratch: &mut Vec<usize>) -> Vec<usize> {
        let mut last_collapsing: Option<usize> = None;
        let mut to_collapse = Vec::new();
        for (i, &gi) in groups.iter().enumerate() {
            let g = &infos[gi];
            if g.collapsed {
                continue;
            }
            if let Some(l) = last_collapsing {
                // a group close to the last collapsing one is not checked unless it is larger
                let last = &infos[l];
                if g.size <= last.size
                    && within_range(last.position_lower, last.position_upper, g.position_lower, g.position_upper, self.max_duplicate_distance * 2)
                {
                    continue;
                }
            }
            if !self.find_proximate_groups(infos, groups, i, scratch) {
                continue;
            }
            to_collapse.push(i);
            last_collapsing = Some(gi);
        }
        to_collapse
    }

    /// SbxDuplicateCollapser.collapseGroups(groups): the largest groups take in their neighbours until nothing
    /// more collapses
    fn collapse_sorted_groups(&self, infos: &mut [GroupInfo], groups: &[usize]) {
        let mut scratch = Vec::new();
        loop {
            let mut to_collapse = self.find_groups_to_collapse(infos, groups, &mut scratch);
            if to_collapse.is_empty() {
                break;
            }
            // largest first; the sort is stable
            to_collapse.sort_by_key(|&i| -infos[groups[i]].size);
            for index in to_collapse {
                let main = groups[index];
                if infos[main].collapsed {
                    continue;
                }
                if !self.find_proximate_groups(infos, groups, index, &mut scratch) {
                    continue;
                }
                for &other in scratch.iter() {
                    if infos[other].collapsed {
                        continue;
                    }
                    self.add_group(infos, main, other);
                    infos[other].collapsed = true;
                }
            }
        }
    }

    /// SbxDuplicateCollapser.collapseGroups(duplicateGroups, singleReads)
    pub fn collapse_groups(&self, fragments: FragmentCoordReads, names: &[String]) -> FragmentCoordReads {
        let mut single_reads: Vec<Option<ReadInfo>> = fragments.single_reads.into_iter().map(Some).collect();
        let mut duplicate_groups: Vec<Option<DuplicateGroup>> = fragments.duplicate_groups.into_iter().map(Some).collect();

        // the lists by orientation of the lower position and primary or supplementary, in the order a
        // java.util.HashMap walks the keys "R", "F", "R_S" and "F_S"
        // R7: Java takes the orientation of the lower position, which for an unpaired read is forward on both
        // strands (the lower position of a reverse read is its 3' end), so forward and reverse reads collapse
        // together. Corrected: the orientation of the read, as in the fragment key.
        let java_orientation = crate::compat::java("R7");
        let key_index = |c: &FragmentCoords| -> usize {
            let forward = if c.unpaired && !java_orientation && !c.read_is_lower { c.upper_forward } else { c.lower_forward };
            match (forward, c.supp.is_some()) {
                (false, false) => 0,
                (true, false) => 1,
                (false, true) => 2,
                (true, true) => 3,
            }
        };
        let mut infos: Vec<GroupInfo> = Vec::with_capacity(single_reads.len() + duplicate_groups.len());
        let mut key_groups: [Vec<usize>; 4] = Default::default();
        let mut add = |infos: &mut Vec<GroupInfo>, source: Source, coords: &FragmentCoords, size: i32| {
            key_groups[key_index(coords)].push(infos.len());
            infos.push(GroupInfo {
                source,
                position_lower: coords.pos_lower,
                position_upper: coords.pos_upper,
                position_lower_min: coords.pos_lower,
                position_lower_max: coords.pos_lower,
                position_upper_min: coords.pos_upper,
                position_upper_max: coords.pos_upper,
                size,
                collapsed: false,
                collapsed_groups: Vec::new(),
            });
        };
        for (i, r) in single_reads.iter().enumerate() {
            add(&mut infos, Source::Single(i), &r.as_ref().unwrap().coords, 1);
        }
        for (i, g) in duplicate_groups.iter().enumerate() {
            let g = g.as_ref().unwrap();
            add(&mut infos, Source::Group(i), &g.coords, g.reads.len() as i32);
        }

        let mut final_duplicate_groups = Vec::new();
        let mut final_single_reads = Vec::new();
        for groups in key_groups.iter_mut() {
            if groups.len() > 1 {
                // by lower then upper position; stable
                groups.sort_by(|&a, &b| {
                    (infos[a].position_lower, infos[a].position_upper).cmp(&(infos[b].position_lower, infos[b].position_upper))
                });
                self.collapse_sorted_groups(&mut infos, groups);
            }
            for &gi in groups.iter() {
                if infos[gi].collapsed {
                    continue;
                }
                let has_collapsed_groups = !infos[gi].collapsed_groups.is_empty();
                let mut duplicate_group = match infos[gi].source {
                    Source::Single(i) => {
                        let read_info = single_reads[i].take().unwrap();
                        if !has_collapsed_groups {
                            final_single_reads.push(read_info);
                            continue;
                        }
                        DuplicateGroup::new(vec![read_info.read], read_info.coords)
                    }
                    Source::Group(i) => duplicate_groups[i].take().unwrap(),
                };
                if has_collapsed_groups {
                    // the reads of the collapsed groups count towards the consensus
                    let collapsed = std::mem::take(&mut infos[gi].collapsed_groups);
                    for other in collapsed {
                        take_all_reads(&mut infos, other, &mut single_reads, &mut duplicate_groups, &mut duplicate_group.reads);
                    }
                }
                let non_consensus = check_supplementary_consistency(&duplicate_group.reads, self.max_duplicate_distance, names);
                if let Some(indices) = &non_consensus {
                    crate::sbx::eager_decode_on_remove(&mut duplicate_group.reads, indices);
                }
                move_non_consensus_reads(&mut duplicate_group, non_consensus);
                final_duplicate_groups.push(duplicate_group);
            }
        }
        FragmentCoordReads { duplicate_groups: final_duplicate_groups, single_reads: final_single_reads }
    }
}

/// GroupInfo.allReads: the group's own reads, then those of the groups collapsed into it
fn take_all_reads(
    infos: &mut [GroupInfo], gi: usize, single_reads: &mut [Option<ReadInfo>], duplicate_groups: &mut [Option<DuplicateGroup>],
    out: &mut Vec<bam::Record>,
) {
    match infos[gi].source {
        Source::Single(i) => out.push(single_reads[i].take().unwrap().read),
        Source::Group(i) => out.append(&mut duplicate_groups[i].take().unwrap().reads),
    }
    let collapsed = std::mem::take(&mut infos[gi].collapsed_groups);
    for other in collapsed {
        take_all_reads(infos, other, single_reads, duplicate_groups, out);
    }
}

fn move_non_consensus_reads(group: &mut DuplicateGroup, non_consensus: Option<Vec<usize>>) {
    let Some(indices) = non_consensus else { return };
    let mut i = 0;
    let reads = std::mem::take(&mut group.reads);
    for read in reads {
        if indices.contains(&i) {
            group.non_consensus_reads.push(read);
        } else {
            group.reads.push(read);
        }
        i += 1;
    }
}

/// SbxDuplicateCollapser.checkSupplementaryConsistency: with supplementary reads in the group, the indices of
/// the reads whose own unclipped ends are not in range of those of the template (the first primary read, or the
/// first read). None when the group has no supplementary read or every read is in range.
pub fn check_supplementary_consistency<R: std::borrow::Borrow<bam::Record>>(
    reads: &[R], max_dup_distance: i32, names: &[String],
) -> Option<Vec<usize>> {
    let supp_count = reads.iter().filter(|r| has((*r).borrow(), FLAG_SUPPLEMENTARY)).count();
    if supp_count == 0 {
        return None;
    }
    let template =
        if supp_count == reads.len() { 0 } else { reads.iter().position(|r| !has(r.borrow(), FLAG_SUPPLEMENTARY)).unwrap_or(0) };
    let template_coords = FragmentCoords::from_read_as(reads[template].borrow(), names, true, false);
    let mut non_consensus: Option<Vec<usize>> = None;
    for (i, read) in reads.iter().enumerate() {
        if i == template {
            continue;
        }
        let coords = FragmentCoords::from_read_as(read.borrow(), names, true, false);
        if !within_range(template_coords.pos_lower, template_coords.pos_upper, coords.pos_lower, coords.pos_upper, max_dup_distance) {
            non_consensus.get_or_insert_with(Vec::new).push(i);
        }
    }
    non_consensus
}
