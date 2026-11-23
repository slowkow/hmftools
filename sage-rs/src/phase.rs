//! Local phasing: CandidateVariantPhaser, PhasedGroupCollection, PhasedVariantGroup, PhasingUtils.
//! Read counters are identified by their index in the tumor counter list; `pos` maps index to position.

use crate::constants::{is_sbx, is_ultima, jround};
use crate::evidence::ReadContextCounter;
use crate::evidence_stage::PhaseSink;
use crate::filters::SageVariant;
use crate::types::*;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct PhasedGroup {
    pub id: i32,
    pub positive: Vec<usize>,
    pub negative: Vec<usize>,
    pub read_count: i32,
    pub allocated: f64,
    pub(crate) pos_min: i32,
    pos_max: i32,
    var_min: i32,
    var_max: i32,
}

fn first_pos(c: &[usize], pos: &[i32]) -> i32 {
    pos[c[0]]
}
fn last_pos(c: &[usize], pos: &[i32]) -> i32 {
    pos[c[c.len() - 1]]
}
fn any_overlap(a: &[usize], b: &[usize]) -> bool {
    a.iter().any(|x| b.contains(x))
}

impl PhasedGroup {
    fn new(id: i32, pos_min: i32, pos_max: i32, positive: Vec<usize>, negative: Vec<usize>, pos: &[i32]) -> Self {
        let (mut var_min, mut var_max) = (pos_min, pos_max);
        if !negative.is_empty() {
            var_min = var_min.min(first_pos(&negative, pos));
            var_max = var_max.max(last_pos(&negative, pos));
        }
        PhasedGroup { id, positive, negative, read_count: 1, allocated: 0.0, pos_min, pos_max, var_min, var_max }
    }
    fn positions_overlap(&self, o: &PhasedGroup) -> bool {
        positions_overlap(self.pos_min, self.pos_max, o.pos_min, o.pos_max)
    }
    fn is_valid(&self) -> bool {
        !self.positive.is_empty() && self.positive.len() + self.negative.len() >= 2
    }
    fn exact_match(&self, pos_min: i32, pos_max: i32, positive: &[usize], negative: &[usize]) -> bool {
        pos_min == self.pos_min
            && pos_max == self.pos_max
            && self.positive.len() == positive.len()
            && self.negative.len() == negative.len()
            && self.positive.iter().all(|x| positive.contains(x))
            && self.negative.iter().all(|x| negative.contains(x))
    }
    fn positives_match(&self, o: &PhasedGroup) -> bool {
        o.pos_min == self.pos_min
            && o.pos_max == self.pos_max
            && self.positive.len() == o.positive.len()
            && self.positive.iter().all(|x| o.positive.contains(x))
    }
    fn is_subset_of(&self, o: &PhasedGroup) -> bool {
        if o.positive.len() < self.positive.len() {
            return false;
        }
        if !self.positive.iter().all(|x| o.positive.contains(x)) {
            return false;
        }
        !(any_overlap(&o.positive, &self.negative) || any_overlap(&self.positive, &o.negative))
    }
    fn populate_common(&self, o: &PhasedGroup, common_pos: &mut Vec<usize>, common_neg: &mut Vec<usize>) -> bool {
        if any_overlap(&o.positive, &self.negative) || any_overlap(&self.positive, &o.negative) {
            return false;
        }
        common_pos.extend(self.positive.iter().filter(|x| o.positive.contains(x)));
        common_neg.extend(self.negative.iter().filter(|x| o.negative.contains(x)));
        !common_pos.is_empty()
    }
    fn has_common_subset(&self, o: &PhasedGroup, common_pos: &[usize], common_neg: &[usize]) -> bool {
        if any_overlap(&o.positive, &self.negative) || any_overlap(&self.positive, &o.negative) {
            return false;
        }
        common_pos.iter().all(|x| self.positive.contains(x)) && common_neg.iter().all(|x| self.negative.contains(x))
    }
    fn merge_negatives(&mut self, negative: &[usize], pos: &[i32]) {
        for x in negative {
            if !self.negative.contains(x) {
                self.negative.push(*x);
            }
        }
        if !self.negative.is_empty() {
            let min = self.negative.iter().map(|&i| pos[i]).min().unwrap();
            let max = self.negative.iter().map(|&i| pos[i]).max().unwrap();
            self.var_min = self.pos_min.min(min);
            self.var_max = self.pos_max.max(max);
        }
    }
    fn merge(&mut self, o: &PhasedGroup, pos: &[i32]) {
        if o.positive.len() == self.positive.len() {
            self.read_count += o.read_count;
        } else {
            self.allocated += o.read_count as f64;
        }
        self.allocated += o.allocated;
        self.merge_negatives(&o.negative, pos);
    }
    fn merge_fraction(&mut self, o: &PhasedGroup, alloc_fraction: f64, pos: &[i32]) {
        self.allocated += (o.read_count as f64 + o.allocated) * alloc_fraction;
        for &rc in &o.positive {
            let mut matched = false;
            let mut index = 0;
            while index < self.positive.len() {
                let c = self.positive[index];
                if c == rc {
                    matched = true;
                    break;
                }
                if pos[c] <= pos[rc] {
                    index += 1;
                } else {
                    break;
                }
            }
            if !matched {
                self.positive.insert(index, rc);
            }
        }
        self.pos_min = first_pos(&self.positive, pos);
        self.pos_max = last_pos(&self.positive, pos);
        self.merge_negatives(&o.negative, pos);
    }
    /// cullReadCounters
    fn cull(&mut self, valid: &[bool], pos: &[i32]) -> bool {
        let localised = |i: usize| valid[i] && position_within(pos[i], self.var_min, self.var_max);
        if !(0..valid.len()).any(localised) {
            return true;
        }
        let invalid_pos: Vec<usize> = self.positive.iter().copied().filter(|&i| !localised(i)).collect();
        let invalid_neg: Vec<usize> = self.negative.iter().copied().filter(|&i| !localised(i)).collect();
        if invalid_pos.is_empty() && invalid_neg.is_empty() {
            return false;
        }
        if invalid_pos.len() == self.positive.len() {
            return true;
        }
        self.positive.retain(|i| !invalid_pos.contains(i));
        self.negative.retain(|i| !invalid_neg.contains(i));
        self.pos_min = first_pos(&self.positive, pos);
        self.pos_max = last_pos(&self.positive, pos);
        self.var_min = self.pos_min;
        self.var_max = self.pos_max;
        if !self.negative.is_empty() {
            self.var_min = self.var_min.min(first_pos(&self.negative, pos));
            self.var_max = self.var_max.max(last_pos(&self.negative, pos));
        }
        true
    }
}

pub(crate) struct GroupCollection {
    /// keyed by the minimum positive position; Java uses a HashMap but only ever sorts by key
    pub(crate) map: Vec<(i32, Vec<PhasedGroup>)>,
    pub(crate) groups: Vec<PhasedGroup>,
    min_position: i32,
    max_position: i32,
    group_count: usize,
    java_bounds: bool,
    java_merge: bool,
}

impl GroupCollection {
    fn new(java_bounds: bool, java_merge: bool) -> Self {
        GroupCollection { map: Vec::new(), groups: Vec::new(), min_position: 0, max_position: 0, group_count: 0, java_bounds, java_merge }
    }
    fn overlaps(&self, min: i32, max: i32) -> bool {
        positions_overlap(self.min_position, self.max_position, min, max)
    }
    fn add(&mut self, pos_min: i32, pos_max: i32, id: i32, positive: &[usize], negative: &[usize], pos: &[i32]) -> bool {
        let idx = match self.map.iter().position(|(k, _)| *k == pos_min) {
            Some(i) => i,
            None => {
                self.map.push((pos_min, Vec::new()));
                self.map.len() - 1
            }
        };
        let groups = &mut self.map[idx].1;
        for g in groups.iter_mut() {
            if g.exact_match(pos_min, pos_max, positive, negative) {
                g.read_count += 1;
                g.merge_negatives(negative, pos);
                return false;
            }
        }
        let g = PhasedGroup::new(id, pos_min, pos_max, positive.to_vec(), negative.to_vec(), pos);
        // S2, Java: the emptiness test is on this key's list, not on the collection, so the first group
        // of every later key replaces the collection's minimum instead of lowering it. Corrected: only
        // the first group of the collection sets it.
        let first = if self.java_bounds { groups.is_empty() } else { self.group_count == 0 };
        self.min_position = if first { g.var_min } else { self.min_position.min(g.var_min) };
        self.max_position = self.max_position.max(g.var_max);
        groups.push(g);
        self.group_count += 1;
        true
    }
    fn merge(&mut self, other: GroupCollection, pos: &[i32]) {
        self.group_count += other.group_count;
        for (k, list) in other.map {
            match self.map.iter_mut().find(|(key, _)| *key == k) {
                // S1, Java: putAll replaces this collection's list for any key the other also has, and
                // those groups are lost. Corrected: the other's groups join the list, and one that
                // matches a group already there exactly is counted into it.
                Some(e) if self.java_merge => e.1 = list,
                Some(e) => {
                    for g in list {
                        match e.1.iter_mut().find(|x| x.exact_match(g.pos_min, g.pos_max, &g.positive, &g.negative)) {
                            Some(x) => {
                                x.read_count += g.read_count;
                                x.merge_negatives(&g.negative, pos);
                            }
                            None => e.1.push(g),
                        }
                    }
                }
                None => self.map.push((k, list)),
            }
        }
        self.min_position = self.min_position.min(other.min_position);
        self.max_position = self.max_position.max(other.max_position);
    }
    pub(crate) fn finalise(&mut self) {
        let mut map = std::mem::take(&mut self.map);
        map.sort_by_key(|(k, _)| *k);
        for (_, list) in map {
            self.groups.extend(list);
        }
    }
}

pub struct Phaser {
    pub(crate) collections: Vec<GroupCollection>,
    next_group_id: i32,
    pub(crate) pos: Vec<i32>,
    /// compat::java("S2") and ("S1"), read once per region
    java_bounds: bool,
    java_merge: bool,
}

impl PhaseSink for Phaser {
    fn register(&mut self, counters: &[ReadContextCounter], positive: &[usize], negative: &[usize]) {
        if positive.is_empty() || positive.len() + negative.len() < 2 {
            return;
        }
        if self.pos.is_empty() {
            self.pos = counters.iter().map(|c| c.position()).collect();
        }
        let pos = &self.pos;
        let pos_min = first_pos(positive, pos);
        let pos_max = last_pos(positive, pos);
        let var_min = if negative.is_empty() { pos_min } else { pos_min.min(first_pos(negative, pos)) };
        let var_max = if negative.is_empty() { pos_max } else { pos_max.max(last_pos(negative, pos)) };
        let mut index = 0;
        while index < self.collections.len() {
            if self.collections[index].overlaps(var_min, var_max) {
                if self.collections[index].add(pos_min, pos_max, self.next_group_id, positive, negative, pos) {
                    self.next_group_id += 1;
                    if index + 1 < self.collections.len() {
                        let (nmin, nmax) = (self.collections[index + 1].min_position, self.collections[index + 1].max_position);
                        if self.collections[index].overlaps(nmin, nmax) {
                            let next = self.collections.remove(index + 1);
                            self.collections[index].merge(next, pos);
                        }
                    }
                }
                return;
            } else if var_max < self.collections[index].min_position {
                break;
            }
            index += 1;
        }
        let mut c = GroupCollection::new(self.java_bounds, self.java_merge);
        c.add(pos_min, pos_max, self.next_group_id, positive, negative, pos);
        self.next_group_id += 1;
        self.collections.push(c);
    }
}

pub(crate) fn merge_matching(groups: &mut Vec<PhasedGroup>, pos: &[i32]) {
    if groups.len() < 2 {
        return;
    }
    let mut i = 0;
    while i < groups.len() {
        let supers: Vec<usize> =
            (i + 1..groups.len()).filter(|&j| groups[i].positions_overlap(&groups[j]) && groups[i].is_subset_of(&groups[j])).collect();
        if supers.is_empty() {
            if groups[i].positive.len() == 1 {
                groups.remove(i);
                continue;
            }
        } else if !supers.iter().any(|&j| groups[j].positive.len() > groups[i].positive.len()) {
            for &j in &supers {
                let other = groups[j].clone();
                groups[i].merge(&other, pos);
            }
            for &j in supers.iter().rev() {
                groups.remove(j);
            }
        }
        i += 1;
    }
}

pub(crate) fn merge_by_extension(groups: &mut Vec<PhasedGroup>, pos: &[i32]) {
    if groups.len() < 2 {
        return;
    }
    let mut modified: Vec<i32> = Vec::new();
    let mut initial = true;
    while initial || !modified.is_empty() {
        let last_modified = std::mem::take(&mut modified);
        let mut i = 0;
        while i < groups.len() {
            if !initial && !last_modified.contains(&groups[i].id) {
                i += 1;
                continue;
            }
            let mut candidates: Vec<usize> = vec![i];
            let mut common_pos: Vec<usize> = Vec::new();
            let mut common_neg: Vec<usize> = Vec::new();
            for j in 0..groups.len() {
                if j == i || !groups[i].positions_overlap(&groups[j]) {
                    continue;
                }
                if common_pos.is_empty() {
                    if !groups[i].populate_common(&groups[j], &mut common_pos, &mut common_neg) {
                        continue;
                    }
                    candidates.push(j);
                } else if candidates.iter().all(|&x| groups[j].has_common_subset(&groups[x], &common_pos, &common_neg)) {
                    candidates.push(j);
                }
            }
            if candidates.len() == 1 {
                i += 1;
                continue;
            }
            common_pos.extend(common_neg.iter());
            let min_subset = common_pos.iter().map(|&x| pos[x]).min().unwrap_or(0);
            let max_subset = common_pos.iter().map(|&x| pos[x]).max().unwrap_or(0);
            let lower: Vec<usize> =
                candidates.iter().copied().filter(|&x| groups[x].var_min < min_subset && groups[x].var_max == max_subset).collect();
            let upper: Vec<usize> =
                candidates.iter().copied().filter(|&x| groups[x].var_min == min_subset && groups[x].var_max > max_subset).collect();
            if lower.len() == 1 || upper.len() == 1 {
                let merged = if lower.len() == 1 { lower[0] } else { upper[0] };
                let mut merging = if lower.len() == 1 { upper.clone() } else { lower.clone() };
                if merging.is_empty() {
                    merging = candidates.iter().copied().filter(|&x| x != merged).collect();
                }
                let total: f64 = merging.iter().map(|&x| groups[x].read_count).sum::<i32>() as f64;
                let merged_group = groups[merged].clone();
                for &x in &merging {
                    let fraction = groups[x].read_count as f64 / total;
                    groups[x].merge_fraction(&merged_group, fraction, pos);
                    if !modified.contains(&groups[x].id) {
                        modified.push(groups[x].id);
                    }
                }
                groups.remove(merged);
                break;
            }
            i += 1;
        }
        initial = false;
    }
}

pub(crate) fn merge_uninformative(groups: &mut Vec<PhasedGroup>, pos: &[i32]) {
    // S3, Java: the test for "a single variant not present as a +ve in another group" looks in the
    // group itself, where it is always found, so such a group is only removed when it is the only one;
    // and the loop index is not adjusted when groups ahead of it are removed, so the group that moves
    // into the current place is skipped. Corrected: the other groups are searched, and the index
    // follows the current group.
    let java = crate::compat::java("S3");
    let mut index = 0;
    while index < groups.len() {
        if groups[index].positive.len() == 1 {
            let counter = groups[index].positive[0];
            let alone = if java {
                groups.len() == 1
            } else {
                !groups.iter().enumerate().any(|(j, g)| j != index && g.positive.contains(&counter))
            };
            if alone {
                groups.remove(index);
                continue;
            }
        }
        let matching: Vec<usize> = (0..groups.len())
            .filter(|&j| j != index && (groups[index].positives_match(&groups[j]) || groups[j].is_subset_of(&groups[index])))
            .collect();
        if !matching.is_empty() {
            for &j in &matching {
                let other = groups[j].clone();
                groups[index].merge(&other, pos);
            }
            for &j in matching.iter().rev() {
                groups.remove(j);
            }
            if !java {
                index -= matching.iter().filter(|&&j| j < index).count();
            }
        }
        index += 1;
    }
}

impl Phaser {
    pub fn new() -> Self {
        Phaser {
            collections: Vec::new(),
            next_group_id: 0,
            pos: Vec::new(),
            java_bounds: crate::compat::java("S2"),
            java_merge: crate::compat::java("S1"),
        }
    }

    /// assignLocalPhaseSets: `passing` and `valid` are indexed like the tumor counters.
    pub fn assign_local_phase_sets(&mut self, passing: &[bool], valid: &[bool], next_lps: &mut i32) -> Vec<LpsAssignment> {
        let mut out = Vec::new();
        if !passing.iter().any(|x| *x) {
            return out;
        }
        let pos = self.pos.clone();
        let mut has_groups = false;
        for c in self.collections.iter_mut() {
            c.finalise();
            let mut threshold = (jround((c.groups.len() as f64).log10()) as i32 - 2).max(1);
            if is_ultima() || is_sbx() {
                // lower for the sequencing technologies with longer reads
                threshold = (threshold - 1).max(1);
            }
            c.groups.retain_mut(|g| {
                !(g.read_count < threshold
                    || !g.positive.iter().any(|&i| passing[i])
                    || (g.cull(valid, &pos) && !g.is_valid()))
            });
            has_groups |= !c.groups.is_empty();
        }
        if !has_groups {
            return out;
        }
        for c in self.collections.iter_mut() {
            if c.groups.is_empty() {
                continue;
            }
            merge_matching(&mut c.groups, &pos);
            merge_by_extension(&mut c.groups, &pos);
            merge_uninformative(&mut c.groups, &pos);
            merge_by_extension(&mut c.groups, &pos);
        }
        for c in &self.collections {
            for g in &c.groups {
                if (g.read_count as f64 + g.allocated) < 2.0 {
                    continue;
                }
                let lps = *next_lps;
                *next_lps += 1;
                out.push(LpsAssignment { lps, counters: g.positive.clone(), read_count: g.read_count, allocated: g.allocated });
            }
        }
        out
    }
}

/// A phase set id to add to each of the listed tumor counters.
pub struct LpsAssignment {
    pub lps: i32,
    pub counters: Vec<usize>,
    pub read_count: i32,
    pub allocated: f64,
}

/// Iteration order of a Java HashMap<Integer, _> / HashSet<Integer> built by inserting `keys` in order.
pub(crate) fn java_int_hash_order(keys: &[i32]) -> Vec<i32> {
    let mut unique: Vec<i32> = Vec::new();
    for k in keys {
        if !unique.contains(k) {
            unique.push(*k);
        }
    }
    let mut capacity = 16usize;
    while unique.len() as f64 > capacity as f64 * 0.75 {
        capacity *= 2;
    }
    let mut items: Vec<(usize, usize, i32)> = unique
        .iter()
        .enumerate()
        .map(|(seq, &k)| {
            let h = k as u32;
            (((h ^ (h >> 16)) as usize) & (capacity - 1), seq, k)
        })
        .collect();
    items.sort();
    items.into_iter().map(|(_, _, k)| k).collect()
}

/// The passing phase set ids in the iteration order of RegionTask.mPassingPhaseSets.
pub fn passing_phase_sets(variants: &[SageVariant]) -> Vec<i32> {
    let mut keys = Vec::new();
    for v in variants {
        if v.is_passing() && v.has_local_phase_sets() {
            keys.extend_from_slice(v.local_phase_sets());
        }
    }
    java_int_hash_order(&keys)
}

/// PhasingUtils.removeUninformativeLps. `variants` are the final variants; `passing_sets` is in
/// Java HashSet iteration order.
pub fn remove_uninformative_lps(variants: &mut [SageVariant], passing_sets: &[i32]) {
    const SUBSET_READ_COUNT_LIMIT: f64 = 0.25;
    // lps -> indices of variants with that lps, in the map's iteration order
    let order = java_int_hash_order(passing_sets);
    let mut lps_variants: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut lps_passing: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut uninformative: Vec<i32> = Vec::new();
    let mut processed: Vec<i32> = Vec::new();
    let mut single_passing: Vec<i32> = Vec::new();
    for &lps in passing_sets {
        let vs: Vec<usize> = (0..variants.len()).filter(|&i| variants[i].has_lps(lps)).collect();
        let passing: Vec<usize> = vs.iter().copied().filter(|&i| variants[i].is_passing()).collect();
        let single = vs.len() == 1 && !passing.is_empty() && variants[passing[0]].local_phase_sets().len() == 1;
        lps_variants.insert(lps, vs);
        if passing.is_empty() {
            uninformative.push(lps);
            processed.push(lps);
            continue;
        }
        if single && !single_passing.contains(&lps) {
            single_passing.push(lps);
        }
        lps_passing.insert(lps, passing);
    }
    for &lps in &order {
        if processed.contains(&lps) {
            continue;
        }
        processed.push(lps);
        let passing = &lps_passing[&lps];
        let mut max_read_count = variants[passing[0]].lps_read_count(lps);
        let mut max_lps = lps;
        let mut matched = vec![lps];
        for &other in &order {
            if other == lps || processed.contains(&other) {
                continue;
            }
            let other_passing = &lps_passing[&other];
            if other_passing.len() == passing.len() && other_passing.iter().all(|x| passing.contains(x)) {
                processed.push(other);
                let other_count = variants[other_passing[0]].lps_read_count(other);
                matched.push(other);
                if other_count > max_read_count {
                    max_read_count = other_count;
                    max_lps = other;
                }
            }
        }
        if matched.len() > 1 {
            uninformative.extend(matched.iter().filter(|&&x| x != max_lps));
            if passing.len() == 1 && !single_passing.contains(&max_lps) {
                single_passing.push(max_lps);
            }
        }
    }
    for &lps in &order {
        if uninformative.contains(&lps) {
            continue;
        }
        let passing = &lps_passing[&lps];
        let read_count = variants[passing[0]].lps_read_count(lps);
        for &other in &order {
            if other == lps || uninformative.contains(&other) {
                continue;
            }
            let other_passing = &lps_passing[&other];
            if other_passing.len() > passing.len() && passing.iter().all(|x| other_passing.contains(x)) {
                let other_count = variants[other_passing[0]].lps_read_count(other);
                if (read_count as f64) < other_count as f64 * SUBSET_READ_COUNT_LIMIT {
                    uninformative.push(lps);
                    if other_passing.len() == 1 && !single_passing.contains(&other) {
                        single_passing.push(other);
                    }
                    break;
                }
            }
        }
    }
    for &lps in &uninformative {
        for &i in &lps_variants[&lps] {
            variants[i].remove_lps(lps);
        }
    }
    // HashSet<Integer> iteration order
    for lps in java_int_hash_order(&single_passing) {
        let vs = &lps_variants[&lps];
        if vs.iter().all(|&i| variants[i].has_local_phase_sets() && variants[i].local_phase_sets().len() == 1) {
            for &i in vs {
                variants[i].remove_lps(lps);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compat::with_java;

    const POS: [i32; 6] = [100, 200, 300, 400, 500, 600];

    fn phaser() -> Phaser {
        let mut p = Phaser::new();
        p.pos = POS.to_vec();
        p
    }

    fn group_count(p: &Phaser) -> usize {
        p.collections.iter().map(|c| c.map.iter().map(|(_, l)| l.len()).sum::<usize>()).sum()
    }

    /// S1: a collection merged into another keeps its groups where both have a group at the same key
    #[test]
    fn merged_collections_keep_all_groups() {
        let run = || {
            let mut p = phaser();
            p.register(&[], &[2, 3], &[]); // 300-400
            p.register(&[], &[1], &[0]); // 100-200, a second collection, key 200
            p.register(&[], &[1, 2], &[]); // 200-300 joins the first collection under key 200 and merges the two
            assert_eq!(p.collections.len(), 1);
            group_count(&p)
        };
        assert_eq!(run(), 3);
        assert_eq!(with_java(true, run), 2);
    }

    /// S1: a group that matches one already there exactly is counted into it
    #[test]
    fn merge_counts_exact_matches_once() {
        let pos = POS.to_vec();
        let mut a = GroupCollection::new(false, false);
        a.add(200, 300, 0, &[1, 2], &[], &pos);
        let mut b = GroupCollection::new(false, false);
        b.add(200, 300, 1, &[1, 2], &[], &pos);
        b.add(200, 400, 2, &[1, 3], &[], &pos);
        a.merge(b, &pos);
        assert_eq!(a.map.len(), 1);
        assert_eq!(a.map[0].1.len(), 2);
        assert_eq!(a.map[0].1[0].read_count, 2);
        assert_eq!(a.group_count, 3);
    }

    /// S2: a group under a new key lowers the collection's minimum, it does not replace it
    #[test]
    fn collection_minimum_is_not_overwritten() {
        let run = || {
            let mut p = phaser();
            p.register(&[], &[0, 2], &[]); // 100-300
            p.register(&[], &[1, 2], &[]); // 200-300, a new key in the same collection
            (p.collections[0].min_position, p.collections[0].max_position)
        };
        assert_eq!(run(), (100, 300));
        assert_eq!(with_java(true, run), (200, 300));
    }

    /// S3: a single-variant group is dropped when no other group has its variant
    #[test]
    fn uninformative_single_variant_group_is_removed() {
        let run = || {
            let pos = POS.to_vec();
            let mut groups =
                vec![PhasedGroup::new(0, 100, 100, vec![0], vec![1], &pos), PhasedGroup::new(1, 200, 300, vec![1, 2], vec![], &pos)];
            merge_uninformative(&mut groups, &pos);
            groups.iter().map(|g| g.id).collect::<Vec<_>>()
        };
        assert_eq!(run(), vec![1]);
        assert_eq!(with_java(true, run), vec![0, 1]);
    }

    /// S3: after groups ahead of the current one are removed, the next group is not skipped
    #[test]
    fn uninformative_merge_visits_every_group() {
        let run = || {
            let pos = POS.to_vec();
            let mut groups = vec![
                PhasedGroup::new(0, 200, 200, vec![1], vec![], &pos),
                PhasedGroup::new(1, 200, 300, vec![1, 2], vec![], &pos),
                PhasedGroup::new(2, 400, 500, vec![3, 4], vec![], &pos),
                PhasedGroup::new(3, 400, 500, vec![3, 4], vec![5], &pos),
            ];
            merge_uninformative(&mut groups, &pos);
            groups.iter().map(|g| (g.id, g.read_count, g.allocated)).collect::<Vec<_>>()
        };
        assert_eq!(run(), vec![(1, 1, 1.0), (2, 2, 0.0)]);
        // Java skips group 2, and group 3 then takes it in
        assert_eq!(with_java(true, run), vec![(1, 1, 1.0), (3, 2, 0.0)]);
    }
}
