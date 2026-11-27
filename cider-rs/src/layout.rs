//! ReadLayout (the consensus of overlapping reads, lined up by their aligned position) and LayoutForest
//! (the tree that assigns reads to layouts along the best supported branches).
//!
//! The forest keeps its nodes in an arena (indices for Java's object references); children and level
//! lists are sorted with stable sorts, as Kotlin's sortByDescending, and maxByOrNull keeps the first
//! maximum.

use std::collections::HashSet;
use std::sync::{Arc, OnceLock};

use crate::constants::*;

const N: u8 = b'N';

/// ReadKey: read name and first of pair
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ReadKey {
    pub name: Arc<str>,
    pub first_of_pair: bool,
}

/// ReadLayout.Read: a read lined up in a layout; `source` is what the read stands for (the candidate and
/// slice of VJReadLayoutBuilder), shared by its copies
#[derive(Clone, Debug)]
pub struct LayoutRead {
    pub key: ReadKey,
    pub sequence: Arc<[u8]>,
    pub base_qualities: Arc<[u8]>,
    pub aligned_position: i32,
    pub source: Option<Arc<crate::layout_builder::VjLayoutSource>>,
}

impl LayoutRead {
    pub fn read_length(&self) -> i32 {
        self.sequence.len() as i32
    }
    /// Read.copy(alignedPosition)
    pub fn copy(&self, aligned_position: i32) -> LayoutRead {
        let mut r = self.clone();
        r.aligned_position = aligned_position;
        r
    }
}

#[derive(Clone, Debug, Default)]
struct BaseSupport {
    likely_base: u8,
    a: i32,
    c: i32,
    g: i32,
    t: i32,
}

impl BaseSupport {
    fn new() -> Self {
        BaseSupport { likely_base: N, ..Default::default() }
    }
    fn count(&self, base: u8) -> i32 {
        match base {
            b'A' => self.a,
            b'C' => self.c,
            b'G' => self.g,
            b'T' => self.t,
            _ => 0,
        }
    }
    /// addToCount: true if the likely base changed
    fn add_to_count(&mut self, base: u8) -> bool {
        match base {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' => self.t += 1,
            _ => {}
        }
        if self.likely_base != base && self.count(base) > self.count(self.likely_base) {
            self.likely_base = base;
            return true;
        }
        false
    }
}

#[derive(Clone, Debug, Default)]
struct SequenceSupport {
    sequence: Vec<u8>,
    support: Vec<BaseSupport>,
}

impl SequenceSupport {
    fn ensure_length(&mut self, l: usize) {
        if l > self.support.len() {
            let n = l - self.support.len();
            for _ in 0..n {
                self.support.push(BaseSupport::new());
            }
            self.sequence.extend(std::iter::repeat_n(N, n));
        }
    }
    fn shift_right_by(&mut self, s: i32) {
        if s <= 0 {
            return;
        }
        let s = s as usize;
        let mut v = Vec::with_capacity(self.support.len() + s);
        v.extend(std::iter::repeat_with(BaseSupport::new).take(s));
        v.append(&mut self.support);
        self.support = v;
        let mut seq = vec![N; s];
        seq.append(&mut self.sequence);
        self.sequence = seq;
    }
    fn update_sequence(&mut self) {
        if self.sequence.len() != self.support.len() {
            self.sequence = vec![0; self.support.len()];
        }
        for (s, b) in self.sequence.iter_mut().zip(&self.support) {
            *s = b.likely_base;
        }
    }
    fn counts(&self) -> Vec<i32> {
        self.sequence.iter().zip(&self.support).map(|(&s, b)| b.count(s)).collect()
    }
}

/// ReadLayout
#[derive(Debug)]
pub struct ReadLayout {
    pub id: String,
    all: SequenceSupport,
    high_qual: SequenceSupport,
    pub aligned_position: i32,
    reads: Vec<LayoutRead>,
    read_set: HashSet<ReadKey>,
    cache: OnceLock<Vec<u8>>,
}

impl Default for ReadLayout {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadLayout {
    pub fn new() -> ReadLayout {
        ReadLayout {
            id: String::new(),
            all: SequenceSupport::default(),
            high_qual: SequenceSupport::default(),
            aligned_position: 0,
            reads: Vec::new(),
            read_set: HashSet::new(),
            cache: OnceLock::new(),
        }
    }

    pub fn reads(&self) -> &[LayoutRead] {
        &self.reads
    }

    pub fn length(&self) -> i32 {
        self.all.support.len() as i32
    }

    /// consensusSequence: the high quality base where there is one, else the base of all reads
    pub fn consensus_sequence(&self) -> &[u8] {
        self.cache.get_or_init(|| {
            (0..self.all.sequence.len()).map(|i| if self.high_qual.sequence[i] == N { self.all.sequence[i] } else { self.high_qual.sequence[i] }).collect()
        })
    }

    pub fn consensus_string(&self) -> String {
        String::from_utf8_lossy(self.consensus_sequence()).into_owned()
    }

    /// highConfidenceSequence: N where fewer than max(reads * fraction, 1) high quality reads support the base
    pub fn high_confidence_sequence(&self, high_qual_read_fraction: f64) -> Vec<u8> {
        let min_count = ((self.reads.len() as f64 * high_qual_read_fraction) as i32).max(1);
        let counts = self.high_qual.counts();
        self.high_qual.sequence.iter().zip(&counts).map(|(&b, &c)| if c >= min_count { b } else { N }).collect()
    }

    pub fn high_qual_support_counts(&self) -> Vec<i32> {
        self.high_qual.counts()
    }

    pub fn high_qual_support_string(&self) -> String {
        counts_to_string(&self.high_qual_support_counts())
    }

    /// getHighQualSequenceSupportAt: (likely base, its count), (N, 0) for none
    pub fn high_qual_support_at(&self, index: usize) -> (u8, i32) {
        let s = &self.high_qual.support[index];
        if s.likely_base == N { (N, 0) } else { (s.likely_base, s.count(s.likely_base)) }
    }

    /// addRead: false if a read with the same key is in the layout
    pub fn add_read(&mut self, read: LayoutRead, min_base_quality: u8) -> bool {
        if !self.read_set.insert(read.key.clone()) {
            return false;
        }
        if self.reads.is_empty() {
            self.aligned_position = read.aligned_position;
        } else if read.aligned_position > self.aligned_position {
            let change = read.aligned_position - self.aligned_position;
            if change > 0 {
                self.all.shift_right_by(change);
                self.high_qual.shift_right_by(change);
            }
            self.aligned_position = read.aligned_position;
        }
        let offset = (self.aligned_position - read.aligned_position) as usize;
        let end = read.sequence.len() + offset;
        self.all.ensure_length(end);
        self.high_qual.ensure_length(end);
        let mut all_changed = false;
        let mut hq_changed = false;
        for (i, (&b, &q)) in read.sequence.iter().zip(read.base_qualities.iter()).enumerate() {
            if b != N {
                all_changed = self.all.support[offset + i].add_to_count(b) | all_changed;
                if q >= min_base_quality {
                    hq_changed = self.high_qual.support[offset + i].add_to_count(b) | hq_changed;
                }
            }
        }
        if all_changed {
            self.all.update_sequence();
        }
        if hq_changed {
            self.high_qual.update_sequence();
        }
        // CI11: Java clears its consensus cache only when a likely base changed, so a read that lengthens or
        // shifts the layout with N bases only leaves the old consensus (shorter, or misaligned) in place.
        // Corrected: cleared after every read added.
        if all_changed || hq_changed || !crate::compat::java("CI11") {
            self.cache = OnceLock::new();
        }
        self.reads.push(read);
        true
    }

    /// mergeIn
    pub fn merge_in(&mut self, layout: &ReadLayout, shift: i32, min_base_quality: u8) {
        for r in &layout.reads {
            if shift == 0 {
                self.add_read(r.clone(), min_base_quality);
            } else {
                self.add_read(r.copy(r.aligned_position + shift), min_base_quality);
            }
        }
    }

    /// ReadLayout.merge
    pub fn merge(l1: &ReadLayout, l2: &ReadLayout, shift1: i32, shift2: i32, min_base_quality: u8) -> ReadLayout {
        let mut m = ReadLayout::new();
        m.merge_in(l1, shift1, min_base_quality);
        m.merge_in(l2, shift2, min_base_quality);
        m
    }
}

/// CiderUtils.countsToString: base 36 digits, '#' for 36 or more
pub fn counts_to_string(counts: &[i32]) -> String {
    counts
        .iter()
        .map(|&c| if c >= 36 { '#' } else { std::char::from_digit(c.max(0) as u32, 36).unwrap() })
        .collect()
}

/// LayoutForest.Read: `source` is the index of what it stands for
#[derive(Clone, Debug)]
pub struct ForestRead {
    pub source: usize,
    pub sequence: Arc<[u8]>,
    pub base_qualities: Arc<[u8]>,
    pub aligned_position: i32,
}

#[derive(Debug)]
struct Node {
    position: i32,
    parent: Option<usize>,
    base: u8,
    high_qual_base: u8,
    count: i32,
    high_quality_count: i32,
    reads: Vec<usize>,
    children: Vec<usize>,
}

impl Node {
    fn support(&self) -> f64 {
        debug_assert!(self.high_quality_count <= self.count);
        self.high_quality_count as f64 + (self.count - self.high_quality_count) as f64 * LOW_QUAL_BASE_FACTOR
    }
}

/// LayoutForest
pub struct LayoutForest {
    min_base_quality: u8,
    min_overlap_length: i32,
    pub node_seal_factor: i32,
    nodes: Vec<Node>,
    roots: Vec<usize>,
    aligned_position: i32,
    level_nodes: Vec<Vec<usize>>,
    reads: Vec<ForestRead>,
}

impl LayoutForest {
    pub fn new(min_base_quality: u8, min_overlap_length: i32, node_seal_factor: i32) -> LayoutForest {
        LayoutForest {
            min_base_quality,
            min_overlap_length,
            node_seal_factor,
            nodes: Vec::new(),
            roots: Vec::new(),
            aligned_position: i32::MIN,
            level_nodes: Vec::new(),
            reads: Vec::new(),
        }
    }

#[allow(dead_code)]
    pub fn num_levels(&self) -> usize {
        self.level_nodes.len()
    }

    fn support(&self, n: usize) -> f64 {
        self.nodes[n].support()
    }

    fn to_layout_position(&self, read_aligned_pos: i32) -> i32 {
        self.aligned_position - read_aligned_pos
    }

    fn matches_node(&self, node: usize, base: u8, q: u8) -> bool {
        let n = &self.nodes[node];
        n.base == base || q < self.min_base_quality || n.high_qual_base == N || base == N
    }

    fn new_node(&mut self, position: i32, parent: Option<usize>) -> usize {
        self.nodes.push(Node { position, parent, base: N, high_qual_base: N, count: 0, high_quality_count: 0, reads: Vec::new(), children: Vec::new() });
        self.nodes.len() - 1
    }

    fn create_child(&mut self, parent: usize, level: i32, base: u8) -> usize {
        assert!(!self.is_node_sealed(parent), "java.lang.IllegalArgumentException: Failed requirement (child of a sealed node)");
        let c = self.new_node(level, Some(parent));
        self.nodes[c].base = base;
        self.nodes[parent].children.push(c);
        self.level_nodes_mut(level).push(c);
        c
    }

    fn level_nodes_mut(&mut self, level: i32) -> &mut Vec<usize> {
        while self.level_nodes.len() <= level as usize {
            self.level_nodes.push(Vec::new());
        }
        &mut self.level_nodes[level as usize]
    }

    fn add_to_node_count(&mut self, node: usize, base: u8, q: u8) {
        let hq_base = if q >= self.min_base_quality { base } else { N };
        let n = &mut self.nodes[node];
        if hq_base != N {
            assert!(n.high_qual_base == N || n.high_qual_base == hq_base, "java.lang.IllegalArgumentException: Failed requirement (node base)");
            n.high_quality_count += 1;
            if n.high_qual_base == N {
                n.base = hq_base;
                n.high_qual_base = hq_base;
            }
        }
        if base != N && n.base == base {
            n.count += 1;
        }
    }

    fn is_node_sealed(&self, node: usize) -> bool {
        self.nodes[node].children.iter().any(|&c| self.nodes[c].high_quality_count >= self.node_seal_factor)
    }

    fn check_update_aligned_position(&mut self, read_aligned_pos: i32) {
        if read_aligned_pos <= self.aligned_position {
            return;
        }
        let diff = read_aligned_pos.wrapping_sub(self.aligned_position);
        if !self.level_nodes.is_empty() {
            for l in &self.level_nodes {
                for &n in l {
                    self.nodes[n].position += diff;
                }
            }
            let mut v: Vec<Vec<usize>> = (0..diff).map(|_| Vec::new()).collect();
            v.append(&mut self.level_nodes);
            self.level_nodes = v;
            crate::warn!("WARNING: adding reads in wrong order. Reads should be added from largest alignedPosition to smallest");
        }
        self.aligned_position = read_aligned_pos;
    }

    /// tryAddRead(read): adds the read, and the reads that sealing nodes on its way takes off the tree
    pub fn try_add_read(&mut self, read: ForestRead) -> bool {
        self.reads.push(read);
        let r = self.reads.len() - 1;
        self.try_add_read_index(r)
    }

    fn try_add_read_index(&mut self, r: usize) -> bool {
        let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
        let result = self.try_add_read_inner(r, &mut queue);
        let mut count = 1;
        while let Some(next) = queue.pop_front() {
            let c = count;
            count += 1;
            if c % ADD_READ_INFINITE_LOOP_CHECK == 0 {
                crate::warn!("Too many reads ({}) reassigned in tryAddRead call. Checking tree for seal nodes", count);
                self.check_and_seal_nodes();
            }
            self.try_add_read_inner(next, &mut queue);
        }
        result
    }

    fn try_add_read_inner(&mut self, r: usize, reassign: &mut std::collections::VecDeque<usize>) -> bool {
        let aligned = self.reads[r].aligned_position;
        self.check_update_aligned_position(aligned);
        let pos = self.to_layout_position(aligned);
        if self.reads[r].sequence.is_empty() {
            return false;
        }
        if self.level_nodes.len() as i64 >= pos as i64 + self.min_overlap_length as i64 {
            let mut best: Option<usize> = None;
            for i in 0..self.level_nodes[pos as usize].len() {
                let branch = self.level_nodes[pos as usize][i];
                if (best.is_none() || self.support(branch) > self.support(best.unwrap())) && self.overlaps_branch(r, branch) {
                    best = Some(branch);
                }
            }
            if let Some(b) = best {
                self.add_read_to_branch(r, b, reassign);
                return true;
            }
        }
        let root = self.new_node(pos, None);
        let (b0, q0) = (self.reads[r].sequence[0], self.reads[r].base_qualities[0]);
        self.nodes[root].base = b0;
        self.nodes[root].high_qual_base = if q0 >= self.min_base_quality { b0 } else { N };
        self.level_nodes_mut(pos).push(root);
        self.add_read_to_branch(r, root, reassign);
        self.roots.push(root);
        true
    }

    /// overlapsBranch: the first min overlap bases match the branch (low quality bases match anything),
    /// and the rest can be added without branching a sealed node
    fn overlaps_branch(&self, r: usize, branch_root: usize) -> bool {
        let read = &self.reads[r];
        let pos = self.to_layout_position(read.aligned_position);
        assert!(pos >= 0 && self.nodes[branch_root].position == pos);
        let mut current = branch_root;
        if !self.matches_node(branch_root, read.sequence[0], read.base_qualities[0]) {
            return false;
        }
        'outer: for i in 1..read.sequence.len() {
            let q = read.base_qualities[i];
            let base = read.sequence[i];
            let hq_base = if q >= self.min_base_quality { base } else { N };
            let children = &self.nodes[current].children;
            if i as i32 >= self.min_overlap_length && children.is_empty() {
                break;
            }
            for &c in children {
                if self.nodes[c].base == base {
                    current = c;
                    continue 'outer;
                }
            }
            if hq_base == N {
                match self.max_support_child(current) {
                    Some(c) => {
                        current = c;
                        continue 'outer;
                    }
                    None => return false,
                }
            }
            for &c in children {
                if self.nodes[c].high_qual_base == N {
                    current = c;
                    continue 'outer;
                }
            }
            if i as i32 >= self.min_overlap_length {
                return !self.is_node_sealed(current);
            }
            return false;
        }
        true
    }

    /// maxByOrNull { support }: the first child of the highest support
    fn max_support_child(&self, node: usize) -> Option<usize> {
        let mut best: Option<usize> = None;
        for &c in &self.nodes[node].children {
            if best.is_none_or(|b| self.support(b) < self.support(c)) {
                best = Some(c);
            }
        }
        best
    }

    fn add_read_to_branch(&mut self, r: usize, branch_root: usize, reassign: &mut std::collections::VecDeque<usize>) {
        let read = self.reads[r].clone();
        let pos = self.to_layout_position(read.aligned_position);
        if read.sequence.is_empty() {
            return;
        }
        self.add_to_node_count(branch_root, read.sequence[0], read.base_qualities[0]);
        if let Some(p) = self.nodes[branch_root].parent {
            self.check_seal_node(p, reassign, Some(branch_root));
        }
        let mut current = branch_root;
        for i in 1..read.sequence.len() {
            let level = pos + i as i32;
            let q = read.base_qualities[i];
            let base = read.sequence[i];
            let hq_base = if q >= self.min_base_quality { base } else { N };
            let mut child = self.nodes[current].children.iter().copied().find(|&c| self.nodes[c].base == base);
            if child.is_none() {
                child = if hq_base == N {
                    self.max_support_child(current)
                } else {
                    self.nodes[current].children.iter().copied().find(|&c| self.nodes[c].high_qual_base == N)
                };
            }
            let child = match child {
                Some(c) => c,
                None => {
                    let c = self.create_child(current, level, base);
                    self.nodes[c].base = base;
                    c
                }
            };
            self.add_to_node_count(child, base, q);
            self.check_seal_node(current, reassign, Some(child));
            current = child;
        }
        self.nodes[current].reads.push(r);
    }

    /// sorts the children by support, descending (stable)
    fn sort_children(&mut self, node: usize) {
        let mut children = std::mem::take(&mut self.nodes[node].children);
        children.sort_by(|&a, &b| self.nodes[b].support().total_cmp(&self.nodes[a].support()));
        self.nodes[node].children = children;
    }

    /// checkSealNode: when the top child's support is at least the seal factor and more than the seal
    /// factor times the next, the other children go and their reads are to be added again
    fn check_seal_node(&mut self, node: usize, reassign: &mut impl Extend<usize>, protected: Option<usize>) -> bool {
        if self.nodes[node].children.len() <= 1 {
            return false;
        }
        self.sort_children(node);
        let children = &self.nodes[node].children;
        if protected.is_some_and(|p| children[0] != p) {
            return false;
        }
        let (s0, s1) = (self.support(children[0]), self.support(children[1]));
        let f = self.node_seal_factor as f64;
        if s0 < f || s0 <= s1 * f {
            return false;
        }
        let removed: Vec<usize> = children[1..].to_vec();
        for c in removed {
            let mut visited = Vec::new();
            self.depth_first_visit(c, &mut |_: &mut LayoutForest, n: usize| visited.push(n), false);
            for n in visited {
                reassign.extend(self.nodes[n].reads.iter().copied());
                self.remove_node(n);
            }
        }
        self.nodes[node].children.truncate(1);
        true
    }

    fn remove_node(&mut self, node: usize) {
        let pos = self.nodes[node].position as usize;
        let l = &mut self.level_nodes[pos];
        let i = l.iter().position(|&x| x == node).expect("java.lang.IllegalArgumentException: Failed requirement (node in its level)");
        l.remove(i);
    }

    /// depthFirstVisit; the visitor may change the children of the node it visits. With `seal` the visitor
    /// is checkSealNode collecting into the forest's reassign list.
    fn depth_first_visit(&mut self, start: usize, visitor: &mut dyn FnMut(&mut LayoutForest, usize), _seal: bool) {
        let mut stack: Vec<usize> = vec![0];
        let mut current = start;
        loop {
            let child_index = *stack.last().unwrap();
            if child_index == 0 {
                visitor(self, current);
            }
            if child_index == self.nodes[current].children.len() {
                if current == start {
                    return;
                }
                stack.pop();
                current = self.nodes[current].parent.expect("parent");
                continue;
            }
            *stack.last_mut().unwrap() += 1;
            stack.push(0);
            current = self.nodes[current].children[child_index];
        }
    }

    fn sort_node_children(&mut self) {
        for i in 0..self.roots.len() {
            let root = self.roots[i];
            self.depth_first_visit(root, &mut |f: &mut LayoutForest, n: usize| f.sort_children(n), false);
        }
    }

    fn sort_level_nodes(&mut self) {
        for l in 0..self.level_nodes.len() {
            let mut v = std::mem::take(&mut self.level_nodes[l]);
            v.sort_by(|&a, &b| self.nodes[b].support().total_cmp(&self.nodes[a].support()));
            self.level_nodes[l] = v;
        }
    }

    /// reassignReads: seals nodes and adds their reads again, with a falling seal factor
    fn reassign_reads(&mut self) {
        let original = self.node_seal_factor;
        for _ in 0..REASSIGN_READ_ROUNDS {
            let any = self.check_and_seal_nodes();
            if !any && self.node_seal_factor <= 1 {
                break;
            }
            self.node_seal_factor = (self.node_seal_factor - 1).max(1);
        }
        self.node_seal_factor = original;
    }

    /// checkAndSealNodes: true if reads were reassigned
    fn check_and_seal_nodes(&mut self) -> bool {
        self.sort_level_nodes();
        let mut reassign: Vec<usize> = Vec::new();
        loop {
            let prev = reassign.len();
            for i in 0..self.roots.len() {
                let root = self.roots[i];
                let mut collected: Vec<usize> = Vec::new();
                self.depth_first_visit(
                    root,
                    &mut |f: &mut LayoutForest, n: usize| {
                        f.check_seal_node(n, &mut collected, None);
                    },
                    true,
                );
                reassign.extend(collected);
            }
            if reassign.len() == prev {
                break;
            }
        }
        if reassign.is_empty() {
            return false;
        }
        crate::debug!("reassigning {} reads", reassign.len());
        // CI10: Java sorts a copy of the list by aligned position and drops it, so the reads go back in the
        // order collected; kept (sorting ascending would add them against the order the forest expects)
        for r in reassign {
            self.try_add_read_index(r);
        }
        true
    }

    /// buildReadLayouts: a layout from each leaf, along the highest supported paths; each read once
    pub fn build_read_layouts(&mut self, create: &dyn Fn(&ForestRead) -> LayoutRead) -> Vec<ReadLayout> {
        self.sort_node_children();
        self.reassign_reads();
        let mut leaves = Vec::new();
        for i in 0..self.roots.len() {
            let root = self.roots[i];
            self.depth_first_visit(
                root,
                &mut |f: &mut LayoutForest, n: usize| {
                    if n != root && f.nodes[n].children.is_empty() {
                        leaves.push(n);
                    }
                },
                false,
            );
        }
        let mut used = vec![false; self.reads.len()];
        let mut layouts = Vec::with_capacity(leaves.len());
        for leaf in leaves {
            layouts.push(self.build_layout_from_leaf(leaf, create, &mut used));
        }
        layouts
    }

    fn build_layout_from_leaf(&self, leaf: usize, create: &dyn Fn(&ForestRead) -> LayoutRead, used: &mut [bool]) -> ReadLayout {
        assert!(self.nodes[leaf].children.is_empty(), "cannot build layout from non leaf node");
        let mut reads = Vec::new();
        let mut current = Some(leaf);
        while let Some(n) = current {
            for &r in &self.nodes[n].reads {
                if used[r] {
                    continue;
                }
                used[r] = true;
                reads.push(create(&self.reads[r]));
            }
            current = self.nodes[n].parent;
        }
        reads.reverse();
        let mut layout = ReadLayout::new();
        for r in reads {
            layout.add_read(r, self.min_base_quality);
        }
        layout
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::read::test_util::phred;

    pub const MIN_BASE_QUALITY: u8 = 30;
    const MIN_OVERLAP_BASES: i32 = 3;
    const MIN_SUPPORT_TO_SEAL_NODE: i32 = 3;

    /// TestUtils.createLayoutRead
    pub fn layout_read(id: &str, seq: &str, quals: &str, aligned_position: i32, first_of_pair: bool) -> LayoutRead {
        LayoutRead {
            key: ReadKey { name: Arc::from(id), first_of_pair },
            sequence: Arc::from(seq.as_bytes()),
            base_qualities: Arc::from(phred(quals).as_slice()),
            aligned_position,
            source: None,
        }
    }

    pub fn names() -> std::cell::RefCell<Vec<String>> {
        std::cell::RefCell::new(Vec::new())
    }

    struct Forest {
        f: LayoutForest,
        ids: Vec<String>,
    }

    impl Forest {
        fn new() -> Forest {
            Self::with_seal(MIN_SUPPORT_TO_SEAL_NODE)
        }
        fn with_seal(seal: i32) -> Forest {
            Forest { f: LayoutForest::new(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, seal), ids: Vec::new() }
        }
        fn add(&mut self, id: &str, seq: &str, quals: &str, aligned: i32) -> bool {
            self.ids.push(id.to_string());
            self.f.try_add_read(ForestRead {
                source: self.ids.len() - 1,
                sequence: Arc::from(seq.as_bytes()),
                base_qualities: Arc::from(phred(quals).as_slice()),
                aligned_position: aligned,
            })
        }
        fn layouts(&mut self) -> Vec<ReadLayout> {
            let ids = self.ids.clone();
            self.f.build_read_layouts(&|r: &ForestRead| LayoutRead {
                key: ReadKey { name: Arc::from(ids[r.source].as_str()), first_of_pair: true },
                sequence: r.sequence.clone(),
                base_qualities: r.base_qualities.clone(),
                aligned_position: r.aligned_position,
                source: None,
            })
        }
    }

    fn sorted(mut l: Vec<ReadLayout>) -> Vec<ReadLayout> {
        l.sort_by_key(|x| x.consensus_string());
        l
    }

    // LayoutForestTest.testAddReadSimple1
    #[test]
    fn add_read_simple() {
        let mut f = Forest::new();
        assert!(f.add("read1", "CAGGT", "FFFFF", 0));
        assert_eq!(f.f.num_levels(), 5);
        let l = f.layouts();
        assert_eq!(l.len(), 1);
        assert_eq!(l[0].consensus_string(), "CAGGT");
        assert_eq!(l[0].high_qual_support_string(), "11111");
        assert!(f.add("read2", "CAGGT", "FF:FF", 0));
        let l = f.layouts();
        assert_eq!(l.len(), 1);
        assert_eq!(l[0].consensus_string(), "CAGGT");
        assert_eq!(l[0].high_qual_support_string(), "22122");
    }

    // testAddReadWrongOrder
    #[test]
    fn add_read_wrong_order() {
        let mut f = Forest::new();
        f.add("r2", "CAGCT", "FFFFF", 4);
        f.add("r1", "GCAGCT", "FFFFFF", 5);
        assert_eq!(f.f.num_levels(), 6);
        let l = f.layouts();
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGCT");
        assert_eq!(l[0].high_qual_support_string(), "11111");
        assert_eq!(l[1].consensus_string(), "GCAGCT");
        assert_eq!(l[1].high_qual_support_string(), "111111");
    }

    // testAddReadSimpleBranch
    #[test]
    fn add_read_simple_branch() {
        let mut f = Forest::new();
        assert!(f.add("read1", "CAGGT", "FFFFF", 0));
        assert!(f.add("read2", "CATGT", "FFFFF", 0));
        let l = sorted(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGGT");
        assert_eq!(l[0].high_qual_support_string(), "11111");
        assert_eq!(l[1].consensus_string(), "CATGT");
        assert_eq!(l[1].high_qual_support_string(), "11111");
    }

    // testAddReadLowQualMismatch
    #[test]
    fn add_read_low_qual_mismatch() {
        let mut f = Forest::new();
        f.add("read1", "CAGGT", "FFFFF", 0);
        f.add("read2", "CATAT", "FF::F", 0);
        let l = sorted(f.layouts());
        assert_eq!(l.len(), 1);
        assert_eq!(l[0].consensus_string(), "CAGGT");
        assert_eq!(l[0].high_qual_support_string(), "22112");
    }

    // testAddReadLowQualMismatch1
    #[test]
    fn add_read_low_qual_mismatch_1() {
        let mut f = Forest::new();
        f.add("read1", "CATAT", "FF::F", 0);
        f.add("read2", "CAGCT", "FFFFF", 0);
        let l = sorted(f.layouts());
        assert_eq!(l.len(), 1);
        assert_eq!(l[0].consensus_string(), "CAGCT");
        assert_eq!(l[0].high_qual_support_string(), "22112");
    }

    // testLowQualBranch
    #[test]
    fn low_qual_branch() {
        let mut f = Forest::new();
        f.add("read1", "CAGCT", "FFFFF", 0);
        f.add("read2", "CATCA", "FF:FF", 0);
        f.add("read3", "CAGCA", "FFFFF", 0);
        let l = sorted(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGCA");
        assert_eq!(l[0].high_qual_support_string(), "22122");
        assert_eq!(l[1].consensus_string(), "CAGCT");
        assert_eq!(l[1].high_qual_support_string(), "11111");
    }

    // testLowQualBranchExtend
    #[test]
    fn low_qual_branch_extend() {
        let mut f = Forest::new();
        f.add("read1", "CAGCT", "FF:FF", 0);
        f.add("read2", "CATCA", "FF:FF", 0);
        f.add("read3", "CAGTAT", "FFF:FF", 0);
        let l = sorted(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGCAT");
        assert_eq!(l[0].high_qual_support_string(), "221121");
        assert_eq!(l[1].consensus_string(), "CAGCT");
        assert_eq!(l[1].high_qual_support_string(), "11011");
    }

    fn by_length(mut l: Vec<ReadLayout>) -> Vec<ReadLayout> {
        l.sort_by(|a, b| b.length().cmp(&a.length()));
        l
    }

    // testAddReadAfterError
    #[test]
    fn add_read_after_error() {
        let mut f = Forest::new();
        f.add("r1", "GCAGCT", "FFFFFF", 5);
        f.add("r2", "CAGCT", "FFFFF", 4);
        f.add("r3", "AGCTA", "FFFFF", 3);
        f.add("r4", "AGCTG", "FFFFF", 3);
        f.add("r5", "CTGAA", "FFFFF", 1);
        let l = by_length(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "GCAGCTGAA");
        assert_eq!(l[0].high_qual_support_string(), "123344211");
        assert_eq!(l[0].reads().len(), 4);
        assert_eq!(l[0].aligned_position, 5);
        assert_eq!(l[1].consensus_string(), "AGCTA");
        assert_eq!(l[1].high_qual_support_string(), "11111");
        assert_eq!(l[1].reads().len(), 1);
        assert_eq!(l[1].aligned_position, 3);
    }

    // testAddNewRootBranch
    #[test]
    fn add_new_root_branch() {
        let mut f = Forest::new();
        assert!(f.add("read1", "CAGCTGGAC", "FFFFFFFFF", 6));
        assert!(f.add("read2", "TCCA", "FF:F", 4));
        let l = by_length(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGCTGGAC");
        assert_eq!(l[0].aligned_position, 6);
        assert_eq!(l[1].consensus_string(), "TCCA");
        assert_eq!(l[1].aligned_position, 4);
        assert!(f.add("read3", "CCATT", "FFFFF", 3));
        let l = by_length(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGCTGGAC");
        assert_eq!(l[0].aligned_position, 6);
        assert_eq!(l[1].consensus_string(), "TCCATT");
        assert_eq!(l[1].aligned_position, 4);
    }

    // testAddReadToHighestSupport
    #[test]
    fn add_read_to_highest_support() {
        let mut f = Forest::new();
        assert!(f.add("read1", "CAGCT", "FF:FF", 0));
        assert!(f.add("read2", "CATCA", "FF:FF", 0));
        assert!(f.add("read3", "CAGTAT", "FFF:FF", 0));
        let l = sorted(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "CAGCAT");
        assert_eq!(l[0].high_qual_support_string(), "221121");
        assert_eq!(l[1].consensus_string(), "CAGCT");
        assert_eq!(l[1].high_qual_support_string(), "11011");
    }

    // testCannotAddToSealedNode
    #[test]
    fn cannot_add_to_sealed_node() {
        let mut f = Forest::with_seal(2);
        f.add("r1", "GCAGCT", "FFFFFF", 5);
        f.add("r2", "CAGCT", "FFFFF", 4);
        f.add("r3", "AGCTGA", "FFFFFF", 3);
        f.add("r4", "AGCTGA", "FFFFFF", 3);
        f.add("r5", "AGCTAA", "FFFFFF", 3);
        assert_eq!(f.f.roots.len(), 2);
        let l = by_length(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "GCAGCTGA");
        assert_eq!(l[0].high_qual_support_string(), "12444422");
        assert_eq!(l[0].reads().len(), 4);
        assert_eq!(l[0].aligned_position, 5);
        assert_eq!(l[1].consensus_string(), "AGCTAA");
        assert_eq!(l[1].high_qual_support_string(), "111111");
        assert_eq!(l[1].reads().len(), 1);
        assert_eq!(l[1].aligned_position, 3);
    }

    // testReassignedReads
    #[test]
    fn reassigned_reads() {
        let mut f = Forest::new();
        f.add("r1", "GCAGCT", "FFFFFF", 4);
        f.add("r2", "AGCTGAA", "FFFFFFF", 2);
        f.add("r3", "AGCTAA", "FFFFFF", 2);
        f.add("r4", "GCTAAC", "FFFFFF", 1);
        f.add("r5", "TAGCTGA", "FFFFFFF", 3);
        assert_eq!(f.f.roots.len(), 2);
        let l = by_length(f.layouts());
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].consensus_string(), "GCAGCTAAC");
        assert_eq!(l[0].reads().len(), 3);
        assert_eq!(l[0].aligned_position, 4);
        assert_eq!(l[1].consensus_string(), "TAGCTGAA");
        assert_eq!(l[1].high_qual_support_string(), "12222221");
        assert_eq!(l[1].reads().len(), 2);
        assert_eq!(l[1].aligned_position, 3);
    }

    #[test]
    fn layout_add_and_merge() {
        let mut l = ReadLayout::new();
        assert!(l.add_read(layout_read("a", "ACGT", "FFFF", 1, true), MIN_BASE_QUALITY));
        assert!(!l.add_read(layout_read("a", "ACGT", "FFFF", 1, true), MIN_BASE_QUALITY));
        assert!(l.add_read(layout_read("b", "TACG", "FF#F", 2, true), MIN_BASE_QUALITY));
        assert_eq!(l.aligned_position, 2);
        assert_eq!(l.consensus_string(), "TACGT");
        assert_eq!(l.high_qual_support_string(), "12121");
        let m = ReadLayout::merge(&l, &l, 0, 1, MIN_BASE_QUALITY);
        // the same keys: the second layout adds nothing
        assert_eq!(m.reads().len(), 2);
        assert_eq!(counts_to_string(&[0, 9, 10, 35, 36, 100]), "09az##");
    }
}
