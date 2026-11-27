//! VDJ sequences: the anchors (VJAnchorByBlosum, VJAnchorByReadMatch), VDJSequence, VdjBuilderUtils and
//! VDJSequenceBuilder (joining V and J layouts, completing one-sided layouts with the BLOSUM search, and
//! merging identical sequences).

use std::collections::HashMap;
use std::sync::Arc;

use anyhow::{bail, Result};
use hmf_common::codon::amino_acid_from_bases;

use crate::blosum::{AnchorSearcher, SearchMode};
use crate::compat;
use crate::constants::*;
use crate::genes::{GeneDatastore, Vj, VjGeneType};
use crate::layout::{LayoutRead, ReadLayout};
use crate::layout_builder::{candidate, LayoutAdaptor};
use crate::screener::MatchMethod;

const N: u8 = b'N';

#[derive(Clone, Debug, PartialEq)]
pub enum AnchorKind {
    /// VJAnchorByBlosum: the template genes (datastore indices) and the similarity score
    Blosum { template_genes: Vec<usize>, similarity_score: i32 },
    /// VJAnchorByReadMatch: the match method of the reads and their number
    ReadMatch { match_method: String, num_reads: i32 },
}

/// VJAnchor
#[derive(Clone, Debug, PartialEq)]
pub struct VjAnchor {
    pub vj: Vj,
    pub gene_type: VjGeneType,
    /// the last base of V + 1, or the first base of J
    pub anchor_boundary: i32,
    pub template_anchor_seq: String,
    pub kind: AnchorKind,
}

impl VjAnchor {
    pub fn match_method(&self) -> &str {
        match &self.kind {
            AnchorKind::Blosum { .. } => "blosum",
            AnchorKind::ReadMatch { match_method, .. } => match_method,
        }
    }
#[allow(dead_code)]
    pub fn num_reads(&self) -> Option<i32> {
        match &self.kind {
            AnchorKind::ReadMatch { num_reads, .. } => Some(*num_reads),
            _ => None,
        }
    }
    fn with_boundary(&self, b: i32) -> VjAnchor {
        let mut a = self.clone();
        a.anchor_boundary = b;
        a
    }
}

/// VDJSequence: a slice of a layout with its anchors (positions relative to the slice)
#[derive(Clone, Debug)]
pub struct VdjSequence {
    pub layout: Arc<ReadLayout>,
    pub layout_slice_start: i32,
    pub layout_slice_end: i32,
    pub v_anchor: Option<VjAnchor>,
    pub j_anchor: Option<VjAnchor>,
}

fn floor_mod(a: i32, m: i32) -> i32 {
    a.rem_euclid(m)
}

impl VdjSequence {
    pub fn new(layout: Arc<ReadLayout>, start: i32, end: i32, v: Option<VjAnchor>, j: Option<VjAnchor>) -> VdjSequence {
        assert!(start >= 0, "java.lang.IllegalArgumentException: Failed requirement (layoutSliceStart >= 0)");
        assert!(end <= layout.length(), "java.lang.IllegalArgumentException: Failed requirement (layoutSliceEnd <= layout.length)");
        VdjSequence { layout, layout_slice_start: start, layout_slice_end: end, v_anchor: v, j_anchor: j }
    }
    pub fn num_reads(&self) -> i32 {
        self.layout.reads().len() as i32
    }
    pub fn length(&self) -> i32 {
        self.layout_slice_end - self.layout_slice_start
    }
    pub fn sequence(&self) -> &[u8] {
        &self.layout.consensus_sequence()[self.layout_slice_start as usize..self.layout_slice_end as usize]
    }
    pub fn v_anchor_boundary(&self) -> Option<i32> {
        self.v_anchor.as_ref().map(|a| a.anchor_boundary)
    }
    pub fn j_anchor_boundary(&self) -> Option<i32> {
        self.j_anchor.as_ref().map(|a| a.anchor_boundary)
    }
    /// aminoAcidSequence
    pub fn amino_acid_sequence(&self) -> String {
        let shift = floor_mod(self.v_anchor_boundary().or(self.j_anchor_boundary()).unwrap_or(0), 3) as usize;
        let s = self.sequence();
        amino_acid_from_bases(&s[shift.min(s.len())..])
    }
    /// aminoAcidSequenceFormatted (for logs)
    pub fn amino_acid_sequence_formatted(&self) -> String {
        let s = self.sequence();
        let sub = |a: i32, b: i32| -> &[u8] { &s[(a.max(0) as usize).min(s.len())..(b.max(0) as usize).min(s.len())] };
        match (&self.v_anchor, &self.j_anchor) {
            (None, None) => String::new(),
            (None, Some(j)) => insert_dashes(&amino_acid_from_bases(sub(floor_mod(j.anchor_boundary, 3), s.len() as i32)), &[j.anchor_boundary / 3]),
            (Some(v), None) => insert_dashes(&amino_acid_from_bases(sub(floor_mod(v.anchor_boundary, 3), s.len() as i32)), &[v.anchor_boundary / 3]),
            (Some(v), Some(j)) => format!(
                "{}-{}",
                insert_dashes(&amino_acid_from_bases(sub(floor_mod(v.anchor_boundary, 3), j.anchor_boundary)), &[v.anchor_boundary / 3]),
                amino_acid_from_bases(self.j_anchor_sequence().unwrap_or_default())
            ),
        }
    }
    #[allow(dead_code)]
    /// sequenceFormatted: dashes at the anchor boundaries
    pub fn sequence_formatted(&self) -> String {
        let mut positions = Vec::new();
        if let Some(v) = self.v_anchor_boundary() {
            positions.push(v);
        }
        if let Some(j) = self.j_anchor_boundary() {
            positions.push(j);
        }
        insert_dashes(&String::from_utf8_lossy(self.sequence()), &positions)
    }
    #[allow(dead_code)]
    /// cdr3SequenceShort: without the C and the W/F
    pub fn cdr3_sequence_short(&self) -> &[u8] {
        let s = self.sequence();
        if let Some(v) = self.v_anchor_boundary() {
            let end = self.j_anchor_boundary().unwrap_or(self.length());
            return &s[v as usize..end as usize];
        }
        if let Some(j) = self.j_anchor_boundary() {
            return &s[floor_mod(j, 3) as usize..j as usize];
        }
        s
    }
    pub fn v_anchor_length(&self) -> i32 {
        self.v_anchor_boundary().unwrap_or(0)
    }
    pub fn j_anchor_length(&self) -> i32 {
        self.length() - self.j_anchor_boundary().unwrap_or(0)
    }
    /// vAnchorSequence: sequence.take(boundary)
    pub fn v_anchor_sequence(&self) -> Option<&[u8]> {
        let b = self.v_anchor_boundary()?;
        let s = self.sequence();
        Some(&s[..(b.max(0) as usize).min(s.len())])
    }
    /// jAnchorSequence: sequence.substring(boundary)
    pub fn j_anchor_sequence(&self) -> Option<&[u8]> {
        let b = self.j_anchor_boundary()?;
        let s = self.sequence();
        Some(&s[(b.max(0) as usize).min(s.len())..])
    }
    pub fn cdr3_start(&self) -> i32 {
        if let Some(v) = self.v_anchor_boundary() {
            return (v - 3).max(0);
        }
        if let Some(j) = self.j_anchor_boundary() {
            return floor_mod(j, 3);
        }
        0
    }
    pub fn cdr3_end(&self) -> i32 {
        match self.j_anchor_boundary() {
            None => self.length(),
            Some(j) => ((j + 3).min(self.length())).max(self.cdr3_start()),
        }
    }
    pub fn cdr3_sequence(&self) -> &[u8] {
        &self.sequence()[self.cdr3_start() as usize..self.cdr3_end() as usize]
    }
    pub fn support_counts(&self) -> Vec<i32> {
        self.layout.high_qual_support_counts()[self.layout_slice_start as usize..self.layout_slice_end as usize].to_vec()
    }
    pub fn is_in_frame(&self) -> bool {
        match (self.v_anchor_boundary(), self.j_anchor_boundary()) {
            (Some(v), Some(j)) => floor_mod(j - v, 3) == 0,
            _ => true,
        }
    }
    pub fn support_at(&self, index: i32) -> (u8, i32) {
        self.layout.high_qual_support_at((self.layout_slice_start + index) as usize)
    }
}

/// CiderUtils.insertDashes: a dash before each position inside the string
pub fn insert_dashes(s: &str, positions: &[i32]) -> String {
    let mut out = s.to_string();
    for &p in positions.iter().rev() {
        if p >= 0 && (p as usize) < s.len() {
            out.insert(p as usize, '-');
        }
    }
    out
}

// ---- VdjBuilderUtils ----

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SequenceOverlap {
    pub seq1_offset: i32,
    pub seq2_offset: i32,
    pub overlap_bases: i32,
    pub high_qual_match_bases: i32,
}

/// calcSequenceWordHashes: 2 bits per base, a hash for each word without N
pub fn calc_sequence_word_hashes(seq: &[u8], word_size: i32) -> Result<Vec<i32>> {
    if word_size > 16 {
        bail!("illegal word size: {word_size} > 16");
    }
    let mask = (0xFFFFFFFFu64 >> (32 - 2 * word_size)) as i64 as i32;
    let mut hashes = Vec::new();
    let mut hash: i32 = 0;
    let mut next_word_end = word_size - 1;
    for (i, &b) in seq.iter().enumerate() {
        if b == N {
            next_word_end += word_size;
            hash = 0;
            continue;
        }
        let h = match b {
            b'A' => 0,
            b'T' => 1,
            b'C' => 2,
            b'G' => 3,
            _ => bail!("unknown base: {}", b as char),
        };
        hash = (hash.wrapping_shl(2)).wrapping_add(h) & mask;
        if i as i32 == next_word_end {
            hashes.push(hash);
            next_word_end += 1;
        }
    }
    Ok(hashes)
}

/// findSequenceOverlap: the first offset (from seq1 starting right of seq2's start to left of it) where
/// the high quality bases agree and more than `min_overlapped_bases` of them match
pub fn find_sequence_overlap(seq1: &[u8], seq2: &[u8], min_overlapped_bases: i32) -> Option<SequenceOverlap> {
    let (n1, n2) = (seq1.len() as i32, seq2.len() as i32);
    for i in -(n1 - min_overlapped_bases)..(n2 - min_overlapped_bases) {
        let mut matched = true;
        let mut hq = 0;
        let overlap = if i > 0 { (n2 - i).min(n1) } else { (n1 + i).min(n2) };
        for j in 0..overlap {
            let (a, b) = if i > 0 { (seq1[j as usize], seq2[(i + j) as usize]) } else { (seq1[(-i + j) as usize], seq2[j as usize]) };
            if a != N && b != N {
                if a == b {
                    hq += 1;
                } else {
                    matched = false;
                    break;
                }
            }
        }
        if matched && hq > min_overlapped_bases {
            return Some(if i > 0 {
                SequenceOverlap { seq1_offset: i, seq2_offset: 0, overlap_bases: overlap, high_qual_match_bases: hq }
            } else {
                SequenceOverlap { seq1_offset: 0, seq2_offset: -i, overlap_bases: overlap, high_qual_match_bases: hq }
            });
        }
    }
    None
}

/// getCompareRanges: the ranges (start, end exclusive) of both sequences lined up at an anchor
pub fn compare_ranges(v1: &VdjSequence, v2: &VdjSequence) -> ((i32, i32), (i32, i32)) {
    assert_eq!(v1.v_anchor.is_some(), v2.v_anchor.is_some());
    assert_eq!(v1.j_anchor.is_some(), v2.j_anchor.is_some());
    let offset = match (v1.v_anchor_boundary(), v2.v_anchor_boundary(), v1.j_anchor_boundary(), v2.j_anchor_boundary()) {
        (Some(a), Some(b), _, _) => a - b,
        (_, _, Some(a), Some(b)) => a - b,
        _ => panic!("java.lang.RuntimeException (no anchors)"),
    };
    let (mut s1, mut s2) = (0, 0);
    if offset >= 0 {
        s1 = offset;
    } else {
        s2 = -offset;
    }
    let n = (v1.length() - s1).min(v2.length() - s2);
    let (mut e1, mut e2) = (s1 + n, s2 + n);
    if v1.v_anchor.is_none() && v2.v_anchor.is_none() {
        s1 = (v1.j_anchor_boundary().unwrap() - PARTIAL_VDJ_UNANCHORED_LENGTH_BASES).max(s1);
        s2 = (v2.j_anchor_boundary().unwrap() - PARTIAL_VDJ_UNANCHORED_LENGTH_BASES).max(s2);
    } else if v1.j_anchor.is_none() && v2.j_anchor.is_none() {
        e1 = (v1.v_anchor_boundary().unwrap() + PARTIAL_VDJ_UNANCHORED_LENGTH_BASES).min(e1);
        e2 = (v2.v_anchor_boundary().unwrap() + PARTIAL_VDJ_UNANCHORED_LENGTH_BASES).min(e2);
    }
    ((s1, e1), (s2, e2))
}

/// vdjSequenceCompare: walks the lined-up bases, calling `diff` where the high quality bases differ;
/// false as soon as it does
pub fn vdj_sequence_compare(v1: &VdjSequence, v2: &VdjSequence, diff: &mut dyn FnMut((u8, i32), (u8, i32)) -> bool) -> bool {
    if let (Some(va), Some(ja), Some(vb), Some(jb)) = (v1.v_anchor_boundary(), v1.j_anchor_boundary(), v2.v_anchor_boundary(), v2.j_anchor_boundary()) {
        if ja - va != jb - vb {
            return false;
        }
    }
    let ((s1, e1), (s2, e2)) = compare_ranges(v1, v2);
    let (mut i1, mut i2) = (s1, s2);
    while i1 < e1 && i2 < e2 {
        let a = v1.support_at(i1);
        let b = v2.support_at(i2);
        if a.0 != b.0 && !diff(a, b) {
            return false;
        }
        i1 += 1;
        i2 += 1;
    }
    true
}

/// mergeAnchors
///
/// CI6: of two anchors matched by reads Java means to keep the one with more reads, but compares the
/// first anchor with itself (maxOf(anchor1, anchor1, ...)), so it always keeps the first. Corrected: the one
/// with more reads (the first on a tie).
pub fn merge_anchors(a1: Option<&VjAnchor>, a2: Option<&VjAnchor>, len1: i32, len2: i32, boundary: i32) -> Option<VjAnchor> {
    let a1 = a1?;
    let a2 = a2.expect("java.lang.IllegalArgumentException: Failed requirement (anchors on both)");
    match (&a1.kind, &a2.kind) {
        (AnchorKind::ReadMatch { num_reads: n1, .. }, AnchorKind::ReadMatch { num_reads: n2, .. }) => {
            let more = if compat::java("CI6") || n1 >= n2 { a1 } else { a2 };
            let mut a = more.with_boundary(boundary);
            if let AnchorKind::ReadMatch { num_reads, .. } = &mut a.kind {
                *num_reads = n1 + n2;
            }
            Some(a)
        }
        (AnchorKind::ReadMatch { .. }, _) => Some(a1.with_boundary(boundary)),
        (_, AnchorKind::ReadMatch { .. }) => Some(a2.with_boundary(boundary)),
        _ => Some(if len1 >= len2 { a1 } else { a2 }.with_boundary(boundary)),
    }
}

/// mergeVDJs: one VDJ from two lined up at their anchors, on the layout of both
pub fn merge_vdjs(v1: &VdjSequence, v2: &VdjSequence, min_base_quality: u8) -> Result<VdjSequence> {
    assert_eq!(v1.v_anchor.is_some(), v2.v_anchor.is_some());
    assert_eq!(v1.j_anchor.is_some(), v2.j_anchor.is_some());
    if let (Some(va), Some(ja)) = (v1.v_anchor_boundary(), v1.j_anchor_boundary()) {
        if ja - va != v2.j_anchor_boundary().unwrap() - v2.v_anchor_boundary().unwrap() {
            crate::error!("cannot merge {} and {}, different lengths between anchors", v1.amino_acid_sequence_formatted(), v2.amino_acid_sequence_formatted());
            bail!("cannot merge VDJs: different lengths between anchors");
        }
    }
    let vb = v1.v_anchor_boundary().unwrap_or(0).max(v2.v_anchor_boundary().unwrap_or(0));
    let jb = v1.j_anchor_boundary().unwrap_or(0).max(v2.j_anchor_boundary().unwrap_or(0));
    let v_anchor = merge_anchors(v1.v_anchor.as_ref(), v2.v_anchor.as_ref(), v1.v_anchor_length(), v2.v_anchor_length(), vb);
    let j_anchor = merge_anchors(v1.j_anchor.as_ref(), v2.j_anchor.as_ref(), v1.j_anchor_length(), v2.j_anchor_length(), jb);
    let j_anchor_length = v1.j_anchor_length().max(v2.j_anchor_length());
    let (p, s) = if v1.num_reads() > v2.num_reads() || (v1.num_reads() == v2.num_reads() && v1.length() > v2.length()) { (v1, v2) } else { (v2, v1) };
    let anchor_ref = |x: &VdjSequence| x.v_anchor_boundary().unwrap_or_else(|| x.j_anchor_boundary().unwrap());
    let p_aligned = p.layout.aligned_position - p.layout_slice_start - anchor_ref(p);
    let s_aligned = s.layout.aligned_position - s.layout_slice_start - anchor_ref(s);
    let shift = p_aligned - s_aligned;
    let mut merged = ReadLayout::merge(&p.layout, &s.layout, 0, shift, min_base_quality);
    merged.id = format!("{};{}", p.layout.id, s.layout.id);
    let boundary_diff = match &v_anchor {
        Some(v) => p.v_anchor_boundary().unwrap() - v.anchor_boundary,
        None => p.j_anchor_boundary().unwrap() - j_anchor.as_ref().unwrap().anchor_boundary,
    };
    let slice_start = p.layout_slice_start + boundary_diff + merged.aligned_position - p.layout.aligned_position;
    assert!(v_anchor.is_some() || j_anchor.is_some());
    let slice_end = match &j_anchor {
        None => {
            let post_v = (p.length() - p.v_anchor_boundary().unwrap()).max(s.length() - s.v_anchor_boundary().unwrap());
            slice_start + v_anchor.as_ref().unwrap().anchor_boundary + post_v
        }
        Some(j) => slice_start + j.anchor_boundary + j_anchor_length,
    };
    if slice_start < 0 {
        crate::info!(
            "unable to merge {} and {}, returning one of the vdj: {}",
            v1.amino_acid_sequence_formatted(),
            v2.amino_acid_sequence_formatted(),
            p.amino_acid_sequence_formatted()
        );
        return Ok(p.clone());
    }
    let slice_end = slice_end.min(merged.length());
    Ok(VdjSequence::new(Arc::new(merged), slice_start, slice_end, v_anchor, j_anchor))
}

// ---- VDJSequenceBuilder ----

pub struct VdjSequenceBuilder<'a> {
    pub adaptor: &'a dyn LayoutAdaptor,
    pub ds: &'a GeneDatastore,
    pub searcher: &'a dyn AnchorSearcher,
    pub min_base_quality: u8,
    pub min_overlapped_bases: i32,
}

impl VdjSequenceBuilder<'_> {
    fn anchor_range(&self, vj: Vj, layout: &ReadLayout) -> Option<(i32, i32)> {
        self.adaptor.anchor_range(vj, layout)
    }

    /// createVJAnchorByReadMatch
    fn anchor_by_read_match(&self, anchor_boundary: i32, layout: &ReadLayout, gene_type: VjGeneType) -> VjAnchor {
        let method = match self.adaptor.anchor_match_method(layout) {
            MatchMethod::Align => "align",
            MatchMethod::Blosum => "blosum",
        };
        VjAnchor {
            vj: gene_type.vj(),
            gene_type,
            anchor_boundary,
            template_anchor_seq: self.adaptor.template_anchor_sequence(layout),
            kind: AnchorKind::ReadMatch { match_method: method.to_string(), num_reads: layout.reads().len() as i32 },
        }
    }

    /// buildVDJSequences: `layouts` holds the layouts of every gene type, in gene type order
    pub fn build_vdj_sequences(
        &self,
        layouts: &[(VjGeneType, Vec<Arc<ReadLayout>>)],
        mate_layout_reads: &HashMap<String, Vec<LayoutRead>>,
        threads: usize,
    ) -> Result<Vec<VdjSequence>> {
        let mut vdjs: Vec<VdjSequence> = Vec::new();
        let mut one_sided: Vec<(VjGeneType, Vec<Arc<ReadLayout>>)> = layouts.to_vec();
        one_sided.sort_by_key(|(t, _)| *t);
        let index = |t: VjGeneType, os: &Vec<(VjGeneType, Vec<Arc<ReadLayout>>)>| os.iter().position(|(x, _)| *x == t);

        // joinVjLayouts
        let v_types: Vec<VjGeneType> = one_sided.iter().map(|(t, _)| *t).filter(|t| t.vj() == Vj::V).collect();
        for v_type in v_types {
            for &j_type in v_type.paired() {
                let vi = index(v_type, &one_sided);
                let ji = index(j_type, &one_sided);
                let mut v_layouts = vi.map(|i| std::mem::take(&mut one_sided[i].1)).unwrap_or_default();
                let mut j_layouts = ji.map(|i| std::mem::take(&mut one_sided[i].1)).unwrap_or_default();
                self.join_by_shared_reads(v_type, j_type, &mut v_layouts, &mut j_layouts, mate_layout_reads, &mut vdjs)?;
                self.join_by_word_hash(v_type, j_type, &mut v_layouts, &mut j_layouts, mate_layout_reads, &mut vdjs)?;
                if let Some(i) = vi {
                    one_sided[i].1 = v_layouts;
                }
                if let Some(i) = ji {
                    one_sided[i].1 = j_layouts;
                }
            }
        }
        crate::info!("built {} vdj by joining V with J layouts", vdjs.len());

        for (gene_type, layouts) in one_sided.iter_mut() {
            let mut kept = Vec::with_capacity(layouts.len());
            for layout in layouts.drain(..) {
                match self.try_complete_layout_with_blosum(*gene_type, &layout)? {
                    Some(v) => vdjs.push(v),
                    None => kept.push(layout),
                }
            }
            *layouts = kept;
        }
        for (gene_type, layouts) in &one_sided {
            for layout in layouts {
                if let Some(v) = self.try_create_one_sided_vdj(*gene_type, layout) {
                    vdjs.push(v);
                }
            }
        }
        vdjs.sort_by(|a, b| b.num_reads().cmp(&a.num_reads()));
        crate::info!("found {} vdj sequences before merge", vdjs.len());

        // grouped by (V gene type, J gene type), merged per group, the groups in key order (null first)
        let mut groups: Vec<((Option<VjGeneType>, Option<VjGeneType>), Vec<VdjSequence>)> = Vec::new();
        for v in vdjs {
            if v.v_anchor.is_none() && v.j_anchor.is_none() {
                continue;
            }
            let key = (v.v_anchor.as_ref().map(|a| a.gene_type), v.j_anchor.as_ref().map(|a| a.gene_type));
            match groups.iter_mut().find(|(k, _)| *k == key) {
                Some((_, l)) => l.push(v),
                None => groups.push((key, vec![v])),
            }
        }
        groups.sort_by(|a, b| a.0.cmp(&b.0));
        let min_bq = self.min_base_quality;
        let merged = parallel_map(groups.into_iter().map(|(_, l)| l).collect(), threads, |l| merge_identical(l, min_bq))?;
        let mut vdjs: Vec<VdjSequence> = merged.into_iter().flatten().collect();
        vdjs.sort_by(|a, b| b.num_reads().cmp(&a.num_reads()));
        crate::info!("{} vdj sequences after merge", vdjs.len());
        Ok(vdjs)
    }

    /// joinVjLayoutsBySharedReads: the J layouts that share the most read names with a V layout first
    fn join_by_shared_reads(
        &self,
        v_type: VjGeneType,
        j_type: VjGeneType,
        v_layouts: &mut Vec<Arc<ReadLayout>>,
        j_layouts: &mut Vec<Arc<ReadLayout>>,
        mates: &HashMap<String, Vec<LayoutRead>>,
        out: &mut Vec<VdjSequence>,
    ) -> Result<()> {
        let mut consumed = vec![false; j_layouts.len()];
        let mut by_read: HashMap<&str, Vec<usize>> = HashMap::new();
        for (ji, j) in j_layouts.iter().enumerate() {
            for r in j.reads() {
                by_read.entry(&r.key.name).or_default().push(ji);
            }
        }
        let mut kept = Vec::with_capacity(v_layouts.len());
        for v in v_layouts.drain(..) {
            let mut counts: HashMap<usize, i32> = HashMap::new();
            for r in v.reads() {
                if let Some(l) = by_read.get(&*r.key.name) {
                    for &ji in l {
                        *counts.entry(ji).or_insert(0) += 1;
                    }
                }
            }
            let joined = self.try_candidates(&v, j_layouts, counts, &consumed, v_type, j_type, mates)?;
            match joined {
                Some((ji, vdj)) => {
                    out.push(vdj);
                    consumed[ji] = true;
                }
                None => kept.push(v),
            }
        }
        *v_layouts = kept;
        let mut i = 0;
        j_layouts.retain(|_| {
            i += 1;
            !consumed[i - 1]
        });
        Ok(())
    }

    /// the candidate J layouts by descending count, then layout id; the first that overlaps
    #[allow(clippy::too_many_arguments)]
    fn try_candidates(
        &self,
        v: &Arc<ReadLayout>,
        j_layouts: &[Arc<ReadLayout>],
        counts: HashMap<usize, i32>,
        consumed: &[bool],
        v_type: VjGeneType,
        j_type: VjGeneType,
        mates: &HashMap<String, Vec<LayoutRead>>,
    ) -> Result<Option<(usize, VdjSequence)>> {
        let mut candidates: Vec<(usize, i32)> = counts.into_iter().filter(|(ji, _)| !consumed[*ji]).collect();
        candidates.sort_by(|a, b| (-a.1).cmp(&-b.1).then_with(|| j_layouts[a.0].id.cmp(&j_layouts[b.0].id)));
        for (ji, _) in candidates {
            if let Some(vdj) = self.try_overlap_vj(v, &j_layouts[ji], v_type, j_type, mates)? {
                return Ok(Some((ji, vdj)));
            }
        }
        Ok(None)
    }

    /// joinVjLayoutsByWordHash: 8-mer hashes of the V layouts after the anchor and the J layouts before it
    fn join_by_word_hash(
        &self,
        v_type: VjGeneType,
        j_type: VjGeneType,
        v_layouts: &mut Vec<Arc<ReadLayout>>,
        j_layouts: &mut Vec<Arc<ReadLayout>>,
        mates: &HashMap<String, Vec<LayoutRead>>,
        out: &mut Vec<VdjSequence>,
    ) -> Result<()> {
        let mut by_hash: HashMap<i32, Vec<usize>> = HashMap::new();
        for (ji, j) in j_layouts.iter().enumerate() {
            let Some(range) = self.anchor_range(Vj::J, j) else { continue };
            if range.1 < 0 {
                continue;
            }
            let seq = j.high_confidence_sequence(MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION);
            let end = (range.0 as usize).min(seq.len());
            for h in calc_sequence_word_hashes(&seq[..end], VJ_JOIN_HASH_WORD_SIZE)? {
                by_hash.entry(h).or_default().push(ji);
            }
        }
        let mut consumed = vec![false; j_layouts.len()];
        let mut kept = Vec::with_capacity(v_layouts.len());
        for v in v_layouts.drain(..) {
            let Some(range) = self.anchor_range(Vj::V, &v) else {
                kept.push(v);
                continue;
            };
            if range.1 < 0 {
                kept.push(v);
                continue;
            }
            let seq = v.high_confidence_sequence(MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION);
            let start = ((range.1 + 1).max(0) as usize).min(seq.len());
            let mut counts: HashMap<usize, i32> = HashMap::new();
            for h in calc_sequence_word_hashes(&seq[start..], VJ_JOIN_HASH_WORD_SIZE)? {
                if let Some(l) = by_hash.get(&h) {
                    for &ji in l {
                        *counts.entry(ji).or_insert(0) += 1;
                    }
                }
            }
            match self.try_candidates(&v, j_layouts, counts, &consumed, v_type, j_type, mates)? {
                Some((ji, vdj)) => {
                    out.push(vdj);
                    consumed[ji] = true;
                }
                None => kept.push(v),
            }
        }
        *v_layouts = kept;
        let mut i = 0;
        j_layouts.retain(|_| {
            i += 1;
            !consumed[i - 1]
        });
        Ok(())
    }

    /// tryOverlapVJ: a V and a J layout whose high confidence sequences overlap, as one VDJ
    pub fn try_overlap_vj(
        &self,
        v_layout: &Arc<ReadLayout>,
        j_layout: &Arc<ReadLayout>,
        v_type: VjGeneType,
        j_type: VjGeneType,
        mates: &HashMap<String, Vec<LayoutRead>>,
    ) -> Result<Option<VdjSequence>> {
        if !v_type.paired().contains(&j_type) {
            return Ok(None);
        }
        let v_seq = v_layout.high_confidence_sequence(MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION);
        let j_seq = j_layout.high_confidence_sequence(MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION);
        let Some(overlap) = find_sequence_overlap(&v_seq, &j_seq, self.min_overlapped_bases) else { return Ok(None) };
        let j_shift = (overlap.seq1_offset + v_layout.aligned_position) - (overlap.seq2_offset + j_layout.aligned_position);
        let mut combined = ReadLayout::merge(v_layout, j_layout, 0, j_shift, self.min_base_quality);
        combined.id = format!("{};{}", v_layout.id, j_layout.id);
        self.try_extend_combined_layout_with_mates(&mut combined, mates, v_type, j_type, j_shift);

        let (mut v0, mut v1) = self.adaptor.extrapolated_anchor_range(Vj::V, v_layout);
        let (mut j0, mut j1) = self.adaptor.extrapolated_anchor_range(Vj::J, j_layout);
        let v_shift = combined.aligned_position - v_layout.aligned_position;
        v0 += v_shift;
        v1 += v_shift;
        let j_anchor_shift = combined.aligned_position - j_layout.aligned_position - j_shift;
        j0 += j_anchor_shift;
        j1 += j_anchor_shift;
        let slice_end = (j1 + 1).min(combined.length());
        let mut v_boundary = v1 + 1;
        let mut j_boundary = j0;
        let len = combined.length();
        if v_boundary < MIN_ANCHOR_LENGTH_BASES || v_boundary > len - MIN_ANCHOR_LENGTH_BASES || j_boundary < MIN_ANCHOR_LENGTH_BASES || j_boundary > len - MIN_ANCHOR_LENGTH_BASES
        {
            return Ok(None);
        }
        let (slice_start, v_anchor) = if v1 < j0 {
            let s = v0.max(0);
            v_boundary -= s;
            j_boundary -= s;
            (s, Some(self.anchor_by_read_match(v_boundary, v_layout, v_type)))
        } else {
            (0, None)
        };
        let j_anchor = self.anchor_by_read_match(j_boundary, j_layout, j_type);
        Ok(Some(VdjSequence::new(Arc::new(combined), slice_start, slice_end, v_anchor, Some(j_anchor))))
    }

    /// tryCompleteLayoutWithBlosum: the missing anchor found by the BLOSUM search
    pub fn try_complete_layout_with_blosum(&self, gene_type: VjGeneType, layout: &Arc<ReadLayout>) -> Result<Option<VdjSequence>> {
        let targets = gene_type.paired();
        let seq = layout.consensus_sequence();
        let Some(range) = self.anchor_range(gene_type.vj(), layout) else { return Ok(None) };
        let (search_start, search_end) = if gene_type.vj() == Vj::V { (range.1 + 1, seq.len() as i32) } else { (0, range.0) };
        let Some(m) = self.searcher.search(seq, targets, SearchMode::DisallowNegSimilarity, search_start, search_end)? else { return Ok(None) };
        let ((mut v0, mut v1), (mut j0, mut j1)) =
            if gene_type.vj() == Vj::V { (range, (m.anchor_start, m.anchor_end - 1)) } else { ((m.anchor_start, m.anchor_end - 1), range) };
        let layout_start = v0.max(0);
        let layout_end = (j1 + 1).min(layout.length());
        v0 -= layout_start;
        v1 -= layout_start;
        j0 -= layout_start;
        j1 -= layout_start;
        let _ = (v0, j1);
        let first = m.template_genes.first().map(|&t| self.ds.template(t).gene_type).expect("java.util.NoSuchElementException: List is empty.");
        let blosum = |vj: Vj, boundary: i32| VjAnchor {
            vj,
            gene_type: first,
            anchor_boundary: boundary,
            template_anchor_seq: m.template_anchor_seq.clone(),
            kind: AnchorKind::Blosum { template_genes: m.template_genes.clone(), similarity_score: m.similarity_score },
        };
        if m.template_genes.iter().any(|&t| self.ds.template(t).gene_type != first) {
            bail!("VJAnchorByBlosum: gene types size != 1");
        }
        let (v_anchor, j_anchor) = if gene_type.vj() == Vj::V {
            (self.anchor_by_read_match(v1 + 1, layout, gene_type), blosum(Vj::J, j0))
        } else {
            (blosum(Vj::V, v1 + 1), self.anchor_by_read_match(j0, layout, gene_type))
        };
        Ok(Some(VdjSequence::new(layout.clone(), layout_start, layout_end, Some(v_anchor), Some(j_anchor))))
    }

    /// tryCreateOneSidedVdj: 60 bases past the V anchor, or before the J anchor
    pub fn try_create_one_sided_vdj(&self, gene_type: VjGeneType, layout: &Arc<ReadLayout>) -> Option<VdjSequence> {
        let range = self.anchor_range(gene_type.vj(), layout)?;
        Some(if gene_type.vj() == Vj::V {
            let start = range.0.max(0);
            let end = layout.length().min(range.1 + 1 + PARTIAL_VDJ_UNANCHORED_LENGTH_BASES);
            let v = self.anchor_by_read_match(range.1 + 1 - start, layout, gene_type);
            VdjSequence::new(layout.clone(), start, end, Some(v), None)
        } else {
            let start = (range.0 - PARTIAL_VDJ_UNANCHORED_LENGTH_BASES).max(0);
            let end = layout.length().min(range.1 + 1);
            let j = self.anchor_by_read_match(range.0 - start, layout, gene_type);
            VdjSequence::new(layout.clone(), start, end, None, Some(j))
        })
    }

    /// tryExtendCombinedLayoutWithMates: the mates of the layout's reads that belong to its V or J gene
    /// type, J reads moved by the J shift
    fn try_extend_combined_layout_with_mates(
        &self,
        combined: &mut ReadLayout,
        mates: &HashMap<String, Vec<LayoutRead>>,
        v_type: VjGeneType,
        j_type: VjGeneType,
        j_shift: i32,
    ) {
        let mut reads = Vec::new();
        for r in combined.reads() {
            let Some(m) = mates.get(&*r.key.name) else { continue };
            for mr in m {
                let mc = candidate(mr);
                if mc.gene_type == v_type || mc.gene_type == j_type {
                    reads.push(if mc.gene_type.vj() == Vj::J { mr.copy(mr.aligned_position + j_shift) } else { mr.clone() });
                }
            }
        }
        self.try_extend_layout_by_sequence(combined, &reads);
    }

    /// tryExtendLayoutBySequence: adds, round by round, the reads that overlap the layout by 20 bases with
    /// at most one mismatch
    fn try_extend_layout_by_sequence(&self, layout: &mut ReadLayout, reads: &[LayoutRead]) {
        let mut remaining: Vec<usize> = (0..reads.len()).collect();
        loop {
            let candidates: Vec<(usize, i32)> = remaining
                .iter()
                .filter_map(|&i| {
                    let r = &reads[i];
                    let start = layout.aligned_position - r.aligned_position;
                    let end = start + r.read_length();
                    let overlap = end.min(layout.length()) - start.max(0);
                    if overlap >= LAYOUT_MIN_READ_OVERLAP_BASES { Some((i, start)) } else { None }
                })
                .collect();
            let mut mismatches = vec![0i32; candidates.len()];
            let consensus = layout.consensus_sequence();
            for (ci, &(i, start)) in candidates.iter().enumerate() {
                let seq = &reads[i].sequence;
                // the consensus positions the read covers
                let from = start.max(0);
                let to = (start + seq.len() as i32).min(consensus.len() as i32);
                for position in from..to {
                    let rb = seq[(position - start) as usize];
                    let cb = consensus[position as usize];
                    if rb != cb && rb != N {
                        mismatches[ci] += 1;
                    }
                }
            }
            let matching: Vec<usize> = candidates.iter().zip(&mismatches).filter(|(_, m)| **m <= LAYOUT_MATE_EXTEND_MISMATCHES_MAX).map(|(c, _)| c.0).collect();
            if matching.is_empty() {
                break;
            }
            for i in matching {
                remaining.retain(|&x| x != i);
                layout.add_read(reads[i].clone(), self.min_base_quality);
            }
        }
    }
}

/// mergeIdentical for one gene type pair: identical sequences merged, then those one base apart (with one
/// read at the difference) dropped
fn merge_identical(list: Vec<VdjSequence>, min_base_quality: u8) -> Result<Vec<VdjSequence>> {
    let mut out = list;
    let (mut i, mut j) = (0usize, 1usize);
    while i < out.len() && j < out.len() {
        let identical = {
            let mut diffs = 0;
            vdj_sequence_compare(&out[i], &out[j], &mut |a, b| {
                if a.1 >= 1 && b.1 >= 1 {
                    diffs += 1;
                }
                diffs == 0
            })
        };
        if identical {
            let merged = merge_vdjs(&out[i], &out[j], min_base_quality)?;
            out[i] = merged;
            out.remove(j);
        } else {
            j += 1;
        }
        if j >= out.len() {
            i += 1;
            j = i + 1;
        }
    }
    let (mut i, mut j) = (0usize, 1usize);
    while i < out.len() && j < out.len() {
        let similar = {
            let mut diffs = 0;
            let mut max_reads = 0;
            vdj_sequence_compare(&out[i], &out[j], &mut |a, b| {
                if a.1 >= 1 && b.1 >= 1 {
                    diffs += 1;
                    max_reads = max_reads.max(a.1.min(b.1));
                }
                diffs <= 1 && max_reads < 2
            })
        };
        if similar {
            out.remove(j);
        } else {
            j += 1;
        }
        if j >= out.len() {
            i += 1;
            j = i + 1;
        }
    }
    Ok(out)
}

/// runs `f` on each item on up to `threads` threads; the results in item order
pub fn parallel_map<T: Send, R: Send>(items: Vec<T>, threads: usize, f: impl Fn(T) -> Result<R> + Sync) -> Result<Vec<R>> {
    let n = items.len();
    if threads <= 1 || n <= 1 {
        return items.into_iter().map(f).collect();
    }
    let queue = std::sync::Mutex::new(items.into_iter().enumerate().collect::<std::collections::VecDeque<_>>());
    let results: std::sync::Mutex<Vec<Option<Result<R>>>> = std::sync::Mutex::new((0..n).map(|_| None).collect());
    std::thread::scope(|s| {
        for _ in 0..threads.min(n) {
            s.spawn(|| loop {
                let next = queue.lock().unwrap().pop_front();
                let Some((i, item)) = next else { break };
                let r = f(item);
                results.lock().unwrap()[i] = Some(r);
            });
        }
    });
    results.into_inner().unwrap().into_iter().map(|r| r.expect("result")).collect()
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::layout::tests::layout_read;

    pub const MIN_BQ: u8 = 30;

    /// TestUtils.createLayout: one read of high quality
    pub fn create_layout(seq: &str, aligned_position: i32) -> ReadLayout {
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(1);
        let id = format!("createLayout::autoReadId::{}", NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed));
        let q = ((MIN_BQ * 2 + 33) as char).to_string().repeat(seq.len());
        let mut l = ReadLayout::new();
        l.add_read(layout_read(&id, seq, &q, aligned_position, true), MIN_BQ);
        l
    }

    fn read_match(vj: Vj, gene_type: VjGeneType, boundary: i32, template: &str, num_reads: i32) -> VjAnchor {
        VjAnchor {
            vj,
            gene_type,
            anchor_boundary: boundary,
            template_anchor_seq: template.into(),
            kind: AnchorKind::ReadMatch { match_method: "exact".into(), num_reads },
        }
    }

    /// TestUtils.createVDJ(layout, v, j, start, end)
    pub fn create_vdj_layout(layout: Arc<ReadLayout>, v: Option<i32>, j: Option<i32>, start: i32, end: i32) -> VdjSequence {
        let cons = layout.consensus_string();
        let n = layout.reads().len() as i32;
        let va = v.map(|b| read_match(Vj::V, VjGeneType::Ighv, b, &cons[start as usize..start as usize + b as usize], n));
        let ja = j.map(|b| read_match(Vj::J, VjGeneType::Ighj, b, &cons[b as usize..], n));
        VdjSequence::new(layout, start, end, va, ja)
    }

    /// TestUtils.createVDJ(sequence, v, j)
    pub fn create_vdj(seq: &str, v: Option<i32>, j: Option<i32>) -> VdjSequence {
        let layout = Arc::new(create_layout(seq, 0));
        let len = layout.length();
        create_vdj_layout(layout, v, j, 0, len)
    }

    // VdjBuilderUtilsTest
    #[test]
    fn word_hashes() {
        let h = calc_sequence_word_hashes(b"ATCG", 2).unwrap();
        // AT 0b0001, TC 0b0110, CG 0b1011
        assert_eq!(h, vec![1, 6, 11]);
        let h = calc_sequence_word_hashes(b"ATNCGA", 2).unwrap();
        assert_eq!(h, vec![1, 11, 12]);
        assert!(calc_sequence_word_hashes(b"A", 17).is_err());
    }

    #[test]
    fn sequence_overlap() {
        // seq2 starts 10 bases into seq1 (i = -10: seq2 offset 10)
        let s1 = b"AAAAACCCCCGGGGGTTTTTACGTACGTAC";
        let s2 = b"GGGGGTTTTTACGTACGTACTTTTTTTTTT";
        let o = find_sequence_overlap(s1, s2, 10).unwrap();
        assert_eq!(o, SequenceOverlap { seq1_offset: 0, seq2_offset: 10, overlap_bases: 20, high_qual_match_bases: 20 });
        // and the other way round
        let o = find_sequence_overlap(s2, s1, 10).unwrap();
        assert_eq!(o, SequenceOverlap { seq1_offset: 10, seq2_offset: 0, overlap_bases: 20, high_qual_match_bases: 20 });
        // more matched bases needed than there are
        assert!(find_sequence_overlap(s1, s2, 20).is_none());
        // N matches anything
        let s3 = b"GGGGGTTTNTACGTACGTACTTTTTTTTTT";
        let o = find_sequence_overlap(s1, s3, 10).unwrap();
        assert_eq!(o.high_qual_match_bases, 19);
    }

    #[test]
    fn merge_anchor_counts() {
        let a = read_match(Vj::V, VjGeneType::Ighv, 10, "A", 2);
        let b = read_match(Vj::V, VjGeneType::Ighv, 12, "B", 5);
        let fixed = compat::with_java(false, || merge_anchors(Some(&a), Some(&b), 10, 12, 12).unwrap());
        assert_eq!((fixed.template_anchor_seq.as_str(), fixed.num_reads(), fixed.anchor_boundary), ("B", Some(7), 12));
        let java = compat::with_java(true, || merge_anchors(Some(&a), Some(&b), 10, 12, 12).unwrap());
        assert_eq!((java.template_anchor_seq.as_str(), java.num_reads()), ("A", Some(7)));
    }

    #[test]
    fn vdj_cdr3() {
        // V anchor of 9 bases, J anchor at 21
        let seq = "TGTGCGAGATTTCCCGGGTTTTGGGGCCAGGG";
        let v = create_vdj(seq, Some(9), Some(21));
        assert_eq!(std::str::from_utf8(v.cdr3_sequence()).unwrap(), &seq[6..24]);
        assert!(v.is_in_frame());
        assert_eq!(v.v_anchor_sequence().unwrap(), &seq.as_bytes()[..9]);
        assert_eq!(v.j_anchor_sequence().unwrap(), &seq.as_bytes()[21..]);
    }
}
