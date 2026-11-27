//! VJReadLayoutBuilder: read candidates to layout reads (trimmed slices, aligned at the anchor boundary)
//! and the layouts of one gene type; ReadSlice.

use std::sync::Arc;

use anyhow::{bail, Result};
use hmf_common::codon::{reverse_complement, swap_dna_base};
use hmf_common::kotlin_random::KotlinRandom;

use crate::compat;
use crate::constants::*;
use crate::genes::{Vj, VjGeneType};
use crate::layout::{ForestRead, LayoutForest, LayoutRead, ReadKey, ReadLayout};
use crate::read::CiderRead;
use crate::screener::{MatchMethod, VjReadCandidate};

/// ReadSlice: a read, reverse complemented if asked, then cut to [start, end)
#[derive(Clone, Debug)]
pub struct ReadSlice {
    pub read: Arc<CiderRead>,
    pub reverse_complement: bool,
    pub slice_start: i32,
    pub slice_end: i32,
}

#[allow(dead_code)]
impl ReadSlice {
    pub fn new(read: Arc<CiderRead>, reverse_complement: bool, slice_start: i32, slice_end: i32) -> ReadSlice {
        assert!(slice_end > slice_start, "slice End <= sliceStart");
        ReadSlice { read, reverse_complement, slice_start, slice_end }
    }
    pub fn read_length(&self) -> i32 {
        self.slice_end - self.slice_start
    }
    pub fn read_string(&self) -> Vec<u8> {
        let s = self.read.read_string();
        let s = if self.reverse_complement { reverse_complement(s) } else { s.to_vec() };
        s[self.slice_start as usize..self.slice_end as usize].to_vec()
    }
    pub fn base_qualities(&self) -> Vec<u8> {
        let mut q = self.read.quals.clone();
        if self.reverse_complement {
            q.reverse();
        }
        q[self.slice_start as usize..self.slice_end as usize].to_vec()
    }
    pub fn base_at(&self, slice_pos: i32) -> u8 {
        let b = self.read.read_string()[self.slice_position_to_read_position(slice_pos) as usize];
        if self.reverse_complement { swap_dna_base(b) } else { b }
    }
    pub fn read_position_to_slice_position(&self, read_pos: i32) -> i32 {
        let pos = if self.reverse_complement { self.read.read_length() - read_pos - 1 } else { read_pos };
        pos - self.slice_start
    }
    pub fn slice_position_to_read_position(&self, slice_pos: i32) -> i32 {
        let pos = slice_pos + self.slice_start;
        if self.reverse_complement { self.read.read_length() - pos - 1 } else { pos }
    }
    /// readRangeToSliceRange: (start, end exclusive)
    pub fn read_range_to_slice_range(&self, start: i32, end_exclusive: i32) -> (i32, i32) {
        let len = self.read.read_length();
        let s = if self.reverse_complement { len - end_exclusive } else { start };
        let e = if self.reverse_complement { len - start } else { end_exclusive };
        (s - self.slice_start, e - self.slice_start)
    }
}

/// VjLayoutRead's own fields: the slice and the candidate it comes from
#[derive(Debug)]
pub struct VjLayoutSource {
    pub slice: ReadSlice,
    pub candidate: Arc<VjReadCandidate>,
}

/// the candidate of a layout read (IVJReadLayoutAdaptor.toReadCandidate)
pub fn candidate(read: &LayoutRead) -> &VjReadCandidate {
    &read.source.as_ref().expect("layout read without a candidate").candidate
}

/// CiderUtils.numTrailingPolyG
pub fn num_trailing_poly_g(seq: &[u8], slice_end: i32) -> i32 {
    for i in 0..slice_end {
        if seq[(slice_end - i - 1) as usize] != b'G' {
            return i;
        }
    }
    slice_end.max(0)
}

/// CiderUtils.numLeadingPolyC
pub fn num_leading_poly_c(seq: &[u8], slice_start: i32) -> i32 {
    for i in slice_start.max(0)..seq.len() as i32 {
        if seq[i as usize] != b'C' {
            return i - slice_start;
        }
    }
    seq.len() as i32 - slice_start
}

/// IVJReadLayoutAdaptor: what the VDJ builder asks of a layout (the tests give their own answers)
pub trait LayoutAdaptor: Sync {
    /// getExtrapolatedAnchorRange: (first, last) inclusive, possibly outside the layout or empty
    fn extrapolated_anchor_range(&self, vj: Vj, layout: &ReadLayout) -> (i32, i32);
    fn anchor_match_method(&self, layout: &ReadLayout) -> MatchMethod;
    fn template_anchor_sequence(&self, layout: &ReadLayout) -> String;
    /// getAnchorRange: the part of the anchor inside the layout, (first, last) inclusive
    fn anchor_range(&self, vj: Vj, layout: &ReadLayout) -> Option<(i32, i32)> {
        let (first, last) = self.extrapolated_anchor_range(vj, layout);
        if first >= layout.length() || last < 0 {
            return None;
        }
        Some((first.max(0), (last + 1).min(layout.length()) - 1))
    }
}

impl LayoutAdaptor for VjReadLayoutBuilder {
    fn extrapolated_anchor_range(&self, vj: Vj, layout: &ReadLayout) -> (i32, i32) {
        VjReadLayoutBuilder::extrapolated_anchor_range(vj, layout)
    }
    fn anchor_match_method(&self, layout: &ReadLayout) -> MatchMethod {
        VjReadLayoutBuilder::anchor_match_method(layout)
    }
    fn template_anchor_sequence(&self, layout: &ReadLayout) -> String {
        VjReadLayoutBuilder::template_anchor_sequence(layout)
    }
}

#[allow(dead_code)]
pub struct LayoutBuildResult {
    pub layouts: Vec<Arc<ReadLayout>>,
    pub num_candidate_reads: usize,
    pub num_high_qual_reads: usize,
    pub num_reads_used: usize,
    pub down_sampled: bool,
}

/// VJReadLayoutBuilder
pub struct VjReadLayoutBuilder {
    pub trim_bases: i32,
    pub min_base_quality: i32,
}

impl VjReadLayoutBuilder {
    pub fn new(trim_bases: i32, min_base_quality: i32) -> Self {
        VjReadLayoutBuilder { trim_bases, min_base_quality }
    }

    /// determineReadSlice: adapter trim, constant trim and poly-G tail trim, then the reverse complement
    fn determine_read_slice(&self, read: &Arc<CiderRead>, use_rev: bool) -> Option<ReadSlice> {
        let len = read.read_length();
        let mut start = read.adapter_trim.0 + self.trim_bases;
        let mut end = len - read.adapter_trim.1 - self.trim_bases;
        let s = read.read_string();
        if !read.negative_strand() {
            let n = num_trailing_poly_g(s, end);
            if n >= MIN_POLY_G_TRIM_COUNT {
                end -= n + POLY_G_TRIM_EXTRA_BASE_COUNT;
            }
        } else {
            let n = num_leading_poly_c(s, start);
            if n >= MIN_POLY_G_TRIM_COUNT {
                start += n + POLY_G_TRIM_EXTRA_BASE_COUNT;
            }
        }
        if end - start < 5 {
            return None;
        }
        if use_rev {
            let t = start;
            start = len - end;
            end = len - t;
        }
        Some(ReadSlice::new(read.clone(), use_rev, start, end))
    }

    /// readCandidateToLayoutRead: aligned at the last base of a V anchor or the first of a J anchor
    pub fn candidate_to_layout_read(&self, c: &Arc<VjReadCandidate>) -> Option<LayoutRead> {
        let slice = self.determine_read_slice(&c.read, c.use_reverse_complement)?;
        let aligned = if c.gene_type.vj() == Vj::V { c.anchor_offset_end - 1 - slice.slice_start } else { c.anchor_offset_start - slice.slice_start };
        let sequence: Arc<[u8]> = Arc::from(slice.read_string().as_slice());
        let quals: Arc<[u8]> = Arc::from(slice.base_qualities().as_slice());
        if sequence.len() != quals.len() {
            panic!("java.lang.IllegalArgumentException: sequence.length != baseQualities.size");
        }
        let key = ReadKey { name: Arc::from(c.read.name.as_str()), first_of_pair: c.read.first_in_pair() };
        Some(LayoutRead { key, sequence, base_qualities: quals, aligned_position: aligned, source: Some(Arc::new(VjLayoutSource { slice, candidate: c.clone() })) })
    }

    /// getAnchorMatchMethod: the first read's
    pub fn anchor_match_method(layout: &ReadLayout) -> MatchMethod {
        candidate(layout.reads().first().expect("read candidate list is empty")).match_method
    }

    /// getTemplateAnchorSequence: of the first read with the highest similarity score
    pub fn template_anchor_sequence(layout: &ReadLayout) -> String {
        let mut best: Option<&VjReadCandidate> = None;
        for r in layout.reads() {
            let c = candidate(r);
            if best.is_none_or(|b| b.similarity_score < c.similarity_score) {
                best = Some(c);
            }
        }
        best.expect("read candidate list is empty").template_anchor_sequence.clone()
    }

    /// getExtrapolatedAnchorRange: (first, last) inclusive, possibly outside the layout or empty
    pub fn extrapolated_anchor_range(vj: Vj, layout: &ReadLayout) -> (i32, i32) {
        let anchor_length = layout.reads().iter().map(|r| {
            let c = candidate(r);
            c.anchor_offset_end - c.anchor_offset_start
        });
        let anchor_length = anchor_length.max().unwrap_or(0);
        match vj {
            Vj::V => (layout.aligned_position - anchor_length + 1, layout.aligned_position),
            Vj::J => (layout.aligned_position, layout.aligned_position + anchor_length - 1),
        }
    }

    /// getAnchorRange: the part of the anchor inside the layout, (first, last) inclusive
#[allow(dead_code)]
    pub fn anchor_range(vj: Vj, layout: &ReadLayout) -> Option<(i32, i32)> {
        let (first, last) = Self::extrapolated_anchor_range(vj, layout);
        if first >= layout.length() || last < 0 {
            return None;
        }
        Some((first.max(0), (last + 1).min(layout.length()) - 1))
    }

    /// buildLayouts for one gene type
    pub fn build_layouts(
        &self,
        gene_type: VjGeneType,
        candidates: &[Arc<VjReadCandidate>],
        min_matched_bases: i32,
        max_read_count_per_gene: i32,
        max_low_qual_base_fraction: f64,
    ) -> Result<LayoutBuildResult> {
        let mut reads: Vec<LayoutRead> = Vec::new();
        for c in candidates {
            let Some(r) = self.candidate_to_layout_read(c) else { continue };
            let low = r.base_qualities.iter().filter(|&&q| (q as i32) < self.min_base_quality).count();
            let fraction = low as f64 / r.read_length() as f64;
            // Doubles.lessOrEqual
            if fraction - max_low_qual_base_fraction < 1e-10 {
                reads.push(r);
            }
        }
        // descending by aligned position, mean quality, end, name, then second of pair before first
        let avg = |r: &LayoutRead| -> f64 {
            let mut sum = 0.0;
            for &q in r.base_qualities.iter() {
                sum += (q as i8) as f64;
            }
            if r.base_qualities.is_empty() { f64::NAN } else { sum / r.base_qualities.len() as f64 }
        };
        let mut keyed: Vec<(f64, LayoutRead)> = reads.into_iter().map(|r| (avg(&r), r)).collect();
        keyed.sort_by(|(qa, a), (qb, b)| {
            b.aligned_position
                .cmp(&a.aligned_position)
                .then_with(|| qb.total_cmp(qa))
                .then_with(|| (b.aligned_position + b.read_length()).cmp(&(a.aligned_position + a.read_length())))
                .then_with(|| b.key.name.as_bytes().cmp(a.key.name.as_bytes()))
                .then_with(|| (!b.key.first_of_pair).cmp(&!a.key.first_of_pair))
        });
        let mut reads: Vec<LayoutRead> = keyed.into_iter().map(|(_, r)| r).collect();
        let num_high_qual_reads = reads.len();
        // CI5: Java downsamples when the count reaches the limit (>=), though nothing is dropped at the limit
        // and the option says reads are downsampled when they exceed it. Corrected: above the limit.
        let over = if compat::java("CI5") { reads.len() >= max_read_count_per_gene as usize } else { reads.len() > max_read_count_per_gene as usize };
        let down_sampled = if over {
            let keep = max_read_count_per_gene as f64 / reads.len() as f64;
            crate::info!(
                "building {} layouts, read count({}) > limit({}), downsampling(frac to keep: {})",
                gene_type.name(),
                reads.len(),
                max_read_count_per_gene,
                hmf_common::fmt::java_fixed(keep, 3)
            );
            let mut random = KotlinRandom::new(0);
            reads.retain(|_| random.next_double() < keep);
            true
        } else {
            false
        };
        crate::info!("building {} layouts from {} reads", gene_type.name(), reads.len());
        let mut forest = LayoutForest::new(self.min_base_quality as u8, min_matched_bases, LAYOUT_MIN_SUPPORT_TO_SEAL_NODE);
        for (i, r) in reads.iter().enumerate() {
            forest.try_add_read(ForestRead { source: i, sequence: r.sequence.clone(), base_qualities: r.base_qualities.clone(), aligned_position: r.aligned_position });
        }
        let mut layouts = forest.build_read_layouts(&|fr: &ForestRead| reads[fr.source].clone());
        // sortedByDescending { reads.size }, stable
        layouts.sort_by(|a, b| b.reads().len().cmp(&a.reads().len()));
        // readCountSanityCheck
        let mut keys: std::collections::HashSet<&ReadKey> = reads.iter().map(|r| &r.key).collect();
        for l in &layouts {
            for r in l.reads() {
                keys.remove(&r.key);
            }
        }
        if !keys.is_empty() {
            crate::error!("built layouts are missing {} input reads", keys.len());
            bail!("built layouts are missing {} input reads", keys.len());
        }
        crate::info!("built {} {} layouts from {} reads", layouts.len(), gene_type.name(), reads.len());
        Ok(LayoutBuildResult {
            layouts: layouts.into_iter().map(Arc::new).collect(),
            num_candidate_reads: candidates.len(),
            num_high_qual_reads,
            num_reads_used: reads.len(),
            down_sampled,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::test_util::test_read;
    use crate::read::{FLAG_FIRST, FLAG_PAIRED, FLAG_REVERSE};

    // TestUtils.createReadCandidate
    pub fn read_candidate(seq: &str, negative: bool, use_rev: bool, vj: Vj, start: i32, end: i32) -> Arc<VjReadCandidate> {
        let mut r = test_read();
        r.flags = FLAG_PAIRED | FLAG_FIRST | if negative { FLAG_REVERSE } else { 0 };
        r.bases = seq.as_bytes().to_vec();
        r.quals = vec![37; seq.len()];
        let gene_type = if vj == Vj::V { VjGeneType::Trav } else { VjGeneType::Traj };
        Arc::new(VjReadCandidate {
            read: Arc::new(r),
            templates: Vec::new(),
            gene_type,
            template_anchor_sequence: "CACGTG".into(),
            match_method: MatchMethod::Align,
            use_reverse_complement: use_rev,
            anchor_offset_start: start,
            anchor_offset_end: end,
            left_soft_clip: 0,
            right_soft_clip: 0,
            similarity_score: i32::MIN,
        })
    }

    fn check(builder: &VjReadLayoutBuilder, seq: &str, negative: bool, use_rev: bool, vj: Vj, start: i32, end: i32) {
        let r = builder.candidate_to_layout_read(&read_candidate(seq, negative, use_rev, vj, start, end)).unwrap();
        assert_eq!(std::str::from_utf8(&r.sequence).unwrap(), "AAGACACGGC");
        // aligned at the last base of a V anchor, the first of a J anchor
        assert_eq!(r.aligned_position, 5);
    }

    // VjReadLayoutBuilderTest
    #[test]
    fn candidate_to_layout_read_pos_strand() {
        let b = VjReadLayoutBuilder::new(0, 30);
        check(&b, "AAGACACGGC", false, false, Vj::V, 2, 6);
        check(&b, "AAGACACGGC", false, false, Vj::J, 5, 8);
    }

    #[test]
    fn candidate_to_layout_read_neg_strand() {
        let b = VjReadLayoutBuilder::new(0, 30);
        check(&b, "GCCGTGTCTT", false, true, Vj::V, 2, 6);
        check(&b, "GCCGTGTCTT", false, true, Vj::J, 5, 8);
    }

    #[test]
    fn candidate_to_layout_read_poly_g() {
        let b = VjReadLayoutBuilder::new(0, 30);
        let seq = format!("AAGACACGGCGTACT{}", "G".repeat(8));
        check(&b, &seq, false, false, Vj::V, 2, 6);
        check(&b, &seq, false, false, Vj::J, 5, 8);
    }

    #[test]
    fn candidate_to_layout_read_poly_c() {
        let b = VjReadLayoutBuilder::new(0, 30);
        let seq = format!("{}GTACTAAGACACGGC", "C".repeat(8));
        check(&b, &seq, true, false, Vj::V, 2 + 13, 6 + 13);
        check(&b, &seq, true, false, Vj::J, 5 + 13, 8 + 13);
    }

    #[test]
    fn candidate_to_layout_read_trim_bases_poly_g() {
        let b = VjReadLayoutBuilder::new(1, 30);
        let seq = format!("TAAGACACGGCGTACT{}A", "G".repeat(8));
        check(&b, &seq, false, false, Vj::V, 2 + 1, 6 + 1);
        check(&b, &seq, false, false, Vj::J, 5 + 1, 8 + 1);
    }

    #[test]
    fn poly_g_counts() {
        assert_eq!(num_trailing_poly_g(b"GGG", 3), 3);
        assert_eq!(num_trailing_poly_g(b"AGG", 3), 2);
        assert_eq!(num_leading_poly_c(b"CCC", 0), 3);
        assert_eq!(num_leading_poly_c(b"ACC", 1), 2);
    }
}
