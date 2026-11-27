//! VdjAnnotator: the VDJ sequences in output order with their filters, read counts and anchor scores,
//! and CiderApplication.passesHardFilters.

use std::collections::HashMap;

use anyhow::{anyhow, Result};
use hmf_common::codon::amino_acid_from_bases;

use crate::alignment::{shm_status, AlignmentAnnotation, AlignmentStatus};
use crate::blosum::similarity_score;
use crate::constants::*;
use crate::genes::{IgTcrLocus, Vj};
use crate::layout_builder::candidate;
use crate::read::{Cigar, CigarOp};
use crate::vdj::{AnchorKind, VdjSequence, VjAnchor};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Filter {
    Pass,
    MatchesRef,
    NoVAnchor,
    NoJAnchor,
    PoorVAnchor,
    PoorJAnchor,
    Partial,
    NoVdjAlignment,
    Duplicate,
    MinLength,
    MaxLength,
    Cdr3Deleted,
    NoHighQualSupport,
    VJOverlap,
}

impl Filter {
    pub fn name(self) -> &'static str {
        match self {
            Filter::Pass => "PASS",
            Filter::MatchesRef => "MATCHES_REF",
            Filter::NoVAnchor => "NO_V_ANCHOR",
            Filter::NoJAnchor => "NO_J_ANCHOR",
            Filter::PoorVAnchor => "POOR_V_ANCHOR",
            Filter::PoorJAnchor => "POOR_J_ANCHOR",
            Filter::Partial => "PARTIAL",
            Filter::NoVdjAlignment => "NO_VDJ_ALIGNMENT",
            Filter::Duplicate => "DUPLICATE",
            Filter::MinLength => "MIN_LENGTH",
            Filter::MaxLength => "MAX_LENGTH",
            Filter::Cdr3Deleted => "CDR3_DELETED",
            Filter::NoHighQualSupport => "NO_HIGH_QUAL_SUPPORT",
            Filter::VJOverlap => "V_J_OVERLAP",
        }
    }
}

/// VdjAnnotation; `vdj` is the index of the sequence, `alignment` that of its alignment annotation
pub struct VdjAnnotation {
    pub vdj: usize,
    pub filters: Vec<Filter>,
    pub v_aligned_reads: i32,
    pub j_aligned_reads: i32,
    pub v_non_split_reads: i32,
    pub j_non_split_reads: i32,
    pub cdr3_support_min: i32,
    pub v_similarity_score: Option<i32>,
    pub j_similarity_score: Option<i32>,
    pub v_primer_match_count: i32,
    pub j_primer_match_count: i32,
    pub alignment: Option<usize>,
    pub somatic_hypermutation_status: Option<&'static str>,
}

impl VdjAnnotation {
    pub fn passes_filter(&self) -> bool {
        self.filters.contains(&Filter::Pass)
    }
}

pub fn locus(vdj: &VdjSequence) -> IgTcrLocus {
    match &vdj.v_anchor {
        Some(v) => v.gene_type.locus(),
        None => vdj.j_anchor.as_ref().expect("an anchor").gene_type.locus(),
    }
}

/// calcCdr3SupportMin: the least high quality support over the whole sequence
pub fn cdr3_support_min(vdj: &VdjSequence) -> i32 {
    vdj.support_counts().into_iter().min().unwrap_or(0)
}

/// CiderUtils.getAdjustedAlignmentBlocks: (read start 1-based, length); blocks run across I, D and P
pub fn adjusted_alignment_blocks(cigar: &Cigar) -> Vec<(i32, i32)> {
    let mut blocks = Vec::new();
    let mut read_base = 1;
    let mut start = -1;
    for &(op, len) in cigar {
        let len = len as i32;
        match op {
            CigarOp::P | CigarOp::D => {}
            CigarOp::I => read_base += len,
            CigarOp::M | CigarOp::Eq | CigarOp::X => {
                if start == -1 {
                    start = read_base;
                }
                read_base += len;
            }
            CigarOp::S | CigarOp::H | CigarOp::N => {
                if start != -1 {
                    blocks.push((start, read_base - start));
                    start = -1;
                }
                if op == CigarOp::S {
                    read_base += len;
                }
            }
        }
    }
    if start != -1 {
        blocks.push((start, read_base - start));
    }
    blocks
}

/// countNonSplitReads: reads with an aligned block over the anchor boundary +- 30 bases
pub fn count_non_split_reads(vdj: &VdjSequence, vj: Vj) -> i32 {
    let aligned = vdj.layout.aligned_position - vdj.layout_slice_start;
    let boundary = match vj {
        Vj::V => match vdj.v_anchor_boundary() {
            Some(b) => b,
            None => return 0,
        },
        Vj::J => match vdj.j_anchor_boundary() {
            Some(b) => b,
            None => return 0,
        },
    };
    let region_start = (boundary - MIN_NON_SPLIT_READ_STRADDLE_LENGTH).max(0);
    let region_end = (boundary + MIN_NON_SPLIT_READ_STRADDLE_LENGTH).min(vdj.length());
    let mut count = 0;
    for read in vdj.layout.reads() {
        let source = read.source.as_ref().expect("layout read source");
        let record = &source.candidate.read;
        let pos_in_vdj = aligned - read.aligned_position;
        let slice_start = region_start - pos_in_vdj;
        let slice_end = region_end - pos_in_vdj;
        if !record.unmapped() {
            for (start, len) in adjusted_alignment_blocks(&record.cigar) {
                let (l, r) = source.slice.read_range_to_slice_range(start - 1, start - 1 + len);
                assert!(l < r, "java.lang.IllegalArgumentException: Failed requirement (align range)");
                if slice_start >= l && slice_end <= r {
                    count += 1;
                }
            }
        }
    }
    count
}

fn anchor_reads(a: &Option<VjAnchor>) -> i32 {
    match a {
        Some(VjAnchor { kind: AnchorKind::ReadMatch { num_reads, .. }, .. }) => *num_reads,
        _ => 0,
    }
}

/// vdjInKnownMatchesRefList
fn in_known_matches_ref_list(vdj: &VdjSequence) -> bool {
    amino_acid_from_bases(vdj.cdr3_sequence()).starts_with(MATCHES_REF_KNOWN_CDR3_AA_PREFIX)
}

fn matches_ref(vdj: &VdjSequence, v_aligned: i32, j_aligned: i32, v_non_split: i32, j_non_split: i32) -> bool {
    let max_non_split = MAX_NONSPLIT_READS.min((vdj.num_reads() / 2).max(1));
    ((v_aligned == 0 || j_aligned == 0) && (j_non_split + v_non_split) >= max_non_split) || in_known_matches_ref_list(vdj)
}

/// vdjMatchesRef(vdj): whether the sequence looks like the reference, so needs no alignment
pub fn vdj_matches_ref(vdj: &VdjSequence) -> bool {
    matches_ref(vdj, anchor_reads(&vdj.v_anchor), anchor_reads(&vdj.j_anchor), count_non_split_reads(vdj, Vj::V), count_non_split_reads(vdj, Vj::J))
}

/// calcAnchorSimilarity
fn anchor_similarity(vdj: &VdjSequence, anchor: &VjAnchor) -> Result<i32> {
    let seq = if anchor.vj == Vj::V { vdj.v_anchor_sequence() } else { vdj.j_anchor_sequence() }.unwrap_or_default();
    similarity_score(anchor.vj, anchor.template_anchor_seq.as_bytes(), seq)
        .map_err(|_| anyhow!("cannot calc similarity score: {} and {}", anchor.template_anchor_seq, String::from_utf8_lossy(seq)))
}

/// annotateFilterReasons
fn filter_reasons(vdj: &VdjSequence, duplicate: bool, matches_ref: bool, support_min: i32, alignment: Option<&AlignmentAnnotation>) -> Vec<Filter> {
    let mut f = Vec::new();
    if matches_ref || alignment.is_some_and(|a| a.status == AlignmentStatus::NoRearrangement) {
        f.push(Filter::MatchesRef);
    }
    match alignment {
        None => {
            if vdj.v_anchor.is_none() {
                f.push(Filter::NoVAnchor);
            }
            if vdj.j_anchor.is_none() {
                f.push(Filter::NoJAnchor);
            }
            if let Some(VjAnchor { kind: AnchorKind::Blosum { similarity_score, .. }, .. }) = &vdj.v_anchor {
                if *similarity_score < 0 {
                    f.push(Filter::PoorVAnchor);
                }
            }
            if let Some(VjAnchor { kind: AnchorKind::Blosum { similarity_score, .. }, .. }) = &vdj.j_anchor {
                if *similarity_score < 0 {
                    f.push(Filter::PoorJAnchor);
                }
            }
        }
        Some(a) => {
            if let (Some(v), Some(j)) = (&a.v_gene, &a.j_gene) {
                if v.alignment.query_end - j.alignment.query_start >= ANNOTATION_V_J_OVERLAP_MAX {
                    f.push(Filter::VJOverlap);
                } else if matches!(a.status, AlignmentStatus::VD | AlignmentStatus::DJ | AlignmentStatus::VOnly | AlignmentStatus::JOnly) {
                    f.push(Filter::Partial);
                } else if a.status == AlignmentStatus::NoVdjAlignment {
                    f.push(Filter::NoVdjAlignment);
                }
            } else if matches!(a.status, AlignmentStatus::VD | AlignmentStatus::DJ | AlignmentStatus::VOnly | AlignmentStatus::JOnly) {
                f.push(Filter::Partial);
            } else if a.status == AlignmentStatus::NoVdjAlignment {
                f.push(Filter::NoVdjAlignment);
            }
        }
    }
    if duplicate {
        f.push(Filter::Duplicate);
    }
    let cdr3_len = vdj.cdr3_sequence().len() as i32;
    if cdr3_len < CDR3_FILTER_AA_MIN_LENGTH * 3 {
        f.push(Filter::MinLength);
    }
    if cdr3_len > CDR3_FILTER_AA_MAX_LENGTH * 3 {
        f.push(Filter::MaxLength);
    }
    if let (Some(v), Some(j)) = (vdj.v_anchor_boundary(), vdj.j_anchor_boundary()) {
        if v > j {
            f.push(Filter::Cdr3Deleted);
        }
    }
    if support_min <= 0 {
        f.push(Filter::NoHighQualSupport);
    }
    if f.is_empty() {
        f.push(Filter::Pass);
    }
    f
}

/// sortAndAnnotateVdjs: by descending CDR3 support, reads, length, CDR3 sequence; a sequence is a
/// duplicate when one before it has the same CDR3
pub fn sort_and_annotate(
    vdjs: &[VdjSequence],
    alignments: &[AlignmentAnnotation],
    primer_counts: &[(i32, i32)],
) -> Result<Vec<VdjAnnotation>> {
    let support: Vec<i32> = vdjs.iter().map(cdr3_support_min).collect();
    let mut order: Vec<usize> = (0..vdjs.len()).collect();
    order.sort_by(|&a, &b| {
        support[b]
            .cmp(&support[a])
            .then(vdjs[b].num_reads().cmp(&vdjs[a].num_reads()))
            .then(vdjs[b].length().cmp(&vdjs[a].length()))
            .then_with(|| vdjs[b].cdr3_sequence().cmp(vdjs[a].cdr3_sequence()))
    });
    let mut first_by_cdr3: HashMap<&[u8], usize> = HashMap::new();
    for &i in &order {
        let cdr3 = vdjs[i].cdr3_sequence();
        match first_by_cdr3.get(cdr3) {
            None => {
                first_by_cdr3.insert(cdr3, i);
            }
            Some(&e) => assert!(support[e] >= support[i], "java.lang.IllegalArgumentException: Failed requirement (duplicate support)"),
        }
    }
    let alignment_of: HashMap<usize, usize> = alignments.iter().enumerate().map(|(k, a)| (a.vdj, k)).collect();
    let mut out = Vec::with_capacity(order.len());
    for &i in &order {
        let vdj = &vdjs[i];
        let duplicate = first_by_cdr3[vdj.cdr3_sequence()] != i;
        let alignment = alignment_of.get(&i).copied();
        let a = alignment.map(|k| &alignments[k]);
        let v_aligned = anchor_reads(&vdj.v_anchor);
        let j_aligned = anchor_reads(&vdj.j_anchor);
        let v_non_split = count_non_split_reads(vdj, Vj::V);
        let j_non_split = count_non_split_reads(vdj, Vj::J);
        let support_min = support[i];
        let mref = matches_ref(vdj, v_aligned, j_aligned, v_non_split, j_non_split);
        let filters = filter_reasons(vdj, duplicate, mref, support_min, a);
        let v_similarity = match &vdj.v_anchor {
            Some(an) => Some(anchor_similarity(vdj, an)?),
            None => None,
        };
        let j_similarity = match &vdj.j_anchor {
            Some(an) => Some(anchor_similarity(vdj, an)?),
            None => None,
        };
        let shm = a.and_then(|a| a.v_gene.as_ref()).and_then(|g| g.comparison.as_ref()).map(|c| shm_status(c.pct_identity));
        let (vp, jp) = primer_counts.get(i).copied().unwrap_or((0, 0));
        out.push(VdjAnnotation {
            vdj: i,
            filters,
            v_aligned_reads: v_aligned,
            j_aligned_reads: j_aligned,
            v_non_split_reads: v_non_split,
            j_non_split_reads: j_non_split,
            cdr3_support_min: support_min,
            v_similarity_score: v_similarity,
            j_similarity_score: j_similarity,
            v_primer_match_count: vp,
            j_primer_match_count: jp,
            alignment,
            somatic_hypermutation_status: shm,
        });
    }
    Ok(out)
}

/// passesHardFilters: matches to the reference go (without -report_match_ref_seq), and short sequences
/// with only a V or a J
pub fn passes_hard_filters(a: &VdjAnnotation, alignment: Option<&AlignmentAnnotation>, report_match_ref_seq: bool) -> bool {
    if !report_match_ref_seq && a.filters.contains(&Filter::MatchesRef) {
        return false;
    }
    if a.filters.contains(&Filter::MinLength) {
        match alignment {
            None => {
                if a.filters.contains(&Filter::NoVAnchor) || a.filters.contains(&Filter::NoJAnchor) {
                    return false;
                }
            }
            Some(al) => {
                if al.status == AlignmentStatus::VOnly || al.status == AlignmentStatus::JOnly {
                    return false;
                }
            }
        }
    }
    true
}

/// the number of layout reads matched by alignment in a layout (for the layout file)
pub fn split_read_counts(layout: &crate::layout::ReadLayout) -> (usize, usize) {
    let mut five = 0;
    let mut ten = 0;
    for r in layout.reads() {
        let c = candidate(r);
        let m = c.left_soft_clip.max(c.right_soft_clip);
        if m >= 5 {
            five += 1;
        }
        if m >= 10 {
            ten += 1;
        }
    }
    (five, ten)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::cigar_from_str;

    #[test]
    fn adjusted_blocks() {
        assert_eq!(adjusted_alignment_blocks(&cigar_from_str("10S20M5I20M3D10M15S")), vec![(11, 55)]);
        assert_eq!(adjusted_alignment_blocks(&cigar_from_str("20M1000N30M")), vec![(1, 20), (21, 30)]);
        assert_eq!(adjusted_alignment_blocks(&cigar_from_str("5H20M")), vec![(1, 20)]);
    }

    #[test]
    fn known_matches_ref() {
        // CTXGPKXELRT... as amino acids: TGT ACT TAA GGT CCT AAA TAG GAA CTT AGA ACT
        let v = crate::vdj::tests::create_vdj("TGTACTTAAGGTCCTAAATAGGAACTTAGAACTTGGGGG", Some(3), Some(33));
        assert!(in_known_matches_ref_list(&v));
    }
}
