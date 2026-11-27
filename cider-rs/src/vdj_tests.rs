//! Ports of the Java tests VDJSequenceBuilderTest, VdjBuilderUtilsTest, VdjAnnotatorTest, CiderUtilsTest,
//! ReadSliceTest and SomaticHypermutationTest (synthetic sequences of the Java tests).

use std::collections::HashMap;
use std::sync::Arc;

use anyhow::Result;

use crate::alignment::{compare_vj_region_to_imgt, Alignment, ImgtSequence, ShmGeneComparison};
use crate::annotate::{adjusted_alignment_blocks, count_non_split_reads};
use crate::blosum::{AnchorBlosumMatch, AnchorSearcher, SearchMode};
use crate::genes::{GeneDatastore, Strand, Vj, VjAnchorTemplate, VjGeneType};
use crate::layout::tests::layout_read;
use crate::layout::{LayoutRead, ReadLayout};
use crate::layout_builder::{LayoutAdaptor, ReadSlice, VjLayoutSource};
use crate::read::test_util::test_read;
use crate::read::{cigar_from_str, FLAG_FIRST, FLAG_PAIRED};
use crate::screener::{MatchMethod, VjReadCandidate};
use crate::vdj::tests::{create_layout, create_vdj, create_vdj_layout, MIN_BQ};
use crate::vdj::{calc_sequence_word_hashes, compare_ranges, find_sequence_overlap, merge_vdjs, AnchorKind, VdjSequence, VdjSequenceBuilder};

/// MockVJReadLayoutAdaptor: the anchor range of each layout as given
struct MockAdaptor {
    ranges: HashMap<usize, (i32, i32)>,
}

impl MockAdaptor {
    fn new() -> Self {
        MockAdaptor { ranges: HashMap::new() }
    }
    /// `first until end` as Kotlin writes it
    fn set(&mut self, layout: &Arc<ReadLayout>, first: i32, end: i32) {
        self.ranges.insert(Arc::as_ptr(layout) as usize, (first, end - 1));
    }
}

impl LayoutAdaptor for MockAdaptor {
    fn extrapolated_anchor_range(&self, _: Vj, layout: &ReadLayout) -> (i32, i32) {
        self.ranges[&(layout as *const ReadLayout as usize)]
    }
    fn anchor_match_method(&self, _: &ReadLayout) -> MatchMethod {
        MatchMethod::Align
    }
    fn template_anchor_sequence(&self, _: &ReadLayout) -> String {
        "TGACCC".into()
    }
}

/// MockAnchorBlosumSearcher: the match it is given
struct MockSearcher {
    m: Option<AnchorBlosumMatch>,
}

impl AnchorSearcher for MockSearcher {
    fn search(&self, _: &[u8], _: &[VjGeneType], _: SearchMode, _: i32, _: i32) -> Result<Option<AnchorBlosumMatch>> {
        Ok(self.m.clone())
    }
}

/// a datastore with TestUtils.ighJ1 (index 0) and ighV1_18 (index 1)
fn test_datastore() -> GeneDatastore {
    GeneDatastore::new(
        vec![
            VjAnchorTemplate::new(VjGeneType::Ighj, "IGHJ1", "01", None, "TGGGGCCAGGGCACCCTGGTCACCGTCTCC", None),
            VjAnchorTemplate::new(VjGeneType::Ighv, "IGHV1-18", "01", None, "AGATCTGACGACACGGCCGTGTATTACTGT", None),
        ],
        Vec::new(),
    )
    .unwrap()
}

fn builder<'a>(adaptor: &'a MockAdaptor, ds: &'a GeneDatastore, searcher: &'a MockSearcher, min_overlap: i32) -> VdjSequenceBuilder<'a> {
    VdjSequenceBuilder { adaptor, ds, searcher, min_base_quality: MIN_BQ, min_overlapped_bases: min_overlap }
}

fn layout(seq: &str) -> Arc<ReadLayout> {
    Arc::new(create_layout(seq, 0))
}

fn str_of(b: Option<&[u8]>) -> &str {
    std::str::from_utf8(b.unwrap()).unwrap()
}

fn check_vdj(v: &VdjSequence, v_type: VjGeneType, j_type: VjGeneType) {
    assert_eq!(v.v_anchor.as_ref().unwrap().gene_type, v_type);
    assert_eq!(v.j_anchor.as_ref().unwrap().gene_type, j_type);
    assert_eq!(str_of(v.v_anchor_sequence()), "GAATACC");
    assert_eq!(str_of(v.j_anchor_sequence()), "GAGTG");
    assert_eq!(std::str::from_utf8(v.cdr3_sequence_short()).unwrap(), "CACATCCTGA");
}

fn is_read_match(a: &crate::vdj::VjAnchor) -> bool {
    matches!(a.kind, AnchorKind::ReadMatch { .. })
}

const SEQ: &str = "TGCGAATACCCACATCCTGAGAGTGTCAGA";

// VDJSequenceBuilderTest.testBuildVdjFromVLayout
#[test]
fn build_vdj_from_v_layout() {
    let ds = test_datastore();
    let mut adaptor = MockAdaptor::new();
    let l = layout(SEQ);
    adaptor.set(&l, 3, 10);
    let searcher = MockSearcher {
        m: Some(AnchorBlosumMatch { anchor_start: 20, anchor_end: 25, template_anchor_seq: "CCCCC".into(), template_genes: vec![0], similarity_score: 10 }),
    };
    let b = builder(&adaptor, &ds, &searcher, crate::constants::MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES);
    let v = b.try_complete_layout_with_blosum(VjGeneType::Ighv, &l).unwrap().unwrap();
    check_vdj(&v, VjGeneType::Ighv, VjGeneType::Ighj);
    assert!(is_read_match(v.v_anchor.as_ref().unwrap()));
    assert!(!is_read_match(v.j_anchor.as_ref().unwrap()));
}

// testBuildVdjFromJLayout
#[test]
fn build_vdj_from_j_layout() {
    let ds = test_datastore();
    let mut adaptor = MockAdaptor::new();
    let l = layout(SEQ);
    adaptor.set(&l, 20, 25);
    let searcher = MockSearcher {
        m: Some(AnchorBlosumMatch { anchor_start: 3, anchor_end: 10, template_anchor_seq: "CCCCC".into(), template_genes: vec![1], similarity_score: 10 }),
    };
    let b = builder(&adaptor, &ds, &searcher, crate::constants::MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES);
    let v = b.try_complete_layout_with_blosum(VjGeneType::Ighj, &l).unwrap().unwrap();
    check_vdj(&v, VjGeneType::Ighv, VjGeneType::Ighj);
    assert!(!is_read_match(v.v_anchor.as_ref().unwrap()));
    assert!(is_read_match(v.j_anchor.as_ref().unwrap()));
}

fn overlap(v_seq: &str, j_seq: &str, v_range: (i32, i32), j_range: (i32, i32), min_overlap: i32) -> VdjSequence {
    let ds = test_datastore();
    let mut adaptor = MockAdaptor::new();
    let vl = layout(v_seq);
    let jl = layout(j_seq);
    adaptor.set(&vl, v_range.0, v_range.1);
    adaptor.set(&jl, j_range.0, j_range.1);
    let searcher = MockSearcher { m: None };
    let b = builder(&adaptor, &ds, &searcher, min_overlap);
    let v = b.try_overlap_vj(&vl, &jl, VjGeneType::Trav, VjGeneType::Traj, &HashMap::new()).unwrap().unwrap();
    check_vdj(&v, VjGeneType::Trav, VjGeneType::Traj);
    assert!(is_read_match(v.v_anchor.as_ref().unwrap()) && is_read_match(v.j_anchor.as_ref().unwrap()));
    v
}

// testBuildVdjFromOverlapLayouts
#[test]
fn build_vdj_from_overlap_layouts() {
    overlap(&SEQ[..21], &SEQ[8..], (3, 10), (12, 17), 10);
}

// testBuildVdjFromOverlapLayouts1: the J layout upstream of the V layout
#[test]
fn build_vdj_from_overlap_layouts_1() {
    let seq = "GATGCGAATACCCACATCCTGAGAGTGTCAGA";
    let v = overlap(&seq[3..], &seq[..28], (2, 9), (22, 27), 10);
    assert_eq!(v.layout.consensus_string(), seq);
}

// testBuildVdjFromOverlapLayouts3: the V layout contains the J layout
#[test]
fn build_vdj_from_overlap_layouts_3() {
    let seq = "GATGCGAATACCCACATCCTGAGAGTGTCAGA";
    let v = overlap(seq, &seq[7..28], (5, 12), (15, 20), 10);
    assert_eq!(v.layout.consensus_string(), seq);
}

// testBuildVdjFromOverlapLayouts4: the J layout contains the V layout
#[test]
fn build_vdj_from_overlap_layouts_4() {
    let seq = "GATGCGAATACCCACATCCTGAGAGTGTCAGA";
    let v = overlap(&seq[4..25], seq, (1, 8), (22, 27), 10);
    assert_eq!(v.layout.consensus_string(), seq);
}

// testOverlapLayoutsNegativeJAnchorPos
#[test]
fn overlap_layouts_negative_j_anchor_pos() {
    let seq = "TGCGAATACCCACATCCTGAGAGTGTCAGATAACT";
    overlap(&seq[..32], &seq[24..], (3, 10), (-4, 1), 5);
}

// testOverlapLayoutsNegativeVAnchorPos
#[test]
fn overlap_layouts_negative_v_anchor_pos() {
    let seq = "ACTTGCGAATACCCACATCCTGAGAGTGTCAGATAACT";
    overlap(&seq[..7], &seq[1..], (6, 13), (22, 27), 5);
}

// testBuildOneSidedVdjVOnly
#[test]
fn build_one_sided_vdj_v_only() {
    let ds = test_datastore();
    let mut adaptor = MockAdaptor::new();
    let l = layout(SEQ);
    adaptor.set(&l, 3, 10);
    let searcher = MockSearcher { m: None };
    let b = builder(&adaptor, &ds, &searcher, crate::constants::MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES);
    let v = b.try_create_one_sided_vdj(VjGeneType::Ighv, &l).unwrap();
    assert!(v.j_anchor.is_none());
    assert_eq!(v.v_anchor.as_ref().unwrap().gene_type, VjGeneType::Ighv);
    assert_eq!(v.sequence_formatted(), "GAATACC-CACATCCTGAGAGTGTCAGA");
    assert_eq!(str_of(v.v_anchor_sequence()), "GAATACC");
    assert_eq!(std::str::from_utf8(v.cdr3_sequence_short()).unwrap(), "CACATCCTGAGAGTGTCAGA");
}

// testBuildOneSidedVdjJOnly
#[test]
fn build_one_sided_vdj_j_only() {
    let ds = test_datastore();
    let mut adaptor = MockAdaptor::new();
    let l = layout("TGCGAATACCCACATCTGAGTGTCAGA");
    adaptor.set(&l, 15, 22);
    let searcher = MockSearcher { m: None };
    let b = builder(&adaptor, &ds, &searcher, crate::constants::MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES);
    let v = b.try_create_one_sided_vdj(VjGeneType::Ighj, &l).unwrap();
    assert!(v.v_anchor.is_none());
    assert_eq!(v.j_anchor.as_ref().unwrap().gene_type, VjGeneType::Ighj);
    assert_eq!(v.sequence_formatted(), "TGCGAATACCCACAT-CTGAGTG");
    assert_eq!(str_of(v.j_anchor_sequence()), "CTGAGTG");
    assert_eq!(std::str::from_utf8(v.cdr3_sequence_short()).unwrap(), "TGCGAATACCCACAT");
}

// testBuildVdjFromOverlapLayoutsWordHash
#[test]
fn build_vdj_from_overlap_layouts_word_hash() {
    // join_by_word_hash is private to the builder; build_vdj_sequences runs it after the shared-read join
    // (no shared reads here) and leaves no one-sided layout
    let ds = test_datastore();
    let mut adaptor = MockAdaptor::new();
    let vl = layout(&SEQ[..21]);
    let jl = layout(&SEQ[8..]);
    adaptor.set(&vl, 3, 10);
    adaptor.set(&jl, 12, 17);
    let searcher = MockSearcher { m: None };
    let b = builder(&adaptor, &ds, &searcher, 10);
    let layouts = vec![(VjGeneType::Trav, vec![vl]), (VjGeneType::Traj, vec![jl])];
    let vdjs = b.build_vdj_sequences(&layouts, &HashMap::new(), 1).unwrap();
    assert_eq!(vdjs.len(), 1);
    check_vdj(&vdjs[0], VjGeneType::Trav, VjGeneType::Traj);
}

// ---- VdjBuilderUtilsTest ----

#[test]
fn calc_sequence_hashes() {
    let h = calc_sequence_word_hashes(b"AAAAATAGAAAAAA", 6).unwrap();
    assert_eq!(h, vec![1, 4, 19, 76, 304, 1216, 768, 3072, 0]);
}

#[test]
fn find_vj_overlap_1() {
    let s1 = b"GACTCAGTCCAGGAAAGTATT";
    let s2 = b"GTCCAGGAAAGTATTATACTCA";
    let o = find_sequence_overlap(s1, s2, 10).unwrap();
    assert_eq!((o.seq1_offset, o.seq2_offset, o.overlap_bases, o.high_qual_match_bases), (0, 6, 15, 15));
    let o = find_sequence_overlap(s2, s1, 10).unwrap();
    assert_eq!((o.seq1_offset, o.seq2_offset, o.overlap_bases, o.high_qual_match_bases), (6, 0, 15, 15));
    assert!(find_sequence_overlap(s1, s2, 16).is_none());
    assert!(find_sequence_overlap(s2, s1, 16).is_none());
}

#[test]
fn find_vj_overlap_long_short() {
    let s1 = b"GACTCAGTCCAGGAAAGTATTATACTCA";
    let s2 = b"GTCCAGGAAAGTATT";
    let o = find_sequence_overlap(s1, s2, 10).unwrap();
    assert_eq!((o.seq1_offset, o.seq2_offset, o.overlap_bases), (0, 6, 15));
    let o = find_sequence_overlap(s2, s1, 10).unwrap();
    assert_eq!((o.seq1_offset, o.seq2_offset, o.overlap_bases), (6, 0, 15));
}

/// ranges as (first, last)
fn ranges(a: &VdjSequence, b: &VdjSequence) -> ((i32, i32), (i32, i32)) {
    let ((s1, e1), (s2, e2)) = compare_ranges(a, b);
    ((s1, e1 - 1), (s2, e2 - 1))
}

#[test]
fn get_compare_ranges() {
    let a = create_vdj("GATGCGAATACC", Some(4), Some(7));
    let b = create_vdj("TGCGAATACCAG", Some(2), Some(5));
    assert_eq!(ranges(&a, &b), ((2, 11), (0, 9)));
}

#[test]
fn get_compare_ranges_v_only() {
    let a = create_vdj("GATGCGAATACC", Some(4), None);
    let b = create_vdj("TGCGAATACCAG", Some(2), None);
    assert_eq!(ranges(&a, &b), ((2, 11), (0, 9)));
}

#[test]
fn get_compare_ranges_v_only_long() {
    let a = create_vdj(&"T".repeat(100), Some(4), None);
    let b = create_vdj(&"T".repeat(100), Some(2), None);
    let p = crate::constants::PARTIAL_VDJ_UNANCHORED_LENGTH_BASES;
    assert_eq!(ranges(&a, &b), ((2, 4 + p - 1), (0, 2 + p - 1)));
}

#[test]
fn get_compare_ranges_j_only() {
    // as the Java test: V boundaries of 7 and 9
    let a = create_vdj("TGCGAATACCAG", Some(7), None);
    let b = create_vdj("GATGCGAATACC", Some(9), None);
    assert_eq!(ranges(&a, &b), ((0, 9), (2, 11)));
}

#[test]
fn get_compare_ranges_j_only_long() {
    let a = create_vdj(&"T".repeat(100), None, Some(90));
    let b = create_vdj(&"T".repeat(100), None, Some(95));
    let p = crate::constants::PARTIAL_VDJ_UNANCHORED_LENGTH_BASES;
    assert_eq!(ranges(&a, &b), ((90 - p, 94), (95 - p, 99)));
}

#[test]
fn merge_vdjs_simple() {
    let seq = "ATGCTGGTGT";
    let mut l1 = ReadLayout::new();
    l1.add_read(layout_read("read1", seq, "FF:FFFF:FF", 4, true), MIN_BQ);
    let mut l2 = ReadLayout::new();
    l2.add_read(layout_read("read2", seq, "FFFF::FFFF", 4, true), MIN_BQ);
    let v1 = create_vdj_layout(Arc::new(l1), Some(3), Some(7), 0, 10);
    let v2 = create_vdj_layout(Arc::new(l2), Some(3), Some(7), 0, 10);
    let m = merge_vdjs(&v1, &v2, MIN_BQ).unwrap();
    assert_eq!(std::str::from_utf8(m.sequence()).unwrap(), seq);
}

fn merge_both_ways(l1: (&str, i32), v1: (Option<i32>, Option<i32>, i32, i32), l2: (&str, i32), v2: (Option<i32>, Option<i32>, i32, i32)) -> Vec<VdjSequence> {
    let a = create_vdj_layout(Arc::new(create_layout(l1.0, l1.1)), v1.0, v1.1, v1.2, v1.3);
    let b = create_vdj_layout(Arc::new(create_layout(l2.0, l2.1)), v2.0, v2.1, v2.2, v2.3);
    vec![merge_vdjs(&a, &b, MIN_BQ).unwrap(), merge_vdjs(&b, &a, MIN_BQ).unwrap()]
}

#[test]
fn merge_vdjs_long_short() {
    let a = create_vdj_layout(Arc::new(create_layout("CCTCAGGTGAT", 5)), Some(3), Some(7), 2, 11);
    assert_eq!(a.sequence_formatted(), "TCA-GGTG-AT");
    let b = create_vdj_layout(Arc::new(create_layout("AGGTGATAA", 3)), Some(1), Some(5), 0, 8);
    assert_eq!(b.sequence_formatted(), "A-GGTG-ATA");
    for m in merge_both_ways(("CCTCAGGTGAT", 5), (Some(3), Some(7), 2, 11), ("AGGTGATAA", 3), (Some(1), Some(5), 0, 8)) {
        assert!(m.v_anchor.is_some() && m.j_anchor.is_some());
        assert_eq!(m.layout.consensus_string(), "CCTCAGGTGATAA");
        assert_eq!(m.sequence_formatted(), "TCA-GGTG-ATA");
        assert_eq!((m.v_anchor_boundary(), m.j_anchor_boundary()), (Some(3), Some(7)));
    }
}

#[test]
fn merge_vdjs_v_only() {
    for m in merge_both_ways(("CCTCAGGTG", 5), (Some(3), None, 2, 9), ("AGGTGATAA", 3), (Some(1), None, 0, 7)) {
        assert!(m.v_anchor.is_some() && m.j_anchor.is_none());
        assert_eq!(m.layout.consensus_string(), "CCTCAGGTGATAA");
        assert_eq!(m.sequence_formatted(), "TCA-GGTGAT");
        assert_eq!(m.v_anchor_boundary(), Some(3));
    }
}

#[test]
fn merge_vdjs_j_only() {
    for m in merge_both_ways(("CCTCAGGTGAT", 5), (None, Some(7), 2, 11), ("AGGTGATAA", 3), (None, Some(5), 0, 8)) {
        assert!(m.v_anchor.is_none() && m.j_anchor.is_some());
        assert_eq!(m.layout.consensus_string(), "CCTCAGGTGATAA");
        assert_eq!(m.sequence_formatted(), "TCAGGTG-ATA");
        assert_eq!(m.j_anchor_boundary(), Some(7));
    }
}

// ---- VdjAnnotatorTest ----

/// a VDJ of one read whose source is the candidate read `read_seq` with `cigar`, sliced as given
fn vdj_with_read(vdj_seq: &str, v: Option<i32>, j: Option<i32>, reads: &[(&str, i32, &str, bool, i32, i32)]) -> VdjSequence {
    // reads: (layout sequence, aligned position, candidate read sequence + cigar, reverse complement, slice start, slice end)
    let mut l = ReadLayout::new();
    for (i, &(seq, aligned, cigar, rev, s0, s1)) in reads.iter().enumerate() {
        let (read_seq, cigar) = cigar.split_once(' ').unwrap();
        let mut r = test_read();
        r.flags = FLAG_PAIRED | FLAG_FIRST;
        r.bases = read_seq.as_bytes().to_vec();
        r.quals = vec![37; read_seq.len()];
        r.cigar = cigar_from_str(cigar);
        r.alignment_start = 1;
        r.tid = 0;
        let r = Arc::new(r);
        let c = Arc::new(VjReadCandidate {
            read: r.clone(),
            templates: Vec::new(),
            gene_type: VjGeneType::Trav,
            template_anchor_sequence: "CACGTG".into(),
            match_method: MatchMethod::Align,
            use_reverse_complement: false,
            anchor_offset_start: 0,
            anchor_offset_end: 30,
            left_soft_clip: 0,
            right_soft_clip: 0,
            similarity_score: i32::MIN,
        });
        let mut lr: LayoutRead = layout_read(&format!("read{i}"), seq, &"F".repeat(seq.len()), aligned, true);
        lr.source = Some(Arc::new(VjLayoutSource { slice: ReadSlice::new(r, rev, s0, s1), candidate: c }));
        l.add_read(lr, MIN_BQ);
    }
    let len = l.length();
    let _ = vdj_seq;
    create_vdj_layout(Arc::new(l), v, j, 0, len)
}

#[test]
fn count_non_split_read() {
    let read_seq = "TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC";
    let vdj_seq = &read_seq[20..90];
    let count = |cigar: &str, rev: bool| {
        let spec = format!("{read_seq} {cigar}");
        let v = vdj_with_read(vdj_seq, Some(10), None, &[(vdj_seq, 0, &spec, rev, 20, 90)]);
        count_non_split_reads(&v, Vj::V)
    };
    assert_eq!(count("100M", false), 1);
    assert_eq!(count("20S80M", false), 1);
    assert_eq!(count("21S79M", false), 0);
    assert_eq!(count("60M40S", false), 1);
    assert_eq!(count("59M41S", false), 0);
    // the slice reverse complemented
    assert_eq!(count("100M", true), 1);
    assert_eq!(count("80M20S", true), 1);
    assert_eq!(count("79M21S", true), 0);
    assert_eq!(count("40S60M", true), 1);
    assert_eq!(count("41S59M", true), 0);
}

#[test]
fn count_non_split_read_aligned_pos() {
    let vdj_seq = "GAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGA";
    let read2 = &vdj_seq[20..];
    let count = |cigar2: &str| {
        // as in the Java test, the first read's record has no CIGAR (no aligned blocks)
        let spec1 = format!("{vdj_seq} ");
        let spec2 = format!("{read2} {cigar2}");
        // the layout's aligned position is the first read's, 0; the second read starts 20 bases on
        let v = vdj_with_read(vdj_seq, Some(10), Some(50), &[(vdj_seq, 0, &spec1, false, 0, vdj_seq.len() as i32), (read2, -20, &spec2, false, 0, read2.len() as i32)]);
        assert_eq!(std::str::from_utf8(v.sequence()).unwrap(), vdj_seq);
        count_non_split_reads(&v, Vj::J)
    };
    assert_eq!(count("50M"), 1);
    assert_eq!(count("49M1S"), 0);
    assert_eq!(count("1S49M"), 0);
}

// ---- CiderUtilsTest ----

#[test]
fn get_adjusted_alignment_blocks() {
    let b = |c: &str| adjusted_alignment_blocks(&cigar_from_str(c));
    assert_eq!(b("10S100M10S"), vec![(11, 100)]);
    assert_eq!(b("100M"), vec![(1, 100)]);
    assert_eq!(b("100M10S"), vec![(1, 100)]);
    assert_eq!(b("10S50M2I50M10S"), vec![(11, 102)]);
    assert_eq!(b("10S50M5D50M10S"), vec![(11, 100)]);
    assert_eq!(b("10S50M1000N50M10S"), vec![(11, 50), (61, 50)]);
    assert_eq!(b("10H50M1000N50M10H"), vec![(1, 50), (51, 50)]);
}

// ---- ReadSliceTest ----

#[test]
fn read_pos_to_slice_pos() {
    let mut r = test_read();
    r.bases = b"AGCTAGCTAG".to_vec();
    r.quals = vec![30; 10];
    let s = ReadSlice::new(Arc::new(r), false, 2, 7);
    assert_eq!(s.read_length(), 5);
    assert_eq!(s.read_string(), b"CTAGC");
    assert_eq!(s.read_position_to_slice_position(3), 1);
    assert_eq!(s.slice_position_to_read_position(1), 3);
}

#[test]
fn read_pos_to_slice_pos_reverse_comp() {
    let mut r = test_read();
    r.bases = b"CTAGCTAGCT".to_vec();
    r.quals = vec![30; 10];
    let s = ReadSlice::new(Arc::new(r), true, 2, 7);
    assert_eq!(s.read_length(), 5);
    assert_eq!(s.read_string(), b"CTAGC");
    assert_eq!(s.read_position_to_slice_position(3), 4);
    assert_eq!(s.slice_position_to_read_position(4), 3);
    assert_eq!(s.base_at(4), b'C');
    assert_eq!(s.read_position_to_slice_position(5), 2);
    assert_eq!(s.slice_position_to_read_position(2), 5);
    assert_eq!(s.base_at(2), b'A');
}

// ---- SomaticHypermutationTest ----

fn alignment(query: &str, q: (i32, i32), r: (i32, i32), strand: Strand, edit: i32, cigar: &str, ref_len: usize) -> Alignment {
    Alignment {
        query_seq: query.into(),
        query_start: q.0,
        query_end: q.1,
        ref_contig: "test".into(),
        ref_start: r.0,
        ref_end: r.1,
        strand,
        alignment_score: 0,
        edit_distance: edit,
        cigar: cigar_from_str(cigar),
        ref_contig_length: ref_len as i32,
    }
}

fn imgt(seq: &str, before: i32, after: i32) -> ImgtSequence {
    ImgtSequence { gene_name: "test".into(), allele: "1".into(), sequence_with_ref: seq.as_bytes().to_vec(), ref_before: before, ref_after: after }
}

fn rc(s: &str) -> String {
    String::from_utf8(hmf_common::codon::reverse_complement(s.as_bytes())).unwrap()
}

/// `a until b` as (a, b - 1)
fn until(a: i32, b: i32) -> (i32, i32) {
    (a, b - 1)
}

#[test]
fn shm_status_from_identity() {
    use crate::alignment::shm_status;
    for (p, s) in [(98.0, "UNMUTATED"), (98.01, "UNMUTATED"), (99.0, "UNMUTATED"), (99.99, "UNMUTATED"), (100.0, "UNMUTATED")] {
        assert_eq!(shm_status(p), s);
    }
    for p in [97.0, 97.01, 97.99] {
        assert_eq!(shm_status(p), "MUTATED_BORDERLINE");
    }
    for p in [96.99, 96.0, 95.0, 90.0] {
        assert_eq!(shm_status(p), "MUTATED");
    }
}

fn shm(seq_length: i32, imgt_length: i32, pct: f64, indel: i32, clip: i32, unaligned: i32) -> ShmGeneComparison {
    ShmGeneComparison { seq_length, imgt_length, pct_identity: pct, indel_bases: indel, clip_bases: clip, imgt_ref_unaligned: unaligned }
}

#[test]
fn compare_vj_region_to_imgt_v_forward() {
    let layout = "GGGTTGGGGGCAAAAACAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT";
    let imgt_seq = "AAAAAGGGGGAAAAAAAAAACCCCCAAAAAAA";
    let q = until(2, 47);
    let a = alignment(&layout[2..47], until(1, 21), until(3, 23), Strand::Forward, 27, "1S20M24S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::V, 20, &imgt(imgt_seq, 5, 7), q, &a);
    assert_eq!(c, Some(shm(15, 15, 100.0 * 13.0 / 15.0, 0, 0, 0)));
}

#[test]
fn compare_vj_region_to_imgt_v_reverse() {
    let layout = "GGGTTGGGGGCAAAAACAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT";
    let imgt_seq = rc("AAAAAGGGGGAAAAAAAAAACCCCCAAAAAAA");
    let n = imgt_seq.len() as i32;
    let a = alignment(&layout[2..47], until(1, 21), until(n - 23, n - 3), Strand::Reverse, 27, "24S20M1S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::V, 20, &imgt(&imgt_seq, 7, 5), until(2, 47), &a);
    assert_eq!(c, Some(shm(15, 15, 100.0 * 13.0 / 15.0, 0, 0, 0)));
}

#[test]
fn compare_vj_region_to_imgt_j_forward() {
    let layout = "GGGGGGGGGGAAAAAAAAAACCCCCCCTCTCAAAAAAAAGTTTTTTTTTT";
    let imgt_seq = "AAAAACCCCCAAAAAAAAAATTTTTAAAAA";
    let a = alignment(&layout[2..48], until(25, 45), until(7, 27), Strand::Forward, 27, "25S20M1S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::J, 30, &imgt(imgt_seq, 5, 5), until(2, 48), &a);
    assert_eq!(c, Some(shm(15, 15, 100.0 * 13.0 / 15.0, 0, 0, 0)));
}

#[test]
fn compare_vj_region_to_imgt_j_reverse() {
    let layout = "GGGGGGGGGGAAAAAAAAAACCCCCCCTCTCAAAAAAAAGTTTTTTTTTT";
    let imgt_seq = rc("AAAAACCCCCAAAAAAAAAATTTTTAAAAA");
    let n = imgt_seq.len() as i32;
    let a = alignment(&layout[2..48], until(25, 45), until(n - 27, n - 7), Strand::Reverse, 27, "1S20M25S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::J, 30, &imgt(&imgt_seq, 5, 5), until(2, 48), &a);
    assert_eq!(c, Some(shm(15, 15, 100.0 * 13.0 / 15.0, 0, 0, 0)));
}

#[test]
fn compare_vj_region_to_imgt_clip_v() {
    let layout = "GGGGGGGGGGAAAAAAAAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT";
    let imgt_seq = "AAAAAGGGGGAAAAAAAAAACCCCCAAAAA";
    let a = alignment(&layout[2..47], until(5, 25), until(7, 27), Strand::Forward, 25, "5S20M20S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::V, 20, &imgt(imgt_seq, 5, 5), until(2, 47), &a);
    assert_eq!(c, Some(shm(13, 13, 100.0, 0, 2, 2)));
}

#[test]
fn compare_vj_region_to_imgt_clip_j() {
    let layout = "GGGGGGGGGGAAAAAAAAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT";
    let imgt_seq = "AAAAACCCCCAAAAAAAAAATTTTTAAAAA";
    let a = alignment(&layout[2..48], until(25, 40), until(7, 22), Strand::Forward, 31, "25S15M6S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::J, 30, &imgt(imgt_seq, 5, 5), until(2, 48), &a);
    assert_eq!(c, Some(shm(12, 12, 100.0, 0, 3, 3)));
}

#[test]
fn compare_vj_region_to_imgt_v_assembly_short() {
    let layout = "TTTTTTTGGGAAAAAAAAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT";
    let imgt_seq = "AAAAAGGGGGGGGGAAAAAAAAAACCCCCAAAAA";
    let a = alignment(&layout[5..47], until(7, 25), until(11, 29), Strand::Forward, 25, "2S18M22S", imgt_seq.len());
    let c = compare_vj_region_to_imgt(layout.as_bytes(), Vj::V, 20, &imgt(imgt_seq, 5, 5), until(5, 47), &a);
    assert_eq!(c, Some(shm(13, 13, 100.0, 0, 2, 6)));
}
