//! Tests of the visualiser's own logic on synthetic input. (The Java vis package has no unit tests.)

use std::cmp::Ordering;

use crate::assembly::sam::{cigar_elements_from_str, cigar_elements_to_str, CigarElement};
use crate::assembly::sequence_diff::{BaseQualType, SequenceDiffInfo, SequenceDiffType};
use crate::assembly::vis::assembly_visualiser::priority_queue_order_for_test;
use crate::assembly::vis::read_view_model::test_access::{correct, map_and_collapse};
use crate::assembly::vis::vis_cigar;
use crate::common::vis::svg::VisCigarOp;

#[test]
fn priority_queue_order_matches_java() {
    // the poll order of java.util.PriorityQueue for these (key1, key2, tag) items, from testdata/VisProbe.java
    let items: [(i32, i32, i32); 13] =
        [(5, 9, 0), (3, 7, 1), (5, 9, 2), (3, 7, 3), (1, 2, 4), (5, 9, 5), (3, 7, 6), (8, 8, 7), (1, 2, 8), (5, 1, 9), (3, 7, 10), (0, 0, 11), (5, 9, 12)];
    let order = priority_queue_order_for_test(&items, |a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let tags: Vec<String> = order.iter().map(|x| x.2.to_string()).collect();
    let expected = include_str!("../../common/vis/testdata/expected.tsv").lines().find_map(|l| l.strip_prefix("pq.order\t")).unwrap();
    assert_eq!(expected, tags.join(" "));
    assert!(priority_queue_order_for_test(&[] as &[i32], |_, _| Ordering::Equal).is_empty());
    assert_eq!(vec![7], priority_queue_order_for_test(&[7], |a: &i32, b: &i32| a.cmp(b)));
}

fn cigar(text: &str) -> Vec<CigarElement> {
    cigar_elements_from_str(text)
}

fn diff(read_index: i32, type_: SequenceDiffType, indel_length: i32) -> SequenceDiffInfo {
    let mut info = SequenceDiffInfo::new(read_index, 0, None, type_, BaseQualType::High);
    info.indel_length = indel_length;
    info
}

#[test]
fn cigar_operators_map_and_collapse() {
    for text in ["10S90M", "5M2D5M", "3M2I5M", "2H4M", "4M3D", "7M"] {
        assert_eq!(text, cigar_elements_to_str(&map_and_collapse(&cigar(text))));
    }
    assert!(map_and_collapse(&[]).is_empty());
}

#[test]
fn indel_differences_are_worked_into_the_cigar() {
    let text = |r: &(Option<Vec<CigarElement>>, i32)| r.0.as_ref().map(|c| cigar_elements_to_str(c));
    // no differences: the read's own cigar is used
    assert_eq!((None, 0), correct(true, &[], &cigar("10M")).unwrap());
    // a base difference changes nothing
    let r = correct(false, &[diff(3, SequenceDiffType::Base, 0)], &cigar("10M")).unwrap();
    assert_eq!((Some("10M".to_string()), 0), (text(&r), r.1));
    // an insert of two bases from read index 4, forwards and backwards
    let r = correct(true, &[diff(4, SequenceDiffType::Insert, 2)], &cigar("10M")).unwrap();
    assert_eq!((Some("4M2I4M".to_string()), 0), (text(&r), r.1));
    let r = correct(false, &[diff(4, SequenceDiffType::Insert, -2)], &cigar("10M")).unwrap();
    assert_eq!((Some("3M2I5M".to_string()), 2), (text(&r), r.1));
    // an insert that runs past the read is cut short
    let r = correct(true, &[diff(9, SequenceDiffType::Insert, 3)], &cigar("10M")).unwrap();
    assert_eq!((Some("9M1I".to_string()), 0), (text(&r), r.1));
    // a deletion of three bases after read index 4
    let r = correct(true, &[diff(4, SequenceDiffType::Delete, 3)], &cigar("10M")).unwrap();
    assert_eq!((Some("5M3D5M".to_string()), 0), (text(&r), r.1));
    let r = correct(false, &[diff(4, SequenceDiffType::Delete, 3)], &cigar("10M")).unwrap();
    assert_eq!((Some("5M3D5M".to_string()), -3), (text(&r), r.1));
    // a deletion at a read index the cigar does not have is Java's NullPointerException
    assert!(correct(true, &[diff(12, SequenceDiffType::Delete, 1)], &cigar("10M")).is_err());
    // a repeat expansion and contraction
    let r = correct(true, &[diff(4, SequenceDiffType::Repeat, 2)], &cigar("10M")).unwrap();
    assert_eq!((Some("5M2I3M".to_string()), 0), (text(&r), r.1));
    let r = correct(false, &[diff(4, SequenceDiffType::Repeat, 2)], &cigar("10M")).unwrap();
    assert_eq!((Some("3M2I5M".to_string()), 2), (text(&r), r.1));
    let r = correct(true, &[diff(4, SequenceDiffType::Repeat, -2)], &cigar("10M")).unwrap();
    assert_eq!((Some("5M2D5M".to_string()), 0), (text(&r), r.1));
    let r = correct(false, &[diff(4, SequenceDiffType::Repeat, -2)], &cigar("10M")).unwrap();
    assert_eq!((Some("4M2D6M".to_string()), -2), (text(&r), r.1));
    // a repeat contraction outside the read is ignored, unlike a deletion
    let r = correct(true, &[diff(12, SequenceDiffType::Repeat, -2)], &cigar("10M")).unwrap();
    assert_eq!((Some("10M".to_string()), 0), (text(&r), r.1));
    // soft clips keep their place; an insert over a clipped base replaces the clip operator
    let r = correct(true, &[diff(1, SequenceDiffType::Insert, 1)], &cigar("3S7M")).unwrap();
    assert_eq!((Some("1S1I1S7M".to_string()), 0), (text(&r), r.1));
    // other kinds of difference are not implemented in Java
    assert!(correct(true, &[diff(4, SequenceDiffType::Other, 0)], &cigar("10M")).is_err());
    assert!(correct(true, &[diff(4, SequenceDiffType::Match, 0)], &cigar("10M")).is_err());
}

#[test]
fn cigar_operators_for_view_models() {
    let ops: Vec<VisCigarOp> = vis_cigar(&cigar("1M1=1X1S1I1D1H1N1P")).into_iter().map(|x| x.1).collect();
    assert_eq!(
        vec![
            VisCigarOp::Match,
            VisCigarOp::Match,
            VisCigarOp::Match,
            VisCigarOp::S,
            VisCigarOp::I,
            VisCigarOp::D,
            VisCigarOp::Skip,
            VisCigarOp::Skip,
            VisCigarOp::Pad
        ],
        ops
    );
}

// corrected behaviour of the visualiser findings (REPORT.md 4.2), each beside the Java behaviour

#[test]
fn clips_at_both_ends_are_taken_off_the_alignment_cigar() {
    use crate::assembly::vis::assembly_visualiser::aligned_cigar_range;
    use crate::common::compat;
    // finding A33: the middle alignment of a chain is clipped at both ends
    let both = cigar("30S50M20S");
    assert_eq!((1, 2), aligned_cigar_range(&both, true));
    assert_eq!((1, 2), aligned_cigar_range(&both, false));
    assert_eq!((0, 1), aligned_cigar_range(&cigar("50M20S"), true));
    assert_eq!((1, 4), aligned_cigar_range(&cigar("30S20M2I28M"), false));
    assert_eq!((0, 1), aligned_cigar_range(&cigar("50M"), true));
    compat::with_java(true, || {
        assert_eq!((0, 2), aligned_cigar_range(&both, true));
        assert_eq!((1, 3), aligned_cigar_range(&both, false));
    });
}

#[test]
fn title_and_file_name_without_a_reference_or_a_tumor_sample() {
    use crate::assembly::vis::assembly_visualiser::{file_sample_id, variant_title};
    use crate::common::compat;
    // finding A32
    let tumor = vec!["SAMPLE_T".to_string()];
    let reference = vec!["SAMPLE_R".to_string()];
    let none: Vec<String> = Vec::new();
    assert_eq!("SAMPLE_T[normal: SAMPLE_R]", variant_title(&tumor, &reference).unwrap());
    assert_eq!("SAMPLE_T", variant_title(&tumor, &none).unwrap());
    assert_eq!("[normal: SAMPLE_R]", variant_title(&none, &reference).unwrap());
    assert_eq!("SAMPLE_T", file_sample_id(&tumor, &none).unwrap());
    assert_eq!("SAMPLE_R", file_sample_id(&none, &reference).unwrap());
    compat::with_java(true, || {
        assert_eq!("SAMPLE_T[normal: SAMPLE_R]", variant_title(&tumor, &reference).unwrap());
        assert!(variant_title(&tumor, &none).is_err());
        assert!(variant_title(&none, &reference).is_err());
        assert!(file_sample_id(&none, &reference).is_err());
    });
}

#[test]
fn drawing_does_not_change_the_lists_for_the_next_drawing() {
    use crate::assembly::tests::{create_alignment, create_read, test_env, TEST_READ_ID};
    use crate::assembly::vis::VisState;
    use crate::common::compat;
    // finding A35: what one drawing reverses is seen by the next in Java
    let _env = test_env();
    let alignment = create_alignment("1", 100, 149, 30, 79, "30S50M20S");
    let read = create_read(TEST_READ_ID, 100, &"ACGT".repeat(10), "10S30M");
    read.m().bam_record_mut().base_qualities = (1..=40).collect();
    let mut state = VisState::default();
    state.align_cigar(alignment).reverse();
    assert_eq!("30S50M20S", cigar_elements_to_str(state.align_cigar(alignment)));
    state.read(&read).quals.reverse();
    state.read(&read).cigar.reverse();
    assert_eq!(1, state.read(&read).quals[0]);
    assert_eq!("10S30M", cigar_elements_to_str(&state.read(&read).cigar));
    compat::with_java(true, || {
        let mut state = VisState::default();
        state.align_cigar(alignment).reverse();
        assert_eq!("20S50M30S", cigar_elements_to_str(state.align_cigar(alignment)));
        state.read(&read).quals.reverse();
        state.read(&read).cigar.reverse();
        assert_eq!(40, state.read(&read).quals[0]);
        assert_eq!("30M10S", cigar_elements_to_str(&state.read(&read).cigar));
    });
}
