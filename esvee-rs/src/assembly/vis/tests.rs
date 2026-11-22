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
