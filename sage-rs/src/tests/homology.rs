//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! common/HomologyTest.java

use super::util::*;
use crate::candidate::CandidateFinder;
use crate::read_context::{find_homology, find_left_homology_shift, Microhomology};
use crate::refseq::RefSequence;
use crate::types::*;

fn create_simple_variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

fn ref_sequence(start: i32, bases: &str) -> RefSequence {
    RefSequence::from_bases(start, bases.as_bytes().to_vec())
}

/// TestUtils.buildCigarString(alignedLength, leftSoftClip, rightSoftClip)
fn build_cigar_string(aligned_length: usize, left_soft_clip: usize, right_soft_clip: usize) -> String {
    let mut s = String::new();
    if left_soft_clip > 0 {
        s.push_str(&format!("{left_soft_clip}S"));
    }
    s.push_str(&format!("{aligned_length}M"));
    if right_soft_clip > 0 {
        s.push_str(&format!("{right_soft_clip}S"));
    }
    s
}

/// HomologyTest.findHomology
fn find_homology_test(variant: &SimpleVariant, read_bases: &[u8], var_read_index: i32) -> Option<Microhomology> {
    find_homology(variant, read_bases, var_read_index, variant.is_insert())
}

// HomologyTest.testHomologyLeftAlignment
#[test]
fn test_homology_left_alignment() {
    //            50        60        70        80
    //            0123456789012345678901234567890
    // ref bases: CTGTCTGTGACTCGGAAAAAAAACTCCCTGA

    // test 1: insert scenario:
    let var = create_simple_variant(70, "A", "AAAAA");

    let read_bases = format!("{}AAAAA{}", &REF_BASES_200[60..70], &REF_BASES_200[71..80]);
    let left_hom_shift = find_left_homology_shift(&var, &ref_sequence_200(), read_bases.as_bytes(), 10);
    assert_eq!(6, left_hom_shift);

    // test 2: dinucleotide repeats
    let ref_bases = "CTGTCTGTGACTCGGATATATATATATATATCCCCTTGCGCTTCCCAGGT";
    let ref_seq = ref_sequence(0, ref_bases);
    //            10        20        30
    //            012345678901234567890123
    // ref bases: CTCGGATATATATATATATATCCC

    let var = create_simple_variant(30, "T", "TATAT");

    let read_bases = format!("{}TATAT{}", &ref_bases[10..30], &ref_bases[31..40]);
    let left_hom_shift = find_left_homology_shift(&var, &ref_seq, read_bases.as_bytes(), 20);
    assert_eq!(16, left_hom_shift);

    // test 3: duplicate scenario
    let ref_bases = "CTGTCTGTGACAAACCCGGGTCGGATCCCGGTAGGTAT";
    //               0         10        20        30
    //               0123456789012345678901234567890123456789
    let ref_seq = ref_sequence(0, ref_bases);

    let inserted_bases = "AAACCCGGG";
    let var = create_simple_variant(19, "G", &format!("G{inserted_bases}"));

    let read_bases = format!("{}G{}{}", &ref_bases[10..19], inserted_bases, &ref_bases[20..35]);

    //          0         10        20        30        40
    //          012345678901234567890123456789012345678901234567890
    // read     CTGTCTGTGACAAACCCGGGAAACCCGGGTCGGATCCCGGTAGGTAT

    let left_hom_shift = find_left_homology_shift(&var, &ref_seq, read_bases.as_bytes(), 9);
    assert_eq!(9, left_hom_shift);

    // test 4: a duplication later on
    let ref_bases = "CTGTCTGTGACAAACCCGGGAAACCCGGGTCGGATCCCGGTAGGTAT";
    //               0         10        20        30        40
    //               01234567890123456789012345678901234567890123456789
    let ref_seq = ref_sequence(0, ref_bases);
    let var = create_simple_variant(28, "G", &format!("G{inserted_bases}"));

    let read_bases = format!("{}G{}{}", &ref_bases[10..28], inserted_bases, &ref_bases[20..35]);

    let left_hom_shift = find_left_homology_shift(&var, &ref_seq, read_bases.as_bytes(), 18);
    assert_eq!(18, left_hom_shift);
}

// HomologyTest.testHomologyLeftAlignmentImperfectHomology
#[test]
fn test_homology_left_alignment_imperfect_homology() {
    //                         10        20        30        40        50        60        70
    //                1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    let ref_bases = "TTTGGTGAGGAATGGGATACTGTCTGAAAGTATCTCTACAAAGGGATAAAGAGTAACTTTACATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAATAGGACAAATCA";
    let read_bases = "CATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAATAGGATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAAT";
    let alt = "ATCTAAGAAACCGACAGACATTAACTTAATCAAATGATCAAGGTGAAAATTACCAAATAGGA";

    let var = create_simple_variant(124, &alt[0..1], alt);

    let ref_seq = ref_sequence(1, ref_bases);

    let left_hom_shift = find_left_homology_shift(&var, &ref_seq, read_bases.as_bytes(), 62);
    assert_eq!(62, left_hom_shift);
}

// HomologyTest.testHomologyLeftAlignmentCandidate
#[test]
fn test_homology_left_alignment_candidate() {
    // test left-aligning a candidate variant found in a right soft-clip
    // the soft-clip must be long enough to have the duplicated section and 12-bases of matching ref

    let ref_bases = "CTGTCTGTGACACGTTGCAGTCGGATCCCGGTAGGTATTGTGACAAACCGTAGCTTGACCAG";
    //               0         10        20        30        40        50        60
    //               01234567890123456789012345678901234567890123456789012345678901
    let ref_seq = ref_sequence(0, ref_bases);

    // up to and including the ref base are duplicated
    let var_position = 30;
    let dup_length = 10;
    let duplicated_bases = &ref_bases[var_position - dup_length + 1..var_position + 1];

    let soft_clip_bases = format!("{}{}", duplicated_bases, &ref_bases[31..46]); // then continues in the ref sequence

    let read_start = 5;
    let ref_aligned_bases = &ref_bases[read_start..var_position + 1];
    let read_bases = format!("{ref_aligned_bases}{soft_clip_bases}");

    // RefContextCache + RefContextConsumer are one CandidateFinder in the port
    let cfg = test_config();
    let ref_data = empty_ref_data();
    let region = BaseRegion::new(0, 200);
    let mut finder = CandidateFinder::new(&cfg, CHR_1, region, &ref_seq, &[], &[], &ref_data);

    let cigar = build_cigar_string(ref_aligned_bases.len(), 0, soft_clip_bases.len());

    let read = build_sam_record(read_start as i32, &cigar, &read_bases);
    let read_clone = build_sam_record(read_start as i32, &cigar, &read_bases);
    finder.process_read(&read);
    finder.process_read(&read_clone);

    let alt_candidates = finder.finish();
    assert_eq!(1, alt_candidates.len());

    let alt_candidate = &alt_candidates[0];
    let read_context = &alt_candidate.ctx;
    assert_eq!(20, read_context.variant.pos);

    // again but this time with an indel in the read bases to confuse where the implied ref position is
    let del_position = 20;
    let read_bases =
        format!("{}{}{}", &ref_aligned_bases[0..del_position], &ref_aligned_bases[del_position + 2..], soft_clip_bases);

    let cigar = "20M2D4M25S";
    let read = build_sam_record(read_start as i32, cigar, &read_bases);
    let read_clone = build_sam_record(read_start as i32, cigar, &read_bases);

    let mut finder = CandidateFinder::new(&cfg, CHR_1, region, &ref_seq, &[], &[], &ref_data);
    finder.process_read(&read);
    finder.process_read(&read_clone);

    let alt_candidates = finder.finish();
    assert_eq!(2, alt_candidates.len());

    let alt_candidate = alt_candidates.iter().find(|x| x.ctx.variant.pos == 18);
    assert!(alt_candidate.is_some());
    let read_context = &alt_candidate.unwrap().ctx;
    assert_eq!(17, read_context.core_position_start);
    assert_eq!(39, read_context.core_position_end);
}

// HomologyTest.testHomology
#[test]
fn test_homology() {
    let var = create_simple_variant(26, "A", "AT");

    let read_bases = format!("{}AT{}", &REF_BASES_200[10..26], &REF_BASES_200[27..40]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "16M1I13M";
    let read = build_sam_record_quals(10, read_cigar, &read_bases, base_quals);

    let homology = find_homology_test(&var, &read.bases, 16);

    assert!(homology.is_some());
    let homology = homology.unwrap();
    assert_eq!(b"T".to_vec(), homology.bases);
    assert_eq!(4, homology.length);

    let var = create_simple_variant(26, "A", "ATT");

    let read_bases = format!("{}ATT{}", &REF_BASES_200[10..26], &REF_BASES_200[27..40]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "16M2I13M";
    let read = build_sam_record_quals(10, read_cigar, &read_bases, base_quals);

    let homology = find_homology_test(&var, &read.bases, 16);

    assert!(homology.is_some());
    let homology = homology.unwrap();
    assert_eq!(b"TT".to_vec(), homology.bases);
    assert_eq!(4, homology.length);

    // delete
    let var = create_simple_variant(64, "GAAA", "G");

    let read_bases = format!("{}{}", &REF_BASES_200[50..65], &REF_BASES_200[69..80]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "15M3D11M";
    let read = build_sam_record_quals(50, read_cigar, &read_bases, base_quals);

    let homology = find_homology_test(&var, &read.bases, 15);

    assert!(homology.is_some());
    let homology = homology.unwrap();
    assert_eq!(b"AAA".to_vec(), homology.bases);
    assert_eq!(3, homology.length);

    // checks read bases not ref bases to determine homology
    let var = create_simple_variant(26, "A", "AAAA");
    let read_bases = format!("{}AAAAAAAAAAAAAAAAAAAAATG", &REF_BASES_200[10..26]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "16M3I20M";
    let read = build_sam_record_quals(10, read_cigar, &read_bases, base_quals);

    let homology = find_homology_test(&var, &read.bases, 16);

    assert!(homology.is_some());
    let homology = homology.unwrap();
    assert_eq!(b"AAA".to_vec(), homology.bases);
    assert_eq!(17, homology.length);

    // multiple copies and then finishes with a partial copy
    // eg G(AACTC)AACTCAACTCAACCCTTT -> GAACTCAACTCAA(CTCAA)CCCTTT

    let var = create_simple_variant(26, "G", "GAACTC");
    let read_bases = format!("{}GAACTCAACTCAACTCAACCCTTT", &REF_BASES_200[10..26]);
    let base_quals = default_quals(read_bases.len());
    let read_cigar = "16M5I18M";
    let read = build_sam_record_quals(10, read_cigar, &read_bases, base_quals);

    let homology = find_homology_test(&var, &read.bases, 16);

    assert!(homology.is_some());
    let homology = homology.unwrap();
    assert_eq!(b"AACTC".to_vec(), homology.bases);
    assert_eq!(13, homology.length);
}
