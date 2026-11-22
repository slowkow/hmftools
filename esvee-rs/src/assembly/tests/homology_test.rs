//! HomologyTest
use super::*;
use crate::assembly::alignment::breakend_builder::BreakendBuilder;
use crate::assembly::alignment::homology_data::HomologyData;
use crate::assembly::alignment::md_tag::{MdTag, MdTagType};
use crate::assembly::types::StructuralVariantType::{Del, Dup, Ins, Inv};

fn build_homology(left_md_tag: &str, right_md_tag: &str, sequence: &str) -> HomologyData {
    let seq_index_end = sequence.len() as i32;
    let pos_start = 1;
    let pos_end = pos_start + sequence.len() as i32 - 1;
    let cigar = make_cigar_string(sequence, 0, 0);
    let left_alignment = create_alignment_full(CHR_1, pos_start, pos_end, false, DEFAULT_MAP_QUAL, sequence.len() as i32, 0, seq_index_end, &cigar, "", left_md_tag);
    let right_alignment = create_alignment_full(CHR_1, pos_start, pos_end, false, DEFAULT_MAP_QUAL, sequence.len() as i32, 0, seq_index_end, &cigar, "", right_md_tag);
    HomologyData::determine_homology(sequence, &left_alignment, &right_alignment).expect("homology")
}

fn check_homology(homology: &HomologyData, bases: &str, exact_start: i32, exact_end: i32, inexact_start: i32, inexact_end: i32) {
    assert_eq!(bases, homology.homology);
    assert_eq!((exact_start, exact_end), (homology.exact_start, homology.exact_end));
    assert_eq!((inexact_start, inexact_end), (homology.inexact_start, homology.inexact_end));
}

#[test]
fn test_indel_homology_exact() {
    let _env = test_env();
    let bases_start = "AACCGG";
    // exact match
    let homology = HomologyData::determine_indel_homology(bases_start, bases_start, bases_start.len() as i32).expect("homology");
    assert_eq!(bases_start, homology.homology);
    assert_eq!((-3, 3), (homology.exact_start, homology.exact_end));
    let homology = HomologyData::determine_indel_homology(bases_start, "AACCTT", bases_start.len() as i32).expect("homology");
    assert_eq!("AACC", homology.homology);
    assert_eq!((-2, 2), (homology.exact_start, homology.exact_end));
}

#[test]
fn test_indel_homology_inexact() {
    let _env = test_env();
    let bases_start = "AACCGTTTAAAA";
    let bases_end = "AACCTTTTACCC";
    let homology = HomologyData::determine_indel_homology(bases_start, bases_end, bases_start.len() as i32).expect("homology");
    check_homology(&homology, "AACCGTTTA", -2, 2, -5, 4);
}

#[test]
fn test_homology() {
    let _env = test_env();
    // test 1: exact match
    let assembly_overlap = "TTCTTCTTCTC";
    check_homology(&build_homology("11", "11", assembly_overlap), assembly_overlap, -6, 5, -6, 5);
    // test 2: first base now no longer matches
    let assembly_overlap = "TTCATCTTCTC";
    check_homology(&build_homology("0T2A7", "11", assembly_overlap), "", 0, 0, 0, 11);
    // test 2b: as above but reversed
    check_homology(&build_homology("11", "7T2C", assembly_overlap), "", 0, 0, -11, 0);
    // test 3: first base matches, range of lowest mismatches is 0-1
    check_homology(&build_homology("1T9", "11", assembly_overlap), "T", -1, 0, -1, 10);
    // test 4: second base has mismatch
    check_homology(&build_homology("2C8", "11", assembly_overlap), "TT", -1, 1, -1, 10);
    // test 5: min mismatches in range 4-6
    check_homology(&build_homology("6T4", "3A7", assembly_overlap), "TC", -1, 1, -5, 6);
    // test 6: ref sequences match but have the same difference from the assembly
    check_homology(&build_homology("5C5", "5C5", assembly_overlap), assembly_overlap, -6, 5, -6, 5);
    // test 7
    check_homology(&build_homology("2C0A0T0C1T0C0T", "11", assembly_overlap), "TT", -1, 1, -1, 10);
    // test 8
    check_homology(&build_homology("2^C8", "10", "TTATATTCTC"), "TT", -1, 1, -1, 9);
    // with an even number of bases
    check_homology(&build_homology("2", "2", "AA"), "AA", -1, 1, -1, 1);
}

#[test]
fn test_position_shift() {
    let _env = test_env();
    let homology = HomologyData::new("", -2, 1, -2, 1);
    assert_eq!(-1, homology.position_adjustment(FORWARD));
    assert_eq!(2, homology.position_adjustment(REVERSE));
}

#[test]
fn test_homology_scenario() {
    let _env = test_env();
    let mut ref_genome = MockRefGenome::new();
    ref_genome.put(CHR_1, REF_BASES_400);
    // first a DEL with 3 bases of exact homology - 'TAG'
    let homology_bases = "TAG";
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 121, FORWARD, CHR_1, 202, REVERSE, "", homology_bases);
    let alignment1 = create_alignment_full(CHR_1, 22, 121, false, DEFAULT_MAP_QUAL, 101, 0, 100, "100M50S", "", "100");
    let alignment2 = create_alignment_full(CHR_1, 202, 302, false, DEFAULT_MAP_QUAL, 101, 97, 197, "50S100M", "", "100");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Del, first.sv_type());
    assert_eq!(120, first.position);
    assert_eq!(FORWARD, first.orient);
    assert_eq!(homology_bases, first.homology.homology);
    assert_eq!((-2, 1), (first.homology.exact_start, first.homology.exact_end));
    let second = assembly_alignment.breakends()[1];
    assert_eq!(204, second.position);
    assert_eq!(REVERSE, second.orient);
    assert_eq!((-2, 1), (second.homology.exact_start, second.homology.exact_end));

    // test 2: positive INV
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 121, FORWARD, CHR_1, 340, FORWARD, "", homology_bases);
    let alignment1 = create_alignment_full(CHR_1, 22, 121, false, DEFAULT_MAP_QUAL, 101, 0, 100, "100M50S", "", "100");
    let alignment2 = create_alignment_full(CHR_1, 241, 340, true, DEFAULT_MAP_QUAL, 100, 0, 100, "100M50S", "", "100");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Inv, first.sv_type());
    assert_eq!(120, first.position);
    assert_eq!(FORWARD, first.orient);
    assert_eq!(homology_bases, first.homology.homology);
    assert_eq!((-2, 1), (first.homology.exact_start, first.homology.exact_end));
    let second = assembly_alignment.breakends()[1];
    assert_eq!(338, second.position);
    assert_eq!(FORWARD, second.orient);
    assert_eq!("CTA", second.homology.homology);
    assert_eq!((-1, 2), (second.homology.exact_start, second.homology.exact_end));

    // test 3: negative INV
    let pos1 = 176;
    let pos2 = 280;
    let junction1 = new_junction(CHR_1, pos1, REVERSE);
    let junction2 = new_junction(CHR_1, pos2, REVERSE);
    let ref_length = 100;
    let extension_length = 50;
    let first_ref_bases = &REF_BASES_400[pos1 as usize..(pos1 + ref_length) as usize]; // starts with TAG
    let second_ref_bases = &REF_BASES_400[pos2 as usize..(pos2 + ref_length) as usize]; // runs into CTA
    let first_extension_bases = reverse_complement_str(&second_ref_bases[3..3 + extension_length]);
    let first_assembly_bases = format!("{first_extension_bases}{first_ref_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = JunctionAssembly::new_test(junction1, first_assembly_bases.as_bytes(), &base_quals, extension_length as i32);
    let second_extension_bases = reverse_complement_str(&first_ref_bases[3..3 + extension_length]);
    let second_assembly_bases = format!("{second_extension_bases}{second_ref_bases}");
    let second_assembly = JunctionAssembly::new_test(junction2, second_assembly_bases.as_bytes(), &base_quals, extension_length as i32);
    let homology_bases = "CTA"; // for the second breakend, which means last 3 ref bases of first breakend reads TAG
    let homology_bases_rev = reverse_complement_str(homology_bases);
    let link = try_assembly_overlap(first_assembly, second_assembly).expect("link");
    assert_eq!(homology_bases_rev, link.overlap_bases());
    let assembly_alignment = create_assembly_alignment_link(link);
    let alignment1 = create_alignment_full(CHR_1, pos1, pos1 + ref_length - 1, true, DEFAULT_MAP_QUAL, 100, 97, 197, "50S100M", "", "100");
    let alignment2 = create_alignment_full(CHR_1, pos2, pos2 + ref_length - 1, false, DEFAULT_MAP_QUAL, 100, 97, 197, "50S100M", "", "100");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Inv, first.sv_type());
    assert_eq!(pos1 + 2, first.position);
    assert_eq!(REVERSE, first.orient);
    assert_eq!(homology_bases_rev, first.homology.homology);
    assert_eq!((-2, 1), (first.homology.exact_start, first.homology.exact_end));
    let second = assembly_alignment.breakends()[1];
    assert_eq!(pos2 + 1, second.position);
    assert_eq!(REVERSE, second.orient);
    assert_eq!(homology_bases, second.homology.homology);
    assert_eq!((-1, 2), (second.homology.exact_start, second.homology.exact_end));
}

#[test]
fn test_md_tag_parsing() {
    let _env = test_env();
    let md_tag = MdTag::new("5G3");
    assert_eq!(3, md_tag.elements().len());
    assert_eq!(9, md_tag.base_length());
    assert_eq!(5, md_tag.elements()[0].length);
    assert_eq!(MdTagType::Match, md_tag.elements()[0].type_);
    assert_eq!(MdTagType::Snv, md_tag.elements()[1].type_);
    assert_eq!(3, md_tag.elements()[2].length);
    assert_eq!(MdTagType::Match, md_tag.elements()[2].type_);

    let md_tag = MdTag::new("0A0A3T");
    assert_eq!(4, md_tag.elements().len());
    assert_eq!(6, md_tag.base_length());
    assert_eq!(MdTagType::Snv, md_tag.elements()[0].type_);
    assert_eq!(MdTagType::Snv, md_tag.elements()[1].type_);
    assert_eq!(3, md_tag.elements()[2].length);
    assert_eq!(MdTagType::Snv, md_tag.elements()[3].type_);

    // deletes
    let md_tag = MdTag::new("5^AA0T5^GGG5");
    assert_eq!(6, md_tag.elements().len());
    assert_eq!(21, md_tag.base_length());
    assert_eq!(16, md_tag.sequence_length());
    let expected = [(MdTagType::Match, Some(5)), (MdTagType::Del, Some(2)), (MdTagType::Snv, None), (MdTagType::Match, Some(5)), (MdTagType::Del, Some(3)), (MdTagType::Match, Some(5))];
    for (i, (type_, length)) in expected.iter().enumerate() {
        assert_eq!(*type_, md_tag.elements()[i].type_);
        if let Some(length) = length {
            assert_eq!(*length, md_tag.elements()[i].length);
        }
    }
    assert_eq!(16, md_tag.extract_sub_sequence(0, 15, false).len());

    // a long mismatch sequence
    let md_tag = MdTag::new("12G16T12G1T8C0A6C1C0C9C0C0T0T7T20A0C0C0C7C0A1C5^G1T8A9G3C9C0T8C3A6T3^CTC14C9G13C679");
    assert_eq!(61, md_tag.elements().len());
    assert_eq!(244, md_tag.extract_sub_sequence(0, 243, false).len());
}

#[test]
fn test_indel_insert_exact_homology() {
    let _env = test_env();
    let mut ref_genome = MockRefGenome::new();
    let initial_bases = &REF_BASES_400[0..100];
    let end_bases = &REF_BASES_400[200..300];
    // bases 100 -> 199 - has a 20-base repeated section 4x from 110-190
    let bases_with_homology = "GGGCTGATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTGGTTGATCGT";
    ref_genome.put(CHR_1, &format!("{initial_bases}{bases_with_homology}{end_bases}"));
    // an insert with a duplicated sequence
    let duplicated_bases = &bases_with_homology[10..50];
    // assumes there are 100 ref bases either side of the junction
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 109, FORWARD, CHR_1, 110, REVERSE, duplicated_bases, "");
    let alignment = create_alignment_full(CHR_1, 10, 209, false, DEFAULT_MAP_QUAL, 200, 0, 220, "100M40I100M", "", "");
    let alignments = vec![alignment];
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&alignments);
    assert_eq!(2, assembly_alignment.breakends().len());
    // converted to DUP with breakends at 110:-1, 149:1, inexact homology = exact homology = 0 at start, 40 at end
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Dup, first.sv_type());
    assert_eq!(110, first.position);
    assert_eq!(REVERSE, first.orient);
    check_homology(&first.homology, duplicated_bases, 0, 40, 0, 40);
    let second = assembly_alignment.breakends()[1];
    assert_eq!(149, second.position);
    assert_eq!(FORWARD, second.orient);

    // now with a mismatch mid-way through (10 bases into the duplicated section) - no longer has full exact homology
    let bases_with_homology = "GGGCTGATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCTTGGGTTTTTAAAAACCCCCGGGGGTTTTTGGTTGATCGT";
    ref_genome.put(CHR_1, &format!("{initial_bases}{bases_with_homology}{end_bases}"));
    assembly_alignment.m().breakends_mut().clear();
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&alignments);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Dup, first.sv_type());
    assert_eq!(110, first.position);
    assert_eq!(REVERSE, first.orient);
    check_homology(&first.homology, "AAAAACCCCC", 0, 10, 0, 10);
    let second = assembly_alignment.breakends()[1];
    assert_eq!(149, second.position);
    assert_eq!(FORWARD, second.orient);

    // converts to a partial DUP since full section is not inserted
    let bases_with_homology = "GGGCTGATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCTTTTGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTGGTTGATCGT";
    ref_genome.put(CHR_1, &format!("{initial_bases}{bases_with_homology}{end_bases}"));
    assembly_alignment.m().breakends_mut().clear();
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&alignments);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Dup, first.sv_type());
    assert_eq!(110, first.position);
    assert_eq!(REVERSE, first.orient);
    assert!(first.homology.is_no_homology());
    let second = assembly_alignment.breakends()[1];
    assert_eq!(139, second.position);
    assert_eq!(FORWARD, second.orient);

    // mismatch at first base, remains as an INS
    let bases_with_homology = "GGGCTGATGCTTAAACCCCCGGGGGTTTTTAAAAACCCCCTTGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTGGTTGATCGT";
    ref_genome.put(CHR_1, &format!("{initial_bases}{bases_with_homology}{end_bases}"));
    assembly_alignment.m().breakends_mut().clear();
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&alignments);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Ins, first.sv_type());
    assert_eq!(109, first.position);
    assert_eq!(FORWARD, first.orient);
    assert!(first.homology.is_no_homology());
    let second = assembly_alignment.breakends()[1];
    assert_eq!(110, second.position);
    assert_eq!(REVERSE, second.orient);
}
