//! AlignmentTest
use super::*;
use crate::assembly::alignment::alignment_filters::filter_alignments;
use crate::assembly::alignment::assembly_aligner::{has_longer_minority_extensions, has_suspect_extension, AssemblyAligner};
use crate::assembly::alignment::breakend::Breakend;
use crate::assembly::alignment::breakend_builder::BreakendBuilder;
use crate::assembly::alignment::bwa::{BwaMemAlignment, IBwaMemAligner};
use crate::assembly::alignment::homology_data::HomologyData;
use crate::assembly::constants::{SSX2_MAX_MAP_QUAL, WEAK_ASSEMBLY_DISC_RATE, WEAK_ASSEMBLY_MIN_INS_LENGTH, WEAK_ASSEMBLY_MIN_REF_REPEAT};
use crate::assembly::globals::set_sample_discordant_rate;
use crate::assembly::types::support_read::{SupP, SupportRead};
use crate::assembly::types::StructuralVariantType::{Bnd, Del, Dup, Sgl};
use crate::common::compat;
use crate::common::java_sort::java_sort;

fn ref_genome() -> MockRefGenome {
    let random = ref_bases_random_100();
    let mut ref_genome = MockRefGenome::new();
    ref_genome.put(CHR_1, &format!("{random}{REF_BASES_400}{random}"));
    ref_genome.put(CHR_2, &format!("{random}{REF_BASES_400}{random}")); // matches chr1
    ref_genome.put(CHR_3, &format!("{random}{REF_BASES_200}{random}")); // has much more repetition
    ref_genome
}

fn add_assembly_read(assembly: AsmP, junction_offset: i32) {
    // a negative junction offset translates to soft-clip length
    let soft_clip = if junction_offset < 0 { junction_offset.abs() } else { 0 };
    let read_length = 50;
    let read_bases = &REF_BASES_400[0..read_length as usize];
    let (read_start, cigar) = if assembly.is_forward_junction() {
        (
            assembly.junction().position - junction_offset - read_length + 1,
            if soft_clip > 0 { format!("{}M{}S", read_length - soft_clip, soft_clip) } else { format!("{read_length}M") },
        )
    } else {
        (
            (assembly.junction().position + junction_offset).max(assembly.junction().position),
            if soft_clip > 0 { format!("{}S{}M", soft_clip, read_length - soft_clip) } else { format!("{read_length}M") },
        )
    };
    let read = create_read(&READ_ID_GENERATOR.next_id(), read_start, read_bases, &cigar);
    assembly.m().add_junction_read(&read);
}

/// new AlignData(new ChrBaseRegion(chr, start, end), seqStart, seqEnd, mapQual, score, flags, cigar, DEFAULT_NM, xaTag, mdTag)
fn align_data(chromosome: &str, pos_start: i32, pos_end: i32, seq_start: i32, seq_end: i32, map_qual: i32, score: i32, cigar: &str, xa_tag: Option<&str>, md_tag: &str) -> AlignP {
    P::new(AlignData::new(chromosome, pos_start, pos_end, seq_start, seq_end, map_qual, score, 0, cigar, DEFAULT_NM, xa_tag, Some(md_tag)))
}

fn run_filter(assembly_alignment: AsmAlignP, alignments: &[AlignP]) -> (Vec<AlignP>, Vec<AlignP>) {
    let mut valid_alignments: Vec<AlignP> = Vec::new();
    let mut low_qual_alignments: Vec<AlignP> = Vec::new();
    filter_alignments(assembly_alignment, alignments, &mut valid_alignments, &mut low_qual_alignments);
    (valid_alignments, low_qual_alignments)
}

fn sort_breakends(assembly_alignment: AsmAlignP) {
    java_sort(assembly_alignment.m().breakends_mut(), |a, b| a.compare_to(b));
}

#[test]
fn test_alignment_sequence() {
    let _env = test_env();
    let ref_bases1 = &REF_BASES_400[100..200];
    let ext_bases1 = &REF_BASES_400[300..350];
    let assembly_bases1 = format!("{ref_bases1}{ext_bases1}");
    let assembly1 = create_assembly(CHR_1, 200, FORWARD, &assembly_bases1, ref_bases1.len() as i32 - 1);
    add_assembly_read(assembly1, -30);
    add_assembly_read(assembly1, 10);
    let assembly_alignment = AssemblyAlignment::from_assembly(assembly1);
    assert_eq!(Some(assembly_bases1.as_str()), assembly_alignment.full_sequence());
    assert_eq!(150, assembly_alignment.full_sequence_length());
    assert_eq!(80, assembly1.support()[0].full_assembly_index_start());
    assert_eq!(40, assembly1.support()[1].full_assembly_index_start());

    let assembly_bases1 = format!("{ext_bases1}{ref_bases1}");
    let assembly1 = create_assembly(CHR_1, 100, REVERSE, &assembly_bases1, ext_bases1.len() as i32);
    let assembly_alignment = AssemblyAlignment::from_assembly(assembly1);
    assert_eq!(Some(assembly_bases1.as_str()), assembly_alignment.full_sequence());
    assert_eq!(150, assembly_alignment.full_sequence_length());

    // from a pair of assemblies
    let ref_bases1 = &REF_BASES_400[101..201];
    let ext_bases1 = &REF_BASES_400[300..350];
    let assembly1 = create_assembly(CHR_1, 200, FORWARD, &format!("{ref_bases1}{ext_bases1}"), ref_bases1.len() as i32 - 1);
    let ref_bases2 = &REF_BASES_400[300..400];
    let ext_bases2 = &REF_BASES_400[151..201];
    let assembly2 = create_assembly(CHR_1, 300, REVERSE, &format!("{ext_bases2}{ref_bases2}"), ext_bases2.len() as i32);
    add_assembly_read(assembly2, -30);
    assert_eq!(30, assembly2.support()[0].junction_read_start_distance());
    add_assembly_read(assembly2, 10);
    assert_eq!(-10, assembly2.support()[1].junction_read_start_distance());
    let check_sequence = |assembly1: AsmP, assembly2: AsmP, inserted_bases: &str, overlap: &str, full_sequence: &str| {
        let assembly_link = AssemblyLink::new(assembly1, assembly2, LinkType::Split, inserted_bases, overlap);
        let assembly_alignment = create_assembly_alignment_link(assembly_link);
        assert_eq!(Some(full_sequence), assembly_alignment.full_sequence());
        assert_eq!(full_sequence.len() as i32, assembly_alignment.full_sequence_length());
    };
    check_sequence(assembly1, assembly2, "", "", &format!("{ref_bases1}{ref_bases2}"));
    assert_eq!(70, assembly2.support()[0].full_assembly_index_start());
    assert_eq!(110, assembly2.support()[1].full_assembly_index_start());
    // with overlap
    check_sequence(assembly1, assembly2, "", &ref_bases2[0..30], &format!("{}{}", ref_bases1, &ref_bases2[30..]));
    // with inserted bases
    let inserted_bases = "TTTTT";
    let assembly1 = create_assembly(CHR_1, 200, FORWARD, &format!("{ref_bases1}{inserted_bases}{ext_bases1}"), ref_bases1.len() as i32 - 1);
    let assembly2 = create_assembly(CHR_1, 300, REVERSE, &format!("{ext_bases2}{inserted_bases}{ref_bases2}"), (ext_bases2.len() + inserted_bases.len()) as i32);
    check_sequence(assembly1, assembly2, inserted_bases, "", &format!("{ref_bases1}{inserted_bases}{ref_bases2}"));
    // same again but with matching orientations
    let assembly1 = create_assembly(CHR_1, 200, FORWARD, &format!("{}{}{}", ref_bases1, inserted_bases, reverse_complement_str(ext_bases1)), ref_bases1.len() as i32 - 1);
    let assembly_bases2 = format!("{}{}{}", ref_bases2, reverse_complement_str(inserted_bases), reverse_complement_str(ext_bases2));
    let assembly2 = create_assembly(CHR_1, 349, FORWARD, &assembly_bases2, ref_bases2.len() as i32 - 1);
    check_sequence(assembly1, assembly2, inserted_bases, "", &format!("{}{}{}", ref_bases1, inserted_bases, reverse_complement_str(ref_bases2)));
    // both negative with an overlap
    let overlap = reverse_complement_str(&ref_bases2[0..10]);
    let ext_bases1 = reverse_complement_str(&ref_bases2[10..60]);
    let assembly1 = create_assembly(CHR_1, 101, REVERSE, &format!("{ext_bases1}{ref_bases1}"), ext_bases1.len() as i32);
    let ext_bases2 = reverse_complement_str(&ref_bases1[10..60]);
    let assembly2 = create_assembly(CHR_1, 349, REVERSE, &format!("{ext_bases2}{ref_bases2}"), ext_bases2.len() as i32);
    check_sequence(assembly1, assembly2, "", &overlap, &format!("{}{}", reverse_complement_str(ref_bases1), &ref_bases2[10..]));
}

#[test]
fn test_basic_sv_types() {
    let _env = test_env();
    let ref_genome = ref_genome();
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 200, FORWARD, CHR_2, 250, REVERSE, "", "");
    let alignment1 = create_alignment(CHR_1, 101, 200, 0, 100, "100M50S");
    let alignment2 = create_alignment(CHR_2, 250, 349, 100, 200, "50S100M");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Bnd, first.sv_type());
    assert_eq!(200, first.position);
    assert_eq!(FORWARD, first.orient);
    assert!(!first.homology.exists());
    let second = assembly_alignment.breakends()[1];
    assert_eq!(250, second.position);
    assert_eq!(REVERSE, second.orient);

    // again for a DUP with inserted bases
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 100, REVERSE, CHR_1, 350, FORWARD, "AAAA", "");
    let alignment1 = create_alignment(CHR_1, 251, 350, 0, 100, "100M50S");
    let alignment2 = create_alignment(CHR_1, 100, 199, 104, 204, "50S100M");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
    assert_eq!(2, assembly_alignment.breakends().len());
    sort_breakends(assembly_alignment);
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Dup, first.sv_type());
    assert_eq!(100, first.position);
    assert_eq!(REVERSE, first.orient);
    assert!(!first.homology.exists());
    let second = assembly_alignment.breakends()[1];
    assert_eq!(350, second.position);
    assert_eq!(FORWARD, second.orient);

    // outer singles (the Java test builds zero-qual alignments with alt alignments here but does not pass them in)
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 200, FORWARD, CHR_1, 250, REVERSE, "", "");
    let alignment1 = create_alignment(CHR_1, 101, 200, 2, 101, "50S100M50S");
    let alignment2 = create_alignment(CHR_1, 250, 349, 101, 199, "50S100M50S");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
    assert_eq!(4, assembly_alignment.breakends().len());
    sort_breakends(assembly_alignment);
    let breakend = assembly_alignment.breakends()[0];
    assert_eq!(Sgl, breakend.sv_type());
    assert_eq!(101, breakend.position);
    assert_eq!(REVERSE, breakend.orient);
    assert_eq!(50, breakend.inserted_bases.len());
    let breakend = assembly_alignment.breakends()[1];
    assert_eq!(Del, breakend.sv_type());
    assert_eq!(200, breakend.position);
    assert_eq!(FORWARD, breakend.orient);
    let breakend = assembly_alignment.breakends()[2];
    assert_eq!(Del, breakend.sv_type());
    assert_eq!(250, breakend.position);
    assert_eq!(REVERSE, breakend.orient);
    let breakend = assembly_alignment.breakends()[3];
    assert_eq!(Sgl, breakend.sv_type());
    assert_eq!(349, breakend.position);
    assert_eq!(FORWARD, breakend.orient);
    assert_eq!(50, breakend.inserted_bases.len());
}

#[test]
fn test_indel_types() {
    let _env = test_env();
    let ref_genome = ref_genome();
    // first a DEL
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_3, 200, FORWARD, CHR_3, 251, REVERSE, "", "");
    let alignment = create_alignment(CHR_3, 101, 350, 0, 200, "100M50D100M");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment]);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Del, first.sv_type());
    assert_eq!(200, first.position);
    assert_eq!(FORWARD, first.orient);
    assert!(!first.homology.exists());
    let second = assembly_alignment.breakends()[1];
    assert_eq!(251, second.position);
    assert_eq!(REVERSE, second.orient);

    // a DEL with homology
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_3, 132, FORWARD, CHR_1, 190, REVERSE, "", "");
    let alignment = create_alignment(CHR_3, 80, 250, 0, 120, "53M57D60M");
    BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment]);
    assert_eq!(2, assembly_alignment.breakends().len());
    let first = assembly_alignment.breakends()[0];
    assert_eq!(Del, first.sv_type());
    assert_eq!(134, first.position);
    assert_eq!(FORWARD, first.orient);
    assert_eq!("CCCC", first.homology.homology);
    assert_eq!((-2, 2), (first.homology.exact_start, first.homology.exact_end));
    let second = assembly_alignment.breakends()[1];
    assert_eq!(192, second.position);
    assert_eq!(REVERSE, second.orient);
}

#[test]
fn test_alignment_qual_calcs() {
    let _env = test_env();
    let ref_genome = ref_genome();
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 200, FORWARD, CHR_2, 200, REVERSE, "", "");
    let alignment1 = align_data(CHR_1, 101, 200, 0, 100, 60, 100, "100M", Some(""), "");
    let alignment2 = align_data(CHR_2, 1, 100, 100, 200, 60, 100, "100M", Some(""), "");
    let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, &[alignment1, alignment2]);
    assert_eq!(0, low_qual_alignments.len());
    assert_eq!(2, valid_alignments.len());
    assert_eq!(100, valid_alignments[0].adjusted_alignment());
    assert_eq!(100, valid_alignments[1].adjusted_alignment());
    assert_close(60.0, valid_alignments[0].modified_map_qual(), 0.1);

    // now with lower alignment score and overlap
    let zero_align = align_data(CHR_1, 101, 200, 0, 1, 0, 1, "1M", Some(""), "");
    let alignment1 = align_data(CHR_1, 101, 200, 0, 80, 60, 80, "110M", Some(""), "");
    let alignment2 = align_data(CHR_2, 201, 300, 60, 140, 100, 80, "100M", Some(""), "");
    let alignment3 = align_data(CHR_3, 1, 100, 130, 200, 60, 70, "100M", Some(""), "");
    let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, &[zero_align, alignment1, zero_align, alignment2, zero_align, alignment3]);
    assert_eq!(4, low_qual_alignments.len());
    assert_eq!(2, valid_alignments.len());
    assert_eq!(60, alignment1.adjusted_alignment());
    assert_eq!(60, alignment3.adjusted_alignment());
    assert_eq!(48, alignment2.adjusted_alignment());
}

#[test]
fn test_alignment_short_adjusted_aligned_length_filter() {
    let _env = test_env();
    let ref_genome = ref_genome();
    let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 200, FORWARD, CHR_2, 200, REVERSE, "", "");
    // exclusion of short adjusted alignment-length segments
    let alt_alignment = "2,+200000,100M,0";
    let alignment1 = align_data(CHR_1, 101, 200, 0, 80, 60, 80, "80M", Some(""), "");
    let alignment2 = align_data(CHR_2, 201, 220, 80, 100, 60, 20, "20M", Some(alt_alignment), "");
    let alignment3 = align_data(CHR_3, 1, 100, 100, 200, 60, 100, "100M", Some(""), "");
    let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, &[alignment1, alignment2, alignment3]);
    assert_eq!(1, low_qual_alignments.len());
    assert_eq!(2, valid_alignments.len());

    // short alignments are permitted on the outside
    let alignment1 = align_data(CHR_1, 101, 120, 0, 20, 60, 80, "20M", Some(alt_alignment), "");
    let alignment2 = align_data(CHR_2, 201, 300, 20, 100, 60, 80, "100M", Some(""), "");
    let alignment3 = align_data(CHR_2, 350, 370, 100, 120, 60, 100, "20M", Some(alt_alignment), "");
    let (valid_alignments, _) = run_filter(assembly_alignment, &[alignment1, alignment2, alignment3]);
    assert_eq!(3, valid_alignments.len());
}

#[test]
fn test_alignment_low_map_qual_filters() {
    let _env = test_env();
    let ref_genome = ref_genome();
    let assembly_alignment = create_assembly_alignment_ext(&ref_genome, CHR_1, 300, FORWARD, CHR_1, 350, REVERSE, "", "", 50, 200);
    // scenario 1
    let cigar = "100M";
    let alt_alignment1 = "2,+200000,100M,0";
    let alt_alignment2 = "1,+2000,100M,0";
    let alignment1 = align_data(CHR_1, 101, 200, 0, 100, 0, 100, cigar, Some(alt_alignment1), "");
    let alignment2 = align_data(CHR_1, 20000, 20100, 101, 200, 0, 100, cigar, Some(alt_alignment2), "");
    let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, &[alignment1, alignment2]);
    assert_eq!(0, low_qual_alignments.len());
    assert_eq!(2, valid_alignments.len());
    assert!(alignment1.has_low_map_qual_alignment());
    assert!(alignment2.has_low_map_qual_alignment());
    assert!(alignment1.selected_alt_alignment().is_none());
    assert_eq!(1, alignment1.unselected_alt_alignments().len());
    assert_eq!(2000, alignment2.selected_alt_alignment().unwrap().position);
    assert!(alignment1.has_low_map_qual_short_sv_link());
    assert!(alignment2.has_low_map_qual_short_sv_link());

    // scenario 2
    let alignment1 = align_data(CHR_1, 101, 200, 0, 100, 60, 100, cigar, Some(alt_alignment1), "");
    let alignment2 = align_data(CHR_1, 201, 300, 101, 200, 60, 100, cigar, Some(alt_alignment2), "");
    let alignment3 = align_data(CHR_2, 101, 200, 201, 300, 0, 100, cigar, Some("2,+200000,100M,0"), "");
    let alignment4 = align_data(CHR_2, 20000, 20100, 301, 400, 0, 100, cigar, Some("1,+2000,100M,0"), "");
    let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, &[alignment1, alignment2, alignment3, alignment4]);
    assert_eq!(0, low_qual_alignments.len());
    assert_eq!(4, valid_alignments.len());
    for alignment in [alignment3, alignment4] {
        assert!(alignment.has_low_map_qual_alignment());
        assert!(alignment.selected_alt_alignment().is_none());
        assert_eq!(1, alignment.unselected_alt_alignments().len());
        assert!(alignment.has_low_map_qual_short_sv_link());
    }
}

#[test]
fn test_keep_adjacent_local_short_adjusted_alignments() {
    let _env = test_env();
    let ref_genome = ref_genome();
    let assembly_alignment = create_assembly_alignment_ext(&ref_genome, CHR_1, 300, FORWARD, CHR_1, 350, REVERSE, "", "", 50, 200);
    let cigar = "100M";
    let middle_cigar = "50M";
    let alignment1 = align_data(CHR_1, 101, 200, 0, 100, 60, 100, cigar, None, "");
    let alignment3 = align_data(CHR_1, 20000, 20100, 201, 300, 60, 100, cigar, None, "");
    let check = |alignments: &[AlignP], low_qual: usize, valid: usize| {
        let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, alignments);
        assert_eq!(low_qual, low_qual_alignments.len());
        assert_eq!(valid, valid_alignments.len());
    };
    let alignment2 = align_data(CHR_1, 1101, 1150, 101, 150, 60, 40, middle_cigar, None, "");
    check(&[alignment1, alignment2, alignment3], 0, 3);
    // rescued by upper alignment
    let alignment2 = align_data(CHR_1, 19101, 19150, 101, 150, 60, 40, middle_cigar, None, "");
    check(&[alignment1, alignment2, alignment3], 0, 3);
    // no longer rescued if remote or too far away
    let alignment2 = align_data(CHR_2, 1101, 1150, 101, 150, 60, 40, middle_cigar, None, "");
    check(&[alignment1, alignment2, alignment3], 1, 2);
    let alignment2 = align_data(CHR_1, 10000, 10049, 101, 150, 60, 40, middle_cigar, None, "");
    check(&[alignment1, alignment2, alignment3], 1, 2);
    // only needs 1 anchoring high-qual alignment
    let alignment2 = align_data(CHR_1, 19101, 19150, 101, 150, 60, 40, middle_cigar, None, "");
    check(&[alignment2, alignment3], 0, 2);
    // not rescued if has too many mismatches
    let alignment2 = align_data(CHR_1, 19101, 19150, 101, 150, 60, 40, "20M10I20M", None, "10A10C10G10T10");
    check(&[alignment2, alignment3], 1, 1);
}

#[test]
fn test_specific_alignment_handling() {
    let _env = test_env();
    let ref_genome = ref_genome();
    let assembly_alignment = create_assembly_alignment_ext(&ref_genome, CHR_1, 300, FORWARD, CHR_1, 350, REVERSE, "", "", 50, 200);
    let cigar = "100M";
    let alignment1 = align_data(CHR_1, 101, 200, 0, 100, 60, 100, cigar, Some(""), "");
    let low_map_qual = SSX2_MAX_MAP_QUAL - 1;
    // the first version 37 SSX2 region of the Java constants, with the gene's (reverse) orientation
    let (ssx2_chromosome, ssx2_start, ssx2_end) = ("X", 52729628, 52731680);
    let alt_alignment2 = format!("{},-{},100M,{}", ssx2_chromosome, ssx2_start + 1, low_map_qual);
    let alignment2 = align_data(CHR_2, 20000, 20100, 101, 200, 0, 100, cigar, Some(&alt_alignment2), "");
    let (valid_alignments, low_qual_alignments) = run_filter(assembly_alignment, &[alignment1, alignment2]);
    assert_eq!(0, low_qual_alignments.len());
    assert_eq!(2, valid_alignments.len());
    let specific_alignment = valid_alignments.iter().find(|x| x.chromosome() == ssx2_chromosome && x.position_start() <= ssx2_end && x.position_end() >= ssx2_start);
    assert!(specific_alignment.is_none());
    assert!(alignment2.has_selected_alt_alignment());
    let ssx2_alt_alignment = alignment2.selected_alt_alignment().unwrap();
    assert_eq!(ssx2_chromosome, ssx2_alt_alignment.chromosome);
    assert!(ssx2_alt_alignment.position >= ssx2_start && ssx2_alt_alignment.position <= ssx2_end);
}

/// hmf-common MockBwaMemAligner.forAnyQuery
struct MockBwaMemAligner {
    alignments: Vec<BwaMemAlignment>,
}

impl IBwaMemAligner for MockBwaMemAligner {
    fn align_sequences(&self, sequences: &[&[u8]]) -> Vec<Vec<BwaMemAlignment>> {
        sequences.iter().map(|_| self.alignments.clone()).collect()
    }
}

fn create_bwa_alignment(chr_ref_id: i32, ref_pos_start: i32, ref_pos_end: i32, flags: i32, seq_start: i32, seq_end: i32, cigar: &str, md_tag: &str, xa_tag: &str, map_qual: i32, n_mismatches: i32, aligner_score: i32) -> BwaMemAlignment {
    BwaMemAlignment {
        sam_flag: flags,
        ref_id: chr_ref_id,
        ref_start: ref_pos_start - 1,
        ref_end: ref_pos_end,
        seq_start,
        seq_end,
        map_qual,
        n_mismatches,
        aligner_score,
        suboptimal_score: 0,
        cigar: cigar.to_string(),
        md_tag: Some(md_tag.to_string()),
        xa_tag: Some(xa_tag.to_string()),
        mate_ref_id: 0,
        mate_ref_start: 0,
        template_len: 0,
    }
}

#[test]
fn test_alignment_requery() {
    let env = test_env();
    let ref_genome = ref_genome();
    // need to ensure the full assembly sequence length matches the alignments returned
    let assembly_alignment = create_assembly_alignment_ext(&ref_genome, CHR_1, 200, FORWARD, CHR_2, 200, REVERSE, "", "", 50, 200);
    // the Java test also lists two initial alignments with a soft-clip, but its mock returns the requery alignment for every query
    let requery_alignments = vec![create_bwa_alignment(0, 210, 259, 0, 0, 49, "50M", "", "", 60, 0, 50)];
    let aligner = MockBwaMemAligner { alignments: requery_alignments };
    let mut assembly_aligner = AssemblyAligner::new(&env.config, &aligner);
    assert!(assembly_aligner.process_assembly(assembly_alignment).is_some());
    AssemblyAligner::allocate_support(&env.config, assembly_alignment, &env.config.combined_sample_ids());
}

#[test]
fn test_longer_minority_extensions() {
    let _env = test_env();
    let read = create_read(&READ_ID_GENERATOR.next_id(), 100, "", "");
    let support = |distance: i32| -> SupP { SupportRead::new(&read, SupportType::Junction, distance, 0, 10) };
    let mut reads: Vec<SupP> = vec![support(10), support(10), support(10), support(24)];
    assert!(!has_longer_minority_extensions(&reads, REVERSE, false)); // ratio too low
    reads = vec![support(10), support(10), support(10), support(35)];
    assert!(has_longer_minority_extensions(&reads, REVERSE, false));
    reads.push(support(35));
    assert!(!has_longer_minority_extensions(&reads, REVERSE, false)); // too many longer extending reads
    for _ in 0..20 {
        reads.push(support(10));
    }
    assert!(has_longer_minority_extensions(&reads, REVERSE, false));
    assert!(has_longer_minority_extensions(&reads, REVERSE, true));
}

#[test]
fn test_suspect_extensions() {
    let _env = test_env();
    let ref_bases1 = format!("{}{}", &REF_BASES_400[100..150], "AT".repeat(WEAK_ASSEMBLY_MIN_REF_REPEAT as usize));
    let ext_bases1 = &REF_BASES_400[300..350];
    let assembly = create_assembly(CHR_1, 200, FORWARD, &format!("{ref_bases1}{ext_bases1}"), ref_bases1.len() as i32 - 1);
    let assembly_alignment = AssemblyAlignment::from_assembly(assembly);
    let position_start = 10;
    let position_end = 100;
    let inserted_bases = "A".repeat(WEAK_ASSEMBLY_MIN_INS_LENGTH as usize);
    let homology_data = HomologyData::new("", 0, 0, 5, 2);
    let add_breakends = |inserted_bases: &str| {
        let breakend_start = Breakend::new(assembly_alignment, CHR_1, position_start, FORWARD, inserted_bases, Some(homology_data.clone()));
        let breakend_end = Breakend::new(assembly_alignment, CHR_2, position_end, FORWARD, inserted_bases, Some(homology_data.clone()));
        breakend_start.m().set_other_breakend(Some(breakend_end));
        breakend_end.m().set_other_breakend(Some(breakend_start));
        assembly_alignment.m().add_breakend(breakend_start);
        assembly_alignment.m().add_breakend(breakend_end);
    };
    add_breakends(&inserted_bases);
    // condition 1 - homology and discordant rate
    set_sample_discordant_rate(WEAK_ASSEMBLY_DISC_RATE + 0.01);
    assert!(has_suspect_extension(assembly_alignment));
    // condition 2 - insert length
    set_sample_discordant_rate(0.0);
    assert!(has_suspect_extension(assembly_alignment));
    // condition 3 - ref repeat
    assembly_alignment.m().breakends_mut().clear();
    add_breakends("");
    assert!(has_suspect_extension(assembly_alignment));
}

/// (chromosome, position, orientation) of each breakend, Java and corrected
type BreakendCoords = Vec<(String, i32, Orientation)>;

fn breakend_coords(assembly_alignment: AsmAlignP) -> BreakendCoords {
    assembly_alignment.breakends().iter().map(|b| (b.chromosome.to_string(), b.position, b.orient)).collect()
}

#[test]
fn a38_pms2_alternative_gives_breakend_coordinates() {
    let _env = test_env();
    let ref_genome = ref_genome();
    // a segment in PMS2CL (version 37, chromosome 7 at 6.76 Mb) with map quality 0 and an alternative in
    // PMS2 at 7:6010001-6010100, beside a segment of chromosome 1; `pms2cl_first` puts it first
    let breakends = |pms2cl_first: bool, alt_strand: char| -> (BreakendCoords, BreakendCoords) {
        let run = || {
            let assembly_alignment = create_assembly_alignment(&ref_genome, CHR_1, 200, FORWARD, CHR_2, 250, REVERSE, "", "");
            let xa = format!("7,{alt_strand}6010001,100M,0;");
            let (alignment1, alignment2) = if pms2cl_first {
                (align_data("7", 6760001, 6760100, 0, 100, 0, 100, "100M100S", Some(&xa), ""), create_alignment(CHR_1, 101, 200, 100, 200, "100S100M"))
            } else {
                (create_alignment(CHR_1, 101, 200, 0, 100, "100M100S"), align_data("7", 6760001, 6760100, 100, 200, 0, 100, "100S100M", Some(&xa), ""))
            };
            BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
            breakend_coords(assembly_alignment)
        };
        (run(), compat::with_java(true, run))
    };
    let pms2 = |position: i32, orient: Orientation| ("7".to_string(), position, orient);
    let chr1 = |position: i32, orient: Orientation| (CHR_1.to_string(), position, orient);
    // second segment, alternative on the + strand: the breakend is at its start, facing reverse
    let (corrected, java) = breakends(false, '+');
    assert_eq!(vec![chr1(200, FORWARD), pms2(6010001, REVERSE)], corrected);
    assert_eq!(vec![chr1(200, FORWARD), pms2(6010001, FORWARD)], java);
    // on the - strand: at its last base, facing forward
    let (corrected, java) = breakends(false, '-');
    assert_eq!(vec![chr1(200, FORWARD), pms2(6010100, FORWARD)], corrected);
    assert_eq!(vec![chr1(200, FORWARD), pms2(6010001, REVERSE)], java);
    // first segment, + strand: at its last base, facing forward
    let (corrected, java) = breakends(true, '+');
    assert_eq!(vec![pms2(6010100, FORWARD), chr1(101, REVERSE)], corrected);
    assert_eq!(vec![pms2(6010001, FORWARD), chr1(101, REVERSE)], java);
    // first segment, - strand: at its start, facing reverse, which Java gets right
    let (corrected, java) = breakends(true, '-');
    assert_eq!(vec![pms2(6010001, REVERSE), chr1(101, REVERSE)], corrected);
    assert_eq!(corrected, java);
}

#[test]
fn a41_indel_anchors_are_the_segment_bases_either_side() {
    let _env = test_env();
    let ref_genome = ref_genome();
    // a second segment at 100-399 of a 400-base sequence, with a 40-base deletion; the first segment
    // (chromosome 1) makes a BND with it. Counts the DEL breakends formed from the cigar.
    let del_breakends = |cigar: &str| {
        let run = || {
            let assembly_alignment = create_assembly_alignment_ext(&ref_genome, CHR_1, 200, FORWARD, CHR_2, 250, REVERSE, "", "", 200, 200);
            assert_eq!(400, assembly_alignment.full_sequence_length());
            let alignment1 = create_alignment(CHR_1, 101, 200, 0, 100, "100M300S");
            let alignment2 = create_alignment(CHR_2, 150, 489, 100, 400, cigar);
            BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment1, alignment2]);
            assembly_alignment.breakends().iter().filter(|x| x.sv_type() == Del).count()
        };
        (run(), compat::with_java(true, run))
    };
    // 240 aligned bases before the deletion and 60 after: Java takes the right anchor as 300 - 340
    assert_eq!((2, 0), del_breakends("100S240M40D60M"));
    // 30 before and 270 after: Java takes the left anchor as 130, from the sequence start
    assert_eq!((0, 2), del_breakends("100S30M40D270M"));
    // 150 either side: called in both
    assert_eq!((2, 2), del_breakends("100S150M40D150M"));
}

#[test]
fn a41_reverse_alignment_insert_is_read_from_its_own_place() {
    let _env = test_env();
    let ref_genome = ref_genome();
    // along the reference strand the query is R1 (100 bases of 1:101-200), an inserted I (40), R2 (160 of
    // 1:201-360); the assembly sequence is its reverse complement, aligned to the reverse strand
    let r1 = ref_genome.base_string(CHR_1, 101, 200);
    let r2 = ref_genome.base_string(CHR_1, 201, 360);
    let insert = &REF_BASES_200[..40];
    let query = format!("{r1}{insert}{r2}");
    let sequence = reverse_complement_str(&query);
    let indel = || {
        let assembly = create_assembly(CHR_1, 200, FORWARD, &sequence, 99);
        let assembly_alignment = AssemblyAlignment::from_assembly(assembly);
        let alignment = create_alignment_full(CHR_1, 101, 360, true, DEFAULT_MAP_QUAL, 300, 0, 300, "100M40I160M", "", "");
        BreakendBuilder::new(Some(&ref_genome), assembly_alignment).form_breakends(&[alignment]);
        let breakends = assembly_alignment.breakends();
        assert_eq!(2, breakends.len());
        (breakends[0].inserted_bases.clone(), breakends[0].segments()[0].indel_sequence_indices().unwrap())
    };
    let (corrected, java) = (indel(), compat::with_java(true, indel));
    // in the assembly sequence the insert is at 160-199, reverse complemented
    assert_eq!((insert.to_string(), [160, 200]), corrected);
    assert_eq!(reverse_complement_str(insert), sequence[160..200]);
    // Java takes 100-139, bases of R2
    assert_eq!([100, 140], java.1);
    assert_ne!(insert, java.0);
}
