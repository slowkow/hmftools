//! NEO's PmNeoEpitopeTest, ported

use super::*;
use crate::cat;
use crate::epitope::*;
use hmf_common::neo::{NeoEpitopeFusion, NeoEpitopeType};

#[test]
fn test_snv_point_mutations() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(120);
    ref_genome.put(CHR_1, &chr1_bases);

    // tests: SNP, delete, insertion

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(85);

    let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    let pm_data = create_point_mutation(42, &sub(&chr1_bases, 42, 43), &"A");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(42, ne_data.position(FS_UP));
    assert_eq!(43, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_DOWN));

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);
    assert_eq!(3, ne_data.exon_rank[FS_UP]);
    assert_eq!(3, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 40, 42), pm_data.alt);
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 43, 51), sub(&chr1_bases, 60, 61));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // SNP requiring all bases plus for phasing
    let pm_data = create_point_mutation(46, &sub(&chr1_bases, 46, 47), &"A");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(46, ne_data.position(FS_UP));
    assert_eq!(47, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_DOWN));

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert_eq!(3, ne_data.exon_rank[FS_UP]);
    assert_eq!(3, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 46));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(pm_data.alt, sub(&chr1_bases, 47, 49));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 49, 51), sub(&chr1_bases, 60, 67));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // point mutation on last base of exon
    let pm_data = create_point_mutation(30, &sub(&chr1_bases, 30, 31), &"A");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(30, ne_data.position(FS_UP));
    assert_eq!(31, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_DOWN));

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Intronic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert_eq!(2, ne_data.exon_rank[FS_UP]);
    assert_eq!(3, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 28);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 28, 30), pm_data.alt);
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = sub(&chr1_bases, 40, 49);
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);


    // test again on the reverse strand
    let coding_start = Some(25);
    let coding_end = Some(65);
    let trans_data_neg_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    let pm_data = create_point_mutation(49, &sub(&chr1_bases, 49, 50), &"A");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(49, ne_data.position(FS_UP));
    assert_eq!(48, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_DOWN));

    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert_eq!(4, ne_data.exon_rank[FS_UP]);
    assert_eq!(4, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&sub(&chr1_bases, 60, 66));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 48, 49), pm_data.alt, sub(&chr1_bases, 50, 51)));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 30, 31), sub(&chr1_bases, 40, 48)));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // point mutation on last base of exon
    let pm_data = create_point_mutation(40, &sub(&chr1_bases, 40, 41), &"A");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(40, ne_data.position(FS_UP));
    assert_eq!(39, ne_data.position(FS_DOWN));
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Intronic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&sub(&chr1_bases, 42, 51));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 30, 31), pm_data.alt, sub(&chr1_bases, 41, 42)));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&sub(&chr1_bases, 25, 30));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);
}

#[test]
fn test_mnv_point_mutations() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(120);

    ref_genome.put(CHR_1, &chr1_bases);

    // tests: SNP, delete, insertion

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(85);

    let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // mutation of 2 bases at first base of codon
    let pm_data = create_point_mutation(43, &sub(&chr1_bases, 43, 45), &"GG");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(43, ne_data.position(FS_UP));
    assert_eq!(45, ne_data.position(FS_DOWN));

    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 43));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(pm_data.alt, sub(&chr1_bases, 45, 46));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 46, 51), sub(&chr1_bases, 60, 64));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // 2-base MNV mutation on last base of exon, making 1 novel codon
    let pm_data = create_point_mutation(30, &sub(&chr1_bases, 30, 32), &"AA");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(30, ne_data.position(FS_UP));
    assert_eq!(32, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 28);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    // FIXME:
    let novel_bases = cat!(sub(&chr1_bases, 28, 30), sub(&pm_data.alt, 0, 1));
    // assertEquals(novelBases, neData.NovelCodonBases);

    let down_bases = sub(&chr1_bases, 40, 49);
    //assertEquals(downBases, neData.CodingBases[FS_DOWN]);


    // 5-base MNV mutation on second base of exon, making 3 novel codons
    let pm_data = create_point_mutation(42, &sub(&chr1_bases, 42, 47), &"AAAAA");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(42, ne_data.position(FS_UP));
    assert_eq!(47, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 40, 42), pm_data.alt, sub(&chr1_bases, 47, 49));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 49, 51), sub(&chr1_bases, 60, 67));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // /*
    // // test again on the reverse strand
    // codingStart = 25;
    // codingEnd = 65;
    // TranscriptData transDataNegStrand = createTransExons(
    // GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");
    // 
    // pmData = createPointMutation(
    // 49, chr1Bases.substring(49, 50), "A", NONSENSE_OR_FRAMESHIFT);
    // 
    // neData = new PmNeoEpitope(pmData);
    // 
    // neData.setTranscriptData(transDataNegStrand, transDataNegStrand);
    // 
    // assertEquals(49, neData.position(FS_UP));
    // assertEquals(48, neData.position(FS_DOWN));
    // assertEquals(ORIENT_REV, neData.orientation(FS_UP));
    // assertEquals(ORIENT_FWD, neData.orientation(FS_DOWN));
    // 
    // assertEquals(CODING, neData.CodingType[FS_UP]);
    // assertEquals(CODING, neData.CodingType[FS_DOWN]);
    // assertEquals(EXONIC, neData.RegionType[FS_UP]);
    // assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
    // assertEquals(PHASE_2, neData.Phases[FS_UP]);
    // assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
    // assertEquals(4, neData.ExonRank[FS_UP]);
    // assertEquals(4, neData.ExonRank[FS_DOWN]);
    // Assert.assertTrue(neData.phaseMatched());
    // 
    // neData.setCodingBases(refGenome, 3);
    // 
    // upBases = reverseStrandBases(chr1Bases.substring(60, 66));
    // assertEquals(upBases, neData.CodingBases[FS_UP]);
    // 
    // novelBases = reverseStrandBases(chr1Bases.substring(48, 49) + pmData.Alt + chr1Bases.substring(50, 51));
    // assertEquals(novelBases, neData.NovelCodonBases);
    // 
    // downBases = reverseStrandBases(chr1Bases.substring(30, 31) + chr1Bases.substring(40, 48));
    // assertEquals(downBases, neData.CodingBases[FS_DOWN]);
    // 
    // // point mutation on last base of exon
    // pmData = createPointMutation(
    // 40, chr1Bases.substring(40, 41), "A", NONSENSE_OR_FRAMESHIFT);
    // 
    // neData = new PmNeoEpitope(pmData);
    // 
    // neData.setTranscriptData(transDataNegStrand, transDataNegStrand);
    // 
    // assertEquals(40, neData.position(FS_UP));
    // assertEquals(39, neData.position(FS_DOWN));
    // assertEquals(ORIENT_REV, neData.orientation(FS_UP));
    // assertEquals(ORIENT_FWD, neData.orientation(FS_DOWN));
    // 
    // assertEquals(CODING, neData.CodingType[FS_UP]);
    // assertEquals(CODING, neData.CodingType[FS_DOWN]);
    // assertEquals(EXONIC, neData.RegionType[FS_UP]);
    // assertEquals(INTRONIC, neData.RegionType[FS_DOWN]);
    // assertEquals(PHASE_2, neData.Phases[FS_UP]);
    // assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
    // assertEquals(4, neData.ExonRank[FS_UP]);
    // assertEquals(3, neData.ExonRank[FS_DOWN]);
    // Assert.assertTrue(neData.phaseMatched());
    // 
    // neData.setCodingBases(refGenome, 3);
    // 
    // upBases = reverseStrandBases(chr1Bases.substring(42, 51));
    // assertEquals(upBases, neData.CodingBases[FS_UP]);
    // 
    // novelBases = reverseStrandBases(chr1Bases.substring(30, 31) + pmData.Alt + chr1Bases.substring(41, 42));
    // assertEquals(novelBases, neData.NovelCodonBases);
    // 
    // downBases = reverseStrandBases(chr1Bases.substring(25, 30));
    // assertEquals(downBases, neData.CodingBases[FS_DOWN]);
    // */
}

#[test]
fn test_deletion_point_mutations() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(120);

    ref_genome.put(CHR_1, &chr1_bases);

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(85);

    let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // delete - of 3 bases, in phase
    let pm_data = create_point_mutation(41, &sub(&chr1_bases, 41, 45), &sub(&chr1_bases, 41, 42));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(41, ne_data.position(FS_UP));
    assert_eq!(45, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_DOWN));

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert_eq!(3, ne_data.exon_rank[FS_UP]);
    assert_eq!(3, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 40, 41), pm_data.alt, sub(&chr1_bases, 45, 46));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 46, 51), sub(&chr1_bases, 60, 64));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // delete - of 2 bases, out of phase
    let pm_data = create_point_mutation(46, &sub(&chr1_bases, 46, 49), &sub(&chr1_bases, 46, 47));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(46, ne_data.position(FS_UP));
    assert_eq!(49, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 46));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 46, 47), sub(&chr1_bases, 49, 51), sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);

    // delete - of 2 bases, out of phase, starting mnid codon
    let pm_data = create_point_mutation(47, &sub(&chr1_bases, 47, 50), &sub(&chr1_bases, 47, 48));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(47, ne_data.position(FS_UP));
    assert_eq!(50, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 46));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 46, 48), sub(&chr1_bases, 50, 51), sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);

    // delete - of 2 bases, out of phase, but deleting a new codon's first 2 bases
    let pm_data = create_point_mutation(42, &sub(&chr1_bases, 42, 45), &sub(&chr1_bases, 42, 43));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(42, ne_data.position(FS_UP));
    assert_eq!(45, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 43));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 45, 51), sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);

    // same again but with enough room to get required upstream bases
    let pm_data = create_point_mutation(45, &sub(&chr1_bases, 45, 48), &sub(&chr1_bases, 45, 46));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(45, ne_data.position(FS_UP));
    assert_eq!(48, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 46));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 48, 51), sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);


    // delete of 3 bases at last base of codon, no novel codon, in phase
    let pm_data = create_point_mutation(42, &sub(&chr1_bases, 42, 46), &sub(&chr1_bases, 42, 43));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(42, ne_data.position(FS_UP));
    assert_eq!(46, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 43));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    assert_eq!("", ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 46, 51), sub(&chr1_bases, 60, 64));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    ne_data.set_amino_acids(&ref_genome, 3).unwrap();
    let up_wildtype_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 51), sub(&chr1_bases, 60, 61));
    let up_wild_aas = get_amino_acids(&up_wildtype_bases, true);
    assert_eq!(up_wild_aas, ne_data.upstream_wild_type_acids);

    // same again but allowing for enough coding bases upstream to satisfy the required count
    let pm_data = create_point_mutation(45, &sub(&chr1_bases, 45, 49), &sub(&chr1_bases, 45, 46));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(45, ne_data.position(FS_UP));
    assert_eq!(49, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 46));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    assert_eq!("", ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 49, 51), sub(&chr1_bases, 60, 67));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);


    // delete on negative strand
    let coding_start = Some(25);
    let coding_end = Some(65);
    let trans_data_neg_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // 3 base DEL - deleting a codon
    let pm_data = create_point_mutation(41, &sub(&chr1_bases, 41, 45), &sub(&chr1_bases, 41, 42));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(45, ne_data.position(FS_UP));
    assert_eq!(41, ne_data.position(FS_DOWN));

    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 45, 51), sub(&chr1_bases, 60, 63)));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    assert_eq!("", ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 42)));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // 3 base DEL - crossing a codon boundary
    let pm_data = create_point_mutation(42, &sub(&chr1_bases, 42, 46), &sub(&chr1_bases, 42, 43));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(46, ne_data.position(FS_UP));
    assert_eq!(42, ne_data.position(FS_DOWN));

    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 48, 51), sub(&chr1_bases, 60, 66)));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 42, 43), sub(&chr1_bases, 46, 48)));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 42)));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    ne_data.set_amino_acids(&ref_genome, 3).unwrap();
    let up_wildtype_bases = cat!(sub(&chr1_bases, 30, 31), sub(&chr1_bases, 40, 51), sub(&chr1_bases, 60, 66));
    let up_wild_aas = get_amino_acids(&reverse_complement_bases(&up_wildtype_bases), true);
    assert_eq!(up_wild_aas, ne_data.upstream_wild_type_acids);

    // again but a base up
    let pm_data = create_point_mutation(43, &sub(&chr1_bases, 43, 47), &sub(&chr1_bases, 43, 44));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(47, ne_data.position(FS_UP));
    assert_eq!(43, ne_data.position(FS_DOWN));

    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 48, 51), sub(&chr1_bases, 60, 66)));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 42, 44), sub(&chr1_bases, 47, 48)));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 42)));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);
}

#[test]
fn test_insert_point_mutations() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(120);

    ref_genome.put(CHR_1, &chr1_bases);

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(85);

    let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // insert of 2 bases at first base of codon - frameshift
    let pm_data = create_point_mutation(28, &sub(&chr1_bases, 28, 29), &"AAA");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(28, ne_data.position(FS_UP));
    assert_eq!(29, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_DOWN));

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);  // whatever it would usually be
    assert_eq!(2, ne_data.exon_rank[FS_UP]);
    assert_eq!(2, ne_data.exon_rank[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 28);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let down_bases = cat!(pm_data.alt, sub(&chr1_bases, 29, 31), sub(&chr1_bases, 40, 51), sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);
    assert_eq!(down_bases, ne_data.novel_codon_bases);  // due to lack of phasing

    // insert of 3 bases at last base of codon - making 1 novel codon, starting after first base and in phase
    let pm_data = create_point_mutation(42, &sub(&chr1_bases, 42, 43), &cat!(sub(&chr1_bases, 42, 43), "AAA"));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(42, ne_data.position(FS_UP));
    assert_eq!(43, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);  // pushed out 3 by the inserted bases, unch
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 43));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = sub(&pm_data.alt, 1, pm_data.alt.len() as i32);
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 43, 51), sub(&chr1_bases, 60, 61));

    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // insert of 3 bases at 2nd base of codon - 2 novel codons, starting after first base and in phase
    let pm_data = create_point_mutation(41, &sub(&chr1_bases, 41, 42), &cat!(sub(&chr1_bases, 41, 42), "AAA"));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(41, ne_data.position(FS_UP));
    assert_eq!(42, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);  // pushed out 6 by the inserted bases, unch
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 40, 41), pm_data.alt, sub(&chr1_bases, 42, 43));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = cat!(sub(&chr1_bases, 43, 51), sub(&chr1_bases, 60, 61));

    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // negative strand

    // insert of 3 bases, in phase
    let coding_start = Some(25);
    let coding_end = Some(65);
    let trans_data_neg_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    let pm_data = create_point_mutation(43, &sub(&chr1_bases, 43, 44), &"AAAA");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(44, ne_data.position(FS_UP));
    assert_eq!(43, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);  // pushed out 2 by the inserted bases
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);  // whatever it would usually be
    assert_eq!(4, ne_data.exon_rank[FS_UP]);
    assert_eq!(4, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 45, 51), sub(&chr1_bases, 60, 63)));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 42, 43), pm_data.alt, sub(&chr1_bases, 44, 45)));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 42)));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // insert of 3 bases, in phase, with inserted bases forming the new codon
    let pm_data = create_point_mutation(44, &sub(&chr1_bases, 44, 45), &cat!(sub(&chr1_bases, 44, 45), "AAA"));

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_neg_strand, trans_data_neg_strand);

    assert_eq!(45, ne_data.position(FS_UP));
    assert_eq!(44, ne_data.position(FS_DOWN));
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);  // pushed out 2 by the inserted bases
    assert_eq!(PHASE_1, ne_data.phases[FS_DOWN]);  // whatever it would usually be
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 45, 51), sub(&chr1_bases, 60, 63)));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 42, 44), pm_data.alt));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = reverse_complement_bases(&cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 42)));
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);
}

#[test]
fn test_stop_lost_point_mutations() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = cat!(generate_random_bases(64), "TGA", generate_random_bases(60));

    ref_genome.put(CHR_1, &chr1_bases);

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(66);

    let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // stop-lost
    let pm_data = create_point_mutation(64, &sub(&chr1_bases, 64, 65), &"A");

    let mut ne_data = NeoEpitope::point_mutation(&pm_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert!(ne_data.phase_matched());

    assert_eq!(NeoEpitopeType::Missense, ne_data.variant_type());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();
    ne_data.set_amino_acids(&ref_genome, 3).unwrap();
    assert_eq!(NeoEpitopeType::StopLost, ne_data.variant_type());

    let up_bases = cat!(sub(&chr1_bases, 46, 51), sub(&chr1_bases, 60, 64));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(pm_data.alt, sub(&chr1_bases, 65, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);
}

