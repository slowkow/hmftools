//! NEO's NeoEpitopeUtilsTest, ported

use super::*;
use crate::cat;
use crate::epitope::*;
use hmf_common::neo::{NeoEpitopeFusion, NeoEpitopeType};

#[test]
fn test_dna_rna_routines() {
    // AminoAcidRna.swapDnaToRna / swapRnaToDna are not used by the port; reverseComplementBases is
    let dna_bases = "AGCTTCGACT";
    let reverse_strand_dna = reverse_complement_bases(&dna_bases);
    assert_eq!(reverse_strand_dna, "AGTCGAAGCT");
    assert_eq!(dna_bases, reverse_complement_bases(&reverse_strand_dna));
}

#[test]
fn test_upstream_coding_bases() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(100);

    ref_genome.put(CHR_1, &chr1_bases);

    // tests: repeated for each strand

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(85);

    let trans_data_up = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // upstream: intronic
    let ne_position = 32;  // intronic
    let ne_orientation = ORIENT_REV;

    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 9).unwrap();

    let act_coding_bases = sub(&chr1_bases, 40, 49);
    assert_eq!(act_coding_bases, coding_bases);

    // spanning 3 exons
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 25).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 40, 51), sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 83));
    assert_eq!(act_coding_bases, coding_bases);

    // upstream: intronic, coding bases insufficient for required bases
    let ne_position = 76;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = sub(&chr1_bases, 80, 86);
    assert_eq!(act_coding_bases, coding_bases);

    // upstream: exonic
    let ne_position = 44;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 44, 51), sub(&chr1_bases, 60, 63));
    assert_eq!(act_coding_bases, coding_bases);

    // upstream: exonic, coding bases insufficient for required bases
    let ne_position = 68;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 68, 71), sub(&chr1_bases, 80, 86));
    assert_eq!(act_coding_bases, coding_bases);

    // upstream: exonic, coding finishes in same exon
    let ne_position = 81;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = sub(&chr1_bases, 81, 86);
    assert_eq!(act_coding_bases, coding_bases);

    // test again with reverse orientation
    let ne_orientation = ORIENT_FWD;
    let ne_position = 55;  // intronic

    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 9).unwrap();

    let act_coding_bases = sub(&chr1_bases, 42, 51);
    assert_eq!(act_coding_bases, coding_bases);

    // spanning 3 exons
    let ne_position = 75;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 25).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 51), sub(&chr1_bases, 60, 71));
    assert_eq!(act_coding_bases, coding_bases);

    // intronic, coding bases insufficient for required bases
    let ne_position = 35;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(act_coding_bases, coding_bases);

    // exonic
    let ne_position = 44;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 26, 31), sub(&chr1_bases, 40, 45));
    assert_eq!(act_coding_bases, coding_bases);

    // exonic, coding bases insufficient for required bases
    let ne_position = 40;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 41));
    assert_eq!(act_coding_bases, coding_bases);

    // exonic, coding finishes in same exon
    let ne_position = 29;
    let coding_bases = get_upstream_coding_bases(&ref_genome, trans_data_up, CHR_1, ne_position, ne_orientation, 10).unwrap();

    let act_coding_bases = sub(&chr1_bases, 25, 30);
    assert_eq!(act_coding_bases, coding_bases);
}

#[test]
fn test_downstream_coding_bases() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(280);

    ref_genome.put(CHR_1, &chr1_bases);

    // tests: repeated for each strand
    // intronic within coding section
    // starting in an exon but skipping it
    // starting in an exon and using it
    // upstream starting at first splice acceptor (2nd exon)
    // getting all remaining bases
    // non-coding - get all bases

    let exon_starts = [100, 120, 140, 160, 180, 200, 220, 240, 260];
    let coding_start = Some(125);
    let coding_end = Some(245);

    let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    // intronic
    let ne_position = 132;
    let ne_orientation = ORIENT_REV;

    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 9, true, true, false).unwrap();

    let act_coding_bases = sub(&chr1_bases, 140, 149);
    assert_eq!(act_coding_bases, coding_bases);

    // not allowed to start in the same exon
    let ne_position = 144;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 15, false, true, false).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 160, 171), sub(&chr1_bases, 180, 184));
    assert_eq!(act_coding_bases, coding_bases);

    // allowed to start in an exon
    let ne_position = 147;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 10, true, true, false).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 47, 51), sub(&chr1_bases, 160, 166));
    assert_eq!(act_coding_bases, coding_bases);

    // upstream
    let ne_position = 80;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 15, true, true, false).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 120, 131), sub(&chr1_bases, 140, 144));
    assert_eq!(act_coding_bases, coding_bases);

    // all remaining bases
    let ne_position = 215;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 1, false, true, true).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 220, 231), sub(&chr1_bases, 240, 251), sub(&chr1_bases, 260, 271));
    assert_eq!(act_coding_bases, coding_bases);

    // non-coding transcript
    let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, None, None, false, "");

    let ne_position = 215;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 1, false, true, true).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 220, 231), sub(&chr1_bases, 240, 251), sub(&chr1_bases, 260, 271));
    assert_eq!(act_coding_bases, coding_bases);


    // repeated for negative strand
    let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    let ne_position = 155;
    let ne_orientation = ORIENT_FWD;

    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 11, false, true, false).unwrap();

    let act_coding_bases = sub(&chr1_bases, 140, 151);
    assert_eq!(act_coding_bases, coding_bases);

    // not allowed to start in the same exon
    let ne_position = 165;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 10, false, true, false).unwrap();

    let act_coding_bases = sub(&chr1_bases, 141, 151);
    assert_eq!(act_coding_bases, coding_bases);

    // allowed to start in an exon
    let ne_position = 165;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 10, true, true, false).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 147, 151), sub(&chr1_bases, 160, 166));
    assert_eq!(act_coding_bases, coding_bases);

    // upstream
    let ne_position = 300;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 15, true, true, false).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 227, 231), sub(&chr1_bases, 240, 251));
    assert_eq!(act_coding_bases, coding_bases);

    // all remaining bases
    let ne_position = 135;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 1, false, true, true).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 100, 111), sub(&chr1_bases, 120, 131));
    assert_eq!(act_coding_bases, coding_bases);

    // non-coding transcript
    let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, None, None, false, "");

    let ne_position = 135;
    let coding_bases = get_downstream_coding_bases(&ref_genome, trans_data, CHR_1, ne_position, ne_orientation, 1, false, true, true).unwrap();

    let act_coding_bases = cat!(sub(&chr1_bases, 100, 111), sub(&chr1_bases, 120, 131));
    assert_eq!(act_coding_bases, coding_bases);
}

#[test]
fn test_skipped_acceptor_donors() {

    let exon_starts = [100, 300, 500, 700, 900];

    let trans_data1 = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 100, None, None, false, "");

    let exon_starts = [200, 500, 800, 1100];

    let trans_data2 = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 100, None, None, false, "");

    let pos_boundaries = [50, 1200];

    let trans_data_list = &[trans_data1, trans_data2];

    let skipped_count = find_skipped_exon_boundaries(trans_data_list, pos_boundaries, true, true);
    assert_eq!(6, skipped_count);  // exon start at 500 is repeated so not double-counted

    let pos_boundaries = [600, 1200];

    let skipped_count = find_skipped_exon_boundaries(trans_data_list, pos_boundaries, true, true);
    assert_eq!(4, skipped_count);  // exon start at 500 is repeated so not double-counted

    let pos_boundaries = [1200, 1400];

    let skipped_count = find_skipped_exon_boundaries(trans_data_list, pos_boundaries, true, true);
    assert_eq!(0, skipped_count);  // exon start at 500 is repeated so not double-counted

    // upstream donor skipping
    let pos_boundaries = [350, 1300];

    let skipped_count = find_skipped_exon_boundaries(trans_data_list, pos_boundaries, false, false);
    assert_eq!(6, skipped_count);  // exon start at 500 is repeated so not double-counted

    let pos_boundaries = [650, 750];

    let skipped_count = find_skipped_exon_boundaries(trans_data_list, pos_boundaries, false, false);
    assert_eq!(0, skipped_count);  // exon start at 500 is repeated so not double-counted
}

#[test]
fn test_peptide_combinations() {

    // standard missense type
    let up_aas = "ABC";
    let novel_aas = "D";
    let down_aas = "EFG";

    let peptides = generate_peptides(up_aas, novel_aas, down_aas, [3, 4], 1);
    assert_eq!(7, peptides.len());

    let peptide_strings: Vec<&str> = peptides.iter().map(|x| x.peptide.as_str()).collect();
    assert!(peptide_strings.contains(&"BCD"));
    assert!(peptide_strings.contains(&"CDE"));
    assert!(peptide_strings.contains(&"DEF"));
    assert!(peptide_strings.contains(&"ABCD"));
    assert!(peptide_strings.contains(&"BCDE"));
    assert!(peptide_strings.contains(&"CDEF"));
    assert!(peptide_strings.contains(&"DEFG"));

    // insertion type
    let novel_aas = "DXX";

    let peptides = generate_peptides(up_aas, novel_aas, down_aas, [3, 4], 0);
    let peptide_strings: Vec<&str> = peptides.iter().map(|x| x.peptide.as_str()).collect();
    assert_eq!(11, peptides.len());
    assert!(peptide_strings.contains(&"BCD"));
    assert!(peptide_strings.contains(&"CDX"));
    assert!(peptide_strings.contains(&"DXX"));
    assert!(peptide_strings.contains(&"XXE"));
    assert!(peptide_strings.contains(&"XEF"));
    assert!(peptide_strings.contains(&"ABCD"));
    assert!(peptide_strings.contains(&"BCDX"));
    assert!(peptide_strings.contains(&"CDXX"));
    assert!(peptide_strings.contains(&"DXXE"));
    assert!(peptide_strings.contains(&"XXEF"));
    assert!(peptide_strings.contains(&"XEFG"));

    // no novel - like a fusion
    let novel_aas = "";
    let peptides = generate_peptides(up_aas, novel_aas, down_aas, [3, 4], 0);
    let peptide_strings: Vec<&str> = peptides.iter().map(|x| x.peptide.as_str()).collect();

    assert_eq!(5, peptides.len());
    assert!(peptide_strings.contains(&"BCE"));
    assert!(peptide_strings.contains(&"CEF"));
    assert!(peptide_strings.contains(&"ABCE"));
    assert!(peptide_strings.contains(&"BCEF"));
    assert!(peptide_strings.contains(&"CEFG"));

    // all downstream bases novel
    let novel_aas = "VWXYZ";
    let down_aas = "";
    let peptides = generate_peptides(up_aas, novel_aas, down_aas, [3, 4], 0);
    let peptide_strings: Vec<&str> = peptides.iter().map(|x| x.peptide.as_str()).collect();

    assert_eq!(10, peptides.len());
    assert!(peptide_strings.contains(&"BCV"));
    assert!(peptide_strings.contains(&"CVW"));
    assert!(peptide_strings.contains(&"VWX"));
    assert!(peptide_strings.contains(&"WXY"));
    assert!(peptide_strings.contains(&"XYZ"));

    assert!(peptide_strings.contains(&"ABCV"));
    assert!(peptide_strings.contains(&"BCVW"));
    assert!(peptide_strings.contains(&"CVWX"));
    assert!(peptide_strings.contains(&"VWXY"));
    assert!(peptide_strings.contains(&"WXYZ"));

    // test flanks specifically
    let up_aas = "ABCDEFG";
    let novel_aas = "";
    let down_aas = "TUVWXYZ";
    let peptides = generate_peptides(up_aas, novel_aas, down_aas, [3, 3], 5);
    assert_eq!(2, peptides.len());
    assert!(peptides.iter().any(|x| x.peptide == "FGT" && x.up_flank == "ABCDE" && x.down_flank == "UVWXY"));
    assert!(peptides.iter().any(|x| x.peptide == "GTU" && x.up_flank == "BCDEF" && x.down_flank == "VWXYZ"));
}

