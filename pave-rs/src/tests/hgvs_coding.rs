//! HgvsCodingTest

use super::*;
use crate::impact::effects::{CodingType, RegionType, VariantEffect};
use crate::impact::{SE_END, SE_START};

use crate::impact::{hgvs_coding, CodingContext};
use crate::impact::utils::{create_right_aligned_variant, find_variant_impacts};
use crate::variant::NO_LOCAL_PHASE_SET;
use super::homology_realign;

    #[test]
    fn test_point_mutations()
    {
        let mut pos = 100;
        let mut ref_ = "A".to_string();
        let mut alt = "C".to_string();
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut cc = CodingContext::default();
        cc.strand = POS_STRAND;

        cc.region_type = RegionType::Intronic;
        cc.coding_type = CodingType::Utr5p;
        cc.coding_base = 10;
        cc.nearest_exon_distance = -50;

        assert_eq!("c.-10-50A>C", hgvs_coding::generate(&var, &cc));

        cc.region_type = RegionType::Intronic;
        cc.coding_type = CodingType::Utr3p;
        cc.coding_base = 20;
        cc.nearest_exon_distance = 10;

        assert_eq!("c.*20+10A>C", hgvs_coding::generate(&var, &cc));

        cc.region_type = RegionType::Intronic;
        cc.coding_type = CodingType::Coding;
        cc.coding_base = 10;
        cc.nearest_exon_distance = -50;

        assert_eq!("c.10-50A>C", hgvs_coding::generate(&var, &cc));

        cc.nearest_exon_distance = 50;

        assert_eq!("c.10+50A>C", hgvs_coding::generate(&var, &cc));

        cc.region_type = RegionType::Exonic;
        assert_eq!("c.10A>C", hgvs_coding::generate(&var, &cc));
    }

    #[test]
    fn test_coding_duplication()
    {
        // pos codons: M 20-22, A 23-25, C 26-28, D 29-31, L 32-34, L 35-37, G 38-40, H 41-43, E 44-46, stop_x 47-49
        // amino acids:           M  A  C  D  L  L  G  H  E  X
        // exon 1:                20        30
        // position:              012345678901234
        let mut ref_coding_bases1 = "ATGGCTTGTGACTTA".to_string();

        // intron:                   40
        // position             567890123456789
        // bases into intron:   12345   7654321
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string();

        // exon 2:                50        60
        // position:              012345678901234
        let mut ref_coding_bases2 = "TTAGGACACGAGTAA".to_string();

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, intronic_bases, ref_coding_bases2, generate_test_bases(20));
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 50], 24, Some(20), Some(64), false, "");

        // duplication of a codon
        let mut pos = 25;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = sub(&ref_bases, pos, pos + 4);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        let mut impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.7_9dupTGT", impact.coding.hgvs);
        assert_eq!("p.Cys3dup", impact.pc().hgvs);

        // duplication of 1 intronic base
        pos = 40;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 2);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::Intronic, impact.top_effect().unwrap());
        assert_eq!(6, impact.coding.nearest_exon_distance);
        assert_eq!("c.15+7dupG", impact.coding.hgvs);

        // duplication of 2 intronic bases
        pos = 40;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 3);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::Intronic, impact.top_effect().unwrap());
        assert_eq!(6, impact.coding.nearest_exon_distance);
        assert_eq!("c.15+7_15+8dupGG", impact.coding.hgvs);

        // test nearer to the next exon
        pos = 47;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 3);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::SpliceAcceptor, impact.top_effect().unwrap());

        assert_eq!(16, impact.coding.coding_base);
        assert_eq!(-3, impact.coding.nearest_exon_distance);

        // the base positions 48-49 are duplicated, so these are intronic bases -1 and -2
        assert_eq!("c.16-2_16-1dupAA", impact.coding.hgvs);

        // also for a single duplicated base
        alt = sub(&ref_bases, pos, pos + 2);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.16-2dupA", impact.coding.hgvs);

        // again on the negative strand
        let neg_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[10, 50], 24, Some(20), Some(64), false, "");

        // coding goes from 64 -> 50 then 34 -> 20

        // first an exonic dup
        pos = 58;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 4);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, neg_trans, &ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.4_6dupCTC", impact.coding.hgvs);

        // dup in intron following an exon
        pos = 43;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 4);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, neg_trans, &ref_genome);
        assert_eq!(VariantEffect::Intronic, impact.top_effect().unwrap());
        assert_eq!(15, impact.coding.coding_base);
        assert_eq!(6, impact.coding.nearest_exon_distance);

        // bases at pos 44-46 are duplicated, these are intronic bases 6-4
        assert_eq!("c.15+4_15+6dupGGG", impact.coding.hgvs);

        // single-base dup in prior intron
        pos = 35;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = cat!(sub(&ref_bases, pos, pos + 1), sub(&ref_bases, pos + 1, pos + 2)); // duplicating the next base
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, neg_trans, &ref_genome);

        assert_eq!(16, impact.coding.coding_base);
        assert_eq!(-2, impact.coding.nearest_exon_distance); // nearest exon ends at 34, using upstream pos of 37

        assert_eq!("c.16-2dupA", impact.coding.hgvs);

        pos = 36;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = cat!(ref_, sub(&ref_bases, pos + 1, pos + 2), sub(&ref_bases, pos + 1, pos + 2));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, neg_trans, &ref_genome);

        assert_eq!(16, impact.coding.coding_base);
        assert_eq!(-3, impact.coding.nearest_exon_distance);

        assert_eq!("c.16-4_16-3dupAA", impact.coding.hgvs);
    }

    #[test]
    fn test_realigned_coding_duplication1()
    {
        // pos codons: M 20-22 (1-3) A 23-25 (4-6) E 26-28 (7-9) N 29-31 (10-12) L 32-34 (13-15) G 35-37 (16-18)
        // amino acids:           M  A  E  N  L  G
        // exon 1:                20        30
        // position:              012345678901234567
        // coding position:       123456789012345678
        let mut ref_coding_bases1 = "ATGGCTGAAAACTTAGGC".to_string();

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, generate_test_bases(30));
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 50], 24, Some(20), Some(64), false, "");

        let gene_data_cache = homology_realign::gene_cache_with(homology_realign::create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, pos_trans.strand, pos_trans.trans_start, pos_trans.trans_end), pos_trans);

        // test 1: realignment for homopolymer expansion
        let mut pos = 26;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, "A");
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        var.set_variant_details(NO_LOCAL_PHASE_SET, "A", "A", 4);
        { let ra = create_right_aligned_variant(&var, &ref_genome).unwrap(); var.set_realigned_variant(ra); }

        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();

        assert!(var.gene_impacts(GENE_NAME_1).is_some());
        let mut impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();

        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());
        assert_eq!("c.11dupA", impact.coding.hgvs);
    }

    #[test]
    fn test_realigned_coding_duplication2()
    {
        // pos codons: M 20-22 (1-3) A 23-25 (4-6) E 26-28 (7-9) E 29-31 (10-12) E 32-34 (13-15) N 35-37 (16-18)
        // amino acids:           M  A  E  E  E  N  L  G
        // exon 1:                20        30
        // position:              012345678901234567
        // coding position:       123456789012345678
        let mut ref_coding_bases1 = "ATGGCTGAAGAAGAACTT".to_string();

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, generate_test_bases(30));
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 50], 34, Some(20), Some(64), false, "");

        let gene_data_cache = homology_realign::gene_cache_with(homology_realign::create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, pos_trans.strand, pos_trans.trans_start, pos_trans.trans_end), pos_trans);

        // test 1: realignment for homopolymer expansion
        let mut pos = 25;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, "GAA");
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        var.set_variant_details(NO_LOCAL_PHASE_SET, "GAA", "GAA", 3);
        { let ra = create_right_aligned_variant(&var, &ref_genome).unwrap(); var.set_realigned_variant(ra); }

        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();

        assert!(var.gene_impacts(GENE_NAME_1).is_some());
        let mut impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();

        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());
        assert_eq!("c.13_15dupGAA", impact.coding.hgvs);
    }


    #[test]
    fn test_amino_acid_right_alignment()
    {
        // while there is no homology, the amino acid impact is right aligned - only relevant for positive-strand genes

        // amino acids:           M   A   P   P   P   P   Q
        // AA index:              1   2   3   4   5   6   7
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        // coding index:          123 456 789 012 345 678 901
        let mut ref_coding_bases1 = "ATG GCT CCA CCC CCG CCT CAG".replace(" ", "");
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string(); // as in previous test

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, intronic_bases);
        let mut ref_genome = MockRefGenome::default();
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 56], 30, Some(20), Some(70), false, "");

        let gene_data_cache = homology_realign::gene_cache_with(homology_realign::create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, pos_trans.strand, pos_trans.trans_start, pos_trans.trans_end), pos_trans);

        // duplication of first P should be converted to a right-aligned AA change
        let mut pos = 25;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = sub(&ref_bases, pos, pos + 4);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();

        assert!(var.gene_impacts(GENE_NAME_1).is_some());
        let mut impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.7_9dupCCA", impact.coding.hgvs);
        assert_eq!("p.Pro6dup", impact.pc().hgvs);
    }

    #[test]
    fn test_amino_acid_right_alignment2()
    {
        // amino acids:         M   Q   L   L   R   G
        // AA index:            1   2   3   4   5   6
        // exon 1:              20           30
        // position:            012 345 678 901
        // coding index:        123 456 789 012 345 678
        let mut ref_coding_bases = "ATG CAA TTA TTA AGG GGA CAC".replace(" ", "");
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string();

        // variant is inserting RV = AGA (R) and GTG (V)
        // ref_: A         AGG
        // alt: A AGA GTG AGG

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases, intronic_bases);
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());

        let gene_data = homology_realign::create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 20, 90);

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 56], 30, Some(20), Some(70), false, "");

        let gene_data_cache = homology_realign::gene_cache_with(gene_data.clone(), pos_trans);

        

        // duplication of first P should be converted to a right-aligned AA change
        let mut pos = 31;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, "AGAGTG");
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();

        assert!(var.gene_impacts(GENE_NAME_1).is_some());
        let mut impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        // assertEquals("c.12insAGAGTG", impact.coding.hgvs);
        assert_eq!("p.Arg5_Gly6insValArg", impact.pc().hgvs);
    }

    #[test]
    fn test_amino_acid_right_alignment3()
    {
        // amino acids:           M   G   L   *
        // AA index:              1   2   3   4
        // exon 1:                20           30
        // position:              012 345 678 901
        // coding index:          123 456 789 012
        let mut ref_coding_bases =  "ATG GGA CTA TGA".replace(" ", "");
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string();

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases, intronic_bases);
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());

        let gene_data = homology_realign::create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, 20, 90);

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 56], 30, Some(20), Some(70), false, "");

        let gene_data_cache = homology_realign::gene_cache_with(gene_data.clone(), pos_trans);

        

        // duplication of first P should be converted to a right-aligned AA change
        let mut pos = 25;
        let mut ref_ = sub(&ref_bases, pos, pos + 2);
        let mut alt = sub(&ref_bases, pos, pos + 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();

        assert!(var.gene_impacts(GENE_NAME_1).is_some());
        let mut impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        assert_eq!("c.7delC", impact.coding.hgvs);
        assert_eq!("p.Leu3fs", impact.pc().hgvs);
    }

    #[ignore]
    #[test]
    fn test_right_alignment_of_coding_of_in_frame_duplication()
    {
        // pos codons: M 20-22, A 23-25, D 26-28, A 29-31, D 32-34, S 35-37, Q 38-40, L 41-43, G 44-46, H 47-49, E 50-52, stop_x 47-49
        // amino acids:           M   A   D   A   D   S   Q
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        let mut ref_coding_bases1 = "ATG GCT GAT GCT GAT TCG CAG".replace(" ", "");

        // intron:              41     50
        // position             123456789012345
        // bases into intron:   12345   7654321
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string();

        // amino acids:           L   G   H   E   X
        // exon 2:                56   60           70
        // position:              678 901 234 567 890
        let mut ref_coding_bases2 = "TTA GGA CAC GAG TAA".replace(" ", "");

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, intronic_bases, ref_coding_bases2, generate_test_bases(20));
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 56], 30, Some(20), Some(70), false, "");

        let gene_data_cache = homology_realign::gene_cache_with(homology_realign::create_ensembl_gene_data(GENE_ID_1, GENE_NAME_1, CHR_1, pos_trans.strand, pos_trans.trans_start, pos_trans.trans_end), pos_trans);

        // duplication of first AD in ADAD should be reported right-maximal
        // ATG [GCT GAT] GCT GAT TCG CAG -> ATG [GCT GAT GCT GAT] GCT GAT TCG CAG = ATG GCT GAT [GCT GAT GCT GAT] TCG CAG
        let mut pos = 23;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = sub(&ref_bases, pos, pos + 7);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        { let ra = create_right_aligned_variant(&var, &ref_genome).unwrap(); var.set_realigned_variant(ra); }

        // TODO: check if use of RA var passes
        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();

        let mut impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();
        // VariantTransImpact impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.9_14dupGCTGAT", impact.coding.hgvs);
        assert_eq!("p.Ala4_Asp5dup", impact.pc().hgvs);

        // duplication of DA in ADAD should give the same result
        // ATG GCT [GAT GCT] GAT TCG CAG -> ATG GCT [GAT GCT GAT GCT] GAT TCG CAG = ATG GCT GAT [GCT GAT GCT GAT] TCG CAG
        pos = 23 + 3;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 7);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        find_variant_impacts(&mut var, gene_data_cache, None, &ref_genome).unwrap();
        impact = var.gene_impacts(GENE_NAME_1).unwrap()[0].clone();

        // impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.9_14dupGCTGAT", impact.coding.hgvs);
        assert_eq!("p.Ala4_Asp5dup", impact.pc().hgvs);

        // duplication of second AD in ADAD should give the same result
        // ATG GCT GAT [GCT GAT] TCG CAG -> ATG GCT GAT [GCT GAT GCT GAT] TCG CAG
        pos = 23 + 6;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = sub(&ref_bases, pos, pos + 7);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.9_14dupGCTGAT", impact.coding.hgvs);
        assert_eq!("p.Ala4_Asp5dup", impact.pc().hgvs);
    }

    #[test]
    fn test_right_alignment_of_coding_of_frameshift_duplication()
    {
        // pos codons: M 20-22, A 23-25, D 26-28, M 29-31, A 32-34, S 35-37, Q 38-40, L 41-43, G 44-46, H 47-49, E 50-52, stop_x 47-49
        // amino acids:           M   A   D   M   A   S   Q
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        // coding base:           123 456 789 012 345 678 901
        let mut ref_coding_bases1 = "ATG GCT GAT ATG GCC TCG CAG".replace(" ", "");
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string(); // same as previous test
        let mut ref_coding_bases2 = "TTAGGACACGAGTAA".to_string(); // ditto

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, intronic_bases, ref_coding_bases2, generate_test_bases(20));
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 56], 30, Some(20), Some(70), false, "");

        let mut pos = 29;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = sub(&ref_bases, pos, pos + 6);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        let mut impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        // ATG GCT GAT A[TG GCC] TCG CAG -> ATG GCT GAT A[TG GCC TG GCC] TCG CAG = ATG GCT GAT AT[G GCC TG GCC T]CG CAG
        assert_eq!("c.11_15dupTGGCC", impact.coding.hgvs);
        // new AA seq is M A D M A W P
        assert_eq!("p.Ser6fs", impact.pc().hgvs);
    }

    #[test]
    fn test_right_alignment_of_coding_of_frameshift_deletion()
    {
        // amino acids:           M   K   E   K   K   K   P
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        let mut ref_coding_bases1 = "ATG AAG GAA AAA AAA AAG CCT".replace(" ", "");
        let mut intronic_bases = "TTTTGGGGCCCCAAA".to_string(); // same as previous test
        let mut ref_coding_bases2 = "TTAGGACACGAGTAA".to_string(); // ditto

        let mut ref_genome = MockRefGenome::default();

        let mut ref_bases = cat!(generate_test_bases(20), ref_coding_bases1, intronic_bases, ref_coding_bases2, generate_test_bases(20));
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 56], 30, Some(20), Some(70), false, "");

        let mut pos = 34;
        let mut ref_ = sub(&ref_bases, pos, pos + 3);
        let mut alt = sub(&ref_bases, pos, pos + 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);

        let mut impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        assert_eq!("c.16_17delAA", impact.coding.hgvs);
        assert_eq!("p.Lys6fs", impact.pc().hgvs);
    }

    #[test]
    fn test_dup_in5_p_utr()
    {
        // random(20) exon1 intron1 exon2 intron2 exon3(coding) intron3 exon4(coding) random(20)
        // all exons and intronns of length 12
        let mut exon1 = "CTAGGACACGAG".to_string();
        let mut intron1 = "CGAGGGCCCAAA".to_string();
        let mut exon2 = "GGACACGAGTAA".to_string();
        let mut intron2 = "GGGCCCAAATTT".to_string();
        let mut exon3 = "ATGAAGGAACCT".to_string();
        let mut intron3 = "AAAGGGCCCTTT".to_string();
        let mut exon4 = "AAGATGGAACCT".to_string();
        let mut ref_bases = cat!(generate_test_bases(20), exon1, intron1, exon2, intron2, exon3, intron3, exon4, generate_test_bases(20));

        let mut ref_genome = MockRefGenome::default();
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        

        let pos_trans = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[20, 44, 68, 92], 11, Some(20), Some(104), false, "");

        let mut var = create_dup_variant(&ref_bases, 20, 3);
        let mut impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.2_4dupTAG", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 31, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.12+1_12+4dupCGAG", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 33, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.12+3_12+6dupAGGG", impact.coding.hgvs);

        //CGAGGGCCCAAA
        var = create_dup_variant(&ref_bases, 34, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.12+4_12+7dupGGGC", impact.coding.hgvs); // should this be 12+4_13-5?

        var = create_dup_variant(&ref_bases, 36, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.12+6_12+9dupGCCC", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 37, 4);
        impact = classify(&var, pos_trans, &ref_genome);
//        assertEquals("c.13-5_13-2dupCCAA", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 38, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.13-5_13-2dupCCAA", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 39, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.13-4_13-1dupCAAA", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 40, 4);
        impact = classify(&var, pos_trans, &ref_genome);
        assert_eq!("c.13-3_13-0dupAAAG", impact.coding.hgvs);

        var = create_dup_variant(&ref_bases, 42, 2);
        impact = classify(&var, pos_trans, &ref_genome);
//                assertEquals("c.13-6_11-3dupCCAA", impact.coding.hgvs);
    }

fn create_dup_variant(ref_bases: &str, pos: i32, length: i32) -> VariantData {
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = sub(&ref_bases, pos, pos + length + 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);
        let mut alt_bases = alt[1..].to_string();
        var.set_variant_details(NO_LOCAL_PHASE_SET, &alt_bases, &alt_bases, 1);
        var
}
