//! ProteinImpactTest

use super::*;
use crate::impact::effects::VariantEffect;

struct Fx {
    ref_genome: MockRefGenome,
    ref_bases: String,
    pos_trans: &'static TranscriptData,
    neg_trans: &'static TranscriptData,
}

impl Fx {
    /// an exon with specific amino acids on each strand
    fn new() -> Fx {
        let pre_gene = 10;
        let pre_post_coding = 10;
        let mut ref_bases = generate_test_bases(pre_gene);
        // pos codons: M 20-22, A 23-25, C 26-28, D 29-31, L 32-34, L 35-37, G 38-40, H 41-43, E 44-46, stop X 47-49
        ref_bases += &generate_test_bases(pre_post_coding);
        let mut coding_bases = get_amino_acid_codon(b'M', 0);
        coding_bases += &get_amino_acids_codons("ACDLLGHE", false);
        coding_bases += "TAA";
        ref_bases += &coding_bases;
        ref_bases += &generate_test_bases(pre_gene);

        let mut trans_start = pre_gene;
        let mut coding_start = trans_start + pre_post_coding;
        let mut coding_end = coding_start + coding_bases.len() as i32 - 1;
        let mut trans_end = coding_end + pre_post_coding;
        let pos_trans =
            create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[trans_start], trans_end - trans_start, Some(coding_start), Some(coding_end), false, "");

        let next_rand = 100 - trans_end;
        ref_bases += &generate_test_bases(next_rand - 1);
        ref_bases += &generate_test_bases(pre_gene);
        ref_bases += &rc(&coding_bases);
        ref_bases += &generate_test_bases(pre_gene);

        trans_start = trans_end + next_rand;
        coding_start = trans_start + pre_post_coding;
        coding_end = coding_start + coding_bases.len() as i32 - 1;
        trans_end = coding_end + pre_post_coding;
        // neg codons: M 139-37, A 136-134, C 133-131, D 130-128, L 127-125, L 124-122, G 121-119, H 118-116, E 115-113, stop X 112-110
        let neg_trans =
            create_trans_exons(GENE_ID_2, TRANS_ID_2, NEG_STRAND, &[trans_start], trans_end - trans_start, Some(coding_start), Some(coding_end), false, "");

        let mut ref_genome = MockRefGenome::default();
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());
        Fx { ref_genome, ref_bases, pos_trans, neg_trans }
    }

    /// checkHgvsStrings: the positive strand, then the same change on the negative strand
    fn check_hgvs_strings(&self, pos: i32, ref_len: i32, alt: &str, effect: VariantEffect, coding_str: &str, protein_str: &str) {
        let r = sub(&self.ref_bases, pos, pos + ref_len);
        let var = VariantData::new(CHR_1, pos, &r, alt);
        let impact = classify(&var, self.pos_trans, &self.ref_genome);
        assert_eq!(Some(effect), impact.top_effect());
        assert_eq!(coding_str, impact.coding.hgvs);
        assert_eq!(protein_str, impact.pc().hgvs);

        let neg_pos_equiv = self.neg_trans.coding_end.unwrap() - (pos - self.pos_trans.coding_start.unwrap());
        let neg_pos = if ref_len == alt.len() as i32 { neg_pos_equiv - (ref_len - 1) } else { neg_pos_equiv - ref_len };
        let neg_ref = sub(&self.ref_bases, neg_pos, neg_pos + ref_len);
        let neg_alt = if var.is_insert() {
            cat!(neg_ref, rc(&alt[1..]))
        } else if var.is_deletion() {
            sub(&neg_ref, 0, 1)
        } else {
            rc(alt)
        };
        let var_neg = VariantData::new(CHR_1, neg_pos, &neg_ref, &neg_alt);
        let impact_neg = classify(&var_neg, self.neg_trans, &self.ref_genome);
        assert_eq!(Some(effect), impact_neg.top_effect());
        assert_eq!(coding_str, impact_neg.coding.hgvs);
        assert_eq!(protein_str, impact_neg.pc().hgvs);
    }
}
use crate::impact::{SE_END, SE_START};

    #[test]
    fn test_synonymous_missense()
    {
        let mut fx = Fx::new();
        // SNV - synonymous on first A, 3rd codon base
        let mut pos = 25;
        let mut ref_ = sub(&fx.ref_bases, pos, pos + 1);
        let mut alt = "C".to_string();
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::Synonymous, impact.top_effect().unwrap());

        assert_eq!("c.6T>C", impact.coding.hgvs);
        assert_eq!("p.Ala2=", impact.pc().hgvs);

        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::Synonymous, "c.6T>C", "p.Ala2=");

        // same codon but missense
        pos = 24;
        alt = "G".to_string();
        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::Missense, "c.5C>G", "p.Ala2Gly");

        // MNV same codon synonymous: L=TTA -> CTG
        pos = 32;
        alt = "CTG".to_string();

        fx.check_hgvs_strings(pos, 3, &alt, VariantEffect::Synonymous, "c.13_15delTTAinsCTG", "p.Leu5=");

        // now missense
        alt = "GGC".to_string();
        fx.check_hgvs_strings(pos, 3, &alt, VariantEffect::Missense, "c.13_15delTTAinsGGC", "p.Leu5Gly");

        // MNV synonymous spanning 2 codons: L L, TTA TTA -> TTG -> CTG
        pos = 34;
        alt = "GC".to_string();
        fx.check_hgvs_strings(pos, 2, &alt, VariantEffect::Synonymous, "c.15_16delATinsGC", "p.Leu5_Leu6=");

        // now missense: L L, TTA TTA -> TCG -> GTA
        pos = 33;
        alt = "CGG".to_string();
        fx.check_hgvs_strings(pos, 3, &alt, VariantEffect::Missense, "c.14_16delTATinsCGG", "p.Leu5_Leu6delinsSerVal");

        // causing stop gained, across 2 codons, L(TTA) L(TTA) -> (C)TGT (X)TAA
        // eg p.Cys495_Val496delinsArg*
        pos = 33;
        alt = "GTTA".to_string();
        fx.check_hgvs_strings(pos, 4, &alt, VariantEffect::StopGained, "c.14_17delTATTinsGTTA", "p.Leu5_Leu6delinsCys*");

        // start lost
        pos = 21;
        alt = "GT".to_string();
        fx.check_hgvs_strings(pos, 2, &alt, VariantEffect::StartLost, "c.2_3delTGinsGT", "p.Met1?");

        // stop lost
        pos = 48;
        alt = "C".to_string();
        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::StopLost, "c.29A>C", "p.Ter10Serext*?");

        // stop gained but in the first codon where an MNV is affecting 2
        pos = 27;
        alt = "AAC".to_string();
        fx.check_hgvs_strings(pos, 3, &alt, VariantEffect::StopGained, "c.8_10delGTGinsAAC", "p.Cys3_Asp4delins*");
    }

    #[test]
    fn test_inframe_deletion()
    {
        let mut fx = Fx::new();
        // p.Lys2del
        // conservative DEL: codon 2 = A
        let mut pos = 22;
        let mut ref_ = sub(&fx.ref_bases, pos, pos + 4);
        let mut alt = sub(&ref_, 0, 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());

        assert_eq!("c.4_6delGCT", impact.coding.hgvs);
        assert_eq!("p.Ala2del", impact.pc().hgvs);

        fx.check_hgvs_strings(pos, 4, &alt, VariantEffect::InframeDeletion, "c.4_6delGCT", "p.Ala2del");

        // 3 codons conservative
        fx.check_hgvs_strings(pos, 10, &alt, VariantEffect::InframeDeletion, "c.4_12delGCTTGTGAC", "p.Ala2_Asp4del");

        // 2 codons non-conservative: deletes part of A and C and makes a G
        pos = 23;
        alt = sub(&fx.ref_bases, pos, pos + 1);

        fx.check_hgvs_strings(pos, 4, &alt, VariantEffect::InframeDeletion, "c.5_7delCTT", "p.Ala2_Cys3delinsGly");

        // causing a stop gained: L(TTA) G(GGA) -> TGA
        pos = 35;
        alt = sub(&fx.ref_bases, pos, pos + 1);
        fx.check_hgvs_strings(pos, 4, &alt, VariantEffect::StopGained, "c.17_19delTAG", "p.Leu6_Gly7delins*");

        // causing a stop lost H(CAC) E(GAG) X(TAA) but report as frameshift since the first AA to change is not a stop
        pos = 41;
        alt = sub(&fx.ref_bases, pos, pos + 1);
        fx.check_hgvs_strings(pos, 7, &alt, VariantEffect::StopLost, "c.23_28delACGAGT", "p.His8fs");

        // causing a stop lost E(GAG) X(TAA) but the inframe DEL makes a synonymous E
        pos = 44;
        alt = sub(&fx.ref_bases, pos, pos + 1);
        fx.check_hgvs_strings(pos, 4, &alt, VariantEffect::StopLost, "c.26_28delAGT", "p.Glu9ext*?");
    }

    #[test]
    fn test_inframe_insertion()
    {
        let mut fx = Fx::new();
        // duplication (single AA) p.Gln8dup
        // duplication (range) p.Gly4_Gln6dup
        // insertion p.Lys2_Leu3insGlnSer
        // insertion (conservative stop) p.Ser81_Val82ins*
        // insertion (non conservative) p.Cys28delinsTrpVal

        // insert a codon in between A1 and C2
        let mut pos = 25;
        let mut ref_ = sub(&fx.ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, get_amino_acid_codon(b'R', 0));
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.6_7insAGG", impact.coding.hgvs);
        assert_eq!("p.Ala2_Cys3insArg", impact.pc().hgvs);

        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::InframeInsertion, "c.6_7insAGG", "p.Ala2_Cys3insArg");

        // insert multiple codons
        alt = cat!(ref_, get_amino_acids_codons("RVS", false));
        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::InframeInsertion, "c.6_7insAGGGTGTCA", "p.Ala2_Cys3insArgValSer");

        // non-conservative insert: CDL -> CEGCL   C + D(GAC) + L(TTA) + L -> C + E(GAA) (G)GGG (C)TGT + L
        pos = 30;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, get_amino_acids_codons("RV", false));
        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::InframeInsertion, "c.11_12insAGGGTG", "p.Asp4delinsGluGlyCys");

        // pos codons: M 20-22, A 23-25, C 26-28, D 29-31, L 32-34, L 35-37, G 38-40, H 41-43, E 44-46, stop_x 47-49
        // M  A  C  D  L  L  G  H  E  X
        // 1  2  3  4  5  6  7  8  9  10
        // ATGGCTTGTGACTTATTAGGACACGAGTAA

        // DLLG -> DCCLLG, or shortened DL -> DCCL, so should be p.D4_L5insCC
        // ref_: L
        //      T       TA
        // alt: C     C   L
        //      T  GT TGT T  TA
        pos = 32;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "GTTGTT"); // get_amino_acids_codons("CC", false)
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());
        assert_eq!("p.Asp4_Leu5insCysCys", impact.pc().hgvs);

        // neg codons: M 139-37, A 136-134, C 133-131, D 130-128, L 127-125, L 124-122, G 121-119, H 118-116, E 115-113, stop_x 112-110
        //   X  E  H  G  L  L  D  C  A  M
        // TTACTCGTGTCCTAATAAGTCACAAGCCAT GATCGATCGA
        // 110                            140

        pos = 123;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, get_amino_acids_codons("GC", true));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.neg_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());
        // assertEquals("p.Leu5_Leu6insCysGly", impact.pc().hgvs);

        // duplications

        // single G repeated
        pos = 37;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "GGA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.18_19insGGA", impact.coding.hgvs);
        assert_eq!("p.Gly7dup", impact.pc().hgvs);

        // repeated more than 1 is no longer a dup
        pos = 37;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "GGAGGAGGA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());
        assert_eq!("p.Gly7_His8insGlyGlyGly", impact.pc().hgvs);

        // range of AAs: L(TTA) G(GGA) H(CAC) E(GAG), T>TAGGACACGA - duplicating GHE
        pos = 36;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AGGACACGA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("c.17_18insAGGACACGA", impact.coding.hgvs);
        assert_eq!("p.Gly7_Glu9dup", impact.pc().hgvs);

        // both Ls are duplicated but treat this is a range rather than repeat of a single
        pos = 31;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "TTATTA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeInsertion, impact.top_effect().unwrap());

        assert_eq!("p.Leu5_Leu6dup", impact.pc().hgvs);
    }

    #[test]
    fn test_frameshifts()
    {
        let mut fx = Fx::new();
        let mut pos = 26;
        let mut ref_ = sub(&fx.ref_bases, pos, pos + 2);
        let mut alt = sub(&ref_, 0, 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        assert_eq!("c.8delG", impact.coding.hgvs);
        assert_eq!("p.Cys3fs", impact.pc().hgvs);

        fx.check_hgvs_strings(pos, 2, &alt, VariantEffect::Frameshift, "c.8delG", "p.Cys3fs");

        // longer DEL still uses first ref_ codon
        fx.check_hgvs_strings(pos, 12, &alt, VariantEffect::Frameshift, "c.8_18delGTGACTTATTA", "p.Cys3fs");

        // from an insert
        alt = cat!(sub(&fx.ref_bases, pos, pos + 1), "CC");
        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::Frameshift, "c.7_8insCC", "p.Cys3fs");

        // causing a stop gained: L(TTA) -> TGA, and no stop-gained identifier even if one is added
        pos = 35;
        ref_ = sub(&fx.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "GA");
        fx.check_hgvs_strings(pos, 1, &alt, VariantEffect::Frameshift, "c.16_17insGA", "p.Leu6*");

        // SNV and deletion combined: ie CT -> A, should be c.5_6delCTinsA instead of c.6delT
        pos = 24;
        ref_ = sub(&fx.ref_bases, pos, pos + 2);
        alt = "A".to_string();
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, fx.pos_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        assert_eq!(5, impact.coding.coding_base);
        assert_eq!(1, impact.coding.deleted_coding_bases);
        assert_eq!("c.5_6delCTinsA", impact.coding.hgvs);
        assert_eq!("p.Ala2fs", impact.pc().hgvs);
    }

    #[test]
    fn test_specific_scenarios()
    {
        let mut fx = Fx::new();
        let mut pre_gene = 10;
        let mut pre_post_coding = 10;
        let mut ref_bases = generate_test_bases(pre_gene);

        // codons: M 20-22, G 23-25, P 26-28, P 29-31, H 32-34, L 35-37, stop_x
        // M  G  P  H  L  X
        // ATGGGACCCCACTTATAA
        // 20

        ref_bases += &generate_test_bases(pre_post_coding);
        let mut coding_bases = get_amino_acid_codon(b'M', 0);
        coding_bases += &"GGACCCCCCCACTTA"; // G, P, P, H, L
        coding_bases += "TAA";

        ref_bases += &coding_bases;
        ref_bases += &generate_test_bases(pre_gene);

        fx.ref_genome.map.insert(CHR_2.to_string(), ref_bases.clone());

        let mut trans_start = pre_gene;
        let mut coding_start = trans_start + pre_post_coding;
        let mut coding_end = coding_start + coding_bases.len() as i32 - 1;
        let mut trans_end = coding_end + pre_post_coding;

        let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[trans_start], trans_end - trans_start, Some(coding_start), Some(coding_end), false, "");

        // frameshift DEL
        let mut pos = 25;
        let mut ref_ = sub(&ref_bases, pos, pos + 2);
        let mut alt = sub(&ref_, 0, 1);
        let mut var = VariantData::new(CHR_2, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data, &fx.ref_genome);
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        assert_eq!("c.7delC", impact.coding.hgvs);
        assert_eq!("p.His5fs", impact.pc().hgvs);

        // conservative inframe deletion with homology - deletes a P, but can be shifted further
        pos = 25;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_2, pos, &ref_, &alt);

        impact = classify(&var, trans_data, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());

        assert_eq!("c.7_9delCCC", impact.coding.hgvs);
        assert_eq!("p.Pro4del", impact.pc().hgvs);
    }

    #[test]
    fn test_specific_scenarios2()
    {
        let mut fx = Fx::new();
        let mut pre_gene = 10;
        let mut pre_post_coding = 10;
        let mut ref_bases = generate_test_bases(pre_gene);

        // codons: X 20-22, K 23-25, K 26-28, K 29-31, M 32-34
        // ATT TTTTTTTGC CAT
        // 20

        ref_bases += &generate_test_bases(pre_post_coding);
        let mut coding_bases = "ATTTTTTTTTGCCAT".to_string(); // M, A, K, K, X

        ref_bases += &coding_bases;
        ref_bases += &generate_test_bases(pre_gene);

        fx.ref_genome.map.insert(CHR_2.to_string(), ref_bases.clone());

        let mut trans_start = pre_gene;
        let mut coding_start = trans_start + pre_post_coding;
        let mut coding_end = coding_start + coding_bases.len() as i32 - 1;
        let mut trans_end = coding_end + pre_post_coding;

        let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[trans_start], trans_end - trans_start, Some(coding_start), Some(coding_end), false, "");

        let mut pos = 28;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = cat!(sub(&ref_, 0, 1), "T");
        let mut var = VariantData::new(CHR_2, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data, &fx.ref_genome);
        assert_eq!(VariantEffect::Frameshift, impact.top_effect().unwrap());

        assert_eq!("c.6_7insA", impact.coding.hgvs);
        // assertEquals("p.Ala5fs", impact.pc().hgvs);

        // test 2: inframe deletion of 2 codons leading up to exon boundary
        ref_bases = generate_test_bases(pre_gene);

        ref_bases += &generate_test_bases(pre_gene);
        ref_bases += &rc("TAA");
        ref_bases += &get_amino_acids_codons("GHC", true);
        ref_bases += &generate_test_bases(pre_gene); // intron
        ref_bases += &get_amino_acids_codons("LRD", true);
        ref_bases += &rc("ATG");
        ref_bases += &generate_test_bases(pre_gene);

        let mut coding_length = 8 * 3;

        trans_start = pre_gene;
        coding_start = trans_start + pre_post_coding;
        coding_end = 53;

        // 10-19 3'UTR, 20-22 X, 23-25 C, 26-28 H 29-31 G, 32-41 intron 42-44 D, 45-47 R, 48-50 L, 51-53 M, 54-63 5'UTR
        // TTACTCGTGTCCTAATAAGTCACAAGCCAT GATCGATCGA
        // 110                            140
        let exon_starts = [10, 42];

        let neg_trans = create_trans_exons(GENE_ID_2, TRANS_ID_2, NEG_STRAND, &exon_starts, 21, Some(coding_start), Some(coding_end), false, "");

        fx.ref_genome.map.insert(CHR_3.to_string(), ref_bases.clone());

        // deleting last 2 amino acids of second coding exon
        pos = 25;
        ref_ = sub(&ref_bases, pos, pos + 7);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_3, pos, &ref_, &alt);

        impact = classify(&var, neg_trans, &fx.ref_genome);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());
        assert_eq!("c.13_18delGGACAC", impact.coding.hgvs);
        assert_eq!("p.Gly5_His6del", impact.pc().hgvs);
    }

    #[test]
    fn test_synonymous_missense_impacts_old()
    {
        let mut fx = Fx::new();
        // inferior version of the above, remove once new tests are complete
        let mut ref_genome = MockRefGenome::default();

        let exon_starts = [0, 100, 200];

        // codons start on at 10, 13, 16 etc
        let coding_start = Some(10);
        let coding_end = Some(250);

        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 80, coding_start, coding_end, false, "");

        

        let mut ref_bases = generate_test_bases(300);

        // set the specific AAs for a region of this mock ref_ genome
        // S: TCA -> TCG - so last base of codon can change
        // I: ATC -> ATT
        // L: TTA -> CTA - first base changes

        //                     40      43      46
        let mut amino_acid_seq = cat!("TCA", "ATC", "TTA");
        ref_bases = cat!(sub(&ref_bases, 0, 40), amino_acid_seq, ref_bases[49..].to_string());
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());

        // test SNVs at the 1st and 3rd codon bases

        // last base of a codon changes
        let mut codon_pos = 40;
        let mut codon = sub(&ref_bases, codon_pos, codon_pos + 3);
        let mut amino_acid = find_amino_acid_for_codon(&codon);

        let mut alt = "G".to_string();
        let mut syn_codon = cat!(sub(&ref_bases, codon_pos, codon_pos + 2), alt);
        assert_eq!(amino_acid, find_amino_acid_for_codon(&syn_codon));

        let mut pos = 42;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos, &ref_genome);
        assert_eq!(VariantEffect::Synonymous, impact.top_effect().unwrap());

        // now missense
        pos = 41;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "T".to_string();
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert_eq!(VariantEffect::Missense, impact.top_effect().unwrap());

        // first base of a codon changes
        codon_pos = 46;
        codon = sub(&ref_bases, codon_pos, codon_pos + 3);
        amino_acid = find_amino_acid_for_codon(&codon);

        alt = "C".to_string();
        syn_codon = cat!(alt, sub(&ref_bases, codon_pos + 1, codon_pos + 3));
        assert_eq!(amino_acid, find_amino_acid_for_codon(&syn_codon));

        pos = 46;
        ref_ = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert_eq!(VariantEffect::Synonymous, impact.top_effect().unwrap());

        // test for an MNV spanning 3 codons - first in the last codon pos at 42 then all the next and 2 into the final one at 46
        pos = 42;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = cat!("G", "ATT", "C");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert_eq!(VariantEffect::Synonymous, impact.top_effect().unwrap());

        // now missense by change middle codon
        alt = cat!("G", "AAT", "C");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert_eq!(VariantEffect::Missense, impact.top_effect().unwrap());

        // test reverse strand
        let trans_data_neg = create_trans_exons(GENE_ID_2, TRANS_ID_2, NEG_STRAND, &exon_starts, 80, coding_start, coding_end, false, "");

        // coding starts at 250 so codons start at 250, 247, 244 etc
        // still change the AA sequence S, I then L - for the range 239-241, 242-244 and 245-247
        let mut amino_acid_seq_rev = rc(&amino_acid_seq);
        ref_bases = cat!(sub(&ref_bases, 0, 239), amino_acid_seq_rev, ref_bases[248..].to_string());
        ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());

        codon_pos = 242;
        codon = sub(&ref_bases, codon_pos, codon_pos + 3);
        amino_acid = find_amino_acid_for_codon(&rc(&codon));

        // change last base of a codon, which is the first base
        alt = "A".to_string();
        syn_codon = cat!(alt, sub(&ref_bases, codon_pos + 1, codon_pos + 3));
        assert_eq!(amino_acid, find_amino_acid_for_codon(&rc(&syn_codon)));

        pos = 242;
        ref_ = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(VariantEffect::Synonymous, impact.top_effect().unwrap());

        // now missense
        pos = 243;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "T".to_string();
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(VariantEffect::Missense, impact.top_effect().unwrap());

        // test a MNV spanning 3 codons as before
        pos = 241;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = rc(&cat!("G", "ATT", "C"));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(VariantEffect::Synonymous, impact.top_effect().unwrap());

        // and missense
        alt = rc(&cat!("G", "AAT", "C"));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(VariantEffect::Missense, impact.top_effect().unwrap());
    }

    #[test]
    fn test_extra_coding_bases()
    {
        let mut fx = Fx::new();
        let mut ref_genome = MockRefGenome::default();

        let chr1 = generate_test_bases(100);
        ref_genome.map.insert(CHR_1.to_string(), chr1.clone());

        let exon_starts = [10, 30, 50, 70];

        // codons start at 10, 13, 16 etc
        let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, Some(15), Some(75), false, "");
        let extra = |exon: &hmf_common::ensembl::ExonData, start: i32, n: i32, up: bool| {
            crate::impact::protein::get_extra_bases_from(trans_data_pos_strand, &ref_genome, CHR_1, exon, start, n, up).unwrap().unwrap()
        };

        let exon = &trans_data_pos_strand.exons[1];
        assert_eq!(sub(&chr1, 31, 32), extra(exon, 32, 1, false));
        assert_eq!(sub(&chr1, 30, 32), extra(exon, 32, 2, false));
        assert_eq!(cat!(sub(&chr1, 20, 21), sub(&chr1, 30, 32)), extra(exon, 32, 3, false));
        assert_eq!(cat!(sub(&chr1, 39, 41), sub(&chr1, 50, 52)), extra(exon, 38, 4, true));

        let _trans_data_neg_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, Some(15), Some(75), false, "");

        let exon = &trans_data_pos_strand.exons[2];
        assert_eq!(cat!(sub(&chr1, 39, 41), sub(&chr1, 50, 52)), extra(exon, 52, 4, false));
        assert_eq!(cat!(sub(&chr1, 59, 61), sub(&chr1, 70, 74)), extra(exon, 58, 6, true));
    }
