//! PhasedVariantsTest: INDELs of a local phase set reclassified together

use super::*;
#[allow(unused_imports)]
use crate::impact::effects::VariantEffect;
use crate::impact::phased::{process_phased_variants, PhasedVariantClassifier};
use crate::impact::VariantTransImpact;
use crate::variant::NO_LOCAL_PHASE_SET;

/// The test class's classifier: variants classified on one transcript, then phased together
struct Harness {
    ref_genome: MockRefGenome,
    ref_bases: String,
}

impl Harness {
    fn new() -> Harness {
        let ref_genome = create_mock_genome(201);
        let ref_bases = ref_genome.bases().to_string();
        Harness { ref_genome, ref_bases }
    }

    /// classifyVariant: the impact, added to the variant
    fn classify(&mut self, var: &mut VariantData, trans: &'static TranscriptData) -> VariantTransImpact {
        let impact = classify(var, trans, &self.ref_genome);
        var.add_impact(GENE_NAME_1, impact.clone());
        impact
    }

    /// processPhasedVariants(NO_LOCAL_PHASE_SET) on the variants classified since the last call, in order
    fn process(&mut self, vars: &[&VariantData]) -> Vec<VariantData> {
        let mut classifier = PhasedVariantClassifier::default();
        for v in vars {
            classifier.add_variant((*v).clone());
        }
        process_phased_variants(&mut classifier, NO_LOCAL_PHASE_SET, &self.ref_genome)
    }

    fn impact(&self, var: &VariantData) -> VariantTransImpact {
        var.gene_impacts(GENE_NAME_1).unwrap()[0].clone()
    }
}

    #[test]
    fn test_inframe_del_pair()
    {
        let mut h = Harness::new();
        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");

        // 2 DELs making inframe DEL
        let mut pos = 20;

        let mut ref_ = sub(&h.ref_bases, pos, pos + 3);
        let mut alt = sub(&h.ref_bases, pos, pos + 1);
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data_pos);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 25;
        ref_ = sub(&h.ref_bases, pos, pos + 2);
        alt = sub(&h.ref_bases, pos, pos + 1);
        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data_pos);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.phased_frameshift());
        assert!(!impact1.coding.is_frameshift);
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.has_effect(VariantEffect::PhasedInframeDeletion));

        // as before but with an extra base so remaining out of frame
        var1.impacts_mut().clear();
        impact1 = h.classify(&mut var1, trans_data_pos);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 25;
        ref_ = sub(&h.ref_bases, pos, pos + 3);
        alt = sub(&h.ref_bases, pos, pos + 1);
        var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        impact2 = h.classify(&mut var2, trans_data_pos);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        // 2 DELs which overlap ref_ codon bases but can form a frameshift
        // codons 13-15, 16-18, 19-21, 22-24 etc
        pos = 16;
        ref_ = sub(&h.ref_bases, pos, pos + 3);
        alt = sub(&h.ref_bases, pos, pos + 1); // deletes 17 & 18;
        var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        impact1 = h.classify(&mut var1, trans_data_pos);
        assert_eq!(16, impact1.pc().ref_coding_base_start());
        assert_eq!(21, impact1.pc().ref_coding_base_end());

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 20;
        ref_ = sub(&h.ref_bases, pos, pos + 2);
        alt = sub(&h.ref_bases, pos, pos + 1); // deletes 21 being the 3rd codon base which they share
        var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        impact2 = h.classify(&mut var2, trans_data_pos);
        assert_eq!(19, impact2.pc().ref_coding_base_start());
        assert_eq!(24, impact2.pc().ref_coding_base_end());

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.phased_frameshift());
        assert!(!impact1.coding.is_frameshift);
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.has_effect(VariantEffect::PhasedInframeDeletion));
    }

    #[test]
    fn test_inframe_del_ins_pair()
    {
        let mut h = Harness::new();
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");

        // DEL and INS making inframe insert
        let mut pos = 20;

        let mut ref_ = sub(&h.ref_bases, pos, pos + 3);
        let mut alt = sub(&h.ref_bases, pos, pos + 1);
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data_neg);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 25;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AAAAAAAA");
        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data_neg);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.phased_frameshift());
        assert!(!impact1.coding.is_frameshift);
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.has_effect(VariantEffect::PhasedInframeInsertion));
        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.has_effect(VariantEffect::PhasedInframeInsertion));

        // again but remaining out-of-frame
        var1.impacts_mut().clear();

        impact1 = h.classify(&mut var1, trans_data_neg);

        pos = 25;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AAAAAAA");
        var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        impact2 = h.classify(&mut var2, trans_data_neg);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));
    }

    #[test]
    fn test_del_snv_overlap_pair()
    {
        let mut h = Harness::new();
        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");

        // DEL and SNV, overlapping in the same codon (ala EGFR)
        let mut pos = 26;

        let mut ref_ = sub(&h.ref_bases, pos, pos + 10);
        let mut alt = sub(&h.ref_bases, pos, pos + 1);

        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data_pos);

        assert!(impact1.has_effect(VariantEffect::InframeDeletion));
        assert_eq!("IDRS", impact1.pc().ref_amino_acids);
        assert_eq!("M", impact1.pc().alt_amino_acids);
        pos = 36;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = "T".to_string();

        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data_pos);

        assert!(impact2.has_effect(VariantEffect::Synonymous));
        assert_eq!("S", impact2.pc().alt_amino_acids);

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.phased_frameshift());
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.has_effect(VariantEffect::PhasedInframeDeletion));

        // again but with the DEL conservative and not overlapping with the SNV despite its ref_ codon bases overlapping
        // 13-18 are deleted, SNV is at 20
        pos = 12;
        ref_ = sub(&h.ref_bases, pos, pos + 7);
        alt = sub(&h.ref_bases, pos, pos + 1);

        var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        impact1 = h.classify(&mut var1, trans_data_pos);

        assert!(impact1.has_effect(VariantEffect::InframeDeletion));
        assert_eq!("SIDR", impact1.pc().ref_amino_acids);
        assert_eq!("SR", impact1.pc().alt_amino_acids);

        pos = 20;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = "T".to_string();

        var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        impact2 = h.classify(&mut var2, trans_data_pos);

        assert!(impact2.has_effect(VariantEffect::Missense));
        assert_eq!("R", impact2.pc().ref_amino_acids);
        assert_eq!("L", impact2.pc().alt_amino_acids);

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(!impact1.phased_frameshift());
    }

    #[test]
    fn test_missense_synonymous_del_ins_pair()
    {
        let mut h = Harness::new();
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");

        // DEL and INS making missense and symonymous effects
        let mut pos = 20;

        let mut ref_ = sub(&h.ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, "AGCT");
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data_neg);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 25;
        ref_ = sub(&h.ref_bases, pos, pos + 5);
        alt = sub(&h.ref_bases, pos, pos + 1);
        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data_neg);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.phased_frameshift());
        assert!(!impact1.coding.is_frameshift);
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.has_effect(VariantEffect::PhasedMissense));
        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.has_effect(VariantEffect::PhasedMissense));
        assert_eq!("p.Ser5_Ile6delinsAlaIle", impact1.pc().hgvs);

        // now synonynous
        pos = 20;
        ref_ = sub(&h.ref_bases, pos, pos + 5);
        alt = sub(&h.ref_bases, pos, pos + 1);
        var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(2, "", "", 0);
        impact1 = h.classify(&mut var1, trans_data_neg);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 28;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "ATCG");
        var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(2, "", "", 0);
        impact2 = h.classify(&mut var2, trans_data_neg);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);

        assert!(impact1.phased_frameshift());
        assert!(!impact1.coding.is_frameshift);
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.has_effect(VariantEffect::PhasedSynonymous));
        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.has_effect(VariantEffect::PhasedSynonymous));
        assert_eq!("p.Arg4_Arg8=", impact1.pc().hgvs);
    }

    #[test]
    fn test_multi_indel()
    {
        let mut h = Harness::new();
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");

        // DEL and INS making missense and symonymous effects
        // SNVs at either end are ignored
        let mut pos = 14;
        let mut ref_ = sub(&h.ref_bases, pos, pos + 1);
        let mut alt = generate_alt(&ref_);
        let mut var_ignore1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var_ignore1.set_variant_details(1, "", "", 0);
        let mut impact_ignore1 = h.classify(&mut var_ignore1, trans_data_neg);
        assert!(impact_ignore1.has_effect(VariantEffect::Missense));

        pos = 20;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AG");
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data_neg);

        assert!(impact1.coding.is_frameshift);
        assert!(impact1.has_effect(VariantEffect::Frameshift));
        assert_eq!("p.Ile48fs", impact1.pc().hgvs);

        pos = 25;
        ref_ = sub(&h.ref_bases, pos, pos + 8);
        alt = sub(&h.ref_bases, pos, pos + 1);
        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data_neg);

        assert!(impact2.coding.is_frameshift);
        assert!(impact2.has_effect(VariantEffect::Frameshift));

        pos = 35;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = cat!(ref_, "GG");
        let mut var3 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var3.set_variant_details(1, "", "", 0);
        let mut impact3 = h.classify(&mut var3, trans_data_neg);

        assert!(impact3.coding.is_frameshift);
        assert!(impact3.has_effect(VariantEffect::Frameshift));
        assert_eq!("p.Arg42fs", impact3.pc().hgvs);

        pos = 40;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = generate_alt(&ref_);
        let mut var_ignore2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var_ignore2.set_variant_details(1, "", "", 0);
        let mut impact_ignore2 = h.classify(&mut var_ignore2, trans_data_neg);
        assert!(impact_ignore2.has_effect(VariantEffect::Synonymous));

        pos = 46;
        ref_ = sub(&h.ref_bases, pos, pos + 1);
        alt = generate_alt(&ref_);
        let mut var_ignore3 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var_ignore3.set_variant_details(1, "", "", 0);
        let mut impact_ignore3 = h.classify(&mut var_ignore3, trans_data_neg);
        assert!(impact_ignore3.has_effect(VariantEffect::Synonymous));

        let out = h.process(&[&var_ignore1, &var1, &var2, &var3, &var_ignore2, &var_ignore3]);
        var_ignore1 = out[0].clone();
        impact_ignore1 = h.impact(&var_ignore1);
        var1 = out[1].clone();
        impact1 = h.impact(&var1);
        var2 = out[2].clone();
        impact2 = h.impact(&var2);
        var3 = out[3].clone();
        impact3 = h.impact(&var3);
        var_ignore2 = out[4].clone();
        impact_ignore2 = h.impact(&var_ignore2);
        var_ignore3 = out[5].clone();
        impact_ignore3 = h.impact(&var_ignore3);

        assert!(!impact1.coding.is_frameshift);
        assert!(!impact1.has_effect(VariantEffect::Frameshift));
        assert!(impact1.phased_frameshift());
        assert!(impact1.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(impact2.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(impact3.has_effect(VariantEffect::PhasedInframeDeletion));

        assert!(!impact_ignore1.phased_frameshift());
        assert!(!impact_ignore2.phased_frameshift());
        assert!(!impact_ignore3.phased_frameshift());
        assert!(!impact_ignore1.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(!impact_ignore2.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(!impact_ignore3.has_effect(VariantEffect::PhasedInframeDeletion));
    }

    #[test]
    fn test_ref_condon_overlap()
    {
        let mut h = Harness::new();
        // variants remain as-is since do not overlap each other's impacted ref_ or alt codons

        let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");

        // positions:      012345678901234567890123     456789012345     67890
        let mut ref_bases = cat!("XCCCCCCCCCCAAAAAAAAAACCC", "GGTGGCGGCGGC", "AAAAATTTTTAAAAATTTTT");

        let mut pos = 26;
        let mut ref_ = sub(&ref_bases, pos, pos + 7);
        let mut alt = sub(&ref_bases, pos, pos + 1);
        let mut var_ignore1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var_ignore1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var_ignore1, trans_data);
        assert!(impact1.has_effect(VariantEffect::InframeDeletion));

        pos = 35;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "T".to_string();
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var1, trans_data);
        assert!(impact2.has_effect(VariantEffect::Missense));

        let out = h.process(&[&var_ignore1, &var1]);
        var_ignore1 = out[0].clone();
        impact1 = h.impact(&var_ignore1);
        var1 = out[1].clone();
        impact2 = h.impact(&var1);

        assert!(!impact1.phased_frameshift());
        assert!(impact1.has_effect(VariantEffect::InframeDeletion));
        assert!(!impact2.phased_frameshift());
        assert!(impact2.has_effect(VariantEffect::Missense));
    }

    #[test]
    fn test_splice_edge_overlap()
    {
        let mut h = Harness::new();
        let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 20], 12, Some(1), Some(25), false, "");

        //                              intron
        //                           10        20        30        40
        // positions:      01234567890123456789012345678901234567890
        let mut ref_bases = "XATGCATGGGCAGTTTTTTTTTTTTTAAAAATTTTTAAAAATTTTT".to_string();

        h.ref_genome.map.insert(CHR_2.to_string(), ref_bases.clone());

        // start
        let mut pos = 8;
        let mut ref_ = sub(&ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, "A");
        let mut var1 = VariantData::new(CHR_2, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 9;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "T".to_string();
        let mut var2 = VariantData::new(CHR_2, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data);
        assert!(impact2.has_effect(VariantEffect::Synonymous));

        pos = 11;
        ref_ = sub(&ref_bases, pos, pos + 9); // AGGTAACTT>A
        alt = sub(&ref_bases, pos, pos + 1);
        let mut var3 = VariantData::new(CHR_2, pos, &ref_, &alt);
        var3.set_variant_details(1, "", "", 0);
        let mut impact3 = h.classify(&mut var3, trans_data);
        assert!(impact3.has_effect(VariantEffect::Frameshift));

        let out = h.process(&[&var1, &var2, &var3]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);
        var3 = out[2].clone();
        impact3 = h.impact(&var3);

        assert_eq!("p.Gln4Ser", impact1.hgvs_protein());
        assert!(impact1.phased_frameshift());
        assert!(impact1.has_effect(VariantEffect::PhasedMissense));
        assert!(impact2.has_effect(VariantEffect::PhasedMissense));
        assert!(impact3.has_effect(VariantEffect::PhasedMissense));
    }

    #[test]
    fn test_mixed_overlaps()
    {
        let mut h = Harness::new();
        let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 100], 80, Some(10), Some(190), false, "");
        // positions:      01234567890123456789012     34567890123456789012345678901234567
        let mut ref_bases = cat!("XCCCCCCCCCCAAAAAAAAAACC", "TTGAGAAGTCTGCCAGCACTGAGAGGAAAATTAAC", "AAAAATTTTTAAAAATTTTT");

        let mut pos = 23;
        let mut ref_ = sub(&ref_bases, pos, pos + 8);
        let mut alt = sub(&ref_bases, pos, pos + 1);
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data);
        assert!(impact1.has_effect(VariantEffect::Frameshift));

        pos = 34;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "T".to_string();
        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data);
        assert!(impact2.has_effect(VariantEffect::Synonymous));

        pos = 35;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = cat!(ref_, "T");
        let mut var3 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var3.set_variant_details(1, "", "", 0);
        let mut impact3 = h.classify(&mut var3, trans_data);
        assert!(impact3.has_effect(VariantEffect::Frameshift));

        pos = 36;
        ref_ = sub(&ref_bases, pos, pos + 13);
        alt = sub(&ref_bases, pos, pos + 1);
        let mut var4 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var4.set_variant_details(1, "", "", 0);
        let mut impact4 = h.classify(&mut var4, trans_data);
        assert!(impact4.has_effect(VariantEffect::InframeDeletion));

        // the last variant doesn't overlap
        pos = 53;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "A".to_string();
        let mut var5 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var5.set_variant_details(1, "", "", 0);
        let mut impact5 = h.classify(&mut var5, trans_data);
        assert!(impact5.has_effect(VariantEffect::Synonymous));

        let out = h.process(&[&var1, &var2, &var3, &var4, &var5]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);
        var3 = out[2].clone();
        impact3 = h.impact(&var3);
        var4 = out[3].clone();
        impact4 = h.impact(&var4);
        var5 = out[4].clone();
        impact5 = h.impact(&var5);

        assert!(!impact2.coding.is_frameshift);
        assert!(!impact2.has_effect(VariantEffect::Frameshift));
        assert!(impact2.phased_frameshift());
        assert!(impact1.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(impact2.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(impact3.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(impact4.has_effect(VariantEffect::PhasedInframeDeletion));
        assert!(impact5.has_effect(VariantEffect::Synonymous));
    }

    #[test]
    fn test_mixed_overlaps2()
    {
        let mut h = Harness::new();
        let trans_data = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[0], 80, Some(5), Some(70), false, "");

        // positions:      0123456789012345678901234567890123456789012345678901234567
        //                 0         10        20        30        40        50
        let mut ref_bases = cat!("XCCCCATGAAAAAAAAAACCCCCCCCCCGGGGGGGAAGCTGATGTTCAGGAGTG", generate_test_bases(40));
        h.ref_genome.map.insert(CHR_1.to_string(), ref_bases.clone());

        // pos(12:50384540-50384542) variant(GA>G)"
        let mut pos = 40;
        let mut ref_ = sub(&ref_bases, pos, pos + 2);
        let mut alt = sub(&ref_, 0, 1);
        let mut var1 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var1.set_variant_details(1, "", "", 0);
        let mut impact1 = h.classify(&mut var1, trans_data);

        // pos(12:50384543) variant(G>C)
        pos = 43;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = "C".to_string();
        let mut var2 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var2.set_variant_details(1, "", "", 0);
        let mut impact2 = h.classify(&mut var2, trans_data);

        // pos(12:50384545-50384546) variant(T>TCTCAAGA)
        pos = 45;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = cat!(ref_, "CTCAAGA");
        let mut var3 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var3.set_variant_details(1, "", "", 0);
        let mut impact3 = h.classify(&mut var3, trans_data);

        assert!(impact3.coding.is_frameshift);
        assert!(impact3.has_effect(VariantEffect::Frameshift));

        // pos(12:50384548-50384549) variant(G>GCTT)
        pos = 48;
        ref_ = sub(&ref_bases, pos, pos + 1);
        alt = cat!(ref_, "CTT");
        let mut var4 = VariantData::new(CHR_1, pos, &ref_, &alt);
        var4.set_variant_details(1, "", "", 0);
        let mut impact4 = h.classify(&mut var4, trans_data);

        let out = h.process(&[&var1, &var2, &var3, &var4]);
        var1 = out[0].clone();
        impact1 = h.impact(&var1);
        var2 = out[1].clone();
        impact2 = h.impact(&var2);
        var3 = out[2].clone();
        impact3 = h.impact(&var3);
        var4 = out[3].clone();
        impact4 = h.impact(&var4);

        assert_eq!("CTGATGTTCAGG", impact4.pc().ref_codon_bases);
        assert_eq!("CTGTCTTCTCAAGACAGCTTG", impact4.pc().alt_codon_bases);
        assert_eq!("p.Pro8_His10delinsGlnAlaValLeuArgArg", impact4.hgvs_protein());

        assert!(!impact3.coding.is_frameshift);
        assert!(!impact3.has_effect(VariantEffect::Frameshift));
        assert!(impact3.phased_frameshift());
        assert!(impact3.has_effect(VariantEffect::PhasedInframeInsertion));
    }


