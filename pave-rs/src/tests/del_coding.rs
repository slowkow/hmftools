//! DelCodingContextTest

use super::*;
use crate::impact::effects::{CodingType, RegionType, VariantEffect};
use crate::impact::{SE_END, SE_START};

    #[test]
    fn test_inframe_dels()
    {
        let ref_genome = create_mock_genome(150);
        let ref_bases = ref_genome.bases().to_string();

        

        let trans_data_pos = create_pos_transcript();

        // DELs crossing codons
        let mut pos = 30;
        let mut ref_codon_bases = sub(&ref_bases, pos, pos + 9);

        let mut ref_ = sub(&ref_codon_bases, 0, 4);
        let mut alt = sub(&ref_codon_bases, 0, 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos, &ref_genome);

        assert!(!impact.coding.is_frameshift);
        assert_eq!(sub(&ref_codon_bases, 0, 6), impact.pc().ref_codon_bases);

        let mut alt_codon_bases = cat!(alt, sub(&ref_codon_bases, 4, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());

        ref_ = sub(&ref_codon_bases, 1, 5);
        alt = sub(&ref_codon_bases, 1, 2);
        var = VariantData::new(CHR_1, 31, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert!(!impact.coding.is_frameshift);
        assert_eq!(sub(&ref_codon_bases, 0, 6), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 1), alt, sub(&ref_codon_bases, 5, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        ref_ = sub(&ref_codon_bases, 2, 6);
        alt = sub(&ref_codon_bases, 2, 3);
        var = VariantData::new(CHR_1, 32, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert!(!impact.coding.is_frameshift);
        assert_eq!(ref_codon_bases, impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 2), alt, sub(&ref_codon_bases, 6, 9));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());

        // repeat for negative strand
        let trans_data_neg = create_neg_transcript();

        // for codons 80-78 and 77-75
        pos = 72;
        ref_codon_bases = sub(&ref_bases, pos, pos + 9);

        ref_ = sub(&ref_codon_bases, 0, 4);
        alt = sub(&ref_codon_bases, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert!(!impact.coding.is_frameshift);
        assert_eq!(sub(&ref_codon_bases, 0, 6), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(alt, sub(&ref_codon_bases, 4, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());

        // shift forward one base each time
        ref_ = sub(&ref_codon_bases, 1, 5);
        alt = sub(&ref_codon_bases, 1, 2);
        var = VariantData::new(CHR_1, { pos += 1; pos }, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert!(!impact.coding.is_frameshift);
        assert_eq!(sub(&ref_codon_bases, 0, 6), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 1), alt, sub(&ref_codon_bases, 5, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        ref_ = sub(&ref_codon_bases, 2, 6);
        alt = sub(&ref_codon_bases, 2, 3);
        var = VariantData::new(CHR_1, { pos += 1; pos }, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert!(!impact.coding.is_frameshift);
        assert_eq!(ref_codon_bases, impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 2), alt, sub(&ref_codon_bases, 6, 9));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());
    }

    #[test]
    fn test_dels_across_splice()
    {
        // DELs crossing exon boundaries
        let ref_genome = create_mock_genome(150);
        let ref_bases = ref_genome.bases().to_string();

        

        let trans_data_pos = create_pos_transcript();

        let mut pos = 20;
        let mut ref_ = sub(&ref_bases, pos, pos + 4);
        let mut alt = sub(&ref_bases, pos, pos + 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(6, impact.coding.coding_base);
        assert!(!impact.coding.is_frameshift);

        pos = 19;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert!(impact.coding.is_frameshift);

        // frameshift since 18 remains while 19-20 are deleted
        pos = 18;
        ref_ = sub(&ref_bases, pos, pos + 6);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert!(impact.coding.is_frameshift);

        // inframe
        pos = 17;
        ref_ = sub(&ref_bases, pos, pos + 6);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert!(!impact.coding.is_frameshift);

        assert_eq!(sub(&ref_bases, 15, 21), impact.pc().ref_codon_bases);

        let mut alt_codon_bases = sub(&ref_bases, 15, 18);
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);
        assert!(impact.has_effect(VariantEffect::InframeDeletion));

        // starting before the next exon
        pos = 28;
        ref_ = sub(&ref_bases, pos, pos + 7);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert!(impact.coding.is_frameshift);

        pos = 29; // deletes first 3 bases of exon
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert!(!impact.coding.is_frameshift);
        assert_eq!(7, impact.coding.coding_base);
        assert_eq!(30, impact.coding.coding_position_range[SE_START]);
        assert_eq!(33, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.7_9delTCG", impact.coding.hgvs);

        // inframe
        pos = 28;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert!(!impact.coding.is_frameshift);

        assert_eq!(sub(&ref_bases, 30, 36), impact.pc().ref_codon_bases);

        alt_codon_bases = sub(&ref_bases, 33, 36);
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);
        assert!(impact.has_effect(VariantEffect::InframeDeletion));

        // negative strand
        let trans_data_neg = create_neg_transcript();

        pos = 88;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert!(impact.coding.is_frameshift);

        pos = 88;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert!(!impact.coding.is_frameshift);

        pos = 78; // deletes last 2 bases of exon, which are first 2 after splice acceptor, index 7 & 8
        ref_ = sub(&ref_bases, pos, pos + 7);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert!(impact.coding.is_frameshift);
        assert_eq!(7, impact.coding.coding_base);
        assert_eq!(78, impact.coding.coding_position_range[SE_START]);
        assert_eq!(80, impact.coding.coding_position_range[SE_END]);

        pos = 76; // deletes last 3 bases of exon but none of the intron, doesn't span splice
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert!(impact.coding.is_frameshift);
        assert_eq!(7, impact.coding.coding_base);
        assert_eq!(76, impact.coding.coding_position_range[SE_START]);
        assert_eq!(80, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.7_10delCGAT", impact.coding.hgvs);

        // deletes splice acceptor covering coding and intronic bases, coding bases 7-9
        pos = 77;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert!(!impact.coding.is_frameshift);
        assert!(impact.coding.spans_splice_junction);
        assert_eq!(3, impact.coding.deleted_coding_bases);
        assert_eq!(-2, impact.coding.nearest_exon_distance);
        assert_eq!("c.7-1_9delTCGA", impact.coding.hgvs);

        pos = 67;
        ref_ = sub(&ref_bases, pos, pos + 6);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert!(!impact.coding.is_frameshift);
        assert!(impact.coding.spans_splice_junction);
        assert_eq!(3, impact.coding.deleted_coding_bases);
        assert_eq!(3, impact.coding.nearest_exon_distance);
        assert_eq!("c.15_17+2delCGATC", impact.coding.hgvs);

        // stop-codon deletion example
        pos = 12;
        ref_ = sub(&ref_bases, pos, pos + 6);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        // assertTrue(impact.coding.is_frameshift);
        assert_eq!(15, impact.coding.coding_position_range[SE_START]);
        assert_eq!(18, impact.coding.coding_position_range[SE_END]);
        assert_eq!(3, impact.coding.deleted_coding_bases);
        assert_eq!(3, impact.coding.nearest_exon_distance);
        assert!(impact.coding.spans_coding_end);
        assert_eq!("c.43_*2delTCGAT", impact.coding.hgvs);
    }

    #[test]
    fn test_utr_variants()
    {
        // DELs crossing exon boundaries
        let ref_genome = create_mock_genome(150);
        let ref_bases = ref_genome.bases().to_string();

        

        let trans_data_neg = create_neg_transcript();

        // DEL in exon before coding begins - 110 = coding-base 6, pos 112 at cb 8, so actually deletes 9-11
        let mut pos = 112;
        let mut ref_ = sub(&ref_bases, pos, pos + 4);
        let mut alt = sub(&ref_, 0, 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(12, impact.coding.coding_base);
        assert_eq!(112, impact.coding.coding_position_range[SE_START]);
        assert_eq!(116, impact.coding.coding_position_range[SE_END]);
        assert_eq!(3, impact.coding.deleted_coding_bases);
        assert_eq!("c.-11_-9delGAT", impact.coding.hgvs);

        // 5'UTR delete of exonic bases spanning into next intron - coding base 6 at 110, and 9 at 113, deleted bases are 110-112
        pos = 109;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(9, impact.coding.coding_base);
        assert_eq!(3, impact.coding.deleted_coding_bases);
        assert_eq!(110, impact.coding.coding_position_range[SE_START]);
        assert_eq!(113, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-8_-6delCGA", impact.coding.hgvs);

        // starting in the next intron
        pos = 107;
        ref_ = sub(&ref_bases, pos, pos + 6);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(9, impact.coding.coding_base);
        assert_eq!(3, impact.coding.deleted_coding_bases);
        assert_eq!(110, impact.coding.coding_position_range[SE_START]);
        assert_eq!(113, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-8_-6+2delCGATC", impact.coding.hgvs);

        // spanning into the next intro
        pos = 118;
        ref_ = sub(&ref_bases, pos, pos + 6);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(16, impact.coding.coding_base);
        assert_eq!(2, impact.coding.deleted_coding_bases);
        assert_eq!(118, impact.coding.coding_position_range[SE_START]);
        assert_eq!(120, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-16-3_-15delGATCG", impact.coding.hgvs);

        // delete exonic bases just before coding
        let trans_data_pos = create_pos_transcript();

        pos = 12;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);
        assert_eq!(1, impact.coding.coding_base);
        assert!(impact.coding.spans_coding_start);
        assert_eq!(2, impact.coding.deleted_coding_bases);
        assert_eq!(15, impact.coding.coding_position_range[SE_START]);
        assert_eq!(17, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-2_2delATCG", impact.coding.hgvs);
    }

    #[test]
    fn test_inframe_indel_impacts()
    {
        let ref_genome = create_mock_genome(300);
        let ref_bases = ref_genome.bases().to_string();

        let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 100, 200], 80, Some(10), Some(250), false, "");

        

        // inframe conservative deletion of 2 codons (ie 3 -> 1)
        let mut pos = 48; // 3rd base of codon
        let mut ref_ = sub(&ref_bases, pos, pos + 7);
        let mut alt = sub(&ref_bases, pos, pos + 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos_strand, &ref_genome);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());
        assert_eq!("c.40_45delATCGAT", impact.coding.hgvs);

        // inframe disruptive (spans first and last codon boundaries) deletion of 3 codons (ie 4 -> 1)
        pos = 47; // first base of codon
        ref_ = sub(&ref_bases, pos, pos + 10);
        alt = sub(&ref_bases, pos, pos + 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos_strand, &ref_genome);
        assert_eq!(VariantEffect::InframeDeletion, impact.top_effect().unwrap());
    }
