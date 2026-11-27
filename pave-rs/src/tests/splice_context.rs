//! SpliceContextTest

use super::*;
use crate::impact::effects::{CodingType, RegionType, VariantEffect};
use crate::impact::{SE_END, SE_START};

use crate::impact::splice::determine_indel_splice_impact;
use crate::impact::SpliceImpactType;
use hmf_common::ensembl::ExonData;

fn ref_alt_splice_bases_match(variant: &VariantData, ref_genome: &MockRefGenome, exon: &ExonData, pos_strand: bool) -> bool {
    determine_indel_splice_impact(variant, ref_genome, exon, pos_strand).unwrap() != SpliceImpactType::BaseShift
}

    #[test]
    fn test_delete_splice_bases_pos_strand()
    {
        let m_ref_genome = create_mock_genome(41);
        let m_ref_bases = m_ref_genome.bases().to_string();
        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[20], 10, Some(20), Some(30), false, "");

        let mut exon = &trans_data_pos.exons[0];

        // test deleting the whole splice region
        let mut pos = 15;
        let mut ref_ = sub(&m_ref_bases, pos, pos + 6);
        let mut alt = sub(&ref_, 0, 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert_eq!(SpliceImpactType::RegionDeleted, determine_indel_splice_impact(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()).unwrap());

        // and on the donor side
        pos = 29;
        ref_ = sub(&m_ref_bases, pos, pos + 7);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert_eq!(SpliceImpactType::RegionDeleted, determine_indel_splice_impact(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()).unwrap());

        // deletes outside region are non-disruptive
        pos = 15;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 19;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 28;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 35;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // acceptor positions - modify A3 only but doesn't become a G so is ok
        pos = 17;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = "A".to_string();
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // acceptor positions - delete A3 and cannot pull exonic bases, so remains disruptive
        pos = 15;
        ref_ = sub(&m_ref_bases, pos, pos + 3);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // deletes A3 and A2
        pos = 16;
        ref_ = sub(&m_ref_bases, pos, pos + 3);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // now with A3 becoming disruptive
        pos = 16;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // delete A2 only
        pos = 17;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // delete A2 and A1
        pos = 17;
        ref_ = sub(&m_ref_bases, pos, pos + 3);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // donor bases - starting with delete of D-1 only
        pos = 29;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 30;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 34;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));
    }

    #[test]
    fn test_delete_splice_bases_neg_strand()
    {
        let m_ref_genome = create_mock_genome(41);
        let m_ref_bases = m_ref_genome.bases().to_string();
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[20], 10, Some(20), Some(30), false, "");

        let mut exon = &trans_data_neg.exons[0];

        // deletes outside region are non-disruptive
        let mut pos = 12;
        let mut ref_ = sub(&m_ref_bases, pos, pos + 3);
        let mut alt = sub(&ref_, 0, 1);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        pos = 20;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        pos = 29;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        pos = 33;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // acceptor positions - delete A3 only but doesn't become a G so is ok
        pos = 32;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // deletes A3 and A2
        pos = 31;
        ref_ = sub(&m_ref_bases, pos, pos + 3);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // now with A3 becoming disruptive
        pos = 32;
        ref_ = sub(&m_ref_bases, pos, pos + 3);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // delete A2 only
        pos = 31;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // delete A2 and A1
        pos = 30;
        ref_ = sub(&m_ref_bases, pos, pos + 3);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // donor bases - starting with delete of D-1 only
        pos = 19;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // D1
        pos = 18;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // D4
        pos = 15;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // D5
        pos = 14;
        ref_ = sub(&m_ref_bases, pos, pos + 2);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));
    }

    #[test]
    fn test_insert_splice_bases_pos_strand()
    {
        let m_ref_genome = create_mock_genome(41);
        let m_ref_bases = m_ref_genome.bases().to_string();
        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[20], 10, Some(20), Some(30), false, "");

        let mut exon = &trans_data_pos.exons[0];

        let mut pos = 17;
        let mut ref_ = sub(&m_ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, generate_alt(&sub(&m_ref_bases, pos + 1, pos + 2)));
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 2);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // now insert bases that match the A1-3 bases - these will be considered disruptive before the exon
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 3);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // starting at A2, inserting matching bases
        pos = 18;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, 17, 20);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // starting at A1 - makes no difference since inserts into exon
        pos = 19;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = "AAAA".to_string();
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // test D-1 through D5 still on positive strand
        pos = 30;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 2);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 30;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 5);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 30;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 6);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // from D2
        pos = 32;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 6);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 32;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = cat!(ref_, generate_alt(&sub(&m_ref_bases, pos + 1, pos + 4)));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        // from D4
        pos = 34;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 2);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));

        pos = 34;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = cat!(ref_, generate_alt(&sub(&m_ref_bases, pos + 1, pos + 2)));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_pos.pos_strand()));
    }

    #[test]
    fn test_insert_splice_bases_neg_strand()
    {
        let m_ref_genome = create_mock_genome(41);
        let m_ref_bases = m_ref_genome.bases().to_string();
        // repeat for negative strand
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[20], 10, Some(20), Some(30), false, "");

        let mut exon = &trans_data_neg.exons[0];

        // inserts outside the range are non-disruptive
        let mut pos = 33;
        let mut ref_ = sub(&m_ref_bases, pos, pos + 1);
        let mut alt = cat!(ref_, "AAAA");
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        pos = 29;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AAAA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        pos = 14;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AAAA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        pos = 20;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = cat!(ref_, "AAAA");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));


        // first the splice acceptors
        pos = 30;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 4);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert_eq!(SpliceImpactType::OutsideRange, determine_indel_splice_impact(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()).unwrap());

        alt = cat!(ref_, generate_alt(&sub(&m_ref_bases, pos + 1, pos + 4)));
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert_eq!(SpliceImpactType::OutsideRange, determine_indel_splice_impact(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()).unwrap());

        // from A2 - first with A3 not changing to the required C
        pos = 32;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 4);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        alt = cat!(ref_, "CCCC");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // then the donors

        // from D5
        pos = 15;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 4);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // ok if covers the full range with homologous bases
        alt = sub(&m_ref_bases, pos, pos + 6);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        alt = cat!(ref_, "A");
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));

        // D1
        pos = 19;
        ref_ = sub(&m_ref_bases, pos, pos + 1);
        alt = sub(&m_ref_bases, pos, pos + 4);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        assert!(!ref_alt_splice_bases_match(&var, &m_ref_genome, exon, trans_data_neg.pos_strand()));
    }
