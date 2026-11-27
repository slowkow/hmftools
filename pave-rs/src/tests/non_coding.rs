//! NonCodingContextTest

use super::*;
use crate::impact::effects::{CodingType, RegionType, VariantEffect};
use crate::impact::{SE_END, SE_START};

fn determine_context(variant: &VariantData, trans: &TranscriptData) -> crate::impact::CodingContext {
    let mut cc = crate::impact::CodingContext::default();
    crate::impact::coding::determine_context(variant, trans, &mut cc).unwrap();
    cc
}

fn create_variant(pos: i32) -> VariantData {
    VariantData::new(CHR_1, pos, "A", "C")
}

    #[test]
    fn test_non_coding()
    {
        let exon_starts = [10, 30, 50];

        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, None, None, false, "");

        // intronic
        let mut pos = 22;
        let mut var = create_variant(pos);

        let mut cc = determine_context(&var, trans_data_pos);

        assert_eq!(CodingType::NonCoding, cc.coding_type);
        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(1, cc.exon_rank);
        assert_eq!(11, cc.coding_base); // exonic bases since start
        assert_eq!(2, cc.nearest_exon_distance); // 2 bases into next intron

        pos = 47;
        var = create_variant(pos);

        cc = determine_context(&var, trans_data_pos);

        assert_eq!(CodingType::NonCoding, cc.coding_type);
        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(3, cc.exon_rank);
        assert_eq!(23, cc.coding_base); // exonic bases since start
        assert_eq!(-3, cc.nearest_exon_distance); // 2 bases into next intron

        // exonic
        pos = 35;
        var = create_variant(pos);

        cc = determine_context(&var, trans_data_pos);

        assert_eq!(CodingType::NonCoding, cc.coding_type);
        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(2, cc.exon_rank);
        assert_eq!(17, cc.coding_base); // exonic bases since start
        assert_eq!(0, cc.nearest_exon_distance); // 2 bases into next intron

        // negative strand
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, None, None, false, "");

        // intronic
        pos = 48;
        var = create_variant(pos);

        cc = determine_context(&var, trans_data_neg);

        assert_eq!(CodingType::NonCoding, cc.coding_type);
        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(1, cc.exon_rank);
        assert_eq!(11, cc.coding_base); // exonic bases since start
        assert_eq!(2, cc.nearest_exon_distance); // 2 bases into next intron

        pos = 23;
        var = create_variant(pos);

        cc = determine_context(&var, trans_data_neg);

        assert_eq!(CodingType::NonCoding, cc.coding_type);
        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(3, cc.exon_rank);
        assert_eq!(23, cc.coding_base); // exonic bases since start
        assert_eq!(-3, cc.nearest_exon_distance); // 2 bases into next intron

        // exonic
        pos = 35;
        var = create_variant(pos);

        cc = determine_context(&var, trans_data_neg);

        assert_eq!(CodingType::NonCoding, cc.coding_type);
        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(2, cc.exon_rank);
        assert_eq!(17, cc.coding_base); // exonic bases since start
        assert_eq!(0, cc.nearest_exon_distance); // 2 bases into next intron
    }

    #[test]
    fn test_pos_strand()
    {
        let exon_starts = [10, 30, 50, 70, 90];

        let coding_start = Some(35);
        let coding_end = Some(55);

        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

        // upstream
        let mut pos = 5;
        let mut var = create_variant(pos);

        let mut cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Upstream, cc.region_type);
        assert_eq!(CodingType::Unknown, cc.coding_type);

        // 5' UTR exonic
        pos = 15;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr5p, cc.coding_type);
        assert_eq!(1, cc.exon_rank);
        assert_eq!(11, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);

        // 5' UTR intronic
        pos = 22;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(CodingType::Utr5p, cc.coding_type);
        assert_eq!(1, cc.exon_rank);
        assert_eq!(6, cc.coding_base);
        assert_eq!(2, cc.nearest_exon_distance);

        // 5' UTR exonic in same exon as coding begins
        pos = 32;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr5p, cc.coding_type);
        assert_eq!(2, cc.exon_rank);
        assert_eq!(3, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);

        // 3'UTR exonic in same exon as coding ends
        pos = 59;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr3p, cc.coding_type);
        assert_eq!(3, cc.exon_rank);
        assert_eq!(4, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);

        // 3'UTR intronic
        pos = 82;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(CodingType::Utr3p, cc.coding_type);
        assert_eq!(4, cc.exon_rank);
        assert_eq!(16, cc.coding_base);
        assert_eq!(2, cc.nearest_exon_distance);

        // 3'UTR exonic
        pos = 95;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_pos);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr3p, cc.coding_type);
        assert_eq!(5, cc.exon_rank);
        assert_eq!(22, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);
    }

    #[test]
    fn test_neg_strand()
    {
        let exon_starts = [10, 30, 50, 70, 90];

        let coding_start = Some(35);
        let coding_end = Some(55);

        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

        // upstream
        let mut pos = 105;
        let mut var = create_variant(pos);

        let mut cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Upstream, cc.region_type);
        assert_eq!(CodingType::Unknown, cc.coding_type);

        // 5' UTR exonic
        pos = 95;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr5p, cc.coding_type);
        assert_eq!(1, cc.exon_rank);
        assert_eq!(22, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);

        // 5' UTR intronic
        pos = 88;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(CodingType::Utr5p, cc.coding_type);
        assert_eq!(1, cc.exon_rank);
        assert_eq!(17, cc.coding_base);
        assert_eq!(2, cc.nearest_exon_distance);

        // 5' UTR exonic in same exon as coding begins
        pos = 58;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr5p, cc.coding_type);
        assert_eq!(3, cc.exon_rank);
        assert_eq!(3, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);

        // 3'UTR exonic in same exon as coding ends
        pos = 30;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr3p, cc.coding_type);
        assert_eq!(4, cc.exon_rank);
        assert_eq!(5, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);

        // 3'UTR intronic
        pos = 22;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Intronic, cc.region_type);
        assert_eq!(CodingType::Utr3p, cc.coding_type);
        assert_eq!(5, cc.exon_rank);
        assert_eq!(6, cc.coding_base);
        assert_eq!(-2, cc.nearest_exon_distance);

        // 3'UTR exonic
        pos = 15;
        var = create_variant(pos);
        cc = determine_context(&var, trans_data_neg);

        assert_eq!(RegionType::Exonic, cc.region_type);
        assert_eq!(CodingType::Utr3p, cc.coding_type);
        assert_eq!(5, cc.exon_rank);
        assert_eq!(11, cc.coding_base);
        assert_eq!(0, cc.nearest_exon_distance);
    }

    #[test]
    fn test_upstream_promotor_variant()
    {
        let ref_genome = create_mock_genome(800);
        let ref_bases = ref_genome.bases().to_string();

        let exon_starts = [100, 200, 300, 400];
        let coding_start = Some(125);
        let coding_end = Some(425);

        let trans_data = create_trans_exons(crate::constants::PROMOTOR_UPSTREAM_GENE_IDS[0], TRANS_ID_1, NEG_STRAND, &exon_starts, 50, coding_start, coding_end, true, "");

        

        // upstream of coding end (start of transcript), within max coding distance range
        let mut pos = 475;
        let mut var = create_snv(pos, &ref_bases);

        let mut impact = classify(&var, trans_data, &ref_genome);
        assert_eq!(VariantEffect::UpstreamGene, impact.top_effect().unwrap());
        assert_eq!(CodingType::Enhancer, impact.coding.coding_type);
        assert_eq!(50, impact.coding.coding_base);
        assert_eq!("c.-50G>C", impact.coding.hgvs);

        pos = 726;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data, &ref_genome);
        assert_eq!(VariantEffect::UpstreamGene, impact.top_effect().unwrap());
        assert_eq!(CodingType::Unknown, impact.coding.coding_type);
        assert_eq!(0, impact.coding.coding_base);
        assert_eq!("", impact.coding.hgvs);

        // a common prod SNV
        pos = 549;
        var =  create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data, &ref_genome);
        assert_eq!(VariantEffect::UpstreamGene, impact.top_effect().unwrap());
        assert_eq!(CodingType::Enhancer, impact.coding.coding_type);
        assert_eq!(124, impact.coding.coding_base);
        assert_eq!("c.-124T>A", impact.coding.hgvs);

        // test again with an MNV which is a common TERT variant
        pos = 549;
        var =  VariantData::new(CHR_1, pos, &sub(&ref_bases, pos, pos + 2), "CC");

        impact = classify(&var, trans_data, &ref_genome);
        assert_eq!(VariantEffect::UpstreamGene, impact.top_effect().unwrap());
        assert_eq!(CodingType::Enhancer, impact.coding.coding_type);
        assert_eq!(125, impact.coding.coding_base);
        assert_eq!("c.-125_-124delATinsGG", impact.coding.hgvs);
    }

    #[test]
    fn test_non_coding_and_upstream_trans()
    {
        let ref_genome = create_mock_genome(500);
        let ref_bases = ref_genome.bases().to_string();

        let exon_starts = [100, 200, 300, 400];
        let coding_start = Some(125);
        let coding_end = Some(425);

        let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 50, coding_start, coding_end, false, "");

        

        // pre-gene
        let mut pos = 50;
        let mut var = create_snv(pos, &ref_bases);

        let mut impact = classify(&var, trans_data_pos_strand, &ref_genome);
        assert_eq!(VariantEffect::UpstreamGene, impact.top_effect().unwrap());

        // 5' UTR
        pos = 120;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_pos_strand, &ref_genome);
        assert_eq!(VariantEffect::FivePrimeUtr, impact.top_effect().unwrap());

        // 3' UTR
        pos = 440;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_pos_strand, &ref_genome);
        assert_eq!(VariantEffect::ThreePrimeUtr, impact.top_effect().unwrap());

        // intronic
        pos = 175;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_pos_strand, &ref_genome);
        assert_eq!(VariantEffect::Intronic, impact.top_effect().unwrap());

        let trans_data_neg_strand = create_trans_exons(GENE_ID_2, TRANS_ID_2, NEG_STRAND, &exon_starts, 50, coding_start, coding_end, false, "");

        // pre-gene
        pos = 490;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_neg_strand, &ref_genome);
        assert_eq!(VariantEffect::UpstreamGene, impact.top_effect().unwrap());

        // intronic
        pos = 375;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_neg_strand, &ref_genome);
        assert_eq!(VariantEffect::Intronic, impact.top_effect().unwrap());


        // non-coding exonic
        let trans_data_non_coding = create_trans_exons(GENE_ID_2, TRANS_ID_2, NEG_STRAND, &exon_starts, 50, None, None, false, "");

        // intronic
        pos = 160;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_non_coding, &ref_genome);
        assert_eq!(VariantEffect::Intronic, impact.top_effect().unwrap());

        // exonic has special classification
        pos = 220;
        var = create_snv(pos, &ref_bases);

        impact = classify(&var, trans_data_non_coding, &ref_genome);
        assert_eq!(VariantEffect::NonCodingTranscript, impact.top_effect().unwrap());

        // coding range for exonic
        pos = 220;
        let mut ref_ = sub(&ref_bases, pos, pos + 4);
        let mut alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_non_coding, &ref_genome);
        assert_eq!(130, impact.coding.coding_base);
        assert_eq!(220, impact.coding.coding_position_range[SE_START]);
        assert_eq!(223, impact.coding.coding_position_range[SE_END]);
        assert_eq!("n.130_133delGATCinsCGAT", impact.coding.hgvs);

        // non-coding spanning exonic boundaries
        pos = 298;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_non_coding, &ref_genome);
        assert_eq!(99, impact.coding.coding_base);
        assert_eq!(300, impact.coding.coding_position_range[SE_START]);
        assert_eq!(303, impact.coding.coding_position_range[SE_END]);
        assert_eq!("n.100_102+1delATCG", impact.coding.hgvs);

        pos = 348;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = sub(&ref_, 0, 1);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_non_coding, &ref_genome);
        assert_eq!(51, impact.coding.coding_base);
        assert_eq!(348, impact.coding.coding_position_range[SE_START]);
        assert_eq!(350, impact.coding.coding_position_range[SE_END]);
        assert_eq!("n.51-1_52delGAT", impact.coding.hgvs);
    }
