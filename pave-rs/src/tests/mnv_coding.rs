//! MnvCodingContextTest

use super::*;
use crate::impact::effects::{CodingType, RegionType, VariantEffect};
use crate::impact::{SE_END, SE_START};

    #[test]
    fn test_non_exonic_positions()
    {
        let ref_genome = create_mock_genome(150);
        let ref_bases = ref_genome.bases().to_string();

        let exon_starts = [0, 30, 60, 90];

        let trans_data_pos = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, Some(35), Some(65), false, "");

        

        // test 5'UTR
        let mut pos = 22;
        let mut ref_ = sub(&ref_bases, pos, pos + 4);
        let mut alt = generate_alt(&ref_);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos, &ref_genome);

        // first check general coding context fields
        assert_eq!(5, impact.coding.coding_base); // since closer
        assert_eq!(CodingType::Utr5p, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(2, impact.coding.exon_rank);
        assert_eq!(-8, impact.coding.nearest_exon_distance);

        pos = 15;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        // first check general coding context fields
        assert_eq!(6, impact.coding.coding_base); // since closer
        assert_eq!(1, impact.coding.exon_rank);
        assert_eq!(5, impact.coding.nearest_exon_distance);

        // coding
        pos = 45;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(6, impact.coding.coding_base); // since closer to previous
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(2, impact.coding.exon_rank);
        assert_eq!(5, impact.coding.nearest_exon_distance);

        pos = 51;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(7, impact.coding.coding_base); // since closer to previous
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(3, impact.coding.exon_rank);
        assert_eq!(-9, impact.coding.nearest_exon_distance);

        // 3'UTR
        pos = 72;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(5, impact.coding.coding_base); // since closer to previous
        assert_eq!(CodingType::Utr3p, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(3, impact.coding.exon_rank);
        assert_eq!(2, impact.coding.nearest_exon_distance);

        pos = 81;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(6, impact.coding.coding_base); // since closer to previous
        assert_eq!(4, impact.coding.exon_rank);
        assert_eq!(-9, impact.coding.nearest_exon_distance);

        // repeat for negative strand
        let trans_data_neg = create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &exon_starts, 10, Some(35), Some(65), false, "");

        // test 5'UTR
        pos = 85;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(6, impact.coding.coding_base); // since closert to upstream exon
        assert_eq!(CodingType::Utr5p, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(1, impact.coding.exon_rank);
        assert_eq!(2, impact.coding.nearest_exon_distance);

        pos = 73;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(5, impact.coding.coding_base);
        assert_eq!(2, impact.coding.exon_rank);
        assert_eq!(-6, impact.coding.nearest_exon_distance);

        // coding
        pos = 51;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(6, impact.coding.coding_base);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(2, impact.coding.exon_rank);
        assert_eq!(6, impact.coding.nearest_exon_distance);

        pos = 45;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(7, impact.coding.coding_base); // since closer to previous
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(3, impact.coding.exon_rank);
        assert_eq!(-8, impact.coding.nearest_exon_distance);

        // 3'UTR
        pos = 22;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(5, impact.coding.coding_base); // since closer to previous
        assert_eq!(CodingType::Utr3p, impact.coding.coding_type);
        assert_eq!(RegionType::Intronic, impact.coding.region_type);
        assert_eq!(3, impact.coding.exon_rank);
        assert_eq!(5, impact.coding.nearest_exon_distance);

        pos = 15;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(6, impact.coding.coding_base); // since closer to previous
        assert_eq!(4, impact.coding.exon_rank);
        assert_eq!(-8, impact.coding.nearest_exon_distance);
    }

    #[test]
    fn test_mnv_basic()
    {
        let ref_genome = create_mock_genome(150);
        let ref_bases = ref_genome.bases().to_string();

        let trans_data_pos = create_pos_transcript();

        

        // MNV spanning 2 codons at all positions - codons 30-32 and 33-35
        let mut pos = 30;
        let mut ref_codon_bases = sub(&ref_bases, pos, pos + 6);

        let mut ref_ = sub(&ref_codon_bases, 0, 4);
        let mut alt = generate_alt(&ref_);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos, &ref_genome);

        // first check general coding context fields
        assert_eq!(7, impact.coding.coding_base);
        assert_eq!(30, impact.coding.coding_position_range[SE_START]);
        assert_eq!(33, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_1, impact.coding.upstream_phase);
        assert_eq!(2, impact.coding.exon_rank);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Exonic, impact.coding.region_type);

        // then the protein context
        assert!(impact.has_protein_context());
        assert_eq!(3, impact.pc().codon_index);
        assert_eq!("p.Ser3_Ile4delinsArgPhe", impact.pc().hgvs);

        let mut alt_codon_bases = cat!(alt, sub(&ref_codon_bases, 4, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // shift up one base
        pos = 31;
        ref_ = sub(&ref_codon_bases, 1, 5);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);
        impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(8, impact.coding.coding_base);
        assert_eq!(31, impact.coding.coding_position_range[SE_START]);
        assert_eq!(34, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_2, impact.coding.upstream_phase);
        assert_eq!(3, impact.pc().codon_index);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 1), alt, sub(&ref_codon_bases, 5, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        pos = 32;
        ref_ = sub(&ref_codon_bases, 2, 6);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);
        impact = classify(&var, trans_data_pos, &ref_genome);

        assert_eq!(9, impact.coding.coding_base);
        assert_eq!(32, impact.coding.coding_position_range[SE_START]);
        assert_eq!(35, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(3, impact.pc().codon_index);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 2), alt);
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // again but on negative strand
        let trans_data_neg = create_neg_transcript();

        // MNV spanning 2 codons at all positions - codons 80-78 and 77-75
        pos = 75;
        ref_codon_bases = sub(&ref_bases, pos, pos + 6);

        ref_ = sub(&ref_codon_bases, 2, 6);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, 77, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(7, impact.coding.coding_base);
        assert_eq!(77, impact.coding.coding_position_range[SE_START]);
        assert_eq!(80, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_1, impact.coding.upstream_phase);
        assert_eq!(4, impact.coding.exon_rank);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Exonic, impact.coding.region_type);

        // then the protein context
        assert!(impact.has_protein_context());
        assert_eq!(3, impact.pc().codon_index);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 2), alt);
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // shift along 1 each time
        ref_ = sub(&ref_codon_bases, 1, 5);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, 76, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(8, impact.coding.coding_base);
        assert_eq!(76, impact.coding.coding_position_range[SE_START]);
        assert_eq!(79, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_2, impact.coding.upstream_phase);
        assert_eq!(3, impact.pc().codon_index);

        alt_codon_bases = cat!(sub(&ref_codon_bases, 0, 1), alt, sub(&ref_codon_bases, 5, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // and again
        ref_ = sub(&ref_codon_bases, 0, 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, 75, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(9, impact.coding.coding_base);
        assert_eq!(75, impact.coding.coding_position_range[SE_START]);
        assert_eq!(78, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(3, impact.pc().codon_index);

        alt_codon_bases = cat!(alt, sub(&ref_codon_bases, 4, 6));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // crossing coding start
        pos = 94;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);
        assert_eq!(1, impact.coding.coding_base);
        assert!(impact.coding.spans_coding_start);
        assert_eq!(94, impact.coding.coding_position_range[SE_START]);
        assert_eq!(95, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-2_2delTCGAinsATCG", impact.coding.hgvs);
    }

    #[test]
    fn test_mnv_across_splice()
    {
        let ref_genome = create_mock_genome(150);

        

        // MNV spanning exon boundary - shortens the effect
        let ref_bases = ref_genome.bases().to_string();

        let trans_data_pos = create_pos_transcript();

        // firstly a 2-lots across the exon
        let mut ref_codon_bases = sub(&ref_bases, 15, 21);

        let mut pos = 20;
        let mut ref_ = sub(&ref_bases, pos, pos + 2);
        let mut alt = generate_alt(&ref_);
        let mut var = VariantData::new(CHR_1, pos, &ref_, &alt);

        let mut impact = classify(&var, trans_data_pos, &ref_genome);

        // first check general coding context fields
        assert_eq!(6, impact.coding.coding_base);
        assert_eq!(20, impact.coding.coding_position_range[SE_START]);
        assert_eq!(20, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(1, impact.coding.exon_rank);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert!(impact.coding.spans_splice_junction);
        assert_eq!(RegionType::Exonic, impact.coding.region_type);

        assert!(impact.has_protein_context());
        assert_eq!(2, impact.pc().codon_index);
        assert_eq!(sub(&ref_codon_bases, 3, 6), impact.pc().ref_codon_bases);

        let mut alt_codon_bases = cat!(sub(&ref_codon_bases, 3, 5), sub(&alt, 0, 1));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // same again but spanning 2 exons with the codons: 36-38, 39-50, 51-53
        ref_codon_bases = cat!(sub(&ref_bases, 36, 41), sub(&ref_bases, 50, 54));

        pos = 38;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        // first check general coding context fields
        assert_eq!(15, impact.coding.coding_base);
        assert_eq!(38, impact.coding.coding_position_range[SE_START]);
        assert_eq!(40, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(2, impact.coding.exon_rank);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Exonic, impact.coding.region_type);
        assert_eq!("c.15_17+1delTCGAinsCGAT", impact.coding.hgvs);
        assert!(impact.has_protein_context());
        assert_eq!(5, impact.pc().codon_index);
        assert_eq!(cat!(sub(&ref_bases, 36, 41), sub(&ref_bases, 50, 51)), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_bases, 36, 38), sub(&alt, 0, 3), sub(&ref_bases, 50, 51));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // and starting before the next exon
        pos = 48;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_pos, &ref_genome);

        // first check general coding context fields
        assert_eq!(18, impact.coding.coding_base);
        assert_eq!(50, impact.coding.coding_position_range[SE_START]);
        assert_eq!(51, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(3, impact.coding.exon_rank);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Exonic, impact.coding.region_type);
        assert_eq!("c.18-2_19delGATCinsATCG", impact.coding.hgvs);

        assert!(impact.has_protein_context());
        assert_eq!(6, impact.pc().codon_index);
        assert_eq!(cat!(sub(&ref_bases, 39, 41), sub(&ref_bases, 50, 54)), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_bases, 39, 41), sub(&alt, 2, 4), sub(&ref_bases, 52, 54));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // for negative strand, is the phase set correctly even though the variant starts in the upstream intron?

        let trans_data_neg = create_neg_transcript();

        // and starting at the end of the last exon
        pos = 88;
        ref_ = sub(&ref_bases, pos, pos + 3);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(6, impact.coding.coding_base);
        assert_eq!(90, impact.coding.coding_position_range[SE_START]);
        assert_eq!(90, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(3, impact.coding.exon_rank);
        assert!(impact.coding.spans_splice_junction);
        assert_eq!(CodingType::Coding, impact.coding.coding_type);
        assert_eq!(RegionType::Exonic, impact.coding.region_type);

        assert!(impact.has_protein_context());
        assert_eq!(2, impact.pc().codon_index);
        assert_eq!(sub(&ref_bases, 90, 93), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&alt, 2, 3), sub(&ref_bases, 91, 93));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // spanning exons 2 and 3 codons: 74-72, 71-60, 59-57
        pos = 68;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        // first check general coding context fields
        assert_eq!(16, impact.coding.coding_base);
        assert_eq!(70, impact.coding.coding_position_range[SE_START]);
        assert_eq!(71, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_1, impact.coding.upstream_phase);
        assert_eq!(4, impact.coding.exon_rank);
        assert!(impact.coding.spans_splice_junction);

        assert!(impact.has_protein_context());
        assert_eq!(6, impact.pc().codon_index);
        assert_eq!(cat!(sub(&ref_bases, 60, 61), sub(&ref_bases, 70, 72)), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_bases, 60, 61), sub(&alt, 2, 4));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        pos = 59;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(18, impact.coding.coding_base);
        assert_eq!(59, impact.coding.coding_position_range[SE_START]);
        assert_eq!(60, impact.coding.coding_position_range[SE_END]);
        assert_eq!(PHASE_0, impact.coding.upstream_phase);
        assert_eq!(5, impact.coding.exon_rank);
        assert!(impact.coding.spans_splice_junction);

        assert!(impact.has_protein_context());
        assert_eq!(6, impact.pc().codon_index);
        assert_eq!(cat!(sub(&ref_bases, 57, 61), sub(&ref_bases, 70, 72)), impact.pc().ref_codon_bases);

        alt_codon_bases = cat!(sub(&ref_bases, 57, 59), sub(&alt, 0, 2), sub(&ref_bases, 70, 72));
        assert_eq!(alt_codon_bases, impact.pc().alt_codon_bases);

        // coding base is the lowest numerically (ie 19 is at pos 59, not the coding base at the lower position (ie 56)
        pos = 56;
        ref_ = sub(&ref_bases, pos, pos + 4);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(19, impact.coding.coding_base);
        assert_eq!(56, impact.coding.coding_position_range[SE_START]);
        assert_eq!(59, impact.coding.coding_position_range[SE_END]);

        // 5'UTR exonic
        pos = 112;
        ref_ = sub(&ref_bases, pos, pos + 3);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(10, impact.coding.coding_base); // being the further exonic base from start of coding
        assert_eq!(112, impact.coding.coding_position_range[SE_START]);
        assert_eq!(114, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-10_-8delATCinsGAT", impact.coding.hgvs);

        // spanning into next lower intron
        pos = 108;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(8, impact.coding.coding_base); // being the further exonic base from start of coding
        assert_eq!(110, impact.coding.coding_position_range[SE_START]);
        assert_eq!(112, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-8_-6+2delCGATCinsTCGAT", impact.coding.hgvs);

        // next upper intron
        pos = 118;
        ref_ = sub(&ref_bases, pos, pos + 5);
        alt = generate_alt(&ref_);
        var = VariantData::new(CHR_1, pos, &ref_, &alt);

        impact = classify(&var, trans_data_neg, &ref_genome);

        assert_eq!(16, impact.coding.coding_base);
        assert_eq!(118, impact.coding.coding_position_range[SE_START]);
        assert_eq!(120, impact.coding.coding_position_range[SE_END]);
        assert_eq!("c.-16-2_-14delATCGAinsGATCG", impact.coding.hgvs);

    }
