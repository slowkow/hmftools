//! InsertCodingContextTest

use super::*;
use crate::impact::effects::{CodingType, RegionType, VariantEffect};
use crate::impact::hgvs_protein::HGVS_SPLICE_UNKNOWN;
use crate::impact::{SE_END, SE_START};

#[test]
fn test_inserts() {
    let ref_genome = create_mock_genome(150);
    let rb = ref_genome.bases().to_string();
    let trans_pos = create_pos_transcript();

    // an insert at the last base of an exon does nothing
    let r = sub(&rb, 20, 21);
    let alt = r.clone() + "AAA";
    let var = VariantData::new(CHR_1, 20, &r, &alt);
    let impact = classify(&var, trans_pos, &ref_genome);
    assert_eq!(CodingType::Coding, impact.coding.coding_type);
    assert_eq!(RegionType::Intronic, impact.coding.region_type);
    assert_eq!(6, impact.coding.coding_base);
    assert_eq!(HGVS_SPLICE_UNKNOWN, impact.protein.as_ref().unwrap().hgvs);

    // frameshift
    let r = sub(&rb, 19, 20);
    let alt = r.clone() + "AAAA";
    let var = VariantData::new(CHR_1, 19, &r, &alt);
    let impact = classify(&var, trans_pos, &ref_genome);
    assert_eq!(5, impact.coding.coding_base);
    assert_eq!(19, impact.coding.coding_position_range[SE_START]);
    assert_eq!(19, impact.coding.coding_position_range[SE_END]);
    assert!(impact.coding.is_frameshift);

    // inframe insert
    let r = sub(&rb, 19, 20);
    let alt = r.clone() + "AAA";
    let var = VariantData::new(CHR_1, 19, &r, &alt);
    let impact = classify(&var, trans_pos, &ref_genome);
    assert_eq!(5, impact.coding.coding_base);
    assert_eq!(19, impact.coding.coding_position_range[SE_START]);
    assert_eq!(19, impact.coding.coding_position_range[SE_END]);
    assert_eq!(PHASE_2, impact.coding.upstream_phase);
    assert_eq!(1, impact.coding.exon_rank);
    assert_eq!(CodingType::Coding, impact.coding.coding_type);
    assert_eq!(RegionType::Exonic, impact.coding.region_type);
    assert!(!impact.coding.is_frameshift);
    let pc = impact.protein.as_ref().unwrap();
    assert_eq!(sub(&rb, 15, 21) + &sub(&rb, 30, 33), pc.ref_codon_bases);
    let alt_codon_bases = sub(&rb, 15, 19) + &alt + &sub(&rb, 20, 21) + &sub(&rb, 30, 33);
    assert_eq!(alt_codon_bases, pc.alt_codon_bases);
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());

    // inframe insert at the exon start
    let r = sub(&rb, 50, 51);
    let alt = r.clone() + "AAA";
    let var = VariantData::new(CHR_1, 50, &r, &alt);
    let impact = classify(&var, trans_pos, &ref_genome);
    assert_eq!(18, impact.coding.coding_base);
    assert_eq!(50, impact.coding.coding_position_range[SE_START]);
    assert_eq!(50, impact.coding.coding_position_range[SE_END]);
    assert_eq!(3, impact.coding.exon_rank);
    assert!(!impact.coding.is_frameshift);
    let pc = impact.protein.as_ref().unwrap();
    assert_eq!(sub(&rb, 39, 41) + &sub(&rb, 50, 54), pc.ref_codon_bases);
    let alt_codon_bases = sub(&rb, 39, 41) + &alt + &sub(&rb, 51, 54);
    assert_eq!(alt_codon_bases, pc.alt_codon_bases);
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());

    // an insert just before the exon start
    let r = sub(&rb, 29, 30);
    let alt = r.clone() + "A";
    let var = VariantData::new(CHR_1, 29, &r, &alt);
    let impact = classify(&var, trans_pos, &ref_genome);
    assert_eq!(7, impact.coding.coding_base);
    assert_eq!(30, impact.coding.coding_position_range[SE_START]);
    assert_eq!(30, impact.coding.coding_position_range[SE_END]);
    assert_eq!(2, impact.coding.exon_rank);
    assert!(impact.coding.is_frameshift);
    assert!(impact.in_splice_region());
    assert_eq!(Some(VariantEffect::Frameshift), impact.top_effect());

    // the negative strand
    let trans_neg = create_neg_transcript();

    // at the exon boundary: the inserted bases go into the coding region
    let r = sub(&rb, 80, 81);
    let alt = r.clone() + "CGG";
    let var = VariantData::new(CHR_1, 80, &r, &alt);
    let impact = classify(&var, trans_neg, &ref_genome);
    assert_eq!(CodingType::Coding, impact.coding.coding_type);
    assert_eq!(RegionType::Exonic, impact.coding.region_type);
    assert_eq!(7, impact.coding.coding_base);
    assert_eq!(80, impact.coding.coding_position_range[SE_START]);
    assert_eq!(80, impact.coding.coding_position_range[SE_END]);
    assert_eq!(4, impact.coding.exon_rank);
    assert!(impact.in_splice_region());
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());

    // frameshift
    let r = sub(&rb, 91, 92);
    let alt = r.clone() + "AAAA";
    let var = VariantData::new(CHR_1, 91, &r, &alt);
    let impact = classify(&var, trans_neg, &ref_genome);
    assert_eq!(4, impact.coding.coding_base);
    assert_eq!(92, impact.coding.coding_position_range[SE_START]);
    assert_eq!(92, impact.coding.coding_position_range[SE_END]);
    assert!(impact.coding.is_frameshift);

    // inframe insert
    let r = sub(&rb, 91, 92);
    let alt = r.clone() + "GGG";
    let var = VariantData::new(CHR_1, 91, &r, &alt);
    let impact = classify(&var, trans_neg, &ref_genome);
    assert_eq!(4, impact.coding.coding_base);
    assert_eq!(92, impact.coding.coding_position_range[SE_START]);
    assert_eq!(92, impact.coding.coding_position_range[SE_END]);
    assert_eq!(PHASE_1, impact.coding.upstream_phase);
    assert_eq!(3, impact.coding.exon_rank);
    assert_eq!(CodingType::Coding, impact.coding.coding_type);
    assert_eq!(RegionType::Exonic, impact.coding.region_type);
    assert!(!impact.coding.is_frameshift);
    let pc = impact.protein.as_ref().unwrap();
    assert_eq!(sub(&rb, 78, 81) + &sub(&rb, 90, 96), pc.ref_codon_bases);
    let alt_codon_bases = sub(&rb, 78, 81) + &sub(&rb, 90, 91) + &alt + &sub(&rb, 92, 96);
    assert_eq!(alt_codon_bases, pc.alt_codon_bases);
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());

    // inframe insert at the exon end
    let r = sub(&rb, 70, 71);
    let alt = r.clone() + "AAA";
    let var = VariantData::new(CHR_1, 70, &r, &alt);
    let impact = classify(&var, trans_neg, &ref_genome);
    assert_eq!(16, impact.coding.coding_base);
    assert_eq!(71, impact.coding.coding_position_range[SE_START]);
    assert_eq!(71, impact.coding.coding_position_range[SE_END]);
    assert_eq!(4, impact.coding.exon_rank);
    assert!(!impact.coding.is_frameshift);
    let pc = impact.protein.as_ref().unwrap();
    assert_eq!(sub(&rb, 57, 61) + &sub(&rb, 70, 75), pc.ref_codon_bases);
    let alt_codon_bases = sub(&rb, 57, 61) + &alt + &sub(&rb, 71, 75);
    assert_eq!(alt_codon_bases, pc.alt_codon_bases);
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());
    assert_eq!("p.Asp6delinsValTyr", pc.hgvs);

    // at phase 0
    let r = sub(&rb, 92, 93);
    let alt = r.clone() + "GGG";
    let var = VariantData::new(CHR_1, 92, &r, &alt);
    let impact = classify(&var, trans_neg, &ref_genome);
    assert_eq!(3, impact.coding.coding_base);
    assert_eq!(93, impact.coding.coding_position_range[SE_START]);
    assert_eq!(93, impact.coding.coding_position_range[SE_END]);
    assert_eq!(PHASE_0, impact.coding.upstream_phase);
    assert_eq!(3, impact.coding.exon_rank);
    assert!(!impact.coding.is_frameshift);
    let pc = impact.protein.as_ref().unwrap();
    assert_eq!(sub(&rb, 90, 96), pc.ref_codon_bases);
    let alt_codon_bases = sub(&rb, 90, 93) + &alt[1..] + &sub(&rb, 93, 96);
    assert_eq!(alt_codon_bases, pc.alt_codon_bases);
    assert_eq!(Some(VariantEffect::InframeInsertion), impact.top_effect());
}
