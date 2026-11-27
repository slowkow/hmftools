//! NEO's SvNeoEpitopeTest, ported

use super::*;
use crate::cat;
use crate::epitope::*;
use hmf_common::neo::{NeoEpitopeFusion, NeoEpitopeType};

#[test]
fn test_sv_fusions() {

    let mut ref_genome = MockRefGenome::default();

    let chr1_bases = generate_random_bases(120);

    ref_genome.put(CHR_1, &chr1_bases);

    // tests: intronic vs exonic, with and without insertions

    let exon_starts = [0, 20, 40, 60, 80, 100];
    let coding_start = Some(25);
    let coding_end = Some(85);

    let trans_data_pos_strand = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 10, coding_start, coding_end, false, "");

    let valid_trans = [trans_data_pos_strand.trans_name.clone(), trans_data_pos_strand.trans_name.clone()];

    // intronic
    let sv_data = create_fusion(35, ORIENT_FWD, 55, ORIENT_REV, "", &valid_trans);

    let mut ne_data = NeoEpitope::fusion(&sv_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(35, ne_data.position(FS_UP));
    assert_eq!(55, ne_data.position(FS_DOWN));
    assert_eq!(ORIENT_FWD, ne_data.orientation(FS_UP));
    assert_eq!(ORIENT_REV, ne_data.orientation(FS_DOWN));

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Intronic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Intronic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert_eq!(2, ne_data.exon_rank[FS_UP]);
    assert_eq!(4, ne_data.exon_rank[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 2).unwrap();

    let up_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));

    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);

    ne_data.set_amino_acids(&ref_genome, 3).unwrap();
    let up_wildtype_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 50));
    let up_wild_aas = get_amino_acids(&up_wildtype_bases, true);
    assert_eq!(up_wild_aas, ne_data.upstream_wild_type_acids);

    // intronic to exonic - skips to next exon
    let sv_data = create_fusion(35, ORIENT_FWD, 45, ORIENT_REV, "", &valid_trans);

    let mut ne_data = NeoEpitope::fusion(&sv_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Intronic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_0, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_2, ne_data.phases[FS_DOWN]);
    assert_eq!(2, ne_data.exon_rank[FS_UP]);
    assert_eq!(4, ne_data.exon_rank[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 2).unwrap();

    let up_bases = sub(&chr1_bases, 25, 31);
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 60, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);

    // exonic to exonic, in phase
    let sv_data = create_fusion(44, ORIENT_FWD, 63, ORIENT_REV, "", &valid_trans);

    let mut ne_data = NeoEpitope::fusion(&sv_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_2, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert_eq!(3, ne_data.exon_rank[FS_UP]);
    assert_eq!(4, ne_data.exon_rank[FS_DOWN]);
    assert!(ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 2).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 28, 31), sub(&chr1_bases, 40, 43));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 43, 45), sub(&chr1_bases, 63, 64));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    let down_bases = sub(&chr1_bases, 64, 70);
    assert_eq!(down_bases, ne_data.coding_bases[FS_DOWN]);

    // exonic to exonic with an insert sequence, out of phase
    let ins_sequence = "AA";

    let sv_data = create_fusion(44, ORIENT_FWD, 63, ORIENT_REV, ins_sequence, &valid_trans);

    let mut ne_data = NeoEpitope::fusion(&sv_data);

    ne_data.set_transcript_data(trans_data_pos_strand, trans_data_pos_strand);

    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_UP]);
    assert_eq!(CodingType::Coding, ne_data.coding_type[FS_DOWN]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_UP]);
    assert_eq!(RegionType::Exonic, ne_data.region_type[FS_DOWN]);
    assert_eq!(PHASE_1, ne_data.phases[FS_UP]);
    assert_eq!(PHASE_0, ne_data.phases[FS_DOWN]);
    assert_eq!(3, ne_data.exon_rank[FS_UP]);
    assert_eq!(4, ne_data.exon_rank[FS_DOWN]);
    assert!(!ne_data.phase_matched());

    ne_data.set_coding_bases(&ref_genome, 3).unwrap();

    let up_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 43));
    assert_eq!(up_bases, ne_data.coding_bases[FS_UP]);

    let novel_bases = cat!(sub(&chr1_bases, 43, 45), ins_sequence, sub(&chr1_bases, 63, 71), sub(&chr1_bases, 80, 91), sub(&chr1_bases, 100, 111));
    assert_eq!(novel_bases, ne_data.novel_codon_bases);

    assert_eq!("", ne_data.coding_bases[FS_DOWN]);

    ne_data.set_amino_acids(&ref_genome, 3).unwrap();
    let up_wildtype_bases = cat!(sub(&chr1_bases, 25, 31), sub(&chr1_bases, 40, 51), sub(&chr1_bases, 60, 61));
    let up_wild_aas = get_amino_acids(&up_wildtype_bases, true);
    assert_eq!(up_wild_aas, ne_data.upstream_wild_type_acids);

    // intronic to upstream, skips to first splice acceptor

}

#[test]
fn test_non_disruptive_sv() {

    let exon_starts = [100, 300, 500, 700, 900];
    let coding_start = Some(150);
    let coding_end = Some(750);

    let trans_data1 = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &exon_starts, 100, coding_start, coding_end, false, "");

    let valid_trans = [trans_data1.trans_name.clone(), trans_data1.trans_name.clone()];

    // intronic DEL
    let sv_positions = [420, 480];

    assert!(sv_is_non_disruptive_in_coding_transcript(&sv_positions, trans_data1));

    let sv_positions = [320, 480];  // exonic
    assert!(!sv_is_non_disruptive_in_coding_transcript(&sv_positions, trans_data1));

    let sv_positions = [420, 680];  // deletes an exon
    assert!(!sv_is_non_disruptive_in_coding_transcript(&sv_positions, trans_data1));

    let sv_positions = [820, 880];  // UTR region
    assert!(!sv_is_non_disruptive_in_coding_transcript(&sv_positions, trans_data1));
}


/// NE14: Java's immune prefixes repeat IG_J in place of IG_V
#[test]
fn ne14_immune_transcripts() {
    let mut t = create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[0, 20], 10, Some(5), Some(25), false, "IG_V_gene").clone();
    assert!(!crate::compat::with_java(true, || is_immune_transcript(&t)));
    assert!(crate::compat::with_java(false, || is_immune_transcript(&t)));
    for (bio, skip) in [("TR_V_gene", true), ("IG_J_gene", true), ("IG_D_gene", true), ("IG_C_gene", false), ("protein_coding", false)] {
        t.bio_type = bio.to_string();
        for java in [true, false] {
            assert_eq!(crate::compat::with_java(java, || is_immune_transcript(&t)), skip, "{bio}");
        }
    }
}
