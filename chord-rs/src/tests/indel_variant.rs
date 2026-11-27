//! IndelVariantTest

use super::*;
use crate::indel::{count_homologous_bases, count_repeat_units, IndelVariant};
use crate::refgenome::RefGenomeSource;
use crate::vcf_file::SmallVariant;

fn ref_genome() -> RefGenomeSource {
    RefGenomeSource::load(&testdata("fasta/indel_regions.fasta")).unwrap()
}

#[test]
fn can_get_deletion_flank_sequences() {
    // based on a deletion in COLO829v003T
    let deletion = IndelVariant::new(SmallVariant::new("1", 21, "CAA", "C"));
    let g = ref_genome();
    let left = deletion.flanking_ref_bases_left(&g, 1).unwrap();
    let right = deletion.flanking_ref_bases_right(&g, 3).unwrap();
    assert_eq!(deletion.indel_length, 2);
    assert_eq!(left, "CC");
    assert_eq!(right, "AAGTCT");
}

#[test]
fn can_get_insertion_flank_sequences() {
    // based on an insertion in COLO829v003T
    let insertion = IndelVariant::new(SmallVariant::new("17", 21, "G", "GGGTGT"));
    let g = ref_genome();
    let left = insertion.flanking_ref_bases_left(&g, 1).unwrap();
    let right = insertion.flanking_ref_bases_right(&g, 3).unwrap();
    assert_eq!(insertion.indel_length, 5);
    assert_eq!(left, "GGTCG");
    assert_eq!(right, "GGTGTGGTGGCTCAC");
}

#[test]
fn can_count_homologous_bases() {
    let indel = "AAGGCACTGGAAAGAATATAGAGA";
    let flank = "AAGGCACTCTCCACTCTGCTTTGTACATTTCTCATTACTGTATTAGCATCTTATAATCATAAATAACACCTT";
    assert_eq!(count_homologous_bases(indel, flank), 8);
}

#[test]
fn can_count_repeat_units() {
    assert_eq!(count_repeat_units("CT", "CTCTCAAA"), 3);
}

#[test]
fn cant_get_indel_sequence() {
    let deletion = IndelVariant::new(SmallVariant::new("chr1", 1, "AGGG", "A"));
    assert_eq!(deletion.indel_sequence(), "GGG");
    let insertion = IndelVariant::new(SmallVariant::new("chr1", 1, "C", "CTTT"));
    assert_eq!(insertion.indel_sequence(), "TTT");
}
