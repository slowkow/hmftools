//! SnvPrepTest

use super::*;
use crate::config::{PURPLE_SOMATIC_VCF_SUFFIX, PURPLE_SV_VCF_SUFFIX};

#[test]
fn can_prep_snvs() {
    let mut c = config();
    c.snv_indel_vcf_file = Some(format!("{}{MINIMAL_SAMPLE}{PURPLE_SOMATIC_VCF_SUFFIX}", input_vcf_dir()));
    c.ref_genome_file = Some(dummy_genome_fasta());
    c.output_dir = Some(tmp_output_dir("can_prep_snvs"));
    let actual = snv_counts(&c, MINIMAL_SAMPLE).unwrap();
    // there are 96 contexts in total, but only test the first few
    let expected = vec![
        mc("A[C>A]A", 10),
        mc("A[C>A]C", 6),
        mc("A[C>A]G", 8),
        mc("A[C>A]T", 2),
        mc("C[C>A]A", 4),
        mc("C[C>A]C", 4),
        mc("C[C>A]G", 6),
        mc("C[C>A]T", 5),
        mc("G[C>A]A", 4),
        mc("G[C>A]C", 4),
        mc("G[C>A]G", 9),
        mc("G[C>A]T", 5),
        mc("T[C>A]A", 6),
        mc("T[C>A]C", 4),
        mc("T[C>A]G", 0),
        mc("T[C>A]T", 6),
        mc("A[C>G]A", 8),
        mc("A[C>G]C", 1),
        mc("A[C>G]G", 8),
        mc("A[C>G]T", 11),
    ];
    assert_eq!(actual.len(), 96);
    assert_eq!(actual[..expected.len()], expected[..]);
}

#[test]
fn can_prep_snvs_with_non_standard_nucleotides() {
    let mut c = config();
    c.snv_indel_vcf_file = Some(format!("{}{NON_STANDARD_NUC_GENOME_SAMPLE}{PURPLE_SOMATIC_VCF_SUFFIX}", input_vcf_dir()));
    c.ref_genome_file = Some(non_standard_nuc_genome_fasta());
    c.output_dir = Some(tmp_output_dir("non_standard"));
    let actual = snv_counts(&c, NON_STANDARD_NUC_GENOME_SAMPLE).unwrap();
    assert!(actual.iter().map(|m| m.count).sum::<i32>() > 0);
}

/// Java expects an IllegalStateException; it comes from the missing reference genome (SnvPrep's
/// constructor), before the SV VCF is read
#[test]
fn providing_wrong_vcf_type_throws_error() {
    let mut c = config();
    c.sample_ids = vec![MINIMAL_SAMPLE.into()];
    c.snv_indel_vcf_file = Some(format!("{}{MINIMAL_SAMPLE}{PURPLE_SV_VCF_SUFFIX}", input_vcf_dir()));
    let e = snv_counts(&c, MINIMAL_SAMPLE).unwrap_err();
    assert!(format!("{e}").contains("IllegalStateException"));
}
