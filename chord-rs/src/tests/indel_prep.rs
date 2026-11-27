//! IndelPrepTest

use super::*;
use crate::config::PURPLE_SOMATIC_VCF_SUFFIX;

#[test]
fn can_prep_indels() {
    let mut c = config();
    c.snv_indel_vcf_file = Some(format!("{}{MINIMAL_SAMPLE}{PURPLE_SOMATIC_VCF_SUFFIX}", input_vcf_dir()));
    c.ref_genome_file = Some(dummy_genome_fasta());
    c.output_dir = Some(tmp_output_dir("can_prep_indels"));
    let actual = indel_counts(&c, MINIMAL_SAMPLE).unwrap();
    let expected = vec![
        mc("del.rep.len.1", 30),
        mc("del.rep.len.2", 3),
        mc("del.rep.len.3", 1),
        mc("del.rep.len.4", 1),
        mc("del.rep.len.5", 0),
        mc("ins.rep.len.1", 25),
        mc("ins.rep.len.2", 2),
        mc("ins.rep.len.3", 1),
        mc("ins.rep.len.4", 1),
        mc("ins.rep.len.5", 0),
        mc("del.mh.bimh.1", 8),
        mc("del.mh.bimh.2", 11),
        mc("del.mh.bimh.3", 0),
        mc("del.mh.bimh.4", 0),
        mc("del.mh.bimh.5", 1),
        mc("ins.mh.bimh.1", 9),
        mc("ins.mh.bimh.2", 12),
        mc("ins.mh.bimh.3", 2),
        mc("ins.mh.bimh.4", 0),
        mc("ins.mh.bimh.5", 0),
        mc("del.none.len.1", 75),
        mc("del.none.len.2", 43),
        mc("del.none.len.3", 26),
        mc("del.none.len.4", 7),
        mc("del.none.len.5", 8),
        mc("ins.none.len.1", 99),
        mc("ins.none.len.2", 55),
        mc("ins.none.len.3", 35),
        mc("ins.none.len.4", 14),
        mc("ins.none.len.5", 13),
    ];
    assert_eq!(actual, expected);
}

#[test]
fn providing_wrong_vcf_type_throws_error() {
    let mut c = config();
    c.sample_ids = vec![MINIMAL_SAMPLE.into()];
    c.snv_indel_vcf_file = Some(format!("{}{MINIMAL_SAMPLE}{PURPLE_SOMATIC_VCF_SUFFIX}", input_vcf_dir()));
    let e = indel_counts(&c, MINIMAL_SAMPLE).unwrap_err();
    assert!(format!("{e}").contains("IllegalStateException"));
}
