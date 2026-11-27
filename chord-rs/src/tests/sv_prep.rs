//! SvPrepTest

use super::*;
use crate::config::{PURPLE_SOMATIC_VCF_SUFFIX, PURPLE_SV_VCF_SUFFIX};

#[test]
fn can_prep_structural_variants() {
    let mut c = config();
    c.sv_vcf_file = Some(format!("{}{MINIMAL_SAMPLE}{PURPLE_SV_VCF_SUFFIX}", input_vcf_dir()));
    c.output_dir = Some(tmp_output_dir("can_prep_svs"));
    let actual = sv_counts(&c, MINIMAL_SAMPLE).unwrap();
    let expected = vec![
        mc("DEL_0e00_1e03_bp", 1),
        mc("DEL_1e03_1e04_bp", 0),
        mc("DEL_1e04_1e05_bp", 1),
        mc("DEL_1e05_1e06_bp", 0),
        mc("DEL_1e06_1e07_bp", 0),
        mc("DEL_1e07_Inf_bp", 0),
        mc("DUP_0e00_1e03_bp", 2),
        mc("DUP_1e03_1e04_bp", 0),
        mc("DUP_1e04_1e05_bp", 0),
        mc("DUP_1e05_1e06_bp", 0),
        mc("DUP_1e06_1e07_bp", 0),
        mc("DUP_1e07_Inf_bp", 0),
        mc("INV_0e00_1e03_bp", 0),
        mc("INV_1e03_1e04_bp", 0),
        mc("INV_1e04_1e05_bp", 0),
        mc("INV_1e05_1e06_bp", 0),
        mc("INV_1e06_1e07_bp", 0),
        mc("INV_1e07_Inf_bp", 1),
        mc("TRA", 1),
    ];
    assert_eq!(actual, expected);
}

#[test]
fn can_prep_structural_variants_from_empty_vcf() {
    let mut c = config();
    c.sample_ids = vec![EMPTY_SAMPLE.into()];
    c.purple_dir = Some(input_vcf_dir());
    c.output_dir = Some(tmp_output_dir("empty_svs"));
    let counts = sv_counts(&c, EMPTY_SAMPLE).unwrap();
    assert_eq!(counts.iter().map(|m| m.count).sum::<i32>(), 0);
}

#[test]
fn providing_wrong_vcf_type_throws_error() {
    let mut c = config();
    c.sample_ids = vec![MINIMAL_SAMPLE.into()];
    c.sv_vcf_file = Some(format!("{}{MINIMAL_SAMPLE}{PURPLE_SOMATIC_VCF_SUFFIX}", input_vcf_dir()));
    let e = sv_counts(&c, MINIMAL_SAMPLE).unwrap_err();
    assert!(format!("{e}").contains("IllegalStateException"));
}
