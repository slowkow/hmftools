//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/dedup/IndelDedupTest.java
//!
//! Java's IndelDeduper reads the reference through MockRefGenome. The port's RefGenome is an indexed
//! FASTA, so the tests write the mock chromosome to a temporary FASTA file. MockRefGenome is zero-based
//! (position p is string index p), so the FASTA holds the Java string without its first character.

use super::util::*;
use super::variant_dedup::{
    add_local_phase_set, create_read_context, create_sage_variant_from_context, create_simple_variant, set_tumor_quality,
};
use crate::constants::{DEFAULT_FLANK_LENGTH, MIN_CORE_DISTANCE};
use crate::dedup::{build_alt_bases, dedup_indels};
use crate::filters::{SageVariant, SoftFilter};
use crate::refseq::RefGenome;
use crate::types::*;
use std::path::PathBuf;

//                                       10        20        30        40        50        60        70
//                             0123456789012345678901234567890123456789012345678901234567890123456789012
const TEST_REF_BASES: &str = "XAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGTAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGT";

/// A directory for temporary test files, next to the test binary (inside the cargo target directory).
pub(super) fn test_tmp_dir() -> PathBuf {
    let exe = std::env::current_exe().unwrap();
    let dir = exe.parent().unwrap().join("sage_rs_test_tmp");
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Stands in for MockRefGenome: `chromosomes` are (name, Java RefGenomeMap string). The first character
/// of each string is dropped so that 1-based FASTA positions equal Java string indices.
pub(super) fn mock_ref_genome(test_name: &str, chromosomes: &[(&str, &str)]) -> RefGenome {
    let path = test_tmp_dir().join(format!("{test_name}.fasta"));
    let _ = std::fs::remove_file(path.with_extension("fasta.fai"));
    let mut fasta = String::new();
    for (name, bases) in chromosomes {
        fasta.push_str(&format!(">{}\n{}\n", name, &bases[1..]));
    }
    std::fs::write(&path, fasta).unwrap();
    RefGenome::open(path.to_str().unwrap()).unwrap()
}

/// IndelDedupTest.createIndelTestVariant
fn create_indel_test_variant(
    position: i32, index: usize, var_read_bases: &str, ref_: &str, alt: &str, local_phase_set: i32,
) -> SageVariant {
    let variant = create_simple_variant(position, ref_, alt);

    let min_core = MIN_CORE_DISTANCE as usize;
    let flank = DEFAULT_FLANK_LENGTH as usize;

    // for these tests, the core is +/- 2 around the index no matter what the indel or variant
    let left_flank = &var_read_bases[index - min_core - flank..index - 2];

    let right_core_start = index + alt.len();
    let core_end_index = right_core_start + min_core;
    let right_flank = &var_read_bases[core_end_index..core_end_index + flank];

    let left_core = &var_read_bases[index - min_core..index];
    let right_core = &var_read_bases[right_core_start..right_core_start + min_core];

    let read_context = create_read_context(variant, left_core, right_core, left_flank, right_flank);

    let mut sage_variant = create_sage_variant_from_context(read_context);

    add_local_phase_set(&mut sage_variant, local_phase_set, 1);
    set_tumor_quality(&mut sage_variant, 5, 1000);

    sage_variant
}

// IndelDedupTest.testVariantAltBuilding
#[test]
fn test_variant_alt_building() {
    //               100       110       120       130
    //               01234567890123456789012345678901
    let ref_bases = "AAAAGGGGCCCCTTTTAAAACCCCGGGGTTTT";

    let variants = vec![
        create_simple_variant(102, "A", "T"),
        create_simple_variant(126, "GGT", "AAA"),
        create_simple_variant(115, "TAAAAC", "T"),
        create_simple_variant(109, "C", "CGGGG"),
    ];
    let variant_refs: Vec<&SimpleVariant> = variants.iter().collect();
    let alt_bases = build_alt_bases(ref_bases.as_bytes(), 100, 131, &variant_refs);

    assert_eq!(std::str::from_utf8(&alt_bases).unwrap(), "AATAGGGGCCGGGGCCTTTTCCCGGAAATTT");
}

// IndelDedupTest.testFilteredIndelRecovery
#[test]
fn test_filtered_indel_recovery() {
    let genome = mock_ref_genome("indel_dedup_filtered_indel_recovery", &[(CHR_1, TEST_REF_BASES)]);

    // an indel is not recovered when the other variant is largely unphased with it
    let combined_read_bases = format!(
        "{}{}A{}", // SNV 1
        &TEST_REF_BASES[0..21],
        &TEST_REF_BASES[30..32],
        &TEST_REF_BASES[33..70]
    );

    let indel_bases = format!("{}{}", &TEST_REF_BASES[0..21], &TEST_REF_BASES[30..70]);

    let del = create_indel_test_variant(20, 20, &indel_bases, &TEST_REF_BASES[20..30], &TEST_REF_BASES[20..21], 1);

    // variant not required by the INDEL
    let var1 = create_indel_test_variant(32, 32, &combined_read_bases, &TEST_REF_BASES[32..33], "A", 1);

    let mut variants = vec![del, var1];
    dedup_indels(&mut variants, &genome, CHR_1);

    assert!(variants[0].is_passing()); // del
    assert!(!variants[1].is_passing()); // var1
}

// IndelDedupTest.testIndelsLocalPhaseSetConditions
#[test]
fn test_indels_local_phase_set_conditions() {
    let genome = mock_ref_genome("indel_dedup_indels_local_phase_set_conditions", &[(CHR_1, TEST_REF_BASES)]);

    // an indel is not recovered when the other variant is largely unphased with it
    let combined_read_bases = format!(
        "{}{}A{}", // SNV 1
        &TEST_REF_BASES[0..21],
        &TEST_REF_BASES[30..32],
        &TEST_REF_BASES[33..70]
    );

    let del =
        create_indel_test_variant(20, 20, &combined_read_bases, &TEST_REF_BASES[20..30], &TEST_REF_BASES[20..21], 1);

    // variant part of the solution
    let var1 = create_indel_test_variant(32, 32, &combined_read_bases, &TEST_REF_BASES[32..33], "A", 1);

    // initially filtered but recoverable - no longer allow non-passing INDELs to be recovered
    // del.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());

    let mut variants = vec![del, var1];
    dedup_indels(&mut variants, &genome, CHR_1);

    assert!(variants[0].is_passing()); // del
    assert!(variants[1].is_passing()); // var1

    // test again but this time with the SNV mostly unphased with the INDEL
    variants[0].filters.clear();
    variants[0].filters.insert(SoftFilter::MinTumorQual);

    variants[1].tumor_counters[0].add_local_phase_set(2, 100, 0.0);

    dedup_indels(&mut variants, &genome, CHR_1);

    assert!(!variants[0].is_passing()); // del
    assert!(variants[1].is_passing()); // var1
}

// IndelDedupTest.testRepeatedInserts
#[test]
fn test_repeated_inserts() {
    // 4 inserts all with the same read context

    //                               10        20
    //                     0123456789012345678901234567890
    let chr2_ref_bases = "XGATCGATCGATCTCTCTCTCGATCGATCGA";
    let _genome = mock_ref_genome("indel_dedup_repeated_inserts", &[(CHR_1, TEST_REF_BASES), ("2", chr2_ref_bases)]);

    // the rest of the Java test, which built four phased CT inserts at positions 12, 14, 16 and 18 and
    // expected only the left-most to pass, is commented out in Java, so there is nothing more to port
}
