//! Tests of the port: the fit, rounding, the file formats, and each corrected finding in both modes.

use crate::compat::with_java;
use crate::fitter::{allocation_lines, round, SignatureAllocation};
use crate::lsq::LeastSquaresFit;
use crate::matrix::load_matrix_lines;
use crate::position_freq::{initialise_position_cache, standard_chromosome_lengths, PositionFrequencies};
use crate::snv_loader::{bucket_names, variant_context};

#[test]
fn bucket_order() {
    let names = bucket_names();
    assert_eq!(names.len(), 96);
    assert_eq!(names[0], "C>A_ACA");
    assert_eq!(names[15], "C>A_TCT");
    assert_eq!(names[16], "C>G_ACA");
    assert_eq!(names[95], "T>G_TTT");
}

#[test]
fn contexts_on_the_pyrimidine_strand() {
    assert_eq!(variant_context("C", "T", "ACG").unwrap(), "C>T_ACG");
    // G>A in AGC is C>T in GCT on the other strand
    assert_eq!(variant_context("G", "A", "AGC").unwrap(), "C>T_GCT");
    assert_eq!(variant_context("A", "C", "TAG").unwrap(), "T>G_CTA");
    assert!(variant_context("A", "C", "").is_err());
}

#[test]
fn java_rounding_units() {
    // Math.pow(10, -5) is 0x3ee4f8b588e368f0 in Java, one unit in the last place below 1e-5
    assert_eq!(hmf_common::java_strictmath::pow(10.0, -5.0).to_bits(), 0x3ee4f8b588e368f0);
    assert_eq!(hmf_common::java_strictmath::pow(10.0, -3.0).to_bits(), 0x3f50624dd2f1a9fc);
    assert_eq!(round(632.5814, 3), 632.581);
    assert_eq!(hmf_common::fmt::java_fixed(round(0.068564, 5), 5), "0.06856");
}

#[test]
fn nnls_exact_fit() {
    // two signatures over four buckets; counts = 3 * s1 + 2 * s2
    let sigs = vec![vec![0.5, 0.0], vec![0.5, 0.25], vec![0.0, 0.25], vec![0.0, 0.5]];
    let counts: Vec<f64> = (0..4).map(|r| 3.0 * sigs[r][0] + 2.0 * sigs[r][1]).collect();
    let mut f = LeastSquaresFit::new(4, 2).unwrap();
    f.initialise(&sigs, &counts);
    f.solve().unwrap();
    assert!((f.contribs()[0] - 3.0).abs() < 1e-12 && (f.contribs()[1] - 2.0).abs() < 1e-12);
}

#[test]
fn nnls_keeps_contributions_non_negative() {
    // the unconstrained fit of b would give s2 a negative contribution
    let sigs = vec![vec![1.0, 1.0], vec![0.0, 1.0], vec![0.0, 0.0]];
    let counts = vec![2.0, -1.0, 0.0];
    let mut f = LeastSquaresFit::new(3, 2).unwrap();
    f.initialise(&sigs, &counts);
    f.solve().unwrap();
    assert!(f.contribs().iter().all(|c| *c >= 0.0));
    assert!((f.contribs()[0] - 2.0).abs() < 1e-12 && f.contribs()[1] == 0.0);
}

#[test]
fn allocation_file_format() {
    let a = vec![
        SignatureAllocation { signature: "Sig1".into(), allocation: round(12.3456, 3), percent: round(0.123456, 5) },
        SignatureAllocation { signature: "MISALLOC".into(), allocation: 0.0, percent: 0.0 },
    ];
    assert_eq!(allocation_lines(&a), "signature\tallocation\tpercent\nSig1\t12.346\t0.12346\nMISALLOC\t0.000\t0.00000\n");
    assert_eq!(allocation_lines(&[]), "signature\tallocation\tpercent\n");
}

#[test]
fn matrix_rows_of_another_length_stay_zero() {
    let lines: Vec<String> = ["A,B", "1,2", "3", "5,"].iter().map(|s| s.to_string()).collect();
    let mut names = Vec::new();
    let m = load_matrix_lines(&lines, &mut names, None, false).unwrap();
    assert_eq!(names, vec!["A", "B"]);
    assert_eq!(m.data, vec![vec![1.0, 2.0], vec![0.0, 0.0], vec![5.0, 0.0]]);
}

#[test]
fn si4_sample_list_shifts_the_ignored_column() {
    let lines: Vec<String> = ["BucketName,S1,S2", "C>A_ACA,1,2"].iter().map(|s| s.to_string()).collect();
    // Java: the configured sample is in the list the header is added to, so the list holds S1 twice and
    // every row has one field too few: all counts stay 0, for three "samples"
    let mut ids = vec!["S1".to_string()];
    let m = load_matrix_lines(&lines, &mut ids, Some(&["BucketName"]), false).unwrap();
    assert_eq!(ids, vec!["S1", "S1", "S2"]);
    assert_eq!(m.data, vec![vec![0.0, 0.0, 0.0]]);
    let mut fresh = Vec::new();
    let m = load_matrix_lines(&lines, &mut fresh, Some(&["BucketName"]), false).unwrap();
    assert_eq!(fresh, vec!["S1", "S2"]);
    assert_eq!(m.data, vec![vec![1.0, 2.0]]);
}

#[test]
fn si1_position_frequencies_kept_only_when_corrected() {
    with_java(true, || {
        let mut pf = PositionFrequencies::new(1_000_000);
        pf.add_position("chr1", 1_500_000);
        assert!(pf.chr_pos_bucket_frequencies.is_empty());
    });
    with_java(false, || {
        let mut pf = PositionFrequencies::new(1_000_000);
        pf.add_position("chr1", 1_500_000);
        pf.add_position("chr1", 1_999_999);
        assert_eq!(pf.chr_pos_bucket_frequencies["chr1"][&1_000_000], 2);
    });
}

#[test]
fn si2_position_cache_starts() {
    let lengths = standard_chromosome_lengths(true);
    // a bucket larger than chr1: one bucket per chromosome
    let (java_index, java_size) = with_java(true, || initialise_position_cache(true, 300_000_000, &lengths));
    assert_eq!(java_size, 24);
    assert_eq!(java_index["chr2"], 0);
    assert_eq!(java_index["chrY"], 0);
    let (index, size) = with_java(false, || initialise_position_cache(true, 300_000_000, &lengths));
    assert_eq!(size, 24);
    assert_eq!(index["chr2"], 1);
    assert_eq!(index["chrY"], 23);
    // the usual bucket sizes agree
    let a = with_java(true, || initialise_position_cache(true, 1_000_000, &lengths));
    let b = with_java(false, || initialise_position_cache(true, 1_000_000, &lengths));
    assert_eq!(a, b);
    assert_eq!(a.0["chr2"], 249);
}
