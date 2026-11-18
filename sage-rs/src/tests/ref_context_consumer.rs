//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/candidate/RefContextConsumerTest.java

use crate::candidate::mnv_length;
use crate::constants::MIN_SOFT_CLIP_MIN_BASE_QUAL;
use crate::evidence::exceeds_soft_clip_low_base_qual;

// RefContextConsumerTest.testMnvLength
#[test]
fn test_mnv_length() {
    assert_eq!(1, mnv_length(0, 0, b"A", b"CC"));
    assert_eq!(1, mnv_length(0, 0, b"AA", b"C"));

    assert_eq!(2, mnv_length(0, 0, b"AAA", b"CC"));
    assert_eq!(2, mnv_length(0, 0, b"AA", b"CCC"));
    assert_eq!(2, mnv_length(0, 0, b"AAC", b"CCC"));

    assert_eq!(3, mnv_length(0, 0, b"AAA", b"CCC"));
    assert_eq!(3, mnv_length(0, 0, b"ACA", b"CCCC"));
    assert_eq!(3, mnv_length(0, 0, b"ACAC", b"CCCC"));

    assert_eq!(3, mnv_length(0, 0, b"AAAA", b"CCCC"));
    assert_eq!(3, mnv_length(1, 0, b"AAAA", b"CCCC"));
    assert_eq!(2, mnv_length(2, 0, b"AAAA", b"CCCC"));
    assert_eq!(1, mnv_length(3, 0, b"AAAA", b"CCCC"));
}

// RefContextConsumerTest.testLowQualSoftClips
#[test]
fn test_low_qual_soft_clips() {
    let mut base_qualities = vec![MIN_SOFT_CLIP_MIN_BASE_QUAL as u8; 30];

    assert!(!exceeds_soft_clip_low_base_qual(&base_qualities, 0, 10));
    assert!(!exceeds_soft_clip_low_base_qual(&base_qualities, 20, 10));

    for q in base_qualities.iter_mut().take(25).skip(5) {
        *q = (MIN_SOFT_CLIP_MIN_BASE_QUAL - 1) as u8;
    }

    assert!(exceeds_soft_clip_low_base_qual(&base_qualities, 0, 10));
    assert!(exceeds_soft_clip_low_base_qual(&base_qualities, 20, 10));

    assert!(!exceeds_soft_clip_low_base_qual(&base_qualities, 0, 7));
    assert!(!exceeds_soft_clip_low_base_qual(&base_qualities, 23, 7));
}
