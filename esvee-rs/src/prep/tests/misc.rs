//! MiscTest

use super::util::*;
use crate::prep::depth::DepthTracker;
use crate::prep::tracker::ReadIdTrimmer;

#[test]
fn test_read_trimmer() {
    let mut trimmer = ReadIdTrimmer::new(false);
    let read_id = "READ_01";
    assert_eq!(read_id.as_bytes(), trimmer.trim(read_id.as_bytes()));

    let mut trimmer = ReadIdTrimmer::new(true);
    let read_id = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
    assert_eq!(b"1:1673:32280:4946", trimmer.trim(read_id.as_bytes()));
    let read_id = "A00260:251:HLYGFDSXY:1:2:3:4";
    assert_eq!(b"1:2:3:4", trimmer.trim(read_id.as_bytes()));

    // disable on a read with different delims or too short
    let read_id = "A0:25:HL:1:2:3:4";
    assert_eq!(read_id.as_bytes(), trimmer.trim(read_id.as_bytes()));
    assert!(!trimmer.enabled());

    let mut trimmer = ReadIdTrimmer::new(true);
    let read_id = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
    assert_eq!(b"1:1673:32280:4946", trimmer.trim(read_id.as_bytes()));
    let read_id = "A00260:251:HLYGFDSX:1:2:3:4"; // diff delim position
    assert_eq!(read_id.as_bytes(), trimmer.trim(read_id.as_bytes()));
    assert!(!trimmer.enabled());

    // other formats are not currently supported
    let mut trimmer = ReadIdTrimmer::new(true);
    let read_id = "011852_2-UGAv3-2-1333458495";
    assert_eq!(read_id.as_bytes(), trimmer.trim(read_id.as_bytes()));
    assert!(!trimmer.enabled());
}

fn assert_close(expected: f64, actual: f64, delta: f64) {
    assert!((expected - actual).abs() <= delta, "expected {expected} actual {actual}");
}

#[test]
fn test_depth_tracker() {
    let ref_bases = generate_random_bases(500);
    let mut depth_tracker = DepthTracker::new(1, 1000, 100);
    for start in [1, 10, 120, 180, 1001, 951] {
        let read = sam_record(&next_id(), CHR_1, start, &ref_bases[..50], "50M");
        depth_tracker.process_read(&read);
    }
    assert_close(1.0, depth_tracker.calc_depth(1), 0.1);
    assert_close(1.0, depth_tracker.calc_depth(100), 0.1);
    assert_close(0.7, depth_tracker.calc_depth(101), 0.1);
    assert_close(0.7, depth_tracker.calc_depth(199), 0.1);
    assert_close(0.0, depth_tracker.calc_depth(301), 0.1);
    assert_close(0.5, depth_tracker.calc_depth(1000), 0.1);
    assert_close(0.0, depth_tracker.calc_depth(1001), 0.1);
}

// not in Java: -low_base_qual_threshold (common::base_qual)

#[test]
fn below_qual_threshold_keeps_java_byte_semantics() {
    use crate::prep::filters::below_qual_threshold;
    // a quality above 127 is a negative Java byte, and so low, at any threshold
    for t in [0u8, 20, 26, 40, 93] {
        for q in 0..=255u8 {
            assert_eq!(q <= t || q >= 128, below_qual_threshold(q, t), "quality {q}, threshold {t}");
        }
    }
}

#[test]
fn low_base_qual_threshold_option_keeps_low_qual_reads() {
    use crate::common::base_qual::set_test_base_qual_thresholds;
    use crate::prep::filters::filter_low_qual_read;
    use crate::prep::filters::ReadFilterConfig;
    let bases = generate_random_bases(100);
    let config = ReadFilterConfig::default();
    let filtered_at = |bands: Option<(u8, u8)>| {
        set_test_base_qual_thresholds(bands);
        let mut rec = sam_record_flags(&next_id(), CHR_1, 1000, &bases, "100M", build_flags(true, false, false), 24);
        let filtered = filter_low_qual_read(&mut rec, &config);
        set_test_base_qual_thresholds(None);
        filtered
    };
    // every base at Q24 is low quality in Java's band (Q26 and below), so the read is dropped
    assert!(filtered_at(None));
    // with the band at Q20 and below it is kept
    assert!(!filtered_at(Some((20, 30))));
}

// not in Java: finding P12, one read name prefix for the run

const RUN_A: &str = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
// another run, prefix of the same length: Java cuts it to the same lane:tile:x:y id as RUN_A
const RUN_B: &str = "A00260:252:HLYGFDSXZ:1:1673:32280:4946";
// another run, prefix one character shorter: Java switches trimming off here
const RUN_C: &str = "A0026:251:HLYGFDSXY:1:1673:32280:4946";

#[test]
fn corrected_trimmer_cuts_only_the_run_prefix() {
    use crate::common::compat::with_java;
    use crate::prep::tracker::trim_prefix;
    with_java(false, || {
        let prefix = trim_prefix(RUN_A.as_bytes());
        assert_eq!(Some(b"A00260:251:HLYGFDSXY:".as_slice()), prefix.as_deref());
        let mut trimmer = ReadIdTrimmer::for_run(true, prefix.as_ref());
        // a name of another run keeps its whole name, so the two fragments stay apart
        assert_eq!(RUN_B.as_bytes(), trimmer.trim(RUN_B.as_bytes()));
        assert_eq!(b"1:1673:32280:4946", trimmer.trim(RUN_A.as_bytes()));
        assert_eq!(RUN_C.as_bytes(), trimmer.trim(RUN_C.as_bytes()));
        // and trimming stays on
        assert!(trimmer.enabled());
        assert_eq!(b"1:1673:32280:4946", trimmer.trim(RUN_A.as_bytes()));
        // the prefix alone is not cut to an empty id
        assert_eq!(b"A00260:251:HLYGFDSXY:", trimmer.trim(b"A00260:251:HLYGFDSXY:"));
        // no prefix (a first name without three ':'): no trimming
        assert_eq!(None, trim_prefix(b"011852_2-UGAv3-2-1333458495"));
        assert_eq!(None, trim_prefix(b"A00260:251:HLYGFDSXY:"));
        let mut trimmer = ReadIdTrimmer::for_run(true, None);
        assert!(!trimmer.enabled());
        assert_eq!(RUN_A.as_bytes(), trimmer.trim(RUN_A.as_bytes()));
    });
}

#[test]
fn java_trimmer_ids_depend_on_the_first_name() {
    use crate::common::compat::with_java;
    use crate::prep::tracker::trim_prefix;
    with_java(true, || {
        // two runs with prefixes of one length: both names are cut to the same id
        let mut trimmer = ReadIdTrimmer::for_run(true, trim_prefix(RUN_A.as_bytes()).as_ref());
        assert_eq!(trimmer.trim(RUN_A.as_bytes()), trimmer.trim(RUN_B.as_bytes()));
        // a shorter prefix switches trimming off, and the mate of a read cut before keeps its whole name
        let mut trimmer = ReadIdTrimmer::for_run(true, None);
        assert_eq!(b"1:1673:32280:4946", trimmer.trim(RUN_A.as_bytes()));
        assert_eq!(RUN_C.as_bytes(), trimmer.trim(RUN_C.as_bytes()));
        assert_eq!(RUN_A.as_bytes(), trimmer.trim(RUN_A.as_bytes()));
        assert!(!trimmer.enabled());
    });
}
