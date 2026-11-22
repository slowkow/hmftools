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
