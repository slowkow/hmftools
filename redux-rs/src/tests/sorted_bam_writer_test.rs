//! Port of SortedBamWriterTest.

use super::util::*;
use crate::out::RawBamOut;
use crate::sorted::SortedBamWriter;
use rust_htslib::bam;

fn create_read(chromosome: &str, position: i32) -> bam::Record {
    create_sam_record(TEST_READ_ID, chromosome, position, "", TEST_READ_CIGAR, CHR_1, 100, false, false, None, false, Some(TEST_READ_CIGAR))
}

/// SortedBamWriterTest.testSortedReadCache
///
/// The Java test passes no writer; the port always writes, so the records go to a scratch file. The last Java
/// assertion, on writeCount(), has no counterpart: the port does not count its write batches (a statistic only).
#[test]
fn test_sorted_read_cache() {
    let path = std::env::temp_dir().join(format!("redux_rs_tests_sorted_{}.bam", std::process::id()));
    let mut out = RawBamOut::create(&path.to_string_lossy(), 1).unwrap();
    let (tid1, tid2) = (chromosome_ordinal(CHR_1), chromosome_ordinal(CHR_2));

    // new SortedBamConfig(20, 5, 1)
    let mut read_cache = SortedBamWriter::with_config(20, 5, 1);

    read_cache.initialise_start_position(tid1, 50);

    read_cache.set_upper_writable_position(100);
    read_cache.set_upper_bound_position(110);

    // first test boundaries
    assert!(read_cache.can_write_record(&create_read(CHR_1, 105)));
    assert!(!read_cache.can_write_record(&create_read(CHR_2, 105))); // wrong chromosome
    assert!(!read_cache.can_write_record(&create_read(CHR_1, 115))); // too high

    read_cache.add_record(create_read(CHR_1, 100), &mut out).unwrap();
    read_cache.add_record(create_read(CHR_1, 110), &mut out).unwrap();

    assert_eq!(0, read_cache.reads_written);
    assert_eq!(2, read_cache.cached());

    // add a record past the buffer position
    read_cache.set_upper_bound_position(170);
    read_cache.set_upper_writable_position(120);
    read_cache.add_record(create_read(CHR_1, 170), &mut out).unwrap();
    assert_eq!(2, read_cache.reads_written);
    assert_eq!(1, read_cache.cached());

    assert!(!read_cache.can_write_record(&create_read(CHR_1, 40))); // too low

    read_cache.add_record(create_read(CHR_1, 180), &mut out).unwrap();
    assert_eq!(2, read_cache.reads_written);
    assert_eq!(2, read_cache.cached());

    // grow capacity for records within the position buffer
    read_cache.set_upper_bound_position(200);

    for _ in 0..10 {
        read_cache.add_record(create_read(CHR_1, 181), &mut out).unwrap();
    }

    assert_eq!(2, read_cache.reads_written);
    assert_eq!(12, read_cache.cached());

    read_cache.set_upper_writable_position(190);
    read_cache.add_record(create_read(CHR_1, 200), &mut out).unwrap();

    assert_eq!(14, read_cache.reads_written);
    assert_eq!(1, read_cache.cached());

    read_cache.flush(&mut out).unwrap();

    assert_eq!(15, read_cache.reads_written);
    assert_eq!(0, read_cache.cached());

    // new chromosome
    read_cache.initialise_start_position(tid2, 1);
    read_cache.set_upper_bound_position(50);

    read_cache.add_record(create_read(CHR_2, 50), &mut out).unwrap();

    assert_eq!(15, read_cache.reads_written);
    assert_eq!(1, read_cache.cached());

    assert_eq!(15, out.written);
    out.close().unwrap();
    let _ = std::fs::remove_file(&path);
}
