//! Ports of ReadCacheTest (top level) and duplicate/ReadCacheTest.

use super::util::*;
use crate::dup::*;
use rust_htslib::bam;

/// new ReadCache(groupSize, maxSoftClipLength, useFragmentOrientation)
fn read_cache(group_size: i32, max_soft_clip: i32) -> ReadCache {
    ReadCache::new(
        group_size, max_soft_clip, false, DEFAULT_POP_DISTANCE_CHECK, DEFAULT_DYNAMIC_READ_COUNT_THRESHOLD, DEFAULT_LOG_READ_COUNT_THRESHOLD,
    )
}

fn create_record(ids: &mut ReadIdGenerator, chromosome: &str, read_start: i32) -> bam::Record {
    create_sam_record(
        &ids.next_id(), chromosome, read_start, &test_read_bases(), TEST_READ_CIGAR, CHR_1, 200, false, false, None, true,
        Some(TEST_READ_CIGAR),
    )
}

fn create_spliced_record(ids: &mut ReadIdGenerator, read_start: i32, cigar: &str) -> bam::Record {
    create_sam_record(&ids.next_id(), CHR_1, read_start, &"A".repeat(100), cigar, CHR_1, 3000, true, false, None, false, Some(TEST_READ_CIGAR))
}

fn process_unspliced_duplicates(read_cache: &mut ReadCache, ids: &mut ReadIdGenerator, names: &[String]) {
    read_cache.process_read(create_record(ids, CHR_1, 700), names);
    read_cache.process_read(create_record(ids, CHR_1, 700), names);
    read_cache.process_read(create_record(ids, CHR_1, 900), names);
    read_cache.process_read(create_record(ids, CHR_1, 900), names);
}

/// ReadCacheTest.testBasics
#[test]
fn test_basics() {
    let names = names();
    let ids = &mut ReadIdGenerator::new();
    let mut read_cache = read_cache(100, 100);

    let read1 = create_record(ids, CHR_1, 50);

    read_cache.process_read(read1, &names);

    assert_eq!(1, read_cache.cached_read_position_groups());
    assert_eq!(1, read_cache.cached_frag_coord_groups());
    assert_eq!(1, read_cache.cached_read_count());

    let mut frag_coord_reads = read_cache.pop_reads();
    assert!(frag_coord_reads.is_none());

    let read2a = create_record(ids, CHR_1, 80);
    let read2b = create_record(ids, CHR_1, 80);

    read_cache.process_read(read2a, &names);
    read_cache.process_read(read2b, &names);

    let read3 = create_record(ids, CHR_1, 120);

    read_cache.process_read(read3, &names);

    assert_eq!(2, read_cache.cached_read_position_groups());
    assert_eq!(4, read_cache.cached_read_count());
    assert_eq!(3, read_cache.cached_frag_coord_groups());

    frag_coord_reads = read_cache.pop_reads();
    assert!(frag_coord_reads.is_none());

    let read4 = create_record(ids, CHR_1, 199);

    read_cache.process_read(read4, &names);

    assert_eq!(2, read_cache.cached_read_position_groups());

    // the next read triggers the first group to be processed
    let read5 = create_record(ids, CHR_1, 201);
    read_cache.process_read(read5, &names);

    frag_coord_reads = read_cache.pop_reads();
    let popped = frag_coord_reads.expect("popped reads");
    assert_eq!(1, popped.single_reads.len());
    assert_eq!(1, popped.duplicate_groups.len());

    assert_eq!(2, read_cache.cached_read_position_groups());
    assert_eq!(3, read_cache.cached_read_count());
    assert_eq!(3, read_cache.cached_frag_coord_groups());

    // next read is on a new chromosome so clears all previous
    let read6 = create_record(ids, CHR_2, 50);
    read_cache.process_read(read6, &names);

    frag_coord_reads = read_cache.pop_reads();
    let popped = frag_coord_reads.expect("popped reads");
    assert_eq!(3, popped.single_reads.len());
    assert_eq!(0, popped.duplicate_groups.len());

    assert_eq!(1, read_cache.cached_read_position_groups());
    assert_eq!(1, read_cache.cached_read_count());
    assert_eq!(1, read_cache.cached_frag_coord_groups());

    // some more duplicates in the next group
    let read7a = create_record(ids, CHR_1, 150);
    let read7b = create_record(ids, CHR_1, 150);
    let read7c = create_record(ids, CHR_1, 150);
    read_cache.process_read(read7a, &names);
    read_cache.process_read(read7b, &names);
    read_cache.process_read(read7c, &names);

    assert_eq!(2, read_cache.cached_read_position_groups());
    assert_eq!(4, read_cache.cached_read_count());
    assert_eq!(2, read_cache.cached_frag_coord_groups());

    frag_coord_reads = read_cache.evict_all();
    let popped = frag_coord_reads.expect("evicted reads");
    assert_eq!(1, popped.single_reads.len());
    assert_eq!(1, popped.duplicate_groups.len());
    assert_eq!(3, popped.duplicate_groups[0].total_read_count());

    assert_eq!(0, read_cache.cached_read_position_groups());
    assert_eq!(0, read_cache.cached_read_count());
    assert_eq!(0, read_cache.cached_frag_coord_groups());
}

/// ReadCacheTest.testMixedPrimarySupplementaries
#[test]
fn test_mixed_primary_supplementaries() {
    let names = names();
    let ids = &mut ReadIdGenerator::new();
    let mut read_cache = read_cache(100, 100);
    let bases = test_read_bases();

    let supp_data1 = SuppData::new(CHR_2, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

    let read1 = create_sam_record(
        &ids.next_id(), CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 200, false, false, Some(&supp_data1), true, Some(TEST_READ_CIGAR),
    );

    let supp_data2 = SuppData::new(CHR_2, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

    let read2 = create_sam_record(
        &ids.next_id(), CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 200, false, true, Some(&supp_data2), true, Some(TEST_READ_CIGAR),
    );

    read_cache.process_read(read1, &names);
    read_cache.process_read(read2, &names);

    assert_eq!(1, read_cache.cached_read_position_groups());
    assert_eq!(2, read_cache.cached_read_count());
    assert_eq!(2, read_cache.cached_frag_coord_groups());
}

/// ReadCacheTest.testFragmentCoordsAreOnlyPoppedOnce
#[test]
fn test_fragment_coords_are_only_popped_once() {
    let names = names();
    let read_cigar = "12S80M51S";
    let mate_cigar = "143M";
    let bases = "A".repeat(143);

    let read1 = create_sam_record("READ_001", CHR_1, 1, &bases, read_cigar, CHR_1, 1, false, false, None, true, Some(mate_cigar));
    let mate1 = create_sam_record("READ_001", CHR_1, 1, &bases, mate_cigar, CHR_1, 1, true, false, None, false, Some(read_cigar));
    let read2 = create_sam_record("READ_002", CHR_1, 1, &bases, read_cigar, CHR_1, 1, false, false, None, true, Some(mate_cigar));
    let mate2 = create_sam_record("READ_002", CHR_1, 1, &bases, mate_cigar, CHR_1, 1, true, false, None, false, Some(read_cigar));

    let coords1 = FragmentCoords::from_read(&read1, &names, false);
    let mate_coords1 = FragmentCoords::from_read(&mate1, &names, false);
    let coords2 = FragmentCoords::from_read(&read2, &names, false);
    let mate_coords2 = FragmentCoords::from_read(&mate2, &names, false);

    // FragmentCoords.equals compares the keys
    assert_eq!(coords1.key, coords2.key);
    assert_eq!(mate_coords1.key, mate_coords2.key);
    assert_ne!(coords1.key, mate_coords1.key);

    // the SAM strings of the Java test: the records are told apart by name and flags
    let id = |r: &bam::Record| (read_name(r), r.flags());
    let mut expected_duplicate_groups = vec![vec![id(&read1), id(&read2)], vec![id(&mate1), id(&mate2)]];

    let mut read_cache = read_cache(DEFAULT_GROUP_SIZE, DEFAULT_MAX_SOFT_CLIP);
    read_cache.process_read(read1, &names);
    read_cache.process_read(mate1, &names);
    read_cache.process_read(read2, &names);
    read_cache.process_read(mate2, &names);

    let fragment_coord_reads = read_cache.evict_all().expect("evicted reads");

    assert_eq!(2, fragment_coord_reads.duplicate_groups.len());
    assert_eq!(0, fragment_coord_reads.single_reads.len());

    let mut actual_duplicate_groups: Vec<Vec<(String, u16)>> = fragment_coord_reads
        .duplicate_groups
        .iter()
        .map(|g| {
            let mut ids: Vec<(String, u16)> = g.reads.iter().map(id).collect();
            ids.sort();
            ids
        })
        .collect();
    actual_duplicate_groups.sort();
    expected_duplicate_groups.iter_mut().for_each(|g| g.sort());
    expected_duplicate_groups.sort();

    assert_eq!(expected_duplicate_groups, actual_duplicate_groups);
}

/// ReadCacheTest.testReverseSplicedReadWithinThresholdIsStillMarkedDuplicate
#[test]
fn test_reverse_spliced_read_within_threshold_is_still_marked_duplicate() {
    let names = names();
    let ids = &mut ReadIdGenerator::new();
    let mut read_cache = read_cache(DEFAULT_GROUP_SIZE, DEFAULT_MAX_SOFT_CLIP);

    let spliced_read = create_spliced_record(ids, 100, "50M2000N50M");
    let spliced_duplicate = create_spliced_record(ids, 600, "50M1500N50M");
    assert_eq!(FragmentCoords::from_read(&spliced_read, &names, false).key, FragmentCoords::from_read(&spliced_duplicate, &names, false).key);

    read_cache.process_read(spliced_read, &names);
    read_cache.process_read(spliced_duplicate, &names);

    let mut frag_coord_reads = read_cache.pop_reads();
    assert!(frag_coord_reads.is_none());

    process_unspliced_duplicates(&mut read_cache, ids, &names);

    frag_coord_reads = read_cache.pop_reads();
    let popped = frag_coord_reads.expect("popped reads");
    assert_eq!(1, popped.duplicate_groups.len());
    assert_eq!(0, popped.single_reads.len());

    // this read takes the read position past the threshold from the second spliced read, releasing both
    let read3 = create_record(ids, CHR_1, 1700);
    read_cache.process_read(read3, &names);

    frag_coord_reads = read_cache.pop_reads();
    let popped = frag_coord_reads.expect("popped reads");
    assert_eq!(2, popped.duplicate_groups.len());
    assert_eq!(0, popped.single_reads.len());
}

/// ReadCacheTest.testReverseSplicedReadBeyondThresholdMissesDuplicate
#[test]
fn test_reverse_spliced_read_beyond_threshold_misses_duplicate() {
    let names = names();
    let ids = &mut ReadIdGenerator::new();
    let mut read_cache = read_cache(DEFAULT_GROUP_SIZE, DEFAULT_MAX_SOFT_CLIP);

    let spliced_read = create_spliced_record(ids, 100, "50M2000N50M");
    let spliced_duplicate = create_spliced_record(ids, 1200, "50M900N50M");
    assert_eq!(FragmentCoords::from_read(&spliced_read, &names, false).key, FragmentCoords::from_read(&spliced_duplicate, &names, false).key);

    read_cache.process_read(spliced_read, &names);

    let mut frag_coord_reads = read_cache.pop_reads();
    assert!(frag_coord_reads.is_none());

    process_unspliced_duplicates(&mut read_cache, ids, &names);

    frag_coord_reads = read_cache.pop_reads();
    let popped = frag_coord_reads.expect("popped reads");
    assert_eq!(1, popped.duplicate_groups.len());
    assert_eq!(0, popped.single_reads.len());

    // this read takes the read position past the threshold, so the spliced read is released on its own
    let read3 = create_record(ids, CHR_1, 1100);
    read_cache.process_read(read3, &names);

    frag_coord_reads = read_cache.pop_reads();
    let popped = frag_coord_reads.expect("popped reads");
    assert_eq!(1, popped.duplicate_groups.len());
    assert_eq!(1, popped.single_reads.len());

    read_cache.process_read(spliced_duplicate, &names);

    let popped = read_cache.evict_all().expect("evicted reads");
    assert_eq!(0, popped.duplicate_groups.len());
    assert_eq!(2, popped.single_reads.len());
}

/// ReadCacheTest.testReadCacheReadPositionWithinGroupBounds
#[test]
fn test_read_cache_read_position_within_group_bounds() {
    let names = names();
    let read_length = 143;
    let alignment_start = 84026397;
    let cigar = "133M10S";
    let read_bases = "A".repeat(read_length);
    let mate_cigar = format!("{read_length}M");

    let read1 = create_sam_record("READ_001", CHR_1, alignment_start, &read_bases, cigar, CHR_2, 1, false, false, None, true, Some(&mate_cigar));
    let read2 = clone_sam_record(&read1, "READ_002");

    let mut read_cache = read_cache(DEFAULT_GROUP_SIZE, DEFAULT_MAX_SOFT_CLIP);
    read_cache.process_read(read1, &names);
    read_cache.process_read(read2, &names);
    let fragment_coord_reads = read_cache.evict_all().expect("evicted reads");

    assert_eq!(0, fragment_coord_reads.single_reads.len());
    assert_eq!(1, fragment_coord_reads.duplicate_groups.len());

    let duplicate_group_reads = &fragment_coord_reads.duplicate_groups[0].reads;
    assert_eq!(2, duplicate_group_reads.len());
    let mut group_names: Vec<String> = duplicate_group_reads.iter().map(read_name).collect();
    group_names.sort();
    assert_eq!(vec!["READ_001", "READ_002"], group_names);
}

/// duplicate/ReadCacheTest.testReadCacheAboveDynamicThreshold
#[test]
fn test_read_cache_above_dynamic_threshold() {
    let names = names();
    let ids = &mut ReadIdGenerator::new();
    let test_read_bases = "A".repeat(10);
    let test_cigar = "10M";

    let mut read_cache = ReadCache::new(100, 50, false, 100, 100, 100);

    // create reads above the dynamic threshold but without standard purging occurring
    let read_pos1 = 100;
    let read_pos2 = 110;
    let read_pos3 = 120;

    for _ in 0..40 {
        for read_pos in [read_pos1, read_pos2, read_pos3] {
            let read =
                create_sam_record(&ids.next_id(), CHR_1, read_pos, &test_read_bases, test_cigar, CHR_1, 200, false, false, None, true, Some(test_cigar));
            read_cache.process_read(read, &names);
        }
    }

    assert_eq!(120, read_cache.cached_read_count());
    assert_eq!(3, read_cache.cached_frag_coord_groups());

    let fragment_coord_reads = read_cache.pop_reads().expect("popped reads");
    assert_eq!(80, fragment_coord_reads.total_read_count());
}
