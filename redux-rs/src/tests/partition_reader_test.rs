//! Port of PartitionReaderTest.

use super::util::*;
use crate::rec::{self, FLAG_MATE_UNMAPPED, FLAG_UNMAPPED};
use rust_htslib::bam;
use std::collections::HashMap;

/// The reference genome of the Java test class: chromosome 1 is REF_BASES repeated 40 times
fn ref_genome() -> MockRefGenome {
    let bases = ref_bases_repeat_40();
    let mut ref_genome = MockRefGenome::new(false);
    ref_genome.chromosome_lengths.insert(CHR_1.to_string(), bases.len() as i32);
    ref_genome.ref_genome_map.insert(CHR_1.to_string(), bases);
    ref_genome
}

/// new ReduxConfig(mRefGenome, false, false, true, new ReadUnmapper(Collections.emptyMap()))
fn consensus_env() -> TestEnv {
    TestEnv::new(ref_genome(), true, read_unmapper(HashMap::new()))
}

/// new ReduxConfig(mRefGenome, false, false, false, READ_UNMAPPER)
fn unmapping_env() -> TestEnv {
    TestEnv::new(ref_genome(), false, test_read_unmapper())
}

fn read_unmapped(read: &bam::Record) -> bool {
    rec::has(read, FLAG_UNMAPPED)
}
fn mate_unmapped(read: &bam::Record) -> bool {
    rec::has(read, FLAG_MATE_UNMAPPED)
}

/// PartitionReaderTest.testMultipleSupplementaries
#[test]
fn test_multiple_supplementaries() {
    // one primary has multiple supplementaries and the other has a different supp
    let env = consensus_env();
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();

    let read_pos = 100;
    let mate_pos = 200;
    let supp_pos = 1800;

    let read1_name = ids.next_id();
    let mut read1 =
        create_sam_record(&read1_name, CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, None, true, Some(TEST_READ_CIGAR));

    let read1_supps = [
        SuppData::new(CHR_1, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
        SuppData::new(CHR_2, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
        SuppData::new(CHR_3, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
    ];

    set_string_attribute(&mut read1, b"SA", &alignments_to_sam_tag(&read1_supps));
    set_string_attribute(&mut read1, b"MC", TEST_READ_CIGAR);

    let read2_name = ids.next_id();
    let mut read2 = create_sam_record(
        &read2_name, CHR_1, read_pos + 2, &bases, "2S98M", CHR_1, mate_pos, false, false,
        Some(&SuppData::new(CHR_2, supp_pos + 1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );
    set_string_attribute(&mut read2, b"MC", TEST_READ_CIGAR);

    partition_reader.setup_region(CHR_1, 1, 1000);

    partition_reader.process_read(&mut read1);
    partition_reader.process_read(&mut read2);
    partition_reader.flush_read_positions();

    let mut mate1 =
        create_sam_record(&read1_name, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false, None, false, Some(TEST_READ_CIGAR));
    set_second_in_pair(&mut mate1);

    let mut mate2 = create_sam_record(
        &read2_name, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos + 2, true, false, None, false, Some(TEST_READ_CIGAR),
    );
    set_second_in_pair(&mut mate2);

    partition_reader.process_read(&mut mate1);
    partition_reader.process_read(&mut mate2);

    partition_reader.setup_region(CHR_1, 1001, 2000);

    let mut read1_supp1 = create_sam_record(
        &read1_name, CHR_1, supp_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, true,
        Some(&SuppData::new(CHR_1, read_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.process_read(&mut read1_supp1);
    partition_reader.post_process_region();

    let mut read1_supp2 = create_sam_record(
        &read1_name, CHR_2, supp_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, true,
        Some(&SuppData::new(CHR_1, read_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    let mut read2_supp1 = create_sam_record(
        &read2_name, CHR_2, supp_pos + 1, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, true,
        Some(&SuppData::new(CHR_1, read_pos + 2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.setup_region(CHR_2, 1001, 2000);
    partition_reader.process_read(&mut read1_supp2);
    partition_reader.process_read(&mut read2_supp1);
    partition_reader.post_process_region();

    let mut read1_supp3 = create_sam_record(
        &read1_name, CHR_3, supp_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, true,
        Some(&SuppData::new(CHR_1, read_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.setup_region(CHR_3, 1001, 2000);
    partition_reader.process_read(&mut read1_supp3);
    partition_reader.post_process_region();

    assert_eq!(8, partition_reader.non_consensus_write_count());
    assert_eq!(1, partition_reader.consensus_write_count()); // since the supplementaries aren't mapped to the same location
}

/// PartitionReaderTest.testInconsistentSupplementaries
#[test]
fn test_inconsistent_supplementaries() {
    // supplementaries arriving out of order and mapped to different locations
    let env = consensus_env();
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();
    let supp = |position: i32| SuppData::new(CHR_1, position, SUPP_POS_STRAND, TEST_READ_CIGAR, 1);

    let supp2_pos1 = 950;
    let read_pos = 1950;
    let supp2_pos2 = 2950;
    let mate_pos = 3500;
    let supp1_pos1 = 4500;
    let supp1_pos2 = 5500;

    let read_id1 = ids.next_id();
    let read_id2 = ids.next_id();

    let mut read2_supp1 = create_sam_record(
        &read_id2, CHR_1, supp2_pos1, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, true, Some(&supp(read_pos)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.setup_region(CHR_1, 1, 999);
    partition_reader.process_read(&mut read2_supp1);
    partition_reader.post_process_region();

    let mut read1 = create_sam_record(
        &read_id1, CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, Some(&supp(supp1_pos2)), true, Some(TEST_READ_CIGAR),
    );
    set_string_attribute(&mut read1, b"MC", TEST_READ_CIGAR);

    let mut read2 = create_sam_record(
        &read_id2, CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, Some(&supp(supp2_pos1)), true, Some(TEST_READ_CIGAR),
    );
    set_string_attribute(&mut read2, b"MC", TEST_READ_CIGAR);

    let mut region_start = 1000;
    partition_reader.setup_region(CHR_1, region_start, region_start + 999);
    partition_reader.process_read(&mut read1);
    partition_reader.process_read(&mut read2);
    partition_reader.post_process_region();

    let mut read2_supp2 = create_sam_record(
        &read_id2, CHR_1, supp2_pos2, &bases, TEST_READ_CIGAR, CHR_1, read_pos, false, true, Some(&supp(mate_pos)), true, Some(TEST_READ_CIGAR),
    );
    set_second_in_pair(&mut read2_supp2);

    region_start += 1000;
    partition_reader.setup_region(CHR_1, region_start, region_start + 999);
    partition_reader.process_read(&mut read2_supp2);
    partition_reader.post_process_region();

    let mut mate1 = create_sam_record(
        &read_id1, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false, Some(&supp(supp1_pos1)), false, Some(TEST_READ_CIGAR),
    );
    set_second_in_pair(&mut mate1);

    let mut mate2 = create_sam_record(
        &read_id2, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false, Some(&supp(supp2_pos2)), false, Some(TEST_READ_CIGAR),
    );
    set_second_in_pair(&mut mate2);

    region_start += 1000;
    partition_reader.setup_region(CHR_1, region_start, region_start + 999);
    partition_reader.process_read(&mut mate1);
    partition_reader.process_read(&mut mate2);
    partition_reader.post_process_region();

    let mut read1_supp1 = create_sam_record(
        &read_id1, CHR_1, supp1_pos2, &bases, TEST_READ_CIGAR, CHR_1, read_pos, false, true, Some(&supp(mate_pos)), true, Some(TEST_READ_CIGAR),
    );
    set_second_in_pair(&mut read1_supp1);

    let mut read1_supp2 = create_sam_record(
        &read_id1, CHR_1, supp1_pos1, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, true, Some(&supp(read_pos)), true, Some(TEST_READ_CIGAR),
    );

    // note last 2 regions are out of order
    region_start = 5000;
    partition_reader.setup_region(CHR_1, region_start, region_start + 999);
    partition_reader.process_read(&mut read1_supp1);
    partition_reader.post_process_region();

    region_start = 4000;
    partition_reader.setup_region(CHR_1, region_start, region_start + 999);
    partition_reader.process_read(&mut read1_supp2);
    partition_reader.post_process_region();

    assert_eq!(8, partition_reader.non_consensus_write_count());
    assert_eq!(2, partition_reader.consensus_write_count());
}

/// PartitionReaderTest.testUnmappedMateRead
#[test]
fn test_unmapped_mate_read() {
    let env = unmapping_env();
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();

    // mate is unmapped at original location by RegionUnmapper - simulated by changing its coords and setting the unmapped attribute
    // it is then processed twice at the new location
    let read_name = ids.next_id();
    let mut read =
        create_sam_record(&read_name, CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_3, 500, false, false, None, true, Some(TEST_READ_CIGAR));

    // in an unmapping region
    let mut mate =
        create_sam_record(&read_name, CHR_3, 500, &bases, TEST_READ_CIGAR, CHR_1, 100, true, false, None, false, Some(TEST_READ_CIGAR));

    let mut mate_original = clone_sam_record(&mate, &read_name);

    // mate processed by RegionUnmapper at original location
    check_transform_read(&mut mate, CHR_3);
    assert!(read_unmapped(&mate));
    assert!(!mate_unmapped(&mate));

    // mate processed again at original location by partition reader
    partition_reader.setup_region(CHR_3, 1, 1000);
    assert!(partition_reader.process_read(&mut mate_original));
    assert!(read_unmapped(&mate_original));

    partition_reader.post_process_region();

    // confirm was unmapped again but not written to BAM (the Java test checks the writer of its test class here,
    // which is not the one this reader writes to)
    assert_eq!(0, partition_reader.non_consensus_write_count());

    // read and now-unmapped mate processed at new location
    partition_reader.setup_region(CHR_1, 1, 1000);
    partition_reader.process_read(&mut read);
    partition_reader.process_read(&mut mate);

    partition_reader.post_process_region();

    // the reader keeps both records, so their state is checked as written
    assert_eq!(2, partition_reader.non_consensus_write_count());
    assert!(mate_unmapped(partition_reader.written(&read_name, 0, FLAG_UNMAPPED)));
    assert!(read_unmapped(partition_reader.written(&read_name, FLAG_UNMAPPED, 0)));
}

/// PartitionReaderTest.testFullyUnmappedReads
#[test]
fn test_fully_unmapped_reads() {
    // scenario 1 - both reads go from mapped to unmapped and are dropped by the RegionUnmapped and the PartitionReader
    let env = unmapping_env();
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();

    let mut read_name = ids.next_id();
    let mut read =
        create_sam_record(&read_name, CHR_3, 500, &bases, TEST_READ_CIGAR, CHR_4, 1000, false, false, None, true, Some(TEST_READ_CIGAR));

    // in an unmapping region
    let mut mate =
        create_sam_record(&read_name, CHR_4, 1000, &bases, TEST_READ_CIGAR, CHR_3, 500, true, false, None, false, Some(TEST_READ_CIGAR));

    let mut read_original = clone_sam_record(&read, &read_name);
    let mut mate_original = clone_sam_record(&mate, &read_name);

    // reads are processed by RegionUnmapper at original location
    check_transform_read(&mut read, CHR_3);
    assert!(read_unmapped(&read));
    assert!(mate_unmapped(&read));

    check_transform_read(&mut mate, CHR_4);
    assert!(read_unmapped(&mate));
    assert!(mate_unmapped(&mate));

    // reads are processed again at original location by partition reader
    partition_reader.setup_region(CHR_3, 1, 1000);
    assert!(partition_reader.process_read(&mut read_original));
    assert!(read_unmapped(&read_original));
    assert!(mate_unmapped(&read_original));
    partition_reader.post_process_region();

    partition_reader.setup_region(CHR_4, 1, 4000);
    assert!(partition_reader.process_read(&mut mate_original));
    assert!(read_unmapped(&mate_original));
    assert!(mate_unmapped(&mate_original));
    partition_reader.post_process_region();

    // confirm was unmapped again but not written to BAM
    assert_eq!(0, partition_reader.non_consensus_write_count());

    // scenario 2 - one read is already unmapped and its mate also becomes unmapped
    read_name = ids.next_id();
    read = create_sam_record(&read_name, CHR_3, 500, &bases, TEST_READ_CIGAR, CHR_4, 1000, false, false, None, false, Some(TEST_READ_CIGAR));
    rec::set_flag(&mut read, FLAG_MATE_UNMAPPED, true);
    assert!(mate_unmapped(&read));

    // in an unmapping region
    mate = create_sam_record(&read_name, CHR_3, 500, &bases, NO_CIGAR, CHR_3, 500, false, false, None, false, Some(TEST_READ_CIGAR));
    rec::set_flag(&mut mate, FLAG_UNMAPPED, true);

    read_original = clone_sam_record(&read, &read_name);
    mate_original = clone_sam_record(&mate, &read_name);

    // reads are processed by RegionUnmapper at original location
    check_transform_read(&mut read, CHR_3);
    assert!(read_unmapped(&read));
    assert!(mate_unmapped(&read));

    check_transform_read(&mut mate, CHR_3);
    assert!(read_unmapped(&mate));
    assert!(mate_unmapped(&mate));

    // reads are processed again at original location by partition reader
    partition_reader.setup_region(CHR_3, 1, 1000);
    assert!(partition_reader.process_read(&mut read_original));
    assert!(partition_reader.process_read(&mut mate_original));
    assert!(read_unmapped(&read_original));
    assert!(mate_unmapped(&read_original));
    assert!(read_unmapped(&mate_original));
    assert!(mate_unmapped(&mate_original));
    partition_reader.post_process_region();

    // confirm was unmapped again but not written to BAM
    assert_eq!(0, partition_reader.non_consensus_write_count());
}

/// Processes one duplicate group of the read names given and returns the CR attribute of its consensus read
fn consensus_read_attribute(env: &TestEnv, read_names: &[&str]) -> String {
    let mut partition_reader = env.partition_reader();
    partition_reader.setup_region(CHR_1, 1, 1000);
    for read_name in read_names {
        let mut read = create_sam_record(
            read_name, CHR_1, 100, &test_read_bases(), TEST_READ_CIGAR, CHR_1, 1000, false, false, None, false, Some(TEST_READ_CIGAR),
        );
        partition_reader.process_read(&mut read);
    }
    partition_reader.post_process_region();

    let consensus_reads: Vec<&bam::Record> = partition_reader.written_records().iter().filter(|r| rec::has_aux(r, b"CR")).collect();
    assert_eq!(1, consensus_reads.len());
    string_attribute(consensus_reads[0], b"CR").unwrap()
}

/// PartitionReaderTest.testIlluminaPCRClusterCount
#[test]
fn test_illumina_pcr_cluster_count() {
    let env = consensus_env();

    // invalid read names
    assert_eq!("2;2", consensus_read_attribute(&env, &["READ_001", "READ_002"]));

    // two optical duplicates
    assert_eq!("2;2;1", consensus_read_attribute(&env, &["A00624:8:HHKYHDSXX:4:2140:1:1", "A00624:8:HHKYHDSXX:4:2141:1:1"]));

    // two pcr duplicates
    assert_eq!("2;2;2", consensus_read_attribute(&env, &["A00624:8:HHKYHDSXX:4:2140:1:1", "A00624:8:HHKYHDSXX:5:2140:1:1"]));

    // three duplicates on different tiles
    assert_eq!(
        "3;3;3",
        consensus_read_attribute(
            &env,
            &["A00624:8:HHKYHDSXX:4:2140:1:1", "A00624:8:HHKYHDSXX:4:2141:1:1", "A00624:8:HHKYHDSXX:4:2142:1:1"],
        )
    );

    // seven duplicates with five pcr clusters
    assert_eq!(
        "7;7;5",
        consensus_read_attribute(
            &env,
            &[
                "A00624:8:HHKYHDSXX:4:2140:1:1",
                "A00624:8:HHKYHDSXX:4:2140:1:1",
                "A00624:8:HHKYHDSXX:4:2141:1:1",
                "A00624:8:HHKYHDSXX:4:2141:3001:1",
                "A00624:8:HHKYHDSXX:4:2142:1:1",
                "A00624:8:HHKYHDSXX:4:2142:1:1",
                "A00624:8:HHKYHDSXX:4:2142:3001:1",
            ],
        )
    );
}

/// Not a Java test. PartitionReader.checkMissingMateCigar: with -fail_supp_no_mate_cigar a paired primary read whose
/// mate is mapped must carry the mate cigar, else the run fails.
#[test]
#[ignore = "fails until /tmp/redux-rs-tests-fixes.patch is applied: the port accepts -fail_supp_no_mate_cigar but does not act on it"]
fn test_missing_mate_cigar_fails() {
    let mut env = TestEnv::new(ref_genome(), false, read_unmapper(HashMap::new()));
    env.config.fail_on_missing_supp_mate_cigar = true;
    let mut partition_reader = env.partition_reader();
    let bases = test_read_bases();
    partition_reader.setup_region(CHR_1, 1, 1000);

    // with a mate cigar, or without a mapped mate, or not a primary read: accepted
    let mut with_cigar = create_sam_record("READ_001", CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 300, false, false, None, true, Some(TEST_READ_CIGAR));
    assert!(partition_reader.reader.process_sam_record(&mut with_cigar).is_ok());
    let mut mate_unmapped = create_sam_record("READ_002", CHR_1, 110, &bases, TEST_READ_CIGAR, NO_CHROMOSOME_NAME, 0, false, false, None, false, None);
    assert!(partition_reader.reader.process_sam_record(&mut mate_unmapped).is_ok());
    let supp = SuppData::new(CHR_1, 500, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
    let mut supplementary = create_sam_record("READ_003", CHR_1, 120, &bases, TEST_READ_CIGAR, CHR_1, 300, false, true, Some(&supp), true, None);
    assert!(partition_reader.reader.process_sam_record(&mut supplementary).is_ok());

    let mut without_cigar = create_sam_record("READ_004", CHR_1, 130, &bases, TEST_READ_CIGAR, CHR_1, 300, false, false, None, true, None);
    let result = partition_reader.reader.process_sam_record(&mut without_cigar);
    assert!(result.is_err());
    assert!(format!("{:#}", result.unwrap_err()).contains("missing mate CIGAR"));
}
