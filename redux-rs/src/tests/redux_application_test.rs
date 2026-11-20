//! Port of ReduxApplicationTest.

use super::util::*;
use crate::rec::{self, FLAG_MATE_UNMAPPED, FLAG_UNMAPPED};
use crate::regions::{split_regions_into_partitions, ChrBaseRegion, SpecificRegions};
use crate::unmap::{UnmappingRegion, UNMAP_MAX_NON_OVERLAPPING_BASES, UNMAP_MIN_HIGH_DEPTH};
use std::collections::HashMap;

/// ReduxApplicationTest.testNonDuplicateFragment
#[test]
fn test_non_duplicate_fragment() {
    let env = TestEnv::test_config();
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();

    partition_reader.setup_region(CHR_1, 1, 1000);

    // test non-duplicate fragments processing reads in various orders
    let read_pos = 100;
    let mate_pos = 200;
    let supp_pos = 2000;

    let read1_name = ids.next_id();
    let mut read1 = create_sam_record(
        &read1_name, CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false,
        Some(&SuppData::new(CHR_1, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.process_read(&mut read1);
    partition_reader.flush_read_positions();

    assert_eq!(1, partition_reader.non_consensus_write_count());

    let mut mate1 =
        create_sam_record(&read1_name, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false, None, false, Some(TEST_READ_CIGAR));
    set_second_in_pair(&mut mate1);

    partition_reader.process_read(&mut mate1);

    let mut supp1 = create_sam_record(
        &read1_name, CHR_1, supp_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, false, true,
        Some(&SuppData::new(CHR_1, read_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.post_process_region();

    partition_reader.setup_region(CHR_1, 1001, 2000);

    partition_reader.process_read(&mut supp1);
    partition_reader.post_process_region();

    assert_eq!(3, partition_reader.non_consensus_write_count());
}

/// ReduxApplicationTest.testSupplementaryFragments
#[test]
fn test_supplementary_fragments() {
    let env = TestEnv::test_config();
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();

    partition_reader.setup_region(CHR_1, 1, 1000);

    // first a supplementary linked to a distant mate but in the same partition as the primary
    let read_pos = 100;
    let mate_pos = 1800;
    let supp_pos = 200;

    let read1_name = ids.next_id();
    let mut read1 =
        create_sam_record(&read1_name, CHR_1, read_pos, &bases, TEST_READ_CIGAR, CHR_1, mate_pos, false, false, None, true, Some(TEST_READ_CIGAR));

    partition_reader.process_read(&mut read1);

    let mut supp1 = create_sam_record(
        &read1_name, CHR_1, supp_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, false, true,
        Some(&SuppData::new(CHR_1, mate_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), true, Some(TEST_READ_CIGAR),
    );

    partition_reader.process_read(&mut supp1);
    partition_reader.post_process_region();

    partition_reader.setup_region(CHR_1, 1001, 2000);

    let mut mate1 = create_sam_record(
        &read1_name, CHR_1, mate_pos, &bases, TEST_READ_CIGAR, CHR_1, read_pos, true, false,
        Some(&SuppData::new(CHR_1, supp_pos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), false, Some(TEST_READ_CIGAR),
    );
    set_second_in_pair(&mut mate1);

    partition_reader.process_read(&mut mate1);
    partition_reader.post_process_region();

    assert_eq!(3, partition_reader.non_consensus_write_count());
}

/// ReduxApplicationTest.testUnmapRegionReads
#[test]
fn test_unmap_region_reads() {
    let mut chr_locations_map = HashMap::new();
    let region1 = UnmappingRegion { start: 550, end: 650, max_depth: UNMAP_MIN_HIGH_DEPTH + 1 };
    let region2 = UnmappingRegion { start: 900, end: 1300, max_depth: UNMAP_MIN_HIGH_DEPTH + 1 };
    chr_locations_map.insert(CHR_1.to_string(), vec![region1, region2]);
    let read_unmapper = read_unmapper(chr_locations_map);

    let bases_repeat = ref_bases_repeat_40();
    let mut ref_genome = MockRefGenome::new(false);
    ref_genome.chromosome_lengths.insert(CHR_1.to_string(), bases_repeat.len() as i32);
    ref_genome.ref_genome_map.insert(CHR_1.to_string(), bases_repeat);

    let env = TestEnv::new(ref_genome, false, read_unmapper);
    let mut partition_reader = env.partition_reader();
    let ids = &mut ReadIdGenerator::new();
    let bases = test_read_bases();

    partition_reader.setup_region(CHR_1, 1, 1999);

    // first fragments have supps on chr 1, primaries in the excluded regions and mates on chr 3

    // a read not overlapping enough with a region to be unmapped but supp is unmapped
    let read1_name = ids.next_id();
    let mut read1 = create_sam_record(
        &read1_name, CHR_1, 550 - UNMAP_MAX_NON_OVERLAPPING_BASES - 2, &bases, TEST_READ_CIGAR, CHR_3, 800, false, false,
        Some(&SuppData::new(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1)), false, Some(TEST_READ_CIGAR),
    );

    // the reader keeps this read, so its state is checked once written, at the end
    assert!(!partition_reader.process_read(&mut read1));

    // read overlaps so that the non overlapping bases of the aligned part does not exceed the maximum
    let mut read2 = create_sam_record(
        &ids.next_id(), CHR_1, 551 + UNMAP_MAX_NON_OVERLAPPING_BASES - 1, &bases, TEST_READ_CIGAR, CHR_3, 800, false, false, None, false,
        Some(TEST_READ_CIGAR),
    );

    partition_reader.process_read(&mut read2);
    assert!(rec::has(&read2, FLAG_UNMAPPED));

    // mate is unmapped
    let read3_name = ids.next_id();
    let mut read3 = create_sam_record(
        &read3_name, CHR_1, 100, &bases, TEST_READ_CIGAR, CHR_1, 550 - UNMAP_MAX_NON_OVERLAPPING_BASES + 1, false, false, None, false,
        Some(TEST_READ_CIGAR),
    );

    assert!(!partition_reader.process_read(&mut read3));

    // unmapped read then has mate unmapped
    let mut read4 =
        create_sam_record(&ids.next_id(), CHR_1, 1000, &bases, NO_CIGAR, CHR_1, 1000, false, false, None, false, Some(TEST_READ_CIGAR));
    rec::set_flag(&mut read4, FLAG_UNMAPPED, true);

    partition_reader.process_read(&mut read4);
    assert!(rec::has(&read4, FLAG_UNMAPPED));
    assert!(rec::has(&read4, FLAG_MATE_UNMAPPED));

    // not in the Java test: the reads kept by the reader are written out to check them
    partition_reader.post_process_region();

    let read1 = partition_reader.written(&read1_name, 0, 0);
    assert_eq!(None, string_attribute(read1, b"SA"));
    assert!(!rec::has(read1, FLAG_UNMAPPED));

    let read3 = partition_reader.written(&read3_name, 0, 0);
    assert!(!rec::has(read3, FLAG_UNMAPPED));
    assert!(rec::has(read3, FLAG_MATE_UNMAPPED));
}

/// ReduxApplicationTest.testRegionPartitionSplits
///
/// The Java test passes no reference genome, which takes the chromosome lengths from the built-in v37 coordinates;
/// here the lengths come from the index of the placeholder reference.
#[test]
fn test_region_partition_splits() {
    let genome = crate::refgenome::RefGenome::open(placeholder_fasta()).unwrap();

    let mut specific_regions = SpecificRegions::default();
    specific_regions.regions.push(ChrBaseRegion::new(CHR_1, 1, 1000000));
    specific_regions.chromosomes.push(CHR_1.to_string());

    let mut partition_regions = split_regions_into_partitions(&specific_regions, 4, false, &genome, true);
    assert_eq!(4, partition_regions.len());

    let mut region = &partition_regions[0][0];
    assert_eq!(1, region.start);
    assert_eq!(250000, region.end);

    region = &partition_regions[1][0];
    assert_eq!(250001, region.start);
    assert_eq!(500000, region.end);

    region = &partition_regions[2][0];
    assert_eq!(500001, region.start);
    assert_eq!(750000, region.end);

    region = &partition_regions[3][0];
    assert_eq!(750001, region.start);
    assert_eq!(1000000, region.end);

    specific_regions = SpecificRegions::default();
    for chromosome in [CHR_1, CHR_2, CHR_3] {
        specific_regions.regions.push(ChrBaseRegion::new(chromosome, 1, 10000));
        specific_regions.chromosomes.push(chromosome.to_string());
    }

    partition_regions = split_regions_into_partitions(&specific_regions, 1, false, &genome, true);
    assert_eq!(1, partition_regions.len());
    assert_eq!(specific_regions.regions[0], partition_regions[0][0]);
    assert_eq!(specific_regions.regions[1], partition_regions[0][1]);
    assert_eq!(specific_regions.regions[2], partition_regions[0][2]);

    // test with human genome
    specific_regions.regions.clear();
    specific_regions.chromosomes.clear();
    specific_regions.chromosomes.push(CHR_1.to_string());
    specific_regions.chromosomes.push(CHR_2.to_string());

    partition_regions = split_regions_into_partitions(&specific_regions, 3, false, &genome, true);
    assert_eq!(3, partition_regions.len());

    region = &partition_regions[0][0];
    assert_eq!(CHR_1, region.chromosome);
    assert_eq!(1, region.start);
    assert_eq!(164149998, region.end);

    region = &partition_regions[1][0];
    assert_eq!(CHR_1, region.chromosome);
    assert_eq!(164149999, region.start);
    assert_eq!(249250621, region.end);

    region = &partition_regions[1][1];
    assert_eq!(CHR_2, region.chromosome);
    assert_eq!(1, region.start);
    assert_eq!(79049375, region.end);

    region = &partition_regions[2][0];
    assert_eq!(CHR_2, region.chromosome);
    assert_eq!(79049376, region.start);
    assert_eq!(243199373, region.end);
}
