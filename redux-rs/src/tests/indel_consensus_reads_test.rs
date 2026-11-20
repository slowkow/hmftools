//! Port of consensus/IndelConsensusReadsTest.

use super::consensus_reads_test::{char_at, create_consensus_input_record, delete_char_at, next_base, next_umi_read_id, set_char_at, UMI_ID_1};
use super::util::*;
use crate::consensus::ConsensusOutcome::*;
use crate::consensus::ConsensusReads;
use crate::rec;
use rust_htslib::bam;

struct Fixture {
    consensus_reads: ConsensusReads,
    ids: ReadIdGenerator,
    ref_bases: String,
}

impl Fixture {
    fn new() -> Self {
        let ref_bases = ref_bases();
        let mut ref_genome = MockRefGenome::new(false);
        ref_genome.ref_genome_map.insert(CHR_1.to_string(), ref_bases.clone());
        Fixture { consensus_reads: consensus_reads(&ref_genome), ids: ReadIdGenerator::new(), ref_bases }
    }
    fn next_read_id(&mut self) -> String {
        next_umi_read_id(UMI_ID_1, &mut self.ids)
    }
    /// IndelConsensusReadsTest.createRecord
    fn record(&mut self, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool) -> bam::Record {
        let id = self.next_read_id();
        create_sam_record(&id, CHR_1, read_start, read_bases, cigar, CHR_1, 5000, is_reversed, false, None, !is_reversed, Some(cigar))
    }
    fn input(&mut self, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool) -> bam::Record {
        let id = self.next_read_id();
        create_consensus_input_record(&id, read_start, read_bases, cigar, is_reversed)
    }
}

/// IndelConsensusReadsTest.testMatchingIndelReads
#[test]
fn test_matching_indel_reads() {
    // indels with no CIGAR differences use the standard consensus routines
    let mut f = Fixture::new();
    let mut reads = Vec::new();

    let pos_start = 13;

    let consensus_bases = format!("{}C{}", &REF_BASES_A[0..5], &REF_BASES_A[5..10]);
    let indel_cigar = "2S3M1I3M2S";
    reads.push(f.input(pos_start, &consensus_bases, indel_cigar, false));
    reads.push(f.input(pos_start, &consensus_bases, indel_cigar, false));
    reads.push(f.input(pos_start, &consensus_bases, indel_cigar, false));

    let read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(IndelMatch, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!(indel_cigar, cigar_string(&read_info.consensus_read));
    assert_eq!(pos_start, rec::start(&read_info.consensus_read));
}

/// IndelConsensusReadsTest.testInsertMismatches
#[test]
fn test_insert_mismatches() {
    // first 2 reads without an insert vs 1 with one ignored
    let mut f = Fixture::new();
    let a = REF_BASES_A;
    let mut reads_reversed = false;
    let mut read1 = f.record(11, &a[0..10], "10M", reads_reversed);

    // matching but without the soft-clip
    let mut read2 = f.record(11, &a[0..10], "10M", reads_reversed);

    let mut indel_bases = format!("{}C{}", &a[0..5], &a[5..10]);
    let mut read3 = f.record(12, &indel_bases, "1S4M1I4M1S", reads_reversed);

    indel_bases = format!("{}GG{}", &a[0..3], &a[3..10]);
    let mut read4 = f.record(11, &indel_bases, "3M2I7M", reads_reversed);

    indel_bases = format!("{}TTT{}TTTT{}", &a[0..4], &a[4..7], &a[7..10]);
    let read5 = f.record(12, &indel_bases, "1S3M3I3M4I2M1S", reads_reversed);

    let mut consensus_bases = a[0..10].to_string();

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3, read4, read5], UMI_ID_1);
    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!("10M", cigar_string(&read_info.consensus_read));
    assert_eq!(11, rec::start(&read_info.consensus_read));

    // now with 2 inserts in the consensus
    consensus_bases = format!("{}TT{}GGG{}", &a[0..3], &a[3..6], &a[6..10]);
    let consensus_cigar = "1S2M2I3M3I3M1S";

    reads_reversed = true;

    // unclipped 11, start = 12, end = 19, unclipped end = 20
    read1 = f.record(12, &consensus_bases, consensus_cigar, reads_reversed);
    read2 = clone_sam_record(&read1, &f.next_read_id());

    // unclipped 10, start = 12, end = 17, unclipped end = 19
    read3 = f.record(12, &a[0..10], "2S6M2S", reads_reversed);

    // has the same indels but shorter
    // unclipped 11, start = 12, end = 19, unclipped end = 20
    indel_bases = format!("{}TTT{}GG{}", &a[0..3], &a[3..6], &a[6..10]);
    read4 = f.record(12, &indel_bases, "1S2M3I3M2I3M1S", reads_reversed);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3, read4], UMI_ID_1);
    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!(consensus_cigar, cigar_string(&read_info.consensus_read));
    assert_eq!(12, rec::start(&read_info.consensus_read));
}

/// IndelConsensusReadsTest.testDeleteMismatches
#[test]
fn test_delete_mismatches() {
    // first 2 reads without a delete vs 1 with one ignored
    let mut f = Fixture::new();
    let a = REF_BASES_A;
    let mut read1 = f.record(11, &a[0..10], "10M", false);

    let mut read2 = f.record(11, &a[0..10], "10M", false);

    let mut indel_bases = format!("{}{}", &a[0..5], &a[6..10]);
    let mut read3 = f.record(12, &indel_bases, "1S4M1D3M1S", false);

    indel_bases = format!("{}{}", &a[0..3], &a[7..10]);
    let read4 = f.record(11, &indel_bases, "3M4D3M", false);

    indel_bases = format!("{}{}{}", &a[0..3], &a[4..6], &a[7..10]);
    let read5 = f.record(11, &indel_bases, "3M1D2M1D3M", false);

    let mut consensus_bases = a[0..10].to_string();

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3, read4, read5], UMI_ID_1);
    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!("10M", cigar_string(&read_info.consensus_read));
    assert_eq!(11, rec::start(&read_info.consensus_read));

    // now test a delete in the consensus read
    let consensus_cigar = "4M3D3M";
    let del_bases = format!("{}{}", &a[0..4], &a[7..10]);
    read1 = f.record(11, &del_bases, consensus_cigar, false);
    set_base_qualities(&mut read1, 11);
    read2 = f.record(11, &del_bases, consensus_cigar, false);
    set_base_qualities(&mut read2, 11);

    indel_bases = format!("{}CCC", &a[0..7]); // no delete and differing aligned bases after the delete
    read3 = f.record(11, &indel_bases, "10M", false);
    set_base_qualities(&mut read3, 35);

    consensus_bases = format!("{}CCC", &a[0..4]);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);
    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!(consensus_cigar, cigar_string(&read_info.consensus_read));

    let calc_base_qual = 11; // lowered from 13

    assert_eq!(calc_base_qual, read_info.consensus_read.qual()[4]);
    assert_eq!(calc_base_qual, read_info.consensus_read.qual()[5]);
    assert_eq!(calc_base_qual, read_info.consensus_read.qual()[6]);

    // a basic test for following bases correctly
    let read_bases = "CTTCGATAATGGCCGGGCCG";
    let shorted_read_bases = &read_bases[0..read_bases.len() - 2];
    let shorter_clipped = "3S7M10D5M3S";
    read1 = f.record(8, shorted_read_bases, shorter_clipped, false); // unclipped is 7
    read2 = f.record(10, read_bases, "5S5M10D5M5S", false); // unclipped is 5

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2], UMI_ID_1);
    assert_eq!(shorter_clipped, cigar_string(&read_info.consensus_read));
    assert_eq!(shorted_read_bases, read_string(&read_info.consensus_read));
    assert_eq!(8, rec::start(&read_info.consensus_read));
}

/// IndelConsensusReadsTest.testIndelMismatchedTemplateSelection
#[test]
fn test_indel_mismatched_template_selection() {
    // favour the read with the least number of soft-clip bases
    let mut f = Fixture::new();
    let mut read1 = f.record(10, "CTTCGATAAT", "7M1D3M", true);
    let mut read2 = f.record(13, "TTCGATAAAT", "2S8M", true);

    // they are marked as same read cause read is reversed and the end are the same
    assert_eq!(rec::end(&read1), rec::end(&read2));

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &[read1.clone(), read2.clone()], UMI_ID_1);
    assert_eq!("7M1D3M", cigar_string(&read_info.consensus_read));
    assert_eq!("CTTCGATAAT", read_string(&read_info.consensus_read));
    assert_eq!(10, rec::start(&read_info.consensus_read));

    // a second supporting read results in the extra base being dropped
    let mut read3 = f.record(13, "TTCGATAAAT", "2S8M", true);
    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);
    assert_eq!("2S8M", cigar_string(&read_info.consensus_read));
    assert_eq!("TTCGATAAAT", read_string(&read_info.consensus_read));
    assert_eq!(13, rec::start(&read_info.consensus_read));

    read1 = f.record(10, "CTTCGATAAAAT", "7M1I4M", true);
    read2 = f.record(13, "TTCGATAAAT", "2S8M", true);

    // they are marked as same read cause read is reversed and the end are the same
    assert_eq!(rec::end(&read1), rec::end(&read2));

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1.clone(), read2.clone()], UMI_ID_1);
    assert_eq!("7M1I4M", cigar_string(&read_info.consensus_read));
    assert_eq!("CTTCGATAAAAT", read_string(&read_info.consensus_read));
    assert_eq!(10, rec::start(&read_info.consensus_read));

    // add one more read so the first C will get chopped
    read3 = f.record(13, "TTCGATAAAT", "2S8M", true);
    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);
    assert_eq!("2S8M", cigar_string(&read_info.consensus_read));
    assert_eq!("TTCGATAAAT", read_string(&read_info.consensus_read));
    assert_eq!(13, rec::start(&read_info.consensus_read));
}

/// IndelConsensusReadsTest.testIndelMismatchesHardClips
#[test]
fn test_indel_mismatches_hard_clips() {
    let mut f = Fixture::new();
    let mut read_bases = "CTTCGATAATGGCCG";
    let mut read1 = f.record(10, read_bases, "3H8M7S", false);
    let mut read2 = f.record(15, read_bases, "5S10M5H", false); // unclipped start = 10
    let mut read3 = f.record(10, read_bases, "2H13M2S", false);

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);
    assert_eq!("2H13M2S", cigar_string(&read_info.consensus_read));
    assert_eq!(read_bases, read_string(&read_info.consensus_read));
    assert_eq!(10, rec::start(&read_info.consensus_read));

    read_bases = "CTTCGATAATGGCCGGGCCG";
    read1 = f.record(10, read_bases, "5H10M10D10M10H", false);
    read2 = f.record(10, read_bases, "5H10M10D10M10H", false);
    read3 = f.record(15, read_bases, "10H5M10D10M10H", false);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);
    assert_eq!("5H10M10D10M10H", cigar_string(&read_info.consensus_read));
    assert_eq!(read_bases, read_string(&read_info.consensus_read));
    assert_eq!(10, rec::start(&read_info.consensus_read));
}

/// IndelConsensusReadsTest.testSoftClipDeletes
#[test]
fn test_soft_clip_deletes() {
    // test combinations of soft-clips vs aligned bases
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let read1 = f.record(40, &ref_bases[40..100], "30S30M", true);

    let read2 = f.record(7, &format!("{}{}", &ref_bases[37..67], &ref_bases[70..100]), "30M3D30M", true);

    let read3 = f.record(40, &ref_bases[40..100], "30S30M", true);
    let read4 = f.record(39, &ref_bases[40..100], "29S31M", true);

    let read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3, read4], UMI_ID_1);
    assert_eq!("30S30M", cigar_string(&read_info.consensus_read)); // was 30M2D31M when favouring non-SCs
    assert_eq!(40, rec::start(&read_info.consensus_read));
    assert_eq!(69, rec::end(&read_info.consensus_read));

    let consensus_bases = &ref_bases[40..100];
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
}

/// IndelConsensusReadsTest.testIndelMismatchNumMutations
#[test]
fn test_indel_mismatch_num_mutations() {
    // no mutations
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let pos_start = 11;
    let read_length = 10;
    let consensus_del_index = 2;
    let consensus_bases = delete_char_at(&ref_bases[pos_start..pos_start + read_length + 1], consensus_del_index);
    let pos = pos_start as i32;

    let mut read1 = f.input(pos, &consensus_bases, "1M1D9M", false);
    let mut read2 = f.input(pos, &consensus_bases, "2M1D8M", false);
    set_second_in_pair(&mut read2);
    let mut read3 = f.input(pos, &consensus_bases, "2M1D8M", false);
    set_second_in_pair(&mut read3);

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);

    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(Some(1), rec::aux_int(&read_info.consensus_read, b"NM"));

    // mismatch
    let mutated_base_index = 8;
    let mutated_bases = set_char_at(&consensus_bases, mutated_base_index, next_base(char_at(&consensus_bases, mutated_base_index)));

    read1 = f.input(pos, &mutated_bases, "1M1D9M", false);
    read2 = f.input(pos, &mutated_bases, "2M1D8M", false);
    set_second_in_pair(&mut read2);
    read3 = f.input(pos, &mutated_bases, "2M1D8M", false);
    set_second_in_pair(&mut read3);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);

    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(Some(2), rec::aux_int(&read_info.consensus_read, b"NM"));
}
