//! Port of consensus/ConsensusReadsTest.

use super::util::*;
use crate::consensus::ConsensusOutcome::*;
use crate::consensus::{form_consensus_read_id, ConsensusReads, BASE_QUAL_MINIMUM, CONSENSUS_MAX_DEPTH};
use crate::rec::{self, FLAG_MATE_UNMAPPED, FLAG_UNMAPPED};
use rust_htslib::bam;

pub const UMI_ID_1: &str = "TAGTAG";

/// ConsensusReadsTest.NEXT_BASE_MAP
pub fn next_base(base: char) -> char {
    match base {
        'G' => 'C',
        'C' => 'A',
        'A' => 'T',
        'T' => 'G',
        other => panic!("no next base for {other}"),
    }
}

/// StringBuilder.setCharAt
pub fn set_char_at(bases: &str, index: usize, base: char) -> String {
    let mut chars: Vec<char> = bases.chars().collect();
    chars[index] = base;
    chars.into_iter().collect()
}

/// StringBuilder.deleteCharAt
pub fn delete_char_at(bases: &str, index: usize) -> String {
    let mut chars: Vec<char> = bases.chars().collect();
    chars.remove(index);
    chars.into_iter().collect()
}

pub fn char_at(bases: &str, index: usize) -> char {
    bases.as_bytes()[index] as char
}

/// ConsensusReadsTest.nextUmiReadId
pub fn next_umi_read_id(umi_id: &str, ids: &mut ReadIdGenerator) -> String {
    format!("ABCD:{}:{}", ids.next_id(), umi_id)
}

/// ConsensusReadsTest.createConsensusInputRecord
pub fn create_consensus_input_record(read_id: &str, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool) -> bam::Record {
    create_sam_record(read_id, CHR_1, read_start, read_bases, cigar, CHR_1, 5000, is_reversed, false, None, true, Some(TEST_READ_CIGAR))
}

struct Fixture {
    ref_genome_one_based: MockRefGenome,
    consensus_reads: ConsensusReads,
    ids: ReadIdGenerator,
    ref_bases: String,
}

impl Fixture {
    fn new() -> Self {
        let ref_bases = ref_bases();
        let mut ref_genome = MockRefGenome::new(false);
        ref_genome.ref_genome_map.insert(CHR_1.to_string(), ref_bases.clone());
        ref_genome.chromosome_lengths.insert(CHR_1.to_string(), ref_bases.len() as i32);
        let consensus_reads = consensus_reads(&ref_genome);

        let mut ref_genome_one_based = MockRefGenome::new(true);
        ref_genome_one_based.ref_genome_map.insert(CHR_1.to_string(), ref_bases.clone());
        ref_genome_one_based.chromosome_lengths.insert(CHR_1.to_string(), ref_bases.len() as i32);

        Fixture { ref_genome_one_based, consensus_reads, ids: ReadIdGenerator::new(), ref_bases }
    }
    fn next_read_id(&mut self) -> String {
        next_umi_read_id(UMI_ID_1, &mut self.ids)
    }
    fn input(&mut self, read_start: i32, read_bases: &str, cigar: &str, is_reversed: bool) -> bam::Record {
        let id = self.next_read_id();
        create_consensus_input_record(&id, read_start, read_bases, cigar, is_reversed)
    }
}

/// ConsensusReadsTest.testBasicConsensusReads
#[test]
fn test_basic_consensus_reads() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let mut pos_start = 11;
    let consensus_bases = &ref_bases[pos_start as usize..21];
    reads.push(f.input(pos_start, consensus_bases, "10M", false));
    reads.push(f.input(pos_start, consensus_bases, "10M", false));
    reads.push(f.input(pos_start, consensus_bases, "10M", false));

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!("10M", cigar_string(&read_info.consensus_read));
    assert_eq!(pos_start, rec::start(&read_info.consensus_read));

    // mismatch at some of the bases and so determined by qual
    reads.clear();

    let mut read1 = f.input(pos_start, consensus_bases, "10M", false);
    set_base_qualities(&mut read1, DEFAULT_QUAL - 1);
    reads.push(read1);

    let bases2 = "AATAATAATA";
    let read2 = f.input(pos_start + 2, bases2, "2S8M", false);
    reads.push(read2);

    let bases3 = "AATAAGAAAA";
    let mut read3 = f.input(pos_start + 4, bases3, "4S6M", false);
    set_base_qualities(&mut read3, DEFAULT_QUAL - 1);
    reads.push(read3);

    // 3rd base is 'T' due to max qual of 2nd and 3rd reads
    // 6th base is 'T' due to max qual of 2nd read
    // 9th base is 'A' due to max qual of 1st and 3rd reads

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(pos_start, rec::start(&read_info.consensus_read));
    assert_eq!("10M", cigar_string(&read_info.consensus_read));

    assert_eq!("AATAATAAAA", read_string(&read_info.consensus_read));
    assert_eq!(11, read_info.consensus_read.qual()[2]); // lowered from 19
    assert_eq!(BASE_QUAL_MINIMUM, read_info.consensus_read.qual()[5] as i32);
    assert_eq!(11, read_info.consensus_read.qual()[8]);

    // test preference for reference base when quals are equal
    reads.clear();
    pos_start = 16;
    let read1 = f.input(pos_start, REF_BASES_A, "10M", false);
    reads.push(read1);

    let read2 = f.input(pos_start, REF_BASES_C, "10M", false);
    reads.push(read2);

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!("AAAAACCCCC", read_string(&read_info.consensus_read));

    // soft-clips at both ends
    pos_start = 11;
    let read1 = f.input(pos_start + 3, REF_BASES_A, "3S6M1S", false);
    let read2 = f.input(pos_start + 3, REF_BASES_A, "3S6M1S", false);

    let read3 = f.input(pos_start + 2, REF_BASES_A, "2S5M3S", false);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);

    assert_eq!(pos_start + 3, rec::start(&read_info.consensus_read));
    assert_eq!("3S6M1S", cigar_string(&read_info.consensus_read));
    assert_eq!(REF_BASES_A, read_string(&read_info.consensus_read));

    // longer reads at the non-fragment start end
    pos_start = 11;
    let longer = format!("{REF_BASES_A}CCC");
    let read1 = f.input(pos_start, REF_BASES_A, "8M2S", false);

    let read2 = f.input(pos_start, &longer, "7M6S", false);
    let read3 = f.input(pos_start, &longer, "7M6S", false);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);

    assert_eq!(pos_start, rec::start(&read_info.consensus_read));
    assert_eq!("7M6S", cigar_string(&read_info.consensus_read));
    assert_eq!(longer, read_string(&read_info.consensus_read));

    // same again on the reverse strand
    pos_start = 11;
    let longer = format!("CCCC{REF_BASES_A}");
    let read1 = f.input(pos_start + 2, REF_BASES_A, "2S6M2S", true);

    let read2 = f.input(pos_start, &longer, "4S7M3S", true);
    let read3 = f.input(pos_start, &longer, "4S7M3S", true);

    read_info = create_consensus_read(&mut f.consensus_reads, &[read1, read2, read3], UMI_ID_1);

    assert_eq!(pos_start, rec::start(&read_info.consensus_read));
    assert_eq!("4S7M3S", cigar_string(&read_info.consensus_read));
    assert_eq!(longer, read_string(&read_info.consensus_read));
}

/// ConsensusReadsTest.testMultipleSoftClipLengths
#[test]
fn test_multiple_soft_clip_lengths() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let pos_start = 1;

    let consensus_bases = &ref_bases[1..];

    let cigar1 = "100M";
    reads.push(f.input(pos_start, consensus_bases, cigar1, true));
    reads.push(f.input(pos_start, consensus_bases, cigar1, true));

    let cigar2 = "17S83M";
    reads.push(f.input(pos_start + 17, consensus_bases, cigar2, true));

    let cigar3 = "68S32M";
    reads.push(f.input(pos_start + 68, consensus_bases, cigar3, true));

    let cigar4 = "94S6M";
    reads.push(f.input(pos_start + 94, consensus_bases, cigar4, true));
    reads.push(f.input(pos_start + 94, consensus_bases, cigar4, true));
    reads.push(f.input(pos_start + 94, consensus_bases, cigar4, true));
    reads.push(f.input(pos_start + 94, consensus_bases, cigar4, true));

    let cigar5 = "99S1M";
    reads.push(f.input(pos_start + 99, consensus_bases, cigar5, true));

    let read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!(cigar4, cigar_string(&read_info.consensus_read));
    assert_eq!(pos_start + 94, rec::start(&read_info.consensus_read));
}

/// ConsensusReadsTest.testMateAlignmentConsensus
#[test]
fn test_mate_alignment_consensus() {
    // variable mate alignments from soft-clips, and if all else equal by sorted read ID
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let pos_start = 1;

    let consensus_bases = &ref_bases[1..];
    let read_cigar = "100M";

    let make = |f: &mut Fixture, mate_start: i32, mate_cigar: &str| {
        let id = f.next_read_id();
        let mut read = create_sam_record_basic(&id, CHR_1, pos_start, consensus_bases, read_cigar, CHR_1, mate_start, false, false, None);
        rec::set_flag(&mut read, rec::FLAG_MATE_REVERSE, false);
        set_string_attribute(&mut read, b"MC", mate_cigar);
        read
    };

    let mate_cigar1 = "100M";
    let read1 = make(&mut f, 1000, mate_cigar1);
    reads.push(read1);

    let mate_cigar2 = "10S90M";
    let read2 = make(&mut f, 1010, mate_cigar2);
    let read2b = make(&mut f, 1010, mate_cigar2);
    reads.push(read2);
    reads.push(read2b);

    let mate_cigar3 = "5S95M";
    let read3 = make(&mut f, 1005, mate_cigar3);
    let read3b = make(&mut f, 1005, mate_cigar3);
    let read3c = make(&mut f, 1005, mate_cigar3);

    // will be sorted alphabetically
    let consensus_read_id = form_consensus_read_id(read3.qname(), None);
    reads.push(read3c);
    reads.push(read3b);
    reads.push(read3);

    let read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(consensus_read_id, read_info.consensus_read.qname());
    assert_eq!(1005, rec::mate_start(&read_info.consensus_read));
    assert_eq!(Some(mate_cigar3.to_string()), string_attribute(&read_info.consensus_read, b"MC"));
}

/// ConsensusReadsTest.testSoftClippedOverChromosomeEnd
#[test]
fn test_soft_clipped_over_chromosome_end() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let pos_start = ref_bases.len() - 5;
    let read_bases1 = format!("{}{}", &ref_bases[pos_start..], "T".repeat(5));
    reads.push(f.input(pos_start as i32, &read_bases1, "5M5S", false));

    let read_bases2 = format!("{}{}", &ref_bases[pos_start..], "A".repeat(5));
    reads.push(f.input(pos_start as i32, &read_bases2, "5M5S", false));

    let mut consensus = consensus_reads(&f.ref_genome_one_based);
    let mut read_info = create_consensus_read(&mut consensus, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);

    consensus = consensus_reads(&f.ref_genome_one_based);
    consensus.set_chromosome_length(f.ref_genome_one_based.chromosome_lengths[CHR_1]);
    read_info = create_consensus_read(&mut consensus, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
}

/// ConsensusReadsTest.testSoftClippedOverChromosomeEndWithDel
#[test]
fn test_soft_clipped_over_chromosome_end_with_del() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    // ref:    CGTTAC
    // read 1: CGT AC
    // read 2: CG TAC
    // net:    CGTTAC

    let pos_start = ref_bases.len() as i32 - 5;
    let read_bases1 = "CGTACTTTTT";
    let read_cigar1 = "3M1D2M5S";
    reads.push(f.input(pos_start, read_bases1, read_cigar1, false));

    let read_bases2 = "CGTACAAAAA";
    let read_cigar2 = "2M1D3M5S";
    reads.push(f.input(pos_start, read_bases2, read_cigar2, false));

    let mut consensus = consensus_reads(&f.ref_genome_one_based);
    let mut read_info = create_consensus_read(&mut consensus, &reads, UMI_ID_1);
    assert_eq!(IndelMismatch, read_info.outcome);

    consensus = consensus_reads(&f.ref_genome_one_based);
    consensus.set_chromosome_length(f.ref_genome_one_based.chromosome_lengths[CHR_1]);
    read_info = create_consensus_read(&mut consensus, &reads, UMI_ID_1);
    assert_eq!(read_cigar2, cigar_string(&read_info.consensus_read));
    assert_eq!(read_bases1, read_string(&read_info.consensus_read)); // takes first if cannot determine ref position
    assert_eq!(IndelMismatch, read_info.outcome);
}

/// ConsensusReadsTest.testSoftClippedOverChromosomeStart
#[test]
fn test_soft_clipped_over_chromosome_start() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let read_bases1 = format!("{}{}", "T".repeat(5), &ref_bases[0..5]);
    reads.push(f.input(1, &read_bases1, "5S5M", false));

    let read_bases2 = format!("{}{}", "A".repeat(5), &ref_bases[0..5]);
    reads.push(f.input(1, &read_bases2, "5S5M", false));

    let mut consensus = consensus_reads(&f.ref_genome_one_based);
    let read_info = create_consensus_read(&mut consensus, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
}

/// ConsensusReadsTest.testSoftClippedOverChromosomeStartWithDel
#[test]
fn test_soft_clipped_over_chromosome_start_with_del() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let read_bases1 = format!("{}{}", "T".repeat(5), &ref_bases[0..5]);
    reads.push(f.input(1, &read_bases1, "5S2M1D3M", false));

    let read_bases2 = format!("{}{}", "A".repeat(5), &ref_bases[0..5]);
    reads.push(f.input(1, &read_bases2, "5S1M1D4M", false));

    let mut consensus = consensus_reads(&f.ref_genome_one_based);
    let read_info = create_consensus_read(&mut consensus, &reads, UMI_ID_1);
    assert_eq!(IndelMismatch, read_info.outcome);
}

/// ConsensusReadsTest.testDualStrandPrefersRefBase
#[test]
fn test_dual_strand_prefers_ref_base() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let pos_start = 11;
    let read_length = 10;
    let cigar = "10M";
    let consensus_bases = &ref_bases[pos_start..pos_start + read_length];
    let read1 = f.input(pos_start as i32, consensus_bases, cigar, false);

    let mutated_base_index = 5;
    let mutated_base = next_base(char_at(consensus_bases, mutated_base_index));
    let mutated_bases = set_char_at(consensus_bases, mutated_base_index, mutated_base);
    let mut read2 = f.input(pos_start as i32, &mutated_bases, cigar, false);
    set_second_in_pair(&mut read2);

    let mut read3 = f.input(pos_start as i32, &mutated_bases, cigar, false);
    set_second_in_pair(&mut read3);

    let reads = vec![read1, read2, read3];
    let read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);

    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!("10M", cigar_string(&read_info.consensus_read));
    assert_eq!(37, read_info.consensus_read.qual()[0] as i32); // non mutated base
    assert_eq!(BASE_QUAL_MINIMUM, read_info.consensus_read.qual()[mutated_base_index] as i32); // mutated base
    assert_eq!(pos_start as i32, rec::start(&read_info.consensus_read));
}

/// ConsensusReadsTest.testDualStrandPrefersRefBaseWithDel
#[test]
fn test_dual_strand_prefers_ref_base_with_del() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let pos_start = 11;
    let read_length = 10;
    let consensus_del_index = 2;
    let consensus_bases = delete_char_at(&ref_bases[pos_start..pos_start + read_length + 1], consensus_del_index);

    let read1 = f.input(pos_start as i32, &consensus_bases, "1M1D9M", false);

    let mutated_base_index = 5;
    let mutated_base = next_base(char_at(&consensus_bases, mutated_base_index));
    let mutated_bases = set_char_at(&consensus_bases, mutated_base_index, mutated_base);
    let mut read2 = f.input(pos_start as i32, &mutated_bases, "2M1D8M", false);
    set_second_in_pair(&mut read2);

    let mut read3 = f.input(pos_start as i32, &mutated_bases, "2M1D8M", false);
    set_second_in_pair(&mut read3);

    let reads = vec![read1, read2, read3];
    let read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);

    assert_eq!(IndelMismatch, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!("2M1D8M", cigar_string(&read_info.consensus_read));
    assert_eq!(pos_start as i32, rec::start(&read_info.consensus_read));
}

/// ConsensusReadsTest.testMaxDepthFiltering
#[test]
fn test_max_depth_filtering() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    let pos_start = 11;
    let read_length = 10;
    let consensus_cigar = "10M";
    let consensus_bases = &ref_bases[pos_start..pos_start + read_length];
    for _ in 0..CONSENSUS_MAX_DEPTH {
        reads.push(f.input(pos_start as i32, consensus_bases, consensus_cigar, false));
    }

    let mutated_base_index = 5;
    let mutated_bases = set_char_at(consensus_bases, mutated_base_index, next_base(char_at(consensus_bases, mutated_base_index)));
    for _ in 0..CONSENSUS_MAX_DEPTH + 1 {
        reads.push(f.input(pos_start as i32, &mutated_bases, consensus_cigar, false));
    }

    let read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(consensus_bases, read_string(&read_info.consensus_read));
    assert_eq!(consensus_cigar, cigar_string(&read_info.consensus_read));
    assert_eq!(pos_start as i32, rec::start(&read_info.consensus_read));
}

fn num_mutations(read: &bam::Record) -> i64 {
    rec::aux_int(read, b"NM").expect("NM attribute")
}

/// ConsensusReadsTest.testNumMutationsAttribute
#[test]
fn test_num_mutations_attribute() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let mut reads = Vec::new();

    // no mutations
    let pos_start = 11;
    let read_bases = &ref_bases[pos_start..pos_start + 10];
    let pos = pos_start as i32;
    reads.push(f.input(pos, read_bases, "10M", false));
    reads.push(f.input(pos, read_bases, "10M", false));

    let mut read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(0, num_mutations(&read_info.consensus_read));

    // indels only
    reads.clear();

    let indel_read_bases = &ref_bases[pos_start..pos_start + 7];
    reads.push(f.input(pos, indel_read_bases, "1M2I1M3D3M", false));
    reads.push(f.input(pos, indel_read_bases, "1M2I1M3D3M", false));

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(IndelMatch, read_info.outcome);
    assert_eq!(5, num_mutations(&read_info.consensus_read));

    // mismatches only
    reads.clear();

    let mut mutated_bases = set_char_at(read_bases, 0, next_base(char_at(read_bases, 0)));
    mutated_bases = set_char_at(&mutated_bases, 1, next_base(char_at(read_bases, 2)));
    mutated_bases = set_char_at(&mutated_bases, 4, next_base(char_at(read_bases, 4)));
    reads.push(f.input(pos, &mutated_bases, "10M", false));
    reads.push(f.input(pos, &mutated_bases, "10M", false));

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(3, num_mutations(&read_info.consensus_read));

    // indels and mismatches
    reads.clear();

    let mut indel_mutated_bases = set_char_at(indel_read_bases, 0, next_base(char_at(read_bases, 0)));
    indel_mutated_bases = set_char_at(&indel_mutated_bases, 1, next_base(char_at(read_bases, 2)));
    indel_mutated_bases = set_char_at(&indel_mutated_bases, 4, next_base(char_at(read_bases, 4)));
    reads.push(f.input(pos, &indel_mutated_bases, "2M1I1M3D3M", false));
    reads.push(f.input(pos, &indel_mutated_bases, "2M1I1M3D3M", false));

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(IndelMatch, read_info.outcome);
    assert_eq!(7, num_mutations(&read_info.consensus_read));

    // left soft clip
    reads.clear();

    let soft_clip_bases = set_char_at(read_bases, 0, next_base(char_at(read_bases, 0)));
    reads.push(f.input(pos + 1, &soft_clip_bases, "1S9M", false));
    reads.push(f.input(pos + 1, &soft_clip_bases, "1S9M", false));

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(0, num_mutations(&read_info.consensus_read));

    // left hard clip
    reads.clear();

    let hard_clip_bases = delete_char_at(read_bases, 0);
    reads.push(f.input(pos + 1, &hard_clip_bases, "1H9M", false));
    reads.push(f.input(pos + 1, &hard_clip_bases, "1H9M", false));

    read_info = create_consensus_read(&mut f.consensus_reads, &reads, UMI_ID_1);
    assert_eq!(AlignmentOnly, read_info.outcome);
    assert_eq!(0, num_mutations(&read_info.consensus_read));
}

/// ConsensusReadsTest.testUnmappedMates
#[test]
fn test_unmapped_mates() {
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();
    let consensus_bases = &ref_bases[0..50];

    let read_cigar1 = "50M";
    let id1 = f.next_read_id();
    let mut read1 =
        create_sam_record(&id1, CHR_1, 1, consensus_bases, read_cigar1, NO_CHROMOSOME_NAME, 1, false, false, None, false, Some(NO_CIGAR));
    rec::set_flag(&mut read1, FLAG_MATE_UNMAPPED, true);
    set_mate_alignment_start(&mut read1, 1);

    let read_cigar2 = "1S49M";
    let id2 = f.next_read_id();
    let mut read2 =
        create_sam_record(&id2, CHR_1, 2, consensus_bases, read_cigar2, NO_CHROMOSOME_NAME, 2, false, false, None, false, Some(NO_CIGAR));
    rec::set_flag(&mut read2, FLAG_MATE_UNMAPPED, true);
    set_mate_alignment_start(&mut read2, 2);

    let coords = create_fragment_coords(&read1);
    let read_consensus_info =
        create_consensus_read_with_coords(&mut f.consensus_reads, &[read1.clone(), read2.clone()], &coords, Some(""));

    assert_eq!(read_consensus_info.template_read, read_name(&read1));
    assert_eq!(read_cigar1, cigar_string(&read_consensus_info.consensus_read));
    assert_eq!(1, rec::start(&read_consensus_info.consensus_read));
    assert_eq!(1, rec::mate_start(&read_consensus_info.consensus_read));
    assert!(rec::has(&read_consensus_info.consensus_read, FLAG_MATE_UNMAPPED));

    let mut mate1 = create_sam_record(
        &read_name(&read1), NO_CHROMOSOME_NAME, 1, consensus_bases, NO_CIGAR, CHR_1, 1, false, false, None, false, Some(read_cigar1),
    );
    rec::set_flag(&mut mate1, FLAG_UNMAPPED, true);
    set_alignment_start(&mut mate1, 1);

    let mut mate2 = create_sam_record(
        &read_name(&read2), NO_CHROMOSOME_NAME, 2, consensus_bases, NO_CIGAR, CHR_1, 2, false, false, None, false, Some(read_cigar2),
    );
    rec::set_flag(&mut mate2, FLAG_UNMAPPED, true);
    set_alignment_start(&mut mate2, 2);

    let coords = create_fragment_coords(&mate1);
    let mate_consensus_info = create_consensus_read_with_coords(&mut f.consensus_reads, &[mate1, mate2], &coords, Some(""));

    assert_eq!(NO_CIGAR, cigar_string(&mate_consensus_info.consensus_read));
    assert_eq!(1, rec::start(&mate_consensus_info.consensus_read));
    assert_eq!(1, rec::mate_start(&mate_consensus_info.consensus_read));
    assert!(rec::has(&mate_consensus_info.consensus_read, FLAG_UNMAPPED));
    assert!(!rec::has(&mate_consensus_info.consensus_read, FLAG_MATE_UNMAPPED));
}

/// ConsensusReadsTest.testDifferentUnclippedPositionStarts
#[test]
fn test_different_unclipped_position_starts() {
    // the template read has no soft-clips and other more outer bases are ignored
    let mut f = Fixture::new();
    let ref_bases = f.ref_bases.clone();

    // pos         0123456789
    // read 1/2:   GTACGTACGT
    // read 3:  AAAGTACGTACGT T
    // read 4:     GTACGTACGT ACT

    let mut read_bases1 = ref_bases[20..30].to_string();
    let mut read_cigar1 = "10M";
    let mut read1 = f.input(10, &read_bases1, read_cigar1, false);

    let mut read2 = clone_sam_record(&read1, &f.next_read_id());

    // as in the Java test, the bases are two short of this cigar
    let mut read_bases3 = format!("AAA{}T", &ref_bases[22..30]); // unclipped pos start = 9
    let mut read_cigar3 = "3S10M1S";
    let mut read3 = f.input(12, &read_bases3, read_cigar3, false);

    let mut read_bases4 = format!("{}TT", &ref_bases[20..32]); // unclipped pos start = 8
    let mut read_cigar4 = "2S10M2S";
    let mut read4 = f.input(10, &read_bases4, read_cigar4, false);

    let coords = create_fragment_coords(&read1);
    let mut read_consensus_info = create_consensus_read_with_coords(
        &mut f.consensus_reads, &[read1.clone(), read2.clone(), read3.clone(), read4.clone()], &coords, Some(""),
    );

    assert_eq!(read_consensus_info.template_read, read_name(&read1));
    assert_eq!(read_cigar1, cigar_string(&read_consensus_info.consensus_read));
    assert_eq!(10, rec::start(&read_consensus_info.consensus_read));
    assert_eq!(read_bases1, read_string(&read_consensus_info.consensus_read));

    // reads without inner bases are skipped
    // read 1/2:  AAAA GTACGTACGT TTTT
    // read 3:       G GTACGTACGT A
    // read 4:          TACGTACGT AT
    // read 5:       G GTACGTACG
    // read 6:      GG GTACGTACGT CC
    // net:       AAAG GTACGTACGT ATTT

    read_bases1 = format!("AAAA{}TTTT", &ref_bases[20..30]);
    read_cigar1 = "4S10M4S";
    read1 = f.input(10, &read_bases1, read_cigar1, false);

    read2 = clone_sam_record(&read1, &f.next_read_id());

    read_bases3 = format!("G{}A", &ref_bases[20..30]); // unclipped pos start = 9
    read_cigar3 = "3S8M1S";
    read3 = f.input(12, &read_bases3, read_cigar3, false);

    read_bases4 = format!("{}AT", &ref_bases[21..30]); // unclipped pos start = 8
    read_cigar4 = "9M2S";
    read4 = f.input(11, &read_bases4, read_cigar4, false);

    let read_bases5 = format!("G{}", &ref_bases[20..29]);
    let read_cigar5 = "10M";
    let read5 = f.input(9, &read_bases5, read_cigar5, false);

    let read_bases6 = format!("GG{}CC", &ref_bases[20..30]);
    let read_cigar6 = "12M2S";
    let read6 = f.input(8, &read_bases6, read_cigar6, false);

    let coords = create_fragment_coords(&read1);
    read_consensus_info = create_consensus_read_with_coords(
        &mut f.consensus_reads, &[read1.clone(), read2, read3, read4, read5, read6], &coords, Some(""),
    );

    assert_eq!(read_consensus_info.template_read, read_name(&read1));
    assert_eq!(read_cigar1, cigar_string(&read_consensus_info.consensus_read));
    assert_eq!(10, rec::start(&read_consensus_info.consensus_read));

    let net_bases = format!("AAAG{}ATTT", &ref_bases[20..30]);
    assert_eq!(net_bases, read_string(&read_consensus_info.consensus_read));
}
