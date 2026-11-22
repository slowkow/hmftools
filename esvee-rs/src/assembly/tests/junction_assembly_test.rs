//! JunctionAssemblyTest
use super::*;
use crate::assembly::assembly_deduper::dedup_proximate_assemblies;
use crate::assembly::assembly_utils::{mismatches_per_comparison_length, DNA_BASE_BYTES};
use crate::assembly::constants::{ASSEMBLY_SPLIT_MIN_READ_SUPPORT, LOW_BASE_QUAL_THRESHOLD};
use crate::assembly::extension_seq_builder::{ExtensionSeqBuilder, ReadJunctionInfo};
use crate::assembly::indel_builder::calc_indel_inferred_unclipped_positions;
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::junction_read_types::JunctionReadTypes;
use crate::assembly::sequence_compare::matched_assembly_sequences;
use crate::assembly::sequence_diff::SequenceDiffType;
use crate::assembly::types::repeat_info::{get_repeat_count_ext, RepeatInfo};

const MISMATCHES_AND_DELETIONS_ATTRIBUTE: &str = "MD";

fn next_id() -> String {
    READ_ID_GENERATOR.next_id()
}
fn read_at(start: i32, bases: &str, cigar: &str) -> ReadP {
    create_read(&next_id(), start, bases, cigar)
}
fn set_base(read: &ReadP, index: usize, base: u8) {
    read.m().bam_record_mut().read_bases[index] = base;
}
fn set_qual(read: &ReadP, index: usize, qual: u8) {
    read.m().bam_record_mut().base_qualities[index] = qual;
}
fn base(read: &ReadP, index: usize) -> u8 {
    read.get_bases()[index]
}
fn process_junction(junction: P<Junction>, reads: &[ReadP]) -> Vec<AsmP> {
    JunctionAssembler::new(junction, None, None).process_junction(reads)
}

#[test]
fn test_misc_assembly_functions() {
    let _env = test_env();
    for (length, expected) in [(1, 0), (15, 1), (100, 1), (101, 2), (201, 3), (399, 3), (500, 4)] {
        assert_eq!(expected, mismatches_per_comparison_length(length));
    }
}

#[test]
fn test_basic_junction_assemblies() {
    let _env = test_env();
    let ext_bases = &REF_BASES_400[300..350];
    let junction = new_junction(CHR_1, 100, FORWARD);
    // assembly has 2 high-qual junction reads
    let read_bases = format!("{}{}", &REF_BASES_400[51..101], ext_bases);
    let read1 = read_at(51, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    let read_bases2 = &read_bases[1..];
    let read2 = read_at(52, &read_bases, &make_cigar_string(read_bases2, 0, ext_bases.len()));
    // the 3rd read doesn't overlap the junction enough to count as a junction read
    // (as in Java its cigar comes from a sequence shorter than the soft-clip, giving a 42M40S cigar for 82 bases)
    let read_bases3 = format!("{}{}", &REF_BASES_400[21..101], &ext_bases[0..2]);
    let read3 = read_at(21, &read_bases3, &format!("{}M{}S", read_bases3.len() as i32 - ext_bases.len() as i32, ext_bases.len()));
    let mut junction_assembler = JunctionAssembler::new(junction, None, None);
    let assemblies = junction_assembler.process_junction(&[read1.clone(), read2, read3]);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    assert_eq!(50, assembly.ref_base_length());
    assert_eq!(51, assembly.ref_base_position());
    assert_eq!(2, assembly.support_count()); // read counts as support, just doesn't extend the ref bases

    // again where the second read is not an initial extension read - has the same bases but a different SC position
    let read_bases = read1.get_bases_string()[1..].to_string();
    let read2 = read_at(52, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len() - 5));
    let assemblies = junction_assembler.process_junction(&[read1, read2]);
    assert_eq!(1, assemblies.len());
}

fn contained_candidate_junction_read(read: &ReadP, candidate_reads: &[ReadJunctionInfo]) -> bool {
    candidate_reads.iter().any(|x| x.read == *read)
}

fn check_junction_read_types(junction: P<Junction>, reads: &[ReadP], ref_genome: &MockRefGenome) {
    let junction_read_types = JunctionReadTypes::new(junction, reads, Some(ref_genome));
    assert_eq!(3, junction_read_types.non_junction_reads().len());
    let candidate_reads = junction_read_types.candidate_junction_reads();
    assert_eq!(1, candidate_reads.iter().filter(|x| x.matches_junction).count());
    assert_eq!(5, candidate_reads.len());
    for read in &reads[0..5] {
        assert!(contained_candidate_junction_read(read, candidate_reads));
    }
}

#[test]
fn test_junction_read_types() {
    let _env = test_env();
    let ref_bases = &REF_BASES_200[10..30];
    let ext_bases = &REF_BASES_200[100..140];
    let junction = new_junction(CHR_1, 29, FORWARD);
    let mut ref_genome = MockRefGenome::new_indexing(false);
    ref_genome.put(CHR_1, REF_BASES_200);
    let read_bases = format!("{ref_bases}{ext_bases}");
    // a read soft-clipping at the junction
    let read1 = read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    // a read soft-clipping before the junction
    let read2 = read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len() + 5));
    // a read soft-clipping after the junction
    let read3 = read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len() - 5));
    // a read without soft-clipping but evidence of SNVs
    let read4 = read_at(10, &read_bases[0..30], "30M");
    read4.m().bam_record_mut().set_integer_attribute(NUM_MUTATONS_ATTRIBUTE, 2);
    read4.m().bam_record_mut().set_string_attribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, "22A1C1A3");
    // an indel read inferring a soft-clip
    let read5 = read_at(10, &read_bases, "40M10I20M");
    calc_indel_inferred_unclipped_positions(read5.m());
    // a read soft-clipping too far from the junction
    let read6a = read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len() - 21));
    let longer_read_bases = format!("{}{}", &REF_BASES_200[10..55], &ext_bases[0..25]);
    let read6b = read_at(10, &longer_read_bases, &make_cigar_string(&longer_read_bases, 0, 25));
    // a read without soft-clipping or SNVs
    let read_ref_bases = &REF_BASES_200[10..35];
    let read7 = read_at(10, read_ref_bases, &make_cigar_string(read_ref_bases, 0, 0));
    read7.m().bam_record_mut().set_integer_attribute(NUM_MUTATONS_ATTRIBUTE, 2); // assume has SNVs elsewhere in the read
    check_junction_read_types(junction, &[read1, read2, read3, read4, read5, read6a, read6b, read7], &ref_genome);

    // repeat for neg orientation
    let junction = new_junction(CHR_1, 40, REVERSE);
    let read_bases = format!("{ext_bases}{ref_bases}");
    let read1 = read_at(40, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    let read2 = read_at(35, &read_bases, &make_cigar_string(&read_bases, ext_bases.len() + 5, 0));
    let read3 = read_at(45, &read_bases, &make_cigar_string(&read_bases, ext_bases.len() - 5, 0));
    let read4 = read_at(30, &read_bases[0..30], "30M");
    read4.m().bam_record_mut().set_integer_attribute(NUM_MUTATONS_ATTRIBUTE, 2);
    read4.m().bam_record_mut().set_string_attribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, "3A1C1A22");
    let read5 = read_at(20, &read_bases, "20M10I40M");
    calc_indel_inferred_unclipped_positions(read5.m());
    let longer_read_bases = format!("{}{}", &ext_bases[0..25], &REF_BASES_200[10..55]);
    let read6a = read_at(61, &longer_read_bases, &make_cigar_string(&longer_read_bases, 25, 0));
    let read6b = read_at(19, &longer_read_bases, &make_cigar_string(&longer_read_bases, 15, 0));
    let read7 = read_at(30, &read_bases, &make_cigar_string(&read_bases, 0, 0));
    check_junction_read_types(junction, &[read1, read2, read3, read4, read5, read6a, read6b, read7], &ref_genome);
}

/// makeCigarString where the soft-clip may exceed the sequence, as some Java tests do: the M length is just the difference
fn cigar_sc_right(read_bases: &str, sc_right: usize) -> String {
    format!("{}M{}S", read_bases.len() as i32 - sc_right as i32, sc_right)
}
fn cigar_sc_left(read_bases: &str, sc_left: usize) -> String {
    format!("{}S{}M", sc_left, read_bases.len() as i32 - sc_left as i32)
}

#[test]
fn test_pos_junction_extension_sequence() {
    let _env = test_env();
    let ref_bases = &REF_BASES_200[10..30];
    let ext_bases = &REF_BASES_200[100..140];
    let junction = new_junction(CHR_1, 29, FORWARD);
    // first a basic assembly with all reads agreeing
    let read_bases = format!("{ref_bases}{ext_bases}");
    let read1 = read_at(10, &read_bases, &cigar_sc_right(&read_bases, ext_bases.len()));
    let read_bases = format!("{}{}", ref_bases, &ext_bases[0..35]);
    let read2 = read_at(10, &read_bases, &cigar_sc_right(&read_bases, ext_bases.len()));
    let read_bases = format!("{}{}", ref_bases, &ext_bases[0..32]);
    let read3 = read_at(10, &read_bases, &cigar_sc_right(&read_bases, ext_bases.len()));
    // a read without soft-clip but mismatches agreeing with the extension
    let aligned_matching_bases = &read_bases[0..22];
    let read4 = read_at(10, aligned_matching_bases, &make_cigar_string(aligned_matching_bases, 0, 0));
    read4.m().bam_record_mut().set_integer_attribute(NUM_MUTATONS_ATTRIBUTE, 2);
    let reads = vec![read1.clone(), read2.clone(), read3.clone(), read4];
    let mut ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    let sequence = format!("{}{}", &ref_bases[19..], ext_bases);
    assert_eq!(sequence, ext_seq_builder.junction_sequence());
    assert_eq!(4, ext_seq_builder.form_assembly_support().len());

    // adding non-extension reads
    let read_bases = read1.get_bases_string();
    let junc_read = read_at(10, &read_bases, &cigar_sc_right(&read_bases, ext_bases.len()));
    let read_info = ext_seq_builder.check_add_junction_read(&junc_read).expect("read state");
    assert!(!read_info.mismatched());

    // low qual mismatches in all 3 three reads but still overall agreement
    set_base(&read2, 30, get_next_base_byte(base(&read1, 30)));
    set_base(&read2, 31, get_next_base_byte(base(&read1, 31)));
    set_qual(&read2, 30, 20);
    set_qual(&read2, 31, 20);
    set_base(&read3, 32, get_next_base_byte(base(&read1, 32)));
    set_base(&read3, 33, get_next_base_byte(base(&read1, 33)));
    set_qual(&read3, 32, 20);
    set_qual(&read3, 33, 20);
    let read_bases = format!("{}{}TTTT{}", ref_bases, &ext_bases[0..29], &ext_bases[33..38]);
    let read4 = read_at(10, &read_bases, &cigar_sc_right(&read_bases, ext_bases.len()));
    let reads = vec![read1.clone(), read2, read3, read4];
    let mut ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    assert_eq!(sequence, ext_seq_builder.junction_sequence());
    assert_eq!(3, ext_seq_builder.form_assembly_support().len());

    // reads with novel INDELs
    let new_ext_bases = format!("{}A{}", &ext_bases[0..10], &ext_bases[10..20]);
    let read_bases = format!("{ref_bases}{new_ext_bases}");
    let junc_read = read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, new_ext_bases.len()));
    let read_info = ext_seq_builder.check_add_junction_read(&junc_read).expect("read state");
    assert!(!read_info.mismatched());
    assert_eq!(SequenceDiffType::Insert, read_info.mismatches()[0].type_);
    let new_ext_bases = format!("{}{}", &ext_bases[0..10], &ext_bases[11..20]);
    let read_bases = format!("{ref_bases}{new_ext_bases}");
    let junc_read = read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, new_ext_bases.len()));
    let read_info = ext_seq_builder.check_add_junction_read(&junc_read).expect("read state");
    assert!(!read_info.mismatched());
    assert_eq!(SequenceDiffType::Delete, read_info.mismatches()[0].type_);

    // a read with most low qual mismatches
    let junc_read = clone_read(&read1, &next_id());
    let read_index_start = junc_read.bases_length() as usize - ext_bases.len() + 5;
    for i in 0..15 {
        let read_index = read_index_start + i;
        set_base(&junc_read, read_index, get_next_base_byte(base(&junc_read, read_index)));
        set_qual(&junc_read, read_index, LOW_BASE_QUAL_THRESHOLD);
    }
    let read_info = ext_seq_builder.check_add_junction_read(&junc_read).expect("read state");
    assert!(!read_info.mismatched());
    assert!(!ext_seq_builder.sufficient_qual_matches(&read_info));
    assert_eq!(17, read_info.mismatch_count(false));
}

#[test]
fn test_neg_junction_extension_sequence() {
    let _env = test_env();
    let ref_bases = &REF_BASES_200[0..20];
    let ext_bases = &REF_BASES_200[100..140];
    let junc_position = 100;
    let junction = new_junction(CHR_1, junc_position, REVERSE);
    // first a basic assembly with all reads agreeing
    let read_bases = format!("{ext_bases}{ref_bases}");
    let read1 = read_at(junc_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    let read_bases = format!("{}{}", &ext_bases[5..], ref_bases);
    let read2 = read_at(junc_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len() - 5, 0));
    let read_bases = format!("{}{}", &ext_bases[8..], ref_bases);
    let read3 = read_at(junc_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len() - 8, 0));
    let aligned_matching_bases = &read_bases[read_bases.len() - 22..];
    let read4 = read_at(junc_position - 2, aligned_matching_bases, &make_cigar_string(aligned_matching_bases, 0, 0));
    let reads = vec![read1.clone(), read2.clone(), read3.clone(), read4];
    let mut ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    let sequence = format!("{}{}", ext_bases, &ref_bases[0..1]);
    assert_eq!(sequence, ext_seq_builder.junction_sequence());
    assert_eq!(4, ext_seq_builder.form_assembly_support().len());

    // adding non-extension reads
    let read_bases = read1.get_bases_string();
    let junc_read = read_at(junc_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    let read_parse_state = ext_seq_builder.check_add_junction_read(&junc_read).expect("read state");
    assert!(!read_parse_state.mismatched());

    // low qual mismatches in all 3 three reads but still overall agreement
    set_base(&read2, 10, get_next_base_byte(base(&read1, 30)));
    set_base(&read2, 11, get_next_base_byte(base(&read1, 31)));
    set_qual(&read2, 10, 20);
    set_qual(&read2, 11, 20);
    set_base(&read3, 12, get_next_base_byte(base(&read1, 32)));
    set_base(&read3, 13, get_next_base_byte(base(&read1, 33)));
    set_qual(&read3, 12, 20);
    set_qual(&read3, 13, 20);
    let read_bases = format!("{}TTTT{}{}", &ext_bases[0..10], &ext_bases[14..], ref_bases);
    let read4 = read_at(junc_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    let reads = vec![read1, read2, read3, read4];
    let ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    assert_eq!(sequence, ext_seq_builder.junction_sequence());
    assert_eq!(3, ext_seq_builder.form_assembly_support().len());
}

#[test]
fn test_high_repeat_extension_sequence() {
    let _env = test_env();
    let ref_bases = &REF_BASES_200[0..20];
    let junction = new_junction(CHR_1, 29, FORWARD);
    let buffer_bases = "CAGTAGGAGCTTTA";
    let consensus_repeat = "AAAAAAAA";
    let consensus_ext_bases = format!("{buffer_bases}{consensus_repeat}{buffer_bases}");
    let read_with = |ext_bases: &str| -> ReadP {
        let read_bases = format!("{ref_bases}{ext_bases}");
        read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()))
    };
    let read1 = read_with(&consensus_ext_bases);
    // reads agreeing to establish a base-line for repeats
    let read1b = clone_read(&read1, &next_id());
    let read1c = clone_read(&read1, &next_id());
    let read2 = read_with(&format!("{buffer_bases}{consensus_repeat}A{buffer_bases}")); // extra A
    let read3 = read_with(&format!("{}{}{}", buffer_bases, &consensus_repeat[1..], buffer_bases)); // 1 less A
    let read4 = read_with(&format!("{}{}{}", buffer_bases, &consensus_repeat[3..], buffer_bases)); // 3 less As, too many and will mismatch
    let read5 = read_with(&format!("{buffer_bases}{consensus_repeat}AA{buffer_bases}")); // 2 extra As
    let read6 = read_with(&format!("{buffer_bases}{consensus_repeat}AAA{buffer_bases}")); // 3 extra As
    let reads = vec![read1, read1b, read1c, read2, read3, read4, read5, read6];
    let ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    assert_eq!(format!("{}{}", &ref_bases[19..], consensus_ext_bases), ext_seq_builder.junction_sequence());
    let support_reads = ext_seq_builder.form_assembly_support();
    assert_eq!(6, support_reads.len());
    assert_eq!(5, support_reads.iter().filter(|x| x.extension_base_mismatches() == 0).count());
}

#[test]
fn test_long_repeat_with_low_quals() {
    let _env = test_env();
    // first test repeat bases with low quals of differing counts
    let test_bases = format!("CAG{}", "T".repeat(20));
    let mut test_bytes = test_bases.into_bytes();
    let mut test_quals = build_base_quals(test_bytes.len(), 37);
    let non_repeat_base = DNA_BASE_BYTES[0];
    test_bytes[6] = non_repeat_base;
    test_quals[6] = 11;
    assert_eq!(20, get_repeat_count_ext(&test_bytes, "T", 3, true, 0, Some(&test_quals)));
    // now with too many mismatches so reduces the observed repeat count
    test_bytes[12] = non_repeat_base;
    test_quals[12] = 11;
    test_bytes[16] = non_repeat_base;
    test_quals[16] = 11;
    assert_eq!(3, get_repeat_count_ext(&test_bytes, "T", 3, true, 0, Some(&test_quals)));

    let ref_bases = &REF_BASES_200[0..20];
    let junction = new_junction(CHR_1, 29, FORWARD);
    let buffer_bases = "CAG";
    let consensus_repeat = "T".repeat(20);
    let consensus_ext_bases = format!("{buffer_bases}{consensus_repeat}{buffer_bases}");
    let read_with = |ext_bases: &str| -> ReadP {
        let read_bases = format!("{ref_bases}{ext_bases}");
        read_at(10, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()))
    };
    let read1 = read_with(&consensus_ext_bases);
    let read1b = clone_read(&read1, &next_id());
    let read1c = clone_read(&read1, &next_id());
    let read2 = read_with(&format!("{}{}{}", buffer_bases, &consensus_repeat[2..], buffer_bases)); // 2 less
    // a read with a repeat matching but for an acceptable number of low-qual mismatches
    let read3 = read_with(&format!("{buffer_bases}{consensus_repeat}A{buffer_bases}")); // extra A
    set_base(&read3, 27, non_repeat_base);
    set_base(&read3, 40, non_repeat_base);
    set_qual(&read3, 27, 37);
    set_qual(&read3, 40, 1);
    // a read with too many low-qual mismatches
    let read4 = read_with(&format!("{buffer_bases}{consensus_repeat}{buffer_bases}"));
    for index in [30, 35, 40] {
        set_base(&read4, index, non_repeat_base);
        set_qual(&read4, index, 1);
    }
    let read5 = read_with(&format!("{}{}{}", buffer_bases, &consensus_repeat[3..], buffer_bases)); // 3 less, 1 low qual is OK
    set_qual(&read5, 30, 1);
    set_base(&read5, 30, non_repeat_base);
    let reads = vec![read1, read1b, read1c, read2, read3, read4, read5];
    let ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    assert_eq!(format!("{}{}", &ref_bases[19..], consensus_ext_bases), ext_seq_builder.junction_sequence());
    assert_eq!(6, ext_seq_builder.form_assembly_support().len());
}

#[test]
fn test_repeat_mismatches_extension_sequence() {
    let _env = test_env();
    let ref_bases = &REF_BASES_200[0..20];
    let junction_position = 50;
    let junction = new_junction(CHR_1, 50, REVERSE);
    let buffer = "GATCGTAGGATC";
    let ca_repeat = "CACACACA";
    let read_with = |ext_bases: &str| -> ReadP {
        let read_bases = format!("{ext_bases}{ref_bases}");
        read_at(junction_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0))
    };
    let read1 = read_with(&format!("{buffer}{ca_repeat}{buffer}"));
    let read2 = read_with(&format!("{buffer}{ca_repeat}CA{buffer}")); // 1 extra
    let read3 = read_with(&format!("{}{}{}", buffer, &ca_repeat[2..], buffer)); // one less CA
    let read4 = read_with(&format!("{buffer}{ca_repeat}CACA{buffer}")); // 2 extra CAs, too many for this number of repeats
    let read1b = clone_read(&read1, &next_id());
    let reads = vec![read1.clone(), read2.clone(), read3.clone(), read4.clone(), read1b];
    let mut ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    let consensus_ext_bases = format!("{buffer}{ca_repeat}{buffer}"); // matches read 1
    assert_eq!(format!("{}{}", consensus_ext_bases, &ref_bases[0..1]), ext_seq_builder.junction_sequence());
    let consensus_repeat = ext_seq_builder.repeats()[0].clone();
    assert_eq!(4, consensus_repeat.count);
    assert_eq!(12, consensus_repeat.index_start()); // moved earlier after trimming
    assert_eq!(19, consensus_repeat.index_end());
    let support_reads = ext_seq_builder.form_assembly_support();
    assert_eq!(5, support_reads.len());
    assert_eq!(4, support_reads.iter().filter(|x| x.extension_base_mismatches() == 0).count());

    // adding further junction reads with differing repeats
    let read1c = clone_read(&read1, &next_id());
    let read2b = clone_read(&read2, &next_id());
    let read3b = clone_read(&read3, &next_id());
    let read4b = clone_read(&read4, &next_id());
    assert!(ext_seq_builder.check_add_junction_read(&read1c).is_some());
    let read_info = ext_seq_builder.check_add_junction_read(&read2b).expect("read state");
    assert_eq!(1, read_info.mismatches().len());
    assert_eq!(SequenceDiffType::Repeat, read_info.mismatches()[0].type_);
    assert_eq!(5, read_info.mismatches()[0].repeat_count);
    let read_info = ext_seq_builder.check_add_junction_read(&read3b).expect("read state");
    assert_eq!(1, read_info.mismatches().len());
    assert_eq!(SequenceDiffType::Repeat, read_info.mismatches()[0].type_);
    assert_eq!(3, read_info.mismatches()[0].repeat_count);
    let read_info = ext_seq_builder.check_add_junction_read(&read4b).expect("read state");
    assert_eq!(2, read_info.mismatches().len());
    assert_eq!(SequenceDiffType::Base, read_info.mismatches()[0].type_); // was too many to be considered a repeat diff
}

fn max_repeat(extension_seq_builder: &ExtensionSeqBuilder) -> Option<RepeatInfo> {
    // the first of the repeats with the highest count (a stable sort by descending count in Java)
    let mut best: Option<&RepeatInfo> = None;
    for repeat in extension_seq_builder.repeats() {
        if best.is_none_or(|b| repeat.count > b.count) {
            best = Some(repeat);
        }
    }
    best.cloned()
}

#[test]
fn tes_ref_base_repeat_expansion_mismatches_extension_sequence() {
    let _env = test_env();
    let junction_position = 50;
    let junction = new_junction(CHR_1, junction_position, REVERSE);
    let ref_bases = format!("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG{}", &REF_BASES_200[0..20]); // has 10 repeats
    let extra_ext_bases = "GATCGTAGGATCGTAGGATCGTAGGATC"; // other non-repeated bases
    let read_with = |ext_bases: &str| -> ReadP {
        let read_bases = format!("{ext_bases}{ref_bases}");
        read_at(junction_position, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0))
    };
    let ext_bases1 = format!("{}{}", extra_ext_bases, "CAG".repeat(4)); // 4 repeats
    let read1 = read_with(&ext_bases1);
    let read1b = clone_read(&read1, &next_id());
    let read3 = read_with(&format!("{}{}", extra_ext_bases, "CAG".repeat(6))); // 6 repeats, 2 more than consensus
    let read4 = read_with(&format!("{}{}", extra_ext_bases, "CAG".repeat(3))); // 1 less CAG
    let read5 = read_with("CAG"); // insufficient overlap bases
    let reads = vec![read1, read1b, read3, read4, read5];
    let ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    let repeat = max_repeat(&ext_seq_builder).expect("repeat");
    assert_eq!(4, repeat.count);
    assert_eq!(10, ext_seq_builder.ref_base_repeat_count());
    assert_eq!(format!("{}{}", ext_bases1, &ref_bases[0..1]), ext_seq_builder.junction_sequence());
    let support_reads = ext_seq_builder.form_assembly_support();
    assert_eq!(4, support_reads.len());
    assert_eq!(3, support_reads.iter().filter(|x| x.extension_base_mismatches() == 0).count());
    assert_eq!("RC=5;EM=2;LQ=1;OR=1;REP:4xCAG", ext_seq_builder.build_information());

    // a forward orientation junction
    let junction = new_junction(CHR_1, junction_position, FORWARD);
    let ref_bases = format!("{}CACACACACACACACACACA", &REF_BASES_200[0..20]); // has 10 repeats
    let read_start_pos = junction_position - ref_bases.len() as i32 + 1;
    let read_with = |ext_bases: &str| -> ReadP {
        let read_bases = format!("{ref_bases}{ext_bases}");
        read_at(read_start_pos, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()))
    };
    let ext_bases1 = format!("{}{}", "CA".repeat(4), extra_ext_bases); // has 4 repeats
    let read1 = read_with(&ext_bases1);
    let read1b = clone_read(&read1, &next_id());
    let read3 = read_with(&format!("{}{}", "CA".repeat(2), ext_bases1)); // 2 extra repeats
    let read4 = read_with(&format!("{}{}", "CA".repeat(3), extra_ext_bases)); // 1 less CA
    let read5 = read_with("CA"); // insufficient overlap bases
    let reads = vec![read1, read1b, read3, read4, read5];
    let ext_seq_builder = ExtensionSeqBuilder::new(junction, &reads);
    let repeat = max_repeat(&ext_seq_builder).expect("repeat");
    assert_eq!(4, repeat.count);
    assert_eq!(10, ext_seq_builder.ref_base_repeat_count());
    assert_eq!(format!("{}{}", &ref_bases[ref_bases.len() - 1..], ext_bases1), ext_seq_builder.junction_sequence());
    let support_reads = ext_seq_builder.form_assembly_support();
    assert_eq!(4, support_reads.len());
    assert_eq!(3, support_reads.iter().filter(|x| x.extension_base_mismatches() == 0).count());
}

fn next_assembly_base(assembly: AsmP, target_index: usize, source: AsmP, source_index: usize) {
    let new_base = get_next_base_byte(source.bases()[source_index]);
    assembly.m().bases_mut()[target_index] = new_base;
}

#[test]
fn test_exact_assembly_dedup() {
    let _env = test_env();
    let assembly_bases = &REF_BASES_200[0..50];
    let pos_junction = new_junction(CHR_1, 60, FORWARD);
    let base_quals = build_default_base_quals(assembly_bases.len());
    let assemblies: Vec<AsmP> = (0..5).map(|_| JunctionAssembly::new_test(pos_junction, assembly_bases.as_bytes(), &base_quals, 21)).collect();
    let (assembly1, assembly2, assembly3, assembly4, assembly5) = (assemblies[0], assemblies[1], assemblies[2], assemblies[3], assemblies[4]);
    let shared_read = create_read(TEST_READ_ID, 40, &assembly_bases[0..30], "30M");
    let shared_read2 = create_read(TEST_READ_ID_2, 40, &assembly_bases[0..30], "30M");
    for assembly in &assemblies {
        assembly.m().add_junction_read(&shared_read);
    }
    for assembly in &assemblies {
        assembly.m().add_junction_read(&shared_read2);
    }
    // add 2 mismatches to each of the other assemblies
    next_assembly_base(assembly2, 2, assembly2, 2);
    next_assembly_base(assembly2, 30, assembly2, 30);
    next_assembly_base(assembly3, 10, assembly2, 10);
    next_assembly_base(assembly3, 20, assembly2, 20);
    // add 5 to the next one
    for index in 20..25 {
        next_assembly_base(assembly4, index, assembly4, index);
    }
    assert!(matched_assembly_sequences(&assembly1, &assembly2));
    assert!(matched_assembly_sequences(&assembly1, &assembly3));
    assert!(!matched_assembly_sequences(&assembly1, &assembly4));
    assert!(matched_assembly_sequences(&assembly1, &assembly5));
}

#[test]
fn test_assembly_no_dedup_with_shared_reads() {
    let _env = test_env();
    let indel_junction_neg = P::new(Junction::new(CHR_1, 130, REVERSE, false, true, false, 0, None));
    let indel_junction_pos = P::new(Junction::new(CHR_1, 131, FORWARD, false, true, false, 0, None));
    let left_ext_bases = &REF_BASES_400[200..240];
    let right_ext_bases = &REF_BASES_400[260..300];
    let left_ref_bases = &REF_BASES_400[100..131];
    let right_ref_bases = &REF_BASES_400[131..171];
    let insert_bases = right_ref_bases; // a duplication
    let indel_read_bases = format!("{left_ext_bases}{left_ref_bases}{insert_bases}{right_ref_bases}{right_ext_bases}");
    // both pairs of junctions share a read which supports them both, so are not deduped
    let indel_read = create_read_mate(&next_id(), CHR_1, 100, &indel_read_bases, "40S31M40I40M40S", CHR_1, 300, true);
    let read_assembly = |junction: P<Junction>, read: &ReadP, junction_index: usize, set_indel: bool| -> AsmP {
        let assembly = JunctionAssembly::new_test(junction, read.get_bases(), read.get_base_quality(), junction_index as i32);
        assembly.m().add_junction_read(read);
        if set_indel {
            assembly.m().set_indel_coords(read.indel_coords());
        }
        assembly
    };
    let indel_assembly_pos = read_assembly(indel_junction_pos, &indel_read, left_ext_bases.len() + left_ref_bases.len(), true);
    let indel_assembly_neg = read_assembly(indel_junction_neg, &indel_read, left_ext_bases.len() + left_ref_bases.len() + 1, true);
    let split_junction_pos = new_junction(CHR_1, 170, FORWARD);
    let split_junction_neg = new_junction(CHR_1, 100, REVERSE);
    let split_assembly_neg = read_assembly(split_junction_neg, &indel_read, left_ext_bases.len(), false);
    let split_assembly_pos = read_assembly(split_junction_pos, &indel_read, indel_read.bases_length() as usize - right_ext_bases.len() - 1, false);
    let mut existing_assemblies = vec![split_assembly_pos, split_assembly_neg];
    let mut new_assemblies = vec![indel_assembly_pos, indel_assembly_neg];
    let mut dedup_indels: Vec<AsmP> = Vec::new();
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(2, existing_assemblies.len());
    assert_eq!(2, new_assemblies.len());

    // again without a shared read and they are deduped
    let split_read = create_read_mate(&next_id(), CHR_1, 100, &indel_read_bases, "40S31M30M70S", CHR_1, 300, true);
    let split_assembly_neg = read_assembly(split_junction_neg, &split_read, left_ext_bases.len(), false);
    let split_assembly_pos = read_assembly(split_junction_pos, &split_read, split_read.bases_length() as usize - right_ext_bases.len() - 1, false);
    let mut existing_assemblies = vec![split_assembly_pos, split_assembly_neg];
    let mut new_assemblies = vec![indel_assembly_pos, indel_assembly_neg];
    dedup_proximate_assemblies(&mut existing_assemblies, &mut new_assemblies, &mut dedup_indels);
    assert_eq!(2, new_assemblies.len());
    assert!(existing_assemblies.is_empty()); // the indel assemblies are kept
}

#[test]
fn test_junction_assembly_splitting() {
    let _env = test_env();
    let ref_sequence = REF_BASES_400;
    let junction = new_junction(CHR_1, 100, FORWARD);
    let junc_read_bases1 = format!("{}{}", &ref_sequence[51..101], &ref_sequence[200..250]);
    let junc_read1 = create_read_mate(&next_id(), CHR_1, 51, &junc_read_bases1, "50M50S", CHR_1, 1, false);
    let mut junction_reads = vec![junc_read1.clone()];
    for _ in 0..9 {
        junction_reads.push(clone_read(&junc_read1, &next_id()));
    }
    junction_reads.push(create_read_mate(&next_id(), CHR_1, 52, &junc_read_bases1[1..], "49M50S", CHR_1, 1, false));
    let junc_read_bases2 = format!("{}{}", &ref_sequence[51..101], &ref_sequence[270..320]);
    let junc_read2 = create_read_mate(&next_id(), CHR_1, 51, &junc_read_bases2, "50M50S", CHR_1, 1, false);
    junction_reads.push(junc_read2.clone());
    for _ in 0..ASSEMBLY_SPLIT_MIN_READ_SUPPORT - 1 {
        junction_reads.push(clone_read(&junc_read2, &next_id()));
    }
    junction_reads.push(create_read_mate(&next_id(), CHR_1, 52, &junc_read_bases2[1..], "49M50S", CHR_1, 1, false));
    let assemblies = process_junction(junction, &junction_reads);
    assert_eq!(2, assemblies.len());
    assert_eq!(11, assemblies[0].support_count());
    assert_eq!(6, assemblies[1].support_count());
}

#[test]
fn test_junction_assembly_ref_side_soft_clip_read_assignment() {
    let _env = test_env();
    // reads which cross the dominant ref-side soft clip are removed from junction support
    let junction = new_junction(CHR_1, 200, FORWARD);
    let ext_bases_shared = "AAGGCCTTAA";
    let ext_bases1 = format!("{}{}", ext_bases_shared, &REF_BASES_400[300..340]);
    let ext_bases2 = format!("{}{}", ext_bases_shared, &REF_BASES_400[350..390]);
    // reads which support the first assembly and which soft-clip to the left
    let other_soft_clip_bases = &REF_BASES_400[0..50];
    let junc_read_bases1 = format!("{}{}{}", other_soft_clip_bases, &REF_BASES_400[151..201], ext_bases1);
    let junc_read1 = create_read_mate(&next_id(), CHR_1, 151, &junc_read_bases1, "50S50M50S", CHR_2, 100, false);
    let junc_read1b = create_read_mate(&next_id(), CHR_1, 152, &junc_read_bases1[1..], "50S49M50S", CHR_2, 100, false);
    let mut junction_reads = vec![junc_read1.clone(), junc_read1b];
    for _ in 0..5 {
        junction_reads.push(clone_read(&junc_read1, &next_id()));
    }
    let junc_read_bases2 = format!("{}{}", &REF_BASES_400[101..201], ext_bases2);
    let junc_read2 = create_read_mate(&next_id(), CHR_1, 101, &junc_read_bases2, "100M50S", CHR_2, 100, false);
    let junc_read2b = create_read_mate(&next_id(), CHR_1, 102, &junc_read_bases2[1..], "99M50S", CHR_2, 100, false);
    junction_reads.push(junc_read2.clone());
    junction_reads.push(junc_read2b);
    for _ in 0..5 {
        junction_reads.push(clone_read(&junc_read2, &next_id()));
    }
    // then a read which could support either extension but which goes past the soft-clip and will be removed from the first assembly
    let junc_read_bases3 = format!("{}{}", &REF_BASES_400[100..201], ext_bases_shared);
    junction_reads.push(create_read_mate(&next_id(), CHR_1, 100, &junc_read_bases3, "101M10S", CHR_2, 100, false));
    let assemblies = process_junction(junction, &junction_reads);
    assert_eq!(2, assemblies.len());
    assert_eq!(7, assemblies[0].support_count());
    assert_eq!(151, assemblies[0].ref_base_position());
    assert_eq!(8, assemblies[1].support_count());
    assert_eq!(100, assemblies[1].ref_base_position());
}
