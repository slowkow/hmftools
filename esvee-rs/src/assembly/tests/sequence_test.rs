//! SequenceTest
use super::*;
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::read_parse_state::{ReadParseState, RpsP};
use crate::assembly::sequence_builder::SequenceBuilder;
use crate::assembly::sequence_compare::compare_sequences;
use crate::assembly::sequence_diff::SequenceDiffType;
use crate::assembly::types::repeat_info::{
    build_trimmed_ref_base_sequence, find_dual_base_repeat, find_dual_dual_repeat, find_multi_base_repeat, find_repeats, find_repeats_min, find_single_base_repeat,
    RepeatInfo,
};

const REF_BUFFER: &str = "CCCCC";
const NO_MISMATCH_LIMIT: i32 = -1;
// hmf-common SbxBamUtils
const SBX_SIMPLEX_QUAL: u8 = 27;
const SBX_DUPLEX_QUAL: u8 = 40;

fn read_100(bases: &str, cigar: &str) -> ReadP {
    create_read(&READ_ID_GENERATOR.next_id(), 100, bases, cigar)
}

fn states(build_forwards: bool, reads: &[&ReadP], start_index: usize) -> Vec<RpsP> {
    reads.iter().map(|r| ReadParseState::new(build_forwards, r, start_index as i32)).collect()
}

fn set_qual(read: &ReadP, index: usize, qual: u8) {
    read.m().bam_record_mut().base_qualities[index] = qual;
}

fn max_read_length(read_states: &[RpsP]) -> i32 {
    read_states.iter().map(|x| x.read().bases_length()).max().unwrap_or(0)
}

#[test]
fn test_seq_build_snv_mismatches() {
    let _env = test_env();
    let read_ext_bases1 = "ACGTAAGTACGT";
    let read_ext_bases3 = "ACGGAAGGACGG";
    let read_bases1 = format!("{REF_BUFFER}{read_ext_bases1}");
    let read_bases3 = format!("{REF_BUFFER}{read_ext_bases3}");
    let cigar = make_cigar_string(&read_bases1, 0, 0);
    let read1 = read_100(&read_bases1, &cigar);
    let read2 = read_100(&read_bases1, &cigar);
    let read3 = read_100(&read_bases3, &cigar);
    let consensus_base_length = (read_bases1.len() - REF_BUFFER.len()) as i32;
    let read_parse_states = states(true, &[&read1, &read2, &read3], REF_BUFFER.len());
    let seq_builder = SequenceBuilder::new_ext(&read_parse_states, true, consensus_base_length, true);
    assert_eq!(read_ext_bases1, seq_builder.base_string());
    assert!(seq_builder.repeats().is_empty());
    assert_eq!(3, read_parse_states[2].mismatches().len());

    // favouring high over medium quals
    set_sbx_sequencing();
    set_qual(&read1, 8, SBX_SIMPLEX_QUAL); // read with medium qual is the tie-breaker
    set_qual(&read2, 8, SBX_DUPLEX_QUAL);
    set_qual(&read3, 8, SBX_DUPLEX_QUAL);
    set_qual(&read1, 12, SBX_SIMPLEX_QUAL);
    set_qual(&read2, 12, SBX_SIMPLEX_QUAL);
    set_qual(&read3, 12, SBX_DUPLEX_QUAL);
    for x in &read_parse_states {
        x.m().reset_all();
    }
    let seq_builder = SequenceBuilder::new(&read_parse_states, true, consensus_base_length);
    assert_eq!("ACGTAAGTACGT", seq_builder.base_string());
    set_illumina_sequencing();

    // building the other direction
    let read_bases1 = format!("{read_ext_bases1}{REF_BUFFER}");
    let read_bases3 = format!("{read_ext_bases3}{REF_BUFFER}");
    let read1 = read_100(&read_bases1, &cigar);
    let read2 = read_100(&read_bases1, &cigar);
    let read3 = read_100(&read_bases3, &cigar);
    let read_parse_states = states(false, &[&read1, &read2, &read3], read_ext_bases1.len() - 1);
    let seq_builder = SequenceBuilder::new(&read_parse_states, false, consensus_base_length);
    assert_eq!(read_ext_bases1, seq_builder.base_string());
    assert!(seq_builder.repeats().is_empty());
}

#[test]
fn test_seq_build_indel_mismatches() {
    let _env = test_env();
    let read_ext_bases1 = "ACGTAAGTACGTGGCC";
    let read_bases1 = format!("{REF_BUFFER}{read_ext_bases1}");
    // with a del then an insert
    let read_ext_bases3 = format!("{}{}A{}", &read_ext_bases1[0..4], &read_ext_bases1[5..10], &read_ext_bases1[10..]);
    let read_bases3 = format!("{REF_BUFFER}{read_ext_bases3}");
    // with 2 inserts then a del
    let read_ext_bases4 = format!("{}A{}A{}{}", &read_ext_bases1[0..3], &read_ext_bases1[3..6], &read_ext_bases1[6..10], &read_ext_bases1[11..]);
    let read_bases4 = format!("{REF_BUFFER}{read_ext_bases4}");
    let cigar = make_cigar_string(&read_bases1, 0, 0);
    let read1 = read_100(&read_bases1, &cigar);
    let read2 = read_100(&read_bases1, &cigar);
    let read3 = read_100(&read_bases3, &cigar);
    let read4 = read_100(&read_bases4, &cigar);
    let read_parse_states = states(true, &[&read1, &read2, &read3, &read4], REF_BUFFER.len());
    let consensus_base_length = max_read_length(&read_parse_states) - REF_BUFFER.len() as i32;
    let seq_builder = SequenceBuilder::new_ext(&read_parse_states, true, consensus_base_length, true);
    assert_eq!(read_ext_bases1, seq_builder.base_string());
    assert!(seq_builder.repeats().is_empty());
    let read_state3 = &read_parse_states[2];
    assert_eq!(3, read_state3.mismatches().len());
    assert_eq!(SequenceDiffType::Delete, read_state3.mismatches()[1].type_);
    assert_eq!(SequenceDiffType::Insert, read_state3.mismatches()[2].type_);
    let read_state4 = &read_parse_states[3];
    assert_eq!(3, read_state4.mismatches().len());
    assert_eq!(SequenceDiffType::Insert, read_state4.mismatches()[0].type_);
    assert_eq!(SequenceDiffType::Insert, read_state4.mismatches()[1].type_);
    assert_eq!(SequenceDiffType::Delete, read_state4.mismatches()[2].type_);
}

#[test]
fn test_seq_build_low_qual_indel() {
    let _env = test_env();
    let read_ext_bases1 = "AACCGGTTAACC";
    let read_bases1 = format!("{REF_BUFFER}{read_ext_bases1}");
    // a low-qual 2-base insert, should not be treated as a low-qual SNV mismatch
    let read_ext_bases3 = format!("{}AT{}", &read_ext_bases1[0..6], &read_ext_bases1[6..9]);
    let read_bases3 = format!("{REF_BUFFER}{read_ext_bases3}");
    let cigar = make_cigar_string(&read_bases1, 0, 0);
    let read1 = read_100(&read_bases1, &cigar);
    let read2 = read_100(&read_bases1, &cigar);
    let read3 = read_100(&read_bases3, &make_cigar_string(&read_bases3, 0, 0));
    set_qual(&read3, 11, 10);
    set_qual(&read3, 14, 10); // only one of the indel bases needs to be low-qual
    let read_parse_states = states(true, &[&read1, &read2, &read3], REF_BUFFER.len());
    let consensus_base_length = max_read_length(&read_parse_states) - REF_BUFFER.len() as i32;
    let seq_builder = SequenceBuilder::new_ext(&read_parse_states, true, consensus_base_length, false);
    assert_eq!(read_ext_bases1, seq_builder.base_string());
    assert!(seq_builder.repeats().is_empty());
    let read_state3 = &read_parse_states[2];
    assert_eq!(1, read_state3.mismatches().len());
    let seq_diff_info = &read_state3.mismatches()[0];
    assert_eq!(SequenceDiffType::Insert, seq_diff_info.type_);
    assert_eq!(2, seq_diff_info.indel_length);
    assert_close(0.0, seq_diff_info.mismatch_penalty, 0.1);

    // again moving the other direction
    let read_bases1 = format!("{read_ext_bases1}{REF_BUFFER}");
    let read_ext_bases3 = format!("{}AT{}", &read_ext_bases1[3..6], &read_ext_bases1[6..]);
    let read_bases3 = format!("{read_ext_bases3}{REF_BUFFER}");
    let read1 = read_100(&read_bases1, &cigar);
    let read2 = read_100(&read_bases1, &cigar);
    let read3 = read_100(&read_bases3, &make_cigar_string(&read_bases3, 0, 0));
    set_qual(&read3, 4, 10);
    let read_parse_states = vec![
        ReadParseState::new(false, &read1, read_ext_bases1.len() as i32 - 1),
        ReadParseState::new(false, &read2, read_ext_bases1.len() as i32 - 1),
        ReadParseState::new(false, &read3, read_ext_bases3.len() as i32 - 1),
    ];
    let consensus_base_length = read_ext_bases1.len() as i32 + 1;
    let seq_builder = SequenceBuilder::new_ext(&read_parse_states, false, consensus_base_length, false);
    assert_eq!(read_ext_bases1, seq_builder.base_string());
    assert!(seq_builder.repeats().is_empty());
    let read_state3 = &read_parse_states[2];
    assert_eq!(1, read_state3.mismatches().len());
    let seq_diff_info = &read_state3.mismatches()[0];
    assert_eq!(SequenceDiffType::Insert, seq_diff_info.type_);
    assert_eq!(2, seq_diff_info.indel_length);
    assert_close(0.0, seq_diff_info.mismatch_penalty, 0.1);
}

#[test]
fn test_seq_build_repeat_mismatch_all_combos() {
    let _env = test_env();
    let read_pre_bases = "GCGTAACCG";
    let read_post_bases = "GGCAGGTG";
    let repeats = ["A", "AT", "AAT", "GGCC", "AATTT"];
    let consensus_repeat_count = 6;
    let read_repeat_counts = [consensus_repeat_count, consensus_repeat_count, 4, 5, 7, 8];
    for build_forwards in [true, false] {
        for repeat in repeats {
            let mut reads: Vec<ReadP> = Vec::new();
            let mut read_states: Vec<RpsP> = Vec::new();
            for repeat_count in read_repeat_counts {
                let mut read_bases = format!("{}{}{}", read_pre_bases, repeat.repeat(repeat_count), read_post_bases);
                if build_forwards {
                    read_bases = format!("{REF_BUFFER}{read_bases}");
                } else {
                    read_bases.push_str(REF_BUFFER);
                }
                let read = read_100(&read_bases, &make_cigar_string(&read_bases, 0, 0));
                let read_index_start = if build_forwards { REF_BUFFER.len() } else { read_bases.len() - REF_BUFFER.len() - 1 };
                read_states.push(ReadParseState::new(build_forwards, &read, read_index_start as i32));
                reads.push(read);
            }
            let consensus_base_length = max_read_length(&read_states) - REF_BUFFER.len() as i32;
            let seq_builder = SequenceBuilder::new(&read_states, build_forwards, consensus_base_length);
            let first_read_bases = reads[0].get_bases_string();
            let consensus_bases = if build_forwards { &first_read_bases[REF_BUFFER.len()..] } else { &first_read_bases[..first_read_bases.len() - REF_BUFFER.len()] };
            assert_eq!(consensus_bases, seq_builder.base_string());
            assert_eq!(1, seq_builder.repeats().len());
            let consensus_repeat = &seq_builder.repeats()[0];
            assert_eq!(consensus_repeat_count as i32, consensus_repeat.count);
            for (r, repeat_count) in read_repeat_counts.iter().enumerate() {
                let read_state = &read_states[r];
                if *repeat_count == consensus_repeat_count {
                    assert!(read_state.mismatches().is_empty());
                } else {
                    assert_eq!(1, read_state.mismatches().len());
                    let seq_diff_info = &read_state.mismatches()[0];
                    let read_repeat_index_start = seq_diff_info.repeat_index(consensus_repeat, build_forwards, true);
                    let read_repeat_index_end = seq_diff_info.repeat_index(consensus_repeat, build_forwards, false);
                    assert_eq!(read_repeat_index_start + (*repeat_count * repeat.len()) as i32 - 1, read_repeat_index_end);
                }
            }
        }
    }
}

#[test]
fn test_seq_build_repeat_misc() {
    let _env = test_env();
    let read_ext_bases1 = "CCGGAAAACCGG"; // 4 As but mismatch in read 3 occurs after only 3, below the min repeat count
    let read_bases1 = format!("{REF_BUFFER}{read_ext_bases1}");
    let read_ext_bases3 = format!("{}{}", &read_ext_bases1[0..7], &read_ext_bases1[8..]);
    let read_bases3 = format!("{REF_BUFFER}{read_ext_bases3}");
    let cigar = make_cigar_string(&read_bases1, 0, 0);
    let read1 = read_100(&read_bases1, &cigar);
    let read2 = read_100(&read_bases1, &cigar);
    let read3 = read_100(&read_bases3, &make_cigar_string(&read_bases3, 0, 0));
    let read_parse_states = states(true, &[&read1, &read2, &read3], REF_BUFFER.len());
    let consensus_base_length = max_read_length(&read_parse_states) - REF_BUFFER.len() as i32;
    let seq_builder = SequenceBuilder::new_ext(&read_parse_states, true, consensus_base_length, true);
    assert_eq!(read_ext_bases1, seq_builder.base_string());
    assert_eq!(1, seq_builder.repeats().len());
    let consensus_repeat = &seq_builder.repeats()[0];
    assert_eq!(4, consensus_repeat.count);
    assert_eq!(4, consensus_repeat.index);
    assert_eq!(7, consensus_repeat.last_index());
    let read_state3 = &read_parse_states[2];
    assert_eq!(1, read_state3.mismatches().len());
    assert_eq!(SequenceDiffType::Repeat, read_state3.mismatches()[0].type_);
    assert_eq!(3, read_state3.mismatches()[0].repeat_count);
}

fn check_repeat(repeat_info: Option<RepeatInfo>, bases: &str, count: i32, index: Option<i32>) {
    let repeat_info = repeat_info.expect("repeat");
    assert_eq!(bases, repeat_info.bases);
    assert_eq!(count, repeat_info.count);
    if let Some(index) = index {
        assert_eq!(index, repeat_info.index);
    }
}

#[test]
fn test_repeat_types() {
    let _env = test_env();
    // single base repeats
    let bases = b"AAACCTTTTT";
    let mut min_repeat_count = 4;
    // first check limits
    check_repeat(find_single_base_repeat(bases, 5, min_repeat_count), "T", 5, None);
    check_repeat(find_single_base_repeat(bases, 6, min_repeat_count), "T", 4, None);
    assert!(find_single_base_repeat(bases, 7, min_repeat_count).is_none());
    let bases = b"AAACCTTTTGAAAAAATGC";
    assert!(find_single_base_repeat(bases, 0, min_repeat_count).is_none());
    check_repeat(find_single_base_repeat(bases, 5, min_repeat_count), "T", 4, Some(5));
    check_repeat(find_single_base_repeat(bases, 10, min_repeat_count), "A", 6, Some(10));
    min_repeat_count = 3;
    let bases = b"ATATGATATATGGAT";
    assert!(find_dual_base_repeat(bases, 0, min_repeat_count).is_none());
    check_repeat(find_dual_base_repeat(bases, 5, min_repeat_count), "AT", 3, Some(5));
    let bases = b"ATGATGGGCCATGATGATGATG";
    check_repeat(find_multi_base_repeat(bases, 0, 3, 2), "ATG", 2, Some(0));
    check_repeat(find_multi_base_repeat(bases, 10, 3, min_repeat_count), "ATG", 4, Some(10));
    let bases = b"AATTAAGGCCTTCCTTCCTTCGTT";
    assert!(find_dual_dual_repeat(bases, 0, min_repeat_count).is_none());
    check_repeat(find_dual_dual_repeat(bases, 8, min_repeat_count), "CCTT", 3, Some(8));
    // five-base repeat
    let bases = b"ACTGTACTGTACTGTACTGTAATTGGTTGGAAGGT";
    check_repeat(find_multi_base_repeat(bases, 0, 5, min_repeat_count), "ACTGT", 4, None);
    // all together
    let bases = b"ATTTTTAACTCTCTAAAGTCGTCGTCAACCTTCCTTCCTTCGTT";
    let repeats = find_repeats_min(bases, min_repeat_count);
    assert_eq!(5, repeats.len());
    for (i, (repeat_bases, count, index)) in [("T", 5, 1), ("CT", 3, 8), ("A", 3, 14), ("GTC", 3, 17), ("CCTT", 3, 28)].iter().enumerate() {
        check_repeat(Some(repeats[i].clone()), repeat_bases, *count, Some(*index));
    }
}

fn compare(first_bases: &str, first_repeats: &[RepeatInfo], second_bases: &str, second_repeats: &[RepeatInfo], first_quals: &[u8], second_quals: &[u8]) -> f64 {
    compare_sequences(
        first_bases.as_bytes(), first_quals, 0, first_quals.len() as i32 - 1, first_repeats, second_bases.as_bytes(), second_quals, 0, second_quals.len() as i32 - 1,
        second_repeats, NO_MISMATCH_LIMIT as f64,
    )
}

#[test]
fn test_sequence_repeat_comparisons() {
    let _env = test_env();
    // the repeat multiple is just sufficient to count for one of the sequences only
    let first_bases = "ACGTA CTCTCT ACGT GAGA ACGT".replace(' ', "");
    let second_bases = "ACGTA CTCT ACGT GAGAGA ACGT".replace(' ', "");
    let first_repeats = find_repeats_min(first_bases.as_bytes(), 3);
    assert_eq!(1, first_repeats.len());
    let first_base_quals = build_default_base_quals(first_bases.len());
    let second_repeats = find_repeats_min(second_bases.as_bytes(), 3);
    assert_eq!(1, second_repeats.len());
    let second_base_quals = build_default_base_quals(second_bases.len());
    for _ in 0..2 {
        assert_close(2.0, compare(&first_bases, &first_repeats, &second_bases, &second_repeats, &first_base_quals, &second_base_quals), 0.1);
    }
}

#[test]
fn test_sequence_comparisons() {
    let _env = test_env();
    let first_bases = "ATTTTT   AACTCTCTCT AGA GG CTGACGTA   TTCC".replace(' ', "");
    let second_bases = "ATTTTTTT   AACTCTCT   AGA  CTGACGTA G TTCC".replace(' ', "");
    let first_repeats = find_repeats_min(first_bases.as_bytes(), 3);
    assert_eq!(2, first_repeats.len());
    let first_base_quals = build_default_base_quals(first_bases.len());
    let second_repeats = find_repeats_min(second_bases.as_bytes(), 3);
    assert_eq!(2, second_repeats.len());
    let second_base_quals = build_default_base_quals(second_bases.len());
    // diffs are extra Ts in the second, the 'CT' repeat in first, then extra 'GG' in first then extra 'G' in second
    for _ in 0..2 {
        assert_close(4.0, compare(&first_bases, &first_repeats, &second_bases, &second_repeats, &first_base_quals, &second_base_quals), 0.1);
    }

    // a single SNV at the start of a homopolymer can look like two 1-based INDELs or a repeat count diff
    // but only count this as one mismatch if it can be explained by a single SNV
    let first_bases = "AACGTTTTTAGCTGA";
    let second_bases = "AACTTTTTTAGCTGA";
    let first_repeats = find_repeats(first_bases.as_bytes());
    let second_repeats = find_repeats(second_bases.as_bytes());
    let base_quals = build_default_base_quals(second_bases.len());
    for _ in 0..2 {
        assert_close(1.0, compare(first_bases, &first_repeats, second_bases, &second_repeats, &base_quals, &base_quals), 0.1);
    }
}

#[test]
fn test_longer_sequence_comparisons() {
    let _env = test_env();
    let first_bases = "TTTTTTGTATTAAGTCTAATA C TTTTTTT  AACTTAAGTGTAGATTTTTTT AAA  TGCTCC A TAA C GGT T TTATTTATA C GATTTTTGTCACTG".replace(' ', "");
    let second_bases = "TTTTTTGTATTAAGTCTAATA G TTTTTTTT AACTTAAGTGTAGATTTTTT  AAAA TGCTCC G TAA T GGT G TTATTTATA T GATTTTTGTCACTG".replace(' ', "");
    let first_repeats = find_repeats(first_bases.as_bytes());
    assert_eq!(5, first_repeats.len());
    let first_base_quals = build_default_base_quals(first_bases.len());
    let second_repeats = find_repeats(second_bases.as_bytes());
    assert_eq!(5, second_repeats.len());
    let second_base_quals = build_default_base_quals(second_bases.len());
    for _ in 0..2 {
        assert_close(8.0, compare(&first_bases, &first_repeats, &second_bases, &second_repeats, &first_base_quals, &second_base_quals), 0.1);
    }
}

#[test]
fn test_assembly_ref_bases_trimmed() {
    let _env = test_env();
    let pos_junction = new_junction(CHR_1, 60, FORWARD);
    let extension_sequence = "ACGTTCGTAAAAAAGGGGGGACGTACGTCCCC";
    let ref_base_sequence = "ACGTAGAGAGAGACGTCCCCACGG";
    let assembly_sequence = format!("{ref_base_sequence}{extension_sequence}");
    let base_quals = build_default_base_quals(assembly_sequence.len());
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 37, &assembly_sequence, "24M32S");
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 38, &assembly_sequence[1..], "23M32S");
    let assembly = JunctionAssembler::new(pos_junction, None, None).process_junction(&[read1, read2])[0];
    assert_eq!(Some("ACGT_AG4_ACGT_C4_ACGG"), assembly.ref_bases_repeated_trimmed()); // 4 + 4 + 4 + 2 + 4
    assert_eq!(18, assembly.ref_base_trim_length());
    assert_eq!("ACGT_AG4_ACGT", build_trimmed_ref_base_sequence(&assembly, 12));

    let neg_junction = new_junction(CHR_1, 60, REVERSE);
    let assembly_sequence = format!("{extension_sequence}{ref_base_sequence}");
    let assembly = JunctionAssembly::new_test(neg_junction, assembly_sequence.as_bytes(), &base_quals, extension_sequence.len() as i32);
    assembly.m().build_repeat_info();
    assert_eq!(Some("ACGT_AG4_ACGT_C4_ACGG"), assembly.ref_bases_repeated_trimmed());
}
