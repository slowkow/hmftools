//! ReadUtilsTest
use super::*;
use crate::assembly::read::read_utils::INVALID_INDEX;
use crate::assembly::read_parse_state::{ReadParseState, RpsP};
use crate::assembly::sam::CigarOp::{D, I, M};

fn base_at(index: usize) -> u8 {
    REF_BASES_200.as_bytes()[index]
}

fn check(read_state: &RpsP, ref_position: i32, read_index: i32, op: CigarOp) {
    assert_eq!(ref_position, read_state.ref_position());
    assert_eq!(read_index, read_state.read_index());
    assert_eq!(Some(op), read_state.operator());
}

fn move_next(read_state: &RpsP, count: usize) {
    for _ in 0..count {
        read_state.m().move_next();
    }
}

/// htsjdk SAMRecord.getReadPositionAtReferencePosition: 1-based read position, 0 if the position is not an aligned base
fn get_read_position_at_reference_position(record: &SamRecord, pos: i32) -> i32 {
    let mut ref_pos = record.alignment_start;
    let mut read_pos = 1;
    for element in &record.cigar {
        match element.op {
            CigarOp::M | CigarOp::EQ | CigarOp::X => {
                if pos >= ref_pos && pos < ref_pos + element.len {
                    return read_pos + (pos - ref_pos);
                }
                ref_pos += element.len;
                read_pos += element.len;
            }
            CigarOp::I | CigarOp::S => read_pos += element.len,
            CigarOp::D | CigarOp::N => ref_pos += element.len,
            _ => {}
        }
    }
    0
}

#[test]
fn test_forward_read_parse_state() {
    let _env = test_env();
    let read = create_read(&READ_ID_GENERATOR.next_id(), 11, &REF_BASES_200[1..41], "10S30M");
    let read_state = ReadParseState::new_cigar(true, &read, 10, true); // from the soft-clip length
    check(&read_state, 11, 10, M);
    assert_eq!(base_at(11), read_state.current_base());
    move_next(&read_state, 29);
    check(&read_state, 40, 39, M);
    assert_eq!(base_at(40), read_state.current_base());
    assert!(!read_state.exhausted());
    move_next(&read_state, 1);
    assert!(read_state.exhausted());

    // again with indels
    let read = create_read(&READ_ID_GENERATOR.next_id(), 1, &REF_BASES_200[1..36], "5M5D5M5I10M10S");
    let read_state = ReadParseState::new_cigar(true, &read, 0, true);
    check(&read_state, 1, 0, M);
    move_next(&read_state, 5);
    check(&read_state, 6, 4, D);
    // move through the delete
    move_next(&read_state, 1);
    check(&read_state, 7, 4, D);
    move_next(&read_state, 3);
    check(&read_state, 10, 4, D);
    move_next(&read_state, 1);
    check(&read_state, 11, 5, M);
    // move onto insert
    move_next(&read_state, 5);
    check(&read_state, 15, 10, I);
    move_next(&read_state, 1);
    check(&read_state, 15, 11, I);
    // through rest of insert
    move_next(&read_state, 4);
    check(&read_state, 16, 15, M);
    assert_eq!(base_at(16), read_state.current_base());
    // stop at final soft-clipping
    move_next(&read_state, 10);
    assert!(read_state.exhausted());
}

#[test]
fn test_reverse_read_state() {
    let _env = test_env();
    let read = create_read(&READ_ID_GENERATOR.next_id(), 1, &REF_BASES_200[1..41], "30M10S");
    let read_state = ReadParseState::new_cigar(false, &read, 29, true);
    check(&read_state, 30, 29, M);
    assert_eq!(base_at(30), read_state.current_base());
    move_next(&read_state, 29);
    check(&read_state, 1, 0, M);
    assert_eq!(base_at(1), read_state.current_base());
    assert!(!read_state.exhausted());
    move_next(&read_state, 1);
    assert!(read_state.exhausted());

    // again with indels
    let read = create_read(&READ_ID_GENERATOR.next_id(), 11, &REF_BASES_200[1..36], "10S10M5D5M5I5M");
    let read_junction_index = read.bases_length() - 1;
    let read_state = ReadParseState::new_cigar(false, &read, read_junction_index, true);
    check(&read_state, 35, 34, M);
    move_next(&read_state, 5);
    check(&read_state, 31, 29, I);
    // move through the insert
    move_next(&read_state, 1);
    check(&read_state, 31, 28, I);
    move_next(&read_state, 4);
    check(&read_state, 30, 24, M);
    // move to start of delete
    move_next(&read_state, 5);
    check(&read_state, 25, 20, D);
    // move past delete
    move_next(&read_state, 5);
    check(&read_state, 20, 19, M);
    // move to last aligned base
    move_next(&read_state, 9);
    check(&read_state, 11, 10, M);
    // stop at final soft-clipping
    move_next(&read_state, 1);
    assert!(read_state.exhausted());
}

#[test]
fn test_hard_clipped_reads() {
    let _env = test_env();
    let read_bases = format!("{}{}", &REF_BASES_200[1..16], &REF_BASES_200[22..53]);
    let read = create_read(&READ_ID_GENERATOR.next_id(), 1, &read_bases, "105H15M6D31M");
    let read_state = ReadParseState::new_cigar(true, &read, 15, true);
    check(&read_state, 22, 15, M);
    assert_eq!(base_at(22), read_state.current_base());
    move_next(&read_state, 1);
    check(&read_state, 23, 16, M);
    assert_eq!(base_at(23), read_state.current_base());

    let read_bases = format!("{}{}", &REF_BASES_200[1..32], &REF_BASES_200[38..53]);
    let read = create_read(&READ_ID_GENERATOR.next_id(), 1, &read_bases, "31M6D15M105H");
    let read_state = ReadParseState::new_cigar(false, &read, 31, true);
    check(&read_state, 38, 31, M);
    assert_eq!(base_at(38), read_state.current_base());
    move_next(&read_state, 1);
    check(&read_state, 37, 31, D);
    move_next(&read_state, 6);
    check(&read_state, 31, 30, M);
    assert_eq!(base_at(31), read_state.current_base());
}

#[test]
fn test_misc_read_functions() {
    let _env = test_env();
    let random = ref_bases_random_100();
    let cigar_str = "5S20M8S";
    let read = create_read(TEST_READ_ID, 20, &random[15..48], cigar_str);
    assert_eq!(15, read.unclipped_start());
    assert_eq!(47, read.unclipped_end());
    assert_eq!(3, read.cigar_elements().len());
    assert_eq!(cigar_str, read.cigar_string());
    assert_eq!(5, read.get_read_index_at_reference_position(20));
    assert_eq!(INVALID_INDEX, read.get_read_index_at_reference_position(19));
    assert_eq!(4, read.get_read_index_at_reference_position_ext(19, true));
    assert_eq!(0, read.get_read_index_at_reference_position_ext(15, true));
    assert_eq!(INVALID_INDEX, read.get_read_index_at_reference_position_ext(14, true));
    assert_eq!(24, read.get_read_index_at_reference_position(39));
    assert_eq!(INVALID_INDEX, read.get_read_index_at_reference_position(40));
    assert_eq!(25, read.get_read_index_at_reference_position_ext(40, true));
    assert_eq!(32, read.get_read_index_at_reference_position_ext(47, true));
    assert_eq!(INVALID_INDEX, read.get_read_index_at_reference_position_ext(48, true));

    // with INDELs
    let cigar_str = "4S20M10D10M5I10M3S";
    let read_bases = format!("{}{}GGGGG{}", &random[16..40], &random[50..60], &random[60..73]);
    let read = create_read(TEST_READ_ID, 20, &read_bases, cigar_str);
    assert_eq!(16, read.unclipped_start());
    assert_eq!(69, read.alignment_end());
    assert_eq!(72, read.unclipped_end());
    assert_eq!(7, read.cigar_elements().len());
    assert_eq!(52, read.bases_length());
    assert_eq!(0, read.get_read_index_at_reference_position_ext(16, true));
    assert_eq!(4, read.get_read_index_at_reference_position_ext(20, false));
    assert_eq!(5, read.get_read_index_at_reference_position_ext(21, false));
    assert_eq!(23, read.get_read_index_at_reference_position_ext(39, false)); // start of DEL
    assert_eq!(23, read.get_read_index_at_reference_position_ext(48, false)); // within DEL
    assert_eq!(24, read.get_read_index_at_reference_position_ext(50, false)); // end of DEL
    assert_eq!(33, read.get_read_index_at_reference_position_ext(59, false)); // end of next M
    assert_eq!(39, read.get_read_index_at_reference_position_ext(60, false)); // base after insert
    assert_eq!(48, read.get_read_index_at_reference_position_ext(69, false)); // final aligned base
    assert_eq!(49, read.get_read_index_at_reference_position_ext(70, true));
    assert_eq!(50, read.get_read_index_at_reference_position_ext(71, true));
    assert_eq!(51, read.get_read_index_at_reference_position_ext(72, true)); // final SC base

    // compare to the SAM record method, which is 1-based
    let record = read.bam_record();
    assert_eq!(5, get_read_position_at_reference_position(record, 20));
    assert_eq!(24, get_read_position_at_reference_position(record, 39)); // start of DEL
    assert_eq!(25, get_read_position_at_reference_position(record, 50)); // end of DEL
    assert_eq!(34, get_read_position_at_reference_position(record, 59)); // end of next M
    assert_eq!(40, get_read_position_at_reference_position(record, 60)); // base after insert
    assert_eq!(49, get_read_position_at_reference_position(record, 69)); // final aligned base
}

// release of reads when a junction group is finished (read::read_release; no Java counterpart)

fn linked_pair(id: &str) -> (ReadP, ReadP) {
    let bases = "A".repeat(50);
    let first = create_read_mate(id, CHR_1, 100, &bases, "50M", CHR_1, 300, true);
    let second = create_read_mate(id, CHR_1, 300, &bases, "50M", CHR_1, 100, false);
    first.m().bam_record_mut().set_flag(FLAG_READ_PAIRED, true);
    first.m().bam_record_mut().set_flag(FLAG_FIRST_OF_PAIR, true);
    second.m().bam_record_mut().set_flag(FLAG_READ_PAIRED, true);
    second.m().bam_record_mut().set_flag(FLAG_FIRST_OF_PAIR, false);
    second.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
    Read::make_read_links(&first, &second);
    assert!(first.mate_read().is_some() && second.mate_read().is_some());
    (first, second)
}

#[test]
fn test_release_unreferenced_reads() {
    use crate::assembly::read::read_release::release_unreferenced_reads;

    // a pair nothing else refers to is freed: the mate links are removed and nothing is retained
    let (first, second) = linked_pair("READ_01");
    let retained = release_unreferenced_reads(vec![first, second]);
    assert!(retained.is_empty());

    // a pair with one read in use: that read is retained and reports its mate from a copy
    let (first, second) = linked_pair("READ_02");
    let holder = first.clone();
    let (mate_chromosome, mate_start, mate_end, mate_fragment_end) = (holder.mate_chromosome(), holder.mate_alignment_start(), holder.mate_alignment_end(), holder.mate_fragment_end());
    let mate_mapped = holder.is_mate_mapped();
    let retained = release_unreferenced_reads(vec![first, second]);
    assert_eq!(1, retained.len());
    assert!(retained[0] == holder);
    assert!(holder.mate_read().is_none());
    assert_eq!(mate_chromosome, holder.mate_chromosome());
    assert_eq!(mate_start, holder.mate_alignment_start());
    assert_eq!(mate_end, holder.mate_alignment_end());
    assert_eq!(mate_fragment_end, holder.mate_fragment_end());
    assert_eq!(mate_mapped, holder.is_mate_mapped());
    drop(retained);
    assert_eq!(1, holder.strong_count());

    // both reads in use: they stay linked until neither is, then both are freed
    let (first, second) = linked_pair("READ_03");
    let (holder1, holder2) = (first.clone(), second.clone());
    let retained = release_unreferenced_reads(vec![first, second]);
    assert_eq!(2, retained.len());
    assert!(holder1.mate_read() == Some(&holder2) && holder2.mate_read() == Some(&holder1));
    drop(holder1);
    drop(holder2);
    let retained = release_unreferenced_reads(retained);
    assert!(retained.is_empty());

    // a read without a mate link
    let single = create_read("READ_04", 100, &"A".repeat(50), "50M");
    let holder = single.clone();
    let retained = release_unreferenced_reads(vec![single]);
    assert_eq!(1, retained.len());
    drop(holder);
    assert!(release_unreferenced_reads(retained).is_empty());
}
