//! Ported from the Java SAGE test suite: evidence/RealignmentTest.java

use super::util::*;
use crate::constants::*;
use crate::evidence::{check_realignment, raw_context, realigned_read_index_position, ReadContextCounter, Realigned};
use crate::matcher::{MatchType, Matcher};
use crate::read_context::{create_context, VariantReadContext};
use crate::refseq::RefSequence;
use crate::types::*;
use std::sync::Arc;

/// Java String.substring(start, end)
fn sub(s: &str, start: usize, end: usize) -> &str {
    &s[start..end]
}

/// VariantUtils.createSimpleVariant(position, ref, alt)
fn create_simple_variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

/// new RefSequence(0, bases)
fn ref_sequence(bases: &str) -> RefSequence {
    RefSequence::from_bases(0, bases.as_bytes().to_vec())
}

/// VariantReadContextBuilder(DEFAULT_FLANK_LENGTH).createContext
fn build_read_context(var: &SimpleVariant, read: &Read, var_index_in_read: i32, ref_seq: &RefSequence) -> VariantReadContext {
    create_context(var, read, var_index_in_read, ref_seq, DEFAULT_FLANK_LENGTH).expect("read context")
}

/// VariantUtils.createReadContextMatcher
fn create_read_context_matcher(ctx: &VariantReadContext) -> Matcher {
    Matcher::new(ctx, true, false)
}

/// VariantUtils.createReadCounter(0, readContext)
fn create_read_counter(ctx: &VariantReadContext, fixture: &Fixture) -> ReadContextCounter {
    ReadContextCounter::new(Arc::new(ctx.clone()), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &fixture.qc(), false)
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

// Java: testInsertRealignment
#[test]
fn test_insert_realignment() {
    let ref_bases = format!(
        "X{}{}{}",
        sub(REF_BASES_200, 0, 100),
        // 0123456789012345678901234567890123456789
        //             ->
        "GTTGTTGTTGTCGTTGTTGTTGTTGGTTTTTCTGAGACAGAGTC",
        sub(REF_BASES_200, 0, 100)
    );

    let ref_sequence = ref_sequence(&ref_bases);

    let inserted_bases = "TTGTTGTTGTTGGTTTTTC";

    let var_build_read_bases =
        format!("{}{}{}{}", sub(&ref_bases, 91, 101), sub(&ref_bases, 101, 114), inserted_bases, sub(&ref_bases, 114, 150));

    let read_cigar = "33M19I36M";

    let var_build_read = build_sam_record(91, read_cigar, &var_build_read_bases);
    let ref_ = sub(&ref_bases, 113, 114);

    let var = create_simple_variant(113, ref_, &format!("{ref_}{inserted_bases}"));
    let read_context = build_read_context(&var, &var_build_read, 22, &ref_sequence);

    assert_eq!(113, read_context.variant.pos);
    assert_eq!(11, read_context.var_index);
    assert_eq!("CGTTGTTGTTGTTGGTTTTTCTTGTTGTTGTTGGTTTTTCTGA", read_context.core_str());

    // variant is at read index 41 when should be 22
    let realigned_start_pos = 72;
    let read_cigar = build_cigar_string(var_build_read_bases.len());
    let realigned_read = build_sam_record(realigned_start_pos, &read_cigar, &var_build_read_bases);

    let matcher = create_read_context_matcher(&read_context);

    let realigned_read_index = realigned_read_index_position(&read_context, &realigned_read);
    let match_type = matcher.determine_read_match(&read_context, &realigned_read, realigned_read_index).kind;

    assert_eq!(MatchType::Full, match_type);
}

// Java: testSoftClipRealignment
#[test]
fn test_soft_clip_realignment() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let ref_bases = format!("{}{}{}", sub(REF_BASES_200, 0, 100), "GGGGATTGGTCTATCTATCTATCTAGG", sub(REF_BASES_200, 0, 100));
    //                                                             0123456789012345678901234567890123456789

    let ref_sequence = ref_sequence(&ref_bases);

    let inserted_bases = "TCTA";

    let var_position = 108;
    let read_pos_start = 89;
    let var_index_in_read = var_position - read_pos_start;
    let var_build_read_bases = format!(
        "{}{}{}",
        sub(&ref_bases, read_pos_start as usize, var_position as usize + 1),
        inserted_bases,
        sub(&ref_bases, var_position as usize + 1, 149)
    );

    let read_cigar = "20M4I40M";

    let var_build_read = build_sam_record(read_pos_start, read_cigar, &var_build_read_bases);
    let ref_ = sub(&ref_bases, var_position as usize, var_position as usize + 1);

    let var = create_simple_variant(var_position, ref_, &format!("{ref_}{inserted_bases}"));
    let read_context = build_read_context(&var, &var_build_read, var_index_in_read, &ref_sequence);

    assert_eq!(11, read_context.var_index);
    assert_eq!("TCTA".as_bytes(), read_context.homology.as_ref().unwrap().bases.as_slice());
    assert_eq!(126, read_context.core_position_end);
    assert_eq!("GGTCTATCTATCTATCTATCTAGG", read_context.core_str());

    let matcher = create_read_context_matcher(&read_context);

    // confirm the initial read is a full match
    let match_type = matcher.determine_read_match(&read_context, &var_build_read, var_index_in_read).kind;

    assert_eq!(MatchType::Full, match_type);

    // position:        SC       94
    // position:        SC       4567890123456789012
    // index:           0123456789012345678901234567
    let read_bases = format!("{}{}", "ATTGGTCTATCTATCTATCTATCTAGG", sub(REF_BASES_200, 0, 20));
    let realigned_start_pos = read_pos_start + 20;
    let read_cigar = "9S38M";
    let realigned_read = build_sam_record(realigned_start_pos, read_cigar, &read_bases);

    let (read_variant_index, _) = raw_context(&var, &realigned_read);
    assert_eq!(8, read_variant_index);

    let realigned_read_index = realigned_read_index_position(&read_context, &realigned_read);
    assert_eq!(4, realigned_read_index);
    let match_type = matcher.determine_read_match(&read_context, &realigned_read, realigned_read_index).kind;

    assert_eq!(MatchType::Full, match_type);

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    read_context_counter.process_read(&realigned_read, 0, None, &qc);

    assert_eq!(1, read_context_counter.counts.realigned);
}

const REF_BASE_START: &str = "ACGTTGCAACGTTGCAGGGG"; // 20 bases

// Java: testInsertRealignment2
#[test]
fn test_insert_realignment2() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let ref_bases = format!("{}{}{}", REF_BASE_START, "CTTTCTTTTCTTTCTTTTTCTTTA", sub(REF_BASES_200, 0, 50));
    // index                                           0123456789012345678901234567890123456789
    // pos                                             20        30        40

    let ref_sequence_ = ref_sequence(&format!("{}{}", ref_bases, sub(REF_BASES_200, 50, 100)));

    let inserted_bases = "TTTCTTTC";

    let var_position = 35;
    let read_pos_start = 20;
    let var_index_in_read = var_position - read_pos_start;
    let var_build_read_bases = format!(
        "{}{}{}",
        sub(&ref_bases, read_pos_start as usize, var_position as usize + 1),
        inserted_bases,
        &ref_bases[var_position as usize + 1..]
    );

    let read_cigar = "16M8I20M";

    let var_build_read = build_sam_record(read_pos_start, read_cigar, &var_build_read_bases);
    let ref_ = sub(&ref_bases, var_position as usize, var_position as usize + 1).to_string();

    let var = create_simple_variant(var_position, &ref_, &format!("{ref_}{inserted_bases}"));
    let read_context = build_read_context(&var, &var_build_read, var_index_in_read, &ref_sequence_);

    assert_eq!(12, read_context.var_index);
    assert_eq!("TTTCTTT".as_bytes(), read_context.homology.as_ref().unwrap().bases.as_slice());
    assert_eq!(44, read_context.core_position_end);
    assert_eq!("CTTTTTCTTTCTTTCTTTAC", read_context.core_str());

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    read_context_counter.process_read(&var_build_read, 0, None, &qc);
    assert_eq!(1, read_context_counter.counts.full);

    let (read_variant_index, _) = raw_context(&var, &var_build_read);
    assert_eq!(15, read_variant_index);

    let matcher = create_read_context_matcher(&read_context);

    // test 1: basic realignment
    let read_bases = format!("{}{}", "CTTTCTTTTTCTTTCTTTCTTTA", sub(REF_BASES_200, 0, 20)); // 17 + 20
    let realigned_start_pos = read_pos_start;
    let read_cigar = "5M1I9M2D22M";
    let realigned_read = build_sam_record(realigned_start_pos, read_cigar, &read_bases);

    let realigned_read_index = realigned_read_index_position(&read_context, &realigned_read);
    let match_type = matcher.determine_read_match(&read_context, &realigned_read, realigned_read_index).kind;

    assert_eq!(MatchType::Full, match_type);

    read_context_counter.process_read(&realigned_read, 0, None, &qc);

    assert_eq!(1, read_context_counter.counts.realigned);

    // test 2: jitter realignment
    let read_bases = format!("{}{}", "CTTTCTTTTTCTTTCTTTCTTTCTTTA", sub(REF_BASES_200, 0, 20));
    let read_cigar = "9M1I10M2I25M";
    let realigned_read = build_sam_record(realigned_start_pos, read_cigar, &read_bases);

    let realigned_read_index = realigned_read_index_position(&read_context, &realigned_read);
    let mut matcher = create_read_context_matcher(&read_context);
    let realigned_type = check_realignment(&read_context, &mut matcher, &realigned_read, 16, realigned_read_index, None, false);

    // works if realignedReadIndex is 6 rather than 10, since the 1xTATT lengthened jitter causes the realignedIndex to be off
    assert!(realigned_type == Realigned::Lengthened);

    let ref_bases = format!("{}{}{}", REF_BASE_START, "TTCTTTTCTTTCTTTTTCTTTA", sub(REF_BASES_200, 0, 60));
    // index                                           0123456789012345678901234567890123456789
    // pos                                             20        30        40

    let ref_sequence_ = ref_sequence(&format!("{}{}", ref_bases, sub(REF_BASES_200, 60, 110)));

    let var_position = 33;
    let ref_ = sub(&ref_bases, var_position as usize, var_position as usize + 1).to_string();
    let var = create_simple_variant(var_position, &ref_, &format!("{ref_}{inserted_bases}"));

    let read_pos_start = 20; // drop back to make the flank long enough
    let var_index_in_read = var_position - read_pos_start;
    let var_build_read_bases = format!(
        "{}{}{}",
        sub(&ref_bases, read_pos_start as usize, var_position as usize + 1),
        inserted_bases,
        &ref_bases[var_position as usize + 1..]
    );
    let read_cigar = "12M8I30M";
    let var_build_read = build_sam_record(read_pos_start, read_cigar, &var_build_read_bases);

    let read_context = build_read_context(&var, &var_build_read, var_index_in_read, &ref_sequence_);

    // 1:33 T>TTTTCTTTC read(TCTTTTCTTTCTTTTTCTTTCTTTCTTTACGCAATATTCG 11M8I20M) pos(21-0) index(10-12-29) repeat(14: TTTC-3) homology(TTTCTTT length(7)) alt(12-20) ref(CTTTTTCTTTAC)
    let read_bases = format!("{}{}", "TCTTTTTCTTTCTTTCTTTA", sub(REF_BASES_200, 0, 20)); // 20 + 20
    let realigned_start_pos = read_pos_start;
    let read_cigar = "12M2D28M";
    let realigned_read = build_sam_record(realigned_start_pos, read_cigar, &read_bases);

    let realigned_read_index = realigned_read_index_position(&read_context, &realigned_read);
    let matcher = create_read_context_matcher(&read_context);
    let match_type = matcher.determine_read_match(&read_context, &realigned_read, realigned_read_index).kind;
    assert_eq!(MatchType::Full, match_type);

    let mut read_context_counter = create_read_counter(&read_context, &fixture);
    read_context_counter.process_read(&realigned_read, 0, None, &qc);
    assert_eq!(1, read_context_counter.counts.realigned);

    // Java does not recompute ref for the new position
    let var_position = 36;
    let var_index_in_read = var_position - read_pos_start;
    let var = create_simple_variant(var_position, &ref_, &format!("{ref_}{inserted_bases}"));
    let read_context = build_read_context(&var, &var_build_read, var_index_in_read, &ref_sequence_);

    let read_bases = format!("{}{}", "CTTTCTTTTTCTTTCTTTCTTTA", sub(REF_BASES_200, 0, 20)); // 23 + 20
    let realigned_start_pos = 26;
    let read_cigar = "6S9M2D8M";
    let realigned_read = build_sam_record(realigned_start_pos, read_cigar, &read_bases);

    let realigned_read_index = realigned_read_index_position(&read_context, &realigned_read);
    let matcher = create_read_context_matcher(&read_context);
    let _match_type = matcher.determine_read_match(&read_context, &realigned_read, realigned_read_index).kind;
    // Java: // assertEquals(FULL, matchType);

    let mut read_context_counter = create_read_counter(&read_context, &fixture);
    read_context_counter.process_read(&realigned_read, 0, None, &qc);
    // Java: //  assertEquals(1, readContextCounter.readCounts().Realigned); // TODO
}

// Java: testDelBeforeUnclippedStart
#[test]
fn test_del_before_unclipped_start() {
    let fixture = Fixture::new();
    let qc = fixture.qc();

    let ref_bases = format!("{}{}{}", sub(REF_BASES_200, 0, 100), "GGGGATTGGTCTATCTATCTATCTAGG", sub(REF_BASES_200, 0, 100));
    //                                                             0123456789012345678901234567890123456789

    let ref_sequence = ref_sequence(&ref_bases);

    let deleted_bases = "TCTATCTA";

    let var_position = 108;
    let read_pos_start = 89;
    let var_index_in_read = var_position - read_pos_start;
    let var_build_read_bases = format!(
        "{}{}",
        sub(&ref_bases, read_pos_start as usize, var_position as usize + 1),
        sub(&ref_bases, var_position as usize + 1 + deleted_bases.len(), 149)
    );

    let read_cigar = "20M8D32M";

    let var_build_read = build_sam_record(read_pos_start, read_cigar, &var_build_read_bases);
    let ref_ = sub(&ref_bases, var_position as usize, var_position as usize + 1);

    let var = create_simple_variant(var_position, &format!("{ref_}{deleted_bases}"), ref_);
    let read_context = build_read_context(&var, &var_build_read, var_index_in_read, &ref_sequence);

    let mut read_context_counter = create_read_counter(&read_context, &fixture);

    let read_bases = format!("{}{}", "TGGTCTATCTAGG", sub(REF_BASES_200, 0, 20));
    let read_cigar = "3S30M";
    let realigned_pos = var_position + deleted_bases.len() as i32 + 1;
    let realigned_read = build_sam_record(realigned_pos, read_cigar, &read_bases);

    // core end position for this read is 126
    // index in read for this should be:
    // 3S   index 0-2   pre position 116
    // 30M  index 3-32  positions 116-145

    read_context_counter.process_read(&realigned_read, 0, None, &qc);
    assert_eq!(1, read_context_counter.counts.realigned);

    // check low-qual base check has been shifted correctly
    let read_bases_with_error = format!("{}A{}", sub(&read_bases, 0, 2), &read_bases[3..]);
    let mut realigned_read = build_sam_record(realigned_pos, read_cigar, &read_bases_with_error);
    realigned_read.quals[2] = 10;

    // should fail since now a low-qual
    read_context_counter.process_read(&realigned_read, 0, None, &qc);
    // Java: // assertEquals(2, readContextCounter.readCounts().Realigned);

    let read_bases_with_error = format!("{}A{}", sub(&read_bases, 0, 1), &read_bases[2..]);
    let mut realigned_read = build_sam_record(realigned_pos, read_cigar, &read_bases_with_error);
    realigned_read.quals[1] = 10;

    // this one should be permitted as a low-qual pass
    read_context_counter.process_read(&realigned_read, 0, None, &qc);
    // Java: // assertEquals(2, readContextCounter.readCounts().Realigned);
}

// Java: testRealignmentIndexOffset
#[test]
fn test_realignment_index_offset() {
    let ref_bases = format!(
        "X{}{}{}",
        sub(REF_BASES_200, 0, 100),
        // 0123456789012345678901234567890123456789
        //             ->
        "GATTTTTTTTTTGA",
        sub(REF_BASES_200, 0, 100)
    );

    let ref_sequence = ref_sequence(&ref_bases);

    let inserted_bases = "C";

    let var_build_read_bases =
        format!("{}{}{}{}", sub(&ref_bases, 91, 101), sub(&ref_bases, 101, 105), inserted_bases, sub(&ref_bases, 105, 150));

    let read_cigar = "14M1I45M";

    let var_build_read = build_sam_record(91, read_cigar, &var_build_read_bases);
    let ref_ = sub(&ref_bases, 104, 105);

    let var = create_simple_variant(104, ref_, &format!("{ref_}{inserted_bases}"));
    let read_context = build_read_context(&var, &var_build_read, 13, &ref_sequence);

    assert_eq!(104, read_context.variant.pos);
    assert_eq!(11, read_context.var_index);
    assert_eq!("TTCTTTTTTTTG", read_context.core_str());

    // read has an extra T at start of poly-T rather than C in middle of poly-T
    let poly_t_var_build_read_bases =
        format!("{}{}{}{}", sub(&ref_bases, 91, 101), sub(&ref_bases, 101, 105), "T", sub(&ref_bases, 105, 150));
    let poly_t_read_cigar = "12M1I47M";
    let mut base_quals = default_quals(poly_t_var_build_read_bases.len());
    base_quals[14] = 11; // the base that was a C in the T>TC insert
    let poly_t_var_build_read = build_sam_record_quals(91, poly_t_read_cigar, &poly_t_var_build_read_bases, base_quals);

    let mut matcher = create_read_context_matcher(&read_context);
    let (read_var_index, _) = raw_context(&var, &poly_t_var_build_read);
    let realigned_read_index = realigned_read_index_position(&read_context, &poly_t_var_build_read);
    let realigned_type =
        check_realignment(&read_context, &mut matcher, &poly_t_var_build_read, read_var_index, realigned_read_index, None, false);

    assert!(realigned_type == Realigned::None);
}
