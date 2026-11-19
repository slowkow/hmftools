//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! evidence/FragmentSyncTest.java. `CombinedSyncData.formFragmentRead` is `sync::form_fragment_read`;
//! FragmentSyncOutcome (CombinedRecord + SyncType) is the `sync::Outcome` enum.

use super::util::*;
use crate::constants::*;
use crate::evidence::{FragmentData, ReadContextCounter};
use crate::read_context::VariantReadContext;
use crate::sync::{form_fragment_read, Outcome};
use crate::types::*;
use std::sync::Arc;

const READ_ID: &str = "READ_01";

/// MockRefGenome.generateRandomBases: not random, cycles through ACGT
fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect()
}

/// MockRefGenome.getNextBase
fn get_next_base(base: &str) -> String {
    match base {
        "A" => "C",
        "C" => "G",
        "G" => "T",
        "T" => "A",
        other => other,
    }
    .to_string()
}

/// FragmentSyncTest.REF_BASES
fn ref_bases() -> String {
    format!("X{}", generate_random_bases(100))
}

fn set_negative_strand(read: &mut Read, negative: bool) {
    if negative {
        read.flags |= FLAG_REVERSE;
    } else {
        read.flags &= !FLAG_REVERSE;
    }
}

/// assertNotNull(formFragmentRead(first, second).CombinedRecord)
fn combined_record(first: &Read, second: &Read) -> Read {
    match form_fragment_read(first, second) {
        Outcome::Combined(read) => read,
        Outcome::NoOverlap => panic!("expected a combined record, got NO_OVERLAP"),
        Outcome::Inversion => panic!("expected a combined record, got INVERSION"),
        Outcome::CigarMismatch => panic!("expected a combined record, got CIGAR_MISMATCH"),
        Outcome::Exception => panic!("expected a combined record, got EXCEPTION"),
    }
}

fn read_string(read: &Read) -> &str {
    std::str::from_utf8(&read.bases).unwrap()
}

const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";

/// VariantUtils.createReadContext(variant, leftCore, rightCore): no repeat or homology, default flanks
fn create_read_context(variant: &SimpleVariant, left_core: &str, right_core: &str) -> VariantReadContext {
    let (left_flank, right_flank) = (TEST_LEFT_FLANK, TEST_RIGHT_FLANK);
    let core_index_start = left_flank.len() as i32;
    let var_read_index = core_index_start + left_core.len() as i32;
    let core_index_end = var_read_index + variant.alt_len() - 1 + right_core.len() as i32;
    let ref_bases = format!("{}{}{}", left_core, variant.ref_str(), right_core);
    let read_bases = format!("{}{}{}{}{}", left_flank, left_core, variant.alt_str(), right_core, right_flank);
    let alignment_start = variant.pos - var_read_index;
    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + right_core.len() as i32 + right_flank.len() as i32 - 1;
    VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.into_bytes(),
        read_cigar: format!("{}M", read_bases.len()),
        read_bases: read_bases.into_bytes(),
        core_index_start,
        var_index: var_read_index,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start: variant.pos - left_core.len() as i32,
        core_position_end: variant.pos + right_core.len() as i32,
        has_indel_in_core: false,
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    }
}

// FragmentSyncTest.testCombinedRecords
#[test]
fn test_combined_records() {
    let r = ref_bases();

    // first a basic match
    let first = create_sam_record(READ_ID, 1, &r[1..21], "20M");

    let mut second = create_sam_record(READ_ID, 5, &r[5..25], "20M");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(1, combined.start);
    assert_eq!(24, combined.end);
    assert_eq!(&r[1..25], read_string(&combined));
    assert_eq!("24M", cigar_str(&combined.cigar));

    // soft-clips extending each end by varying amounts
    let first = create_sam_record(READ_ID, 6, &r[1..21], "5S10M5S");

    let mut second = create_sam_record(READ_ID, 12, &r[10..30], "2S16M2S");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(6, combined.start);
    assert_eq!(27, combined.end);
    assert_eq!(&r[1..30], read_string(&combined));
    assert_eq!("5S22M2S", cigar_str(&combined.cigar));

    // a delete
    let first_bases = format!("{}{}", &r[1..11], &r[12..22]);
    let first = create_sam_record(READ_ID, 1, &first_bases, "10M1D10M");

    let second_bases = format!("{}{}", &r[6..11], &r[12..27]);
    let mut second = create_sam_record(READ_ID, 6, &second_bases, "5M1D15M");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(1, combined.start);
    assert_eq!(26, combined.end);
    let combined_bases = format!("{}{}", &r[1..11], &r[12..27]);
    assert_eq!(combined_bases, read_string(&combined));
    assert_eq!("10M1D15M", cigar_str(&combined.cigar));

    // with a longer delete
    let first_bases = format!("{}{}", &r[1..11], &r[16..26]);
    let first = create_sam_record(READ_ID, 1, &first_bases, "10M5D10M");

    let second_bases = format!("{}{}", &r[6..11], &r[16..31]);
    let mut second = create_sam_record(READ_ID, 6, &second_bases, "5M5D15M");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(1, combined.start);
    assert_eq!(30, combined.end);
    let combined_bases = format!("{}{}", &r[1..11], &r[16..31]);
    assert_eq!(combined_bases, read_string(&combined));
    assert_eq!("10M5D15M", cigar_str(&combined.cigar));

    // multiple deletes
    let first_bases = format!("{}{}{}{}", &r[1..11], &r[16..26], &r[28..38], &r[41..46]);
    let first = create_sam_record(READ_ID, 1, &first_bases, "10M5D10M2D10M3D5M");

    let second_bases = format!("{}{}{}{}", &r[6..11], &r[16..26], &r[28..38], &r[41..51]);
    let mut second = create_sam_record(READ_ID, 6, &second_bases, "5M5D10M2D10M3D10M");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(1, combined.start);
    assert_eq!(50, combined.end);
    let combined_bases = format!("{}{}{}{}", &r[1..11], &r[16..26], &r[28..38], &r[41..51]);
    assert_eq!(combined_bases, read_string(&combined));
    assert_eq!("10M5D10M2D10M3D10M", cigar_str(&combined.cigar));

    // with an insert
    let first_bases = format!("{}CCC{}", &r[1..11], &r[11..21]);
    let first = create_sam_record(READ_ID, 1, &first_bases, "10M3I10M");

    let second_bases = format!("{}CCC{}", &r[6..11], &r[11..26]);
    let mut second = create_sam_record(READ_ID, 6, &second_bases, "5M3I15M");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(1, combined.start);
    assert_eq!(25, combined.end);
    let combined_bases = format!("{}CCC{}", &r[1..11], &r[11..26]);
    assert_eq!(combined_bases, read_string(&combined));
    assert_eq!("10M3I15M", cigar_str(&combined.cigar));

    // more complicated example
    let first_bases = format!("{}CCC{}{}AA{}", &r[1..11], &r[11..21], &r[26..36], &r[36..46]);

    let first = create_sam_record(READ_ID, 1, &first_bases, "10M3I10M5D10M2I5M5S");

    let second_bases = format!("{}CCC{}{}AA{}", &r[6..11], &r[11..21], &r[26..36], &r[36..51]);
    let mut second = create_sam_record(READ_ID, 6, &second_bases, "5M3I10M5D10M2I15M");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(1, combined.start);
    assert_eq!(50, combined.end);
    let combined_bases = format!("{}CCC{}{}AA{}", &r[1..11], &r[11..21], &r[26..36], &r[36..51]);

    assert_eq!(combined_bases, read_string(&combined));
    assert_eq!("10M3I10M5D10M2I15M", cigar_str(&combined.cigar));
}

// FragmentSyncTest.testCombinedRecords2
#[test]
fn test_combined_records2() {
    let r = ref_bases();

    let first = create_sam_record(READ_ID, 10, &format!("{}{}", &r[10..25], &r[35..55]), "15M10D20M");

    let mut second = create_sam_record(READ_ID, 45, &r[45..75], "20M10S");
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(10, combined.start);
    assert_eq!(64, combined.end);
    assert_eq!("15M10D30M10S", cigar_str(&combined.cigar));
    assert_eq!(format!("{}{}", &r[10..25], &r[35..75]), read_string(&combined));
}

// FragmentSyncTest.testSoftClips
#[test]
fn test_soft_clips() {
    let r = ref_bases();

    let mut first = create_sam_record(READ_ID, 20, &r[20..35], "10M5S");
    set_negative_strand(&mut first, true);

    let second = create_sam_record(READ_ID, 20, &r[17..30], "3S10M");

    let combined = combined_record(&first, &second);
    assert_eq!(20, combined.start);
    assert_eq!(29, combined.end);
    assert_eq!("3S10M5S", cigar_str(&combined.cigar));
    assert_eq!(&r[17..35], read_string(&combined));
}

// FragmentSyncTest.testReadStrandMismatches
#[test]
fn test_read_strand_mismatches() {
    // one read supports the variant, one doesn't and this factors into the read strand determination
    let r = ref_bases();
    let position = 20usize;

    let ref_base = &r[position..position + 1];
    let alt_base = get_next_base(ref_base);

    let variant = SimpleVariant::new(position as i32, ref_base.as_bytes(), alt_base.as_bytes());

    let read_context = create_read_context(&variant, &r[position - 2..position], &r[position + 1..position + 3]);

    // VariantUtils.createReadCounter(readContext)
    let fixture = Fixture::new();
    let qc = fixture.qc();
    let mut read_context_counter =
        ReadContextCounter::new(Arc::new(read_context), Tier::LowConfidence, DEFAULT_MAX_READ_DEPTH, 1, &qc, false);

    let mut first = create_sam_record(READ_ID, 10, &r[10..40], "30M");
    first.quals[10] = 11;

    let mut second =
        create_sam_record(READ_ID, 10, &format!("{}{}{}", &r[10..position], alt_base, &r[position + 1..40]), "30M");
    set_negative_strand(&mut second, true);

    let consensus_read = combined_record(&first, &second);

    for _ in 0..3 {
        let fragment = FragmentData { first: first.clone(), second: second.clone() };
        read_context_counter.process_read(&consensus_read, 1, Some(&fragment), &qc);
    }

    // swap and re-process
    set_negative_strand(&mut first, true);
    set_negative_strand(&mut second, false);

    let fragment = FragmentData { first: first.clone(), second: second.clone() };
    read_context_counter.process_read(&consensus_read, 1, Some(&fragment), &qc);

    assert_eq!(4, read_context_counter.alt_read_strand_bias.depth());
    assert!((0.25 - read_context_counter.alt_read_strand_bias.bias()).abs() <= 0.01);
}

// FragmentSyncTest.testSplitReads
#[test]
fn test_split_reads() {
    let r = ref_bases();

    // first a basic match
    let mut first = create_sam_record(READ_ID, 20, &format!("{}{}", &r[20..30], &r[60..80]), "10M30N20M");
    set_negative_strand(&mut first, true);

    let second = create_sam_record(READ_ID, 10, &format!("{}{}", &r[10..30], &r[60..70]), "20M30N10M");

    let combined = combined_record(&first, &second);
    assert_eq!(10, combined.start);
    assert_eq!(79, combined.end);
    assert_eq!(format!("{}{}", &r[10..30], &r[60..80]), read_string(&combined));
    assert_eq!("20M30N20M", cigar_str(&combined.cigar));
}

// FragmentSyncTest.testMismatches
#[test]
fn test_mismatches() {
    let r = ref_bases();

    // first a basic match
    let first = create_sam_record(READ_ID, 1, &format!("{}C{}", &r[1..11], &r[11..21]), "10M1I10M");

    let mut second = create_sam_record(READ_ID, 5, &r[5..25], "20M");
    set_negative_strand(&mut second, true);

    assert!(matches!(form_fragment_read(&first, &second), Outcome::CigarMismatch));

    // an inversion
    let first = create_sam_record(READ_ID, 10, &r[10..40], "30M");

    let second = create_sam_record(READ_ID, 15, &r[15..45], "30M");

    assert!(matches!(form_fragment_read(&first, &second), Outcome::Inversion));

    // off by 1
    let first = create_sam_record(READ_ID, 1, &format!("{}C{}", &r[1..11], &r[11..21]), "10M1I10M");

    let mut second = create_sam_record(READ_ID, 2, &format!("{}C{}", &r[1..12], &r[12..21]), "10M1I10M");
    set_negative_strand(&mut second, true);

    assert!(matches!(form_fragment_read(&first, &second), Outcome::CigarMismatch));

    // (the Java "too many mismatches" / BASE_MISMATCH case is commented out in the Java test)

    // INDEL vs aligned
    let first = create_sam_record(READ_ID, 10, &format!("{}{}", &r[10..20], &r[25..40]), "10M5D15M");

    let mut second = create_sam_record(READ_ID, 15, &r[15..45], "30M");
    set_negative_strand(&mut second, true);

    assert!(matches!(form_fragment_read(&first, &second), Outcome::CigarMismatch));

    // different INDELs
    let first = create_sam_record(READ_ID, 10, &format!("{}C{}", &r[10..20], &r[20..40]), "10M1I20M");

    let mut second = create_sam_record(READ_ID, 10, &format!("{}CC{}", &r[10..20], &r[20..40]), "10M2I20M");
    set_negative_strand(&mut second, true);

    assert!(matches!(form_fragment_read(&first, &second), Outcome::CigarMismatch));

    // non-overlapping but different INDELs
    let first = create_sam_record(READ_ID, 1, &format!("{}C{}", &r[1..6], &r[6..41]), "5M1I35M");

    let mut second = create_sam_record(READ_ID, 30, &format!("{}{}", &r[30..40], &r[45..75]), "10M5D30M");
    set_negative_strand(&mut second, true);

    assert!(matches!(form_fragment_read(&first, &second), Outcome::CigarMismatch));
}

// FragmentSyncTest.testTruncatedFragments
#[test]
fn test_truncated_fragments() {
    let r = ref_bases();

    let mut first = create_sam_record(READ_ID, 5, &r[5..35], "30M");

    first.insert_size = 25;
    set_negative_strand(&mut first, true);

    let mut second = create_sam_record(READ_ID, 10, &r[10..40], "30M");
    second.insert_size = -25;

    // the Java test reads CombinedRecord before asserting SyncType == COMBINED; both are the Combined arm here
    let combined = combined_record(&first, &second);
    assert_eq!(25, combined.insert_size);
    assert_eq!(10, combined.start);
    assert_eq!(34, combined.end);
    assert_eq!("25M", cigar_str(&combined.cigar));
    assert_eq!(25, combined.quals.len());

    assert_eq!(&r[10..35], read_string(&combined));

    // test again with soft-clips
    // actual alignment is 11 -> 42 = 22 bases
    let frag_length = 22;
    let mut first = create_sam_record(READ_ID, 11, &r[9..50], "2S35M4S");
    first.insert_size = frag_length;

    let mut second = create_sam_record(READ_ID, 11, &r[5..45], "6S32M2S");

    second.insert_size = -frag_length;
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(frag_length, combined.insert_size);
    assert_eq!(11, combined.start);
    assert_eq!(42, combined.end);
    assert_eq!("2S32M2S", cigar_str(&combined.cigar));
    assert_eq!(36, combined.quals.len());

    assert_eq!(&r[9..45], read_string(&combined));

    // test again with soft-clips on the 3' ends and a DEL

    let frag_length = 43 - 8 + 1;
    let mut first = create_sam_record(READ_ID, 8, &format!("{}{}", &r[8..23], &r[28..45]), "15M5D15M2S");
    first.insert_size = frag_length;

    let mut second = create_sam_record(READ_ID, 11, &format!("{}{}", &r[5..23], &r[28..43]), "6S12M5D15M");

    second.insert_size = -frag_length;
    set_negative_strand(&mut second, true);

    let combined = combined_record(&first, &second);
    assert_eq!(frag_length, combined.insert_size);
    assert_eq!(8, combined.start);
    assert_eq!(42, combined.end);
    assert_eq!("15M5D15M", cigar_str(&combined.cigar));
    assert_eq!(30, combined.quals.len());

    assert_eq!(format!("{}{}", &r[8..23], &r[28..43]), read_string(&combined));
}

// FragmentSyncTest.testTruncatedFragmentIndelMismatch
#[test]
fn test_truncated_fragment_indel_mismatch() {
    let r = ref_bases();

    // 3' read has a delete instead of a soft-clip at the start and end
    let mut first = create_sam_record(READ_ID, 4, &format!("{}{}", &r[4..8], &r[10..40]), "4M2D30M");

    first.insert_size = 30;
    set_negative_strand(&mut first, true);

    let mut second = create_sam_record(READ_ID, 10, &format!("{}{}", &r[10..40], &r[42..46]), "30M2D2M2S");
    second.insert_size = -30;

    // SyncType == COMBINED is the Combined arm
    let combined = combined_record(&first, &second);
    assert_eq!(30, combined.insert_size);
    assert_eq!(10, combined.start);
    assert_eq!(39, combined.end);
    assert_eq!("30M", cigar_str(&combined.cigar));
    assert_eq!(30, combined.quals.len());

    assert_eq!(&r[10..40], read_string(&combined));
}

/// Not in the Java suite: the fast path for two plain-match reads gives exactly what the general
/// base-by-base walk gives, for every arrangement of the two reads.
#[test]
fn test_plain_match_fast_path_equals_general_path() {
    use crate::sync::form_fragment_read_with;
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut next = |bound: u64| -> u64 {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) % bound
    };
    let mut combined = 0;
    for _ in 0..20000 {
        let mut make = |start: i32, next: &mut dyn FnMut(u64) -> u64| -> Read {
            let len = 1 + next(40) as usize;
            let bases: String = (0..len).map(|_| ['A', 'C', 'G', 'T'][next(4) as usize]).collect();
            // a few quality values, so that equal and differing quals both occur
            let quals: Vec<u8> = (0..len).map(|_| [1u8, 11, 25, 37, 12, 30][next(6) as usize]).collect();
            let mut read = build_sam_record_quals(start, &format!("{len}M"), &bases, quals);
            read.name = READ_ID.to_string();
            read
        };
        let mut first = make(100 + next(30) as i32, &mut next);
        let mut second = make(100 + next(30) as i32, &mut next);
        // mostly matching bases over the overlap, as real mates have
        if next(4) != 0 {
            for pos in first.start.max(second.start)..=first.end.min(second.end) {
                if next(10) != 0 {
                    second.bases[(pos - second.start) as usize] = first.bases[(pos - first.start) as usize];
                }
            }
        }
        set_negative_strand(&mut first, next(2) == 0);
        set_negative_strand(&mut second, next(8) != 0 && !first.negative_strand());
        let fast = form_fragment_read_with(&first, &second, crate::bamio::blank_read(), true);
        let general = form_fragment_read_with(&first, &second, crate::bamio::blank_read(), false);
        match (fast, general) {
            (Outcome::Combined(a), Outcome::Combined(b)) => {
                combined += 1;
                assert_eq!(a.cigar, b.cigar);
                assert_eq!(a.bases, b.bases);
                assert_eq!(a.quals, b.quals);
                assert_eq!((a.start, a.end, a.flags, a.mapq), (b.start, b.end, b.flags, b.mapq));
            }
            (Outcome::NoOverlap, Outcome::NoOverlap) | (Outcome::Inversion, Outcome::Inversion) => {}
            (Outcome::CigarMismatch, Outcome::CigarMismatch) | (Outcome::Exception, Outcome::Exception) => {}
            _ => panic!("fast and general paths disagree on the outcome"),
        }
    }
    assert!(combined > 2000);
}
