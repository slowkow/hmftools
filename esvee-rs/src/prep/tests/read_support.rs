//! ReadSupportTest

use super::util::*;
use crate::prep::filters::*;
use crate::prep::tracker::{has_exact_junction_support, has_other_junction_support};
use crate::prep::types::Orientation::{Forward, Reverse};
use crate::prep::types::*;

fn ref_bases() -> String {
    let unit = "AAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCC";
    format!("{}{}", unit.repeat(8), generate_random_bases(200))
}

fn with_base(bases: &str, changes: &[(usize, char)]) -> String {
    let mut b: Vec<u8> = bases.as_bytes().to_vec();
    for &(i, c) in changes {
        b[i] = c as u8;
    }
    String::from_utf8(b).unwrap()
}

#[test]
fn test_exact_supporting_reads() {
    let ref_bases = ref_bases();
    let junction_read = prep_read(sam_record(&next_id(), CHR_1, 230, &ref_bases[..100], "30S70M"));
    let junction_data = JunctionData::new(230, Reverse, Some(0));
    let supports = |read: &PrepRead, junction: &JunctionData, junc_read: &PrepRead| has_exact_junction_support(read, junction, Some(junc_read));

    // exact position match
    let support_read = prep_read(sam_record(&next_id(), CHR_1, 230, &ref_bases[29..100], "1S70M"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    // offset to an earlier base
    let support_read = prep_read(sam_record(&next_id(), CHR_1, 225, &ref_bases[22..100], "3S72M"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    // offset to an later base
    let id = next_id();
    let bases = ref_bases[24..100].to_string();
    let support_read = prep_read(sam_record(&id, CHR_1, 234, &bases, "10S66M"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    // testing mismatches - 1 is allowed
    let bases = with_base(&bases, &[(0, 'A')]);
    let support_read = prep_read(sam_record(&id, CHR_1, 234, &bases, "10S66M"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    let bases = with_base(&bases, &[(1, 'A')]);
    let support_read = prep_read(sam_record(&id, CHR_1, 234, &bases, "10S66M"));
    assert!(!supports(&support_read, &junction_data, &junction_read));

    let mut rec = sam_record(&id, CHR_1, 234, &bases, "10S66M");
    let mut quals = rec.qual().to_vec();
    quals[0] = 10;
    let cigar = rec.cigar().take();
    rec.set(id.as_bytes(), Some(&cigar), bases.as_bytes(), &quals);
    let support_read = prep_read(rec);
    assert!(supports(&support_read, &junction_data, &junction_read));

    // for positive orientation
    let junction_read = prep_read(sam_record("01", CHR_1, 200, &ref_bases[..80], "50M30S"));
    let junction_data = JunctionData::new(249, Forward, Some(0));

    // exact position match
    let support_read = prep_read(sam_record("02", CHR_1, 200, &ref_bases[..53], "50M3S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    let support_read = prep_read(sam_record("02", CHR_1, 200, &ref_bases[..79], "49M30S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    let support_read = prep_read(sam_record("02", CHR_1, 200, &ref_bases[..71], "51M20S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    // offset to an earlier base
    let support_read = prep_read(sam_record("02", CHR_1, 200, &ref_bases[..60], "45M15S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    let support_read = prep_read(sam_record("02", CHR_1, 210, &ref_bases[10..60], "35M15S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    // offset to an later base
    let bases = ref_bases[..65].to_string();
    let support_read = prep_read(sam_record("02", CHR_1, 200, &bases, "55M10S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    // testing mismatches
    let bases = with_base(&bases, &[(58, 'T'), (62, 'T')]);
    let support_read = prep_read(sam_record("02", CHR_1, 200, &bases, "55M10S"));
    assert!(supports(&support_read, &junction_data, &junction_read));

    let bases = with_base(&bases, &[(61, 'T')]);
    let support_read = prep_read(sam_record("02", CHR_1, 200, &bases, "55M10S"));
    assert!(!supports(&support_read, &junction_data, &junction_read));

    let mut rec = sam_record("02", CHR_1, 200, &bases, "55M10S");
    let mut quals = rec.qual().to_vec();
    quals[58] = 10;
    quals[62] = 10;
    let cigar = rec.cigar().take();
    rec.set(b"02", Some(&cigar), bases.as_bytes(), &quals);
    let support_read = prep_read(rec);
    assert!(supports(&support_read, &junction_data, &junction_read));
}

#[test]
fn test_read_repeat_trimming() {
    let ref_bases = ref_bases();
    let read_bases = format!("ACGTACGTAA{}{}GGGGGGGGGG{}", &ref_bases[10..40], &ref_bases[50..70], &ref_bases[70..100]);
    let read = prep_read(sam_record(&next_id(), CHR_1, 10, &read_bases, "10S30M10D20M10I30M"));
    assert!(above_repeat_trimmed_alignment_threshold(&read, MIN_ALIGNMENT_BASES, true));

    let (repeat1, repeat2) = ("ACGTT", "AAGG");
    let read_bases = format!("CGT{}TGAT{}AAT", repeat1.repeat(6), repeat2.repeat(8));
    let cigar = format!("{}M", read_bases.len());
    let read = prep_read(sam_record(&next_id(), CHR_1, 10, &read_bases, &cigar));
    assert!(!above_repeat_trimmed_alignment_threshold(&read, MIN_ALIGNMENT_BASES, true));
}

fn create_read(read_id: &str, read_start: i32, first_in_pair: bool, reversed: bool) -> PrepRead {
    prep_read(sam_record_flags(read_id, CHR_1, read_start, "", "100M", build_flags(first_in_pair, reversed, false), DEFAULT_BASE_QUAL))
}

#[test]
fn test_distant_supporting_reads() {
    let ref_bases = ref_bases();
    let filters = ReadFilterConfig::default();

    // first negative orientation
    let _junction_read = prep_read(sam_record(&next_id(), CHR_1, 230, &ref_bases[..100], "30S70M"));
    let junction_data = JunctionData::new(230, Reverse, Some(0));

    // distant within range, correct orientation
    let support_read = create_read(&next_id(), 500, false, true);
    assert!(has_other_junction_support(&support_read, &junction_data, &filters));

    // too far
    let support_read = create_read(&next_id(), 1500, false, true);
    assert!(!has_other_junction_support(&support_read, &junction_data, &filters));

    // wrong orientation
    let support_read = create_read(&next_id(), 500, false, false);
    assert!(!has_other_junction_support(&support_read, &junction_data, &filters));

    // wrong side of junction
    let support_read = create_read(&next_id(), 220, false, true);
    assert!(!has_other_junction_support(&support_read, &junction_data, &filters));

    // other junction side
    let _junction_read = prep_read(sam_record(&next_id(), CHR_1, 1500, &ref_bases[..100], "70M30S"));
    let junction_data = JunctionData::new(1500, Forward, Some(0));

    let support_read = prep_read(sam_record_mate(&next_id(), CHR_1, 1000, CHR_3, 100, true, false, None));
    assert!(has_other_junction_support(&support_read, &junction_data, &filters));

    // too far
    let support_read = create_read(&next_id(), 100, true, false);
    assert!(!has_other_junction_support(&support_read, &junction_data, &filters));

    // wrong orientation
    let support_read = create_read(&next_id(), 100, true, true);
    assert!(!has_other_junction_support(&support_read, &junction_data, &filters));

    // wrong side of junction
    let support_read = create_read(&next_id(), 1510, true, false);
    assert!(!has_other_junction_support(&support_read, &junction_data, &filters));
}

#[test]
fn test_read_filter_flags() {
    let is_set = |filters: i32, f: i32| filters & f != 0;
    let mut filters = 0;
    filters |= F_SOFT_CLIP_LENGTH;
    assert!(is_set(filters, F_SOFT_CLIP_LENGTH));
    assert!(!is_set(filters, F_INSERT_MAP_OVERLAP));

    filters |= F_INSERT_MAP_OVERLAP;
    assert!(is_set(filters, F_SOFT_CLIP_LENGTH));
    assert!(is_set(filters, F_INSERT_MAP_OVERLAP));

    filters &= !F_SOFT_CLIP_LENGTH;
    assert!(!is_set(filters, F_SOFT_CLIP_LENGTH));
    assert!(is_set(filters, F_INSERT_MAP_OVERLAP));

    filters &= !F_INSERT_MAP_OVERLAP;
    assert!(!is_set(filters, F_SOFT_CLIP_LENGTH));
    assert!(!is_set(filters, F_INSERT_MAP_OVERLAP));
}

#[test]
fn test_repetitive_breaks() {
    let bases = generate_random_bases(30);
    assert!(!is_repetitive_section_break(bases.as_bytes(), true, 10));
    assert!(!is_repetitive_section_break(bases.as_bytes(), false, 10));

    let bases = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    assert!(is_repetitive_section_break(bases.as_bytes(), true, 10));
    assert!(is_repetitive_section_break(bases.as_bytes(), false, 10));

    // 2-base repeats
    let bases = "ATATATATATATATATATATATATATATAT";
    assert!(is_repetitive_section_break(bases.as_bytes(), true, 10));
    assert!(is_repetitive_section_break(bases.as_bytes(), false, 10));

    // with an error
    let bases = "ATATATAGATATATATATATAGATATATAT";
    assert!(!is_repetitive_section_break(bases.as_bytes(), true, 10));
    assert!(!is_repetitive_section_break(bases.as_bytes(), false, 10));

    // 3-base repeats
    let bases = "ATCATCATCATCATCATCATCATCATCATCATC";
    assert!(is_repetitive_section_break(bases.as_bytes(), true, 10));
    assert!(is_repetitive_section_break(bases.as_bytes(), false, 10));

    // with an error
    let bases = "ATCATCATGATCATCATCATCATCGTCATCATC";
    assert!(!is_repetitive_section_break(bases.as_bytes(), true, 10));
    assert!(!is_repetitive_section_break(bases.as_bytes(), false, 10));
}
