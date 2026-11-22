//! ReadAdjustmentsTest
use super::*;
use crate::assembly::constants::{BASE_QUAL_MINIMUM, LINE_BASE_A, LINE_BASE_T, POLY_G_TRIM_LENGTH};
use crate::assembly::line_utils::{find_ref_base_repeat_count, has_line_tail};
use crate::assembly::read::read_adjustments::{find_low_base_qual_trim_count_ext, mark_line_soft_clips, trim_low_qual_soft_clip_bases, trim_poly_g_sequences};
use crate::assembly::remote_region_finder::filter_low_qual_read;
use crate::assembly::seq_tech_utils::trim_illumina_adapter_bases;

const LOW_QUAL_BASE: u8 = 10;

fn set_quals(read: &ReadP, base_qualities: &[u8]) {
    read.m().bam_record_mut().base_qualities = base_qualities.to_vec();
}

fn set_reversed(read: &ReadP) {
    read.m().bam_record_mut().set_flag(FLAG_READ_REVERSE_STRAND, true);
}

#[test]
fn test_poly_g_trimming() {
    let _env = test_env();
    let poly_g_section = "GGGGGGGG";
    let poly_c_section = "CCCCCCCC";
    let other_bases = "AAAACCCCGAGATTTT";
    // wrong end
    let read_bases = format!("{poly_g_section}{other_bases}");
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, 0, 0));
    assert!(!trim_poly_g_sequences(read.m()));
    // too few Gs
    let read_bases = format!("{}{}", other_bases, &poly_g_section[0..3]);
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, 0, 0));
    assert!(!trim_poly_g_sequences(read.m()));
    // trimmed on the 3' end
    let read_bases = format!("{}{}", other_bases, &poly_g_section[0..POLY_G_TRIM_LENGTH as usize]);
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, 0, 0));
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(other_bases.len() as i32, read.bases_length());
    assert_eq!(read.alignment_start() + other_bases.len() as i32 - 1, read.alignment_end());
    assert_eq!(other_bases, read.get_bases_string());
    assert_eq!(make_cigar_string(other_bases, 0, 0), read.cigar_string());
    // again with a soft-clipping at the end
    let read_bases = format!("{}{}", other_bases, &poly_g_section[0..8]);
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, 0, 4));
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(other_bases.len() as i32, read.bases_length());
    assert_eq!(read.alignment_start() + other_bases.len() as i32 - 1, read.alignment_end());
    assert_eq!(other_bases, read.get_bases_string());
    assert_eq!(make_cigar_string(other_bases, 0, 0), read.cigar_string());
    // negative strand - trimmed at the start
    let read_bases = format!("{}{}", &poly_c_section[0..6], other_bases);
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, 3, 0));
    set_reversed(&read);
    assert_eq!(97, read.unclipped_start());
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(103, read.alignment_start());
    assert_eq!(103, read.unclipped_start());
    assert_eq!(103 + other_bases.len() as i32 - 1, read.alignment_end());
    assert_eq!(other_bases, read.get_bases_string());
    assert_eq!(make_cigar_string(other_bases, 0, 0), read.cigar_string());
    // trimming into delete
    let read_bases = format!("{}{}", other_bases, &poly_g_section[0..5]);
    let read = create_read(TEST_READ_ID, 100, &read_bases, "16M2D5M");
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(115, read.alignment_end());
    assert_eq!("16M", read.cigar_string());
    let read_bases = format!("{}{}", other_bases, &poly_g_section[0..6]);
    let read = create_read(TEST_READ_ID, 100, &read_bases, "16M2D5M");
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(114, read.alignment_end());
    assert_eq!("15M", read.cigar_string());
    // again from the other end
    let read_bases = format!("{}{}", &poly_c_section[0..5], other_bases);
    let read = create_read(TEST_READ_ID, 100, &read_bases, "5M2D16M");
    set_reversed(&read);
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(107, read.alignment_start());
    assert_eq!("16M", read.cigar_string());
    let read_bases = format!("{}{}", &poly_c_section[0..6], other_bases);
    let read = create_read(TEST_READ_ID, 100, &read_bases, "5M2D16M");
    set_reversed(&read);
    assert!(trim_poly_g_sequences(read.m()));
    assert_eq!(108, read.alignment_start());
    assert_eq!("15M", read.cigar_string());
}

#[test]
fn test_low_base_filtering() {
    let _env = test_env();
    let read_bases = ref_bases_random_100();
    let mut base_qualities = build_default_base_quals(read_bases.len());
    for q in base_qualities.iter_mut().take(50) {
        *q = LOW_QUAL_BASE;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    set_quals(&read, &base_qualities);
    assert!(!filter_low_qual_read(read.m().bam_record_mut()));
    base_qualities[50] = LOW_QUAL_BASE;
    set_quals(&read, &base_qualities);
    assert!(filter_low_qual_read(read.m().bam_record_mut()));
}

#[test]
fn test_soft_clip_low_base_qual_trimming() {
    let _env = test_env();
    let read_bases = ref_bases_random_100();
    let mut base_qualities = build_default_base_quals(read_bases.len());
    let trim_length = 8;
    let n = base_qualities.len();
    for i in 0..trim_length {
        base_qualities[i] = LOW_QUAL_BASE;
        base_qualities[n - i - 1] = LOW_QUAL_BASE;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    set_quals(&read, &base_qualities);
    assert!(!trim_low_qual_soft_clip_bases(read.m())); // nothing without soft-clips

    let read = create_read(TEST_READ_ID, 110, &read_bases, &make_cigar_string(&read_bases, 10, 10));
    set_quals(&read, &base_qualities);
    assert!(trim_low_qual_soft_clip_bases(read.m()));
    assert_eq!(110, read.alignment_start());
    assert_eq!(189, read.alignment_end());
    assert_eq!(191, read.unclipped_end());
    assert_eq!("10S80M2S", read.cigar_string());

    let read = create_read(TEST_READ_ID, 110, &read_bases, &make_cigar_string(&read_bases, 10, 10));
    for q in base_qualities.iter_mut().skip(90) {
        *q = DEFAULT_BASE_QUAL;
    }
    base_qualities[95] = LOW_QUAL_BASE;
    base_qualities[98] = LOW_QUAL_BASE;
    set_quals(&read, &base_qualities);
    assert!(trim_low_qual_soft_clip_bases(read.m()));
    assert_eq!(110, read.alignment_start());
    assert_eq!(189, read.alignment_end());
    assert_eq!(197, read.unclipped_end());
    assert_eq!("10S80M8S", read.cigar_string());

    let read = create_read(TEST_READ_ID, 110, &read_bases, &make_cigar_string(&read_bases, 10, 10));
    set_reversed(&read);
    set_quals(&read, &base_qualities);
    assert!(trim_low_qual_soft_clip_bases(read.m()));
    assert_eq!(110, read.alignment_start());
    assert_eq!(108, read.unclipped_start());
    assert_eq!("2S80M10S", read.cigar_string());
}

#[test]
fn test_low_base_qual_trimming() {
    let _env = test_env();
    let read_bases = ref_bases_random_100();
    let mut base_qualities = build_default_base_quals(read_bases.len());
    // low qual will be 70% of outer bases
    let n = base_qualities.len();
    for i in 3..10 {
        base_qualities[i] = LOW_QUAL_BASE;
        base_qualities[n - i - 1] = LOW_QUAL_BASE;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    set_quals(&read, &base_qualities);
    Read::trim_low_qual_bases(&read);
    assert!(read.low_qual_trimmed());
    assert_eq!("90M", read.cigar_string());
    assert_eq!(10, read.trim_count());
    assert_eq!(189, read.alignment_end());

    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    set_quals(&read, &base_qualities);
    set_reversed(&read);
    Read::trim_low_qual_bases(&read);
    assert!(read.low_qual_trimmed());
    assert_eq!("90M", read.cigar_string());
    assert_eq!(10, read.trim_count());
    assert_eq!(110, read.alignment_start());
}

#[test]
fn test_uncertain_base_qual_trimming() {
    let _env = test_env();
    let read_bases = ref_bases_random_100();
    let mut base_qualities = build_default_base_quals(read_bases.len());
    let n = base_qualities.len();
    for i in 3..10 {
        base_qualities[i] = BASE_QUAL_MINIMUM;
        base_qualities[n - i - 1] = BASE_QUAL_MINIMUM;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, TEST_CIGAR_100);
    set_quals(&read, &base_qualities);
    assert_eq!(10, find_low_base_qual_trim_count_ext(&read, 0, read.bases_length() - 1, true, true));
    assert_eq!(10, find_low_base_qual_trim_count_ext(&read, 0, read.bases_length() - 1, false, true));
}

#[test]
fn test_illumina_adapter() {
    let _env = test_env();
    let read_bases = ref_bases_random_100();
    let read = create_read(TEST_READ_ID, 110, &read_bases, "20S60M20S");
    read.m().bam_record_mut().set_flag(FLAG_MATE_REVERSE_STRAND, true);
    let mate_read = create_read(TEST_READ_ID, 100, &read_bases, "20S60M20S");
    set_reversed(&mate_read);
    read.m().set_mate_read(&mate_read);
    mate_read.m().set_mate_read(&read);
    trim_illumina_adapter_bases(&read);
    trim_illumina_adapter_bases(&mate_read);
    assert_eq!(10, read.trim_count());
    assert_eq!("20S60M10S", read.cigar_string());
    assert_eq!(179, read.unclipped_end());
    assert_eq!(10, mate_read.trim_count());
    assert_eq!("10S60M20S", mate_read.cigar_string());
    assert_eq!(90, mate_read.unclipped_start());
    // the reads refer to each other: break the cycle so they are released
    read.m().clear_mate_read();
    mate_read.m().clear_mate_read();
}

#[test]
fn test_mark_line_tail_read() {
    let _env = test_env();
    let random = ref_bases_random_100();
    let non_line_sequence = "GCTGCTGTCGTGTCC";
    let left_tail = |line_sequence: &str, local_repeat: &str, default_quals: &[u8]| -> bool {
        let soft_clip_bases = format!("{non_line_sequence}{line_sequence}");
        let read_bases = format!("{soft_clip_bases}{local_repeat}{random}");
        has_line_tail(read_bases.as_bytes(), soft_clip_bases.len() as i32 - 1, true, LINE_BASE_A, Some(default_quals))
    };
    // as in Java, one quality array sized for the first read is used throughout
    let default_quals = build_default_base_quals(non_line_sequence.len() + 16 + random.len());
    assert!(left_tail("AAAAAAAAAAAAAAAA", "", &default_quals));
    assert!(left_tail("AAAAAAGAAAAAATAAAA", "", &default_quals));
    // too many non-line bases
    assert!(!left_tail("AAAAAAGAAACAATAAAA", "", &default_quals));
    // an extension of a local repeat
    assert!(!left_tail("AAAAAAAAAAAAAAAA", "AAAAAAAAAA", &default_quals));
    // reverse side
    let soft_clip_bases = format!("TTTTTCTTTTTTTGTTTT{non_line_sequence}");
    let read_bases = format!("{random}{soft_clip_bases}");
    assert!(has_line_tail(read_bases.as_bytes(), random.len() as i32, false, LINE_BASE_T, Some(&default_quals)));
    let read_bases = format!("{random}TTTTTTTTTT{soft_clip_bases}");
    assert!(!has_line_tail(read_bases.as_bytes(), random.len() as i32 + 10, false, LINE_BASE_T, Some(&default_quals)));

    // ref base repeats
    for (ref_bases, expected) in
        [("AAAAAAAAAA", 10), ("AAAAAGAAAA", 9), ("AAAAAGAAGA", 5), ("AAAAAAAAAAGAAAA", 14), ("AAAAAAAAAAGAAAAAAAAAAG", 20), ("AAAAAAAAAAGAAAAG", 10), ("AAAAAACGT", 6)]
    {
        assert_eq!(expected, find_ref_base_repeat_count(ref_bases.as_bytes(), 0, true, LINE_BASE_A), "{ref_bases}");
    }
}

#[test]
fn test_line_low_qual_trimming() {
    let _env = test_env();
    let random = ref_bases_random_100();
    let line_sequence = "AAAAAAAAAAAAAAAA";
    let soft_clip_bases = format!("GCTGCTGTCGTGTCC{line_sequence}");
    let read_bases = format!("{soft_clip_bases}{random}");
    let mut base_qualities = build_default_base_quals(read_bases.len());
    let low_qual_length = read_bases.len() - 2;
    for q in base_qualities.iter_mut().take(low_qual_length) {
        *q = LOW_QUAL_BASE;
    }
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, soft_clip_bases.len(), 0));
    set_quals(&read, &base_qualities);
    set_reversed(&read);
    mark_line_soft_clips(read.m());
    assert!(trim_low_qual_soft_clip_bases(read.m()));
    assert_eq!(84, read.unclipped_start());
    assert_eq!("16S100M", read.cigar_string());

    // keep 2 extra bases since the non-line bases are further in
    let soft_clip_bases = format!("GCTGCTGTCGTGT{line_sequence}CC");
    let read_bases = format!("{soft_clip_bases}{random}");
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, soft_clip_bases.len(), 0));
    set_quals(&read, &base_qualities);
    set_reversed(&read);
    mark_line_soft_clips(read.m());
    assert!(trim_low_qual_soft_clip_bases(read.m()));
    assert_eq!(82, read.unclipped_start());
    assert_eq!("18S100M", read.cigar_string());

    // no LINE trimming if an extension of a matching repeat (as in Java, the qualities are those of the shorter read)
    let read_bases = format!("{soft_clip_bases}{line_sequence}{random}");
    let read = create_read(TEST_READ_ID, 100, &read_bases, &make_cigar_string(&read_bases, soft_clip_bases.len(), 0));
    set_quals(&read, &base_qualities);
    set_reversed(&read);
    mark_line_soft_clips(read.m());
    assert!(!read.has_line_tail_start());
    assert!(!read.has_line_tail_end());

    // the other orientation
    let line_sequence = reverse_complement_str(line_sequence);
    let soft_clip_bases = format!("CC{line_sequence}GCTGCTGTCGTGTCC");
    let read_bases = format!("{random}{soft_clip_bases}");
    let base_qualities = vec![LOW_QUAL_BASE; read_bases.len()];
    let read = create_read(TEST_READ_ID, 101, &read_bases, &make_cigar_string(&read_bases, 0, soft_clip_bases.len()));
    set_quals(&read, &base_qualities);
    mark_line_soft_clips(read.m());
    assert!(trim_low_qual_soft_clip_bases(read.m()));
    assert_eq!(200, read.alignment_end());
    assert_eq!(218, read.unclipped_end());
    assert_eq!("100M18S", read.cigar_string());

    // expansion of local repeat
    let read_bases = format!("{random}{line_sequence}{soft_clip_bases}");
    let read = create_read(TEST_READ_ID, 101, &read_bases, &make_cigar_string(&read_bases, 0, soft_clip_bases.len()));
    set_quals(&read, &build_default_base_quals(read_bases.len()));
    mark_line_soft_clips(read.m());
    assert!(!read.has_line_tail_start());
    assert!(!read.has_line_tail_end());
}

#[test]
fn test_low_qual_line_tail_read() {
    let _env = test_env();
    // uncertain qual bases can still rule out a LINE sequence
    let line_sequence = "AAAACAAAACAAAACAAAAC";
    let non_line_sequence = "GCTGCTGTCGTGTCC";
    let soft_clip_bases = format!("{non_line_sequence}{line_sequence}");
    let read_bases = format!("{}{}", soft_clip_bases, ref_bases_random_100());
    let mut default_quals = build_default_base_quals(read_bases.len());
    default_quals[4] = 1;
    default_quals[9] = 1;
    default_quals[14] = 1;
    assert!(!has_line_tail(read_bases.as_bytes(), soft_clip_bases.len() as i32 - 1, true, LINE_BASE_A, Some(&default_quals)));
}
