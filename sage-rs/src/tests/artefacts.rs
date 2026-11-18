//! Ported from the Java SAGE test suite: seqtech/ArtefactsTest.java (the Illumina artefact context)

use super::util::*;
use crate::read_context::{has_indel_in_core, ArtefactContext, VariantReadContext, ARTEFACT_NO_INDEX};
use crate::types::*;

const FLANK_BASES: &str = "AAAAAGGGGG";
const HOMOPOLYMER_SEQ: &str = "TTTTTTTT";

//            10         20         30
// 0123456789 01234567 89012345 67890123
// AAAAAGGGGG ACGTTGCA TTTTTTTT ACGTTGCA AAAAAGGGGG
fn ref_bases() -> String {
    format!("{FLANK_BASES}ACGTTGCA{HOMOPOLYMER_SEQ}ACGTTGCA{FLANK_BASES}")
}

const INVALID_BASE_QUAL: i32 = -1;

/// VariantUtils.createReadContext(variant, readBases, leftCoreIndex, varIndex, rightCoreIndex): a context
/// with no repeat or homology
fn create_read_context(
    variant: &SimpleVariant, read_bases: &str, left_core_index: usize, var_index: usize, right_core_index: usize,
) -> VariantReadContext {
    let left_flank = &read_bases[0..left_core_index];
    let left_core = &read_bases[left_core_index..var_index];
    let right_core = &read_bases[var_index + 1..right_core_index + 1];
    let right_flank = &read_bases[right_core_index + 1..];

    let core_index_start = left_flank.len() as i32;
    let var_read_index = core_index_start + left_core.len() as i32;
    let core_index_end = var_read_index + variant.alt_len() - 1 + right_core.len() as i32;
    let ref_bases = format!("{}{}{}", left_core, variant.ref_str(), right_core);
    let read_bases = format!("{}{}{}{}{}", left_flank, left_core, variant.alt_str(), right_core, right_flank);

    let alignment_start = variant.pos - var_read_index;
    let core_right_pos_start = variant.pos + variant.ref_len().min(variant.alt_len());
    let alignment_end = core_right_pos_start + right_core.len() as i32 + right_flank.len() as i32 - 1;
    let read_cigar = vec![CigarEl::new(read_bases.len() as i32, Op::M)];

    VariantReadContext {
        variant: variant.clone(),
        alignment_start,
        alignment_end,
        ref_bases: ref_bases.into_bytes(),
        read_bases: read_bases.into_bytes(),
        core_index_start,
        var_index: var_read_index,
        core_index_end,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start: variant.pos - left_core.len() as i32,
        core_position_end: variant.pos + right_core.len() as i32,
        read_cigar: cigar_str(&read_cigar),
        has_indel_in_core: has_indel_in_core(&read_cigar, core_index_start, core_index_end),
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    }
}

/// IlluminaArtefactContext.buildContext, asserted not null
fn build_artefact_context(read_context: &VariantReadContext) -> ArtefactContext {
    ArtefactContext::build(read_context).expect("no exception").expect("artefact context")
}

fn variant(position: i32, ref_: &str, alt: &str) -> SimpleVariant {
    SimpleVariant::new(position, ref_.as_bytes(), alt.as_bytes())
}

/// buildSamRecord(alignmentStart, cigar, readBases, qualities) on the negative strand
fn negative_strand_record(alignment_start: i32, cigar: &str, read_bases: &str, quals: Vec<u8>) -> Read {
    let mut record = build_sam_record_quals(alignment_start, cigar, read_bases, quals);
    record.flags |= FLAG_REVERSE;
    record
}

// ArtefactsTest.testHomopolymerArtefacts
#[test]
fn test_homopolymer_artefacts() {
    let ref_bases = ref_bases();
    let pos = 100;

    // A>T where SNV is immediately next to (left of) upstream homopolymer (HP)
    let mut var = variant(pos, "A", "T");

    let mut var_index = 17usize;
    let mut read_context_bases = format!("{}T{}", &ref_bases[0..var_index], &ref_bases[var_index + 1..]);

    let mut read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    let mut artefact_context = build_artefact_context(&read_context);
    assert!(artefact_context.offset[SE_START] == ARTEFACT_NO_INDEX);
    assert_eq!(0, artefact_context.offset[SE_END]);
    assert_eq!(b'T' as i8, artefact_context.base[SE_END]);

    let read_length = read_context_bases.len();
    let mut cigar = format!("{read_length}M");
    let mut read_qualities = default_quals(read_length);
    let low_qual_base = 10u8;
    read_qualities[var_index] = low_qual_base;

    let mut record = build_sam_record_quals(90, &cigar, &read_context_bases, read_qualities);

    let mut adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(adjusted_base_qual, INVALID_BASE_QUAL);

    record.flags |= FLAG_REVERSE;

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);

    // single base DEL on left of HP
    var = variant(pos, "CA", "C");
    var_index = 16;
    read_context_bases = format!("{}{}", &ref_bases[0..var_index + 1], &ref_bases[var_index + 2..]);

    read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    artefact_context = build_artefact_context(&read_context);
    assert_eq!(1, artefact_context.offset[SE_END]);
    assert_eq!(b'T' as i8, artefact_context.base[SE_END]);

    read_qualities = default_quals(read_length);
    read_qualities[var_index + 1] = low_qual_base;

    record = negative_strand_record(20, &cigar, &read_context_bases, read_qualities);

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);

    // similar before but where the variant is left aligned so not immediately next to the HP
    var = variant(pos, "TG", "T");
    var_index = 15;

    //             0123456 7890
    // ACGTTTGG -> ACGTTTG TTTT

    //            10             20         30
    // 0123456789 012345     678901234 567890123
    // AAAAAGGGGG ACGTTT [G] GTTTTTTTT ACGTTGCA AAAAAGGGGG
    // low-qual               X
    read_context_bases = format!("{FLANK_BASES}ACGTTTG{}", &ref_bases[18..]);

    read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    artefact_context = build_artefact_context(&read_context);
    assert_eq!(2, artefact_context.offset[SE_END]);
    assert_eq!(b'T' as i8, artefact_context.base[SE_END]);

    read_qualities = default_quals(read_length);
    read_qualities[17] = low_qual_base;

    record = negative_strand_record(20, &cigar, &read_context_bases, read_qualities);

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);

    // similar to above but the delete is 2-bases left-aligned, so not immediately left of the HP
    //            10              20         30
    // 0123456789 01234567      89012345 67890123
    // AAAAAGGGGG AACGTAGT [GT] TTTTTTTT ACGTTGCA AAAAAGGGGG
    // low-qual                 X
    var = variant(pos, "AGT", "A");
    var_index = 15;

    read_context_bases = format!("{FLANK_BASES}AACGTAGT{HOMOPOLYMER_SEQ}{}", &ref_bases[26..]);

    read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    artefact_context = build_artefact_context(&read_context);
    assert_eq!(3, artefact_context.offset[SE_END]);
    assert_eq!(b'T' as i8, artefact_context.base[SE_END]);

    read_qualities = default_quals(read_length);
    read_qualities[18] = low_qual_base;

    record = negative_strand_record(20, &cigar, &read_context_bases, read_qualities);

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);

    // with germline insertion which needs to be ignored when looking for the ref base

    // AAAAAGGGGG ACGTTGC insT [T] TTTTTTTT ACGTTGCA AAAAAGGGGG   - read
    // AAAAAGGGGG ACGTTGC      [A] TTTTTTTT ACGTTGCA AAAAAGGGGG   - ref

    // an SNV which follows a germline T insertion, which should be ignored when finding the first ref base
    //            10                20         30
    // 0123456789 0123456    7  8  9  01234567 890123
    // AAAAAGGGGG AACGTAG insT [T] TTTTTTTT ACGTTGCA AAAAAGGGGG
    // low-qual                 X
    var = variant(pos, "A", "T");
    var_index = 18;

    read_context_bases = format!("{FLANK_BASES}AACGTAGTT{HOMOPOLYMER_SEQ}{}", &ref_bases[26..]);

    read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    artefact_context = build_artefact_context(&read_context);
    assert_eq!(-1, artefact_context.offset[SE_END]);
    assert_eq!(b'T' as i8, artefact_context.base[SE_END]);

    read_qualities = default_quals(read_length);
    read_qualities[18] = low_qual_base;

    cigar = "17M1I27M".to_string();
    record = negative_strand_record(20, &cigar, &read_context_bases, read_qualities);

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);
}

// ArtefactsTest.testHomopolymerArtefactsPosStrand
#[test]
fn test_homopolymer_artefacts_pos_strand() {
    let ref_bases = ref_bases();
    let pos = 100;

    let mut var = variant(pos, "A", "T");

    let mut var_index = 26usize;
    let mut read_context_bases = format!("{}T{}", &ref_bases[0..var_index], &ref_bases[var_index + 1..]);

    let mut read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    let mut artefact_context = build_artefact_context(&read_context);
    assert!(artefact_context.offset[SE_END] == ARTEFACT_NO_INDEX);
    assert_eq!(0, artefact_context.offset[SE_START]);

    let read_length = read_context_bases.len();
    let cigar = format!("{read_length}M");
    let mut read_qualities = default_quals(read_length);
    let low_qual_base = 10u8;
    read_qualities[var_index] = low_qual_base;

    let mut record = build_sam_record_quals(20, &cigar, &read_context_bases, read_qualities);

    let mut adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);

    // ignore -ve strand for downstream HPs
    record.flags |= FLAG_REVERSE;

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(INVALID_BASE_QUAL, adjusted_base_qual);

    // single base DEL on left of HP
    var = variant(pos, "CA", "C");
    var_index = 16;
    read_context_bases = format!("{}{}", &ref_bases[0..var_index + 1], &ref_bases[var_index + 2..]);

    read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    artefact_context = build_artefact_context(&read_context);
    assert_eq!(1, artefact_context.offset[SE_END]);

    read_qualities = default_quals(read_length);
    read_qualities[var_index + 1] = low_qual_base;

    record = negative_strand_record(20, &cigar, &read_context_bases, read_qualities);

    adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);
}

// ArtefactsTest.testSnvWithAdjacentInsert
#[test]
fn test_snv_with_adjacent_insert() {
    let pos = 115;
    let left_flank = "AAAATCAGGA";
    let right_flank = "AAAGGCAGAA";

    let var = variant(pos, "G", "T");

    let var_index = 15usize;
    let read_context_bases = format!("{left_flank}CTGTTTTTTTTTTTTTTTTTTA{right_flank}");

    let read_context = create_read_context(&var, &read_context_bases, var_index - 2, var_index, var_index + 10);

    let artefact_context = build_artefact_context(&read_context);

    let read_length = read_context_bases.len();
    let cigar = format!("{read_length}M");
    let mut read_qualities = default_quals(read_length);
    let low_qual_base = 10u8;
    read_qualities[var_index - 1] = low_qual_base;

    let record = negative_strand_record(20, &cigar, &read_context_bases, read_qualities);

    let adjusted_base_qual = artefact_context.find_applicable_base_qual(&record, var_index as i32);
    assert_eq!(low_qual_base as i32, adjusted_base_qual);
}
