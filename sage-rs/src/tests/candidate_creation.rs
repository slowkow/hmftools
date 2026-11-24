//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/candidate/CandidateCreationTest.java
//!
//! Java's RefContextConsumer + RefContextCache pair is `CandidateFinder` here: `process_read` is
//! RefContextConsumer.processRead and `finish` is RefContextCache.altCandidates.

use super::util::*;
use crate::candidate::{soft_clip_insert, CandidateFinder, ReadContextCandidate};
use crate::constants::{set_test_base_qual_thresholds, LOW_BASE_QUAL_THRESHOLD};
use crate::refseq::RefSequence;
use crate::types::*;

/// VariantUtils.TEST_LEFT_FLANK / TEST_RIGHT_FLANK
pub(super) const TEST_LEFT_FLANK: &str = "ACCGCTGACT";
pub(super) const TEST_RIGHT_FLANK: &str = "CTGAGACTCA";

/// MockRefGenome.generateRandomBases: not random, cycles through A, C, G, T
pub(super) fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| ['A', 'C', 'G', 'T'][i % 4]).collect()
}

/// TestUtils.buildCigarString(alignedLength)
fn build_cigar_string(aligned_length: usize) -> String {
    format!("{aligned_length}M")
}

fn find_candidate<'a>(candidates: &'a [ReadContextCandidate], ref_: &str, alt: &str) -> Option<&'a ReadContextCandidate> {
    candidates.iter().find(|x| x.ctx.variant.ref_str() == ref_ && x.ctx.variant.alt_str() == alt)
}

// CandidateCreationTest.testSoftClipInsert
#[test]
fn test_soft_clip_insert() {
    let ref_base_str = generate_random_bases(100);

    let ref_bases = RefSequence::from_bases(100, ref_base_str.as_bytes().to_vec());

    let insert_bases = "AAAAA";

    // first an insert on the right
    let mut sc_start_ref_index = 51;
    let mut sc_bases = format!("{}{}", insert_bases, &ref_base_str[sc_start_ref_index..71]);
    let mut read_bases = format!("{}{}", &ref_base_str[20..sc_start_ref_index], sc_bases);
    let mut sc_read_index = 31;

    let alt_read =
        soft_clip_insert(120, 150, read_bases.as_bytes(), sc_bases.len() as i32, sc_read_index, &ref_bases, false);

    assert!(alt_read.is_some());
    let (ref_, alt) = alt_read.unwrap();
    assert_eq!("G", std::str::from_utf8(&ref_).unwrap());
    assert_eq!(format!("G{insert_bases}"), std::str::from_utf8(&alt).unwrap());

    // then on the left
    sc_start_ref_index = 0;
    sc_bases = format!("{}{}", &ref_base_str[sc_start_ref_index..20], insert_bases);
    read_bases = format!("{}{}", sc_bases, &ref_base_str[20..51]);
    sc_read_index = 0;

    let alt_read =
        soft_clip_insert(120, 150, read_bases.as_bytes(), sc_bases.len() as i32, sc_read_index, &ref_bases, true);

    assert!(alt_read.is_some());
    let (ref_, alt) = alt_read.unwrap();
    assert_eq!("T", std::str::from_utf8(&ref_).unwrap());
    assert_eq!(format!("T{insert_bases}"), std::str::from_utf8(&alt).unwrap());
}

// CandidateCreationTest.testSnvBeforeInsert
#[test]
fn test_snv_before_insert() {
    // create an SNV directly before a 1-base insert - confirm impact on read contexts and how reads are handled

    let ref_bases = format!("{}{}{}{}", &REF_BASES_200[0..100], TEST_LEFT_FLANK, "ACGTTCCAACCTTGCA", &REF_BASES_200[0..100]);
    let ref_sequence = RefSequence::from_bases(0, ref_bases.as_bytes().to_vec());

    //            110               120         130
    // 0123456789 0123456 7      8 9012345  6789012345
    // test flank ACGTTCC A>G insT ACCTTGCA AAAAAGGGGG

    let position = 117;
    let variant = SimpleVariant::new(position, b"A", b"G");

    let region = BaseRegion::new(0, 200);

    let cfg = test_config();
    let ref_data = empty_ref_data();
    let mut ref_context_consumer = CandidateFinder::new(&cfg, CHR_1, region, &ref_sequence, &[], &[], &ref_data);

    // send through a read with the SNV, then another with the SNV also immediately followed by an insert
    let mut read_bases = format!("{}{}{}", &ref_bases[100..117], variant.alt_str(), &ref_bases[118..135]);

    let mut cigar = build_cigar_string(read_bases.len());

    // the read's alignment start with the first base of the read context
    let mut read = build_sam_record(100, &cigar, &read_bases);
    read.mapq = 60;
    let mut read_clone = build_sam_record(100, &cigar, &read_bases);
    read_clone.mapq = 60;

    ref_context_consumer.process_read(&read);
    ref_context_consumer.process_read(&read_clone); // repeat to add support

    read_bases = format!("{}{}T{}", &ref_bases[100..117], variant.alt_str(), &ref_bases[118..135]);

    cigar = "18M1I17M".to_string();

    read = build_sam_record(100, &cigar, &read_bases);
    read.mapq = 60;
    read_clone = build_sam_record(100, &cigar, &read_bases);
    read_clone.mapq = 60;

    ref_context_consumer.process_read(&read);
    ref_context_consumer.process_read(&read_clone);

    let alt_contexts = ref_context_consumer.finish();
    assert_eq!(2, alt_contexts.len());

    let snv = find_candidate(&alt_contexts, variant.ref_str(), variant.alt_str());
    assert!(snv.is_some());

    let snv_insert = find_candidate(&alt_contexts, variant.ref_str(), "GT");
    assert!(snv_insert.is_some());
}

// CandidateCreationTest.testIlluminaLowQualCoreVariants
#[test]
fn test_illumina_low_qual_core_variants() {
    let ref_buffer = "A".repeat(10);

    // index / pos                                   0123456789012345
    let ref_bases = format!("{}{}{}{}{}", ref_buffer, TEST_LEFT_FLANK, "ACGTACGCAGCTGTCA", TEST_RIGHT_FLANK, ref_buffer);
    let ref_sequence = RefSequence::from_bases(0, ref_bases.as_bytes().to_vec());

    let mut variant1 = SimpleVariant::new(25, b"C", b"A");
    let mut variant2 = SimpleVariant::new(30, b"C", b"A");

    let region = BaseRegion::new(0, 200);
    let cfg = test_config();
    let ref_data = empty_ref_data();
    let mut ref_context_consumer = CandidateFinder::new(&cfg, CHR_1, region, &ref_sequence, &[], &[], &ref_data);

    // send through a read with the SNV, then another with the SNV also immediately followed by an insert
    let mut read_bases =
        format!("{}{}{}{}{}", &ref_bases[1..25], variant1.alt_str(), &ref_bases[26..30], variant2.alt_str(), &ref_bases[31..]);
    let mut cigar = build_cigar_string(read_bases.len());

    // the read's alignment start with the first base of the read context
    let mut read = build_sam_record(1, &cigar, &read_bases);
    read.quals[24] = LOW_BASE_QUAL_THRESHOLD;

    let mut read2 = build_sam_record(1, &cigar, &read_bases);
    read2.quals[24] = LOW_BASE_QUAL_THRESHOLD;

    ref_context_consumer.process_read(&read);
    ref_context_consumer.process_read(&read2);

    let mut alt_contexts = ref_context_consumer.finish();
    assert_eq!(1, alt_contexts.len());
    assert_eq!(30, alt_contexts[0].ctx.variant.pos);

    // refContextCache.clear(): `finish` consumes the finder and altCandidates() has already evicted every
    // ref context, so the cleared cache and its consumer are a new finder here
    let mut ref_context_consumer = CandidateFinder::new(&cfg, CHR_1, region, &ref_sequence, &[], &[], &ref_data);

    // test indels
    variant1 = SimpleVariant::new(25, b"C", b"CA");
    variant2 = SimpleVariant::new(30, b"CT", b"C");

    read_bases =
        format!("{}{}{}{}{}", &ref_bases[1..25], variant1.alt_str(), &ref_bases[26..30], variant2.alt_str(), &ref_bases[32..52]);
    cigar = "25M1I5M1D20M".to_string();

    // the read's alignment start with the first base of the read context
    read = build_sam_record(1, &cigar, &read_bases);
    read.quals[24] = LOW_BASE_QUAL_THRESHOLD;

    read2 = build_sam_record(1, &cigar, &read_bases);
    read2.quals[25] = LOW_BASE_QUAL_THRESHOLD;

    ref_context_consumer.process_read(&read);
    ref_context_consumer.process_read(&read2);

    alt_contexts = ref_context_consumer.finish();
    assert_eq!(1, alt_contexts.len());
    assert_eq!(30, alt_contexts[0].ctx.variant.pos);
}

// not in Java: -low_base_qual_threshold moves the band of bases too low to create a candidate on their own
#[test]
fn test_low_base_qual_threshold_option() {
    let ref_buffer = "A".repeat(10);
    let ref_bases = format!("{}{}{}{}{}", ref_buffer, TEST_LEFT_FLANK, "ACGTACGCAGCTGTCA", TEST_RIGHT_FLANK, ref_buffer);
    let ref_sequence = RefSequence::from_bases(0, ref_bases.as_bytes().to_vec());
    let variant1 = SimpleVariant::new(25, b"C", b"A");
    let variant2 = SimpleVariant::new(30, b"C", b"A");
    let region = BaseRegion::new(0, 200);
    let cfg = test_config();
    let ref_data = empty_ref_data();
    let read_bases =
        format!("{}{}{}{}{}", &ref_bases[1..25], variant1.alt_str(), &ref_bases[26..30], variant2.alt_str(), &ref_bases[31..]);
    let cigar = build_cigar_string(read_bases.len());
    let mut read = build_sam_record(1, &cigar, &read_bases);
    read.quals[24] = LOW_BASE_QUAL_THRESHOLD;
    let mut read2 = build_sam_record(1, &cigar, &read_bases);
    read2.quals[24] = LOW_BASE_QUAL_THRESHOLD;

    let candidates_at = |low_high: Option<(u8, u8)>| -> Vec<i32> {
        set_test_base_qual_thresholds(low_high);
        let mut finder = CandidateFinder::new(&cfg, CHR_1, region, &ref_sequence, &[], &[], &ref_data);
        finder.process_read(&read);
        finder.process_read(&read2);
        let mut positions: Vec<i32> = finder.finish().iter().map(|x| x.ctx.variant.pos).collect();
        set_test_base_qual_thresholds(None);
        positions.sort();
        positions
    };
    // Java's band (Q26 and below): only the second SNV, whose base is not low quality, is a candidate
    assert_eq!(candidates_at(None), [30]);
    assert_eq!(candidates_at(Some((LOW_BASE_QUAL_THRESHOLD, 30))), [30]);
    // with the band at Q20 and below the Q26 base no longer counts as low quality
    assert_eq!(candidates_at(Some((20, 30))), [25, 30]);
}
