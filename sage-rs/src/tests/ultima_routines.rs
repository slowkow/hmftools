//! Ported from the Java SAGE test suite: seqtech/UltimaRoutinesTest.java
//!
//! Java's RefContextConsumer + RefContextCache pair is `CandidateFinder` here; a cleared cache is a new
//! finder (see candidate_creation.rs).

use super::candidate_creation::TEST_LEFT_FLANK;
use super::util::*;
use crate::candidate::CandidateFinder;
use crate::constants::*;
use crate::refseq::RefSequence;
use crate::types::*;

/// refContextConsumer.processRead over each read with a new consumer and cache, then altCandidates():
/// whether the only candidate is the variant, or there are none
fn alt_candidates(ref_bases: &str, reads: &[&Read]) -> Vec<SimpleVariant> {
    let ref_sequence = RefSequence::from_bases(1, ref_bases.as_bytes().to_vec());
    let region = BaseRegion::new(1, 200);
    let cfg = test_config();
    let ref_data = empty_ref_data();
    let mut ref_context_consumer = CandidateFinder::new(&cfg, CHR_1, region, &ref_sequence, &[], &[], &ref_data);
    for read in reads {
        ref_context_consumer.process_read(read);
    }
    ref_context_consumer.finish().iter().map(|c| c.ctx.variant.clone()).collect()
}

/// SamRecordTestUtils.cloneSamRecord
fn clone_sam_record(read: &Read) -> Read {
    let mut clone = read.clone();
    clone.name = next_read_id();
    clone
}

// UltimaRoutinesTest.testLowQualBasesInCore
#[test]
fn test_low_qual_bases_in_core() {
    set_test_sequencing_type(Some(SequencingType::Ultima));
    // low-qual bases in core invalidates candidate evidence unless the variant has a long homopolymer
    let mut ref_bases = format!("{}{}{}{}", &REF_BASES_200[0..10], TEST_LEFT_FLANK, "AACTT", &REF_BASES_200[0..10]);

    let position = 23;
    let variant = SimpleVariant::new(position, b"C", b"G");

    // first a valid read without low-qual bases
    let mut read_bases = format!("{}{}{}", &ref_bases[0..22], variant.alt_str(), &ref_bases[23..]);
    let mut cigar = format!("{}M", read_bases.len());

    // the read's alignment start with the first base of the read context
    let mut read = build_sam_record(1, &cigar, &read_bases);
    let mut read2 = clone_sam_record(&read);

    let mut alt_contexts = alt_candidates(&ref_bases, &[&read, &read2]);
    assert_eq!(1, alt_contexts.len());
    assert!(alt_contexts.iter().any(|x| *x == variant));

    // repeat with low-qual bases
    let mut low_core_quals = "20=10-29";
    read2.low_qual_tag = Some(low_core_quals.to_string());
    let read3 = clone_sam_record(&read2);

    alt_contexts = alt_candidates(&ref_bases, &[&read2, &read3]);
    assert!(alt_contexts.is_empty());

    // now a high-qual read followed by a low-qual, still passes the min raw alt threshold
    read = build_sam_record(1, &cigar, &read_bases);
    read2 = clone_sam_record(&read);
    read2.low_qual_tag = Some(low_core_quals.to_string());

    alt_contexts = alt_candidates(&ref_bases, &[&read, &read2]);
    assert_eq!(1, alt_contexts.len());
    assert!(alt_contexts.iter().any(|x| *x == variant));

    // test again for a variant with a long HP repeat

    // max repeat will be the GA, but the 9xTs will cause the low-qual check to be skipped
    ref_bases = format!("{}{}{}{}{}", &REF_BASES_200[0..10], TEST_LEFT_FLANK, "AACGAGAGATTTT", "TTTTTACGTA", &REF_BASES_200[0..10]);

    read_bases = format!("{}{}{}", &ref_bases[0..22], variant.alt_str(), &ref_bases[23..]);
    cigar = format!("{}M", read_bases.len());

    read = build_sam_record(1, &cigar, &read_bases);

    low_core_quals = "10=25-34";
    read.low_qual_tag = Some(low_core_quals.to_string());
    read2 = clone_sam_record(&read);
    read2.low_qual_tag = Some(low_core_quals.to_string());

    alt_contexts = alt_candidates(&ref_bases, &[&read, &read2]);
    assert_eq!(1, alt_contexts.len());
    assert!(alt_contexts.iter().any(|x| *x == variant));
    set_test_sequencing_type(None);
}
