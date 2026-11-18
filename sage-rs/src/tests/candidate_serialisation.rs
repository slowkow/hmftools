//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/candidate/CandidateSerialisationTest.java
//!
//! Java's vcf/CandidateSerialisation rebuilds a Candidate and its read context from a VCF record
//! (`toCandidate`), which is what append mode reads its input with. The port does not implement append
//! mode and has no equivalent of `toCandidate`, so:
//!
//! - testCandidateCreation is only partly ported: building the three read contexts and checking they are
//!   valid is kept; the `toContext` -> `toCandidate` round trip and its
//!   `recreatedContext.matches(readContext)` assertions are skipped.
//! - testCreateFromVariantContext is skipped entirely: it only calls `toCandidate` on a hand-built
//!   VariantContext with an RC_INFO attribute.

use super::util::*;
use crate::constants::DEFAULT_FLANK_LENGTH;
use crate::read_context::create_context;
use crate::refseq::RefSequence;
use crate::types::*;

//             10        20        30        40        50
//   0123456789012345678901234567890123456789012345678901
const REF_BASES: &str = "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

fn ref_sequence() -> RefSequence {
    RefSequence::from_bases(1, REF_BASES.as_bytes().to_vec())
}

// CandidateSerialisationTest.testCandidateCreation (read context creation only, see the module comment)
#[test]
fn test_candidate_creation() {
    let ref_sequence = ref_sequence();

    // test 1: a simple SNV
    let mut position: usize = 21;
    let mut read_bases = format!("{}C{}", &REF_BASES[1..position], &REF_BASES[position + 1..50]);

    let mut read_cigar = format!("{}M", read_bases.len());
    let mut read = build_sam_record(1, &read_cigar, &read_bases);

    let mut var = SimpleVariant::new(position as i32, b"A", b"C");

    let read_context = create_context(&var, &read, position as i32 - 1, &ref_sequence, DEFAULT_FLANK_LENGTH);

    assert!(read_context.is_some_and(|x| x.is_valid()));

    // skipped: Candidate -> CandidateSerialisation.toContext -> toCandidate, recreatedContext.matches(readContext)

    // a delete
    position = 17;
    var = SimpleVariant::new(position as i32, b"ACCC", b"A");

    read_bases = format!("{}{}", &REF_BASES[1..position + 1], &REF_BASES[position + 4..50]);
    read_cigar = "17M3D29M".to_string();
    read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = create_context(&var, &read, position as i32 - 1, &ref_sequence, DEFAULT_FLANK_LENGTH);

    assert!(read_context.is_some_and(|x| x.is_valid()));

    // skipped: the serialisation round trip, as above

    // an insert
    position = 21;
    var = SimpleVariant::new(position as i32, b"A", b"ACCC");

    read_bases = REF_BASES[1..position + 1].to_string(); // up to and including the ref base
    read_bases += &format!("CCC{}", &REF_BASES[position + 2..50]);
    read_cigar = "21M3I27M".to_string();
    read = build_sam_record(1, &read_cigar, &read_bases);

    let read_context = create_context(&var, &read, position as i32 - 1, &ref_sequence, DEFAULT_FLANK_LENGTH);

    assert!(read_context.is_some_and(|x| x.is_valid()));

    // skipped: the serialisation round trip, as above
}

// CandidateSerialisationTest.testCreateFromVariantContext: skipped, CandidateSerialisation.toCandidate
// (candidate file / append mode input) is not implemented in the port.
