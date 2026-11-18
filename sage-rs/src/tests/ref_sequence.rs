//! Ported from the Java SAGE test suite; see the Java class named in the module list in mod.rs.
//! common/RefSequenceTest.java

use crate::refseq::RefSequence;

/// MockRefGenome.generateRandomBases: not random, cycles through A, C, G, T
fn generate_random_bases(length: usize) -> String {
    const DNA_BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    (0..length).map(|i| DNA_BASES[i % 4]).collect()
}

// RefSequenceTest.testRefSequence
#[test]
fn test_ref_sequence() {
    let ref_bases = generate_random_bases(101);
    let ref_sequence = RefSequence::from_bases(100, ref_bases.as_bytes().to_vec());

    assert_eq!(200, ref_sequence.end);
    assert_eq!(10, ref_sequence.index(110));
    // RefSequence.length(): the port has no such method, the bases are a public field
    assert_eq!(101, ref_sequence.bases.len());
    assert_eq!(Some(ref_bases[5..16].as_bytes()), ref_sequence.position_bases(105, 115));
    // RefSequence.trinucleotideContext(110) is not ported (the port only takes trinucleotides from a read
    // context's ref bases), so the same three bases are read with base_range
    assert_eq!(Some(ref_bases[9..12].as_bytes()), ref_sequence.base_range(109, 111));
}
