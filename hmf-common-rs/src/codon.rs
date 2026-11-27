//! hmf-common `codon.Codons` and `codon.Nucleotides`, and htsjdk's `SequenceUtil` complement: amino acids
//! from DNA, base complements and reverse complements.

/// Codons.STOP_AMINO_ACID
pub const STOP_AMINO_ACID: u8 = b'X';
/// Codons.UNKNOWN: a codon with a base other than A, C, G or T
pub const UNKNOWN: u8 = b'.';

// the lookup of Codons, indexed by base (A 0, T 1, C 2, G 3) * 16 + * 4 + ...
const CODON_LOOKUP: &[u8; 64] = b"KNNKIIIMTTTTRSSRXYYXLFFLSSSSXCCWQHHQLLLLPPPPRRRREDDEVVVVAAAAGGGG";

#[inline]
fn base_to_int(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'T' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Codons.codonToAminoAcid(dna, index): UNKNOWN for a base other than A, C, G, T (upper case)
#[inline]
pub fn codon_to_amino_acid(dna: &[u8], index: usize) -> u8 {
    match (base_to_int(dna[index]), base_to_int(dna[index + 1]), base_to_int(dna[index + 2])) {
        (Some(a), Some(b), Some(c)) => CODON_LOOKUP[a * 16 + b * 4 + c],
        _ => UNKNOWN,
    }
}

/// Codons.aminoAcidFromBases: one amino acid per complete codon from the start
pub fn amino_acid_from_bases(dna: &[u8]) -> String {
    let mut s = String::with_capacity(dna.len() / 3);
    let mut i = 0;
    while i + 2 < dna.len() {
        s.push(codon_to_amino_acid(dna, i) as char);
        i += 3;
    }
    s
}

/// As `amino_acid_from_bases`, into a byte buffer (cleared first).
pub fn amino_acids_into(dna: &[u8], out: &mut Vec<u8>) {
    out.clear();
    let mut i = 0;
    while i + 2 < dna.len() {
        out.push(codon_to_amino_acid(dna, i));
        i += 3;
    }
}

/// Nucleotides.swapDnaBase / complement: A-T and C-G (upper case); anything else unchanged
#[inline]
pub fn swap_dna_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b,
    }
}

/// htsjdk SequenceUtil.complement: A-T and C-G in either case; anything else unchanged
#[inline]
pub fn htsjdk_complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b,
    }
}

/// htsjdk SequenceUtil.reverseComplement
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| htsjdk_complement(b)).collect()
}

/// Nucleotides.reverseComplementBases
pub fn reverse_complement_bases(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| swap_dna_base(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn amino_acids() {
        assert_eq!(amino_acid_from_bases(b"ATGTAATGGTGTNNA"), "MXWC.");
        assert_eq!(amino_acid_from_bases(b"AT"), "");
        assert_eq!(amino_acid_from_bases(b"ATGC"), "M");
        // some codons of the standard table
        for (codon, aa) in [("AAA", 'K'), ("TTT", 'F'), ("CCC", 'P'), ("GGG", 'G'), ("TAG", 'X'), ("TGA", 'X'), ("TGC", 'C'), ("CAT", 'H')] {
            assert_eq!(codon_to_amino_acid(codon.as_bytes(), 0) as char, aa, "{codon}");
        }
    }

    #[test]
    fn complements() {
        assert_eq!(reverse_complement(b"ACGTNacgt"), b"acgtNACGT".to_vec());
        assert_eq!(reverse_complement_bases(b"ACGTNa"), b"aNACGT".to_vec());
    }
}
