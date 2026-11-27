//! hmf-common `aminoacid.BlosumMapping`: BLOSUM62 scores between amino acids, with the stop codon 'X'
//! given a penalty against every other amino acid (-10 by default) and 0 against itself.

const AMINO_ACIDS: &[u8; 21] = b"ACDEFGHIKLMNPQRSTVWYX";

const BLOSUM62_COLUMNS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
const BLOSUM62: [[i8; 20]; 20] = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
];

/// BlosumMapping.DEFAULT_STOP_CODON_PENALTY
pub const DEFAULT_STOP_CODON_PENALTY: i32 = -10;

/// BlosumMapping, as a table over byte pairs; `None` where Java throws IllegalArgumentException ("invalid
/// amino acid"). Java takes the characters case-insensitively (Character.toUpperCase).
pub struct BlosumMapping {
    table: Vec<i16>,
}

const INVALID: i16 = i16::MIN;

impl Default for BlosumMapping {
    fn default() -> Self {
        Self::new(DEFAULT_STOP_CODON_PENALTY)
    }
}

impl BlosumMapping {
    pub fn new(stop_codon_penalty: i32) -> BlosumMapping {
        assert!((-128..=0).contains(&stop_codon_penalty));
        let mut table = vec![INVALID; 256 * 256];
        let index = |aa: u8| AMINO_ACIDS.iter().position(|&a| a == aa).unwrap();
        let mut m = [[0i16; 21]; 21];
        for (r, &row_aa) in BLOSUM62_COLUMNS.iter().enumerate() {
            for (c, &col_aa) in BLOSUM62_COLUMNS.iter().enumerate() {
                // Java stores a half matrix: the value set last for a pair wins (the matrix is symmetric)
                m[index(row_aa)][index(col_aa)] = BLOSUM62[r][c] as i16;
                m[index(col_aa)][index(row_aa)] = BLOSUM62[r][c] as i16;
            }
        }
        let stop = index(b'X');
        for i in 0..21 {
            m[i][stop] = stop_codon_penalty as i16;
            m[stop][i] = stop_codon_penalty as i16;
        }
        m[stop][stop] = 0;
        for (i, &a) in AMINO_ACIDS.iter().enumerate() {
            for (j, &b) in AMINO_ACIDS.iter().enumerate() {
                for a2 in [a, a.to_ascii_lowercase()] {
                    for b2 in [b, b.to_ascii_lowercase()] {
                        table[a2 as usize * 256 + b2 as usize] = m[i][j];
                    }
                }
            }
        }
        BlosumMapping { table }
    }

    /// map(aa1, aa2); None for an invalid amino acid
    #[inline]
    pub fn map(&self, aa1: u8, aa2: u8) -> Option<i32> {
        let v = self.table[aa1 as usize * 256 + aa2 as usize];
        if v == INVALID { None } else { Some(v as i32) }
    }

    /// selfMapping(aa)
    #[inline]
    pub fn self_mapping(&self, aa: u8) -> Option<i32> {
        self.map(aa, aa)
    }

    /// calcSequenceSum(sequence): the sum of self mappings
    pub fn calc_sequence_sum(&self, seq: &[u8]) -> Option<i32> {
        let mut total = 0;
        for &aa in seq {
            total += self.self_mapping(aa)?;
        }
        Some(total)
    }

    /// calcSequenceSum(seq1, seq2): 0 for different lengths, the self sum for equal sequences
    pub fn calc_sequence_sum_pair(&self, seq1: &[u8], seq2: &[u8]) -> Option<i32> {
        if seq1.len() != seq2.len() {
            return Some(0);
        }
        if seq1 == seq2 {
            return self.calc_sequence_sum(seq1);
        }
        let mut total = 0;
        for (&a, &b) in seq1.iter().zip(seq2) {
            total += self.map(a, b)?;
        }
        Some(total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values() {
        let b = BlosumMapping::default();
        assert_eq!(b.map(b'A', b'A'), Some(4));
        assert_eq!(b.map(b'W', b'W'), Some(11));
        assert_eq!(b.map(b'W', b'F'), Some(1));
        assert_eq!(b.map(b'F', b'W'), Some(1));
        assert_eq!(b.map(b'C', b'X'), Some(-10));
        assert_eq!(b.map(b'X', b'X'), Some(0));
        assert_eq!(b.map(b'.', b'A'), None);
        assert_eq!(b.map(b'B', b'A'), None);
        assert_eq!(b.map(b'a', b'r'), Some(-1));
        assert_eq!(b.calc_sequence_sum(b"CAR"), Some(9 + 4 + 5));
        assert_eq!(b.calc_sequence_sum_pair(b"CAR", b"CAK"), Some(9 + 4 + 2));
        assert_eq!(b.calc_sequence_sum_pair(b"CAR", b"CA"), Some(0));
    }
}
