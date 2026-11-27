//! hmf-common LocalSequenceAligner (Smith-Waterman) and TEAL's TelomereMatcher.
//!
//! The aligner keeps Java's tie rules exactly: the diagonal wins ties, then the left move (a deletion), then
//! the up move (an insertion); a gap is extended when the neighbouring cell's traceback is the same move;
//! scores of 0 or less end the alignment; the first cell (in row-major order) with the highest score is
//! where the traceback starts.

use crate::telo::{CANONICAL_TELOMERE_SEQ, CANONICAL_TELOMERE_SEQ_REV, MIN_TELOMERE_MATCH_BASES};

const TRACEBACK_END: u8 = 0;
const TRACEBACK_DIAG: u8 = 1;
const TRACEBACK_LEFT: u8 = 2;
const TRACEBACK_UP: u8 = 3;

/// The alignment of the first sequence: its aligned range [start, end) and the number of matched bases
/// on the traceback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Alignment {
    pub first_start: usize,
    pub first_end: usize,
    pub num_matches: usize,
    pub score: i32,
}

pub struct LocalSequenceAligner {
    match_score: i32,
    mismatch_score: i32,
    gap_opening_score: i32,
    gap_extension_score: i32,
}

thread_local! {
    static WORK: std::cell::RefCell<(Vec<i32>, Vec<u8>)> = const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

impl LocalSequenceAligner {
    pub const fn new(match_score: i32, mismatch_score: i32, gap_opening_score: i32, gap_extension_score: i32) -> Self {
        LocalSequenceAligner { match_score, mismatch_score, gap_opening_score, gap_extension_score }
    }

    /// alignSequence(seq, refSeq)
    pub fn align(&self, seq: &[u8], ref_seq: &[u8]) -> Alignment {
        WORK.with(|w| {
            let (scores, tb) = &mut *w.borrow_mut();
            self.align_with(seq, ref_seq, scores, tb)
        })
    }

    fn align_with(&self, seq: &[u8], ref_seq: &[u8], scores: &mut Vec<i32>, tb: &mut Vec<u8>) -> Alignment {
        let rows = seq.len() + 1;
        let cols = ref_seq.len() + 1;
        scores.clear();
        scores.resize(rows * cols, 0);
        tb.clear();
        tb.resize(rows * cols, TRACEBACK_END);
        // the first row and column are 0 with traceback END
        let (mut hx, mut hy, mut highest) = (0usize, 0usize, 0i32);
        for x in 1..rows {
            let sx = seq[x - 1];
            let row = x * cols;
            let prev = (x - 1) * cols;
            for y in 1..cols {
                let mut diag = scores[prev + y - 1];
                diag += if sx == ref_seq[y - 1] { self.match_score } else { self.mismatch_score };
                let left_i = row + y - 1;
                let left = scores[left_i] + if tb[left_i] == TRACEBACK_LEFT { self.gap_extension_score } else { self.gap_opening_score };
                let up_i = prev + y;
                let up = scores[up_i] + if tb[up_i] == TRACEBACK_UP { self.gap_extension_score } else { self.gap_opening_score };
                let (score, t) = if diag >= left {
                    if diag >= up { (diag, TRACEBACK_DIAG) } else { (up, TRACEBACK_UP) }
                } else if left >= up {
                    (left, TRACEBACK_LEFT)
                } else {
                    (up, TRACEBACK_UP)
                };
                let i = row + y;
                if score <= 0 {
                    scores[i] = 0;
                    tb[i] = TRACEBACK_END;
                } else {
                    scores[i] = score;
                    tb[i] = t;
                    if score > highest {
                        highest = score;
                        hx = x;
                        hy = y;
                    }
                }
            }
        }
        // the traceback from the highest cell
        let (mut x, mut y) = (hx, hy);
        let mut num_matches = 0usize;
        loop {
            match tb[x * cols + y] {
                TRACEBACK_DIAG => {
                    if seq[x - 1] == ref_seq[y - 1] {
                        num_matches += 1;
                    }
                    x -= 1;
                    y -= 1;
                }
                TRACEBACK_LEFT => y -= 1,
                TRACEBACK_UP => x -= 1,
                _ => break,
            }
        }
        Alignment { first_start: x, first_end: hx, num_matches, score: highest }
    }
}

/// TelomereMatcher's aligner: the same score for a mismatch and a gap, for aligning to a TTAGGG template
static ALIGNER: LocalSequenceAligner = LocalSequenceAligner::new(1, -1, -1, -1);

/// TelomereMatcher.TelomereMatch
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TelomereMatch {
    pub match_start: usize,
    pub match_end: usize,
    pub num_matched_bases: usize,
    pub matched_sequence: Vec<u8>,
}

fn telomere_template(length: usize, g_rich: bool) -> Vec<u8> {
    let hexamer = if g_rich { CANONICAL_TELOMERE_SEQ } else { CANONICAL_TELOMERE_SEQ_REV }.as_bytes();
    (0..length).map(|i| hexamer[i % 6]).collect()
}

/// TelomereMatcher.findTelomereSegmentImpl. The match threshold Java takes is not used (TE3).
pub fn find_telomere_segment(seq: &[u8], g_rich: bool) -> Option<TelomereMatch> {
    // some insertion / deletion plus a full hexamer
    let template_length = (seq.len() as f64 * 1.2) as usize + 6;
    let a = ALIGNER.align(seq, &telomere_template(template_length, g_rich));
    if a.num_matches < MIN_TELOMERE_MATCH_BASES {
        return None;
    }
    Some(TelomereMatch {
        match_start: a.first_start,
        match_end: a.first_end,
        num_matched_bases: a.num_matches,
        matched_sequence: seq[a.first_start..a.first_end].to_vec(),
    })
}

#[cfg(test)]
pub fn find_g_telomere_segment(seq: &[u8]) -> Option<TelomereMatch> {
    find_telomere_segment(seq, true)
}

#[cfg(test)]
pub fn find_c_telomere_segment(seq: &[u8]) -> Option<TelomereMatch> {
    find_telomere_segment(seq, false)
}

/// TelomereMatcher.matchesTelomere (the threshold is not used, TE3)
pub fn matches_telomere(seq: &[u8], min_telomere_match_length: usize, g_rich: bool) -> bool {
    let triple: &[u8] = if g_rich { b"GGG" } else { b"CCC" };
    if !seq.windows(3).any(|w| w == triple) {
        return false;
    }
    if seq.len() < min_telomere_match_length {
        // shorter than the length: a full match is required
        return find_telomere_segment(seq, g_rich).is_some();
    }
    match find_telomere_segment(seq, g_rich) {
        Some(m) => m.num_matched_bases >= min_telomere_match_length,
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// TelomereMatcherTest.testCalcTelomereMatch
    #[test]
    fn calc_telomere_match() {
        let seq = b"ACTACCATTAGGGTTAGGGACTA";
        let m = find_g_telomere_segment(seq).unwrap();
        assert_eq!(m.num_matched_bases, 12);
        assert_eq!(m.match_start, 7);
        assert_eq!(m.match_end, 19);
        assert_eq!(m.matched_sequence, b"TTAGGGTTAGGG".to_vec());
        assert!(matches_telomere(seq, 12, true));
    }

    /// TelomereMatcherTest.testCalcTelomereNoMatch
    #[test]
    fn calc_telomere_no_match() {
        let seq = b"TGCAGCTTAACTGAGAGCCGCTCCTCTTCTCT";
        assert!(find_g_telomere_segment(seq).is_none());
        assert!(!matches_telomere(seq, 12, true));
    }

    /// TelomereMatcherTest.testCalcTelomereMatch2
    #[test]
    fn calc_telomere_match2() {
        let seq = b"ACTATTAGGGCTTAGGGACTA";
        let m = find_g_telomere_segment(seq).unwrap();
        assert_eq!(m.num_matched_bases, 12);
        assert_eq!(m.match_start, 4);
        assert_eq!(m.match_end, 17);
        assert_eq!(m.matched_sequence, b"TTAGGGCTTAGGG".to_vec());
        assert!(matches_telomere(seq, 12, true));
    }

    /// TelomereMatcherTest.testTelomereMatch3
    #[test]
    fn telomere_match3() {
        let seq = b"TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGGTAGGGTGAGGGTTAGGGTGAGGGTTAGGGTGAGGGTGAGGGTTAGGGTGAGGGTTAGGGTGAGGGAGAGGGTTAGGG";
        let m = find_g_telomere_segment(seq).unwrap();
        assert_eq!(m.num_matched_bases, 110);
        assert_eq!(m.match_start, 0);
        assert_eq!(m.match_end, 120);
        assert_eq!(m.matched_sequence, seq.to_vec());
        assert!(matches_telomere(seq, 12, true));
    }

    /// TelomereMatcherTest.testCalcTelomereMatch4
    #[test]
    fn calc_telomere_match4() {
        let seq = b"GTACTATTAGGGCTTAGGGACTA";
        let m = find_g_telomere_segment(seq).unwrap();
        assert_eq!(m.num_matched_bases, 12);
        assert_eq!(m.match_start, 6);
        assert_eq!(m.match_end, 19);
        assert_eq!(m.matched_sequence, b"TTAGGGCTTAGGG".to_vec());
        assert!(matches_telomere(seq, 12, true));
    }

    /// The tie rules: an empty alignment, and a gap opened then extended
    #[test]
    fn aligner_basics() {
        let a = LocalSequenceAligner::new(1, -4, -6, -1);
        let r = a.align(b"", b"ACGT");
        assert_eq!((r.first_start, r.first_end, r.num_matches, r.score), (0, 0, 0, 0));
        let r = a.align(b"ACGTACGT", b"TTACGTACGTTT");
        assert_eq!((r.first_start, r.first_end, r.num_matches, r.score), (0, 8, 8, 8));
        // an insertion of 3 bases costs 8 (open 6, extend 1 twice): 16 + 16 - 8 beats either half alone
        let half = b"ACGTTGCAACGTAGCT";
        let seq = [&half[..], b"GGG", &half[..]].concat();
        let r = a.align(&seq, &[&half[..], &half[..]].concat());
        assert_eq!((r.first_start, r.first_end, r.num_matches, r.score), (0, 35, 32, 24));
    }
}
