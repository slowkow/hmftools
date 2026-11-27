//! BlosumSimilarityCalc and AnchorBlosumSearcher: how much a stretch of DNA looks like a V or J anchor,
//! by the BLOSUM62 similarity of its amino acids to the anchor templates of the gene file.
//!
//! The search tries every template at every position, as Java does, and keeps the first best score. The
//! score of a window without N comes from per-codon tables (the same sums as Java's string code); a
//! window with an N takes Java's string code (`similarity_score`).

use anyhow::{bail, Result};
use hmf_common::blosum::BlosumMapping;
use hmf_common::codon::{amino_acid_from_bases, codon_to_amino_acid};

use crate::compat;
use crate::constants::*;
use crate::genes::{GeneDatastore, Vj, VjGeneType};

fn round_up_to_multiple(v: i32, m: i32) -> i32 {
    ((v as f64 / m as f64).ceil() * m as f64) as i32
}

pub fn round_down_to_multiple(v: i32, m: i32) -> i32 {
    ((v as f64 / m as f64).floor() * m as f64) as i32
}

fn blosum() -> &'static BlosumMapping {
    static B: std::sync::OnceLock<BlosumMapping> = std::sync::OnceLock::new();
    B.get_or_init(BlosumMapping::default)
}

/// BlosumSimilarityCalc.calcBlosumDistance
fn blosum_distance(ref_aa: &[u8], seq_aa: &[u8]) -> Result<i32> {
    if ref_aa.len() != seq_aa.len() {
        return Ok(i32::MAX);
    }
    let b = blosum();
    match (b.calc_sequence_sum(ref_aa), b.calc_sequence_sum_pair(ref_aa, seq_aa)) {
        (Some(a), Some(c)) => Ok(a - c),
        _ => bail!("cannot calc blosum for seq: {} and {}", String::from_utf8_lossy(ref_aa), String::from_utf8_lossy(seq_aa)),
    }
}

/// BlosumSimilarityCalc.calcSimilarityScore(vj, templateDnaSeq, dnaSeq)
pub fn similarity_score(vj: Vj, template: &[u8], seq: &[u8]) -> Result<i32> {
    let mut t = template.to_vec();
    let mut s = seq.to_vec();
    if template.len() != seq.len() {
        let trim = round_down_to_multiple(template.len().min(seq.len()) as i32, 3) as usize;
        if vj == Vj::V {
            t = t[t.len() - trim..].to_vec();
            s = s[s.len() - trim..].to_vec();
        } else {
            t.truncate(trim);
            s.truncate(trim);
        }
    }
    let mut unknown = 0;
    if t.contains(&b'N') || s.contains(&b'N') {
        let (t0, s0) = (t.clone(), s.clone());
        for i in 0..t0.len() {
            if t0[i] == b'N' {
                t[i] = s0[i];
                unknown += 1;
            }
            if s0[i] == b'N' {
                s[i] = t0[i];
                unknown += 1;
            }
        }
        if t.contains(&b'N') || s.contains(&b'N') {
            bail!(
                "calcSimilarityScore: both ({} and {}) contain N at same location",
                String::from_utf8_lossy(template),
                String::from_utf8_lossy(seq)
            );
        }
    }
    let t_aa = amino_acid_from_bases(&t);
    let s_aa = amino_acid_from_bases(&s);
    let distance = blosum_distance(t_aa.as_bytes(), s_aa.as_bytes())?;
    let similarity = (MAX_BLOSUM_DIFF_PER_AA * t_aa.len() as i32).wrapping_sub(BLOSUM_SIMILARITY_SCORE_CONSTANT).wrapping_sub(distance);
    Ok(similarity.wrapping_sub(unknown * BLOSUM_UNKNOWN_BASE_PENALTY))
}

/// AnchorBlosumMatch; the template genes are indices into the datastore's templates
#[derive(Clone, Debug, PartialEq)]
pub struct AnchorBlosumMatch {
    pub anchor_start: i32,
    pub anchor_end: i32,
    pub template_anchor_seq: String,
    pub template_genes: Vec<usize>,
    pub similarity_score: i32,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SearchMode {
    AllowNegSimilarity,
    DisallowNegSimilarity,
}

/// IAnchorBlosumSearcher
pub trait AnchorSearcher: Sync {
    /// searchForAnchor(sequence, targetAnchorGeneTypes, mode, startOffset, endOffset)
    fn search(&self, seq: &[u8], types: &[VjGeneType], mode: SearchMode, start: i32, end: i32) -> Result<Option<AnchorBlosumMatch>>;
    /// searchForAnchor(readString, mode): all gene types over the whole sequence
    fn search_all(&self, seq: &[u8], mode: SearchMode) -> Result<Option<AnchorBlosumMatch>> {
        self.search(seq, &VjGeneType::ALL, mode, 0, seq.len() as i32)
    }
}

/// a template sequence with its amino acids
struct TemplateSeq {
    seq: String,
    aa: Vec<u8>,
    has_n: bool,
}

/// AnchorBlosumSearcher
pub struct AnchorBlosumSearcher<'a> {
    datastore: &'a GeneDatastore,
    min_partial_anchor_base_length: i32,
    /// per gene type, the sorted anchor sequences
    templates: Vec<Vec<TemplateSeq>>,
    /// self(t) - map(t, s) for amino acid bytes t, s; i16::MIN where Java throws
    distance: Vec<i16>,
}

const INVALID: i16 = i16::MIN;

impl<'a> AnchorBlosumSearcher<'a> {
    pub fn new(datastore: &'a GeneDatastore, min_partial_anchor_aa_length: i32) -> Self {
        let templates = VjGeneType::ALL
            .iter()
            .map(|&t| {
                datastore
                    .anchor_sequences(t)
                    .iter()
                    .map(|s| TemplateSeq { seq: s.clone(), aa: amino_acid_from_bases(s.as_bytes()).into_bytes(), has_n: s.contains('N') })
                    .collect()
            })
            .collect();
        let b = blosum();
        let mut distance = vec![INVALID; 256 * 256];
        for t in 0..256usize {
            for s in 0..256usize {
                if let (Some(st), Some(m)) = (b.self_mapping(t as u8), b.map(t as u8, s as u8)) {
                    distance[t * 256 + s] = (st - m) as i16;
                }
            }
        }
        AnchorBlosumSearcher { datastore, min_partial_anchor_base_length: min_partial_anchor_aa_length * 3, templates, distance }
    }
}

/// the read's amino acid at every codon start, and the running count of N
struct SeqInfo {
    aa: Vec<u8>,
    n_prefix: Vec<u32>,
}

impl SeqInfo {
    fn new(seq: &[u8]) -> SeqInfo {
        let mut aa = vec![b'.'; seq.len()];
        for (p, a) in aa.iter_mut().enumerate() {
            if p + 2 < seq.len() {
                *a = codon_to_amino_acid(seq, p);
            }
        }
        let mut n_prefix = Vec::with_capacity(seq.len() + 1);
        n_prefix.push(0);
        let mut n = 0;
        for &b in seq {
            if b == b'N' {
                n += 1;
            }
            n_prefix.push(n);
        }
        SeqInfo { aa, n_prefix }
    }
    fn has_n(&self, start: usize, end: usize) -> bool {
        self.n_prefix[end] != self.n_prefix[start]
    }
}

impl AnchorBlosumSearcher<'_> {
    /// tryMatchWithBlosum: the score of the template at `input_anchor_start`, None where Java returns null
    #[allow(clippy::too_many_arguments)]
    fn try_match(
        &self,
        gene_type: VjGeneType,
        seq: &[u8],
        info: &SeqInfo,
        input_anchor_start: i32,
        template: &TemplateSeq,
        mode: SearchMode,
    ) -> Result<Option<(i32, i32, i32)>> {
        let t_len = template.seq.len() as i32;
        let len = seq.len() as i32;
        let mut anchor_start = input_anchor_start;
        let mut anchor_end = input_anchor_start + t_len;
        let mut trim_left = 0i32;
        let mut trim_right = 0i32;
        if anchor_start < 0 {
            if gene_type.vj() == Vj::J {
                return Ok(None);
            }
            trim_left = round_up_to_multiple(-anchor_start, 3);
            anchor_start += trim_left;
        }
        // CI3: Java takes an anchor that ends at the last base of the sequence as running past its end
        // (>= in place of >), which rejects a V anchor there. Corrected: only anchors past the end.
        let past_end = if compat::java("CI3") { anchor_end >= len } else { anchor_end > len };
        if past_end {
            if gene_type.vj() == Vj::V {
                return Ok(None);
            }
            trim_right = round_up_to_multiple(anchor_end - len, 3);
            anchor_end -= trim_right;
        }
        // String.substring and dropLast as Java has them
        let trimmed_len = t_len - trim_left - trim_right;
        if trim_left > t_len || trimmed_len < 0 {
            bail!("java.lang.StringIndexOutOfBoundsException (template trim)");
        }
        if anchor_start > anchor_end || anchor_end > len {
            bail!("java.lang.StringIndexOutOfBoundsException (begin {anchor_start}, end {anchor_end}, length {len})");
        }
        let potential_len = anchor_end - anchor_start;
        if trimmed_len < t_len && potential_len < self.min_partial_anchor_base_length {
            return Ok(None);
        }
        if anchor_start < 0 || anchor_end > len {
            return Ok(None);
        }
        let (a0, a1) = (anchor_start as usize, anchor_end as usize);
        let score = if template.has_n || info.has_n(a0, a1) || trimmed_len != potential_len || trim_left % 3 != 0 {
            let t = &template.seq.as_bytes()[trim_left as usize..(t_len - trim_right) as usize];
            similarity_score(gene_type.vj(), t, &seq[a0..a1])?
        } else {
            // equal lengths, no N: 3 * codons - 6 - the BLOSUM distance of the amino acids
            let codons = (trimmed_len / 3) as usize;
            let t_aa = &template.aa[(trim_left / 3) as usize..(trim_left / 3) as usize + codons];
            let mut distance = 0i32;
            for (c, &ta) in t_aa.iter().enumerate() {
                let sa = info.aa[a0 + 3 * c];
                let d = self.distance[ta as usize * 256 + sa as usize];
                if d == INVALID {
                    // Java throws; take its exact path for the message
                    let t = &template.seq.as_bytes()[trim_left as usize..(t_len - trim_right) as usize];
                    return similarity_score(gene_type.vj(), t, &seq[a0..a1]).map(|s| Some((anchor_start, anchor_end, s)));
                }
                distance += d as i32;
            }
            MAX_BLOSUM_DIFF_PER_AA * codons as i32 - BLOSUM_SIMILARITY_SCORE_CONSTANT - distance
        };
        if score >= 0 || mode == SearchMode::AllowNegSimilarity { Ok(Some((anchor_start, anchor_end, score))) } else { Ok(None) }
    }
}

impl AnchorSearcher for AnchorBlosumSearcher<'_> {
    fn search(&self, seq: &[u8], types: &[VjGeneType], mode: SearchMode, start: i32, end: i32) -> Result<Option<AnchorBlosumMatch>> {
        let info = SeqInfo::new(seq);
        // (gene type, template index, start, end, score)
        let mut best: Option<(VjGeneType, usize, i32, i32, i32)> = None;
        for &gene_type in types {
            let templates = &self.templates[gene_type.ordinal()];
            for i in start..end {
                for (ti, template) in templates.iter().enumerate() {
                    let anchor_pos = if gene_type.vj() == Vj::V { i - template.seq.len() as i32 + 1 } else { i };
                    if let Some((s, e, score)) = self.try_match(gene_type, seq, &info, anchor_pos, template, mode)? {
                        if best.is_none_or(|b| score > b.4) {
                            best = Some((gene_type, ti, s, e, score));
                        }
                    }
                }
            }
        }
        Ok(best.map(|(gene_type, ti, s, e, score)| {
            let seq = &self.templates[gene_type.ordinal()][ti].seq;
            AnchorBlosumMatch {
                anchor_start: s,
                anchor_end: e,
                template_anchor_seq: seq.clone(),
                template_genes: self.datastore.by_type_anchor_sequence(gene_type, seq).to_vec(),
                similarity_score: score,
            }
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genes::VjAnchorTemplate;

    // BlosumSimilarityCalcTest
    #[test]
    fn similarity_calc_simple() {
        let template = "CAACCTGAAGATTTTGCAACTTACTACTGT";
        let test = "CAACCTGAAGATGTTGGATTTTACTTCTGT";
        assert_eq!(similarity_score(Vj::V, template.as_bytes(), test.as_bytes()).unwrap(), 2);
        // one longer than the other: V aligns right, J aligns left
        assert_eq!(similarity_score(Vj::V, format!("TGTTGT{template}").as_bytes(), test.as_bytes()).unwrap(), 2);
        assert_eq!(similarity_score(Vj::J, format!("{template}TGTTGT").as_bytes(), test.as_bytes()).unwrap(), 2);
    }

    #[test]
    fn similarity_calc_1() {
        assert_eq!(similarity_score(Vj::V, b"GAGGCTGAGGATGTTGGGGTTTATTACTGC", b"CAACCTGAAGATGTTGGATTTTACTTCTGT").unwrap(), 7);
    }

    #[test]
    fn similarity_with_n() {
        // an N is replaced by the other base and costs 2; N at the same place in both is an error
        assert_eq!(similarity_score(Vj::V, b"TGT", b"TNT").unwrap(), 3 - 6 - 2);
        assert!(similarity_score(Vj::V, b"TNT", b"TNT").is_err());
    }

    // AnchorBlosumSearcherTest
    fn datastore() -> GeneDatastore {
        let t = |gt, name: &str, seq: &str| VjAnchorTemplate::new(gt, name, "01", None, seq, None);
        GeneDatastore::new(
            vec![
                t(VjGeneType::Ighj, "IGHJ1", "TGGGGCCAGGGCACCCTGGTCACCGTCTCC"),
                t(VjGeneType::Ighj, "IGHJ6", "TGGGGGCAAGGGACCACGGTCACCGTCTCC"),
                t(VjGeneType::Ighv, "IGHV1-18", "AGATCTGACGACACGGCCGTGTATTACTGT"),
                t(VjGeneType::Ighv, "IGHV3-7", "AGAGCCGAGGACACGGCTGTGTATTACTGT"),
            ],
            Vec::new(),
        )
        .unwrap()
    }

    const V_ANCHOR: &str = "CGAGTCGAAGACACGGCTGTGTATTACTGT";
    const J_ANCHOR: &str = "TGGGGCCAAGGGACCACGGTCACCGTCTCC";
    const CDR3: &str = "GCGAGACCGAAATTTTATAGTAATGGCTTGGCGGGTATGGACGTC";

    fn search(searcher: &AnchorBlosumSearcher, seq: &str, gt: VjGeneType) -> AnchorBlosumMatch {
        searcher.search(seq.as_bytes(), &[gt], SearchMode::AllowNegSimilarity, 0, seq.len() as i32).unwrap().unwrap()
    }

    #[test]
    fn find_j_anchor() {
        let d = datastore();
        let searcher = AnchorBlosumSearcher::new(&d, 8);
        let full = format!("TCACTGTCTCTGCAAATGAATGACCTG{V_ANCHOR}{CDR3}{J_ANCHOR}CCAGCATAAAACACGATGGAAGTGATA");
        let v_end = full.find(V_ANCHOR).unwrap() + V_ANCHOR.len();
        let test = &full[v_end..];
        let m = search(&searcher, test, VjGeneType::Ighj);
        assert_eq!(m.anchor_start, 45);
        assert_eq!(&test[m.anchor_start as usize..m.anchor_end as usize], J_ANCHOR);
    }

    #[test]
    fn find_j_anchor_partial() {
        let d = datastore();
        let searcher = AnchorBlosumSearcher::new(&d, 7);
        let j = &J_ANCHOR[..J_ANCHOR.len() - 8];
        let full = format!("TCACTGTCTCTGCAAATGAATGACCTG{V_ANCHOR}{CDR3}{j}");
        let v_end = full.find(V_ANCHOR).unwrap() + V_ANCHOR.len();
        let test = &full[v_end..];
        let m = search(&searcher, test, VjGeneType::Ighj);
        assert_eq!(m.anchor_start, 45);
        assert_eq!(&test[m.anchor_start as usize..m.anchor_end as usize], &j[..j.len() - 1]);
    }

    #[test]
    fn find_v_anchor() {
        let d = datastore();
        let searcher = AnchorBlosumSearcher::new(&d, 8);
        let full = format!("TCACTGTCTCTGCAAATGAATGACCTG{V_ANCHOR}{CDR3}{J_ANCHOR}CCAGCATAAAACACGATGGAAGTGATA");
        let test = &full[..full.find(J_ANCHOR).unwrap()];
        let m = search(&searcher, test, VjGeneType::Ighv);
        assert_eq!(m.anchor_start, 27);
        assert_eq!(&test[m.anchor_start as usize..m.anchor_end as usize], V_ANCHOR);
    }

    #[test]
    fn find_v_anchor_partial() {
        let d = datastore();
        let searcher = AnchorBlosumSearcher::new(&d, 7);
        let v = &V_ANCHOR[8..];
        let full = format!("{v}{CDR3}{J_ANCHOR}CCAGCATAAAACACGATGGAAGTGATA");
        let test = &full[..full.find(J_ANCHOR).unwrap()];
        let m = search(&searcher, test, VjGeneType::Ighv);
        assert_eq!(m.anchor_start, 1);
        assert_eq!(&test[m.anchor_start as usize..m.anchor_end as usize], &v[1..]);
    }

    #[test]
    fn v_anchor_at_sequence_end() {
        // CI3: a V anchor that ends at the last base of the sequence
        let d = datastore();
        let searcher = AnchorBlosumSearcher::new(&d, 8);
        let seq = format!("{CDR3}{V_ANCHOR}");
        let end = seq.len() as i32;
        let run = |java| {
            compat::with_java(java, || searcher.search(seq.as_bytes(), &[VjGeneType::Ighv], SearchMode::DisallowNegSimilarity, 0, end).unwrap())
        };
        let fixed = run(false).unwrap();
        assert_eq!((fixed.anchor_start, fixed.anchor_end), (end - 30, end));
        assert!(run(true).is_none_or(|m| m.anchor_end < end));
    }

    #[test]
    fn fast_path_equals_string_path() {
        // every placement of every template on sequences with and without N: the codon tables give the
        // score of Java's string code
        let d = datastore();
        let searcher = AnchorBlosumSearcher::new(&d, 3);
        let seqs = [format!("TCACTG{V_ANCHOR}{CDR3}{J_ANCHOR}CCAGCA"), format!("TCANTG{V_ANCHOR}{CDR3}{J_ANCHOR}CCNGCA").replace("GACACG", "GANACG")];
        for s in &seqs {
            let info = SeqInfo::new(s.as_bytes());
            for gt in [VjGeneType::Ighv, VjGeneType::Ighj] {
                for t in &searcher.templates[gt.ordinal()] {
                    for p in -(t.seq.len() as i32)..(s.len() as i32) {
                        if let Some((a, b, score)) = searcher.try_match(gt, s.as_bytes(), &info, p, t, SearchMode::AllowNegSimilarity).unwrap() {
                            let tt = &t.seq.as_bytes()[(a - p) as usize..(b - p) as usize];
                            assert_eq!(score, similarity_score(gt.vj(), tt, &s.as_bytes()[a as usize..b as usize]).unwrap(), "{s} {p}");
                        }
                    }
                }
            }
        }
    }
}
