//! HlaSequenceLoci, HlaSequenceFile, SequenceCount, SequenceMatchType and HlaExonSequences.

use crate::hla::{AId, Alleles};
use crate::tok::{self, Tok, DEL, WILD};
use hmf_common::java_hash::{bucket_index, grown_capacity};

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum MatchType {
    Full,
    Wild,
    NoLoci,
    Mismatch,
}

impl MatchType {
    pub fn rank(self) -> i32 {
        match self {
            MatchType::Full => 0,
            MatchType::Wild => 1,
            MatchType::NoLoci => 1,
            MatchType::Mismatch => 3,
        }
    }
    pub fn is_better(self, other: MatchType) -> bool {
        self.rank() < other.rank()
    }
}

/// HlaSequenceLoci: an allele and its sequence, one token per locus
#[derive(Clone, Debug)]
pub struct SeqLoci {
    pub allele: AId,
    pub seq: Vec<Tok>,
    pub has_deletes: bool,
    pub has_inserts: bool,
    pub wildcard_count: u32,
    pub has_exon_boundary_wildcards: bool,
}

impl SeqLoci {
    /// new HlaSequenceLoci(allele, sequences): marks the allele as having wildcards if this sequence has any
    pub fn new(alleles: &mut Alleles, allele: AId, seq: Vec<Tok>) -> SeqLoci {
        let has_deletes = seq.contains(&DEL);
        let has_inserts = seq.iter().any(|t| tok::len(*t) > 1);
        let wildcard_count = seq.iter().filter(|t| **t == WILD).count() as u32;
        if wildcard_count > 0 {
            alleles.set_has_wildcards(allele);
        }
        SeqLoci { allele, seq, has_deletes, has_inserts, wildcard_count, has_exon_boundary_wildcards: false }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn has_wildcards(&self) -> bool {
        self.wildcard_count > 0
    }
    pub fn has_indels(&self) -> bool {
        self.has_deletes || self.has_inserts
    }
    #[inline]
    pub fn at(&self, locus: usize) -> Tok {
        self.seq[locus]
    }

    /// setExonBoundaryWildcards
    pub fn set_exon_boundary_wildcards(&mut self, boundaries: &[i32]) {
        for &locus in boundaries {
            if locus as usize >= self.seq.len() {
                return;
            }
            if self.seq[locus as usize] == WILD {
                self.has_exon_boundary_wildcards = true;
                return;
            }
        }
    }

    /// sequence(): all loci joined, deletes removed
    pub fn sequence_str(&self) -> String {
        let mut s = String::new();
        for t in &self.seq {
            if *t == DEL {
                continue;
            }
            if tok::is_single(*t) {
                tok::push(*t, &mut s);
            } else {
                s.push_str(&tok::to_string(*t).replace('.', ""));
            }
        }
        s
    }

    /// sequence(startLocus, endLocus)
    pub fn sequence_range(&self, start: usize, end: usize) -> String {
        tok::join(self.seq[start..=end].iter().copied())
    }

    /// sequence(indices): loci past the end give a wildcard
    pub fn sequence_at(&self, indices: &[i32]) -> Vec<u8> {
        let mut s = Vec::with_capacity(indices.len());
        for &i in indices {
            if (i as usize) < self.seq.len() {
                tok::push_bytes(self.seq[i as usize], &mut s);
            } else {
                s.push(b'*');
            }
        }
        s
    }

    /// determineMatchType(List<String> targetSequences, loci): locus by locus
    pub fn match_type_toks(&self, targets: &[Tok], loci: &[i32]) -> MatchType {
        if loci.is_empty() {
            return MatchType::Mismatch;
        }
        let mut wild = false;
        for (i, &locus) in loci.iter().enumerate() {
            if locus as usize >= self.seq.len() {
                return MatchType::Mismatch;
            }
            let s = self.seq[locus as usize];
            if s == WILD {
                wild = true;
            } else if s != targets[i] {
                return MatchType::Mismatch;
            }
        }
        if wild { MatchType::Wild } else { MatchType::Full }
    }

    /// determineMatchType(String targetSequence, indices): on the joined strings
    pub fn match_type_str(&self, target: &[u8], indices: &[i32]) -> MatchType {
        if indices.is_empty() {
            return MatchType::Mismatch;
        }
        let hla = self.sequence_at(indices);
        if hla == target {
            return MatchType::Full;
        }
        if hla.len() != target.len() {
            return MatchType::Mismatch;
        }
        let mut wild = 0;
        for i in 0..target.len() {
            let h = hla[i];
            if h != b'*' && h != target[i] {
                return MatchType::Mismatch;
            }
            if h == b'*' {
                wild += 1;
            }
        }
        if wild > 0 { MatchType::Wild } else { MatchType::Full }
    }

    pub fn consistent_with(&self, target: &[u8], indices: &[i32]) -> bool {
        self.match_type_str(target, indices) != MatchType::Mismatch
    }

    pub fn consistent_with_any<S: AsRef<[u8]>>(&self, targets: &[S], indices: &[i32]) -> bool {
        targets.iter().any(|t| self.match_type_str(t.as_ref(), indices) != MatchType::Mismatch)
    }
}

/// HlaSequenceFile.createFromReference: one locus per character, inserts between '|' characters
pub fn parse_reference_sequence(s: &str) -> Vec<Tok> {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut in_multi = false;
    let mut cur: Vec<u8> = Vec::new();
    for &c in b {
        let is_multi = c == b'|';
        if in_multi || is_multi {
            if in_multi && is_multi {
                in_multi = false;
                out.push(tok::intern(std::str::from_utf8(&cur).unwrap()));
                cur.clear();
            } else if is_multi {
                in_multi = true;
            } else {
                cur.push(c);
            }
        } else {
            out.push(tok::of_byte(c));
        }
    }
    out
}

/// HlaSequenceLoci.create(allele, sequence, reference): an aligned sequence against its reference
pub fn create_aligned(sequence: &str, reference: &str) -> Vec<Tok> {
    const IDENTICAL: u8 = b'-';
    const DELETION: u8 = b'.';
    const EXON_BOUNDARY: u8 = b'|';
    let s = sequence.as_bytes();
    let r = reference.as_bytes();
    let mut seqs: Vec<String> = Vec::new();
    let mut ins_len = 0usize;
    for i in 0..s.len() {
        let c = s[i];
        let rc = if i < r.len() { r[i] } else { IDENTICAL };
        let inserted = c != DELETION && (i >= r.len() || rc == DELETION);
        let ignored = (c == DELETION && rc == DELETION) || c == EXON_BOUNDARY;
        if ins_len > 0 && !inserted {
            let ins = &sequence[i - ins_len..i];
            if let Some(last) = seqs.last_mut() {
                last.push_str(ins);
            } else {
                // Java's sequences.set(-1, ...) throws here
                panic!("insert before the first locus");
            }
            ins_len = 0;
        }
        if inserted {
            ins_len += 1;
        } else if !ignored {
            if c == IDENTICAL {
                seqs.push((rc as char).to_string());
            } else {
                seqs.push((c as char).to_string());
            }
        }
    }
    seqs.iter().map(|x| tok::intern(x)).collect()
}

/// Codons.aminoAcidFromBases on a base string
pub fn amino_acids_from_bases(dna: &[u8]) -> String {
    let mut s = String::with_capacity(dna.len() / 3);
    let mut i = 0usize;
    while i + 2 < dna.len() {
        s.push(codon_to_amino_acid(dna[i], dna[i + 1], dna[i + 2]));
        i += 3;
    }
    s
}

const CODON_LOOKUP: &[u8; 64] = b"KNNKIIIMTTTTRSSRXYYXLFFLSSSSXCCWQHHQLLLLPPPPRRRREDDEVVVVAAAAGGGG";

fn base_to_int(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'T' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Codons.codonToAminoAcid: '.' for a codon with an unknown base
pub fn codon_to_amino_acid(a: u8, b: u8, c: u8) -> char {
    match (base_to_int(a), base_to_int(b), base_to_int(c)) {
        (Some(x), Some(y), Some(z)) => CODON_LOOKUP[x * 16 + y * 4 + z] as char,
        _ => '.',
    }
}

/// Codons.aminoAcidToCodon: the first codon in A, C, G, T order; None where Java throws
pub fn amino_acid_to_codon(aa: u8) -> Option<[u8; 3]> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    for &x in &BASES {
        for &y in &BASES {
            for &z in &BASES {
                if codon_to_amino_acid(x, y, z) as u8 == aa {
                    return Some([x, y, z]);
                }
            }
        }
    }
    None
}

/// HlaSequenceLoci.buildAminoAcidSequenceFromNucleotides (the HLA-Y protein sequences when the protein file
/// has none)
pub fn amino_acids_from_nucleotides(nuc: &SeqLoci, template: &SeqLoci) -> Vec<Tok> {
    let mut seq = String::new();
    let mut i = 0usize;
    while i + 2 < nuc.len() {
        let (a, b, c) = (nuc.seq[i], nuc.seq[i + 1], nuc.seq[i + 2]);
        if a == DEL || b == DEL || c == DEL {
            seq.push('.');
        } else {
            let s = tok::concat3(a, b, c);
            seq.push_str(&amino_acids_from_bases(s.as_bytes()));
        }
        i += 3;
    }
    create_aligned(&seq, &template.sequence_str())
}

/// The counts of one locus: tokens in first-insertion order (a Guava HashMultiset on a java.util.HashMap)
#[derive(Clone, Debug, Default)]
pub struct LocusCounts {
    pub entries: Vec<(Tok, u32)>,
    pub total: u32,
}

impl LocusCounts {
    #[inline]
    pub fn add(&mut self, t: Tok) {
        self.total += 1;
        for e in self.entries.iter_mut() {
            if e.0 == t {
                e.1 += 1;
                return;
            }
        }
        self.entries.push((t, 1));
    }
    pub fn count(&self, t: Tok) -> u32 {
        self.entries.iter().find(|e| e.0 == t).map(|e| e.1).unwrap_or(0)
    }
    pub fn contains(&self, t: Tok) -> bool {
        self.entries.iter().any(|e| e.0 == t)
    }
    /// the entries in the HashMap's iteration order
    pub fn java_order(&self) -> Vec<(Tok, u32)> {
        if self.entries.len() <= 1 {
            return self.entries.clone();
        }
        let cap = grown_capacity(16, self.entries.len());
        let mut v: Vec<(usize, (Tok, u32))> = self.entries.iter().map(|e| (bucket_index(tok::java_hash(e.0), cap), *e)).collect();
        v.sort_by_key(|x| x.0);
        v.into_iter().map(|x| x.1).collect()
    }
}

/// SequenceCount: per-locus counts of the nucleotides or amino acids of a set of fragments
#[derive(Clone, Debug, Default)]
pub struct SequenceCount {
    pub min_evidence_factor: f64,
    pub min_evidence_support: i32,
    /// by locus; `None` where no fragment has the locus
    pub loci: Vec<Option<LocusCounts>>,
}

static EMPTY_COUNTS: LocusCounts = LocusCounts { entries: Vec::new(), total: 0 };

impl SequenceCount {
    pub fn new(min_evidence_factor: f64, min_evidence_support: i32) -> Self {
        SequenceCount { min_evidence_factor, min_evidence_support, loci: Vec::new() }
    }

    /// the @VisibleForTesting constructor: counts per locus (min evidence factor 0), each locus's entries in
    /// the order given
    #[cfg(test)]
    pub fn from_counts(counts: &[Vec<(&str, u32)>]) -> Self {
        let mut c = SequenceCount::new(0.0, 1);
        for (l, m) in counts.iter().enumerate() {
            let mut lc = LocusCounts::default();
            for (t, n) in m {
                lc.entries.push((tok::intern(t), *n));
                lc.total += *n;
            }
            if c.loci.len() <= l {
                c.loci.resize(l + 1, None);
            }
            c.loci[l] = Some(lc);
        }
        c
    }

    #[inline]
    pub fn add(&mut self, locus: i32, t: Tok) {
        let l = locus as usize;
        if l >= self.loci.len() {
            self.loci.resize(l + 1, None);
        }
        self.loci[l].get_or_insert_with(LocusCounts::default).add(t);
    }

    pub fn has(&self, locus: i32) -> bool {
        locus >= 0 && (locus as usize) < self.loci.len() && self.loci[locus as usize].is_some()
    }

    /// get(locus): the counts, empty where none
    pub fn get(&self, locus: i32) -> &LocusCounts {
        if locus < 0 {
            return &EMPTY_COUNTS;
        }
        self.loci.get(locus as usize).and_then(|x| x.as_ref()).unwrap_or(&EMPTY_COUNTS)
    }

    /// seqCountsByLoci().keySet(), ascending
    pub fn loci_keys(&self) -> impl Iterator<Item = i32> + '_ {
        self.loci.iter().enumerate().filter(|(_, x)| x.is_some()).map(|(i, _)| i as i32)
    }

    pub fn depth(&self, locus: i32) -> u32 {
        self.get(locus).total
    }

    /// max(MIN_EVIDENCE_SUPPORT, (int) ceil(support * minEvidenceFactor))
    fn min_count(&self, total: u32) -> i64 {
        let c = (total as f64 * self.min_evidence_factor).ceil() as i32;
        self.min_evidence_support.max(c) as i64
    }

    /// getMinEvidenceSequences: tokens with enough support, in the multiset's iteration order
    pub fn min_evidence_sequences(&self, locus: i32) -> Vec<Tok> {
        let c = self.get(locus);
        if c.total == 0 {
            return Vec::new();
        }
        let min = self.min_count(c.total);
        c.java_order().into_iter().filter(|e| e.1 as i64 >= min).map(|e| e.0).collect()
    }

    /// the number of tokens with enough support (order-free)
    pub fn min_evidence_count(&self, locus: i32) -> usize {
        let c = self.get(locus);
        if c.total == 0 {
            return 0;
        }
        let min = self.min_count(c.total);
        c.entries.iter().filter(|e| e.1 as i64 >= min).count()
    }

    /// whether `t` has enough support at the locus
    pub fn has_min_evidence(&self, locus: i32, t: Tok) -> bool {
        let c = self.get(locus);
        if c.total == 0 {
            return false;
        }
        let min = self.min_count(c.total);
        c.entries.iter().any(|e| e.0 == t && e.1 as i64 >= min)
    }

    pub fn heterozygous_loci(&self) -> Vec<i32> {
        self.loci_keys().filter(|l| self.min_evidence_count(*l) > 1).collect()
    }
    pub fn homozygous_loci(&self) -> Vec<i32> {
        self.loci_keys().filter(|l| self.min_evidence_count(*l) == 1).collect()
    }

    /// getMaxCountSequence: "-" where the locus has no counts; ties go to the first in iteration order
    pub fn max_count_sequence(&self, locus: i32) -> Tok {
        let c = self.get(locus);
        if c.total == 0 {
            return b'-' as Tok;
        }
        let mut best: Option<(Tok, u32)> = None;
        for e in c.java_order() {
            match best {
                Some(b) if b.1 >= e.1 => {}
                _ => best = Some(e),
            }
        }
        best.unwrap().0
    }

    /// writeVertically: per locus, the five most frequent tokens and counts (raw depth first if given)
    pub fn vertical_text(&self, raw: Option<&SequenceCount>) -> String {
        let mut out = String::new();
        for locus in self.loci_keys() {
            out.push_str(&locus.to_string());
            if let Some(r) = raw {
                out.push('\t');
                out.push_str(&r.get(locus).total.to_string());
            }
            let mut v = self.get(locus).java_order();
            // stable sort by count, descending
            v.sort_by(|a, b| b.1.cmp(&a.1));
            for e in v.into_iter().take(5) {
                out.push('\t');
                out.push_str(&tok::to_string(e.0));
                out.push('\t');
                out.push_str(&e.1.to_string());
            }
            out.push('\n');
        }
        out
    }

    /// the private extractHeterozygousLociSequences: per locus, the supported tokens plus those of the extra
    /// sequences, where there are two or more
    pub fn het_loci_sequences(&self, extra: &[&SeqLoci]) -> std::collections::BTreeMap<i32, Vec<Tok>> {
        let mut out = std::collections::BTreeMap::new();
        for locus in self.loci_keys() {
            let mut set = self.min_evidence_sequences(locus);
            for s in extra {
                if locus as usize >= s.len() {
                    continue;
                }
                let t = s.at(locus as usize);
                if !set.contains(&t) {
                    set.push(t);
                }
            }
            if set.len() > 1 {
                out.insert(locus, set);
            }
        }
        out
    }
}

/// HlaExonSequences.ExonSequence
#[derive(Clone, Debug)]
pub struct ExonSequence {
    pub seq: Vec<Tok>,
    pub has_wildcards: bool,
}

const EXON_CHUNK_SIZE: usize = 46;

fn chunk(acids: &[Tok], out: &mut Vec<ExonSequence>) {
    let n = acids.len();
    if n == 0 {
        return;
    }
    let count = n.div_ceil(EXON_CHUNK_SIZE);
    let size = n.div_ceil(count);
    let mut start = 0;
    while start < n {
        let end = (start + size).min(n);
        let s = acids[start..end].to_vec();
        let w = s.contains(&WILD);
        out.push(ExonSequence { seq: s, has_wildcards: w });
        start = end;
    }
}

/// HlaExonSequences.create: the protein sequence split at the exon boundaries, and each exon into chunks of
/// at most 46 amino acids
pub fn exon_sequences(boundaries: &[i32], seq: &SeqLoci) -> Vec<ExonSequence> {
    let acids = &seq.seq;
    let mut out = Vec::new();
    let mut index = 0usize;
    for &b in boundaries {
        if index >= acids.len() {
            break;
        }
        let to = ((b + 1) as usize).min(acids.len());
        if to > index {
            chunk(&acids[index..to], &mut out);
        }
        index = (b + 1) as usize;
    }
    if index < acids.len() {
        chunk(&acids[index..], &mut out);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_parsing() {
        let s = parse_reference_sequence("AB|CDE|F.");
        assert_eq!(tok::join(s.iter().copied()), "ABCDEF.");
        assert_eq!(s.len(), 5);
        assert_eq!(tok::to_string(s[2]), "CDE");
        let a = create_aligned("--A.-", "ABCDE");
        assert_eq!(tok::join(a.iter().copied()), "ABA.E");
        let ins = create_aligned("A-BXY-", "A.B..C");
        // the inserted X and Y attach to B; the '.' against '.' is ignored
        assert_eq!(ins.iter().map(|t| tok::to_string(*t)).collect::<Vec<_>>(), vec!["A-", "BXY", "C"]);
    }

    #[test]
    fn codons() {
        assert_eq!(amino_acids_from_bases(b"ATGTAAGC"), "MX");
        assert_eq!(amino_acids_from_bases(b"ANG"), ".");
        assert_eq!(&amino_acid_to_codon(b'M').unwrap(), b"ATG");
        assert_eq!(&amino_acid_to_codon(b'X').unwrap(), b"TAA");
        assert!(amino_acid_to_codon(b'*').is_none());
    }

    #[test]
    fn counts_order() {
        let mut c = SequenceCount::new(0.0, 1);
        for b in [b'G', b'A', b'Q', b'A'] {
            c.add(0, b as Tok);
        }
        // buckets of 16: A=65->1, Q=81->1, G=71->7
        assert_eq!(c.min_evidence_sequences(0), vec![b'A' as Tok, b'Q' as Tok, b'G' as Tok]);
        assert_eq!(c.max_count_sequence(0), b'A' as Tok);
        assert_eq!(c.heterozygous_loci(), vec![0]);
    }
}
