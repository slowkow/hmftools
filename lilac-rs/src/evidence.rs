//! PhasedEvidence, CombineEvidence, ExtendEvidence, PhasedEvidenceFactory, NucleotideFiltering and
//! Candidates: the phasing of amino acids and the candidate alleles it allows.

use std::collections::HashSet;

use hmf_common::java_hash::{bucket_index, grown_capacity};

use crate::config::{LilacConfig, Params};
use crate::fragment::{count_amino_acids, Fragment};
use crate::hla::{AId, Alleles, HlaGene};
use crate::refdata::{ExpectedAlleles, HlaContext};
use crate::seq::{SeqLoci, SequenceCount};
use crate::tok::WILD;

/// PhasedEvidence: the amino acid sequences seen over a set of loci, with their fragment counts
#[derive(Clone, Debug)]
pub struct PhasedEvidence {
    pub loci: Vec<i32>,
    /// evidence sequence to fragment count (the map's order never reaches the output)
    pub evidence: Vec<(String, i32)>,
    /// List.hashCode of the loci at construction (Java caches it; LI5 changes the loci later)
    pub hash: i32,
}

/// List<Integer>.hashCode
pub fn list_hash(v: &[i32]) -> i32 {
    let mut h: i32 = 1;
    for x in v {
        h = h.wrapping_mul(31).wrapping_add(*x);
    }
    h
}

impl PhasedEvidence {
    pub fn new(loci: Vec<i32>, evidence: Vec<(String, i32)>) -> PhasedEvidence {
        let hash = list_hash(&loci);
        PhasedEvidence { loci, evidence, hash }
    }

    /// PhasedEvidence.evidence(fragments, indices)
    pub fn from_fragments<'a>(frags: impl IntoIterator<Item = &'a Fragment>, indices: &[i32]) -> PhasedEvidence {
        let mut ev: Vec<(String, i32)> = Vec::new();
        let mut idx: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        for f in frags {
            if !f.contains_aa_loci(indices) {
                continue;
            }
            let s = f.aas_str(indices);
            match idx.get(&s) {
                Some(i) => ev[*i].1 += 1,
                None => {
                    idx.insert(s.clone(), ev.len());
                    ev.push((s, 1));
                }
            }
        }
        PhasedEvidence::new(indices.to_vec(), ev)
    }

    pub fn total(&self) -> i32 {
        self.evidence.iter().map(|x| x.1).sum()
    }
    pub fn min_evidence(&self) -> i32 {
        self.evidence.iter().map(|x| x.1).min().unwrap_or(0)
    }

    /// removeSingles(minTotalEvidence)
    pub fn remove_singles(self, min_total: i32) -> PhasedEvidence {
        if self.total() >= min_total {
            let ev = self.evidence.into_iter().filter(|x| x.1 > 1).collect();
            return PhasedEvidence { loci: self.loci, evidence: ev, hash: self.hash };
        }
        self
    }

    fn distinct_count<'a>(it: impl Iterator<Item = &'a str>) -> usize {
        it.collect::<HashSet<&str>>().len()
    }

    fn unambiguous_tail_length(&self) -> usize {
        let end = self.loci.len();
        for i in 0..self.loci.len() {
            let start = end - i - 1;
            let n = Self::distinct_count(self.evidence.iter().map(|x| substring(&x.0, start, end)));
            if n == self.evidence.len() {
                return i + 1;
            }
        }
        self.loci.len()
    }

    fn unambiguous_head_length(&self) -> usize {
        for i in 0..self.loci.len() {
            let len = i + 1;
            let n = Self::distinct_count(self.evidence.iter().map(|x| substring(&x.0, 0, len)));
            if n == self.evidence.len() {
                return len;
            }
        }
        self.loci.len()
    }

    pub fn unambiguous_head_indices(&self) -> Vec<i32> {
        let n = self.unambiguous_head_length().min(self.loci.len());
        self.loci[..n].to_vec()
    }
    pub fn unambiguous_tail_indices(&self) -> Vec<i32> {
        let min = self.loci.len() - self.unambiguous_tail_length();
        self.loci[min..].to_vec()
    }

    /// compareTo: by total evidence, descending, then by the smallest count, ascending
    pub fn compare(&self, other: &PhasedEvidence) -> std::cmp::Ordering {
        let t = other.total() - self.total();
        if t != 0 {
            return t.cmp(&0);
        }
        let m = self.min_evidence() - other.min_evidence();
        m.cmp(&0)
    }

    /// equals as Java has it: equal when the loci lists have the same size and differ at every position
    /// (LI11), or are the same object
    pub fn java_equals(&self, other: &PhasedEvidence) -> bool {
        if self.loci.len() != other.loci.len() {
            return false;
        }
        for i in 0..self.loci.len() {
            if self.loci[i] == other.loci[i] {
                return false;
            }
        }
        true
    }

    /// inconsistentEvidence: the evidence no candidate is consistent with
    pub fn inconsistent_evidence(&self, candidates: &[&SeqLoci]) -> PhasedEvidence {
        let ev = self.evidence.iter().filter(|x| !candidates.iter().any(|c| c.consistent_with(x.0.as_bytes(), &self.loci))).cloned().collect();
        PhasedEvidence { loci: self.loci.clone(), evidence: ev, hash: self.hash }
    }

    pub fn display(&self) -> String {
        format!(
            "PhasedEvidence(head={} tail={} loci={} types={} indices={:?}, total={})",
            self.unambiguous_head_length(),
            self.unambiguous_tail_length(),
            self.loci.len(),
            self.evidence.len(),
            self.loci,
            self.total()
        )
    }
}

/// String.substring(begin, end) on character positions; Java throws when out of range
fn substring(s: &str, begin: usize, end: usize) -> &str {
    let b = s.as_bytes();
    if begin > end || end > b.len() {
        panic!("java.lang.StringIndexOutOfBoundsException: begin {begin}, end {end}, length {}", b.len());
    }
    &s[begin..end]
}

fn sorted(mut v: Vec<i32>) -> Vec<i32> {
    v.sort();
    v
}

/// CombineEvidence.canCombine(left, right)
pub fn can_combine(left: &PhasedEvidence, right: &PhasedEvidence) -> bool {
    let inter = sorted(left.loci.iter().copied().filter(|x| right.loci.contains(x)).collect());
    if inter.is_empty() {
        return false;
    }
    let ul: Vec<i32> = left.loci.iter().copied().filter(|x| !inter.contains(x)).collect();
    let ur: Vec<i32> = right.loci.iter().copied().filter(|x| !inter.contains(x)).collect();
    let right_min = ur.iter().copied().min().unwrap_or(0);
    let left_min = ul.iter().copied().min().unwrap_or(0);
    if !ur.is_empty() && ul.iter().any(|x| *x > right_min) {
        return false;
    }
    if !ul.is_empty() && ur.iter().any(|x| *x < left_min) {
        return false;
    }
    let lstart = left.loci.len() - inter.len();
    let lc: HashSet<&str> = left.evidence.iter().map(|x| substring(&x.0, lstart, x.0.len())).collect();
    let rc: HashSet<&str> = right.evidence.iter().map(|x| substring(&x.0, 0, inter.len())).collect();
    lc.len() == rc.len() && lc.iter().all(|x| rc.contains(x))
}

pub fn can_combine3(left: &PhasedEvidence, common: &PhasedEvidence, right: &PhasedEvidence) -> bool {
    can_combine(left, common) && can_combine(common, right)
}

/// CombineEvidence.combineOverlapping
pub fn combine_overlapping(left: &PhasedEvidence, right: &PhasedEvidence) -> PhasedEvidence {
    let mut union: Vec<i32> = left.loci.clone();
    for x in &right.loci {
        if !left.loci.contains(x) {
            union.push(*x);
        }
    }
    union.sort();
    let inter = sorted(left.loci.iter().copied().filter(|x| right.loci.contains(x)).collect());
    let ul = left.loci.iter().filter(|x| !inter.contains(x)).count();
    let mut results: Vec<(String, i32)> = Vec::new();
    for (lk, lv) in &left.evidence {
        let lu = substring(lk, 0, ul);
        let lcom = substring(lk, ul, lk.len());
        for (rk, rv) in &right.evidence {
            let rcom = substring(rk, 0, inter.len());
            if lcom == rcom {
                let combined = format!("{lu}{rk}");
                let v = (*lv).min(*rv);
                match results.iter_mut().find(|x| x.0 == combined) {
                    Some(x) => x.1 = v,
                    None => results.push((combined, v)),
                }
            }
        }
    }
    PhasedEvidence::new(union, results)
}

pub fn combine(left: &PhasedEvidence, common: &PhasedEvidence, right: &PhasedEvidence) -> PhasedEvidence {
    let lc = combine_overlapping(left, common);
    combine_overlapping(&lc, right)
}

/// ExtendEvidence
pub struct ExtendEvidence<'a> {
    pub config: &'a LilacConfig,
    pub het_loci: Vec<i32>,
    pub frags: &'a [Fragment],
    pub expected: &'a ExpectedAlleles,
}

/// An evidence in the arena of the phasing loop (index = object identity)
type EvId = usize;

impl<'a> ExtendEvidence<'a> {
    fn min_total_fragments(&self, indices: &[i32]) -> i32 {
        self.expected.expected(indices) * self.config.min_fragments_per_allele
    }

    /// pairedEvidence
    pub fn paired_evidence(&self) -> Vec<PhasedEvidence> {
        let mut out = Vec::new();
        for i in 0..self.het_loci.len().saturating_sub(1) {
            let indices = vec![self.het_loci[i], self.het_loci[i + 1]];
            if !self.frags.iter().any(|f| f.contains_aa_loci(&indices)) {
                continue;
            }
            let min_total = self.min_total_fragments(&indices);
            let left = PhasedEvidence::from_fragments(self.frags.iter(), &indices[..1]);
            let right = PhasedEvidence::from_fragments(self.frags.iter(), &indices[1..]);
            let combined = PhasedEvidence::from_fragments(self.frags.iter(), &indices).remove_singles(self.config.min_fragments_to_remove_single);
            if combined.total() >= min_total && can_combine3(&left, &combined, &right) {
                if self.config.debug_phasing {
                    crate::info!("  paired Evidence: {}", combined.display());
                }
                out.push(combined);
            } else if self.config.debug_phasing {
                crate::info!("  failed paired Evidence: {}", combined.display());
            }
        }
        out.sort_by(|a, b| a.compare(b));
        out
    }

    /// merge(current, left, right): the combined evidence and the two it replaces, or nothing
    fn merge3(&self, arena: &[PhasedEvidence], current: EvId, left: EvId, right: EvId) -> Option<(PhasedEvidence, [EvId; 2])> {
        let (l, r) = (&arena[left], &arena[right]);
        let tail = l.unambiguous_tail_indices();
        let head = r.unambiguous_head_indices();
        let mut merge: Vec<i32> = tail.clone();
        for x in &head {
            if !tail.contains(x) {
                merge.push(*x);
            }
        }
        merge.sort();
        let filtered: Vec<&Fragment> = self.frags.iter().filter(|f| f.contains_aa_loci(&merge)).collect();
        let _ = current;
        if filtered.is_empty() {
            return None;
        }
        let min_total = self.min_total_fragments(&merge);
        let me = PhasedEvidence::from_fragments(filtered.into_iter(), &merge).remove_singles(self.config.min_fragments_to_remove_single);
        if can_combine3(l, &me, r) {
            let c = combine(l, &me, r);
            if c.total() >= min_total {
                return Some((c, [left, right]));
            }
        }
        None
    }

    /// merge(current, others)
    fn merge(&self, arena: &[PhasedEvidence], current: EvId, others: &[EvId]) -> Option<(PhasedEvidence, [EvId; 2])> {
        let cur = &arena[current];
        let min_e = cur.loci.iter().copied().min().unwrap_or(0);
        let max_e = cur.loci.iter().copied().max().unwrap_or(0);
        let with_max: Vec<EvId> = others.iter().copied().filter(|x| *x != current && arena[*x].loci.contains(&max_e)).collect();
        let with_min: Vec<EvId> = others.iter().copied().filter(|x| *x != current && arena[*x].loci.contains(&min_e)).collect();
        if !with_min.is_empty() && !with_max.is_empty() {
            let (left, right) = (with_min[0], with_max[0]);
            if arena[left].total() > arena[right].total() {
                let r = self.merge3(arena, current, left, current);
                if r.is_none() {
                    return self.merge3(arena, current, current, right);
                }
                return r;
            }
            let r = self.merge3(arena, current, current, right);
            if r.is_none() {
                return self.merge3(arena, current, left, current);
            }
            return r;
        }
        if !with_min.is_empty() {
            return self.merge3(arena, current, with_min[0], current);
        }
        if !with_max.is_empty() {
            return self.merge3(arena, current, current, with_max[0]);
        }
        None
    }
}

/// Membership in Java's HashSet<PhasedEvidence>: the same hash and the same object or equals (LI11); with
/// the correction, the same object only
fn java_set_contains(arena: &[PhasedEvidence], set: &[EvId], x: EvId) -> bool {
    set.iter().any(|y| {
        if *y == x {
            return true;
        }
        let eq = arena[*y].hash == arena[x].hash && arena[x].java_equals(&arena[*y]);
        if eq {
            crate::debug!("LI11: phased evidence {:?} taken as equal to {:?}", arena[x].loci, arena[*y].loci);
        }
        eq && crate::compat::java("LI11")
    })
}

/// the elements of a HashSet<PhasedEvidence> built by adding these in order, in its iteration order
fn java_evidence_set(arena: &[PhasedEvidence], adds: impl IntoIterator<Item = EvId>) -> Vec<EvId> {
    let mut set: Vec<EvId> = Vec::new();
    for x in adds {
        if !java_set_contains(arena, &set, x) {
            set.push(x);
        }
    }
    let cap = grown_capacity(16, set.len());
    let mut keyed: Vec<(usize, EvId)> = set.into_iter().map(|e| (bucket_index(arena[e].hash, cap), e)).collect();
    keyed.sort_by_key(|x| x.0);
    keyed.into_iter().map(|x| x.1).collect()
}

/// PhasedEvidenceFactory.evidence(expectedAlleles, fragments)
pub fn phased_evidence(config: &LilacConfig, params: &Params, gene: HlaGene, expected: &ExpectedAlleles, frags: &[Fragment]) -> Vec<PhasedEvidence> {
    let counts = count_amino_acids(params, params.min_evidence_factor, frags.iter());
    let het = counts.heterozygous_loci();
    if config.debug_phasing {
        crate::debug!("phasing {} records:", gene.long_name());
        crate::debug!("  heterozygous Indices: {het:?}");
    }
    let ext = ExtendEvidence { config, het_loci: het, frags, expected };
    let mut arena: Vec<PhasedEvidence> = Vec::new();
    let mut finalised: Vec<EvId> = Vec::new();
    let mut unprocessed: Vec<EvId> = Vec::new();
    for e in ext.paired_evidence() {
        arena.push(e);
        unprocessed.push(arena.len() - 1);
    }
    while !unprocessed.is_empty() {
        let top = unprocessed.remove(0);
        let others = java_evidence_set(&arena, finalised.iter().copied().chain(unprocessed.iter().copied()));
        match ext.merge(&arena, top, &others) {
            Some((parent, children)) => {
                // removeAll(Sets.newHashSet(left, right))
                let child_set = {
                    let mut s: Vec<EvId> = Vec::new();
                    for c in children {
                        if !java_set_contains(&arena, &s, c) {
                            s.push(c);
                        }
                    }
                    s
                };
                finalised.retain(|x| !java_set_contains(&arena, &child_set, *x));
                unprocessed.retain(|x| !java_set_contains(&arena, &child_set, *x));
                arena.push(parent);
                unprocessed.push(arena.len() - 1);
                if config.debug_phasing {
                    crate::debug!("  Produced child: {}", arena.last().unwrap().display());
                }
            }
            None => finalised.push(top),
        }
        // Collections.sort: stable, with a consistent comparator
        unprocessed.sort_by(|a, b| arena[*a].compare(&arena[*b]));
    }
    let mut out: Vec<PhasedEvidence> = finalised.into_iter().map(|i| arena[i].clone()).collect();
    out.sort_by_key(|e| e.loci[0]);
    out
}

/// NucleotideFiltering.filterCandidatesOnAminoAcidBoundaries
pub fn filter_on_amino_acid_boundaries<'s>(params: &Params, boundaries: &[i32], gene: HlaGene, candidates: Vec<&'s SeqLoci>, frags: &[Fragment]) -> Vec<&'s SeqLoci> {
    let mut results = candidates;
    for &b in boundaries {
        let ns = b * 3;
        let span = frags
            .iter()
            .filter(|f| (gene == HlaGene::Drb3 || gene == HlaGene::Drb4 || f.read_gene == gene) && f.contains_all_nucleotides(&[ns, ns + 1, ns + 2]))
            .count();
        if (span as i64) < params.min_depth_filter as i64 {
            continue;
        }
        let seqs = |loci: &[i32]| -> Vec<Vec<u8>> {
            let mut set: Vec<Vec<u8>> = Vec::new();
            for f in frags {
                if !f.contains_all_nucleotides(loci) {
                    continue;
                }
                let s = f.nucleotides_str(loci.iter().copied());
                if !set.contains(&s) {
                    set.push(s);
                }
            }
            set
        };
        let starts = seqs(&[ns]);
        let ends = seqs(&[ns + 1, ns + 2]);
        results.retain(|s| s.consistent_with_any(&starts, &[ns]) && s.consistent_with_any(&ends, &[ns + 1, ns + 2]));
    }
    results
}

/// NucleotideFiltering.calcNucleotideHeterogygousLoci: per gene, the heterozygous nucleotide loci at its
/// exon boundaries
pub fn nucleotide_het_loci(gc: &crate::refdata::GeneCache, het: &[i32]) -> Vec<(HlaGene, Vec<i32>)> {
    gc.gene_names
        .iter()
        .map(|g| {
            let mut nb = Vec::new();
            for b in gc.aa_boundaries(*g) {
                nb.extend([b * 3, b * 3 + 1, b * 3 + 2]);
            }
            (*g, het.iter().copied().filter(|x| nb.contains(x)).collect())
        })
        .collect()
}

/// Candidates.filterSequencesByMinSupport
pub fn filter_by_min_support<'s>(params: &Params, candidates: Vec<&'s SeqLoci>, counts: &SequenceCount, boundaries: &[i32], raw: Option<&SequenceCount>) -> Vec<&'s SeqLoci> {
    let mut cands = candidates;
    for locus in counts.loci_keys().collect::<Vec<_>>() {
        if boundaries.contains(&locus) {
            continue;
        }
        if let Some(r) = raw {
            if (r.get(locus).total as i64) < params.min_depth_filter as i64 {
                continue;
            }
        }
        let expected = counts.min_evidence_sequences(locus);
        if expected.is_empty() {
            continue;
        }
        cands.retain(|s| {
            if (locus as usize) < s.len() {
                let t = s.at(locus as usize);
                return t == WILD || expected.contains(&t);
            }
            false
        });
    }
    cands
}

/// Candidates.unphasedCandidates
pub fn unphased_candidates(params: &Params, alleles: &Alleles, ctx: &HlaContext, frags: &[Fragment], aa_seqs: &[SeqLoci], nuc_seqs: &[SeqLoci], common: &[AId], raw_aa: Option<&SequenceCount>) -> Vec<AId> {
    let gene = ctx.gene;
    crate::debug!("gene({}) determining un-phased candidates from frags({})", gene.long_name(), frags.len());
    let counts = count_amino_acids(params, params.min_evidence_factor, frags.iter());
    let gene_cands: Vec<&SeqLoci> = aa_seqs.iter().filter(|s| alleles.gene(s.allele) == gene).collect();
    let mut bset = ctx.aa_boundaries.clone();
    bset.sort();
    bset.dedup();
    let aa_cands = filter_by_min_support(params, gene_cands, &counts, &bset, raw_aa);
    let aa_alleles: Vec<AId> = aa_cands.iter().map(|s| s.allele).collect();
    let specific: HashSet<AId> = aa_alleles.iter().map(|a| alleles.four_digit(*a)).collect();
    if specific.is_empty() {
        crate::warn!("gene({}) no candidates after amino acid filtering - reverting to common allele gene candidates", gene.long_name());
        return common.iter().copied().filter(|a| alleles.gene(*a) == gene).collect();
    }
    crate::debug!("gene({}) {} candidates after amino acid filtering", gene.long_name(), aa_cands.len());
    let nuc_cands: Vec<&SeqLoci> = nuc_seqs.iter().filter(|s| specific.contains(&alleles.four_digit(s.allele))).collect();
    let nuc_specific = filter_on_amino_acid_boundaries(params, &ctx.aa_boundaries, gene, nuc_cands, frags);
    let mut out: Vec<AId> = Vec::new();
    for s in nuc_specific {
        let f = alleles.four_digit(s.allele);
        if !out.contains(&f) {
            out.push(f);
        }
    }
    if out.is_empty() {
        crate::warn!("gene({}) 0 candidates after exon boundary filtering - reverting to amino acid candidates", gene.long_name());
        return aa_alleles;
    }
    crate::debug!("gene({}) {} candidates after exon boundary filtering", gene.long_name(), out.len());
    out
}

/// Candidates.filterCandidates: drops candidates inconsistent with any phased evidence. Low-depth loci are
/// dropped from the evidence; Java drops them from the evidence's own loci list (LI5), which the haplotype
/// QC later reads against the full sequences.
pub fn filter_candidates<'s>(params: &Params, initial: Vec<&'s SeqLoci>, evidence: &mut [PhasedEvidence], raw: Option<&SequenceCount>) -> Vec<&'s SeqLoci> {
    let mut cands = initial;
    for ev in evidence.iter_mut() {
        let mut target_loci: Vec<i32> = ev.loci.clone();
        let mut targets: Vec<Vec<u8>> = ev.evidence.iter().map(|x| x.0.as_bytes().to_vec()).collect();
        let mut low: Vec<usize> = Vec::new();
        for (j, l) in target_loci.iter().enumerate() {
            if let Some(r) = raw {
                if (r.get(*l).total as i64) < params.min_depth_filter as i64 {
                    low.push(j);
                }
            }
        }
        if low.len() as i64 >= target_loci.len() as i64 - 1 {
            continue;
        }
        if !low.is_empty() {
            crate::debug!("LI5: phased evidence at loci {:?}: {} low-depth loci dropped", ev.loci, low.len());
            low.sort_by(|a, b| b.cmp(a));
            for &j in &low {
                target_loci.remove(j);
                for t in targets.iter_mut() {
                    if j >= t.len() {
                        panic!("java.lang.StringIndexOutOfBoundsException: index {j}, length {}", t.len());
                    }
                    t.remove(j);
                }
            }
            if crate::compat::java("LI5") {
                ev.loci = target_loci.clone();
            }
        }
        cands.retain(|s| s.consistent_with_any(&targets, &target_loci));
    }
    cands
}
