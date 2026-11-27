//! Fragment, FragmentUtils, NucleotideFragmentFactory, the nucleotide and amino acid enrichments, the
//! gene enrichment of NucleotideGeneEnrichment, AminoAcidFragmentPipeline and SequenceCountDiff.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::Arc;

use anyhow::{bail, Result};
use hmf_common::java_hash::{bucket_index, grown_capacity, java_string_hash};

use crate::bam::{is_uncertain_base_qual, Read, ReadCtx, Rec};
use crate::config::Params;
use crate::constants::*;
use crate::hla::HlaGene;
use crate::refdata::{GeneCache, GeneEnrichment, HlaContext};
use crate::seq::{self, SeqLoci, SequenceCount};
use crate::tok::{self, Tok, DEL};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Nuc {
    pub bases: Tok,
    pub low_qual: bool,
    pub qual: i8,
}

impl Nuc {
    pub fn high_qual(bases: Tok) -> Nuc {
        Nuc { bases, low_qual: false, qual: MISSING_BASE_QUAL }
    }
}

/// isLowBaseQual on a Java byte (signed: quals of 128 and above are negative and so low)
pub fn is_low_qual_byte(q: i8) -> bool {
    (q as i32) <= LOW_BASE_QUAL_THRESHOLD as i32
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Scope {
    BaseQualFiltered,
    HlaY,
    NoHetLoci,
    UnmatchedAminoAcid,
    WildOnly,
    Candidate,
    Solution,
    Unset,
}

impl Scope {
    pub fn is_unmatched(self) -> bool {
        matches!(self, Scope::Candidate | Scope::UnmatchedAminoAcid | Scope::WildOnly)
    }
    pub fn name(self) -> &'static str {
        match self {
            Scope::BaseQualFiltered => "BASE_QUAL_FILTERED",
            Scope::HlaY => "HLA_Y",
            Scope::NoHetLoci => "NO_HET_LOCI",
            Scope::UnmatchedAminoAcid => "UNMATCHED_AMINO_ACID",
            Scope::WildOnly => "WILD_ONLY",
            Scope::Candidate => "CANDIDATE",
            Scope::Solution => "SOLUTION",
            Scope::Unset => "UNSET",
        }
    }
}

/// A sorted map from locus to value (Java's TreeMap<Integer, ...>)
#[derive(Clone, Debug, Default)]
pub struct LociMap<T: Copy> {
    pub v: Vec<(i32, T)>,
}

impl<T: Copy> LociMap<T> {
    pub fn new() -> Self {
        LociMap { v: Vec::new() }
    }
    /// from (locus, value) pairs; a repeated locus keeps the last value (TreeMap.put)
    pub fn from_pairs(pairs: impl IntoIterator<Item = (i32, T)>) -> Self {
        let mut v: Vec<(i32, T)> = pairs.into_iter().collect();
        let sorted = v.windows(2).all(|w| w[0].0 < w[1].0);
        if !sorted {
            let mut m = std::collections::BTreeMap::new();
            for (k, x) in v {
                m.insert(k, x);
            }
            v = m.into_iter().collect();
        }
        LociMap { v }
    }
    #[inline]
    fn find(&self, k: i32) -> Result<usize, usize> {
        self.v.binary_search_by_key(&k, |x| x.0)
    }
    #[inline]
    pub fn get(&self, k: i32) -> Option<T> {
        self.find(k).ok().map(|i| self.v[i].1)
    }
    #[inline]
    pub fn contains(&self, k: i32) -> bool {
        self.find(k).is_ok()
    }
    pub fn insert(&mut self, k: i32, x: T) {
        match self.find(k) {
            Ok(i) => self.v[i].1 = x,
            Err(i) => self.v.insert(i, (k, x)),
        }
    }
    pub fn len(&self) -> usize {
        self.v.len()
    }
    pub fn is_empty(&self) -> bool {
        self.v.is_empty()
    }
    pub fn first_key(&self) -> Option<i32> {
        self.v.first().map(|x| x.0)
    }
    pub fn last_key(&self) -> Option<i32> {
        self.v.last().map(|x| x.0)
    }
    pub fn keys(&self) -> impl Iterator<Item = i32> + '_ {
        self.v.iter().map(|x| x.0)
    }
    pub fn retain(&mut self, mut f: impl FnMut(i32, &T) -> bool) {
        self.v.retain(|x| f(x.0, &x.1));
    }
    pub fn clear(&mut self) {
        self.v.clear();
    }
}

static NEXT_UID: AtomicU64 = AtomicU64::new(1);

fn next_uid() -> u64 {
    NEXT_UID.fetch_add(1, AtomicOrdering::Relaxed)
}

/// fragment.Fragment
#[derive(Clone, Debug)]
pub struct Fragment {
    /// object identity (Java compares fragments by reference in a few places)
    pub uid: u64,
    pub reads: Vec<Arc<Rec>>,
    pub read_gene: HlaGene,
    /// in insertion order (a HashSet of enum constants: insertion order under -XX:hashCode=2)
    pub genes: Vec<HlaGene>,
    pub raw: LociMap<Nuc>,
    pub nucs: LociMap<Nuc>,
    pub qual_filtered: bool,
    pub aa_conversions: u32,
    pub aas: LociMap<Tok>,
    pub scope: Scope,
}

impl Fragment {
    pub fn new(read: Arc<Rec>, read_gene: HlaGene, genes: Vec<HlaGene>, nucs: impl IntoIterator<Item = (i32, Nuc)>) -> Fragment {
        let m = LociMap::from_pairs(nucs);
        Fragment {
            uid: next_uid(),
            reads: vec![read],
            read_gene,
            genes,
            raw: m.clone(),
            nucs: m,
            qual_filtered: false,
            aa_conversions: 0,
            aas: LociMap::new(),
            scope: Scope::Unset,
        }
    }

    pub fn id(&self) -> &str {
        &self.reads[0].name
    }
    pub fn read_info(&self) -> String {
        self.reads[0].read_info()
    }

    pub fn add_read(&mut self, r: &Arc<Rec>) {
        if self.reads.iter().any(|x| x.flags == r.flags) {
            return;
        }
        self.reads.push(r.clone());
    }
    pub fn add_reads(&mut self, other: &Fragment) {
        for r in &other.reads {
            self.add_read(r);
        }
    }
    pub fn contains_gene(&self, g: HlaGene) -> bool {
        self.genes.contains(&g)
    }
    pub fn add_gene(&mut self, g: HlaGene) {
        if !self.genes.contains(&g) {
            self.genes.push(g);
        }
    }

    /// addNucleotide
    pub fn add_nucleotide(&mut self, locus: i32, n: Nuc) {
        if !self.nucs.contains(locus) {
            self.nucs.insert(locus, n);
            if !self.qual_filtered {
                self.raw.insert(locus, n);
                return;
            }
        }
        if !self.raw.contains(locus) {
            self.raw.insert(locus, n);
        }
    }

    pub fn has_nucleotides(&self) -> bool {
        !self.nucs.is_empty()
    }
    pub fn contains_indel(&self) -> bool {
        self.nucs.v.iter().any(|x| x.1.bases == DEL || tok::len(x.1.bases) > 1)
    }
    #[inline]
    pub fn contains_nucleotide(&self, locus: i32) -> bool {
        self.nucs.contains(locus)
    }
    pub fn contains_all_nucleotides(&self, loci: &[i32]) -> bool {
        loci.iter().all(|l| self.nucs.contains(*l))
    }
    /// nucleotide(locus): "" when missing
    pub fn nucleotide(&self, locus: i32) -> Tok {
        self.nucs.get(locus).map(|n| n.bases).unwrap_or(tok::EMPTY)
    }
    pub fn nucleotides_str(&self, loci: impl IntoIterator<Item = i32>) -> Vec<u8> {
        let mut s = Vec::new();
        for l in loci {
            tok::push_bytes(self.nucleotide(l), &mut s);
        }
        s
    }
    pub fn min_nucleotide_locus(&self) -> i32 {
        self.nucs.first_key().unwrap_or(-1)
    }
    pub fn max_nucleotide_locus(&self) -> i32 {
        self.nucs.last_key().unwrap_or(-1)
    }
    pub fn min_aa_locus(&self) -> i32 {
        self.aas.first_key().unwrap_or(-1)
    }
    pub fn max_aa_locus(&self) -> i32 {
        self.aas.last_key().unwrap_or(-1)
    }

    /// removeLowQualBases
    pub fn remove_low_qual_bases(&mut self) {
        if self.qual_filtered {
            return;
        }
        self.qual_filtered = true;
        self.nucs.retain(|_, n| !n.low_qual);
    }

    /// buildAminoAcids: an amino acid for each complete codon
    pub fn build_amino_acids(&mut self) {
        self.aa_conversions += 1;
        self.aas.clear();
        let mut out: Vec<(i32, Tok)> = Vec::new();
        let v = &self.nucs.v;
        for (i, &(locus, _)) in v.iter().enumerate() {
            if locus % 3 != 0 {
                continue;
            }
            // the next two loci, which for sorted unique keys are the next two entries
            if i + 2 < v.len() && v[i + 1].0 == locus + 1 && v[i + 2].0 == locus + 2 {
                let aa = codon_amino_acid(v[i].1.bases, v[i + 1].1.bases, v[i + 2].1.bases);
                out.push((locus / 3, aa));
            }
        }
        self.aas = LociMap { v: out };
    }

    pub fn contains_aa(&self, locus: i32) -> bool {
        self.aas.contains(locus)
    }
    pub fn contains_aa_loci(&self, loci: &[i32]) -> bool {
        loci.iter().all(|l| self.aas.contains(*l))
    }
    /// aminoAcid(locus): "" when missing
    pub fn aa(&self, locus: i32) -> Tok {
        self.aas.get(locus).unwrap_or(tok::EMPTY)
    }
    pub fn aas_str(&self, loci: &[i32]) -> String {
        let mut s = String::new();
        for l in loci {
            tok::push(self.aa(*l), &mut s);
        }
        s
    }
    pub fn filter_aa_on_loci(&mut self, keep: &[i32]) {
        self.aas.retain(|l, _| keep.binary_search(&l).is_ok());
    }

    /// getRawNucleotide(locus, checkUncertainBaseQual)
    pub fn raw_nucleotide(&self, locus: i32, check_uncertain: bool, params: &Params) -> Tok {
        let Some(n) = self.raw.get(locus) else { return tok::EMPTY };
        if check_uncertain && !params.is_ultima() && n.qual != MISSING_BASE_QUAL && is_uncertain_base_qual(n.qual) {
            return tok::EMPTY;
        }
        n.bases
    }

    /// getLowQualAminoAcid(locus, checkUncertainBaseQual): the amino acid of the raw bases
    pub fn low_qual_aa(&self, locus: i32, check_uncertain: bool, params: &Params) -> Tok {
        let s = locus * 3;
        let (Some(a), Some(b), Some(c)) = (self.raw.get(s), self.raw.get(s + 1), self.raw.get(s + 2)) else {
            return tok::EMPTY;
        };
        if check_uncertain && !params.is_ultima() {
            for n in [a, b, c] {
                if n.qual != MISSING_BASE_QUAL && is_uncertain_base_qual(n.qual) {
                    return tok::EMPTY;
                }
            }
        }
        codon_amino_acid(a.bases, b.bases, c.bases)
    }

    pub fn is_scope_set(&self) -> bool {
        self.scope != Scope::Unset
    }
    pub fn clear_scope(&mut self) {
        self.scope = Scope::Unset;
    }
    /// setScope(scope, override)
    pub fn set_scope(&mut self, s: Scope, over: bool) {
        if self.scope != Scope::Unset {
            if over {
                if s != Scope::HlaY {
                    crate::debug!("frag({}: {}) overriding existing scope: {} -> {}", self.id(), self.read_info(), self.scope.name(), s.name());
                }
                self.scope = s;
            }
        } else {
            self.scope = s;
        }
    }

    /// FragmentUtils.validateFragment
    pub fn validate(&self) -> bool {
        if self.genes.is_empty() {
            crate::warn!("{} {} has no genes", self.id(), self.read_info());
            return false;
        }
        if self.nucs.is_empty() || self.raw.is_empty() {
            crate::warn!("{} {} has no bases", self.id(), self.read_info());
            return false;
        }
        if self.aa_conversions > 0 {
            for (i, (_, a)) in self.aas.v.iter().enumerate() {
                if *a == tok::EMPTY {
                    crate::warn!("{} {} empty amino-acid, index({i})", self.id(), self.read_info());
                    return false;
                }
            }
            if self.aa_conversions > 1 {
                crate::warn!("{} {} amino-acid conversion repeated({})", self.id(), self.read_info(), self.aa_conversions);
                return false;
            }
        }
        true
    }

    /// FragmentUtils.copyNucleotideFragment: a new fragment from the raw bases and the reads
    pub fn copy_nucleotides(&self) -> Fragment {
        let mut f = Fragment {
            uid: next_uid(),
            reads: vec![self.reads[0].clone()],
            read_gene: self.read_gene,
            genes: self.genes.clone(),
            raw: self.raw.clone(),
            nucs: self.raw.clone(),
            qual_filtered: false,
            aa_conversions: 0,
            aas: LociMap::new(),
            scope: Scope::Unset,
        };
        for r in &self.reads[1..] {
            f.add_read(r);
        }
        f
    }
}

/// FragmentUtils.formCodonAminoAcid on three loci's bases: "." for three deletes
pub fn codon_amino_acid(a: Tok, b: Tok, c: Tok) -> Tok {
    if a == DEL && b == DEL && c == DEL {
        return DEL;
    }
    if tok::is_single(a) && tok::is_single(b) && tok::is_single(c) {
        return seq::codon_to_amino_acid(a as u8, b as u8, c as u8) as u8 as Tok;
    }
    let s = tok::concat3(a, b, c);
    tok::intern(&seq::amino_acids_from_bases(s.as_bytes()))
}

/// FragmentUtils.mergeFragments: frag2's genes and the nucleotides at loci frag1 lacks go into frag1
pub fn merge_fragments(frag1: &mut Fragment, frag2: &Fragment) {
    for g in &frag2.genes {
        frag1.add_gene(*g);
    }
    for &(locus, n) in &frag2.nucs.v {
        if frag1.contains_nucleotide(locus) {
            continue;
        }
        frag1.add_nucleotide(locus, n);
        frag1.add_reads(frag2);
    }
}

/// The order Java iterates a HashMap<String, ...> keyed by read id, built by puts in this order
pub fn read_id_map_order(ids: &[&str]) -> Vec<usize> {
    let cap = grown_capacity(16, ids.len());
    let mut v: Vec<(usize, usize)> = ids.iter().enumerate().map(|(i, s)| (bucket_index(java_string_hash(s), cap), i)).collect();
    v.sort_by_key(|x| x.0);
    v.into_iter().map(|x| x.1).collect()
}

/// FragmentUtils.mergeFragmentsById: fragments with one read id merged, in the HashMap's order of ids
pub fn merge_fragments_by_id(fragments: Vec<Fragment>) -> Vec<Fragment> {
    let mut groups: Vec<Vec<Fragment>> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    for f in fragments {
        let id = f.id().to_string();
        match index.get(&id) {
            Some(i) => groups[*i].push(f),
            None => {
                index.insert(id, groups.len());
                groups.push(vec![f]);
            }
        }
    }
    let ids: Vec<String> = groups.iter().map(|g| g[0].id().to_string()).collect();
    let order = read_id_map_order(&ids.iter().map(|s| s.as_str()).collect::<Vec<_>>());
    let mut groups: Vec<Option<Vec<Fragment>>> = groups.into_iter().map(Some).collect();
    let mut out = Vec::new();
    for i in order {
        let mut g = groups[i].take().unwrap().into_iter();
        let mut first = g.next().unwrap();
        for f in g {
            merge_fragments(&mut first, &f);
        }
        out.push(first);
    }
    out
}

/// FragmentUtils.calcAminoAcidIndices
pub fn calc_amino_acid_indices(nuc_start: i32, nuc_end: i32) -> Vec<i32> {
    let start = nuc_start / 3 + if nuc_start % 3 == 0 { 0 } else { 1 };
    let end = (nuc_end + 1) / 3 - 1;
    if end <= start {
        return vec![start, end];
    }
    (start..=end).collect()
}

/// FragmentUtils.expandIndices
pub fn expand_indices(aa_loci: &[i32]) -> Vec<i32> {
    aa_loci.iter().flat_map(|x| [3 * x, 3 * x + 1, 3 * x + 2]).collect()
}

/// SequenceCount.buildFromNucleotides
pub fn count_nucleotides<'a>(params: &Params, factor: f64, frags: impl IntoIterator<Item = &'a Fragment>) -> SequenceCount {
    let mut c = SequenceCount::new(factor, params.min_evidence_support);
    for f in frags {
        for &(l, n) in &f.nucs.v {
            c.add(l, n.bases);
        }
    }
    c
}

/// SequenceCount.buildFromAminoAcids
pub fn count_amino_acids<'a>(params: &Params, factor: f64, frags: impl IntoIterator<Item = &'a Fragment>) -> SequenceCount {
    let mut c = SequenceCount::new(factor, params.min_evidence_support);
    for f in frags {
        for &(l, a) in &f.aas.v {
            c.add(l, a);
        }
    }
    c
}

/// An amino acid sequence with inserts or deletes, as a string to search (SuffixTree.indices)
pub struct IndelSequence {
    pub seq: SeqLoci,
    pub text: Vec<u8>,
}

impl IndelSequence {
    /// SuffixTree(text).indices(pattern): the start of every occurrence, ascending
    pub fn indices(&self, pattern: &[u8]) -> Vec<i32> {
        let t = &self.text;
        let n = pattern.len();
        let mut out = Vec::new();
        if n == 0 || n > t.len() {
            return out;
        }
        for i in 0..=t.len() - n {
            if t[i] == pattern[0] && &t[i..i + n] == pattern {
                out.push(i as i32);
            }
        }
        out
    }
}

/// NucleotideFragmentFactory
pub struct Factory<'a> {
    pub gene_cache: &'a GeneCache,
    pub inserts: Vec<IndelSequence>,
    pub deletes: Vec<IndelSequence>,
    pub read_ctx: ReadCtx<'a>,
}

impl<'a> Factory<'a> {
    pub fn new(gene_cache: &'a GeneCache, inserts: &[SeqLoci], deletes: &[SeqLoci], read_ctx: ReadCtx<'a>) -> Factory<'a> {
        let mk = |s: &SeqLoci| IndelSequence { text: s.sequence_str().into_bytes(), seq: s.clone() };
        Factory { gene_cache, inserts: inserts.iter().map(mk).collect(), deletes: deletes.iter().map(mk).collect(), read_ctx }
    }

    /// createFragment: None where Java returns null, Err where it throws
    pub fn create_fragment(&self, read: &Read, gene: HlaGene, gene_strand: i8) -> Result<Option<Fragment>> {
        if read.read_index_start < 0 || read.read_index_end < read.read_index_start || read.position_end < read.position_start {
            return Ok(None);
        }
        let reverse = gene_strand == -1;
        let start_loci = self.gene_cache.nucleotide_locus(read.position_start);
        let end_loci = self.gene_cache.nucleotide_locus(read.position_end);
        let sam_start = if !reverse { start_loci } else { end_loci };
        let sam_end = if !reverse { end_loci } else { start_loci };
        let read_len = (read.read_index_end - read.read_index_start + 1) as usize;
        if read.read_index_end as usize >= read.record.read_len() {
            bail!("java.lang.StringIndexOutOfBoundsException: index {}, length {}", read.read_index_end, read.record.read_len());
        }
        let (bases, low, quals) = read.coding_region(self.read_ctx.params, reverse);

        if read.contains_valid_indel() || read.contains_soft_clip() {
            let aa_indices = calc_amino_acid_indices(sam_start, sam_end);
            let first_aa = aa_indices[0];
            let nuc_start = first_aa * 3;
            let start = nuc_start - sam_start;
            if start < 0 || start as usize >= bases.len() {
                return Ok(None);
            }
            let aas = seq::amino_acids_from_bases(&bases[start as usize..]);
            if !aas.is_empty() {
                let lo = first_aa - read.soft_clipped_start / 3 - read.max_valid_indel_size();
                let hi = first_aa + read.max_valid_indel_size() + read.soft_clipped_end / 3;
                if let Some(f) = self.check_matched_indel(read, gene, aas.as_bytes(), lo, hi, &self.inserts)? {
                    return Ok(Some(f));
                }
                if let Some(f) = self.check_matched_indel(read, gene, aas.as_bytes(), lo, hi, &self.deletes)? {
                    return Ok(Some(f));
                }
            }
            if read.contains_valid_indel() {
                return Ok(None);
            }
        }
        if sam_start < 0 || sam_end < 0 {
            return Ok(None);
        }
        let mut nucs: Vec<(i32, Nuc)> = Vec::with_capacity((sam_end - sam_start + 1).max(0) as usize);
        let mut queue: Vec<&crate::indel::Indel> = read.ignored_indels.iter().collect();
        // a PriorityQueue by read index: taken from the smallest
        queue.sort_by_key(|x| x.read_index);
        let mut q = 0usize;
        let inc: i32 = if reverse { -1 } else { 1 };
        let mut i: i64 = if reverse { read_len as i64 - 1 } else { 0 };
        let mut loci = if reverse { sam_end } else { sam_start };
        let mut read_index = read.read_index_start;
        while i >= 0 && (i as usize) < read_len {
            if q < queue.len() && queue[q].read_index < read_index {
                q += 1;
            }
            let mut current: Option<&crate::indel::Indel> = None;
            if q < queue.len() && queue[q].read_index == read_index {
                current = Some(queue[q]);
                q += 1;
            }
            let iu = i as usize;
            let nuc = Nuc { bases: tok::of_byte(bases[iu]), low_qual: low[iu], qual: quals[iu] };
            match current {
                None => {
                    nucs.push((loci, nuc));
                    loci += inc;
                    read_index += 1;
                    i += inc as i64;
                }
                Some(ind) if ind.is_insert => {
                    nucs.push((loci, nuc));
                    loci += inc;
                    read_index += 1 + ind.length;
                    i += inc as i64 * (1 + ind.length) as i64;
                }
                Some(ind) => {
                    nucs.push((loci, nuc));
                    loci += inc;
                    read_index += 1;
                    i += inc as i64;
                    let rb = ind.ref_.as_bytes();
                    for &b in rb.iter().skip(1) {
                        if loci < sam_start || loci > sam_end {
                            break;
                        }
                        let b = if reverse { crate::bam::swap_dna_base(b) } else { b };
                        nucs.push((loci, Nuc::high_qual(tok::of_byte(b))));
                        loci += inc;
                    }
                }
            }
        }
        if reverse {
            nucs.reverse();
        }
        Ok(Some(Fragment::new(read.record.clone(), gene, vec![gene], nucs)))
    }

    fn check_matched_indel(&self, read: &Read, gene: HlaGene, aas: &[u8], lo: i32, hi: i32, seqs: &[IndelSequence]) -> Result<Option<Fragment>> {
        let mut matched: Vec<(usize, Vec<i32>)> = Vec::new();
        for (k, s) in seqs.iter().enumerate() {
            let idx: Vec<i32> = s.indices(aas).into_iter().filter(|x| *x >= lo && *x <= hi).collect();
            if !idx.is_empty() {
                matched.push((k, idx));
            }
        }
        for (k, idx) in matched {
            let f = create_indel_fragment(read, gene, idx[0], aas, &seqs[k].seq)?;
            if !f.nucs.is_empty() {
                return Ok(Some(f));
            }
        }
        Ok(None)
    }

    /// createAlignmentFragments: one fragment per overlapping alignment block, all built from the whole read
    /// and merged (so the blocks only decide whether there is a fragment)
    pub fn create_alignment_fragments(&self, read: &Read, gene: HlaGene, strand: i8) -> Result<Option<Fragment>> {
        let blocks = read.alignments_only(&self.read_ctx)?;
        let mut frags = Vec::new();
        for _ in &blocks {
            if let Some(f) = self.create_fragment(read, gene, strand)? {
                frags.push(f);
            }
        }
        if frags.is_empty() {
            return Ok(None);
        }
        if frags.len() == 1 {
            return Ok(frags.pop());
        }
        Ok(merge_fragments_by_id(frags).into_iter().next())
    }
}

/// createIndelFragment: the nucleotides of the indel sequence's codons from its matching start
fn create_indel_fragment(read: &Read, gene: HlaGene, start: i32, bam_seq: &[u8], s: &SeqLoci) -> Result<Fragment> {
    // endLoci: the last locus whose sequence still prefixes the read's amino acids
    let mut end = s.len() as i32 - 1;
    let mut tmp: Vec<u8> = Vec::new();
    for locus in start..s.len() as i32 {
        tok::push_bytes(s.at(locus as usize), &mut tmp);
        if !bam_seq.starts_with(&tmp) {
            end = locus - 1;
            break;
        }
    }
    let aa_loci: Vec<i32> = if end >= start { (start..=end).collect() } else { Vec::new() };
    let nuc_loci = expand_indices(&aa_loci);
    let mut nucs: Vec<Tok> = Vec::with_capacity(nuc_loci.len());
    for &l in &aa_loci {
        nucs.extend(nucleotides_from_amino_acid(s.at(l as usize))?);
    }
    let pairs: Vec<(i32, Nuc)> = nuc_loci.into_iter().zip(nucs).map(|(l, b)| (l, Nuc { bases: b, low_qual: false, qual: MISSING_BASE_QUAL })).collect();
    Ok(Fragment::new(read.record.clone(), gene, vec![gene], pairs))
}

/// NucleotideFragmentFactory.createNucleotidesFromAminoAcid: the first codon of each amino acid, as three
/// loci (the third holding the rest for an insert)
pub fn nucleotides_from_amino_acid(aa: Tok) -> Result<[Tok; 3]> {
    if aa == DEL {
        return Ok([DEL, DEL, DEL]);
    }
    let aa_str = tok::to_string(aa);
    let mut codons = Vec::new();
    for c in aa_str.bytes() {
        match seq::amino_acid_to_codon(c) {
            Some(x) => codons.extend_from_slice(&x),
            None => bail!("java.lang.IllegalArgumentException: Unknown amino acid {}", c as char),
        }
    }
    Ok([tok::of_byte(codons[0]), tok::of_byte(codons[1]), tok::intern(std::str::from_utf8(&codons[2..]).unwrap())])
}

/// NucleotideFragmentFactory.calculatePercentileBaseQuality
pub fn percentile_base_quality<'a>(frags: impl IntoIterator<Item = &'a Fragment>, percentile: f64) -> i8 {
    let max_qual = LOW_BASE_QUAL_THRESHOLD as usize * 2;
    let mut freq = vec![0u64; max_qual + 1];
    let mut total: i64 = 0;
    for f in frags {
        for (_, n) in &f.raw.v {
            if n.qual == MISSING_BASE_QUAL {
                continue;
            }
            total += 1;
            // min(baseQual, maxBaseQual) on a signed byte
            let q = (n.qual as i32).min(max_qual as i32);
            if q < 0 {
                // Java indexes the array with a negative value and throws
                panic!("java.lang.ArrayIndexOutOfBoundsException: Index {q} out of bounds for length {}", max_qual + 1);
            }
            freq[q as usize] += 1;
        }
    }
    let entry = hmf_common::java_math::java_round(total as f64 * percentile);
    let mut cum: i64 = 0;
    for (i, c) in freq.iter().enumerate() {
        cum += *c as i64;
        if cum >= entry {
            return i as i8;
        }
    }
    LOW_BASE_QUAL_THRESHOLD as i8
}

/// NucleotideFragmentFactory.calculateGeneCoverage: per gene, the raw depth at each coding base
pub fn gene_coverage<'a>(gc: &GeneCache, frags: impl IntoIterator<Item = &'a Fragment>) -> Vec<(HlaGene, Vec<i32>)> {
    let mut out: Vec<(HlaGene, Vec<i32>)> = gc.gene_names.iter().map(|g| (*g, vec![0; gc.nuc_length(*g).max(0) as usize])).collect();
    for f in frags {
        let Some(d) = out.iter_mut().find(|x| x.0 == f.read_gene) else {
            panic!("java.lang.NullPointerException (no depth for gene {})", f.read_gene.long_name());
        };
        for l in f.raw.keys() {
            if l >= 0 && (l as usize) < d.1.len() {
                d.1[l as usize] += 1;
            }
        }
    }
    out
}

/// NucleotideGeneEnrichment.checkAdditionalGenes: add the genes a fragment cannot be told apart from
pub fn check_additional_genes(e: &GeneEnrichment, gc: &GeneCache, f: &mut Fragment) {
    if f.contains_indel() {
        return;
    }
    let max_locus = f.max_nucleotide_locus();
    for &gene in &e.genes {
        if f.contains_gene(gene) || max_locus > gc.nuc_length(gene) {
            continue;
        }
        let mut add = false;
        for &other in &e.genes {
            if gene == other {
                continue;
            }
            if !f.contains_gene(other) {
                continue;
            }
            let b = e.filter(gene, other);
            if b.is_none() || max_locus < b.unwrap() * 3 {
                add = true;
                break;
            }
        }
        if add {
            f.add_gene(gene);
        }
    }
}

/// NucleotideFragmentQualEnrichment.applyQualityFilter: the fragment itself if all its nucleotides pass,
/// else a new fragment with those that do
pub fn apply_quality_filter(f: &Fragment, high: &SequenceCount, raw: &SequenceCount, raw_nuc: Option<&SequenceCount>, params: &Params) -> Option<Fragment> {
    let mut all = true;
    let mut kept: Vec<(i32, Nuc)> = Vec::with_capacity(f.nucs.len());
    for &(locus, n) in &f.nucs.v {
        if let Some(r) = raw_nuc {
            if (r.get(locus).total as i64) < params.min_depth_filter as i64 {
                kept.push((locus, n));
                continue;
            }
        }
        if high.has_min_evidence(locus, n.bases) && raw.has_min_evidence(locus, n.bases) {
            kept.push((locus, n));
        } else {
            all = false;
        }
    }
    if all {
        return None;
    }
    let mut nf = Fragment::new(f.reads[0].clone(), f.read_gene, f.genes.clone(), kept);
    nf.add_reads(f);
    Some(nf)
}

/// NucleotideSpliceEnrichment: the bases to add at exon boundaries where the high-qual reference is homozygous
pub struct SpliceInfo {
    counts: SequenceCount,
    hom_starts: Vec<i32>,
    hom_ends: Vec<i32>,
}

impl SpliceInfo {
    pub fn new(high_qual: &[Fragment], boundaries: &[i32], params: &Params) -> SpliceInfo {
        let counts = count_nucleotides(params, params.min_evidence_factor, high_qual.iter());
        let hom: Vec<i32> = counts.homozygous_loci();
        let is_hom = |l: i32| hom.binary_search(&l).is_ok();
        // the sets' iteration order does not matter: each fragment gets every addition
        let starts: Vec<i32> = boundaries.iter().map(|x| x * 3).collect();
        let hom_starts: Vec<i32> = starts.iter().copied().filter(|x| is_hom(*x)).collect();
        let hom_ends: Vec<i32> = starts.iter().copied().filter(|x| is_hom(x + 1) && is_hom(x + 2)).collect();
        SpliceInfo { counts, hom_starts, hom_ends }
    }

    /// applySpliceInfo on one fragment, in place
    pub fn apply(&self, f: &mut Fragment) {
        for &s in &self.hom_starts {
            if !f.contains_nucleotide(s) && f.contains_nucleotide(s + 1) && f.contains_nucleotide(s + 2) {
                let b = self.counts.min_evidence_sequences(s)[0];
                f.add_nucleotide(s, Nuc::high_qual(b));
            }
        }
        for &e in &self.hom_ends {
            if f.contains_nucleotide(e) && !f.contains_nucleotide(e + 1) && !f.contains_nucleotide(e + 2) {
                let b1 = self.counts.min_evidence_sequences(e + 1)[0];
                f.add_nucleotide(e + 1, Nuc::high_qual(b1));
                let b2 = self.counts.min_evidence_sequences(e + 2)[0];
                f.add_nucleotide(e + 2, Nuc::high_qual(b2));
            }
        }
    }
}

/// AminoAcidQualEnrichment.applyQualFilter: keep the amino acids with enough support (or low raw depth)
pub fn apply_aa_qual_filter(f: &mut Fragment, count: &SequenceCount, raw_aa: Option<&SequenceCount>, params: &Params) {
    let mut keep: Vec<i32> = Vec::new();
    for &(locus, a) in &f.aas.v {
        if let Some(r) = raw_aa {
            if (r.get(locus).total as i64) < params.min_depth_filter as i64 {
                keep.push(locus);
                continue;
            }
        }
        if count.has_min_evidence(locus, a) {
            keep.push(locus);
        }
    }
    f.filter_aa_on_loci(&keep);
}

/// SequenceCountDiff
#[derive(Clone, Debug)]
pub struct SequenceCountDiff {
    pub locus: i32,
    pub seq: Tok,
    pub ref_count: u32,
    pub ref_depth: u32,
    pub tumor_count: u32,
    pub tumor_depth: u32,
}

/// SequenceCountDiff.create: sequences seen in one sample only
pub fn sequence_count_diffs(r: &SequenceCount, t: &SequenceCount) -> Vec<SequenceCountDiff> {
    let mut out = Vec::new();
    let n = r.loci.len().max(t.loci.len());
    for l in 0..n as i32 {
        let in_ref = r.has(l);
        let in_tum = t.has(l);
        if !in_ref && !in_tum {
            continue;
        }
        let rc = r.get(l);
        let tc = t.get(l);
        if in_ref {
            for (s, c) in rc.java_order() {
                if !tc.contains(s) {
                    out.push(SequenceCountDiff { locus: l, seq: s, ref_count: c, ref_depth: rc.total, tumor_count: 0, tumor_depth: 0 });
                }
            }
        }
        if in_tum {
            for (s, c) in tc.java_order() {
                if !rc.contains(s) {
                    out.push(SequenceCountDiff { locus: l, seq: s, ref_count: 0, ref_depth: 0, tumor_count: c, tumor_depth: tc.total });
                }
            }
        }
    }
    out
}

/// AminoAcidFragmentPipeline
pub struct AaPipeline {
    /// copies of the reference fragments, for phasing
    pub original: Vec<Fragment>,
    /// mHighQualRefAminoAcidFragments: indices into the reference fragments
    pub high_qual_ref: Vec<usize>,
    /// per gene, in the order the genes were set
    pub ref_nuc_counts: Vec<(HlaGene, SequenceCount)>,
    pub ref_aa_counts: Vec<(HlaGene, SequenceCount)>,
}

/// createHighQualAminoAcidFragments, in place: drops low-qual bases, builds amino acids; returns the
/// indices of the fragments with nucleotides left (the others get BASE_QUAL_FILTERED)
pub fn high_qual_amino_acid_fragments(frags: &mut [Fragment], subset: impl IntoIterator<Item = usize>) -> Vec<usize> {
    let mut out = Vec::new();
    for i in subset {
        let f = &mut frags[i];
        f.remove_low_qual_bases();
        if !f.has_nucleotides() {
            f.set_scope(Scope::BaseQualFiltered, false);
            continue;
        }
        f.build_amino_acids();
        out.push(i);
    }
    out
}

/// The raw per-gene counts of the reference (AminoAcidFragmentPipeline.RAW_REF_NUCLEOTIDE_COUNTS and
/// RAW_REF_AMINO_ACID_COUNTS, static in Java)
#[derive(Default)]
pub struct RawCounts {
    pub nuc: HashMap<HlaGene, SequenceCount>,
    pub aa: HashMap<HlaGene, SequenceCount>,
}

impl AaPipeline {
    /// new AminoAcidFragmentPipeline(referenceFragments)
    pub fn new(ref_frags: &mut [Fragment]) -> AaPipeline {
        let original: Vec<Fragment> = ref_frags.iter().map(|f| f.copy_nucleotides()).collect();
        let n = ref_frags.len();
        let high_qual_ref = high_qual_amino_acid_fragments(ref_frags, 0..n);
        AaPipeline { original, high_qual_ref, ref_nuc_counts: Vec::new(), ref_aa_counts: Vec::new() }
    }

    /// referencePhasingFragments(context): the fragments of a gene filtered and enriched for phasing
    pub fn reference_phasing_fragments(&mut self, ctx: &HlaContext, gc: &GeneCache, raw: &mut RawCounts, params: &Params, copy_originals: bool) -> Vec<Fragment> {
        let gene = ctx.gene;
        if gene != HlaGene::Drb3 && gene != HlaGene::Drb4 {
            let raw_gene: Vec<&Fragment> = self.original.iter().filter(|f| f.read_gene == gene).collect();
            let raw_nuc = count_nucleotides(params, params.min_evidence_factor, raw_gene.iter().copied());
            let mut acid_frags: Vec<Fragment> = raw_gene.iter().map(|f| f.copy_nucleotides()).collect();
            for f in acid_frags.iter_mut() {
                f.build_amino_acids();
            }
            let raw_aa = count_amino_acids(params, params.min_evidence_factor, acid_frags.iter());
            raw.nuc.insert(gene, raw_nuc);
            raw.aa.insert(gene, raw_aa);
        }
        let gene_idx: Vec<usize> = (0..self.original.len()).filter(|i| self.original[*i].contains_gene(gene)).collect();
        if gene_idx.is_empty() {
            return Vec::new();
        }
        let mut high_qual: Vec<Fragment> = gene_idx.iter().map(|i| self.original[*i].copy_nucleotides()).collect();
        for f in high_qual.iter_mut() {
            f.remove_low_qual_bases();
        }
        high_qual.retain(|f| f.has_nucleotides());

        // NucleotideFragmentQualEnrichment.qualityFilterFragments
        let high_counts = count_nucleotides(params, params.min_high_qual_evidence_factor, high_qual.iter());
        let raw_counts = count_nucleotides(params, params.min_evidence_factor, gene_idx.iter().map(|i| &self.original[*i]));
        let raw_nuc = raw.nuc.get(&gene);
        // each fragment that passes in full is the original object (Java then enriches it in place, LI3);
        // others are new fragments
        let filtered: Vec<Option<Fragment>> = gene_idx.iter().map(|i| apply_quality_filter(&self.original[*i], &high_counts, &raw_counts, raw_nuc, params)).collect();
        let max_common = gc.max_common_aa_boundary;
        let mut boundaries: Vec<i32> = ctx.aa_boundaries.iter().copied().filter(|x| *x <= max_common).collect();
        boundaries.sort();
        boundaries.dedup();
        let splice = SpliceInfo::new(&high_qual, &boundaries, params);
        let mut enriched: Vec<Fragment> = Vec::with_capacity(gene_idx.len());
        for (f, i) in filtered.into_iter().zip(gene_idx.iter()) {
            match f {
                Some(mut nf) => {
                    splice.apply(&mut nf);
                    enriched.push(nf);
                }
                None if copy_originals => {
                    let mut c = self.original[*i].clone();
                    c.uid = next_uid();
                    splice.apply(&mut c);
                    enriched.push(c);
                }
                None => {
                    // LI3: Java enriches the original object, which later genes' phasing then sees
                    splice.apply(&mut self.original[*i]);
                    enriched.push(self.original[*i].clone());
                }
            }
        }

        // AminoAcidQualEnrichment.qualityFilterAminoAcidFragments
        let mut qf: Vec<Fragment> = enriched.iter().map(|f| f.copy_nucleotides()).collect();
        for f in qf.iter_mut() {
            f.build_amino_acids();
        }
        let hq_aa = count_amino_acids(params, params.min_evidence_factor, qf.iter());
        let mut out: Vec<Fragment> = enriched.iter().filter(|f| f.has_nucleotides()).map(|f| f.copy_nucleotides()).collect();
        let raw_aa = raw.aa.get(&gene);
        for f in out.iter_mut() {
            f.build_amino_acids();
            apply_aa_qual_filter(f, &hq_aa, raw_aa, params);
        }
        let nuc_counts = count_nucleotides(params, params.min_evidence_factor, out.iter());
        let aa_counts = count_amino_acids(params, params.min_evidence_factor, out.iter());
        self.set_counts(gene, nuc_counts, aa_counts);
        out
    }

    fn set_counts(&mut self, gene: HlaGene, n: SequenceCount, a: SequenceCount) {
        match self.ref_nuc_counts.iter_mut().find(|x| x.0 == gene) {
            Some(x) => x.1 = n,
            None => self.ref_nuc_counts.push((gene, n)),
        }
        match self.ref_aa_counts.iter_mut().find(|x| x.0 == gene) {
            Some(x) => x.1 = a,
            None => self.ref_aa_counts.push((gene, a)),
        }
    }

    /// getReferenceNucleotides: the union over genes of the nucleotides counted at each locus, as a list by
    /// position among the loci with counts (not by locus: LI4)
    pub fn reference_nucleotides(&self) -> Vec<Vec<Tok>> {
        let mut m: std::collections::BTreeMap<i32, Vec<Tok>> = std::collections::BTreeMap::new();
        for (_, c) in &self.ref_nuc_counts {
            for l in c.loci_keys() {
                let e = m.entry(l).or_default();
                for (t, _) in &c.get(l).entries {
                    if !e.contains(t) {
                        e.push(*t);
                    }
                }
            }
        }
        let n = m.keys().last().map(|x| *x + 1).unwrap_or(0) as usize;
        if n != m.len() {
            crate::debug!("LI4: reference nucleotide loci: {} with counts of {} (by position, Java shifts {} loci)", m.len(), n, n - m.len());
        }
        if crate::compat::java("LI4") {
            m.into_values().collect()
        } else {
            let n = m.keys().last().map(|x| *x + 1).unwrap_or(0) as usize;
            let mut v = vec![Vec::new(); n];
            for (l, s) in m {
                v[l as usize] = s;
            }
            v
        }
    }

    /// calcComparisonCoverageFragments: the tumor (or RNA) fragments after quality filtering, without those
    /// holding a nucleotide or amino acid the reference lacks
    pub fn comparison_coverage_fragments(&self, ref_frags: &[Fragment], frags: &mut [Fragment], params: &Params) -> Vec<usize> {
        let n = frags.len();
        let high = high_qual_amino_acid_fragments(frags, 0..n);
        if high.is_empty() {
            return Vec::new();
        }
        let ref_nuc = count_nucleotides(params, params.min_evidence_factor, self.high_qual_ref.iter().map(|i| &ref_frags[*i]));
        let ref_aa = count_amino_acids(params, params.min_evidence_factor, self.high_qual_ref.iter().map(|i| &ref_frags[*i]));
        let t_nuc = count_nucleotides(params, params.min_evidence_factor, high.iter().map(|i| &frags[*i]));
        let t_aa = count_amino_acids(params, params.min_evidence_factor, high.iter().map(|i| &frags[*i]));
        let nd: Vec<SequenceCountDiff> = sequence_count_diffs(&ref_nuc, &t_nuc).into_iter().filter(|x| x.tumor_count > 0).collect();
        let ad: Vec<SequenceCountDiff> = sequence_count_diffs(&ref_aa, &t_aa).into_iter().filter(|x| x.tumor_count > 0).collect();
        high.into_iter()
            .filter(|i| {
                let f = &frags[*i];
                let nv = nd.iter().any(|d| f.contains_nucleotide(d.locus) && f.nucleotide(d.locus) == d.seq);
                let av = ad.iter().any(|d| f.contains_aa(d.locus) && f.aa(d.locus) == d.seq);
                !(nv || av)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn amino_acid_indices() {
        assert_eq!(calc_amino_acid_indices(0, 8), vec![0, 1, 2]);
        assert_eq!(calc_amino_acid_indices(1, 8), vec![1, 2]);
        assert_eq!(calc_amino_acid_indices(1, 7), vec![1, 1]);
        assert_eq!(calc_amino_acid_indices(-1, 5), vec![1, 1]);
    }

    #[test]
    fn suffix_search() {
        let s = IndelSequence { seq: SeqLoci { allele: 0, seq: vec![], has_deletes: false, has_inserts: false, wildcard_count: 0, has_exon_boundary_wildcards: false }, text: b"ABCABCAB".to_vec() };
        assert_eq!(s.indices(b"AB"), vec![0, 3, 6]);
        assert_eq!(s.indices(b"CABX"), Vec::<i32>::new());
        assert_eq!(s.indices(b"ABCABCAB"), vec![0]);
    }
}

#[cfg(test)]
mod suffix_tree_tests {
    use super::*;

    /// hmf-common's SuffixTree, ported as written, to check that its indices are every occurrence
    struct Node {
        text: Vec<u8>,
        position: i32,
        children: Vec<usize>,
    }
    struct Tree {
        nodes: Vec<Node>,
    }
    impl Tree {
        fn new(text: &[u8]) -> Tree {
            let mut t = Tree { nodes: vec![Node { text: Vec::new(), position: -1, children: Vec::new() }] };
            for i in 0..text.len() {
                let mut s = text[i..].to_vec();
                s.push(b'$');
                t.add_suffix(&s, i as i32);
            }
            t
        }
        fn path(&self, pattern: &[u8], start: usize, partial: bool) -> Vec<Option<usize>> {
            let mut nodes = Vec::new();
            for &c in &self.nodes[start].children {
                let nt = &self.nodes[c].text;
                if pattern[0] == nt[0] {
                    if partial && pattern.len() <= nt.len() {
                        nodes.push(Some(c));
                        return nodes;
                    }
                    let cl = nt.len().min(pattern.len());
                    for j in 1..cl {
                        if pattern[j] != nt[j] {
                            if partial {
                                nodes.push(Some(c));
                            }
                            return nodes;
                        }
                    }
                    nodes.push(Some(c));
                    if pattern.len() > cl {
                        let n2 = self.path(&pattern[cl..], c, partial);
                        if !n2.is_empty() {
                            nodes.extend(n2);
                        } else if !partial {
                            nodes.push(None);
                        }
                    }
                    return nodes;
                }
            }
            nodes
        }
        fn add_suffix(&mut self, suffix: &[u8], pos: i32) {
            let mut nodes: Vec<usize> = self.path(suffix, 0, true).into_iter().map(|x| x.unwrap()).collect();
            if nodes.is_empty() {
                self.nodes.push(Node { text: suffix.to_vec(), position: pos, children: Vec::new() });
                let id = self.nodes.len() - 1;
                self.nodes[0].children.push(id);
                return;
            }
            let last = nodes.pop().unwrap();
            let prefix_len: usize = nodes.iter().map(|n| self.nodes[*n].text.len()).sum();
            let new_text = suffix[prefix_len..].to_vec();
            let cur = self.nodes[last].text.clone();
            let mut common = 0;
            while common < cur.len().min(new_text.len()) && cur[common] == new_text[common] {
                common += 1;
            }
            if common != cur.len() {
                let child = Node { text: cur[common..].to_vec(), position: self.nodes[last].position, children: std::mem::take(&mut self.nodes[last].children) };
                self.nodes.push(child);
                let cid = self.nodes.len() - 1;
                self.nodes[last].children.push(cid);
                self.nodes[last].text = cur[..common].to_vec();
                self.nodes[last].position = -1;
            }
            self.nodes.push(Node { text: new_text[common..].to_vec(), position: pos, children: Vec::new() });
            let id = self.nodes.len() - 1;
            self.nodes[last].children.push(id);
        }
        fn positions(&self, n: usize, out: &mut Vec<i32>) {
            if self.nodes[n].text.ends_with(b"$") {
                out.push(self.nodes[n].position);
            }
            for &c in &self.nodes[n].children {
                self.positions(c, out);
            }
        }
        fn indices(&self, pattern: &[u8]) -> Vec<i32> {
            let p = self.path(pattern, 0, false);
            match p.last() {
                Some(Some(n)) => {
                    let mut v = Vec::new();
                    self.positions(*n, &mut v);
                    v.sort();
                    v
                }
                _ => Vec::new(),
            }
        }
    }

    #[test]
    fn suffix_tree_indices_are_all_occurrences() {
        let mut x: u64 = 7;
        let mut rnd = |n: u64| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (x >> 33) % n
        };
        for _ in 0..300 {
            let len = 5 + rnd(60) as usize;
            let text: Vec<u8> = (0..len).map(|_| b"ACDE"[rnd(4) as usize]).collect();
            let tree = Tree::new(&text);
            let s = IndelSequence { seq: SeqLoci { allele: 0, seq: vec![], has_deletes: false, has_inserts: false, wildcard_count: 0, has_exon_boundary_wildcards: false }, text: text.clone() };
            for _ in 0..20 {
                let plen = 1 + rnd(8) as usize;
                let pattern: Vec<u8> = if rnd(2) == 0 && plen < len {
                    let st = rnd((len - plen) as u64) as usize;
                    text[st..st + plen].to_vec()
                } else {
                    (0..plen).map(|_| b"ACDE"[rnd(4) as usize]).collect()
                };
                assert_eq!(tree.indices(&pattern), s.indices(&pattern), "{:?} in {:?}", String::from_utf8_lossy(&pattern), String::from_utf8_lossy(&text));
            }
        }
    }
}
