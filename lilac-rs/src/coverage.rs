//! FragmentAlleles, AlleleCoverage, ComplexCoverage, FragmentAlleleMapper and ComplexBuilder.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};

use hmf_common::java_hash::{bucket_index, grown_capacity, jdk19_collection_capacity};

use crate::config::Params;
use crate::constants::*;
use crate::fragment::{Fragment, Scope};
use crate::hla::{AId, Alleles, HlaGene};
use crate::refdata::{GeneCache, RefData};
use crate::seq::{MatchType, SeqLoci};
use crate::tok::{Tok, WILD};

/// The Java iteration order of a HashSet / HashMap of alleles with `cap` buckets, the keys given in
/// first-insertion order
pub fn allele_set_order(alleles: &Alleles, keys: &[AId], cap: usize) -> Vec<AId> {
    let mut v: Vec<(usize, AId)> = keys.iter().map(|a| (bucket_index(alleles.hash(*a), cap), *a)).collect();
    v.sort_by_key(|x| x.0);
    v.into_iter().map(|x| x.1).collect()
}

fn dedup(v: &[AId]) -> Vec<AId> {
    let mut out: Vec<AId> = Vec::with_capacity(v.len());
    for a in v {
        if !out.contains(a) {
            out.push(*a);
        }
    }
    out
}

/// FragmentAlleles: a fragment's full and wildcard matches (HashSets: kept in Java's iteration order)
#[derive(Clone, Debug)]
pub struct FragAlleles {
    /// index into the fragment list the alleles were built from
    pub fi: usize,
    pub full: Vec<AId>,
    pub wild: Vec<AId>,
    full_sorted: Vec<AId>,
    wild_sorted: Vec<AId>,
}

impl FragAlleles {
    /// new FragmentAlleles(fragment, full, wild): Sets.newHashSet(list) of each
    pub fn new(alleles: &Alleles, fi: usize, full: &[AId], wild: &[AId]) -> FragAlleles {
        let f = dedup(full);
        let w = dedup(wild);
        let full = allele_set_order(alleles, &f, jdk19_collection_capacity(full.len()));
        let wild = allele_set_order(alleles, &w, jdk19_collection_capacity(wild.len()));
        let mut full_sorted = full.clone();
        full_sorted.sort();
        let mut wild_sorted = wild.clone();
        wild_sorted.sort();
        FragAlleles { fi, full, wild, full_sorted, wild_sorted }
    }
    #[inline]
    pub fn has_full(&self, a: AId) -> bool {
        self.full_sorted.binary_search(&a).is_ok()
    }
    #[inline]
    pub fn has_wild(&self, a: AId) -> bool {
        self.wild_sorted.binary_search(&a).is_ok()
    }
    #[inline]
    pub fn contains(&self, a: AId) -> bool {
        self.has_full(a) || self.has_wild(a)
    }
    /// removeAll from both sets (the order of the rest is kept)
    pub fn remove_all(&mut self, remove: &[AId]) {
        self.full.retain(|a| !remove.contains(a));
        self.wild.retain(|a| !remove.contains(a));
        self.full_sorted.retain(|a| !remove.contains(a));
        self.wild_sorted.retain(|a| !remove.contains(a));
    }
}

/// FragmentAlleles.filter: the fragments with any of the alleles, restricted to them
pub fn filter_frag_alleles(alleles: &Alleles, fas: &[FragAlleles], keep: &[AId]) -> Vec<FragAlleles> {
    let mut out = Vec::new();
    for fa in fas {
        if keep.iter().any(|a| fa.contains(*a)) {
            let full: Vec<AId> = keep.iter().copied().filter(|a| fa.has_full(*a)).collect();
            let wild: Vec<AId> = keep.iter().copied().filter(|a| fa.has_wild(*a)).collect();
            out.push(FragAlleles::new(alleles, fa.fi, &full, &wild));
        }
    }
    out
}

/// AlleleCoverage
#[derive(Clone, Debug)]
pub struct AlleleCoverage {
    pub allele: AId,
    pub unique: i32,
    pub shared: f64,
    pub wild: f64,
    pub total: f64,
}

impl AlleleCoverage {
    pub fn new(allele: AId, unique: i32, shared: f64, wild: f64) -> AlleleCoverage {
        AlleleCoverage { allele, unique, shared, wild, total: unique as f64 + shared + wild }
    }
    /// compareTo: by unique, then shared coverage
    pub fn compare(&self, o: &AlleleCoverage) -> Ordering {
        if self.unique != o.unique {
            return if self.unique > o.unique { Ordering::Greater } else { Ordering::Less };
        }
        if self.shared != o.shared {
            return if self.shared > o.shared { Ordering::Greater } else { Ordering::Less };
        }
        Ordering::Equal
    }
    /// toString: A*01:01[total,unique,shared,wild]
    pub fn display(&self, alleles: &Alleles) -> String {
        format!(
            "{}[{},{},{},{}]",
            alleles.name(self.allele),
            hmf_common::fmt::java_fixed(self.total, 0),
            self.unique,
            hmf_common::fmt::java_fixed(self.shared, 0),
            hmf_common::fmt::java_fixed(self.wild, 0)
        )
    }
}

pub fn coverage_list_str(alleles: &Alleles, v: &[AlleleCoverage]) -> String {
    v.iter().map(|x| x.display(alleles)).collect::<Vec<_>>().join(", ")
}

/// An insertion-ordered map of alleles to values, iterated in Java HashMap order
struct AlleleMap<V> {
    keys: Vec<AId>,
    vals: HashMap<AId, V>,
}

impl<V: Copy> AlleleMap<V> {
    fn new() -> Self {
        AlleleMap { keys: Vec::new(), vals: HashMap::new() }
    }
    fn get(&self, a: AId) -> Option<V> {
        self.vals.get(&a).copied()
    }
    fn update(&mut self, a: AId, init: V, f: impl FnOnce(V) -> V) {
        match self.vals.get_mut(&a) {
            Some(v) => *v = f(*v),
            None => {
                self.keys.push(a);
                self.vals.insert(a, init);
            }
        }
    }
    fn java_keys(&self, alleles: &Alleles) -> Vec<AId> {
        allele_set_order(alleles, &self.keys, grown_capacity(16, self.keys.len()))
    }
}

/// AlleleCoverage.buildCoverage (proteinCoverage, or groupCoverage with `as_group`)
pub fn build_coverage(alleles: &Alleles, fas: &[FragAlleles], as_group: bool) -> Vec<AlleleCoverage> {
    let mut unique: AlleleMap<i32> = AlleleMap::new();
    let mut combined: AlleleMap<f64> = AlleleMap::new();
    let mut wildm: AlleleMap<f64> = AlleleMap::new();
    for fa in fas {
        let (full, wild): (Vec<AId>, Vec<AId>) = if as_group {
            // Collectors.toCollection(Sets::newHashSet): a default HashSet filled in the source's order
            let g = |v: &[AId]| -> Vec<AId> {
                let d = dedup(&v.iter().map(|a| alleles.group(*a)).collect::<Vec<_>>());
                allele_set_order(alleles, &d, grown_capacity(16, d.len()))
            };
            (g(&fa.full), g(&fa.wild))
        } else {
            (fa.full.clone(), fa.wild.clone())
        };
        if full.len() == 1 && wild.is_empty() {
            unique.update(full[0], 1, |v| v + 1);
        } else {
            let c = 1.0 / (full.len() + wild.len()) as f64;
            for a in &full {
                combined.update(*a, c, |v| v + c);
            }
            for a in &wild {
                wildm.update(*a, c, |v| v + c);
            }
        }
    }
    let mut list: Vec<AId> = unique.java_keys(alleles);
    for a in combined.java_keys(alleles) {
        if !list.contains(&a) {
            list.push(a);
        }
    }
    let mut out: Vec<AlleleCoverage> = list
        .into_iter()
        .map(|a| AlleleCoverage::new(a, unique.get(a).unwrap_or(0), combined.get(a).unwrap_or(0.0), wildm.get(a).unwrap_or(0.0)))
        .collect();
    // Collections.sort(results, Collections.reverseOrder()): stable
    out.sort_by(|x, y| y.compare(x));
    out
}

/// ComplexCoverage
#[derive(Clone, Debug)]
pub struct ComplexCoverage {
    pub total: i32,
    pub unique: i32,
    pub shared: i32,
    pub wild: i32,
    pub alleles: Vec<AlleleCoverage>,
    pub cohort_frequency_penalty: f64,
    pub cohort_frequency_total: f64,
    pub score: f64,
    pub complexity_penalty: f64,
    pub complexity: f64,
    pub homozygous_count: i32,
    pub recovered_count: i32,
    pub wildcard_count: i32,
}

impl ComplexCoverage {
    /// ComplexCoverage.create: the sums (shared and wild rounded), the alleles sorted by name
    pub fn create(gc: &GeneCache, alleles: &Alleles, list: Vec<AlleleCoverage>) -> ComplexCoverage {
        let mut unique = 0;
        let mut shared = 0.0;
        let mut wild = 0.0;
        for c in &list {
            unique += c.unique;
            shared += c.shared;
            wild += c.wild;
        }
        let mut sorted = list;
        sorted.sort_by(|a, b| hmf_common::java_text::compare(alleles.name(a.allele), alleles.name(b.allele)));
        Self::with_counts(gc, unique, crate::math::round_int(shared), crate::math::round_int(wild), sorted)
    }

    pub fn with_counts(gc: &GeneCache, unique: i32, shared: i32, wild: i32, alleles: Vec<AlleleCoverage>) -> ComplexCoverage {
        ComplexCoverage {
            total: unique + shared + wild,
            unique,
            shared,
            wild,
            homozygous_count: gc.expect_allele_count as i32 - alleles.len() as i32,
            alleles,
            cohort_frequency_penalty: 0.0,
            cohort_frequency_total: 0.0,
            score: 0.0,
            complexity_penalty: 0.0,
            complexity: 0.0,
            recovered_count: 0,
            wildcard_count: 0,
        }
    }

    pub fn empty(gc: &GeneCache) -> ComplexCoverage {
        Self::with_counts(gc, 0, 0, 0, Vec::new())
    }

    pub fn allele_ids(&self) -> Vec<AId> {
        self.alleles.iter().map(|x| x.allele).collect()
    }

    pub fn is_homozygous(&self, alleles: &Alleles, a: AId) -> bool {
        let g = alleles.gene(a);
        self.alleles.iter().filter(|x| alleles.gene(x.allele) == g).count() == 1
    }

    /// expandToSixAlleles: a homozygous gene's coverage split over two copies of its allele. Java sets the
    /// second copy's wildcard coverage to zero (LI12); the correction gives it the remainder.
    pub fn expand_to_six(&mut self, gc: &GeneCache, alleles: &Alleles) {
        if self.alleles.len() == gc.expect_allele_count {
            return;
        }
        let existing = std::mem::take(&mut self.alleles);
        for g in &gc.gene_names {
            let gene_cov: Vec<AlleleCoverage> = existing.iter().filter(|x| alleles.gene(x.allele) == *g).cloned().collect();
            if gene_cov.len() == 2 || gene_cov.len() != 1 {
                self.alleles.extend(gene_cov);
                continue;
            }
            let s = &gene_cov[0];
            if s.wild > 0.0 {
                crate::debug!("LI12: homozygous {} has wildcard coverage {}", alleles.name(s.allele), s.wild);
            }
            let first = AlleleCoverage::new(s.allele, s.unique / 2, s.shared / 2.0, s.wild / 2.0);
            let rest_wild = if crate::compat::java("LI12") { 0.0 } else { s.wild - first.wild };
            let rest = AlleleCoverage::new(s.allele, s.unique - first.unique, s.shared - first.shared, rest_wild);
            self.alleles.push(first);
            self.alleles.push(rest);
        }
    }

    /// populateMissingCoverage: one entry per allele of the solution, zero where there was no coverage
    pub fn populate_missing(&mut self, gc: &GeneCache, solution: &[AId]) {
        if self.alleles.len() == gc.expect_allele_count {
            return;
        }
        let existing = std::mem::take(&mut self.alleles);
        let mut seen: Vec<AId> = Vec::new();
        for a in solution {
            if seen.contains(a) {
                continue;
            }
            seen.push(*a);
            match existing.iter().find(|x| x.allele == *a) {
                Some(c) => self.alleles.push(c.clone()),
                None => self.alleles.push(AlleleCoverage::new(*a, 0, 0.0, 0.0)),
            }
        }
    }

    /// HlaComplexFile.asString
    pub fn as_line(&self, alleles: &Alleles, genes: Option<&str>) -> String {
        use hmf_common::fmt::java_fixed;
        let mut v: Vec<String> = Vec::new();
        if let Some(g) = genes {
            v.push(g.to_string());
        }
        v.push(java_fixed(self.score, 2));
        v.push(java_fixed(self.complexity_penalty, 2));
        v.push(java_fixed(self.complexity, 2));
        v.push(self.homozygous_count.to_string());
        v.push(java_fixed(self.cohort_frequency_penalty, 2));
        v.push(java_fixed(self.cohort_frequency_total, 2));
        v.push(self.recovered_count.to_string());
        v.push(self.wildcard_count.to_string());
        v.push(self.total.to_string());
        v.push(self.unique.to_string());
        v.push(self.shared.to_string());
        v.push(self.wild.to_string());
        if genes.is_some() {
            v.push(self.alleles.iter().map(|x| x.display(alleles)).collect::<Vec<_>>().join(";"));
        } else {
            for x in &self.alleles {
                v.push(x.display(alleles));
            }
        }
        v.join("\t")
    }
}

/// ComplexBuilder.calcProteinCoverage
pub fn calc_protein_coverage(gc: &GeneCache, alleles: &Alleles, fas: &[FragAlleles], keep: &[AId]) -> ComplexCoverage {
    let f = filter_frag_alleles(alleles, fas, keep);
    ComplexCoverage::create(gc, alleles, build_coverage(alleles, &f, false))
}

fn calc_group_coverage(gc: &GeneCache, alleles: &Alleles, fas: &[FragAlleles], keep: &[AId]) -> ComplexCoverage {
    let f = filter_frag_alleles(alleles, fas, keep);
    ComplexCoverage::create(gc, alleles, build_coverage(alleles, &f, true))
}

/// FragmentAlleleMapper
pub struct Mapper {
    pub aa_het_loci: Vec<(HlaGene, BTreeMap<i32, Vec<Tok>>)>,
    pub nuc_het_loci: Vec<(HlaGene, Vec<i32>)>,
    pub ref_nucleotides: Vec<Vec<Tok>>,
    /// known stop-loss allele and the uids of its fragments
    pub stop_loss: Vec<(AId, Vec<u64>)>,
    pub threads: usize,
}

fn het_loci_contains(v: &[i32], x: i32) -> bool {
    v.contains(&x)
}

impl Mapper {
    /// createFragmentAlleles: the fragments with full support; the scope of the others is set
    pub fn create_fragment_alleles(&self, rd: &RefData, params: &Params, frags: &mut [Fragment], subset: &[usize], aa_seqs: &[&SeqLoci], nuc_seqs: &[&SeqLoci]) -> Vec<FragAlleles> {
        crate::info!(
            "building frag-alleles from aminoAcids(frags={} candSeq={}) nucFrags(hetLoci={} candSeq={} nucs={}) knownIndels({})",
            subset.len(),
            aa_seqs.len(),
            self.nuc_het_loci.len(),
            nuc_seqs.len(),
            self.ref_nucleotides.len(),
            self.stop_loss.iter().map(|x| x.1.len()).sum::<usize>()
        );
        // each fragment is mapped on its own (in parallel), then the scopes are set in order
        let todo: Vec<usize> = subset.iter().copied().filter(|fi| frags[*fi].scope != Scope::HlaY).collect();
        let frags_ro: &[Fragment] = frags;
        let map_one = |fi: usize| -> FragAlleles {
            match self.stop_loss.iter().find(|x| x.1.contains(&frags_ro[fi].uid)) {
                Some(x) => FragAlleles::new(&rd.alleles, fi, &[x.0], &[]),
                None => self.map_fragment(rd, params, &frags_ro[fi], fi, aa_seqs, nuc_seqs),
            }
        };
        let mapped: Vec<FragAlleles> = if self.threads > 1 && todo.len() >= 64 {
            let chunk = todo.len().div_ceil(self.threads);
            std::thread::scope(|s| {
                let handles: Vec<_> = todo.chunks(chunk).map(|c| s.spawn(move || c.iter().map(|fi| map_one(*fi)).collect::<Vec<_>>())).collect();
                handles.into_iter().flat_map(|h| h.join().unwrap()).collect()
            })
        } else {
            todo.iter().map(|fi| map_one(*fi)).collect()
        };
        let mut out = Vec::new();
        for fa in mapped {
            if !fa.full.is_empty() {
                out.push(fa);
            } else {
                let fi = fa.fi;
                self.set_fragment_scope(&mut frags[fi], &fa);
            }
        }
        out
    }

    fn set_fragment_scope(&self, f: &mut Fragment, fa: &FragAlleles) {
        if f.scope == Scope::HlaY {
            return;
        }
        f.clear_scope();
        if !fa.full.is_empty() {
            return;
        }
        if !fa.wild.is_empty() {
            f.set_scope(Scope::WildOnly, false);
            return;
        }
        let mut has_het = false;
        for (g, m) in &self.aa_het_loci {
            if !f.contains_gene(*g) {
                continue;
            }
            if f.aas.keys().any(|l| m.contains_key(&l)) {
                has_het = true;
                break;
            }
        }
        f.set_scope(if has_het { Scope::UnmatchedAminoAcid } else { Scope::NoHetLoci }, false);
    }

    fn map_fragment(&self, rd: &RefData, params: &Params, f: &Fragment, fi: usize, aa_seqs: &[&SeqLoci], nuc_seqs: &[&SeqLoci]) -> FragAlleles {
        let alleles = &rd.alleles;
        let nm = self.nucleotide_matches(rd, params, f, nuc_seqs);
        let nm_order = nm.java_keys(alleles);
        let full_nuc: Vec<AId> = nm_order.iter().copied().filter(|a| nm.get(*a) == Some(MatchType::Full)).collect();
        let wild_nuc: Vec<AId> = nm_order.iter().copied().filter(|a| nm.get(*a) == Some(MatchType::Wild)).collect();
        let am = self.amino_acid_matches(rd, params, f, aa_seqs);
        let am_order = am.java_keys(alleles);
        let full_aa: Vec<AId> = am_order.iter().copied().filter(|a| am.get(*a) == Some(MatchType::Full)).collect();
        let mut wild_aa: Vec<AId> = am_order.iter().copied().filter(|a| am.get(*a) == Some(MatchType::Wild)).collect();
        let hom_aa: Vec<AId> = am_order.iter().copied().filter(|a| am.get(*a) == Some(MatchType::NoLoci)).collect();
        if full_nuc.is_empty() && wild_nuc.is_empty() {
            if full_aa.is_empty() && !wild_aa.is_empty() && !hom_aa.is_empty() {
                return FragAlleles::new(alleles, fi, &[], &[]);
            }
            return FragAlleles::new(alleles, fi, &full_aa, &wild_aa);
        }
        let full_set: HashSet<AId> = full_aa.iter().copied().collect();
        let hom_set: HashSet<AId> = hom_aa.iter().copied().collect();
        let consistent: Vec<AId> = full_nuc.iter().copied().filter(|a| full_set.contains(a) || hom_set.contains(a)).collect();
        for a in &full_aa {
            if !wild_aa.contains(a) && wild_nuc.contains(a) {
                wild_aa.push(*a);
            }
        }
        FragAlleles::new(alleles, fi, &consistent, &wild_aa)
    }

    fn nucleotide_matches(&self, rd: &RefData, params: &Params, f: &Fragment, seqs: &[&SeqLoci]) -> AlleleMap<MatchType> {
        let alleles = &rd.alleles;
        let mut matches: AlleleMap<MatchType> = AlleleMap::new();
        let mut gene_loci: Vec<(HlaGene, Vec<i32>, Vec<Tok>)> = Vec::new();
        let mut missed: HashMap<i32, Tok> = HashMap::new();
        for (g, ref_loci) in &self.nuc_het_loci {
            let mut matched: Vec<i32> = f.nucs.keys().filter(|l| het_loci_contains(ref_loci, *l)).collect();
            for &m in ref_loci.iter() {
                if f.nucs.contains(m) {
                    continue;
                }
                if m < 0 || m as usize >= self.ref_nucleotides.len() {
                    continue;
                }
                let lq = f.raw_nucleotide(m, true, params);
                if lq == crate::tok::EMPTY {
                    continue;
                }
                if self.ref_nucleotides[m as usize].contains(&lq) {
                    matched.push(m);
                    missed.insert(m, lq);
                }
            }
            matched.sort();
            let nucs: Vec<Tok> = matched.iter().map(|l| missed.get(l).copied().unwrap_or_else(|| f.nucleotide(*l))).collect();
            gene_loci.push((*g, matched, nucs));
        }
        if gene_loci.iter().all(|x| x.1.is_empty()) {
            return matches;
        }
        for s in seqs {
            let a = s.allele;
            let protein = alleles.four_digit(a);
            let existing = matches.get(protein);
            if existing == Some(MatchType::Full) {
                continue;
            }
            let gene = alleles.gene(a);
            if !f.contains_gene(gene) {
                continue;
            }
            let Some(gl) = gene_loci.iter().find(|x| x.0 == gene) else { continue };
            if gl.1.is_empty() {
                continue;
            }
            let mt = if s.has_exon_boundary_wildcards {
                let bounds = rd.gene_cache.nuc_boundaries(gene);
                let mut loci = Vec::with_capacity(gl.1.len());
                let mut nucs = Vec::with_capacity(gl.1.len());
                for (k, &l) in gl.1.iter().enumerate() {
                    let wb = bounds.contains(&l) && (l as usize) < s.len() && s.at(l as usize) == WILD;
                    if !wb {
                        loci.push(l);
                        nucs.push(gl.2[k]);
                    }
                }
                if loci.is_empty() { MatchType::Full } else { s.match_type_toks(&nucs, &loci) }
            } else {
                s.match_type_toks(&gl.2, &gl.1)
            };
            if mt == MatchType::Mismatch {
                continue;
            }
            if existing.is_none() || mt.is_better(existing.unwrap()) {
                matches.update(protein, mt, |_| mt);
            }
        }
        matches
    }

    fn amino_acid_matches(&self, rd: &RefData, params: &Params, f: &Fragment, seqs: &[&SeqLoci]) -> AlleleMap<MatchType> {
        let alleles = &rd.alleles;
        let mut matches: AlleleMap<MatchType> = AlleleMap::new();
        let mut gene_loci: Vec<(HlaGene, Vec<i32>, Vec<Tok>)> = Vec::new();
        for (g, het) in &self.aa_het_loci {
            let mut loci: Vec<i32> = f.aas.keys().filter(|l| het.contains_key(l)).collect();
            let mut missed: HashMap<i32, Tok> = HashMap::new();
            for (&m, cands) in het.iter() {
                if loci.binary_search(&m).is_ok() {
                    continue;
                }
                let lq = f.low_qual_aa(m, true, params);
                if lq == crate::tok::EMPTY {
                    continue;
                }
                if cands.contains(&lq) {
                    missed.insert(m, lq);
                }
            }
            if !missed.is_empty() {
                loci.extend(missed.keys().copied());
                loci.sort();
            }
            let aas: Vec<Tok> = loci.iter().map(|l| missed.get(l).copied().unwrap_or_else(|| f.aa(*l))).collect();
            gene_loci.push((*g, loci, aas));
        }
        for s in seqs {
            let a = s.allele;
            let gene = alleles.gene(a);
            if !f.contains_gene(gene) {
                matches.update(a, MatchType::NoLoci, |_| MatchType::NoLoci);
                continue;
            }
            let gl = gene_loci.iter().find(|x| x.0 == gene);
            let Some(gl) = gl.filter(|x| !x.1.is_empty()) else {
                matches.update(a, MatchType::NoLoci, |_| MatchType::NoLoci);
                continue;
            };
            let mt = if s.has_exon_boundary_wildcards {
                let bounds = rd.gene_cache.aa_boundaries(gene);
                let mut loci = Vec::with_capacity(gl.1.len());
                let mut aas = Vec::with_capacity(gl.1.len());
                for (k, &l) in gl.1.iter().enumerate() {
                    let wb = bounds.contains(&l) && (l as usize) < s.len() && s.at(l as usize) == WILD;
                    if !wb {
                        loci.push(l);
                        aas.push(gl.2[k]);
                    }
                }
                s.match_type_toks(&aas, &loci)
            } else {
                s.match_type_toks(&gl.2, &gl.1)
            };
            if mt == MatchType::Mismatch {
                continue;
            }
            matches.update(a, mt, |_| mt);
        }
        matches
    }
}

/// findWildcardAlleles: the wildcard alleles of the fragments, in a HashSet's order
pub fn find_wildcard_alleles(alleles: &Alleles, fas: &[FragAlleles]) -> Vec<AId> {
    let mut keys: Vec<AId> = Vec::new();
    let mut seen: HashSet<AId> = HashSet::new();
    for fa in fas {
        for a in fa.full.iter().chain(fa.wild.iter()) {
            if alleles.has_wildcards(*a) && seen.insert(*a) {
                keys.push(*a);
            }
        }
    }
    allele_set_order(alleles, &keys, grown_capacity(16, keys.len()))
}

/// findUnsupportedWildcards: wildcard alleles with fewer than two fragments that support no other allele
pub fn find_unsupported_wildcards(alleles: &Alleles, fas: &[FragAlleles], wildcards: &[AId]) -> Vec<AId> {
    let mut supported = Vec::new();
    for &a in wildcards {
        let mut n = 0;
        for fa in fas {
            if !fa.has_full(a) {
                continue;
            }
            if fa.full.iter().any(|x| !alleles.has_wildcards(*x)) {
                continue;
            }
            n += 1;
        }
        if n >= MIN_WILDCARD_FRAGMENTS {
            supported.push(a);
        }
    }
    if !supported.is_empty() {
        crate::info!("filtered wildcards({} -> {}) supported alleles: {}", wildcards.len(), supported.len(), alleles.join(&supported));
    } else if !wildcards.is_empty() {
        crate::info!("removed all {} wildcard candidates", wildcards.len());
    }
    wildcards.iter().copied().filter(|a| !supported.contains(a)).collect()
}

/// filterUnsupportedWildcardFragments: the alleles removed from every fragment; fragments left with none are
/// dropped and get WILD_ONLY
pub fn filter_unsupported_wildcard_fragments(fas: &mut Vec<FragAlleles>, frags: &mut [Fragment], unsupported: &[AId]) {
    let mut removed = 0;
    let mut i = 0;
    while i < fas.len() {
        fas[i].remove_all(unsupported);
        if fas[i].full.is_empty() && fas[i].wild.is_empty() {
            frags[fas[i].fi].set_scope(Scope::WildOnly, false);
            fas.remove(i);
            removed += 1;
        } else {
            i += 1;
        }
    }
    if removed > 0 {
        crate::debug!("removed frags({removed}) with wildcard alleles removed");
    }
}

/// HlaComplex: alleles in gene order
pub type Complex = Vec<AId>;

/// ComplexBuilder
pub struct ComplexBuilder {
    pub unique_group_alleles: Vec<AId>,
    pub unique_protein_alleles: Vec<AId>,
    pub confirmed_recovered_alleles: Vec<AId>,
    pub gene_count: usize,
}

fn take_n(v: Vec<AId>, n: usize) -> Vec<AId> {
    v.into_iter().take(n).collect()
}

impl ComplexBuilder {
    pub fn new(gene_count: usize) -> ComplexBuilder {
        ComplexBuilder { unique_group_alleles: Vec::new(), unique_protein_alleles: Vec::new(), confirmed_recovered_alleles: Vec::new(), gene_count }
    }

    /// findUnique: per gene the two best uniquely supported alleles (one per group where groups are given)
    fn find_unique(rd: &RefData, cov: &ComplexCoverage, confirmed_groups: &[AId], total_frags: usize) -> Vec<AlleleCoverage> {
        let alleles = &rd.alleles;
        let is_group = confirmed_groups.is_empty();
        let required = if is_group { total_frags as f64 * MIN_HI_CONF_UNIQUE_GROUP_COVERAGE } else { total_frags as f64 * MIN_HI_CONF_UNIQUE_PROTEIN_COVERAGE };
        let mut unique: Vec<AlleleCoverage> = cov.alleles.iter().filter(|x| x.unique as f64 >= required).cloned().collect();
        unique.sort_by(|x, y| y.compare(x));
        let mut results: Vec<AlleleCoverage> = Vec::new();
        for g in &rd.gene_cache.gene_names {
            let gene_cov: Vec<&AlleleCoverage> = unique.iter().filter(|x| alleles.gene(x.allele) == *g).collect();
            if gene_cov.is_empty() {
                continue;
            }
            let gene_groups: Vec<AId> = confirmed_groups.iter().copied().filter(|x| alleles.gene(*x) == *g).collect();
            let mut count = 0;
            for c in gene_cov {
                if gene_groups.len() > 1 {
                    let grp = alleles.group(c.allele);
                    let matched = results.iter().filter(|x| gene_groups.contains(&alleles.group(x.allele)) && alleles.group(x.allele) == grp).count();
                    if matched >= 1 {
                        continue;
                    }
                }
                results.push(c.clone());
                count += 1;
                if count >= 2 {
                    break;
                }
            }
        }
        results.sort_by(|x, y| y.compare(x));
        results
    }

    /// filterWithUniqueGroups
    fn filter_with_unique_groups(rd: &RefData, cands: &[AId], groups: &[AId], recovered: &[AId]) -> Vec<AId> {
        let alleles = &rd.alleles;
        let mut out = Vec::new();
        for &a in cands {
            let g = alleles.gene(a);
            let gene_groups: Vec<AId> = groups.iter().copied().filter(|x| alleles.gene(*x) == g).collect();
            let grp = alleles.group(a);
            if recovered.contains(&a) && !gene_groups.contains(&grp) {
                continue;
            }
            if gene_groups.len() < 2 || gene_groups.contains(&grp) {
                out.push(a);
            }
        }
        out
    }

    /// filterCandidates
    pub fn filter_candidates(&mut self, rd: &RefData, fas: &[FragAlleles], candidates: &[AId], recovered: &[AId]) {
        let alleles = &rd.alleles;
        crate::info!("filtering candidates from fragAlleles({}) candidates({}) recovered({})", fas.len(), candidates.len(), recovered.len());
        let group_cov = calc_group_coverage(&rd.gene_cache, alleles, fas, candidates);
        let total = fas.len();
        let unique_groups = Self::find_unique(rd, &group_cov, &[], total);
        let low_unique = total as f64 * MIN_LOW_CONF_GROUP_UNIQUE_COVERAGE;
        let low_total = total as f64 * MIN_LOW_CONF_GROUP_TOTAL_COVERAGE;
        let mut low_conf: Vec<AlleleCoverage> = Vec::new();
        let mut discarded: Vec<AlleleCoverage> = Vec::new();
        for c in &group_cov.alleles {
            if unique_groups.iter().any(|u| u.allele == c.allele) {
                continue;
            }
            if c.unique as f64 >= low_unique || c.total >= low_total {
                low_conf.push(c.clone());
            } else {
                discarded.push(c.clone());
            }
        }
        crate::debug!("  confirmed {} unique allele groups{}{}", unique_groups.len(), if unique_groups.is_empty() { "" } else { ": " }, coverage_list_str(alleles, &unique_groups));
        crate::debug!("  found {} insufficiently unique allele groups{}{}", low_conf.len(), if low_conf.is_empty() { "" } else { ": " }, coverage_list_str(alleles, &low_conf));
        crate::debug!("  discarded {} allele groups{}{}", discarded.len(), if discarded.is_empty() { "" } else { ": " }, coverage_list_str(alleles, &discarded));
        let unique_group_alleles: Vec<AId> = unique_groups.iter().map(|x| x.allele).collect();
        let mut after = Self::filter_with_unique_groups(rd, candidates, &unique_group_alleles, recovered);

        // getTopLowConfGroups
        let mut per_gene: HashMap<HlaGene, i32> = HashMap::new();
        for u in &unique_groups {
            *per_gene.entry(alleles.gene(u.allele)).or_insert(0) += 1;
        }
        let mut top_low: Vec<AId> = Vec::new();
        for c in &low_conf {
            let g = alleles.gene(c.allele);
            let n = per_gene.entry(g).or_insert(0);
            if *n < 2 {
                top_low.push(c.allele);
            }
            *n += 1;
        }
        let common_low: Vec<AId> = rd.common_alleles.iter().copied().filter(|x| !after.contains(x) && top_low.contains(&alleles.group(*x))).collect();
        after.extend(common_low.iter().copied());
        crate::debug!("  keeping {} common allele(s) from insufficiently unique allele groups{}{}", common_low.len(), if common_low.is_empty() { "" } else { ": " }, alleles.join(&common_low));

        let rare: Vec<AId> = after.iter().copied().filter(|x| rd.freqs.frequency(alleles.gene(*x), alleles.name(*x)) <= RARE_ALLELES_FREQ_CUTOFF).collect();
        let dubious_groups: HashSet<AId> = after.iter().copied().filter(|x| alleles.has_wildcards(*x) || rare.contains(x)).map(|x| alleles.group(x)).collect();
        let common_dubious: Vec<AId> = rd.common_alleles.iter().copied().filter(|x| !after.contains(x) && dubious_groups.contains(&alleles.group(*x))).collect();
        after.extend(common_dubious.iter().copied());
        crate::debug!(
            "  keeping {} common allele(s) in the same 2-digit group as wildcard or rare allele candidates{}{}",
            common_dubious.len(),
            if common_dubious.is_empty() { "" } else { ": " },
            alleles.join(&common_dubious)
        );
        let stop_loss: Vec<AId> = recovered.iter().copied().filter(|x| rd.stop_loss_alleles.iter().any(|s| s.1 == *x) && !after.contains(x)).collect();
        after.extend(stop_loss.iter().copied());
        crate::debug!("  keeping {} allele(s) with stop loss indel{}{}", stop_loss.len(), if stop_loss.is_empty() { "" } else { ": " }, alleles.join(&stop_loss));
        for r in recovered {
            if after.contains(r) {
                self.confirmed_recovered_alleles.push(*r);
            }
        }
        if !self.confirmed_recovered_alleles.is_empty() {
            hmf_common::java_sort::java_sort(&mut self.confirmed_recovered_alleles, |a, b| alleles.compare(*a, *b));
            crate::debug!("  keeping {} recovered alleles from insufficiently unique groups: {}", self.confirmed_recovered_alleles.len(), alleles.join(&self.confirmed_recovered_alleles));
        } else if !recovered.is_empty() {
            crate::info!("  no recovered alleles kept from insufficiently unique groups");
        }
        let protein_cov = calc_protein_coverage(&rd.gene_cache, alleles, fas, &after);
        let unique_proteins: Vec<AlleleCoverage> = Self::find_unique(rd, &protein_cov, &unique_group_alleles, total).into_iter().filter(|x| !recovered.contains(&x.allele)).collect();
        let mut insufficient: Vec<AlleleCoverage> = protein_cov.alleles.iter().filter(|x| x.unique > 0 && !unique_proteins.iter().any(|u| u.allele == x.allele)).cloned().collect();
        insufficient.sort_by(|x, y| y.compare(x));
        crate::info!("  confirmed {} unique proteins{}{}", unique_proteins.len(), if unique_proteins.is_empty() { "" } else { ": " }, coverage_list_str(alleles, &unique_proteins));
        crate::info!("  found {} insufficiently unique proteins{}{}", insufficient.len(), if insufficient.is_empty() { "" } else { ": " }, coverage_list_str(alleles, &insufficient));
        // unique protein filtering is no longer applied: filterWithUniqueProteins with no confirmed proteins
        self.unique_protein_alleles = Self::filter_with_unique_groups(rd, &after, &[], &[]);
        self.unique_group_alleles = unique_group_alleles;
    }

    /// buildComplexesByGene
    pub fn build_complexes_by_gene(alleles: &Alleles, gene: HlaGene, groups: &[AId], candidates: &[AId]) -> Vec<Complex> {
        let confirmed = take_n(groups.iter().copied().filter(|x| alleles.gene(*x) == gene).collect(), 2);
        let cands: Vec<AId> = candidates.iter().copied().filter(|x| alleles.gene(*x) == gene).collect();
        // Sets.newHashSet(candidates): a HashSet from a collection
        let cand_set = {
            let d = dedup(&cands);
            allele_set_order(alleles, &d, jdk19_collection_capacity(cands.len()))
        };
        let group_set = |g: AId| -> Vec<AId> {
            let d = dedup(&cands.iter().copied().filter(|x| alleles.group(*x) == g).collect::<Vec<_>>());
            allele_set_order(alleles, &d, grown_capacity(16, d.len()))
        };
        if confirmed.len() == 2 {
            return combine_alleles(alleles, &group_set(confirmed[0]), &group_set(confirmed[1]));
        }
        if confirmed.len() == 1 {
            let first = group_set(confirmed[0]);
            let mut out: Vec<Complex> = first.iter().map(|a| vec![*a]).collect();
            out.extend(combine_alleles(alleles, &first, &cand_set));
            return out;
        }
        let mut out: Vec<Complex> = cands.iter().map(|a| vec![*a]).collect();
        out.extend(combine_alleles(alleles, &cand_set, &cand_set));
        out
    }

    /// buildAlleleComplexes: the per-gene complexes combined over the genes, in nested order
    pub fn build_allele_complexes(rd: &RefData, groups: &[AId], candidates: &[AId]) -> Vec<Vec<Complex>> {
        rd.gene_cache.gene_names.iter().filter(|g| !g.is_pseudo()).map(|g| Self::build_complexes_by_gene(&rd.alleles, *g, groups, candidates)).collect()
    }
}

/// cartesianAlleleProduct + combineAlleles: the pairs in a HashSet<Pair>'s order (Pair hash: left ^ right)
fn combine_alleles(alleles: &Alleles, first: &[AId], second: &[AId]) -> Vec<Complex> {
    let mut pairs: Vec<(AId, AId)> = Vec::new();
    let mut seen: HashSet<(AId, AId)> = HashSet::new();
    for &a1 in first {
        for &a2 in second {
            if a1 == a2 {
                continue;
            }
            let p = if alleles.compare(a1, a2) != Ordering::Greater { (a1, a2) } else { (a2, a1) };
            if seen.insert(p) {
                pairs.push(p);
            }
        }
    }
    let cap = grown_capacity(16, pairs.len());
    let mut keyed: Vec<(usize, (AId, AId))> = pairs.into_iter().map(|p| (bucket_index(alleles.hash(p.0) ^ alleles.hash(p.1), cap), p)).collect();
    keyed.sort_by_key(|x| x.0);
    keyed.into_iter().map(|x| vec![x.1 .0, x.1 .1]).collect()
}
