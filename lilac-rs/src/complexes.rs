//! FragmentAlleleMatrix, ComplexCoverageCalculator with its CoverageCalcTasks, and ComplexCoverageRanking.
//!
//! Java computes, for every candidate complex (up to about a million), the coverage of its alleles over
//! every fragment, sums the shared and wildcard fractions in floating point and rounds them. Only the
//! rounded totals are used to rank the complexes, so the port computes them exactly with integers (each
//! fraction 1/k, k <= 12 for six alleles, is a whole number of 1/27720 units). The sum in floating point
//! rounds the same way except when the exact value is half way between two integers; those complexes are
//! recomputed as Java does. Fragments supported by the alleles of one or two genes only are summed once per
//! gene complex or pair of gene complexes, so a complex costs a handful of additions plus the fragments
//! supported by three or more genes. The complexes Java keeps get their per-allele coverage computed as
//! Java does (in floating point, fragment by fragment).

use std::cmp::Ordering;
use std::collections::HashMap;

use crate::coverage::{AlleleCoverage, Complex, ComplexCoverage, FragAlleles};
use crate::hla::{AId, Alleles, HlaGene};
use crate::refdata::RefData;
use crate::seq::{ExonSequence, SeqLoci};
use crate::tok::WILD;

const NONE: u8 = 0;
const FULL: u8 = 1;
const WILDB: u8 = 2;
/// 1/k for k = 1..12 in whole units
const UNITS: i64 = 27720;

/// FragmentAlleleMatrix: the support of each fragment for each allele
pub struct Matrix {
    pub alleles: Vec<AId>,
    index: HashMap<AId, usize>,
    pub frag_count: usize,
    /// [allele][fragment]: NONE, FULL, WILD or both
    dense: Vec<Vec<u8>>,
}

impl Matrix {
    pub fn new(fas: &[FragAlleles], alleles: &[AId]) -> Matrix {
        let mut index = HashMap::new();
        for (i, a) in alleles.iter().enumerate() {
            index.insert(*a, i);
        }
        let mut dense = vec![vec![NONE; fas.len()]; alleles.len()];
        for (f, fa) in fas.iter().enumerate() {
            for a in &fa.full {
                if let Some(i) = index.get(a) {
                    dense[*i][f] = FULL;
                }
            }
            for a in &fa.wild {
                if let Some(i) = index.get(a) {
                    dense[*i][f] |= WILDB;
                }
            }
        }
        Matrix { alleles: alleles.to_vec(), index, frag_count: fas.len(), dense }
    }

    pub fn index_of(&self, a: AId) -> Option<usize> {
        self.index.get(&a).copied()
    }

    /// FragmentAlleleMatrix.create: the alleles' coverage in Java's floating point
    pub fn create(&self, complex: &[AId]) -> Vec<AlleleCoverage> {
        let mut idx = Vec::with_capacity(complex.len());
        for a in complex {
            match self.index.get(a) {
                Some(i) => idx.push(*i),
                None => return Vec::new(),
            }
        }
        let n = idx.len();
        let mut unique = vec![0i32; n];
        let mut shared = vec![0f64; n];
        let mut wild = vec![0f64; n];
        let mut fl = vec![false; n];
        let mut wl = vec![false; n];
        for f in 0..self.frag_count {
            let mut fc = 0;
            let mut wc = 0;
            let mut fi = 0usize;
            for k in 0..n {
                let s = self.dense[idx[k]][f];
                fl[k] = s & FULL != 0;
                wl[k] = s & WILDB != 0;
                if fl[k] {
                    fc += 1;
                    fi = k;
                }
                if wl[k] {
                    wc += 1;
                }
            }
            if fc == 1 && wc == 0 {
                unique[fi] += 1;
            } else if fc > 0 || wc > 0 {
                let c = 1.0 / (fc + wc) as f64;
                for k in 0..n {
                    if fl[k] {
                        shared[k] += c;
                    }
                    if wl[k] {
                        wild[k] += c;
                    }
                }
            }
        }
        (0..n).map(|k| AlleleCoverage::new(complex[k], unique[k], shared[k], wild[k])).collect()
    }
}

/// CoverageCalcTask.calcTotalFragments
pub fn total_fragments(list: &[AlleleCoverage]) -> i32 {
    let mut u = 0;
    let mut s = 0.0;
    let mut w = 0.0;
    for c in list {
        u += c.unique;
        s += c.shared;
        w += c.wild;
    }
    u + crate::math::round_int(s) + crate::math::round_int(w)
}

/// The rounded totals of a complex: unique, shared and wildcard coverage
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct Totals {
    pub unique: i32,
    pub shared: i32,
    pub wild: i32,
}

impl Totals {
    pub fn total(&self) -> i32 {
        self.unique + self.shared + self.wild
    }
}

/// half-up rounding of n / UNITS; None when exactly half way (floating point decides)
#[inline]
fn round_units(n: i64) -> Option<i32> {
    if n % UNITS == UNITS / 2 {
        return None;
    }
    Some(((n + UNITS / 2) / UNITS) as i32)
}

#[derive(Copy, Clone, Default)]
struct Acc {
    unique: i64,
    shared: i64,
    wild: i64,
}

impl Acc {
    fn add(&mut self, o: &Acc) {
        self.unique += o.unique;
        self.shared += o.shared;
        self.wild += o.wild;
    }
    #[inline]
    fn fragment(&mut self, fc: i64, wc: i64) {
        if fc == 1 && wc == 0 {
            self.unique += 1;
        } else if fc + wc > 0 {
            let c = UNITS / (fc + wc);
            self.shared += fc * c;
            self.wild += wc * c;
        }
    }
}

/// The candidate complexes as Java builds them: the product of per-gene complexes, in nested order
pub struct Product {
    pub per_gene: Vec<Vec<Complex>>,
}

impl Product {
    pub fn len(&self) -> usize {
        if self.per_gene.is_empty() {
            return 0;
        }
        self.per_gene.iter().map(|x| x.len()).product()
    }
    /// the per-gene complex indices of complex `i`
    pub fn decode(&self, mut i: usize, out: &mut [usize]) {
        for g in (0..self.per_gene.len()).rev() {
            let n = self.per_gene[g].len();
            out[g] = i % n;
            i /= n;
        }
    }
    pub fn alleles(&self, i: usize) -> Vec<AId> {
        let mut idx = vec![0; self.per_gene.len()];
        self.decode(i, &mut idx);
        let mut v = Vec::with_capacity(6);
        for (g, x) in idx.iter().enumerate() {
            v.extend_from_slice(&self.per_gene[g][*x]);
        }
        v
    }
}

/// Exact totals for every complex of a product
struct FastTotals<'a> {
    product: &'a Product,
    matrix: &'a Matrix,
    /// per gene, per gene complex: matrix allele indices
    gene_idx: Vec<Vec<Vec<usize>>>,
    /// per single gene: the totals of its pure fragments for each gene complex
    single: Vec<Vec<Acc>>,
    /// per pair of genes (g1 < g2): totals of the fragments supported by exactly those two genes
    pairs: Vec<((usize, usize), Vec<Acc>)>,
    /// the fragments supported by three or more genes
    multi: Vec<usize>,
}

impl<'a> FastTotals<'a> {
    fn new(product: &'a Product, matrix: &'a Matrix, alleles: &Alleles, genes: &[HlaGene]) -> FastTotals<'a> {
        let ng = product.per_gene.len();
        let gene_idx: Vec<Vec<Vec<usize>>> = product
            .per_gene
            .iter()
            .map(|cs| cs.iter().map(|c| c.iter().map(|a| matrix.index_of(*a).unwrap()).collect()).collect())
            .collect();
        // the gene of each matrix allele (position in `genes`)
        let gene_pos: Vec<usize> = matrix.alleles.iter().map(|a| genes.iter().position(|g| *g == alleles.gene(*a)).unwrap_or(usize::MAX)).collect();
        let mut masks = vec![0u32; matrix.frag_count];
        for (ai, col) in matrix.dense.iter().enumerate() {
            let gp = gene_pos[ai];
            if gp == usize::MAX {
                continue;
            }
            for (f, s) in col.iter().enumerate() {
                if *s != NONE {
                    masks[f] |= 1 << gp;
                }
            }
        }
        let counts = |f: usize, idx: &[usize]| -> (i64, i64) {
            let mut fc = 0;
            let mut wc = 0;
            for &i in idx {
                let s = matrix.dense[i][f];
                if s & FULL != 0 {
                    fc += 1;
                }
                if s & WILDB != 0 {
                    wc += 1;
                }
            }
            (fc, wc)
        };
        let mut single: Vec<Vec<Acc>> = vec![Vec::new(); ng];
        for g in 0..ng {
            let frags: Vec<usize> = (0..matrix.frag_count).filter(|f| masks[*f] == 1 << g).collect();
            single[g] = gene_idx[g]
                .iter()
                .map(|idx| {
                    let mut a = Acc::default();
                    for &f in &frags {
                        let (fc, wc) = counts(f, idx);
                        a.fragment(fc, wc);
                    }
                    a
                })
                .collect();
        }
        let mut pairs = Vec::new();
        for g1 in 0..ng {
            for g2 in g1 + 1..ng {
                let m = (1u32 << g1) | (1u32 << g2);
                let frags: Vec<usize> = (0..matrix.frag_count).filter(|f| masks[*f] == m).collect();
                let n2 = gene_idx[g2].len();
                let mut t = vec![Acc::default(); gene_idx[g1].len() * n2];
                if !frags.is_empty() {
                    for (i1, idx1) in gene_idx[g1].iter().enumerate() {
                        for (i2, idx2) in gene_idx[g2].iter().enumerate() {
                            let mut a = Acc::default();
                            for &f in &frags {
                                let (fc1, wc1) = counts(f, idx1);
                                let (fc2, wc2) = counts(f, idx2);
                                a.fragment(fc1 + fc2, wc1 + wc2);
                            }
                            t[i1 * n2 + i2] = a;
                        }
                    }
                }
                pairs.push(((g1, g2), t));
            }
        }
        let multi: Vec<usize> = (0..matrix.frag_count).filter(|f| masks[*f].count_ones() >= 3).collect();
        crate::debug!(
            "complex coverage: fragments({}) single-gene({}) multi-gene({})",
            matrix.frag_count,
            (0..matrix.frag_count).filter(|f| masks[*f].count_ones() == 1).count(),
            multi.len()
        );
        FastTotals { product, matrix, gene_idx, single, pairs, multi }
    }

    /// the totals of complex `i`, None where a floating point sum decides the rounding
    fn totals(&self, i: usize, idx: &mut [usize]) -> Option<Totals> {
        self.product.decode(i, idx);
        let mut acc = Acc::default();
        for (g, x) in idx.iter().enumerate() {
            acc.add(&self.single[g][*x]);
        }
        for ((g1, g2), t) in &self.pairs {
            if !t.is_empty() {
                acc.add(&t[idx[*g1] * self.gene_idx[*g2].len() + idx[*g2]]);
            }
        }
        if !self.multi.is_empty() {
            for &f in &self.multi {
                let mut fc = 0;
                let mut wc = 0;
                for (g, x) in idx.iter().enumerate() {
                    for &a in &self.gene_idx[g][*x] {
                        let s = self.matrix.dense[a][f];
                        if s & FULL != 0 {
                            fc += 1;
                        }
                        if s & WILDB != 0 {
                            wc += 1;
                        }
                    }
                }
                acc.fragment(fc, wc);
            }
        }
        Some(Totals { unique: acc.unique as i32, shared: round_units(acc.shared)?, wild: round_units(acc.wild)? })
    }
}

/// A scored complex: an index into its complex list and its totals
#[derive(Clone, Debug)]
pub struct Scored {
    pub index: u32,
    pub totals: Totals,
    pub score: f64,
    pub cohort_total: f64,
    pub cohort_penalty: f64,
    pub complexity: f64,
    pub complexity_penalty: f64,
    pub recovered: i32,
    pub wildcards: i32,
}

/// CoverageCalcTask.canCull with its running maximum
struct Culler {
    max: i32,
    perc: f64,
}

impl Culler {
    fn cull(&mut self, total: i32) -> bool {
        if total > self.max {
            self.max = total;
            return false;
        }
        if self.max - total < 40 {
            return false;
        }
        if total as f64 > self.max as f64 * (1.0 - self.perc) {
            return false;
        }
        true
    }
}

const CULL_COMPLEX_COUNT: usize = 100000;

/// ComplexCoverageCalculator.calculateComplexCoverages on a product: the totals of the complexes in the
/// order Java collects them, without those its tasks cull. With more than one thread and 10,000 complexes
/// or more Java deals the complexes out to the threads in turn, culls within each thread's share and joins
/// the shares; the correction (LI6) keeps the single-thread order and culling at any thread count.
/// LI7: Java builds a thread pool of -threads threads for 10,000 complexes or more (unless -threads is 1),
/// which throws for 0 or fewer
pub fn check_thread_pool(n: usize, java_threads: i32) -> anyhow::Result<()> {
    if crate::compat::java("LI7") && java_threads < 1 && n >= 10000 {
        anyhow::bail!("java.lang.IllegalArgumentException (Executors.newFixedThreadPool({java_threads}) in ComplexCoverageCalculator)");
    }
    Ok(())
}

pub fn calculate_product(product: &Product, fas: &[FragAlleles], alleles: &Alleles, genes: &[HlaGene], threads: usize, java_threads: i32, top_score_threshold: f64) -> anyhow::Result<Vec<(u32, Totals)>> {
    check_thread_pool(product.len(), java_threads)?;
    let mut all: Vec<AId> = Vec::new();
    for cs in &product.per_gene {
        for c in cs {
            for a in c {
                if !all.contains(a) {
                    all.push(*a);
                }
            }
        }
    }
    let matrix = Matrix::new(fas, &all);
    let n = product.len();
    let fast = FastTotals::new(product, &matrix, alleles, genes);
    // totals of all complexes, in parallel chunks
    let threads = threads.max(1);
    let chunk = n.div_ceil(threads).max(1);
    let mut totals: Vec<Totals> = vec![Totals::default(); n];
    std::thread::scope(|s| {
        for (k, part) in totals.chunks_mut(chunk).enumerate() {
            let fast = &fast;
            let matrix = &matrix;
            s.spawn(move || {
                let mut idx = vec![0usize; product.per_gene.len()];
                for (j, t) in part.iter_mut().enumerate() {
                    let i = k * chunk + j;
                    *t = match fast.totals(i, &mut idx) {
                        Some(t) => t,
                        None => {
                            let list = matrix.create(&product.alleles(i));
                            exact_totals(&list)
                        }
                    };
                }
            });
        }
    });
    let perc = crate::math::java_min(top_score_threshold * 5.0, 0.99);
    let java_split = crate::compat::java("LI6") && java_threads != 1 && n >= 10000;
    let tasks: Vec<Vec<usize>> = if java_split {
        let t = java_threads.max(1) as usize;
        let mut v: Vec<Vec<usize>> = vec![Vec::new(); t];
        for i in 0..n {
            v[i % t].push(i);
        }
        v
    } else {
        vec![(0..n).collect()]
    };
    let mut out = Vec::with_capacity(n);
    for task in tasks {
        let check = task.len() >= CULL_COMPLEX_COUNT;
        let mut culler = Culler { max: 0, perc };
        let mut culled = 0usize;
        for (k, &i) in task.iter().enumerate() {
            if check && k > 0 && k % CULL_COMPLEX_COUNT == 0 {
                crate::debug!("complexes({k}) processed, discard({culled}, {}%)", hmf_common::fmt::java_fixed(100.0 * culled as f64 / k as f64, 2));
            }
            let t = totals[i];
            if check && culler.cull(t.total()) {
                culled += 1;
                continue;
            }
            out.push((i as u32, t));
        }
        if check {
            crate::debug!("LI6: complex task of {} culled {}", task.len(), culled);
        }
    }
    Ok(out)
}

/// The totals ComplexCoverage.create gives from Java's floating point coverage
pub fn exact_totals(list: &[AlleleCoverage]) -> Totals {
    let mut u = 0;
    let mut s = 0.0;
    let mut w = 0.0;
    for c in list {
        u += c.unique;
        s += c.shared;
        w += c.wild;
    }
    Totals { unique: u, shared: crate::math::round_int(s), wild: crate::math::round_int(w) }
}

/// ComplexCoverageCalculator.calculateComplexCoverages on a list of complexes (Java's floating point)
pub fn calculate_list(rd: &RefData, complexes: &[Complex], fas: &[FragAlleles], java_threads: i32, top_score_threshold: f64) -> anyhow::Result<Vec<ComplexCoverage>> {
    check_thread_pool(complexes.len(), java_threads)?;
    let mut all: Vec<AId> = Vec::new();
    for c in complexes {
        for a in c {
            if !all.contains(a) {
                all.push(*a);
            }
        }
    }
    let matrix = Matrix::new(fas, &all);
    let n = complexes.len();
    let perc = crate::math::java_min(top_score_threshold * 5.0, 0.99);
    let java_split = crate::compat::java("LI6") && java_threads != 1 && n >= 10000;
    let tasks: Vec<Vec<usize>> = if java_split {
        let t = java_threads.max(1) as usize;
        let mut v: Vec<Vec<usize>> = vec![Vec::new(); t];
        for i in 0..n {
            v[i % t].push(i);
        }
        v
    } else {
        vec![(0..n).collect()]
    };
    let mut out = Vec::new();
    for task in tasks {
        let check = task.len() >= CULL_COMPLEX_COUNT;
        let mut culler = Culler { max: 0, perc };
        for &i in &task {
            let list = matrix.create(&complexes[i]);
            if check && culler.cull(total_fragments(&list)) {
                continue;
            }
            out.push(ComplexCoverage::create(&rd.gene_cache, &rd.alleles, list));
        }
    }
    Ok(out)
}

/// ComplexCoverageRanking
pub struct Ranker<'a> {
    pub rd: &'a RefData,
    pub max_score_diff: f64,
    pub gene_count: usize,
    pub freq_score_penalty: f64,
    /// log10(max(frequency, MIN_POPULATION_FREQUENCY)) per allele
    penalties: HashMap<AId, f64>,
    complexity: Complexity,
}

impl<'a> Ranker<'a> {
    pub fn new(rd: &'a RefData, max_score_diff: f64, gene_count: usize, freq_score_penalty: f64) -> Ranker<'a> {
        Ranker { rd, max_score_diff, gene_count, freq_score_penalty, penalties: HashMap::new(), complexity: Complexity::new() }
    }

    fn penalty(&mut self, a: AId) -> f64 {
        if let Some(p) = self.penalties.get(&a) {
            return *p;
        }
        let al = &self.rd.alleles;
        let f = self.rd.freqs.frequency(al.gene(a), al.name(a));
        let p = crate::math::log10(crate::math::java_max(f, crate::constants::MIN_POPULATION_FREQUENCY));
        self.penalties.insert(a, p);
        p
    }

    /// calcScore for a complex with these alleles (sorted by name) and total coverage
    pub fn score(&mut self, sorted: &[AId], total: i32, recovered: &[AId], sequences: &[&SeqLoci]) -> Scored {
        let al = &self.rd.alleles;
        let mut cohort = 0.0;
        for &a in sorted {
            let p = self.penalty(a);
            cohort += p;
            let g = self.rd.alleles.gene(a);
            if sorted.iter().filter(|x| self.rd.alleles.gene(**x) == g).count() == 1 {
                cohort += p;
            }
        }
        let recovered_count = sorted.iter().filter(|a| recovered.contains(a)).count() as i32;
        let mut wildcards = 0;
        if !sequences.is_empty() {
            for &a in sorted {
                if !al.has_wildcards(a) {
                    continue;
                }
                if let Some(s) = sequences.iter().find(|s| al.name(s.allele) == al.name(a)) {
                    wildcards += s.wildcard_count as i32;
                }
            }
        }
        let total_f = total as f64;
        let penalty_factor = 3.0 / self.gene_count as f64;
        let complexity = penalty_factor * complexity_count(&mut self.complexity, &self.rd.exon_sequences, &self.rd.alleles, sorted) as f64;
        let complexity_penalty = -complexity * crate::constants::SOLUTION_COMPLEXITY_PENALTY_WEIGHT * total_f;
        let cohort_penalty = penalty_factor * cohort * self.freq_score_penalty * total_f;
        let score = total_f + cohort_penalty + complexity_penalty - recovered_count as f64 * crate::constants::RECOVERY_SCORE_PENALTY * total_f;
        Scored {
            index: 0,
            totals: Totals::default(),
            score,
            cohort_total: cohort,
            cohort_penalty,
            complexity,
            complexity_penalty,
            recovered: recovered_count,
            wildcards,
        }
    }

    /// rankCandidates on scored items: the indices of the items kept, in rank order. `alleles_of` gives each
    /// item's alleles sorted by name.
    pub fn rank(&self, items: &[Scored], alleles_of: &dyn Fn(&Scored) -> Vec<AId>) -> Vec<usize> {
        if items.is_empty() {
            return Vec::new();
        }
        let al = &self.rd.alleles;
        let cmp = |a: &usize, b: &usize| -> Ordering {
            let (x, y) = (&items[*a], &items[*b]);
            if (x.score - y.score).abs() > 0.0001 {
                return if x.score < y.score { Ordering::Greater } else { Ordering::Less };
            }
            let ax = alleles_of(x);
            let ay = alleles_of(y);
            for i in 0..ax.len().min(ay.len()) {
                let c = al.compare(ax[i], ay[i]);
                if c != Ordering::Equal {
                    return c;
                }
            }
            Ordering::Equal
        };
        if self.max_score_diff == 0.0 {
            let mut idx: Vec<usize> = (0..items.len()).collect();
            hmf_common::java_sort::java_sort(&mut idx, cmp);
            return idx;
        }
        let top_score = items.iter().map(|x| x.score).fold(f64::NEG_INFINITY, crate::math::java_max);
        let top_cov = items.iter().map(|x| x.totals.total()).max().unwrap_or(0);
        let inclusion = top_score - self.max_score_diff * top_cov as f64;
        let base = top_score - 0.25 * top_cov as f64;
        let mut idx: Vec<usize> = (0..items.len()).filter(|i| items[*i].score >= base).collect();
        hmf_common::java_sort::java_sort(&mut idx, cmp);
        let mut out = Vec::new();
        for i in idx {
            if items[i].score >= inclusion || out.len() < 2 {
                out.push(i);
            } else {
                break;
            }
        }
        out
    }
}

/// Copies a score onto a ComplexCoverage
pub fn apply_score(c: &mut ComplexCoverage, s: &Scored) {
    c.score = s.score;
    c.cohort_frequency_total = s.cohort_total;
    c.cohort_frequency_penalty = s.cohort_penalty;
    c.complexity = s.complexity;
    c.complexity_penalty = s.complexity_penalty;
    c.recovered_count = s.recovered;
    c.wildcard_count = s.wildcards;
}

/// Scores and ranks ComplexCoverages (rankCandidates)
pub fn rank_coverages(ranker: &mut Ranker, mut covs: Vec<ComplexCoverage>, recovered: &[AId], sequences: &[&SeqLoci]) -> Vec<ComplexCoverage> {
    let mut items = Vec::with_capacity(covs.len());
    for (i, c) in covs.iter().enumerate() {
        let ids = c.allele_ids();
        let mut s = ranker.score(&ids, c.total, recovered, sequences);
        s.index = i as u32;
        s.totals = Totals { unique: c.unique, shared: c.shared, wild: c.wild };
        items.push(s);
    }
    for (c, s) in covs.iter_mut().zip(items.iter()) {
        apply_score(c, s);
    }
    let ids: Vec<Vec<AId>> = covs.iter().map(|c| c.allele_ids()).collect();
    let order = ranker.rank(&items, &|s: &Scored| ids[s.index as usize].clone());
    let mut slots: Vec<Option<ComplexCoverage>> = covs.into_iter().map(Some).collect();
    order.into_iter().map(|i| slots[i].take().unwrap()).collect()
}

/// solutionComplexity, with the exon chunks of each allele interned and their matches cached
pub struct Complexity {
    /// per allele: chunk id per exon index
    chunks: HashMap<AId, Vec<u32>>,
    /// per chunk id: sequence and wildcard flag
    seqs: Vec<ExonSequence>,
    index: HashMap<(usize, Vec<crate::tok::Tok>), u32>,
    matches: HashMap<(u32, u32), bool>,
}

impl Complexity {
    pub fn new() -> Complexity {
        Complexity { chunks: HashMap::new(), seqs: Vec::new(), index: HashMap::new(), matches: HashMap::new() }
    }

    fn of_allele(&mut self, exons: &HashMap<AId, Vec<ExonSequence>>, alleles: &Alleles, a: AId) -> Vec<u32> {
        if let Some(c) = self.chunks.get(&a) {
            return c.clone();
        }
        let f = alleles.four_digit(a);
        let ex = exons.get(&f).unwrap_or_else(|| panic!("java.lang.NullPointerException (no exon sequences for {})", alleles.name(f)));
        let mut ids = Vec::with_capacity(ex.len());
        for (i, e) in ex.iter().enumerate() {
            let key = (i, e.seq.clone());
            let id = match self.index.get(&key) {
                Some(id) => *id,
                None => {
                    let id = self.seqs.len() as u32;
                    self.seqs.push(e.clone());
                    self.index.insert(key, id);
                    id
                }
            };
            ids.push(id);
        }
        self.chunks.insert(a, ids.clone());
        ids
    }

    fn matches(&mut self, a: u32, b: u32) -> bool {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(m) = self.matches.get(&key) {
            return *m;
        }
        let (x, y) = (&self.seqs[a as usize].seq, &self.seqs[b as usize].seq);
        let m = x.len() == y.len() && x.iter().zip(y.iter()).all(|(p, q)| *p == WILD || *q == WILD || p == q);
        self.matches.insert(key, m);
        m
    }

}

/// ComplexCoverageRanking.solutionComplexity: the exon chunks of the alleles that differ from the others
pub fn complexity_count(cx: &mut Complexity, exons: &HashMap<AId, Vec<ExonSequence>>, alleles: &Alleles, sorted: &[AId]) -> i32 {
    if exons.is_empty() {
        return 0;
    }
    let per: Vec<Vec<u32>> = sorted.iter().map(|a| cx.of_allele(exons, alleles, *a)).collect();
    let max_exons = per.iter().map(|v| v.len()).max().unwrap_or(0);
    let mut count = 0;
    for e in 0..max_exons {
        let mut acids: Vec<u32> = per.iter().filter(|v| e < v.len()).map(|v| v[e]).collect();
        count += 1;
        // stable sort: the chunks with wildcards first
        acids.sort_by_key(|id| if cx.seqs[*id as usize].has_wildcards { 0 } else { 1 });
        for i in 1..acids.len() {
            let mut uniq = true;
            for j in 0..acids.len() {
                if i == j {
                    continue;
                }
                if j > i && !cx.seqs[acids[j] as usize].has_wildcards {
                    break;
                }
                if cx.seqs[acids[i] as usize].seq.len() != cx.seqs[acids[j] as usize].seq.len() {
                    continue;
                }
                if cx.matches(acids[i], acids[j]) {
                    uniq = false;
                    break;
                }
            }
            if uniq {
                count += 1;
            }
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hla::Alleles;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// the exact integer totals equal Java's floating point totals for every complex, and the half-way cases
    /// are the ones left to floating point
    #[test]
    fn fast_totals_match_floating_point() {
        let mut al = Alleles::new();
        let genes = [HlaGene::A, HlaGene::B, HlaGene::C];
        let mut per_gene_alleles: Vec<Vec<AId>> = Vec::new();
        for (gi, g) in ["A", "B", "C"].iter().enumerate() {
            let v: Vec<AId> = (1..=5).map(|k| al.request_four_digit(&format!("{g}*0{}:0{k}", gi + 1)).unwrap()).collect();
            per_gene_alleles.push(v);
        }
        for seed in 1..30u64 {
            let mut rng = Rng(seed);
            let mut fas = Vec::new();
            for f in 0..200usize {
                // mostly one gene, some two or three
                let ng = match rng.below(10) {
                    0 => 3,
                    1 | 2 => 2,
                    _ => 1,
                };
                let mut full = Vec::new();
                let mut wild = Vec::new();
                let start = rng.below(3) as usize;
                for k in 0..ng {
                    let g = (start + k) % 3;
                    for a in &per_gene_alleles[g] {
                        match rng.below(4) {
                            0 => full.push(*a),
                            1 => wild.push(*a),
                            2 if rng.below(3) == 0 => {
                                full.push(*a);
                                wild.push(*a);
                            }
                            _ => {}
                        }
                    }
                }
                fas.push(FragAlleles::new(&al, f, &full, &wild));
            }
            let per_gene: Vec<Vec<Complex>> = per_gene_alleles
                .iter()
                .map(|v| {
                    let mut c: Vec<Complex> = v.iter().map(|a| vec![*a]).collect();
                    for i in 0..v.len() {
                        for j in i + 1..v.len() {
                            c.push(vec![v[i], v[j]]);
                        }
                    }
                    c
                })
                .collect();
            let product = Product { per_gene };
            let all: Vec<AId> = per_gene_alleles.iter().flatten().copied().collect();
            let matrix = Matrix::new(&fas, &all);
            let fast = FastTotals::new(&product, &matrix, &al, &genes);
            let mut idx = vec![0; 3];
            for i in 0..product.len() {
                let exact = exact_totals(&matrix.create(&product.alleles(i)));
                match fast.totals(i, &mut idx) {
                    Some(t) => assert_eq!(t, exact, "seed {seed} complex {i}"),
                    None => {}
                }
            }
        }
    }
}
