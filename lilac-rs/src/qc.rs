//! The QC classes (AminoAcidQC, BamQC, CoverageQC, Haplotype, HaplotypeQC, SomaticVariantQC, LilacQC),
//! SomaticCodingCount and SolutionSummary.

use hmf_common::fmt::java_fixed;

use crate::config::Params;
use crate::constants::*;
use crate::coverage::{AlleleCoverage, ComplexCoverage};
use crate::evidence::PhasedEvidence;
use crate::fragment::{Fragment, Scope};
use crate::hla::{AId, Alleles, GeneSelector, HlaGene};
use crate::refdata::GeneCache;
use crate::seq::{SeqLoci, SequenceCount};
use crate::tok::{self, WILD};

pub struct AminoAcidQC {
    pub unused: i32,
    pub max_frags: i32,
}

impl AminoAcidQC {
    /// AminoAcidQC.create
    pub fn create(winners: &[&SeqLoci], hlay: Option<&[SeqLoci]>, counts: &SequenceCount, unmatched: &[Haplotype]) -> AminoAcidQC {
        let mut unused = 0;
        let mut largest = 0;
        for locus in counts.heterozygous_loci() {
            if unmatched.iter().any(|h| locus >= h.start && locus <= h.end) {
                continue;
            }
            let actual: Vec<tok::Tok> = winners.iter().filter(|s| (locus as usize) < s.len()).map(|s| s.at(locus as usize)).collect();
            let hy: Option<Vec<tok::Tok>> = hlay.map(|v| v.iter().filter(|s| (locus as usize) < s.len()).map(|s| s.at(locus as usize)).collect());
            if actual.contains(&WILD) {
                continue;
            }
            for (aa, count) in counts.get(locus).java_order() {
                if actual.contains(&aa) {
                    continue;
                }
                if hy.as_ref().is_some_and(|h| h.contains(&aa)) {
                    continue;
                }
                if (count as i32) < LOG_UNMATCHED_HAPLOTYPE_SUPPORT {
                    continue;
                }
                crate::warn!("  UNMATCHED_AMINO_ACID - amino acid sequence({} count={} locus={}) not in winning solution", tok::to_string(aa), count, locus);
                unused += 1;
                largest = largest.max(count as i32);
            }
        }
        AminoAcidQC { unused, max_frags: largest }
    }
}

pub struct BamQC {
    pub discarded_alignment_fragments: i32,
    pub discarded_indels: i32,
    pub discarded_indel_max_frags: i32,
    pub low_coverage: Vec<(HlaGene, i32)>,
}

impl BamQC {
    pub fn new(gc: &GeneCache, discarded_alignment_fragments: i32, discarded_indels: i32, discarded_indel_max_frags: i32, depth: &[(HlaGene, Vec<i32>)]) -> BamQC {
        let low_coverage = gc
            .gene_names
            .iter()
            .map(|g| {
                let d = depth.iter().find(|x| x.0 == *g).map(|x| x.1.as_slice()).unwrap_or(&[]);
                (*g, d.iter().filter(|x| (**x as f64) < WARN_LOW_COVERAGE_DEPTH).count() as i32)
            })
            .collect();
        BamQC { discarded_alignment_fragments, discarded_indels, discarded_indel_max_frags, low_coverage }
    }
    pub fn total_low_coverage(&self) -> i32 {
        self.low_coverage.iter().map(|x| x.1).sum()
    }
    fn body(&self, gc: &GeneCache) -> Vec<String> {
        let low: Vec<String> = gc
            .gene_names
            .iter()
            .filter(|g| !g.is_pseudo())
            .map(|g| format!("{}={}", g.short_name(), self.low_coverage.iter().find(|x| x.0 == *g).map(|x| x.1).unwrap_or(0)))
            .collect();
        vec![self.discarded_indels.to_string(), self.discarded_indel_max_frags.to_string(), self.discarded_alignment_fragments.to_string(), low.join(";")]
    }
}

pub struct CoverageQC {
    pub counts_by_gene: Vec<(HlaGene, i32)>,
    pub total_fragments: i32,
    pub fitted: i32,
    pub unmatched: i32,
    pub uninformative: i32,
    pub hla_y: i32,
    pub unique: i32,
    pub shared: i32,
    pub wildcard: i32,
    pub percent_unique: f64,
    pub percent_shared: f64,
    pub percent_wildcard: f64,
}

impl CoverageQC {
    #[allow(clippy::too_many_arguments)]
    pub fn new(counts_by_gene: Vec<(HlaGene, i32)>, total: i32, unmatched: i32, uninformative: i32, hla_y: i32, unique: i32, shared: i32, wildcard: i32) -> CoverageQC {
        let fitted = unique + shared + wildcard;
        CoverageQC {
            counts_by_gene,
            total_fragments: total,
            fitted,
            unmatched,
            uninformative,
            hla_y,
            unique,
            shared,
            wildcard,
            percent_unique: 1.0 * unique as f64 / fitted as f64,
            percent_shared: 1.0 * shared as f64 / fitted as f64,
            percent_wildcard: 1.0 * wildcard as f64 / fitted as f64,
        }
    }

    /// CoverageQC.create: counts of the scopes of the reference fragments
    pub fn create(gc: &GeneCache, alleles: &Alleles, frags: &[Fragment], subset: &[usize], winner: &ComplexCoverage) -> CoverageQC {
        let ids = winner.allele_ids();
        let mut counts = Vec::new();
        for g in &gc.gene_names {
            if g.is_pseudo() {
                continue;
            }
            let mut set: Vec<AId> = ids.iter().copied().filter(|a| alleles.gene(*a) == *g).collect();
            set.sort();
            set.dedup();
            counts.push((*g, set.len() as i32));
        }
        if counts.iter().any(|x| x.1 == 0) {
            let msg: Vec<String> = counts.iter().map(|x| format!("{} {} alleles", x.1, x.0.short_name())).collect();
            crate::warn!("  UNMATCHED ALLELE: {}", msg.join(", "));
        }
        let total = subset.len() as i32;
        let (mut sol, mut hy, mut uninf, mut unmatched) = (0, 0, 0, 0);
        for &i in subset {
            let f = &frags[i];
            match f.scope {
                Scope::HlaY => hy += 1,
                Scope::Solution => sol += 1,
                Scope::NoHetLoci => uninf += 1,
                s if s.is_unmatched() => unmatched += 1,
                s => crate::warn!("frag({} : {}) scope({}) unassigned", f.id(), f.read_info(), s.name()),
            }
        }
        let unassigned = total - sol - hy - uninf - unmatched;
        if unassigned != 0 {
            crate::warn!("fragment count inconsistency: total({total}) sol({sol}) unmatch({unmatched}) uninform({uninf}) hlaY({hy}) unassigned({unassigned})");
        }
        CoverageQC::new(counts, total, unmatched, uninf, hy, winner.unique, winner.shared, winner.wild)
    }

    fn body(&self, gc: &GeneCache) -> Vec<String> {
        let types: Vec<String> = gc
            .gene_names
            .iter()
            .filter(|g| !g.is_pseudo())
            .map(|g| format!("{}={}", g.short_name(), self.counts_by_gene.iter().find(|x| x.0 == *g).map(|x| x.1).unwrap_or(0)))
            .collect();
        vec![
            types.join(";"),
            self.total_fragments.to_string(),
            self.fitted.to_string(),
            self.unmatched.to_string(),
            self.uninformative.to_string(),
            self.hla_y.to_string(),
            java_fixed(self.percent_unique, 3),
            java_fixed(self.percent_shared, 3),
            java_fixed(self.percent_wildcard, 3),
        ]
    }
}

#[derive(Clone, Debug)]
pub struct Haplotype {
    pub start: i32,
    pub end: i32,
    pub phased_fragments: i32,
    pub haplotype: Vec<u8>,
    pub matching_fragments: i32,
}

impl Haplotype {
    /// Haplotype.fromString: startLocus,haplotype
    pub fn from_csv(line: &str) -> Haplotype {
        let items = hmf_common::java_text::split_all(line, ',');
        let start = hmf_common::java_text::parse_int(items[0]).unwrap();
        let h = items[1].as_bytes().to_vec();
        Haplotype { start, end: start + h.len() as i32 - 1, phased_fragments: 0, haplotype: h, matching_fragments: 0 }
    }

    /// Haplotype.create: the evidence sequence over its loci, filled between them with the most common
    /// amino acid
    pub fn create(loci: &[i32], seq: &str, count: i32, counts: &SequenceCount) -> Haplotype {
        let start = loci.iter().copied().min().unwrap_or(0);
        let end = loci.iter().copied().max().unwrap_or(0);
        let sparse = seq.as_bytes();
        let mut h: Vec<u8> = Vec::new();
        for locus in start..=end {
            if let Some(i) = loci.iter().position(|x| *x == locus) {
                if i >= sparse.len() {
                    panic!("java.lang.StringIndexOutOfBoundsException: index {i}, length {}", sparse.len());
                }
                h.push(sparse[i]);
            } else {
                tok::push_bytes(counts.max_count_sequence(locus), &mut h);
            }
        }
        Haplotype { start, end, phased_fragments: count, haplotype: h, matching_fragments: 0 }
    }

    /// contains: the other lies within this and agrees at every locus
    pub fn contains(&self, o: &Haplotype) -> bool {
        if o.start >= self.start && o.end <= self.end {
            for locus in o.start.max(self.start)..=o.end.min(self.end) {
                if self.char_at(locus) != o.char_at(locus) {
                    return false;
                }
            }
            return true;
        }
        false
    }

    fn char_at(&self, locus: i32) -> u8 {
        let i = (locus - self.start) as usize;
        if i >= self.haplotype.len() {
            panic!("java.lang.StringIndexOutOfBoundsException: index {i}, length {}", self.haplotype.len());
        }
        self.haplotype[i]
    }

    /// sequence(locus): "" outside, or where the haplotype is shorter than its range
    fn sequence(&self, locus: i32) -> Option<u8> {
        if locus < self.start || locus > self.end {
            return None;
        }
        let i = (locus - self.start) as usize;
        if i >= self.haplotype.len() {
            crate::error!("haplotype({}) invalid sequence request at locus({locus})", self.display());
            return None;
        }
        Some(self.haplotype[i])
    }

    pub fn display(&self) -> String {
        format!(
            "{} loci({} - {}) frags(supporting={} phased={})",
            String::from_utf8_lossy(&self.haplotype),
            self.start,
            self.end,
            self.matching_fragments,
            self.phased_fragments
        )
    }
}

/// HaplotypeQC.unmatchedHaplotype: the evidence sequences no winning (or HLA-Y) sequence explains
pub fn unmatched_haplotype(e: &PhasedEvidence, min_evidence: i32, winners: &[&SeqLoci], counts: &SequenceCount, hlay: Option<&[SeqLoci]>) -> Vec<Haplotype> {
    let mut out = Vec::new();
    for (seq, count) in &e.evidence {
        let b = seq.as_bytes();
        if winners.iter().any(|w| w.consistent_with(b, &e.loci)) {
            continue;
        }
        if hlay.is_some_and(|h| h.iter().any(|w| w.consistent_with(b, &e.loci))) {
            continue;
        }
        if *count < min_evidence {
            continue;
        }
        let h = Haplotype::create(&e.loci, seq, *count, counts);
        if (h.haplotype.len() as i32) < h.end - h.start + 1 {
            let loci: Vec<String> = e.loci.iter().map(|x| x.to_string()).collect();
            crate::error!("invalid haplotype({}) created from evidence({seq}) and amino-acid indices: {}", h.display(), loci.join(";"));
        }
        out.push(h);
    }
    out
}

pub struct HaplotypeQC {
    pub unused: i32,
    pub max_frags: i32,
    pub pon: i32,
    pub unmatched: Vec<Haplotype>,
}

fn pon_haplotypes() -> Vec<Haplotype> {
    hmf_common::java_text::lines_of(include_str!("../resources/haplotypes.csv")).iter().map(|l| Haplotype::from_csv(l)).collect()
}

impl HaplotypeQC {
    /// HaplotypeQC.create
    pub fn create(params: &Params, winners: &[&SeqLoci], hlay: Option<&[SeqLoci]>, evidence: &[PhasedEvidence], counts: &SequenceCount, unmatched_frags: &[&Fragment]) -> HaplotypeQC {
        let pon = pon_haplotypes();
        let mut all: Vec<Haplotype> = Vec::new();
        for e in evidence {
            all.extend(unmatched_haplotype(e, LOG_UNMATCHED_HAPLOTYPE_SUPPORT, winners, counts, hlay));
        }
        // sorted by matching fragments (all zero here: the order is kept), distinct by start, end and sequence
        let mut keys: std::collections::HashSet<(i32, i32, Vec<u8>)> = std::collections::HashSet::new();
        let mut distinct: Vec<Haplotype> = Vec::new();
        for h in all {
            if keys.insert((h.start, h.end, h.haplotype.clone())) {
                distinct.push(h);
            }
        }
        for h in distinct.iter_mut() {
            for f in unmatched_frags {
                if !(f.min_aa_locus() <= h.end && f.max_aa_locus() >= h.start) {
                    continue;
                }
                let mut matched = 0;
                let mut low = 0;
                let mut mismatch = false;
                for locus in h.start..=h.end {
                    let aa = f.aas.get(locus);
                    let faa = match aa {
                        Some(a) => a,
                        None => {
                            let l = f.low_qual_aa(locus, false, params);
                            if l == tok::EMPTY {
                                continue;
                            }
                            l
                        }
                    };
                    let fs = tok::to_string(faa);
                    let hs = h.sequence(locus).map(|c| (c as char).to_string()).unwrap_or_default();
                    if fs != hs {
                        mismatch = true;
                        break;
                    }
                    if aa.is_some() {
                        matched += 1;
                    } else {
                        low += 1;
                    }
                }
                if mismatch {
                    continue;
                }
                if matched + low > 0 {
                    h.matching_fragments += 1;
                }
            }
        }
        distinct.sort_by_key(|h| h.start);
        crate::debug!(
            "haplotype QC: {} distinct unmatched haplotypes, {} with {} or more fragments",
            distinct.len(),
            distinct.iter().filter(|h| h.matching_fragments >= LOG_UNMATCHED_HAPLOTYPE_SUPPORT).count(),
            LOG_UNMATCHED_HAPLOTYPE_SUPPORT
        );
        let (mut pon_count, mut unused, mut max_support) = (0, 0, 0);
        for h in &distinct {
            if h.matching_fragments < LOG_UNMATCHED_HAPLOTYPE_SUPPORT {
                continue;
            }
            if pon.iter().any(|p| p.contains(h)) {
                pon_count += 1;
                crate::info!("  UNMATCHED_PON_HAPLOTYPE - {}", h.display());
            } else {
                max_support = max_support.max(h.matching_fragments);
                unused += 1;
                crate::info!("  UNMATCHED_HAPLOTYPE {}", h.display());
            }
        }
        HaplotypeQC { unused, max_frags: max_support, pon: pon_count, unmatched: distinct }
    }
}

/// SomaticCodingCount
#[derive(Clone, Debug)]
pub struct SomaticCodingCount {
    pub allele: AId,
    pub inframe_indel: f64,
    pub missense: f64,
    pub nonsense: f64,
    pub splice: f64,
    pub synonymous: f64,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum CodingEffect {
    Nonsense,
    Splice,
    Missense,
    Synonymous,
    NoneEffect,
    Undefined,
}

impl CodingEffect {
    pub fn parse(s: &str) -> Option<CodingEffect> {
        Some(match s {
            "NONSENSE_OR_FRAMESHIFT" => CodingEffect::Nonsense,
            "SPLICE" => CodingEffect::Splice,
            "MISSENSE" => CodingEffect::Missense,
            "SYNONYMOUS" => CodingEffect::Synonymous,
            "NONE" => CodingEffect::NoneEffect,
            "UNDEFINED" => CodingEffect::Undefined,
            _ => return None,
        })
    }
    pub fn name(self) -> &'static str {
        match self {
            CodingEffect::Nonsense => "NONSENSE_OR_FRAMESHIFT",
            CodingEffect::Splice => "SPLICE",
            CodingEffect::Missense => "MISSENSE",
            CodingEffect::Synonymous => "SYNONYMOUS",
            CodingEffect::NoneEffect => "NONE",
            CodingEffect::Undefined => "UNDEFINED",
        }
    }
}

impl SomaticCodingCount {
    pub fn new(a: AId) -> SomaticCodingCount {
        SomaticCodingCount { allele: a, inframe_indel: 0.0, missense: 0.0, nonsense: 0.0, splice: 0.0, synonymous: 0.0 }
    }
    pub fn total(&self) -> f64 {
        self.inframe_indel + self.missense + self.nonsense + self.splice + self.synonymous
    }
    fn add(&mut self, indel: bool, e: CodingEffect, c: f64) {
        if indel && e == CodingEffect::Missense {
            self.inframe_indel += c;
            return;
        }
        match e {
            CodingEffect::Missense => self.missense += c,
            CodingEffect::Nonsense => self.nonsense += c,
            CodingEffect::Splice => self.splice += c,
            CodingEffect::Synonymous => self.synonymous += c,
            _ => {}
        }
    }
    /// SomaticCodingCount.addVariant: each variant allele's share to its first matching count
    pub fn add_variant(counts: &mut [SomaticCodingCount], alleles: &Alleles, indel: bool, effect: CodingEffect, variant_alleles: &[AId]) {
        let c = 1.0 / variant_alleles.len() as f64;
        for a in variant_alleles {
            if let Some(m) = counts.iter_mut().find(|x| alleles.name(x.allele) == alleles.name(*a)) {
                m.add(indel, effect, c);
            }
        }
    }
}

pub struct SomaticVariantQC {
    pub variant_count: i32,
    pub allele_count: f64,
}

impl SomaticVariantQC {
    pub fn matched(&self) -> i32 {
        crate::math::round_int(crate::math::java_min(self.allele_count, self.variant_count as f64))
    }
    pub fn unmatched(&self) -> i32 {
        self.variant_count - self.matched()
    }
    pub fn create(variant_count: i32, counts: &[SomaticCodingCount]) -> SomaticVariantQC {
        let total = hmf_common::java_math::double_stream_sum(counts.iter().map(|x| x.total()));
        let r = SomaticVariantQC { variant_count, allele_count: total };
        if r.unmatched() > 0 {
            crate::warn!("  UNASSIGNED_VARIANT - {} variants found but {} assigned", variant_count, hmf_common::fmt::java_double_to_string(total));
        }
        r
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum QcStatus {
    Pass,
    Fail,
    WarnUnmatchedAllele,
    WarnUnmatchedSomaticVariant,
    WarnUnmatchedHaplotype,
    WarnUnmatchedAminoAcid,
    WarnLowCoverage,
    WarnLowBaseQual,
    WarnUnmatchedIndel,
}

impl QcStatus {
    pub fn name(self) -> &'static str {
        match self {
            QcStatus::Pass => "PASS",
            QcStatus::Fail => "FAIL",
            QcStatus::WarnUnmatchedAllele => "WARN_UNMATCHED_ALLELE",
            QcStatus::WarnUnmatchedSomaticVariant => "WARN_UNMATCHED_SOMATIC_VARIANT",
            QcStatus::WarnUnmatchedHaplotype => "WARN_UNMATCHED_HAPLOTYPE",
            QcStatus::WarnUnmatchedAminoAcid => "WARN_UNMATCHED_AMINO_ACID",
            QcStatus::WarnLowCoverage => "WARN_LOW_COVERAGE",
            QcStatus::WarnLowBaseQual => "WARN_LOW_BASE_QUAL",
            QcStatus::WarnUnmatchedIndel => "WARN_UNMATCHED_INDEL",
        }
    }
}

pub struct LilacQC {
    pub genes: GeneSelector,
    /// in insertion order: a HashSet of enum constants iterates so under -XX:hashCode=2 (LI8)
    pub status: Vec<QcStatus>,
    pub score_margin: f64,
    pub next_solution: String,
    pub hla_y_allele: Option<String>,
    pub median_base_qual: Option<i32>,
    pub aa_qc: AminoAcidQC,
    pub bam_qc: BamQC,
    pub coverage_qc: CoverageQC,
    pub haplotype_qc: HaplotypeQC,
    pub somatic_qc: SomaticVariantQC,
}

impl LilacQC {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        params: &Params,
        genes: GeneSelector,
        score_margin: f64,
        next_solution: String,
        median_base_qual: Option<i32>,
        hla_y_allele: Option<String>,
        aa_qc: AminoAcidQC,
        bam_qc: BamQC,
        coverage_qc: CoverageQC,
        haplotype_qc: HaplotypeQC,
        somatic_qc: SomaticVariantQC,
    ) -> LilacQC {
        let mut q = LilacQC { genes, status: Vec::new(), score_margin, next_solution, hla_y_allele, median_base_qual, aa_qc, bam_qc, coverage_qc, haplotype_qc, somatic_qc };
        q.populate_status(params);
        q
    }

    fn populate_status(&mut self, p: &Params) {
        let add = |s: QcStatus, v: &mut Vec<QcStatus>| {
            if !v.contains(&s) {
                v.push(s);
            }
        };
        let mut st = Vec::new();
        if self.coverage_qc.total_fragments == 0 || self.bam_qc.total_low_coverage() >= p.fail_low_coverage_threshold {
            add(QcStatus::Fail, &mut st);
            self.status = st;
            return;
        }
        if let Some(m) = self.median_base_qual {
            if (m as f64) < WARN_LOW_BASE_QUAL_THRESHOLD {
                add(QcStatus::WarnLowBaseQual, &mut st);
            }
        }
        if self.bam_qc.total_low_coverage() as f64 >= p.warn_low_coverage_threshold {
            add(QcStatus::WarnLowCoverage, &mut st);
        }
        let hap_warn = self.coverage_qc.total_fragments as f64 * p.warn_unmatched_haplotype_support;
        if self.haplotype_qc.max_frags as f64 >= hap_warn {
            add(QcStatus::WarnUnmatchedHaplotype, &mut st);
        }
        if self.aa_qc.max_frags as f64 >= hap_warn {
            add(QcStatus::WarnUnmatchedAminoAcid, &mut st);
        }
        if self.coverage_qc.counts_by_gene.iter().any(|x| x.1 == 0) {
            add(QcStatus::WarnUnmatchedAllele, &mut st);
        }
        let indel_warn = self.coverage_qc.total_fragments as f64 * p.warn_indel_threshold;
        if self.bam_qc.discarded_indel_max_frags as f64 >= indel_warn {
            add(QcStatus::WarnUnmatchedIndel, &mut st);
        }
        if self.somatic_qc.unmatched() > 0 {
            add(QcStatus::WarnUnmatchedSomaticVariant, &mut st);
        }
        if st.is_empty() {
            add(QcStatus::Pass, &mut st);
        }
        self.status = st;
    }

    pub fn header() -> Vec<&'static str> {
        vec![
            "Genes",
            "Status",
            "ScoreMargin",
            "NextSolutionAlleles",
            "MedianBaseQuality",
            "HlaYAllele",
            "DiscardedIndels",
            "DiscardedIndelMaxFrags",
            "DiscardedAlignmentFragments",
            "LowCoverageBases",
            "geneTypes",
            "TotalFragments",
            "FittedFragments",
            "UnmatchedFragments",
            "UninformativeFragments",
            "HlaYFragments",
            "PercentUnique",
            "PercentShared",
            "PercentWildcard",
            "UnusedAminoAcids",
            "UnusedAminoAcidMaxFrags",
            "UnusedHaplotypes",
            "UnusedHaplotypeMaxFrags",
            "SomaticVariantsMatched",
            "SomaticVariantsUnmatched",
        ]
    }

    pub fn body(&self, gc: &GeneCache) -> Vec<String> {
        let mut v = vec![
            self.genes.name().to_string(),
            self.status.iter().map(|s| s.name()).collect::<Vec<_>>().join(";"),
            java_fixed(self.score_margin, 3),
            self.next_solution.clone(),
            self.median_base_qual.map(|x| x.to_string()).unwrap_or_else(|| "-".to_string()),
            self.hla_y_allele.clone().unwrap_or_else(|| "NONE".to_string()),
        ];
        v.extend(self.bam_qc.body(gc));
        v.extend(self.coverage_qc.body(gc));
        v.push(self.aa_qc.unused.to_string());
        v.push(self.aa_qc.max_frags.to_string());
        v.push(self.haplotype_qc.unused.to_string());
        v.push(self.haplotype_qc.max_frags.to_string());
        v.push(self.somatic_qc.matched().to_string());
        v.push(self.somatic_qc.unmatched().to_string());
        v
    }

    pub fn log(&self, gc: &GeneCache, sample: &str) {
        let h = Self::header();
        let b = self.body(gc);
        let items: Vec<String> = h.iter().zip(b.iter()).map(|(x, y)| format!("{x}={y}")).collect();
        crate::info!("{sample} QC Stats: {}", items.join(","));
    }
}

/// SolutionSummary
pub struct SolutionSummary {
    pub genes: GeneSelector,
    pub reference: Option<ComplexCoverage>,
    pub tumor: Option<ComplexCoverage>,
    pub tumor_copy_number: Vec<f64>,
    pub coding: Vec<SomaticCodingCount>,
    pub rna: Option<ComplexCoverage>,
}

pub const LILAC_ALLELE_HEADER: &str = "Genes\tAllele\tRefTotal\tRefUnique\tRefShared\tRefWild\tTumorTotal\tTumorUnique\tTumorShared\tTumorWild\tRnaTotal\tRnaUnique\tRnaShared\tRnaWild\tTumorCopyNumber\tSomaticMissense\tSomaticNonsenseOrFrameshift\tSomaticSplice\tSomaticSynonymous\tSomaticInframeIndel";

impl SolutionSummary {
    /// SolutionSummary.create: the coding counts sorted by HlaAllele.compareTo
    pub fn create(alleles: &Alleles, genes: GeneSelector, reference: ComplexCoverage, tumor: ComplexCoverage, tumor_copy_number: Vec<f64>, coding: &[SomaticCodingCount], rna: ComplexCoverage) -> SolutionSummary {
        let mut sorted = coding.to_vec();
        hmf_common::java_sort::java_sort(&mut sorted, |a, b| alleles.compare(a.allele, b.allele));
        SolutionSummary { genes, reference: Some(reference), tumor: Some(tumor), tumor_copy_number, coding: sorted, rna: Some(rna) }
    }

    /// the lines of lilac.tsv for this gene group
    pub fn lines(&self, alleles: &Alleles) -> Vec<String> {
        let Some(r) = &self.reference else { return Vec::new() };
        let t = self.tumor.as_ref();
        let n = if !r.alleles.is_empty() {
            r.alleles.len()
        } else if let Some(t) = t.filter(|t| !t.alleles.is_empty()) {
            t.alleles.len()
        } else {
            0
        };
        let mut out = Vec::new();
        for i in 0..n {
            let ref_allele = if !r.alleles.is_empty() { r.alleles[i].allele } else { t.unwrap().alleles[i].allele };
            let none = AlleleCoverage::new(ref_allele, 0, 0.0, 0.0);
            let rc = if !r.alleles.is_empty() { &r.alleles[i] } else { &none };
            let tc = match t {
                Some(t) if !t.alleles.is_empty() => &t.alleles[i],
                _ => &none,
            };
            let rna = match &self.rna {
                Some(x) if !x.alleles.is_empty() => &x.alleles[i],
                _ => &none,
            };
            let cn = self.tumor_copy_number[i];
            // LI9: Java pairs the i-th coding count (sorted by HlaAllele.compareTo) with the i-th allele (sorted
            // by name); the correction takes the allele's own count
            if alleles.name(self.coding[i].allele) != alleles.name(ref_allele) {
                crate::debug!("LI9: row {} allele {} paired with the coding counts of {}", i, alleles.name(ref_allele), alleles.name(self.coding[i].allele));
            }
            let cc = if crate::compat::java("LI9") {
                &self.coding[i]
            } else {
                self.coding.iter().find(|c| alleles.name(c.allele) == alleles.name(ref_allele)).unwrap_or(&self.coding[i])
            };
            let r = |x: f64| crate::math::round_int(x).to_string();
            out.push(
                [
                    self.genes.name().to_string(),
                    alleles.name(ref_allele).to_string(),
                    r(rc.total),
                    rc.unique.to_string(),
                    r(rc.shared),
                    r(rc.wild),
                    r(tc.total),
                    tc.unique.to_string(),
                    r(tc.shared),
                    r(tc.wild),
                    r(rna.total),
                    rna.unique.to_string(),
                    r(rna.shared),
                    r(rna.wild),
                    java_fixed(cn, 2),
                    java_fixed(cc.missense, 1),
                    java_fixed(cc.nonsense, 1),
                    java_fixed(cc.splice, 1),
                    java_fixed(cc.synonymous, 1),
                    java_fixed(cc.inframe_indel, 1),
                ]
                .join("\t"),
            );
        }
        out
    }
}
