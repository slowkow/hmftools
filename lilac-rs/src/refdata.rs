//! ReferenceData, GeneCache, CohortFrequency, HlaContextFactory, ExpectedAlleles and the gene-pair filters of
//! NucleotideGeneEnrichment: everything loaded per gene group.

use std::collections::HashMap;

use anyhow::{bail, Context, Result};
use hmf_common::java_hash::{bucket_index, grown_capacity};
use hmf_common::java_text::{parse_double, parse_int, split};

use crate::config::LilacConfig;
use crate::hla::{AId, Alleles, GeneSelector, HlaGene};
use crate::indel::Indel;
use crate::seq::{self, ExonSequence, SeqLoci};

#[derive(Clone, Debug)]
pub struct Exon {
    pub start: i32,
    pub end: i32,
    pub rank: i32,
    pub phase: i32,
    pub end_phase: i32,
}

#[derive(Clone, Debug)]
pub struct Transcript {
    pub trans_id: i32,
    pub trans_name: String,
    pub gene_id: String,
    pub strand: i8,
    pub trans_start: i32,
    pub trans_end: i32,
    pub coding_start: i32,
    pub coding_end: i32,
    pub exons: Vec<Exon>,
}

impl Transcript {
    pub fn pos_strand(&self) -> bool {
        self.strand == 1
    }

    /// TranscriptUtils.calcCodingBases(...).CodingBases for a position within the coding region
    pub fn coding_bases(&self, position: i32) -> i32 {
        let coding_start = self.coding_start;
        let coding_end = self.coding_end;
        let pos_in_coding = position >= coding_start && position <= coding_end;
        let mut in_coding = false;
        let mut ended = false;
        let mut is_exonic = false;
        let mut total = 0;
        let mut coding_bases = 0;
        for exon in &self.exons {
            if !is_exonic && position >= exon.start && position <= exon.end {
                is_exonic = true;
            }
            if !in_coding {
                if exon.end >= coding_start {
                    in_coding = true;
                    total += exon.end - coding_start + 1;
                    if position >= coding_start {
                        if position < exon.end {
                            coding_bases += position - coding_start + 1;
                        } else {
                            coding_bases += exon.end - coding_start + 1;
                        }
                    }
                }
            } else if !ended {
                if exon.start > coding_end {
                    ended = true;
                } else if exon.end >= coding_end {
                    ended = true;
                    total += coding_end - exon.start + 1;
                    if position >= exon.start {
                        if position < coding_end {
                            coding_bases += position - exon.start + 1;
                        } else {
                            coding_bases += coding_end - exon.start + 1;
                        }
                    }
                } else {
                    total += exon.end - exon.start + 1;
                    if position >= exon.start {
                        if position < exon.end {
                            coding_bases += position - exon.start + 1;
                        } else {
                            coding_bases += exon.end - exon.start + 1;
                        }
                    }
                }
                if ended {
                    break;
                }
            }
        }
        if self.strand == -1 {
            if pos_in_coding && is_exonic {
                coding_bases = total - coding_bases + 1;
            } else {
                coding_bases = total - coding_bases;
            }
        }
        coding_bases
    }
}

/// ReferenceData.loadHlaTranscripts: the transcripts in the jar, in file order (Java puts them in a HashMap
/// and sorts the genes, so the order does not matter)
pub fn load_hla_transcripts(v38: bool, genes: Option<GeneSelector>) -> Vec<(HlaGene, Transcript)> {
    let text = if v38 {
        include_str!("../resources/hla_transcripts_v38.csv")
    } else {
        include_str!("../resources/hla_transcripts_v37.csv")
    };
    let lines = hmf_common::java_text::lines_of(text);
    let idx = hmf_common::java_text::fields_index_map(&lines[0], ',');
    let f = |n: &str| idx[n];
    let trans_name_index = idx.get("TransName").or(idx.get("Trans")).copied().unwrap();
    let mut out: Vec<(HlaGene, Transcript)> = Vec::new();
    let mut current_gene = String::new();
    let mut current: Option<(Option<HlaGene>, Transcript)> = None;
    let finish = |c: Option<(Option<HlaGene>, Transcript)>, out: &mut Vec<(HlaGene, Transcript)>| {
        if let Some((Some(g), t)) = c {
            if genes.is_none_or(|s| s.contains(g)) {
                out.push((g, t));
            }
        }
    };
    for line in &lines[1..] {
        let items = split(line, ',');
        let gene_id = items[f("GeneId")];
        let trans_id = parse_int(items[f("TransId")]).unwrap();
        if gene_id != current_gene {
            finish(current.take(), &mut out);
            current_gene = gene_id.to_string();
            let gene = HlaGene::from_string(items[f("GeneName")]);
            current = Some((
                gene,
                Transcript {
                    trans_id,
                    trans_name: items[trans_name_index].to_string(),
                    gene_id: gene_id.to_string(),
                    strand: parse_int(items[f("Strand")]).unwrap() as i8,
                    trans_start: parse_int(items[f("TransStart")]).unwrap(),
                    trans_end: parse_int(items[f("TransEnd")]).unwrap(),
                    coding_start: parse_int(items[f("CodingStart")]).unwrap(),
                    coding_end: parse_int(items[f("CodingEnd")]).unwrap(),
                    exons: Vec::new(),
                },
            ));
        }
        current.as_mut().unwrap().1.exons.push(Exon {
            start: parse_int(items[f("ExonStart")]).unwrap(),
            end: parse_int(items[f("ExonEnd")]).unwrap(),
            rank: parse_int(items[f("ExonRank")]).unwrap(),
            phase: parse_int(items[f("ExonPhase")]).unwrap(),
            end_phase: parse_int(items[f("ExonEndPhase")]).unwrap(),
        });
    }
    finish(current.take(), &mut out);
    out
}

/// ReferenceData.trimDrb1Transcripts: HLA-DRB1 restricted to its last two exons, with a coding start that
/// leaves a multiple of three coding bases
fn trim_drb1(t: &Transcript) -> Transcript {
    let exons: Vec<Exon> = t.exons[t.exons.len() - 2..].to_vec();
    let mut nuc_count = 0;
    for e in &exons {
        let start = t.coding_start.max(e.start);
        let end = t.coding_end.min(e.end);
        nuc_count += end - start + 1;
    }
    let mut trimmed = t.clone();
    trimmed.coding_start = exons[0].start + (nuc_count % 3);
    trimmed.exons = exons;
    trimmed
}

/// GeneCache
#[derive(Clone, Debug)]
pub struct GeneCache {
    /// sorted by enum order
    pub gene_names: Vec<HlaGene>,
    /// in gene_names order
    pub transcripts: Vec<Transcript>,
    pub expect_allele_count: usize,
    pub aa_boundaries: HashMap<HlaGene, Vec<i32>>,
    pub nuc_boundaries: HashMap<HlaGene, Vec<i32>>,
    pub nuc_lengths: HashMap<HlaGene, i32>,
    pub max_common_aa_boundary: i32,
}

impl GeneCache {
    pub fn new(map: &[(HlaGene, Transcript)]) -> GeneCache {
        let mut gene_names: Vec<HlaGene> = map.iter().map(|x| x.0).collect();
        gene_names.sort();
        gene_names.dedup();
        let get = |g: HlaGene| map.iter().rev().find(|x| x.0 == g).map(|x| x.1.clone()).unwrap();
        let transcripts: Vec<Transcript> = gene_names.iter().map(|g| get(*g)).collect();
        let mut aa_boundaries = HashMap::new();
        let mut nuc_boundaries = HashMap::new();
        let mut nuc_lengths = HashMap::new();
        // in insertion order for findMaxCommonAminoAcidBoundary (any order gives the same answer)
        let mut aa_ordered: Vec<Vec<i32>> = Vec::new();
        for (g, t) in gene_names.iter().zip(transcripts.iter()) {
            let (aa, nuc, len) = exon_boundaries(t);
            aa_ordered.push(aa.clone());
            aa_boundaries.insert(*g, aa);
            nuc_boundaries.insert(*g, nuc);
            nuc_lengths.insert(*g, len);
        }
        let max_common_aa_boundary = if aa_ordered.len() == 1 {
            *aa_ordered[0].last().unwrap_or(&-1)
        } else {
            let mut m = -1;
            for b in &aa_ordered[0] {
                if aa_ordered[1..].iter().all(|s| s.contains(b)) {
                    m = *b;
                }
            }
            m
        };
        GeneCache {
            expect_allele_count: gene_names.len() * 2,
            gene_names,
            transcripts,
            aa_boundaries,
            nuc_boundaries,
            nuc_lengths,
            max_common_aa_boundary,
        }
    }

    pub fn transcript(&self, g: HlaGene) -> Option<&Transcript> {
        self.gene_names.iter().position(|x| *x == g).map(|i| &self.transcripts[i])
    }

    pub fn aa_boundaries(&self, g: HlaGene) -> &[i32] {
        self.aa_boundaries.get(&g).map(|v| v.as_slice()).unwrap_or(&[])
    }
    pub fn nuc_boundaries(&self, g: HlaGene) -> &[i32] {
        self.nuc_boundaries.get(&g).map(|v| v.as_slice()).unwrap_or(&[])
    }
    pub fn nuc_length(&self, g: HlaGene) -> i32 {
        *self.nuc_lengths.get(&g).unwrap_or(&0)
    }

    /// LilacUtils.calcNucelotideLocus: the zero-based coding base index of a position in the first transcript
    /// whose coding region holds it, -1 for none
    pub fn nucleotide_locus(&self, position: i32) -> i32 {
        nucleotide_locus(self.transcripts.iter(), position)
    }
}

pub fn nucleotide_locus<'a>(transcripts: impl Iterator<Item = &'a Transcript>, position: i32) -> i32 {
    for t in transcripts {
        if position < t.coding_start || position > t.coding_end {
            continue;
        }
        return t.coding_bases(position) - 1;
    }
    -1
}

/// GeneCache.setExonBoundaryValues: amino acid exon boundaries, nucleotide boundaries, coding length
fn exon_boundaries(t: &Transcript) -> (Vec<i32>, Vec<i32>, i32) {
    let fwd = t.pos_strand();
    let n = t.exons.len() as i32;
    let mut idx: i32 = if fwd { 0 } else { n - 1 };
    let mut total = 0;
    let mut aa = Vec::new();
    let mut nuc = Vec::new();
    while idx >= 0 && idx < n {
        let e = &t.exons[idx as usize];
        let mut within = false;
        if fwd {
            if e.start > t.coding_end {
                break;
            }
            if e.end >= t.coding_start {
                within = true;
            }
        } else {
            if e.end < t.coding_start {
                break;
            }
            if e.start <= t.coding_end {
                within = true;
            }
        }
        if within {
            total += e.end.min(t.coding_end) - e.start.max(t.coding_start) + 1;
            let boundary = if fwd { e.end <= t.coding_end } else { e.start >= t.coding_start };
            if boundary {
                let a = total / 3;
                aa.push(a);
                nuc.push(a * 3);
                nuc.push(a * 3 + 1);
                nuc.push(a * 3 + 2);
            }
        }
        idx += if fwd { 1 } else { -1 };
    }
    (aa, nuc, total)
}

/// CohortFrequency: allele name to cohort frequency, kept in the HashMap's insertion order
#[derive(Clone, Debug, Default)]
pub struct CohortFrequency {
    /// (name, frequency, removed), first-insertion order
    entries: Vec<(String, f64, bool)>,
    index: HashMap<String, usize>,
}

impl CohortFrequency {
    pub fn load(genes: GeneSelector, path: &str) -> Result<CohortFrequency> {
        let mut cf = CohortFrequency::default();
        if path.is_empty() {
            return Ok(cf);
        }
        let lines = match hmf_common::java_text::read_lines(path) {
            Ok(l) => l,
            Err(e) => {
                crate::error!("failed to read cohort allele frequency file({path}): {e}");
                return Ok(cf);
            }
        };
        for line in &lines {
            if line.starts_with("Allele") {
                continue;
            }
            let items = split(line, ',');
            if items.len() != 2 {
                crate::error!("invalid allele frequency record: {line}");
                return Ok(cf);
            }
            let Some(freq) = parse_double(items[1]) else {
                bail!("java.lang.NumberFormatException: For input string: \"{}\"", items[1]);
            };
            let Some(p) = crate::hla::parse_allele(items[0]) else {
                bail!("invalid allele in frequency file: {}", items[0]);
            };
            if !genes.contains(p.gene) {
                continue;
            }
            let name = crate::hla::allele_name(p.gene, &p.group, &p.protein, &p.synonymous, &p.synonymous_non_coding);
            match cf.index.get(&name) {
                Some(i) => cf.entries[*i].1 = freq,
                None => {
                    cf.index.insert(name.clone(), cf.entries.len());
                    cf.entries.push((name, freq, false));
                }
            }
        }
        crate::info!("loaded {} allele frequencies from file({path})", cf.len());
        Ok(cf)
    }

    pub fn len(&self) -> usize {
        self.entries.iter().filter(|e| !e.2).count()
    }

    /// put(allele, frequency)
    pub fn put(&mut self, name: &str, freq: f64) {
        match self.index.get(name) {
            Some(i) => {
                self.entries[*i].1 = freq;
                self.entries[*i].2 = false;
            }
            None => {
                self.index.insert(name.to_string(), self.entries.len());
                self.entries.push((name.to_string(), freq, false));
            }
        }
    }

    pub fn remove(&mut self, name: &str) {
        if let Some(i) = self.index.remove(name) {
            self.entries[i].2 = true;
        }
    }

    pub fn contains(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }

    /// getAlleleFrequency
    pub fn frequency(&self, gene: HlaGene, name: &str) -> f64 {
        if !gene.has_frequencies() {
            return crate::constants::RARE_ALLELES_FREQ_CUTOFF;
        }
        self.index.get(name).map(|i| self.entries[*i].1).unwrap_or(0.0)
    }

    /// the entries in the HashMap's iteration order (the table keeps the size it grew to)
    pub fn java_order(&self) -> Vec<(&str, f64)> {
        let cap = grown_capacity(16, self.entries.len());
        let mut v: Vec<(usize, &str, f64)> = self
            .entries
            .iter()
            .filter(|e| !e.2)
            .map(|e| (bucket_index(hmf_common::java_hash::java_string_hash(&e.0), cap), e.0.as_str(), e.1))
            .collect();
        v.sort_by_key(|x| x.0);
        v.into_iter().map(|x| (x.1, x.2)).collect()
    }
}

/// ExpectedAlleles
#[derive(Clone, Debug)]
pub struct ExpectedAlleles {
    pub default_count: i32,
    pub counts: Option<Vec<i32>>,
}

impl ExpectedAlleles {
    fn at(&self, locus: i32) -> i32 {
        match &self.counts {
            Some(c) if locus >= 0 && (locus as usize) < c.len() => c[locus as usize],
            _ => self.default_count,
        }
    }
    /// expectedAlleles(loci): the minimum over the loci, 0 for none
    pub fn expected(&self, loci: &[i32]) -> i32 {
        loci.iter().map(|l| self.at(*l)).min().unwrap_or(0)
    }

    /// ExpectedAlleles.expectedAlleles(geneCount, filters). Index i of the array holds the count for locus
    /// i + 1 but is looked up for locus i, as in Java (an apparent off-by-one, kept: intent unclear).
    pub fn create(gene_count: i32, filters: &std::collections::BTreeMap<i32, i32>) -> ExpectedAlleles {
        if filters.is_empty() {
            return ExpectedAlleles { default_count: 2 * gene_count, counts: None };
        }
        let max = *filters.keys().last().unwrap();
        let size = (max - 1).max(0) as usize;
        if max - 1 < 0 {
            // new int[negative] throws in Java
            panic!("java.lang.NegativeArraySizeException: {}", max - 1);
        }
        let mut counts = vec![0; size];
        for (i, c) in counts.iter_mut().enumerate() {
            let locus = i as i32 + 1;
            let n: i32 = filters.range(..=locus).map(|x| *x.1).sum();
            *c = 2 * (gene_count - n);
        }
        let n: i32 = filters.values().sum();
        ExpectedAlleles { default_count: 2 * (gene_count - n), counts: Some(counts) }
    }
}

/// NucleotideGeneEnrichment's gene-pair filters: for each pair of genes, the first amino acid exon boundary
/// that only one of them has
#[derive(Clone, Debug)]
pub struct GeneEnrichment {
    pub genes: Vec<HlaGene>,
    pairs: Vec<(HlaGene, HlaGene, Option<i32>)>,
}

impl GeneEnrichment {
    pub fn create(boundaries: &HashMap<HlaGene, Vec<i32>>) -> Option<GeneEnrichment> {
        let genes: Vec<HlaGene> = if boundaries.contains_key(&HlaGene::A) {
            GeneSelector::MhcClass1.genes().to_vec()
        } else if boundaries.contains_key(&HlaGene::Drb3) {
            GeneSelector::HlaDrb.genes().to_vec()
        } else {
            return None;
        };
        let empty: Vec<i32> = Vec::new();
        let mut pairs = Vec::new();
        for i in 0..genes.len().saturating_sub(1) {
            for j in i + 1..genes.len() {
                let b1 = boundaries.get(&genes[i]).unwrap_or(&empty);
                let b2 = boundaries.get(&genes[j]).unwrap_or(&empty);
                let unique: Vec<i32> = b1.iter().filter(|x| !b2.contains(x)).chain(b2.iter().filter(|x| !b1.contains(x))).copied().collect();
                pairs.push((genes[i], genes[j], unique.iter().min().copied()));
            }
        }
        Some(GeneEnrichment { genes, pairs })
    }

    /// getFilter(gene1, gene2)
    pub fn filter(&self, g1: HlaGene, g2: HlaGene) -> Option<i32> {
        self.pairs.iter().find(|p| (p.0 == g1 && p.1 == g2) || (p.0 == g2 && p.1 == g1)).and_then(|p| p.2)
    }

    /// getFilters(gene): boundary to the number of other genes it separates from this one
    pub fn filters(&self, gene: HlaGene) -> std::collections::BTreeMap<i32, i32> {
        let mut m = std::collections::BTreeMap::new();
        for other in &self.genes {
            if *other == gene {
                continue;
            }
            if let Some(f) = self.filter(gene, *other) {
                *m.entry(f).or_insert(0) += 1;
            }
        }
        m
    }
}

/// HlaContext
#[derive(Clone, Debug)]
pub struct HlaContext {
    pub gene: HlaGene,
    pub aa_boundaries: Vec<i32>,
    pub expected: ExpectedAlleles,
}

/// The two sequence files, read once for all gene groups (Java reads them again for each group)
pub struct SequenceFiles {
    pub nucleotides: Vec<String>,
    pub amino_acids: Vec<String>,
}

impl SequenceFiles {
    pub fn load(resource_dir: &str) -> Result<SequenceFiles> {
        let n = format!("{resource_dir}hla_ref_nucleotide_sequences.csv");
        let a = format!("{resource_dir}hla_ref_aminoacid_sequences.csv");
        let nucleotides = hmf_common::java_text::read_lines(&n).with_context(|| format!("failed to load ref sequence data from file({n})"))?;
        let amino_acids = hmf_common::java_text::read_lines(&a).with_context(|| format!("failed to load ref sequence data from file({a})"))?;
        Ok(SequenceFiles { nucleotides, amino_acids })
    }
}

pub const CLASS_1_EXCLUDED_ALLELES: &[&str] =
    &["A*31:135", "A*33:191", "A*02:783", "B*07:282", "A*30:205", "A*30:207", "A*30:225", "A*30:228", "A*01:81", "A*01:237", "B*40:278"];
pub const HLA_DRB1_EXCLUDED_ALLELES: &[&str] = &["DRB1*14:242", "DRB1*14:141", "DRB1*14:262", "DRB1*12:57"];
pub const STOP_LOSS_ON_C_ALLELE: &str = "C*04:09N";
pub const DEFLATE_TEMPLATE: &str = "A*01:01";

pub const V37_HLA_REGION: (&str, i32, i32) = ("6", 29908037, 33055978);
pub const V38_HLA_REGION: (&str, i32, i32) = ("chr6", 29940260, 33090696);

/// The reference bases of the HLA region (ReferenceData.HLA_REF_BASES)
#[derive(Clone, Debug)]
pub struct RefBases {
    pub start: i32,
    pub bases: Vec<u8>,
}

impl RefBases {
    pub fn load(ref_genome: &str, v38: bool) -> Result<RefBases> {
        let (chr, start, end) = if v38 { V38_HLA_REGION } else { V37_HLA_REGION };
        let fa = hmf_common::refgenome::RefGenome::open(ref_genome).with_context(|| format!("Failed to load ref genome: {ref_genome}"))?;
        let bases = fa.bases(chr, start, end)?;
        Ok(RefBases { start, bases })
    }
    /// refBases(startPos, endPos)
    pub fn get(&self, start: i32, end: i32) -> Result<String> {
        let s = start - self.start;
        let e = end - self.start + 1;
        if s < 0 || e < s || e as usize > self.bases.len() {
            bail!("java.lang.StringIndexOutOfBoundsException: begin {s}, end {e}, length {}", self.bases.len());
        }
        Ok(String::from_utf8_lossy(&self.bases[s as usize..e as usize]).into_owned())
    }
}

/// ReferenceData of one gene group
pub struct RefData {
    pub genes: GeneSelector,
    pub alleles: Alleles,
    pub gene_cache: GeneCache,
    pub nuc_seqs: Vec<SeqLoci>,
    pub aa_seqs: Vec<SeqLoci>,
    /// AminoAcidSequencesWithInserts and WithDeletes (HLA-Y included)
    pub aa_with_inserts: Vec<SeqLoci>,
    pub aa_with_deletes: Vec<SeqLoci>,
    pub hlay_nuc_seqs: Option<Vec<SeqLoci>>,
    pub hlay_aa_seqs: Option<Vec<SeqLoci>>,
    /// four-digit allele to its exon sequences
    pub exon_sequences: HashMap<AId, Vec<ExonSequence>>,
    pub common_alleles: Vec<AId>,
    /// KnownStopLossIndelAlleles
    pub stop_loss_alleles: Vec<(Indel, AId)>,
    /// STOP_LOSS_ON_C_INDEL
    pub stop_loss_indel: Option<Indel>,
    pub freqs: CohortFrequency,
    pub pon_indels: Vec<Indel>,
    pub excluded: Vec<&'static str>,
    pub gene_exon_boundaries: HashMap<HlaGene, Vec<i32>>,
    pub gene_enrichment: Option<GeneEnrichment>,
    pub contexts: Vec<HlaContext>,
    /// the config's actual and restricted alleles, as cached four-digit alleles after loading
    pub actual_alleles: Vec<AId>,
    pub restricted_alleles: Vec<AId>,
    pub transcripts: Vec<(HlaGene, Transcript)>,
}

impl RefData {
    /// new ReferenceData(resourceDir, config) for the gene group
    pub fn new(config: &LilacConfig, genes: GeneSelector) -> Result<RefData> {
        let mut transcripts = load_hla_transcripts(config.v38, Some(genes));
        if genes == GeneSelector::HlaDrb1 {
            for t in transcripts.iter_mut() {
                if t.0 == HlaGene::Drb1 {
                    t.1 = trim_drb1(&t.1);
                }
            }
        }
        let gene_cache = GeneCache::new(&transcripts);
        let mut excluded = Vec::new();
        if genes.covers_mhc_class1() {
            excluded.extend_from_slice(CLASS_1_EXCLUDED_ALLELES);
        }
        if genes.contains(HlaGene::Drb1) {
            excluded.extend_from_slice(HLA_DRB1_EXCLUDED_ALLELES);
        }
        let mut gene_exon_boundaries = HashMap::new();
        for g in &gene_cache.gene_names {
            gene_exon_boundaries.insert(*g, gene_cache.aa_boundaries(*g).to_vec());
        }
        let gene_enrichment = GeneEnrichment::create(&gene_exon_boundaries);
        let gene_count = genes.genes().len() as i32;
        let contexts = genes
            .genes()
            .iter()
            .map(|g| {
                let filters = gene_enrichment.as_ref().map(|e| e.filters(*g)).unwrap_or_default();
                HlaContext {
                    gene: *g,
                    aa_boundaries: gene_exon_boundaries.get(g).cloned().unwrap_or_default(),
                    expected: ExpectedAlleles::create(gene_count, &filters),
                }
            })
            .collect();
        let freq_file = if config.resource_dir.is_empty() { String::new() } else { format!("{}lilac_allele_frequencies.csv", config.resource_dir) };
        let freqs = CohortFrequency::load(genes, &freq_file)?;
        Ok(RefData {
            genes,
            alleles: Alleles::new(),
            gene_cache,
            nuc_seqs: Vec::new(),
            aa_seqs: Vec::new(),
            aa_with_inserts: Vec::new(),
            aa_with_deletes: Vec::new(),
            hlay_nuc_seqs: if genes.covers_mhc_class1() { Some(Vec::new()) } else { None },
            hlay_aa_seqs: if genes.covers_mhc_class1() { Some(Vec::new()) } else { None },
            exon_sequences: HashMap::new(),
            common_alleles: Vec::new(),
            stop_loss_alleles: Vec::new(),
            stop_loss_indel: None,
            freqs,
            pon_indels: crate::indel::pon_indels(config.v38),
            excluded,
            gene_exon_boundaries,
            gene_enrichment,
            contexts,
            actual_alleles: Vec::new(),
            restricted_alleles: Vec::new(),
            transcripts,
        })
    }

    /// excludeAllele, with the config's actual and restricted allele names as given
    fn exclude_allele(&self, a: AId, actual: &[String], restricted: &[String]) -> bool {
        if self.alleles.gene(a) == HlaGene::H {
            return true;
        }
        let f = self.alleles.four_digit(a);
        let name = self.alleles.name(f);
        if self.excluded.contains(&name) {
            return true;
        }
        if !restricted.is_empty() {
            if actual.iter().any(|x| x == name) {
                return false;
            }
            return !restricted.iter().any(|x| x == name);
        }
        false
    }

    /// loadSequenceFile(lines, out, isProtein)
    fn load_sequences(&mut self, lines: &[String], protein: bool, actual: &[String], restricted: &[String], template: &mut Option<SeqLoci>) -> Vec<SeqLoci> {
        self.load_sequence_lines(&lines[1.min(lines.len())..], protein, actual, restricted, template)
    }

    /// loadSequenceFile(lines, sequences, isProteinFile): lines without the header
    pub fn load_sequence_lines(&mut self, lines: &[String], protein: bool, actual: &[String], restricted: &[String], template: &mut Option<SeqLoci>) -> Vec<SeqLoci> {
        let mut out = Vec::new();
        let all_genes = self.genes.all_genes();
        for line in lines.iter() {
            let items = split(line, ',');
            if items.len() != 2 {
                break;
            }
            let allele_str = items[0];
            // the gene first: Java creates the allele for every line, which only fills its cache
            let Some(p) = crate::hla::parse_allele(allele_str) else { continue };
            if !all_genes.contains(&p.gene) {
                continue;
            }
            let a = if protein { self.alleles.request_four_digit(allele_str).unwrap() } else { self.alleles.request(allele_str).unwrap() };
            let is_template = protein && template.is_none() && self.alleles.name(a) == DEFLATE_TEMPLATE;
            let exclude = self.exclude_allele(a, actual, restricted);
            if !is_template && exclude {
                continue;
            }
            let toks = seq::parse_reference_sequence(items[1]);
            let mut s = SeqLoci::new(&mut self.alleles, a, toks);
            if p.gene == HlaGene::Drb1 {
                let nuc_count = self.gene_cache.nuc_length(HlaGene::Drb1);
                let trim = if protein { nuc_count / 3 } else { nuc_count } as usize;
                let trimmed = s.seq[..trim.min(s.len())].to_vec();
                s = SeqLoci::new(&mut self.alleles, a, trimmed);
            }
            if is_template {
                *template = Some(s.clone());
            }
            if !exclude {
                out.push(s);
            }
        }
        out
    }

    /// ReferenceData.load
    pub fn load(&mut self, files: &SequenceFiles, config: &LilacConfig, actual: &[String], restricted: &[String]) -> Result<()> {
        let mut template: Option<SeqLoci> = None;
        crate::info!("reading nucleotide file: {}hla_ref_nucleotide_sequences.csv", config.resource_dir);
        self.nuc_seqs = self.load_sequences(&files.nucleotides, false, actual, restricted, &mut template);
        crate::info!("loaded {} sequences from file {}hla_ref_nucleotide_sequences.csv", self.nuc_seqs.len(), config.resource_dir);
        crate::info!("reading protein file: {}hla_ref_aminoacid_sequences.csv", config.resource_dir);
        self.aa_seqs = self.load_sequences(&files.amino_acids, true, actual, restricted, &mut template);
        crate::info!("loaded {} sequences from file {}hla_ref_aminoacid_sequences.csv", self.aa_seqs.len(), config.resource_dir);

        // populateAminoAcidSequenceLookup
        for s in &self.aa_seqs {
            let a = s.allele;
            if self.alleles.four_digit(a) != a {
                bail!("allele({}) is not four-digit", self.alleles.name(a));
            }
            if self.alleles.gene(a).is_pseudo() {
                continue;
            }
            if !self.exon_sequences.contains_key(&a) {
                let ex = seq::exon_sequences(self.gene_cache.aa_boundaries(self.alleles.gene(a)), s);
                self.exon_sequences.insert(a, ex);
            }
        }
        // drop frequencies of alleles without sequences
        let names: Vec<String> = self.freqs.java_order().iter().map(|x| x.0.to_string()).collect();
        for name in names {
            let known = self.alleles.find_four_digit(&name).is_some_and(|a| self.exon_sequences.contains_key(&a));
            if known {
                continue;
            }
            crate::warn!("allele({name}) with cohort frequency has no loaded sequences, dropping from allele frequencies");
            self.freqs.remove(&name);
        }
        // rebuildProteinAlleles of the actual and restricted alleles
        self.actual_alleles = actual.iter().filter_map(|x| four_digit_name(x).and_then(|n| self.alleles.request_four_digit(&n))).collect();
        self.restricted_alleles = restricted.iter().filter_map(|x| four_digit_name(x).and_then(|n| self.alleles.request_four_digit(&n))).collect();

        // loadCommonAlleles
        let common: Vec<String> = self
            .freqs
            .java_order()
            .into_iter()
            .filter(|x| x.1 >= crate::constants::COMMON_ALLELES_FREQ_CUTOFF)
            .map(|x| x.0.to_string())
            .collect();
        for n in common {
            let a = self.alleles.request_four_digit(&n).unwrap();
            self.common_alleles.push(a);
        }
        if !self.common_alleles.is_empty() {
            crate::info!("loaded {} common alleles", self.common_alleles.len());
        }
        // loadStopLossRecoveryAllele
        if self.genes.covers_mhc_class1() {
            let chr = if config.v38 { "chr6" } else { "6" };
            let indel = if config.v38 { Indel::of(chr, 31269338, "CN", "C") } else { Indel::of(chr, 31237115, "CN", "C") };
            let a = self.alleles.request_four_digit(STOP_LOSS_ON_C_ALLELE).unwrap();
            self.stop_loss_indel = Some(indel.clone());
            self.stop_loss_alleles.push((indel, a));
        }
        if let Some(hy) = self.hlay_nuc_seqs.as_mut() {
            let (y, rest): (Vec<SeqLoci>, Vec<SeqLoci>) = std::mem::take(&mut self.nuc_seqs).into_iter().partition(|s| self.alleles.gene(s.allele) == HlaGene::Y);
            *hy = y;
            self.nuc_seqs = rest;
        }
        // markExonBoundaryWildcards
        for s in self.aa_seqs.iter_mut() {
            let g = self.alleles.gene(s.allele);
            if g.is_pseudo() {
                continue;
            }
            s.set_exon_boundary_wildcards(self.gene_cache.aa_boundaries.get(&g).map(|v| v.as_slice()).unwrap_or(&[]));
        }
        for s in self.nuc_seqs.iter_mut() {
            let g = self.alleles.gene(s.allele);
            if g.is_pseudo() {
                continue;
            }
            s.set_exon_boundary_wildcards(self.gene_cache.nuc_boundaries.get(&g).map(|v| v.as_slice()).unwrap_or(&[]));
        }
        // the sequences with inserts and deletes, taken while the HLA-Y protein sequences are still in the
        // list (Java fills these lists before it moves HLA-Y out, so they keep its sequences)
        for s in self.aa_seqs.iter() {
            if s.has_inserts {
                self.aa_with_inserts.push(s.clone());
            } else if s.has_deletes || self.stop_loss_alleles.iter().any(|x| x.1 == s.allele) {
                self.aa_with_deletes.push(s.clone());
            }
        }
        // buildHlaYAminoAcidSequences
        if self.hlay_aa_seqs.is_some() {
            let (y, rest): (Vec<SeqLoci>, Vec<SeqLoci>) = std::mem::take(&mut self.aa_seqs).into_iter().partition(|s| self.alleles.gene(s.allele) == HlaGene::Y);
            self.aa_seqs = rest;
            if !y.is_empty() {
                self.hlay_aa_seqs = Some(y);
            } else {
                let t = template.as_ref().context("no amino acid template sequence")?;
                let mut built = Vec::new();
                for n in self.hlay_nuc_seqs.as_ref().unwrap().clone() {
                    let toks = seq::amino_acids_from_nucleotides(&n, t);
                    let a = self.alleles.four_digit(n.allele);
                    built.push(SeqLoci::new(&mut self.alleles, a, toks));
                }
                self.hlay_aa_seqs = Some(built);
            }
        }
        Ok(())
    }
}

/// the four-digit form of an allele name (HlaAllele.fromString(x).asFourDigit().toString())
pub fn four_digit_name(s: &str) -> Option<String> {
    let p = crate::hla::parse_allele(s)?;
    Some(crate::hla::allele_name(p.gene, &p.group, &p.protein, "", ""))
}
