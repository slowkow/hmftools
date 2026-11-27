//! The `expression` package: fragment category counts, the expected counts per gene collection, the
//! expectation-maximisation fit of transcript abundances, TPMs, the cohort percentiles and the gene and
//! transcript results.

use crate::common::*;
use crate::compat;
use crate::config::{FragmentSize, IsofoxConfig, HIGH_EXPRESSION_FOLD_CHANGE_HIGH, HIGH_EXPRESSION_FOLD_CHANGE_LOW, HIGH_EXPRESSION_PERCENTILE_HIGH, MAX_GENE_PERC_CONTRIBUTION};
use crate::gene::{find_exon_region, GeneCollection};
use anyhow::{anyhow, Result};
use hmf_common::driver_panel::DriverGene;
use hmf_common::ensembl::{GeneData, TranscriptData};
use hmf_common::fmt::{java_fixed, java_sci};
use hmf_common::java_hash::{JavaHashMap, JavaHashSet};
use hmf_common::java_seq_hash::SeqHashMap;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU32, Ordering};

/// CategoryCountsData: the fragments that could support a set of transcripts and/or unspliced genes
#[derive(Clone, Debug)]
pub struct CategoryCountsData {
    pub transcripts: Vec<i32>,
    pub unspliced_genes: Vec<String>,
    pub fragment_count: f64,
    pub multi_mapped: i32,
    pub counts_by_gc: Option<Vec<f64>>,
    pub counts_by_length: Option<Vec<f64>>,
    pub combined_key: String,
}

impl CategoryCountsData {
    pub fn new(transcripts: Vec<i32>, unspliced_genes: Vec<String>) -> Self {
        let combined_key = form_transcript_ids(&transcripts, &unspliced_genes);
        CategoryCountsData { transcripts, unspliced_genes, fragment_count: 0.0, multi_mapped: 0, counts_by_gc: None, counts_by_length: None, combined_key }
    }
    /// from an expected counts file category
    pub fn from_key(key: &str, frag_lengths: usize) -> Self {
        let mut c = CategoryCountsData {
            transcripts: Vec::new(),
            unspliced_genes: Vec::new(),
            fragment_count: 0.0,
            multi_mapped: 0,
            counts_by_gc: None,
            counts_by_length: Some(vec![0.0; frag_lengths]),
            combined_key: key.to_string(),
        };
        for item in hmf_common::java_text::split(key, '-') {
            if item.contains("ENSG") {
                c.unspliced_genes.push(item.to_string());
            } else {
                // Integer.parseInt: Java does not catch a bad value
                c.transcripts.push(item.parse().unwrap_or(0));
            }
        }
        c
    }
    pub fn matches(&self, transcripts: &[i32], genes: &[&str]) -> bool {
        if self.transcripts.len() != transcripts.len() || self.unspliced_genes.len() != genes.len() {
            return false;
        }
        transcripts.iter().all(|t| self.transcripts.contains(t)) && genes.iter().all(|g| self.unspliced_genes.iter().any(|x| x == g))
    }
    pub fn add_counts(&mut self, count: f64, multi_mapped: bool) {
        self.fragment_count += count;
        if multi_mapped {
            self.multi_mapped += 1;
        }
    }
    pub fn add_gc_ratio_counts(&mut self, count: f64, multi_mapped: bool, idx: &[i32; 2], counts: &[f64; 2]) {
        self.fragment_count += count;
        if multi_mapped {
            self.multi_mapped += 1;
        }
        let gc = self.counts_by_gc.as_mut().expect("GC ratio counts");
        for i in 0..2 {
            if idx[i] >= 0 {
                gc[idx[i] as usize] += counts[i];
            }
        }
    }
    pub fn apply_gc_adjustments(&mut self, adj: &[f64]) {
        self.fragment_count = 0.0;
        let gc = self.counts_by_gc.as_mut().expect("GC ratio counts");
        for i in 0..gc.len() {
            gc[i] *= adj[i];
            self.fragment_count += gc[i];
        }
    }
    pub fn apply_frequencies(&mut self, rates: &[f64]) {
        let c = self.counts_by_length.as_mut().unwrap();
        if rates.len() != c.len() {
            return;
        }
        self.fragment_count = 0.0;
        for i in 0..c.len() {
            c[i] *= rates[i];
            self.fragment_count += c[i];
        }
    }
}

fn form_transcript_ids(transcripts: &[i32], genes: &[String]) -> String {
    let mut ids: Vec<i32> = Vec::with_capacity(transcripts.len());
    for t in transcripts {
        let mut index = 0;
        while index < ids.len() {
            if *t < ids[index] {
                break;
            }
            index += 1;
        }
        ids.insert(index, *t);
    }
    let mut items: Vec<String> = ids.iter().map(|x| x.to_string()).collect();
    items.extend(genes.iter().cloned());
    items.join("-")
}

pub fn has_gene_identifier(s: &str) -> bool {
    s.starts_with("ENSG")
}

// ---- ExpectedCountsCache ----

struct GeneSet {
    categories: Vec<Box<str>>,
    /// rows x fragment lengths
    counts: Vec<i32>,
    gene_ids: Vec<Box<str>>,
    /// how often the fragment length frequencies were applied to this set (Java changes the cached counts)
    uses: AtomicU32,
}

pub struct ExpectedCountsCache {
    index: SeqHashMap<String, usize>,
    sets: Vec<GeneSet>,
    frag_lengths: usize,
    pub valid: bool,
}

impl ExpectedCountsCache {
    pub fn new(config: &mut IsofoxConfig) -> Result<ExpectedCountsCache> {
        let mut cache = ExpectedCountsCache { index: SeqHashMap::new(), sets: Vec::new(), frag_lengths: config.fragment_size_data.len(), valid: true };
        if let Some(f) = config.exp_counts_file.clone() {
            if std::path::Path::new(&f).exists() {
                cache.valid = cache.load(&f, config)?;
            }
        }
        Ok(cache)
    }

    fn load(&mut self, filename: &str, config: &mut IsofoxConfig) -> Result<bool> {
        let delim = if filename.ends_with(".tsv") || filename.ends_with(".tsv.gz") { '\t' } else { ',' };
        let mut reader = hmf_common::gz::LineReader::open(filename)?;
        let Some(header) = reader.next_line()? else {
            hmf_common::error!("empty calculated expected counts file({})", filename);
            return Ok(false);
        };
        let header = std::str::from_utf8(header)?.to_string();
        let items: Vec<&str> = hmf_common::java_text::split_all(&header, delim);
        let file_count = items.len() as i32 - 2;
        if config.fragment_size_data.is_empty() {
            let mut lens = Vec::new();
            for it in items.iter().skip(3) {
                let l: i32 = it.replace("Length_", "").parse().map_err(|_| anyhow!("invalid header {it}"))?;
                config.fragment_size_data.push(FragmentSize { length: l, frequency: 1 });
                lens.push(l.to_string());
            }
            hmf_common::info!("fragment lengths({}) set from expected transcript counts file", lens.join(", "));
        } else if config.fragment_size_data.len() as i32 != file_count {
            hmf_common::error!("expected counts file has {} fragment lengths vs configuredCount({})", file_count, config.fragment_size_data.len());
            return Ok(false);
        }
        let n = config.fragment_size_data.len();
        self.frag_lengths = n;
        let fields = hmf_common::java_text::fields_index_map(&header, delim);
        let gs_idx = *fields.get("GeneSetId").ok_or_else(|| anyhow!("no GeneSetId column"))?;
        let cat_idx = *fields.get("Category").ok_or_else(|| anyhow!("no Category column"))?;
        let mut current = String::new();
        let mut cur_set: Option<usize> = None;
        while let Some(line) = reader.next_line()? {
            let line = std::str::from_utf8(line)?;
            let values: Vec<&str> = hmf_common::java_text::split_all(line, delim);
            let gs = values[gs_idx];
            let cat = values[cat_idx];
            if gs != current || cur_set.is_none() {
                current = gs.to_string();
                self.sets.push(GeneSet { categories: Vec::new(), counts: Vec::new(), gene_ids: Vec::new(), uses: AtomicU32::new(0) });
                let si = self.sets.len() - 1;
                self.index.insert(current.clone(), si);
                cur_set = Some(si);
            }
            let set = &mut self.sets[cur_set.unwrap()];
            set.categories.push(cat.into());
            for item in hmf_common::java_text::split(cat, '-') {
                if item.contains("ENSG") && !set.gene_ids.iter().any(|g| &**g == item) {
                    set.gene_ids.push(item.into());
                }
            }
            for i in 0..n {
                let v = values.get(cat_idx + i + 1).ok_or_else(|| anyhow!("short expected counts line"))?;
                set.counts.push(v.parse().map_err(|_| anyhow!("invalid expected count: {v}"))?);
            }
        }
        hmf_common::info!("loaded {} gene expected counts from file({})", self.index.len(), filename);
        Ok(true)
    }

    fn set_matches(&self, si: usize, gene_ids: &[&str]) -> bool {
        let s = &self.sets[si];
        gene_ids.iter().all(|g| s.gene_ids.iter().any(|x| &**x == *g))
    }

    /// getGeneExpectedRatesData with the fragment length frequencies applied (applyFragmentLengthDistributionToExpectedCounts)
    pub fn gene_expected_rates_data(&self, chr_id: &str, gene_ids: &[&str], rates: Option<&[f64]>) -> Option<Vec<CategoryCountsData>> {
        let mut si = self.index.get(&chr_id.to_string()).copied();
        if si.is_none() || !self.set_matches(si.unwrap(), gene_ids) {
            si = None;
            for (_, i) in self.index.iter() {
                if self.set_matches(*i, gene_ids) {
                    si = Some(*i);
                    break;
                }
            }
        }
        let si = si?;
        let set = &self.sets[si];
        let n = self.frag_lengths;
        // the frequencies are applied to the cached counts in Java, so a set used again has them applied
        // again (IF1); corrected, each use starts from the file's counts
        let applications = match rates {
            None => 0,
            Some(_) => {
                let prior = set.uses.fetch_add(1, Ordering::SeqCst);
                if compat::java("IF1") { prior + 1 } else { 1 }
            }
        };
        let mut out = Vec::with_capacity(set.categories.len());
        for (ci, cat) in set.categories.iter().enumerate() {
            let mut c = CategoryCountsData::from_key(cat, n);
            for i in 0..n {
                let v = set.counts[ci * n + i];
                c.fragment_count += v as f64;
                c.counts_by_length.as_mut().unwrap()[i] += v as f64;
            }
            if let Some(r) = rates {
                for _ in 0..applications {
                    c.apply_frequencies(r);
                }
            }
            out.push(c);
        }
        Some(out)
    }
}

/// ExpectedRatesData: categories (rows) and transcripts or unspliced genes (columns) of a gene collection
#[derive(Clone, Debug, Default)]
pub struct ExpectedRatesData {
    pub categories: Vec<String>,
    pub transcript_ids: Vec<String>,
    /// column-major: definitions[col * rows + row]
    pub definitions: Vec<f64>,
    pub has_definitions: bool,
}

impl ExpectedRatesData {
    pub fn valid_data(&self) -> bool {
        !self.categories.is_empty() && self.has_definitions && !self.transcript_ids.is_empty()
    }
    pub fn category_index(&self, c: &str) -> i32 {
        self.categories.iter().position(|x| x == c).map(|x| x as i32).unwrap_or(-1)
    }
    fn add_category(&mut self, c: &str) {
        if !self.categories.iter().any(|x| x == c) {
            self.categories.push(c.to_string());
        }
    }
}

/// ExpectedRatesCommon.formTranscriptDefinitions
pub fn form_transcript_definitions(data: &[CategoryCountsData]) -> ExpectedRatesData {
    let mut e = ExpectedRatesData::default();
    for c in data {
        if c.fragment_count > 0.0 || c.transcripts.is_empty() {
            e.add_category(&c.combined_key);
        }
    }
    // createTransComboDataMap
    let mut map: SeqHashMap<String, Vec<usize>> = SeqHashMap::new();
    for (ci, c) in data.iter().enumerate() {
        let mut names: JavaHashSet<String> = JavaHashSet::new();
        for g in &c.unspliced_genes {
            names.insert(g.clone());
        }
        for t in &c.transcripts {
            names.insert(t.to_string());
        }
        for n in names.iter() {
            match map.get_mut(n) {
                Some(l) => l.push(ci),
                None => {
                    map.insert(n.clone(), vec![ci]);
                }
            }
        }
    }
    let rows = e.categories.len();
    for (k, _) in map.iter() {
        e.transcript_ids.push(k.clone());
    }
    if rows == 0 || e.transcript_ids.is_empty() {
        return e;
    }
    e.definitions = vec![0.0; rows * e.transcript_ids.len()];
    e.has_definitions = true;
    for ti in 0..e.transcript_ids.len() {
        let mut counts = vec![0.0; rows];
        let list = map.get(&e.transcript_ids[ti]).unwrap();
        for ci in list {
            let c = &data[*ci];
            if c.fragment_count > 0.0 {
                let cat = e.category_index(&c.combined_key);
                if cat < 0 {
                    hmf_common::error!("invalid category index from transKey({})", c.combined_key);
                    return e;
                }
                counts[cat as usize] = c.fragment_count;
            }
        }
        crate::adjusts::convert_to_percentages(&mut counts);
        e.definitions[ti * rows..(ti + 1) * rows].copy_from_slice(&counts);
    }
    e
}

fn sum_vector(v: &[f64]) -> f64 {
    let mut t = 0.0;
    for x in v {
        t += x;
    }
    t
}

/// SigUtils.calculateFittedCounts
pub fn calculate_fitted_counts(defs: &[f64], rows: usize, cols: usize, allocs: &[f64]) -> Vec<f64> {
    let mut fitted = vec![0.0; rows];
    for t in 0..cols {
        let a = allocs[t];
        let col = &defs[t * rows..(t + 1) * rows];
        for c in 0..rows {
            fitted[c] += a * col[c];
        }
    }
    fitted
}

/// SigUtils.calcResiduals: (total, percent)
pub fn calc_residuals(counts: &[f64], fitted: &[f64], total_counts: f64) -> (f64, f64) {
    let mut total = 0.0;
    for c in 0..counts.len() {
        let diff = fitted[c] - counts[c];
        total += diff.abs();
    }
    (total, total / total_counts)
}

/// ExpectationMaxFit.performFit
pub fn perform_fit(counts: &[f64], defs: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    const MAX_ITERATIONS: usize = 20;
    const RESIDUALS_EXIT_PERC: f64 = 0.001;
    let total = sum_vector(counts);
    let init = 1.0 / cols as f64;
    let mut allocs = vec![init; cols];
    let mut new_allocs = vec![0.0; cols];
    let mut iteration = 0;
    while iteration < MAX_ITERATIONS {
        let mut factors = vec![0.0; rows];
        for t in 0..cols {
            new_allocs[t] = 0.0;
            let col = &defs[t * rows..(t + 1) * rows];
            for c in 0..rows {
                factors[c] += allocs[t] * col[c];
            }
        }
        for t in 0..cols {
            let col = &defs[t * rows..(t + 1) * rows];
            let ta = allocs[t];
            for c in 0..rows {
                if factors[c] == 0.0 {
                    continue;
                }
                new_allocs[t] += counts[c] * col[c] * ta / factors[c];
            }
        }
        let fitted = calculate_fitted_counts(defs, rows, cols, &new_allocs);
        let (_, perc) = calc_residuals(counts, &fitted, total);
        if perc < RESIDUALS_EXIT_PERC {
            break;
        }
        for t in 0..cols {
            allocs[t] = new_allocs[t] / total;
        }
        iteration += 1;
    }
    new_allocs
}

/// TranscriptResult
#[derive(Clone, Debug)]
pub struct TranscriptResult {
    pub trans: &'static TranscriptData,
    pub exonic_bases: i32,
    pub exonic_bases_covered: i32,
    pub sj_supported: i32,
    pub unique_sj_supported: i32,
    pub unique_sj_fragments: i32,
    pub unique_non_sj_fragments: i32,
    pub discordant_fragments: i32,
    pub effective_length: f64,
    pub fit_allocation: f64,
    pub raw_fit_allocation: f64,
    pub raw_tpm: f64,
    pub adjusted_tpm: f64,
    pub multi_mapped: f64,
}

pub fn calc_effective_length(trans_length: i32, data: &[FragmentSize]) -> f64 {
    if data.is_empty() {
        return trans_length as f64;
    }
    let mut freq_total: i64 = 0;
    let mut bases_total: i64 = 0;
    for fl in data {
        if fl.length >= trans_length {
            continue;
        }
        let possible = (trans_length - fl.length) as i64;
        freq_total += fl.frequency as i64;
        bases_total += possible * fl.frequency as i64;
    }
    if freq_total == 0 {
        return 0.0;
    }
    bases_total as f64 / freq_total as f64
}

fn is_splice_junction_unique(trans_name: &str, trans_list: &[&TranscriptData], exon_end: i32, exon_start: i32) -> bool {
    for t in trans_list {
        if t.trans_name == trans_name {
            continue;
        }
        for i in 1..t.exons.len() {
            if t.exons[i - 1].end == exon_end && t.exons[i].start == exon_start {
                return false;
            }
        }
    }
    true
}

impl TranscriptResult {
    pub fn new(gc: &GeneCollection, gi: usize, trans: &'static TranscriptData, frag_lengths: &[FragmentSize]) -> Self {
        let gene = &gc.genes[gi];
        let mut sj_supported = 0;
        let mut exonic_bases = 0;
        let mut covered = 0;
        let mut unique_sj = 0;
        for (i, exon) in trans.exons.iter().enumerate() {
            let Some(ri) = find_exon_region(&gc.exon_regions, Some(&gene.exon_regions), exon.start, exon.end) else { continue };
            let r = &gc.exon_regions[ri];
            covered += r.base_coverage(1);
            exonic_bases += exon.end - exon.start + 1;
            if i > 0 {
                let sj = r.transcript_junction_match_count(trans.trans_id, SE_START);
                let prev = &trans.exons[i - 1];
                let unique = is_splice_junction_unique(&trans.trans_name, &gene.transcripts, prev.end, exon.start);
                if sj[TRANS_COUNT] > 0 {
                    sj_supported += 1;
                    if unique {
                        unique_sj += 1;
                    }
                }
            }
        }
        let frags = gc.transcript_read_count(trans.trans_id);
        TranscriptResult {
            trans,
            exonic_bases,
            exonic_bases_covered: covered,
            sj_supported,
            unique_sj_supported: unique_sj,
            unique_sj_fragments: frags[FragmentMatchType::Spliced as usize][UNIQUE_TRANS_COUNT],
            unique_non_sj_fragments: frags[FragmentMatchType::Long as usize][UNIQUE_TRANS_COUNT],
            discordant_fragments: frags[FragmentMatchType::Discordant as usize][TRANS_COUNT],
            effective_length: calc_effective_length(exonic_bases, frag_lengths),
            fit_allocation: 0.0,
            raw_fit_allocation: 0.0,
            raw_tpm: 0.0,
            adjusted_tpm: 0.0,
            multi_mapped: 0.0,
        }
    }
    pub fn fragments_per_kb(&self) -> f64 {
        if self.effective_length > 0.0 { self.fit_allocation / (self.effective_length / 1000.0) } else { 0.0 }
    }
    pub fn header() -> String {
        [
            "GeneId",
            "GeneName",
            "TransId",
            "TransName",
            "ExonCount",
            "TranscriptLength",
            "EffectiveLength",
            "FittedFragments",
            "RawFittedFragments",
            "AdjTPM",
            "RawTPM",
            "TranscriptBasesCovered",
            "SJSupported",
            "UniqueSJSupported",
            "UniqueSJFragments",
            "UniqueNonSJFragments",
            "DiscordantFragments",
            if crate::version::v201() { "LowMapQualFrags" } else { "MultiMappedFragments" },
        ]
        .join("\t")
    }
    pub fn to_line(&self, gene: &GeneData) -> String {
        [
            gene.gene_id.clone(),
            gene.gene_name.clone(),
            self.trans.trans_id.to_string(),
            self.trans.trans_name.clone(),
            self.trans.exons.len().to_string(),
            self.exonic_bases.to_string(),
            java_fixed(self.effective_length, 0),
            java_fixed(self.fit_allocation, 1),
            java_fixed(self.raw_fit_allocation, 1),
            java_sci(self.adjusted_tpm, 3),
            java_sci(self.raw_tpm, 3),
            self.exonic_bases_covered.to_string(),
            self.sj_supported.to_string(),
            self.unique_sj_supported.to_string(),
            self.unique_sj_fragments.to_string(),
            self.unique_non_sj_fragments.to_string(),
            self.discordant_fragments.to_string(),
            java_fixed(self.multi_mapped, 1),
        ]
        .join("\t")
    }
}

/// GeneResult
#[derive(Clone, Debug)]
pub struct GeneResult {
    pub gene: &'static GeneData,
    pub collection_id: String,
    pub intronic_length: i32,
    pub trans_count: i32,
    pub spliced_alloc: f64,
    pub unspliced_alloc: f64,
    pub raw_tpm: f64,
    pub adjusted_tpm: f64,
    pub fit_residuals: f64,
    pub multi_mapped: f64,
    pub median_tpm_cohort: f64,
    pub percentile_tpm_cohort: f64,
    pub median_tpm_cancer: f64,
    pub percentile_tpm_cancer: f64,
    pub reported: bool,
}

impl GeneResult {
    pub fn new(gc: &GeneCollection, gi: usize) -> Self {
        let g = &gc.genes[gi];
        let exonic = gc.exonic_region_length(gi);
        // IF2: the gene length and its exonic regions are measured as end - start
        let intronic = if compat::java("IF2") {
            (g.gene.length() as i64 - exonic) as i32
        } else {
            let merged = crate::gene::generate_common_exonic_regions(&gc.exon_regions, &g.exon_regions);
            let exonic_bases: i64 = merged.iter().map(|x| (x[1] - x[0] + 1) as i64).sum();
            ((g.gene.length() + 1) as i64 - exonic_bases) as i32
        };
        GeneResult {
            gene: g.gene,
            collection_id: gc.chr_id(),
            intronic_length: intronic,
            trans_count: g.transcripts.len() as i32,
            spliced_alloc: 0.0,
            unspliced_alloc: 0.0,
            raw_tpm: 0.0,
            adjusted_tpm: 0.0,
            fit_residuals: 0.0,
            multi_mapped: 0.0,
            median_tpm_cohort: 0.0,
            percentile_tpm_cohort: 0.0,
            median_tpm_cancer: 0.0,
            percentile_tpm_cancer: 0.0,
            reported: false,
        }
    }
    pub fn header() -> String {
        [
            "GeneId",
            "GeneName",
            "Chromosome",
            "GeneLength",
            "IntronicLength",
            "TranscriptCount",
            "GeneSetId",
            "SplicedFragments",
            "UnsplicedFragments",
            "AdjTPM",
            "RawTPM",
            "FitResiduals",
            if crate::version::v201() { "LowMapQualFrags" } else { "MultiMappedFragments" },
            "MedianTpmCancer",
            "PercentileCancer",
            "MedianTpmCohort",
            "PercentileCohort",
            "ReportedStatus",
        ]
        .join("\t")
    }
    pub fn to_line(&self) -> String {
        [
            self.gene.gene_id.clone(),
            self.gene.gene_name.clone(),
            self.gene.chromosome.clone(),
            self.gene.length().to_string(),
            self.intronic_length.to_string(),
            self.trans_count.to_string(),
            self.collection_id.clone(),
            java_fixed(self.spliced_alloc, 0),
            java_fixed(self.unspliced_alloc, 0),
            java_sci(self.adjusted_tpm, 3),
            java_sci(self.raw_tpm, 3),
            java_fixed(self.fit_residuals, 1),
            java_fixed(self.multi_mapped, 1),
            java_sci(self.median_tpm_cancer, 3),
            java_fixed(self.percentile_tpm_cancer, 3),
            java_sci(self.median_tpm_cohort, 3),
            java_fixed(self.percentile_tpm_cohort, 3),
            if self.reported { "REPORTED".to_string() } else { "NOT_REPORTED".to_string() },
        ]
        .join("\t")
    }
}

/// GeneCollectionSummary
pub struct GeneCollectionSummary {
    pub chr_id: String,
    pub gene_ids: Vec<&'static str>,
    pub gene_names: String,
    pub trans_category_counts: Vec<CategoryCountsData>,
    pub gene_results: Vec<GeneResult>,
    pub transcript_results: Vec<TranscriptResult>,
    total_multi_mapped: i32,
    pub fit_allocations: JavaHashMap<String, f64>,
    pub fit_residuals: f64,
    /// the expected rates of the first fit, kept for the GC adjustment re-fit
    pub expected_rates: Option<ExpectedRatesData>,
    /// transcript id -> name of the collection's transcripts
    pub trans_id_names: HashMap<i32, &'static str>,
}

impl GeneCollectionSummary {
    pub fn new(chr_id: String, gene_ids: Vec<&'static str>, gene_names: String, counts: Vec<CategoryCountsData>) -> Self {
        let total_multi_mapped = counts.iter().map(|c| c.multi_mapped).sum();
        GeneCollectionSummary {
            chr_id,
            gene_ids,
            gene_names,
            trans_category_counts: counts,
            gene_results: Vec::new(),
            transcript_results: Vec::new(),
            total_multi_mapped,
            fit_allocations: JavaHashMap::new(),
            fit_residuals: 0.0,
            expected_rates: None,
            trans_id_names: HashMap::new(),
        }
    }
    pub fn splice_genes_count(&self) -> i32 {
        self.gene_results.iter().filter(|g| g.spliced_alloc > 0.0).count() as i32
    }
    pub fn fit_allocation(&self, name: &str) -> f64 {
        self.fit_allocations.get(&name.to_string()).copied().unwrap_or(0.0)
    }
    pub fn set_fit_allocations(&mut self) {
        let mut gene_totals: HashMap<&str, f64> = HashMap::new();
        for i in 0..self.transcript_results.len() {
            let name = self.transcript_results[i].trans.trans_name.as_str();
            let a = self.fit_allocation(name);
            self.transcript_results[i].fit_allocation = a;
            let gid = self.transcript_results[i].trans.gene_id.as_str();
            match gene_totals.get_mut(gid) {
                Some(t) => *t += a,
                None => {
                    gene_totals.insert(gid, a);
                }
            }
        }
        for i in 0..self.gene_results.len() {
            let gid = self.gene_results[i].gene.gene_id.as_str();
            let spliced = gene_totals.get(gid).copied().unwrap_or(0.0);
            let unspliced = self.fit_allocation(gid);
            self.gene_results[i].spliced_alloc = spliced;
            self.gene_results[i].unspliced_alloc = unspliced;
        }
    }
    pub fn allocate_residuals_to_genes(&mut self) {
        if self.gene_results.len() == 1 {
            self.gene_results[0].fit_residuals = self.fit_residuals;
        } else {
            let total: i64 = self.gene_results.iter().map(|g| g.gene.length() as i64).sum();
            for g in &mut self.gene_results {
                g.fit_residuals = g.gene.length() as f64 / total as f64 * self.fit_residuals;
            }
        }
    }
    pub fn apply_gc_adjustments(&mut self, adj: &[f64]) {
        let mut orig = 0.0;
        let mut new_total = 0.0;
        for c in &mut self.trans_category_counts {
            orig += c.fragment_count;
            c.apply_gc_adjustments(adj);
            new_total += c.fragment_count;
        }
        let factor = orig / new_total;
        for c in &mut self.trans_category_counts {
            c.fragment_count *= factor;
        }
    }
    pub fn assign_multi_mapped_fragments(&mut self) {
        if self.total_multi_mapped == 0 {
            return;
        }
        let mut total_trans = 0.0;
        for t in &self.transcript_results {
            total_trans += t.fit_allocation;
        }
        let mut total_unspliced = 0.0;
        for g in &self.gene_results {
            total_unspliced += g.unspliced_alloc;
        }
        let total = total_trans + total_unspliced;
        if total == 0.0 {
            return;
        }
        let mm = self.total_multi_mapped as f64;
        let spliced_mm = mm * total_trans / total;
        let unspliced_mm = mm * total_unspliced / total;
        for t in &mut self.transcript_results {
            t.multi_mapped = if total_trans > 0.0 { t.fit_allocation / total_trans * spliced_mm } else { 0.0 };
        }
        for g in &mut self.gene_results {
            let s = if total_trans > 0.0 { g.spliced_alloc / total_trans * spliced_mm } else { 0.0 };
            let u = if total_unspliced > 0.0 { g.unspliced_alloc / total_unspliced * unspliced_mm } else { 0.0 };
            g.multi_mapped = s + u;
        }
    }
}

/// TranscriptExpression: the fit of one gene collection
pub struct TranscriptExpression<'a> {
    pub config: &'a IsofoxConfig,
    pub cache: &'a ExpectedCountsCache,
    pub frequency_rates: Vec<f64>,
}

impl<'a> TranscriptExpression<'a> {
    pub fn new(config: &'a IsofoxConfig, cache: &'a ExpectedCountsCache) -> Self {
        let mut total = 0.0;
        for f in &config.fragment_size_data {
            total += f.frequency as f64;
        }
        let mut rates = Vec::new();
        if total > 0.0 {
            for f in &config.fragment_size_data {
                rates.push(f.frequency as f64 / total);
            }
        }
        TranscriptExpression { config, cache, frequency_rates: rates }
    }

    /// runTranscriptEstimation; `check_cached` for the re-fit after the GC adjustment
    pub fn run_transcript_estimation(&self, summary: &mut GeneCollectionSummary, check_cached: bool) {
        let exp: ExpectedRatesData = if check_cached {
            match &summary.expected_rates {
                Some(e) => e.clone(),
                None => {
                    hmf_common::debug!("gene({}) invalid expected rates or actuals data", summary.gene_names);
                    return;
                }
            }
        } else {
            let rates = if self.config.apply_fragment_length_adjust { Some(self.frequency_rates.as_slice()) } else { None };
            match self.cache.gene_expected_rates_data(&summary.chr_id, &summary.gene_ids, rates) {
                None => {
                    hmf_common::warn!("genes({}: {}) expected counts data not loaded", summary.chr_id, summary.gene_ids.join(";"));
                    return;
                }
                Some(data) => form_transcript_definitions(&data),
            }
        };
        if !exp.valid_data() {
            hmf_common::debug!("gene({}) invalid expected rates or actuals data", summary.gene_names);
            return;
        }
        if !check_cached && self.config.apply_gc_bias_adjust() {
            summary.expected_rates = Some(exp.clone());
        }
        let rows = exp.categories.len();
        let cols = exp.transcript_ids.len();
        let mut counts = vec![0.0; rows];
        for tc in &summary.trans_category_counts {
            if tc.fragment_count > 0.0 {
                let ci = exp.category_index(&tc.combined_key);
                if ci >= 0 {
                    counts[ci as usize] = tc.fragment_count;
                }
            }
        }
        let total = sum_vector(&counts);
        if total == 0.0 {
            return;
        }
        let allocs = perform_fit(&counts, &exp.definitions, rows, cols);
        let fitted = calculate_fitted_counts(&exp.definitions, rows, cols, &allocs);
        let (res_total, _) = calc_residuals(&counts, &fitted, total);
        summary.fit_residuals = res_total;
        for (ti, tid) in exp.transcript_ids.iter().enumerate() {
            let name: Option<String> = if has_gene_identifier(tid) {
                Some(tid.clone())
            } else {
                tid.parse::<i32>().ok().and_then(|id| summary.trans_id_names.get(&id)).map(|s| s.to_string())
            };
            if let Some(n) = name {
                summary.fit_allocations.insert(n, allocs[ti]);
            }
        }
    }
}

/// TranscriptExpression.calcTpmFactors: raw and adjusted (capping any gene at 1% of the total)
pub fn calc_tpm_factors(summaries: &[&GeneCollectionSummary], enriched_gene_ids: &[String]) -> [f64; 2] {
    let mut set: Vec<f64> = Vec::with_capacity(200_000);
    let mut max_fpkb = 0.0f64;
    let mut total = 0.0;
    let mut raw_total = 0.0;
    for s in summaries {
        for t in &s.transcript_results {
            let f = t.fragments_per_kb();
            raw_total += f;
            // ISOFOX 2.0.1 leaves the enriched genes out of the adjusted TPM total
            if enriched_gene_ids.iter().any(|x| *x == t.trans.gene_id) {
                continue;
            }
            if f < 0.5 {
                total += f;
            } else {
                set.push(f);
                total += f;
                max_fpkb = max_fpkb.max(f);
            }
        }
    }
    let raw = raw_total / 1_000_000.0;
    if max_fpkb > MAX_GENE_PERC_CONTRIBUTION * total {
        let mut iterations = 0;
        let mut max_total = total;
        let mut next_total = total * 0.5;
        let mut min_total = 0.0f64;
        while iterations < 10 {
            let max_allowed = next_total * MAX_GENE_PERC_CONTRIBUTION;
            let mut max_value = 0.0f64;
            total = 0.0;
            for f in &set {
                let v = f.min(max_allowed);
                max_value = max_value.max(v);
                total += v;
            }
            let perc = max_value / total;
            if (perc - MAX_GENE_PERC_CONTRIBUTION).abs() < 0.001 {
                break;
            }
            if perc > MAX_GENE_PERC_CONTRIBUTION {
                max_total = total;
                next_total = (total * 0.5).max(min_total);
            } else {
                min_total = total;
                next_total = (total + max_total) * 0.5;
            }
            iterations += 1;
        }
    }
    [raw, total / 1_000_000.0]
}

pub fn set_transcripts_per_million(summary: &mut GeneCollectionSummary, factors: [f64; 2]) {
    let mut gene_tpms: HashMap<&str, [f64; 2]> = HashMap::new();
    for t in &mut summary.transcript_results {
        let f = t.fragments_per_kb();
        let raw = if factors[0] > 0.0 { f / factors[0] } else { 0.0 };
        let adj = if factors[1] > 0.0 { f / factors[1] } else { 0.0 };
        t.raw_tpm = raw;
        t.adjusted_tpm = adj;
        match gene_tpms.get_mut(t.trans.gene_id.as_str()) {
            None => {
                gene_tpms.insert(t.trans.gene_id.as_str(), [raw, adj]);
            }
            Some(g) => {
                g[0] += raw;
                g[1] += adj;
            }
        }
    }
    for g in &mut summary.gene_results {
        // Java takes the gene's entry unchecked: a gene without transcripts cannot reach here
        let v = gene_tpms.get(g.gene.gene_id.as_str()).copied().unwrap_or([0.0, 0.0]);
        g.raw_tpm = v[0];
        g.adjusted_tpm = v[1];
    }
}

// ---- CohortGenePercentiles ----

pub const PAN_CANCER: &str = "ALL";
pub const PERCENTILE_COUNT: usize = 101;
const INVALID_VALUE: f64 = -1.0;

struct GeneCohort {
    medians: Vec<(Box<str>, f64)>,
    percentiles: Vec<(Box<str>, Vec<f64>)>,
}

/// CohortGenePercentiles, keeping the percentiles of the pan-cancer and the configured cancer type only
pub struct CohortGenePercentiles {
    genes: HashMap<String, GeneCohort>,
    max_medians: HashMap<String, f64>,
    cancer_types: std::collections::HashSet<String>,
}

impl CohortGenePercentiles {
    pub fn load(filename: &str, cancer_type: Option<&str>) -> Result<CohortGenePercentiles> {
        let mut c = CohortGenePercentiles { genes: HashMap::new(), max_medians: HashMap::new(), cancer_types: Default::default() };
        if !std::path::Path::new(filename).exists() {
            hmf_common::error!("invalid gene cohort distribution file({})", filename);
            return Ok(c);
        }
        let delim = if filename.ends_with(".tsv") || filename.ends_with(".tsv.gz") { '\t' } else { ',' };
        let mut reader = hmf_common::gz::LineReader::open(filename)?;
        let header = match reader.next_line()? {
            Some(h) => std::str::from_utf8(h)?.to_string(),
            None => return Ok(c),
        };
        let fields = hmf_common::java_text::fields_index_map(&header, delim);
        let gi = *fields.get("GeneId").ok_or_else(|| anyhow!("no GeneId"))?;
        let ci = *fields.get("CancerType").ok_or_else(|| anyhow!("no CancerType"))?;
        let mi = *fields.get("Median").ok_or_else(|| anyhow!("no Median"))?;
        let pi = *fields.get("Pct_0").ok_or_else(|| anyhow!("no Pct_0"))?;
        let mut current = String::new();
        while let Some(line) = reader.next_line()? {
            let line = std::str::from_utf8(line)?;
            let values: Vec<&str> = hmf_common::java_text::split_all(line, delim);
            let gene = values[gi];
            let cancer = values[ci];
            if !c.cancer_types.contains(cancer) {
                c.cancer_types.insert(cancer.to_string());
            }
            let median = hmf_common::java_text::parse_double(values[mi]).ok_or_else(|| anyhow!("invalid median {}", values[mi]))?;
            if current != gene {
                current = gene.to_string();
                c.genes.insert(current.clone(), GeneCohort { medians: Vec::new(), percentiles: Vec::new() });
            }
            let g = c.genes.get_mut(&current).unwrap();
            match g.medians.iter_mut().find(|(k, _)| &**k == cancer) {
                Some(e) => e.1 = median,
                None => g.medians.push((cancer.into(), median)),
            }
            if cancer == PAN_CANCER || Some(cancer) == cancer_type {
                let mut pct = vec![0.0; PERCENTILE_COUNT];
                for i in pi..values.len() {
                    if i - pi >= PERCENTILE_COUNT {
                        anyhow::bail!("too many percentile values for gene {gene}");
                    }
                    pct[i - pi] = hmf_common::java_text::parse_double(values[i]).ok_or_else(|| anyhow!("invalid percentile {}", values[i]))?;
                }
                match g.percentiles.iter_mut().find(|(k, _)| &**k == cancer) {
                    Some(e) => e.1 = pct,
                    None => g.percentiles.push((cancer.into(), pct)),
                }
            }
        }
        hmf_common::info!("loaded distribution for {} genes from file({})", c.genes.len(), filename);
        for (k, g) in &c.genes {
            let mut m = f64::NEG_INFINITY;
            for (_, v) in &g.medians {
                m = m.max(*v);
            }
            c.max_medians.insert(k.clone(), if g.medians.is_empty() { 0.0 } else { m });
        }
        Ok(c)
    }
    pub fn has_cancer_type(&self, t: &str) -> bool {
        self.cancer_types.contains(t)
    }
    pub fn tpm_median_across_cancer_types(&self, gene_id: &str) -> f64 {
        self.max_medians.get(gene_id).copied().unwrap_or(0.0)
    }
    pub fn tpm_median(&self, gene_id: &str, cancer: &str) -> f64 {
        self.genes.get(gene_id).and_then(|g| g.medians.iter().find(|(k, _)| &**k == cancer)).map(|x| x.1).unwrap_or(INVALID_VALUE)
    }
    pub fn tpm_percentile(&self, gene_id: &str, cancer: &str, tpm: f64) -> f64 {
        match self.genes.get(gene_id).and_then(|g| g.percentiles.iter().find(|(k, _)| &**k == cancer)) {
            None => INVALID_VALUE,
            Some((_, p)) => get_percentile(p, tpm),
        }
    }
}

/// Percentiles.getPercentile(values, value)
pub fn get_percentile(p: &[f64], value: f64) -> f64 {
    if p.len() != PERCENTILE_COUNT {
        return INVALID_VALUE;
    }
    if value < p[0] {
        return 0.0;
    } else if value > p[p.len() - 1] {
        return (p.len() - 1) as f64 * 0.01;
    }
    let mut same_start: i32 = -1;
    let mut same_end: i32 = 0;
    let n = p.len();
    for i in 0..n - 1 {
        if value == p[i] && p[i] == p[i + 1] {
            if same_start == -1 {
                same_start = i as i32;
            }
            same_end = i as i32 + 1;
            if i < n - 2 {
                continue;
            }
        }
        if value >= p[i] && value <= p[i + 1] {
            if value > p[i] && value == p[i + 1] && i < n - 2 {
                continue;
            }
            if value == p[i] && same_start >= 0 && same_start < same_end {
                return (0.5 * same_start as f64 + 0.5 * same_end as f64) * 0.01;
            }
            if p[i + 1] == p[i] {
                return i as f64 * 0.01;
            }
            let upper = (value - p[i]) / (p[i + 1] - p[i]);
            return (upper * (i + 1) as f64 + (1.0 - upper) * i as f64) * 0.01;
        }
    }
    INVALID_VALUE
}

/// setCohortDistributionValues
pub fn set_cohort_distribution_values(
    summary: &mut GeneCollectionSummary,
    cohort: Option<&CohortGenePercentiles>,
    cancer_type: Option<&str>,
    driver_genes: &HashMap<String, DriverGene>,
) {
    let Some(cohort) = cohort else { return };
    let use_cancer = cancer_type.map(|c| cohort.has_cancer_type(c)).unwrap_or(false);
    for g in &mut summary.gene_results {
        let gid = g.gene.gene_id.as_str();
        let tpm = g.adjusted_tpm;
        let med = cohort.tpm_median(gid, PAN_CANCER);
        let pct = cohort.tpm_percentile(gid, PAN_CANCER, tpm);
        let (med_c, pct_c) = if use_cancer {
            let ct = cancer_type.unwrap();
            (cohort.tpm_median(gid, ct), cohort.tpm_percentile(gid, ct, tpm))
        } else {
            (-1.0, -1.0)
        };
        g.median_tpm_cohort = med;
        g.percentile_tpm_cohort = pct;
        g.median_tpm_cancer = med_c;
        g.percentile_tpm_cancer = pct_c;
        if report_high_expression(g, use_cancer, cohort, driver_genes) {
            g.reported = true;
        }
    }
}

fn report_high_expression(g: &GeneResult, use_cancer: bool, cohort: &CohortGenePercentiles, drivers: &HashMap<String, DriverGene>) -> bool {
    let Some(d) = drivers.get(&g.gene.gene_name) else { return false };
    if !d.report_high_expression {
        return false;
    }
    let (median, pct) = if use_cancer {
        (g.median_tpm_cancer, g.percentile_tpm_cancer)
    } else {
        (cohort.tpm_median_across_cancer_types(&g.gene.gene_id), g.percentile_tpm_cohort)
    };
    if median <= 0.0 || pct <= 0.0 || g.adjusted_tpm <= median {
        return false;
    }
    let lfc = (g.adjusted_tpm / median).ln();
    if lfc > HIGH_EXPRESSION_FOLD_CHANGE_HIGH {
        return true;
    }
    lfc > HIGH_EXPRESSION_FOLD_CHANGE_LOW && pct > HIGH_EXPRESSION_PERCENTILE_HIGH
}

/// PanelTpmNormaliser
pub struct PanelTpmNormaliser {
    factors: HashMap<String, f64>,
    /// the gene ids in Java HashMap order
    pub gene_ids: std::collections::HashSet<String>,
}

impl PanelTpmNormaliser {
    pub fn new(filename: Option<&str>) -> Result<PanelTpmNormaliser> {
        let mut p = PanelTpmNormaliser { factors: HashMap::new(), gene_ids: Default::default() };
        let Some(f) = filename else { return Ok(p) };
        let lines = hmf_common::java_text::read_lines(f)?;
        if lines.is_empty() {
            return Ok(p);
        }
        let delim = if f.ends_with(".tsv") { '\t' } else { ',' };
        let fields = hmf_common::java_text::fields_index_map(&lines[0], delim);
        let gi = *fields.get("GeneId").ok_or_else(|| anyhow!("no GeneId"))?;
        let ai = *fields.get("TpmAdjustFactor").ok_or_else(|| anyhow!("no TpmAdjustFactor"))?;
        for l in &lines[1..] {
            let v: Vec<&str> = hmf_common::java_text::split(l, delim);
            let factor = hmf_common::java_text::parse_double(v[ai]).ok_or_else(|| anyhow!("invalid factor {}", v[ai]))?;
            p.factors.insert(v[gi].to_string(), factor);
            p.gene_ids.insert(v[gi].to_string());
        }
        hmf_common::info!("loaded {} panel gene TPM adjustment factors from file({})", p.factors.len(), f);
        Ok(p)
    }
    pub fn enabled(&self) -> bool {
        !self.factors.is_empty()
    }
    pub fn gene_count(&self) -> usize {
        self.factors.len()
    }
    pub fn apply(&self, summary: &mut GeneCollectionSummary) {
        if !self.enabled() {
            return;
        }
        for g in &mut summary.gene_results {
            let Some(f) = self.factors.get(&g.gene.gene_id) else { continue };
            if g.adjusted_tpm == 0.0 {
                continue;
            }
            if *f <= 0.0 {
                hmf_common::warn!("gene({}:{}) missing or invalid panel TPM adjustment factor({})", g.gene.gene_id, g.gene.gene_name, f);
            } else {
                g.adjusted_tpm /= f;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn category_keys() {
        let c = CategoryCountsData::new(vec![30, 10, 20], vec!["ENSG01".into()]);
        assert_eq!(c.combined_key, "10-20-30-ENSG01");
        let f = CategoryCountsData::from_key("10-20-ENSG02", 3);
        assert_eq!(f.transcripts, vec![10, 20]);
        assert_eq!(f.unspliced_genes, vec!["ENSG02".to_string()]);
        assert!(c.matches(&[10, 20, 30], &["ENSG01"]));
        assert!(!c.matches(&[10, 20], &["ENSG01"]));
    }

    // ExpectationMaxFit on a definition matrix with a unique category per transcript recovers the counts
    #[test]
    fn em_fit() {
        // categories: A only, B only, shared; transcripts A and B
        let defs = vec![0.5, 0.0, 0.5, 0.0, 0.5, 0.5];
        let counts = vec![50.0, 30.0, 80.0];
        let a = perform_fit(&counts, &defs, 3, 2);
        assert!((a[0] - 100.0).abs() < 1.0, "{a:?}");
        assert!((a[1] - 60.0).abs() < 1.0, "{a:?}");
    }

    #[test]
    fn effective_length() {
        let fl = vec![FragmentSize { length: 100, frequency: 1 }, FragmentSize { length: 200, frequency: 3 }];
        assert_eq!(calc_effective_length(150, &fl), 50.0);
        assert_eq!(calc_effective_length(300, &fl), (200.0 + 300.0) / 4.0);
        assert_eq!(calc_effective_length(50, &fl), 0.0);
    }

    // PercentilesTest-style checks of getPercentile
    #[test]
    fn percentiles() {
        let mut p = vec![0.0; 101];
        for (i, x) in p.iter_mut().enumerate() {
            *x = i as f64;
        }
        assert_eq!(get_percentile(&p, 50.0), 0.5);
        assert!((get_percentile(&p, 50.5) - 0.505).abs() < 1e-12);
        assert_eq!(get_percentile(&p, -1.0), 0.0);
        assert_eq!(get_percentile(&p, 101.0), 1.0);
        let z = vec![0.0; 101];
        assert_eq!(get_percentile(&z, 0.0), 0.5);
    }
}
