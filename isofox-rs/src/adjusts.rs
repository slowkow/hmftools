//! The `adjusts` package: GC ratio counts, the expected GC counts of transcripts (GcTranscriptCalculator)
//! and the fragment length distribution (FragmentSizeCalcs).

use crate::bam::{rec_end, rec_mate_start, rec_start, record_cigar, BamReader, SliceFilter};
use crate::common::*;
use crate::config::{FragmentSize, IsofoxConfig, WriteType};
use crate::gene::find_next_overlapping_genes;
use crate::read::{CigarElem, OP_M, OP_N, OP_S};
use anyhow::Result;
use hmf_common::ensembl::{EnsemblDataCache, GeneData, TranscriptData};
use hmf_common::java_hash::JavaHashMap;
use rust_htslib::bam;
use std::collections::HashMap;

pub const GC_RATIO_BUCKET: f64 = 0.01;

/// GcRatioCounts: counts by GC ratio bucket, 0.00 to 1.00
#[derive(Clone, Debug)]
pub struct GcRatioCounts {
    pub ratios: Vec<f64>,
    pub counts: Vec<f64>,
}

impl Default for GcRatioCounts {
    fn default() -> Self {
        Self::new()
    }
}

impl GcRatioCounts {
    pub fn new() -> Self {
        let n = (1.0 / GC_RATIO_BUCKET).floor() as usize + 1;
        let mut ratios = vec![0.0; n];
        let mut g = 0.0;
        for r in ratios.iter_mut() {
            *r = g;
            g += GC_RATIO_BUCKET;
        }
        GcRatioCounts { ratios, counts: vec![0.0; n] }
    }
    pub fn size(&self) -> usize {
        self.ratios.len()
    }
    pub fn counts_total(&self) -> f64 {
        let mut t = 0.0;
        for c in &self.counts {
            t += c;
        }
        t
    }
    pub fn clear_counts(&mut self) {
        self.counts.iter_mut().for_each(|c| *c = 0.0);
    }
    fn ratio_index(&self, gc_ratio: f64) -> i32 {
        let raw = gc_ratio * (self.ratios.len() - 1) as f64;
        let r = java_round(raw).max(0);
        r.min(self.ratios.len() as i64 - 1) as i32
    }
    /// determineRatioData: the one or two buckets of a ratio and their shares
    pub fn determine_ratio_data(&self, gc_ratio: f64, idx: &mut [i32; 2], counts: &mut [f64; 2]) {
        let lower = (gc_ratio / GC_RATIO_BUCKET).floor() * GC_RATIO_BUCKET;
        if lower >= 1.0 {
            idx[0] = self.ratio_index(lower);
            idx[1] = -1;
            counts[0] = 1.0;
            return;
        }
        idx[0] = self.ratio_index(lower);
        idx[1] = idx[0] + 1;
        counts[1] = (gc_ratio - lower) / GC_RATIO_BUCKET;
        counts[0] = 1.0 - counts[1];
    }
    pub fn add_gc_ratio_count(&mut self, idx: i32, count: f64) {
        if idx < 0 || idx as usize >= self.counts.len() {
            return;
        }
        self.counts[idx as usize] += count;
    }
    /// mergeRatioCounts (VectorUtils.sumVectors(other, counts))
    pub fn merge_ratio_counts(&mut self, other: &[f64]) {
        if other.len() != self.counts.len() {
            return;
        }
        for i in 0..other.len() {
            self.counts[i] += other[i];
        }
    }
    pub fn percentile_ratio(&self, percentile: f64) -> f64 {
        let total = self.counts_total();
        let mut current = 0.0;
        let mut prev = 0.0;
        for i in 0..self.ratios.len() {
            let r = self.ratios[i];
            let f = self.counts[i];
            let next = (current + f) / total;
            if next >= percentile {
                return if prev > 0.0 { (prev + r) * 0.5 } else { r };
            }
            current += f;
            prev = r;
        }
        0.0
    }
}

/// Math.round(double)
pub fn java_round(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    (x + 0.5).floor() as i64
}

pub fn calc_gc_count(bases: &[u8]) -> i32 {
    bases.iter().filter(|b| **b == b'G' || **b == b'C').count() as i32
}

#[allow(dead_code)]
pub fn calc_gc_ratio(bases: &[u8]) -> f64 {
    calc_gc_count(bases) as f64 / bases.len() as f64
}

/// SigUtils.convertToPercentages
pub fn convert_to_percentages(v: &mut [f64]) {
    let mut total = 0.0;
    for x in v.iter() {
        total += x;
    }
    if total <= 0.0 {
        return;
    }
    for x in v.iter_mut() {
        *x /= total;
    }
}

/// GcTranscriptCalculator
pub struct GcTranscriptCalculator {
    cache: HashMap<String, Vec<f64>>,
    pub transcript_fit_gc_counts: GcRatioCounts,
    pub gc_ratio_adjustments: Vec<f64>,
}

const MAX_ADJUST_FACTOR: f64 = 3.0;
const ADJUST_LOWER_BOUND: usize = 20;
const ADJUST_UPPER_BOUND: usize = 80;

impl GcTranscriptCalculator {
    pub fn new(config: &IsofoxConfig) -> Result<Self> {
        let counts = GcRatioCounts::new();
        let n = counts.size();
        let mut calc = GcTranscriptCalculator { cache: HashMap::new(), transcript_fit_gc_counts: counts, gc_ratio_adjustments: vec![0.0; n] };
        if let Some(f) = &config.exp_gc_ratios_file {
            calc.load(f)?;
        }
        Ok(calc)
    }

    fn load(&mut self, filename: &str) -> Result<()> {
        if !std::path::Path::new(filename).exists() {
            hmf_common::error!("invalid expected GC ratios file");
            return Ok(());
        }
        let delim = if filename.ends_with(".tsv") || filename.ends_with(".tsv.gz") { '\t' } else { ',' };
        let mut reader = hmf_common::gz::LineReader::open(filename)?;
        let mut first = true;
        let expected = 1 + GcRatioCounts::new().size();
        while let Some(line) = reader.next_line()? {
            let line = std::str::from_utf8(line)?;
            if first {
                first = false;
                continue;
            }
            let values: Vec<&str> = hmf_common::java_text::split_all(line, delim);
            if values.len() != expected {
                hmf_common::error!("invalid exp GC ratio data length({}) vs expected({}): {}", values.len(), expected, line);
                return Ok(());
            }
            let mut freqs = vec![0.0; expected - 1];
            for i in 1..values.len() {
                freqs[i - 1] = hmf_common::java_text::parse_double(values[i]).ok_or_else(|| anyhow::anyhow!("invalid GC ratio value: {}", values[i]))?;
            }
            self.cache.insert(values[0].to_string(), freqs);
        }
        hmf_common::info!("loaded {} transcript expected GC ratios from file({})", self.cache.len(), filename);
        Ok(())
    }

    /// generateGcCountsFromFit for one chromosome's gene collections
    pub fn generate_gc_counts_from_fit(&mut self, fit_allocations: &[(&str, &JavaHashMap<String, f64>)]) {
        let freqs = &mut self.transcript_fit_gc_counts.counts;
        for (gene_names, allocs) in fit_allocations {
            for (trans_name, fit) in allocs.iter() {
                let Some(tf) = self.cache.get(trans_name) else {
                    hmf_common::warn!("genes({}) transcript({}) missing expected GC ratio counts from cache", gene_names, trans_name);
                    // IF4: Java stops here, leaving the chromosome's later gene collections out of the expected
                    // GC counts; corrected, only this transcript is left out
                    if crate::compat::java("IF4") {
                        return;
                    }
                    continue;
                };
                for i in 0..freqs.len() {
                    freqs[i] += tf[i] * fit;
                }
            }
        }
    }

    pub fn calc_gc_ratio_adjustments(&mut self, global: &GcRatioCounts) {
        let expected = &self.transcript_fit_gc_counts.counts;
        let actual = &global.counts;
        let actual_total = global.counts_total();
        let expected_total = self.transcript_fit_gc_counts.counts_total();
        if expected_total == 0.0 || actual_total == 0.0 {
            hmf_common::error!(
                "invalid expected({}) or actual({}) totals",
                hmf_common::fmt::java_fixed(expected_total, 0),
                hmf_common::fmt::java_fixed(actual_total, 0)
            );
            return;
        }
        for i in 0..actual.len() {
            if !(ADJUST_LOWER_BOUND..=ADJUST_UPPER_BOUND).contains(&i) {
                self.gc_ratio_adjustments[i] = 1.0;
                continue;
            }
            let ap = actual[i] / actual_total;
            let ep = expected[i] / expected_total;
            self.gc_ratio_adjustments[i] = if ap > 0.0 && ep > 0.0 {
                (ep / ap).min(MAX_ADJUST_FACTOR).max(1.0 / MAX_ADJUST_FACTOR)
            } else if ap == 0.0 && ep == 0.0 {
                1.0
            } else if ap == 0.0 {
                MAX_ADJUST_FACTOR
            } else {
                1.0 / MAX_ADJUST_FACTOR
            };
        }
    }
}

// ---- FragmentSizeCalcs ----

const MIN_GENE_LENGTH: i32 = 1000;
const MAX_GENE_LENGTH: i32 = 1_000_000;
const MAX_GENE_TRANS: usize = 50;
const MAX_GENE_FRAGMENT_COUNT: i32 = 5000;
const MAX_FRAGMENT_LENGTH: i32 = 5000;

/// the fragment lengths sampled on one chromosome
pub struct FragmentSizeResult {
    pub fragment_lengths: Vec<FragmentSize>,
    pub max_read_length: i32,
    /// FRAG_LENGTH_BY_GENE lines
    pub gene_lines: Vec<String>,
}

fn has_op(c: &[CigarElem], op: u8) -> bool {
    c.iter().any(|e| e.op == op)
}

fn length_bucket(len: i32) -> i32 {
    if len < 1000 {
        return len;
    }
    if len < 3000 {
        return 10 * java_round(len as f64 / 10.0) as i32;
    }
    100 * java_round(len as f64 / 100.0) as i32
}

fn add_fragment_length(len: i32, lengths: &mut Vec<FragmentSize>) -> bool {
    let fl = length_bucket(len);
    if fl == 0 {
        return false;
    }
    let mut index = 0;
    let mut exists = false;
    while index < lengths.len() {
        if lengths[index].length < fl {
            index += 1;
            continue;
        }
        if lengths[index].length == fl {
            lengths[index].frequency += 1;
            exists = true;
        }
        break;
    }
    if !exists {
        lengths.insert(index, FragmentSize { length: fl, frequency: 1 });
    }
    true
}

/// FragmentSizeCalcs.calcSampleFragmentSize for one chromosome
pub fn calc_fragment_sizes(
    config: &IsofoxConfig,
    cache: &'static EnsemblDataCache,
    chromosome: &str,
    genes: &[&'static GeneData],
    required_frag_count: i32,
) -> Result<FragmentSizeResult> {
    let mut result = FragmentSizeResult { fragment_lengths: Vec::new(), max_read_length: 0, gene_lines: Vec::new() };
    if genes.is_empty() || required_frag_count == 0 {
        hmf_common::error!("chromosome({}) fragment size uninitialised", chromosome);
        return Ok(result);
    }
    let Some(bam_file) = &config.bam_file else { return Ok(result) };
    let mut reader = BamReader::open(bam_file, config.ref_genome_file.as_deref())?;
    let filter = SliceFilter { min_map_quality: config.single_map_quality, keep_duplicates: false, keep_supplementaries: false, keep_secondaries: false };
    let has_specific = config.filters.has_specific_regions();
    let by_gene = config.write_type(WriteType::FragLengthByGene);

    let mut overlapping: Vec<&'static GeneData> = Vec::new();
    let mut current = 0;
    let mut processed_fragments = 0;
    let mut trans_list: Vec<&'static TranscriptData> = Vec::new();
    // the read cache (FragmentTracker) of the current genes: read name -> present
    let mut tracker: HashMap<Vec<u8>, ()> = HashMap::new();

    while current < genes.len() {
        current = find_next_overlapping_genes(genes, current, &mut overlapping);
        let v201 = crate::version::v201();
        // ISOFOX 2.0.1 skips the enriched genes
        if v201 && overlapping.iter().any(|g| config.filters.enriched_gene_ids.iter().any(|x| *x == g.gene_id)) {
            continue;
        }
        trans_list.clear();
        tracker.clear();
        let mut range = [0, 0];
        for (i, g) in overlapping.iter().enumerate() {
            range[0] = if i == 0 { g.gene_start } else { g.gene_start.min(range[0]) };
            range[1] = if i == 0 { g.gene_end } else { g.gene_end.max(range[1]) };
            if let Some(t) = cache.transcripts(&g.gene_id) {
                trans_list.extend(t.iter());
            }
        }
        if (v201 || !has_specific) && (trans_list.is_empty() || trans_list.len() > MAX_GENE_TRANS) {
            continue;
        }
        let gene_length = range[1] - range[0];
        if !(MIN_GENE_LENGTH..=MAX_GENE_LENGTH).contains(&gene_length) {
            continue;
        }
        // ISOFOX 2.0.1 skips genes overlapping an excluded region
        if v201 && !config.filters.find_excluded_regions(chromosome, range[0], range[1]).is_empty() {
            continue;
        }

        let mut current_count = 0;
        let mut gene_lengths: Vec<FragmentSize> = Vec::new();
        let max_read_length = &mut result.max_read_length;
        let lengths = &mut result.fragment_lengths;
        reader.slice(chromosome, range[0], range[1], &filter, |rec: &bam::Record| {
            if !is_candidate_record(rec, config, *max_read_length, range, &trans_list) {
                return true;
            }
            current_count += 1;
            if current_count >= MAX_GENE_FRAGMENT_COUNT {
                return false;
            }
            *max_read_length = (*max_read_length).max(rec.seq_len() as i32);
            if rec.is_paired() {
                let name = rec.qname().to_vec();
                if tracker.remove(&name).is_none() {
                    tracker.insert(name, ());
                    return true;
                }
            }
            let isize = inferred_insert_size_abs(rec);
            if add_fragment_length(isize, lengths) {
                processed_fragments += 1;
            }
            if by_gene && add_fragment_length(isize, &mut gene_lengths) {
                processed_fragments += 1;
            }
            true
        })?;

        if by_gene && !gene_lengths.is_empty() {
            let mut names = overlapping[0].gene_name.clone();
            for g in overlapping.iter().take(10).skip(1) {
                names.push(';');
                names.push_str(&g.gene_name);
            }
            for fl in &gene_lengths {
                result.gene_lines.push(format!("{},{},{},{},{},{},{}", names, overlapping.len(), chromosome, range[0], range[1], fl.length, fl.frequency));
            }
        }
        if processed_fragments >= required_frag_count {
            break;
        }
    }
    Ok(result)
}

/// SamRecordUtils.inferredInsertSizeAbs
fn inferred_insert_size_abs(rec: &bam::Record) -> i32 {
    if rec.is_paired() { (rec.insert_size() as i32).wrapping_abs() } else { (rec_end(rec) - rec_start(rec) + 1).wrapping_abs() }
}

fn is_candidate_record(rec: &bam::Record, config: &IsofoxConfig, max_read_length: i32, range: [i32; 2], trans: &[&TranscriptData]) -> bool {
    let paired = rec.is_paired();
    // ISOFOX 2.0.1 compares the mate's chromosome even for an unpaired read ("*"), so never takes one
    if !paired && crate::version::v201() {
        return false;
    }
    let isize = inferred_insert_size_abs(rec);
    if isize > MAX_FRAGMENT_LENGTH {
        return false;
    }
    if paired && (rec.mtid() != rec.tid() || (rec.is_mate_reverse() == rec.is_reverse())) {
        // mateNegativeStrand is false for an unpaired read; the read is paired here
        return false;
    }
    let cigar = record_cigar(rec);
    if has_op(&cigar, OP_N) || !has_op(&cigar, OP_M) {
        return false;
    }
    if !paired && has_op(&cigar, OP_S) {
        return false;
    }
    let read_length = max_read_length.max(config.read_length);
    if read_length > 0 && isize > read_length && has_op(&cigar, OP_S) && cigar.iter().any(|e| e.op == OP_S && e.len > 2) {
        return false;
    }
    if paired && !position_within(rec_mate_start(rec), range[0], range[1]) {
        return false;
    }
    let (s, e) = (rec_start(rec), rec_end(rec));
    for t in trans {
        if t.exons.iter().any(|x| positions_overlap(s, e, x.start, x.end)) {
            return false;
        }
    }
    true
}

/// FragmentSizeCalcs.mergeData
pub fn merge_fragment_lengths(all: &mut Vec<FragmentSize>, other: &[FragmentSize]) {
    for o in other {
        let mut index = 0;
        let mut exists = false;
        while index < all.len() {
            if all[index].length < o.length {
                index += 1;
                continue;
            }
            if all[index].length == o.length {
                all[index].frequency += o.frequency;
                exists = true;
            }
            break;
        }
        if !exists {
            all.insert(index, FragmentSize { length: o.length, frequency: o.frequency });
        }
    }
}

/// setConfigFragmentLengthData: the observed frequencies of the expected-rate fragment lengths
pub fn set_config_fragment_length_data(config: &mut IsofoxConfig, lengths: &[FragmentSize]) {
    let n = config.fragment_size_data.len();
    let mut range_max = 0;
    for i in 0..n {
        let range_min = if i == 0 { 0 } else { range_max + 1 };
        range_max = if i == n - 1 {
            config.max_fragment_length - 1
        } else {
            (config.fragment_size_data[i].length + config.fragment_size_data[i + 1].length) / 2
        };
        let mut count = 0;
        for fl in lengths {
            if fl.length >= range_min && fl.length <= range_max {
                count += fl.frequency;
            }
        }
        config.fragment_size_data[i].frequency = count;
        hmf_common::debug!("fragmentLength({}) frequency({})", config.fragment_size_data[i].length, count);
    }
}

/// calcPercentileData
pub fn calc_percentile_data(lengths: &[FragmentSize], percentiles: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    let total: i64 = lengths.iter().map(|x| x.frequency as i64).sum();
    let total = total as f64;
    for p in percentiles {
        let mut current: i32 = 0;
        let mut prev = 0;
        for fl in lengths {
            let next = (current + fl.frequency) as f64 / total;
            if next >= *p {
                let v = if prev > 0 { (prev + fl.length) as f64 * 0.5 } else { fl.length as f64 };
                out.push(v);
                break;
            }
            current += fl.frequency;
            prev = fl.length;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_ratio_buckets() {
        let g = GcRatioCounts::new();
        assert_eq!(g.size(), 101);
        let mut idx = [-1, -1];
        let mut c = [0.0, 0.0];
        g.determine_ratio_data(0.425, &mut idx, &mut c);
        assert_eq!(idx, [42, 43]);
        assert!((c[1] - 0.5).abs() < 1e-9);
        g.determine_ratio_data(1.0, &mut idx, &mut c);
        assert_eq!(idx, [100, -1]);
    }

    #[test]
    fn fragment_length_percentiles() {
        let mut l = Vec::new();
        for len in [100, 150, 150, 200, 1004, 1006] {
            add_fragment_length(len, &mut l);
        }
        assert_eq!(l.iter().map(|x| (x.length, x.frequency)).collect::<Vec<_>>(), vec![(100, 1), (150, 2), (200, 1), (1000, 1), (1010, 1)]);
        let p = calc_percentile_data(&l, &[0.05, 0.5, 0.95]);
        assert_eq!(p, vec![100.0, 125.0, 1005.0]);
    }
}
