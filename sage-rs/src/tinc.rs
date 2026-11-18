//! TINC (tumor in normal contamination): TincAnalyser, TincCalculator, VariantCache, TincVcfWriter.
//!
//! Runs over a finished SAGE VCF: estimates the fraction of tumor DNA in the reference sample from
//! the reference alt support of likely-somatic SNVs, then, if that level is above zero, re-tests the
//! germline soft filters with the implied tumor contribution backed out and recovers variants that
//! now pass. The VCF is rewritten with gnomAD / PON annotations. Variants outside chr1-22,X,Y are
//! dropped from the rewritten file, as in Java (it only queries those chromosomes).

use crate::constants::{jround, LONG_INSERT_LENGTH};
use crate::filters::{
    above_max_germline_relative_qual, above_max_germline_vaf, above_max_mnv_indel_germline_alt_support,
    include_ref_jitter_in_msi_indel,
    is_panel_or_long_indel_repeat_variant, tier_config,
    ALL_FILTERS,
};
use crate::fmt::{java_fixed, vcf_double};
use crate::stats::{binomial_cdf, binomial_probability};
use crate::types::{chromosome_rank, is_human_chromosome, Tier};
use crate::vcfout::{index_vcf, sorted_header_lines};
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use rust_htslib::bgzf;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const TINC_MIN_VARIANTS: usize = 200;
const TINC_GERMLINE_ABQ_MIN: i32 = 30;
const TINC_MAX_FITTING_VARIANTS: usize = 15_000;
const TINC_GERMLINE_DEPTH_LOW: f64 = 0.5;
const TINC_GERMLINE_DEPTH_HIGH: f64 = 1.0 + TINC_GERMLINE_DEPTH_LOW;
const TINC_MQF_LIMIT: f64 = 25.0;
const TINC_TUMOR_AF_LOWER_LIMIT: f64 = 50.0;
const TINC_TUMOR_AF_UPPER_LIMIT: f64 = 300.0;
const TINC_TUMOR_AF_UPPER_TEST_MIN: f64 = 0.011;
pub const TINC_GERMLINE_MAX_AD: i32 = 10;
const TINC_RECOVERY_FACTOR: f64 = 2.5;
const TINC_RECOVERY_MIN: f64 = 0.03;
const TINC_RECOVERY_GERMLINE_AF_PROB: f64 = 0.001;
const DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT: f64 = 3.0;

const F_REL_QUAL: &str = "maxGermlineRelQual";
const F_VAF: &str = "maxGermlineVAF";
const F_ALT_SUPPORT: &str = "maxGermlineAltSupport";
const RECOVERY_FILTERS: [&str; 4] = [F_REL_QUAL, F_VAF, F_ALT_SUPPORT, "PASS"];

pub const TINC_HEADER_LINES: [&str; 6] = [
    "##INFO=<ID=TINC_RECOVERED,Number=0,Type=Flag,Description=\"Variant recovered from germline filters by TINC detection\">",
    "##INFO=<ID=PON_COUNT,Number=1,Type=Integer,Description=\"Cohort frequency for variant\">",
    "##INFO=<ID=PON_MAX,Number=1,Type=Integer,Description=\"Max read depth in any sample with variant\">",
    "##INFO=<ID=PON_READS,Number=1,Type=Integer,Description=\"Average reads in sample with variant\">",
    "##INFO=<ID=PON_MULTI,Number=1,Type=String,Description=\"Status of the variant across PONs\">",
    "##INFO=<ID=GND_FREQ,Number=1,Type=Float,Description=\"Gnomad variant frequency\">",
];

#[derive(Clone, Debug)]
pub struct TincConfig {
    pub v38: bool,
    pub input_vcf: String,
    pub tumor_id: String,
    pub reference_id: String,
    pub pon_file: Option<String>,
    pub gnomad_dir: Option<String>,
    pub gnomad_file: Option<String>,
    pub threads: usize,
    pub write_fit_variants: bool,
    pub fit_variants_file: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FilterReason {
    None,
    Indel,
    Filtered,
    Gnomad,
    Pon,
    LowAbq,
    GermlineAf,
    Mqf,
    NearIndel,
    RefDepth,
}

const FILTER_REASONS: [(FilterReason, &str); 10] = [
    (FilterReason::None, "NONE"),
    (FilterReason::Indel, "INDEL"),
    (FilterReason::Filtered, "FILTERED"),
    (FilterReason::Gnomad, "GNOMAD"),
    (FilterReason::Pon, "PON"),
    (FilterReason::LowAbq, "LOW_ABQ"),
    (FilterReason::GermlineAf, "GERMLINE_AF"),
    (FilterReason::Mqf, "MQF"),
    (FilterReason::NearIndel, "NEAR_INDEL"),
    (FilterReason::RefDepth, "REF_DEPTH"),
];

/// tinc.VariantData: the parts of a VCF record TINC uses, plus its annotations and recovery state.
#[derive(Clone, Debug, Default)]
pub struct VariantData {
    pub position: i32,
    pub ref_: String,
    pub alt: String,
    pub reference_depth: i32,
    pub reference_alt_frags: i32,
    pub tumor_depth: i32,
    pub tumor_alt_frags: i32,
    /// false when loaded from a fit-variants file (no VCF record behind it)
    pub has_context: bool,
    pub filters: Vec<String>,
    pub tier: Option<Tier>,
    /// QUAL and RC_REPC, for the indel repeat test in the SBX germline thresholds
    pub qual: f64,
    pub repeat_count: i32,
    /// the length of RC_REPS, for the Ultima reference jitter test
    pub repeat_bases_length: i32,
    pub map_qual_factor: f64,
    pub nearby_indel: bool,
    pub has_lps: bool,
    /// ARCBQ[1] of each sample
    pub ref_abq: i32,
    pub tumor_abq: i32,
    /// RC_CNT / RC_QUAL: Full, PartialCore, Core, Realigned, Ref, Total
    pub ref_counts: [i32; 6],
    pub ref_quals: [i32; 6],
    pub tumor_counts: [i32; 6],
    pub tumor_quals: [i32; 6],
    pub ref_simple_alt: i32,
    pub ref_jitter: [i32; 2],

    pub pon_filtered: bool,
    pub pon_sample_count: i32,
    pub pon_max_read_count: i32,
    pub pon_mean_read_count: i32,
    pub gnomad_frequency: Option<f64>,

    pub reference_alt_frags_reduction: f64,
    pub new_filters: Option<Vec<String>>,
}

fn int_list<const N: usize>(s: Option<&str>) -> [i32; N] {
    let mut out = [0; N];
    if let Some(s) = s {
        for (i, v) in s.split(',').enumerate().take(N) {
            out[i] = v.parse().unwrap_or(0);
        }
    }
    out
}

fn alt_support(c: &[i32; 6]) -> i32 {
    c[0] + c[1] + c[2] + c[3]
}

impl VariantData {
    /// the fit-variants-file constructor
    pub fn simple(position: i32, ref_: &str, alt: &str, ref_depth: i32, ref_alt: i32, tumor_depth: i32, tumor_alt: i32) -> Self {
        VariantData {
            position,
            ref_: ref_.to_string(),
            alt: alt.to_string(),
            reference_depth: ref_depth,
            reference_alt_frags: ref_alt,
            tumor_depth,
            tumor_alt_frags: tumor_alt,
            ..Default::default()
        }
    }

    /// From the tab-split fields of a VCF data line; `ref_col` / `tumor_col` are the sample column indices.
    pub fn from_fields(f: &[&str], ref_col: usize, tumor_col: usize) -> Result<Self> {
        if f.len() <= ref_col.max(tumor_col) {
            bail!("VCF line has too few columns");
        }
        let keys: Vec<&str> = f[8].split(':').collect();
        let geno = |col: usize, key: &str| -> Option<&str> {
            let idx = keys.iter().position(|k| *k == key)?;
            f[col].split(':').nth(idx).filter(|v| *v != ".")
        };
        let dp = |col: usize| geno(col, "DP").and_then(|v| v.parse().ok()).unwrap_or(-1);
        let ad1 = |col: usize| int_list::<2>(geno(col, "AD"))[1];
        let mut v = VariantData {
            position: f[1].parse().context("VCF POS")?,
            ref_: f[3].to_string(),
            alt: f[4].split(',').next().unwrap_or("").to_string(),
            reference_depth: dp(ref_col),
            reference_alt_frags: ad1(ref_col),
            tumor_depth: dp(tumor_col),
            tumor_alt_frags: ad1(tumor_col),
            has_context: true,
            qual: f[5].parse().unwrap_or(0.0),
            filters: f[6].split(';').filter(|x| !x.is_empty() && *x != "PASS" && *x != ".").map(|x| x.to_string()).collect(),
            ref_abq: int_list::<2>(geno(ref_col, "ARCBQ"))[1],
            tumor_abq: int_list::<2>(geno(tumor_col, "ARCBQ"))[1],
            ref_counts: int_list(geno(ref_col, "RC_CNT")),
            ref_quals: int_list(geno(ref_col, "RC_QUAL")),
            tumor_counts: int_list(geno(tumor_col, "RC_CNT")),
            tumor_quals: int_list(geno(tumor_col, "RC_QUAL")),
            ref_simple_alt: geno(ref_col, "SAC").and_then(|v| v.parse().ok()).unwrap_or(0),
            ref_jitter: int_list(geno(ref_col, "RC_JIT")),
            ..Default::default()
        };
        for item in f[7].split(';') {
            let (key, value) = item.split_once('=').unwrap_or((item, ""));
            match key {
                "TIER" => v.tier = tier_from_str(value),
                "MQF" => v.map_qual_factor = value.parse().unwrap_or(0.0),
                "RC_REPC" => v.repeat_count = value.parse().unwrap_or(0),
                "RC_REPS" => v.repeat_bases_length = value.len() as i32,
                "NEARBY_INDEL" => v.nearby_indel = true,
                "LPS" => v.has_lps = true,
                _ => {}
            }
        }
        Ok(v)
    }

    pub fn is_indel(&self) -> bool {
        self.alt.len() != self.ref_.len()
    }
    pub fn is_insert(&self) -> bool {
        self.alt.len() > self.ref_.len()
    }
    /// As Java: requires alt to equal ref, so this is never true for a real variant.
    pub fn is_mnv(&self) -> bool {
        self.ref_.len() > 1 && self.alt == self.ref_
    }
    pub fn is_long_indel(&self) -> bool {
        self.ref_.len().max(self.alt.len()) as i32 >= LONG_INSERT_LENGTH
    }
    /// VariantTier.UNKNOWN is treated as the low confidence tier wherever TINC uses the tier
    pub fn tier(&self) -> Tier {
        self.tier.unwrap_or(Tier::LowConfidence)
    }
    pub fn is_passing(&self) -> bool {
        !self.has_context || self.filters.is_empty()
    }
    pub fn only_germline_filtered(&self) -> bool {
        self.has_context && self.filters.iter().all(|x| RECOVERY_FILTERS.contains(&x.as_str()))
    }
    pub fn recovered(&self) -> bool {
        matches!(&self.new_filters, Some(f) if f.is_empty())
    }
    pub fn tumor_af(&self) -> f64 {
        if self.tumor_depth > 0 {
            self.tumor_alt_frags as f64 / self.tumor_depth as f64
        } else {
            0.0
        }
    }
    pub fn reference_af(&self) -> f64 {
        if self.reference_depth > 0 {
            self.reference_alt_frags as f64 / self.reference_depth as f64
        } else {
            0.0
        }
    }
    pub fn calc_reduced_alt_value(&self, alt_value: f64) -> f64 {
        let reference_af = self.reference_af();
        if reference_af <= 0.0 {
            return 0.0;
        }
        let reduced_af = (reference_af - self.reference_alt_frags_reduction).max(0.0);
        alt_value * (reduced_af / reference_af)
    }
    pub fn calc_reduced_alt_count(&self, alt_count: i32) -> i32 {
        jround(self.calc_reduced_alt_value(alt_count as f64)) as i32
    }
    fn build_new_filters(&mut self) {
        self.new_filters =
            Some(self.filters.iter().filter(|x| ALL_FILTERS.iter().any(|(_, n, _)| n == x)).cloned().collect());
    }
    fn has_new_filter(&self, name: &str) -> bool {
        self.new_filters.as_ref().is_some_and(|f| f.iter().any(|x| x == name))
    }
    fn remove_new_filter(&mut self, name: &str) {
        if let Some(f) = self.new_filters.as_mut() {
            f.retain(|x| x != name);
        }
    }
}

fn tier_from_str(s: &str) -> Option<Tier> {
    match s {
        "HOTSPOT" => Some(Tier::Hotspot),
        "PANEL" => Some(Tier::Panel),
        "HIGH_CONFIDENCE" => Some(Tier::HighConfidence),
        "LOW_CONFIDENCE" => Some(Tier::LowConfidence),
        _ => None,
    }
}

/// VariantCache.ChromosomeTask.checkFilters: is this variant usable for the TINC fit
pub fn check_filters(v: &VariantData) -> FilterReason {
    let mut reason = FilterReason::None;
    if v.is_indel() {
        reason = FilterReason::Indel;
    } else if !v.is_passing() && !v.only_germline_filtered() {
        reason = FilterReason::Filtered;
    } else if v.reference_alt_frags > 0 && v.ref_abq < TINC_GERMLINE_ABQ_MIN {
        reason = FilterReason::LowAbq;
    }
    if reason == FilterReason::None {
        if v.map_qual_factor < TINC_MQF_LIMIT {
            reason = FilterReason::Mqf;
        } else if v.nearby_indel {
            reason = FilterReason::NearIndel;
        }
    }
    if reason == FilterReason::None {
        if v.gnomad_frequency.is_some() {
            reason = FilterReason::Gnomad;
        } else if v.pon_filtered {
            reason = FilterReason::Pon;
        }
    }
    if reason == FilterReason::None
        && v.reference_alt_frags > 0
        && binomial_cdf(v.reference_depth, 0.5, v.reference_alt_frags) >= TINC_RECOVERY_GERMLINE_AF_PROB
    {
        reason = FilterReason::GermlineAf;
    }
    reason
}

// ---------------------------------------------------------------------------------------------
// TincCalculator

/// TincCalculator.populateDefaultLevels, with the same floating point accumulation
pub fn populate_default_levels() -> Vec<f64> {
    let rounded = |level: f64, increment: f64| (jround(level / increment) as i32) as f64 * increment;
    let mut levels = Vec::new();
    let mut level = 0.0;
    while level < 0.03 {
        levels.push(rounded(level, 0.005));
        level += 0.005;
    }
    level = 0.03;
    while level < 0.15 {
        levels.push(rounded(level, 0.01));
        level += 0.01;
    }
    level = 0.15;
    while level <= 1.0 {
        levels.push(level);
        level += 0.05;
    }
    levels
}

fn compute_contamination_score(variants: &[&VariantData], test_level: f64) -> f64 {
    let mut fit_penalty = 0.0;
    for test_ad in 0..=TINC_GERMLINE_MAX_AD {
        let mut probability_total = 0.0;
        let mut variant_count = 0;
        for v in variants {
            let expected_af = test_level * v.tumor_af();
            probability_total += binomial_probability(v.reference_depth, expected_af, test_ad);
            if v.reference_alt_frags == test_ad
                || (test_ad == TINC_GERMLINE_MAX_AD && v.reference_alt_frags > TINC_GERMLINE_MAX_AD)
            {
                variant_count += 1;
            }
        }
        fit_penalty += (variant_count as f64 - probability_total).abs();
    }
    fit_penalty
}

/// TincCalculator.calculate. Java's early exit never triggers (its counter is never incremented),
/// so every level is scored.
pub fn calculate(variants: &[&VariantData], test_levels: &[f64]) -> f64 {
    if variants.len() < TINC_MIN_VARIANTS {
        return 0.0;
    }
    let mut tumor_af_total = 0.0;
    for v in variants {
        tumor_af_total += v.tumor_af();
    }
    if tumor_af_total < TINC_TUMOR_AF_LOWER_LIMIT {
        return 0.0;
    }
    let exclude_lower_levels = tumor_af_total < TINC_TUMOR_AF_UPPER_LIMIT;
    let mut lowest_score = -1.0;
    let mut lowest_level = -1.0;
    for &level in test_levels {
        if exclude_lower_levels && level > 0.0 && level < TINC_TUMOR_AF_UPPER_TEST_MIN {
            continue;
        }
        let score = compute_contamination_score(variants, level);
        if lowest_score < 0.0 || score < lowest_score {
            lowest_score = score;
            lowest_level = level;
        }
    }
    eprintln!(
        "tinc level({}) tumorAFTotal({}) minPenalty({})",
        java_fixed(lowest_level, 3),
        java_fixed(tumor_af_total, 0),
        java_fixed(lowest_score, 0)
    );
    lowest_level
}

// ---------------------------------------------------------------------------------------------
// TincAnalyser.recoverVariants

fn recheck_max_germline_vaf(v: &mut VariantData) {
    if !v.has_new_filter(F_VAF) {
        return;
    }
    let mut adjusted_ref_alt = alt_support(&v.ref_counts) + v.ref_simple_alt;
    let indel_length_abs = (v.ref_.len() as i32 - v.alt.len() as i32).abs();
    if v.is_long_indel() || include_ref_jitter_in_msi_indel(v.repeat_count, v.repeat_bases_length, indel_length_abs) {
        adjusted_ref_alt += v.ref_jitter[0] + v.ref_jitter[1];
    }
    adjusted_ref_alt = v.calc_reduced_alt_count(adjusted_ref_alt);
    let max_germline_vaf = tier_config(v.tier()).max_germline_vaf;
    let valid_indel_repeat = is_panel_or_long_indel_repeat_variant(v.tier(), v.qual, v.is_indel(), v.repeat_count);
    if above_max_germline_vaf(v.tier(), valid_indel_repeat, v.tumor_af(), adjusted_ref_alt, v.ref_counts[5], max_germline_vaf) {
        return;
    }
    v.remove_new_filter(F_VAF);
}

fn recheck_max_germline_relative_qual(v: &mut VariantData) {
    if !v.has_new_filter(F_REL_QUAL) {
        return;
    }
    let tumor_qual = (v.tumor_quals[0] + v.tumor_quals[1] + v.tumor_quals[3]) as f64;
    let ref_qual = v.calc_reduced_alt_value((v.ref_quals[0] + v.ref_quals[1] + v.ref_quals[3]) as f64);
    if above_max_germline_relative_qual(
        v.tier(),
        tumor_qual,
        v.tumor_af(),
        v.tumor_counts[5],
        v.tumor_abq as f64,
        ref_qual,
        v.ref_counts[5],
        alt_support(&v.ref_counts),
        v.ref_abq as f64,
        is_panel_or_long_indel_repeat_variant(v.tier(), v.qual, v.is_indel(), v.repeat_count),
    ) {
        return;
    }
    v.remove_new_filter(F_REL_QUAL);
}

fn recheck_max_germline_alt_support(v: &mut VariantData) {
    if !v.has_new_filter(F_ALT_SUPPORT) {
        return;
    }
    let is_long_insert = v.is_insert() && v.is_long_indel();
    let new_alt_support = v.calc_reduced_alt_count(alt_support(&v.ref_counts));
    if above_max_mnv_indel_germline_alt_support(v.tier(), v.is_mnv(), is_long_insert, v.ref_counts[5], new_alt_support) {
        return;
    }
    v.remove_new_filter(F_ALT_SUPPORT);
}

/// Returns the number of variants recovered.
pub fn recover_variants(variants: &mut [VariantData], tinc_level: f64) -> usize {
    let recovery_tinc = TINC_RECOVERY_MIN + tinc_level * TINC_RECOVERY_FACTOR;
    let mut recovered = 0;
    for v in variants.iter_mut() {
        if v.is_passing() {
            continue;
        }
        if !v.filters.iter().any(|x| RECOVERY_FILTERS.contains(&x.as_str())) {
            continue;
        }
        if v.pon_filtered {
            continue;
        }
        v.build_new_filters();
        v.reference_alt_frags_reduction = (v.tumor_af() * recovery_tinc).min(1.0);
        recheck_max_germline_vaf(v);
        recheck_max_germline_relative_qual(v);
        recheck_max_germline_alt_support(v);
        if v.recovered() {
            recovered += 1;
        }
    }
    recovered
}

/// TincAnalyser.filterOutVariant: drop variants kept only because TINC raised the germline alt support limit
fn filter_out_variant(v: &VariantData) -> bool {
    let reduced = v.calc_reduced_alt_value(v.reference_alt_frags as f64);
    if reduced <= DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT {
        return false;
    }
    !v.has_lps
}

// ---------------------------------------------------------------------------------------------
// gnomAD and PON resources, loaded only at the positions present in the VCF

struct PonEntry {
    ref_: String,
    alt: String,
    samples: i32,
    max_reads: i32,
    total_reads: i32,
}
type Positions = HashMap<String, HashSet<i32>>;
type PonMap = HashMap<String, HashMap<i32, Vec<PonEntry>>>;
type GnomadMap = HashMap<String, HashMap<i32, Vec<(String, String, f64)>>>;

fn open_text(path: &str) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path).with_context(|| format!("opening {path}"))?;
    let reader: Box<dyn Read + Send> = if path.ends_with(".gz") { Box::new(MultiGzDecoder::new(file)) } else { Box::new(file) };
    Ok(Box::new(BufReader::with_capacity(1 << 20, reader)))
}

fn parse_i32(b: &[u8]) -> Option<i32> {
    std::str::from_utf8(b).ok()?.parse().ok()
}

fn trim_eol(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}

fn load_pon(path: &str, wanted: &Positions) -> Result<PonMap> {
    let mut reader = open_text(path)?;
    let mut line = Vec::new();
    reader.read_until(b'\n', &mut line)?;
    let header = String::from_utf8_lossy(trim_eol(&line)).to_string();
    let cols: Vec<&str> = header.split('\t').collect();
    if cols.len() < 4 {
        bail!("pon file({path}) has insufficient column count({})", cols.len());
    }
    let col = |names: &[&str]| names.iter().find_map(|n| cols.iter().position(|c| c == n));
    let need = |name: &str| col(&[name]).with_context(|| format!("PON file {path} has no {name} column"));
    let (chr_i, pos_i, ref_i, alt_i) = (need("Chromosome")?, need("Position")?, need("Ref")?, need("Alt")?);
    let samples_i = col(&["SampleCount", "SamplesCount"]);
    let max_i = col(&["MaxSampleReadCount", "MaxSampleReads"]);
    let total_i = col(&["TotalReadCount", "TotalReads"]);

    let mut map: PonMap = HashMap::new();
    let mut current_chr: Vec<u8> = Vec::new();
    let mut current_wanted: Option<&HashSet<i32>> = None;
    let mut current_name = String::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let f: Vec<&[u8]> = trim_eol(&line).split(|b| *b == b'\t').collect();
        if f.len() <= chr_i.max(pos_i).max(ref_i).max(alt_i) {
            continue;
        }
        if f[chr_i] != current_chr.as_slice() {
            current_chr = f[chr_i].to_vec();
            current_name = String::from_utf8_lossy(&current_chr).to_string();
            current_wanted = wanted.get(&current_name);
        }
        let Some(set) = current_wanted else { continue };
        let Some(pos) = parse_i32(f[pos_i]) else { bail!("invalid PON position in {path}") };
        if !set.contains(&pos) {
            continue;
        }
        let int = |i: Option<usize>| i.and_then(|i| f.get(i)).and_then(|b| parse_i32(b)).unwrap_or(0);
        map.entry(current_name.clone()).or_default().entry(pos).or_default().push(PonEntry {
            ref_: String::from_utf8_lossy(f[ref_i]).to_string(),
            alt: String::from_utf8_lossy(f[alt_i]).to_string(),
            samples: int(samples_i),
            max_reads: int(max_i),
            total_reads: int(total_i),
        });
    }
    Ok(map)
}

/// One gnomAD csv; `file_chromosome` is set for per-chromosome files, which have no chromosome column.
fn load_gnomad_file(path: &str, file_chromosome: Option<&str>, wanted: &Positions, out: &mut GnomadMap) -> Result<()> {
    let mut reader = open_text(path)?;
    let mut line = Vec::new();
    reader.read_until(b'\n', &mut line)?; // header
    let offset = if file_chromosome.is_some() { 0 } else { 1 };
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let f: Vec<&[u8]> = trim_eol(&line).split(|b| *b == b',').collect();
        if f.len() < offset + 4 {
            continue;
        }
        let chr = match file_chromosome {
            Some(c) => c.to_string(),
            None => String::from_utf8_lossy(f[0]).to_string(),
        };
        let Some(set) = wanted.get(&chr) else { continue };
        let Some(pos) = parse_i32(f[offset]) else { bail!("invalid gnomAD position in {path}") };
        if !set.contains(&pos) {
            continue;
        }
        let freq: f64 = std::str::from_utf8(f[offset + 3]).ok().and_then(|s| s.parse().ok()).context("invalid gnomAD frequency")?;
        out.entry(chr).or_default().entry(pos).or_default().push((
            String::from_utf8_lossy(f[offset + 1]).to_string(),
            String::from_utf8_lossy(f[offset + 2]).to_string(),
            freq,
        ));
    }
    Ok(())
}

/// GnomadCache.loadAllFrequencyFiles: the file for a chromosome, e.g. gnomad_variants_chr10_v38.csv.gz
fn gnomad_chromosome_files(dir: &str, chromosomes: &[String]) -> Result<Vec<(String, String)>> {
    let mut files: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(dir).with_context(|| format!("reading gnomAD directory {dir}"))? {
        let path = entry?.path();
        if !path.is_dir() && path.to_string_lossy().contains("gnomad_variants") {
            files.push(path.to_string_lossy().to_string());
        }
    }
    files.sort();
    let mut out = Vec::new();
    for chr in chromosomes {
        let bare = chr.strip_prefix("chr").unwrap_or(chr);
        let no_id = format!("gnomad_variants_chr{bare}.csv");
        let with_id = format!("gnomad_variants_chr{bare}_");
        match files.iter().find(|f| f.ends_with(&no_id) || f.contains(&with_id)) {
            Some(f) => out.push((chr.clone(), f.clone())),
            None => bail!("missing Gnomad chromosome({chr}) file"),
        }
    }
    Ok(out)
}

fn load_resources(cfg: &TincConfig, wanted: &Positions, chromosomes: &[String]) -> Result<(PonMap, GnomadMap)> {
    let gnomad_jobs: Vec<(Option<String>, String)> = if let Some(file) = &cfg.gnomad_file {
        vec![(None, file.clone())]
    } else if let Some(dir) = &cfg.gnomad_dir {
        gnomad_chromosome_files(dir, chromosomes)?.into_iter().map(|(c, f)| (Some(c), f)).collect()
    } else {
        Vec::new()
    };
    let next = AtomicUsize::new(0);
    let gnomad: Mutex<GnomadMap> = Mutex::new(HashMap::new());
    let mut pon: Result<PonMap> = Ok(HashMap::new());
    let errors: Mutex<Vec<anyhow::Error>> = Mutex::new(Vec::new());
    std::thread::scope(|scope| {
        let pon_handle = cfg.pon_file.as_ref().map(|path| scope.spawn(move || load_pon(path, wanted)));
        for _ in 0..cfg.threads.max(1).min(gnomad_jobs.len()) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                let Some((chr, file)) = gnomad_jobs.get(i) else { break };
                let mut local: GnomadMap = HashMap::new();
                match load_gnomad_file(file, chr.as_deref(), wanted, &mut local) {
                    Ok(()) => gnomad.lock().unwrap().extend(local),
                    Err(e) => errors.lock().unwrap().push(e),
                }
            });
        }
        if let Some(h) = pon_handle {
            pon = h.join().unwrap();
        }
    });
    if let Some(e) = errors.into_inner().unwrap().pop() {
        return Err(e);
    }
    Ok((pon?, gnomad.into_inner().unwrap()))
}

/// PonCache.filterOnTierCriteria with PON_FILTERS_V37 / V38
fn pon_filter_on_tier(v38: bool, tier: Option<Tier>, samples: i32, max_reads: i32) -> bool {
    let (req_samples, req_reads) = match (v38, tier) {
        (true, Some(Tier::Hotspot)) => (6, 5),
        (true, Some(Tier::Panel)) => (3, 3),
        (true, _) => (3, 0),
        (false, Some(Tier::Hotspot)) => (10, 5),
        (false, Some(Tier::Panel)) => (6, 5),
        (false, _) => (6, 0),
    };
    samples >= req_samples && max_reads >= req_reads
}

// ---------------------------------------------------------------------------------------------
// VariantCache

fn downsample_fitting(fitting: &mut Vec<usize>) {
    let nth = fitting.len() / TINC_MAX_FITTING_VARIANTS;
    if nth < 2 {
        return;
    }
    let start = fitting.len();
    let mut kept = Vec::with_capacity(start / nth + 1);
    let mut counter = 0;
    for &i in fitting.iter() {
        counter += 1;
        if counter >= nth {
            kept.push(i);
            counter = 0;
        }
    }
    *fitting = kept;
    eprintln!("down-sampled variants({start} -> {})", fitting.len());
}

fn filter_by_average_reference_depth(variants: &[VariantData], fitting: &mut Vec<usize>) -> usize {
    let depth_total: i64 = fitting.iter().map(|&i| variants[i].reference_depth as i64).sum();
    if depth_total == 0 {
        return 0;
    }
    let average = jround(depth_total as f64 / fitting.len() as f64) as i32;
    let before = fitting.len();
    fitting.retain(|&i| {
        let d = variants[i].reference_depth as f64;
        !(d < TINC_GERMLINE_DEPTH_LOW * average as f64 || d > TINC_GERMLINE_DEPTH_HIGH * average as f64)
    });
    let removed = before - fitting.len();
    if removed > 0 {
        eprintln!("filtered {removed} variants vs average ref depth({average})");
    }
    removed
}

struct VcfInput {
    header: Vec<String>,
    column_line: String,
}

/// Calls `on_record(ref_col, tumor_col, fields, line)` for each data line, in file order.
fn read_vcf(cfg: &TincConfig, mut on_record: impl FnMut(usize, usize, &[&str], &str) -> Result<()>) -> Result<VcfInput> {
    let bgzf_reader = bgzf::Reader::from_path(&cfg.input_vcf).with_context(|| format!("opening {}", cfg.input_vcf))?;
    let mut reader = BufReader::with_capacity(1 << 20, bgzf_reader);
    let mut header = Vec::new();
    let mut input: Option<VcfInput> = None;
    let mut cols = (0, 0);
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let text = line.trim_end_matches(['\n', '\r']);
        if text.starts_with("##") {
            header.push(text.to_string());
        } else if text.starts_with('#') {
            let names: Vec<&str> = text.split('\t').collect();
            let find = |id: &str| {
                names
                    .iter()
                    .skip(9)
                    .position(|c| *c == id)
                    .map(|i| i + 9)
                    .with_context(|| format!("sample {id} is not in {}", cfg.input_vcf))
            };
            cols = (find(&cfg.reference_id)?, find(&cfg.tumor_id)?);
            input = Some(VcfInput { header: std::mem::take(&mut header), column_line: text.to_string() });
        } else if !text.is_empty() {
            if input.is_none() {
                bail!("{} has no header", cfg.input_vcf);
            }
            let fields: Vec<&str> = text.split('\t').collect();
            on_record(cols.0, cols.1, &fields, text)?;
        }
    }
    input.with_context(|| format!("{} has no header", cfg.input_vcf))
}

pub fn generate_tinc_vcf_filename(input_vcf: &str) -> String {
    input_vcf.replace(".vcf.gz", ".tinc.vcf.gz")
}

pub struct TincAnalyser {
    cfg: TincConfig,
    /// the variants that will be rewritten, in file order
    variants: Vec<VariantData>,
    chromosomes: Vec<String>,
    /// index into `chromosomes` for each variant
    chromosome_index: Vec<u32>,
    /// Java only queries chr1-22,X,Y, except when the fit variants come from a file
    human_only: bool,
    pub level: f64,
}

impl TincAnalyser {
    pub fn new(cfg: TincConfig) -> Self {
        TincAnalyser { cfg, variants: Vec::new(), chromosomes: Vec::new(), chromosome_index: Vec::new(), human_only: true, level: 0.0 }
    }

    pub fn variants(&self) -> &[VariantData] {
        &self.variants
    }

    fn load_vcf_variants(&mut self) -> Result<()> {
        let human_only = self.human_only;
        let path = self.cfg.input_vcf.clone();
        let mut variants = Vec::new();
        let mut chromosomes: Vec<String> = Vec::new();
        let mut chromosome_index = Vec::new();
        read_vcf(&self.cfg, |ref_col, tumor_col, fields, _| {
            let chr = fields[0];
            if human_only && !is_human_chromosome(chr) {
                return Ok(());
            }
            if chromosomes.last().map(|c| c.as_str()) != Some(chr) {
                if let Some(last) = chromosomes.last() {
                    if human_only && chromosome_rank(last) >= chromosome_rank(chr) {
                        bail!("{path} is not sorted in chromosome order ({last} before {chr})");
                    }
                }
                chromosomes.push(chr.to_string());
            }
            variants.push(VariantData::from_fields(fields, ref_col, tumor_col)?);
            chromosome_index.push(chromosomes.len() as u32 - 1);
            Ok(())
        })?;
        self.variants = variants;
        self.chromosomes = chromosomes;
        self.chromosome_index = chromosome_index;
        Ok(())
    }

    fn annotate(&mut self) -> Result<()> {
        let mut wanted: Positions = HashMap::new();
        for (v, &c) in self.variants.iter().zip(&self.chromosome_index) {
            wanted.entry(self.chromosomes[c as usize].clone()).or_default().insert(v.position);
        }
        let (pon, gnomad) = load_resources(&self.cfg, &wanted, &self.chromosomes)?;
        let v38 = self.cfg.v38;
        for (v, &c) in self.variants.iter_mut().zip(&self.chromosome_index) {
            let chr = &self.chromosomes[c as usize];
            if let Some(list) = gnomad.get(chr).and_then(|m| m.get(&v.position)) {
                if let Some((_, _, freq)) = list.iter().find(|(r, a, _)| *r == v.ref_ && *a == v.alt) {
                    v.gnomad_frequency = Some(*freq);
                    v.pon_filtered = true;
                }
            }
            if let Some(list) = pon.get(chr).and_then(|m| m.get(&v.position)) {
                if let Some(e) = list.iter().find(|e| e.ref_ == v.ref_ && e.alt == v.alt) {
                    v.pon_sample_count = e.samples;
                    v.pon_max_read_count = e.max_reads;
                    v.pon_mean_read_count =
                        if e.samples > 0 { jround(e.total_reads as f64 / e.samples as f64) as i32 } else { 0 };
                    if pon_filter_on_tier(v38, v.tier, e.samples, e.max_reads) {
                        v.pon_filtered = true;
                    }
                }
            }
        }
        Ok(())
    }

    /// TincAnalyser.run: `rewrite` says whether the variants are needed for writing a VCF afterwards.
    pub fn run(&mut self, apply_recovery: bool, rewrite: bool) -> Result<()> {
        eprintln!("running TINC analyser");
        if let Some(file) = self.cfg.fit_variants_file.clone() {
            let fit = load_fit_variants(&file)?;
            eprintln!("loaded {} fit variants from file: {file}", fit.len());
            if rewrite {
                self.human_only = false;
                self.load_vcf_variants()?;
            }
            self.level = calculate(&fit.iter().collect::<Vec<_>>(), &populate_default_levels());
        } else {
            self.load_vcf_variants()?;
            self.annotate()?;
            let mut counts = [0usize; FILTER_REASONS.len()];
            let mut fitting: Vec<usize> = Vec::new();
            for (i, v) in self.variants.iter().enumerate() {
                let reason = check_filters(v);
                if reason == FilterReason::None {
                    fitting.push(i);
                } else {
                    counts[FILTER_REASONS.iter().position(|(r, _)| *r == reason).unwrap()] += 1;
                }
            }
            if fitting.len() >= TINC_MAX_FITTING_VARIANTS {
                downsample_fitting(&mut fitting);
            }
            counts[9] += filter_by_average_reference_depth(&self.variants, &mut fitting);
            let reasons: Vec<String> = FILTER_REASONS.iter().zip(counts).map(|((_, n), c)| format!("{n}={c}")).collect();
            eprintln!("variants({}) fitCount({}) reasons: {}", self.variants.len(), fitting.len(), reasons.join(","));
            if self.cfg.write_fit_variants {
                self.write_fit_variants(&fitting)?;
            }
            let fit: Vec<&VariantData> = fitting.iter().map(|&i| &self.variants[i]).collect();
            self.level = calculate(&fit, &populate_default_levels());
        }
        if apply_recovery && self.level > 0.0 {
            let recovered = recover_variants(&mut self.variants, self.level);
            if recovered > 0 {
                eprintln!("recovered {recovered} variants");
            }
        }
        Ok(())
    }

    fn write_fit_variants(&self, fitting: &[usize]) -> Result<()> {
        let path = self.cfg.input_vcf.replace(".vcf.gz", ".fit_variant.tsv");
        let mut out = std::io::BufWriter::new(File::create(&path).with_context(|| format!("creating {path}"))?);
        writeln!(out, "Chromosome\tPosition\tRef\tAlt\tRefDepth\tRefFrags\tTumorDepth\tTumorFrags")?;
        for &i in fitting {
            let v = &self.variants[i];
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.chromosomes[self.chromosome_index[i] as usize],
                v.position,
                v.ref_,
                v.alt,
                v.reference_depth,
                v.reference_alt_frags,
                v.tumor_depth,
                v.tumor_alt_frags
            )?;
        }
        Ok(())
    }

    /// TincVcfWriter.writeVcf: rewrite the input VCF (possibly in place) with annotations and recovered filters.
    pub fn write_vcf(&self, output_vcf: &str) -> Result<()> {
        let tmp = format!("{output_vcf}.tmp");
        let mut out = bgzf::Writer::from_path(&tmp).with_context(|| format!("creating {tmp}"))?;
        let mut body: Vec<u8> = Vec::with_capacity(1 << 20);
        let mut index = 0;
        let human_only = self.human_only;
        let input = read_vcf(&self.cfg, |_, _, fields, text| {
            if human_only && !is_human_chromosome(fields[0]) {
                return Ok(());
            }
            let v = self.variants.get(index).context("VCF changed while TINC was running")?;
            index += 1;
            if let Some(line) = output_line(v, fields, text) {
                body.extend_from_slice(line.as_bytes());
                body.push(b'\n');
            }
            Ok(())
        })?;
        // the header is only complete once the input has been read, and the input may be the file being replaced
        let mut lines = input.header;
        if self.level > 0.0 {
            lines.push(format!("##tincLevel={}", java_fixed(self.level, 3)));
        }
        lines.extend(TINC_HEADER_LINES.iter().map(|x| x.to_string()));
        for line in sorted_header_lines(lines) {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
        out.write_all(input.column_line.as_bytes())?;
        out.write_all(b"\n")?;
        out.write_all(&body)?;
        out.flush()?;
        drop(out);
        std::fs::rename(&tmp, output_vcf).with_context(|| format!("replacing {output_vcf}"))?;
        index_vcf(output_vcf)
    }
}

/// The rewritten record, or None if TINC drops it.
fn output_line(v: &VariantData, fields: &[&str], text: &str) -> Option<String> {
    if !v.is_passing() && !v.recovered() && filter_out_variant(v) {
        return None;
    }
    let has_pon = v.pon_sample_count > 0 || v.pon_max_read_count > 0;
    let filters_changed = matches!(&v.new_filters, Some(f) if f.len() != v.filters.len());
    if v.gnomad_frequency.is_none() && !has_pon && !filters_changed {
        return Some(text.to_string());
    }
    let mut info: Vec<(&str, String)> = if fields[7] == "." {
        Vec::new()
    } else {
        fields[7].split(';').map(|item| (item.split_once('=').map_or(item, |(k, _)| k), item.to_string())).collect()
    };
    if let Some(freq) = v.gnomad_frequency {
        info.push(("GND_FREQ", format!("GND_FREQ={}", vcf_double(freq))));
    }
    if has_pon {
        info.push(("PON_COUNT", format!("PON_COUNT={}", v.pon_sample_count)));
        info.push(("PON_MAX", format!("PON_MAX={}", v.pon_max_read_count)));
        info.push(("PON_READS", format!("PON_READS={}", v.pon_mean_read_count)));
    }
    let mut filter = fields[6].to_string();
    if filters_changed {
        if v.recovered() {
            info.push(("TINC_RECOVERED", "TINC_RECOVERED".to_string()));
            filter = "PASS".to_string();
        } else {
            let mut names = v.new_filters.clone().unwrap();
            names.sort();
            filter = names.join(";");
        }
    }
    info.sort_by(|a, b| a.0.cmp(b.0));
    let info: Vec<String> = info.into_iter().map(|(_, item)| item).collect();
    let mut out: Vec<&str> = fields.to_vec();
    let info = info.join(";");
    out[6] = &filter;
    out[7] = &info;
    Some(out.join("\t"))
}

fn load_fit_variants(path: &str) -> Result<Vec<VariantData>> {
    let mut reader = open_text(path)?;
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let cols: Vec<&str> = header.trim_end().split('\t').collect();
    let col = |name: &str| cols.iter().position(|c| *c == name).with_context(|| format!("{path} has no {name} column"));
    let (pos_i, ref_i, alt_i) = (col("Position")?, col("Ref")?, col("Alt")?);
    let (rd_i, rf_i, td_i, tf_i) = (col("RefDepth")?, col("RefFrags")?, col("TumorDepth")?, col("TumorFrags")?);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        out.push(VariantData::simple(
            f[pos_i].parse()?,
            f[ref_i],
            f[alt_i],
            f[rd_i].parse()?,
            f[rf_i].parse()?,
            f[td_i].parse()?,
            f[tf_i].parse()?,
        ));
    }
    Ok(out)
}

/// The standalone entry point, as TincAnalyser.main: `sage-rs tinc -input_vcf ... -tumor ... -reference ...`
pub fn run_standalone(args: &[String]) -> Result<()> {
    const TINC_FLAGS: [&str; 3] = ["write_tinc_vcf", "write_fit_variants", "log_debug"];
    let mut vals: HashMap<String, String> = HashMap::new();
    let mut i = 0;
    while i < args.len() {
        let Some(key) = args[i].strip_prefix('-') else { bail!("unexpected argument: {}", args[i]) };
        if TINC_FLAGS.contains(&key) {
            vals.insert(key.to_string(), "true".to_string());
            i += 1;
        } else {
            let v = args.get(i + 1).with_context(|| format!("missing value for -{key}"))?;
            vals.insert(key.to_string(), v.clone());
            i += 2;
        }
    }
    let req = |k: &str| vals.get(k).cloned().with_context(|| format!("missing required option -{k}"));
    let fit_variants_file = vals.get("fit_variants_file").cloned();
    let cfg = TincConfig {
        v38: vals.get("ref_genome_version").map(|v| v.contains("38")).unwrap_or(false),
        input_vcf: req("input_vcf")?,
        tumor_id: req("tumor")?,
        reference_id: req("reference")?,
        pon_file: vals.get("pon_file").cloned(),
        gnomad_dir: vals.get("gnomad_freq_dir").cloned(),
        gnomad_file: vals.get("gnomad_freq_file").cloned(),
        threads: vals.get("threads").map(|v| v.parse()).transpose()?.unwrap_or(1),
        write_fit_variants: fit_variants_file.is_none() && vals.contains_key("write_fit_variants"),
        fit_variants_file,
    };
    let rewrite = vals.contains_key("write_tinc_vcf");
    let output = generate_tinc_vcf_filename(&cfg.input_vcf);
    let mut analyser = TincAnalyser::new(cfg);
    analyser.run(true, rewrite)?;
    if rewrite {
        analyser.write_vcf(&output)?;
    }
    eprintln!("TINC complete");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // TincTest.createTincVariantFromContext
    fn variant_from_context(ref_depth: i32, ref_alt: i32, tumor_depth: i32, tumor_alt: i32, mqf: f64) -> VariantData {
        let geno = |depth: i32, ad0: i32, alt: i32| {
            let counts = format!("{alt},0,0,0,{},{depth}", depth - alt);
            format!("./.:{ad0},{alt}:30,30:{depth}:{counts}:{counts}")
        };
        let ref_geno = geno(ref_depth, 0, ref_alt);
        let tumor_geno = geno(tumor_depth, ref_alt, tumor_alt);
        let info = format!("MQF={mqf}");
        let fields = ["chr1", "100", ".", "A", "C", ".", ".", &info, "GT:AD:ARCBQ:DP:RC_CNT:RC_QUAL", &ref_geno, &tumor_geno];
        VariantData::from_fields(&fields, 9, 10).unwrap()
    }

    #[test]
    fn tinc_calculations() {
        // TincTest.testTincCalculations
        let mut variants: Vec<VariantData> = Vec::new();
        assert_eq!(calculate(&variants.iter().collect::<Vec<_>>(), &populate_default_levels()), 0.0);
        for _ in 0..20 {
            for j in 0..=TINC_GERMLINE_MAX_AD {
                variants.push(VariantData::simple(100, "A", "C", 30, j, 100, 30));
            }
        }
        let level = calculate(&variants.iter().collect::<Vec<_>>(), &populate_default_levels());
        assert!((level - 0.75).abs() < 0.001, "{level}");
    }

    #[test]
    fn tinc_recovery() {
        // TincTest.testTincRecovery
        let mut v = variant_from_context(30, 5, 100, 30, 30.0);
        v.filters = vec![F_REL_QUAL.to_string(), F_VAF.to_string(), F_ALT_SUPPORT.to_string()];
        let mut variants = vec![v];
        recover_variants(&mut variants, 0.2);
        assert!(variants[0].new_filters.as_ref().unwrap().is_empty());
        recover_variants(&mut variants, 0.001);
        let filters = variants[0].new_filters.clone().unwrap();
        assert_eq!(filters.len(), 2);
        assert!(filters.iter().any(|x| x == F_REL_QUAL));
        assert!(filters.iter().any(|x| x == F_VAF));
    }

    #[test]
    fn tinc_filters() {
        // TincTest.testTincFilters
        let mqf = 30.0;
        assert_eq!(check_filters(&variant_from_context(30, 10, 100, 50, mqf)), FilterReason::GermlineAf);

        let mut v = variant_from_context(30, 1, 100, 50, mqf);
        v.pon_filtered = true;
        assert_eq!(check_filters(&v), FilterReason::Pon);

        let mut v = variant_from_context(30, 1, 100, 50, mqf);
        v.gnomad_frequency = Some(0.01);
        assert_eq!(check_filters(&v), FilterReason::Gnomad);

        let mut v = variant_from_context(30, 1, 100, 50, mqf);
        v.ref_abq = 25;
        assert_eq!(check_filters(&v), FilterReason::LowAbq);

        assert_eq!(check_filters(&variant_from_context(30, 1, 100, 50, 10.0)), FilterReason::Mqf);
    }

    #[test]
    fn default_levels() {
        let levels = populate_default_levels();
        assert_eq!(levels[0], 0.0);
        assert_eq!(levels[1], 0.005);
        // the last level is 0.95, not 1.0: the accumulated value overshoots 1.0, as in Java
        assert_eq!(levels.len(), 35);
        assert_eq!(*levels.last().unwrap(), 0.9500000000000003);
        assert_eq!(levels[23], 0.39999999999999997);
    }

    #[test]
    fn downsampling_keeps_every_nth() {
        let mut fitting: Vec<usize> = (0..45_000).collect();
        downsample_fitting(&mut fitting);
        assert_eq!(fitting.len(), 15_000);
        assert_eq!(&fitting[..3], &[2, 5, 8]);
    }
}
