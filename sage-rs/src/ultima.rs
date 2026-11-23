//! Ultima sequencing support: the seqtech/Ultima* classes, Homopolymer, MergedHomopolymers, RefMask
//! and the parts of hmf-common UltimaBamUtils that SAGE uses.
//!
//! Java signals several conditions here with exceptions (array index out of bounds, a null reference
//! genome result, IllegalArgumentException). Those are `Err(())` (`Ex<T>`); each caller does what the
//! Java catch block that would have received the exception does.

use crate::constants::*;
use crate::evidence::{read_index_from_position, ReadContextCounter};
use crate::read_context::{ReadCigarInfo, VariantReadContext};
use crate::refseq::RefSequence;
use crate::types::*;
use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::sync::atomic::{AtomicI32, AtomicI8, Ordering};
use std::sync::OnceLock;

/// The result of code that can throw in Java.
pub type Ex<T> = std::result::Result<T, ()>;

// UltimaBamUtils
pub const ULTIMA_MAX_QUAL: i8 = 35;
pub const ULTIMA_MAX_HP_LEN: i32 = 20;
pub const ULTIMA_INVALID_QUAL: i8 = -1;
pub const HALF_PHRED_SCORE_SCALING: i8 = 3;
const CYCLE_BASES: [u8; 4] = *b"TGCA";

// UltimaUtils
const MAX_HOMOPOLYMER: i32 = ULTIMA_MAX_HP_LEN;
/// INVALID_BASE (-1 as a Java byte)
const INVALID_BASE: u8 = 0xFF;
const TP_ZERO_BASE_QUAL: i8 = 0;
const T0_EXPECTED_QUAL_THRESHOLD_SIMPLE: f64 = 18.0;
const T0_EXPECTED_QUAL_THRESHOLD_COMPLEX: f64 = 24.0;
const MSI_SAMPLE_REPEAT_LENGTH_MAX: i32 = 8;
pub const GERMLINE_VAF_INDEL_REPEAT_THRESHOLD: f64 = 2.5;
pub const GERMLINE_VAF_INDEL_REPEAT_MIN_THRESHOLD_FACTOR: f64 = 4.0;
pub const GERMLINE_VAF_REL_QUAL_RATIO_THRESHOLD: f64 = 0.17;
const INVALID_BASE_QUAL: i8 = -1;

// FilterConfig
pub const ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN: i32 = 8;
/// FilterConfig.ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD (-ultima_cand_min_high_bq), a static in Java too
static ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD: AtomicI32 = AtomicI32::new(1);

pub fn set_candidate_min_high_bq_threshold(threshold: i32) {
    ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD.store(threshold, Ordering::Relaxed);
}
pub fn candidate_min_high_bq_threshold() -> i32 {
    ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD.load(Ordering::Relaxed)
}

fn at(bases: &[u8], index: i32) -> Ex<u8> {
    if index < 0 || index >= bases.len() as i32 {
        Err(())
    } else {
        Ok(bases[index as usize])
    }
}

/// Nucleotides.swapDnaBase
fn swap_dna_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        o => o,
    }
}

fn is_cycle_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

// --- UltimaBamUtils ---

/// UltimaBamUtils.extractLowQualIndices over the UQ tag ("count=index,start-end,..."); an absent tag
/// gives an empty list
pub fn extract_low_qual_indices(read: &Read) -> Vec<i32> {
    let Some(tag) = &read.low_qual_tag else { return Vec::new() };
    let Some((_, items)) = tag.split_once('=') else { return Vec::new() };
    let mut values = Vec::new();
    for item in items.split(',') {
        if let Some((start, end)) = item.split_once('-') {
            if let (Ok(start), Ok(end)) = (start.parse::<i32>(), end.parse::<i32>()) {
                values.extend(start..=end);
            }
        } else if let Ok(v) = item.parse::<i32>() {
            values.push(v);
        }
    }
    values
}

/// UltimaBamUtils.cycleCount
pub fn cycle_count(first_base: u8, last_base: u8, inner_base: u8) -> i32 {
    if !is_cycle_base(first_base) || !is_cycle_base(last_base) || !is_cycle_base(inner_base) {
        return -1;
    }
    let bases = [first_base, inner_base, last_base];
    let mut base_index = 0;
    let mut cycle_index = 0;
    let mut cycle_count = 1;
    while base_index < bases.len() {
        if CYCLE_BASES[cycle_index] == bases[base_index] {
            base_index += 1;
        } else {
            cycle_index += 1;
        }
        if cycle_index == CYCLE_BASES.len() {
            cycle_index = 0;
            cycle_count += 1;
        }
    }
    cycle_count
}

/// UltimaBamUtils.isBaseInCycle
pub fn is_base_in_cycle(start_base: u8, end_base: u8, test_base: u8) -> bool {
    if start_base == end_base || start_base == test_base {
        return false;
    }
    if !is_cycle_base(start_base) || !is_cycle_base(end_base) {
        return false;
    }
    let mut started = false;
    let mut index = 0;
    // Java loops until it returns; a test base outside TGCA ends at the start or end base
    loop {
        if !started {
            if CYCLE_BASES[index] == start_base {
                started = true;
            }
        } else if CYCLE_BASES[index] == end_base || CYCLE_BASES[index] == start_base {
            return false;
        } else if CYCLE_BASES[index] == test_base {
            return true;
        }
        index = if index == CYCLE_BASES.len() - 1 { 0 } else { index + 1 };
    }
}

// --- UltimaQualRecalibration ---

const RECALIBRATED_QUAL_MAX_HP_LENGTH: i32 = 5;
const MAX_TP_T0_RECALIBRATION: i8 = 50;
const T0_OUT_OF_CYCLE: &str = "OUT_OF_CYCLE";
const DEFAULT_RECALIBRATION: &str = include_str!("../resources/ultima_qual_recalibration.tsv");

/// UltimaQualRecalibration: recalibrated TP quals by (homopolymer length, base, tp is zero) and T0 quals
/// by (bases either side, base)
pub struct UltimaQualRecalibration {
    tp_qual_map: HashMap<(i32, u8, bool), i8>,
    t0_qual_map: HashMap<([u8; 2], u8), i8>,
    out_of_cycle_t0_qual: i8,
}

impl UltimaQualRecalibration {
    fn empty() -> Self {
        UltimaQualRecalibration {
            tp_qual_map: HashMap::new(),
            t0_qual_map: HashMap::new(),
            out_of_cycle_t0_qual: MAX_TP_T0_RECALIBRATION,
        }
    }

    /// loadRecalibrationData. Keys that formKey cannot produce are never looked up, so they are dropped.
    pub fn from_lines<'a>(lines: impl Iterator<Item = &'a str>) -> Result<Self> {
        let mut lines = lines;
        let header = lines.next().context("empty Ultima base-qual recalibration file")?;
        let col = |name: &str| header.split('\t').position(|c| c == name).with_context(|| format!("missing column {name}"));
        let (type_index, key_index, rq_index) = (col("Type")?, col("Key")?, col("RecalibratedQual")?);
        let mut recal = Self::empty();
        for line in lines {
            let values: Vec<&str> = line.split('\t').collect();
            if values.len() <= type_index.max(key_index).max(rq_index) {
                bail!("invalid Ultima base-qual recalibration line: {line}");
            }
            let key = values[key_index];
            let qual = values[rq_index].parse::<i8>().with_context(|| format!("invalid recalibrated qual: {line}"))?;
            let qual = qual.min(MAX_TP_T0_RECALIBRATION);
            match values[type_index] {
                "TP_RECAL" => {
                    // 1_A_FALSE
                    let parts: Vec<&str> = key.split('_').collect();
                    if let [length, base, tp_is_zero] = parts[..] {
                        if let (Ok(length), 1, true) =
                            (length.parse::<i32>(), base.len(), tp_is_zero == "TRUE" || tp_is_zero == "FALSE")
                        {
                            recal.tp_qual_map.insert((length, base.as_bytes()[0], tp_is_zero == "TRUE"), qual);
                        }
                    }
                }
                "T0_RECAL" => {
                    // AC_T
                    if key == T0_OUT_OF_CYCLE {
                        recal.out_of_cycle_t0_qual = qual;
                    } else if let Some((tnc, base)) = key.split_once('_') {
                        if tnc.len() == 2 && base.len() == 1 {
                            let tnc = tnc.as_bytes();
                            recal.t0_qual_map.insert(([tnc[0], tnc[1]], base.as_bytes()[0]), qual);
                        }
                    }
                }
                other => bail!("invalid Ultima base-qual recalibration type: {other}"),
            }
        }
        Ok(recal)
    }

    pub fn out_of_cycle_t0_qual(&self) -> i8 {
        self.out_of_cycle_t0_qual
    }

    /// calcTpRecalibratedQual
    pub fn calc_tp_recalibrated_qual(
        &self, read_qual: i8, homopolymer_length: i32, base: u8, tp_is_zero: bool, is_reverse: bool,
    ) -> i8 {
        if read_qual < max_raw_qual() || homopolymer_length > RECALIBRATED_QUAL_MAX_HP_LENGTH {
            return read_qual;
        }
        let bqr_qual = self.get_tp_recalibrated_qual(homopolymer_length, base, tp_is_zero, is_reverse);
        if bqr_qual != INVALID_BASE_QUAL {
            bqr_qual
        } else {
            read_qual
        }
    }

    /// getTpRecalibratedQual
    pub fn get_tp_recalibrated_qual(&self, homopolymer_length: i32, base: u8, tp_is_zero: bool, is_reverse: bool) -> i8 {
        let base_to_lookup = if is_reverse { swap_dna_base(base) } else { base };
        self.tp_qual_map.get(&(homopolymer_length, base_to_lookup, tp_is_zero)).copied().unwrap_or(INVALID_BASE_QUAL)
    }

    /// calcT0RecalibratedQual
    pub fn calc_t0_recalibrated_qual(&self, read_qual: i8, read: &Read, var_read_index: i32) -> i8 {
        if read_qual < max_raw_qual() {
            return read_qual;
        }
        let bqr_qual = self.get_t0_recalibrated_qual(read, var_read_index);
        if bqr_qual != INVALID_BASE_QUAL {
            bqr_qual
        } else {
            read_qual
        }
    }

    /// getT0RecalibratedQual(record, varReadIndex)
    pub fn get_t0_recalibrated_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        if var_read_index <= 0 || var_read_index >= read.len() - 1 {
            return ULTIMA_INVALID_QUAL;
        }
        let i = var_read_index as usize;
        let is_reverse = read.negative_strand();
        let variant_base = read.bases[i];
        let base_to_lookup = if is_reverse { swap_dna_base(variant_base) } else { variant_base };
        let tnc = [read.bases[i - 1], read.bases[i + 1]];
        let tnc_to_lookup = if is_reverse { [swap_dna_base(tnc[1]), swap_dna_base(tnc[0])] } else { tnc };
        self.t0_qual_map.get(&(tnc_to_lookup, base_to_lookup)).copied().unwrap_or(self.out_of_cycle_t0_qual)
    }
}

/// UltimaUtils.BQR_CACHE. Java's is a mutable static; here the table is set once and the maximum raw
/// qual, which SageApplication sets later from the sample BQR files, is kept beside it.
static BQR_CACHE: OnceLock<UltimaQualRecalibration> = OnceLock::new();
static MAX_RAW_QUAL: AtomicI8 = AtomicI8::new(ULTIMA_MAX_QUAL);

pub fn bqr_cache() -> &'static UltimaQualRecalibration {
    BQR_CACHE.get_or_init(UltimaQualRecalibration::empty)
}

/// UltimaUtils.loadBqrCache: from the file, or the table shipped with SAGE. The first table loaded is
/// kept for the life of the process.
pub fn load_bqr_cache(filename: Option<&str>) -> Result<()> {
    let recal = match filename {
        Some(f) => {
            let text = std::fs::read_to_string(f)
                .with_context(|| format!("failed to read Ultima base-qual recalibration file({f})"))?;
            UltimaQualRecalibration::from_lines(text.lines())?
        }
        None => UltimaQualRecalibration::from_lines(DEFAULT_RECALIBRATION.lines())?,
    };
    let _ = BQR_CACHE.set(recal);
    Ok(())
}

/// UltimaUtils.setMaxRawQual
pub fn set_max_raw_qual(qual: i8) {
    MAX_RAW_QUAL.store(qual, Ordering::Relaxed);
}
/// UltimaUtils.maxRawQual
pub fn max_raw_qual() -> i8 {
    MAX_RAW_QUAL.load(Ordering::Relaxed)
}

// --- UltimaUtils ---

/// UltimaUtils.safeQualLookup
fn safe_qual_lookup(quals: &[i8], index: i32) -> i8 {
    if index < 0 || index >= quals.len() as i32 {
        ULTIMA_INVALID_QUAL
    } else {
        quals[index as usize]
    }
}

/// UltimaUtils.calcTpBaseQual
pub fn calc_tp_base_qual(read: &Read, index_start: i32, index_end: i32, tp_search_value: i32) -> i8 {
    if index_start < 0 || index_start >= read.len() {
        return ULTIMA_INVALID_QUAL;
    }
    // extractTpValues: Java fails on a read without the tag
    assert!(!read.tp.is_empty(), "read {} has no Ultima tp tag", read.name);
    let tp_values = &read.tp;
    let mut qual_value1: i8 = -1;
    let mut qual_value2: i8 = -1;
    let mut i = index_start;
    while i <= index_end.min(tp_values.len() as i32 - 1) {
        if tp_values[i as usize] as i32 == tp_search_value {
            if qual_value1 < 0 {
                qual_value1 = read.quals[i as usize] as i8;
            } else {
                qual_value2 = read.quals[i as usize] as i8;
                break;
            }
        }
        i += 1;
    }
    let homopolymer_length = index_end - index_start + 1;
    let homopolymer_base = read.bases[index_start as usize];
    let is_reverse = read.negative_strand();
    let bqr = bqr_cache();

    qual_value1 = bqr.calc_tp_recalibrated_qual(qual_value1, homopolymer_length, homopolymer_base, false, is_reverse);
    qual_value2 = bqr.calc_tp_recalibrated_qual(qual_value2, homopolymer_length, homopolymer_base, false, is_reverse);

    if qual_value1 < 0 {
        // if bases aren't found, use the middle base as: min(40, Q + 6 * abs(required_tp - T))
        let middle_index = index_start + (index_end - index_start) / 2;
        if middle_index < 0 || middle_index >= read.quals.len() as i32 {
            return ULTIMA_INVALID_QUAL;
        }
        qual_value1 = read.quals[middle_index as usize] as i8;
        let tp_value = tp_values[middle_index as usize];
        if tp_value == TP_ZERO_BASE_QUAL {
            return bqr.calc_tp_recalibrated_qual(qual_value1, homopolymer_length, homopolymer_base, true, is_reverse);
        }
        qual_value1 = bqr.calc_tp_recalibrated_qual(qual_value1, homopolymer_length, homopolymer_base, false, is_reverse);
        let bqr_value =
            bqr.calc_tp_recalibrated_qual(max_raw_qual(), homopolymer_length, homopolymer_base, false, is_reverse);
        let adjusted =
            qual_value1 as i32 + HALF_PHRED_SCORE_SCALING as i32 * 2 * (tp_search_value - tp_value as i32).abs();
        return (bqr_value as i32).min(adjusted) as i8;
    }
    if qual_value2 < 0 {
        return qual_value1;
    }
    (qual_value1 as i32 - HALF_PHRED_SCORE_SCALING as i32) as i8
}

/// UltimaUtils.findHomopolymerLength
fn find_homopolymer_length(ref_bases: &[u8], compare_base: u8, start_index: i32, search_up: bool) -> i32 {
    let mut repeat_count = 0;
    let mut i = start_index;
    while i >= 0 && i < ref_bases.len() as i32 {
        if ref_bases[i as usize] != compare_base {
            return repeat_count;
        }
        repeat_count += 1;
        i += if search_up { 1 } else { -1 };
    }
    repeat_count
}

/// UltimaUtils.isCleanSnv
#[allow(dead_code)] // used by the ported tests
pub fn is_clean_snv(ctx: &VariantReadContext) -> bool {
    if !ctx.variant.is_snv() {
        return false;
    }
    if ctx.core_length() != ctx.ref_bases.len() as i32 {
        return false;
    }
    for i in ctx.core_index_start..=ctx.core_index_end {
        if i == ctx.var_index {
            continue;
        }
        if ctx.read_bases[i as usize] != ctx.ref_bases[(i - ctx.core_index_start) as usize] {
            return false;
        }
    }
    true
}

/// UltimaUtils.isMsiIndelOfType
#[allow(dead_code)] // used by the ported tests
pub fn is_msi_indel_of_type(variant: &SimpleVariant, units: &[&str]) -> bool {
    if !variant.is_indel() {
        return false;
    }
    let indel_bases: &[u8] = if variant.is_insert() { &variant.alt[1..] } else { &variant.ref_[1..] };
    for unit in units {
        let num_repeats = indel_bases.len() / unit.len();
        if indel_bases == unit.repeat(num_repeats).as_bytes() {
            return true;
        }
    }
    false
}

#[allow(dead_code)] // used by the ported tests
fn is_dna_base_homopolymer(bases: &[u8], start_index: i32, end_index: i32) -> bool {
    if !is_cycle_base(bases[start_index as usize]) {
        return false;
    }
    ((start_index + 1)..=end_index).all(|i| bases[(i - 1) as usize] == bases[i as usize])
}

/// UltimaUtils.isAdjacentToLongHomopolymer
#[allow(dead_code)] // used by the ported tests
pub fn is_adjacent_to_long_homopolymer(read: &Read, var_index_in_read: i32, long_length: i32) -> bool {
    let read_bases = &read.bases;
    if var_index_in_read + long_length < read.len()
        && is_dna_base_homopolymer(read_bases, var_index_in_read + 1, var_index_in_read + long_length)
    {
        return true;
    }
    var_index_in_read >= long_length
        && is_dna_base_homopolymer(read_bases, var_index_in_read - long_length, var_index_in_read - 1)
}

/// UltimaUtils.hasLongHomopolymerInContext
pub fn has_long_homopolymer_in_context(variant: &SimpleVariant, read: &Read, var_read_index: i32) -> Ex<bool> {
    let standard_length = MIN_CORE_DISTANCE + DEFAULT_FLANK_LENGTH;
    let core_start_index = (var_read_index - standard_length + if variant.is_indel() { 1 } else { 0 }).max(0);
    let core_end_index = (var_read_index + (variant.alt_len() - 1) + standard_length).min(read.len() - 1);
    has_long_homopolymer_in_range(core_start_index, core_end_index, ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN, read)
}

/// UltimaUtils.hasLongHomopolymerInRange
pub fn has_long_homopolymer_in_range(start_index: i32, end_index: i32, long_repeat_length: i32, read: &Read) -> Ex<bool> {
    let mut last_base = at(&read.bases, start_index)?;
    let mut repeat_count = 1;
    for i in start_index + 1..=end_index {
        let next_base = at(&read.bases, i)?;
        if next_base == last_base {
            repeat_count += 1;
            if repeat_count >= long_repeat_length {
                return Ok(true);
            }
        } else {
            last_base = next_base;
            repeat_count = 1;
        }
    }
    Ok(false)
}

/// UltimaUtils.lowQualInReadCore
pub fn low_qual_in_read_core(variant: &SimpleVariant, variant_read_index: i32, low_qual_indices: Option<&[i32]>) -> bool {
    let Some(low_qual_indices) = low_qual_indices else { return false };
    if low_qual_indices.is_empty() {
        return false;
    }
    let core_index_start = variant_read_index - MIN_CORE_DISTANCE + if variant.is_indel() { 1 } else { 0 };
    let core_index_end = variant_read_index + variant.alt_len() - 1 + MIN_CORE_DISTANCE;
    (core_index_start..=core_index_end).any(|index| low_qual_indices.contains(&index))
}

/// UltimaUtils.lowQualInReadContextCore
pub fn low_qual_in_read_context_core(variant_read_index: i32, low_qual_indices: &[i32], ctx: &VariantReadContext) -> bool {
    if low_qual_indices.is_empty() {
        return false;
    }
    // exclude from checking if the variant has a long homopolymer
    if let Some(r) = &ctx.max_repeat {
        if r.repeat_length() == 1 && r.count >= ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN {
            return false;
        }
    }
    // check for a sufficiently long single-base repeat anywhere in the read bases
    let mut last_base = ctx.read_bases[0];
    let mut repeat_count = 1;
    for &next_base in &ctx.read_bases[1..] {
        if next_base == last_base {
            repeat_count += 1;
            if repeat_count >= ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN {
                return false;
            }
        } else {
            last_base = next_base;
            repeat_count = 1;
        }
    }
    let core_index_start = variant_read_index - ctx.left_core_length();
    let core_index_end = variant_read_index + ctx.right_core_length();
    (core_index_start..=core_index_end).any(|index| low_qual_indices.contains(&index))
}

// --- Homopolymer, RefMask, MergedHomopolymers ---

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Homopolymer {
    pub base: u8,
    pub length: i32,
}

impl Homopolymer {
    pub fn new(base: u8, length: i32) -> Self {
        Homopolymer { base, length }
    }
    #[allow(dead_code)] // used by the ported tests
    pub fn expand(&self) -> String {
        (self.base as char).to_string().repeat(self.length as usize)
    }
}

/// Homopolymer.getHomopolymers
pub fn get_homopolymers(bases: &[u8], start_index: i32, end_index: i32) -> Vec<Homopolymer> {
    let mut homopolymers = Vec::new();
    if start_index < 0 || end_index >= bases.len() as i32 || end_index < start_index {
        return homopolymers;
    }
    let mut current_base = bases[start_index as usize];
    let mut current_length = 1;
    for i in start_index + 1..=end_index {
        let base = bases[i as usize];
        if base == current_base {
            current_length += 1;
            continue;
        }
        homopolymers.push(Homopolymer::new(current_base, current_length));
        current_base = base;
        current_length = 1;
    }
    homopolymers.push(Homopolymer::new(current_base, current_length));
    homopolymers
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RefMask {
    pub pos_start: i32,
    pub pos_end: i32,
    pub base_mask: u8,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MergedHomopolymers {
    pub ref_homopolymers: Vec<Homopolymer>,
    pub read_homopolymers: Vec<Homopolymer>,
    pub variant_in_merged_homopolymers: bool,
    pub ref_masks: Vec<RefMask>,
}

// --- qual models ---

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UltimaModelType {
    HomopolymerAdjustment,
    HomopolymerDeletion,
    HomopolymerTransition,
    Snv,
    BaseShift,
    Other,
}

/// UltimaHomopolymerAdjustment
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UltimaHomopolymerAdjustment {
    pub hp_start_index: i32,
    pub hp_end_index: i32,
    /// bases to adjust to get back to the ref
    pub ref_adjust_count: i32,
}

impl UltimaHomopolymerAdjustment {
    pub fn new(hp_start_index: i32, hp_end_index: i32, ref_adjust_count: i32) -> Self {
        UltimaHomopolymerAdjustment { hp_start_index, hp_end_index, ref_adjust_count }
    }
    pub fn calculate_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        calc_tp_base_qual(read, var_read_index + self.hp_start_index, var_read_index + self.hp_end_index, self.ref_adjust_count)
    }
}

/// UltimaHomopolymerDeletion
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UltimaHomopolymerDeletion {
    // the index and read base straddling the variant's deleted HP
    straddle_index_start: i32,
    straddle_index_end: i32,
    straddle_base_start: u8,
    straddle_base_end: u8,
    // whether the delete is in-cycle when read from each direction
    in_cycle_pos_strand: bool,
    in_cycle_neg_strand: bool,
}

impl UltimaHomopolymerDeletion {
    /// fromDelete
    pub fn from_delete(deleted_base: u8, ref_base: u8, post_delete_ref_base: u8) -> Self {
        UltimaHomopolymerDeletion {
            straddle_index_start: 0,
            straddle_index_end: 1,
            straddle_base_start: ref_base,
            straddle_base_end: post_delete_ref_base,
            in_cycle_pos_strand: is_base_in_cycle(ref_base, post_delete_ref_base, deleted_base),
            in_cycle_neg_strand: is_base_in_cycle(
                swap_dna_base(post_delete_ref_base),
                swap_dna_base(ref_base),
                swap_dna_base(deleted_base),
            ),
        }
    }

    /// fromMnv: for an SNV the deleted ref is the ref of the SNV and the new alt is the alt
    pub fn from_mnv(straddle_index: i32, deleted_ref_base: u8, new_alt_base: u8, read_base_before: u8, read_base_after: u8) -> Self {
        let pos_ref = cycle_count(read_base_before, read_base_after, deleted_ref_base);
        let pos_alt = cycle_count(read_base_before, read_base_after, new_alt_base);
        let (before, after) = (swap_dna_base(read_base_before), swap_dna_base(read_base_after));
        let neg_ref = cycle_count(after, before, swap_dna_base(deleted_ref_base));
        let neg_alt = cycle_count(after, before, swap_dna_base(new_alt_base));
        UltimaHomopolymerDeletion {
            straddle_index_start: straddle_index,
            straddle_index_end: straddle_index,
            straddle_base_start: read_base_before,
            straddle_base_end: read_base_after,
            in_cycle_pos_strand: pos_ref == pos_alt,
            in_cycle_neg_strand: neg_ref == neg_alt,
        }
    }

    pub fn calculate_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        let bqr = bqr_cache();
        let in_cycle = if read.negative_strand() { self.in_cycle_neg_strand } else { self.in_cycle_pos_strand };
        if !in_cycle {
            return bqr.out_of_cycle_t0_qual();
        }
        // extractT0Values: Java fails on a read without the tag
        assert!(!read.t0.is_empty(), "read {} has no Ultima t0 tag", read.name);
        let t0_values = &read.t0;
        let index1 = var_read_index + self.straddle_index_start;
        let qual1 = bqr.calc_t0_recalibrated_qual(safe_qual_lookup(t0_values, index1), read, index1);
        if qual1 == ULTIMA_INVALID_QUAL {
            return ULTIMA_INVALID_QUAL;
        }
        let qual2 = if self.straddle_index_end == self.straddle_index_start {
            qual1
        } else {
            let index2 = var_read_index + self.straddle_index_end;
            let qual2 = bqr.calc_t0_recalibrated_qual(safe_qual_lookup(t0_values, index2), read, index2);
            if qual2 == ULTIMA_INVALID_QUAL {
                return ULTIMA_INVALID_QUAL;
            }
            qual2
        };
        qual1.max(qual2)
    }
}

/// UltimaHomopolymerTransitionDeletion
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UltimaHomopolymerTransitionDeletion {
    lower_hp_start_index: i32,
    lower_hp_end_index: i32,
    lower_ref_adjust_count: i32,
    upper_hp_start_index: i32,
    upper_hp_end_index: i32,
    upper_ref_adjust_count: i32,
}

impl UltimaHomopolymerTransitionDeletion {
    fn new(variant: &SimpleVariant, homopolymer_transition_index: i32, lower_hp_length: i32, upper_hp_length: i32) -> Self {
        let del_length = variant.indel_length_abs();
        let lower_del_length = homopolymer_transition_index - 1;
        let upper_del_length = del_length - lower_del_length;
        // the ref base always matches the lower HP base, and the ref + 1 base the upper HP base
        let lower_hp_end_index = 0;
        let upper_hp_start_index = lower_hp_end_index + 1;
        UltimaHomopolymerTransitionDeletion {
            lower_hp_start_index: lower_hp_end_index - lower_hp_length + 1 + lower_del_length,
            lower_hp_end_index,
            lower_ref_adjust_count: lower_del_length,
            upper_hp_start_index,
            upper_hp_end_index: upper_hp_start_index + upper_hp_length - 1 - upper_del_length,
            upper_ref_adjust_count: upper_del_length,
        }
    }

    /// calculateQual: reading the homopolymer bases can run off the read, which throws in Java
    pub fn calculate_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        let lower_qual = calc_tp_base_qual(
            read, var_read_index + self.lower_hp_start_index, var_read_index + self.lower_hp_end_index,
            self.lower_ref_adjust_count,
        );
        if lower_qual == ULTIMA_INVALID_QUAL {
            return ULTIMA_INVALID_QUAL;
        }
        let upper_qual = calc_tp_base_qual(
            read, var_read_index + self.upper_hp_start_index, var_read_index + self.upper_hp_end_index,
            self.upper_ref_adjust_count,
        );
        if upper_qual == ULTIMA_INVALID_QUAL {
            return ULTIMA_INVALID_QUAL;
        }
        // cap at the higher of the BQR quals
        let bqr = bqr_cache();
        let is_reverse = read.negative_strand();
        let lower_hp_length = self.lower_hp_end_index - self.lower_hp_start_index + 1;
        let lower_hp_base = read.bases[(var_read_index + self.lower_hp_start_index) as usize];
        let lower_bqr_qual = bqr.get_tp_recalibrated_qual(lower_hp_length, lower_hp_base, false, is_reverse) as i32;
        let upper_hp_length = self.upper_hp_end_index - self.upper_hp_start_index + 1;
        let upper_hp_base = read.bases[(var_read_index + self.upper_hp_start_index) as usize];
        let upper_bqr_qual = bqr.get_tp_recalibrated_qual(upper_hp_length, upper_hp_base, false, is_reverse) as i32;
        let max_bqr_qual = lower_bqr_qual.max(upper_bqr_qual);
        (lower_qual as i32 + upper_qual as i32).min(max_bqr_qual) as i8
    }
}

/// The left and right components shared by UltimaSnv and UltimaBaseShift.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Adjustments {
    pub left_adjust: Option<UltimaHomopolymerAdjustment>,
    pub right_adjust: Option<UltimaHomopolymerAdjustment>,
    pub left_deletion: Option<UltimaHomopolymerDeletion>,
    pub right_deletion: Option<UltimaHomopolymerDeletion>,
}

impl Adjustments {
    /// UltimaBaseShift.canCompute / UltimaSnv.canCompute
    pub fn can_compute(&self) -> bool {
        (self.left_adjust.is_some() || self.left_deletion.is_some())
            && (self.right_adjust.is_some() || self.right_deletion.is_some())
    }

    fn calculate_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        if !self.can_compute() {
            return ULTIMA_INVALID_QUAL;
        }
        let left_qual = match (&self.left_adjust, &self.left_deletion) {
            (Some(adjust), _) => adjust.calculate_qual(read, var_read_index),
            (None, Some(deletion)) => deletion.calculate_qual(read, var_read_index),
            (None, None) => ULTIMA_INVALID_QUAL,
        };
        if left_qual == ULTIMA_INVALID_QUAL {
            return ULTIMA_INVALID_QUAL;
        }
        let right_qual = match (&self.right_adjust, &self.right_deletion) {
            (Some(adjust), _) => adjust.calculate_qual(read, var_read_index),
            (None, Some(deletion)) => deletion.calculate_qual(read, var_read_index),
            (None, None) => ULTIMA_INVALID_QUAL,
        };
        if right_qual == ULTIMA_INVALID_QUAL {
            return ULTIMA_INVALID_QUAL;
        }
        left_qual.max(right_qual)
    }
}

/// Reference bases by position on the chromosome being processed (RefGenomeInterface.getBases).
/// None stands for the cases where Java returns null or throws.
pub trait RefBases {
    fn get_bases(&self, pos_start: i32, pos_end: i32) -> Option<Vec<u8>>;
}

impl RefBases for RefSequence {
    fn get_bases(&self, pos_start: i32, pos_end: i32) -> Option<Vec<u8>> {
        if pos_start < 1 || pos_end < pos_start {
            return None;
        }
        self.base_range(pos_start, pos_end).map(|b| b.to_vec())
    }
}

/// UltimaSnv constructor
fn build_snv(
    variant: &SimpleVariant, ref_bases: &[u8], ref_var_index: i32, left_read_base: u8, right_read_base: u8,
    ref_genome: &dyn RefBases,
) -> Ex<Adjustments> {
    // the ref bases start 1 base prior to the variant
    let ref_base = at(ref_bases, ref_var_index)?;
    let alt_base = variant.alt[0];
    let left_base = at(ref_bases, ref_var_index - 1)?;
    let right_base = at(ref_bases, ref_var_index + 1)?;

    let mut adjustments = Adjustments::default();
    if left_base == right_base && (left_base == ref_base || left_base == alt_base) {
        // both straddling bases match either the ref or alt, so no need to check for HP adjustment
        return Ok(adjustments);
    }

    // one side will be a delete (either partial or complete), the other a 1-base insert
    let left_matches_ref = left_base == ref_base;
    let left_matches_alt = left_base == alt_base;
    let right_matches_ref = right_base == ref_base;
    let right_matches_alt = right_base == alt_base;
    let left_side_delete = (left_matches_ref && !left_matches_alt) || right_matches_alt;

    let mut lower_hp_length = 0;
    if left_matches_alt || left_matches_ref {
        let lower_ref_base_end = if left_matches_ref { variant.pos } else { variant.pos - 1 };
        let lower_ref_base_start = lower_ref_base_end - MAX_HOMOPOLYMER;
        let lower_ref_bases = ref_genome.get_bases(lower_ref_base_start, lower_ref_base_end).ok_or(())?;
        let last = lower_ref_bases.len() as i32 - 1;
        let hp_base = at(&lower_ref_bases, last)?;
        lower_hp_length = find_homopolymer_length(&lower_ref_bases, hp_base, last, false);
    }

    if left_matches_alt || !left_side_delete {
        // a HP expansion of 1 base on the left
        let new_hp_length = lower_hp_length + 1;
        let left_hp_end_index = 0;
        let left_hp_start_index = left_hp_end_index - new_hp_length + 1;
        adjustments.left_adjust = Some(UltimaHomopolymerAdjustment::new(left_hp_start_index, left_hp_end_index, -1));
    } else if left_matches_ref {
        // a HP contraction of 1 base on the left, starting with the variant's position
        let new_lower_hp_length = lower_hp_length - 1;
        let left_hp_end_index = -1;
        let left_hp_start_index = left_hp_end_index - new_lower_hp_length + 1;
        adjustments.left_adjust = Some(UltimaHomopolymerAdjustment::new(left_hp_start_index, left_hp_end_index, 1));
    } else {
        // HP deletion
        adjustments.left_deletion =
            Some(UltimaHomopolymerDeletion::from_mnv(0, at(ref_bases, 1)?, alt_base, left_read_base, right_read_base));
    }

    let mut upper_hp_length = 0;
    if right_matches_alt || right_matches_ref {
        // 1 being the variant's base, 2 being the one after
        let compare_start_index = if right_matches_ref { 1 } else { 2 };
        let hp_base = at(ref_bases, 2)?;
        upper_hp_length = find_homopolymer_length(ref_bases, hp_base, compare_start_index, true);
    }

    if right_matches_alt || left_side_delete {
        // a HP expansion of 1 base on the right
        let new_hp_length = upper_hp_length + 1;
        let right_hp_start_index = 0;
        let right_hp_end_index = right_hp_start_index + new_hp_length - 1;
        adjustments.right_adjust = Some(UltimaHomopolymerAdjustment::new(right_hp_start_index, right_hp_end_index, -1));
    } else if right_matches_ref {
        // a HP contraction of 1 base on the right, starting with the variant's position
        let new_hp_length = upper_hp_length - 1;
        let right_hp_start_index = 1;
        let right_hp_end_index = right_hp_start_index + new_hp_length - 1;
        adjustments.right_adjust = Some(UltimaHomopolymerAdjustment::new(right_hp_start_index, right_hp_end_index, 1));
    } else {
        adjustments.right_deletion =
            Some(UltimaHomopolymerDeletion::from_mnv(0, at(ref_bases, 1)?, alt_base, left_read_base, right_read_base));
    }
    Ok(adjustments)
}

/// UltimaBaseShift.getAdjustments
#[allow(clippy::too_many_arguments)]
fn base_shift_adjustments(
    variant: &SimpleVariant, ref_bases: &[u8], ref_var_index: i32, left_read_base: u8, right_read_base: u8,
    ref_genome: &dyn RefBases, mnv_length: i32, valid_shift_left: bool,
) -> Ex<Adjustments> {
    // one side will be a delete (either partial or complete), the other a 1-base insert
    let mut adjustments = Adjustments::default();
    let left_ref_base = at(ref_bases, ref_var_index)?;
    let right_ref_base = at(ref_bases, ref_var_index + mnv_length - 1)?;
    let left_alt_base = variant.alt[0];
    let right_alt_base = variant.alt[(mnv_length - 1) as usize];
    let next_left_ref_base = at(ref_bases, ref_var_index - 1)?;
    let next_right_ref_base = at(ref_bases, ref_var_index + mnv_length)?;

    if valid_shift_left {
        // find the deletion or contraction on the left
        if left_ref_base == next_left_ref_base {
            let lower_ref_base_end = variant.pos;
            let lower_ref_base_start = lower_ref_base_end - MAX_HOMOPOLYMER;
            let lower_ref_bases = ref_genome.get_bases(lower_ref_base_start, lower_ref_base_end).ok_or(())?;
            let lower_hp_length =
                find_homopolymer_length(&lower_ref_bases, left_ref_base, lower_ref_bases.len() as i32 - 1, false);
            let new_lower_hp_length = lower_hp_length - 1;
            let left_hp_end_index = -1;
            let left_hp_start_index = left_hp_end_index - new_lower_hp_length + 1;
            adjustments.left_adjust = Some(UltimaHomopolymerAdjustment::new(left_hp_start_index, left_hp_end_index, 1));
        } else {
            let right_base = variant.alt[1];
            adjustments.left_deletion =
                Some(UltimaHomopolymerDeletion::from_mnv(0, left_ref_base, left_alt_base, left_read_base, right_base));
        }

        // check if the insertion on the right is an expansion or a novel insertion
        let right_hp_start_index = mnv_length - 1;
        if right_alt_base == next_right_ref_base {
            // the alt has expanded an existing homopolymer
            let upper_ref_base_start = variant.position_end() + 1;
            let upper_ref_base_end = upper_ref_base_start + MAX_HOMOPOLYMER;
            let upper_ref_bases = ref_genome.get_bases(upper_ref_base_start, upper_ref_base_end).ok_or(())?;
            let upper_hp_length = find_homopolymer_length(&upper_ref_bases, next_right_ref_base, 0, true);
            let new_upper_hp_length = upper_hp_length + 1;
            let right_hp_end_index = right_hp_start_index + new_upper_hp_length - 1;
            adjustments.right_adjust =
                Some(UltimaHomopolymerAdjustment::new(right_hp_start_index, right_hp_end_index, -1));
        } else {
            // a novel insertion, so a HP expansion of 1 base on the right
            adjustments.right_adjust =
                Some(UltimaHomopolymerAdjustment::new(right_hp_start_index, right_hp_start_index, -1));
        }
    } else {
        if right_ref_base == next_right_ref_base {
            let upper_ref_base_start = variant.position_end();
            let upper_ref_base_end = upper_ref_base_start + MAX_HOMOPOLYMER;
            let upper_ref_bases = ref_genome.get_bases(upper_ref_base_start, upper_ref_base_end).ok_or(())?;
            let upper_hp_length = find_homopolymer_length(&upper_ref_bases, right_ref_base, 0, true);
            let new_upper_hp_length = upper_hp_length - 1;
            let right_hp_start_index = mnv_length;
            let right_hp_end_index = right_hp_start_index + new_upper_hp_length - 1;
            adjustments.right_adjust =
                Some(UltimaHomopolymerAdjustment::new(right_hp_start_index, right_hp_end_index, 1));
        } else {
            let left_base = variant.alt[(mnv_length - 2) as usize];
            adjustments.right_deletion = Some(UltimaHomopolymerDeletion::from_mnv(
                mnv_length - 1, right_ref_base, right_alt_base, left_base, right_read_base,
            ));
        }

        let left_hp_end_index = 0;
        if left_alt_base == next_left_ref_base {
            let lower_ref_base_end = variant.pos - 1;
            let lower_ref_base_start = lower_ref_base_end - MAX_HOMOPOLYMER;
            let lower_ref_bases = ref_genome.get_bases(lower_ref_base_start, lower_ref_base_end).ok_or(())?;
            let lower_hp_length =
                find_homopolymer_length(&lower_ref_bases, next_left_ref_base, lower_ref_bases.len() as i32 - 1, false);
            let new_hp_length = lower_hp_length + 1;
            let left_hp_start_index = left_hp_end_index - new_hp_length + 1;
            adjustments.left_adjust = Some(UltimaHomopolymerAdjustment::new(left_hp_start_index, left_hp_end_index, -1));
        } else {
            adjustments.left_adjust = Some(UltimaHomopolymerAdjustment::new(left_hp_end_index, left_hp_end_index, -1));
        }
    }
    Ok(adjustments)
}

/// UltimaBaseShift.refHpLen
fn ref_hp_len(adjustment: &Option<UltimaHomopolymerAdjustment>) -> i32 {
    adjustment.map_or(1, |a| a.hp_end_index - a.hp_start_index + 1 + a.ref_adjust_count)
}

/// UltimaBaseShift.preferOtherAdjustments
fn prefer_other_adjustments(adjustments: &Adjustments, other: &Adjustments) -> bool {
    if !other.can_compute() {
        return false;
    } else if !adjustments.can_compute() {
        return true;
    }
    // if the shift direction is ambiguous, prefer the interpretation that modifies the longer homopolymer
    let max_adj_hp_len = ref_hp_len(&adjustments.left_adjust).max(ref_hp_len(&adjustments.right_adjust));
    let max_other_adj_hp_len = ref_hp_len(&other.left_adjust).max(ref_hp_len(&other.right_adjust));
    max_other_adj_hp_len > max_adj_hp_len
}

/// UltimaBaseShift constructor
fn build_base_shift(
    variant: &SimpleVariant, ref_bases: &[u8], ref_var_index: i32, left_read_base: u8, right_read_base: u8,
    ref_genome: &dyn RefBases,
) -> Ex<Adjustments> {
    let mnv_length = variant.ref_len();
    // a base shift is a deletion or contraction on one side, the remaining bases shifting and the last
    // MNV base an insert; ie ACG -> CGT for a left shift, ACG -> TAC for a right shift
    let n = (mnv_length - 1) as usize;
    let valid_shift_left = (0..n).all(|i| variant.ref_[i + 1] == variant.alt[i]);
    let valid_shift_right = (0..n).all(|i| variant.ref_[i] == variant.alt[i + 1]);
    if !valid_shift_left && !valid_shift_right {
        return Ok(Adjustments::default());
    }
    let mut adjustments = base_shift_adjustments(
        variant, ref_bases, ref_var_index, left_read_base, right_read_base, ref_genome, mnv_length, valid_shift_left,
    )?;
    if valid_shift_left && valid_shift_right {
        // if both are viable options, take the higher prior shift direction
        let other = base_shift_adjustments(
            variant, ref_bases, ref_var_index, left_read_base, right_read_base, ref_genome, mnv_length, false,
        )?;
        if prefer_other_adjustments(&adjustments, &other) {
            adjustments = other;
        }
    }
    Ok(adjustments)
}

/// UltimaQualModel and its subclasses
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UltimaQualModel {
    HomopolymerAdjustment(UltimaHomopolymerAdjustment),
    HomopolymerDeletion(UltimaHomopolymerDeletion),
    HomopolymerTransition(UltimaHomopolymerTransitionDeletion),
    Snv(Adjustments),
    BaseShift(Adjustments),
    /// UltimaQualModelBuilder.OtherVariant: the fall-back for no valid model
    Other,
}

impl UltimaQualModel {
    pub fn model_type(&self) -> UltimaModelType {
        match self {
            UltimaQualModel::HomopolymerAdjustment(_) => UltimaModelType::HomopolymerAdjustment,
            UltimaQualModel::HomopolymerDeletion(_) => UltimaModelType::HomopolymerDeletion,
            UltimaQualModel::HomopolymerTransition(_) => UltimaModelType::HomopolymerTransition,
            UltimaQualModel::Snv(_) => UltimaModelType::Snv,
            UltimaQualModel::BaseShift(_) => UltimaModelType::BaseShift,
            UltimaQualModel::Other => UltimaModelType::Other,
        }
    }
    pub fn can_compute(&self) -> bool {
        match self {
            UltimaQualModel::Snv(a) | UltimaQualModel::BaseShift(a) => a.can_compute(),
            UltimaQualModel::Other => false,
            _ => true,
        }
    }
    pub fn calculate_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        match self {
            UltimaQualModel::HomopolymerAdjustment(m) => m.calculate_qual(read, var_read_index),
            UltimaQualModel::HomopolymerDeletion(m) => m.calculate_qual(read, var_read_index),
            UltimaQualModel::HomopolymerTransition(m) => m.calculate_qual(read, var_read_index),
            UltimaQualModel::Snv(a) | UltimaQualModel::BaseShift(a) => a.calculate_qual(read, var_read_index),
            UltimaQualModel::Other => INVALID_BASE_QUAL,
        }
    }
    /// the model as an UltimaBaseShift / UltimaSnv
    #[allow(dead_code)] // used by the ported tests
    pub fn adjustments(&self) -> Option<&Adjustments> {
        match self {
            UltimaQualModel::Snv(a) | UltimaQualModel::BaseShift(a) => Some(a),
            _ => None,
        }
    }
    /// the model as an UltimaHomopolymerAdjustment
    #[allow(dead_code)] // used by the ported tests
    pub fn homopolymer_adjustment(&self) -> Option<&UltimaHomopolymerAdjustment> {
        match self {
            UltimaQualModel::HomopolymerAdjustment(m) => Some(m),
            _ => None,
        }
    }
}

impl Adjustments {
    /// UltimaBaseShift.isLeftShift
    #[allow(dead_code)] // used by the ported tests
    pub fn is_left_shift(&self) -> bool {
        self.left_deletion.is_some() || self.left_adjust.map_or(false, |a| a.ref_adjust_count > 0)
    }
    /// UltimaBaseShift.isRightShift
    #[allow(dead_code)] // used by the ported tests
    pub fn is_right_shift(&self) -> bool {
        self.right_deletion.is_some() || self.right_adjust.map_or(false, |a| a.ref_adjust_count > 0)
    }
}

// --- UltimaQualModelBuilder ---

/// UltimaQualModelBuilder.getStraddlingReadBases
pub fn get_straddling_read_bases(variant: &SimpleVariant, read_bases: &[u8], var_read_index: i32) -> Ex<[u8; 2]> {
    if variant.is_insert() {
        return Ok([0, 0]); // unused
    }
    if variant.is_delete() {
        Ok([at(read_bases, var_read_index)?, at(read_bases, var_read_index + 1)?])
    } else {
        Ok([at(read_bases, var_read_index - 1)?, at(read_bases, var_read_index + variant.alt_len())?])
    }
}

/// UltimaQualModelBuilder.canSkipRealignedModels
pub fn can_skip_realigned_models(ctx: &VariantReadContext) -> Ex<bool> {
    let indel_length = ctx.variant.indel_length();
    let core_length = ctx.core_length();
    if core_length != ctx.ref_bases.len() as i32 + indel_length {
        return Ok(false);
    }
    let mut read_index = ctx.core_index_start;
    if ctx.variant.is_indel() {
        let mut ref_index = 0;
        while read_index < ctx.read_bases.len() as i32 && ref_index < ctx.ref_bases.len() as i32 {
            if read_index == ctx.var_index + 1 {
                if ctx.variant.is_delete() {
                    ref_index += indel_length.abs();
                } else {
                    read_index += indel_length;
                }
            }
            if at(&ctx.read_bases, read_index)? != at(&ctx.ref_bases, ref_index)? {
                return Ok(false);
            }
            read_index += 1;
            ref_index += 1;
        }
    } else {
        let alt_length = ctx.variant.alt_len();
        for i in 0..core_length {
            if !(read_index >= ctx.var_index && read_index <= ctx.var_index + (alt_length - 1))
                && at(&ctx.read_bases, read_index)? != at(&ctx.ref_bases, i)?
            {
                return Ok(false);
            }
            read_index += 1;
        }
    }
    Ok(true)
}

/// findHomopolymerTransitionCandidate: the index within the deleted ref bases of a transition, otherwise 0
fn find_homopolymer_transition_candidate(variant: &SimpleVariant) -> i32 {
    if !variant.is_delete() || variant.ref_.len() < 3 {
        return 0;
    }
    let first_del_base = variant.ref_[1];
    (2..variant.ref_.len()).find(|&i| variant.ref_[i] != first_del_base).map_or(0, |i| i as i32)
}

/// isHomopolymerDeletion
fn is_homopolymer_deletion(variant: &SimpleVariant, ref_bases: &[u8]) -> Ex<bool> {
    let ref_base = at(ref_bases, 1)?;
    let deleted_base = at(ref_bases, 2)?;
    let del_length = variant.indel_length_abs();
    let post_del_ref_base = at(ref_bases, 1 + variant.ref_len())?;
    if ref_base == deleted_base || post_del_ref_base == deleted_base {
        return Ok(false);
    }
    // for a DEL of 2 bases this checks the 3rd base, being the second deleted base, and then stops
    for i in 3..=1 + del_length {
        if at(ref_bases, i)? != deleted_base {
            return Ok(false);
        }
    }
    Ok(true)
}

fn single_distinct_base(bases: &[u8]) -> bool {
    bases.windows(2).all(|w| w[0] == w[1])
}

/// UltimaQualModelBuilder.buildContext. The straddling read bases are the read bases around the variant.
pub fn build_qual_model(
    ref_genome: &dyn RefBases, variant: &SimpleVariant, straddling_read_bases: [u8; 2], ref_masks: Option<&[RefMask]>,
) -> Ex<UltimaQualModel> {
    let max_homopolymer_length = variant.ref_len().max(MAX_HOMOPOLYMER);
    let ref_base_end = variant.pos + max_homopolymer_length + 1;

    // extract sufficient ref bases to set the context for most scenarios (not for homopolymer transition)
    let orig_ref_bases = ref_genome.get_bases(variant.pos - 1, ref_base_end).ok_or(())?;
    let mut ref_bases = orig_ref_bases.clone();
    if let Some(ref_masks) = ref_masks {
        for ref_mask in ref_masks {
            let start_index = (ref_mask.pos_start - variant.pos + 1).max(0);
            let end_index = (ref_mask.pos_end - variant.pos + 1).min(ref_bases.len() as i32 - 1);
            for i in start_index..=end_index {
                ref_bases[i as usize] = ref_mask.base_mask;
            }
        }
    }
    let ref_var_index = 1; // since the ref bases above start 1 base before the variant

    if variant.is_indel() {
        if variant.is_delete() && is_homopolymer_deletion(variant, &orig_ref_bases)? {
            // the HP base is the first ref base after the variant's position, ie the first base deleted;
            // the straddling read bases are the variant's ref and the first ref base after the delete
            return Ok(UltimaQualModel::HomopolymerDeletion(UltimaHomopolymerDeletion::from_delete(
                at(&orig_ref_bases, 2)?, straddling_read_bases[0], straddling_read_bases[1],
            )));
        } else if variant.is_delete() && is_homopolymer_deletion(variant, &ref_bases)? {
            return Ok(UltimaQualModel::HomopolymerDeletion(UltimaHomopolymerDeletion::from_delete(
                at(&ref_bases, 2)?, straddling_read_bases[0], straddling_read_bases[1],
            )));
        }

        let homopolymer_transition_index = find_homopolymer_transition_candidate(variant);
        if homopolymer_transition_index > 0 {
            let lower_hp_end = variant.pos + homopolymer_transition_index - 1;
            let upper_hp_start = lower_hp_end + 1;
            let lower_ref_base_start = lower_hp_end - MAX_HOMOPOLYMER;
            let upper_ref_base_end = upper_hp_start + MAX_HOMOPOLYMER;
            let lower_ref_bases = ref_genome.get_bases(lower_ref_base_start, lower_hp_end).ok_or(())?;
            let upper_ref_bases = ref_genome.get_bases(upper_hp_start, upper_ref_base_end).ok_or(())?;

            let lower_last = lower_ref_bases.len() as i32 - 1;
            let lower_hp_base = at(&lower_ref_bases, lower_last)?;
            let lower_hp_length = find_homopolymer_length(&lower_ref_bases, lower_hp_base, lower_last, false);
            let upper_hp_base = at(&upper_ref_bases, 0)?;
            let upper_hp_length = find_homopolymer_length(&upper_ref_bases, upper_hp_base, 0, true);

            if lower_hp_length > 2 && upper_hp_length > 2 && variant.ref_len() + variant.pos - 1 <= upper_ref_base_end {
                return Ok(UltimaQualModel::HomopolymerTransition(UltimaHomopolymerTransitionDeletion::new(
                    variant, homopolymer_transition_index, lower_hp_length, upper_hp_length,
                )));
            }
        }

        // the HP search starts 1 base after the variant's ref position
        let hp_base = if variant.is_insert() { variant.alt[1] } else { at(&ref_bases, ref_var_index + 1)? };
        let homopolymer_length = find_homopolymer_length(&ref_bases, hp_base, ref_var_index + 1, true);

        let single_ref_base = single_distinct_base(&variant.ref_[1..]);
        let single_alt_base = single_distinct_base(&variant.alt[1..]);
        if homopolymer_length <= MAX_HOMOPOLYMER && single_ref_base && single_alt_base {
            let homopolymer_start_index = 1;
            // adjust the HP length by the diff observed in the variant
            let homopolymer_end_index = homopolymer_start_index + homopolymer_length - 1 + variant.indel_length();
            let ref_adjust_count = -variant.indel_length();
            return Ok(UltimaQualModel::HomopolymerAdjustment(UltimaHomopolymerAdjustment::new(
                homopolymer_start_index, homopolymer_end_index, ref_adjust_count,
            )));
        }
    } else if variant.is_snv() {
        let model = UltimaQualModel::Snv(build_snv(
            variant, &ref_bases, ref_var_index, straddling_read_bases[0], straddling_read_bases[1], ref_genome,
        )?);
        if model.can_compute() {
            return Ok(model);
        }
    } else {
        // MNV case which may be a base shift
        let model = UltimaQualModel::BaseShift(build_base_shift(
            variant, &ref_bases, ref_var_index, straddling_read_bases[0], straddling_read_bases[1], ref_genome,
        )?);
        if model.can_compute() {
            return Ok(model);
        }
    }
    Ok(UltimaQualModel::Other)
}

// --- UltimaRealignedQualModel(s) and their builder ---

/// UltimaRealignedQualModel
#[derive(Clone, Debug)]
#[allow(dead_code)] // the var index and variant ref index are only for Java's unused MSI jitter qual cache
pub struct UltimaRealignedQualModel {
    pub variant: SimpleVariant,
    pub qual_model: UltimaQualModel,
    pub var_read_index_offset: i32,
    pub var_index: i32,
    pub variant_ref_index: i32,
}

impl UltimaRealignedQualModel {
    pub fn calculate_qual(&self, read: &Read, var_read_index: i32) -> i8 {
        self.qual_model.calculate_qual(read, var_read_index + self.var_read_index_offset)
    }
    #[allow(dead_code)] // used by the ported tests
    pub fn model_type(&self) -> UltimaModelType {
        self.qual_model.model_type()
    }
}

/// UltimaRealignedQualModels
#[derive(Clone, Debug)]
pub struct UltimaRealignedQualModels {
    pub original_qual_model: UltimaQualModel,
    pub realigned_qual_models: Vec<UltimaRealignedQualModel>,
}

impl UltimaRealignedQualModels {
    /// calculateQual: the minimum qual across the models
    pub fn calculate_qual(
        &self, uses_msi_indel_error_qual: bool, msi_indel_error_qual: i8, is_indel: bool, read_index: i32, read: &Read,
    ) -> f64 {
        let mut model_qual: Option<i8> = None;
        if self.original_qual_model.model_type() == UltimaModelType::Other && uses_msi_indel_error_qual {
            model_qual = Some(msi_indel_error_qual);
        } else if self.original_qual_model.can_compute() {
            model_qual = Some(self.original_qual_model.calculate_qual(read, read_index));
        }
        if matches!(model_qual, Some(q) if q < 0) {
            return INVALID_BASE_QUAL as f64;
        }
        for realigned in &self.realigned_qual_models {
            if !realigned.qual_model.can_compute() {
                continue;
            }
            let realigned_model_qual = realigned.calculate_qual(read, read_index);
            if realigned_model_qual < 0 {
                return INVALID_BASE_QUAL as f64;
            }
            model_qual = Some(model_qual.map_or(realigned_model_qual, |q| q.min(realigned_model_qual)));
        }
        // the default is the max raw qual, allowing for the fixed penalty deducted later
        let model_qual =
            model_qual.unwrap_or((max_raw_qual() as i32 + DEFAULT_BASE_QUAL_FIXED_PENALTY) as i8) as f64;
        if is_indel {
            // floor the base qual due to realignment issues from nearby MSI sites dominating indel errors
            // at quals above 40: probabilityToPhredQual(MSI_JITTER_DEFAULT_ERROR_RATE)
            jmin(model_qual, 40.0)
        } else {
            model_qual
        }
    }
}

/// UltimaRealignedQualModelBuilder.mergeSandwichedHomopolymers
pub fn merge_sandwiched_homopolymers(
    ctx: &VariantReadContext, ref_homopolymers0: &[Homopolymer], read_homopolymers0: &[Homopolymer],
) -> MergedHomopolymers {
    let mut ref_homopolymers = ref_homopolymers0.to_vec();
    let mut read_homopolymers = read_homopolymers0.to_vec();
    let mut merged = MergedHomopolymers::default();

    let mut ref_index = 0;
    let mut read_index = 0;
    let mut read_bases_consumed = 0;
    let mut ref_bases_consumed = 0;
    let var_index_in_core = ctx.var_index - ctx.core_index_start;
    let is_indel = ctx.variant.is_indel();

    while ref_index < ref_homopolymers.len() && read_index < read_homopolymers.len() {
        let ref_homopolymer = ref_homopolymers[ref_index];
        let read_homopolymer = read_homopolymers[read_index];
        let ref_remaining = ref_homopolymers.len() - ref_index - 1;
        let read_remaining = read_homopolymers.len() - read_index - 1;

        if ref_homopolymer.base == read_homopolymer.base {
            let net_length = read_homopolymer.length - ref_homopolymer.length;
            let mut contracted = false;
            if net_length >= 2 {
                // look forward in ref
                if ref_index < ref_homopolymers.len() - 1 {
                    let mut extra_ref_index = 1;
                    let mut forward_ref_homopolymer = ref_homopolymers[ref_index + extra_ref_index];
                    let mut extra_ref_bases_consumed = forward_ref_homopolymer.length;
                    loop {
                        if forward_ref_homopolymer.base == ref_homopolymer.base {
                            let mask_position_start = ctx.core_position_start + ref_bases_consumed + ref_homopolymer.length;
                            let ref_mask_length: i32 =
                                ref_homopolymers[ref_index + 1..ref_index + extra_ref_index].iter().map(|h| h.length).sum();
                            let merged_length = ref_mask_length
                                + ref_homopolymers[ref_index].length
                                + ref_homopolymers[ref_index + extra_ref_index].length;
                            merged.ref_masks.push(RefMask {
                                pos_start: mask_position_start,
                                pos_end: mask_position_start + ref_mask_length - 1,
                                base_mask: ref_homopolymer.base,
                            });
                            ref_homopolymers[ref_index] = Homopolymer::new(ref_homopolymer.base, merged_length);
                            ref_homopolymers.drain(ref_index + 1..=ref_index + extra_ref_index);
                            contracted = true;
                            if !is_indel
                                && var_index_in_core >= read_bases_consumed
                                && var_index_in_core < read_bases_consumed + read_homopolymer.length
                            {
                                merged.variant_in_merged_homopolymers = true;
                            }
                            break;
                        } else if extra_ref_bases_consumed >= net_length {
                            break;
                        }
                        extra_ref_index += 1;
                        if ref_index + extra_ref_index >= ref_homopolymers.len() {
                            break;
                        }
                        forward_ref_homopolymer = ref_homopolymers[ref_index + extra_ref_index];
                        extra_ref_bases_consumed += forward_ref_homopolymer.length;
                    }
                }
            } else if net_length <= -2 {
                // look forward in read
                if read_index < read_homopolymers.len() - 1 {
                    let mut extra_read_index = 1;
                    let mut forward_read_homopolymer = read_homopolymers[read_index + extra_read_index];
                    let mut extra_read_bases_consumed = forward_read_homopolymer.length;
                    loop {
                        if forward_read_homopolymer.base == ref_homopolymer.base {
                            let merged_length: i32 =
                                read_homopolymers[read_index..=read_index + extra_read_index].iter().map(|h| h.length).sum();
                            read_homopolymers[read_index] = Homopolymer::new(ref_homopolymer.base, merged_length);
                            read_homopolymers.drain(read_index + 1..=read_index + extra_read_index);
                            contracted = true;
                            if !is_indel
                                && var_index_in_core >= read_bases_consumed
                                && var_index_in_core < read_bases_consumed + merged_length
                            {
                                merged.variant_in_merged_homopolymers = true;
                            }
                            break;
                        } else if extra_read_bases_consumed >= -net_length {
                            break;
                        }
                        extra_read_index += 1;
                        if read_index + extra_read_index >= read_homopolymers.len() {
                            break;
                        }
                        forward_read_homopolymer = read_homopolymers[read_index + extra_read_index];
                        extra_read_bases_consumed += forward_read_homopolymer.length;
                    }
                }
            }
            if !contracted {
                merged.ref_homopolymers.push(ref_homopolymer);
                merged.read_homopolymers.push(read_homopolymer);
                ref_index += 1;
                ref_bases_consumed += ref_homopolymer.length;
                read_index += 1;
                read_bases_consumed += read_homopolymer.length;
            }
        } else if read_remaining <= ref_remaining {
            // novel delete
            merged.ref_homopolymers.push(ref_homopolymer);
            ref_index += 1;
            ref_bases_consumed += ref_homopolymer.length;
        } else {
            // novel insert
            merged.read_homopolymers.push(read_homopolymer);
            read_index += 1;
            read_bases_consumed += read_homopolymer.length;
        }
    }
    merged.ref_homopolymers.extend_from_slice(&ref_homopolymers[ref_index..]);
    merged.read_homopolymers.extend_from_slice(&read_homopolymers[read_index..]);
    merged
}

#[derive(Clone, Copy)]
struct HomopolymerVariant {
    base: u8,
    length: i32,
    variant_pos: i32,
    read_core_index: i32,
    prior_base: u8,
}

/// The matched-homopolymer state that createRealignedQualModels receives.
#[derive(Clone, Copy)]
struct LastMatched {
    read_core_index: i32,
    ref_pos: i32,
    base: u8,
    length: i32,
    prior_ref_base: u8,
    prior_read_base: u8,
}

/// UltimaRealignedQualModelBuilder.createRealignedQualModels
fn create_realigned_qual_models(
    realigned_qual_models: &mut Vec<UltimaRealignedQualModel>, ref_genome: &dyn RefBases, ctx: &VariantReadContext,
    ref_masks: &[RefMask], del_homopolymers: &[HomopolymerVariant], ins_homopolymers: &[HomopolymerVariant],
    last: &LastMatched,
) -> Ex<()> {
    if del_homopolymers.is_empty() && ins_homopolymers.is_empty() {
        return Ok(());
    }
    if last.ref_pos == -1 || last.read_core_index == -1 || last.length == -1 {
        // IllegalArgumentException: the first core homopolymers of the read and ref do not have the same base
        return Err(());
    }
    let orig_read_core_index = ctx.var_index - ctx.core_index_start;

    for (i, del_homopolymer) in del_homopolymers.iter().enumerate() {
        let mut var_read_index_offset = last.read_core_index - orig_read_core_index;
        let mut first_ref_base = del_homopolymer.prior_base;
        // in GGGGAG -> GGGG, don't attach the deleted G to the start of the poly-G
        let attached = del_homopolymer.base == last.base && i == 0;
        let mut first_alt_base = if attached { last.prior_read_base } else { last.base };
        if first_ref_base == INVALID_BASE || first_alt_base == INVALID_BASE {
            first_ref_base = at(&ctx.read_bases, ctx.core_index_start - 1)?;
            first_alt_base = first_ref_base;
        }
        if attached {
            var_read_index_offset -= last.length;
        }
        let mut ref_ = vec![first_ref_base];
        ref_.extend(std::iter::repeat(del_homopolymer.base).take(del_homopolymer.length as usize));
        let variant = SimpleVariant::new(del_homopolymer.variant_pos, &ref_, &[first_alt_base]);

        let adjusted_var_read_index = ctx.var_index + var_read_index_offset;
        let straddling_bases = get_straddling_read_bases(&variant, &ctx.read_bases, adjusted_var_read_index)?;
        let qual_model = build_qual_model(ref_genome, &variant, straddling_bases, Some(ref_masks))?;
        realigned_qual_models.push(UltimaRealignedQualModel {
            variant,
            qual_model,
            var_read_index_offset,
            var_index: adjusted_var_read_index,
            variant_ref_index: del_homopolymer.variant_pos - ctx.core_position_start,
        });
    }

    for (i, ins_homopolymer) in ins_homopolymers.iter().enumerate() {
        let var_read_index_offset = ins_homopolymer.read_core_index - orig_read_core_index;
        // in GGGG -> GGGGAG, don't attach the inserted G to the start of the poly-G
        let mut first_ref_base = if ins_homopolymer.base == last.base && i == 0 { last.prior_ref_base } else { last.base };
        let mut first_alt_base = ins_homopolymer.prior_base;
        if first_ref_base == INVALID_BASE || first_alt_base == INVALID_BASE {
            first_ref_base = at(&ctx.read_bases, ctx.core_index_start - 1)?;
            first_alt_base = first_ref_base;
        }
        let mut alt = vec![first_alt_base];
        alt.extend(std::iter::repeat(ins_homopolymer.base).take(ins_homopolymer.length as usize));
        let variant = SimpleVariant::new(ins_homopolymer.variant_pos, &[first_ref_base], &alt);

        let adjusted_var_read_index = ctx.var_index + var_read_index_offset;
        // S17, Java: the offset is passed where the read index is meant. The bases are not read for an
        // insert, so nothing changes.
        let straddling_index = if crate::compat::java("S17") { var_read_index_offset } else { adjusted_var_read_index };
        let straddling_bases = get_straddling_read_bases(&variant, &ctx.read_bases, straddling_index)?;
        let qual_model = build_qual_model(ref_genome, &variant, straddling_bases, Some(ref_masks))?;
        realigned_qual_models.push(UltimaRealignedQualModel {
            variant,
            qual_model,
            var_read_index_offset,
            var_index: adjusted_var_read_index,
            variant_ref_index: ins_homopolymer.variant_pos - ctx.core_position_start,
        });
    }
    Ok(())
}

/// extendHomopolymerVariants
fn extend_homopolymer_variants(variants: &mut Vec<HomopolymerVariant>, variant: HomopolymerVariant) {
    if let Some(last_variant) = variants.last_mut() {
        if last_variant.base == variant.base {
            last_variant.length += variant.length;
            return;
        }
    }
    variants.push(variant);
}

/// UltimaRealignedQualModelBuilder.getRealignedVariants
pub fn get_realigned_variants(
    ctx: &VariantReadContext, ref_genome: &dyn RefBases, ref_homopolymers: &[Homopolymer],
    read_homopolymers: &[Homopolymer], ref_masks: &[RefMask],
) -> Ex<Vec<UltimaRealignedQualModel>> {
    let mut realigned_variants = Vec::new();
    let mut last = LastMatched {
        read_core_index: -1,
        ref_pos: -1,
        base: 0,
        length: -1,
        prior_ref_base: INVALID_BASE,
        prior_read_base: INVALID_BASE,
    };
    let mut read_core_index = 0;
    let mut ref_index = 0;
    let mut read_index = 0;
    let mut ref_pos = ctx.core_position_start;
    let mut del_homopolymers: Vec<HomopolymerVariant> = Vec::new();
    let mut ins_homopolymers: Vec<HomopolymerVariant> = Vec::new();
    let prior = |list: &[Homopolymer], index: usize| if index == 0 { INVALID_BASE } else { list[index - 1].base };

    while ref_index < ref_homopolymers.len() && read_index < read_homopolymers.len() {
        let ref_homopolymer = ref_homopolymers[ref_index];
        let read_homopolymer = read_homopolymers[read_index];

        if ref_homopolymer.base == read_homopolymer.base {
            create_realigned_qual_models(
                &mut realigned_variants, ref_genome, ctx, ref_masks, &del_homopolymers, &ins_homopolymers, &last,
            )?;
            let matched_length = ref_homopolymer.length.min(read_homopolymer.length);
            last.length = matched_length;
            last.prior_ref_base = prior(ref_homopolymers, ref_index);
            last.prior_read_base = prior(read_homopolymers, read_index);
            last.base = ref_homopolymer.base;
            ref_index += 1;
            read_index += 1;
            last.ref_pos = ref_pos + matched_length - 1;
            ref_pos += ref_homopolymer.length;
            last.read_core_index = read_core_index + matched_length - 1;
            read_core_index += read_homopolymer.length;
            del_homopolymers.clear();
            ins_homopolymers.clear();

            if ref_homopolymer.length < read_homopolymer.length {
                ins_homopolymers.push(HomopolymerVariant {
                    base: read_homopolymer.base,
                    length: read_homopolymer.length - ref_homopolymer.length,
                    variant_pos: ref_pos - ref_homopolymer.length - 1,
                    read_core_index: read_core_index - read_homopolymer.length - 1,
                    prior_base: if read_index <= 1 { INVALID_BASE } else { read_homopolymers[read_index - 2].base },
                });
            } else if ref_homopolymer.length > read_homopolymer.length {
                del_homopolymers.push(HomopolymerVariant {
                    base: ref_homopolymer.base,
                    length: ref_homopolymer.length - read_homopolymer.length,
                    variant_pos: ref_pos - ref_homopolymer.length - 1,
                    read_core_index: read_core_index - read_homopolymer.length - 1,
                    prior_base: if ref_index <= 1 { INVALID_BASE } else { ref_homopolymers[ref_index - 2].base },
                });
            }
            continue;
        }

        let ref_homopolymers_left = ref_homopolymers.len() - ref_index - 1;
        let read_homopolymers_left = read_homopolymers.len() - read_index - 1;
        let del_homopolymer = HomopolymerVariant {
            base: ref_homopolymer.base,
            length: ref_homopolymer.length,
            variant_pos: ref_pos - 1,
            read_core_index: read_core_index - 1,
            prior_base: prior(ref_homopolymers, ref_index),
        };
        let ins_homopolymer = HomopolymerVariant {
            base: read_homopolymer.base,
            length: read_homopolymer.length,
            variant_pos: ref_pos - 1,
            read_core_index: read_core_index - 1,
            prior_base: prior(read_homopolymers, read_index),
        };
        if ref_homopolymers_left == read_homopolymers_left {
            extend_homopolymer_variants(&mut del_homopolymers, del_homopolymer);
            extend_homopolymer_variants(&mut ins_homopolymers, ins_homopolymer);
            ref_index += 1;
            read_index += 1;
            ref_pos += ref_homopolymer.length;
            read_core_index += read_homopolymer.length;
        } else if ref_homopolymers_left > read_homopolymers_left {
            extend_homopolymer_variants(&mut del_homopolymers, del_homopolymer);
            ref_index += 1;
            ref_pos += ref_homopolymer.length;
        } else {
            extend_homopolymer_variants(&mut ins_homopolymers, ins_homopolymer);
            read_index += 1;
            read_core_index += read_homopolymer.length;
        }
    }

    while read_index < read_homopolymers.len() {
        let read_homopolymer = read_homopolymers[read_index];
        extend_homopolymer_variants(
            &mut ins_homopolymers,
            HomopolymerVariant {
                base: read_homopolymer.base,
                length: read_homopolymer.length,
                variant_pos: ref_pos - 1,
                read_core_index: read_core_index - 1,
                prior_base: prior(read_homopolymers, read_index),
            },
        );
        read_index += 1;
        read_core_index += read_homopolymer.length;
    }
    while ref_index < ref_homopolymers.len() {
        let ref_homopolymer = ref_homopolymers[ref_index];
        extend_homopolymer_variants(
            &mut del_homopolymers,
            HomopolymerVariant {
                base: ref_homopolymer.base,
                length: ref_homopolymer.length,
                variant_pos: ref_pos - 1,
                read_core_index: read_core_index - 1,
                prior_base: prior(ref_homopolymers, ref_index),
            },
        );
        ref_index += 1;
        ref_pos += ref_homopolymer.length;
    }

    create_realigned_qual_models(
        &mut realigned_variants, ref_genome, ctx, ref_masks, &del_homopolymers, &ins_homopolymers, &last,
    )?;
    Ok(realigned_variants)
}

/// UltimaRealignedQualModelBuilder.buildUltimaRealignedQualModels
pub fn build_ultima_realigned_qual_models(
    ctx: &VariantReadContext, base_qual_model: UltimaQualModel, ref_genome: &dyn RefBases, skip_sandwich_masking: bool,
) -> Ex<UltimaRealignedQualModels> {
    let ref_homopolymers = get_homopolymers(&ctx.ref_bases, 0, ctx.ref_bases.len() as i32 - 1);
    let read_homopolymers = get_homopolymers(&ctx.read_bases, ctx.core_index_start, ctx.core_index_end);
    let merged = if skip_sandwich_masking {
        MergedHomopolymers { ref_homopolymers, read_homopolymers, ..Default::default() }
    } else {
        merge_sandwiched_homopolymers(ctx, &ref_homopolymers, &read_homopolymers)
    };
    let mut realigned_qual_models =
        get_realigned_variants(ctx, ref_genome, &merged.ref_homopolymers, &merged.read_homopolymers, &merged.ref_masks)?;
    realigned_qual_models.sort_by_key(|m| m.variant.pos);

    if realigned_qual_models.is_empty() && !merged.variant_in_merged_homopolymers {
        if !skip_sandwich_masking {
            return build_ultima_realigned_qual_models(ctx, base_qual_model, ref_genome, true);
        }
        // expected to have realigned variants, but none were found
        return Ok(UltimaRealignedQualModels { original_qual_model: base_qual_model, realigned_qual_models: Vec::new() });
    }
    Ok(UltimaRealignedQualModels { original_qual_model: base_qual_model, realigned_qual_models })
}

/// UltimaQualModelBuilder.setReadContextUltimaModels, returning the models
pub fn build_read_context_ultima_models(ref_genome: &dyn RefBases, ctx: &VariantReadContext) -> Ex<UltimaRealignedQualModels> {
    let straddling_read_bases = get_straddling_read_bases(&ctx.variant, &ctx.read_bases, ctx.var_index)?;
    let qual_model = build_qual_model(ref_genome, &ctx.variant, straddling_read_bases, None)?;
    if can_skip_realigned_models(ctx)? {
        Ok(UltimaRealignedQualModels { original_qual_model: qual_model, realigned_qual_models: Vec::new() })
    } else {
        build_ultima_realigned_qual_models(ctx, qual_model, ref_genome, false)
    }
}

// --- UltimaVariantData ---

/// UltimaVariantData: per-variant homopolymer and T0 qualities over the reads with full support
#[derive(Clone, Debug)]
pub struct UltimaVariantData {
    pub qual_models: Option<UltimaRealignedQualModels>,
    homopolymer_lengths: Vec<i32>,
    homopolymer_left_padding: i32,
    homopolymer_right_padding: i32,
    pub homopolymer_avg_quals: Vec<f64>,
    pub t0_avg_quals: Vec<f64>,
    homopolymer_avg_quals_count: i32,
}

impl UltimaVariantData {
    pub fn new(ctx: &VariantReadContext) -> Self {
        let homopolymers = get_homopolymers(&ctx.read_bases, ctx.core_index_start, ctx.core_index_end);
        // HP extensions from the flanks to the outer homopolymers
        let left_padding = homopolymers.first().map_or(0, |first| {
            ctx.read_bases[..ctx.core_index_start as usize].iter().rev().take_while(|&&b| b == first.base).count()
        });
        let right_padding = homopolymers.last().map_or(0, |last| {
            ctx.read_bases[(ctx.core_index_end + 1) as usize..].iter().take_while(|&&b| b == last.base).count()
        });
        let n = homopolymers.len();
        UltimaVariantData {
            qual_models: None,
            homopolymer_lengths: homopolymers.iter().map(|h| h.length).collect(),
            homopolymer_left_padding: left_padding as i32,
            homopolymer_right_padding: right_padding as i32,
            homopolymer_avg_quals: vec![0.0; n],
            t0_avg_quals: vec![0.0; n],
            homopolymer_avg_quals_count: 0,
        }
    }

    /// paddedHomopolymerLengths
    pub fn padded_homopolymer_lengths(&self) -> Vec<i32> {
        let n = self.homopolymer_lengths.len();
        let mut padded = self.homopolymer_lengths.clone();
        if n > 0 {
            padded[0] += self.homopolymer_left_padding;
            padded[n - 1] += self.homopolymer_right_padding;
        }
        padded
    }

    /// addReadSupportInfo: an exact match of the core, since quals are not applicable
    pub fn add_read_support_info(&mut self, ctx: &VariantReadContext, read: &Read, read_var_index: i32) {
        let start = read_var_index - ctx.left_core_length();
        let end = read_var_index + ctx.right_core_length();
        if start < 0 || end >= read.len() {
            return;
        }
        if &read.bases[start as usize..=end as usize] != ctx.core() {
            return;
        }
        // registerHomopolymerQuals and registerT0Quals; Java fails on a read without the tags
        assert!(!read.tp.is_empty() && !read.t0.is_empty(), "read {} has no Ultima tp or t0 tag", read.name);
        let n = self.homopolymer_avg_quals_count as f64;
        let mut read_index = start;
        for (i, &hp_length) in self.homopolymer_lengths.iter().enumerate() {
            let lookup_index = (if i == 0 { read_index + hp_length - 1 } else { read_index }) as usize;
            let mut homopolymer_qual = read.quals[lookup_index] as i8;
            if read.tp[lookup_index] != 0 && hp_length > 1 {
                homopolymer_qual = (homopolymer_qual as i32 - HALF_PHRED_SCORE_SCALING as i32) as i8;
            }
            let t0_qual = read.t0[lookup_index];
            self.homopolymer_avg_quals[i] = (self.homopolymer_avg_quals[i] * n + homopolymer_qual as f64) / (n + 1.0);
            self.t0_avg_quals[i] = (self.t0_avg_quals[i] * n + t0_qual as f64) / (n + 1.0);
            read_index += hp_length;
        }
        self.homopolymer_avg_quals_count += 1;
    }
}

// --- variant filters ---

/// MIN_HP_QUAL_MAP
fn min_hp_qual(length: i32) -> f64 {
    match length {
        1 => 21.5,
        2 => 20.0,
        3 => 18.5,
        4 => 17.0,
        5 => 15.5,
        6 => 14.0,
        7 => 12.0,
        _ => 10.0,
    }
}

/// PANEL_REPEAT_HP_QUAL_MAP
fn panel_repeat_hp_qual(length: i32) -> f64 {
    match length {
        1 => 19.5,
        2 => 17.5,
        3 => 15.5,
        4 => 14.0,
        5 => 13.5,
        6 => 10.0,
        _ => 6.0,
    }
}

const HP_MAX_LENGTHS: i32 = 5;
const HP_THRESHOLD_NON_INDEL_DEDUCTION: f64 = 2.5;
const GERMLINE_HET_MIN_EXPECTED_VAF: f64 = 0.4;
const GERMLINE_STRAND_BIAS_ADJUST_FACTOR: f64 = 50.0;
const GERMLINE_STRAND_BIAS_ADJUST: f64 = 5.0;

/// UltimaUtils.belowExpectedHpQuals. `is_panel_indel_repeat` is VariantFilters.isPanelIndelRepeatVariant
/// for the variant, which only applies when not in germline mode.
pub fn below_expected_hp_quals(primary_tumor: &ReadContextCounter, is_germline: bool, is_panel_indel_repeat: bool) -> bool {
    let ctx = &primary_tumor.ctx;
    if ctx.variant.is_snv() && !ctx.has_indel_in_core {
        return false;
    }
    let ultima_data = primary_tumor.ultima.as_ref().expect("Ultima variant data");
    let homopolymer_lengths = ultima_data.padded_homopolymer_lengths();
    if primary_tumor.is_long_indel() && *homopolymer_lengths.iter().max().expect("core homopolymers") < HP_MAX_LENGTHS {
        return false;
    }
    if !is_germline
        && primary_tumor.strong_alt_support() as f64 / primary_tumor.depth() as f64 >= GERMLINE_HET_MIN_EXPECTED_VAF
    {
        return false;
    }
    let repeat_length = ctx.max_repeat_count();
    let is_panel_indel_repeat = !is_germline && is_panel_indel_repeat;

    // for germline, add an offset to the threshold based on the variant's strand bias: it starts at 5 and
    // increases by 5 per 0.1 the alt read strand bias is below 40% or above 60%, unless the non-alt read
    // strand bias is on the same side
    let mut strand_bias_threshold_adjustment = 0.0;
    if is_germline {
        let alt_strand_bias = primary_tumor.alt_read_strand_bias.bias();
        let non_alt_strand_bias = primary_tumor.non_alt_read_strand_bias.bias();
        let lower_vaf_limit = GERMLINE_HET_MIN_EXPECTED_VAF;
        let upper_vaf_limit = 1.0 - GERMLINE_HET_MIN_EXPECTED_VAF;
        strand_bias_threshold_adjustment = GERMLINE_STRAND_BIAS_ADJUST;
        if alt_strand_bias < lower_vaf_limit {
            if non_alt_strand_bias >= lower_vaf_limit {
                strand_bias_threshold_adjustment += (lower_vaf_limit - alt_strand_bias) * GERMLINE_STRAND_BIAS_ADJUST_FACTOR;
            }
        } else if alt_strand_bias > upper_vaf_limit && non_alt_strand_bias <= upper_vaf_limit {
            strand_bias_threshold_adjustment += (alt_strand_bias - upper_vaf_limit) * GERMLINE_STRAND_BIAS_ADJUST_FACTOR;
        }
    }

    let msi_long_homopolymer = primary_tumor.qual_cache.is_msi_sample_and_variant
        && ctx.max_repeat.as_ref().map_or(false, |r| r.repeat_length() == 1)
        && repeat_length >= MSI_SAMPLE_REPEAT_LENGTH_MAX
        && ctx.variant.indel_length_abs() > 1;

    for (i, &length) in homopolymer_lengths.iter().enumerate() {
        let mut threshold = if is_panel_indel_repeat {
            panel_repeat_hp_qual(length)
        } else {
            let mut threshold = min_hp_qual(length);
            if !primary_tumor.is_indel {
                threshold -= HP_THRESHOLD_NON_INDEL_DEDUCTION;
            }
            threshold
        };
        threshold += strand_bias_threshold_adjustment;
        if msi_long_homopolymer && length >= MSI_SAMPLE_REPEAT_LENGTH_MAX {
            continue;
        }
        if ultima_data.homopolymer_avg_quals[i] < threshold {
            return true;
        }
    }
    false
}

/// UltimaUtils.belowExpectedT0Quals. `alt_all_one_side` is STRAND_BIAS_CALCS.allOneSide(readStrandBiasAlt).
pub fn below_expected_t0_quals(primary_tumor: &ReadContextCounter, alt_all_one_side: bool) -> bool {
    let ultima_data = primary_tumor.ultima.as_ref().expect("Ultima variant data");
    // the same call succeeded when the qual models were built
    let mut is_simple = can_skip_realigned_models(&primary_tumor.ctx).unwrap_or(false);
    if primary_tumor.ctx.variant.is_mnv() {
        let original_can_compute =
            ultima_data.qual_models.as_ref().expect("Ultima qual models").original_qual_model.can_compute();
        is_simple = !alt_all_one_side && !original_can_compute;
    }
    let threshold = if is_simple { T0_EXPECTED_QUAL_THRESHOLD_SIMPLE } else { T0_EXPECTED_QUAL_THRESHOLD_COMPLEX };
    let min_t0 = ultima_data.t0_avg_quals.iter().copied().fold(f64::INFINITY, f64::min);
    min_t0 < threshold
}

// --- UltimaCoreExtender (and hmf-common ReadCigarState) ---

#[derive(Clone, Copy, Debug)]
struct ReadCigarState {
    ref_position: i32,
    read_index: i32,
    cigar_index: i32,
    element_index: i32,
    element: CigarEl,
    valid: bool,
}

impl ReadCigarState {
    /// ReadCigarState.initialise: past any initial soft-clip
    fn initialise(ref_position: i32, cigar: &[CigarEl]) -> Ex<Self> {
        let first = *cigar.first().ok_or(())?;
        if first.op == Op::S {
            let second = *cigar.get(1).ok_or(())?;
            Ok(ReadCigarState {
                ref_position,
                read_index: first.len,
                cigar_index: 1,
                element_index: 0,
                element: second,
                valid: true,
            })
        } else {
            Ok(ReadCigarState { ref_position, read_index: 0, cigar_index: 0, element_index: 0, element: first, valid: true })
        }
    }

    /// the copy constructor
    fn copy(&self) -> Self {
        // S18, Java: the copy constructor sets the valid flag to true whatever the source holds.
        // Corrected: the flag is copied. Both callers copy a state that is valid, so nothing changes.
        ReadCigarState { valid: crate::compat::java("S18") || self.valid, ..*self }
    }

    fn op(&self) -> Op {
        self.element.op
    }

    /// ReadCigarState.moveState: stops at the end of the read's cigar
    fn move_state(&mut self, cigar: &[CigarEl], move_up: bool) {
        let last_operator = self.op();
        if move_up {
            self.element_index += 1;
            if self.element_index >= self.element.len {
                self.cigar_index += 1;
                if self.cigar_index >= cigar.len() as i32 {
                    self.valid = false;
                    return;
                }
                self.element = cigar[self.cigar_index as usize];
                self.element_index = 0;
            }
            // a soft-clip will have remained at the read start
            if last_operator != Op::S && self.op().consumes_ref() {
                self.ref_position += 1;
            }
            if self.op().consumes_read() {
                self.read_index += 1;
            }
        } else {
            self.element_index -= 1;
            if self.element_index < 0 {
                self.cigar_index -= 1;
                if self.cigar_index < 0 {
                    self.valid = false;
                    return;
                }
                self.element = cigar[self.cigar_index as usize];
                self.element_index = self.element.len - 1;
            }
            if last_operator != Op::S && self.op().consumes_ref() {
                self.ref_position -= 1;
            }
            if self.op().consumes_read() {
                self.read_index -= 1;
            }
        }
    }

    fn off_cigar(&self, cigar: &[CigarEl]) -> bool {
        self.cigar_index < 0 || self.cigar_index >= cigar.len() as i32
    }

    /// The moves that stay within the current element, taken at once: up to `limit` of them where
    /// each brings the target closer, otherwise all of them. Returns whether any were taken.
    fn move_within_element(&mut self, move_up: bool, limit: Option<i32>) -> bool {
        let room = if move_up { self.element.len - 1 - self.element_index } else { self.element_index };
        let steps = limit.map_or(room, |l| room.min(l));
        if steps <= 0 {
            return false;
        }
        let direction = if move_up { 1 } else { -1 };
        self.element_index += direction * steps;
        if self.op().consumes_ref() {
            self.ref_position += direction * steps;
        }
        if self.op().consumes_read() {
            self.read_index += direction * steps;
        }
        true
    }

    /// ReadCigarState.moveToRefPosition
    fn move_to_ref_position(&mut self, cigar: &[CigarEl], target_ref_position: i32) {
        let move_up = self.ref_position < target_ref_position;
        while self.ref_position != target_ref_position {
            let limit = self.op().consumes_ref().then(|| (target_ref_position - self.ref_position).abs());
            if self.move_within_element(move_up, limit) {
                continue;
            }
            self.move_state(cigar, move_up);
            if !self.valid {
                break;
            }
        }
        if self.ref_position != target_ref_position {
            self.valid = false;
        }
    }

    /// ReadCigarState.moveToIndex
    fn move_to_index(&mut self, cigar: &[CigarEl], target_index: i32) {
        while self.read_index != target_index {
            let limit = self.op().consumes_read().then(|| (target_index - self.read_index).abs());
            if self.move_within_element(self.read_index < target_index, limit) {
                continue;
            }
            self.move_state(cigar, self.read_index < target_index);
            if !self.valid {
                break;
            }
        }
        if self.read_index != target_index {
            self.valid = false;
        }
    }

    /// moveToRefPosition as Java has it, one base at a time
    #[cfg(test)]
    fn move_to_ref_position_by_base(&mut self, cigar: &[CigarEl], target_ref_position: i32) {
        let move_up = self.ref_position < target_ref_position;
        while self.ref_position != target_ref_position {
            self.move_state(cigar, move_up);
            if !self.valid {
                break;
            }
        }
        if self.ref_position != target_ref_position {
            self.valid = false;
        }
    }

    /// moveToIndex as Java has it, one base at a time
    #[cfg(test)]
    fn move_to_index_by_base(&mut self, cigar: &[CigarEl], target_index: i32) {
        while self.read_index != target_index {
            self.move_state(cigar, self.read_index < target_index);
            if !self.valid {
                break;
            }
        }
        if self.read_index != target_index {
            self.valid = false;
        }
    }
}

/// Not in the Java suite: moving a cigar state a whole element at a time ends in the state that
/// moving it a base at a time ends in.
#[cfg(test)]
#[test]
fn test_cigar_state_moves_equal_moves_by_base() {
    let mut seed: u64 = 99;
    let mut next = |bound: u64| -> u64 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) % bound
    };
    let describe = |s: &ReadCigarState| format!("{s:?}");
    for _ in 0..3000 {
        let mut cigar: Vec<CigarEl> = Vec::new();
        if next(3) == 0 {
            cigar.push(CigarEl::new(1 + next(6) as i32, Op::S));
        }
        cigar.push(CigarEl::new(1 + next(20) as i32, Op::M));
        for _ in 0..next(5) {
            let op = [Op::I, Op::D, Op::N, Op::I, Op::D][next(5) as usize];
            cigar.push(CigarEl::new(1 + next(4) as i32, op));
            cigar.push(CigarEl::new(1 + next(20) as i32, Op::M));
        }
        if next(3) == 0 {
            cigar.push(CigarEl::new(1 + next(6) as i32, Op::S));
        }
        let read_length: i32 = cigar.iter().filter(|e| e.op.consumes_read()).map(|e| e.len).sum();
        let ref_length: i32 = cigar.iter().filter(|e| e.op.consumes_ref()).map(|e| e.len).sum();
        let start = 1000;
        let mut fast = ReadCigarState::initialise(start, &cigar).unwrap();
        let mut by_base = fast;
        // a sequence of moves, each from where the last one ended, some to targets off the read
        for _ in 0..6 {
            if !fast.valid {
                break;
            }
            if next(2) == 0 {
                let target = next(read_length as u64 + 6) as i32 - 3;
                fast.move_to_index(&cigar, target);
                by_base.move_to_index_by_base(&cigar, target);
            } else {
                let target = start + next(ref_length as u64 + 6) as i32 - 3;
                fast.move_to_ref_position(&cigar, target);
                by_base.move_to_ref_position_by_base(&cigar, target);
            }
            assert_eq!(describe(&fast), describe(&by_base), "cigar {cigar:?}");
        }
    }
}

/// UltimaCoreExtender.moveState: for this routine the adjusted core and flanks cannot be within a soft-clip
fn extender_move_state(state: &mut ReadCigarState, cigar: &[CigarEl], move_up: bool, skip_deletes: bool) {
    state.move_state(cigar, move_up);
    // Java has no cigar bounds test here; a cigar does not start or end with a delete
    while skip_deletes && state.op() == Op::D && !state.off_cigar(cigar) {
        state.move_state(cigar, move_up);
    }
    if state.op() == Op::S {
        state.valid = false;
    }
}

fn ref_base(ref_sequence: &RefSequence, position: i32) -> Ex<u8> {
    at(&ref_sequence.bases, ref_sequence.index(position))
}

/// UltimaCoreExtender.requiresCoreExtension
fn requires_core_extension(
    read_bases: &[u8], ref_sequence: &RefSequence, cigar: &[CigarEl], state: &mut ReadCigarState, search_up: bool,
) -> Ex<bool> {
    if !state.valid {
        return Ok(true); // if unclear return an invalid state
    }
    let n = read_bases.len() as i32;
    if state.read_index < 0 || state.read_index >= n {
        return Ok(true);
    }
    let read_base = read_bases[state.read_index as usize];
    let mut ref_position = state.ref_position;
    let ref_base_ = ref_base(ref_sequence, state.ref_position)?;
    if read_base != ref_base_ {
        return Ok(true);
    }
    extender_move_state(state, cigar, search_up, true);
    if state.read_index < 0 || state.read_index >= n {
        return Ok(true);
    }
    let next_read_base = read_bases[state.read_index as usize];
    ref_position += if search_up { 1 } else { -1 };
    let next_ref_base = ref_base(ref_sequence, ref_position)?;
    Ok(read_base == next_read_base || ref_base_ == next_ref_base)
}

/// UltimaCoreExtender.findCoreExtension: extend the core while any homopolymer exists in either the ref or
/// the read, then exit when the ref matches the read at an M element
fn find_core_extension(
    read_bases: &[u8], ref_sequence: &RefSequence, cigar: &[CigarEl], state: &mut ReadCigarState, search_up: bool,
) -> Ex<bool> {
    let mut read_base = at(read_bases, state.read_index)?;
    let mut ref_position = state.ref_position;
    let mut ref_base_ = ref_base(ref_sequence, ref_position)?;
    let mut check_homopolymers = true;
    let mut required_extension = false;
    loop {
        extender_move_state(state, cigar, search_up, true);
        ref_position += if search_up { 1 } else { -1 };
        if !state.valid {
            break;
        }
        let next_read_base = at(read_bases, state.read_index)?;
        let Ok(mut next_ref_base) = ref_base(ref_sequence, ref_position) else {
            // !refSequence.containsPosition(refPosition)
            state.valid = false;
            return Ok(false);
        };
        let mut has_homopolymer = false;
        if check_homopolymers {
            has_homopolymer = read_base == next_read_base || ref_base_ == next_ref_base;
            if !has_homopolymer {
                check_homopolymers = false;
            }
        }
        if !required_extension && !has_homopolymer && read_base == ref_base_ {
            // immediate exit if no extension was required
            return Ok(false);
        }
        required_extension = true;
        if !has_homopolymer {
            next_ref_base = ref_base(ref_sequence, state.ref_position)?;
            if next_read_base == next_ref_base && matches!(state.op(), Op::M | Op::S) {
                break;
            }
        }
        read_base = next_read_base;
        ref_base_ = next_ref_base;
    }
    Ok(required_extension)
}

/// UltimaCoreExtender.findFlankIndex
fn find_flank_index(state: &mut ReadCigarState, cigar: &[CigarEl], required_flank_length: i32, move_up: bool) {
    if state.op() == Op::S {
        return;
    }
    let mut flank_length = 0;
    while flank_length < required_flank_length || state.op() == Op::I {
        extender_move_state(state, cigar, move_up, true);
        flank_length += 1;
        if !state.valid {
            break;
        }
    }
}

/// UltimaCoreExtender.buildReadCigar
fn build_read_cigar(mut state: ReadCigarState, cigar: &[CigarEl], flank_index_end: i32) -> Vec<CigarEl> {
    let mut elements = Vec::new();
    let mut cigar_length = 1;
    let mut cigar_operator = state.op();
    while state.read_index < flank_index_end {
        extender_move_state(&mut state, cigar, true, false);
        if state.off_cigar(cigar) {
            break; // Java would not terminate; the end index is always within the read
        }
        if state.op() == cigar_operator {
            cigar_length += 1;
        } else {
            elements.push(CigarEl::new(cigar_length, cigar_operator));
            cigar_length = 1;
            cigar_operator = state.op();
        }
    }
    elements.push(CigarEl::new(cigar_length, cigar_operator));
    elements
}

/// UltimaCoreExtender.extendUltimaCore: expand the core so that homopolymers are not cut off in the read
/// or the ref. Ok(None) is Java's null (no valid context); an unchanged core returns a copy of the input.
pub fn extend_ultima_core(
    read_bases: &[u8], ref_sequence: &RefSequence, read_alignment_start: i32, cigar: &[CigarEl],
    read_cigar_info: &ReadCigarInfo, flank_size: i32, in_append_mode: bool,
) -> Ex<Option<ReadCigarInfo>> {
    let mut core_position_start = read_cigar_info.core_position_start;
    let mut core_position_end = read_cigar_info.core_position_end;
    let mut flank_position_start = read_cigar_info.flank_position_start;
    let mut flank_index_start = read_cigar_info.flank_index_start;
    let mut flank_position_end = read_cigar_info.flank_position_end;
    let mut flank_index_end = read_cigar_info.flank_index_end;

    // move to the start of the core
    let mut lower_state = ReadCigarState::initialise(read_alignment_start, cigar)?;
    let mut upper_state = lower_state.copy();
    let read_alignment_end =
        read_alignment_start + cigar.iter().filter(|e| e.op.consumes_ref()).map(|e| e.len).sum::<i32>() - 1;

    let mut extended_start = false;
    let mut extended_end = false;
    let lower_in_soft_clip = read_cigar_info.flank_position_start < read_alignment_start;
    let upper_in_soft_clip = read_cigar_info.flank_position_end > read_alignment_end;

    lower_state.move_to_ref_position(cigar, read_cigar_info.core_position_start);
    upper_state.move_to_ref_position(cigar, read_cigar_info.core_position_end);

    if !lower_in_soft_clip {
        if lower_state.valid {
            extended_start = find_core_extension(read_bases, ref_sequence, cigar, &mut lower_state, false)?;
        }
        if extended_start {
            core_position_start = lower_state.ref_position;
            find_flank_index(&mut lower_state, cigar, flank_size, false);
            if !lower_state.valid && in_append_mode {
                lower_state = ReadCigarState::initialise(read_alignment_start, cigar)?;
                lower_state.move_to_index(cigar, read_cigar_info.flank_index_start);
            }
            flank_position_start = lower_state.ref_position;
            flank_index_start = lower_state.read_index;
        } else {
            // keep the existing flank position
            lower_state = ReadCigarState::initialise(read_alignment_start, cigar)?;
            lower_state.move_to_ref_position(cigar, read_cigar_info.flank_position_start);
        }
    } else {
        // check whether extension would be required and return an invalid context if so
        if requires_core_extension(read_bases, ref_sequence, cigar, &mut lower_state, false)? {
            return Ok(None);
        }
        // move to the existing flank start
        lower_state.move_to_index(cigar, read_cigar_info.flank_index_start);
    }

    if !upper_in_soft_clip {
        if upper_state.valid {
            extended_end = find_core_extension(read_bases, ref_sequence, cigar, &mut upper_state, true)?;
        }
        if extended_end {
            core_position_end = upper_state.ref_position;
            find_flank_index(&mut upper_state, cigar, flank_size, true);
            if !upper_state.valid && in_append_mode {
                upper_state = ReadCigarState::initialise(read_alignment_start, cigar)?;
                upper_state.move_to_index(cigar, read_cigar_info.flank_index_end);
            }
            flank_position_end = upper_state.ref_position;
            flank_index_end = upper_state.read_index;
        } else {
            upper_state = ReadCigarState::initialise(read_alignment_start, cigar)?;
            upper_state.move_to_ref_position(cigar, read_cigar_info.flank_position_end);
        }
    } else {
        if requires_core_extension(read_bases, ref_sequence, cigar, &mut upper_state, true)? {
            return Ok(None);
        }
        upper_state.move_to_index(cigar, read_cigar_info.flank_index_end);
    }

    if !extended_start && !extended_end {
        return Ok(Some(read_cigar_info.clone()));
    }
    if !lower_state.valid || !upper_state.valid {
        return Ok(None);
    }

    // build a new cigar between the 2 flank boundaries
    let new_cigar = build_read_cigar(lower_state.copy(), cigar, upper_state.read_index);
    Ok(Some(ReadCigarInfo {
        cigar: new_cigar,
        read_alignment_start,
        flank_position_start,
        flank_position_end,
        core_position_start,
        core_position_end,
        flank_index_start,
        flank_index_end,
    }))
}

/// The core read indices implied by an extended ReadCigarInfo (VariantReadContextBuilder.buildContext).
pub fn extended_core_indices(info: &ReadCigarInfo) -> (i32, i32) {
    let max_pos_start = info.flank_position_start.max(info.read_alignment_start);
    (
        info.flank_index_start + read_index_from_position(max_pos_start, &info.cigar, info.core_position_start),
        info.flank_index_start + read_index_from_position(max_pos_start, &info.cigar, info.core_position_end),
    )
}
