//! SageConstants (subset in use).
#![allow(dead_code)]

pub const DEFAULT_MIN_MAP_QUALITY: i32 = 10;
pub const DEFAULT_MAX_READ_DEPTH: i32 = 1000;
pub const DEFAULT_MAX_READ_DEPTH_PANEL: i32 = 100_000;
pub const DEFAULT_SLICE_SIZE: i32 = 100_000;
pub const REGION_BLOCK_SIZE: i32 = 100;
pub const DEFAULT_MAX_PARTITION_SLICES: i32 = 10;
pub const DEFAULT_READ_LENGTH: i32 = 151;
pub const DEFAULT_FLANK_LENGTH: i32 = 10;
pub const MIN_CORE_DISTANCE: i32 = 2;
pub const MAX_REPEAT_LENGTH: i32 = 5;
pub const MIN_REPEAT_COUNT: i32 = 3;
pub const OUTER_MIN_REPEAT_COUNT: i32 = 6;
pub const MIN_SECOND_CANDIDATE_FULL_READS: i32 = 3;
pub const MIN_SECOND_CANDIDATE_FULL_READS_PERC: f64 = 0.25;
pub const MATCHING_BASE_QUALITY: u8 = 20;
pub const CORE_LOW_QUAL_MISMATCH_FACTOR: i32 = 8;
pub const FLANK_LOW_QUAL_MISMATCHES: i32 = 3;
pub const SC_READ_EVENTS_FACTOR: f64 = 12.0;
pub const SC_INSERT_REF_TEST_LENGTH: i32 = 12;
pub const SC_INSERT_MIN_LENGTH: i32 = 5;
pub const MIN_INSERT_ALIGNMENT_OVERLAP: i32 = 5;
pub const LONG_INSERT_LENGTH: i32 = 11;
pub const LONG_GERMLINE_INSERT_READ_VS_REF_DIFF: i32 = 2;
pub const DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT: i32 = 2;
pub const DEFAULT_MAP_QUAL_FIXED_PENALTY: i32 = 0;
pub const DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY: i32 = 15;
pub const DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY: f64 = 7.0 * DEFAULT_READ_LENGTH as f64;

pub const LOW_BASE_QUAL_THRESHOLD: u8 = 26;
pub const BASE_QUAL_MINIMUM: u8 = 1;
pub const ILLUMINA_HIGH_BASE_QUAL: u8 = 30;

pub fn is_low_base_qual(q: u8) -> bool {
    q <= LOW_BASE_QUAL_THRESHOLD
}
pub fn is_uncertain_base_qual(q: u8) -> bool {
    q <= BASE_QUAL_MINIMUM
}

/// Java Math.round: floor(x + 0.5), with NaN -> 0.
pub fn jround(x: f64) -> i64 {
    if x.is_nan() {
        0
    } else {
        (x + 0.5).floor() as i64
    }
}

pub const MAX_SOFT_CLIP_LOW_QUAL_COUNT: i32 = 5;
pub const MIN_SOFT_CLIP_HIGH_QUAL_PERC: f64 = 0.75;
pub const MIN_SOFT_CLIP_MIN_BASE_QUAL: i32 = 25;
pub const DEFAULT_DOUBLE_JITTER_REPEAT_COUNT: i32 = 11;
pub const MSI_DOUBLE_JITTER_REPEAT_COUNT: i32 = 21;
pub const MSI_JITTER_MAX_REPEAT_CHANGE: i32 = 5;
pub const MSI_JITTER_DEFAULT_ERROR_RATE: f64 = 0.0001;
pub const MSI_JITTER_MIN_TRINUC_ERROR_RATE: f64 = 0.04;
pub const MSI_JITTER_NOISE_RATE: f64 = 0.00025;
pub const MSI_JITTER_HARD_FILTER_NOISE_RATE: f64 = 0.05;
pub const MSI_JITTER_MIN_RATIO: f64 = 2.0;
pub const MSI_JITTER_MIN_RATIO_HOTSPOT: f64 = 1.5;
pub const MSI_JITTER_RATIO_HOTSPOT_VAF: f64 = 0.03;
pub const JITTER_QUAL_BOOST_MAX_PERC: f64 = 1.3;
pub const INDEL_UNCERTAIN_BASE_REPEAT_MIN: i32 = 3;
pub const REQUIRED_UNIQUE_FRAG_COORDS_2: i32 = 3;
pub const EVIDENCE_MIN_MAP_QUAL: i32 = 1;
pub const HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY: i32 = 10;
pub const DEFAULT_BASE_QUAL_FIXED_PENALTY: i32 = 12;
pub const READ_EDGE_PENALTY_0: i32 = 15;
pub const READ_EDGE_PENALTY_1: i32 = 5;
pub const DEFAULT_HARD_MIN_TUMOR_QUAL: i32 = 50;
pub const DEFAULT_HARD_MIN_TUMOR_VAF: f64 = 0.002;

/// Java Math.max / Math.min for doubles: NaN propagates.
pub fn jmax(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a > b {
        a
    } else {
        b
    }
}
pub fn jmin(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a < b {
        a
    } else {
        b
    }
}

// --- sequencing technology (SageConfig.SEQUENCING_TYPE, a static in Java too) ---

use std::sync::atomic::{AtomicU8, Ordering};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequencingType {
    Illumina = 0,
    Ultima = 1,
    Sbx = 2,
}

static SEQUENCING_TYPE: AtomicU8 = AtomicU8::new(0);

pub fn set_sequencing_type(t: SequencingType) {
    SEQUENCING_TYPE.store(t as u8, Ordering::Relaxed);
}

#[cfg(test)]
thread_local! {
    /// Tests run in parallel threads, so a test selects its sequencing type for its own thread only.
    static TEST_SEQUENCING_TYPE: std::cell::Cell<Option<u8>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
pub fn set_test_sequencing_type(t: Option<SequencingType>) {
    TEST_SEQUENCING_TYPE.with(|c| c.set(t.map(|t| t as u8)));
}

fn sequencing_type() -> u8 {
    #[cfg(test)]
    if let Some(t) = TEST_SEQUENCING_TYPE.with(|c| c.get()) {
        return t;
    }
    SEQUENCING_TYPE.load(Ordering::Relaxed)
}

pub fn is_illumina() -> bool {
    sequencing_type() == SequencingType::Illumina as u8
}
pub fn is_sbx() -> bool {
    sequencing_type() == SequencingType::Sbx as u8
}
pub fn is_ultima() -> bool {
    sequencing_type() == SequencingType::Ultima as u8
}

// --- append mode (SageConfig.AppendMode, a static in Java too) ---

static APPEND_MODE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

pub fn set_append_mode(on: bool) {
    APPEND_MODE.store(on, Ordering::Relaxed);
}

#[cfg(test)]
thread_local! {
    static TEST_APPEND_MODE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
pub fn set_test_append_mode(on: Option<bool>) {
    TEST_APPEND_MODE.with(|c| c.set(on));
}

pub fn is_append_mode() -> bool {
    #[cfg(test)]
    if let Some(on) = TEST_APPEND_MODE.with(|c| c.get()) {
        return on;
    }
    APPEND_MODE.load(Ordering::Relaxed)
}

pub const NON_ILLUMINA_MAX_READ_LENGTH: i32 = 800;
pub const MATCHING_BASE_QUALITY_SBX: u8 = 10;
// SbxBamUtils: quals assigned by REDUX
const SBX_MEDIUM_QUAL_LOWER: i8 = 27;
const SBX_MEDIUM_QUAL_UPPER: i8 = 29;

/// The base quality below which a mismatch against the read context is tolerated.
pub fn matching_base_quality() -> u8 {
    if is_sbx() {
        MATCHING_BASE_QUALITY_SBX
    } else {
        MATCHING_BASE_QUALITY
    }
}

/// BaseQualAdjustment.isMediumBaseQual: only SBX has a medium band (simplex bases)
pub fn is_medium_base_qual(q: i8) -> bool {
    is_sbx() && q >= SBX_MEDIUM_QUAL_LOWER && q <= SBX_MEDIUM_QUAL_UPPER
}

/// BaseQualAdjustment.isHighBaseQual
pub fn is_high_base_qual_byte(q: i8) -> bool {
    if is_sbx() {
        q > SBX_MEDIUM_QUAL_UPPER
    } else if is_ultima() {
        q >= crate::ultima::ULTIMA_MAX_QUAL
    } else {
        q >= ILLUMINA_HIGH_BASE_QUAL as i8
    }
}

/// SageConstants.VIS_VARIANT_BUFFER: the region read around a specific or visualised variant
pub const VIS_VARIANT_BUFFER: i32 = 200;
