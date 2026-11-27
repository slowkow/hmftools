//! CiderConstants and the defaults of CiderParams.

/// for partially rearranged VDJs, the sequence kept is 60 bases before the J anchor or after the V anchor
pub const PARTIAL_VDJ_UNANCHORED_LENGTH_BASES: i32 = 60;

// BLOSUM search
pub const MAX_BLOSUM_DIFF_PER_AA: i32 = 3;
pub const BLOSUM_SIMILARITY_SCORE_CONSTANT: i32 = 6;
pub const BLOSUM_UNKNOWN_BASE_PENALTY: i32 = 2;

pub const CANDIDATE_MIN_PARTIAL_ANCHOR_AA_LENGTH: i32 = 10;
pub const VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH: i32 = 3;

/// the bases a read must overlap the CDR3 (past the anchor) by to be used
pub const MIN_READ_CDR3_OVERLAP: i32 = 20;

/// the bases a read must overlap a layout by to be added
pub const LAYOUT_MIN_READ_OVERLAP_BASES: i32 = 20;
pub const MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES: i32 = 20;

/// mismatches allowed between the consensus and a mate read extending it
pub const LAYOUT_MATE_EXTEND_MISMATCHES_MAX: i32 = 1;

pub const VJ_JOIN_HASH_WORD_SIZE: i32 = 8;
pub const MIN_ANCHOR_LENGTH_BASES: i32 = 3;
pub const LAYOUT_MIN_SUPPORT_TO_SEAL_NODE: i32 = 5;

/// the fraction of reads that must support a base of a layout
pub const MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION: f64 = 0.02;
pub const MIN_POLY_G_TRIM_COUNT: i32 = 2;
pub const POLY_G_TRIM_EXTRA_BASE_COUNT: i32 = 5;

pub const MIN_NON_SPLIT_READ_STRADDLE_LENGTH: i32 = 30;

/// sequences per BWA-MEM call
pub const BWAMEM_BATCH_SIZE: usize = 20000;

pub const ANNOTATION_ALIGN_SCORE_MIN: i32 = 15;
pub const ANNOTATION_VDJ_FLANK_BASES: i32 = 300;
pub const ANNOTATION_VJ_IDENTITY_MIN: f64 = 0.9;
pub const ANNOTATION_MATCH_REF_IDENTITY: f64 = 0.95;
pub const ANNOTATION_ANCHOR_BASE_OVERLAP: i32 = 5;
pub const ANNOTATION_V_J_OVERLAP_MAX: i32 = 30;

// aligner scores
pub const ALIGNER_MATCH_SCORE: i32 = 1;
pub const ALIGNER_MISMATCH_SCORE: i32 = -4;
pub const ALIGNER_GAP_OPENING_SCORE: i32 = -5;
pub const ALIGNER_GAP_EXTEND_SCORE: i32 = -2;
pub const ALIGNER_WORD_SIZE: i32 = 8;

/// MATCHES_REF_KNOWN_CDR3_AA: "CTXGPKXELRT.*"
pub const MATCHES_REF_KNOWN_CDR3_AA_PREFIX: &str = "CTXGPKXELRT";

// VdjAnnotator
pub const CDR3_FILTER_AA_MIN_LENGTH: i32 = 5;
pub const CDR3_FILTER_AA_MAX_LENGTH: i32 = 40;
pub const MAX_NONSPLIT_READS: i32 = 2;

// LayoutForest
pub const LOW_QUAL_BASE_FACTOR: f64 = 0.1;
pub const REASSIGN_READ_ROUNDS: i32 = 6;
pub const ADD_READ_INFINITE_LOOP_CHECK: i32 = 10_000;

// CiderParams defaults
#[allow(dead_code)]
pub const DEFAULT_MAX_FRAGMENT_LENGTH: i32 = 1000;
#[allow(dead_code)]
pub const DEFAULT_MAX_READ_COUNT_PER_GENE: i32 = 100_000;
#[allow(dead_code)]
pub const DEFAULT_MAX_LOW_QUAL_BASES_FRACTION: f64 = 0.1;
#[allow(dead_code)]
pub const DEFAULT_MIN_BASE_QUALITY: i32 = 25;
