//! LilacConstants that do not change with the config (those that do are in `config::Params`).

pub const LOW_BASE_TRIM_PERC: f64 = 0.35;
pub const MAX_LOW_BASE_PERC: f64 = 0.5;

pub const BASE_QUAL_PERCENTILE: f64 = 0.25;
/// BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD: quals at and below are low
pub const LOW_BASE_QUAL_THRESHOLD: u8 = 26;
/// BaseQualAdjustment.BASE_QUAL_MINIMUM: quals at and below are uncertain
pub const BASE_QUAL_MINIMUM: u8 = 1;
/// Nucleotide.MISSING_BASE_QUAL
pub const MISSING_BASE_QUAL: i8 = -1;

pub const COMMON_ALLELES_FREQ_CUTOFF: f64 = 0.005;
pub const RARE_ALLELES_FREQ_CUTOFF: f64 = 0.0001;

pub const MIN_HI_CONF_UNIQUE_GROUP_COVERAGE: f64 = 0.02;
pub const MIN_HI_CONF_UNIQUE_PROTEIN_COVERAGE: f64 = 0.01;
pub const MIN_LOW_CONF_GROUP_UNIQUE_COVERAGE: f64 = 0.001;
pub const MIN_LOW_CONF_GROUP_TOTAL_COVERAGE: f64 = 0.035;
pub const MIN_WILDCARD_FRAGMENTS: usize = 2;

pub const MIN_POPULATION_FREQUENCY: f64 = 0.0001;
pub const SOLUTION_COMPLEXITY_PENALTY_WEIGHT: f64 = 0.002;
pub const RECOVERY_SCORE_PENALTY: f64 = 0.0;

pub const WARN_LOW_BASE_QUAL_THRESHOLD: f64 = 25.0;
pub const LOG_UNMATCHED_HAPLOTYPE_SUPPORT: i32 = 3;
pub const WARN_LOW_COVERAGE_DEPTH: f64 = 10.0;

pub const SPLICE_VARIANT_BUFFER: i32 = 5;
pub const COMPLEX_PERMS_THRESHOLD: i64 = 100000;

/// BamRecordReader.MAX_DISTANCE
pub const MAX_DISTANCE: i32 = 1000;

/// SBX medium base quals (SbxBamUtils)
pub const SBX_MEDIUM_QUAL_LOWER: u8 = 27;
pub const SBX_MEDIUM_QUAL_UPPER: u8 = 29;

pub const LILAC_FILE_CANDIDATE_COVERAGE: &str = "candidates.coverage.tsv";
pub const LILAC_FILE_CANDIDATE_FRAGS: &str = "candidates.fragments.tsv";
pub const LILAC_FILE_FRAGMENTS: &str = "fragments.tsv";
pub const LILAC_FILE_CANDIDATE_AA: &str = "candidates.aminoacids.txt";
pub const LILAC_FILE_CANDIDATE_NUC: &str = "candidates.nucleotides.txt";
pub const LILAC_FILE_HLA_Y_COVERAGE: &str = "HLA_Y.coverage.tsv";
pub const LILAC_FILE_HLA_Y_FRAGMENTS: &str = "HLA_Y.fragments.tsv";
