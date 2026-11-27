//! hmf-common `metrics`: the BamTools metrics files, written by BamTools and read by PURPLE, ORANGE and
//! qsee. BamMetricSummary (`.bam_metric.summary.tsv`), BamMetricCoverage (`.coverage.tsv`),
//! BamMetricFragmentLength (`.frag_length.tsv`), the flag counts file (`.flag_counts.tsv`, written by
//! BamTools' MetricsWriter, read by BamFlagStats) and GeneDepthFile (`.gene_coverage.tsv`,
//! `.exon_coverage.tsv`). Numbers are formatted as Java's String.format does it (`fmt::java_fixed`).

use anyhow::{anyhow, bail, Context, Result};

use crate::fmt::java_fixed;
use crate::java_text::{self, parse_double, parse_int, parse_long};

pub const BAM_METRICS_FILE_ID: &str = ".bam_metric";
pub const COUNT_COLUMN: &str = "Count";

/// FileWriterUtils.checkAddDirSeparator
pub fn check_add_dir_separator(path: &str) -> String {
    if path.is_empty() || path.ends_with('/') { path.to_string() } else { format!("{path}/") }
}

fn metrics_filename(base_path: &str, sample_id: &str, extension: &str) -> String {
    format!("{}{sample_id}{BAM_METRICS_FILE_ID}{extension}", check_add_dir_separator(base_path))
}

/// Files.write(path, lines): each line ended by a newline
fn write_lines(filename: &str, lines: &[String]) -> Result<()> {
    let mut s = String::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    std::fs::write(filename, s).with_context(|| format!("writing {filename}"))
}

fn read_all_lines(filename: &str) -> Result<Vec<String>> {
    java_text::read_lines(filename).with_context(|| format!("reading {filename}"))
}

fn column(f: &std::collections::HashMap<String, usize>, name: &str) -> Result<usize> {
    f.get(name).copied().ok_or_else(|| anyhow!("missing column {name}"))
}

fn value<'a>(values: &[&'a str], i: usize) -> Result<&'a str> {
    values.get(i).copied().ok_or_else(|| anyhow!("missing value in column {}", i + 1))
}

/// ValueFrequency
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ValueFrequency {
    pub value: i32,
    pub count: i64,
}

impl ValueFrequency {
    pub fn new(value: i32, count: i64) -> Self {
        ValueFrequency { value, count }
    }
}

// ---- summary ----

pub const TOTAL_REGION_COLUMN: &str = "TotalRegionBases";
pub const TOTAL_READS_COLUMN: &str = "TotalReads";
pub const DUPLICATES_COLUMN: &str = "DuplicateReads";
pub const DUAL_STRAND_COLUMN: &str = "DualStrandReads";
pub const OFF_TARGET_COLUMN: &str = "OffTargetReads";
pub const MEAN_COVERAGE_COLUMN: &str = "MeanCoverage";
pub const SD_COVERAGE_COLUMN: &str = "StdDevCoverage";
pub const MEDIAN_COVERAGE_COLUMN: &str = "MedianCoverage";
pub const MAD_COVERAGE_COLUMN: &str = "MadCoverage";
pub const LOW_MAPQ_COLUMN: &str = "LowMapQualPercent";
pub const DUPLICATE_PCT_COLUMN: &str = "DuplicatePercent";
pub const UNMAPPED_COLUMN: &str = "UnmappedPercent";
pub const LOW_BASEQ_COLUMN: &str = "LowBaseQualPercent";
pub const OVERLAP_READ_COLUMN: &str = "OverlappingReadPercent";
pub const CAPPED_COVERAGE_COLUMN: &str = "CappedCoverage";
const UNPAIRED_COLUMN: &str = "UnpairedPercent"; // deprecated in v1.5
const DEPTH_COVERAGE: &str = "DepthCoverage";

/// BamMetricSummary
#[derive(Clone, Debug, PartialEq)]
pub struct BamMetricSummary {
    pub total_region_bases: i64,
    pub total_reads: i64,
    pub duplicate_reads: i64,
    pub dual_strand_reads: i64,
    pub off_target_reads: i64,
    pub mean_coverage: f64,
    pub sd_coverage: f64,
    pub median_coverage: i32,
    pub mad_coverage: i32,
    pub low_map_qual_percent: f64,
    pub duplicate_percent: f64,
    pub unmapped_percent: f64,
    pub low_base_qual_percent: f64,
    pub overlapping_read_percent: f64,
    pub capped_coverage_percent: f64,
    pub coverage_levels: Vec<i32>,
    pub coverage_percents: Vec<f64>,
}

impl BamMetricSummary {
    pub fn generate_filename(base_path: &str, sample_id: &str) -> String {
        metrics_filename(base_path, sample_id, ".summary.tsv")
    }

    pub fn total_filtered_percent(&self) -> f64 {
        self.low_map_qual_percent
            + self.duplicate_percent
            + self.unmapped_percent
            + self.low_base_qual_percent
            + self.overlapping_read_percent
            + self.capped_coverage_percent
    }

    /// getCoveragePercent: the percentage for that level, if the file has it
    pub fn get_coverage_percent(&self, level: i32) -> Option<f64> {
        self.coverage_levels.iter().position(|l| *l == level).map(|i| self.coverage_percents[i])
    }

    /// coveragePercent: 0 for a level the file does not have
    pub fn coverage_percent(&self, level: i32) -> f64 {
        self.get_coverage_percent(level).unwrap_or(0.0)
    }

    /// the header and the values line of write()
    pub fn lines(&self) -> Vec<String> {
        let mut header: Vec<String> = [
            TOTAL_REGION_COLUMN,
            TOTAL_READS_COLUMN,
            DUPLICATES_COLUMN,
            DUAL_STRAND_COLUMN,
            OFF_TARGET_COLUMN,
            MEAN_COVERAGE_COLUMN,
            SD_COVERAGE_COLUMN,
            MEDIAN_COVERAGE_COLUMN,
            MAD_COVERAGE_COLUMN,
            LOW_MAPQ_COLUMN,
            DUPLICATE_PCT_COLUMN,
            UNMAPPED_COLUMN,
            LOW_BASEQ_COLUMN,
            OVERLAP_READ_COLUMN,
            CAPPED_COVERAGE_COLUMN,
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        for l in &self.coverage_levels {
            header.push(format!("{DEPTH_COVERAGE}_{l}"));
        }
        let mut values: Vec<String> = vec![
            self.total_region_bases.to_string(),
            self.total_reads.to_string(),
            self.duplicate_reads.to_string(),
            self.dual_strand_reads.to_string(),
            self.off_target_reads.to_string(),
            java_fixed(self.mean_coverage, 3),
            java_fixed(self.sd_coverage, 3),
            self.median_coverage.to_string(),
            self.mad_coverage.to_string(),
            java_fixed(self.low_map_qual_percent, 5),
            java_fixed(self.duplicate_percent, 5),
            java_fixed(self.unmapped_percent, 5),
            java_fixed(self.low_base_qual_percent, 5),
            java_fixed(self.overlapping_read_percent, 5),
            java_fixed(self.capped_coverage_percent, 5),
        ];
        for p in &self.coverage_percents {
            values.push(java_fixed(*p, 5));
        }
        vec![header.join("\t"), values.join("\t")]
    }

    pub fn write(&self, filename: &str) -> Result<()> {
        write_lines(filename, &self.lines())
    }

    /// BamMetricSummary.read: the DepthCoverage_<n> columns in file order; UnpairedPercent for an older
    /// file without UnmappedPercent; OffTargetReads 0 when absent
    pub fn read(filename: &str) -> Result<Self> {
        let lines = read_all_lines(filename)?;
        Self::from_lines(&lines).with_context(|| format!("parsing {filename}"))
    }

    pub fn from_lines(lines: &[String]) -> Result<Self> {
        if lines.len() < 2 {
            bail!("expected a header and a values line");
        }
        let f = java_text::fields_index_map(&lines[0], '\t');
        let values = java_text::split_all(&lines[1], '\t');
        let columns = java_text::split_all(&lines[0], '\t');
        let mut levels = Vec::new();
        let mut percents = Vec::new();
        for (i, c) in columns.iter().enumerate() {
            if c.starts_with(DEPTH_COVERAGE) {
                let level = c.splitn(2, '_').nth(1).and_then(parse_int).ok_or_else(|| anyhow!("invalid column {c}"))?;
                levels.push(level);
                percents.push(parse_double(value(&values, i)?).ok_or_else(|| anyhow!("invalid value for {c}"))?);
            }
        }
        let unmapped_i = match f.get(UNMAPPED_COLUMN) {
            Some(i) => *i,
            None => column(&f, UNPAIRED_COLUMN)?,
        };
        let long = |name: &str| -> Result<i64> { parse_long(value(&values, column(&f, name)?)?).ok_or_else(|| anyhow!("invalid {name}")) };
        let dbl_at = |i: usize, name: &str| -> Result<f64> { parse_double(value(&values, i)?).ok_or_else(|| anyhow!("invalid {name}")) };
        let dbl = |name: &str| -> Result<f64> { dbl_at(column(&f, name)?, name) };
        Ok(BamMetricSummary {
            total_region_bases: long(TOTAL_REGION_COLUMN)?,
            total_reads: long(TOTAL_READS_COLUMN)?,
            duplicate_reads: long(DUPLICATES_COLUMN)?,
            dual_strand_reads: long(DUAL_STRAND_COLUMN)?,
            off_target_reads: if f.contains_key(OFF_TARGET_COLUMN) { long(OFF_TARGET_COLUMN)? } else { 0 },
            mean_coverage: dbl(MEAN_COVERAGE_COLUMN)?,
            sd_coverage: dbl(SD_COVERAGE_COLUMN)?,
            // (int)Double.parseDouble
            median_coverage: java_double_to_int(dbl(MEDIAN_COVERAGE_COLUMN)?),
            mad_coverage: java_double_to_int(dbl(MAD_COVERAGE_COLUMN)?),
            low_map_qual_percent: dbl(LOW_MAPQ_COLUMN)?,
            duplicate_percent: dbl(DUPLICATE_PCT_COLUMN)?,
            unmapped_percent: dbl_at(unmapped_i, UNMAPPED_COLUMN)?,
            low_base_qual_percent: dbl(LOW_BASEQ_COLUMN)?,
            overlapping_read_percent: dbl(OVERLAP_READ_COLUMN)?,
            capped_coverage_percent: dbl(CAPPED_COVERAGE_COLUMN)?,
            coverage_levels: levels,
            coverage_percents: percents,
        })
    }
}

/// a Java (int) cast of a double: towards zero, saturating, NaN to 0 (as Rust's `as`)
fn java_double_to_int(d: f64) -> i32 {
    d as i32
}

// ---- coverage and fragment length distributions ----

pub const COVERAGE_COLUMN: &str = "Coverage";
pub const FRAG_LENGTH_COLUMN: &str = "FragmentLength";

fn write_value_frequencies(filename: &str, value_column: &str, rows: &[ValueFrequency]) -> Result<()> {
    let mut lines = Vec::with_capacity(rows.len() + 1);
    lines.push(format!("{value_column}\t{COUNT_COLUMN}"));
    for r in rows {
        lines.push(format!("{}\t{}", r.value, r.count));
    }
    write_lines(filename, &lines)
}

fn read_value_frequencies(filename: &str, value_column: &str) -> Result<Vec<ValueFrequency>> {
    let lines = read_all_lines(filename)?;
    let header = lines.first().ok_or_else(|| anyhow!("{filename}: empty file"))?;
    let f = java_text::fields_index_map(header, '\t');
    let vi = column(&f, value_column)?;
    let ci = column(&f, COUNT_COLUMN)?;
    let mut rows = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let v = java_text::split_all(line, '\t');
        let value_ = parse_int(value(&v, vi)?).ok_or_else(|| anyhow!("{filename}: invalid value in line: {line}"))?;
        let count = parse_long(value(&v, ci)?).ok_or_else(|| anyhow!("{filename}: invalid count in line: {line}"))?;
        rows.push(ValueFrequency::new(value_, count));
    }
    Ok(rows)
}

/// BamMetricCoverage: bases per coverage level
pub struct BamMetricCoverage;

impl BamMetricCoverage {
    pub fn generate_filename(base_path: &str, sample_id: &str) -> String {
        metrics_filename(base_path, sample_id, ".coverage.tsv")
    }
    pub fn write(filename: &str, coverage_levels: &[ValueFrequency]) -> Result<()> {
        write_value_frequencies(filename, COVERAGE_COLUMN, coverage_levels)
    }
    pub fn read(filename: &str) -> Result<Vec<ValueFrequency>> {
        read_value_frequencies(filename, COVERAGE_COLUMN)
    }
}

/// BamMetricFragmentLength: fragments per length
pub struct BamMetricFragmentLength;

impl BamMetricFragmentLength {
    pub fn generate_filename(base_path: &str, sample_id: &str) -> String {
        metrics_filename(base_path, sample_id, ".frag_length.tsv")
    }
    pub fn write(filename: &str, lengths: &[ValueFrequency]) -> Result<()> {
        write_value_frequencies(filename, FRAG_LENGTH_COLUMN, lengths)
    }
    pub fn read(filename: &str) -> Result<Vec<ValueFrequency>> {
        read_value_frequencies(filename, FRAG_LENGTH_COLUMN)
    }
}

// ---- flag counts ----

pub const FLAGSTAT_TOTAL: &str = "in total (QC-passed reads + QC-failed reads)";
pub const FLAGSTAT_TOTAL_SHORT: &str = "total";
pub const FLAGSTAT_PRIMARY: &str = "primary";
pub const FLAGSTAT_SECONDARY: &str = "secondary";
pub const FLAGSTAT_SUPPLEMENTARY: &str = "supplementary";
pub const FLAGSTAT_DUPLICATE: &str = "duplicates";
pub const FLAGSTAT_PRIMARY_DUPLICATE: &str = "primary duplicates";
pub const FLAGSTAT_MAPPED: &str = "mapped";
pub const FLAGSTAT_PRIMARY_MAPPED: &str = "primary mapped";
pub const FLAGSTAT_PAIRED: &str = "paired in sequencing";
pub const FLAGSTAT_READ1: &str = "read1";
pub const FLAGSTAT_READ2: &str = "read2";
pub const FLAGSTAT_PROPER_PAIR: &str = "properly paired";
pub const FLAGSTAT_BOTH_MAPPED: &str = "with itself and mate mapped";
pub const FLAGSTAT_SINGLE: &str = "singletons";
pub const FLAGSTAT_MATE_DIFF_CHR: &str = "with mate mapped to a different chr";
pub const FLAGSTAT_MATE_DIFF_CHR_MAPQ_5: &str = "with mate mapped to a different chr (mapQ>=5)";

/// BamTools FlagStatType, in its order
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub enum FlagStatType {
    Total,
    Primary,
    Secondary,
    Supplementary,
    Duplicate,
    PrimaryDuplicate,
    Mapped,
    PrimaryMapped,
    Paired,
    Read1,
    Read2,
    ProperlyPaired,
    PairMapped,
    Singleton,
    InterChrPairMapped,
    InterChrPairMapQualGe5,
}

pub const FLAG_STAT_TYPES: usize = 16;

/// BamTools FlagQCStats: QC-passed and QC-failed counts
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FlagQcStats {
    pub passed: i64,
    pub failed: i64,
}

impl FlagQcStats {
    /// "%d + %d"
    pub fn as_string(&self) -> String {
        format!("{} + {}", self.passed, self.failed)
    }
}

/// FlagQCStats.flagStatsPercentages: "(p% : f%)", each computed in float as 100.0f * n / d, "N/A" for a
/// zero denominator
pub fn flag_stats_percentages(numerator: &FlagQcStats, denominator: &FlagQcStats) -> String {
    let pct = |n: i64, d: i64| -> String {
        if d == 0 {
            "N/A".to_string()
        } else {
            let v: f32 = 100.0f32 * (n as f32) / (d as f32);
            format!("{}%", java_fixed(v as f64, 2))
        }
    };
    format!("({} : {})", pct(numerator.passed, denominator.passed), pct(numerator.failed, denominator.failed))
}

/// The flag counts file, as BamTools' MetricsWriter.writeFlagCounts writes it (in the layout of samtools
/// flagstat)
pub struct BamFlagCounts;

impl BamFlagCounts {
    pub fn generate_filename(base_path: &str, sample_id: &str) -> String {
        metrics_filename(base_path, sample_id, BamFlagStats::FILE_EXTENSION)
    }

    pub fn lines(stats: &[FlagQcStats; FLAG_STAT_TYPES]) -> Vec<String> {
        use FlagStatType as T;
        let s = |t: T| stats[t as usize].as_string();
        let pct = |n: T, d: T| flag_stats_percentages(&stats[n as usize], &stats[d as usize]);
        vec![
            format!("{} {FLAGSTAT_TOTAL}", s(T::Total)),
            format!("{} {FLAGSTAT_PRIMARY}", s(T::Primary)),
            format!("{} {FLAGSTAT_SECONDARY}", s(T::Secondary)),
            format!("{} {FLAGSTAT_SUPPLEMENTARY}", s(T::Supplementary)),
            format!("{} {FLAGSTAT_DUPLICATE}", s(T::Duplicate)),
            format!("{} {FLAGSTAT_PRIMARY_DUPLICATE}", s(T::PrimaryDuplicate)),
            format!("{} {FLAGSTAT_MAPPED} {}", s(T::Mapped), pct(T::Mapped, T::Total)),
            format!("{} {FLAGSTAT_PRIMARY_MAPPED} {}", s(T::PrimaryMapped), pct(T::PrimaryMapped, T::Primary)),
            format!("{} {FLAGSTAT_PAIRED}", s(T::Paired)),
            format!("{} {FLAGSTAT_READ1}", s(T::Read1)),
            format!("{} {FLAGSTAT_READ2}", s(T::Read2)),
            format!("{} {FLAGSTAT_PROPER_PAIR} {}", s(T::ProperlyPaired), pct(T::ProperlyPaired, T::Paired)),
            format!("{} {FLAGSTAT_BOTH_MAPPED}", s(T::PairMapped)),
            format!("{} {FLAGSTAT_SINGLE} {}", s(T::Singleton), pct(T::Singleton, T::Paired)),
            format!("{} {FLAGSTAT_MATE_DIFF_CHR}", s(T::InterChrPairMapped)),
            format!("{} {FLAGSTAT_MATE_DIFF_CHR_MAPQ_5}", s(T::InterChrPairMapQualGe5)),
        ]
    }

    pub fn write(filename: &str, stats: &[FlagQcStats; FLAG_STAT_TYPES]) -> Result<()> {
        write_lines(filename, &Self::lines(stats))
    }
}

/// BamFlagStats: the proportions ORANGE and qsee take from the flag counts file
#[derive(Clone, Debug, PartialEq)]
pub struct BamFlagStats {
    pub unique_read_count: i64,
    pub secondary_count: i64,
    pub supplementary_count: i64,
    pub duplicate_proportion: f64,
    pub mapped_proportion: f64,
    pub paired_in_sequencing_proportion: f64,
    pub properly_paired_proportion: f64,
    pub with_itself_and_mate_mapped_proportion: f64,
    pub singleton_proportion: f64,
}

impl BamFlagStats {
    pub const FILE_EXTENSION: &'static str = ".flag_counts.tsv";

    pub fn generate_filename(base_path: &str, sample_id: &str) -> String {
        metrics_filename(base_path, sample_id, Self::FILE_EXTENSION)
    }

    /// BamFlagStats.read: the QC-passed count of the first line containing each label
    pub fn read(filename: &str) -> Result<Self> {
        let lines = read_all_lines(filename)?;
        Self::from_lines(&lines).with_context(|| format!("parsing {filename}"))
    }

    pub fn from_lines(lines: &[String]) -> Result<Self> {
        let by = |sub: &str| -> Option<String> { lines.iter().find(|l| l.contains(sub)).map(|l| java_text::split(l, ' ').first().copied().unwrap_or("").to_string()) };
        let (Some(total), Some(secondary), Some(supplementary), Some(duplicates), Some(mapped), Some(paired), Some(proper), Some(both), Some(single)) = (
            by(FLAGSTAT_TOTAL_SHORT),
            by(FLAGSTAT_SECONDARY),
            by(FLAGSTAT_SUPPLEMENTARY),
            by(FLAGSTAT_DUPLICATE),
            by("mapped ("),
            by(FLAGSTAT_PAIRED),
            by(FLAGSTAT_PROPER_PAIR),
            by(FLAGSTAT_BOTH_MAPPED),
            by(FLAGSTAT_SINGLE),
        ) else {
            bail!("Unable to parse flagstat file correctly");
        };
        let num = |s: &str| parse_long(s).ok_or_else(|| anyhow!("NumberFormatException: {s}"));
        let total_n = num(&total)?;
        let secondary_n = num(&secondary)?;
        let supplementary_n = num(&supplementary)?;
        let unique = total_n - secondary_n - supplementary_n;
        let div = |a: &str, b: i64| -> Result<f64> { round_to_six_decimals(num(a)? as f64 / b as f64) };
        Ok(BamFlagStats {
            unique_read_count: unique,
            secondary_count: secondary_n,
            supplementary_count: supplementary_n,
            duplicate_proportion: div(&duplicates, total_n)?,
            mapped_proportion: div(&mapped, total_n)?,
            paired_in_sequencing_proportion: div(&paired, unique)?,
            properly_paired_proportion: div(&proper, unique)?,
            with_itself_and_mate_mapped_proportion: div(&both, total_n)?,
            singleton_proportion: div(&single, total_n)?,
        })
    }
}

/// BigDecimal.valueOf(value).setScale(6, HALF_UP).doubleValue(): Double.toString's digits rounded half up
/// (an infinite value makes BigDecimal throw)
fn round_to_six_decimals(v: f64) -> Result<f64> {
    if v.is_nan() {
        return Ok(f64::NAN);
    }
    if v.is_infinite() {
        bail!("NumberFormatException: Infinite or NaN");
    }
    Ok(java_fixed(v, 6).parse::<f64>().unwrap())
}

// ---- gene and exon depth ----

pub const FLD_GENE_NAME: &str = "GeneName";
pub const FLD_CHROMOSOME: &str = "Chromosome";
pub const FLD_POS_START: &str = "PosStart";
pub const FLD_POS_END: &str = "PosEnd";
pub const FLD_EXON_RANK: &str = "ExonRank";
pub const FLD_MEDIAN_DEPTH: &str = "MedianDepth";
pub const FLD_MEAN_DEPTH: &str = "MeanDepth";
const FLD_DEPTH_RANGE_PREFIX: &str = "DR_";

/// DepthRange
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct DepthRange {
    pub lower: i32,
    pub upper: i32,
}

impl DepthRange {
    pub fn midpoint(&self) -> i32 {
        (self.lower + self.upper + 1) / 2
    }
}

/// GeneDepthFile.depthRangeFieldName: DR_<n> or DR_<lower>_<upper>
pub fn depth_range_field_name(r: &DepthRange) -> String {
    if r.lower == r.upper { format!("{FLD_DEPTH_RANGE_PREFIX}{}", r.lower) } else { format!("{FLD_DEPTH_RANGE_PREFIX}{}_{}", r.lower, r.upper) }
}

/// GeneDepth; the depth counts in column order
#[derive(Clone, Debug, PartialEq)]
pub struct GeneDepth {
    pub gene: String,
    pub chromosome: String,
    pub pos_start: i32,
    pub pos_end: i32,
    pub depth_counts: Vec<(DepthRange, i32)>,
}

/// ExonDepth
#[derive(Clone, Debug, PartialEq)]
pub struct ExonDepth {
    pub depth: GeneDepth,
    pub exon_rank: i32,
    pub median_depth: f64,
    pub mean_depth: f64,
}

fn coverage_header(ranges: &[(DepthRange, i32)], for_exons: bool) -> String {
    let mut h: Vec<String> = vec![FLD_GENE_NAME.into(), FLD_CHROMOSOME.into(), FLD_POS_START.into(), FLD_POS_END.into()];
    if for_exons {
        h.extend([FLD_EXON_RANK.to_string(), FLD_MEDIAN_DEPTH.to_string(), FLD_MEAN_DEPTH.to_string()]);
    }
    h.extend(ranges.iter().map(|(r, _)| depth_range_field_name(r)));
    h.join("\t")
}

fn depth_line(d: &GeneDepth, exon: Option<(i32, f64, f64)>) -> String {
    let mut s = format!("{}\t{}\t{}\t{}", d.gene, d.chromosome, d.pos_start, d.pos_end);
    if let Some((rank, median, mean)) = exon {
        s.push_str(&format!("\t{rank}\t{}\t{}", java_fixed(median, 0), java_fixed(mean, 2)));
    }
    for (_, c) in &d.depth_counts {
        s.push('\t');
        s.push_str(&c.to_string());
    }
    s
}

/// GeneDepthFile
pub struct GeneDepthFile;

impl GeneDepthFile {
    pub fn generate_gene_coverage_filename(base_path: &str, sample: &str) -> String {
        metrics_filename(base_path, sample, ".gene_coverage.tsv")
    }
    pub fn generate_exon_coverage_filename(base_path: &str, sample: &str) -> String {
        metrics_filename(base_path, sample, ".exon_coverage.tsv")
    }

    /// writeGeneCoverage: sorted by gene name (stable); an empty file for no genes; the depth columns of
    /// the first gene
    pub fn write_gene_coverage(filename: &str, depths: &mut [GeneDepth]) -> Result<()> {
        if depths.is_empty() {
            return write_lines(filename, &[]);
        }
        depths.sort_by(|a, b| java_text::compare(&a.gene, &b.gene));
        let mut lines = vec![coverage_header(&depths[0].depth_counts, false)];
        lines.extend(depths.iter().map(|d| depth_line(d, None)));
        write_lines(filename, &lines)
    }

    /// writeExonCoverage: sorted by gene name then exon rank (stable); median as %.0f, mean as %.2f
    pub fn write_exon_coverage(filename: &str, depths: &mut [ExonDepth]) -> Result<()> {
        if depths.is_empty() {
            return write_lines(filename, &[]);
        }
        depths.sort_by(|a, b| java_text::compare(&a.depth.gene, &b.depth.gene).then(a.exon_rank.cmp(&b.exon_rank)));
        let mut lines = vec![coverage_header(&depths[0].depth.depth_counts, true)];
        lines.extend(depths.iter().map(|d| depth_line(&d.depth, Some((d.exon_rank, d.median_depth, d.mean_depth)))));
        write_lines(filename, &lines)
    }

    pub fn read_gene_coverage(filename: &str) -> Result<Vec<GeneDepth>> {
        Ok(read_depths(filename, false)?.into_iter().map(|e| e.depth).collect())
    }

    pub fn read_exon_coverage(filename: &str) -> Result<Vec<ExonDepth>> {
        read_depths(filename, true)
    }
}

/// GeneDepthFile.readDepths: the depth columns are those from the first DR_ column to the end
fn read_depths(filename: &str, for_exons: bool) -> Result<Vec<ExonDepth>> {
    let lines = read_all_lines(filename)?;
    let header = lines.first().ok_or_else(|| anyhow!("{filename}: empty file"))?;
    let f = java_text::fields_index_map(header, '\t');
    let chr_i = column(&f, FLD_CHROMOSOME)?;
    let gene_i = column(&f, FLD_GENE_NAME)?;
    let start_i = column(&f, FLD_POS_START)?;
    let end_i = column(&f, FLD_POS_END)?;
    let exon_cols = if for_exons { Some((column(&f, FLD_EXON_RANK)?, column(&f, FLD_MEDIAN_DEPTH)?, column(&f, FLD_MEAN_DEPTH)?)) } else { None };
    let columns = java_text::split_all(header, '\t');
    // Java: findFirst().orElse(-1), and the depth ranges from there (none parse from -1: substring fails)
    let depth_start = columns.iter().position(|c| c.starts_with(FLD_DEPTH_RANGE_PREFIX)).ok_or_else(|| anyhow!("{filename}: no depth range columns"))?;
    let mut ranges = Vec::new();
    for c in &columns[depth_start..] {
        let bounds: Vec<&str> = java_text::split(&c[FLD_DEPTH_RANGE_PREFIX.len().min(c.len())..], '_');
        let lower = bounds.first().and_then(|s| parse_int(s)).ok_or_else(|| anyhow!("{filename}: invalid depth column {c}"))?;
        let upper = bounds.last().and_then(|s| parse_int(s)).ok_or_else(|| anyhow!("{filename}: invalid depth column {c}"))?;
        ranges.push(DepthRange { lower, upper });
    }
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let v = java_text::split_all(line, '\t');
        let bad = || anyhow!("{filename}: invalid line: {line}");
        let mut counts = Vec::with_capacity(ranges.len());
        for (i, r) in ranges.iter().enumerate() {
            counts.push((*r, parse_int(value(&v, depth_start + i)?).ok_or_else(bad)?));
        }
        let depth = GeneDepth {
            gene: value(&v, gene_i)?.to_string(),
            chromosome: value(&v, chr_i)?.to_string(),
            pos_start: parse_int(value(&v, start_i)?).ok_or_else(bad)?,
            pos_end: parse_int(value(&v, end_i)?).ok_or_else(bad)?,
            depth_counts: counts,
        };
        let (rank, median, mean) = match exon_cols {
            Some((ri, mi, ai)) => (
                parse_int(value(&v, ri)?).ok_or_else(bad)?,
                parse_double(value(&v, mi)?).ok_or_else(bad)?,
                parse_double(value(&v, ai)?).ok_or_else(bad)?,
            ),
            None => (0, 0.0, 0.0),
        };
        out.push(ExonDepth { depth, exon_rank: rank, median_depth: median, mean_depth: mean });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> String {
        let dir = std::env::temp_dir().join(format!("hmf_common_metrics_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name).to_str().unwrap().to_string()
    }

    #[test]
    fn summary_round_trip() {
        let s = BamMetricSummary {
            total_region_bases: 2_900_000_000,
            total_reads: 123,
            duplicate_reads: -2,
            dual_strand_reads: 0,
            off_target_reads: 7,
            mean_coverage: 30.12345,
            sd_coverage: 0.0005,
            median_coverage: 30,
            mad_coverage: 4,
            low_map_qual_percent: 0.123456,
            duplicate_percent: 0.000005,
            unmapped_percent: 0.0,
            low_base_qual_percent: 1.0,
            overlapping_read_percent: 0.5,
            capped_coverage_percent: 0.25,
            coverage_levels: vec![1, 5],
            coverage_percents: vec![0.999995, 0.5],
        };
        let lines = s.lines();
        assert_eq!(lines[0].split('\t').count(), 17);
        assert!(lines[0].ends_with("CappedCoverage\tDepthCoverage_1\tDepthCoverage_5"));
        // %.3f and %.5f round half up on the shortest decimal digits
        assert_eq!(lines[1], "2900000000\t123\t-2\t0\t7\t30.123\t0.001\t30\t4\t0.12346\t0.00001\t0.00000\t1.00000\t0.50000\t0.25000\t1.00000\t0.50000");
        let p = tmp("s.summary.tsv");
        s.write(&p).unwrap();
        let r = BamMetricSummary::read(&p).unwrap();
        assert_eq!(r.total_region_bases, 2_900_000_000);
        assert_eq!(r.coverage_levels, vec![1, 5]);
        assert_eq!(r.mean_coverage, 30.123);
        assert_eq!(r.coverage_percent(5), 0.5);
        assert_eq!(r.coverage_percent(30), 0.0);
        // an older file: UnpairedPercent, no OffTargetReads
        let old = vec![
            "TotalRegionBases\tTotalReads\tDuplicateReads\tDualStrandReads\tMeanCoverage\tStdDevCoverage\tMedianCoverage\tMadCoverage\tLowMapQualPercent\tDuplicatePercent\tUnpairedPercent\tLowBaseQualPercent\tOverlappingReadPercent\tCappedCoverage\tDepthCoverage_10".to_string(),
            "10\t2\t0\t0\t1.5\t0.5\t1.0\t0\t0.1\t0.2\t0.3\t0.4\t0.5\t0.6\t0.7".to_string(),
        ];
        let o = BamMetricSummary::from_lines(&old).unwrap();
        assert_eq!((o.off_target_reads, o.unmapped_percent, o.median_coverage), (0, 0.3, 1));
        assert!((o.total_filtered_percent() - 2.1).abs() < 1e-12);
    }

    #[test]
    fn distributions() {
        let p = tmp("c.coverage.tsv");
        let rows = vec![ValueFrequency::new(0, 10), ValueFrequency::new(1, 3)];
        BamMetricCoverage::write(&p, &rows).unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), "Coverage\tCount\n0\t10\n1\t3\n");
        assert_eq!(BamMetricCoverage::read(&p).unwrap(), rows);
        let f = tmp("f.frag_length.tsv");
        BamMetricFragmentLength::write(&f, &rows).unwrap();
        assert!(std::fs::read_to_string(&f).unwrap().starts_with("FragmentLength\tCount\n"));
        assert_eq!(BamMetricFragmentLength::read(&f).unwrap(), rows);
        assert_eq!(BamMetricCoverage::generate_filename("dir", "S"), "dir/S.bam_metric.coverage.tsv");
    }

    #[test]
    fn flag_counts() {
        let mut s = [FlagQcStats::default(); FLAG_STAT_TYPES];
        s[FlagStatType::Total as usize] = FlagQcStats { passed: 1000, failed: 3 };
        s[FlagStatType::Primary as usize] = FlagQcStats { passed: 990, failed: 0 };
        s[FlagStatType::Secondary as usize] = FlagQcStats { passed: 2, failed: 0 };
        s[FlagStatType::Supplementary as usize] = FlagQcStats { passed: 8, failed: 0 };
        s[FlagStatType::Duplicate as usize] = FlagQcStats { passed: 100, failed: 0 };
        s[FlagStatType::Mapped as usize] = FlagQcStats { passed: 996, failed: 1 };
        s[FlagStatType::PrimaryMapped as usize] = FlagQcStats { passed: 986, failed: 0 };
        s[FlagStatType::Paired as usize] = FlagQcStats { passed: 990, failed: 0 };
        s[FlagStatType::ProperlyPaired as usize] = FlagQcStats { passed: 2, failed: 0 };
        s[FlagStatType::PairMapped as usize] = FlagQcStats { passed: 980, failed: 0 };
        s[FlagStatType::Singleton as usize] = FlagQcStats { passed: 6, failed: 0 };
        let lines = BamFlagCounts::lines(&s);
        assert_eq!(lines[0], "1000 + 3 in total (QC-passed reads + QC-failed reads)");
        assert_eq!(lines[6], "996 + 1 mapped (99.60% : 33.33%)");
        assert_eq!(lines[7], "986 + 0 primary mapped (99.60% : N/A)");
        assert_eq!(lines[11], "2 + 0 properly paired (0.20% : N/A)");
        let p = tmp("x.flag_counts.tsv");
        BamFlagCounts::write(&p, &s).unwrap();
        let r = BamFlagStats::read(&p).unwrap();
        assert_eq!(r.unique_read_count, 990);
        assert_eq!(r.duplicate_proportion, 0.1);
        assert_eq!(r.mapped_proportion, 0.996);
        assert_eq!(r.singleton_proportion, 0.006);
        assert_eq!(r.paired_in_sequencing_proportion, 1.0);
        assert_eq!(r.properly_paired_proportion, 0.00202);
    }

    #[test]
    fn percentages_in_float() {
        // 100f * 2 / 3 in float is 66.66667 (66.66666412353516 as a double)
        assert_eq!(flag_stats_percentages(&FlagQcStats { passed: 2, failed: 0 }, &FlagQcStats { passed: 3, failed: 0 }), "(66.67% : N/A)");
        // 1/8 exactly: %.2f rounds 12.5 to 12.50
        assert_eq!(flag_stats_percentages(&FlagQcStats { passed: 1, failed: 1 }, &FlagQcStats { passed: 8, failed: 800 }), "(12.50% : 0.13%)");
    }

    #[test]
    fn gene_and_exon_depths() {
        let ranges = vec![DepthRange { lower: 0, upper: 0 }, DepthRange { lower: 30, upper: 39 }, DepthRange { lower: 10000, upper: 10000 }];
        let gd = |gene: &str, counts: [i32; 3]| GeneDepth {
            gene: gene.into(),
            chromosome: "chr1".into(),
            pos_start: 10,
            pos_end: 20,
            depth_counts: ranges.iter().copied().zip(counts).collect(),
        };
        let mut genes = vec![gd("B", [1, 2, 3]), gd("A", [4, 5, 6])];
        let p = tmp("g.gene_coverage.tsv");
        GeneDepthFile::write_gene_coverage(&p, &mut genes).unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), "GeneName\tChromosome\tPosStart\tPosEnd\tDR_0\tDR_30_39\tDR_10000\nA\tchr1\t10\t20\t4\t5\t6\nB\tchr1\t10\t20\t1\t2\t3\n");
        assert_eq!(GeneDepthFile::read_gene_coverage(&p).unwrap(), vec![gd("A", [4, 5, 6]), gd("B", [1, 2, 3])]);
        let mut exons = vec![
            ExonDepth { depth: gd("B", [1, 0, 0]), exon_rank: 2, median_depth: 12.0, mean_depth: 12.345 },
            ExonDepth { depth: gd("B", [0, 1, 0]), exon_rank: 1, median_depth: 0.0, mean_depth: 0.125 },
        ];
        let e = tmp("e.exon_coverage.tsv");
        GeneDepthFile::write_exon_coverage(&e, &mut exons).unwrap();
        let text = std::fs::read_to_string(&e).unwrap();
        assert_eq!(text.lines().next().unwrap(), "GeneName\tChromosome\tPosStart\tPosEnd\tExonRank\tMedianDepth\tMeanDepth\tDR_0\tDR_30_39\tDR_10000");
        // 0.125 to two places: Java rounds half up, 0.13
        assert_eq!(text.lines().nth(1).unwrap(), "B\tchr1\t10\t20\t1\t0\t0.13\t0\t1\t0");
        assert_eq!(text.lines().nth(2).unwrap(), "B\tchr1\t10\t20\t2\t12\t12.35\t1\t0\t0");
        let r = GeneDepthFile::read_exon_coverage(&e).unwrap();
        assert_eq!((r[0].exon_rank, r[0].mean_depth, r[1].median_depth), (1, 0.13, 12.0));
        let empty = tmp("n.gene_coverage.tsv");
        GeneDepthFile::write_gene_coverage(&empty, &mut []).unwrap();
        assert_eq!(std::fs::read_to_string(&empty).unwrap(), "");
    }
}
