//! Jitter count tables and the output files (Java common.redux JitterCountsTable, JitterTableRow, JitterCountsTableFile,
//! JitterModelParams, JitterModelParamsFile, MsiModelPrediction and utils.file.DelimFileWriter).

use super::read::ConsensusType;
use crate::fmt::java_fixed;
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::io::Write;

/// JitterCountsTableFile.MAX_UNITS: rows above this repeat count are not written.
const MAX_UNITS: i32 = 20;

/// JitterTableRow: the reads of one reference repeat count, by jitter. Counts are Java ints.
#[derive(Clone, Debug)]
pub struct JitterTableRow {
    pub ref_num_units: i32,
    pub total_read_count: i32,
    pub jitter_counts: BTreeMap<i32, i32>,
}

impl JitterTableRow {
    pub fn new(ref_num_units: i32) -> Self {
        JitterTableRow { ref_num_units, total_read_count: 0, jitter_counts: BTreeMap::new() }
    }
    pub fn add_reads(&mut self, jitter: i32, num_reads: i32) {
        self.total_read_count = self.total_read_count.wrapping_add(num_reads);
        let c = self.jitter_counts.entry(jitter).or_insert(0);
        *c = c.wrapping_add(num_reads);
    }
    pub fn jitter_read_count(&self, jitter: i32) -> i32 {
        self.jitter_counts.get(&jitter).copied().unwrap_or(0)
    }
}

/// JitterCountsTable: the rows of one repeat unit group and consensus type, sorted by repeat count.
#[derive(Clone, Debug)]
pub struct JitterCountsTable {
    pub repeat_unit: String,
    pub consensus_type: ConsensusType,
    pub rows: Vec<JitterTableRow>,
}

impl JitterCountsTable {
    pub fn new(repeat_unit: &str, consensus_type: ConsensusType) -> Self {
        JitterCountsTable { repeat_unit: repeat_unit.to_string(), consensus_type, rows: Vec::new() }
    }
    /// JitterCountsTable.repeatUnitLength
    pub fn repeat_unit_length(&self) -> usize {
        match self.repeat_unit.split_once('/') {
            Some((first, _)) => first.len(),
            None => 3,
        }
    }
    pub fn get_or_create_row(&mut self, ref_num_units: i32) -> &mut JitterTableRow {
        let i = self.rows.partition_point(|r| r.ref_num_units < ref_num_units);
        if i >= self.rows.len() || self.rows[i].ref_num_units != ref_num_units {
            self.rows.insert(i, JitterTableRow::new(ref_num_units));
        }
        &mut self.rows[i]
    }
    pub fn row(&self, ref_num_units: i32) -> Option<&JitterTableRow> {
        self.rows.iter().find(|r| r.ref_num_units == ref_num_units)
    }
    pub fn read_count(&self, ref_num_units: i32) -> i32 {
        self.row(ref_num_units).map(|r| r.total_read_count).unwrap_or(0)
    }
    /// JitterCountsTable.mergeCounts with the rows of another table.
    pub fn merge_table(&mut self, other: &JitterCountsTable) {
        for row in &other.rows {
            let target = self.get_or_create_row(row.ref_num_units);
            for (&jitter, &count) in &row.jitter_counts {
                target.add_reads(jitter, count);
            }
        }
    }
}

/// JitterModelParams
#[derive(Clone, Debug)]
pub struct JitterModelParams {
    pub repeat_unit: String,
    pub consensus_type: ConsensusType,
    pub optimal_scale_repeat4: f64,
    pub optimal_scale_repeat5: f64,
    pub optimal_scale_repeat6: f64,
    pub scale_fit_gradient: f64,
    pub scale_fit_intercept: f64,
    pub microsatellite_skew: f64,
}

/// java.text.DecimalFormat("#.####") as DelimFileWriter applies it to doubles: at most four decimals, rounded
/// half-even on the exact binary value, no trailing zeros.
pub fn decimal_format_4(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "∞".into() } else { "-∞".into() };
    }
    let mut s = format!("{v:.4}");
    while s.ends_with('0') {
        s.pop();
    }
    if s.ends_with('.') {
        s.pop();
    }
    s
}

/// A text file writer that gzips when the name ends with .gz (FileWriterUtils.createBufferedWriter).
pub fn create_writer(filename: &str) -> Result<Box<dyn Write>> {
    let file = std::fs::File::create(filename).with_context(|| format!("failed to create {filename}"))?;
    if filename.ends_with(".gz") {
        Ok(Box::new(std::io::BufWriter::new(flate2::write::GzEncoder::new(file, flate2::Compression::default()))))
    } else {
        Ok(Box::new(std::io::BufWriter::new(file)))
    }
}

/// JitterCountsTableFile.generateFilename. `output_dir` ends with '/', and Java adds File.separator again.
pub fn counts_table_filename(output_dir: &str, sample: &str) -> String {
    format!("{output_dir}/{sample}.redux.ms_table.tsv.gz")
}

/// JitterModelParamsFile.generateFilename
pub fn model_params_filename(output_dir: &str, sample: &str) -> String {
    format!("{output_dir}/{sample}.redux.jitter_params.tsv")
}

/// MsiModelPrediction.generateFilename
pub fn msi_prediction_filename(output_dir: &str, sample: &str) -> String {
    format!("{output_dir}/{sample}.redux.msi_prediction.tsv")
}

/// The count-10 .. count+10 column names shared by the table and the site file.
pub fn count_columns() -> Vec<String> {
    let mut cols = Vec::new();
    for j in (1..=10).rev() {
        cols.push(format!("count-{j}"));
    }
    for j in 0..=10 {
        cols.push(format!("count+{j}"));
    }
    cols
}

/// JitterCountsTableFile.write
pub fn write_counts_tables(out: &mut dyn Write, tables: &[JitterCountsTable]) -> Result<()> {
    let mut columns: Vec<String> = ["unit", "consensusType", "numUnits", "readCount"].iter().map(|s| s.to_string()).collect();
    columns.extend(count_columns());
    writeln!(out, "{}", columns.join("\t"))?;
    for table in tables {
        for row in &table.rows {
            if row.ref_num_units > MAX_UNITS {
                continue;
            }
            write!(out, "{}\t{}\t{}\t{}", table.repeat_unit, table.consensus_type.name(), row.ref_num_units, row.total_read_count)?;
            for jitter in -10..=10 {
                write!(out, "\t{}", row.jitter_read_count(jitter))?;
            }
            writeln!(out)?;
        }
    }
    Ok(())
}

/// JitterModelParamsFile.write
pub fn write_model_params(out: &mut dyn Write, params: &[JitterModelParams]) -> Result<()> {
    writeln!(
        out,
        "unit\tconsensusType\toptimalScaleRepeat4\toptimalScaleRepeat5\toptimalScaleRepeat6\tscaleFitGradient\tscaleFitIntercept\tmicrosatelliteSkew"
    )?;
    for p in params {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p.repeat_unit,
            p.consensus_type.name(),
            decimal_format_4(p.optimal_scale_repeat4),
            decimal_format_4(p.optimal_scale_repeat5),
            decimal_format_4(p.optimal_scale_repeat6),
            decimal_format_4(p.scale_fit_gradient),
            decimal_format_4(p.scale_fit_intercept),
            decimal_format_4(p.microsatellite_skew)
        )?;
    }
    Ok(())
}

/// MsiModelPrediction.write: a header line and the value as %.4f.
pub fn write_msi_prediction(filename: &str, predicted: f64) -> Result<()> {
    let text = format!("PredictedValue\n{}\n", java_fixed(predicted, 4));
    std::fs::write(filename, text).with_context(|| format!("failed to write {filename}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimal_format() {
        assert_eq!(decimal_format_4(0.13), "0.13");
        assert_eq!(decimal_format_4(1.0), "1");
        assert_eq!(decimal_format_4(-0.16), "-0.16");
        assert_eq!(decimal_format_4(0.060000000000000005), "0.06");
        assert_eq!(decimal_format_4(0.69451), "0.6945");
        // exact ties round to even
        assert_eq!(decimal_format_4(0.03125), "0.0312");
        assert_eq!(decimal_format_4(0.09375), "0.0938");
        assert_eq!(decimal_format_4(0.0), "0");
        assert_eq!(decimal_format_4(-0.00001), "-0");
        assert_eq!(decimal_format_4(12345.67891), "12345.6789");
    }

    #[test]
    fn rows_stay_sorted() {
        let mut t = JitterCountsTable::new("A/T", ConsensusType::None);
        t.get_or_create_row(7).add_reads(0, 5);
        t.get_or_create_row(4).add_reads(-1, 2);
        t.get_or_create_row(7).add_reads(1, 1);
        t.get_or_create_row(5);
        let units: Vec<i32> = t.rows.iter().map(|r| r.ref_num_units).collect();
        assert_eq!(units, vec![4, 5, 7]);
        assert_eq!(t.read_count(7), 6);
        assert_eq!(t.read_count(5), 0);
        assert_eq!(t.read_count(9), 0);
        assert_eq!(t.repeat_unit_length(), 1);
        assert_eq!(JitterCountsTable::new("3-5bp repeat", ConsensusType::None).repeat_unit_length(), 3);
        assert_eq!(JitterCountsTable::new("AC/CA/GT/TG", ConsensusType::None).repeat_unit_length(), 2);
    }

    #[test]
    fn table_file_layout() {
        let mut t = JitterCountsTable::new("A/T", ConsensusType::Dual);
        t.get_or_create_row(4).add_reads(0, 7);
        t.get_or_create_row(4).add_reads(-11, 1);
        t.get_or_create_row(21).add_reads(0, 3);
        let mut out = Vec::new();
        write_counts_tables(&mut out, &[t]).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("unit\tconsensusType\tnumUnits\treadCount\tcount-10\tcount-9"));
        assert!(lines[0].ends_with("count-1\tcount+0\tcount+1\tcount+2\tcount+3\tcount+4\tcount+5\tcount+6\tcount+7\tcount+8\tcount+9\tcount+10"));
        assert_eq!(lines[1], "A/T\tDUAL\t4\t8\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t7\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0");
    }
}
