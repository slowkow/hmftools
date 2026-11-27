//! hmf-common `redux` files as qsee reads them (added 2026-11-27): BqrFile (`.redux.bqr.tsv`),
//! DuplicateFrequency (`.redux.duplicate_freq.tsv`) and JitterCountsTableFile (`.redux.ms_table.tsv.gz`),
//! with BqrKey, BqrRecord, JitterCountsTable and JitterTableRow. ConsensusType is kept as its ordinal
//! order (NONE, SINGLE, DUAL).

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;

use crate::delim_file::DelimFileReader;
use crate::java_text::{fields_index_map, parse_double, parse_int, parse_long, read_lines, split_all};

pub const REDUX_FILE_ID: &str = ".redux";

/// ConsensusType, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConsensusType {
    None,
    Single,
    Dual,
}

impl ConsensusType {
    pub fn name(self) -> &'static str {
        match self {
            ConsensusType::None => "NONE",
            ConsensusType::Single => "SINGLE",
            ConsensusType::Dual => "DUAL",
        }
    }
    /// ConsensusType.valueOf
    pub fn value_of(s: &str) -> Result<ConsensusType> {
        match s {
            "NONE" => Ok(ConsensusType::None),
            "SINGLE" => Ok(ConsensusType::Single),
            "DUAL" => Ok(ConsensusType::Dual),
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.bam.ConsensusType.{s}"),
        }
    }
}

// ---- BQR ----

/// BqrKey with its BqrRecord's count and recalibrated quality
#[derive(Clone, Debug, PartialEq)]
pub struct BqrRecord {
    pub ref_base: u8,
    pub alt: u8,
    pub trinucleotide_context: Vec<u8>,
    /// a Java byte: (byte)Integer.parseInt(value), so values above 127 wrap
    pub quality: i8,
    pub read_type: ConsensusType,
    pub count: i32,
    pub recalibrated_quality: f64,
}

/// BqrFile.generateFilename: basePath + File.separator + sample + ".redux.bqr.tsv"
pub fn bqr_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}{REDUX_FILE_ID}.bqr.tsv")
}

/// FileReaderUtils.getColumnIndex: the column whose name equals `column` ignoring case (Java takes the
/// first such entry in HashMap order; with one such column, as in REDUX's files, that is the column)
fn column_ignore_case(f: &HashMap<String, usize>, column: &str) -> Option<usize> {
    let mut found: Vec<(&String, usize)> = f.iter().filter(|(k, _)| k.eq_ignore_ascii_case(column)).map(|(k, v)| (k, *v)).collect();
    found.sort();
    found.first().map(|(_, v)| *v)
}

/// BqrFile.read. Java catches every exception here, logs it and returns null; this returns the error,
/// and the caller decides what Java's null leads to.
pub fn read_bqr_file(path: &str) -> Result<Vec<BqrRecord>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    bqr_records_from_lines(&lines).with_context(|| format!("failed to read BQR file({path})"))
}

/// BqrFile.read on lines already read
pub fn bqr_records_from_lines(lines: &[String]) -> Result<Vec<BqrRecord>> {
    let header = lines.first().ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: empty file"))?;
    let f = fields_index_map(header, '\t');
    let need = |c: &str| column_ignore_case(&f, c).ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {c}"));
    let ref_i = need("Ref")?;
    let alt_i = need("Alt")?;
    let tnc_i = need("TrinucleotideContext")?;
    let count_i = need("Count")?;
    let orig_i = need("OriginalQual")?;
    let recal_i = need("RecalibratedQual")?;
    let read_type_i = column_ignore_case(&f, "ReadType");
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let v = split_all(line, '\t');
        let get = |i: usize| -> Result<&str> { v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: index {i}")) };
        let first_byte = |i: usize| -> Result<u8> { get(i)?.as_bytes().first().copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0")) };
        let alt = first_byte(alt_i)?;
        let ref_base = first_byte(ref_i)?;
        let context = get(tnc_i)?.as_bytes().to_vec();
        let read_type = match read_type_i {
            Some(i) => ConsensusType::value_of(get(i)?)?,
            None => ConsensusType::None,
        };
        let num = |s: &str| parse_int(s).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""));
        let count = num(get(count_i)?)?;
        let quality = num(get(orig_i)?)? as i8;
        let recal = get(recal_i)?;
        let recalibrated_quality = parse_double(recal).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{recal}\""))?;
        out.push(BqrRecord { ref_base, alt, trinucleotide_context: context, quality, read_type, count, recalibrated_quality });
    }
    Ok(out)
}

// ---- duplicate frequencies ----

pub const NO_DUAL_STRAND: i64 = -1;

/// DuplicateFrequency
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DuplicateFrequency {
    pub read_count: i32,
    pub count: i64,
    pub dual_strand_count: i64,
}

/// DuplicateFrequency.generateFilename: checkAddDirSeparator(basePath) + sample + ".redux.duplicate_freq.tsv"
pub fn duplicate_freq_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}{REDUX_FILE_ID}.duplicate_freq.tsv", crate::metrics::check_add_dir_separator(base_path))
}

/// DuplicateFrequency.read
pub fn read_duplicate_frequencies(path: &str) -> Result<Vec<DuplicateFrequency>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    duplicate_frequencies_from_lines(&lines).with_context(|| format!("reading {path}"))
}

/// DuplicateFrequency.read on lines already read
pub fn duplicate_frequencies_from_lines(lines: &[String]) -> Result<Vec<DuplicateFrequency>> {
    let header = lines.first().ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: empty file"))?;
    let f = fields_index_map(header, '\t');
    let col = |c: &str| f.get(c).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {c}"));
    let read_i = col("DuplicateReadCount")?;
    let freq_i = col("Frequency")?;
    let dual_i = f.get("DualStrandFrequency").copied();
    let mut out = Vec::with_capacity(lines.len());
    for line in &lines[1..] {
        let v = split_all(line, '\t');
        let get = |i: usize| -> Result<&str> { v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: index {i}")) };
        let bad = |s: &str| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\"");
        let r = get(read_i)?;
        let c = get(freq_i)?;
        let read_count = parse_int(r).ok_or_else(|| bad(r))?;
        let count = parse_long(c).ok_or_else(|| bad(c))?;
        let dual_strand_count = match dual_i {
            Some(i) => {
                let d = get(i)?;
                parse_long(d).ok_or_else(|| bad(d))?
            }
            None => NO_DUAL_STRAND,
        };
        out.push(DuplicateFrequency { read_count, count, dual_strand_count });
    }
    Ok(out)
}

// ---- jitter counts tables ----

/// JitterTableRow: the jitter counts in a map by jitter (a Java HashMap; the readers here only look
/// counts up and sum them, so their order does not matter)
#[derive(Clone, Debug, PartialEq)]
pub struct JitterTableRow {
    pub ref_num_units: i32,
    pub repeat_unit: String,
    pub consensus_type: ConsensusType,
    pub total_read_count: i32,
    pub jitter_counts: Vec<(i32, i32)>,
}

impl JitterTableRow {
    pub fn new(ref_num_units: i32, repeat_unit: &str, consensus_type: ConsensusType) -> JitterTableRow {
        JitterTableRow { ref_num_units, repeat_unit: repeat_unit.to_string(), consensus_type, total_read_count: 0, jitter_counts: Vec::new() }
    }
    /// getJitterReadCount: 0 for a jitter not set
    pub fn jitter_read_count(&self, jitter: i32) -> i32 {
        self.jitter_counts.iter().find(|(j, _)| *j == jitter).map(|(_, c)| *c).unwrap_or(0)
    }
    /// setJitterReadCount: puts (or replaces) the count
    pub fn set_jitter_read_count(&mut self, jitter: i32, count: i32) {
        match self.jitter_counts.iter_mut().find(|(j, _)| *j == jitter) {
            Some(e) => e.1 = count,
            None => self.jitter_counts.push((jitter, count)),
        }
    }
}

/// JitterCountsTable: the rows of one repeat unit and consensus type, sorted by the reference number of
/// units
#[derive(Clone, Debug, PartialEq)]
pub struct JitterCountsTable {
    pub repeat_unit: String,
    pub consensus_type: ConsensusType,
    pub rows: Vec<JitterTableRow>,
}

impl JitterCountsTable {
    /// getOrCreateRow: the row of that number of units, inserted in order if new
    pub fn get_or_create_row(&mut self, ref_num_units: i32) -> &mut JitterTableRow {
        let mut at = self.rows.len();
        for (i, r) in self.rows.iter().enumerate() {
            if r.ref_num_units == ref_num_units {
                return &mut self.rows[i];
            }
            if r.ref_num_units > ref_num_units {
                at = i;
                break;
            }
        }
        self.rows.insert(at, JitterTableRow::new(ref_num_units, &self.repeat_unit, self.consensus_type));
        &mut self.rows[at]
    }
}

/// JitterCountsTableFile.generateFilename: basePath + File.separator + sample + ".redux.ms_table.tsv.gz"
pub fn jitter_counts_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}{REDUX_FILE_ID}.ms_table.tsv.gz")
}

/// JitterCountsTableFile.read: the tables in the order their unit and consensus type first appear
pub fn read_jitter_counts_tables(path: &str) -> Result<Vec<JitterCountsTable>> {
    let reader = DelimFileReader::open(path)?;
    let has_consensus = reader.has_column("consensusType");
    let mut tables: Vec<JitterCountsTable> = Vec::new();
    for row in reader.rows() {
        let unit = row.get("unit")?;
        let consensus = if has_consensus { ConsensusType::value_of(row.get("consensusType")?)? } else { ConsensusType::None };
        let ti = match tables.iter().position(|t| t.repeat_unit == unit && t.consensus_type == consensus) {
            Some(i) => i,
            None => {
                tables.push(JitterCountsTable { repeat_unit: unit.to_string(), consensus_type: consensus, rows: Vec::new() });
                tables.len() - 1
            }
        };
        let num_units = row.get_int("numUnits")?;
        // the values are read in Java's order: readCount, +0, +10 .. +1, -10 .. -1
        let read_count = row.get_int("readCount")?;
        let mut counts = vec![(0, row.get_int("count+0")?)];
        for j in (1..=10).rev() {
            counts.push((j, row.get_int(&format!("count+{j}"))?));
        }
        for j in (1..=10).rev() {
            counts.push((-j, row.get_int(&format!("count-{j}"))?));
        }
        let r = tables[ti].get_or_create_row(num_units);
        r.total_read_count = read_count;
        for (j, c) in counts {
            r.set_jitter_read_count(j, c);
        }
    }
    Ok(tables)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bqr_lines() {
        let lines: Vec<String> = ["Alt\tRef\tTrinucleotideContext\tReadType\tCount\tOriginalQual\tRecalibratedQual", "T\tC\tACA\tDUAL\t12\t200\t36.50"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let r = bqr_records_from_lines(&lines).unwrap();
        assert_eq!((r[0].alt, r[0].ref_base, r[0].read_type, r[0].count), (b'T', b'C', ConsensusType::Dual, 12));
        // (byte)200
        assert_eq!(r[0].quality, -56);
        assert_eq!(r[0].recalibrated_quality, 36.5);
        let no_type: Vec<String> = ["alt\tref\ttrinucleotidecontext\tcount\toriginalqual\trecalibratedqual", "T\tC\tACA\t1\t30\t30"].iter().map(|s| s.to_string()).collect();
        assert_eq!(bqr_records_from_lines(&no_type).unwrap()[0].read_type, ConsensusType::None);
        let bad: Vec<String> = [lines[0].clone(), "\tC\tACA\tNONE\t1\t30\t30".to_string()].to_vec();
        assert!(bqr_records_from_lines(&bad).is_err());
    }

    #[test]
    fn duplicate_frequency_lines() {
        let lines: Vec<String> = ["DuplicateReadCount\tFrequency", "1\t100", "2\t7"].iter().map(|s| s.to_string()).collect();
        let d = duplicate_frequencies_from_lines(&lines).unwrap();
        assert_eq!(d[1], DuplicateFrequency { read_count: 2, count: 7, dual_strand_count: NO_DUAL_STRAND });
        assert_eq!(duplicate_freq_filename("", "S"), "S.redux.duplicate_freq.tsv");
        assert_eq!(bqr_filename("", "S"), "/S.redux.bqr.tsv");
    }

    #[test]
    fn jitter_rows_in_order() {
        let mut t = JitterCountsTable { repeat_unit: "A/T".into(), consensus_type: ConsensusType::None, rows: Vec::new() };
        t.get_or_create_row(5).total_read_count = 5;
        t.get_or_create_row(3).total_read_count = 3;
        t.get_or_create_row(4).total_read_count = 4;
        t.get_or_create_row(5).set_jitter_read_count(1, 2);
        assert_eq!(t.rows.iter().map(|r| r.ref_num_units).collect::<Vec<_>>(), vec![3, 4, 5]);
        assert_eq!((t.rows[2].jitter_read_count(1), t.rows[2].jitter_read_count(-1)), (2, 0));
    }
}
