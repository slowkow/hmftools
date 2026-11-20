//! The MSI prediction written with the jitter output (Java redux.ms_model: MsModelCalculator, MsModelCoefficients,
//! MsModelErrorRates, MsModelParams, RepeatUnitData). The training application (MsModelBuilder) is not part of this.

use super::fit::{java_max, java_pow};
use super::read::ConsensusType;
use super::table::{JitterCountsTable, JitterTableRow};
use crate::config::SequencingType;
use anyhow::{anyhow, bail, Context, Result};
use std::collections::{HashMap, HashSet};

const REPEAT_3_5_GROUP: &str = "3-5bp repeat";
const MULTI_BASE_REPEAT: &str = "bp repeat";
const TRAINING_REPEAT_UNIT_MIN_READ_COUNT: i32 = 500;
pub const INVALID_VALUE: f64 = -1.0;

/// MsModelParams: the fields the prediction uses.
#[derive(Clone, Debug)]
pub struct MsModelParams {
    pub repeat_types: Vec<&'static str>,
    pub repeat_count_min: i32,
    pub repeat_count_max: i32,
}

impl MsModelParams {
    /// MsModelParams.sequencingTypeParams
    pub fn for_sequencing_type(sequencing_type: SequencingType) -> Self {
        match sequencing_type {
            SequencingType::Ultima => {
                MsModelParams { repeat_types: vec!["AC/CA/GT/TG", "AG/GA/CT/TC"], repeat_count_min: 11, repeat_count_max: 15 }
            }
            SequencingType::Sbx => MsModelParams { repeat_types: vec!["A/T"], repeat_count_min: 11, repeat_count_max: 15 },
            SequencingType::Illumina => MsModelParams { repeat_types: vec!["A/T"], repeat_count_min: 10, repeat_count_max: 20 },
        }
    }
}

/// RepeatUnitData.repeatUnitAndCountKey, e.g. "A/T_10"
fn repeat_unit_and_count_key(repeat_unit: &str, repeat_count: i32) -> String {
    format!("{repeat_unit}_{repeat_count}")
}

/// RepeatUnitData
#[derive(Clone, Debug)]
pub struct RepeatUnitData {
    pub repeat_unit: String,
    pub consensus_type: ConsensusType,
    pub jitter_row: JitterTableRow,
    pub ruc_key: String,
    error_rate: f64,
    adjusted_error_rate: f64,
}

impl RepeatUnitData {
    fn new(repeat_unit: &str, consensus_type: ConsensusType, row: &JitterTableRow) -> Self {
        RepeatUnitData {
            repeat_unit: repeat_unit.to_string(),
            consensus_type,
            jitter_row: row.clone(),
            ruc_key: repeat_unit_and_count_key(repeat_unit, row.ref_num_units),
            error_rate: 0.0,
            adjusted_error_rate: 0.0,
        }
    }
    /// RepeatUnitData.computeErrorRates: jitter-weighted reads over jitter-weighted reads plus unchanged reads, in ints.
    pub fn compute_error_rates(&mut self) {
        let mut weighted_error_count: i32 = 0;
        let mut weighted_total: i32 = 0;
        for (&jitter, &reads) in &self.jitter_row.jitter_counts {
            let count = jitter.wrapping_abs();
            if count == 0 {
                weighted_total = weighted_total.wrapping_add(reads);
            } else {
                weighted_error_count = weighted_error_count.wrapping_add(count.wrapping_mul(reads));
            }
        }
        weighted_total = weighted_total.wrapping_add(weighted_error_count);
        self.error_rate = if weighted_total > 0 { weighted_error_count as f64 / weighted_total as f64 } else { 0.0 };
    }
    pub fn set_adjusted_error_rate(&mut self, noise_adjustment: f64) {
        self.adjusted_error_rate = java_max(self.error_rate - noise_adjustment, 0.0);
    }
    #[cfg(test)]
    pub fn error_rate(&self) -> f64 {
        self.error_rate
    }
    pub fn adjusted_error_rate(&self) -> f64 {
        self.adjusted_error_rate
    }
}

/// The rows of a tab-separated file with a header, as (column index by name, split lines). Java reads all lines
/// and fails on one with too few values; so does the caller here.
fn read_tsv(filename: &str) -> Result<(HashMap<String, usize>, Vec<Vec<String>>)> {
    let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read {filename}"))?;
    let mut lines = text.lines();
    let header = lines.next().ok_or_else(|| anyhow!("empty file: {filename}"))?;
    let index: HashMap<String, usize> = header.split('\t').enumerate().map(|(i, c)| (c.to_string(), i)).collect();
    let rows = lines.map(|l| l.split('\t').map(str::to_string).collect()).collect();
    Ok((index, rows))
}

fn column(index: &HashMap<String, usize>, name: &str, filename: &str) -> Result<usize> {
    index.get(name).copied().ok_or_else(|| anyhow!("missing column {name} in {filename}"))
}

fn field<'a>(values: &'a [String], i: usize, filename: &str) -> Result<&'a str> {
    values.get(i).map(|s| s.as_str()).ok_or_else(|| anyhow!("too few values in a line of {filename}"))
}

/// MsModelCalculator as used for one sample: coefficients and error rates by repeat unit and count.
pub struct MsModelCalculator {
    params: MsModelParams,
    coefficients: HashMap<String, [f64; 2]>,
    error_rates: HashMap<String, f64>,
}

impl MsModelCalculator {
    pub fn new(params: MsModelParams) -> Self {
        MsModelCalculator { params, coefficients: HashMap::new(), error_rates: HashMap::new() }
    }

    /// The constructor for external use: MsModelCoefficients.read and MsModelErrorRates.read.
    pub fn from_files(params: MsModelParams, coefficients_file: &str, error_rates_file: &str) -> Result<Self> {
        let mut calc = MsModelCalculator::new(params);
        let (index, rows) = read_tsv(coefficients_file)?;
        let (unit, num, c1, c2) = (
            column(&index, "unit", coefficients_file)?,
            column(&index, "numUnits", coefficients_file)?,
            column(&index, "coefficient1", coefficients_file)?,
            column(&index, "coefficient2", coefficients_file)?,
        );
        for values in &rows {
            let f = |i| field(values, i, coefficients_file);
            calc.add_coefficients(f(unit)?, f(num)?.parse()?, f(c1)?.parse()?, f(c2)?.parse()?);
        }
        let (index, rows) = read_tsv(error_rates_file)?;
        let (unit, num, rate) =
            (column(&index, "unit", error_rates_file)?, column(&index, "numUnits", error_rates_file)?, column(&index, "errorRate", error_rates_file)?);
        for values in &rows {
            let f = |i| field(values, i, error_rates_file);
            calc.add_error_rate(f(unit)?, f(num)?.parse()?, f(rate)?.parse()?);
        }
        Ok(calc)
    }

    pub fn add_coefficients(&mut self, repeat_unit: &str, repeat_count: i32, coefficient1: f64, coefficient2: f64) {
        self.coefficients.insert(repeat_unit_and_count_key(repeat_unit, repeat_count), [coefficient1, coefficient2]);
    }

    pub fn add_error_rate(&mut self, repeat_unit: &str, repeat_count: i32, error_rate: f64) {
        self.error_rates.insert(repeat_unit_and_count_key(repeat_unit, repeat_count), error_rate);
    }

    /// MsModelCalculator.calcPredictedMsiRate
    pub fn calc_predicted_msi_rate(&self, tables: &[JitterCountsTable]) -> Result<f64> {
        let mut data = self.build_repeat_unit_data(tables);
        let values = self.calc_ms_indel_per_mb_values(&mut data)?;
        Ok(calc_ms_indel_per_mb(&values))
    }

    /// MsModelCalculator.buildRepeatUnitData: the rows of the model's repeat units within its repeat count range, then
    /// for each unit and count the first consensus type with the highest read count, dropping those under 500 reads.
    pub fn build_repeat_unit_data(&self, tables: &[JitterCountsTable]) -> Vec<RepeatUnitData> {
        let mut list: Vec<RepeatUnitData> = Vec::new();
        for table in tables {
            if !self.params.repeat_types.iter().any(|t| *t == table.repeat_unit) {
                continue;
            }
            let is_multi_base_repeat = table.repeat_unit.contains(MULTI_BASE_REPEAT);
            for row in &table.rows {
                if row.ref_num_units < self.params.repeat_count_min || row.ref_num_units > self.params.repeat_count_max {
                    continue;
                }
                if is_multi_base_repeat {
                    // consolidate the 3-5 base units. Java merges into the first table's own row object; the tables
                    // are not used again afterwards, so a copy gives the same result
                    let existing = list.iter_mut().find(|x| {
                        x.repeat_unit == REPEAT_3_5_GROUP && x.jitter_row.ref_num_units == row.ref_num_units && x.consensus_type == table.consensus_type
                    });
                    match existing {
                        Some(data) => {
                            for (&jitter, &count) in &row.jitter_counts {
                                data.jitter_row.add_reads(jitter, count);
                            }
                        }
                        None => list.push(RepeatUnitData::new(REPEAT_3_5_GROUP, table.consensus_type, row)),
                    }
                } else {
                    list.push(RepeatUnitData::new(&table.repeat_unit, table.consensus_type, row));
                }
            }
        }

        let mut processed: HashSet<&str> = HashSet::new();
        let mut filtered = Vec::new();
        for data in &list {
            if !processed.insert(&data.ruc_key) {
                continue;
            }
            let max_count = list.iter().filter(|x| x.ruc_key == data.ruc_key).map(|x| x.jitter_row.total_read_count).max().unwrap_or(0);
            if max_count < TRAINING_REPEAT_UNIT_MIN_READ_COUNT {
                continue;
            }
            if let Some(best) = list.iter().find(|x| x.ruc_key == data.ruc_key && x.jitter_row.total_read_count == max_count) {
                filtered.push(best.clone());
            }
        }
        filtered
    }

    /// MsModelCalculator.calcMsIndelPerMbValues. Java unboxes the error rate before it looks at the coefficients,
    /// so a unit and count without an error rate is a NullPointerException there and an error here.
    pub fn calc_ms_indel_per_mb_values(&self, data: &mut [RepeatUnitData]) -> Result<Vec<f64>> {
        let mut values = Vec::new();
        for d in data.iter_mut() {
            d.compute_error_rates();
            let Some(&combined_error_rate) = self.error_rates.get(&d.ruc_key) else {
                bail!("java.lang.NullPointerException: no MSI model error rate for {}", d.ruc_key);
            };
            let coefficients = match self.coefficients.get(&d.ruc_key) {
                Some(c) if c[0] != 0.0 && c[1] != 0.0 => c,
                _ => {
                    values.push(INVALID_VALUE);
                    continue;
                }
            };
            d.set_adjusted_error_rate(combined_error_rate);
            let x = d.adjusted_error_rate();
            values.push(java_max(coefficients[0] * x + coefficients[1] * java_pow(x, 3.0), 0.0));
        }
        Ok(values)
    }
}

/// MsModelCalculator.calcMsIndelPerMb: the mean of the valid values, or -1.
pub fn calc_ms_indel_per_mb(values: &[f64]) -> f64 {
    let valid: Vec<f64> = values.iter().copied().filter(|v| *v != INVALID_VALUE).collect();
    if valid.is_empty() {
        return INVALID_VALUE;
    }
    let mut sum = 0.0;
    for v in &valid {
        sum += v;
    }
    sum / valid.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(unit: &str, ct: ConsensusType, rows: &[(i32, &[(i32, i32)])]) -> JitterCountsTable {
        let mut t = JitterCountsTable::new(unit, ct);
        for (units, counts) in rows {
            let row = t.get_or_create_row(*units);
            for (j, c) in counts.iter() {
                row.add_reads(*j, *c);
            }
        }
        t
    }

    fn illumina() -> MsModelCalculator {
        let mut calc = MsModelCalculator::new(MsModelParams::for_sequencing_type(SequencingType::Illumina));
        for n in 10..=20 {
            calc.add_coefficients("A/T", n, 100.0 + n as f64, 2000.0);
            calc.add_error_rate("A/T", n, 0.01);
        }
        calc
    }

    #[test]
    fn error_rate_is_jitter_weighted() {
        let t = table("A/T", ConsensusType::None, &[(10, &[(0, 900), (-1, 50), (2, 25)])]);
        let mut d = RepeatUnitData::new("A/T", ConsensusType::None, &t.rows[0]);
        d.compute_error_rates();
        assert_eq!(d.error_rate(), 100.0 / 1000.0);
        d.set_adjusted_error_rate(0.04);
        assert_eq!(d.adjusted_error_rate(), 0.1 - 0.04);
        d.set_adjusted_error_rate(0.4);
        assert_eq!(d.adjusted_error_rate(), 0.0);
        assert_eq!(d.ruc_key, "A/T_10");
    }

    #[test]
    fn selects_units_counts_and_best_consensus_type() {
        let calc = illumina();
        let tables = vec![
            table("A/T", ConsensusType::None, &[(9, &[(0, 5000)]), (10, &[(0, 600)]), (11, &[(0, 499)]), (12, &[(0, 700)]), (21, &[(0, 9000)])]),
            table("C/G", ConsensusType::None, &[(10, &[(0, 5000)])]),
            table("A/T", ConsensusType::Dual, &[(10, &[(0, 800)]), (12, &[(0, 700)])]),
        ];
        let data = calc.build_repeat_unit_data(&tables);
        let got: Vec<(String, ConsensusType, i32)> = data.iter().map(|d| (d.ruc_key.clone(), d.consensus_type, d.jitter_row.total_read_count)).collect();
        assert_eq!(
            got,
            vec![("A/T_10".to_string(), ConsensusType::Dual, 800), ("A/T_12".to_string(), ConsensusType::None, 700)]
        );
    }

    #[test]
    fn predicted_value() {
        let calc = illumina();
        let tables = vec![table("A/T", ConsensusType::None, &[(10, &[(0, 900), (-1, 100)]), (11, &[(0, 995), (1, 5)])])];
        // 10: x = 0.1 - 0.01; 11: error rate 0.005 is under the noise rate, x = 0
        let x: f64 = 0.1 - 0.01;
        let expected = (110.0 * x + 2000.0 * java_pow(x, 3.0) + 0.0) / 2.0;
        assert_eq!(calc.calc_predicted_msi_rate(&tables).unwrap(), expected);
        // nothing with enough reads
        assert_eq!(calc.calc_predicted_msi_rate(&[]).unwrap(), INVALID_VALUE);
    }

    #[test]
    fn missing_coefficients_and_error_rates() {
        let mut calc = MsModelCalculator::new(MsModelParams::for_sequencing_type(SequencingType::Sbx));
        calc.add_error_rate("A/T", 11, 0.01);
        calc.add_error_rate("A/T", 12, 0.01);
        calc.add_coefficients("A/T", 12, 0.0, 5.0);
        let tables = vec![table("A/T", ConsensusType::None, &[(10, &[(0, 900)]), (11, &[(0, 900)]), (12, &[(0, 900)])])];
        // 10 is outside the SBX range; 11 has no coefficients and 12 a zero one: both invalid
        assert_eq!(calc.calc_predicted_msi_rate(&tables).unwrap(), INVALID_VALUE);
        // a row in range without an error rate fails as Java does
        let tables = vec![table("A/T", ConsensusType::None, &[(13, &[(0, 900)])])];
        assert!(calc.calc_predicted_msi_rate(&tables).is_err());
    }

    #[test]
    fn multi_base_units_are_merged() {
        let mut calc = MsModelCalculator::new(MsModelParams { repeat_types: vec!["3bp repeat", "4bp repeat"], repeat_count_min: 4, repeat_count_max: 6 });
        calc.add_error_rate("3-5bp repeat", 4, 0.0);
        calc.add_coefficients("3-5bp repeat", 4, 10.0, 1.0);
        let tables = vec![
            table("3bp repeat", ConsensusType::None, &[(4, &[(0, 400), (1, 50)])]),
            table("4bp repeat", ConsensusType::None, &[(4, &[(0, 500), (1, 50)])]),
            table("5bp repeat", ConsensusType::None, &[(4, &[(0, 500)])]),
        ];
        let data = calc.build_repeat_unit_data(&tables);
        assert_eq!(data.len(), 1);
        assert_eq!(data[0].ruc_key, "3-5bp repeat_4");
        assert_eq!(data[0].jitter_row.total_read_count, 1000);
        assert_eq!(data[0].jitter_row.jitter_read_count(1), 100);
    }

    #[test]
    fn reads_model_files() {
        let dir = std::env::temp_dir().join(format!("redux_rs_ms_model_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let coef = dir.join("coef.tsv");
        let rates = dir.join("rates.tsv");
        std::fs::write(&coef, "unit\tnumUnits\tcoefficient1\tcoefficient2\nA/T\t10\t1031.5\t37724.25\n").unwrap();
        std::fs::write(&rates, "unit\tnumUnits\terrorRate\nA/T\t10\t0.125\n").unwrap();
        let calc =
            MsModelCalculator::from_files(MsModelParams::for_sequencing_type(SequencingType::Illumina), coef.to_str().unwrap(), rates.to_str().unwrap())
                .unwrap();
        assert_eq!(calc.coefficients["A/T_10"], [1031.5, 37724.25]);
        assert_eq!(calc.error_rates["A/T_10"], 0.125);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
