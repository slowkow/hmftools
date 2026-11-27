//! qsee `cohort`: FeatureMatrix, PercentileTransformer, FeaturePercentiles, CohortPercentiles and its
//! file, CohortPercentilesTrainer (with TrainConfig) and CohortFeaturesWriter.

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;
use std::io::Write;

use hmf_common::delim_file::{write_delim_file, DelimFileReader, JavaWriter};
use hmf_common::java_math::{double_stream_sum, sort_doubles};

use crate::common::{decimal_8, generate_cohort_filename, interpolate_1d, interpolate_2d, SampleType, COL_FEATURE_NAME, COL_FEATURE_TYPE, COL_SAMPLE_TYPE, COL_SOURCE_TOOL};
use crate::compat;
use crate::config::{PrepCategory, QseeConfig};
use crate::feature::{Feature, FeatureKey, FeatureType, SourceTool};
use crate::prep::{run_category_task, PrepError};

// ---- FeatureMatrix ----

/// FeatureMatrix: features by sample (rows) and key (columns, in the order first added); a sample without
/// a feature has a NaN feature of the column's first key
#[derive(Clone, Debug, Default)]
pub struct FeatureMatrix {
    pub sample_ids: Vec<String>,
    pub keys: Vec<FeatureKey>,
    columns: HashMap<FeatureKey, Vec<Feature>>,
    num_samples: usize,
}

impl FeatureMatrix {
    /// new FeatureMatrix(map, expectedSampleIds)
    pub fn with_samples(sample_ids: &[String]) -> FeatureMatrix {
        FeatureMatrix { sample_ids: sample_ids.to_vec(), keys: Vec::new(), columns: HashMap::new(), num_samples: sample_ids.len() }
    }
    /// new FeatureMatrix(map, numSamples)
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn with_num_samples(n: usize) -> FeatureMatrix {
        FeatureMatrix { num_samples: n, ..Default::default() }
    }

    /// addRow: a row index past the expected number of samples is an error (ArrayIndexOutOfBounds)
    pub fn add_row(&mut self, row_id: &str, features: Vec<Feature>) -> Result<()> {
        if !self.sample_ids.iter().any(|s| s == row_id) {
            self.sample_ids.push(row_id.to_string());
        }
        let row = self.sample_ids.iter().position(|s| s == row_id).unwrap();
        for f in features {
            if !self.columns.contains_key(&f.key) {
                self.keys.push(f.key.clone());
                self.columns.insert(f.key.clone(), vec![Feature::new(f.key.clone(), f64::NAN); self.num_samples]);
            }
            let col = self.columns.get_mut(&f.key).unwrap();
            if row >= col.len() {
                bail!("java.lang.ArrayIndexOutOfBoundsException: Index {row} out of bounds for length {}", col.len());
            }
            col[row] = f;
        }
        Ok(())
    }

    /// sortFeatureKeys: Collections.sort, stable
    pub fn sort_feature_keys(&mut self) {
        self.keys.sort_by(|a, b| a.compare(b));
    }

    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn num_features(&self) -> usize {
        self.keys.len()
    }

    pub fn row(&self, index: usize) -> Vec<Feature> {
        self.keys.iter().map(|k| self.columns[k][index].clone()).collect()
    }

    /// getRow(rowId)
    pub fn row_of(&self, row_id: &str) -> Result<Vec<Feature>> {
        let i = self.sample_ids.iter().position(|s| s == row_id).ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: RowId({row_id}) not found"))?;
        Ok(self.row(i))
    }

    pub fn column(&self, key: &FeatureKey) -> &[Feature] {
        &self.columns[key]
    }

    /// getFeatureValues: the values by row
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn feature_values(&self) -> Vec<Vec<f64>> {
        (0..self.num_samples).map(|i| self.row(i).iter().map(|f| f.value).collect()).collect()
    }
}

// ---- PercentileTransformer ----

/// PercentileTransformer
#[derive(Clone, Debug)]
pub struct PercentileTransformer {
    pub percentiles: Vec<f64>,
    pub ref_values: Vec<f64>,
    pub percentiles_deduped: Vec<f64>,
    pub ref_values_deduped: Vec<f64>,
    fitted: bool,
    fitted_with_all_nan: bool,
}

impl PercentileTransformer {
    fn unfitted(percentiles: Vec<f64>) -> PercentileTransformer {
        let n = percentiles.len();
        PercentileTransformer { percentiles, ref_values: vec![f64::NAN; n], percentiles_deduped: Vec::new(), ref_values_deduped: Vec::new(), fitted: false, fitted_with_all_nan: false }
    }

    /// withNumPercentiles: n percentiles from 0 to 100 in steps of 100 / (n - 1)
    pub fn with_num_percentiles(n: i32) -> Result<PercentileTransformer> {
        if n <= 0 {
            bail!("java.lang.IllegalArgumentException: Number of percentiles should be greater than 0");
        }
        let increment = 100.0 / (n - 1) as f64;
        Ok(Self::unfitted((0..n).map(|i| i as f64 * increment).collect()))
    }

    /// withIncrement
    pub fn with_increment(increment: f64) -> Result<PercentileTransformer> {
        if increment <= 0.0 || increment >= 100.0 {
            bail!("java.lang.IllegalArgumentException: Percentile increment should be between 0 and 100");
        }
        // QS12: Java's floating-point remainder rejects increments that divide 100 but are not exact in
        // binary, such as 0.1 (100 % 0.1 is 0.09999999999999445), the example of its own message
        let divides = if compat::java("QS12") { 100.0 % increment == 0.0 } else { ((100.0 / increment).round() * increment - 100.0).abs() <= 1e-9 * 100.0 };
        if !divides {
            bail!("java.lang.IllegalArgumentException: Percentile increment should evenly divide 100 (e.g. 0.1, 1, 10, 25, etc)");
        }
        let n = (100.0 / increment).ceil() as i32 + 1;
        Self::with_num_percentiles(n)
    }

    /// fromPrefitData
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn from_prefit(percentiles: &[f64], ref_values: &[f64]) -> Result<PercentileTransformer> {
        if percentiles.len() != ref_values.len() {
            bail!("java.lang.IllegalArgumentException: Percentiles and ref values should be the same length");
        }
        let not_sorted = |v: &[f64]| (1..v.len()).any(|i| v[i] < v[i - 1]);
        if not_sorted(percentiles) || not_sorted(ref_values) {
            bail!("java.lang.IllegalArgumentException: Percentiles and ref values should be sorted in ascending order");
        }
        if percentiles.iter().any(|p| *p < 0.0 || *p > 100.0) {
            bail!("java.lang.IllegalArgumentException: Percentiles should be between 0 and 100");
        }
        if percentiles.len() < 2 || percentiles[0] != 0.0 || percentiles[percentiles.len() - 1] != 100.0 {
            bail!("java.lang.IllegalArgumentException: Data for at least percentiles 0 and 100 should be provided");
        }
        let mut t = Self::unfitted(percentiles.to_vec());
        t.ref_values = ref_values.to_vec();
        t.fitted = true;
        t.dedup();
        Ok(t)
    }

    /// fit: the reference value at each percentile of the finite values (linear between neighbours)
    pub fn fit(&mut self, values: &[f64]) -> Result<()> {
        if self.fitted {
            bail!("java.lang.IllegalStateException: PercentileTransformer already fitted");
        }
        if values.is_empty() {
            bail!("java.lang.IllegalArgumentException: Cohort values must not be empty");
        }
        let mut finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
        if finite.is_empty() {
            self.fitted = true;
            self.fitted_with_all_nan = true;
            return Ok(());
        }
        sort_doubles(&mut finite);
        for i in 0..self.percentiles.len() {
            self.ref_values[i] = ref_value_at_percentile(self.percentiles[i], &finite);
        }
        self.dedup();
        self.fitted = true;
        Ok(())
    }

    /// dedupWithRefValues: each distinct reference value (in order, a LinkedHashMap keyed by the boxed
    /// Double, so NaN equals NaN and 0.0 differs from -0.0) with the mean of its percentiles
    fn dedup(&mut self) {
        let mut groups: Vec<(f64, Vec<f64>)> = Vec::new();
        let mut index: HashMap<i64, usize> = HashMap::new();
        for (r, p) in self.ref_values.iter().zip(&self.percentiles) {
            let bits = hmf_common::java_math::double_to_long_bits(*r);
            match index.get(&bits) {
                Some(&i) => groups[i].1.push(*p),
                None => {
                    index.insert(bits, groups.len());
                    groups.push((*r, vec![*p]));
                }
            }
        }
        self.ref_values_deduped = groups.iter().map(|(r, _)| *r).collect();
        self.percentiles_deduped = groups.iter().map(|(_, ps)| if ps.len() == 1 { ps[0] } else { double_stream_sum(ps.iter().copied()) / ps.len() as f64 }).collect();
    }

    /// featureValueToPercentile
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn feature_value_to_percentile(&self, v: f64) -> Result<f64> {
        if !self.fitted {
            bail!("java.lang.IllegalStateException: PercentileTransformer not fitted");
        }
        if self.fitted_with_all_nan || v.is_nan() {
            return Ok(f64::NAN);
        }
        Ok(interpolate_2d(v, &self.ref_values_deduped, &self.percentiles_deduped))
    }

    /// percentileToFeatureValue
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn percentile_to_feature_value(&self, p: f64) -> Result<f64> {
        if !self.fitted {
            bail!("java.lang.IllegalStateException: PercentileTransformer not fitted");
        }
        if p.is_nan() {
            bail!("java.lang.IllegalArgumentException: Percentile cannot be NaN");
        }
        if self.fitted_with_all_nan {
            return Ok(f64::NAN);
        }
        Ok(interpolate_2d(p, &self.percentiles, &self.ref_values))
    }
}

fn ref_value_at_percentile(percentile: f64, sorted: &[f64]) -> f64 {
    let position = (percentile / 100.0) * (sorted.len() - 1) as f64;
    let lower = position.floor() as i32;
    let upper = position.ceil() as i32;
    if lower == upper {
        sorted[lower as usize]
    } else {
        interpolate_1d(sorted[lower as usize], sorted[upper as usize], position - lower as f64)
    }
}

// ---- cohort percentiles ----

/// FeaturePercentiles
#[derive(Clone, Debug)]
pub struct FeaturePercentiles {
    pub sample_type: SampleType,
    pub key: FeatureKey,
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub percentiles: Vec<f64>,
    pub ref_values: Vec<f64>,
}

/// CohortPercentiles: by sample type, the percentiles of each feature
#[derive(Clone, Debug, Default)]
pub struct CohortPercentiles {
    pub data: Vec<(SampleType, HashMap<FeatureKey, FeaturePercentiles>)>,
}

impl CohortPercentiles {
    pub fn add(&mut self, fp: FeaturePercentiles) {
        let i = match self.data.iter().position(|(t, _)| *t == fp.sample_type) {
            Some(i) => i,
            None => {
                self.data.push((fp.sample_type, HashMap::new()));
                self.data.len() - 1
            }
        };
        self.data[i].1.insert(fp.key.clone(), fp);
    }

    /// warnIfMissing (getFeaturePercentiles): a warning for a feature the cohort does not have. Java
    /// looks the sample type up first and fails (NullPointerException) when the file has no row of that
    /// type (QS3).
    pub fn warn_if_missing(&self, sample_type: SampleType, key: &FeatureKey) -> Result<()> {
        let by_type = self.data.iter().find(|(t, _)| *t == sample_type).map(|(_, m)| m);
        let found = match by_type {
            Some(m) => m.contains_key(key),
            None => {
                if compat::java("QS3") {
                    bail!("java.lang.NullPointerException: no cohort percentiles of sample type {}", sample_type.name());
                }
                false
            }
        };
        if !found {
            hmf_common::warn!("Missing cohort percentiles for sampleType({}) feature({})", sample_type.name(), key.display());
        }
        Ok(())
    }
}

/// CohortPercentilesFile.read: any failure is logged and ends the run
pub fn read_cohort_percentiles(path: &str) -> Result<CohortPercentiles> {
    read_cohort_percentiles_inner(path).map_err(|e| anyhow!("Failed to load cohort percentiles file: {path}: {e:#}"))
}

fn read_cohort_percentiles_inner(path: &str) -> Result<CohortPercentiles> {
    let reader = DelimFileReader::open(path)?;
    let need = |c: &str| reader.column_index(c).ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {c}"));
    let (i_type, i_ftype, i_name, i_tool) = (need(COL_SAMPLE_TYPE)?, need(COL_FEATURE_TYPE)?, need(COL_FEATURE_NAME)?, need(COL_SOURCE_TOOL)?);
    let mut pct_columns = Vec::new();
    let mut percentiles = Vec::new();
    for (i, c) in reader.column_names().iter().enumerate() {
        if c.starts_with("Pct") {
            let p = c.splitn(2, '_').nth(1).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: column {c}"))?;
            percentiles.push(hmf_common::java_text::parse_double(p).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{p}\""))?);
            pct_columns.push(i);
        }
    }
    let mut cohort = CohortPercentiles::default();
    for row in reader.rows() {
        let sample_type = SampleType::value_of(row.get_at(i_type)?)?;
        let name = row.get_at(i_name)?;
        let feature_type = FeatureType::value_of(row.get_at(i_ftype)?)?;
        let tool = SourceTool::value_of(row.get_at(i_tool)?)?;
        let key = FeatureKey::new(name, feature_type, tool);
        let mut values = Vec::with_capacity(pct_columns.len());
        for i in &pct_columns {
            values.push(row.get_double_at(*i)?);
        }
        cohort.add(FeaturePercentiles { sample_type, key, percentiles: percentiles.clone(), ref_values: values });
    }
    Ok(cohort)
}

/// CohortPercentilesFile.write
pub fn write_cohort_percentiles(path: &str, percentiles: &[f64], cohort: &[FeaturePercentiles]) -> Result<()> {
    let names: Vec<String> = percentiles.iter().map(|p| format!("Pct_{}", decimal_8(*p))).collect();
    let mut columns: Vec<String> = [COL_SAMPLE_TYPE, COL_SOURCE_TOOL, COL_FEATURE_TYPE, COL_FEATURE_NAME].iter().map(|s| s.to_string()).collect();
    columns.extend(names.iter().cloned());
    write_delim_file(path, &columns, cohort, |fp, row| {
        row.set(COL_SAMPLE_TYPE, fp.sample_type.name())?;
        row.set(COL_SOURCE_TOOL, fp.key.source_tool.map(|t| t.name()).unwrap_or("null"))?;
        row.set(COL_FEATURE_TYPE, fp.key.feature_type.map(|t| t.name()).unwrap_or("null"))?;
        row.set(COL_FEATURE_NAME, fp.key.name.as_str())?;
        for i in 0..names.len() {
            row.set_at(4 + i, decimal_8(fp.ref_values[i]));
        }
        Ok(())
    })
}

// ---- training ----

/// FeaturePrep.prepCohortCategory: the category for every sample of the type, added to the matrix in
/// sample order. Java runs the samples on -threads threads and adds each when it is done; the port adds
/// them in sample order whatever the threads (QS8: Java's key order can depend on which sample finishes
/// first). Returns the number of samples with a missing input.
pub fn prep_cohort_category(config: &QseeConfig, category: PrepCategory, sample_type: SampleType, matrix: &mut FeatureMatrix) -> Result<usize, PrepError> {
    let ids = config.sample_ids(sample_type).to_vec();
    let missing = std::sync::atomic::AtomicUsize::new(0);
    let n = ids.len();
    let total = n;
    let log_progress = |i: usize, id: &str| {
        if total > 100 {
            let interval = if total > 10000 { 10000 } else if total > 1000 { 1000 } else { 100 };
            if (i + 1) % interval == 0 || i == total - 1 {
                hmf_common::debug!("category({}) - Progress: {}/{} - current sample: {}", category.name(), i + 1, total, id);
            }
        }
    };
    let threads = (config.threads.max(1) as usize).min(n.max(1));
    let mut results: Vec<Option<Result<Vec<Feature>, PrepError>>> = (0..n).map(|_| None).collect();
    if threads <= 1 {
        for (i, id) in ids.iter().enumerate() {
            log_progress(i, id);
            results[i] = Some(run_category_task(config, category, id, sample_type, Some(&missing)));
            if matches!(results[i], Some(Err(_))) {
                break;
            }
        }
    } else {
        let next = std::sync::atomic::AtomicUsize::new(0);
        let slots: Vec<std::sync::Mutex<Option<Result<Vec<Feature>, PrepError>>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
        std::thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(|| loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= n {
                        break;
                    }
                    log_progress(i, &ids[i]);
                    let r = run_category_task(config, category, &ids[i], sample_type, Some(&missing));
                    *slots[i].lock().unwrap() = Some(r);
                });
            }
        });
        for (i, s) in slots.into_iter().enumerate() {
            results[i] = s.into_inner().unwrap();
        }
    }
    for (i, r) in results.into_iter().enumerate() {
        match r {
            Some(Ok(f)) => matrix.add_row(&ids[i], f).map_err(PrepError::Fatal)?,
            Some(Err(e)) => return Err(e),
            None => {}
        }
    }
    Ok(missing.into_inner())
}

/// TrainConfig on top of QseeConfig
#[derive(Clone, Debug)]
pub struct TrainConfig {
    pub num_percentiles: i32,
    pub percentile_increment: f64,
    pub write_cohort_features: bool,
}

fn create_transformer(t: &TrainConfig) -> Result<PercentileTransformer> {
    if t.percentile_increment.is_nan() { PercentileTransformer::with_num_percentiles(t.num_percentiles) } else { PercentileTransformer::with_increment(t.percentile_increment) }
}

/// CohortFeaturesWriter: `qsee.cohort.features.<tumor|normal>[.<id>].tsv.gz`, a row per feature and a
/// column per sample
struct CohortFeaturesWriter {
    writer: JavaWriter,
}

impl CohortFeaturesWriter {
    fn new(config: &QseeConfig, sample_type: SampleType) -> Result<CohortFeaturesWriter> {
        let path = generate_cohort_filename(&config.output_dir, &format!("features.{}", sample_type.name().to_lowercase()), config.output_id.as_deref(), "tsv.gz");
        let mut writer = JavaWriter::create(&path)?;
        let mut header = vec![COL_SOURCE_TOOL.to_string(), COL_FEATURE_TYPE.to_string(), COL_FEATURE_NAME.to_string()];
        header.extend(config.sample_ids(sample_type).iter().cloned());
        writer.write_all(header.join("\t").as_bytes())?;
        writer.write_all(b"\n")?;
        Ok(CohortFeaturesWriter { writer })
    }

    fn write_category(&mut self, m: &FeatureMatrix) -> Result<()> {
        for key in &m.keys {
            let mut line = vec![key.source_tool.map(|t| t.name()).unwrap_or("null").to_string(), key.feature_type.map(|t| t.name()).unwrap_or("null").to_string(), key.name.clone()];
            line.extend(m.column(key).iter().map(|f| decimal_8(f.value)));
            self.writer.write_all(line.join("\t").as_bytes())?;
            self.writer.write_all(b"\n")?;
        }
        Ok(())
    }
}

fn run_for(config: &QseeConfig, train: &TrainConfig, sample_type: SampleType, percentiles: &[f64]) -> Result<Vec<FeaturePercentiles>, PrepError> {
    let ids = config.sample_ids(sample_type);
    if ids.is_empty() {
        hmf_common::info!("Skipping sampleType({}) as no samples provided", sample_type.name());
        return Ok(Vec::new());
    }
    let categories = crate::app::selected_categories(config);
    let mut out = Vec::new();
    let mut writer = if train.write_cohort_features { Some(CohortFeaturesWriter::new(config, sample_type)?) } else { None };
    for category in categories {
        hmf_common::info!("Extracting cohort features - sampleType({}) category({})", sample_type.name(), category.name());
        let mut m = FeatureMatrix::with_samples(ids);
        let missing = prep_cohort_category(config, category, sample_type, &mut m)?;
        if missing == ids.len() {
            hmf_common::error!("failed prep as no samples had data for sampleType({}) category({})", sample_type.name(), category.name());
            return Err(PrepError::Exit);
        } else if missing > 0 {
            hmf_common::warn!("sampleType({}) category({}) - {}/{} samples had missing input files", sample_type.name(), category.name(), missing, ids.len());
        }
        m.sort_feature_keys();
        if let Some(w) = writer.as_mut() {
            hmf_common::info!("Writing cohort features - sampleType({}) category({})", sample_type.name(), category.name());
            w.write_category(&m)?;
        }
        hmf_common::info!("Calculating percentiles - sampleType({}) category({})", sample_type.name(), category.name());
        for key in &m.keys {
            let values: Vec<f64> = m.column(key).iter().map(|f| f.value).collect();
            let mut t = create_transformer(train)?;
            t.fit(&values)?;
            out.push(FeaturePercentiles { sample_type, key: key.clone(), percentiles: percentiles.to_vec(), ref_values: t.ref_values.clone() });
        }
    }
    if let Some(w) = writer {
        w.writer.finish()?;
    }
    Ok(out)
}

/// CohortPercentilesTrainer.run
pub fn train(config: &QseeConfig, train: &TrainConfig) -> Result<(), PrepError> {
    let percentiles = create_transformer(train)?.percentiles;
    let desc = if percentiles.len() <= 10 {
        percentiles.iter().map(|p| decimal_8(*p)).collect::<Vec<_>>().join(", ")
    } else {
        let n = percentiles.len();
        format!(
            "{}, {}, {}, ..., {}, {}, {}",
            decimal_8(percentiles[0]),
            decimal_8(percentiles[1]),
            decimal_8(percentiles[2]),
            decimal_8(percentiles[n - 3]),
            decimal_8(percentiles[n - 2]),
            decimal_8(percentiles[n - 1])
        )
    };
    hmf_common::info!("Running CohortPercentilesTrainer with percentiles: {desc}");
    let mut all = run_for(config, train, SampleType::Tumor, &percentiles)?;
    all.extend(run_for(config, train, SampleType::Normal, &percentiles)?);
    let path = generate_cohort_filename(&config.output_dir, "percentiles", config.output_id.as_deref(), "tsv.gz");
    hmf_common::info!("Writing cohort percentiles file: {path}");
    write_cohort_percentiles(&path, &percentiles, &all)?;
    Ok(())
}
