//! qsee `common`: sample types, file names and number formats (QseeFileCommon), BinnedFrequencies,
//! LinearInterpolation, MultiFieldStringBuilder and SampleIdsLoader.

use anyhow::{anyhow, bail, Result};

use hmf_common::fmt::decimal_format;
use hmf_common::java_text;
use hmf_common::metrics::check_add_dir_separator;

/// SampleType, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SampleType {
    Tumor,
    Normal,
}

impl SampleType {
    pub const ALL: [SampleType; 2] = [SampleType::Tumor, SampleType::Normal];

    pub fn name(self) -> &'static str {
        match self {
            SampleType::Tumor => "TUMOR",
            SampleType::Normal => "NORMAL",
        }
    }
    /// SampleType.valueOf
    pub fn value_of(s: &str) -> Result<SampleType> {
        match s {
            "TUMOR" => Ok(SampleType::Tumor),
            "NORMAL" => Ok(SampleType::Normal),
            _ => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.qsee.common.SampleType.{s}"),
        }
    }
}

pub const QSEE_FILE_ID: &str = "qsee";
pub const MULTISAMPLE_SAMPLE_ID: &str = "multisample";

pub const COL_SAMPLE_ID: &str = "SampleId";
pub const COL_SAMPLE_TYPE: &str = "SampleType";
pub const COL_FEATURE_TYPE: &str = "FeatureType";
pub const COL_FEATURE_NAME: &str = "FeatureName";
pub const COL_SOURCE_TOOL: &str = "SourceTool";
pub const COL_FEATURE_VALUE: &str = "FeatureValue";

/// QseeFileCommon.DECIMAL_FORMAT: DecimalFormat("0.########") with English symbols and "Inf" for infinity
pub fn decimal_8(v: f64) -> String {
    if v.is_infinite() {
        return if v > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    decimal_format(v, 8)
}

/// DelimFileWriter's default number format for a double: DecimalFormat("#.####") with English symbols
pub fn decimal_4(v: f64) -> String {
    decimal_format(v, 4)
}

/// QseeFileCommon.generateFilename: <dir>/<sample>.qsee.<fileId>[.<outputId>].<extension>
pub fn generate_filename(base_path: &str, sample_id: &str, file_id: &str, output_id: Option<&str>, extension: &str) -> String {
    let mut parts = vec![sample_id, QSEE_FILE_ID, file_id];
    if let Some(o) = output_id {
        parts.push(o);
    }
    parts.push(extension);
    format!("{}{}", check_add_dir_separator(base_path), parts.join("."))
}

/// QseeFileCommon.generateCohortFilename: <dir>/qsee.cohort.<fileId>[.<outputId>].<extension>
pub fn generate_cohort_filename(base_path: &str, file_id: &str, output_id: Option<&str>, extension: &str) -> String {
    let mut parts = vec![QSEE_FILE_ID, "cohort", file_id];
    if let Some(o) = output_id {
        parts.push(o);
    }
    parts.push(extension);
    format!("{}{}", check_add_dir_separator(base_path), parts.join("."))
}

// ---- MultiFieldStringBuilder ----

pub const FIELD_KEY_VALUE_SEPARATOR: &str = ":";
pub const FIELD_SEPARATOR: &str = ";";

/// MultiFieldStringBuilder.formSingleField
pub fn form_single_field(key: &str, value: &str) -> String {
    format!("{key}{FIELD_KEY_VALUE_SEPARATOR}{value}")
}

/// MultiFieldStringBuilder.formMultiField / toString: key:value pairs joined by ';', pairs with an empty
/// value left out
pub fn form_multi_field(pairs: &[(&str, &str)]) -> String {
    pairs.iter().filter(|(_, v)| !v.is_empty()).map(|(k, v)| form_single_field(k, v)).collect::<Vec<_>>().join(FIELD_SEPARATOR)
}

/// MultiFieldStringBuilder.parseMultiFieldAsPairs: String.split(";") then each field split(":"); a field
/// without a value is an error (Java's ArrayIndexOutOfBoundsException)
pub fn parse_multi_field_as_pairs(fields: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for field in java_text::split(fields, ';') {
        let kv = java_text::split(field, ':');
        if kv.len() < 2 {
            bail!("java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length {}", kv.len());
        }
        out.push((kv[0].to_string(), kv[1].to_string()));
    }
    Ok(out)
}

/// MultiFieldStringBuilder.parseMultiFieldAsArray: String.split(";|:")
pub fn parse_multi_field_as_array(fields: &str) -> Vec<String> {
    if fields.is_empty() {
        return vec![String::new()];
    }
    let mut parts: Vec<String> = fields.split([';', ':']).map(String::from).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

// ---- BinnedFrequencies ----

/// BinnedFrequencies: bins from their starts (right-open, the last bin open-ended) and their counts
#[derive(Clone, Debug)]
pub struct BinnedFrequencies {
    pub bin_starts: Vec<i64>,
    pub bin_ends: Vec<i64>,
    pub frequencies: Vec<f64>,
}

impl BinnedFrequencies {
    /// new BinnedFrequencies(binStarts, frequencies). With no bins Java's System.arraycopy is given a
    /// length of -1 and throws; that is an error here too.
    pub fn new(bin_starts: Vec<i64>, frequencies: Vec<f64>) -> Result<BinnedFrequencies> {
        if bin_starts.is_empty() {
            bail!("java.lang.ArrayIndexOutOfBoundsException: arraycopy: length -1 is negative");
        }
        let mut bin_ends: Vec<i64> = bin_starts[1..].to_vec();
        bin_ends.push(i64::MAX);
        Ok(BinnedFrequencies { bin_starts, bin_ends, frequencies })
    }

    /// BinnedFrequencies.fromValueFrequencies
    pub fn from_value_frequencies(values: &[hmf_common::metrics::ValueFrequency]) -> Result<BinnedFrequencies> {
        Self::new(values.iter().map(|v| v.value as i64).collect(), values.iter().map(|v| v.count as f64).collect())
    }

    /// calcFrequencyDensities: each count divided by its bin's width, the last bin's count as it is
    pub fn frequency_densities(&self) -> Vec<f64> {
        let n = self.frequencies.len();
        (0..n)
            .map(|i| {
                if i == n - 1 {
                    self.frequencies[i]
                } else {
                    let width = self.bin_ends[i].wrapping_sub(self.bin_starts[i]) as f64;
                    self.frequencies[i] / width
                }
            })
            .collect()
    }

    /// calcProportionalDensities: the densities divided by the total count (Arrays.stream(..).sum())
    pub fn proportional_densities(&self) -> Vec<f64> {
        let total = hmf_common::java_math::double_stream_sum(self.frequencies.iter().copied());
        self.frequency_densities().into_iter().map(|d| d / total).collect()
    }

    /// formProportionalDensityFeatures: "<binTitle>:<bin start>" and the proportional density
    pub fn proportional_density_features(&self, bin_title: &str) -> Vec<(String, f64)> {
        self.proportional_densities().into_iter().enumerate().map(|(i, d)| (form_single_field(bin_title, &self.bin_starts[i].to_string()), d)).collect()
    }
}

// ---- LinearInterpolation ----

/// LinearInterpolation.interpolate1D
pub fn interpolate_1d(lower: f64, upper: f64, fraction: f64) -> f64 {
    lower + fraction * (upper - lower)
}

#[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
/// Arrays.binarySearch(double[], double): Java's search, with Double.doubleToLongBits for equal values
pub fn java_binary_search(a: &[f64], key: f64) -> i64 {
    let (mut low, mut high) = (0i64, a.len() as i64 - 1);
    while low <= high {
        let mid = ((low + high) as u64 >> 1) as i64;
        let mid_val = a[mid as usize];
        if mid_val < key {
            low = mid + 1;
        } else if mid_val > key {
            high = mid - 1;
        } else {
            let mid_bits = hmf_common::java_math::double_to_long_bits(mid_val);
            let key_bits = hmf_common::java_math::double_to_long_bits(key);
            if mid_bits == key_bits {
                return mid;
            } else if mid_bits < key_bits {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
    }
    -(low + 1)
}

#[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
/// LinearInterpolation.interpolate2D: -Infinity below the first x, +Infinity above the last, the y of an
/// equal x, else linear between the neighbours
pub fn interpolate_2d(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    if x < xs[0] {
        return f64::NEG_INFINITY;
    }
    if x > xs[xs.len() - 1] {
        return f64::INFINITY;
    }
    let m = java_binary_search(xs, x);
    if m >= 0 {
        return ys[m as usize];
    }
    let upper = (-m - 1) as usize;
    let lower = upper - 1;
    let (lx, ux, ly, uy) = (xs[lower], xs[upper], ys[lower], ys[upper]);
    if lx == x {
        return ys[lower];
    }
    let fraction = (x - lx) / (ux - lx);
    interpolate_1d(ly, uy, fraction)
}

// ---- SampleIdsLoader ----

pub const SUPPORTED_SAMPLE_MODES_MSG: &str = "Only tumor-only and tumor-normal sample modes are supported";

/// What ends a run in SampleIdsLoader: an exception (a crash in Java, exit 1 with a stack trace) or a
/// logged error with System.exit(1); both give exit status 1.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct SampleIds {
    pub tumor_ids: Vec<String>,
    pub reference_ids: Vec<String>,
}

impl SampleIds {
    /// SampleIdsLoader.fromStrings
    pub fn from_strings(tumor: Option<&str>, reference: Option<&str>) -> Result<SampleIds> {
        let Some(tumor) = tumor else {
            bail!("java.lang.IllegalStateException: Tumor sample IDs must be provided. {SUPPORTED_SAMPLE_MODES_MSG}");
        };
        let mut ids = SampleIds { tumor_ids: java_text::split(tumor, ',').into_iter().map(String::from).collect(), reference_ids: Vec::new() };
        if let Some(r) = reference {
            ids.reference_ids = java_text::split(r, ',').into_iter().map(String::from).collect();
        }
        ids.check_matching_counts()?;
        Ok(ids)
    }

    /// SampleIdsLoader.fromLines: columns TumorId and, optionally, ReferenceId
    pub fn from_lines(lines: &[String]) -> Result<SampleIds> {
        let header = lines.first().ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
        let f = java_text::fields_index_map(header, '\t');
        let Some(&tumor_i) = f.get("TumorId") else {
            bail!("java.lang.IllegalStateException: Sample IDs file must contain column TumorId. {SUPPORTED_SAMPLE_MODES_MSG}");
        };
        let ref_i = f.get("ReferenceId").copied();
        if lines.len() == 1 {
            bail!("No sample IDs found in sample IDs file");
        }
        let mut ids = SampleIds::default();
        for line in &lines[1..] {
            let v = java_text::split_all(line, '\t');
            let get = |i: usize| -> Result<&str> { v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", v.len())) };
            let tumor = get(tumor_i)?;
            let reference = match ref_i {
                Some(i) => get(i)?,
                None => "",
            };
            if tumor.is_empty() && !reference.is_empty() {
                bail!("java.lang.IllegalStateException: Reference sample({reference}) has no tumor sample ID. {SUPPORTED_SAMPLE_MODES_MSG}");
            }
            if !tumor.is_empty() {
                ids.tumor_ids.push(tumor.to_string());
            }
            if !reference.is_empty() {
                ids.reference_ids.push(reference.to_string());
            }
        }
        ids.check_matching_counts()?;
        hmf_common::debug!("Loaded {} tumor and {} reference sample IDs from file", ids.tumor_ids.len(), ids.reference_ids.len());
        Ok(ids)
    }

    fn check_matching_counts(&self) -> Result<()> {
        if !self.reference_ids.is_empty() && self.reference_ids.len() != self.tumor_ids.len() {
            bail!("java.lang.IllegalStateException: If reference samples provided, no. of reference samples must match no. of tumor samples");
        }
        Ok(())
    }

    /// SampleIdsLoader.fromConfig
    pub fn from_config(cb: &hmf_common::config_builder::ConfigBuilder) -> Result<SampleIds> {
        if !cb.has_value("sample_id_file") && !cb.has_value("tumor") && !cb.has_value("reference") {
            bail!("Sample IDs must be provided to 1) -sample_id_file or to 2) -tumor and/or -reference");
        }
        if cb.has_value("sample_id_file") {
            let path = cb.value("sample_id_file").unwrap();
            let lines = java_text::read_lines(&path).map_err(|e| anyhow!("Failed to load sample IDs file: {e}"))?;
            Self::from_lines(&lines)
        } else {
            Self::from_strings(cb.value("tumor").as_deref(), cb.value("reference").as_deref())
        }
    }

    pub fn get(&self, sample_type: SampleType) -> &[String] {
        match sample_type {
            SampleType::Tumor => &self.tumor_ids,
            SampleType::Normal => &self.reference_ids,
        }
    }

    /// QseeConfig.isSinglePatient
    pub fn is_single_patient(&self) -> bool {
        self.tumor_ids.len() <= 1 && self.reference_ids.len() <= 1
    }
}
