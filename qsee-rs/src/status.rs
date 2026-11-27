//! qsee `status`: ComparisonOperator, QcStatusType, QcStatus, ThresholdKey, QcThreshold,
//! ThresholdRegistry (the default thresholds), ThresholdOverridesFile, PurpleQCStatusConverter,
//! QcStatusFile and DefaultThresholdsWriter.

use anyhow::{anyhow, bail, Result};

use hmf_common::delim_file::{write_delim_file, DelimFileReader};
use hmf_common::fmt::java_double_to_string;

use crate::common::{decimal_4, SampleType, COL_FEATURE_NAME, COL_FEATURE_TYPE, COL_FEATURE_VALUE, COL_SAMPLE_ID, COL_SAMPLE_TYPE, COL_SOURCE_TOOL};
use crate::compat;
use crate::feature::{FeatureType, NumberFormat, SampleFeatures, SummaryTableFeature};

/// ComparisonOperator, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComparisonOperator {
    LessThan,
    LessThanOrEqual,
    GreaterThanOrEqual,
    GreaterThan,
}

impl ComparisonOperator {
    const ALL: [ComparisonOperator; 4] =
        [ComparisonOperator::LessThan, ComparisonOperator::LessThanOrEqual, ComparisonOperator::GreaterThanOrEqual, ComparisonOperator::GreaterThan];

    pub fn name(self) -> &'static str {
        match self {
            ComparisonOperator::LessThan => "LESS_THAN",
            ComparisonOperator::LessThanOrEqual => "LESS_THAN_OR_EQUAL",
            ComparisonOperator::GreaterThanOrEqual => "GREATER_THAN_OR_EQUAL",
            ComparisonOperator::GreaterThan => "GREATER_THAN",
        }
    }
    pub fn operator_string(self) -> &'static str {
        match self {
            ComparisonOperator::LessThan => "<",
            ComparisonOperator::LessThanOrEqual => "<=",
            ComparisonOperator::GreaterThanOrEqual => ">=",
            ComparisonOperator::GreaterThan => ">",
        }
    }
    /// ComparisonOperator.fromString: the name or the symbol
    pub fn from_string(s: &str) -> Result<ComparisonOperator> {
        Self::ALL
            .into_iter()
            .find(|o| o.name() == s || o.operator_string() == s)
            .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: Unknown comparison operator: {s}"))
    }
}

/// QcStatusType, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum QcStatusType {
    None,
    Pass,
    Warn,
    Fail,
}

impl QcStatusType {
    pub fn name(self) -> &'static str {
        match self {
            QcStatusType::None => "NONE",
            QcStatusType::Pass => "PASS",
            QcStatusType::Warn => "WARN",
            QcStatusType::Fail => "FAIL",
        }
    }
    /// displayString: "" for NONE
    pub fn display_string(self) -> &'static str {
        if self == QcStatusType::None { "" } else { self.name() }
    }
    pub fn value_of(s: &str) -> Result<QcStatusType> {
        [QcStatusType::None, QcStatusType::Pass, QcStatusType::Warn, QcStatusType::Fail]
            .into_iter()
            .find(|t| t.name() == s)
            .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.qsee.status.QcStatusType.{s}"))
    }
}

/// QcStatus: a status with the threshold it was judged against (for display)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QcStatus {
    pub status_type: QcStatusType,
    pub operator: Option<ComparisonOperator>,
    pub threshold: f64,
}

impl QcStatus {
    /// QcStatus.createEmpty
    pub fn empty() -> QcStatus {
        QcStatus { status_type: QcStatusType::None, operator: None, threshold: f64::NAN }
    }
    /// QcStatus.createPass: PASS with the same threshold
    pub fn pass(other: &QcStatus) -> QcStatus {
        QcStatus { status_type: QcStatusType::Pass, operator: other.operator, threshold: other.threshold }
    }
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
    pub fn display_string(&self, number_format: NumberFormat) -> Result<String> {
        try_form_display_string(self.status_type, self.operator, self.threshold, number_format)
    }
}

/// QcStatus.formDiplayString: "" for NONE or no threshold; else the operator and the threshold as
/// BigDecimal.valueOf(threshold) (times 100 and with "%" for percentages), stripTrailingZeros,
/// toPlainString. Java fails on an infinite threshold (BigDecimal.valueOf throws
/// NumberFormatException) and on a missing operator (NullPointerException) (QS6).
pub fn try_form_display_string(status_type: QcStatusType, operator: Option<ComparisonOperator>, threshold: f64, number_format: NumberFormat) -> Result<String> {
    if status_type == QcStatusType::None || threshold.is_nan() {
        return Ok(String::new());
    }
    let percent = number_format == NumberFormat::Percent;
    if threshold.is_infinite() {
        if compat::java("QS6") {
            bail!("java.lang.NumberFormatException: Infinite or NaN");
        }
        let op = operator.map(|o| o.operator_string()).unwrap_or("");
        let inf = if threshold > 0.0 { "Inf" } else { "-Inf" };
        return Ok(format!("{op}{inf}{}", if percent { "%" } else { "" }));
    }
    let Some(op) = operator else {
        bail!("java.lang.NullPointerException: no comparison operator for threshold {threshold}");
    };
    let mut value = BigDecimal::value_of(threshold);
    if percent {
        value = value.multiply_100();
    }
    let mut s = value.strip_trailing_zeros().to_plain_string();
    if percent {
        s.push('%');
    }
    Ok(format!("{}{s}", op.operator_string()))
}

/// The part of java.math.BigDecimal that formDiplayString uses: an unscaled decimal integer (digits and a
/// sign) and a scale
#[derive(Clone, Debug, PartialEq)]
struct BigDecimal {
    negative: bool,
    digits: Vec<u8>,
    scale: i64,
}

impl BigDecimal {
    /// BigDecimal.valueOf(double): new BigDecimal(Double.toString(d))
    fn value_of(d: f64) -> BigDecimal {
        let s = java_double_to_string(d);
        let (negative, s) = match s.strip_prefix('-') {
            Some(r) => (true, r.to_string()),
            None => (false, s),
        };
        let (mantissa, exponent) = match s.split_once('E') {
            Some((m, e)) => (m.to_string(), e.parse::<i64>().unwrap()),
            None => (s.clone(), 0),
        };
        let (int_part, frac_part) = mantissa.split_once('.').unwrap_or((&mantissa, ""));
        let mut digits: Vec<u8> = int_part.bytes().chain(frac_part.bytes()).map(|b| b - b'0').collect();
        while digits.len() > 1 && digits[0] == 0 {
            digits.remove(0);
        }
        BigDecimal { negative, digits, scale: frac_part.len() as i64 - exponent }
    }

    /// multiply(BigDecimal.valueOf(100)): the unscaled value times 100, the same scale
    fn multiply_100(mut self) -> BigDecimal {
        if !(self.digits.len() == 1 && self.digits[0] == 0) {
            self.digits.push(0);
            self.digits.push(0);
        }
        self
    }

    /// stripTrailingZeros: a zero value becomes 0 with scale 0
    fn strip_trailing_zeros(mut self) -> BigDecimal {
        if self.digits.iter().all(|d| *d == 0) {
            return BigDecimal { negative: false, digits: vec![0], scale: 0 };
        }
        while self.digits.len() > 1 && *self.digits.last().unwrap() == 0 {
            self.digits.pop();
            self.scale -= 1;
        }
        self
    }

    /// toPlainString
    fn to_plain_string(&self) -> String {
        let mut s = String::new();
        if self.negative {
            s.push('-');
        }
        let digits: String = self.digits.iter().map(|d| (d + b'0') as char).collect();
        if self.scale <= 0 {
            s.push_str(&digits);
            s.extend(std::iter::repeat_n('0', (-self.scale) as usize));
        } else {
            let scale = self.scale as usize;
            if digits.len() > scale {
                let (a, b) = digits.split_at(digits.len() - scale);
                s.push_str(a);
                s.push('.');
                s.push_str(b);
            } else {
                s.push_str("0.");
                s.extend(std::iter::repeat_n('0', scale - digits.len()));
                s.push_str(&digits);
            }
        }
        s
    }
}

/// ThresholdKey
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ThresholdKey {
    pub sample_type: SampleType,
    pub feature_type: FeatureType,
    pub feature_name: Option<String>,
    pub qc_status_type: QcStatusType,
}

impl ThresholdKey {
    pub fn new(sample_type: SampleType, feature_type: FeatureType, feature_name: &str, qc_status_type: QcStatusType) -> ThresholdKey {
        ThresholdKey { sample_type, feature_type, feature_name: Some(feature_name.to_string()), qc_status_type }
    }
    /// toString
    pub fn display(&self) -> String {
        format!(
            "sampleType({}) featureType({}) featureName({}) qcStatusType({})",
            self.sample_type.name(),
            self.feature_type.name(),
            self.feature_name.as_deref().unwrap_or("null"),
            self.qc_status_type.name()
        )
    }
}

/// QcThreshold
#[derive(Clone, Debug)]
pub struct QcThreshold {
    pub key: ThresholdKey,
    pub operator: Option<ComparisonOperator>,
    pub threshold: f64,
    pub determined_elsewhere: bool,
}

impl PartialEq for QcThreshold {
    /// equals: the key, the operator and the threshold (NaN equal to NaN)
    fn eq(&self, other: &QcThreshold) -> bool {
        let values_equal = self.threshold == other.threshold || (self.threshold.is_nan() && other.threshold.is_nan());
        self.key == other.key && self.operator == other.operator && values_equal
    }
}

impl QcThreshold {
    /// getQcStatus: the key's status type with this threshold
    pub fn qc_status(&self) -> QcStatus {
        QcStatus { status_type: self.key.qc_status_type, operator: self.operator, threshold: self.threshold }
    }

    /// getSampleQcStatus: the key's status if the value meets the condition, else PASS (with the
    /// threshold); a NaN threshold never meets it. Java has no operator to switch on for a threshold
    /// without one (NullPointerException); no feature asks for those.
    pub fn sample_qc_status(&self, value: f64) -> QcStatus {
        if self.key.qc_status_type == QcStatusType::None {
            return QcStatus::empty();
        }
        let fails = match self.operator.expect("QcThreshold without a comparison operator") {
            ComparisonOperator::LessThan => value < self.threshold,
            ComparisonOperator::LessThanOrEqual => value <= self.threshold,
            ComparisonOperator::GreaterThanOrEqual => value >= self.threshold,
            ComparisonOperator::GreaterThan => value > self.threshold,
        };
        if fails {
            QcStatus { status_type: self.key.qc_status_type, operator: self.operator, threshold: self.threshold }
        } else {
            QcStatus { status_type: QcStatusType::Pass, operator: self.operator, threshold: self.threshold }
        }
    }
}

// PurpleQCStatus thresholds (hmf-common PurpleQCStatus)
pub const PURPLE_MIN_PURITY: f64 = 0.2;
pub const PURPLE_MAX_DELETED_GENES: f64 = 280.0;
pub const PURPLE_MAX_UNSUPPORTED_SEGMENTS: f64 = 220.0;
pub const PURPLE_CONTAMINATION_WARN_LEVEL: f64 = 0.02;
pub const PURPLE_CONTAMINATION_FAIL_LEVEL: f64 = 0.1;
pub const PURPLE_TINC_WARN_LEVEL: f64 = 0.02;
pub const PURPLE_TINC_FAIL_LEVEL: f64 = 0.25;

/// ThresholdRegistry: the thresholds in insertion order (a LinkedHashMap: a key set again keeps its place)
#[derive(Clone, Debug, Default)]
pub struct ThresholdRegistry {
    thresholds: Vec<QcThreshold>,
}

impl ThresholdRegistry {
    /// createDefault(targetedMode)
    pub fn create_default(targeted_mode: bool) -> ThresholdRegistry {
        let mut r = ThresholdRegistry::default();
        r.set_defaults(targeted_mode);
        r
    }

    /// createWithoutThresholds: the defaults (WGS) with every threshold NaN
    #[allow(dead_code)] // ported from Java; used by the tests or kept as Java has it
        pub fn create_without_thresholds() -> ThresholdRegistry {
        let mut r = ThresholdRegistry::default();
        r.set_defaults(false);
        for t in &mut r.thresholds {
            t.threshold = f64::NAN;
        }
        r
    }

    fn put(&mut self, t: QcThreshold) {
        match self.thresholds.iter_mut().find(|x| x.key == t.key) {
            Some(x) => *x = t,
            None => self.thresholds.push(t),
        }
    }

    fn set(&mut self, sample_type: SampleType, name: &str, qc: QcStatusType, op: Option<ComparisonOperator>, value: f64, elsewhere: bool) {
        self.put(QcThreshold { key: ThresholdKey::new(sample_type, FeatureType::SummaryTable, name, qc), operator: op, threshold: value, determined_elsewhere: elsewhere });
    }

    fn set_common(&mut self, name: &str, qc: QcStatusType, op: ComparisonOperator, value: f64) {
        self.set(SampleType::Tumor, name, qc, Some(op), value, false);
        self.set(SampleType::Normal, name, qc, Some(op), value, false);
    }

    /// setDefaults
    pub fn set_defaults(&mut self, targeted: bool) {
        use ComparisonOperator::*;
        use QcStatusType::{Fail, Warn};
        use SampleType::{Normal, Tumor};
        use SummaryTableFeature as F;
        hmf_common::info!("Using default QC thresholds for {} mode", if targeted { "targeted" } else { "WGS" });
        let choose = |wgs: f64, tgt: f64| if targeted { tgt } else { wgs };
        let nan = f64::NAN;
        self.set_common(F::MappedProportion.name(), Fail, LessThan, 0.95);
        self.set_common(F::LowMapQual.name(), Warn, GreaterThan, 0.05);
        self.set_common(F::LowBaseQual.name(), Warn, GreaterThan, 0.05);
        self.set_common(F::DuplicateReads.name(), Warn, GreaterThan, choose(0.3, 0.99));
        self.set_common(F::DualStrandReads.name(), Warn, GreaterThan, choose(0.5, nan));

        // PURPLE QC thresholds are not handled by qsee; the contamination operators follow PURPLE's (QS4)
        let contamination_op = if compat::java("QS4") { GreaterThan } else { GreaterThanOrEqual };
        self.set(Tumor, F::Purity.name(), Warn, Some(LessThan), PURPLE_MIN_PURITY, true);
        self.set(Tumor, F::DeletedGenes.name(), Warn, Some(GreaterThan), PURPLE_MAX_DELETED_GENES, true);
        self.set(Tumor, F::UnsupportedCnSegments.name(), Warn, Some(GreaterThan), PURPLE_MAX_UNSUPPORTED_SEGMENTS, true);
        self.set(Tumor, F::Contamination.name(), Warn, Some(contamination_op), PURPLE_CONTAMINATION_WARN_LEVEL, true);
        self.set(Tumor, F::Contamination.name(), Fail, Some(contamination_op), PURPLE_CONTAMINATION_FAIL_LEVEL, true);
        self.set(Tumor, F::Tinc.name(), Warn, Some(GreaterThan), PURPLE_TINC_WARN_LEVEL, true);
        self.set(Tumor, F::Tinc.name(), Fail, Some(GreaterThanOrEqual), PURPLE_TINC_FAIL_LEVEL, true);
        self.set(Tumor, "FAIL_NO_TUMOR", Fail, None, nan, true);
        self.set(Tumor, "WARN_GENDER_MISMATCH", Warn, None, nan, true);

        self.set(Tumor, F::MeanCoverage.name(), Warn, Some(LessThan), 70.0, false);
        self.set(Tumor, F::CoverageAbove10.name(), Warn, Some(LessThan), choose(0.9, 0.99), false);
        self.set(Tumor, F::CoverageAbove20.name(), Warn, Some(LessThan), choose(0.9, 0.99), false);
        self.set(Tumor, F::CoverageAbove30.name(), Warn, Some(LessThan), choose(0.9, 0.99), false);
        self.set(Tumor, F::CoverageAbove60.name(), Warn, Some(LessThan), choose(0.8, 0.99), false);
        self.set(Tumor, F::CoverageAbove100.name(), Warn, Some(LessThan), choose(0.1, 0.99), false);
        self.set(Tumor, F::CoverageAbove250.name(), Warn, Some(LessThan), choose(nan, 0.99), false);

        self.set(Normal, F::MeanCoverage.name(), Warn, Some(LessThan), 20.0, false);
        self.set(Normal, F::CoverageAbove10.name(), Warn, Some(LessThan), 0.9, false);
        self.set(Normal, F::CoverageAbove20.name(), Warn, Some(LessThan), 0.8, false);
        self.set(Normal, F::CoverageAbove30.name(), Warn, Some(LessThan), 0.4, false);
        self.set(Normal, F::CoverageAbove60.name(), Warn, Some(LessThan), nan, false);
        self.set(Normal, F::CoverageAbove100.name(), Warn, Some(LessThan), nan, false);
        self.set(Normal, F::CoverageAbove250.name(), Warn, Some(LessThan), nan, false);
    }

    pub fn contains_key(&self, key: &ThresholdKey) -> bool {
        self.thresholds.iter().any(|t| &t.key == key)
    }

    /// getThreshold: an error (NoSuchElementException) for a key not registered
    pub fn get(&self, key: &ThresholdKey) -> Result<&QcThreshold> {
        self.thresholds.iter().find(|t| &t.key == key).ok_or_else(|| anyhow!("java.util.NoSuchElementException: No threshold defined for {}", key.display()))
    }

    pub fn get_summary(&self, sample_type: SampleType, feature: SummaryTableFeature, qc: QcStatusType) -> Result<&QcThreshold> {
        self.get(&ThresholdKey::new(sample_type, FeatureType::SummaryTable, feature.name(), qc))
    }

    pub fn get_all(&self) -> &[QcThreshold] {
        &self.thresholds
    }
}

// ---- ThresholdOverridesFile ----

/// ThresholdOverridesFile.read: the defaults, with each row's threshold put over them. A key that is not
/// a default threshold is logged; the rows after it are not applied and, once all are read, the run ends
/// (RuntimeException). A key given twice keeps the first.
pub fn read_threshold_overrides(path: &str, targeted_mode: bool) -> Result<ThresholdRegistry> {
    let reader = DelimFileReader::open(path)?;
    hmf_common::info!("Reading threshold overrides from: {path}");
    let mut thresholds = ThresholdRegistry::default();
    thresholds.set_defaults(targeted_mode);
    let mut unsupported = false;
    let mut seen: Vec<ThresholdKey> = Vec::new();
    for row in reader.rows() {
        let sample_type = SampleType::value_of(row.get(COL_SAMPLE_TYPE)?)?;
        let feature_type = FeatureType::value_of(row.get(COL_FEATURE_TYPE)?)?;
        let feature_name = row.get(COL_FEATURE_NAME)?;
        let qc = QcStatusType::value_of(row.get("QcStatusType")?)?;
        let op = ComparisonOperator::from_string(row.get("ComparisonOperator")?)?;
        let value = row.get_double("Threshold")?;
        let key = ThresholdKey::new(sample_type, feature_type, feature_name, qc);
        if !thresholds.contains_key(&key) {
            unsupported = true;
            hmf_common::error!("  Unsupported threshold: {}", key.display());
        }
        if unsupported {
            continue;
        }
        if seen.contains(&key) {
            hmf_common::warn!("  Skipping duplicate threshold: {}", key.display());
            continue;
        }
        let existing = thresholds.get(&key)?.clone();
        // QS9: an override of a threshold PURPLE decides changes only the threshold qsee displays
        if existing.determined_elsewhere && !compat::java("QS9") {
            hmf_common::warn!("  Ignoring threshold override determined by PURPLE: {}", key.display());
            seen.push(key);
            continue;
        }
        hmf_common::debug!("  Overriding threshold: {} default({}) override({})", key.display(), java_double_to_string(existing.threshold), java_double_to_string(value));
        thresholds.put(QcThreshold { key: key.clone(), operator: Some(op), threshold: value, determined_elsewhere: false });
        seen.push(key);
    }
    if unsupported {
        bail!("java.lang.RuntimeException: Failed to read threshold overrides file - unsupported threshold entries found");
    }
    Ok(thresholds)
}

/// ThresholdOverridesFile.write (DefaultThresholdsWriter's file)
pub fn write_thresholds(path: &str, thresholds: &[QcThreshold]) -> Result<()> {
    let columns: Vec<String> =
        [COL_SAMPLE_TYPE, COL_FEATURE_TYPE, COL_FEATURE_NAME, "QcStatusType", "ComparisonOperator", "Threshold", "Overridable"].iter().map(|s| s.to_string()).collect();
    write_delim_file(path, &columns, thresholds, |t, row| {
        row.set(COL_SAMPLE_TYPE, t.key.sample_type.name())?;
        row.set(COL_FEATURE_TYPE, t.key.feature_type.name())?;
        match &t.key.feature_name {
            Some(n) => row.set(COL_FEATURE_NAME, n.as_str())?,
            None => row.set_null(COL_FEATURE_NAME)?,
        }
        row.set("QcStatusType", t.key.qc_status_type.name())?;
        let op = if t.determined_elsewhere {
            "NA".to_string()
        } else {
            t.operator.ok_or_else(|| anyhow!("java.lang.NullPointerException: no comparison operator"))?.name().to_string()
        };
        row.set("ComparisonOperator", op)?;
        row.set("Threshold", decimal_4(t.threshold))?;
        row.set("Overridable", if t.determined_elsewhere { "false" } else { "true" })
    })
}

/// DefaultThresholdsWriter.main: `qsee.thresholds.default.tsv` in -output_dir (or the working directory)
pub fn write_default_thresholds(output_dir: &str, overridable_only: bool, targeted_mode: bool) -> Result<()> {
    let path = format!("{}qsee.thresholds.default.tsv", hmf_common::metrics::check_add_dir_separator(output_dir));
    hmf_common::info!("Writing default thresholds to: {path}");
    let registry = ThresholdRegistry::create_default(targeted_mode);
    let mut list: Vec<QcThreshold> = registry.get_all().to_vec();
    // a stable sort by sample type
    list.sort_by_key(|t| t.key.sample_type);
    if overridable_only {
        hmf_common::info!("Skipping thresholds that cannot be overridden");
        list.retain(|t| !t.determined_elsewhere);
    }
    write_thresholds(&path, &list)
}

// ---- PurpleQCStatusConverter ----

/// PurpleQCStatus, in enum order
pub const PURPLE_QC_STATUSES: [&str; 10] = hmf_common::purple_files::PURPLE_QC_STATUSES;

/// PurpleQCStatusConverter.toQseeQcStatus
pub fn purple_to_qsee_status(thresholds: &ThresholdRegistry, purple_status: &str) -> Result<QcStatus> {
    use QcStatusType::{Fail, Warn};
    let get = |name: &str, qc: QcStatusType| -> Result<QcStatus> {
        Ok(thresholds.get(&ThresholdKey::new(SampleType::Tumor, FeatureType::SummaryTable, name, qc))?.qc_status())
    };
    match purple_status {
        "PASS" => Ok(QcStatus::empty()),
        "WARN_DELETED_GENES" => get("DELETED_GENES", Warn),
        "WARN_HIGH_COPY_NUMBER_NOISE" => get("UNSUPPORTED_CN_SEGMENTS", Warn),
        "WARN_LOW_PURITY" => get("PURITY", Warn),
        "WARN_TINC" => get("TINC", Warn),
        "FAIL_TINC" => get("TINC", Fail),
        "WARN_CONTAMINATION" => get("CONTAMINATION", Warn),
        "FAIL_CONTAMINATION" => get("CONTAMINATION", Fail),
        "FAIL_NO_TUMOR" => get("FAIL_NO_TUMOR", Fail),
        "WARN_GENDER_MISMATCH" => get("WARN_GENDER_MISMATCH", Warn),
        other => bail!("java.lang.IllegalArgumentException: QC threshold not defined for PurpleQCStatus({other})"),
    }
}

// ---- QcStatusFile ----

/// QcStatusFile.write: a row for every feature with a QC status (sample order, then feature order)
pub fn write_qc_status_file(path: &str, samples: &[SampleFeatures]) -> Result<()> {
    let mut rows = Vec::new();
    for s in samples {
        for f in &s.features {
            if f.metadata.qc_status.status_type != QcStatusType::None {
                rows.push((s, f));
            }
        }
    }
    let columns: Vec<String> = [COL_SAMPLE_ID, COL_SAMPLE_TYPE, COL_SOURCE_TOOL, COL_FEATURE_TYPE, COL_FEATURE_NAME, COL_FEATURE_VALUE, "QcStatus", "FailCondition"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    write_delim_file(path, &columns, &rows, |(s, f), row| {
        row.set(COL_SAMPLE_ID, s.sample_id.as_str())?;
        row.set(COL_SAMPLE_TYPE, s.sample_type.name())?;
        row.set(COL_SOURCE_TOOL, f.key.source_tool.map(|t| t.name()).unwrap_or("null"))?;
        row.set(COL_FEATURE_TYPE, f.key.feature_type.map(|t| t.name()).unwrap_or("null"))?;
        row.set(COL_FEATURE_NAME, f.key.name.as_str())?;
        row.set(COL_FEATURE_VALUE, decimal_4(f.value))?;
        let q = &f.metadata.qc_status;
        row.set("QcStatus", q.status_type.name())?;
        row.set("FailCondition", try_form_display_string(q.status_type, q.operator, q.threshold, NumberFormat::Number)?)
    })
}
