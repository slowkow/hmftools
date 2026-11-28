//! qsee `feature`: FeatureType, SourceTool, NumberFormat, FeatureKey (with its Comparable order),
//! Feature, FeatureMetadata, SampleFeatures, and the category enums the key order depends on
//! (SummaryTableFeature, DiscordantFragGroup, BaseQualBin, RepeatType, ConsensusType).

use std::cmp::Ordering;

use hmf_common::java_text;
use hmf_common::redux_files::ConsensusType;

use crate::common::{parse_multi_field_as_array, parse_multi_field_as_pairs, SampleType, FIELD_KEY_VALUE_SEPARATOR};
use crate::status::QcStatus;

/// FeatureType, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FeatureType {
    SummaryTable,
    CoverageDistribution,
    FragLengthDistribution,
    GcBias,
    DiscordantFragFreq,
    DuplicateFreq,
    MissedVariantLikelihood,
    /// upstream bdceb2302a: the key of the per-exon LoD VAF threshold (no feature has this type)
    ExonVariantLod,
    BqrPerSnv96Context,
    BqrPerOrigQual,
    MsIndelErrorRates,
    MsIndelErrorBias,
}

impl FeatureType {
    pub const ALL: [FeatureType; 12] = [
        FeatureType::SummaryTable,
        FeatureType::CoverageDistribution,
        FeatureType::FragLengthDistribution,
        FeatureType::GcBias,
        FeatureType::DiscordantFragFreq,
        FeatureType::DuplicateFreq,
        FeatureType::MissedVariantLikelihood,
        FeatureType::ExonVariantLod,
        FeatureType::BqrPerSnv96Context,
        FeatureType::BqrPerOrigQual,
        FeatureType::MsIndelErrorRates,
        FeatureType::MsIndelErrorBias,
    ];

    pub fn name(self) -> &'static str {
        match self {
            FeatureType::SummaryTable => "SUMMARY_TABLE",
            FeatureType::CoverageDistribution => "COVERAGE_DISTRIBUTION",
            FeatureType::FragLengthDistribution => "FRAG_LENGTH_DISTRIBUTION",
            FeatureType::GcBias => "GC_BIAS",
            FeatureType::DiscordantFragFreq => "DISCORDANT_FRAG_FREQ",
            FeatureType::DuplicateFreq => "DUPLICATE_FREQ",
            FeatureType::MissedVariantLikelihood => "MISSED_VARIANT_LIKELIHOOD",
            FeatureType::ExonVariantLod => "EXON_VARIANT_LOD",
            FeatureType::BqrPerSnv96Context => "BQR_PER_SNV96_CONTEXT",
            FeatureType::BqrPerOrigQual => "BQR_PER_ORIG_QUAL",
            FeatureType::MsIndelErrorRates => "MS_INDEL_ERROR_RATES",
            FeatureType::MsIndelErrorBias => "MS_INDEL_ERROR_BIAS",
        }
    }
    /// FeatureType.valueOf
    pub fn value_of(s: &str) -> anyhow::Result<FeatureType> {
        FeatureType::ALL
            .into_iter()
            .find(|t| t.name() == s)
            .ok_or_else(|| anyhow::anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.qsee.feature.FeatureType.{s}"))
    }
}

/// SourceTool, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SourceTool {
    BamMetrics,
    Cobalt,
    Esvee,
    Purple,
    Redux,
}

impl SourceTool {
    pub const ALL: [SourceTool; 5] = [SourceTool::BamMetrics, SourceTool::Cobalt, SourceTool::Esvee, SourceTool::Purple, SourceTool::Redux];

    pub fn name(self) -> &'static str {
        match self {
            SourceTool::BamMetrics => "BAM_METRICS",
            SourceTool::Cobalt => "COBALT",
            SourceTool::Esvee => "ESVEE",
            SourceTool::Purple => "PURPLE",
            SourceTool::Redux => "REDUX",
        }
    }
    /// SourceTool.valueOf
    pub fn value_of(s: &str) -> anyhow::Result<SourceTool> {
        SourceTool::ALL
            .into_iter()
            .find(|t| t.name() == s)
            .ok_or_else(|| anyhow::anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.qsee.feature.SourceTool.{s}"))
    }
}

/// NumberFormat
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NumberFormat {
    None,
    Number,
    Percent,
    Log10,
}

impl NumberFormat {
    /// displayString: "" for NONE, else the name
    pub fn display_string(self) -> &'static str {
        match self {
            NumberFormat::None => "",
            NumberFormat::Number => "NUMBER",
            NumberFormat::Percent => "PERCENT",
            NumberFormat::Log10 => "LOG10",
        }
    }
}

/// FeatureKey: equal (and hashed) by name and type; the source tool takes part in the order only
#[derive(Clone, Debug)]
pub struct FeatureKey {
    pub name: String,
    pub feature_type: Option<FeatureType>,
    pub source_tool: Option<SourceTool>,
}

impl PartialEq for FeatureKey {
    fn eq(&self, other: &FeatureKey) -> bool {
        self.name == other.name && self.feature_type == other.feature_type
    }
}
impl Eq for FeatureKey {}
impl std::hash::Hash for FeatureKey {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        self.feature_type.hash(state);
    }
}

impl FeatureKey {
    pub fn new(name: impl Into<String>, feature_type: FeatureType, source_tool: SourceTool) -> FeatureKey {
        FeatureKey { name: name.into(), feature_type: Some(feature_type), source_tool: Some(source_tool) }
    }

    /// toString, as the log lines show it
    pub fn display(&self) -> String {
        match self.feature_type {
            None => self.name.clone(),
            Some(t) => format!("name({}) type({})", self.name, t.name()),
        }
    }

    /// FeatureKey.compareTo: by type, then source tool (nulls last), then by the name in the type's own
    /// order. Java throws for a name its type cannot parse; the keys made here always parse, and such a
    /// pair compares as equal.
    pub fn compare(&self, other: &FeatureKey) -> Ordering {
        let d = nulls_last(self.feature_type, other.feature_type);
        if d != Ordering::Equal {
            return d;
        }
        let d = nulls_last(self.source_tool, other.source_tool);
        if d != Ordering::Equal {
            return d;
        }
        match self.feature_type {
            Some(FeatureType::SummaryTable) => {
                let (a, b) = (SummaryTableFeature::value_of(&self.name), SummaryTableFeature::value_of(&other.name));
                match (a, b) {
                    (Some(a), Some(b)) => a.cmp(&b),
                    _ => Ordering::Equal,
                }
            }
            Some(FeatureType::DiscordantFragFreq) => {
                let group = |n: &str| java_text::split(n, ':').get(1).and_then(|g| DiscordantFragGroup::value_of(g));
                match (group(&self.name), group(&other.name)) {
                    (Some(a), Some(b)) => a.cmp(&b),
                    _ => Ordering::Equal,
                }
            }
            Some(FeatureType::BqrPerOrigQual) | Some(FeatureType::BqrPerSnv96Context) => {
                let bin = |n: &str| -> Option<BaseQualBin> {
                    let pairs = parse_multi_field_as_pairs(n).ok()?;
                    let last = pairs.last()?;
                    BaseQualBin::value_of(java_text::split(&last.1, ' ').first().copied().unwrap_or(""))
                };
                match (bin(&self.name), bin(&other.name)) {
                    (Some(a), Some(b)) => a.cmp(&b),
                    _ => Ordering::Equal,
                }
            }
            Some(FeatureType::MsIndelErrorRates) | Some(FeatureType::MsIndelErrorBias) => compare_ms_indels(&self.name, &other.name),
            _ => {
                if self.name.contains(FIELD_KEY_VALUE_SEPARATOR) && other.name.contains(FIELD_KEY_VALUE_SEPARATOR) {
                    compare_with_numerical_values(&self.name, &other.name)
                } else {
                    java_text::compare(&self.name, &other.name)
                }
            }
        }
    }
}

fn nulls_last<T: Ord>(a: Option<T>, b: Option<T>) -> Ordering {
    match (a, b) {
        (Some(a), Some(b)) => a.cmp(&b),
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Greater,
        (Some(_), None) => Ordering::Less,
    }
}

fn compare_ms_indels(a: &str, b: &str) -> Ordering {
    let (Ok(pa), Ok(pb)) = (parse_multi_field_as_pairs(a), parse_multi_field_as_pairs(b)) else {
        return Ordering::Equal;
    };
    if pa.len() < 3 || pb.len() < 3 {
        return Ordering::Equal;
    }
    let (Ok(ca), Ok(cb)) = (ConsensusType::value_of(&pa[0].1), ConsensusType::value_of(&pb[0].1)) else {
        return Ordering::Equal;
    };
    let d = ca.cmp(&cb);
    if d != Ordering::Equal {
        return d;
    }
    let (Some(ra), Some(rb)) = (RepeatType::from_display_name(&pa[1].1), RepeatType::from_display_name(&pb[1].1)) else {
        return Ordering::Equal;
    };
    let d = ra.cmp(&rb);
    if d != Ordering::Equal {
        return d;
    }
    let (Some(na), Some(nb)) = (java_text::parse_double(&pa[2].1), java_text::parse_double(&pb[2].1)) else {
        return Ordering::Equal;
    };
    hmf_common::java_math::double_compare(na, nb)
}

/// compareNameWithNumericalValues: the names split at ';' and ':', compared part by part, as numbers
/// where both parts parse (Double.compare), else as strings
fn compare_with_numerical_values(a: &str, b: &str) -> Ordering {
    let pa = parse_multi_field_as_array(a);
    let pb = parse_multi_field_as_array(b);
    for i in 0..pa.len().min(pb.len()) {
        match (java_text::parse_double(&pa[i]), java_text::parse_double(&pb[i])) {
            (Some(x), Some(y)) => {
                let d = hmf_common::java_math::double_compare(x, y);
                if d != Ordering::Equal {
                    return d;
                }
            }
            _ => {
                let d = java_text::compare(&pa[i], &pb[i]);
                if d != Ordering::Equal {
                    return d;
                }
            }
        }
    }
    pa.len().cmp(&pb.len())
}

/// FeatureMetadata
#[derive(Clone, Debug)]
pub struct FeatureMetadata {
    pub feature_group: String,
    pub display_name: String,
    pub number_format: NumberFormat,
    pub qc_status: QcStatus,
}

impl Default for FeatureMetadata {
    /// FeatureMetadata.createEmpty
    fn default() -> FeatureMetadata {
        FeatureMetadata { feature_group: String::new(), display_name: String::new(), number_format: NumberFormat::None, qc_status: QcStatus::empty() }
    }
}

impl FeatureMetadata {
    /// displayString: the vis data file's PlotMetadata column (an error where Java's threshold display
    /// fails)
    pub fn display_string(&self) -> anyhow::Result<String> {
        let q = &self.qc_status;
        let threshold = crate::status::try_form_display_string(q.status_type, q.operator, q.threshold, self.number_format)?;
        Ok(crate::common::form_multi_field(&[
            ("FeatureGroup", &self.feature_group),
            ("DisplayName", &self.display_name),
            ("NumberFormat", self.number_format.display_string()),
            ("QcStatus", self.qc_status.status_type.display_string()),
            ("QcThreshold", &threshold),
        ]))
    }
}

/// Feature
#[derive(Clone, Debug)]
pub struct Feature {
    pub key: FeatureKey,
    pub value: f64,
    pub metadata: FeatureMetadata,
}

impl Feature {
    pub fn new(key: FeatureKey, value: f64) -> Feature {
        Feature { key, value, metadata: FeatureMetadata::default() }
    }
    pub fn with_metadata(key: FeatureKey, value: f64, metadata: FeatureMetadata) -> Feature {
        Feature { key, value, metadata }
    }
}

/// SampleFeatures
#[derive(Clone, Debug)]
pub struct SampleFeatures {
    pub sample_id: String,
    pub sample_type: SampleType,
    pub features: Vec<Feature>,
}

// ---- the category enums the key order and the features use ----

/// SummaryTableGroup
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SummaryTableGroup {
    Tmb,
    Contamination,
    CopyNumber,
    Mapping,
}

impl SummaryTableGroup {
    pub fn display_name(self) -> &'static str {
        match self {
            SummaryTableGroup::Tmb => "Mutational burden",
            SummaryTableGroup::Contamination => "Contamination",
            SummaryTableGroup::CopyNumber => "Copy number",
            SummaryTableGroup::Mapping => "Mapping",
        }
    }
}

/// SummaryTableFeature, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SummaryTableFeature {
    MeanCoverage,
    MappedProportion,
    CoverageAbove10,
    CoverageAbove20,
    CoverageAbove30,
    CoverageAbove60,
    CoverageAbove100,
    CoverageAbove250,
    LowBaseQual,
    LowMapQual,
    DuplicateReads,
    DualStrandReads,
    /// upstream bdceb2302a: the proportion of assessed exons of reportable genes with a LoD VAF WARN
    LowCoverageExons,
    Purity,
    Ploidy,
    LohPercent,
    DeletedGenes,
    UnsupportedCnSegments,
    Tinc,
    Contamination,
    TmbSmallVariants,
    TmbMsIndels,
    TmbStructuralVariants,
}

impl SummaryTableFeature {
    pub const ALL: [SummaryTableFeature; 23] = [
        SummaryTableFeature::MeanCoverage,
        SummaryTableFeature::MappedProportion,
        SummaryTableFeature::CoverageAbove10,
        SummaryTableFeature::CoverageAbove20,
        SummaryTableFeature::CoverageAbove30,
        SummaryTableFeature::CoverageAbove60,
        SummaryTableFeature::CoverageAbove100,
        SummaryTableFeature::CoverageAbove250,
        SummaryTableFeature::LowBaseQual,
        SummaryTableFeature::LowMapQual,
        SummaryTableFeature::DuplicateReads,
        SummaryTableFeature::DualStrandReads,
        SummaryTableFeature::LowCoverageExons,
        SummaryTableFeature::Purity,
        SummaryTableFeature::Ploidy,
        SummaryTableFeature::LohPercent,
        SummaryTableFeature::DeletedGenes,
        SummaryTableFeature::UnsupportedCnSegments,
        SummaryTableFeature::Tinc,
        SummaryTableFeature::Contamination,
        SummaryTableFeature::TmbSmallVariants,
        SummaryTableFeature::TmbMsIndels,
        SummaryTableFeature::TmbStructuralVariants,
    ];

    pub fn name(self) -> &'static str {
        use SummaryTableFeature::*;
        match self {
            MeanCoverage => "MEAN_COVERAGE",
            MappedProportion => "MAPPED_PROPORTION",
            CoverageAbove10 => "COVERAGE_ABOVE_10",
            CoverageAbove20 => "COVERAGE_ABOVE_20",
            CoverageAbove30 => "COVERAGE_ABOVE_30",
            CoverageAbove60 => "COVERAGE_ABOVE_60",
            CoverageAbove100 => "COVERAGE_ABOVE_100",
            CoverageAbove250 => "COVERAGE_ABOVE_250",
            LowBaseQual => "LOW_BASE_QUAL",
            LowMapQual => "LOW_MAP_QUAL",
            DuplicateReads => "DUPLICATE_READS",
            DualStrandReads => "DUAL_STRAND_READS",
            LowCoverageExons => "LOW_COVERAGE_EXONS",
            Purity => "PURITY",
            Ploidy => "PLOIDY",
            LohPercent => "LOH_PERCENT",
            DeletedGenes => "DELETED_GENES",
            UnsupportedCnSegments => "UNSUPPORTED_CN_SEGMENTS",
            Tinc => "TINC",
            Contamination => "CONTAMINATION",
            TmbSmallVariants => "TMB_SMALL_VARIANTS",
            TmbMsIndels => "TMB_MS_INDELS",
            TmbStructuralVariants => "TMB_STRUCTURAL_VARIANTS",
        }
    }

    pub fn value_of(s: &str) -> Option<SummaryTableFeature> {
        SummaryTableFeature::ALL.into_iter().find(|f| f.name() == s)
    }

    pub fn source_tool(self) -> SourceTool {
        if self <= SummaryTableFeature::LowCoverageExons { SourceTool::BamMetrics } else { SourceTool::Purple }
    }

    pub fn group(self) -> SummaryTableGroup {
        use SummaryTableFeature::*;
        match self {
            Purity | Ploidy | LohPercent | DeletedGenes | UnsupportedCnSegments => SummaryTableGroup::CopyNumber,
            Tinc | Contamination => SummaryTableGroup::Contamination,
            TmbSmallVariants | TmbMsIndels | TmbStructuralVariants => SummaryTableGroup::Tmb,
            _ => SummaryTableGroup::Mapping,
        }
    }

    pub fn display_name(self) -> &'static str {
        use SummaryTableFeature::*;
        match self {
            MeanCoverage => "Mean coverage",
            MappedProportion => "Mapped proportion",
            CoverageAbove10 => "Coverage >10x",
            CoverageAbove20 => "Coverage >20x",
            CoverageAbove30 => "Coverage >30x",
            CoverageAbove60 => "Coverage >60x",
            CoverageAbove100 => "Coverage >100x",
            CoverageAbove250 => "Coverage >250x",
            LowBaseQual => "Low qual. bases",
            LowMapQual => "Low map qual. reads",
            DuplicateReads => "Duplicate reads",
            DualStrandReads => "Dual-strand reads",
            LowCoverageExons => "Low coverage exons",
            Purity => "Purity",
            Ploidy => "Ploidy",
            LohPercent => "LOH",
            DeletedGenes => "Deleted genes",
            UnsupportedCnSegments => "Unsupported CN segs",
            Tinc => "Tumor in normal",
            Contamination => "Other DNA",
            TmbSmallVariants => "SNVs/indels per MB",
            TmbMsIndels => "MS indels per MB",
            TmbStructuralVariants => "SVs per MB",
        }
    }

    pub fn number_format(self) -> NumberFormat {
        use SummaryTableFeature::*;
        match self {
            MeanCoverage | Ploidy => NumberFormat::Number,
            DeletedGenes | UnsupportedCnSegments | TmbSmallVariants | TmbMsIndels | TmbStructuralVariants => NumberFormat::Log10,
            _ => NumberFormat::Percent,
        }
    }

    /// SummaryTableFeature.putFeature: the feature with its metadata and QC status
    pub fn feature(self, value: f64, qc_status: QcStatus) -> Feature {
        let key = FeatureKey::new(self.name(), FeatureType::SummaryTable, self.source_tool());
        let metadata = FeatureMetadata {
            feature_group: self.group().display_name().to_string(),
            display_name: self.display_name().to_string(),
            number_format: self.number_format(),
            qc_status,
        };
        Feature::with_metadata(key, value, metadata)
    }
}

/// DiscordantFragGroup, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum DiscordantFragGroup {
    DelDupShort,
    DelDupMedium,
    DelDupLong,
    InvShort,
    InvMedium,
    InvLong,
    Translocation,
}

impl DiscordantFragGroup {
    pub const ALL: [DiscordantFragGroup; 7] = [
        DiscordantFragGroup::DelDupShort,
        DiscordantFragGroup::DelDupMedium,
        DiscordantFragGroup::DelDupLong,
        DiscordantFragGroup::InvShort,
        DiscordantFragGroup::InvMedium,
        DiscordantFragGroup::InvLong,
        DiscordantFragGroup::Translocation,
    ];

    pub fn name(self) -> &'static str {
        match self {
            DiscordantFragGroup::DelDupShort => "DEL_DUP_SHORT",
            DiscordantFragGroup::DelDupMedium => "DEL_DUP_MEDIUM",
            DiscordantFragGroup::DelDupLong => "DEL_DUP_LONG",
            DiscordantFragGroup::InvShort => "INV_SHORT",
            DiscordantFragGroup::InvMedium => "INV_MEDIUM",
            DiscordantFragGroup::InvLong => "INV_LONG",
            DiscordantFragGroup::Translocation => "TRANSLOCATION",
        }
    }
    pub fn value_of(s: &str) -> Option<DiscordantFragGroup> {
        DiscordantFragGroup::ALL.into_iter().find(|g| g.name() == s)
    }
    pub fn display_name(self) -> &'static str {
        match self {
            DiscordantFragGroup::DelDupShort => "Del/Dup 1K-5K",
            DiscordantFragGroup::DelDupMedium => "Del/Dup 5K-100K",
            DiscordantFragGroup::DelDupLong => "Del/Dup >100K",
            DiscordantFragGroup::InvShort => "Inv <=1K",
            DiscordantFragGroup::InvMedium => "Inv 1K-100K",
            DiscordantFragGroup::InvLong => "Inv >100K",
            DiscordantFragGroup::Translocation => "Translocation",
        }
    }
    /// DiscordantFragGroup.fromType, by the DiscordantFragType ordinal (Translocation, InvLt1K,
    /// Inv1To5K, Inv5To100K, InvGt100K, Del1To5K, Del5To100K, DelGt100K, Dup1To5K, Dup5To100K, DupGt100K)
    pub fn from_type_index(i: usize) -> DiscordantFragGroup {
        match i {
            0 => DiscordantFragGroup::Translocation,
            1 => DiscordantFragGroup::InvShort,
            2 | 3 => DiscordantFragGroup::InvMedium,
            4 => DiscordantFragGroup::InvLong,
            5 | 8 => DiscordantFragGroup::DelDupShort,
            6 | 9 => DiscordantFragGroup::DelDupMedium,
            _ => DiscordantFragGroup::DelDupLong,
        }
    }
}

/// BaseQualBin, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BaseQualBin {
    Low,
    Medium,
    High,
}

impl BaseQualBin {
    pub fn name(self) -> &'static str {
        match self {
            BaseQualBin::Low => "LOW",
            BaseQualBin::Medium => "MEDIUM",
            BaseQualBin::High => "HIGH",
        }
    }
    pub fn value_of(s: &str) -> Option<BaseQualBin> {
        match s {
            "LOW" => Some(BaseQualBin::Low),
            "MEDIUM" => Some(BaseQualBin::Medium),
            "HIGH" => Some(BaseQualBin::High),
            _ => None,
        }
    }
}

/// RepeatType, in enum order
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum RepeatType {
    RepeatAT,
    RepeatCG,
    Repeat2Bp,
    Repeat3BpOrLonger,
}

impl RepeatType {
    pub fn display_name(self) -> &'static str {
        match self {
            RepeatType::RepeatAT => "A/T repeat",
            RepeatType::RepeatCG => "C/G repeat",
            RepeatType::Repeat2Bp => "2bp repeat",
            RepeatType::Repeat3BpOrLonger => ">=3bp repeat",
        }
    }
    pub fn from_display_name(s: &str) -> Option<RepeatType> {
        [RepeatType::RepeatAT, RepeatType::RepeatCG, RepeatType::Repeat2Bp, RepeatType::Repeat3BpOrLonger].into_iter().find(|r| r.display_name() == s)
    }
    /// RepeatType.fromRepeatUnit: "A/T", "C/G", two word characters and a '/' ("^\w{2}/.*"), or
    /// "<digits>bp repeat"; anything else is an error (IllegalArgumentException)
    pub fn from_repeat_unit(unit: &str) -> anyhow::Result<RepeatType> {
        if unit == "A/T" {
            return Ok(RepeatType::RepeatAT);
        }
        if unit == "C/G" {
            return Ok(RepeatType::RepeatCG);
        }
        let word = |c: char| c.is_ascii_alphanumeric() || c == '_';
        let chars: Vec<char> = unit.chars().collect();
        // String.matches("^\\w{2}/.*"): '.' does not match a line terminator
        if chars.len() >= 3 && word(chars[0]) && word(chars[1]) && chars[2] == '/' && !chars[3..].iter().any(|&c| is_line_terminator(c)) {
            return Ok(RepeatType::Repeat2Bp);
        }
        if let Some(digits) = unit.strip_suffix("bp repeat") {
            if !digits.is_empty() && digits.chars().all(|c| c.is_ascii_digit()) {
                return Ok(RepeatType::Repeat3BpOrLonger);
            }
        }
        anyhow::bail!("java.lang.IllegalArgumentException: Unexpected repeat unit: {unit}")
    }
}

/// a line terminator as java.util.regex's '.' excludes it
fn is_line_terminator(c: char) -> bool {
    matches!(c, '\n' | '\r' | '\u{85}' | '\u{2028}' | '\u{2029}')
}
