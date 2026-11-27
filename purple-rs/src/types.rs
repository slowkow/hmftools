//! Small types of hmf-common used throughout: the purple and driver enums, Doubles, chromosome
//! comparison as GenomePosition and GenomeRegion compare.

use std::cmp::Ordering;

pub use hmf_common::segmentation::HumanChromosome as Chromosome;

/// An enum with Java's name() and valueOf()
macro_rules! java_enum {
    ($name:ident { $($variant:ident = $text:literal),* $(,)? }) => {
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
        #[allow(dead_code)]
        pub enum $name { $($variant),* }
        #[allow(dead_code)]
        impl $name {
            pub const ALL: &'static [$name] = &[$($name::$variant),*];
            pub fn name(self) -> &'static str {
                match self { $($name::$variant => $text),* }
            }
            pub fn from_name(s: &str) -> Option<$name> {
                match s { $($text => Some($name::$variant),)* _ => None }
            }
            pub fn ordinal(self) -> usize {
                self as usize
            }
        }
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str(self.name())
            }
        }
    };
}
pub(crate) use java_enum;

java_enum!(GermlineStatus {
    HomDeletion = "HOM_DELETION",
    HetDeletion = "HET_DELETION",
    Amplification = "AMPLIFICATION",
    Noise = "NOISE",
    Diploid = "DIPLOID",
    LikelyDiploid = "LIKELY_DIPLOID",
    Excluded = "EXCLUDED",
    Centrometic = "CENTROMETIC",
    Unknown = "UNKNOWN",
});

java_enum!(Gender { Male = "MALE", Female = "FEMALE" });

java_enum!(SegmentSupport {
    Bnd = "BND",
    Inv = "INV",
    Del = "DEL",
    Sgl = "SGL",
    Ins = "INS",
    Inf = "INF",
    Dup = "DUP",
    Multiple = "MULTIPLE",
    Centromere = "CENTROMERE",
    Telomere = "TELOMERE",
    None = "NONE",
    Excl = "EXCL",
    Unknown = "UNKNOWN",
});

impl SegmentSupport {
    pub fn is_sv(self) -> bool {
        (self as usize) <= SegmentSupport::Multiple as usize
    }
    pub fn from_variant(t: SvType) -> SegmentSupport {
        match t {
            SvType::Inf => SegmentSupport::Inf,
            SvType::Bnd => SegmentSupport::Bnd,
            SvType::Del => SegmentSupport::Del,
            SvType::Ins => SegmentSupport::Ins,
            SvType::Dup => SegmentSupport::Dup,
            SvType::Inv => SegmentSupport::Inv,
            SvType::Sgl => SegmentSupport::Sgl,
        }
    }
}

java_enum!(SvType { Bnd = "BND", Del = "DEL", Dup = "DUP", Inf = "INF", Ins = "INS", Inv = "INV", Sgl = "SGL" });

java_enum!(CopyNumberMethod {
    Unknown = "UNKNOWN",
    BafWeighted = "BAF_WEIGHTED",
    StructuralVariant = "STRUCTURAL_VARIANT",
    LongArm = "LONG_ARM",
    HomozygousDeletion = "HOMOZYGOUS_DELETION",
    GermlineHomDeletion = "GERMLINE_HOM_DELETION",
    GermlineHet2HomDeletion = "GERMLINE_HET2HOM_DELETION",
    NonDiploid = "NON_DIPLOID",
});

java_enum!(FittedPurityMethod { Normal = "NORMAL", HighlyDiploid = "HIGHLY_DIPLOID", Somatic = "SOMATIC", NoTumor = "NO_TUMOR" });

java_enum!(RunMode { TumorGermline = "TUMOR_GERMLINE", Tumor = "TUMOR", Germline = "GERMLINE" });

java_enum!(ReportedStatus { None = "NONE", NotReported = "NOT_REPORTED", Reported = "REPORTED" });

java_enum!(GermlineDetectionMethod { Segment = "SEGMENT", ExonCoverage = "EXON_COVERAGE" });

java_enum!(MicrosatelliteStatus { Msi = "MSI", Mss = "MSS", Unknown = "UNKNOWN" });

impl MicrosatelliteStatus {
    pub fn from_indels_per_mb(v: f64) -> MicrosatelliteStatus {
        if greater_or_equal(v, 4.0) { MicrosatelliteStatus::Msi } else { MicrosatelliteStatus::Mss }
    }
}

java_enum!(TumorMutationalStatus { High = "HIGH", Low = "LOW", Unknown = "UNKNOWN" });

impl TumorMutationalStatus {
    pub fn from_burden_per_mb(v: f64) -> TumorMutationalStatus {
        if greater_or_equal(v, 10.0) { TumorMutationalStatus::High } else { TumorMutationalStatus::Low }
    }
    pub fn from_load(load: f64) -> TumorMutationalStatus {
        if load >= 140.0 { TumorMutationalStatus::High } else { TumorMutationalStatus::Low }
    }
}

java_enum!(GermlineAberration {
    None = "NONE",
    MosaicX = "MOSAIC_X",
    Klinefelter = "KLINEFELTER",
    Xyy = "XYY",
    TrisomyX = "TRISOMY_X",
    Trisomy13 = "TRISOMY_13",
    Trisomy15 = "TRISOMY_15",
    Trisomy18 = "TRISOMY_18",
    Trisomy21 = "TRISOMY_21",
    Tertrasomy9 = "TERTRASOMY_9",
});

java_enum!(DriverType {
    Amp = "AMP",
    PartialAmp = "PARTIAL_AMP",
    Del = "DEL",
    HetDel = "HET_DEL",
    Loh = "LOH",
    HomDupDisruption = "HOM_DUP_DISRUPTION",
    HomDelDisruption = "HOM_DEL_DISRUPTION",
    Disruption = "DISRUPTION",
    Mutation = "MUTATION",
    GermlineMutation = "GERMLINE_MUTATION",
    GermlineDeletion = "GERMLINE_DELETION",
    GermlineAmp = "GERMLINE_AMP",
    GermlineDisruption = "GERMLINE_DISRUPTION",
    Unknown = "UNKNOWN",
});

java_enum!(DriverCategory { Onco = "ONCO", Tsg = "TSG" });

java_enum!(LikelihoodMethod {
    Biallelic = "BIALLELIC",
    Hotspot = "HOTSPOT",
    Inframe = "INFRAME",
    Dnds = "DNDS",
    SpliceRegion = "SPLICE_REGION",
    Amp = "AMP",
    Del = "DEL",
    Disruption = "DISRUPTION",
    Germline = "GERMLINE",
    None = "NONE",
});

java_enum!(DriverInterpretation { High = "HIGH", Medium = "MEDIUM", Low = "LOW" });

impl DriverInterpretation {
    pub fn interpret(likelihood: f64) -> DriverInterpretation {
        if likelihood >= 0.8 {
            DriverInterpretation::High
        } else if likelihood >= 0.2 {
            DriverInterpretation::Medium
        } else {
            DriverInterpretation::Low
        }
    }
}

java_enum!(DriverImpact {
    Missense = "MISSENSE",
    Nonsense = "NONSENSE",
    Splice = "SPLICE",
    Inframe = "INFRAME",
    Frameshift = "FRAMESHIFT",
    Unknown = "UNKNOWN",
});

java_enum!(CodingEffect {
    NonsenseOrFrameshift = "NONSENSE_OR_FRAMESHIFT",
    Splice = "SPLICE",
    Missense = "MISSENSE",
    Synonymous = "SYNONYMOUS",
    None = "NONE",
    Undefined = "UNDEFINED",
});

impl CodingEffect {
    pub fn has_protein_impact(self) -> bool {
        matches!(self, CodingEffect::NonsenseOrFrameshift | CodingEffect::Splice | CodingEffect::Missense)
    }
}

java_enum!(VariantType { Mnp = "MNP", Snp = "SNP", Indel = "INDEL", Undefined = "UNDEFINED" });

java_enum!(VariantTier { Unknown = "UNKNOWN", LowConfidence = "LOW_CONFIDENCE", HighConfidence = "HIGH_CONFIDENCE", Panel = "PANEL", Hotspot = "HOTSPOT" });

java_enum!(Pathogenicity {
    Pathogenic = "PATHOGENIC",
    LikelyPathogenic = "LIKELY_PATHOGENIC",
    Benign = "BENIGN",
    LikelyBenign = "LIKELY_BENIGN",
    BenignBlacklist = "BENIGN_BLACKLIST",
    Conflicting = "CONFLICTING",
    Unknown = "UNKNOWN",
});

impl Pathogenicity {
    pub fn is_pathogenic(self) -> bool {
        matches!(self, Pathogenicity::Pathogenic | Pathogenicity::LikelyPathogenic)
    }
}

java_enum!(GenotypeStatus { HomRef = "HOM_REF", Het = "HET", HomAlt = "HOM_ALT", Unknown = "UNKNOWN" });

java_enum!(SomaticLikelihood { Unknown = "UNKNOWN", High = "HIGH", Medium = "MEDIUM", Low = "LOW" });

java_enum!(GermlineReporting { Any = "ANY", VariantNotLost = "VARIANT_NOT_LOST", WildtypeLost = "WILDTYPE_LOST", None = "NONE" });

// ---------------------------------------------------------------------------------------------
// Doubles

pub const EPSILON: f64 = 1e-10;

pub fn equal(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}
pub fn is_zero(v: f64) -> bool {
    equal(v, 0.0)
}
pub fn less_than(v: f64, r: f64) -> bool {
    v - r < -EPSILON
}
pub fn less_or_equal(v: f64, r: f64) -> bool {
    v - r < EPSILON
}
pub fn greater_than(v: f64, r: f64) -> bool {
    v - r > EPSILON
}
pub fn greater_or_equal(v: f64, r: f64) -> bool {
    v - r > -EPSILON
}
pub fn positive(v: f64) -> bool {
    greater_than(v, 0.0)
}
pub fn positive_or_zero(v: f64) -> bool {
    greater_or_equal(v, 0.0)
}
/// Doubles.absDistanceFromInteger
pub fn abs_distance_from_integer(v: f64) -> f64 {
    (v - hmf_common::java_math::java_round(v) as f64).abs()
}
/// Doubles.round(value, places): Math.round(value * 10^places) / 10^places
pub fn round(v: f64, places: i32) -> f64 {
    let log_value = hmf_common::java_math::java_round(jpow(10.0, places as f64)) as f64;
    hmf_common::java_math::java_round(v * log_value) as f64 / log_value
}
/// Doubles.median(List<Double>): sorted (Double.compareTo), the mean of the two middle values when even, 0 if empty
pub fn median(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| hmf_common::java_math::double_compare(*a, *b));
    let n = v.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 0 { (v[n / 2] + v[n / 2 - 1]) / 2.0 } else { v[n / 2] }
}
/// Doubles.median(double[]): NaN for an empty array
pub fn median_array(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    median(values)
}
/// Doubles.clamp
pub fn clamp(v: f64, min: f64, max: f64) -> f64 {
    jmax(jmin(v, max), min)
}

/// Math.round(double)
pub fn jround(v: f64) -> i64 {
    hmf_common::java_math::java_round(v)
}
/// Math.max(double, double)
pub fn jmax(a: f64, b: f64) -> f64 {
    hmf_common::java_math::java_max(a, b)
}
/// Math.min(double, double): NaN if either is NaN, -0.0 below 0.0
pub fn jmin(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && b.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a <= b { a } else { b }
}
pub fn jexp(x: f64) -> f64 {
    hmf_common::java_strictmath::exp(x)
}
pub fn jlog(x: f64) -> f64 {
    hmf_common::java_strictmath::log(x)
}
pub fn jlog10(x: f64) -> f64 {
    hmf_common::java_strictmath::log10(x)
}
pub fn jpow(x: f64, y: f64) -> f64 {
    hmf_common::java_strictmath::pow(x, y)
}
/// Double.compare
pub fn dcmp(a: f64, b: f64) -> Ordering {
    hmf_common::java_math::double_compare(a, b)
}

// ---------------------------------------------------------------------------------------------
// chromosomes

/// HumanChromosome.chromosomeRank
pub fn chromosome_rank(chr: &str) -> i32 {
    hmf_common::regions::chromosome_rank_java(chr)
}

/// ContigComparator: by rank only (equal ranks are equal)
pub fn contig_compare(a: &str, b: &str) -> Ordering {
    chromosome_rank(a).cmp(&chromosome_rank(b))
}

/// GenomePosition.compare
pub fn position_compare(c1: &str, p1: i32, c2: &str, p2: i32) -> Ordering {
    match contig_compare(c1, c2) {
        Ordering::Equal => p1.cmp(&p2),
        o => o,
    }
}

/// GenomeRegion.compareTo: same chromosome name by start only, else ContigComparator
pub fn region_compare(c1: &str, s1: i32, c2: &str, s2: i32) -> Ordering {
    if c1 == c2 { s1.cmp(&s2) } else { contig_compare(c1, c2) }
}

/// HumanChromosome.fromString, panicking as Java throws
pub fn chr(name: &str) -> Chromosome {
    Chromosome::parse(name).unwrap_or_else(|e| panic!("{e}"))
}

pub fn is_autosome(c: Chromosome) -> bool {
    c.ordinal() < 22
}

/// HumanChromosome.isDiploid(gender)
pub fn is_diploid_chromosome(c: Chromosome, gender: Gender) -> bool {
    is_autosome(c) || (gender != Gender::Male && c.short_name() == "X")
}

pub fn is_x(c: Chromosome) -> bool {
    c.ordinal() == 22
}
pub fn is_y(c: Chromosome) -> bool {
    c.ordinal() == 23
}

/// Genomic window helper (genome.region.Window)
#[derive(Clone, Copy)]
pub struct Window {
    pub size: i32,
}
impl Window {
    pub fn start(&self, position: i32) -> i32 {
        (position - 1) / self.size * self.size + 1
    }
    pub fn end(&self, position: i32) -> i32 {
        self.start(position) + self.size - 1
    }
}

/// BaseRegion.positionsWithin(innerStart, innerEnd, outerStart, outerEnd)
pub fn positions_within(inner_start: i32, inner_end: i32, outer_start: i32, outer_end: i32) -> bool {
    inner_start <= inner_end && inner_start >= outer_start && inner_end <= outer_end
}
/// BaseRegion.positionWithin
pub fn position_within(pos: i32, start: i32, end: i32) -> bool {
    pos >= start && pos <= end
}
/// BaseRegion.positionsOverlap
pub fn positions_overlap(s1: i32, e1: i32, s2: i32, e2: i32) -> bool {
    s1 <= e2 && e1 >= s2
}

/// A region on a chromosome (ChrBaseRegion)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChrRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}
impl ChrRegion {
    pub fn new(chromosome: &str, start: i32, end: i32) -> ChrRegion {
        ChrRegion { chromosome: chromosome.to_string(), start, end }
    }
    pub fn contains_position(&self, chromosome: &str, pos: i32) -> bool {
        self.chromosome == chromosome && pos >= self.start && pos <= self.end
    }
    pub fn overlaps(&self, chromosome: &str, start: i32, end: i32) -> bool {
        self.chromosome == chromosome && positions_overlap(self.start, self.end, start, end)
    }
}
