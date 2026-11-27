//! The data types shared by the steps: sites, the evidence counted at each site, tumor BAFs and the
//! LOH state of a tumor-only site. Names follow the Java classes (hmf-common `amber` package and
//! `com.hartwig.hmftools.amber`).

use std::cmp::Ordering;

/// AmberBase. The Java enum's order (G, A, T, C, N) is its ordinal order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AmberBase {
    G,
    A,
    T,
    C,
    N,
}

impl AmberBase {
    /// AmberBase.valueOf for a one-letter name; None where Java throws IllegalArgumentException.
    pub fn from_char(c: u8) -> Option<AmberBase> {
        match c {
            b'G' => Some(AmberBase::G),
            b'A' => Some(AmberBase::A),
            b'T' => Some(AmberBase::T),
            b'C' => Some(AmberBase::C),
            b'N' => Some(AmberBase::N),
            _ => None,
        }
    }
    pub fn from_str(s: &str) -> Option<AmberBase> {
        if s.len() == 1 { AmberBase::from_char(s.as_bytes()[0]) } else { None }
    }
    pub fn as_str(self) -> &'static str {
        match self {
            AmberBase::G => "G",
            AmberBase::A => "A",
            AmberBase::T => "T",
            AmberBase::C => "C",
            AmberBase::N => "N",
        }
    }
    pub fn ordinal(self) -> i32 {
        self as i32
    }
    pub fn complement(self) -> AmberBase {
        match self {
            AmberBase::A => AmberBase::T,
            AmberBase::T => AmberBase::A,
            AmberBase::C => AmberBase::G,
            AmberBase::G => AmberBase::C,
            AmberBase::N => AmberBase::N,
        }
    }
}

/// ContigComparator, as used by GenomePosition.compare: the chromosome ranks only (1-22, X, Y, MT, then
/// every other contig with one rank), so two contigs of one rank compare equal and a stable sort keeps
/// their order.
pub fn compare_chromosomes(a: &str, b: &str) -> Ordering {
    hmf_common::regions::chromosome_rank(a).cmp(&hmf_common::regions::chromosome_rank(b))
}

/// GenomePosition.compare: chromosome, then position.
pub fn compare_positions(chr_a: &str, pos_a: i32, chr_b: &str, pos_b: i32) -> Ordering {
    compare_chromosomes(chr_a, chr_b).then(pos_a.cmp(&pos_b))
}

/// AmberSite: one line of the loci file.
#[derive(Clone, Debug)]
pub struct AmberSite {
    pub chromosome: String,
    pub position: i32,
    pub ref_: String,
    pub alt: String,
    pub snp_check: bool,
    /// GnomadFrequency, 0.5 when the file has none
    pub gnomad_frequency: f64,
}

impl AmberSite {
    pub fn new(chromosome: &str, position: i32, ref_: &str, alt: &str, snp_check: bool) -> AmberSite {
        AmberSite::with_frequency(chromosome, position, ref_, alt, snp_check, 0.5)
    }
    pub fn with_frequency(chromosome: &str, position: i32, ref_: &str, alt: &str, snp_check: bool, gnomad_frequency: f64) -> AmberSite {
        AmberSite { chromosome: chromosome.to_string(), position, ref_: ref_.to_string(), alt: alt.to_string(), snp_check, gnomad_frequency }
    }
    /// AmberSite.matches (equals ignores the SNP-check flag and the frequency)
    pub fn matches(&self, chromosome: &str, position: i32, ref_: &str, alt: &str) -> bool {
        self.chromosome == chromosome && self.position == position && self.ref_ == ref_ && self.alt == alt
    }
    /// AmberSite.hashCode (djb2-style over chromosome, position, ref, alt and the SNP-check flag)
    pub fn java_hash(&self) -> i32 {
        use hmf_common::java_hash::java_string_hash;
        let mut h: i32 = 5381;
        h = h.wrapping_add((h << 5).wrapping_add(java_string_hash(&self.chromosome)));
        h = h.wrapping_add((h << 5).wrapping_add(self.position));
        h = h.wrapping_add((h << 5).wrapping_add(java_string_hash(&self.ref_)));
        h = h.wrapping_add((h << 5).wrapping_add(java_string_hash(&self.alt)));
        h = h.wrapping_add((h << 5).wrapping_add(if self.snp_check { 1231 } else { 1237 }));
        h
    }
}

/// BafLohData: the tumor-only LOH state of a site (written by LohIdentifier, read by the BAF filters).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BafLohData {
    pub valid_baf_high: bool,
    pub valid_baf_medium: bool,
    pub valid_baf_low: bool,
    pub valid_rate_high: f64,
    pub valid_rate_medium: f64,
    pub valid_rate_low: f64,
    pub loh_candidate: bool,
    pub loh_recovered: bool,
}

/// PositionEvidence: the reads counted at one site of one sample.
#[derive(Clone, Debug)]
pub struct PositionEvidence {
    pub chromosome: String,
    pub position: i32,
    pub ref_: AmberBase,
    pub alt: AmberBase,
    pub read_depth: i32,
    pub indel_count: i32,
    pub ref_support: i32,
    pub alt_support: i32,
    pub base_qual_filtered: i32,
    pub map_qual_filtered: i32,
    pub seq_tech_filtered: i32,
    pub loh_data: Option<BafLohData>,
}

impl PositionEvidence {
    /// Java checks Ref != Alt (Preconditions.checkArgument); the caller makes sure of it.
    pub fn new(chromosome: &str, position: i32, ref_: AmberBase, alt: AmberBase) -> PositionEvidence {
        PositionEvidence {
            chromosome: chromosome.to_string(),
            position,
            ref_,
            alt,
            read_depth: 0,
            indel_count: 0,
            ref_support: 0,
            alt_support: 0,
            base_qual_filtered: 0,
            map_qual_filtered: 0,
            seq_tech_filtered: 0,
            loh_data: None,
        }
    }
    /// PositionEvidence.copy: same site, no counts
    pub fn copy_site(&self) -> PositionEvidence {
        PositionEvidence::new(&self.chromosome, self.position, self.ref_, self.alt)
    }
    pub fn is_valid(&self) -> bool {
        self.indel_count == 0
    }
    /// vaf(): NaN when there is no depth; note the denominator is ref + alt, not the depth
    pub fn vaf(&self) -> f64 {
        if self.read_depth == 0 {
            return f64::NAN;
        }
        self.alt_support as f64 / (self.alt_support + self.ref_support) as f64
    }
    pub fn symmetric_vaf(&self) -> f64 {
        let vaf_depth = self.alt_support + self.ref_support;
        if vaf_depth == 0 {
            return f64::NAN;
        }
        let mut k = self.alt_support;
        if k > vaf_depth / 2 {
            k = self.ref_support;
        }
        k as f64 / vaf_depth as f64
    }
    pub fn loh_recovered(&self) -> bool {
        self.loh_data.as_ref().is_some_and(|d| d.loh_recovered)
    }
    pub fn compare(&self, other: &PositionEvidence) -> Ordering {
        compare_positions(&self.chromosome, self.position, &other.chromosome, other.position)
    }
}

/// TumorBAF: tumor evidence at a germline heterozygous site, with the normal's counts.
#[derive(Clone, Debug)]
pub struct TumorBaf {
    pub normal_read_depth: i32,
    pub normal_ref_support: i32,
    pub normal_alt_support: i32,
    pub tumor: PositionEvidence,
}

impl TumorBaf {
    pub fn from_normal(normal: &PositionEvidence) -> TumorBaf {
        TumorBaf {
            normal_read_depth: normal.read_depth,
            normal_ref_support: normal.ref_support,
            normal_alt_support: normal.alt_support,
            tumor: normal.copy_site(),
        }
    }
    pub fn ref_frequency(&self) -> f64 {
        self.tumor.ref_support as f64 / self.tumor.read_depth as f64
    }
    pub fn alt_frequency(&self) -> f64 {
        self.tumor.alt_support as f64 / self.tumor.read_depth as f64
    }
}

/// AmberBAF: one line of the BAF file.
#[derive(Clone, Debug)]
pub struct AmberBaf {
    pub chromosome: String,
    pub position: i32,
    pub tumor_baf: f64,
    pub tumor_depth: i32,
    pub normal_baf: f64,
    pub normal_depth: i32,
    pub loh_recovered: bool,
}

impl AmberBaf {
    pub fn tumor_modified_baf(&self) -> f64 {
        0.5 + (self.tumor_baf - 0.5).abs()
    }
    pub fn normal_modified_baf(&self) -> f64 {
        0.5 + (self.normal_baf - 0.5).abs()
    }
}
