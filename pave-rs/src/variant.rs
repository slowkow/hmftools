//! VariantData: one variant with its derived positions, the SAGE details PAVE uses (local phase set,
//! microhomology, repeat), its transcript impacts by gene and its annotations.

use hmf_common::java_hash::JavaHashMap;
use hmf_common::regions::positions_within;
use hmf_common::vcf::VcfRecord;

use crate::impact::VariantTransImpact;

pub const NO_LOCAL_PHASE_SET: i32 = -1;

// SageVcfTags
pub const LOCAL_PHASE_SET: &str = "LPS";
pub const REPEAT_COUNT: &str = "REP_C";
pub const REPEAT_SEQUENCE: &str = "REP_S";
pub const READ_CONTEXT_REPEAT_COUNT: &str = "RC_REPC";
pub const READ_CONTEXT_REPEAT_SEQUENCE: &str = "RC_REPS";
pub const MICROHOMOLOGY: &str = "MH";

/// MultiPonStatus
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MultiPonStatus {
    Base,
    Artefact,
    Multi,
}

impl MultiPonStatus {
    /// MultiPonStatus.valueOf
    pub fn value_of(s: &str) -> anyhow::Result<MultiPonStatus> {
        match s {
            "BASE" => Ok(MultiPonStatus::Base),
            "ARTEFACT" => Ok(MultiPonStatus::Artefact),
            "MULTI" => Ok(MultiPonStatus::Multi),
            _ => anyhow::bail!("No enum constant com.hartwig.hmftools.common.variant.pon.MultiPonStatus.{s}"),
        }
    }
}

/// VariantTier
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum VariantTier {
    Unknown,
    LowConfidence,
    HighConfidence,
    Panel,
    Hotspot,
}

impl VariantTier {
    pub const TIER: &'static str = "TIER";

    pub fn from_string(s: &str) -> VariantTier {
        match s {
            "HOTSPOT" => VariantTier::Hotspot,
            "PANEL" => VariantTier::Panel,
            "HIGH_CONFIDENCE" => VariantTier::HighConfidence,
            "LOW_CONFIDENCE" => VariantTier::LowConfidence,
            _ => VariantTier::Unknown,
        }
    }

    /// VariantTier.fromContext
    pub fn from_context(rec: &VcfRecord) -> VariantTier {
        VariantTier::from_string(&rec.attr_string(Self::TIER, "UNKNOWN"))
    }
}

/// VariantType
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum VariantType {
    Snp,
    Mnp,
    Indel,
}

impl VariantType {
    pub fn name(self) -> &'static str {
        match self {
            VariantType::Snp => "SNP",
            VariantType::Mnp => "MNP",
            VariantType::Indel => "INDEL",
        }
    }
}

#[derive(Clone)]
pub struct VariantData {
    pub chromosome: String,
    pub position: i32,
    /// pos for an SNV, pos + alt length - 1 for an MNV, pos + deleted bases + 1 for a DEL (the first
    /// base after it), pos + 1 for an INS
    pub end_position: i32,
    pub ref_: String,
    pub alt: String,

    local_phase_set: i32,
    microhomology: String,
    repeat_count: i32,
    repeat_sequence: String,

    indel_base_diff: i32,
    /// the bases changed: none for an insert, the deleted ones for a DEL, all for SNVs and MNVs
    alt_positions: Vec<i32>,

    gene_impacts: JavaHashMap<String, Vec<VariantTransImpact>>,

    realigned_variant: Option<Box<VariantData>>,
    is_realigned_variant: bool,

    context: Option<VcfRecord>,

    pub filters: Vec<String>,
    pub pon_sample_count: i32,
    pub pon_max_read_count: i32,
    pub pon_mean_read_count: i32,
    pub multi_pon_status: MultiPonStatus,
    gnomad_frequency: Option<f64>,
    /// the tier, from the context (read once)
    tier: VariantTier,
}

#[allow(dead_code)]
impl VariantData {
    pub fn new(chromosome: &str, position: i32, ref_: &str, alt: &str) -> VariantData {
        let indel_base_diff = alt.len() as i32 - ref_.len() as i32;
        let mut alt_positions = Vec::new();
        let end_position;
        if indel_base_diff != 0 {
            // rare case where an INDEL has the ref base change as well, eg CT -> A
            let ref_base_changed = alt.as_bytes()[0] != ref_.as_bytes()[0];
            if indel_base_diff > 0 {
                if ref_base_changed {
                    alt_positions.push(position);
                }
                end_position = position + 1;
            } else {
                let count = indel_base_diff.abs();
                let start = if ref_base_changed { 0 } else { 1 };
                for i in start..ref_.len() as i32 {
                    alt_positions.push(position + i);
                }
                end_position = position + count + 1;
            }
        } else {
            for i in 0..ref_.len() as i32 {
                alt_positions.push(position + i);
            }
            end_position = position + ref_.len() as i32 - 1;
        }
        VariantData {
            chromosome: chromosome.to_string(),
            position,
            end_position,
            ref_: ref_.to_string(),
            alt: alt.to_string(),
            local_phase_set: NO_LOCAL_PHASE_SET,
            microhomology: String::new(),
            repeat_count: 0,
            repeat_sequence: String::new(),
            indel_base_diff,
            alt_positions,
            gene_impacts: JavaHashMap::new(),
            realigned_variant: None,
            is_realigned_variant: false,
            context: None,
            filters: Vec::new(),
            pon_sample_count: 0,
            pon_max_read_count: 0,
            pon_mean_read_count: 0,
            multi_pon_status: MultiPonStatus::Base,
            gnomad_frequency: None,
            tier: VariantTier::Unknown,
        }
    }

    /// VariantData.fromContext
    pub fn from_context(rec: VcfRecord) -> VariantData {
        let ref_ = rec.ref_allele.clone();
        // only the first of multiple alts is supported
        let mut alt = if rec.alts.is_empty() { ref_.clone() } else { rec.alts[0].clone() };
        if alt == "*" || alt == "N" {
            alt = ref_.clone();
        }
        let mut variant = VariantData::new(&rec.chrom, rec.pos, &ref_, &alt);

        let local_phase_sets = rec.attr_int_list(LOCAL_PHASE_SET, NO_LOCAL_PHASE_SET);
        let mut repeat_sequence = rec.attr_string(REPEAT_SEQUENCE, "");
        let mut repeat_count = rec.attr_int(REPEAT_COUNT, 0);

        // use the read context repeat if the ref's isn't set
        if repeat_count == 0 {
            let rc_repeat_sequence = rec.attr_string(READ_CONTEXT_REPEAT_SEQUENCE, "");
            let rc_repeat_count = rec.attr_int(READ_CONTEXT_REPEAT_COUNT, 0);
            let indel_bases = if variant.is_insert() { &alt[1..] } else { &ref_[1..] };
            if rc_repeat_count > 0 && indel_bases.starts_with(rc_repeat_sequence.as_str()) {
                repeat_count = rc_repeat_count - variant.base_diff() / rc_repeat_sequence.len() as i32;
                repeat_sequence = rc_repeat_sequence;
            }
        }

        let microhomology = rec.attr_string(MICROHOMOLOGY, "");
        variant.set_variant_details(
            local_phase_sets.first().copied().unwrap_or(NO_LOCAL_PHASE_SET),
            &microhomology,
            &repeat_sequence,
            repeat_count,
        );
        variant.tier = VariantTier::from_context(&rec);
        variant.context = Some(rec);
        variant
    }

    pub fn variant_type(&self) -> VariantType {
        if self.indel_base_diff == 0 {
            return if self.ref_.len() == 1 { VariantType::Snp } else { VariantType::Mnp };
        }
        VariantType::Indel
    }

    pub fn base_diff(&self) -> i32 {
        self.indel_base_diff
    }
    pub fn is_base_change(&self) -> bool {
        self.indel_base_diff == 0
    }
    pub fn is_indel(&self) -> bool {
        self.indel_base_diff != 0
    }
    pub fn is_insert(&self) -> bool {
        self.indel_base_diff > 0
    }
    pub fn is_deletion(&self) -> bool {
        self.indel_base_diff < 0
    }
    pub fn is_mixed(&self) -> bool {
        self.indel_base_diff != 0 && self.alt.as_bytes()[0] != self.ref_.as_bytes()[0]
    }
    pub fn is_snv(&self) -> bool {
        self.indel_base_diff == 0 && self.ref_.len() == 1
    }
    pub fn is_mnv(&self) -> bool {
        self.indel_base_diff == 0 && self.ref_.len() > 1
    }

    pub fn alt_positions(&self) -> &[i32] {
        &self.alt_positions
    }

    pub fn realigned_variant(&self) -> Option<&VariantData> {
        self.realigned_variant.as_deref()
    }
    pub fn realigned_variant_mut(&mut self) -> Option<&mut VariantData> {
        self.realigned_variant.as_deref_mut()
    }
    pub fn set_realigned_variant(&mut self, variant: Option<VariantData>) {
        if let Some(mut v) = variant {
            v.is_realigned_variant = true;
            self.realigned_variant = Some(Box::new(v));
        }
    }
    pub fn set_realigned_variant_none(&mut self) {
        self.realigned_variant = None;
    }
    pub fn is_realigned_variant(&self) -> bool {
        self.is_realigned_variant
    }
    pub fn mark_realigned_variant(&mut self) {
        self.is_realigned_variant = true;
    }

    pub fn context(&self) -> &VcfRecord {
        self.context.as_ref().expect("variant without a VCF record")
    }
    pub fn context_mut(&mut self) -> &mut VcfRecord {
        self.context.as_mut().expect("variant without a VCF record")
    }
    pub fn has_context(&self) -> bool {
        self.context.is_some()
    }
    pub fn set_context(&mut self, rec: VcfRecord) {
        self.tier = VariantTier::from_context(&rec);
        self.context = Some(rec);
    }
    pub fn take_context(&mut self) -> Option<VcfRecord> {
        self.context.take()
    }

    pub fn set_variant_details(&mut self, local_phase_set: i32, microhomology: &str, repeat_sequence: &str, repeat_count: i32) {
        self.local_phase_set = local_phase_set;
        self.microhomology = microhomology.to_string();
        self.repeat_sequence = repeat_sequence.to_string();
        self.repeat_count = repeat_count;
    }

    pub fn has_local_phase_set(&self) -> bool {
        self.local_phase_set != NO_LOCAL_PHASE_SET
    }
    pub fn local_phase_set(&self) -> i32 {
        self.local_phase_set
    }
    pub fn microhomology(&self) -> &str {
        &self.microhomology
    }
    pub fn repeat_sequence(&self) -> &str {
        &self.repeat_sequence
    }
    pub fn repeat_count(&self) -> i32 {
        self.repeat_count
    }

    pub fn impacts(&self) -> &JavaHashMap<String, Vec<VariantTransImpact>> {
        &self.gene_impacts
    }
    pub fn impacts_mut(&mut self) -> &mut JavaHashMap<String, Vec<VariantTransImpact>> {
        &mut self.gene_impacts
    }
    pub fn gene_impacts(&self, gene_name: &str) -> Option<&Vec<VariantTransImpact>> {
        self.gene_impacts.get(&gene_name.to_string())
    }

    /// getCanonicalTransImpacts
    pub fn canonical_trans_impact(&self, gene_name: &str) -> Option<&VariantTransImpact> {
        self.gene_impacts(gene_name)?.iter().find(|x| x.trans.is_canonical)
    }

    /// addImpact: ignored without effects or when the gene already has an impact for the transcript
    pub fn add_impact(&mut self, gene_name: &str, impact: VariantTransImpact) {
        if impact.effects().is_empty() {
            return;
        }
        let key = gene_name.to_string();
        match self.gene_impacts.get_mut(&key) {
            None => {
                self.gene_impacts.insert(key, vec![impact]);
            }
            Some(list) => {
                if list.iter().any(|x| x.trans.trans_id == impact.trans.trans_id) {
                    return;
                }
                list.push(impact);
            }
        }
    }

    /// getRealignedImpact: the realigned variant's impact on the same transcript
    pub fn realigned_impact(&self, gene_name: &str, ti: &VariantTransImpact) -> Option<&VariantTransImpact> {
        self.realigned_variant.as_ref()?.gene_impacts(gene_name)?.iter().find(|x| x.trans.trans_name == ti.trans.trans_name)
    }

    pub fn alt_bases_above(&self, position: i32) -> bool {
        self.alt_bases_outside_position(position, true)
    }
    pub fn alt_bases_below(&self, position: i32) -> bool {
        self.alt_bases_outside_position(position, false)
    }

    /// whether all altered positions are at or inside the limits
    pub fn alt_positions_within(&self, lower: i32, upper: i32) -> bool {
        if self.is_insert() {
            return positions_within(self.position, self.end_position, lower, upper);
        }
        self.alt_bases_above(lower - 1) && self.alt_bases_below(upper + 1)
    }

    pub fn alt_positions_overlap(&self, lower: i32, upper: i32) -> bool {
        !self.alt_bases_below(lower) && !self.alt_bases_above(upper)
    }

    /// whether an altered base is beyond a position, above or below it
    pub fn alt_bases_outside_position(&self, position: i32, is_higher: bool) -> bool {
        if self.indel_base_diff == 0 {
            return if is_higher { self.position > position } else { self.end_position < position };
        }
        if self.is_insert() {
            if is_higher { self.end_position > position } else { self.position < position }
        } else {
            // the first base of a DEL is the ref
            if is_higher { self.position + 1 > position } else { self.end_position - 1 < position }
        }
    }

    pub fn tier(&self) -> VariantTier {
        if self.context.is_none() { VariantTier::Unknown } else { self.tier }
    }

    /// the sample's AF, 0 without the sample or the field
    pub fn sample_vaf(&self, sample_id: &str) -> f64 {
        match self.context().genotype(sample_id) {
            Some(g) => g.attr_double("AF", 0.0),
            None => 0.0,
        }
    }

    pub fn add_filter(&mut self, filter: &str) {
        if !self.filters.iter().any(|f| f == filter) {
            self.filters.push(filter.to_string());
        }
    }

    pub fn set_pon_frequency(&mut self, sample_count: i32, max_read_count: i32, mean_read_count: i32, status: MultiPonStatus) {
        self.pon_sample_count = sample_count;
        self.pon_max_read_count = max_read_count;
        self.pon_mean_read_count = mean_read_count;
        self.multi_pon_status = status;
    }

    pub fn gnomad_frequency(&self) -> Option<f64> {
        self.gnomad_frequency
    }
    pub fn set_gnomad_frequency(&mut self, f: Option<f64>) {
        self.gnomad_frequency = f;
    }

    pub fn tsv_header() -> &'static str {
        "Chromosome\tPosition\tType\tRef\tAlt\tLPS"
    }

    pub fn tsv_data(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}", self.chromosome, self.position, self.variant_type().name(), self.ref_, self.alt, self.local_phase_set)
    }

    pub fn extra_data_header() -> &'static str {
        "Microhomology\tRepeatSequence\tRepeatCount\tLocalPhaseSetId\tFilter\tTier\tQual\tAlleleReadCount\tTotalReadCount"
    }

    /// extraDataTsv: with VariantContextDecorator's filter, tier, qual and allelic depth of the sample
    pub fn extra_data_tsv(&self, sample_id: &str) -> String {
        let rec = self.context();
        // VariantContextDecorator.filter: PASS when unfiltered, else the filters joined by ';'
        let filter = if rec.filters.is_empty() { "PASS".to_string() } else { rec.filters.join(";") };
        let (allele_read_count, total_read_count) = allelic_depth(rec, sample_id);
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.microhomology,
            self.repeat_sequence,
            self.repeat_count,
            self.local_phase_set,
            filter,
            tier_name(VariantTier::from_context(rec)),
            hmf_common::fmt::java_double_to_string(rec.phred_scaled_qual()),
            allele_read_count,
            total_read_count
        )
    }
}

pub fn tier_name(t: VariantTier) -> &'static str {
    match t {
        VariantTier::Unknown => "UNKNOWN",
        VariantTier::LowConfidence => "LOW_CONFIDENCE",
        VariantTier::HighConfidence => "HIGH_CONFIDENCE",
        VariantTier::Panel => "PANEL",
        VariantTier::Hotspot => "HOTSPOT",
    }
}

/// VariantContextDecorator.allelicDepth: AD[1] and the total (the larger of DP and the sum of AD, or the
/// sum without DP); no sample gives zeros, a sample without AD fails as Java's precondition does
fn allelic_depth(rec: &VcfRecord, sample_id: &str) -> (i32, i32) {
    let Some(g) = rec.genotype(sample_id) else { return (0, 0) };
    let Some(ad_text) = g.ad.as_deref() else { panic!("IllegalArgumentException: genotype of {sample_id} without AD") };
    let ad: Vec<i32> = ad_text.split(',').map(|x| x.parse().unwrap_or(0)).collect();
    let alt = ad.get(1).copied().unwrap_or(0);
    let sum: i32 = ad.iter().sum();
    let total = match g.dp.as_deref().and_then(|d| d.parse::<i32>().ok()) {
        Some(dp) => dp.max(sum),
        None => sum,
    };
    (alt, total)
}

impl std::fmt::Display for VariantData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.indel_base_diff == 0 && self.ref_.len() == 1 {
            write!(f, "pos({}:{}) variant({}>{})", self.chromosome, self.position, self.ref_, self.alt)
        } else {
            write!(f, "pos({}:{}-{}) variant({}>{})", self.chromosome, self.position, self.end_position, self.ref_, self.alt)
        }
    }
}
