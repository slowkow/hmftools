//! The impact package: a variant's effect on each transcript it touches (coding context, protein
//! context, HGVS strings, splice effects), the phasing of nearby variants, and the variant's summary
//! impact (canonical and worst).

pub mod builder;
pub mod classifier;
pub mod coding;
pub mod effects;
pub mod hgvs_coding;
pub mod hgvs_protein;
pub mod phased;
pub mod protein;
pub mod special_splice;
pub mod splice;
pub mod utils;

use anyhow::{bail, Result};

use hmf_common::ensembl::TranscriptData;

use effects::{effects_to_string, CodingType, RegionType, VariantEffect};

/// RefGenomeInterface.getBaseString: None (Java's null) for a start before 1 or an end before the start
pub trait RefGenome: Sync {
    fn base_string(&self, chromosome: &str, start: i32, end: i32) -> Result<Option<String>>;
}

/// The bases as a String, Java's "null" for a null String (as `+=` appends it)
pub fn or_null(s: Option<String>) -> String {
    s.unwrap_or_else(|| "null".to_string())
}

/// String.substring(begin, end), with Java's exception for an invalid range
pub fn jsub(s: &str, begin: i32, end: i32) -> Result<&str> {
    if begin < 0 || end > s.len() as i32 || begin > end {
        bail!("StringIndexOutOfBoundsException: begin {begin}, end {end}, length {}", s.len());
    }
    Ok(&s[begin as usize..end as usize])
}

/// String.substring(begin)
pub fn jsub_from(s: &str, begin: i32) -> Result<&str> {
    jsub(s, begin, s.len() as i32)
}

/// String.charAt
pub fn char_at(s: &str, i: i32) -> Result<u8> {
    if i < 0 || i >= s.len() as i32 {
        bail!("StringIndexOutOfBoundsException: index {i}, length {}", s.len());
    }
    Ok(s.as_bytes()[i as usize])
}

pub const PHASE_NONE: i32 = -1;
pub const PHASE_0: i32 = 0;
pub const PHASE_1: i32 = 1;
pub const PHASE_2: i32 = 2;

pub const SE_START: usize = 0;
pub const SE_END: usize = 1;

pub const POS_STRAND: i8 = 1;
pub const NEG_STRAND: i8 = -1;

/// SpliceImpactType
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SpliceImpactType {
    Unknown,
    HomologyShift,
    BaseShift,
    BaseChange,
    RegionDeleted,
    OutsideRange,
}

impl SpliceImpactType {
    pub fn is_disruptive(self) -> bool {
        matches!(self, SpliceImpactType::BaseChange | SpliceImpactType::BaseShift | SpliceImpactType::RegionDeleted | SpliceImpactType::HomologyShift)
    }
    pub fn name(self) -> &'static str {
        match self {
            SpliceImpactType::Unknown => "UNKNOWN",
            SpliceImpactType::HomologyShift => "HOMOLOGY_SHIFT",
            SpliceImpactType::BaseShift => "BASE_SHIFT",
            SpliceImpactType::BaseChange => "BASE_CHANGE",
            SpliceImpactType::RegionDeleted => "REGION_DELETED",
            SpliceImpactType::OutsideRange => "OUTSIDE_RANGE",
        }
    }
}

/// CodingContext
#[derive(Clone, Debug)]
pub struct CodingContext {
    /// favours the more impactful type if the variant spans several
    pub region_type: RegionType,
    pub coding_type: CodingType,
    pub exon_rank: i32,
    pub strand: i8,
    pub upstream_phase: i32,
    /// indexed at 1; for non-coding and UTR regions, the exonic bases from the coding region
    pub coding_base: i32,
    /// the alt bases of an MNV, otherwise the start and end ref bases (enclosing the DEL or INS);
    /// likewise the exonic base range for non-coding and UTR regions
    pub coding_position_range: [i32; 2],
    pub nearest_exon_distance: i32,
    pub spans_splice_junction: bool,
    pub is_frameshift: bool,
    pub spans_coding_start: bool,
    pub spans_coding_end: bool,
    pub deleted_coding_bases: i32,
    pub splice_donor_acceptor_bases: String,
    /// for UTR positions in the intron immediately before or after the coding region
    pub coding_ends_on_exon_boundary: bool,
    pub hgvs: String,
}

impl Default for CodingContext {
    fn default() -> Self {
        CodingContext {
            region_type: RegionType::Unknown,
            coding_type: CodingType::Unknown,
            exon_rank: 0,
            strand: 0,
            upstream_phase: PHASE_NONE,
            coding_base: 0,
            coding_position_range: [0, 0],
            nearest_exon_distance: 0,
            spans_splice_junction: false,
            is_frameshift: false,
            spans_coding_start: false,
            spans_coding_end: false,
            deleted_coding_bases: 0,
            splice_donor_acceptor_bases: String::new(),
            coding_ends_on_exon_boundary: false,
            hgvs: String::new(),
        }
    }
}

impl CodingContext {
    pub fn is_coding(&self) -> bool {
        self.coding_type == CodingType::Coding && self.region_type == RegionType::Exonic
    }

    pub fn coding_ref(&self, variant: &crate::variant::VariantData) -> Result<String> {
        if variant.is_insert() {
            return Ok(variant.ref_.clone());
        }
        let end = if variant.is_deletion() { variant.end_position - 1 } else { variant.end_position };
        self.trim_to_coding(&variant.ref_, variant.position, end)
    }

    pub fn coding_alt(&self, variant: &crate::variant::VariantData) -> Result<String> {
        if variant.is_insert() {
            return Ok(variant.alt.clone());
        }
        if variant.is_deletion() {
            // trimmed to empty if outside the coding region
            if hmf_common::regions::position_within(variant.position, self.coding_position_range[SE_START], self.coding_position_range[SE_END]) {
                return Ok(variant.alt.clone());
            }
            return Ok(String::new());
        }
        self.trim_to_coding(&variant.alt, variant.position, variant.end_position)
    }

    fn trim_to_coding(&self, bases: &str, pos_start: i32, pos_end: i32) -> Result<String> {
        let r = &self.coding_position_range;
        let start_diff = if r[SE_START] > pos_start { r[SE_START] - pos_start } else { 0 };
        let end_diff = if r[SE_END] > 0 && r[SE_END] < pos_end { pos_end - r[SE_END] } else { 0 };
        Ok(jsub(bases, start_diff, bases.len() as i32 - end_diff)?.to_string())
    }

    pub fn tsv_header() -> &'static str {
        "HgvsCoding\tRegionType\tCodingType\tExonRank\tCodingBase\tCodingPosRange\tUpstreamPhase\tSpansSplice\tSpliceDonorAcceptorBases\tNearestExonDistance"
    }

    pub fn to_tsv(&self) -> String {
        let r = &self.coding_position_range;
        let range = if r[SE_START] == r[SE_END] { r[SE_START].to_string() } else { format!("{}-{}", r[SE_START], r[SE_END]) };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.hgvs,
            self.region_type.name(),
            self.coding_type.name(),
            self.exon_rank,
            self.coding_base,
            range,
            self.upstream_phase,
            self.spans_splice_junction,
            self.splice_donor_acceptor_bases,
            self.nearest_exon_distance
        )
    }
}

/// ProteinContext
#[derive(Clone, Debug, Default)]
pub struct ProteinContext {
    /// coding bases expanded to whole codons
    pub ref_codon_bases: String,
    /// as above with the alt in place of the ref
    pub alt_codon_bases: String,
    /// the ranges of the ref codon bases, used for phasing variants
    pub ref_codons_ranges: Vec<[i32; 2]>,
    /// alt codons plus any downstream refs to complete a codon, also extended for inframe
    pub alt_codon_bases_complete: String,
    /// ref codons extended for inframe variants until a different amino acid is found
    pub ref_codon_bases_extended: String,
    /// amino acid index of the ref codon
    pub codon_index: i32,
    pub ref_amino_acids: String,
    pub alt_amino_acids: String,
    pub extra_upstream_codon: bool,
    /// the ref and alt amino acids without the matching ones at either end
    pub net_ref_amino_acids: String,
    pub net_alt_amino_acids: String,
    pub net_codon_index_range: [i32; 2],
    pub is_duplication: bool,
    pub is_phased: bool,
    pub hgvs: String,
}

impl ProteinContext {
    pub fn has_coding_bases(&self) -> bool {
        !self.ref_codon_bases.is_empty()
    }
    pub fn has_amino_acids(&self) -> bool {
        !self.ref_amino_acids.is_empty()
    }
    pub fn has_protein_change(&self) -> bool {
        self.ref_amino_acids != self.alt_amino_acids
    }
    pub fn valid_ref_codon(&self) -> bool {
        !self.ref_codon_bases.is_empty() && self.ref_codon_bases.len() % 3 == 0
    }
    pub fn valid_alt_codon(&self) -> bool {
        self.alt_codon_bases_complete.len() % 3 == 0
    }
    pub fn ref_coding_base_start(&self) -> i32 {
        self.ref_coding_base_position(SE_START)
    }
    pub fn ref_coding_base_end(&self) -> i32 {
        self.ref_coding_base_position(SE_END)
    }
    pub fn ref_coding_base_position(&self, se: usize) -> i32 {
        if self.ref_codons_ranges.is_empty() {
            return -1;
        }
        if se == SE_START { self.ref_codons_ranges[0][SE_START] } else { self.ref_codons_ranges[self.ref_codons_ranges.len() - 1][SE_END] }
    }

    pub fn tsv_header() -> &'static str {
        "HgvsProtein\tRefCodonBases\tAltCodonBases\tCodonIndex\tRefAA\tAltAA\tNetAARange\tNetRefAA\tNetAltAA\tIsDup"
    }

    pub fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}-{}\t{}\t{}\t{}",
            self.hgvs,
            self.ref_codon_bases,
            self.alt_codon_bases,
            self.codon_index,
            self.ref_amino_acids,
            self.alt_amino_acids,
            self.net_codon_index_range[SE_START],
            self.net_codon_index_range[SE_END],
            self.net_ref_amino_acids,
            self.net_alt_amino_acids,
            self.is_duplication
        )
    }

    pub fn empty_tsv() -> &'static str {
        ",,,0,,,0-0,,,"
    }
}

/// VariantTransImpact: a variant's impact on one transcript
#[derive(Clone, Debug)]
pub struct VariantTransImpact {
    pub trans: &'static TranscriptData,
    effects: Vec<VariantEffect>,
    pub coding: CodingContext,
    pub protein: Option<ProteinContext>,
    in_splice_region: bool,
    realigned: bool,
    splice_impact_type: SpliceImpactType,
    phased_frameshift: bool,
}

impl VariantTransImpact {
    pub fn new(trans: &'static TranscriptData) -> VariantTransImpact {
        VariantTransImpact {
            trans,
            effects: Vec::new(),
            coding: CodingContext::default(),
            protein: None,
            in_splice_region: false,
            realigned: false,
            splice_impact_type: SpliceImpactType::Unknown,
            phased_frameshift: false,
        }
    }

    /// addEffect: once each, kept in order of rank (a new effect goes after those of the same rank)
    pub fn add_effect(&mut self, effect: VariantEffect) {
        if self.effects.contains(&effect) {
            return;
        }
        let mut index = 0;
        while index < self.effects.len() {
            if effect.rank() > self.effects[index].rank() {
                break;
            }
            index += 1;
        }
        self.effects.insert(index, effect);
    }

    pub fn effects(&self) -> &[VariantEffect] {
        &self.effects
    }
    pub fn remove_effect(&mut self, effect: VariantEffect) {
        if let Some(i) = self.effects.iter().position(|e| *e == effect) {
            self.effects.remove(i);
        }
    }
    pub fn clear_effects(&mut self) {
        self.effects.clear();
    }

    pub fn top_effect(&self) -> Option<VariantEffect> {
        self.effects.first().copied()
    }
    pub fn has_effect(&self, effect: VariantEffect) -> bool {
        self.effects.contains(&effect)
    }
    pub fn top_rank(&self) -> i32 {
        self.effects.iter().map(|e| e.rank()).max().unwrap_or(-1)
    }

    pub fn is_exonic(&self) -> bool {
        self.coding.region_type == RegionType::Exonic
    }

    pub fn has_coding_bases(&self) -> bool {
        self.protein.as_ref().is_some_and(|p| p.has_coding_bases())
    }
    pub fn has_amino_acids(&self) -> bool {
        self.protein.as_ref().is_some_and(|p| p.has_amino_acids())
    }
    pub fn has_protein_context(&self) -> bool {
        self.protein.is_some()
    }

    /// the protein context, for tests
    #[cfg(test)]
    pub fn pc(&self) -> &ProteinContext {
        self.protein.as_ref().expect("no protein context")
    }

    pub fn mark_splice_region(&mut self) {
        self.in_splice_region = true;
    }
    pub fn in_splice_region(&self) -> bool {
        self.in_splice_region
    }
    pub fn mark_phased_frameshift(&mut self) {
        self.phased_frameshift = true;
    }
    #[allow(dead_code)]
    pub fn phased_frameshift(&self) -> bool {
        self.phased_frameshift
    }
    pub fn mark_realigned(&mut self) {
        self.realigned = true;
    }
    pub fn realigned(&self) -> bool {
        self.realigned
    }
    pub fn set_splice_impact_type(&mut self, t: SpliceImpactType) {
        self.splice_impact_type = t;
    }
    pub fn splice_impact_type(&self) -> SpliceImpactType {
        self.splice_impact_type
    }

    pub fn hgvs_coding(&self) -> &str {
        &self.coding.hgvs
    }
    pub fn hgvs_protein(&self) -> &str {
        self.protein.as_ref().map(|p| p.hgvs.as_str()).unwrap_or("")
    }

    pub fn effects_str(&self) -> String {
        effects_to_string(&self.effects, "&")
    }
    pub fn effects_to_csv(&self) -> String {
        effects_to_string(&self.effects, ";")
    }

    pub fn tsv_header() -> &'static str {
        "TransId\tCanonical\tIsCoding\tStrand\tSpliceRegion\tSpliceImpact\tPhasedInframe\tRealigned\tEffects"
    }

    pub fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.trans.trans_name,
            self.trans.is_canonical,
            !self.trans.non_coding(),
            self.trans.strand,
            self.in_splice_region,
            self.splice_impact_type.name(),
            self.phased_frameshift,
            self.realigned,
            self.effects_to_csv()
        )
    }
}

/// TranscriptData.CodingStart / CodingEnd, which callers use only for coding transcripts (Java unboxes
/// the null Integer otherwise)
pub fn coding_start(t: &TranscriptData) -> Result<i32> {
    t.coding_start.ok_or_else(|| anyhow::anyhow!("NullPointerException: coding start of non-coding transcript {}", t.trans_name))
}
pub fn coding_end(t: &TranscriptData) -> Result<i32> {
    t.coding_end.ok_or_else(|| anyhow::anyhow!("NullPointerException: coding end of non-coding transcript {}", t.trans_name))
}
