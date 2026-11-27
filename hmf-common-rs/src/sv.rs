//! hmf-common `sv`: structural variants read from an ESVEE or PURPLE SV VCF, as LINX (and PURPLE) read
//! them: StructuralVariantType, StructuralVariantFactory (breakend pairs by MATEID, single breakends,
//! the legs' fields), EnrichedStructuralVariantFactory (the PURPLE_* INFO fields) and
//! StructuralVariantData (the flat record of both legs) with its `convertSvData`.
//!
//! The VCF is read as text (plain or gzip/BGZF) and its INFO and FORMAT values are taken the way htsjdk's
//! VariantContext hands them out: `getAttributeAsString` of a value with commas is the list printed as
//! `[a, b]`, `getAttributeAsDouble` / `getAttributeAsInt` parse with Java's rules, "PASS" is an empty
//! filter set, the QUAL "." is -10.

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;

use crate::java_text::{parse_double, parse_int};

/// StructuralVariantType, in declaration order (the ordinal is used for type counts)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StructuralVariantType {
    BND,
    DEL,
    DUP,
    INF,
    INS,
    INV,
    SGL,
}

impl StructuralVariantType {
    pub const VALUES: [StructuralVariantType; 7] = [Self::BND, Self::DEL, Self::DUP, Self::INF, Self::INS, Self::INV, Self::SGL];

    pub fn ordinal(self) -> usize {
        self as usize
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::BND => "BND",
            Self::DEL => "DEL",
            Self::DUP => "DUP",
            Self::INF => "INF",
            Self::INS => "INS",
            Self::INV => "INV",
            Self::SGL => "SGL",
        }
    }
    /// valueOf: None where Java throws IllegalArgumentException
    pub fn from_name(s: &str) -> Option<Self> {
        Self::VALUES.iter().copied().find(|t| t.name() == s)
    }
}

impl std::fmt::Display for StructuralVariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// One VCF record, with the parts the SV readers use
#[derive(Clone, Debug)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: i32,
    pub id: String,
    pub ref_bases: String,
    pub alt: String,
    /// getPhredScaledQual: -10 for "."
    pub qual: f64,
    /// getFilters: empty for PASS and "."
    pub filters: Vec<String>,
    info: Vec<(String, Option<String>)>,
    format: Vec<String>,
    samples: Vec<Vec<String>>,
}

impl VcfRecord {
    /// parses a data line of a VCF (tab-separated, eight fixed columns, then FORMAT and the samples)
    pub fn parse(line: &str) -> Result<VcfRecord> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            bail!("VCF line with fewer than 8 columns: {line}");
        }
        let pos = parse_int(cols[1]).ok_or_else(|| anyhow!("invalid VCF position: {}", cols[1]))?;
        let qual = if cols[5] == "." { -10.0 } else { parse_double(cols[5]).ok_or_else(|| anyhow!("invalid VCF QUAL: {}", cols[5]))? };
        let filters: Vec<String> =
            if cols[6] == "." || cols[6] == "PASS" { Vec::new() } else { cols[6].split(';').map(|s| s.to_string()).collect() };
        let mut info = Vec::new();
        if cols[7] != "." && !cols[7].is_empty() {
            for item in cols[7].split(';') {
                if item.is_empty() {
                    continue;
                }
                match item.split_once('=') {
                    Some((k, v)) => info.push((k.to_string(), Some(v.to_string()))),
                    None => info.push((item.to_string(), None)),
                }
            }
        }
        let format: Vec<String> = if cols.len() > 8 { cols[8].split(':').map(|s| s.to_string()).collect() } else { Vec::new() };
        let samples: Vec<Vec<String>> = cols.iter().skip(9).map(|c| c.split(':').map(|s| s.to_string()).collect()).collect();
        Ok(VcfRecord {
            chrom: cols[0].to_string(),
            pos,
            id: cols[2].to_string(),
            ref_bases: cols[3].to_string(),
            alt: cols[4].to_string(),
            qual,
            filters,
            info,
            format,
            samples,
        })
    }

    fn raw(&self, key: &str) -> Option<&Option<String>> {
        self.info.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    pub fn has_attribute(&self, key: &str) -> bool {
        self.raw(key).is_some()
    }

    /// getAttributeAsString(key, default): a list value (with commas) as Java's List.toString prints it
    pub fn attribute_as_string(&self, key: &str, default: &str) -> String {
        match self.raw(key) {
            None => default.to_string(),
            // a flag's value is Boolean.TRUE
            Some(None) => "true".to_string(),
            Some(Some(v)) => {
                if v.contains(',') {
                    format!("[{}]", v.split(',').collect::<Vec<_>>().join(", "))
                } else {
                    v.clone()
                }
            }
        }
    }

    /// getAttributeAsString(key, null)
    pub fn attribute_as_string_opt(&self, key: &str) -> Option<String> {
        if self.has_attribute(key) { Some(self.attribute_as_string(key, "")) } else { None }
    }

    /// getAttributeAsDouble(key, default)
    pub fn attribute_as_double(&self, key: &str, default: f64) -> Result<f64> {
        match self.raw(key) {
            None => Ok(default),
            Some(None) => bail!("attribute {key} is a flag, not a double"),
            Some(Some(v)) => {
                let first = v.split(',').next().unwrap_or("");
                parse_double(first).ok_or_else(|| anyhow!("NumberFormatException: attribute {key}={v}"))
            }
        }
    }

    /// getAttributeAsInt(key, default)
    pub fn attribute_as_int(&self, key: &str, default: i32) -> Result<i32> {
        match self.raw(key) {
            None => Ok(default),
            Some(None) => bail!("attribute {key} is a flag, not an integer"),
            Some(Some(v)) => {
                let first = v.split(',').next().unwrap_or("");
                parse_int(first).ok_or_else(|| anyhow!("NumberFormatException: attribute {key}={v}"))
            }
        }
    }

    /// getAttributeAsDoubleList(key, default): "." gives the default
    pub fn attribute_as_double_list(&self, key: &str, default: f64) -> Result<Vec<f64>> {
        match self.raw(key) {
            None | Some(None) => Ok(Vec::new()),
            Some(Some(v)) => v
                .split(',')
                .map(|x| if x == "." { Ok(default) } else { parse_double(x).ok_or_else(|| anyhow!("NumberFormatException: {key}={v}")) })
                .collect(),
        }
    }

    /// getAttributeAsIntList(key, default)
    pub fn attribute_as_int_list(&self, key: &str, default: i32) -> Result<Vec<i32>> {
        match self.raw(key) {
            None | Some(None) => Ok(Vec::new()),
            Some(Some(v)) => v
                .split(',')
                .map(|x| if x == "." { Ok(default) } else { parse_int(x).ok_or_else(|| anyhow!("NumberFormatException: {key}={v}")) })
                .collect(),
        }
    }

    /// getAttributeAsBoolean(key, default)
    pub fn attribute_as_bool(&self, key: &str, default: bool) -> bool {
        match self.raw(key) {
            None => default,
            Some(None) => true,
            Some(Some(v)) => v == "1" || v.eq_ignore_ascii_case("true"),
        }
    }

    pub fn genotype_count(&self) -> usize {
        self.samples.len()
    }

    /// Genotype.getExtendedAttribute(key) of the sample at `ordinal`, as its text
    pub fn genotype_attribute(&self, ordinal: usize, key: &str) -> Option<&str> {
        let sample = self.samples.get(ordinal)?;
        let i = self.format.iter().position(|f| f == key)?;
        let v = sample.get(i)?;
        // htsjdk leaves out missing values
        if v == "." || v.is_empty() { None } else { Some(v.as_str()) }
    }

    /// CommonVcfTags.getGenotypeAttributeAsInt
    pub fn genotype_attribute_as_int(&self, ordinal: usize, key: &str, default: i32) -> Result<i32> {
        match self.genotype_attribute(ordinal, key) {
            None => Ok(default),
            Some(v) => parse_int(v).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{v}\"")),
        }
    }

    /// CommonVcfTags.getGenotypeAttributeAsDouble
    pub fn genotype_attribute_as_double(&self, ordinal: usize, key: &str, default: f64) -> Result<f64> {
        match self.genotype_attribute(ordinal, key) {
            None => Ok(default),
            Some(v) => parse_double(v).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{v}\"")),
        }
    }
}

/// A VCF file: the sample names of the #CHROM line and the records
pub struct VcfFile {
    pub samples: Vec<String>,
    pub records: Vec<VcfRecord>,
}

pub fn read_vcf(path: &str) -> Result<VcfFile> {
    let mut samples = Vec::new();
    let mut records = Vec::new();
    crate::java_text::for_each_line(path, |line| {
        if line.starts_with("##") {
            return Ok(());
        }
        if let Some(h) = line.strip_prefix('#') {
            samples = h.split('\t').skip(9).map(|s| s.to_string()).collect();
            return Ok(());
        }
        if line.is_empty() {
            return Ok(());
        }
        records.push(VcfRecord::parse(line)?);
        Ok(())
    })?;
    Ok(VcfFile { samples, records })
}

// SvVcfTags
pub const SV_TYPE: &str = "SVTYPE";
pub const MATE_ID: &str = "MATEID";
pub const CIPOS: &str = "CIPOS";
pub const HOMSEQ: &str = "HOMSEQ";
pub const IHOMPOS: &str = "IHOMPOS";
pub const INSALN: &str = "INSALN";
pub const SEG_ALIGN_LENGTH: &str = "SEGALEN";
pub const ASM_LINKS: &str = "ASMLNKS";
pub const TOTAL_FRAGS: &str = "VF";
pub const REF_DEPTH: &str = "REF";
pub const REF_DEPTH_PAIR: &str = "REFPAIR";
pub const ALLELE_FRACTION: &str = "AF";
pub const PON_FILTER_PON: &str = "PON";
pub const PON_COUNT: &str = "PON_COUNT";
pub const HOTSPOT: &str = "HOTSPOT";
pub const REPEAT_MASK_REPEAT_CLASS: &str = "INSRMRC";
pub const REPEAT_MASK_REPEAT_TYPE: &str = "INSRMRT";
pub const REPEAT_MASK_ORIENTATION: &str = "INSRMRO";
pub const REPEAT_MASK_COVERAGE: &str = "INSRMP";
pub const INFERRED: &str = "INFERRED";
pub const PASS_FILTER: &str = "PASS";
// PURPLE's annotations
pub const PURPLE_AF: &str = "PURPLE_AF";
pub const PURPLE_CN: &str = "PURPLE_CN";
pub const PURPLE_CN_CHANGE: &str = "PURPLE_CN_CHANGE";
pub const PURPLE_JUNCTION_COPY_NUMBER: &str = "PURPLE_JCN";
pub const PURPLE_PLOIDY_INFO: &str = "PURPLE_PLOIDY";

/// HumanChromosome.contains: 1-22, X or Y, with or without the chr prefix
pub fn is_human_chromosome(chromosome: &str) -> bool {
    let c = crate::regions::strip_chr_prefix(chromosome);
    if !c.is_empty() && c.bytes().all(|b| b.is_ascii_digit()) {
        return c.parse::<i64>().is_ok_and(|v| (1..=22).contains(&v));
    }
    c == "X" || c == "Y"
}

/// StructuralVariantLeg: one breakend's fields
#[derive(Clone, Debug, Default)]
pub struct StructuralVariantLeg {
    pub chromosome: String,
    pub position: i32,
    pub orientation: i8,
    pub homology: String,
    pub allele_frequency: Option<f64>,
    pub tumor_variant_fragment_count: Option<i32>,
    pub tumor_reference_fragment_count: Option<i32>,
    pub normal_variant_fragment_count: Option<i32>,
    pub normal_reference_fragment_count: Option<i32>,
    pub start_offset: i32,
    pub end_offset: i32,
    pub inexact_homology_offset_start: i32,
    pub inexact_homology_offset_end: i32,
    pub anchoring_support_distance: i32,
    // EnrichedStructuralVariantLeg
    pub adjusted_allele_frequency: Option<f64>,
    pub adjusted_copy_number: Option<f64>,
    pub adjusted_copy_number_change: Option<f64>,
}

/// StructuralVariant (with the enriched fields): its legs and VCF records
#[derive(Clone, Debug)]
pub struct StructuralVariant {
    pub id: String,
    pub mate_id: Option<String>,
    pub start: StructuralVariantLeg,
    pub end: Option<StructuralVariantLeg>,
    pub insert_sequence: String,
    pub sv_type: StructuralVariantType,
    pub filter: String,
    pub hotspot: bool,
    pub quality_score: f64,
    pub event: String,
    pub start_linked_by: String,
    pub end_linked_by: String,
    pub insert_sequence_alignments: String,
    pub insert_sequence_repeat_class: Option<String>,
    pub insert_sequence_repeat_type: Option<String>,
    pub insert_sequence_repeat_orientation: Option<i8>,
    pub insert_sequence_repeat_coverage: Option<f64>,
    pub junction_copy_number: Option<f64>,
    pub start_context: VcfRecord,
    pub end_context: Option<VcfRecord>,
}

impl StructuralVariant {
    pub fn chromosome(&self, is_start: bool) -> &str {
        if is_start { &self.start.chromosome } else { &self.end.as_ref().unwrap().chromosome }
    }
}

/// StructuralVariantData: the flat record LINX works with
#[derive(Clone, Debug, Default)]
pub struct StructuralVariantData {
    pub id: i32,
    pub vcf_id_start: String,
    pub vcf_id_end: String,
    pub start_chromosome: String,
    pub end_chromosome: String,
    pub start_position: i32,
    pub end_position: i32,
    pub start_orientation: i8,
    pub end_orientation: i8,
    pub start_homology_sequence: String,
    pub end_homology_sequence: String,
    pub start_af: f64,
    pub end_af: f64,
    pub junction_copy_number: f64,
    pub adjusted_start_af: f64,
    pub adjusted_end_af: f64,
    pub adjusted_start_copy_number: f64,
    pub adjusted_end_copy_number: f64,
    pub adjusted_start_copy_number_change: f64,
    pub adjusted_end_copy_number_change: f64,
    pub insert_sequence: String,
    pub sv_type: Option<StructuralVariantType>,
    pub filter: String,
    pub quality_score: f64,
    pub event: String,
    pub start_tumor_variant_fragment_count: i32,
    pub start_tumor_reference_fragment_count: i32,
    pub start_normal_variant_fragment_count: i32,
    pub start_normal_reference_fragment_count: i32,
    pub end_tumor_variant_fragment_count: i32,
    pub end_tumor_reference_fragment_count: i32,
    pub end_normal_variant_fragment_count: i32,
    pub end_normal_reference_fragment_count: i32,
    pub start_interval_offset_start: i32,
    pub start_interval_offset_end: i32,
    pub end_interval_offset_start: i32,
    pub end_interval_offset_end: i32,
    pub inexact_homology_offset_start: i32,
    pub inexact_homology_offset_end: i32,
    pub start_linked_by: String,
    pub end_linked_by: String,
    pub insert_sequence_alignments: String,
    pub insert_sequence_repeat_class: String,
    pub insert_sequence_repeat_type: String,
    pub insert_sequence_repeat_orientation: i8,
    pub insert_sequence_repeat_coverage: f64,
    pub start_anchoring_support_distance: i32,
    pub end_anchoring_support_distance: i32,
    pub pon_count: i32,
}

impl StructuralVariantData {
    pub fn sv_type(&self) -> StructuralVariantType {
        self.sv_type.unwrap_or(StructuralVariantType::BND)
    }

    /// StructuralVariantData.convertSvData(EnrichedStructuralVariant, svId)
    pub fn convert(var: &StructuralVariant, sv_id: i32) -> StructuralVariantData {
        let s = &var.start;
        let e = var.end.as_ref();
        let nn = |v: Option<f64>| v.unwrap_or(0.0);
        let ni = |v: Option<i32>| v.unwrap_or(0);
        StructuralVariantData {
            id: sv_id,
            vcf_id_start: var.id.clone(),
            vcf_id_end: var.mate_id.clone().unwrap_or_default(),
            start_chromosome: s.chromosome.clone(),
            end_chromosome: e.map(|l| l.chromosome.clone()).unwrap_or_else(|| "0".to_string()),
            start_position: s.position,
            end_position: e.map(|l| l.position).unwrap_or(-1),
            start_orientation: s.orientation,
            end_orientation: e.map(|l| l.orientation).unwrap_or(0),
            start_homology_sequence: s.homology.clone(),
            end_homology_sequence: e.map(|l| l.homology.clone()).unwrap_or_default(),
            junction_copy_number: nn(var.junction_copy_number),
            start_af: nn(s.allele_frequency),
            end_af: e.map(|l| nn(l.allele_frequency)).unwrap_or(0.0),
            adjusted_start_af: nn(s.adjusted_allele_frequency),
            adjusted_end_af: e.map(|l| nn(l.adjusted_allele_frequency)).unwrap_or(0.0),
            adjusted_start_copy_number: nn(s.adjusted_copy_number),
            adjusted_end_copy_number: e.map(|l| nn(l.adjusted_copy_number)).unwrap_or(0.0),
            adjusted_start_copy_number_change: nn(s.adjusted_copy_number_change),
            adjusted_end_copy_number_change: e.map(|l| nn(l.adjusted_copy_number_change)).unwrap_or(0.0),
            insert_sequence: var.insert_sequence.clone(),
            sv_type: Some(var.sv_type),
            filter: var.filter.clone(),
            quality_score: var.quality_score,
            event: var.event.clone(),
            start_tumor_variant_fragment_count: ni(s.tumor_variant_fragment_count),
            start_tumor_reference_fragment_count: ni(s.tumor_reference_fragment_count),
            start_normal_variant_fragment_count: ni(s.normal_variant_fragment_count),
            start_normal_reference_fragment_count: ni(s.normal_reference_fragment_count),
            end_tumor_variant_fragment_count: e.map(|l| ni(l.tumor_variant_fragment_count)).unwrap_or(0),
            end_tumor_reference_fragment_count: e.map(|l| ni(l.tumor_reference_fragment_count)).unwrap_or(0),
            end_normal_variant_fragment_count: e.map(|l| ni(l.normal_variant_fragment_count)).unwrap_or(0),
            end_normal_reference_fragment_count: e.map(|l| ni(l.normal_reference_fragment_count)).unwrap_or(0),
            start_interval_offset_start: s.start_offset,
            start_interval_offset_end: s.end_offset,
            end_interval_offset_start: e.map(|l| l.start_offset).unwrap_or(0),
            end_interval_offset_end: e.map(|l| l.end_offset).unwrap_or(0),
            inexact_homology_offset_start: s.inexact_homology_offset_start,
            inexact_homology_offset_end: s.inexact_homology_offset_end,
            start_linked_by: var.start_linked_by.clone(),
            end_linked_by: var.end_linked_by.clone(),
            insert_sequence_alignments: var.insert_sequence_alignments.clone(),
            insert_sequence_repeat_class: var.insert_sequence_repeat_class.clone().unwrap_or_default(),
            insert_sequence_repeat_type: var.insert_sequence_repeat_type.clone().unwrap_or_default(),
            insert_sequence_repeat_orientation: var.insert_sequence_repeat_orientation.unwrap_or(0),
            insert_sequence_repeat_coverage: var.insert_sequence_repeat_coverage.unwrap_or(0.0),
            start_anchoring_support_distance: s.anchoring_support_distance,
            end_anchoring_support_distance: e.map(|l| l.anchoring_support_distance).unwrap_or(0),
            pon_count: 0,
        }
    }
}

/// SvUtils.formSvType
pub fn form_sv_type(
    chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: i8, orient_end: i8, has_inserted_bases: bool,
) -> StructuralVariantType {
    use StructuralVariantType::*;
    if chr_start != chr_end {
        return BND;
    }
    if orient_start != orient_end {
        let pos_diff = (pos_start - pos_end).abs();
        if pos_diff == 1 && has_inserted_bases {
            return INS;
        }
        if pos_diff == 0 {
            return DUP;
        }
        let first_is_lower = pos_start < pos_end;
        if first_is_lower == (orient_start == 1) { DEL } else { DUP }
    } else {
        INV
    }
}

pub fn is_indel(t: StructuralVariantType) -> bool {
    matches!(t, StructuralVariantType::DEL | StructuralVariantType::DUP | StructuralVariantType::INS)
}

/// SvUtils.SHORT_INDEL_LENGTH
pub const SHORT_INDEL_LENGTH: i32 = 1000;

/// StructuralVariantFactory.trimStringListValue
pub fn trim_string_list_value(v: &str) -> String {
    if v.is_empty() {
        return String::new();
    }
    v.replace(['[', ']', ' '], "")
}

/// The breakend notation of an ALT: BREAKEND_REGEX `^(.*)([\[\]])(.+)[\[\]](.*)$` (greedy), as its four groups
pub fn match_breakend_alt(alt: &str) -> Option<(String, String, String, String)> {
    // group 4 is the text after the last bracket; group 2 the bracket before group 3, where group 1 is greedy:
    // the second-to-last bracket overall, with group 3 non-empty
    let chars: Vec<char> = alt.chars().collect();
    let brackets: Vec<usize> = chars.iter().enumerate().filter(|(_, c)| **c == '[' || **c == ']').map(|(i, _)| i).collect();
    if brackets.len() < 2 {
        return None;
    }
    let last = *brackets.last().unwrap();
    // group 1 greedy: the latest bracket before `last` such that group 3 (between) is non-empty
    for &b in brackets.iter().rev().skip(1) {
        if last > b + 1 {
            let s = |a: usize, z: usize| chars[a..z].iter().collect::<String>();
            return Some((s(0, b), chars[b].to_string(), s(b + 1, last), s(last + 1, chars.len())));
        }
    }
    None
}

/// The record filter for the VCF reader (VariantContextFilter)
pub type RecordFilter = fn(&VcfRecord) -> bool;

/// StructuralVariantFactory: pairs breakends by MATEID in file order
pub struct StructuralVariantFactory {
    filter: RecordFilter,
    unmatched: HashMap<String, VcfRecord>,
    complete: Vec<StructuralVariant>,
    ordinals_set: bool,
    reference_ordinal: i32,
    tumor_ordinal: i32,
    /// EnrichedStructuralVariantFactory: the PURPLE fields are read for every variant
    enrich: bool,
}

impl StructuralVariantFactory {
    /// StructuralVariantFactory.build(filter): human chromosomes only, then the filter
    pub fn build(filter: RecordFilter) -> Self {
        StructuralVariantFactory {
            filter,
            unmatched: HashMap::new(),
            complete: Vec::new(),
            ordinals_set: false,
            reference_ordinal: -1,
            tumor_ordinal: -1,
            enrich: false,
        }
    }

    /// with EnrichedStructuralVariantFactory.enrich applied to the results
    pub fn set_enrich(&mut self, enrich: bool) {
        self.enrich = enrich;
    }

    pub fn set_genotype_ordinals(&mut self, reference_ordinal: i32, tumor_ordinal: i32) {
        self.ordinals_set = true;
        self.reference_ordinal = reference_ordinal;
        self.tumor_ordinal = tumor_ordinal;
    }

    pub fn results(self) -> Vec<StructuralVariant> {
        self.complete
    }

    pub fn add_variant_context(&mut self, context: VcfRecord) -> Result<()> {
        if !is_human_chromosome(&context.chrom) || !(self.filter)(&context) {
            return Ok(());
        }
        if is_single_breakend(&context.alt) {
            let sv = self.create_single_breakend(context)?;
            self.complete.push(sv);
        } else {
            let mate = context.attribute_as_string_opt(MATE_ID);
            if let Some(start) = mate.as_ref().and_then(|m| self.unmatched.remove(m)) {
                let sv = self.create_sv(start, context)?;
                self.complete.push(sv);
            } else {
                self.unmatched.insert(context.id.clone(), context);
            }
        }
        Ok(())
    }

    fn create_sv(&self, context_start: VcfRecord, context_end: VcfRecord) -> Result<StructuralVariant> {
        let alt = &context_start.alt;
        let Some((g1, g2, _g3, g4)) = match_breakend_alt(alt) else {
            bail!("ALT {alt} is not in breakend notation");
        };
        let mut start_orientation: i8 = if !g1.is_empty() { 1 } else { -1 };
        let mut end_orientation: i8 = if g2 == "]" { 1 } else { -1 };
        let inserted_sequence = if !g1.is_empty() { g1[1..].to_string() } else { g4[..g4.len() - 1].to_string() };
        let mut sv_type = if context_start.has_attribute(SV_TYPE) {
            let t = context_start.attribute_as_string(SV_TYPE, "");
            Some(StructuralVariantType::from_name(&t).ok_or_else(|| anyhow!("No enum constant StructuralVariantType.{t}"))?)
        } else {
            None
        };
        let (mut leg_start_ctx, mut leg_end_ctx) = (&context_start, &context_end);
        if context_start.chrom == context_end.chrom
            && start_orientation != end_orientation
            && context_start.pos == context_end.pos
            && start_orientation == 1
        {
            leg_start_ctx = &context_end;
            leg_end_ctx = &context_start;
            start_orientation = -1;
            end_orientation = 1;
        }
        if sv_type.is_none() {
            sv_type = Some(form_sv_type(
                &leg_start_ctx.chrom,
                &leg_end_ctx.chrom,
                leg_start_ctx.pos,
                leg_end_ctx.pos,
                start_orientation,
                end_orientation,
                !inserted_sequence.is_empty(),
            ));
        }
        let sv_type = sv_type.unwrap();
        let indel_length = if is_indel(sv_type) { (context_start.pos - context_end.pos).abs() } else { 0 };
        let is_small_del_dup = is_indel(sv_type) && indel_length <= SHORT_INDEL_LENGTH;
        let mut start_leg = self.set_leg_common(leg_start_ctx, is_small_del_dup, start_orientation)?;
        start_leg.homology = leg_start_ctx.attribute_as_string(HOMSEQ, "");
        let mut end_leg = self.set_leg_common(leg_end_ctx, is_small_del_dup, end_orientation)?;
        end_leg.homology = leg_end_ctx.attribute_as_string(HOMSEQ, "");
        let mut sv = self.set_common(&context_start, start_leg)?;
        sv.end = Some(end_leg);
        sv.mate_id = Some(leg_end_ctx.id.clone());
        sv.insert_sequence = inserted_sequence;
        sv.sv_type = sv_type;
        sv.filter = filters(leg_start_ctx, Some(leg_end_ctx));
        sv.start_linked_by = trim_string_list_value(&leg_start_ctx.attribute_as_string(ASM_LINKS, ""));
        sv.end_linked_by = trim_string_list_value(&leg_end_ctx.attribute_as_string(ASM_LINKS, ""));
        let (s, e) = (leg_start_ctx.clone(), leg_end_ctx.clone());
        sv.start_context = s;
        sv.end_context = Some(e);
        if self.enrich {
            enrich(&mut sv)?;
        }
        Ok(sv)
    }

    fn create_single_breakend(&self, context: VcfRecord) -> Result<StructuralVariant> {
        let alt = &context.alt;
        let orientation: i8 = if alt.starts_with('.') { -1 } else { 1 };
        let ref_length = context.ref_bases.len();
        let inserted_sequence =
            if orientation == -1 { alt[1..alt.len() - ref_length].to_string() } else { alt[ref_length..alt.len() - 1].to_string() };
        let mut start_leg = self.set_leg_common(&context, false, orientation)?;
        start_leg.homology = String::new();
        let mut sv = self.set_common(&context, start_leg)?;
        sv.insert_sequence = inserted_sequence;
        sv.sv_type = if context.has_attribute(INFERRED) { StructuralVariantType::INF } else { StructuralVariantType::SGL };
        sv.filter = filters(&context, None);
        sv.start_linked_by = trim_string_list_value(&context.attribute_as_string(ASM_LINKS, ""));
        sv.start_context = context;
        if self.enrich {
            enrich(&mut sv)?;
        }
        Ok(sv)
    }

    fn set_common(&self, context: &VcfRecord, start: StructuralVariantLeg) -> Result<StructuralVariant> {
        let mut sv = StructuralVariant {
            id: context.id.clone(),
            mate_id: None,
            start,
            end: None,
            insert_sequence: String::new(),
            sv_type: StructuralVariantType::SGL,
            filter: String::new(),
            hotspot: context.attribute_as_bool(HOTSPOT, false),
            quality_score: context.qual,
            event: String::new(),
            start_linked_by: String::new(),
            end_linked_by: String::new(),
            insert_sequence_alignments: trim_string_list_value(&context.attribute_as_string(INSALN, "")),
            insert_sequence_repeat_class: None,
            insert_sequence_repeat_type: None,
            insert_sequence_repeat_orientation: None,
            insert_sequence_repeat_coverage: None,
            junction_copy_number: None,
            start_context: context.clone(),
            end_context: None,
        };
        if context.has_attribute(REPEAT_MASK_REPEAT_CLASS) {
            sv.insert_sequence_repeat_class = Some(context.attribute_as_string(REPEAT_MASK_REPEAT_CLASS, ""));
            sv.insert_sequence_repeat_type = Some(context.attribute_as_string(REPEAT_MASK_REPEAT_TYPE, ""));
            sv.insert_sequence_repeat_orientation =
                Some(if context.attribute_as_string(REPEAT_MASK_ORIENTATION, "+") == "+" { 1 } else { -1 });
            sv.insert_sequence_repeat_coverage = Some(context.attribute_as_double(REPEAT_MASK_COVERAGE, 0.0)?);
        }
        Ok(sv)
    }

    fn set_leg_common(&self, context: &VcfRecord, ignore_refpair: bool, orientation: i8) -> Result<StructuralVariantLeg> {
        let mut leg = StructuralVariantLeg {
            chromosome: context.chrom.clone(),
            position: context.pos,
            orientation,
            ..Default::default()
        };
        if context.has_attribute(CIPOS) {
            let cipos = context.attribute_as_int_list(CIPOS, 0)?;
            if cipos.len() == 2 {
                leg.start_offset = cipos[0];
                leg.end_offset = cipos[1];
            }
        }
        let seg_lengths = context.attribute_as_int_list(SEG_ALIGN_LENGTH, 0)?;
        leg.anchoring_support_distance = seg_lengths.iter().copied().max().unwrap_or(0);
        if context.has_attribute(IHOMPOS) {
            let ihompos = context.attribute_as_int_list(IHOMPOS, 0)?;
            if ihompos.len() == 2 {
                leg.inexact_homology_offset_start = ihompos[0];
                leg.inexact_homology_offset_end = ihompos[1];
            }
        }
        let (reference_ordinal, tumor_ordinal) = if self.ordinals_set {
            (self.reference_ordinal, self.tumor_ordinal)
        } else if context.genotype_count() == 1 {
            (-1, 0)
        } else {
            (0, 1)
        };
        if reference_ordinal >= 0 && (reference_ordinal as usize) < context.genotype_count() {
            let o = reference_ordinal as usize;
            let total_frags = context.genotype_attribute_as_int(o, TOTAL_FRAGS, 0)?;
            let ref_frags = context.genotype_attribute_as_int(o, REF_DEPTH, 0)?;
            let ref_pair_frags = context.genotype_attribute_as_int(o, REF_DEPTH_PAIR, 0)?;
            leg.normal_variant_fragment_count = Some(total_frags);
            leg.normal_reference_fragment_count = Some(ref_frags + if ignore_refpair { 0 } else { ref_pair_frags });
        }
        if tumor_ordinal >= 0 && (tumor_ordinal as usize) < context.genotype_count() {
            let o = tumor_ordinal as usize;
            let total_frags = context.genotype_attribute_as_int(o, TOTAL_FRAGS, 0)?;
            let ref_frags = context.genotype_attribute_as_int(o, REF_DEPTH, 0)?;
            let ref_pair_frags = context.genotype_attribute_as_int(o, REF_DEPTH_PAIR, 0)?;
            leg.tumor_variant_fragment_count = Some(total_frags);
            leg.tumor_reference_fragment_count = Some(ref_frags + if ignore_refpair { 0 } else { ref_pair_frags });
            leg.allele_frequency = Some(context.genotype_attribute_as_double(o, ALLELE_FRACTION, 0.0)?);
        }
        Ok(leg)
    }
}

/// StructuralVariantFactory.isSingleBreakend
pub fn is_single_breakend(alt: &str) -> bool {
    let b = alt.as_bytes();
    !b.is_empty() && (b[0] == b'.' || b[b.len() - 1] == b'.')
}

/// StructuralVariantFactory.filters: the union of both records' filters, PASS dropped if others are
/// present, PASS if none, sorted and joined by ';'
fn filters(context: &VcfRecord, paired: Option<&VcfRecord>) -> String {
    let mut set: Vec<String> = Vec::new();
    for f in context.filters.iter().chain(paired.map(|p| p.filters.iter()).into_iter().flatten()) {
        if !set.contains(f) {
            set.push(f.clone());
        }
    }
    if set.len() > 1 {
        set.retain(|f| f != PASS_FILTER);
    }
    if set.is_empty() {
        set.push(PASS_FILTER.to_string());
    }
    set.sort_by(|a, b| crate::java_text::compare(a, b));
    set.join(";")
}

/// EnrichedStructuralVariantFactory.enrich for one variant
fn enrich(sv: &mut StructuralVariant) -> Result<()> {
    let ctx = &sv.start_context;
    sv.junction_copy_number = if ctx.has_attribute(PURPLE_JUNCTION_COPY_NUMBER) {
        Some(ctx.attribute_as_double(PURPLE_JUNCTION_COPY_NUMBER, 0.0)?)
    } else if ctx.has_attribute(PURPLE_PLOIDY_INFO) {
        Some(ctx.attribute_as_double(PURPLE_PLOIDY_INFO, 0.0)?)
    } else {
        None
    };
    enrich_leg(&mut sv.start, &sv.start_context)?;
    if let (Some(end), Some(ec)) = (sv.end.as_mut(), sv.end_context.as_ref()) {
        enrich_leg(end, ec)?;
    }
    Ok(())
}

fn enrich_leg(leg: &mut StructuralVariantLeg, ctx: &VcfRecord) -> Result<()> {
    let first = |key: &str| -> Result<Option<f64>> {
        if !ctx.has_attribute(key) {
            return Ok(None);
        }
        Ok(ctx.attribute_as_double_list(key, 0.0)?.first().copied())
    };
    leg.adjusted_allele_frequency = first(PURPLE_AF)?;
    leg.adjusted_copy_number = first(PURPLE_CN)?;
    leg.adjusted_copy_number_change = first(PURPLE_CN_CHANGE)?;
    Ok(())
}

/// AlwaysPassFilter
pub fn always_pass(_: &VcfRecord) -> bool {
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breakend_alts() {
        assert_eq!(match_breakend_alt("G[chr1:5023347["), Some(("G".into(), "[".into(), "chr1:5023347".into(), "".into())));
        assert_eq!(match_breakend_alt("]chr1:5023299]A"), Some(("".into(), "]".into(), "chr1:5023299".into(), "A".into())));
        assert_eq!(match_breakend_alt("GTT]chr5:100]"), Some(("GTT".into(), "]".into(), "chr5:100".into(), "".into())));
        assert_eq!(match_breakend_alt("G."), None);
        assert!(is_single_breakend("G.") && is_single_breakend(".TTG") && !is_single_breakend("G[1:2["));
    }

    #[test]
    fn types_and_filters() {
        use StructuralVariantType::*;
        assert_eq!(form_sv_type("1", "2", 1, 2, 1, 1, false), BND);
        assert_eq!(form_sv_type("1", "1", 100, 200, 1, -1, false), DEL);
        assert_eq!(form_sv_type("1", "1", 100, 200, -1, 1, false), DUP);
        assert_eq!(form_sv_type("1", "1", 100, 101, 1, -1, true), INS);
        assert_eq!(form_sv_type("1", "1", 100, 200, 1, 1, false), INV);
        let r = VcfRecord::parse("chr1\t10\t1\tA\tA.\t.\tPASS\tINSALN=a,b;LINE\tGT:AF\t./.:0.5\t./.:.").unwrap();
        assert_eq!(r.qual, -10.0);
        assert!(r.filters.is_empty());
        assert_eq!(r.attribute_as_string(INSALN, ""), "[a, b]");
        assert_eq!(trim_string_list_value(&r.attribute_as_string(INSALN, "")), "a,b");
        assert!(r.has_attribute("LINE"));
        assert_eq!(r.genotype_attribute_as_double(0, "AF", 0.0).unwrap(), 0.5);
        assert_eq!(r.genotype_attribute_as_double(1, "AF", 0.0).unwrap(), 0.0);
        let p = VcfRecord::parse("chr1\t10\t2\tA\tA.\t5\tPON;LOW\t.").unwrap();
        assert_eq!(filters(&r, Some(&p)), "LOW;PON");
        assert_eq!(filters(&r, None), "PASS");
        assert!(is_human_chromosome("chr22") && is_human_chromosome("X") && !is_human_chromosome("chr23") && !is_human_chromosome("chrM"));
    }
}
