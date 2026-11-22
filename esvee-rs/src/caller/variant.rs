//! Variants and their breakends (Variant, Breakend, SvDataCache and the parts of hmf-common's
//! StructuralVariantFactory that the caller relies on).
//!
//! Variants live in one list; a breakend is addressed by the index of its variant and whether it is
//! the start or the end.

use super::types::*;
use super::vcf::{Genotype, VcfRecord};
use hashbrown::HashMap;

pub const SE_START: usize = 0;
pub const SE_END: usize = 1;

// SvVcfTags
pub const MATE_ID: &str = "MATEID";
pub const SV_TYPE: &str = "SVTYPE";
pub const CIPOS: &str = "CIPOS";
pub const HOMSEQ: &str = "HOMSEQ";
pub const IHOMPOS: &str = "IHOMPOS";
pub const INSALN: &str = "INSALN";
pub const ASM_INFO: &str = "ASMINFO";
pub const ASM_LENGTH: &str = "ASMLEN";
pub const SEG_REPEAT_LENGTH: &str = "SEGRL";
pub const ASM_LINKS: &str = "ASMLNKS";
pub const TOTAL_FRAGS: &str = "VF";
pub const AVG_FRAG_LENGTH: &str = "AVGLEN";
pub const LINE_SITE: &str = "LINE";
pub const UNIQUE_FRAG_POSITIONS: &str = "UFP";
pub const THREE_PRIME_RANGE: &str = "TPR";
pub const MAX_LOCAL_REPEAT: &str = "MLR";
pub const PROX_JUNC_READ_RATIO: &str = "PJRR";
pub const SUPP_REMOTE_REGION_RATIO: &str = "SRRR";
pub const SAGA_VARIANT: &str = "SAGA";
pub const SPLIT_FRAGS: &str = "SF";
pub const STRAND_BIAS: &str = "SB";
pub const REF_DEPTH: &str = "REF";
pub const REF_DEPTH_PAIR: &str = "REFPAIR";
pub const ALLELE_FRACTION: &str = "AF";

/// SvUtils.SHORT_INDEL_LENGTH
pub const SHORT_INDEL_LENGTH: i32 = 1000;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Interval {
    pub start: i32,
    pub end: i32,
}

impl Interval {
    pub fn length(&self) -> i32 {
        (self.end - self.start).abs()
    }
}

/// The index of a variant and the breakend within it
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BreakendRef {
    pub var: u32,
    pub se: u8,
}

impl BreakendRef {
    pub fn new(var: usize, se: usize) -> BreakendRef {
        BreakendRef { var: var as u32, se: se as u8 }
    }
    pub fn other(self) -> BreakendRef {
        BreakendRef { var: self.var, se: 1 - self.se }
    }
}

#[derive(Clone, Debug)]
pub struct Breakend {
    pub context: VcfRecord,
    pub chromosome: String,
    pub position: i32,
    pub orient: Orient,
    pub is_start: bool,
    pub insert_sequence: String,
    pub confidence_interval: Interval,
    pub inexact_homology: Interval,
    /// index of the reference and tumor genotypes in the record, if the sample is in the run
    pub ref_genotype: Option<usize>,
    pub tumor_genotype: Option<usize>,
    pub is_line_insertion: bool,
    pub line_site_breakend: Option<BreakendRef>,
    pub line_remote_source_site: bool,
}

/// VariantAltInsertCoords.fromRefAlt, the insert sequence only
fn insert_sequence_from_alt(alt: &str, ref_allele: &str) -> String {
    let ref_len = ref_allele.len();
    if alt.starts_with('.') {
        return java_substring(alt, ref_len, alt.len().wrapping_sub(1)).to_string();
    } else if alt.ends_with('.') {
        return java_substring(alt, 1, alt.len().wrapping_sub(ref_len)).to_string();
    }
    if let Some((before, _, _, after)) = match_breakend(alt) {
        if !before.is_empty() {
            return java_substring(before, ref_len, before.len()).to_string();
        } else if !after.is_empty() {
            return java_substring(after, 0, after.len().wrapping_sub(ref_len)).to_string();
        }
    }
    String::new()
}

/// String.substring, with its bounds failure
fn java_substring(s: &str, begin: usize, end: usize) -> &str {
    if begin > end || end > s.len() {
        panic!("StringIndexOutOfBoundsException: begin {begin}, end {end}, length {}", s.len());
    }
    &s[begin..end]
}

/// BREAKEND_REGEX `^(.*)([\[\]])(.+)[\[\]](.*)$` with greedy groups: (before, bracket, mate, after)
pub fn match_breakend(alt: &str) -> Option<(&str, char, &str, &str)> {
    let is_bracket = |c: u8| c == b'[' || c == b']';
    let bytes = alt.as_bytes();
    let last = bytes.iter().rposition(|&c| is_bracket(c))?;
    // group 1 is greedy: the first bracket is the last one that still leaves a non-empty group 3
    let first = bytes[..last].iter().enumerate().rev().find(|&(i, &c)| is_bracket(c) && i + 1 < last).map(|(i, _)| i)?;
    Some((&alt[..first], bytes[first] as char, &alt[first + 1..last], &alt[last + 1..]))
}

/// StructuralVariantFactory.isSingleBreakend
pub fn is_single_breakend(rec: &VcfRecord) -> bool {
    let alt = rec.alt().as_bytes();
    !alt.is_empty() && (alt[0] == b'.' || alt[alt.len() - 1] == b'.')
}

impl Breakend {
    fn new(is_start: bool, context: VcfRecord, position: i32, orient: Orient, reference_ordinal: i32, tumor_ordinal: i32) -> Breakend {
        let cipos = context.attr_int_list(CIPOS, 0);
        let confidence_interval = if cipos.len() == 2 { Interval { start: cipos[0], end: cipos[1] } } else { Interval::default() };
        let insert_sequence = insert_sequence_from_alt(context.alt(), &context.ref_allele);
        let inexact_homology = if context.has_attr(IHOMPOS) {
            let ihompos = context.attr_int_list(IHOMPOS, 0);
            if ihompos.len() < 2 {
                panic!("IndexOutOfBoundsException: IHOMPOS with {} value(s)", ihompos.len());
            }
            Interval { start: ihompos[0], end: ihompos[1] }
        } else {
            Interval::default()
        };
        // VariantContext.getGenotype(int) is null beyond the genotypes of the record
        let genotype_index = |ordinal: i32| if ordinal >= 0 && (ordinal as usize) < context.genotypes.len() { Some(ordinal as usize) } else { None };
        Breakend {
            chromosome: context.chrom.clone(),
            position,
            orient,
            is_start,
            insert_sequence,
            confidence_interval,
            inexact_homology,
            ref_genotype: genotype_index(reference_ordinal),
            tumor_genotype: genotype_index(tumor_ordinal),
            is_line_insertion: context.has_attr(LINE_SITE),
            line_site_breakend: None,
            line_remote_source_site: false,
            context,
        }
    }

    pub fn vcf_id(&self) -> &str {
        &self.context.id
    }
    pub fn is_line(&self) -> bool {
        self.is_line_insertion
    }
    pub fn in_chained_assembly(&self) -> bool {
        self.context.has_attr(ASM_LINKS)
    }
    pub fn anchor_length(&self) -> i32 {
        self.context.attr_int(SEG_REPEAT_LENGTH, 0)
    }
    pub fn unique_fragment_positions(&self) -> i32 {
        self.context.attr_int(UNIQUE_FRAG_POSITIONS, 0)
    }
    fn tumor(&self) -> Option<&Genotype> {
        self.tumor_genotype.map(|i| &self.context.genotypes[i])
    }
    fn reference(&self) -> Option<&Genotype> {
        self.ref_genotype.map(|i| &self.context.genotypes[i])
    }

    /// Breakend.fragmentCount(genotype, tag)
    pub fn fragment_count_tag(genotype: Option<&Genotype>, tag: &str) -> i32 {
        genotype.map(|g| g.attr_int(tag, 0)).unwrap_or(0)
    }
    pub fn fragment_count_of(genotype: &Genotype) -> i32 {
        genotype.attr_int(TOTAL_FRAGS, 0)
    }
    /// Breakend.fragmentCount(): tumor plus reference
    pub fn fragment_count(&self) -> i32 {
        Self::fragment_count_tag(self.tumor(), TOTAL_FRAGS) + Self::fragment_count_tag(self.reference(), TOTAL_FRAGS)
    }

    pub fn strand_bias(&self) -> f64 {
        let genotype = self.tumor().or(self.reference());
        genotype.map(|g| g.attr_double(STRAND_BIAS, 0.5)).unwrap_or(0.5)
    }

    pub fn min_orientation_position_range(&self) -> i32 {
        let ranges = self.context.attr_int_list(THREE_PRIME_RANGE, 0);
        if ranges.len() != 2 {
            return -1;
        }
        let (five, three) = (ranges[0], ranges[1]);
        if five >= 0 && three >= 0 {
            return five.min(three);
        }
        if five >= 0 { five } else { three }
    }
}

/// An original junction of the assembly, from the ASMINFO tag
#[derive(Clone, Debug, PartialEq)]
pub struct Junction {
    pub chromosome: String,
    pub position: i32,
    pub orient: Orient,
    pub discordant_only: bool,
    pub indel_based: bool,
}

impl Junction {
    pub fn soft_clip_based(&self) -> bool {
        !self.indel_based && !self.discordant_only
    }
}

#[derive(Clone, Debug)]
pub struct RepeatMaskAnnotation {
    pub class_type: String,
    pub repeat: String,
    pub coverage: f64,
}

#[derive(Clone, Debug)]
pub struct Variant {
    pub sv_type: SvType,
    pub breakends: [Option<Breakend>; 2],
    pub qual: f64,
    pub insert_sequence: String,
    pub is_short_local_indel: bool,
    pub pon_count: i32,
    pub is_hotspot: bool,
    pub germline: bool,
    pub rm_annotation: Option<RepeatMaskAnnotation>,
    pub filters: FilterSet,
    pub is_saga_matched: bool,
    pub original_junctions: Vec<Junction>,
}

/// Variant.originalJunctions
fn parse_original_junctions(context: &VcfRecord) -> Vec<Junction> {
    let info = context.attr_string(ASM_INFO, "");
    if info.is_empty() {
        return Vec::new();
    }
    let mut junctions = Vec::new();
    for assembly in info.split('_') {
        // Java splits without a limit, so trailing empty items are dropped
        let mut items: Vec<&str> = assembly.split(':').collect();
        while items.last() == Some(&"") {
            items.pop();
        }
        if items.len() < 3 {
            continue;
        }
        let position =
            super::vcf::java_parse_int(items[1]).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{}\"", items[1]));
        let orient = Orient::from_byte_str(items[2]);
        let (mut indel_based, mut discordant_only) = (false, false);
        if items.len() >= 4 {
            indel_based = items[3] == "I";
            discordant_only = items[3] == "D";
        }
        junctions.push(Junction { chromosome: items[0].to_string(), position, orient, discordant_only, indel_based });
    }
    junctions
}

impl Variant {
    /// StructuralVariantFactory.createSingleBreakend followed by the Variant constructor
    pub fn new_sgl(context: VcfRecord, reference_ordinal: i32, tumor_ordinal: i32) -> Variant {
        let alt = context.alt().to_string();
        let orient = if alt.starts_with('.') { Orient::Reverse } else { Orient::Forward };
        let ref_length = context.ref_allele.len();
        let insert_sequence = if orient == Orient::Reverse {
            java_substring(&alt, 1, alt.len().wrapping_sub(ref_length)).to_string()
        } else {
            java_substring(&alt, ref_length, alt.len() - 1).to_string()
        };
        let sv_type = if context.has_attr("INFERRED") { SvType::Inf } else { SvType::Sgl };
        let qual = context.phred_scaled_qual();
        let is_saga_matched = context.has_attr(SAGA_VARIANT);
        let position = context.pos;
        let start = Breakend::new(true, context, position, orient, reference_ordinal, tumor_ordinal);
        let original_junctions = parse_original_junctions(&start.context);
        Variant {
            sv_type,
            breakends: [Some(start), None],
            qual,
            insert_sequence,
            is_short_local_indel: false,
            pon_count: 0,
            is_hotspot: false,
            germline: false,
            rm_annotation: None,
            filters: FilterSet::default(),
            is_saga_matched,
            original_junctions,
        }
    }

    /// StructuralVariantFactory.createSV followed by the Variant constructor. `context_start` is the
    /// breakend that was read first.
    pub fn new_sv(context_start: VcfRecord, context_end: VcfRecord, reference_ordinal: i32, tumor_ordinal: i32) -> Variant {
        let alt = context_start.alt().to_string();
        let Some((before, bracket, _, after)) = match_breakend(&alt) else {
            panic!("IllegalArgumentException: ALT {alt} is not in breakend notation");
        };
        let mut start_orient = if !before.is_empty() { Orient::Forward } else { Orient::Reverse };
        let mut end_orient = if bracket == ']' { Orient::Forward } else { Orient::Reverse };
        let insert_sequence =
            if !before.is_empty() { java_substring(before, 1, before.len()).to_string() } else { java_substring(after, 0, after.len().wrapping_sub(1)).to_string() };

        let type_from_context = if context_start.has_attr(SV_TYPE) { Some(SvType::value_of(&context_start.attr_string(SV_TYPE, ""))) } else { None };
        let qual = context_start.phred_scaled_qual();

        // a DUP of one base has both breakends at the same position: the lower one is the reverse one
        let (leg_start, leg_end) = if context_start.chrom == context_end.chrom
            && start_orient != end_orient
            && context_start.pos == context_end.pos
            && start_orient == Orient::Forward
        {
            start_orient = Orient::Reverse;
            end_orient = Orient::Forward;
            (context_end, context_start)
        } else {
            (context_start, context_end)
        };

        let sv_type = type_from_context.unwrap_or_else(|| {
            form_sv_type(&leg_start.chrom, &leg_end.chrom, leg_start.pos, leg_end.pos, start_orient, end_orient, !insert_sequence.is_empty())
        });

        let is_saga_matched = leg_start.has_attr(SAGA_VARIANT) || leg_end.has_attr(SAGA_VARIANT);
        let (pos_start, pos_end) = (leg_start.pos, leg_end.pos);
        let is_short_local_indel = sv_type.is_indel() && (pos_end - pos_start).abs() <= SHORT_INDEL_LENGTH;

        let start = Breakend::new(true, leg_start, pos_start, start_orient, reference_ordinal, tumor_ordinal);
        let end = Breakend::new(false, leg_end, pos_end, end_orient, reference_ordinal, tumor_ordinal);
        let original_junctions = parse_original_junctions(&start.context);

        Variant {
            sv_type,
            breakends: [Some(start), Some(end)],
            qual,
            insert_sequence,
            is_short_local_indel,
            pon_count: 0,
            is_hotspot: false,
            germline: false,
            rm_annotation: None,
            filters: FilterSet::default(),
            is_saga_matched,
            original_junctions,
        }
    }

    pub fn is_sgl(&self) -> bool {
        self.sv_type == SvType::Sgl
    }
    pub fn start(&self) -> &Breakend {
        self.breakends[SE_START].as_ref().unwrap()
    }
    pub fn end(&self) -> Option<&Breakend> {
        self.breakends[SE_END].as_ref()
    }
    pub fn breakend(&self, se: usize) -> &Breakend {
        self.breakends[se].as_ref().unwrap()
    }
    pub fn breakend_mut(&mut self, se: usize) -> &mut Breakend {
        self.breakends[se].as_mut().unwrap()
    }
    pub fn chromosome_start(&self) -> &str {
        &self.start().chromosome
    }
    pub fn chromosome_end(&self) -> &str {
        if !self.is_sgl() { &self.breakend(SE_END).chromosome } else { "" }
    }
    pub fn pos_start(&self) -> i32 {
        self.start().position
    }
    pub fn pos_end(&self) -> i32 {
        if !self.is_sgl() { self.breakend(SE_END).position } else { -1 }
    }
    pub fn orient_start(&self) -> Orient {
        self.start().orient
    }
    pub fn orient_end(&self) -> Orient {
        self.breakend(SE_END).orient
    }
    pub fn is_short_local(&self) -> bool {
        self.is_short_local_indel
    }
    pub fn sv_length(&self) -> i32 {
        if self.sv_type.has_length() { (self.pos_end() - self.pos_start()).abs() } else { 0 }
    }
    pub fn adjusted_length(&self) -> i32 {
        if !self.sv_type.has_length() {
            return 0;
        }
        if self.sv_type == SvType::Ins {
            return self.insert_sequence.len() as i32;
        }
        let mut position_length = self.pos_end() - self.pos_start();
        if self.sv_type == SvType::Dup {
            position_length += 1;
        }
        position_length + self.insert_sequence.len() as i32
    }
    pub fn average_fragment_length(&self) -> i32 {
        self.start().context.attr_int(AVG_FRAG_LENGTH, 0)
    }
    pub fn split_fragment_count(&self) -> i32 {
        self.start().context.attr_int(TOTAL_FRAGS, 0)
    }
    pub fn is_line_site(&self) -> bool {
        self.breakends.iter().flatten().any(|b| b.is_line() || b.line_site_breakend.is_some())
    }
    pub fn in_chained_assembly(&self) -> bool {
        self.breakends.iter().flatten().any(|b| b.in_chained_assembly())
    }
    pub fn max_unique_fragment_positions(&self) -> i32 {
        self.breakends.iter().flatten().map(|b| b.unique_fragment_positions()).fold(0, i32::max)
    }
    pub fn add_filter(&mut self, filter: FilterType) {
        self.filters.add(filter);
    }
    pub fn is_pass(&self) -> bool {
        self.filters.is_empty()
    }
    pub fn is_filtered(&self) -> bool {
        !self.filters.is_empty()
    }

    /// Breakend.calcAllelicFrequency(genotype)
    pub fn calc_allelic_frequency(&self, genotype: &Genotype) -> f64 {
        if genotype.has_attr(ALLELE_FRACTION) {
            return genotype.attr_double(ALLELE_FRACTION, 0.0);
        }
        let read_pair_support = if self.is_sgl() || !self.is_short_local() { genotype.attr_int(REF_DEPTH_PAIR, 0) } else { 0 };
        let ref_support = genotype.attr_int(REF_DEPTH, 0);
        let fragment_count = Breakend::fragment_count_of(genotype);
        let total_support = (fragment_count + ref_support + read_pair_support) as f64;
        if total_support > 0.0 { fragment_count as f64 / total_support } else { 0.0 }
    }

    /// Breakend.calcAllelicFrequency(): of the tumor if there is one, else of the reference
    pub fn breakend_allelic_frequency(&self, se: usize) -> f64 {
        let breakend = self.breakend(se);
        match breakend.tumor().or(breakend.reference()) {
            Some(g) => self.calc_allelic_frequency(g),
            None => panic!("NullPointerException: breakend without a genotype"),
        }
    }

    /// Breakend.closeToOriginalJunction
    pub fn close_to_original_junction(&self, se: usize, permitted_distance: i32) -> bool {
        let breakend = self.breakend(se);
        self.original_junctions
            .iter()
            .any(|j| breakend.chromosome == j.chromosome && breakend.orient == j.orient && (breakend.position - j.position).abs() <= permitted_distance)
    }
}

/// Breakend.toString
pub fn describe_breakend(variants: &[Variant], b: BreakendRef) -> String {
    let var = &variants[b.var as usize];
    let breakend = var.breakend(b.se as usize);
    format!("{}:{} pos({}:{}:{})", breakend.vcf_id(), var.sv_type.name(), breakend.chromosome, breakend.position, breakend.orient.as_byte())
}

/// The breakends of each chromosome in position order. `order` gives the chromosomes in the order a
/// Java HashMap iterates over them, which the steps that walk the whole map follow.
#[derive(Default)]
pub struct BreakendMap {
    pub chromosomes: Vec<String>,
    pub lists: Vec<Vec<BreakendRef>>,
    pub order: Vec<usize>,
}

impl BreakendMap {
    /// SvDataCache.buildBreakendMap
    pub fn build(variants: &[Variant]) -> BreakendMap {
        let mut index: HashMap<String, usize> = HashMap::new();
        let mut map = BreakendMap::default();
        for (v, var) in variants.iter().enumerate() {
            for se in SE_START..=SE_END {
                let Some(breakend) = &var.breakends[se] else { continue };
                let i = match index.get(&breakend.chromosome) {
                    Some(i) => *i,
                    None => {
                        index.insert(breakend.chromosome.clone(), map.chromosomes.len());
                        map.chromosomes.push(breakend.chromosome.clone());
                        map.lists.push(Vec::new());
                        map.chromosomes.len() - 1
                    }
                };
                map.lists[i].push(BreakendRef::new(v, se));
            }
        }
        for list in map.lists.iter_mut() {
            // BreakendPositionComparator, in a stable sort
            list.sort_by(|a, b| {
                let first = variants[a.var as usize].breakend(a.se as usize);
                let second = variants[b.var as usize].breakend(b.se as usize);
                first.position.cmp(&second.position).then_with(|| {
                    if first.orient == second.orient {
                        std::cmp::Ordering::Equal
                    } else if first.orient.is_forward() {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Greater
                    }
                })
            });
        }
        map.order = java_hash_map_order(&map.chromosomes);
        map
    }

    pub fn get(&self, chromosome: &str) -> Option<&Vec<BreakendRef>> {
        self.chromosomes.iter().position(|c| c == chromosome).map(|i| &self.lists[i])
    }

    /// the lists in HashMap iteration order
    pub fn iter(&self) -> impl Iterator<Item = (&String, &Vec<BreakendRef>)> {
        self.order.iter().map(|&i| (&self.chromosomes[i], &self.lists[i]))
    }
}

/// GenotypeIds.fromVcfHeader: (reference ordinal, tumor ordinal, reference id, tumor id as in the VCF)
pub struct GenotypeIds {
    pub reference_ordinal: i32,
    pub tumor_ordinal: i32,
    pub reference_id: String,
    pub tumor_id: String,
}

impl GenotypeIds {
    pub fn from_vcf_header(samples: &[String], reference_id: Option<&str>, tumor_id: Option<&str>) -> GenotypeIds {
        let mut ids = GenotypeIds { reference_ordinal: -1, tumor_ordinal: -1, reference_id: String::new(), tumor_id: String::new() };
        for check_exact in [true, false] {
            for (i, name) in samples.iter().enumerate() {
                let matches = |id: Option<&str>| id.is_some_and(|id| name == id || (!check_exact && name.contains(id)));
                if ids.tumor_ordinal < 0 && matches(tumor_id) {
                    ids.tumor_id = name.clone();
                    ids.tumor_ordinal = i as i32;
                } else if ids.reference_ordinal < 0 && matches(reference_id) {
                    ids.reference_id = name.clone();
                    ids.reference_ordinal = i as i32;
                }
            }
        }
        ids
    }
    pub fn has_reference(&self) -> bool {
        self.reference_ordinal >= 0
    }
    pub fn has_tumor(&self) -> bool {
        self.tumor_ordinal >= 0
    }
}
