//! Small types shared by the caller: orientation, SV type, filter types, regions, alternative
//! alignments, and the Java HashMap iteration order that some steps depend on.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Orient {
    Forward,
    Reverse,
}

impl Orient {
    pub fn is_forward(self) -> bool {
        self == Orient::Forward
    }
    pub fn as_byte(self) -> i8 {
        if self.is_forward() { 1 } else { -1 }
    }
    /// Orientation.fromByte
    pub fn from_byte(b: i8) -> Orient {
        match b {
            1 => Orient::Forward,
            -1 => Orient::Reverse,
            _ => panic!("Invalid orientation: {b}"),
        }
    }
    /// Orientation.fromByteStr
    pub fn from_byte_str(s: &str) -> Orient {
        Orient::from_byte(s.parse::<i8>().unwrap_or_else(|_| panic!("NumberFormatException: For input string: \"{s}\"")))
    }
    /// Orientation.fromChar
    pub fn from_char(c: char) -> Orient {
        match c {
            '+' => Orient::Forward,
            '-' => Orient::Reverse,
            _ => panic!("Invalid orientation: {c}"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SvType {
    Bnd,
    Del,
    Dup,
    Inf,
    Ins,
    Inv,
    Sgl,
}

impl SvType {
    pub fn name(self) -> &'static str {
        match self {
            SvType::Bnd => "BND",
            SvType::Del => "DEL",
            SvType::Dup => "DUP",
            SvType::Inf => "INF",
            SvType::Ins => "INS",
            SvType::Inv => "INV",
            SvType::Sgl => "SGL",
        }
    }
    /// StructuralVariantType.valueOf
    pub fn value_of(s: &str) -> SvType {
        match s {
            "BND" => SvType::Bnd,
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            "INF" => SvType::Inf,
            "INS" => SvType::Ins,
            "INV" => SvType::Inv,
            "SGL" => SvType::Sgl,
            _ => panic!("IllegalArgumentException: No enum constant StructuralVariantType.{s}"),
        }
    }
    /// SvUtils.isIndel
    pub fn is_indel(self) -> bool {
        matches!(self, SvType::Del | SvType::Dup | SvType::Ins)
    }
    /// Variant.hasLength
    pub fn has_length(self) -> bool {
        matches!(self, SvType::Ins | SvType::Inv | SvType::Del | SvType::Dup)
    }
}

/// SvUtils.formSvType
pub fn form_sv_type(chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: Orient, orient_end: Orient, has_inserted_bases: bool) -> SvType {
    if chr_start != chr_end {
        return SvType::Bnd;
    }
    if orient_start != orient_end {
        let pos_diff = (pos_start - pos_end).abs();
        if pos_diff == 1 && has_inserted_bases {
            return SvType::Ins;
        }
        if pos_diff == 0 {
            return SvType::Dup;
        }
        let first_is_lower = pos_start < pos_end;
        if first_is_lower == orient_start.is_forward() { SvType::Del } else { SvType::Dup }
    } else {
        SvType::Inv
    }
}

/// FilterType, in the order of the Java enum
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum FilterType {
    MinAf,
    MinSupport,
    MinQuality,
    StrandBias,
    SglStrandBias,
    MinLength,
    MinAnchorLength,
    ShortFragLength,
    Duplicate,
    Sgl,
    InvShortLowVafHom,
    InvShortFragLowVaf,
    InvShortIsolated,
    DelShortLowVaf,
    UnpairedThreePrimeRange,
    LineSource,
    SbxStrandBias,
    SbxArtefact,
    SbxInvZeroMotif,
    OffTarget,
    WeakJuncMinAf,
    Pon,
}

pub const FILTER_TYPES: [FilterType; 22] = [
    FilterType::MinAf,
    FilterType::MinSupport,
    FilterType::MinQuality,
    FilterType::StrandBias,
    FilterType::SglStrandBias,
    FilterType::MinLength,
    FilterType::MinAnchorLength,
    FilterType::ShortFragLength,
    FilterType::Duplicate,
    FilterType::Sgl,
    FilterType::InvShortLowVafHom,
    FilterType::InvShortFragLowVaf,
    FilterType::InvShortIsolated,
    FilterType::DelShortLowVaf,
    FilterType::UnpairedThreePrimeRange,
    FilterType::LineSource,
    FilterType::SbxStrandBias,
    FilterType::SbxArtefact,
    FilterType::SbxInvZeroMotif,
    FilterType::OffTarget,
    FilterType::WeakJuncMinAf,
    FilterType::Pon,
];

impl FilterType {
    pub fn vcf_tag(self) -> &'static str {
        self.tag_desc().0
    }
    pub fn vcf_desc(self) -> &'static str {
        self.tag_desc().1
    }
    fn tag_desc(self) -> (&'static str, &'static str) {
        match self {
            FilterType::MinAf => ("minAF", "Variant allele fraction too low"),
            FilterType::MinSupport => ("minSupport", "Below minimum fragment support"),
            FilterType::MinQuality => ("minQual", "Insufficient quality"),
            FilterType::StrandBias => ("strandBias", "Outside valid strand bias range"),
            FilterType::SglStrandBias => ("sglStrandBias", "Single breakend with excessive strand bias and not a mobile line insertion"),
            FilterType::MinLength => ("minLength", "Variant is too short"),
            FilterType::MinAnchorLength => ("minAnchorLength", "Trimmed anchor length is too short"),
            FilterType::ShortFragLength => ("shortFrags", "Average variant fragment length is too short"),
            FilterType::Duplicate => ("dedup", "Event is duplicate of another"),
            FilterType::Sgl => ("sgl", "SGLs filtered entirely"),
            FilterType::InvShortLowVafHom => ("invShortLowVafHom", "Short (<3K) low-VAF INV artefacts with homology"),
            FilterType::InvShortFragLowVaf => ("invShortFragLowVaf", "Short (<300b) low-VAF INV artefacts"),
            FilterType::InvShortIsolated => ("invShortIsolated", "Short (<100b) isolated INV artefacts"),
            FilterType::DelShortLowVaf => ("delShortLowVaf", "Short low-VAF DEL artefacts with homology"),
            FilterType::UnpairedThreePrimeRange => ("unpairedPositionRange", "Unpaired position range"),
            FilterType::LineSource => ("lineSource", "Remote LINE source"),
            FilterType::SbxStrandBias => ("sbxStrandBias", "SBX strand-bias"),
            FilterType::SbxArtefact => ("sbxArtefact", "SBX artefact"),
            FilterType::SbxInvZeroMotif => ("sbxInvZeroMotif", "SBX artefact, zero-length INV with known motif"),
            FilterType::OffTarget => ("offTarget", "Outside target panel"),
            FilterType::WeakJuncMinAf => ("weakJuncMinAF", "Low allele fraction at weak junction"),
            FilterType::Pon => ("PON", "Found in panel of normals"),
        }
    }
    /// the Java enum constant name, as logged
    pub fn java_name(self) -> &'static str {
        match self {
            FilterType::MinAf => "MIN_AF",
            FilterType::MinSupport => "MIN_SUPPORT",
            FilterType::MinQuality => "MIN_QUALITY",
            FilterType::StrandBias => "STRAND_BIAS",
            FilterType::SglStrandBias => "SGL_STRAND_BIAS",
            FilterType::MinLength => "MIN_LENGTH",
            FilterType::MinAnchorLength => "MIN_ANCHOR_LENGTH",
            FilterType::ShortFragLength => "SHORT_FRAG_LENGTH",
            FilterType::Duplicate => "DUPLICATE",
            FilterType::Sgl => "SGL",
            FilterType::InvShortLowVafHom => "INV_SHORT_LOW_VAF_HOM",
            FilterType::InvShortFragLowVaf => "INV_SHORT_FRAG_LOW_VAF",
            FilterType::InvShortIsolated => "INV_SHORT_ISOLATED",
            FilterType::DelShortLowVaf => "DEL_SHORT_LOW_VAF",
            FilterType::UnpairedThreePrimeRange => "UNPAIRED_THREE_PRIME_RANGE",
            FilterType::LineSource => "LINE_SOURCE",
            FilterType::SbxStrandBias => "SBX_STRAND_BIAS",
            FilterType::SbxArtefact => "SBX_ARTEFACT",
            FilterType::SbxInvZeroMotif => "SBX_INV_ZERO_MOTIF",
            FilterType::OffTarget => "OFF_TARGET",
            FilterType::WeakJuncMinAf => "WEAK_JUNC_MIN_AF",
            FilterType::Pon => "PON",
        }
    }
    /// only PON is specific to one of the two output files
    pub fn germline_only(self) -> bool {
        self == FilterType::Pon
    }
    pub fn from_vcf_tag(tag: &str) -> Option<FilterType> {
        FILTER_TYPES.iter().copied().find(|f| f.vcf_tag() == tag)
    }
}

/// A set of filter types. Java keeps a HashSet of enum values, whose iteration order is not defined;
/// here the order is that of the enum.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FilterSet(u32);

impl FilterSet {
    pub fn add(&mut self, f: FilterType) {
        self.0 |= 1 << (f as u8);
    }
    pub fn remove(&mut self, f: FilterType) {
        self.0 &= !(1 << (f as u8));
    }
    pub fn contains(&self, f: FilterType) -> bool {
        self.0 & (1 << (f as u8)) != 0
    }
    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }
    pub fn len(&self) -> usize {
        self.0.count_ones() as usize
    }
    pub fn clear(&mut self) {
        self.0 = 0;
    }
    pub fn iter(&self) -> impl Iterator<Item = FilterType> + '_ {
        FILTER_TYPES.iter().copied().filter(|f| self.contains(*f))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefGenomeVersion {
    V37,
    V38,
}

pub const CHR_PREFIX: &str = "chr";

/// RefGenomeFunctions.stripChrPrefix
pub fn strip_chr_prefix(chromosome: &str) -> &str {
    chromosome.strip_prefix(CHR_PREFIX).unwrap_or(chromosome)
}

impl RefGenomeVersion {
    /// RefGenomeVersion.from(String)
    pub fn from_str(version: &str) -> Option<RefGenomeVersion> {
        match version {
            "V37" | "37" | "HG37" => Some(RefGenomeVersion::V37),
            "V38" | "38" | "HG38" => Some(RefGenomeVersion::V38),
            _ => None,
        }
    }
    pub fn is37(self) -> bool {
        self == RefGenomeVersion::V37
    }
    pub fn versioned_chromosome(self, chromosome: &str) -> String {
        match self {
            RefGenomeVersion::V38 => {
                if chromosome.starts_with(CHR_PREFIX) { chromosome.to_string() } else { format!("{CHR_PREFIX}{chromosome}") }
            }
            RefGenomeVersion::V37 => strip_chr_prefix(chromosome).to_string(),
        }
    }
}

/// HumanChromosome.values() as strings: 1-22, X, Y
pub fn human_chromosomes() -> Vec<String> {
    let mut out: Vec<String> = (1..=22).map(|i| i.to_string()).collect();
    out.push("X".to_string());
    out.push("Y".to_string());
    out
}

/// HumanChromosome.contains
pub fn is_human_chromosome(chromosome: &str) -> bool {
    let trimmed = strip_chr_prefix(chromosome);
    if trimmed.chars().all(|c| c.is_ascii_digit()) {
        return match trimmed.parse::<i32>() {
            Ok(value) => (1..=22).contains(&value),
            // Java (C16): an empty name counts as numeric and the parse throws, as does a number that overflows
            Err(_) if crate::common::compat::java("C16") => panic!("NumberFormatException: For input string: \"{trimmed}\""),
            Err(_) => false,
        };
    }
    trimmed == "X" || trimmed == "Y"
}

/// HumanChromosome.fromString(...).toString(): the name without the prefix, whatever its case
pub fn human_chromosome_name(chromosome: &str) -> String {
    if chromosome.to_lowercase().starts_with(CHR_PREFIX) { chromosome[3..].to_string() } else { chromosome.to_string() }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChrBaseRegion {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

impl ChrBaseRegion {
    pub fn new(chromosome: &str, start: i32, end: i32) -> ChrBaseRegion {
        ChrBaseRegion { chromosome: chromosome.to_string(), start, end }
    }
    pub fn contains_position(&self, chromosome: &str, position: i32) -> bool {
        self.chromosome == chromosome && position >= self.start && position <= self.end
    }
}

/// BaseRegion.positionsOverlap
pub fn positions_overlap(start1: i32, end1: i32, start2: i32, end2: i32) -> bool {
    !(start1 > end2 || end1 < start2)
}

/// The order in which a java.util.HashMap with String keys iterates, given the keys in insertion
/// order: by bucket of the final table, and by insertion within a bucket. (Buckets only turn into
/// trees from eight entries in a table of 64, which the maps here never reach.)
pub fn java_hash_map_order(keys: &[String]) -> Vec<usize> {
    let mut table: usize = 16;
    while keys.len() as f64 > table as f64 * 0.75 {
        table *= 2;
    }
    let mut order: Vec<usize> = (0..keys.len()).collect();
    order.sort_by_key(|&i| {
        let mut h: i32 = 0;
        for unit in keys[i].encode_utf16() {
            h = h.wrapping_mul(31).wrapping_add(unit as i32);
        }
        let h = h as u32;
        ((h ^ (h >> 16)) as usize) & (table - 1)
    });
    order
}

/// StructuralVariantFactory.trimStringListValue
pub fn trim_string_list_value(value: &str) -> String {
    value.chars().filter(|c| !matches!(c, '[' | ']' | ' ')).collect()
}

/// An alternative alignment of an insert sequence, from the INSALN tag
#[derive(Clone, Debug, PartialEq)]
pub struct AlternativeAlignment {
    pub chromosome: String,
    pub position: i32,
    pub orient: Orient,
    pub cigar: String,
    pub num_mutations: i32,
}

/// Java String.split(regex) without a limit: trailing empty strings are dropped
fn java_split<'a>(s: &'a str, delim: char) -> Vec<&'a str> {
    let mut parts: Vec<&str> = s.split(delim).collect();
    while parts.len() > 1 && parts.last() == Some(&"") {
        parts.pop();
    }
    if parts.len() == 1 && parts[0].is_empty() && !s.is_empty() {
        parts.clear();
    }
    parts
}

impl AlternativeAlignment {
    /// AlternativeAlignment.fromVcfTag: `chr:pos|orient|cigar|mutations` items separated by commas
    pub fn from_vcf_tag(vcf_tag: &str) -> Vec<AlternativeAlignment> {
        let trimmed = trim_string_list_value(vcf_tag);
        if trimmed.is_empty() {
            return Vec::new();
        }
        trimmed.split(',').filter_map(Self::from_vcf_string).collect()
    }

    /// AlternativeAlignment.vcfString
    #[cfg(test)]
    pub fn vcf_string(&self) -> String {
        format!("{}:{}|{}|{}|{}", self.chromosome, self.position, if self.orient.is_forward() { '+' } else { '-' }, self.cigar, self.num_mutations)
    }

    /// AlternativeAlignment.toVcfTag
    #[cfg(test)]
    pub fn to_vcf_tag(alignments: &[AlternativeAlignment]) -> String {
        alignments.iter().map(|a| a.vcf_string()).collect::<Vec<_>>().join(",")
    }

    /// AlternativeAlignment.fromLocationTag: the aligner's `chr,±pos,cigar,mutations;` list
    #[cfg(test)]
    pub fn from_location_tag(location_tag: &str) -> Vec<AlternativeAlignment> {
        if location_tag.is_empty() {
            return Vec::new();
        }
        location_tag
            .split(';')
            .filter_map(|mapping| {
                let items: Vec<&str> = mapping.splitn(4, ',').collect();
                if items.len() != 4 || !is_human_chromosome(items[0]) {
                    return None;
                }
                Some(AlternativeAlignment {
                    chromosome: items[0].to_string(),
                    orient: Orient::from_char(items[1].chars().next()?),
                    position: items[1][1..].parse().ok()?,
                    cigar: items[2].to_string(),
                    num_mutations: items[3].parse().ok()?,
                })
            })
            .collect()
    }

    fn from_vcf_string(s: &str) -> Option<AlternativeAlignment> {
        let values: Vec<&str> = s.split('|').collect();
        if values.len() != 4 {
            return None;
        }
        let location = java_split(values[0], ':');
        if location.len() != 2 {
            return None;
        }
        if values[1].is_empty() {
            return None;
        }
        let chromosome = location[0];
        if !is_human_chromosome(chromosome) {
            return None;
        }
        let position = crate::caller::vcf::java_parse_int(location[1])
            .unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{}\"", location[1]));
        let orient = Orient::from_char(values[1].chars().next().unwrap());
        let num_mutations =
            crate::caller::vcf::java_parse_int(values[3]).unwrap_or_else(|| panic!("NumberFormatException: For input string: \"{}\"", values[3]));
        Some(AlternativeAlignment { chromosome: chromosome.to_string(), position, orient, cigar: values[2].to_string(), num_mutations })
    }
}

/// CigarUtils.calcCigarAlignedLength
pub fn calc_cigar_aligned_length(cigar: &str) -> i32 {
    let mut base_length = 0;
    let mut current = 0i32;
    for c in cigar.chars() {
        if matches!(c, 'D' | 'M' | 'N') {
            base_length += current;
            current = 0;
            continue;
        }
        let digit = c as i32 - '0' as i32;
        if (0..=9).contains(&digit) {
            current = current.wrapping_mul(10).wrapping_add(digit);
        } else {
            current = 0;
        }
    }
    base_length
}

/// CigarUtils.cigarElementsFromStr: (length, operator)
pub fn cigar_elements(cigar: &str) -> Vec<(i32, char)> {
    if cigar == "*" {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut length = 0i32;
    for c in cigar.chars() {
        let digit = c as i32 - '0' as i32;
        if (0..=9).contains(&digit) {
            length = length.wrapping_mul(10).wrapping_add(digit);
        } else {
            if !"MIDNSHP=X".contains(c) {
                panic!("IllegalArgumentException: Unrecognized CigarOperator: {c}");
            }
            out.push((length, c));
            length = 0;
        }
    }
    out
}

/// Nucleotides.reverseComplementBases
pub fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}
