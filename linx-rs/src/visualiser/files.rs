//! Package `visualiser.file` (readers) and the CNV inputs of SampleData: the vis_*.tsv files LINX writes,
//! AMBER's BAF file, COBALT's ratio file, PURPLE's segment file and somatic driver catalog.

use anyhow::{bail, Context, Result};
use hmf_common::java_text::{fields_index_map, parse_bool, parse_double, parse_int, split};

use super::genome::{human_chromosome_contains, GenomeRegion};

/// SegmentTerminal
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegmentTerminal {
    TELOMERE,
    CENTROMERE,
    NONE,
}

impl SegmentTerminal {
    pub fn from_str(s: &str) -> SegmentTerminal {
        match s {
            "C" => SegmentTerminal::CENTROMERE,
            "T" => SegmentTerminal::TELOMERE,
            _ => SegmentTerminal::NONE,
        }
    }
}

/// VisGeneAnnotationType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnnotationType {
    FUSION,
    DRIVER,
    DISRUPTION,
    PSEUDOGENE,
    EXON_LOST,
}

impl AnnotationType {
    fn parse(s: &str) -> Result<AnnotationType> {
        Ok(match s {
            "FUSION" => AnnotationType::FUSION,
            "DRIVER" => AnnotationType::DRIVER,
            "DISRUPTION" => AnnotationType::DISRUPTION,
            "PSEUDOGENE" => AnnotationType::PSEUDOGENE,
            "EXON_LOST" => AnnotationType::EXON_LOST,
            _ => bail!("No enum constant com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.{s}"),
        })
    }
}

const SV_TYPES: [&str; 7] = ["BND", "DEL", "DUP", "INF", "INS", "INV", "SGL"];
const RESOLVED_TYPES: [&str; 31] = [
    "NONE", "DUP_BE", "LOW_VAF", "PAIR_INF", "SGL_MAPPED_INF", "DEL", "DUP", "INF", "SGL", "INV", "INS", "UNBAL_TRANS", "SIMPLE_GRP",
    "RECIP_TRANS", "RECIP_TRANS_DEL_DUP", "RECIP_TRANS_DUPS", "RECIP_INV", "RECIP_INV_DUPS", "RECIP_INV_DEL_DUP", "DUP_TI", "DEL_TI",
    "UNBAL_TRANS_TI", "PAIR_OTHER", "RESOLVED_FOLDBACK", "FB_INV_PAIR", "SGL_PAIR_DEL", "SGL_PAIR_DUP", "SGL_PAIR_INS", "LINE", "COMPLEX",
    "DOUBLE_MINUTE",
];

/// VisSvData
#[derive(Clone, Debug)]
pub struct VisSvData {
    pub cluster_id: i32,
    pub chain_id: i32,
    pub sv_id: i32,
    pub sv_type: &'static str,
    pub resolved_type: &'static str,
    pub is_synthetic: bool,
    pub chr_start: String,
    pub chr_end: String,
    pub pos_start: i32,
    pub pos_end: i32,
    pub info_start: String,
    pub info_end: String,
    pub jcn: f64,
    pub in_double_minute: bool,
    pub frame: i32,
}

impl VisSvData {
    /// ResolvedType.isSimple() && !IsSynthetic
    pub fn is_simple_sv(&self) -> bool {
        matches!(self.resolved_type, "DEL" | "DUP" | "INS" | "SGL_PAIR_DEL" | "SGL_PAIR_DUP" | "SGL_PAIR_INS") && !self.is_synthetic
    }
    pub fn is_line_element(&self) -> bool {
        self.resolved_type == "LINE"
    }
    pub fn connectors_only(&self, show_simple_sv_segments: bool) -> bool {
        (!show_simple_sv_segments && self.is_simple_sv()) || self.is_line_element()
    }
    pub fn is_valid_start(&self) -> bool {
        human_chromosome_contains(&self.chr_start)
    }
    pub fn is_valid_end(&self) -> bool {
        human_chromosome_contains(&self.chr_end)
    }
}

/// VisSegment
#[derive(Clone, Debug)]
pub struct VisSegment {
    pub cluster_id: i32,
    pub chain_id: i32,
    pub chromosome: String,
    pub pos_start: String,
    pub pos_end: String,
    pub link_ploidy: f64,
    pub in_double_minute: bool,
    pub frame: i32,
    pub track: i32,
    pub start: i32,
    pub end: i32,
    pub terminal_start: SegmentTerminal,
    pub terminal_end: SegmentTerminal,
}

impl VisSegment {
    /// The constructor: a terminal start takes the end's position and vice versa.
    pub fn new(cluster_id: i32, chain_id: i32, chromosome: &str, pos_start: &str, pos_end: &str, link_ploidy: f64, in_dm: bool) -> Result<VisSegment> {
        let ts = SegmentTerminal::from_str(pos_start);
        let te = SegmentTerminal::from_str(pos_end);
        let int = |s: &str| parse_int(s).with_context(|| format!("For input string: \"{s}\""));
        let start = if ts == SegmentTerminal::NONE { int(pos_start)? } else { int(pos_end)? };
        let end = if te == SegmentTerminal::NONE { int(pos_end)? } else { int(pos_start)? };
        Ok(VisSegment {
            cluster_id,
            chain_id,
            chromosome: chromosome.to_string(),
            pos_start: pos_start.to_string(),
            pos_end: pos_end.to_string(),
            link_ploidy,
            in_double_minute: in_dm,
            frame: 0,
            track: 0,
            start,
            end,
            terminal_start: ts,
            terminal_end: te,
        })
    }
}

impl GenomeRegion for VisSegment {
    fn chromosome(&self) -> &str {
        &self.chromosome
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
}

/// VisCopyNumber
#[derive(Clone, Debug)]
pub struct VisCopyNumber {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub copy_number: f64,
    pub baf: f64,
    pub truncated: bool,
}

impl VisCopyNumber {
    pub fn minor_allele_copy_number(&self) -> f64 {
        f64::max(0.0, (1.0 - self.baf) * self.copy_number)
    }
}

impl GenomeRegion for VisCopyNumber {
    fn chromosome(&self) -> &str {
        &self.chromosome
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
}

/// VisGeneExon
#[derive(Clone, Debug)]
pub struct VisGeneExon {
    pub cluster_id: i32,
    pub gene: String,
    pub transcript: String,
    pub chromosome: String,
    pub annotation_type: AnnotationType,
    pub exon_rank: i32,
    pub exon_start: i32,
    pub exon_end: i32,
}

impl GenomeRegion for VisGeneExon {
    fn chromosome(&self) -> &str {
        &self.chromosome
    }
    fn start(&self) -> i32 {
        self.exon_start
    }
    fn end(&self) -> i32 {
        self.exon_end
    }
}

pub const PD_FIVE_PRIME_UTR: &str = "5-Prime UTR";
pub const PD_THREE_PRIME_UTR: &str = "3-Prime UTR";
pub const PD_NON_CODING: &str = "Non Coding";
pub const PROTEIN_DOMAIN_UTR: &str = "UTR/Non-coding";

/// VisProteinDomain
#[derive(Clone, Debug)]
pub struct VisProteinDomain {
    pub cluster_id: i32,
    pub transcript: String,
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub info: String,
}

impl VisProteinDomain {
    /// name(): UTR and non-coding regions are one category. Java tests the 5' UTR twice and never the 3'
    /// UTR (LV1), so 3' UTRs are drawn as a protein domain of their own.
    pub fn name(&self) -> &str {
        let three_prime = !crate::compat::java("LV1") && self.info == PD_THREE_PRIME_UTR;
        if self.info == PD_FIVE_PRIME_UTR || three_prime || self.info == PD_NON_CODING {
            PROTEIN_DOMAIN_UTR
        } else {
            &self.info
        }
    }
}

impl GenomeRegion for VisProteinDomain {
    fn chromosome(&self) -> &str {
        &self.chromosome
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
}

/// VisFusion
#[derive(Clone, Debug)]
pub struct VisFusion {
    pub cluster_id: i32,
    pub reportable: bool,
    pub gene_name_up: String,
    pub transcript_up: String,
    pub chr_up: String,
    pub pos_up: i32,
    pub strand_up: i32,
    pub region_type_up: String,
    pub fused_exon_up: i32,
    pub gene_name_down: String,
    pub transcript_down: String,
    pub chr_down: String,
    pub pos_down: i32,
    pub strand_down: i32,
    pub region_type_down: String,
    pub fused_exon_down: i32,
}

impl VisFusion {
    pub fn name(&self) -> String {
        format!("{}_{}", self.gene_name_up, self.gene_name_down)
    }
}

// ---------------------------------------------------------------------------------------------------------
// readers: Files.readAllLines, the header's field map, String.split(TSV_DELIM) (trailing empty fields dropped)

struct Table {
    fields: std::collections::HashMap<String, usize>,
    rows: Vec<String>,
}

fn read_table(path: &str) -> Result<Table> {
    let text = std::fs::read_to_string(path).with_context(|| path.to_string())?;
    let mut lines = text.lines();
    let header = lines.next().with_context(|| format!("{path}: no header"))?;
    let fields = fields_index_map(header, '\t');
    Ok(Table { fields, rows: lines.map(|l| l.to_string()).collect() })
}

struct Row<'a> {
    t: &'a Table,
    v: Vec<&'a str>,
}

impl<'a> Row<'a> {
    fn raw(&self, name: &str) -> Result<&'a str> {
        let i = *self.t.fields.get(name).with_context(|| format!("no column {name}"))?;
        self.v.get(i).copied().with_context(|| format!("Index {i} out of bounds for length {}", self.v.len()))
    }
    /// getValue: the field, or the default when the column is missing
    fn str_or(&self, name: &str, default: &'a str) -> Result<&'a str> {
        if self.t.fields.contains_key(name) { self.raw(name) } else { Ok(default) }
    }
    fn int_or(&self, name: &str, default: i32) -> Result<i32> {
        if !self.t.fields.contains_key(name) {
            return Ok(default);
        }
        let s = self.raw(name)?;
        parse_int(s).with_context(|| format!("For input string: \"{s}\""))
    }
    fn double_or(&self, name: &str, default: f64) -> Result<f64> {
        if !self.t.fields.contains_key(name) {
            return Ok(default);
        }
        let s = self.raw(name)?;
        parse_double(s).with_context(|| format!("For input string: \"{s}\""))
    }
    fn bool_or(&self, name: &str, default: bool) -> Result<bool> {
        if !self.t.fields.contains_key(name) {
            return Ok(default);
        }
        Ok(parse_bool(self.raw(name)?))
    }
    fn byte(&self, name: &str) -> Result<i8> {
        let s = self.raw(name)?;
        hmf_common::java_text::parse_byte(s).with_context(|| format!("For input string: \"{s}\""))
    }
}

fn rows(t: &Table) -> impl Iterator<Item = Row<'_>> {
    t.rows.iter().map(move |l| Row { t, v: split(l, '\t') })
}

fn enum_of(value: &str, names: &[&'static str], enum_name: &str) -> Result<&'static str> {
    names.iter().find(|n| **n == value).copied().with_context(|| format!("No enum constant {enum_name}.{value}"))
}

pub fn read_vis_sv_data(path: &str) -> Result<Vec<VisSvData>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(VisSvData {
            cluster_id: r.int_or("ClusterId", 0)?,
            chain_id: r.int_or("ChainId", 0)?,
            sv_id: r.int_or("SvId", 0)?,
            sv_type: enum_of(r.raw("Type")?, &SV_TYPES, "com.hartwig.hmftools.common.sv.StructuralVariantType")?,
            resolved_type: enum_of(r.raw("ResolvedType")?, &RESOLVED_TYPES, "com.hartwig.hmftools.linx.types.ResolvedType")?,
            is_synthetic: r.bool_or("IsSynthetic", false)?,
            chr_start: r.str_or("ChrStart", "")?.to_string(),
            chr_end: r.str_or("ChrEnd", "")?.to_string(),
            pos_start: r.int_or("PosStart", 0)?,
            pos_end: r.int_or("PosEnd", 0)?,
            info_start: {
                r.byte("OrientStart")?;
                r.byte("OrientEnd")?;
                r.str_or("InfoStart", "")?.to_string()
            },
            info_end: r.str_or("InfoEnd", "")?.to_string(),
            jcn: r.double_or("JunctionCopyNumber", 0.0)?,
            in_double_minute: r.bool_or("InDoubleMinute", false)?,
            frame: 0,
        });
    }
    Ok(out)
}

pub fn read_vis_segments(path: &str) -> Result<Vec<VisSegment>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(VisSegment::new(
            r.int_or("ClusterId", 0)?,
            r.int_or("ChainId", 0)?,
            r.str_or("Chromosome", "")?,
            r.str_or("PosStart", "")?,
            r.str_or("PosEnd", "")?,
            r.double_or("LinkPloidy", 0.0)?,
            r.bool_or("InDoubleMinute", false)?,
        )?);
    }
    Ok(out)
}

pub fn read_vis_copy_numbers(path: &str) -> Result<Vec<VisCopyNumber>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(VisCopyNumber {
            chromosome: r.str_or("Chromosome", "")?.to_string(),
            start: r.int_or("Start", 0)?,
            end: r.int_or("End", 0)?,
            copy_number: r.double_or("CopyNumber", 0.0)?,
            baf: r.double_or("BAF", 0.0)?,
            truncated: false,
        });
    }
    Ok(out)
}

pub fn read_vis_gene_exons(path: &str) -> Result<Vec<VisGeneExon>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(VisGeneExon {
            cluster_id: r.int_or("ClusterId", 0)?,
            gene: r.str_or("Gene", "")?.to_string(),
            transcript: r.str_or("Transcript", "")?.to_string(),
            chromosome: r.str_or("Chromosome", "")?.to_string(),
            annotation_type: AnnotationType::parse(r.raw("AnnotationType")?)?,
            exon_rank: r.int_or("ExonRank", 0)?,
            exon_start: r.int_or("ExonStart", 0)?,
            exon_end: r.int_or("ExonEnd", 0)?,
        });
    }
    Ok(out)
}

pub fn read_vis_protein_domains(path: &str) -> Result<Vec<VisProteinDomain>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(VisProteinDomain {
            cluster_id: r.int_or("ClusterId", 0)?,
            transcript: r.str_or("Transcript", "")?.to_string(),
            chromosome: r.str_or("Chromosome", "")?.to_string(),
            start: r.int_or("Start", 0)?,
            end: r.int_or("End", 0)?,
            info: r.str_or("Info", "")?.to_string(),
        });
    }
    Ok(out)
}

pub fn read_vis_fusions(path: &str) -> Result<Vec<VisFusion>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(VisFusion {
            cluster_id: r.int_or("ClusterId", 0)?,
            reportable: r.raw("Reportable")?.eq_ignore_ascii_case("true"),
            gene_name_up: r.str_or("GeneNameUp", "")?.to_string(),
            transcript_up: r.str_or("TranscriptUp", "")?.to_string(),
            chr_up: r.str_or("ChrUp", "")?.to_string(),
            pos_up: r.int_or("PosUp", 0)?,
            strand_up: r.int_or("StrandUp", 0)?,
            region_type_up: r.str_or("RegionTypeUp", "")?.to_string(),
            fused_exon_up: r.int_or("FusedExonUp", 0)?,
            gene_name_down: r.str_or("GeneNameDown", "")?.to_string(),
            transcript_down: r.str_or("TranscriptDown", "")?.to_string(),
            chr_down: r.str_or("ChrDown", "")?.to_string(),
            pos_down: r.int_or("PosDown", 0)?,
            strand_down: r.int_or("StrandDown", 0)?,
            region_type_down: r.str_or("RegionTypeDown", "")?.to_string(),
            fused_exon_down: r.int_or("FusedExonDown", 0)?,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------------------
// CNV inputs, kept per chromosome in file order (position order within a chromosome in AMBER's and
// COBALT's files; the visualiser filters them by region, which keeps that order)

/// AMBER's BAF points of one chromosome: (position, tumorBAF).
#[derive(Clone, Debug, Default)]
pub struct ChrPoints {
    pub chromosome: String,
    pub points: Vec<(i32, f64)>,
    /// whether the positions ascend (then a region's points are found by binary search)
    pub sorted: bool,
}

impl ChrPoints {
    /// The points with start <= position <= end, in file order.
    pub fn in_range(&self, start: i32, end: i32) -> Vec<(i32, f64)> {
        if self.sorted {
            let a = self.points.partition_point(|p| p.0 < start);
            let b = self.points.partition_point(|p| p.0 <= end);
            if a >= b {
                return Vec::new();
            }
            return self.points[a..b].to_vec();
        }
        self.points.iter().filter(|p| start <= p.0 && end >= p.0).copied().collect()
    }
}

fn finish(chrs: &mut [ChrPoints]) {
    for c in chrs.iter_mut() {
        c.sorted = c.points.windows(2).all(|w| w[0].0 <= w[1].0);
    }
}

fn line_reader(path: &str) -> Result<hmf_common::gz::LineReader> {
    hmf_common::gz::LineReader::open(path).with_context(|| path.to_string())
}

/// AmberBAFFile.read(file, true): chromosome, position and tumorBAF of each point.
pub fn read_amber_bafs(path: &str) -> Result<Vec<ChrPoints>> {
    let mut r = line_reader(path)?;
    let header = match r.next_line()? {
        Some(l) => String::from_utf8_lossy(l).to_string(),
        None => return Ok(Vec::new()),
    };
    let f = fields_index_map(&header, '\t');
    let col = |n: &str| f.get(n).copied().with_context(|| format!("{path}: no column {n}"));
    let (ci, pi, bi) = (col("chromosome")?, col("position")?, col("tumorBAF")?);
    let mut out: Vec<ChrPoints> = Vec::new();
    while let Some(line) = r.next_line()? {
        let line = std::str::from_utf8(line)?;
        let v: Vec<&str> = line.split('\t').collect();
        let chr = v[ci];
        let pos = parse_int(v[pi]).with_context(|| format!("For input string: \"{}\"", v[pi]))?;
        let baf = parse_double(v[bi]).with_context(|| format!("For input string: \"{}\"", v[bi]))?;
        if out.last().is_none_or(|c| c.chromosome != chr) {
            match out.iter().position(|c| c.chromosome == chr) {
                Some(i) => {
                    let c = out.remove(i);
                    out.push(c);
                }
                None => out.push(ChrPoints { chromosome: chr.to_string(), points: Vec::new(), sorted: true }),
            }
        }
        out.last_mut().unwrap().points.push((pos, baf));
    }
    finish(&mut out);
    Ok(out)
}

/// CobaltRatioFile.read: chromosome, position and tumorGCRatio, without the windows of ratio -1 (SampleData
/// drops them).
pub fn read_cobalt_ratios(path: &str) -> Result<Vec<ChrPoints>> {
    let mut r = line_reader(path)?;
    let header = match r.next_line()? {
        Some(l) => String::from_utf8_lossy(l).to_string(),
        None => return Ok(Vec::new()),
    };
    let f = fields_index_map(&header, '\t');
    let col = |n: &str| f.get(n).copied().with_context(|| format!("{path}: no column {n}"));
    let (ci, pi, ti) = (col("chromosome")?, col("position")?, col("tumorGCRatio")?);
    let mut out: Vec<ChrPoints> = Vec::new();
    while let Some(line) = r.next_line()? {
        let line = std::str::from_utf8(line)?;
        let v: Vec<&str> = line.split('\t').collect();
        let chr = v[ci];
        if out.last().is_none_or(|c| c.chromosome != chr) {
            match out.iter().position(|c| c.chromosome == chr) {
                Some(i) => {
                    let c = out.remove(i);
                    out.push(c);
                }
                None => out.push(ChrPoints { chromosome: chr.to_string(), points: Vec::new(), sorted: true }),
            }
        }
        let ratio = parse_double(v[ti]).with_context(|| format!("For input string: \"{}\"", v[ti]))?;
        if ratio == -1.0 {
            continue;
        }
        let pos = parse_int(v[pi]).with_context(|| format!("For input string: \"{}\"", v[pi]))?;
        out.last_mut().unwrap().points.push((pos, ratio));
    }
    finish(&mut out);
    Ok(out)
}

/// The fields of a PurpleSegment the visualiser uses.
#[derive(Clone, Debug)]
pub struct PurpleSegment {
    pub chromosome: String,
    pub pos_start: i32,
    pub pos_end: i32,
    pub germline_diploid: bool,
    pub observed_baf: f64,
    pub observed_tumor_ratio: f64,
    pub max_start: i32,
}

pub fn read_purple_segments(path: &str) -> Result<Vec<PurpleSegment>> {
    let t = read_table(path)?;
    let mut out = Vec::new();
    for r in rows(&t) {
        out.push(PurpleSegment {
            chromosome: r.raw("chromosome")?.to_string(),
            pos_start: r.int_or("start", 0)?,
            pos_end: r.int_or("end", 0)?,
            germline_diploid: r.raw("germlineStatus")? == "DIPLOID",
            observed_baf: r.double_or("observedBAF", 0.0)?,
            observed_tumor_ratio: r.double_or("observedTumorRatio", 0.0)?,
            max_start: r.int_or("maxStart", 0)?,
        });
    }
    Ok(out)
}

/// The purple somatic driver catalog's AMP, PARTIAL_AMP, DEL and HET_DEL entries, the first of each gene and
/// driver type (SampleData.PurpleGeneCNVs).
pub fn read_purple_gene_cnvs(path: &str) -> Result<Vec<(String, String)>> {
    let catalog = hmf_common::purple_files::read_driver_catalog(path)?;
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for d in catalog {
        let driver = d.driver.to_string();
        if !matches!(driver.as_str(), "AMP" | "PARTIAL_AMP" | "DEL" | "HET_DEL") {
            continue;
        }
        if seen.insert(format!("{}-{}", d.gene, driver)) {
            out.push((d.gene.clone(), driver));
        }
    }
    Ok(out)
}

/// For the chromosome check of SvVisualiser.submitChromosome
pub fn all_human(chromosomes: &[String]) -> bool {
    chromosomes.iter().all(|c| human_chromosome_contains(c))
}
