//! The visualiser: one HTML file per variant showing the reads behind it. Ports of `vis/VisConfig`,
//! `VariantVis`, `ReadEvidenceRecord`, `ReadTableColumn`, `RenderUtils`, `VisFileBuilder` and of
//! `hmf-common` `vis/HtmlUtil`. The drawing is in `vis_svg`, the CSS and element writers in `vis_html`
//! and the amino acid track in `vis_aa`.

use crate::constants::{jmax, jmin, jround};
use crate::evidence::{FragmentData, ReadContextCounter};
use crate::filters::SageVariant;
use crate::fmt::java_fixed;
use crate::matcher::MatchType;
use crate::quality::QualityScores;
use crate::read_context::VariantReadContext;
use crate::refseq::RefGenome;
use crate::types::*;
use crate::vis_aa::{self, AminoAcidEvent, Transcript};
use crate::vis_html::*;
use crate::vis_svg::*;
use anyhow::{bail, Context, Result};
use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

// SageVisConstants
pub const READ_HEIGHT_PX: f64 = 12.0;
pub const READ_EXTEND_LENGTH: i32 = 75;
const DISPLAY_EVERY_NTH_COORD: i32 = 10;
const MAX_MAPQ_SHADING_CUTTOFF: i32 = 60;
const INSERT_SIZE_CUTTOFF: i32 = 1000;
pub const MAX_READ_UPPER_LIMIT: usize = 1000;

static VIS_ENABLED: AtomicBool = AtomicBool::new(false);

/// Whether the visualiser is on: alignment reading then also keeps the tags only it shows
pub fn is_enabled() -> bool {
    VIS_ENABLED.load(Ordering::Relaxed)
}

pub fn set_enabled(on: bool) {
    VIS_ENABLED.store(on, Ordering::Relaxed);
}

/// VariantVis.VARIANT_INDEXED_BASES_MAP: the read contexts seen for each variant, process-wide
static VARIANT_INDEXED_BASES: Mutex<Option<HashMap<String, BTreeSet<String>>>> = Mutex::new(None);

const MATCH_TYPE_COUNT: usize = 7;
const MATCH_TYPE_NAMES: [&str; MATCH_TYPE_COUNT] = ["NONE", "CORE", "PARTIAL_CORE", "FULL", "REALIGNED", "REF", "SIMPLE_ALT"];
/// SageVisConstants.MAX_READS_PER_TYPE, by ReadContextMatch ordinal
const MAX_READS_PER_TYPE: [usize; MATCH_TYPE_COUNT] = [5, 10, 10, 40, 20, 20, 5];
/// VisFileBuilder.SORTED_MATCH_TYPES: FULL, PARTIAL_CORE, REALIGNED, CORE, REF, NONE, SIMPLE_ALT
const SORTED_MATCH_TYPES: [usize; MATCH_TYPE_COUNT] = [3, 2, 4, 1, 5, 0, 6];

/// VariantImpact, as far as the visualiser reads it
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VariantImpact {
    pub gene_name: String,
    pub canonical_transcript: String,
    pub canonical_effect: String,
    pub canonical_hgvs_protein: String,
}

/// VariantImpactSerialiser.fromAttributeValues
pub fn impact_from_values(values: &[String]) -> Result<VariantImpact, String> {
    if values.len() != 10 {
        return Ok(VariantImpact::default());
    }
    const CODING_EFFECTS: [&str; 6] = ["NONSENSE_OR_FRAMESHIFT", "SPLICE", "MISSENSE", "SYNONYMOUS", "NONE", "UNDEFINED"];
    for i in [3, 8] {
        if !CODING_EFFECTS.contains(&values[i].as_str()) {
            return Err(format!("No enum constant CodingEffect.{}", values[i]));
        }
    }
    values[9].parse::<i32>().map_err(|_| format!("For input string: \"{}\"", values[9]))?;
    Ok(VariantImpact {
        gene_name: values[0].clone(),
        canonical_transcript: values[1].clone(),
        canonical_effect: values[2].clone(),
        canonical_hgvs_protein: values[6].clone(),
    })
}

/// A reportable variant of the PURPLE VCF with its IMPACT values
#[derive(Clone, Debug)]
pub struct ReportableVariant {
    pub chromosome: String,
    pub variant: SimpleVariant,
    pub impact_values: Vec<String>,
}

/// VisConfig
#[derive(Clone, Debug, Default)]
pub struct VisConfig {
    pub pass_only: bool,
    pub max_support_reads: i32,
    pub output_dir: String,
    pub specific_variants: Vec<(String, SimpleVariant)>,
    pub purple_vcf: Option<String>,
    pub reportable_variants: Vec<ReportableVariant>,
    /// the transcripts the reportable variants name, by transcript name
    pub transcripts: HashMap<String, Transcript>,
}

/// SimpleVariant.fromConfig: chromosome:position:ref:alt separated by ';'
pub fn parse_vis_variants(config: &str) -> Result<Vec<(String, SimpleVariant)>> {
    let mut variants = Vec::new();
    for item in config.split(';') {
        let parts: Vec<&str> = item.splitn(4, ':').collect();
        if parts.len() != 4 {
            bail!("invalid vis_variants entry: {item}");
        }
        let position: i32 = parts[1].parse().with_context(|| format!("invalid vis_variants position: {item}"))?;
        variants.push((parts[0].to_string(), SimpleVariant::new(position, parts[2].as_bytes(), parts[3].as_bytes())));
    }
    Ok(variants)
}

/// The records of a PURPLE VCF with the REPORTED flag
pub fn load_reportable_variants(path: &str) -> Result<Vec<ReportableVariant>> {
    let mut out = Vec::new();
    for line in crate::refdata::read_lines(path)? {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }
        let mut reported = false;
        let mut impact_values: Vec<String> = Vec::new();
        for item in fields[7].split(';') {
            match item.split_once('=') {
                None if item == "REPORTED" => reported = true,
                Some(("REPORTED", value)) => reported = value.eq_ignore_ascii_case("true"),
                Some(("IMPACT", value)) => impact_values = value.split(',').map(|v| v.to_string()).collect(),
                _ => {}
            }
        }
        if !reported {
            continue;
        }
        // getAltAlleleWithHighestAlleleCount: the first alt unless allele counts say otherwise (not read here)
        let alt = fields[4].split(',').next().unwrap_or("");
        out.push(ReportableVariant {
            chromosome: fields[0].to_string(),
            variant: SimpleVariant::new(fields[1].parse()?, fields[3].as_bytes(), alt.as_bytes()),
            impact_values,
        });
    }
    Ok(out)
}

impl VisConfig {
    /// VisConfig(configBuilder, outputDir)
    pub fn from_options(
        variants: Option<&str>, pass_only: bool, max_support_reads: i32, output_dir: Option<&str>, purple_vcf: Option<&str>,
        output_vcf: &str,
    ) -> Result<VisConfig> {
        let mut cfg = VisConfig {
            pass_only,
            max_support_reads,
            purple_vcf: purple_vcf.map(|s| s.to_string()),
            ..Default::default()
        };
        if let Some(v) = variants {
            cfg.specific_variants = parse_vis_variants(v)?;
            eprintln!("loaded {} specific variants for vis plot generation", cfg.specific_variants.len());
        } else if let Some(path) = purple_vcf {
            cfg.reportable_variants =
                load_reportable_variants(path).with_context(|| format!("failed to read vis_purple_vcf {path}"))?;
            cfg.specific_variants = cfg.reportable_variants.iter().map(|r| (r.chromosome.clone(), r.variant.clone())).collect();
            eprintln!("loaded {} reportable variants for vis plot generation", cfg.specific_variants.len());
        }
        cfg.output_dir = match output_dir {
            Some(dir) => dir.to_string(),
            None => {
                let parent = std::path::Path::new(output_vcf).parent().map(|p| p.to_string_lossy().to_string()).unwrap_or_default();
                if parent.is_empty() { "vis".to_string() } else { format!("{}/vis", parent.trim_end_matches('/')) }
            }
        };
        if cfg.enabled() {
            std::fs::create_dir_all(&cfg.output_dir).with_context(|| format!("failed to create {}", cfg.output_dir))?;
        }
        Ok(cfg)
    }
    pub fn enabled(&self) -> bool {
        self.purple_vcf.is_some() || self.pass_only || !self.specific_variants.is_empty()
    }
    pub fn has_variants(&self) -> bool {
        !self.specific_variants.is_empty()
    }
    pub fn visualiser_only_mode(&self) -> bool {
        self.purple_vcf.is_some()
    }
    /// processVariant
    pub fn process_variant(&self, chr: &str, variant: &SimpleVariant) -> bool {
        self.specific_variants.is_empty() || self.specific_variants.iter().any(|(c, v)| c == chr && v == variant)
    }
    fn reportable(&self, chr: &str, variant: &SimpleVariant) -> Option<&ReportableVariant> {
        self.reportable_variants.iter().find(|r| r.chromosome == chr && &r.variant == variant)
    }
    /// The Ensembl data ReferenceData loads for the amino acid track, limited to the transcripts the
    /// reportable variants name
    pub fn load_transcripts(&mut self, ensembl_dir: Option<&str>) -> Result<()> {
        if self.purple_vcf.is_none() || self.reportable_variants.is_empty() {
            return Ok(());
        }
        let Some(dir) = ensembl_dir else { return Ok(()) };
        let names: BTreeSet<String> =
            self.reportable_variants.iter().filter(|r| r.impact_values.len() == 10).map(|r| r.impact_values[1].clone()).collect();
        self.transcripts = load_transcripts(dir, &names)?;
        Ok(())
    }
}

/// TranscriptAminoAcids joined with the TranscriptData of the same gene and transcript name
pub fn load_transcripts(dir: &str, names: &BTreeSet<String>) -> Result<HashMap<String, Transcript>> {
    let dir = dir.trim_end_matches('/');
    let mut out: HashMap<String, Transcript> = HashMap::new();
    if names.is_empty() {
        return Ok(out);
    }
    let lines = crate::refdata::read_lines(&format!("{dir}/ensembl_trans_amino_acids.csv"))?;
    let column = |header: &str, name: &str| header.split(',').position(|f| f == name).with_context(|| format!("missing column {name}"));
    let header = lines.first().context("empty amino acid file")?;
    let (gi, ni, ti, ai) = (column(header, "GeneId")?, column(header, "GeneName")?, column(header, "TransName")?, column(header, "AminoAcids")?);
    for line in &lines[1..] {
        let values: Vec<&str> = line.split(',').collect();
        if values.len() <= ai || !names.contains(values[ti]) {
            continue;
        }
        out.insert(
            values[ti].to_string(),
            Transcript {
                gene_id: values[gi].to_string(),
                gene_name: values[ni].to_string(),
                trans_name: values[ti].to_string(),
                pos_strand: true,
                coding_start: None,
                coding_end: None,
                exons: Vec::new(),
                amino_acids: values[ai].as_bytes().to_vec(),
            },
        );
    }
    let lines = crate::refdata::read_lines(&format!("{dir}/ensembl_trans_exon_data.csv"))?;
    let header = lines.first().context("empty transcript file")?;
    let (gi, si, tni, esi, eei, epi, csi, cei) = (
        column(header, "GeneId")?,
        column(header, "Strand")?,
        column(header, "TransName")?,
        column(header, "ExonStart")?,
        column(header, "ExonEnd")?,
        column(header, "ExonPhase")?,
        column(header, "CodingStart")?,
        column(header, "CodingEnd")?,
    );
    let mut seen: BTreeSet<String> = BTreeSet::new();
    for line in &lines[1..] {
        let values: Vec<&str> = line.split(',').collect();
        if values.len() <= gi.max(si).max(tni).max(esi).max(eei).max(epi).max(csi).max(cei) {
            continue;
        }
        let Some(t) = out.get_mut(values[tni]) else { continue };
        if t.gene_id != values[gi] {
            continue;
        }
        let parse_opt = |s: &str| -> Option<i32> { if s.eq_ignore_ascii_case("NULL") { None } else { s.parse().ok() } };
        t.pos_strand = values[si].parse::<i32>()? == 1;
        t.coding_start = parse_opt(values[csi]);
        t.coding_end = parse_opt(values[cei]);
        t.exons.push((values[esi].parse()?, values[eei].parse()?, values[epi].parse()?));
        seen.insert(t.trans_name.clone());
    }
    // a transcript without exon data is one EnsemblDataCache does not hold
    out.retain(|name, _| seen.contains(name));
    Ok(out)
}

/// A read as the evidence stage saw it, kept once however many variants show it
pub struct VisRead {
    pub read: Read,
    /// the two mates a fragment read was formed from
    pub fragment: Option<(Read, Read)>,
}

thread_local! {
    static CURRENT_READ: RefCell<Option<(usize, Arc<VisRead>)>> = const { RefCell::new(None) };
}

/// Called before each read of the evidence stage when the visualiser is on
pub fn begin_read() {
    CURRENT_READ.with(|c| *c.borrow_mut() = None);
}

fn shared_read(read: &Read, fragment: Option<&FragmentData>) -> Arc<VisRead> {
    let key = read as *const Read as usize;
    CURRENT_READ.with(|c| {
        let mut slot = c.borrow_mut();
        if let Some((k, shared)) = slot.as_ref() {
            if *k == key && shared.read.name == read.name && shared.read.start == read.start && shared.read.flags == read.flags {
                return shared.clone();
            }
        }
        let shared = Arc::new(VisRead { read: read.clone(), fragment: fragment.map(|f| (f.first.clone(), f.second.clone())) });
        *slot = Some((key, shared.clone()));
        shared
    })
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Shown {
    /// the read itself, which has no fragment
    Read,
    /// one mate of a fragment, the other being out of view
    First,
    Second,
    /// the fragment read with both mates
    Fragment,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orientation {
    Forward,
    Reverse,
    Mixed,
}

impl Orientation {
    fn sort_key(self) -> i32 {
        match self {
            Orientation::Forward => -1,
            Orientation::Reverse => 1,
            Orientation::Mixed => 0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct VisQualities {
    pub combined_quality: f64,
    pub final_base_quality: f64,
    pub final_map_quality: i32,
    pub seq_tech_base_quality: f64,
}

/// ReadEvidenceRecord
#[derive(Clone)]
pub struct ReadEvidenceRecord {
    shared: Arc<VisRead>,
    shown: Shown,
    pub qualities: Option<VisQualities>,
    variant_position: i32,
}

fn orientation_of(read: &Read) -> Orientation {
    if read.negative_strand() { Orientation::Reverse } else { Orientation::Forward }
}

impl ReadEvidenceRecord {
    /// ReadEvidenceRecord.Read
    pub fn read(&self) -> &Read {
        match (self.shown, self.shared.fragment.as_ref()) {
            (Shown::First, Some(f)) => &f.0,
            (Shown::Second, Some(f)) => &f.1,
            _ => &self.shared.read,
        }
    }
    /// ReadEvidenceRecord.Fragment
    pub fn fragment(&self) -> Option<&(Read, Read)> {
        if self.shown == Shown::Fragment { self.shared.fragment.as_ref() } else { None }
    }
    pub fn orientation(&self) -> Orientation {
        let Some((first, second)) = self.fragment() else { return orientation_of(self.read()) };
        let first_overlaps = position_within(self.variant_position, first.start, first.end);
        let second_overlaps = position_within(self.variant_position, second.start, second.end);
        let (first_orientation, second_orientation) = (orientation_of(first), orientation_of(second));
        match (first_overlaps, second_overlaps) {
            (true, false) => first_orientation,
            (false, true) => second_orientation,
            (true, true) => {
                if first_orientation == second_orientation { first_orientation } else { Orientation::Mixed }
            }
            (false, false) => {
                if first.flag(FLAG_FIRST) { first_orientation } else { second_orientation }
            }
        }
    }
    /// compareTo
    fn compare(&self, other: &ReadEvidenceRecord) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        let by_orientation = self.orientation().sort_key().cmp(&other.orientation().sort_key());
        if by_orientation != Equal {
            return by_orientation;
        }
        let by_map_quality = || other.read().mapq.cmp(&self.read().mapq);
        match (self.qualities.as_ref(), other.qualities.as_ref()) {
            (None, None) => by_map_quality(),
            (None, Some(_)) => Greater,
            (Some(_), None) => Less,
            (Some(a), Some(b)) => {
                if a.combined_quality == b.combined_quality {
                    by_map_quality()
                } else if a.combined_quality > b.combined_quality {
                    Less
                } else {
                    Greater
                }
            }
        }
    }
}

/// VariantVis: what one sample's counter collects for one variant
pub struct VariantVis {
    pub sample_id: String,
    pub chromosome: String,
    pub tier: Tier,
    pub ctx: Arc<VariantReadContext>,
    pub view_region: BaseRegion,
    pub context_borders: Vec<(i32, Vec<(BorderLocation, Color)>)>,
    pub ref_view_model: BaseSeqViewModel,
    pub context_view_model: BaseSeqViewModel,
    pub variant_key: String,
    pub indexed_bases_key: String,
    /// reference sequence names of the sample's alignment file, for the mate's chromosome
    pub sequence_names: Arc<Vec<String>>,
    pub read_count: i32,
    pub records_by_type: [Vec<ReadEvidenceRecord>; MATCH_TYPE_COUNT],
    pub read_count_by_type: [i32; MATCH_TYPE_COUNT],
}

/// VisFileBuilder.contextViewModelFromVariant
pub fn context_view_model_from_variant(ctx: &VariantReadContext) -> BaseSeqViewModel {
    let variant = &ctx.variant;
    let raw = &ctx.read_bases;
    let pos_start = variant.pos - ctx.var_index;
    if variant.ref_len() == variant.alt_len() {
        return BaseSeqViewModel::from_bases(raw, pos_start);
    }
    let var_index = ctx.var_index as usize;
    let mut bases: Vec<BaseViewModel> = raw[..=var_index].iter().map(|b| BaseViewModel::base(*b)).collect();
    if variant.is_delete() {
        let del_len = variant.ref_len() - variant.alt_len();
        bases.extend((0..del_len).map(|_| BaseViewModel::del(false)));
        bases.extend(raw[var_index + 1..].iter().map(|b| BaseViewModel::base(*b)));
    } else {
        let ins_len = (variant.alt_len() - variant.ref_len()) as usize;
        bases.last_mut().unwrap().right_insert_count += ins_len as i32;
        if var_index + ins_len + 1 < raw.len() {
            bases.extend(raw[var_index + ins_len + 1..].iter().map(|b| BaseViewModel::base(*b)));
        }
    }
    BaseSeqViewModel::new(bases, pos_start, None, None)
}

/// the key VARIANT_INDEXED_BASES_MAP holds for a read context
pub fn indexed_bases_key(ctx: &VariantReadContext) -> String {
    format!(
        "{}_{}_{}_{}_{}",
        ctx.var_index,
        ctx.core_index_start,
        ctx.core_index_end,
        ctx.left_flank_length(),
        String::from_utf8_lossy(&ctx.read_bases)
    )
}

pub fn variant_key(chr: &str, variant: &SimpleVariant) -> String {
    format!("{}_{}_{}_{}", chr, variant.pos, variant.ref_str(), variant.alt_str())
}

fn register_indexed_bases(variant_key: &str, indexed_bases_key: &str) {
    let mut map = VARIANT_INDEXED_BASES.lock().unwrap();
    map.get_or_insert_with(HashMap::new).entry(variant_key.to_string()).or_default().insert(indexed_bases_key.to_string());
}

impl VariantVis {
    /// new VariantVis(config, sample, variant, readContext, tier). None for a chromosome
    /// RefGenomeCoordinates does not list, where Java has no reference bases to show.
    pub fn new(
        chr: &str, sample_id: &str, ctx: Arc<VariantReadContext>, tier: Tier, genome: &RefGenome, sequence_names: Arc<Vec<String>>,
    ) -> Result<Option<Box<VariantVis>>> {
        if !is_human_chromosome(chr) {
            return Ok(None);
        }
        let variant = &ctx.variant;
        let core_start = ctx.core_position_start;
        let core_end = ctx.core_position_end;
        let flank_start = core_start - ctx.left_flank_length();
        let flank_end = core_end + ctx.right_flank_length();
        let view_region = BaseRegion::new(core_start - READ_EXTEND_LENGTH, core_end + READ_EXTEND_LENGTH);
        let context_borders = vec![
            (variant.pos, vec![(BorderLocation::Left, BLACK), (BorderLocation::Right, BLACK)]),
            (core_start, vec![(BorderLocation::Left, BLUE)]),
            (core_end, vec![(BorderLocation::Right, BLUE)]),
            (flank_start, vec![(BorderLocation::Left, DARK_BLUE)]),
            (flank_end, vec![(BorderLocation::Right, DARK_BLUE)]),
        ];
        let ref_pos_start = view_region.start.max(1);
        let ref_pos_end = view_region.end.min(genome.chromosome_length(chr));
        let ref_bases = genome.bases(chr, ref_pos_start, ref_pos_end)?;
        let variant_key = variant_key(chr, variant);
        let indexed_bases_key = indexed_bases_key(&ctx);
        register_indexed_bases(&variant_key, &indexed_bases_key);
        Ok(Some(Box::new(VariantVis {
            sample_id: sample_id.to_string(),
            chromosome: chr.to_string(),
            tier,
            view_region,
            context_borders,
            ref_view_model: BaseSeqViewModel::from_bases(&ref_bases, ref_pos_start),
            context_view_model: context_view_model_from_variant(&ctx),
            variant_key,
            indexed_bases_key,
            sequence_names,
            read_count: 0,
            records_by_type: Default::default(),
            read_count_by_type: [0; MATCH_TYPE_COUNT],
            ctx,
        })))
    }

    /// addEvidence
    pub fn add_evidence(
        &mut self, read: &Read, fragment: Option<&FragmentData>, match_type: MatchType, qualities: Option<&QualityScores>,
    ) {
        let index = match_type as usize;
        self.read_count += 1;
        self.read_count_by_type[index] += 1;
        let shown = match fragment {
            None => Shown::Read,
            Some(f) => {
                let visible =
                    |r: &Read| positions_overlap(self.view_region.start, self.view_region.end, unclipped_start(r), unclipped_end(r));
                match (visible(&f.first), visible(&f.second)) {
                    (true, false) => Shown::First,
                    (false, true) => Shown::Second,
                    _ => Shown::Fragment,
                }
            }
        };
        self.records_by_type[index].push(ReadEvidenceRecord {
            shared: shared_read(read, fragment),
            shown,
            qualities: qualities.map(|q| VisQualities {
                combined_quality: q.combined_quality,
                final_base_quality: q.final_base_quality,
                final_map_quality: q.final_map_quality,
                seq_tech_base_quality: q.seq_tech_base_quality,
            }),
            variant_position: self.ctx.variant.pos,
        });
    }
}

/// The number of reads shown for a type holding `count` (downsampleReadEvidenceRecords)
pub fn max_reads_shown(match_type: usize, max_support_reads: i32) -> usize {
    if max_support_reads < 0 {
        MAX_READ_UPPER_LIMIT
    } else if max_support_reads > 0 {
        (max_support_reads as usize).min(MAX_READ_UPPER_LIMIT)
    } else {
        MAX_READS_PER_TYPE[match_type]
    }
}

/// Which of `count` records are shown when at most `max_reads` may be: all of them, or else ones
/// spread evenly over the order the reads arrived in. Java draws them at random.
pub fn downsample_indices(count: usize, max_reads: usize) -> Vec<usize> {
    if count <= max_reads {
        return (0..count).collect();
    }
    (0..max_reads).map(|k| k * count / max_reads).collect()
}

fn sorted_records(vis: &VariantVis, match_type: usize, max_support_reads: i32) -> Vec<&ReadEvidenceRecord> {
    let records = &vis.records_by_type[match_type];
    let mut shown: Vec<&ReadEvidenceRecord> =
        downsample_indices(records.len(), max_reads_shown(match_type, max_support_reads)).into_iter().map(|i| &records[i]).collect();
    shown.sort_by(|a, b| a.compare(b));
    shown
}

/// SamRecordUtils.getMateAlignmentEnd: 0 without a mate cigar
fn mate_alignment_end(read: &Read) -> i32 {
    let Some(cigar) = read.mate_cigar.as_deref() else { return 0 };
    if cigar == "*" {
        return 0;
    }
    // CigarUtils.getReadBoundaryPosition(readStart, cigar, false, false)
    let mut current = read.mate_start;
    let mut element_length = 0;
    for c in cigar.bytes() {
        if matches!(c, b'D' | b'M' | b'S' | b'N') {
            if c == b'S' && read.mate_start == current {
                // the left clip
            } else if c == b'S' {
                break;
            } else {
                current += element_length;
            }
            element_length = 0;
            continue;
        }
        element_length = if c.is_ascii_digit() { element_length * 10 + (c - b'0') as i32 } else { 0 };
    }
    current - 1
}

fn svg_box_cell(color: Color) -> String {
    tag("td", &[], &render_colored_box(READ_HEIGHT_PX, color).svg_element())
}

/// ReadTableColumn: the seven cells before the read, as (rendered `td` without style, column style)
fn read_table_cells(record: &ReadEvidenceRecord) -> Vec<(String, CssBuilder)> {
    let empty = CssBuilder::empty();
    let mut cells = Vec::with_capacity(7);
    // Mate
    let read = match record.fragment() {
        Some((first, _)) => first,
        None => record.read(),
    };
    let mate = if !read.paired() {
        (String::new(), empty.clone())
    } else if read.flag(FLAG_MATE_UNMAPPED) {
        ("UN".to_string(), empty.clone())
    } else if read.tid != read.mate_tid {
        ("TRL".to_string(), empty.color(PURPLE))
    } else if read.negative_strand() == read.flag(FLAG_MATE_REVERSE) {
        ("INV".to_string(), empty.color(BLUE))
    } else {
        let is_dup = if read.negative_strand() {
            read.end < read.mate_start
        } else {
            let mate_end = mate_alignment_end(read);
            mate_end != 0 && read.start > mate_end
        };
        if is_dup {
            ("DUP".to_string(), empty.color(DARK_GREEN))
        } else if read.insert_size.abs() >= INSERT_SIZE_CUTTOFF {
            ("DEL".to_string(), empty.color(RED))
        } else {
            (String::new(), empty.clone())
        }
    };
    cells.push((td(&mate.0), mate.1));
    // MAPQ
    let shade = 1.0 * record.read().mapq as f64 / MAX_MAPQ_SHADING_CUTTOFF as f64;
    cells.push((svg_box_cell(interpolate_colors(RED, GREEN, shade)), empty.clone()));
    // combinedQ, finalBQ, finalMQ, seqTechBQ
    let int_style = empty.text_align("right").padding_right(CssSize::em(0.5f32 as f64));
    let q = record.qualities.as_ref();
    let values = [
        q.map(|q| jround(q.combined_quality) as i32),
        q.map(|q| jround(q.final_base_quality) as i32),
        q.map(|q| q.final_map_quality),
        q.map(|q| jround(q.seq_tech_base_quality) as i32),
    ];
    for value in values {
        cells.push((td(&value.map(|v| v.to_string()).unwrap_or_default()), int_style.clone()));
    }
    // Orientation
    let color = match record.orientation() {
        Orientation::Forward => FORWARD_STRAND_COLOR,
        Orientation::Reverse => REVERSE_STRAND_COLOR,
        Orientation::Mixed => overlapping_fragment_border_color(),
    };
    cells.push((svg_box_cell(color), empty));
    cells
}

/// re-open a rendered `<td>...` to give it a style attribute (`TdTag.withStyle`)
fn with_style(td: &str, style: &CssBuilder) -> String {
    let rest = td.strip_prefix("<td").expect("td");
    format!("<td style=\"{}\"{}", escape(&style.to_string()), rest)
}

fn render_bases(vis: &VariantVis, bases: &BaseSeqViewModel, shade_quals: bool, compare_to_ref: bool) -> String {
    render_base_seq(
        READ_HEIGHT_PX,
        vis.view_region,
        bases,
        shade_quals,
        &vis.context_borders,
        if compare_to_ref { Some(&vis.ref_view_model) } else { None },
    )
    .svg_element()
}

/// SamRecordUtils.getOrientationString
fn orientation_string(read: &Read) -> String {
    if !read.paired() || read.flag(FLAG_UNMAPPED) || read.flag(FLAG_MATE_UNMAPPED) {
        return String::new();
    }
    let letter = |negative: bool| if negative { "R" } else { "F" };
    let (first, second) = if read.flag(FLAG_FIRST) {
        (letter(read.negative_strand()), letter(read.flag(FLAG_MATE_REVERSE)))
    } else {
        (letter(read.flag(FLAG_MATE_REVERSE)), letter(read.negative_strand()))
    };
    format!("{first}1{second}2")
}

fn info_row(key: &str, value: &str) -> String {
    tag("tr", &[], &(td(key) + &td(value)))
}

/// HtmlUtil.renderReadInfoTable(firstRead, secondRead)
pub fn render_read_info_table(first: &Read, second: Option<&Read>, chromosome: &str, sequence_names: &[String]) -> String {
    let mut rows = String::new();
    rows.push_str(&info_row("Read name:", &first.name));
    let mate_alignment = if first.paired() && !first.flag(FLAG_MATE_UNMAPPED) {
        let mate_chromosome =
            if first.mate_tid < 0 { "*" } else { sequence_names.get(first.mate_tid as usize).map_or("*", |s| s.as_str()) };
        let mate_end = mate_alignment_end(first);
        let mate_end = if mate_end == 0 { "?".to_string() } else { mate_end.to_string() };
        format!("{}:{}-{}", mate_chromosome, first.mate_start, mate_end)
    } else {
        "unmapped".to_string()
    };
    rows.push_str(&info_row("Alignment:", &format!("{}:{}-{}, {}", chromosome, first.start, first.end, mate_alignment)));
    let cigar = if first.cigar.is_empty() { "*".to_string() } else { cigar_str(&first.cigar) };
    rows.push_str(&info_row("Cigar:", &format!("{}, {}", cigar, first.mate_cigar.as_deref().unwrap_or("missing"))));
    let insert_size = if first.paired() { first.insert_size } else { first.end - first.start + 1 };
    rows.push_str(&info_row("Insert size:", &insert_size.abs().to_string()));
    rows.push_str(&info_row("Orientation:", &orientation_string(first)));
    let both = |a: i32, b: Option<i32>| match b {
        Some(b) => format!("{a}, {b}"),
        None => a.to_string(),
    };
    rows.push_str(&info_row("MapQ:", &both(first.mapq, second.map(|s| s.mapq))));
    rows.push_str(&info_row("NM:", &both(first.nm.unwrap_or(0), second.map(|s| s.nm.unwrap_or(0)))));
    if first.has_consensus_tag {
        let name = match first.consensus {
            ConsensusType::None => "NONE",
            ConsensusType::Single => "SINGLE",
            ConsensusType::Dual => "DUAL",
        };
        rows.push_str(&info_row("Dup type:", name));
    }
    let dup_count = first.consensus_read_tag.as_deref().map_or("0", |tag| tag.split(';').next().unwrap_or(""));
    rows.push_str(&info_row("Dup count:", dup_count));
    let style = CssBuilder::empty().padding(CssSize::ZERO).margin(CssSize::ZERO).display("none");
    tag("div", &[("class", "read-info"), ("style", &style.to_string())], &styled_table(&rows, &CssBuilder::empty()))
}

/// RenderUtils.renderRead
fn render_read(vis: &VariantVis, record: &ReadEvidenceRecord) -> String {
    let base_div_style = CssBuilder::empty().padding(CssSize::ZERO).margin(CssSize::ZERO);
    let base_div = base_div_style.to_string();
    let mut inner = String::new();
    match record.fragment() {
        None => {
            let model = BaseSeqViewModel::from_read(record.read());
            push_tag(&mut inner, "div", &[("class", "read-svg"), ("style", &base_div)], &render_bases(vis, &model, true, true));
            inner.push_str(&render_read_info_table(record.read(), None, &vis.chromosome, &vis.sequence_names));
        }
        Some((first, second)) => {
            let first_model = BaseSeqViewModel::from_read(first);
            let second_model = BaseSeqViewModel::from_read(second);
            let model = BaseSeqViewModel::from_consensus_fragment(&record.shared.read, Some(&first_model), Some(&second_model));
            push_tag(&mut inner, "div", &[("class", "read-svg"), ("style", &base_div)], &render_bases(vis, &model, true, true));
            let hidden = base_div_style.display("none").to_string();
            for mate in [&first_model, &second_model] {
                push_tag(&mut inner, "div", &[("class", "read-of-fragment-sgv"), ("style", &hidden)], &render_bases(vis, mate, true, true));
            }
            let (first_read, second_read) = if first.flag(FLAG_FIRST) { (first, second) } else { (second, first) };
            inner.push_str(&render_read_info_table(first_read, Some(second_read), &vis.chromosome, &vis.sequence_names));
        }
    }
    tag("div", &[("style", &base_div)], &inner)
}

/// AminoAcidElements
pub struct AminoAcidElements {
    pub gene_region_label: String,
    pub ref_svg: String,
    pub alt_svg: String,
    pub impact: VariantImpact,
}

/// RenderUtils.renderReads: the table rows of one sample
fn render_reads(vis: &VariantVis, is_tumor: bool, aa: Option<&AminoAcidElements>, max_support_reads: i32, rows: &mut String) {
    let empty = CssBuilder::empty();
    let light_gray_bg = empty.background_color(LIGHT_GRAY);
    let vertical_header = light_gray_bg.writing_mode("vertical-rl");
    let header_style = light_gray_bg.text_align("right").to_string();
    let match_type_border = empty.border_bottom(CssSize::px(2.0), "solid", BLACK);
    let match_type_style = vertical_header.merge(&match_type_border).text_align("center").padding(CssSize::px(4.0));
    let font_pt = jround(2.0 / 3.0 * READ_HEIGHT_PX) as i32;
    let table_info_cell = empty.font_size_pt(font_pt).font_weight("bold");
    let spacer_div = empty.height(CssSize::em(1.0)).padding(CssSize::ZERO).margin(CssSize::ZERO);
    const COLUMNS: [&str; 7] = ["Mate", "MAPQ", "combinedQ", "finalBQ", "finalMQ", "seqTechBQ", "Orientation"];
    let full_span = (COLUMNS.len() + 2).to_string();
    let label_span = (COLUMNS.len() + 1).to_string();

    let spacer = tag("div", &[("style", &spacer_div.to_string())], "");
    push_tag(rows, "tr", &[], &tag("td", &[("colspan", &full_span)], &spacer));
    let sample = format!("{} {}", vis.sample_id, if is_tumor { "(tumor)" } else { "(reference)" });
    push_tag(rows, "tr", &[], &tag("td", &[("colspan", &full_span)], &escape(&sample)));

    let mut header = td_styled("Type", &vertical_header);
    for column in COLUMNS {
        header.push_str(&td_styled(column, &vertical_header));
    }
    let coords = render_coords(READ_HEIGHT_PX, vis.view_region, vis.ctx.variant.pos, DISPLAY_EVERY_NTH_COORD).svg_element();
    header.push_str(&td_styled(&coords, &light_gray_bg));
    push_tag(rows, "tr", &[], &header);

    let labelled_row = |rows: &mut String, label: &str, cell: &str| {
        let label_cell = tag("td", &[("colspan", &label_span), ("style", &header_style)], label);
        push_tag(rows, "tr", &[], &(label_cell + cell));
    };
    if let Some(aa) = aa {
        labelled_row(rows, "ref", &tag("td", &[], &aa.ref_svg));
        labelled_row(rows, "predicted", &tag("td", &[], &aa.alt_svg));
        let gene_name_style = empty.font_size_pt(font_pt).font_weight("bold").text_align("center");
        labelled_row(rows, "gene", &td_styled(&escape(&aa.gene_region_label), &gene_name_style));
    }
    labelled_row(rows, "ref", &tag("td", &[], &render_bases(vis, &vis.ref_view_model, false, false)));
    labelled_row(rows, "context", &tag("td", &[], &render_bases(vis, &vis.context_view_model, false, true)));

    for match_type in SORTED_MATCH_TYPES {
        let records = sorted_records(vis, match_type, max_support_reads);
        if records.is_empty() {
            continue;
        }
        let type_content = format!("{}<br>({}/{})", MATCH_TYPE_NAMES[match_type], records.len(), vis.read_count_by_type[match_type]);
        for (i, record) in records.iter().enumerate() {
            let is_last = i == records.len() - 1;
            let mut cols = String::new();
            if i == 0 {
                push_tag(
                    &mut cols,
                    "td",
                    &[("rowspan", &records.len().to_string()), ("style", &match_type_style.to_string())],
                    &type_content,
                );
            }
            for (cell, cell_style) in read_table_cells(record) {
                let mut style = table_info_cell.merge(&cell_style);
                if is_last {
                    style = style.merge(&match_type_border);
                }
                cols.push_str(&with_style(&cell, &style));
            }
            let read = render_read(vis, record);
            if is_last {
                cols.push_str(&td_styled(&read, &match_type_border));
            } else {
                push_tag(&mut cols, "td", &[], &read);
            }
            push_tag(rows, "tr", &[], &cols);
        }
    }
}

/// VisFileBuilder.getSampleInfoRowElements: the cell texts of a row
fn info_cells(key: &str, tumor: Option<String>, reference: Option<String>) -> Vec<String> {
    let mut cells = vec![key.to_string()];
    cells.extend(tumor);
    cells.extend(reference);
    cells
}

fn first_with_vis(counters: &[ReadContextCounter]) -> Option<&ReadContextCounter> {
    counters.iter().find(|c| c.vis.is_some())
}

/// VisFileBuilder.renderSampleSummaryAndCountsTable
fn render_sample_summary_and_counts_table(v: &SageVariant) -> String {
    let empty = CssBuilder::empty();
    let border_style = empty.border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
    let table_style = border_style.margin_right(CssSize::px(10.0));
    let header_style = empty.font_weight("bold").text_align("center").background_color(LIGHT_GRAY).to_string();
    let cell_style = empty.padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5));
    let tumor = first_with_vis(&v.tumor_counters);
    let reference = first_with_vis(&v.ref_counters);
    let col_count = 1 + tumor.is_some() as usize + reference.is_some() as usize;
    let both = |f: &dyn Fn(&ReadContextCounter) -> String| (tumor.map(f), reference.map(f));

    let mut content: Vec<Vec<String>> = Vec::new();
    let (t, r) = both(&|c| c.alt_support().to_string());
    content.push(info_cells("AD", t, r));
    let (t, r) = both(&|c| c.depth().to_string());
    content.push(info_cells("DP", t, r));
    let (t, r) = both(&|c| java_fixed(c.vaf(), 3));
    content.push(info_cells("AF", t, r));
    for match_type in SORTED_MATCH_TYPES {
        let (t, r) = both(&|c| c.vis.as_ref().unwrap().read_count_by_type[match_type].to_string());
        content.push(info_cells(MATCH_TYPE_NAMES[match_type], t, r));
    }

    let bordered_cell = cell_style.merge(&border_style);
    let value_cell = empty.text_align("right").merge(&bordered_cell);
    let mut rows = String::new();
    let header: String = info_cells("Info", tumor.map(|_| "Tumor".to_string()), reference.map(|_| "Normal".to_string()))
        .iter()
        .map(|c| td_styled(&escape(c), &bordered_cell))
        .collect();
    push_tag(&mut rows, "tr", &[("style", &header_style)], &header);
    for (i, cells) in content.iter().enumerate() {
        if i == 3 {
            let spacer = tag(
                "td",
                &[("style", &cell_style.no_border().to_string()), ("colspan", &col_count.to_string())],
                "Support By Type",
            );
            push_tag(&mut rows, "tr", &[("style", &header_style)], &spacer);
        }
        let mut row = td_styled(&escape(&cells[0]), &bordered_cell);
        for cell in &cells[1..] {
            row.push_str(&td_styled(&escape(cell), &value_cell));
        }
        push_tag(&mut rows, "tr", &[], &row);
    }
    tag("div", &[], &styled_table(&rows, &table_style))
}

/// VisFileBuilder.renderSampleQualAndSiteInfo
fn render_sample_qual_and_site_info(v: &SageVariant, first_counter: &ReadContextCounter) -> String {
    let empty = CssBuilder::empty();
    let border_style = empty.border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
    let table_style = border_style.margin_right(CssSize::px(10.0));
    let header_style = empty.font_weight("bold").text_align("center").background_color(LIGHT_GRAY).to_string();
    let cell_style = empty.padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&border_style);
    let tumor = first_with_vis(&v.tumor_counters);
    let reference = first_with_vis(&v.ref_counters);
    let both = |f: &dyn Fn(&ReadContextCounter) -> String| (tumor.map(f), reference.map(f));

    let mut content: Vec<Vec<String>> = Vec::new();
    let mut add = |key: &str, values: (Option<String>, Option<String>)| content.push(info_cells(key, values.0, values.1));
    add("RAW_QUAL", both(&|c| java_double(c.tumor_quality())));
    add("ARCBQ", both(&|c| (c.average_alt_recalibrated_base_quality() as i32).to_string()));
    add("ASBQ", both(&|c| (c.average_alt_seq_tech_base_quality() as i32).to_string()));
    add(
        "ARMQ",
        both(&|c| {
            let alt_support = c.alt_support();
            let avg = if alt_support > 0 { jround(c.qual_counters.alt_map_quality_total as f64 / alt_support as f64) as i32 } else { 0 };
            avg.to_string()
        }),
    );
    add("FSB", both(&|c| java_fixed(c.alt_fragment_strand_bias.bias(), 2)));
    add("RSB", both(&|c| java_fixed(c.alt_read_strand_bias.bias(), 2)));
    add("JIT", both(&|c| format!("{}-{}", c.jitter.shortened, c.jitter.lengthened)));
    add("MUC", both(&|c| c.fragment_coords.min_count().to_string()));
    add(
        "UMI_CNT",
        both(&|c| match c.consensus_type_counts {
            Some(counts) => counts.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "),
            None => "null".to_string(),
        }),
    );
    // Java shows the first counter's two averages, one under each sample
    add(
        "AED",
        (
            tumor.map(|_| java_fixed(first_counter.read_edge_distance.avg_distance_from_edge(), 2)),
            reference.map(|_| java_fixed(first_counter.read_edge_distance.avg_alt_distance_from_edge(), 2)),
        ),
    );

    let value_cell = cell_style.text_align("right");
    let mut rows = String::new();
    let header: String = info_cells("Info", tumor.map(|_| "Tumor".to_string()), reference.map(|_| "Normal".to_string()))
        .iter()
        .map(|c| td_styled(&escape(c), &cell_style))
        .collect();
    push_tag(&mut rows, "tr", &[("style", &header_style)], &header);
    for cells in &content {
        let mut row = td_styled(&escape(&cells[0]), &cell_style);
        for cell in &cells[1..] {
            row.push_str(&td_styled(&escape(cell), &value_cell));
        }
        push_tag(&mut rows, "tr", &[], &row);
    }
    tag("div", &[], &styled_table(&rows, &table_style))
}

/// The filter names in the order Java's `HashSet<String>` iterates them. Names that share a bucket
/// come out in the order of the variant's own filter set, which in Java is a `HashSet` of enum
/// values ordered by identity hash; here they keep the order the filters are declared in.
pub fn filter_names_in_java_set_order(names: &[&str]) -> Vec<String> {
    let mut capacity = 16usize;
    while names.len() as f64 > capacity as f64 * 0.75 {
        capacity *= 2;
    }
    let mut keyed: Vec<(usize, usize, &str)> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let h = crate::sync::java_string_hash(name);
            let spread = (h ^ ((h as u32) >> 16) as i32) as u32 as usize;
            (spread & (capacity - 1), i, *name)
        })
        .collect();
    keyed.sort();
    let mut out: Vec<String> = Vec::new();
    for (_, _, name) in keyed {
        if !out.iter().any(|n| n == name) {
            out.push(name.to_string());
        }
    }
    out
}

/// RenderUtils.renderVariantInfoTable
fn render_variant_info_table(vis: &VariantVis, v: &SageVariant, first_counter: &ReadContextCounter, aa: Option<&AminoAcidElements>) -> String {
    let empty = CssBuilder::empty();
    let border_style = empty.border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
    let header_style = empty.font_weight("bold").text_align("center").background_color(LIGHT_GRAY);
    let cell_style = empty.padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&border_style);
    let span = |text: &str| tag("span", &[], &escape(text));
    let mut values: Vec<(&str, String)> = Vec::new();
    if let Some(aa) = aa {
        values.push(("Gene", span(&aa.impact.gene_name)));
        values.push(("Transcript", span(&aa.impact.canonical_transcript)));
        values.push(("HGVS", span(&aa.impact.canonical_hgvs_protein)));
        values.push(("Coding impact", span(&aa.impact.canonical_effect)));
    }
    let filter_str = if v.filters.is_empty() {
        "PASS".to_string()
    } else {
        let declared: Vec<&str> = crate::filters::ALL_FILTERS.iter().filter(|(f, _, _)| v.filters.contains(f)).map(|(_, name, _)| *name).collect();
        filter_names_in_java_set_order(&declared).join(",")
    };
    let ctx = &vis.ctx;
    let context = span(ctx.left_flank_str())
        + &tag("span", &[("style", &empty.font_weight("bold").to_string())], &escape(ctx.core_str()))
        + &span(ctx.right_flank_str());
    let repeat = match ctx.max_repeat.as_ref() {
        Some(r) => format!("{}x{}", r.count, r.bases_str()),
        None => "NO REPEAT".to_string(),
    };
    // (int)(-10 * logTqp) and round(round(mapQualFactor * 10) / 10)
    let tqp = jmin(jmax(first_counter.tumor_qual_probability, 0.0), 1.0);
    let total_tumor_quality = (-10.0 * jmax(tqp, 1e-20).log10()) as i32;
    let map_qual_factor = jround(jround(first_counter.map_qual_factor * 10.0) as f64 / 10.0) as f64;
    values.push(("Qual", span(&total_tumor_quality.to_string())));
    values.push(("Filter", span(&filter_str)));
    values.push(("Tier", span(vis.tier.name())));
    values.push(("CONTEXT", context));
    values.push(("NEARBY_INDEL", span(if v.near_indel { "true" } else { "false" })));
    values.push(("REPEAT", span(&repeat)));
    values.push(("MQF", span(&java_double(map_qual_factor))));

    let key_style = cell_style.merge(&header_style);
    let mut rows = String::new();
    for (key, value) in values {
        push_tag(&mut rows, "tr", &[], &(td_styled(key, &key_style) + &td_styled(&value, &cell_style)));
    }
    tag("div", &[], &styled_table(&rows, &border_style))
}

/// VisFileBuilder.getAminoAcidsElements
fn amino_acid_elements(
    vis_cfg: &VisConfig, chr: &str, view_region: BaseRegion, variant: &SimpleVariant, genome: &RefGenome,
) -> Result<Option<AminoAcidElements>, String> {
    if vis_cfg.purple_vcf.is_none() {
        return Ok(None);
    }
    let Some(reportable) = vis_cfg.reportable(chr, variant) else { return Ok(None) };
    let impact = impact_from_values(&reportable.impact_values)?;
    let hgvs = impact.canonical_hgvs_protein.as_str();
    let events: Vec<AminoAcidEvent> =
        if hgvs == "p.?" || hgvs == "unknown" || hgvs.is_empty() { Vec::new() } else { vis_aa::parse_amino_acid_events(hgvs)? };
    let transcript = vis_cfg
        .transcripts
        .get(&impact.canonical_transcript)
        .ok_or_else(|| format!("no amino acid or transcript data for transcript({})", impact.canonical_transcript))?;
    let gene_regions = vis_aa::gene_regions(transcript)?;
    let ref_models = vis_aa::ref_gene_region_view_models(transcript, &gene_regions);
    let ref_bases = |start: i32, end: i32| genome.bases(chr, start, end).map_err(|e| e.to_string());
    let alt_context = vis_aa::AltContext { variant, chromosome_length: genome.chromosome_length(chr), ref_bases: &ref_bases };
    let alt_models = vis_aa::alt_gene_region_view_models(transcript, &gene_regions, &events, view_region, &alt_context)?;
    let (ref_canvas, alt_canvas) = render_gene_data(READ_HEIGHT_PX, view_region, transcript.pos_strand, &ref_models, &alt_models);
    Ok(Some(AminoAcidElements {
        gene_region_label: format!("{} {}", transcript.gene_name, vis_aa::gene_region_label(transcript, variant.pos)),
        ref_svg: ref_canvas.svg_element(),
        alt_svg: alt_canvas.svg_element(),
        impact,
    }))
}

/// VisFileBuilder.formFilename, given the read contexts registered for the variant
pub fn form_filename(
    chr: &str, variant: &SimpleVariant, sample_id: &str, gene_name: Option<&str>, variant_type: Option<&str>,
    indexed_bases_keys: &BTreeSet<String>, indexed_bases_key: &str,
) -> String {
    let variant_prefix = match gene_name {
        // MiscToolFiles.generateSageVisFilePrefix
        Some(gene) => {
            let chr_prefixed = if chr.starts_with("chr") { chr.to_string() } else { format!("chr{chr}") };
            format!(".sage.{}_{}_{}_{}_{}", gene, chr_prefixed, variant.pos, variant.ref_str(), variant.alt_str())
        }
        None => variant_key(chr, variant),
    };
    let type_suffix = match variant_type {
        Some(t) => {
            let sanitised: String = t.replace('&', "_").chars().filter(|c| !matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')).collect();
            format!("_{}", sanitised.replace("_variant", ""))
        }
        None => String::new(),
    };
    let filename = format!("{sample_id}{variant_prefix}");
    if indexed_bases_keys.len() <= 1 {
        return format!("{filename}{type_suffix}.html.gz");
    }
    let file_num = indexed_bases_keys.iter().take_while(|k| k.as_str() < indexed_bases_key).count() + 1;
    let width = indexed_bases_keys.len().to_string().len();
    format!("{filename}_{file_num:0width$}{type_suffix}.html.gz")
}

fn javascript() -> String {
    include_str!("../resources/sagevis.js").lines().collect::<Vec<_>>().join("\n")
}

/// The file name and the HTML of a variant; None when Java writes no file
fn build_html(
    tumor_ids: &[String], reference_ids: &[String], vis_cfg: &VisConfig, chr: &str, v: &SageVariant, genome: &RefGenome,
) -> Result<Option<(String, String)>, String> {
    if vis_cfg.pass_only && !v.is_passing() {
        return Ok(None);
    }
    let tumor_vis: Vec<&VariantVis> = v.tumor_counters.iter().filter_map(|c| c.vis.as_deref()).collect();
    let ref_vis: Vec<&VariantVis> = v.ref_counters.iter().filter_map(|c| c.vis.as_deref()).collect();
    let Some(first_vis) = tumor_vis.first().or(ref_vis.first()).copied() else { return Ok(None) };
    let first_counter = v.tumor_counters.first().or(v.ref_counters.first()).ok_or("no counters")?;
    let variant = v.variant();
    let aa = amino_acid_elements(vis_cfg, chr, first_vis.view_region, variant, genome)?;
    let (gene_name, variant_type) = match aa.as_ref() {
        Some(aa) => (Some(aa.impact.gene_name.as_str()), Some(aa.impact.canonical_effect.as_str())),
        None => (None, None),
    };
    let sample_id = tumor_ids.first().or(reference_ids.first()).ok_or("no samples")?;
    let filename = {
        let map = VARIANT_INDEXED_BASES.lock().unwrap();
        let keys = map.as_ref().and_then(|m| m.get(&first_vis.variant_key)).ok_or("variant has no indexed bases")?;
        form_filename(chr, variant, sample_id, gene_name, variant_type, keys, &first_vis.indexed_bases_key)
    };
    let total_read_count: i32 = tumor_vis.iter().chain(ref_vis.iter()).map(|x| x.read_count).sum();
    if total_read_count == 0 {
        return Ok(None);
    }

    let mut read_rows = String::new();
    for vis in &tumor_vis {
        render_reads(vis, true, aa.as_ref(), vis_cfg.max_support_reads, &mut read_rows);
    }
    for vis in &ref_vis {
        render_reads(vis, false, aa.as_ref(), vis_cfg.max_support_reads, &mut read_rows);
    }
    let empty = CssBuilder::empty();
    let read_table = tag("div", &[], &styled_table(&read_rows, &empty.border_spacing(CssSize::ZERO)));
    let vertical_spacer = tag("div", &[("style", &empty.height(CssSize::em(1.0)).to_string())], "");

    // RenderUtils.renderVariantInfo
    let samples = match (tumor_ids.first(), reference_ids.first()) {
        (None, Some(r)) => r.clone(),
        (Some(t), None) => t.clone(),
        (Some(t), Some(r)) => format!("{t}[normal: {r}]"),
        (None, None) => String::new(),
    };
    let info = format!("{} {}:{} {}>{}", samples, chr, variant.pos, variant.ref_str(), variant.alt_str());
    let variant_info = tag("div", &[("style", &empty.font_weight("bold").to_string())], &escape(&info));

    let top = empty.vertical_align("top");
    let sample_info_row = td_styled(&render_sample_summary_and_counts_table(v), &top)
        + &td_styled(&render_sample_qual_and_site_info(v, first_counter), &top)
        + &td_styled(&render_variant_info_table(first_vis, v, first_counter, aa.as_ref()), &top);
    let sample_info = styled_table(&tag("tr", &[], &sample_info_row), &empty);

    let script = tag("script", &[("type", "text/javascript")], &javascript());
    let body = tag(
        "body",
        &[("style", &base_font_style().to_string())],
        &(variant_info + &vertical_spacer + &sample_info + &read_table + &script),
    );
    Ok(Some((filename, tag("html", &[], &body))))
}

/// VisFileBuilder.writeToHtmlFile: a failure is reported and the variant skipped, as in Java
pub fn write_to_html_file(
    tumor_ids: &[String], reference_ids: &[String], vis_cfg: &VisConfig, chr: &str, v: &SageVariant, genome: &RefGenome,
) {
    let describe = || format!("{}:{} {}>{}", chr, v.position(), v.variant().ref_str(), v.variant().alt_str());
    let (filename, html) = match build_html(tumor_ids, reference_ids, vis_cfg, chr, v, genome) {
        Ok(Some(x)) => x,
        Ok(None) => return,
        Err(e) => {
            eprintln!("[WARN ] failed to generate variant({}) visualisation:: {}", describe(), e);
            return;
        }
    };
    let path = std::path::Path::new(&vis_cfg.output_dir).join(filename);
    let write = || -> std::io::Result<()> {
        let file = std::io::BufWriter::new(std::fs::File::create(&path)?);
        let mut w = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        w.write_all(html.as_bytes())?;
        w.write_all(b"\n")?;
        w.finish()?.flush()
    };
    if let Err(e) = write() {
        eprintln!("[WARN ] failed to generate variant({}) visualisation:: {}", describe(), e);
    }
}
