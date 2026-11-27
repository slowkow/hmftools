//! IsofoxConfig, FusionConfig and GeneRegionFilters: the command line as Java ISOFOX 2.1 registers and reads it,
//! or ISOFOX 2.0.1 in its mode (version.rs).

use anyhow::{anyhow, bail, Context, Result};
use hmf_common::config_builder::ConfigBuilder;
use hmf_common::driver_panel::{read_driver_genes, DriverGene};
use hmf_common::known_fusion::KnownFusionCache;
use hmf_common::regions::{ChrBaseRegion, RefGenomeVersion, SpecificRegions};
use hmf_common::segmentation::HumanChromosome;
use hmf_common::{info, warn};
use std::collections::HashMap;

pub const APP_NAME: &str = "Isofox";

pub const DEFAULT_MAX_FRAGMENT_SIZE: i32 = 550;
pub const DEFAULT_FRAG_LENGTH_MIN_COUNT: i32 = 1_000_000;
pub const DEFAULT_SINGLE_MAP_QUALITY: i32 = 60;
/// ISOFOX 2.0.1's default: STAR gives uniquely mapped reads a map quality of 255
pub const DEFAULT_SINGLE_MAP_QUALITY_201: i32 = 255;
/// ISOFOX 2.0.1: multi-mapped fragments are given map quals of 3 or lower
pub const MULTI_MAP_QUALITY_THRESHOLD: u8 = 3;
/// ISOFOX 2.0.1: the buffer around each enriched gene excluded from the BAM slices
pub const ENRICHED_GENE_BUFFER: i32 = 1000;
pub const DEFAULT_EXPECTED_RATE_LENGTHS: &str = "50-1;75-1;100-1;125-1;150-1;200-1;250-1;300-1;400-1;550-1";
pub const DEFAULT_HARD_FILTER_MIN_FRAGS: i32 = 2;

pub const MAX_NOVEL_SJ_DISTANCE: i32 = 500_000;
pub const MAX_GENE_PERC_CONTRIBUTION: f64 = 0.01;
pub const PANEL_LOW_COVERAGE_FACTOR: f64 = 5.0;
pub const HIGH_EXPRESSION_FOLD_CHANGE_LOW: f64 = 1.0;
pub const HIGH_EXPRESSION_FOLD_CHANGE_HIGH: f64 = 2.3;
pub const HIGH_EXPRESSION_PERCENTILE_HIGH: f64 = 0.98;

/// chromosomes processed first, then the rest in HumanChromosome order
pub const PRIORITISED_CHROMOSOMES: [&str; 8] = ["14", "22", "2", "3", "6", "9", "16", "17"];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IsofoxFunction {
    TranscriptCounts,
    AltSpliceJunctions,
    RetainedIntrons,
    Fusions,
    Statistics,
    NeoEpitopes,
    /// ISOFOX 2.0.1 only (BamReadCounter); not ported
    ReadCounts,
}

impl IsofoxFunction {
    pub fn name(self) -> &'static str {
        match self {
            IsofoxFunction::TranscriptCounts => "TRANSCRIPT_COUNTS",
            IsofoxFunction::AltSpliceJunctions => "ALT_SPLICE_JUNCTIONS",
            IsofoxFunction::RetainedIntrons => "RETAINED_INTRONS",
            IsofoxFunction::Fusions => "FUSIONS",
            IsofoxFunction::Statistics => "STATISTICS",
            IsofoxFunction::NeoEpitopes => "NEO_EPITOPES",
            IsofoxFunction::ReadCounts => "READ_COUNTS",
        }
    }
    pub fn value_of(s: &str) -> Result<Self> {
        if s == "READ_COUNTS" && crate::version::v201() {
            return Ok(IsofoxFunction::ReadCounts);
        }
        Ok(match s {
            "TRANSCRIPT_COUNTS" => IsofoxFunction::TranscriptCounts,
            "ALT_SPLICE_JUNCTIONS" => IsofoxFunction::AltSpliceJunctions,
            "RETAINED_INTRONS" => IsofoxFunction::RetainedIntrons,
            "FUSIONS" => IsofoxFunction::Fusions,
            "STATISTICS" => IsofoxFunction::Statistics,
            "NEO_EPITOPES" => IsofoxFunction::NeoEpitopes,
            _ => bail!("No enum constant com.hartwig.hmftools.isofox.IsofoxFunction.{s}"),
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WriteType {
    Exon,
    SpliceJunc,
    FragLength,
    FragLengthByGene,
    Read,
    FusionFragment,
    ChimericPositionData,
    SpliceSite,
    TransCombo,
    GcRatio,
    MultiMapLoci,
    /// ISOFOX 2.0.1 only
    ChimericRead,
}

impl WriteType {
    pub fn value_of(s: &str) -> Result<Self> {
        if crate::version::v201() {
            // ISOFOX 2.0.1's WriteType: CHIMERIC_READ in place of FUSION_FRAGMENT, and no MULTI_MAP_LOCI
            match s {
                "CHIMERIC_READ" => return Ok(WriteType::ChimericRead),
                "FUSION_FRAGMENT" | "MULTI_MAP_LOCI" => bail!("No enum constant com.hartwig.hmftools.isofox.WriteType.{s}"),
                _ => {}
            }
        }
        Ok(match s {
            "EXON" => WriteType::Exon,
            "SPLICE_JUNC" => WriteType::SpliceJunc,
            "FRAG_LENGTH" => WriteType::FragLength,
            "FRAG_LENGTH_BY_GENE" => WriteType::FragLengthByGene,
            "READ" => WriteType::Read,
            "FUSION_FRAGMENT" => WriteType::FusionFragment,
            "CHIMERIC_POSITION_DATA" => WriteType::ChimericPositionData,
            "SPLICE_SITE" => WriteType::SpliceSite,
            "TRANS_COMBO" => WriteType::TransCombo,
            "GC_RATIO" => WriteType::GcRatio,
            "MULTI_MAP_LOCI" => WriteType::MultiMapLoci,
            _ => bail!("No enum constant com.hartwig.hmftools.isofox.WriteType.{s}"),
        })
    }
    pub const ALL: [&'static str; 11] = [
        "EXON",
        "SPLICE_JUNC",
        "FRAG_LENGTH",
        "FRAG_LENGTH_BY_GENE",
        "READ",
        "FUSION_FRAGMENT",
        "CHIMERIC_POSITION_DATA",
        "SPLICE_SITE",
        "TRANS_COMBO",
        "GC_RATIO",
        "MULTI_MAP_LOCI",
    ];
    pub const ALL_201: [&'static str; 10] = [
        "EXON",
        "SPLICE_JUNC",
        "FRAG_LENGTH",
        "FRAG_LENGTH_BY_GENE",
        "READ",
        "CHIMERIC_READ",
        "CHIMERIC_POSITION_DATA",
        "SPLICE_SITE",
        "TRANS_COMBO",
        "GC_RATIO",
    ];
}

/// FragmentSize: a fragment length and its frequency
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FragmentSize {
    pub length: i32,
    pub frequency: i32,
}

/// IG regions and the HLA region (IsofoxConstants.populateImmuneRegions)
pub fn immune_regions(version: RefGenomeVersion) -> Vec<ChrBaseRegion> {
    if version.is_37() {
        vec![
            ChrBaseRegion::new("2", 89156674, 90538397),
            ChrBaseRegion::new("14", 106032614, 107288051),
            ChrBaseRegion::new("22", 22380474, 23265085),
            ChrBaseRegion::new("6", 29690552, 33111102),
        ]
    } else {
        vec![
            ChrBaseRegion::new("chr2", 88857161, 90315836),
            ChrBaseRegion::new("chr14", 105586437, 106879844),
            ChrBaseRegion::new("chr22", 22026076, 22922913),
            ChrBaseRegion::new("chr6", 29722775, 33143325),
        ]
    }
}

/// the RN7SL and RN7SK genes ISOFOX 2.0.1 treats as enriched (IsofoxConstants.populateEnrichedGeneIds)
pub fn default_enriched_gene_ids(version: RefGenomeVersion) -> Vec<String> {
    let ids: [&str; 6] = if version == RefGenomeVersion::V38 {
        ["ENSG00000276168", "ENSG00000274012", "ENSG00000278771", "ENSG00000263740", "ENSG00000283293", "ENSG00000265735"]
    } else {
        ["ENSG00000265150", "ENSG00000258486", "ENSG00000202198", "ENSG00000266037", "ENSG00000263740", "ENSG00000265735"]
    };
    ids.iter().map(|s| s.to_string()).collect()
}

/// GeneRegionFilters
#[derive(Clone, Debug, Default)]
pub struct GeneRegionFilters {
    pub specific_regions: SpecificRegions,
    pub restricted_gene_ids: Vec<String>,
    pub immune_gene_regions: Vec<ChrBaseRegion>,
    has_specific_regions: bool,
    /// ISOFOX 2.0.1: genes counted but not otherwise processed (EnrichedGeneIds)
    pub enriched_gene_ids: Vec<String>,
    /// ISOFOX 2.0.1: regions whose reads are skipped, from `-excluded_regions` plus each enriched gene with a
    /// 1 kb buffer (ExcludedRegions), sorted by start per chromosome
    pub excluded_regions: HashMap<String, Vec<hmf_common::regions::BaseRegion>>,
}

const READ_END_BUFFER: i32 = 150;

impl GeneRegionFilters {
    pub fn has_specific_regions(&self) -> bool {
        self.has_specific_regions
    }
    /// skipRead(chromosome, readStart): with the read end buffer
    pub fn skip_read_start(&self, chromosome: &str, read_start: i32) -> bool {
        self.skip_read(chromosome, read_start, read_start + READ_END_BUFFER, false)
    }
    pub fn skip_read(&self, chromosome: &str, read_start: i32, read_end: i32, is_mate_or_supp: bool) -> bool {
        if !HumanChromosome::contains(chromosome) {
            return true;
        }
        if crate::version::v201() {
            // ISOFOX 2.0.1 applies the specific regions to mates and supplementaries too, and skips reads
            // overlapping an excluded region
            if self.has_specific_regions {
                if self.specific_regions.exclude_chromosome(chromosome) {
                    return true;
                }
                if !self.specific_regions.regions.iter().any(|r| r.overlaps(chromosome, read_start, read_end)) {
                    return true;
                }
            }
            if let Some(ex) = self.excluded_regions.get(chromosome) {
                if ex.iter().any(|x| positions_overlap(x.start, x.end, read_start, read_end)) {
                    return true;
                }
            }
            return false;
        }
        if !is_mate_or_supp && self.has_specific_regions {
            if self.specific_regions.exclude_chromosome(chromosome) {
                return true;
            }
            if !self.specific_regions.regions.iter().any(|r| r.overlaps(chromosome, read_start, read_end)) {
                return true;
            }
        }
        false
    }

    /// findExcludedRegions (2.0.1): the excluded regions overlapping a region, in sorted order
    pub fn find_excluded_regions(&self, chromosome: &str, start: i32, end: i32) -> Vec<hmf_common::regions::BaseRegion> {
        match self.excluded_regions.get(chromosome) {
            Some(l) => l.iter().copied().filter(|x| positions_overlap(x.start, x.end, start, end)).collect(),
            None => Vec::new(),
        }
    }

    /// buildGeneRegions (2.0.1): each enriched gene, with a buffer, becomes an excluded region
    pub fn build_gene_regions(&mut self, cache: &hmf_common::ensembl::EnsemblDataCache) {
        for id in self.enriched_gene_ids.clone() {
            let Some(g) = cache.gene_data_by_id(&id) else {
                warn!("enriched gene ID({}) missing from Ensembl cache", id);
                continue;
            };
            self.excluded_regions
                .entry(g.chromosome.clone())
                .or_default()
                .push(hmf_common::regions::BaseRegion::new(g.gene_start - ENRICHED_GENE_BUFFER, g.gene_end + ENRICHED_GENE_BUFFER));
        }
        // Collections.sort: by start (BaseRegion.compareTo), stable
        for l in self.excluded_regions.values_mut() {
            l.sort_by_key(|r| r.start);
        }
    }
}

fn positions_overlap(s1: i32, e1: i32, s2: i32, e2: i32) -> bool {
    s1 <= e2 && e1 >= s2
}

/// FusionConfig
#[derive(Clone)]
#[allow(dead_code)]
pub struct FusionConfig {
    pub cohort_file: Option<String>,
    pub cache_fragments: bool,
    pub skip_non_genic: bool,
    pub write_chimeric_only: bool,
    pub run_perf_checks: bool,
    pub min_hard_filter_frags: i32,
    pub known_fusions: KnownFusionCache,
}

impl FusionConfig {
    /// FusionConfig() for functions other than FUSIONS and ALT_SPLICE_JUNCTIONS
    pub fn empty() -> Self {
        FusionConfig {
            cohort_file: None,
            cache_fragments: true,
            skip_non_genic: false,
            write_chimeric_only: false,
            run_perf_checks: false,
            min_hard_filter_frags: 0,
            known_fusions: KnownFusionCache::new(),
        }
    }
}

#[allow(dead_code)]
pub struct IsofoxConfig {
    pub sample_id: String,
    pub dna_tumor_id: Option<String>,
    pub output_dir: String,
    pub output_id: Option<String>,
    pub functions: Vec<IsofoxFunction>,
    pub canonical_only: bool,
    pub bam_file: Option<String>,
    pub ref_genome_file: Option<String>,
    pub ref_genome_version: RefGenomeVersion,
    pub gene_distribution_file: Option<String>,
    pub alt_sj_cohort_file: Option<String>,
    pub filters: GeneRegionFilters,
    pub cancer_type: Option<String>,
    pub driver_genes: HashMap<String, DriverGene>,
    pub read_length: i32,
    pub max_fragment_length: i32,
    pub drop_duplicates: bool,
    pub exp_counts_file: Option<String>,
    pub exp_gc_ratios_file: Option<String>,
    pub panel_tpm_norm_file: Option<String>,
    pub neo_dir: Option<String>,
    pub apply_fragment_length_adjust: bool,
    pub fragment_size_data: Vec<FragmentSize>,
    pub write_types: Vec<WriteType>,
    pub fragment_length_sampling_count: i32,
    pub fusions: FusionConfig,
    pub perf_debug: bool,
    pub gene_read_limit: i32,
    pub run_validations: bool,
    pub run_perf_checks: bool,
    pub threads: usize,
    pub log_read_ids: Vec<String>,
    pub single_map_quality: i32,
    pub ensembl_data_dir: Option<String>,
    pub ensembl_given: bool,
}

pub fn register(cb: &mut ConfigBuilder) {
    cb.add_str("sample", true, "Sample ID");
    cb.add_path("ensembl_data_dir", false, "Ensembl data file directory");
    cb.add_str("output_dir", false, "Output directory");
    cb.add_str("output_id", false, "Output file suffix");
    cb.add_logging_options();
    cb.add_path("driver_gene_panel", false, "Path to driver gene panel");
    cb.add_str("functions", false, "List of functional routines to run (see documentation)");
    cb.add_flag("canonical_only", "Check all transcripts, not just canonical");
    cb.add_integer("gene_read_limit", "Per-gene limit on max reads processed (0 = not applied)", 0);
    cb.add_str("ref_genome_version", true, "Ref genome version, 37 or 38");
    cb.add_path("ref_genome", true, "Path to the reference genome fasta file");
    cb.add_integer("long_frag_limit", "Max RNA fragment size", DEFAULT_MAX_FRAGMENT_SIZE);
    cb.add_flag("drop_dups", "Include duplicate fragments in expression calculations");
    cb.add_path("gene_distribution_file", false, "Gene distribution for medians and percentile data");
    cb.add_path("alt_sj_cohort_file", false, "Alternate splice junction cohort file");
    cb.add_str("cancer_type", false, "Cancer type for gene cohort distribution");
    cb.add_integer("frag_length_min_count", "Fragment length measurement - min read fragments required", DEFAULT_FRAG_LENGTH_MIN_COUNT);
    cb.add_flag("frag_length_by_gene", "Write fragment lengths by gene");
    cb.add_path("bam_file", true, "RNA BAM file location");
    let write_types = if crate::version::v201() { WriteType::ALL_201.join(", ") } else { WriteType::ALL.join(", ") };
    cb.add_str_default("write_type", false, &format!("Write types: {write_types}"), "");
    cb.add_path("exp_counts_file", false, "File with generated expected expression rates per transcript");
    cb.add_path("exp_gc_ratios_file", false, "File with generated expected GC ratios per transcript");
    cb.add_path("neo_dir", false, "Path to Neo output");
    cb.add_str("tumor", false, "DNA tumor sample ID, used to read neoepitope file for annotation");
    cb.add_path("panel_tpm_norm_file", false, "Panel TPM normalisation file");
    cb.add_integer("read_length", "Sample sequencing read length, if 0 then is inferred from reads", 0);
    let v201 = crate::version::v201();
    cb.add_integer("single_map_qual", "Map quality for reads mapped to a single location", if v201 { DEFAULT_SINGLE_MAP_QUALITY_201 } else { DEFAULT_SINGLE_MAP_QUALITY });
    cb.add_str_default(
        "exp_rate_frag_lengths",
        false,
        "Fragment sizes and weights for expected transcript calcs (format: length1-freq1;length3-freq2 eg 100-10;150-20) in integer terms",
        DEFAULT_EXPECTED_RATE_LENGTHS,
    );
    cb.add_flag("validate", "Run auto-validations");
    cb.add_flag("run_perf_checks", "Run performance logging routines");
    cb.add_path("filter_reads_file", false, "Only process reads in this file");
    // GeneRegionFilters
    cb.add_path("gene_id_file", false, "Restricted set of Gene IDs in CSV file");
    cb.add_str("specific_chr", false, "Optional: list of chromosomes separated by ;");
    cb.add_str("specific_regions", false, "Optional: list of regions in format chromosome:positionStart-positionEnd, separated by ;");
    if v201 {
        cb.add_str("enriched_gene_ids", false, "List of geneIds to treat as enriched");
        cb.add_path("excluded_regions", false, "List of excluded regions");
    }
    // FusionConfig
    if v201 {
        cb.add_flag("write_chimeric_reads", "Write chimeric read data");
        cb.add_flag("write_chimeric_frags", "Write chimeric fragment data");
    }
    cb.add_flag("fusion_skip_non_genic", "Skip non-genic fusion fragments");
    cb.add_flag("fusion_write_chimeric_only", "Write chimeric reads but no other fusion processing");
    cb.add_path("known_fusion_file", false, "Known fusion file");
    cb.add_path("fusion_cohort_file", false, "Cohort file previously generated");
    cb.add_integer("fusion_min_frags_filter", "Hard filter chimeric translocations", DEFAULT_HARD_FILTER_MIN_FRAGS);
    cb.add_flag("run_fusion_perfs", "Write chimeric fragment data");
    cb.add_thread_options(1);
    cb.add_flag("perf_debug", "Detailed performance tracking and logging");
    cb.add_str("log_read_ids", false, "Log specific read IDs, separated by ';'");
}

/// a CSV file with a GeneId column (loadGeneIdsFile)
fn load_gene_ids_file(filename: &str) -> Result<Vec<String>> {
    let lines = hmf_common::java_text::read_lines(filename).with_context(|| format!("failed to load gene ID file({filename})"))?;
    let mut ids = Vec::new();
    if lines.is_empty() {
        return Ok(ids);
    }
    let header = &lines[0];
    let fields = hmf_common::java_text::fields_index_map(header, ',');
    let idx = *fields.get("GeneId").ok_or_else(|| anyhow!("gene ID file({filename}) has no GeneId column"))?;
    for line in &lines[1..] {
        let values: Vec<&str> = line.split(',').collect();
        if let Some(v) = values.get(idx) {
            ids.push(v.to_string());
        }
    }
    Ok(ids)
}

pub fn parse_fragment_sizes(value: Option<&str>) -> Result<Vec<FragmentSize>> {
    let mut sizes = Vec::new();
    if let Some(v) = value {
        for item in hmf_common::java_text::split(v, ';') {
            let parts: Vec<&str> = hmf_common::java_text::split(item, '-');
            let length = hmf_common::java_text::parse_int(parts.first().copied().unwrap_or("")).ok_or_else(|| anyhow!("invalid fragment length: {item}"))?;
            let freq = hmf_common::java_text::parse_int(parts.get(1).copied().unwrap_or("")).ok_or_else(|| anyhow!("invalid fragment length: {item}"))?;
            sizes.push(FragmentSize { length, frequency: freq.max(1) });
        }
    }
    Ok(sizes)
}

impl IsofoxConfig {
    pub fn from(cb: &ConfigBuilder) -> Result<IsofoxConfig> {
        let sample_id = cb.value("sample").unwrap_or_default();
        let mut functions = Vec::new();
        if cb.has_value("functions") {
            for f in cb.value("functions").unwrap().split(';') {
                functions.push(IsofoxFunction::value_of(f)?);
            }
        } else {
            functions = vec![IsofoxFunction::TranscriptCounts, IsofoxFunction::AltSpliceJunctions, IsofoxFunction::Fusions];
        }
        info!("running function(s): {}", functions.iter().map(|f| f.name()).collect::<Vec<_>>().join(","));

        let output_dir = cb.value("output_dir").map(|d| hmf_common::config_builder::check_add_dir_separator(&d));
        let ref_genome_version = RefGenomeVersion::from_java(&cb.value("ref_genome_version").unwrap_or_default())?;

        let mut driver_genes = HashMap::new();
        if cb.has_value("driver_gene_panel") {
            let file = cb.value("driver_gene_panel").unwrap();
            let genes = read_driver_genes(&file).with_context(|| format!("failed to read driver gene file: {file}"))?;
            for g in genes {
                driver_genes.insert(g.gene.clone(), g);
            }
        }

        let mut filters = GeneRegionFilters { immune_gene_regions: immune_regions(ref_genome_version), ..Default::default() };
        if cb.has_value("gene_id_file") {
            let file = cb.value("gene_id_file").unwrap();
            filters.restricted_gene_ids = load_gene_ids_file(&file)?;
            if !filters.restricted_gene_ids.is_empty() {
                info!("file({}) loaded {} restricted genes", file, filters.restricted_gene_ids.len());
            }
        }
        if crate::version::v201() {
            filters.enriched_gene_ids = match cb.value("enriched_gene_ids") {
                Some(v) if cb.has_value("enriched_gene_ids") => v.split(';').map(|s| s.to_string()).collect(),
                _ => default_enriched_gene_ids(ref_genome_version),
            };
            if let Some(file) = cb.value("excluded_regions") {
                let regions = hmf_common::regions::load_chr_base_regions(&file, false)?;
                for (chr, list) in regions.iter_unordered() {
                    filters.excluded_regions.insert(chr.clone(), list.clone());
                }
                info!("file({}) loaded {} excluded regions", file, filters.excluded_regions.len());
            }
            if cb.has_flag("write_chimeric_reads") || cb.has_flag("write_chimeric_frags") {
                bail!("-write_chimeric_reads and -write_chimeric_frags (ISOFOX 2.0.1 debug output) are not ported");
            }
        }
        filters.specific_regions = SpecificRegions::from_java(cb.value("specific_regions").as_deref(), cb.value("specific_chr").as_deref())?;
        filters.has_specific_regions = filters.specific_regions.has_filters();
        if !filters.specific_regions.regions.is_empty() {
            for r in &filters.specific_regions.regions {
                info!("filtering for specific region: {}:{}-{}", r.chromosome, r.start, r.end);
            }
        } else if !filters.specific_regions.chromosomes.is_empty() {
            info!("filtering for specific chromosomes: {}", filters.specific_regions.chromosomes.join(";"));
        }

        let write_types = match cb.value("write_type") {
            Some(v) if !v.is_empty() => v.split(';').map(WriteType::value_of).collect::<Result<Vec<_>>>()?,
            _ => Vec::new(),
        };
        let threads = cb.integer("threads").max(1) as usize;
        let (exp_counts_file, exp_gc_ratios_file) = if functions.contains(&IsofoxFunction::TranscriptCounts) {
            (cb.value("exp_counts_file"), cb.value("exp_gc_ratios_file"))
        } else {
            (None, None)
        };
        let apply_fragment_length_adjust = exp_counts_file.is_some();
        let default_sampling = if apply_fragment_length_adjust { DEFAULT_FRAG_LENGTH_MIN_COUNT } else { 0 };
        let fragment_length_sampling_count = if cb.has_value("frag_length_min_count") { cb.integer("frag_length_min_count") } else { default_sampling };

        let fusions = if functions.contains(&IsofoxFunction::Fusions) || functions.contains(&IsofoxFunction::AltSpliceJunctions) {
            let mut known = KnownFusionCache::new();
            if let Some(f) = cb.value("known_fusion_file") {
                if !known.load_from_file(&f) {
                    warn!("failed to load known fusion file({f})");
                }
            }
            FusionConfig {
                cohort_file: cb.value("fusion_cohort_file"),
                cache_fragments: false,
                skip_non_genic: cb.has_flag("fusion_skip_non_genic"),
                write_chimeric_only: cb.has_flag("fusion_write_chimeric_only"),
                run_perf_checks: cb.has_flag("run_fusion_perfs"),
                min_hard_filter_frags: cb.integer("fusion_min_frags_filter"),
                known_fusions: known,
            }
        } else {
            FusionConfig::empty()
        };

        let log_read_ids = match cb.value("log_read_ids") {
            Some(v) if cb.has_value("log_read_ids") => v.split(';').map(|s| s.to_string()).collect(),
            _ => Vec::new(),
        };

        Ok(IsofoxConfig {
            sample_id,
            dna_tumor_id: cb.value("tumor"),
            output_dir: output_dir.unwrap_or_default(),
            output_id: cb.value("output_id"),
            functions,
            canonical_only: cb.has_flag("canonical_only"),
            bam_file: cb.value("bam_file"),
            ref_genome_file: cb.value("ref_genome"),
            ref_genome_version,
            gene_distribution_file: cb.value("gene_distribution_file"),
            alt_sj_cohort_file: cb.value("alt_sj_cohort_file"),
            filters,
            cancer_type: cb.value("cancer_type"),
            driver_genes,
            read_length: cb.integer("read_length"),
            max_fragment_length: cb.integer("long_frag_limit"),
            drop_duplicates: cb.has_flag("drop_dups"),
            exp_counts_file,
            exp_gc_ratios_file,
            panel_tpm_norm_file: cb.value("panel_tpm_norm_file"),
            neo_dir: cb.value("neo_dir"),
            apply_fragment_length_adjust,
            fragment_size_data: parse_fragment_sizes(cb.value("exp_rate_frag_lengths").as_deref())?,
            write_types,
            fragment_length_sampling_count,
            fusions,
            perf_debug: cb.has_flag("perf_debug"),
            gene_read_limit: cb.integer("gene_read_limit"),
            run_validations: cb.has_flag("validate"),
            run_perf_checks: cb.has_flag("run_perf_checks"),
            threads,
            log_read_ids,
            single_map_quality: cb.integer("single_map_qual") as i16 as i32,
            ensembl_data_dir: cb.value("ensembl_data_dir"),
            ensembl_given: cb.has_value("ensembl_data_dir"),
        })
    }

    pub fn is_valid(&self) -> bool {
        if self.output_dir.is_empty() {
            hmf_common::error!("no output directory specified");
            return false;
        }
        if self.fragment_size_data.iter().any(|x| x.length <= 0 || x.frequency < 0) {
            hmf_common::error!("invalid fragment lengths");
            return false;
        }
        let bam_ok = self.bam_file.as_deref().map(|b| !b.is_empty() && std::path::Path::new(b).exists()).unwrap_or(false);
        if !bam_ok && !self.run_fusions_only() {
            hmf_common::error!("BAM file({}) missing or not found", self.bam_file.as_deref().unwrap_or("null"));
            return false;
        }
        if self.sample_id.is_empty() {
            hmf_common::error!("sampleId missing");
            return false;
        }
        if self.ref_genome_file.is_none() {
            hmf_common::error!("ref genome missing");
            return false;
        }
        if self.filters.restricted_gene_ids.is_empty()
            && !self.filters.has_specific_regions
            && (self.write_type(WriteType::Exon) || self.write_type(WriteType::Read))
        {
            warn!("writing exon and/or read data for all transcripts may be slow and generate large output files");
        }
        true
    }

    pub fn run_function(&self, f: IsofoxFunction) -> bool {
        self.functions.contains(&f)
    }
    pub fn run_fusions_only(&self) -> bool {
        self.functions.contains(&IsofoxFunction::Fusions) && self.functions.len() == 1
    }
    pub fn run_statistics_only(&self) -> bool {
        self.functions.contains(&IsofoxFunction::Statistics) && self.functions.len() == 1
    }
    pub fn require_fragment_length_calcs(&self) -> bool {
        self.write_types.contains(&WriteType::FragLength) || self.fragment_length_sampling_count > 0
    }
    pub fn write_type(&self, w: WriteType) -> bool {
        self.write_types.contains(&w)
    }
    pub fn require_gc_ratio_calcs(&self) -> bool {
        self.write_types.contains(&WriteType::GcRatio) || self.exp_gc_ratios_file.is_some()
    }
    pub fn apply_gc_bias_adjust(&self) -> bool {
        self.exp_gc_ratios_file.is_some()
    }
    pub fn form_output_file(&self, file_id: &str) -> String {
        match &self.output_id {
            Some(id) => format!("{}{}.isf.{}.{}", self.output_dir, self.sample_id, id, file_id),
            None => format!("{}{}.isf.{}", self.output_dir, self.sample_id, file_id),
        }
    }
}

#[cfg(test)]
impl IsofoxConfig {
    /// IsofoxConfig(refGenome): the configuration the Java tests construct
    pub fn for_tests() -> IsofoxConfig {
        IsofoxConfig {
            sample_id: "TEST".into(),
            dna_tumor_id: None,
            output_dir: String::new(),
            output_id: None,
            functions: vec![IsofoxFunction::TranscriptCounts, IsofoxFunction::AltSpliceJunctions, IsofoxFunction::RetainedIntrons],
            canonical_only: false,
            bam_file: None,
            ref_genome_file: None,
            ref_genome_version: RefGenomeVersion::V37,
            gene_distribution_file: None,
            alt_sj_cohort_file: None,
            filters: GeneRegionFilters { immune_gene_regions: immune_regions(RefGenomeVersion::V37), ..Default::default() },
            cancer_type: None,
            driver_genes: HashMap::new(),
            read_length: 0,
            max_fragment_length: DEFAULT_MAX_FRAGMENT_SIZE,
            drop_duplicates: false,
            exp_counts_file: None,
            exp_gc_ratios_file: None,
            panel_tpm_norm_file: None,
            neo_dir: None,
            apply_fragment_length_adjust: false,
            fragment_size_data: Vec::new(),
            write_types: Vec::new(),
            fragment_length_sampling_count: 0,
            fusions: FusionConfig::empty(),
            perf_debug: false,
            gene_read_limit: 0,
            run_validations: true,
            run_perf_checks: false,
            threads: 1,
            log_read_ids: Vec::new(),
            single_map_quality: DEFAULT_SINGLE_MAP_QUALITY,
            ensembl_data_dir: None,
            ensembl_given: false,
        }
    }
}
