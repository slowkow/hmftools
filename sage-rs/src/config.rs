//! Command line configuration. Option names and defaults follow SageConfig / SageCallConfig.

use crate::constants::*;
use anyhow::{bail, Context, Result};
use std::collections::HashMap;

const FLAGS: &[&str] = &[
    "include_mt", "panel_only", "no_fragment_sync", "skip_bqr", "skip_msi_jitter", "high_depth_mode",
    "log_candidates", "perf_summary", "log_debug", "germline", "disable_soft_filter", "disable_hard_filter", "write_frag_lengths",
    "bqr_write_plot", "log_evidence", "run_tinc", "write_tinc_vcf", "write_fit_variants", "require_gene", "vis_pass_only",
];

#[derive(Clone, Debug)]
pub struct Config {
    pub tumor_ids: Vec<String>,
    pub tumor_bams: Vec<String>,
    pub reference_ids: Vec<String>,
    pub reference_bams: Vec<String>,
    pub ref_genome: String,
    pub v38: bool,
    pub output_vcf: String,
    pub hotspots: Option<String>,
    pub panel_bed: Option<String>,
    pub driver_gene_panel: Option<String>,
    pub high_confidence_bed: Option<String>,
    pub ensembl_data_dir: Option<String>,
    pub jitter_bqr_dir: String,
    pub threads: usize,
    /// not a Java option: the most tumor record data (MB) kept per thread and sample between the
    /// candidate and evidence stages; 0 reads the file twice, as Java does
    pub read_cache_mb: usize,
    pub include_mt: bool,
    pub panel_only: bool,
    pub sync_fragments: bool,
    pub log_candidates: bool,
    /// not a Java option: report where the wall-clock time went
    pub perf_summary: bool,
    pub specific_chr: Vec<String>,
    /// (chromosome, start, end)
    pub specific_regions: Vec<(String, i32, i32)>,
    pub read_length: i32,
    pub min_map_quality: i32,
    pub max_read_depth: i32,
    pub max_read_depth_panel: i32,
    pub region_slice_size: i32,
    pub max_partition_slices: i32,
    pub ref_sample_count: usize,
    pub fixed_map_qual_penalty: i32,
    pub improper_pair_penalty: i32,
    pub read_events_penalty: f64,
    pub hard_min_tumor_raw_alt_support: i32,
    pub hard_min_tumor_qual: i32,
    pub hard_min_tumor_vaf: f64,
    pub skip_bqr: bool,
    pub skip_msi_jitter: bool,
    pub log_evidence: bool,
    pub min_avg_base_qual: i32,
    pub min_avg_base_qual_hotspot: i32,
    pub disable_soft_filter: bool,
    pub disable_hard_filter: bool,
    pub filtered_max_germline_alt_support: i32,
    /// -germline: germline filters are not applied
    pub germline: bool,
    pub high_depth_mode: bool,
    /// SageCallConfig.RunTinc: -run_tinc with exactly one reference sample
    pub run_tinc: bool,
    pub write_tinc_vcf: bool,
    pub write_fit_variants: bool,
    pub pon_file: Option<String>,
    pub gnomad_freq_dir: Option<String>,
    pub gnomad_freq_file: Option<String>,
    /// SageAppendConfig.InputVcf and FilterToGenes (-require_gene): append mode only
    pub input_vcf: Option<String>,
    pub require_gene: bool,
    /// -write_frag_lengths: write <output>.frag_lengths.tsv.gz
    pub write_frag_lengths: bool,
    /// SageConfig.SpecificVariants (-specific_var_file): only these candidates are kept
    pub specific_variants: Vec<(String, crate::types::SimpleVariant)>,
    /// SageConfig.Visualiser (-vis_variants, -vis_pass_only, -vis_max_support_reads, -vis_output_dir,
    /// -vis_purple_vcf); never enabled in append mode
    pub vis: crate::vis::VisConfig,
}

fn list(v: Option<&String>) -> Vec<String> {
    v.map(|s| s.split(',').map(|x| x.to_string()).collect()).unwrap_or_default()
}

/// SpecificRegions.addRegion: a region which overlaps an existing one widens the first such region
/// (regions that come to overlap each other as a result are not merged)
pub fn add_specific_region(chromosomes: &mut Vec<String>, regions: &mut Vec<(String, i32, i32)>, chr: &str, start: i32, end: i32) {
    if !chromosomes.iter().any(|c| c == chr) {
        chromosomes.push(chr.to_string());
    }
    for (c, s, e) in regions.iter_mut() {
        if c == chr && start <= *e && end >= *s {
            *s = (*s).min(start);
            *e = (*e).max(end);
            return;
        }
    }
    regions.push((chr.to_string(), start, end));
}

/// SimpleVariant.loadSimpleVariants: a delimited file with Chromosome, Position, Ref and Alt columns
pub fn load_simple_variants(path: &str) -> Result<Vec<(String, crate::types::SimpleVariant)>> {
    let mut lines = crate::refdata::read_lines(path)?.into_iter();
    let header = lines.next().context("empty file")?;
    let delim = if header.contains(',') { ',' } else { '\t' };
    let col = |name: &str| header.split(delim).position(|f| f == name).with_context(|| format!("missing column {name}"));
    let (chr_i, pos_i, ref_i, alt_i) = (col("Chromosome")?, col("Position")?, col("Ref")?, col("Alt")?);
    let mut out = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split(delim).collect();
        let get = |i: usize| values.get(i).copied().with_context(|| format!("short line: {line}"));
        out.push((
            get(chr_i)?.to_string(),
            crate::types::SimpleVariant::new(get(pos_i)?.parse()?, get(ref_i)?.as_bytes(), get(alt_i)?.as_bytes()),
        ));
    }
    Ok(out)
}

impl Config {
    pub fn from_args(args: &[String]) -> Result<Config> {
        Self::parse(args, false)
    }

    /// SageAppendConfig: the common options plus -input_vcf and -require_gene; there are no tumor samples
    pub fn from_append_args(args: &[String]) -> Result<Config> {
        Self::parse(args, true)
    }

    fn parse(args: &[String], append: bool) -> Result<Config> {
        let mut vals: HashMap<String, String> = HashMap::new();
        let mut i = 0;
        while i < args.len() {
            let a = &args[i];
            let Some(key) = a.strip_prefix('-') else { bail!("unexpected argument: {a}") };
            if FLAGS.contains(&key) {
                vals.insert(key.to_string(), "true".to_string());
                i += 1;
            } else {
                let v = args.get(i + 1).with_context(|| format!("missing value for -{key}"))?;
                vals.insert(key.to_string(), v.clone());
                i += 2;
            }
        }
        let get = |k: &str| vals.get(k);
        let req = |k: &str| vals.get(k).cloned().with_context(|| format!("missing required option -{k}"));
        let int = |k: &str, d: i32| -> Result<i32> { Ok(get(k).map(|v| v.parse()).transpose()?.unwrap_or(d)) };
        let flag = |k: &str| vals.contains_key(k);

        if let Some(seq) = get("sequencing_type") {
            if seq.eq_ignore_ascii_case("SBX") {
                set_sequencing_type(SequencingType::Sbx);
            } else if seq.eq_ignore_ascii_case("ULTIMA") {
                set_sequencing_type(SequencingType::Ultima);
            } else if !seq.eq_ignore_ascii_case("ILLUMINA") {
                bail!("sequencing_type {seq} is not supported; ILLUMINA, SBX and ULTIMA are implemented");
            }
        }
        if is_ultima() {
            // UltimaQualRecalibration.CFG_FILENAME and FilterConfig.ULTIMA_CANDIDATE_MIN_HIGH_BQ
            crate::ultima::load_bqr_cache(get("ultima_bqr_file").map(|s| s.as_str()))?;
            if let Some(v) = get("ultima_cand_min_high_bq") {
                crate::ultima::set_candidate_min_high_bq_threshold(v.parse()?);
            }
        }
        let tumor_bams = list(get("tumor_bam"));
        let reference_bams = list(get("reference_bam"));
        let dir_of = |p: &String| {
            std::path::Path::new(p).parent().map(|d| d.to_string_lossy().to_string()).unwrap_or_else(|| ".".into())
        };
        let jitter_bqr_dir = match get("jitter_bqr_dir") {
            Some(d) => d.clone(),
            None => reference_bams.first().or(tumor_bams.first()).map(dir_of).unwrap_or_else(|| ".".into()),
        };
        let version = req("ref_genome_version")?;
        let mut specific_regions = Vec::new();
        if let Some(s) = get("specific_regions") {
            for item in s.split(';') {
                let (chr, range) = item.split_once(':').context("invalid specific_regions")?;
                let (a, b) = range.split_once('-').context("invalid specific_regions")?;
                specific_regions.push((chr.to_string(), a.parse()?, b.parse()?));
            }
        }
        let mut specific_chr: Vec<String> = get("specific_chr").map(|s| s.split(';').map(|x| x.to_string()).collect()).unwrap_or_default();
        for (c, _, _) in &specific_regions {
            if !specific_chr.contains(c) {
                specific_chr.push(c.clone());
            }
        }
        // SageConfig: the variants file is only read when no region or chromosome was given, and each
        // variant adds a region around itself
        let vis = if append {
            crate::vis::VisConfig::default()
        } else {
            crate::vis::VisConfig::from_options(
                get("vis_variants").map(|s| s.as_str()),
                flag("vis_pass_only"),
                int("vis_max_support_reads", 0)?,
                get("vis_output_dir").map(|s| s.as_str()),
                get("vis_purple_vcf").map(|s| s.as_str()),
                &req("output_vcf")?,
            )?
        };
        let mut specific_variants = Vec::new();
        if !vis.specific_variants.is_empty() {
            // each visualised variant adds a region around itself, and the variants file is then not read
            for (chr, v) in &vis.specific_variants {
                add_specific_region(&mut specific_chr, &mut specific_regions, chr, v.pos - VIS_VARIANT_BUFFER, v.pos + VIS_VARIANT_BUFFER);
            }
        } else if specific_regions.is_empty() && specific_chr.is_empty() {
            if let Some(path) = get("specific_var_file") {
                specific_variants =
                    load_simple_variants(path).with_context(|| format!("failed to load specific variants file: {path}"))?;
                for (chr, v) in &specific_variants {
                    add_specific_region(
                        &mut specific_chr,
                        &mut specific_regions,
                        chr,
                        v.pos - VIS_VARIANT_BUFFER,
                        v.pos + VIS_VARIANT_BUFFER,
                    );
                }
            }
        }
        Ok(Config {
            tumor_ids: if append { Vec::new() } else { list(Some(&req("tumor")?)) },
            tumor_bams,
            reference_ids: list(get("reference")),
            reference_bams,
            ref_genome: req("ref_genome")?,
            v38: version.contains("38"),
            output_vcf: req("output_vcf")?,
            hotspots: get("hotspots").cloned(),
            panel_bed: get("panel_bed").cloned(),
            driver_gene_panel: get("driver_gene_panel").cloned(),
            high_confidence_bed: get("high_confidence_bed").cloned(),
            ensembl_data_dir: get("ensembl_data_dir").cloned(),
            jitter_bqr_dir,
            threads: int("threads", 1)? as usize,
            read_cache_mb: int("read_cache_mb", crate::bamio::DEFAULT_READ_CACHE_MB)?.max(0) as usize,
            include_mt: flag("include_mt"),
            panel_only: flag("panel_only"),
            sync_fragments: !flag("no_fragment_sync"),
            log_candidates: flag("log_candidates"),
            perf_summary: flag("perf_summary"),
            specific_chr,
            specific_regions,
            // SBX and Ultima reads are long and variable, so the length is not sampled from the BAM
            read_length: int("read_length", if is_sbx() || is_ultima() { NON_ILLUMINA_MAX_READ_LENGTH } else { 0 })?,
            min_map_quality: int("min_map_quality", DEFAULT_MIN_MAP_QUALITY)?,
            // when append runs in high depth mode the panel depth limit applies to every tier
            max_read_depth: if append && flag("high_depth_mode") && get("max_read_depth").is_none() {
                int("max_read_depth_panel", DEFAULT_MAX_READ_DEPTH_PANEL)?
            } else {
                int("max_read_depth", DEFAULT_MAX_READ_DEPTH)?
            },
            max_read_depth_panel: int("max_read_depth_panel", DEFAULT_MAX_READ_DEPTH_PANEL)?,
            region_slice_size: int("slice_size", DEFAULT_SLICE_SIZE)?,
            max_partition_slices: int("max_partition_slices", DEFAULT_MAX_PARTITION_SLICES)?,
            ref_sample_count: int("ref_sample_count", 1)? as usize,
            fixed_map_qual_penalty: int("fixed_qual_penalty", DEFAULT_MAP_QUAL_FIXED_PENALTY)?,
            improper_pair_penalty: int("improper_pair_qual_penalty", DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY)?,
            read_events_penalty: get("read_events_qual_penalty")
                .map(|v| v.parse())
                .transpose()?
                .unwrap_or(DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY),
            hard_min_tumor_raw_alt_support: int("hard_min_tumor_raw_alt_support", DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT)?,
            hard_min_tumor_qual: int("hard_min_tumor_qual", DEFAULT_HARD_MIN_TUMOR_QUAL)?,
            hard_min_tumor_vaf: get("hard_min_tumor_vaf").map(|v| v.parse()).transpose()?.unwrap_or(DEFAULT_HARD_MIN_TUMOR_VAF),
            skip_bqr: flag("skip_bqr"),
            skip_msi_jitter: flag("skip_msi_jitter"),
            log_evidence: flag("log_evidence"),
            min_avg_base_qual: int("min_avg_base_qual", 25)?,
            min_avg_base_qual_hotspot: int("min_avg_base_qual_hotspot", 18)?,
            disable_soft_filter: flag("disable_soft_filter"),
            disable_hard_filter: flag("disable_hard_filter"),
            // FilterConfig: the limit is relaxed when TINC will run, and TINC re-applies the default afterwards
            filtered_max_germline_alt_support: int("filtered_max_germline_alt_support", if flag("run_tinc") { 10 } else { 3 })?,
            germline: flag("germline"),
            high_depth_mode: flag("high_depth_mode"),
            run_tinc: flag("run_tinc") && list(get("reference")).len() == 1,
            write_tinc_vcf: flag("write_tinc_vcf"),
            write_fit_variants: flag("write_fit_variants"),
            pon_file: get("pon_file").cloned(),
            gnomad_freq_dir: get("gnomad_freq_dir").cloned(),
            gnomad_freq_file: get("gnomad_freq_file").cloned(),
            input_vcf: if append { Some(req("input_vcf")?) } else { None },
            require_gene: flag("require_gene"),
            write_frag_lengths: flag("write_frag_lengths"),
            specific_variants,
            vis,
        })
    }

    /// SageConfig.processChromosome
    pub fn process_chromosome(&self, chr: &str) -> bool {
        if !self.specific_chr.is_empty() && !self.specific_chr.iter().any(|c| c == chr) {
            return false;
        }
        if crate::types::is_human_chromosome(chr) {
            return true;
        }
        self.include_mt && crate::types::is_mito(chr)
    }
}
