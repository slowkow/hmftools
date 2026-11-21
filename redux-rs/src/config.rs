//! Command line configuration. Option names and defaults follow ReduxConfig.

use crate::args::Args;
use crate::regions::SpecificRegions;
use anyhow::{bail, Result};

pub const FLAGS: &[&str] = &[
    "form_consensus", "drop_duplicates", "bqr_jitter_msi_only", "fail_supp_no_mate_cigar", "std_chr_only", "skip_fully_unmapped",
    "skip_duplicate_marking", "skip_unmapping", "unmap_mt", "unmap_alt_decoy", "keep_interim_bams", "no_write_bam", "run_checks",
    "log_read_cache", "parallel_concat", "process_read_only", "log_debug", "bqr_disable", "bqr_write_plot", "bqr_use_all_regions",
    "write_msi_plots", "write_msi_site_file", "msi_no_bq_filter", "umi_enabled", "umi_duplex", "umi_variable_length",
    "umi_base_diff_stats", "ultima_drop_qual_tags",
];

pub const ITEMS: &[&str] = &[
    "sample", "input_bam", "output_bam", "ref_genome", "ref_genome_version", "sequencing_type", "read_length", "read_output",
    "unmap_regions", "bamtool", "bam_validation", "threads", "partition_ratio", "output_dir", "output_id", "log_level",
    "specific_regions", "specific_chr", "log_read_ids", "perf_log_time", "specific_region_filter", "log_dup_group_size",
    "write_read_base_length", "bqr_debug_tnc", "bqr_debug_qual", "bqr_debug_position", "bqr_debug_consensus", "ref_genome_msi_file",
    "msi_min_map_quality", "msi_max_sites_per_type", "max_site_alt_contribution", "msi_model_coefficients", "msi_model_error_rates",
    "umi_type", "umi_duplex_delim", "umi_defined_ids", "sbx_max_dup_dist", "compress_level",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequencingType {
    Illumina,
    Sbx,
    Ultima,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FilterReadsType {
    None,
    Read,
    Mate,
    MateAndSupp,
}

pub struct Config {
    pub args: Args,
    pub specific_region_filter: FilterReadsType,
    pub sample_id: String,
    pub bam_files: Vec<String>,
    pub ref_genome_file: String,
    pub v38: bool,
    pub sequencing_type: SequencingType,
    pub output_bam: Option<String>,
    /// ends with '/'
    pub output_dir: String,
    pub output_id: Option<String>,
    pub form_consensus: bool,
    pub bqr_jitter_msi_only: bool,
    pub skip_unmapping: bool,
    pub skip_fully_unmapped_reads: bool,
    pub fail_on_missing_supp_mate_cigar: bool,
    pub skip_duplicate_marking: bool,
    pub drop_duplicates: bool,
    pub unmap_regions_file: Option<String>,
    pub unmap_mt: bool,
    pub unmap_alt_decoys: bool,
    pub std_chr_only: bool,
    pub read_length: i32,
    pub bamtool: Option<String>,
    pub parallel_concat: bool,
    pub threads: usize,
    pub partition_thread_ratio: usize,
    pub specific: SpecificRegions,
    pub write_bam: bool,
    pub keep_interim_bams: bool,
    pub log_debug: bool,
    pub umi_enabled: bool,
    pub umi_duplex: bool,
    /// -ultima_drop_qual_tags, read only for Ultima (ReduxConfig: UltimaRoutines.setUltimaConfig)
    pub ultima_drop_qual_tags: bool,
    /// DuplicatesConfig.SbxMaxDuplicateDistance (-sbx_max_dup_dist)
    pub sbx_max_dup_dist: i32,
    pub umi: crate::umi::UmiConfig,
    /// not a Java option: the deflate level of the output BAM (htslib levels 0-9). The default, 4, gives a file
    /// within about 0.5% of Java's size (zlib level 5); 5 is about 2.5% smaller and 7% slower, 3 about 5% larger
    /// and 9% faster, 1 about 7% larger and 21% faster.
    pub compress_level: u32,
}

impl Config {
    pub fn from_args(argv: &[String]) -> Result<Config> {
        let args = Args::parse(argv, FLAGS, ITEMS)?;
        let Some(sample_id) = args.value("sample").map(str::to_string) else { bail!("missing -sample") };
        let Some(bams) = args.value("input_bam") else { bail!("missing -input_bam") };
        let bam_files: Vec<String> = bams.split(',').map(str::to_string).collect();
        let Some(ref_genome_file) = args.value("ref_genome").map(str::to_string) else { bail!("missing -ref_genome") };
        let v38 = match args.value_or("ref_genome_version", "37") {
            "38" | "V38" | "HG38" => true,
            "37" | "V37" | "HG19" => false,
            other => bail!("invalid ref_genome_version: {other}"),
        };
        let sequencing_type = match args.value_or("sequencing_type", "ILLUMINA") {
            "ILLUMINA" => SequencingType::Illumina,
            "SBX" => SequencingType::Sbx,
            "ULTIMA" => SequencingType::Ultima,
            other => bail!("invalid sequencing_type: {other}"),
        };
        let output_bam = args.value("output_bam").map(str::to_string);
        if let Some(o) = &output_bam {
            if bam_files.iter().any(|b| b == o) {
                bail!("output BAM({o}) matches input BAM filename");
            }
        }
        let dir_of = |f: &str| match f.rfind('/') {
            Some(i) => f[..=i].to_string(),
            None => String::new(),
        };
        let output_dir = if let Some(d) = args.value("output_dir") {
            let d = if d.ends_with('/') { d.to_string() } else { format!("{d}/") };
            std::fs::create_dir_all(&d)?;
            d
        } else if let Some(o) = &output_bam {
            dir_of(o)
        } else {
            dir_of(&bam_files[0])
        };
        let msi_only = args.has_flag("bqr_jitter_msi_only");
        let skip_unmapping = msi_only || args.has_flag("skip_unmapping");
        let threads = args.integer("threads", 1)?.max(1) as usize;
        let umi = crate::umi::UmiConfig::from_args(&args)?;
        let (umi_enabled, umi_duplex) = (umi.enabled, umi.duplex);
        let specific = SpecificRegions::from(args.value("specific_regions"), args.value("specific_chr"))?;
        let specific_region_filter = if specific.has_filters() {
            match args.value_or("specific_region_filter", "READ") {
                "NONE" => FilterReadsType::None,
                "READ" => FilterReadsType::Read,
                "MATE" => FilterReadsType::Mate,
                "MATE_AND_SUPP" => FilterReadsType::MateAndSupp,
                other => bail!("invalid specific_region_filter: {other}"),
            }
        } else {
            FilterReadsType::None
        };
        Ok(Config {
            specific_region_filter,
            sample_id,
            bam_files,
            ref_genome_file,
            v38,
            sequencing_type,
            output_bam,
            output_dir,
            output_id: args.value("output_id").map(str::to_string),
            form_consensus: umi_enabled || args.has_flag("form_consensus"),
            bqr_jitter_msi_only: msi_only,
            skip_unmapping,
            skip_fully_unmapped_reads: msi_only || skip_unmapping || args.has_flag("skip_fully_unmapped"),
            fail_on_missing_supp_mate_cigar: !msi_only && args.has_flag("fail_supp_no_mate_cigar"),
            skip_duplicate_marking: msi_only || args.has_flag("skip_duplicate_marking"),
            drop_duplicates: args.has_flag("drop_duplicates"),
            unmap_regions_file: if skip_unmapping { None } else { args.value("unmap_regions").map(str::to_string) },
            unmap_mt: args.has_flag("unmap_mt"),
            unmap_alt_decoys: args.has_flag("unmap_alt_decoy"),
            std_chr_only: args.has_flag("std_chr_only"),
            read_length: args.integer("read_length", 0)?,
            bamtool: args.value("bamtool").map(str::to_string),
            parallel_concat: args.has_flag("parallel_concat"),
            threads,
            partition_thread_ratio: if threads <= 1 { 1 } else { args.integer("partition_ratio", 2)? as usize },
            specific,
            write_bam: !args.has_flag("no_write_bam") && !msi_only,
            keep_interim_bams: args.has_flag("keep_interim_bams"),
            log_debug: args.has_flag("log_debug") || args.value("log_level") == Some("DEBUG"),
            umi_enabled,
            umi_duplex,
            ultima_drop_qual_tags: sequencing_type == SequencingType::Ultima && args.has_flag("ultima_drop_qual_tags"),
            sbx_max_dup_dist: args.integer("sbx_max_dup_dist", crate::sbx::DEFAULT_SBX_MAX_DUPLICATE_DISTANCE)?,
            umi,
            compress_level: args.integer("compress_level", 4)?.clamp(0, 9) as u32,
            args,
        })
    }

    /// ReduxConfig.formFilename: <output_dir><sample>.redux.<type>[.<output_id>].tsv
    pub fn form_filename(&self, file_type: &str) -> String {
        let mut f = format!("{}{}.redux.{}", self.output_dir, self.sample_id, file_type);
        if let Some(id) = &self.output_id {
            f.push('.');
            f.push_str(id);
        }
        f.push_str(".tsv");
        f
    }

    pub fn is_illumina(&self) -> bool {
        self.sequencing_type == SequencingType::Illumina
    }
    pub fn is_sbx(&self) -> bool {
        self.sequencing_type == SequencingType::Sbx
    }
    pub fn is_ultima(&self) -> bool {
        self.sequencing_type == SequencingType::Ultima
    }
}
