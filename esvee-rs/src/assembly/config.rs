//! AssemblyConfig: the command line options of the assembly step.

use std::collections::HashSet;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::assembly::constants::*;
use crate::assembly::globals::*;
use crate::assembly::ptr::P;
use crate::assembly::ref_genome::{RefGenomeRef, RefGenomeSource};
use crate::assembly::seq_tech_utils::set_seq_tech_specifics;
use crate::assembly::types::junction::{JuncP, Junction};
use crate::common::args::Args;
use crate::common::java_sort::java_sort;
use crate::common::regions::SpecificRegions;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum WriteType {
    JuncAssembly,
    AssemblyRead,
    Breakend,
    Vcf,
    PhasedAssembly,
    Alignment,
    DecoyMatches,
    PhaseGroupBuilding,
}

impl WriteType {
    pub const ALL: [WriteType; 8] = [
        WriteType::JuncAssembly,
        WriteType::AssemblyRead,
        WriteType::Breakend,
        WriteType::Vcf,
        WriteType::PhasedAssembly,
        WriteType::Alignment,
        WriteType::DecoyMatches,
        WriteType::PhaseGroupBuilding,
    ];
    pub fn file_id(self) -> &'static str {
        match self {
            WriteType::JuncAssembly => "assembly.tsv",
            WriteType::AssemblyRead => "assembly_read.tsv",
            WriteType::Breakend => "breakend.tsv",
            WriteType::Vcf => "raw.vcf.gz",
            WriteType::PhasedAssembly => "phased_assembly.tsv",
            WriteType::Alignment => "alignment.tsv",
            WriteType::DecoyMatches => "decoy_match.tsv",
            WriteType::PhaseGroupBuilding => "phase_group_building.tsv",
        }
    }
    fn from_name(name: &str) -> Option<Option<WriteType>> {
        // Some(None): a prep write type, accepted and ignored by the assembly step
        Some(Some(match name {
            "JUNC_ASSEMBLY" => WriteType::JuncAssembly,
            "ASSEMBLY_READ" => WriteType::AssemblyRead,
            "BREAKEND" => WriteType::Breakend,
            "VCF" => WriteType::Vcf,
            "PHASED_ASSEMBLY" => WriteType::PhasedAssembly,
            "ALIGNMENT" => WriteType::Alignment,
            "DECOY_MATCHES" => WriteType::DecoyMatches,
            "PHASE_GROUP_BUILDING" => WriteType::PhaseGroupBuilding,
            "PREP_JUNCTION" | "PREP_READ" | "PREP_BAM" | "UNSORTED_BAM" | "CACHE_BAM" | "FRAGMENT_LENGTH_DIST" | "DISCORDANT_STATS" => return Some(None),
            _ => return None,
        }))
    }

    /// WriteType.parseAssemblyTypes
    pub fn parse_assembly_types(config_str: Option<&str>) -> Result<HashSet<WriteType>> {
        let standard = [WriteType::Vcf, WriteType::Breakend, WriteType::JuncAssembly];
        let Some(config_str) = config_str.filter(|s| !s.is_empty()) else { return Ok(standard.into_iter().collect()) };
        let mut types: HashSet<WriteType> = HashSet::new();
        if config_str.contains("ASSEMBLY_STANDARD") {
            types.extend(standard);
        }
        if config_str.contains("ASSEMBLY_ALL") {
            types.extend(WriteType::ALL);
        }
        for item in config_str.split(';') {
            if item == "PREP_STANDARD" || item == "ASSEMBLY_STANDARD" || item == "ASSEMBLY_ALL" {
                continue;
            }
            match WriteType::from_name(item) {
                Some(Some(t)) => {
                    types.insert(t);
                }
                Some(None) => {}
                None => bail!("No enum constant com.hartwig.hmftools.esvee.common.WriteType.{item}"),
            }
        }
        Ok(types)
    }
}

pub struct AssemblyConfig {
    pub tumor_ids: Vec<String>,
    pub reference_ids: Vec<String>,
    pub tumor_bams: Vec<String>,
    pub reference_bams: Vec<String>,
    pub junction_file: Option<String>,
    pub prep_dir: String,
    pub ref_gen_version_v38: bool,
    pub ref_genome_file: Option<String>,
    pub ref_genome: Option<RefGenomeRef>,
    pub ref_genome_image_file: Option<String>,
    pub decoy_genome: Option<String>,
    pub saga_fasta_file: Option<String>,
    pub write_types: HashSet<WriteType>,
    pub output_dir: Option<String>,
    pub output_id: Option<String>,
    pub bam_tool_path: Option<String>,
    pub threads: i32,
    pub purple_dir: Option<String>,
    pub run_visualiser: bool,
    pub assembly_map_qual_threshold: i32,
    pub discordant_only_disabled: bool,
    pub discordant_rate_increment: f64,
    pub specific_chr_regions: SpecificRegions,
    pub specific_junctions: Vec<JuncP>,
    pub log_read_ids: Vec<String>,
    pub assembly_ref_base_write_max: i32,
    pub assembly_detailed_tsv: bool,
    pub apply_remote_phasing_read_check_threshold: bool,
    chromosome_lengths: Vec<(String, i32)>,
}

const FLAGS: &[&str] = &[
    "perf_debug", "disc_only_disabled", "remote_phase_read_check_threshold", "write_candidate_reads", "asm_build_debug", "asm_log_all_reads", "assembly_detailed_tsv",
    "run_assembly_vis", "log_debug",
];
const ITEMS: &[&str] = &[
    "tumor", "tumor_bam", "reference", "reference_bam", "junction_file", "esvee_prep_dir", "sequencing_type", "write_types", "ref_genome", "ref_genome_version",
    "decoy_genome", "saga_fasta", "bwa_lib", "log_read_ids", "specific_junctions", "perf_log_time", "asm_ref_base_write_max", "phase_process_limit",
    "asm_map_qual_threshold", "disc_rate_increment", "purple_dir", "bamtool", "specific_chr", "specific_regions", "output_dir", "output_id", "log_level", "threads",
];

fn split_list(args: &Args, name: &str) -> Vec<String> {
    match args.value(name) {
        Some(v) => v.split(',').map(|s| s.to_string()).collect(),
        None => Vec::new(),
    }
}

/// FileWriterUtils.checkAddDirSeparator
fn check_add_dir_separator(dir: &str) -> String {
    if dir.is_empty() || dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

/// FileWriterUtils.pathFromFile
fn path_from_file(file_path: &str) -> String {
    if file_path.contains('/') {
        let parent = Path::new(file_path).parent().map(|p| p.to_string_lossy().into_owned()).unwrap_or_default();
        format!("{parent}/")
    } else {
        "./".to_string()
    }
}

/// FileCommon.formOutputFile
pub fn form_output_file(output_dir: &str, sample_id: &str, app_stage: &str, file_type: &str, output_id: Option<&str>) -> String {
    let mut filename = format!("{output_dir}{sample_id}.{app_stage}");
    if let Some(id) = output_id {
        filename.push('.');
        filename.push_str(id);
    }
    filename.push('.');
    filename.push_str(file_type);
    filename
}

/// FileCommon.formPrepInputFilename
pub fn form_prep_input_filename(output_dir: &str, sample_id: &str, file_id: &str, output_id: Option<&str>) -> String {
    let with_id = form_output_file(output_dir, sample_id, "esvee.prep", file_id, output_id);
    if Path::new(&with_id).exists() {
        return with_id;
    }
    form_output_file(output_dir, sample_id, "esvee.prep", file_id, None)
}

fn form_prep_bam_filenames(output_dir: &str, sample_ids: &[String], output_id: Option<&str>) -> Vec<String> {
    let mut bam_files = Vec::with_capacity(sample_ids.len());
    for sample_id in sample_ids {
        let mut bam_file = format!("{output_dir}{sample_id}.esvee.prep");
        if let Some(id) = output_id {
            let with_id = format!("{bam_file}.{id}.bam");
            if Path::new(&with_id).exists() {
                bam_files.push(with_id);
                continue;
            }
        }
        bam_file.push_str(".bam");
        if !Path::new(&bam_file).exists() {
            return Vec::new();
        }
        bam_files.push(bam_file);
    }
    bam_files
}

fn load_chromosome_lengths(v38: bool) -> Vec<(String, i32)> {
    let text = if v38 { include_str!("resources/lengths.38.tsv") } else { include_str!("resources/lengths.37.tsv") };
    text.lines()
        .filter_map(|l| {
            let (c, len) = l.split_once('\t')?;
            Some((c.trim_start_matches("chr").to_string(), len.trim().parse().ok()?))
        })
        .collect()
}

impl AssemblyConfig {
    pub fn from_args(argv: &[String]) -> Result<AssemblyConfig> {
        let args = Args::parse(argv, FLAGS, ITEMS)?;
        if args.has_flag("log_debug") || matches!(args.value("log_level"), Some("DEBUG") | Some("TRACE")) {
            crate::common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
        }

        let mut prep_bam_dir: Option<String> = None;
        if args.has_value("tumor_bam") {
            prep_bam_dir = Some(path_from_file(&split_list(&args, "tumor_bam")[0]));
        } else if args.has_value("reference_bam") {
            prep_bam_dir = Some(path_from_file(&split_list(&args, "reference_bam")[0]));
        }
        let output_dir: Option<String> =
            if !args.has_value("output_dir") && prep_bam_dir.is_some() { prep_bam_dir.clone() } else { args.value("output_dir").map(check_add_dir_separator) };
        let prep_dir = if let Some(d) = args.value("esvee_prep_dir") {
            check_add_dir_separator(d)
        } else if let Some(d) = &prep_bam_dir {
            d.clone()
        } else {
            output_dir.clone().context("no output directory or prep directory given")?
        };
        let output_id = args.value("output_id").map(|s| s.to_string());

        let mut tumor_ids = split_list(&args, "tumor");
        let mut tumor_bams: Vec<String> = Vec::new();
        let mut reference_bams: Vec<String> = Vec::new();
        let mut reference_ids: Vec<String> = Vec::new();
        let prep_tumor_bams = form_prep_bam_filenames(&prep_dir, &tumor_ids, output_id.as_deref());
        if prep_tumor_bams.len() == tumor_ids.len() {
            tumor_bams.extend(prep_tumor_bams);
        } else {
            tumor_bams.extend(split_list(&args, "tumor_bam"));
        }
        if args.has_value("reference") {
            reference_ids.extend(split_list(&args, "reference"));
            let prep_ref_bams = form_prep_bam_filenames(&prep_dir, &reference_ids, output_id.as_deref());
            if prep_ref_bams.len() == reference_ids.len() {
                reference_bams.extend(prep_ref_bams);
            } else {
                reference_bams.extend(split_list(&args, "reference_bam"));
            }
        }
        if tumor_ids.is_empty() && reference_ids.is_empty() {
            bail!("no tumor or reference IDs provided");
        }
        if tumor_ids.len() != tumor_bams.len() || reference_ids.len() != reference_bams.len() {
            bail!("tumor and reference IDs must match BAM files");
        }
        // in germline mode, the reference values become the tumor ones
        if tumor_ids.is_empty() && !reference_ids.is_empty() {
            tumor_ids.append(&mut reference_ids);
            tumor_bams.append(&mut reference_bams);
        }

        let junction_file = match args.value("junction_file") {
            Some(f) => Some(f.to_string()),
            None => {
                let f = form_prep_input_filename(&prep_dir, &tumor_ids[0], "junction.tsv", output_id.as_deref());
                if Path::new(&f).exists() { Some(f) } else { None }
            }
        };

        let ref_gen_version_v38 = match args.value("ref_genome_version") {
            None | Some("37") | Some("V37") | Some("HG19") => false,
            Some("38") | Some("V38") | Some("HG38") => true,
            Some(v) => bail!("invalid ref genome version: {v}"),
        };
        let ref_genome_file = args.value("ref_genome").map(|s| s.to_string());
        let ref_genome: Option<RefGenomeRef> = match &ref_genome_file {
            Some(f) => Some(std::sync::Arc::new(RefGenomeSource::open(f).with_context(|| format!("failed to load ref genome: {f}"))?)),
            None => bail!("missing required config: ref_genome"),
        };
        let ref_genome_image_file = ref_genome_file.as_ref().map(|f| format!("{f}.img"));

        let write_types = WriteType::parse_assembly_types(args.value("write_types"))?;

        if let Some(t) = args.value("sequencing_type") {
            match SequencingType::parse(t) {
                Some(st) => set_sequencing_type(st),
                None => bail!("invalid sequencing type: {t}"),
            }
        }
        set_seq_tech_specifics();

        let specific_chr_regions = SpecificRegions::from(args.value("specific_regions"), args.value("specific_chr"))?;
        let mut specific_junctions: Vec<JuncP> = Vec::new();
        if let Some(s) = args.value("specific_junctions") {
            let delim = if s.contains('_') { '_' } else { ';' };
            for item in s.split(delim) {
                if let Some(j) = Junction::from_config_str(item) {
                    specific_junctions.push(P::new(j));
                }
            }
            crate::debug!("loaded {} specific junctions", specific_junctions.len());
            java_sort(&mut specific_junctions, |a, b| a.compare_to(b));
        }
        let purple_dir = args.value("purple_dir").map(|s| s.to_string());
        let run_visualiser = args.has_flag("run_assembly_vis") && !specific_junctions.is_empty();
        if run_visualiser || purple_dir.is_some() {
            set_clear_cached_read(false);
        }
        let has_filters = specific_chr_regions.has_filters() || !specific_junctions.is_empty();
        if write_types.contains(&WriteType::AssemblyRead) && !has_filters {
            crate::warn!("writing assembly reads to TSV without region filtering may result in large output files & impact performance");
        }
        let assembly_detailed_tsv = args.has_flag("assembly_detailed_tsv") || has_filters;
        let log_read_ids: Vec<String> = args.value("log_read_ids").map(|v| v.split(';').map(|s| s.to_string()).collect()).unwrap_or_default();

        let perf_log_time = args.decimal("perf_log_time", 0.0)?;
        set_perf_log_time(perf_log_time);
        set_perf_debug(args.has_flag("perf_debug") || perf_log_time > 0.0);
        set_assembly_build_debug(args.has_flag("asm_build_debug"));
        // as in Java, set from the build debug flag
        set_assembly_log_all_reads(args.has_flag("asm_build_debug"));
        set_write_candidate_reads(args.has_flag("write_candidate_reads"));
        set_phase_processing_limit(args.integer("phase_process_limit", 0)?);

        let assembly_ref_base_write_max = if has_filters { 0 } else { args.integer("asm_ref_base_write_max", DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX)? };
        READ_ID_TRIMMER.reset(!has_filters && tumor_bams.len() == 1);

        Ok(AssemblyConfig {
            tumor_ids,
            reference_ids,
            tumor_bams,
            reference_bams,
            junction_file,
            prep_dir,
            ref_gen_version_v38,
            ref_genome_file,
            ref_genome,
            ref_genome_image_file,
            decoy_genome: args.value("decoy_genome").map(|s| s.to_string()),
            saga_fasta_file: args.value("saga_fasta").map(|s| s.to_string()),
            write_types,
            output_dir,
            output_id,
            bam_tool_path: args.value("bamtool").map(|s| s.to_string()),
            threads: args.integer("threads", 1)?,
            purple_dir,
            run_visualiser,
            assembly_map_qual_threshold: args.integer("asm_map_qual_threshold", DEFAULT_ASSEMBLY_MAP_QUAL_THRESHOLD)?,
            discordant_only_disabled: args.has_flag("disc_only_disabled"),
            discordant_rate_increment: args.decimal("disc_rate_increment", DEFAULT_DISC_RATE_INCREMENT)?,
            specific_chr_regions,
            specific_junctions,
            log_read_ids,
            assembly_ref_base_write_max,
            assembly_detailed_tsv,
            apply_remote_phasing_read_check_threshold: args.has_flag("remote_phase_read_check_threshold"),
            chromosome_lengths: load_chromosome_lengths(ref_gen_version_v38),
        })
    }

    /// the test constructor AssemblyConfig()
    pub fn new_test() -> AssemblyConfig {
        set_perf_debug(false);
        set_perf_log_time(0.0);
        set_phase_processing_limit(0);
        set_assembly_build_debug(false);
        set_write_candidate_reads(false);
        READ_ID_TRIMMER.reset(false);
        AssemblyConfig {
            tumor_ids: Vec::new(),
            reference_ids: Vec::new(),
            tumor_bams: Vec::new(),
            reference_bams: Vec::new(),
            junction_file: None,
            prep_dir: String::new(),
            ref_gen_version_v38: true,
            ref_genome_file: None,
            ref_genome: None,
            ref_genome_image_file: None,
            decoy_genome: None,
            saga_fasta_file: None,
            write_types: HashSet::new(),
            output_dir: None,
            output_id: None,
            bam_tool_path: None,
            threads: 0,
            purple_dir: None,
            run_visualiser: false,
            assembly_map_qual_threshold: -1,
            discordant_only_disabled: false,
            discordant_rate_increment: DEFAULT_DISC_RATE_INCREMENT,
            specific_chr_regions: SpecificRegions::default(),
            specific_junctions: Vec::new(),
            log_read_ids: Vec::new(),
            assembly_ref_base_write_max: 0,
            assembly_detailed_tsv: false,
            apply_remote_phasing_read_check_threshold: false,
            chromosome_lengths: load_chromosome_lengths(true),
        }
    }

    pub fn combined_sample_ids(&self) -> Vec<String> {
        self.tumor_ids.iter().chain(self.reference_ids.iter()).cloned().collect()
    }
    pub fn combined_bam_files(&self) -> Vec<String> {
        self.tumor_bams.iter().chain(self.reference_bams.iter()).cloned().collect()
    }
    pub fn sample_id(&self) -> &str {
        &self.tumor_ids[0]
    }
    pub fn output_filename(&self, write_type: WriteType) -> String {
        form_output_file(self.output_dir.as_deref().unwrap_or(""), self.sample_id(), "esvee", write_type.file_id(), self.output_id.as_deref())
    }
    /// RefGenomeCoordinates.length: zero for anything but chromosomes 1-22, X and Y
    pub fn chromosome_length(&self, chromosome: &str) -> i32 {
        if !crate::assembly::common_utils::human_chromosome_contains(chromosome) {
            return 0;
        }
        let bare = chromosome.strip_prefix("chr").unwrap_or(chromosome);
        self.chromosome_lengths.iter().find(|(c, _)| c == bare).map(|(_, l)| *l).unwrap_or(0)
    }
    pub fn ref_genome(&self) -> Option<&dyn crate::assembly::ref_genome::RefGenomeInterface> {
        self.ref_genome.as_deref()
    }
}
