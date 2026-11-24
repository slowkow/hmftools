//! PrepConfig: the command line of the Java PrepApplication.

use super::filters::*;
use super::types::*;
use crate::common::args::Args;
use anyhow::{bail, Result};

pub const DEFAULT_CHR_PARTITION_SIZE: i32 = 1_000_000;

const FLAGS: &[&str] = &["no_cache_bam", "no_trim_read_id", "no_clean_up", "perf_debug", "log_debug"];
const ITEMS: &[&str] = &[
    "sample",
    "bam_file",
    "ref_genome",
    "ref_genome_version",
    "known_hotspot_file",
    "saga_fasta",
    "partition_size",
    "sequencing_type",
    "write_types",
    "specific_chr",
    "specific_regions",
    "log_read_ids",
    "max_frag_length_override",
    "bam_validation",
    "min_align_bases",
    "min_junction_frags",
    // not Java options: see common::base_qual
    "low_base_qual_threshold",
    "high_base_qual_threshold",
    "bamtool",
    "threads",
    "output_dir",
    "output_id",
    "log_level",
];

/// The options of the step and whether each is a flag (for the tests of the combined mode).
#[cfg(test)]
pub(crate) fn accepted_options() -> Vec<(&'static str, bool)> {
    FLAGS.iter().map(|x| (*x, true)).chain(ITEMS.iter().map(|x| (*x, false))).collect()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum WriteType {
    PrepJunction,
    PrepRead,
    PrepBam,
    UnsortedBam,
    CacheBam,
    FragmentLengthDist,
    DiscordantStats,
}

impl WriteType {
    pub fn file_id(self) -> &'static str {
        match self {
            WriteType::PrepJunction => "junction.tsv",
            WriteType::PrepRead => "read.tsv",
            WriteType::PrepBam => "bam",
            WriteType::UnsortedBam => "unsorted.bam",
            WriteType::CacheBam => "cache",
            WriteType::FragmentLengthDist => "fragment_length.tsv",
            WriteType::DiscordantStats => "disc_stats.tsv",
        }
    }

    fn parse(name: &str) -> Result<Option<WriteType>> {
        Ok(Some(match name {
            "PREP_JUNCTION" => WriteType::PrepJunction,
            "PREP_READ" => WriteType::PrepRead,
            "PREP_BAM" => WriteType::PrepBam,
            "UNSORTED_BAM" => WriteType::UnsortedBam,
            "CACHE_BAM" => WriteType::CacheBam,
            "FRAGMENT_LENGTH_DIST" => WriteType::FragmentLengthDist,
            "DISCORDANT_STATS" => WriteType::DiscordantStats,
            // the assembly step's types are accepted and left out
            "JUNC_ASSEMBLY" | "ASSEMBLY_READ" | "BREAKEND" | "VCF" | "PHASED_ASSEMBLY" | "ALIGNMENT" | "DECOY_MATCHES" | "PHASE_GROUP_BUILDING" => return Ok(None),
            _ => bail!("invalid write type: {name}"),
        }))
    }

    /// WriteType.parsePrepTypes
    pub fn parse_prep_types(config: Option<&str>) -> Result<Vec<WriteType>> {
        let standard = [WriteType::PrepJunction, WriteType::PrepBam, WriteType::FragmentLengthDist, WriteType::DiscordantStats];
        let Some(config) = config.filter(|c| !c.is_empty()) else { return Ok(standard.to_vec()) };
        let mut types: Vec<WriteType> = Vec::new();
        if config.contains("PREP_STANDARD") {
            types.extend(standard);
        }
        for item in config.split(';') {
            if matches!(item, "PREP_STANDARD" | "ASSEMBLY_STANDARD" | "ASSEMBLY_ALL") {
                continue;
            }
            if let Some(t) = WriteType::parse(item)? {
                if !types.contains(&t) {
                    types.push(t);
                }
            }
        }
        Ok(types)
    }
}

/// SpecificRegions: `-specific_regions chr:start-end;...` (or chr:start:end) or `-specific_chr a;b`
#[derive(Clone, Debug, Default)]
pub struct SpecificRegions {
    /// (chromosome name, start, end), sorted by chromosome rank and start
    pub regions: Vec<(String, i32, i32)>,
    pub chromosomes: Vec<String>,
}

impl SpecificRegions {
    pub fn from(regions: Option<&str>, chromosomes: Option<&str>) -> Result<Self> {
        let mut s = SpecificRegions::default();
        if let Some(r) = regions.filter(|r| !r.is_empty()) {
            for item in r.split(';') {
                let parts: Vec<&str> = item.split(':').collect();
                let (chr, start, end) = if parts.len() == 3 {
                    (parts[0], parts[1].parse::<i32>(), parts[2].parse::<i32>())
                } else if parts.len() == 2 && parts[1].contains('-') {
                    let (a, b) = parts[1].split_once('-').unwrap();
                    (parts[0], a.parse::<i32>(), b.parse::<i32>())
                } else {
                    bail!("invalid specific region: {item}");
                };
                let (Ok(start), Ok(end)) = (start, end) else { bail!("invalid specific region: {item}") };
                if start < 1 || end < start {
                    bail!("invalid specific region: {item}");
                }
                s.regions.push((chr.to_string(), start, end));
            }
            // ChrBaseRegion.compareTo: chromosome rank, then start; the sort is stable
            s.regions.sort_by(|a, b| {
                if a.0 == b.0 {
                    a.1.cmp(&b.1)
                } else {
                    let (ra, rb) = (chromosome_rank(&a.0), chromosome_rank(&b.0));
                    if ra == rb { a.0.cmp(&b.0) } else { ra.cmp(&rb) }
                }
            });
            for r in &s.regions {
                if !s.chromosomes.contains(&r.0) {
                    s.chromosomes.push(r.0.clone());
                }
            }
        } else if let Some(c) = chromosomes.filter(|c| !c.is_empty()) {
            s.chromosomes = c.split(';').map(|x| x.to_string()).collect();
            // String.split drops trailing empty items
            while s.chromosomes.last().is_some_and(|x| x.is_empty()) {
                s.chromosomes.pop();
            }
        }
        Ok(s)
    }
    pub fn has_filters(&self) -> bool {
        !self.regions.is_empty() || !self.chromosomes.is_empty()
    }
    pub fn exclude_chromosome(&self, chr: &str) -> bool {
        !(self.chromosomes.is_empty() || self.chromosomes.iter().any(|c| c == chr))
    }
}

pub struct PrepConfig {
    pub sample_ids: Vec<String>,
    pub bam_files: Vec<String>,
    pub ref_genome_file: String,
    pub v38: bool,
    pub known_hotspot_file: Option<String>,
    pub saga_fasta_file: Option<String>,
    pub partition_size: i32,
    pub filter: ReadFilterConfig,
    pub output_dir: String,
    pub output_id: Option<String>,
    pub write_types: Vec<WriteType>,
    pub bam_tool_path: Option<String>,
    pub threads: usize,
    pub use_cache_bam: bool,
    pub trim_read_id: bool,
    pub specific: SpecificRegions,
    pub log_read_ids: Vec<String>,
    pub perf_debug: bool,
    pub max_fragment_length_override: i32,
    pub no_clean_up: bool,
}

/// RefGenomeVersion.versionedChromosome
pub fn versioned_chromosome(v38: bool, chr: &str) -> String {
    let bare = chr.strip_prefix("chr").unwrap_or(chr);
    if v38 { format!("chr{bare}") } else { bare.to_string() }
}

pub const HUMAN_CHROMOSOMES: [&str; 24] =
    ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"];

const LENGTHS_37: [i32; 24] = [
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878,
    107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566,
];
const LENGTHS_38: [i32; 24] = [
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328,
    107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415,
];

/// RefGenomeCoordinates: the length of a human chromosome by its index in `HUMAN_CHROMOSOMES`
pub fn chromosome_length(v38: bool, index: usize) -> i32 {
    if v38 { LENGTHS_38[index] } else { LENGTHS_37[index] }
}

impl PrepConfig {
    pub fn parse(argv: &[String]) -> Result<PrepConfig> {
        let args = Args::parse(argv, FLAGS, ITEMS)?;
        crate::common::base_qual::set_from_args(&args)?;
        if args.has_flag("log_debug") || args.value("log_level") == Some("DEBUG") || args.value("log_level") == Some("TRACE") {
            crate::common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        let list = |name: &str| -> Vec<String> { args.value(name).map(|v| v.split(',').map(|s| s.to_string()).collect()).unwrap_or_default() };
        let sample_ids = list("sample");
        let bam_files = list("bam_file");
        if sample_ids.is_empty() || sample_ids.len() != bam_files.len() {
            bail!("samples({}) and bam({}) not matched or empty", sample_ids.len(), bam_files.len());
        }
        let Some(ref_genome_file) = args.value("ref_genome").map(|s| s.to_string()) else { bail!("missing config: refGenome(false) outputDir(true)") };
        let output_dir = match args.value("output_dir") {
            Some(d) if d.is_empty() || d.ends_with('/') => d.to_string(),
            Some(d) => format!("{d}/"),
            None => match bam_files[0].rfind('/') {
                Some(i) => format!("{}/", &bam_files[0][..i]),
                None => "./".to_string(),
            },
        };
        let output_id = args.value("output_id").map(|s| s.to_string());
        let v38 = match args.value("ref_genome_version") {
            None | Some("V37") | Some("37") | Some("HG37") => false,
            Some("V38") | Some("38") | Some("HG38") => true,
            Some(v) => bail!("cannot parse ref genome version {v}"),
        };
        for name in ["bam_file", "ref_genome", "known_hotspot_file", "saga_fasta", "bamtool"] {
            if let Some(v) = args.value(name) {
                for path in v.split(',') {
                    if !std::path::Path::new(path).exists() {
                        bail!("config item {name}: invalid path = {path}");
                    }
                }
            }
        }
        let mut filter = ReadFilterConfig::new(args.integer("min_align_bases", MIN_ALIGNMENT_BASES)?, args.integer("min_junction_frags", MIN_JUNCTION_SUPPORT)?);
        if let Some(t) = args.value("sequencing_type") {
            filter.sequencing = SequencingType::parse(t).ok_or_else(|| anyhow::anyhow!("invalid sequencing type: {t}"))?;
        }
        let write_types = WriteType::parse_prep_types(args.value("write_types"))?;
        let specific = SpecificRegions::from(args.value("specific_regions"), args.value("specific_chr"))
            .map_err(|_| anyhow::anyhow!("invalid specific config: regions({}) chromosomes({})", args.value_or("specific_regions", "null"), args.value_or("specific_chr", "null")))?;
        let log_read_ids = args.value("log_read_ids").map(|v| v.split(';').map(|s| s.to_string()).collect()).unwrap_or_default();
        let has_filters = specific.has_filters();
        Ok(PrepConfig {
            trim_read_id: !args.has_flag("no_trim_read_id") && !has_filters && bam_files.len() == 1,
            use_cache_bam: !args.has_flag("no_cache_bam") && !has_filters,
            sample_ids,
            bam_files,
            ref_genome_file,
            v38,
            known_hotspot_file: args.value("known_hotspot_file").map(|s| s.to_string()),
            saga_fasta_file: args.value("saga_fasta").map(|s| s.to_string()),
            partition_size: args.integer("partition_size", DEFAULT_CHR_PARTITION_SIZE)?,
            filter,
            output_dir,
            output_id,
            write_types,
            bam_tool_path: args.value("bamtool").map(|s| s.to_string()),
            threads: args.integer("threads", 1)?.max(1) as usize,
            specific,
            log_read_ids,
            perf_debug: args.has_flag("perf_debug"),
            max_fragment_length_override: args.integer("max_frag_length_override", -1)?,
            no_clean_up: args.has_flag("no_clean_up"),
        })
    }

    pub fn sample_id(&self) -> &str {
        &self.sample_ids[0]
    }

    /// formFilename: `<dir><sample>.esvee.prep[.<output id>].<type>`
    pub fn form_filename(&self, write_type: WriteType, sample_id: &str) -> String {
        let mut name = format!("{}{}.esvee.prep", self.output_dir, sample_id);
        if let Some(id) = &self.output_id {
            name.push('.');
            name.push_str(id);
        }
        name.push('.');
        name.push_str(write_type.file_id());
        name
    }

    pub fn has_write_type(&self, t: WriteType) -> bool {
        self.write_types.contains(&t)
    }

    /// writeReads: whether read groups are completed across partitions and written
    pub fn write_reads(&self) -> bool {
        self.has_write_type(WriteType::PrepBam) || self.has_write_type(WriteType::PrepRead)
    }

    /// A config for tests, as the Java PrepConfig(int partitionSize)
    #[cfg(test)]
    pub fn for_tests(partition_size: i32) -> PrepConfig {
        PrepConfig {
            sample_ids: Vec::new(),
            bam_files: Vec::new(),
            ref_genome_file: String::new(),
            v38: false,
            known_hotspot_file: None,
            saga_fasta_file: None,
            partition_size,
            filter: ReadFilterConfig::default(),
            output_dir: String::new(),
            output_id: None,
            write_types: Vec::new(),
            bam_tool_path: None,
            threads: 1,
            use_cache_bam: false,
            trim_read_id: false,
            specific: SpecificRegions::default(),
            log_read_ids: Vec::new(),
            perf_debug: false,
            max_fragment_length_override: -1,
            no_clean_up: false,
        }
    }
}
