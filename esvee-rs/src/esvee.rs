//! The four steps in one call, as the Java class EsveeApplication runs them: `esvee-rs <options>`, with no
//! step name. This is how nf-core/oncoanalyser runs ESVEE.
//!
//! EsveeApplication registers the options of all four steps on one ConfigBuilder, reads the command line
//! once, and runs prep, assembly, depth annotation and calling in turn, each step reading its options from
//! that shared ConfigBuilder. An option registered by several steps keeps its first registration, so it has
//! one value and one default for every step that reads it. Here the command line is checked as that
//! ConfigBuilder checks it, and then each step is given the command line under which it resolves its inputs
//! as its Java config does in this mode:
//!
//! - prep and depth take `-sample` and `-bam_file` when given, else the tumor ids followed by the reference
//!   ids, and the tumor BAMs followed by the reference BAMs (FileCommon.parseSampleIds and parseBamFiles).
//!   Both read the original BAMs.
//! - assembly runs as a sub-routine (`new AssemblyApplication(configBuilder, true)`): `-tumor_bam` and
//!   `-reference_bam` are not read. Its output directory is `-output_dir`, its prep directory
//!   `-esvee_prep_dir` or else the output directory, and its BAMs the prep BAMs found there.
//! - depth reads `-input_vcf` when given, else `<output dir><first sample>.esvee.raw.vcf.gz`. The caller
//!   reads `-input_vcf` when given (the same file as depth), else
//!   `<output dir><tumor, or reference>.esvee.ref_depth.vcf.gz`, and its prep directory is
//!   `-esvee_prep_dir` or the output directory.
//! - `-write_types` is read by prep and by assembly, each taking its own types from it.
//! - any other option goes to every step that registers it.

use anyhow::{bail, Result};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Str,
    Path,
    Paths,
    Integer,
    Decimal,
    Flag,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Str => "STRING",
            Kind::Path => "PATH",
            Kind::Paths => "PATHS",
            Kind::Integer => "INTEGER",
            Kind::Decimal => "DECIMAL",
            Kind::Flag => "FLAG",
        }
    }
}

// the steps that are given an option
const PREP: u8 = 1;
const ASSEMBLY: u8 = 2;
const DEPTH: u8 = 4;
const CALLER: u8 = 8;
/// given to no step as it is: `-stages`, and the sample and BAM options from which prep and depth are
/// given their lists
const NONE: u8 = 0;

struct Opt {
    name: &'static str,
    kind: Kind,
    required: bool,
    default: Option<&'static str>,
    desc: &'static str,
    steps: u8,
}

const fn opt(name: &'static str, kind: Kind, default: Option<&'static str>, desc: &'static str, steps: u8) -> Opt {
    Opt { name, kind, required: false, default, desc, steps }
}

const fn required(name: &'static str, kind: Kind, desc: &'static str, steps: u8) -> Opt {
    Opt { name, kind, required: true, default: None, desc, steps }
}

use Kind::*;

/// The options of EsveeApplication in the order Java registers them (its own `-stages`, then the
/// registerConfig of PrepConfig, AssemblyConfig, DepthConfig and CallerConfig), each with its first
/// registration, and the steps that are given it.
#[rustfmt::skip]
static OPTIONS: &[Opt] = &[
    opt("stages", Str, None, "Stages: PREP;ASSEMBLY;REF_DEPTH;CALLING separated by ';'", NONE),
    // PrepConfig
    opt("sample", Str, None, "List of samples separated by ','", NONE),
    opt("bam_file", Paths, None, "BAM file paths separated by ','", NONE),
    required("ref_genome_version", Str, "Ref genome version, 37 or 38", PREP | ASSEMBLY | DEPTH | CALLER),
    required("ref_genome", Path, "Path to reference genome fasta files", PREP | ASSEMBLY | DEPTH),
    opt("known_hotspot_file", Path, None, "Known fusion hotspot BED file", PREP | CALLER),
    opt("saga_fasta", Path, None, "SAGA resource FASTA file", PREP | ASSEMBLY),
    opt("partition_size", Integer, Some("1000000"), "Partition size", PREP),
    opt("sequencing_type", Str, Some("ILLUMINA"), "Sequencing types: ILLUMINA, ULTIMA, SBX", PREP | ASSEMBLY | CALLER),
    opt("write_types", Str, None, "Write types: PREP_JUNCTION;PREP_READ;PREP_BAM;UNSORTED_BAM;CACHE_BAM;FRAGMENT_LENGTH_DIST;DISCORDANT_STATS;JUNC_ASSEMBLY;ASSEMBLY_READ;BREAKEND;VCF;PHASED_ASSEMBLY;ALIGNMENT;DECOY_MATCHES;PHASE_GROUP_BUILDING separated by ';'", PREP | ASSEMBLY),
    opt("specific_chr", Str, None, "Restrict to chromosome(s) separated by ';'", PREP | ASSEMBLY | DEPTH | CALLER),
    opt("specific_regions", Str, None, "Restrict to regions(s) separated by ';' in format Chr:PosStart:PosEnd or Chr:PosStart-PosEnd", PREP | ASSEMBLY | DEPTH | CALLER),
    opt("log_read_ids", Str, None, "Log specific read IDs, separated by ';'", PREP | ASSEMBLY),
    opt("no_cache_bam", Flag, None, "Write a BAM to cache candidate reads", PREP),
    opt("no_trim_read_id", Flag, None, "Disable use of a shortened readId internally", PREP),
    opt("no_clean_up", Flag, None, "Keep candidate cache BAM files", PREP),
    opt("perf_debug", Flag, None, "Detailed performance tracking and logging", PREP | ASSEMBLY),
    opt("max_frag_length_override", Integer, Some("-1"), "Set max fragment length instead of calculating", PREP),
    opt("bam_validation", Str, Some("SILENT"), "BAM validation: STRICT (default), LENIENT or SILENT", PREP | DEPTH),
    opt("min_align_bases", Integer, Some("50"), "Junction fragment min aligned bases", PREP),
    opt("min_junction_frags", Integer, Some("2"), "Required fragments to keep a junction", PREP),
    opt("bamtool", Path, None, "Path to sambamba or samtools for sort/merge", PREP | ASSEMBLY),
    opt("threads", Integer, Some("1"), "Number of threads, 0 or 1 not multi-threaded", PREP | ASSEMBLY | DEPTH),
    opt("output_dir", Str, None, "Output directory", PREP | ASSEMBLY | DEPTH | CALLER),
    opt("output_id", Str, None, "Output file suffix", PREP | ASSEMBLY | DEPTH | CALLER),
    opt("log_debug", Flag, None, "Log at DEBUG level", PREP | ASSEMBLY | DEPTH | CALLER),
    opt("log_level", Str, None, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE", PREP | ASSEMBLY | DEPTH | CALLER),
    // AssemblyConfig
    opt("tumor", Str, None, "Tumor ID(s), separated by ','", ASSEMBLY | CALLER),
    opt("tumor_bam", Str, None, "Tumor BAM file(s), separated by ','", NONE),
    opt("reference", Str, None, "Reference ID(s), separated by ','", ASSEMBLY | CALLER),
    opt("reference_bam", Str, None, "Reference BAM file(s), separated by ','", NONE),
    opt("junction_file", Paths, None, "Esvee Prep junction file, default is to match by sample name", ASSEMBLY),
    opt("esvee_prep_dir", Paths, None, "Esvee prep input directory", ASSEMBLY | CALLER),
    opt("decoy_genome", Path, None, "Decoy genome image file", ASSEMBLY),
    opt("bwa_lib", Path, None, "Path to BWA library", ASSEMBLY),
    opt("specific_junctions", Str, None, "Specific junctions: format: Chromosome:Position:Orientation:Type (I, D or none), separated by ';'", ASSEMBLY),
    opt("perf_log_time", Decimal, Some("0.0"), "Log performance data for routine exceeding specified time (0 = disabled)", ASSEMBLY | DEPTH),
    opt("asm_ref_base_write_max", Integer, Some("200"), "Cap assembly ref bases in TSV and VCF, use zero to write all", ASSEMBLY),
    opt("phase_process_limit", Integer, Some("0"), "Exclude phase groups above this size from extension and phase sets", ASSEMBLY),
    opt("disc_only_disabled", Flag, None, "Disable discordant-only junctions", ASSEMBLY),
    opt("asm_map_qual_threshold", Integer, Some("10"), "Realign and test assemblies with average map-qual below this threshold", ASSEMBLY),
    opt("remote_phase_read_check_threshold", Flag, None, "Apply remote phase building max read check threshold", ASSEMBLY),
    opt("write_candidate_reads", Flag, None, "Write assembly candidate reads regardless of whether used", ASSEMBLY),
    opt("asm_build_debug", Flag, None, "Log assembly building working", ASSEMBLY),
    opt("asm_log_all_reads", Flag, None, "Log assembly all reads even if not used in support", ASSEMBLY),
    opt("assembly_detailed_tsv", Flag, None, "Log extra assembly info to TSV", ASSEMBLY),
    opt("disc_rate_increment", Decimal, Some("0.01"), "Discordant rate increment", ASSEMBLY),
    opt("purple_dir", Path, None, "Purple data directory for event visualisation", ASSEMBLY),
    opt("run_assembly_vis", Flag, None, "Run visualiser on specific junctions", ASSEMBLY),
    // DepthConfig
    opt("input_vcf", Path, None, "Input VCF", DEPTH | CALLER),
    opt("vcf_tag_prefix", Str, None, "VCF tag prefix for testing & comparison", DEPTH),
    opt("unmap_regions", Path, None, "Genome regions to unmap reads", DEPTH),
    opt("vaf_cap", Decimal, Some("0.001"), "Ref support depth limit as function of variant fragments", DEPTH),
    opt("proximity_distance", Integer, Some("2000"), "Proximity distance to group variants", DEPTH),
    // CallerConfig
    opt("manual_ref_depth", Integer, Some("0"), "Manually set ref depth for testing", CALLER),
    opt("pon_sv_file", Path, None, "PON for SV positions", CALLER),
    opt("pon_sgl_file", Path, None, "PON for SGL positions", CALLER),
    opt("pon_margin", Integer, Some("3"), "PON permitted matching position margin", CALLER),
    opt("pon_sgl_margin", Integer, Some("10"), "PON permitted matching position margin for SGLs", CALLER),
    opt("artefact_pon_sv_file", Path, None, "Additional artefact SV PON file", CALLER),
    opt("artefact_pon_sgl_file", Path, None, "Additional artefact SGL PON file", CALLER),
    opt("write_breakend_tsv", Flag, None, "Rewrite a breakend TSV for with additional caller annotations", CALLER),
    opt("filter_sgls", Flag, None, "Filter SGLs from VCF, intended for tumor-only mode, default=true in target panel", CALLER),
    opt("min_qual", Integer, Some("30"), "Min qual, panel default: 60", CALLER),
    opt("hotspot_min_qual", Integer, Some("30"), "Hotspot min qual", CALLER),
    opt("min_support", Integer, Some("4"), "Min support", CALLER),
    opt("hotspot_min_support", Integer, Some("2"), "Hotspot min support", CALLER),
    opt("sgl_min_support", Integer, Some("6"), "SGL min support", CALLER),
    opt("repeat_mask_file", Path, None, "Repeat mask definitions file", CALLER),
    opt("target_regions_bed", Path, None, "Target regions BED file", CALLER),
];

/// EsveeApplication.Stage, in the order the steps run
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Stage {
    Prep,
    Assembly,
    RefDepth,
    Calling,
}

impl Stage {
    fn from_name(name: &str) -> Option<Stage> {
        match name {
            "PREP" => Some(Stage::Prep),
            "ASSEMBLY" => Some(Stage::Assembly),
            "REF_DEPTH" => Some(Stage::RefDepth),
            "CALLING" => Some(Stage::Calling),
            _ => None,
        }
    }
}

/// Why the command line was refused, with the lines Java logs at ERROR.
#[derive(Debug, PartialEq)]
enum ParseError {
    /// an argument that is not a registered option, or a value where none fits: Java also logs the options
    IncorrectArgument(String),
    /// missing required options, paths that do not exist, numbers that do not parse
    Invalid(Vec<String>),
}

/// The values given on the command line, by the index of the option in `OPTIONS`.
struct Options {
    values: Vec<Option<String>>,
}

/// String.split(delim) without a limit: one item if the delimiter does not occur, else trailing empty
/// items are dropped
fn java_split(s: &str, delim: char) -> Vec<String> {
    if !s.contains(delim) {
        return vec![s.to_string()];
    }
    let mut parts: Vec<String> = s.split(delim).map(str::to_string).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// Files.exists(Paths.get(path)): the empty path is the current directory
fn path_exists(path: &str) -> bool {
    path.is_empty() || std::path::Path::new(path).exists()
}

impl Options {
    fn index(name: &str) -> usize {
        OPTIONS.iter().position(|o| o.name == name).unwrap_or_else(|| panic!("unregistered option {name}"))
    }

    fn value(&self, name: &str) -> Option<&str> {
        self.values[Self::index(name)].as_deref()
    }

    fn has(&self, name: &str) -> bool {
        self.value(name).is_some()
    }

    /// ConfigBuilder.parseCommandLine and isValid
    fn parse(args: &[String]) -> std::result::Result<Options, ParseError> {
        let mut values: Vec<Option<String>> = vec![None; OPTIONS.len()];
        let mut matched: Option<usize> = None;
        for argument in args {
            if let Some(index) = matched.take() {
                let o = &OPTIONS[index];
                // a leading '-' is taken as a value only by a number, which may be negative
                if argument.starts_with('-') && o.kind != Kind::Decimal && o.kind != Kind::Integer {
                    return Err(ParseError::IncorrectArgument(format!("config item({}) has invalid argument: {}", o.name, argument)));
                }
                values[index] = Some(argument.clone());
            } else if let Some(name) = argument.strip_prefix('-') {
                let Some(index) = OPTIONS.iter().position(|o| o.name == name) else {
                    return Err(ParseError::IncorrectArgument(format!("unregistered config item: {argument}")));
                };
                if OPTIONS[index].kind == Kind::Flag {
                    values[index] = Some("true".to_string());
                } else {
                    matched = Some(index);
                }
            } else {
                return Err(ParseError::IncorrectArgument(format!("expecting config name but encountered invalid argument: {argument}")));
            }
        }
        // an option at the end without its value is left unset, as in Java

        let mut errors = Vec::new();
        for (o, value) in OPTIONS.iter().zip(&values) {
            if o.required && value.is_none() {
                errors.push(format!("missing required config: {}: {}", o.name, o.desc));
            } else if o.kind == Kind::Path || o.kind == Kind::Paths {
                let Some(value) = value else { continue };
                let paths: Vec<&str> = if o.kind == Kind::Paths { value.split(',').collect() } else { vec![value.as_str()] };
                if !paths.iter().all(|p| p.contains('*') || p.contains('$') || path_exists(p)) {
                    errors.push(format!("invalid path for config: {} = {}", o.name, value));
                }
            } else if o.kind == Kind::Integer || o.kind == Kind::Decimal {
                let text = value.as_deref().or(o.default);
                let ok = match (o.kind, text) {
                    (Kind::Integer, Some(t)) => crate::depth::parse_java_int(t).is_some(),
                    (_, Some(t)) => crate::depth::parse_java_double(t).is_some(),
                    _ => false,
                };
                if !ok {
                    errors.push(format!("invalid type for config: {}: {} = {}", o.name, o.kind.name(), text.unwrap_or("null")));
                }
            }
        }
        if errors.is_empty() { Ok(Options { values }) } else { Err(ParseError::Invalid(errors)) }
    }

    /// EsveeApplication: every stage unless `-stages` lists some
    fn stages(&self) -> Result<Vec<Stage>> {
        let Some(value) = self.value("stages") else {
            return Ok(vec![Stage::Prep, Stage::Assembly, Stage::RefDepth, Stage::Calling]);
        };
        // split(ITEM_DELIM, -1): empty items are kept, and are not a stage
        let mut stages = Vec::new();
        for name in value.split(';') {
            let Some(stage) = Stage::from_name(name) else {
                bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.esvee.EsveeApplication.Stage.{name}");
            };
            stages.push(stage);
        }
        Ok(stages)
    }

    /// FileCommon.parseSampleBamLists of one option
    fn list(&self, name: &str) -> Vec<String> {
        self.value(name).map(|v| java_split(v, ',')).unwrap_or_default()
    }

    /// FileCommon.parseSampleIds (`sample`, tumor, reference) and parseBamFiles (`bam_file`, tumor_bam,
    /// reference_bam), as prep and depth read them
    fn combined_list(&self, own: &str, tumor: &str, reference: &str) -> Vec<String> {
        if self.has(own) {
            return self.list(own);
        }
        let mut items = self.list(tumor);
        items.extend(self.list(reference));
        items
    }

    fn sample_ids(&self) -> Vec<String> {
        self.combined_list("sample", "tumor", "reference")
    }

    fn bam_files(&self) -> Vec<String> {
        self.combined_list("bam_file", "tumor_bam", "reference_bam")
    }

    /// the command line of one step: the options given to it that have a value, in the order of `OPTIONS`
    fn step_args(&self, step: u8) -> Vec<String> {
        let mut args = Vec::new();
        if step == PREP || step == DEPTH {
            for (name, list) in [("sample", self.sample_ids()), ("bam_file", self.bam_files())] {
                // an empty list is left to the step to report
                if !list.is_empty() {
                    args.push(format!("-{name}"));
                    args.push(list.join(","));
                }
            }
        }
        for (o, value) in OPTIONS.iter().zip(&self.values) {
            let Some(value) = value else { continue };
            if o.steps & step == 0 {
                continue;
            }
            args.push(format!("-{}", o.name));
            if o.kind != Kind::Flag {
                args.push(value.clone());
            }
        }
        args
    }
}

/// ConfigBuilder.printItems: the options, one line each, sorted
fn item_lines() -> Vec<String> {
    let mut lines: Vec<String> = OPTIONS
        .iter()
        .map(|o| {
            let mut s = format!(" -{}:", o.name);
            if o.kind != Kind::Str {
                s.push_str(&format!(" type({})", o.kind.name()));
            }
            if let Some(d) = o.default.filter(|_| o.kind != Kind::Flag) {
                s.push_str(&format!(" default({d})"));
            }
            if o.required {
                s.push_str(" REQUIRED");
            }
            s.push_str(&format!(" desc({})", o.desc));
            s
        })
        .collect();
    lines.sort();
    lines.insert(0, "registered config items:".to_string());
    lines
}

/// What Java prints for `-version`; `-help` adds the options.
fn version_line() -> String {
    format!("Esvee version {}", crate::VERSION)
}

/// Java prints the version only when `-version` is the sole argument, and refuses it among other options;
/// here it is honoured anywhere.
fn wants_version(argv: &[String]) -> bool {
    argv.iter().any(|a| a == "-version")
}

/// ConfigUtils.setLogLevel: Level.valueOf of `-log_level`; a name log4j does not know is an exception
fn check_log_level(options: &Options) -> Result<()> {
    if options.has("log_debug") {
        crate::common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    } else if let Some(level) = options.value("log_level") {
        let name = level.trim().to_uppercase();
        let Some(rank) = crate::depth::LEVEL_NAMES.iter().position(|x| *x == name) else {
            bail!("java.lang.IllegalArgumentException: Unknown level constant [{name}].");
        };
        if rank as u8 >= crate::depth::LEVEL_DEBUG {
            crate::common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
        }
    }
    Ok(())
}

/// EsveeApplication.main and run
pub fn run(argv: &[String]) -> Result<()> {
    if wants_version(argv) {
        println!("{}", version_line());
        std::process::exit(0);
    }
    if argv.len() == 1 && argv[0] == "-help" {
        println!("{}", version_line());
        for line in item_lines() {
            println!("{line}");
        }
        std::process::exit(0);
    }
    let options = match Options::parse(argv) {
        Ok(o) => o,
        Err(ParseError::IncorrectArgument(message)) => {
            crate::error!("{message}");
            for line in item_lines() {
                crate::info!("{line}");
            }
            std::process::exit(1);
        }
        Err(ParseError::Invalid(messages)) => {
            for message in messages {
                crate::error!("{message}");
            }
            std::process::exit(1);
        }
    };
    check_log_level(&options)?;
    if crate::version_log_due() {
        crate::info!("{}", version_line());
    }
    let stages = options.stages()?;
    let start = std::time::Instant::now();

    if stages.contains(&Stage::Prep) {
        crate::prep::run(&options.step_args(PREP))?;
    }
    if stages.contains(&Stage::Assembly) {
        crate::assembly::run(&options.step_args(ASSEMBLY))?;
    }
    if stages.contains(&Stage::RefDepth) {
        crate::depth::run(&options.step_args(DEPTH))?;
    }
    if stages.contains(&Stage::Calling) {
        crate::caller::run(&options.step_args(CALLER))?;
    }
    crate::info!("Esvee complete, mins({:.3})", start.elapsed().as_millis() as f64 / 60000.0);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| s.to_string()).collect()
    }

    /// a file and a folder that exist wherever the tests run
    fn existing_file() -> String {
        format!("{}/Cargo.toml", env!("CARGO_MANIFEST_DIR"))
    }
    fn existing_dir() -> String {
        format!("{}/src/", env!("CARGO_MANIFEST_DIR"))
    }

    /// the command line of nf-core/oncoanalyser (modules/local/esvee), with every path an existing one;
    /// `with_reference` false is its tumor-only form
    fn oncoanalyser_args(with_reference: bool) -> Vec<String> {
        let f = existing_file();
        let d = existing_dir();
        let mut args = strings(&["-tumor", "TUMOR_ID", "-tumor_bam", "tumor.bam"]);
        if with_reference {
            args.extend(strings(&["-reference", "NORMAL_ID", "-reference_bam", "normal.bam"]));
        }
        args.extend(strings(&["-esvee_prep_dir", &d, "-ref_genome", &f, "-ref_genome_version", "38"]));
        for name in ["-known_hotspot_file", "-pon_sgl_file", "-pon_sv_file", "-repeat_mask_file", "-unmap_regions", "-decoy_genome", "-saga_fasta"] {
            args.extend(strings(&[name, &f]));
        }
        args.extend(strings(&["-sequencing_type", "ILLUMINA", "-bamtool", &f]));
        args.extend(strings(&["-write_types", "PREP_JUNCTION;PREP_BAM;FRAGMENT_LENGTH_DIST;JUNC_ASSEMBLY;PHASED_ASSEMBLY;ALIGNMENT;BREAKEND;VCF"]));
        args.extend(strings(&["-threads", "8", "-output_dir", &d]));
        args
    }

    fn parse(args: &[String]) -> Options {
        match Options::parse(args) {
            Ok(o) => o,
            Err(e) => panic!("{e:?}"),
        }
    }

    /// the value that follows `-name` in a step's command line
    fn arg<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
        args.iter().position(|a| a == &format!("-{name}")).and_then(|i| args.get(i + 1)).map(String::as_str)
    }

    fn names(args: &[String]) -> Vec<&str> {
        args.iter().filter_map(|a| a.strip_prefix('-')).filter(|a| !a.is_empty() && !a.starts_with(|c: char| c.is_ascii_digit())).collect()
    }

    #[test]
    fn options_are_those_of_java() {
        // the output of `java -cp esvee.jar com.hartwig.hmftools.esvee.EsveeApplication -help` (ESVEE 2.1)
        let expected = include_str!("testdata/esvee_application_help.txt");
        let mut text = version_line();
        text.push('\n');
        for line in item_lines() {
            text.push_str(&line);
            text.push('\n');
        }
        assert_eq!(text, expected);
    }

    #[test]
    fn each_step_accepts_what_it_is_given() {
        let accepted = [
            (PREP, crate::prep::config::accepted_options()),
            (ASSEMBLY, crate::assembly::config::accepted_options()),
            (DEPTH, crate::depth::accepted_options()),
            (CALLER, crate::caller::accepted_options()),
        ];
        for (step, options) in accepted {
            for o in OPTIONS.iter().filter(|o| o.steps & step != 0) {
                let found = options.iter().find(|(name, _)| *name == o.name);
                assert!(found.is_some(), "step {step} does not accept -{}", o.name);
                assert_eq!(found.unwrap().1, o.kind == Kind::Flag, "-{} is a flag in one and not the other (step {step})", o.name);
            }
            // the sample and BAM lists
            if step == PREP || step == DEPTH {
                for name in ["sample", "bam_file"] {
                    assert!(options.iter().any(|(n, flag)| *n == name && !flag), "step {step} does not accept -{name}");
                }
            }
        }
    }

    #[test]
    fn unknown_and_misplaced_arguments() {
        let base = oncoanalyser_args(true);
        let with = |extra: &[&str]| {
            let mut a = base.clone();
            a.extend(strings(extra));
            Options::parse(&a).err()
        };
        assert_eq!(with(&["-no_such_option", "1"]), Some(ParseError::IncorrectArgument("unregistered config item: -no_such_option".into())));
        // the Rust-only options are taken off the command line in main before this
        assert!(matches!(with(&["-java_compat"]), Some(ParseError::IncorrectArgument(_))));
        assert_eq!(with(&["stray"]), Some(ParseError::IncorrectArgument("expecting config name but encountered invalid argument: stray".into())));
        assert_eq!(with(&["-output_id", "-x"]), Some(ParseError::IncorrectArgument("config item(output_id) has invalid argument: -x".into())));
        // a flag takes no value, so a word after it is out of place
        assert!(matches!(with(&["-no_cache_bam", "true"]), Some(ParseError::IncorrectArgument(_))));
        // a number may be negative
        assert_eq!(with(&["-max_frag_length_override", "-1"]), None);
        assert_eq!(with(&["-perf_log_time", "-0.5"]), None);
        // an option at the end without a value is ignored
        assert_eq!(with(&["-output_id"]), None);
        assert!(!parse(&{
            let mut a = base.clone();
            a.push("-output_id".into());
            a
        })
        .has("output_id"));
    }

    #[test]
    fn values_are_checked_as_java_checks_them() {
        assert_eq!(
            Options::parse(&strings(&["-tumor", "T"])).err(),
            Some(ParseError::Invalid(vec![
                "missing required config: ref_genome_version: Ref genome version, 37 or 38".into(),
                "missing required config: ref_genome: Path to reference genome fasta files".into(),
            ]))
        );
        let base = oncoanalyser_args(true);
        let with = |extra: &[&str]| {
            let mut a = base.clone();
            a.extend(strings(extra));
            Options::parse(&a).err()
        };
        assert_eq!(with(&["-pon_sv_file", "/no/such/file"]), Some(ParseError::Invalid(vec!["invalid path for config: pon_sv_file = /no/such/file".into()])));
        assert_eq!(with(&["-threads", "8x"]), Some(ParseError::Invalid(vec!["invalid type for config: threads: INTEGER = 8x".into()])));
        assert_eq!(with(&["-vaf_cap", "x"]), Some(ParseError::Invalid(vec!["invalid type for config: vaf_cap: DECIMAL = x".into()])));
        // every path of a list must exist; a wildcard is not checked
        let two = format!("{},/no/such/file", existing_file());
        assert!(matches!(with(&["-junction_file", &two]), Some(ParseError::Invalid(_))));
        assert_eq!(with(&["-input_vcf", "/no/such/*.vcf.gz"]), None);
        // the BAM options of -tumor and -reference are text, not checked as paths (prep finds a missing BAM)
        assert_eq!(with(&["-tumor_bam", "/no/such/file"]), None);
        // the prep folder must exist before the run starts
        assert!(matches!(with(&["-esvee_prep_dir", "/no/such/dir/"]), Some(ParseError::Invalid(_))));
    }

    #[test]
    fn version_anywhere() {
        assert!(wants_version(&strings(&["-version"])));
        assert!(wants_version(&strings(&["-tumor", "T", "-version"])));
        assert!(!wants_version(&strings(&["-tumor", "version"])));
        assert_eq!(version_line(), "Esvee version 2.1");
    }

    #[test]
    fn stages() {
        let base = oncoanalyser_args(true);
        assert_eq!(parse(&base).stages().unwrap(), vec![Stage::Prep, Stage::Assembly, Stage::RefDepth, Stage::Calling]);
        let with = |value: &str| {
            let mut a = base.clone();
            a.extend(strings(&["-stages", value]));
            parse(&a).stages()
        };
        assert_eq!(with("CALLING;PREP").unwrap(), vec![Stage::Calling, Stage::Prep]);
        assert_eq!(with("REF_DEPTH").unwrap(), vec![Stage::RefDepth]);
        assert!(with("PREP;").is_err());
        assert!(with("prep").is_err());
        assert!(with("DEPTH").is_err());
    }

    #[test]
    fn log_level() {
        let base = oncoanalyser_args(true);
        let with = |value: &str| {
            let mut a = base.clone();
            a.extend(strings(&["-log_level", value]));
            check_log_level(&parse(&a))
        };
        assert!(with("info").is_ok());
        assert!(with(" WARN ").is_ok());
        assert_eq!(with("VERBOSE").unwrap_err().to_string(), "java.lang.IllegalArgumentException: Unknown level constant [VERBOSE].");
    }

    #[test]
    fn samples_and_bams() {
        let o = parse(&oncoanalyser_args(true));
        assert_eq!(o.sample_ids(), ["TUMOR_ID", "NORMAL_ID"]);
        assert_eq!(o.bam_files(), ["tumor.bam", "normal.bam"]);
        let o = parse(&oncoanalyser_args(false));
        assert_eq!(o.sample_ids(), ["TUMOR_ID"]);
        assert_eq!(o.bam_files(), ["tumor.bam"]);
        // -sample and -bam_file take precedence for prep and depth
        let mut a = oncoanalyser_args(true);
        a.extend(strings(&["-sample", "A,B,", "-bam_file", &format!("{0},{0}", existing_file())]));
        let o = parse(&a);
        assert_eq!(o.sample_ids(), ["A", "B"]);
        assert_eq!(o.bam_files().len(), 2);
        // germline only
        let o = parse(&strings(&["-reference", "N", "-reference_bam", "n.bam", "-ref_genome", &existing_file(), "-ref_genome_version", "38"]));
        assert_eq!(o.sample_ids(), ["N"]);
        assert_eq!(o.bam_files(), ["n.bam"]);
        assert_eq!(java_split("", ','), [""]);
        assert!(java_split(",", ',').is_empty());
        assert_eq!(java_split(",a,,", ','), ["", "a"]);
    }

    #[test]
    fn prep_command_line() {
        let o = parse(&oncoanalyser_args(true));
        let args = o.step_args(PREP);
        assert_eq!(arg(&args, "sample"), Some("TUMOR_ID,NORMAL_ID"));
        assert_eq!(arg(&args, "bam_file"), Some("tumor.bam,normal.bam"));
        assert_eq!(
            names(&args),
            ["sample", "bam_file", "ref_genome_version", "ref_genome", "known_hotspot_file", "saga_fasta", "sequencing_type", "write_types", "bamtool", "threads", "output_dir"]
        );
        let o = parse(&oncoanalyser_args(false));
        let args = o.step_args(PREP);
        assert_eq!(arg(&args, "sample"), Some("TUMOR_ID"));
        assert_eq!(arg(&args, "bam_file"), Some("tumor.bam"));
    }

    #[test]
    fn assembly_command_line() {
        let o = parse(&oncoanalyser_args(true));
        let args = o.step_args(ASSEMBLY);
        // as a sub-routine: the input BAMs are not read, the prep BAMs are found in the prep folder
        assert_eq!(
            names(&args),
            [
                "ref_genome_version", "ref_genome", "saga_fasta", "sequencing_type", "write_types", "bamtool", "threads", "output_dir", "tumor", "reference",
                "esvee_prep_dir", "decoy_genome"
            ]
        );
        assert_eq!(arg(&args, "tumor"), Some("TUMOR_ID"));
        assert_eq!(arg(&args, "reference"), Some("NORMAL_ID"));
        let o = parse(&oncoanalyser_args(false));
        assert!(arg(&o.step_args(ASSEMBLY), "reference").is_none());
    }

    #[test]
    fn depth_command_line() {
        let o = parse(&oncoanalyser_args(true));
        let args = o.step_args(DEPTH);
        // the input BAMs, not the prep BAMs; the input VCF is found from the output folder and first sample
        assert_eq!(names(&args), ["sample", "bam_file", "ref_genome_version", "ref_genome", "threads", "output_dir", "unmap_regions"]);
        assert_eq!(arg(&args, "bam_file"), Some("tumor.bam,normal.bam"));
    }

    #[test]
    fn caller_command_line() {
        let o = parse(&oncoanalyser_args(true));
        let args = o.step_args(CALLER);
        assert_eq!(
            names(&args),
            [
                "ref_genome_version", "known_hotspot_file", "sequencing_type", "output_dir", "tumor", "reference", "esvee_prep_dir", "pon_sv_file", "pon_sgl_file",
                "repeat_mask_file"
            ]
        );
        // -input_vcf, when given, is read by depth and by the caller alike
        let mut a = oncoanalyser_args(true);
        a.extend(strings(&["-input_vcf", &existing_file()]));
        let o = parse(&a);
        assert_eq!(arg(&o.step_args(DEPTH), "input_vcf"), arg(&o.step_args(CALLER), "input_vcf"));
        assert!(arg(&o.step_args(CALLER), "input_vcf").is_some());
    }

    #[test]
    fn shared_options_reach_every_step_that_reads_them() {
        let mut a = oncoanalyser_args(true);
        a.extend(strings(&["-specific_chr", "chr1;chr2", "-log_level", "DEBUG", "-output_id", "run1", "-perf_log_time", "2.5", "-no_cache_bam"]));
        let o = parse(&a);
        for step in [PREP, ASSEMBLY, DEPTH, CALLER] {
            let args = o.step_args(step);
            assert_eq!(arg(&args, "specific_chr"), Some("chr1;chr2"));
            assert_eq!(arg(&args, "log_level"), Some("DEBUG"));
            assert_eq!(arg(&args, "output_id"), Some("run1"));
            assert_eq!(arg(&args, "perf_log_time").is_some(), step == ASSEMBLY || step == DEPTH);
            assert_eq!(args.iter().any(|a| a == "-no_cache_bam"), step == PREP);
            assert_eq!(arg(&args, "threads").is_some(), step != CALLER);
        }
    }
}
