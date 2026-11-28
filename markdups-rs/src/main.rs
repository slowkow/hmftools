//! markdups-rs: Picard's MarkDuplicates as GATK 4.6.1.0 runs it (Picard 3.3.0, htsjdk 4.1.3), with the
//! same command line (`markdups-rs MarkDuplicates --INPUT ... --OUTPUT ... --METRICS_FILE ...`) and the same
//! outputs under -java_compat; Picard's faults corrected by default (NOTES.md, MD1...). The `gatk`
//! wrapper of the Rust image calls it for `gatk [--java-options ...] MarkDuplicates ...`.

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod args;
mod bai;
mod bam;
mod bgzf;
mod compat;
mod input;
mod markdups;
mod metrics;
mod optical;
#[cfg(test)]
mod picard_tests;
mod sort;
mod validate;

use anyhow::{Result, bail};

use crate::bam::SortOrder;
use crate::markdups::{Config, Scoring, Tagging};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// The GATK whose MarkDuplicates this reproduces, and its Picard
pub const GATK_VERSION: &str = "4.6.1.0";
pub const PICARD_VERSION: &str = "3.3.0";

#[macro_export]
macro_rules! log_info {
    ($($arg:tt)*) => {
        eprintln!("INFO  MarkDuplicates - {}", format!($($arg)*))
    };
}

fn usage() -> String {
    format!(
        "markdups-rs {VERSION}: Picard MarkDuplicates as GATK {GATK_VERSION} runs it (Picard {PICARD_VERSION}), in Rust\n\n\
         usage: markdups-rs [--java-options '...'] MarkDuplicates --INPUT in.bam [--INPUT ...] --OUTPUT out.bam --METRICS_FILE out.metrics [GATK options]\n\
         \x20      [--threads N] [-java_compat | -java_compat_except MD1,...]\n\n\
         Threads: --threads N, else MARKDUPS_THREADS, else the CPUs available. -java_compat (or MARKDUPS_JAVA_COMPAT=1)\n\
         reproduces GATK {GATK_VERSION}'s output; without it Picard's faults are corrected (see NOTES.md)."
    )
}

/// Splits off what is not MarkDuplicates' own: --java-options (its -Dsamjdk.compression_level sets the
/// default COMPRESSION_LEVEL, as the gatk launcher passes it to Java), --threads, the tool name.
struct Outer {
    tool: Option<String>,
    threads: Option<usize>,
    compression_default: i32,
    rest: Vec<String>,
    version: bool,
}

fn split_outer(argv: &[String]) -> Result<Outer> {
    let mut o = Outer { tool: None, threads: None, compression_default: 2, rest: Vec::new(), version: false };
    let mut i = 0;
    while i < argv.len() {
        let a = &argv[i];
        match a.as_str() {
            "--java-options" => {
                let Some(v) = argv.get(i + 1) else { bail!("--java-options requires a value") };
                for opt in v.split_whitespace() {
                    if let Some(l) = opt.strip_prefix("-Dsamjdk.compression_level=") {
                        o.compression_default = l.parse().map_err(|_| anyhow::anyhow!("bad -Dsamjdk.compression_level"))?;
                    }
                }
                i += 2;
                continue;
            }
            "--threads" | "-threads" => {
                let Some(v) = argv.get(i + 1) else { bail!("--threads requires a value") };
                o.threads = Some(v.parse().map_err(|_| anyhow::anyhow!("bad --threads value {v}"))?);
                i += 2;
                continue;
            }
            "--version" if o.tool.is_none() => {
                o.version = true;
            }
            _ => {
                if o.tool.is_none() && !a.starts_with('-') {
                    o.tool = Some(a.clone());
                } else {
                    o.rest.push(a.clone());
                }
            }
        }
        i += 1;
    }
    Ok(o)
}

fn threads(given: Option<usize>) -> usize {
    if let Some(t) = given {
        return t.max(1);
    }
    if let Some(t) = std::env::var("MARKDUPS_THREADS").ok().and_then(|v| v.parse::<usize>().ok()) {
        return t.max(1);
    }
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

fn config(p: &args::Parsed, threads: usize) -> Result<Config> {
    p.check_supported()?;
    let scoring = match p.get("DUPLICATE_SCORING_STRATEGY").as_deref() {
        Some("TOTAL_MAPPED_REFERENCE_LENGTH") => Scoring::TotalMappedReferenceLength,
        Some("RANDOM") => Scoring::Random,
        _ => Scoring::SumOfBaseQualities,
    };
    let tagging = match p.get("TAGGING_POLICY").as_deref() {
        Some("OpticalOnly") => Tagging::OpticalOnly,
        Some("All") => Tagging::All,
        _ => Tagging::DontTag,
    };
    let mut assume = p.get("ASSUME_SORT_ORDER").map(|s| SortOrder::parse(s.as_bytes()).unwrap());
    if p.flag("ASSUME_SORTED") && assume.is_none() {
        assume = Some(SortOrder::Coordinate);
    }
    let compression = p.int("COMPRESSION_LEVEL") as i32;
    if !(0..=9).contains(&compression) {
        bail!("COMPRESSION_LEVEL must be between 0 and 9");
    }
    Ok(Config {
        inputs: p.list("INPUT"),
        output: p.get("OUTPUT").unwrap(),
        metrics_file: p.get("METRICS_FILE").unwrap(),
        create_index: p.flag("CREATE_INDEX"),
        remove_duplicates: p.flag("REMOVE_DUPLICATES"),
        remove_sequencing_duplicates: p.flag("REMOVE_SEQUENCING_DUPLICATES"),
        tagging,
        clear_dt: p.flag("CLEAR_DT"),
        add_pg_tag_to_reads: p.flag("ADD_PG_TAG_TO_READS"),
        program_record_id: p.get("PROGRAM_RECORD_ID"),
        program_group_version: p.get("PROGRAM_GROUP_VERSION"),
        program_group_command_line: p.get("PROGRAM_GROUP_COMMAND_LINE"),
        program_group_name: p.get("PROGRAM_GROUP_NAME"),
        comments: p.list("COMMENT"),
        read_name_regex: p.get("READ_NAME_REGEX"),
        pixel_distance: p.int("OPTICAL_DUPLICATE_PIXEL_DISTANCE") as i32,
        max_optical_set_size: p.int("MAX_OPTICAL_DUPLICATE_SET_SIZE"),
        scoring,
        assume_sort_order: assume,
        compression,
        threads,
        stringency: validate::Stringency::parse(p.get("VALIDATION_STRINGENCY").as_deref().unwrap_or("STRICT")),
        command_line: p.command_line(),
        // GATK's @PG VN under -java_compat; corrected output says what wrote it
        version: if compat::java_all() {
            format!("Version:{GATK_VERSION}")
        } else {
            format!("markdups-rs {VERSION} (MarkDuplicates of GATK {GATK_VERSION}, corrected)")
        },
    })
}

fn run() -> Result<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let argv = compat::take_options(&argv)?;
    let outer = split_outer(&argv)?;
    if outer.version {
        println!("markdups-rs {VERSION} (reproduces GATK {GATK_VERSION} MarkDuplicates, Picard {PICARD_VERSION})");
        return Ok(());
    }
    let one_step = match outer.tool.as_deref() {
        Some("MarkDuplicates") => false,
        Some("SortAndMarkDuplicates") => true,
        None => {
            eprintln!("{}", usage());
            bail!("no tool given");
        }
        Some(t) => bail!("markdups-rs runs MarkDuplicates (and SortAndMarkDuplicates), not {t}"),
    };
    if outer.rest.iter().any(|a| a == "--help" || a == "-h") {
        println!("{}", usage());
        return Ok(());
    }
    let (rest, sort_opts) = if one_step { take_sort_options(&outer.rest)? } else { (outer.rest.clone(), SortOptions::default()) };
    let mut parsed = args::parse(&rest, outer.compression_default)?;
    let threads = threads(outer.threads);
    let real_inputs = parsed.list("INPUT");
    if one_step {
        if real_inputs.len() != 1 {
            bail!("SortAndMarkDuplicates takes one --INPUT (the unsorted BAM)");
        }
        // the command line recorded as that of MarkDuplicates on the sorted BAM samtools would have written
        if let Some(name) = &sort_opts.sorted_input_name {
            parsed.given[args::index_of("INPUT")] = Some(args::Value::List(vec![name.clone()]));
        }
    }
    let mut cfg = config(&parsed, threads)?;
    cfg.inputs = real_inputs;
    if one_step {
        return run_one_step(&cfg, &parsed, &sort_opts);
    }
    let start = std::time::Instant::now();
    log_info!("markdups-rs {VERSION} ({} mode), {threads} threads", if compat::java_all() { "-java_compat" } else { "corrected" });
    let inputs = input::Inputs { paths: cfg.inputs.clone(), assume_sort_order: cfg.assume_sort_order };
    let pass1_threads = if threads <= 1 { 0 } else { threads - 1 };
    let mut dups = markdups::find_duplicates(&cfg, &inputs, pass1_threads)?;
    log_info!(
        "Read {} records. {} pairs never matched. Marking {} records as duplicates ({:.1} s).",
        dups.records,
        dups.unmatched_pairs,
        dups.num_duplicate_indices,
        start.elapsed().as_secs_f64()
    );
    let (rt, ct) = markdups::pass2_threads(threads);
    markdups::write_output(&cfg, &inputs, &mut dups, rt, ct)?;
    log_info!("Done in {:.1} s", start.elapsed().as_secs_f64());
    Ok(())
}

#[derive(Default)]
struct SortOptions {
    sorted_input_name: Option<String>,
    samtools_command_line: Option<String>,
    samtools_version: Option<String>,
    memory_mb: Option<usize>,
}

/// The options of SortAndMarkDuplicates that are not MarkDuplicates'.
fn take_sort_options(argv: &[String]) -> Result<(Vec<String>, SortOptions)> {
    let mut o = SortOptions::default();
    let mut rest = Vec::new();
    let mut i = 0;
    while i < argv.len() {
        let slot = match argv[i].as_str() {
            "--SORTED_INPUT_NAME" => Some(&mut o.sorted_input_name),
            "--SAMTOOLS_SORT_COMMAND_LINE" => Some(&mut o.samtools_command_line),
            "--SAMTOOLS_VERSION" => Some(&mut o.samtools_version),
            _ => None,
        };
        if let Some(slot) = slot {
            let Some(v) = argv.get(i + 1) else { bail!("{} requires a value", argv[i]) };
            *slot = Some(v.clone());
            i += 2;
            continue;
        }
        if argv[i] == "--SORT_MEMORY_MB" {
            let Some(v) = argv.get(i + 1) else { bail!("--SORT_MEMORY_MB requires a value") };
            o.memory_mb = Some(v.parse().map_err(|_| anyhow::anyhow!("bad --SORT_MEMORY_MB {v}"))?);
            i += 2;
            continue;
        }
        rest.push(argv[i].clone());
        i += 1;
    }
    Ok((rest, o))
}

/// samtools sort and MarkDuplicates in one: sort into runs, then MarkDuplicates' two passes over them.
fn run_one_step(cfg: &Config, parsed: &args::Parsed, o: &SortOptions) -> Result<()> {
    let start = std::time::Instant::now();
    let threads = cfg.threads;
    log_info!("markdups-rs {VERSION} SortAndMarkDuplicates ({} mode), {threads} threads", if compat::java_all() { "-java_compat" } else { "corrected" });
    let tmp = parsed.list("TMP_DIR").first().cloned().unwrap_or_else(|| std::env::temp_dir().to_string_lossy().to_string());
    let memory = o.memory_mb.unwrap_or(2048) << 20;
    let (rt, ct) = markdups::pass2_threads(threads);
    let sorted = sort::sort_input(&cfg.inputs[0], &tmp, memory, rt, ct, o.samtools_command_line.as_deref(), o.samtools_version.as_deref().unwrap_or("1.21"))?;
    log_info!("Sorted in {:.1} s", start.elapsed().as_secs_f64());
    let pass1_threads = if threads <= 1 { 0 } else { threads - 1 };
    let mut dups = markdups::find_duplicates(cfg, &sorted, pass1_threads)?;
    log_info!("Marking {} records as duplicates ({:.1} s).", dups.num_duplicate_indices, start.elapsed().as_secs_f64());
    markdups::write_output(cfg, &sorted, &mut dups, rt, ct)?;
    log_info!("Done in {:.1} s", start.elapsed().as_secs_f64());
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("ERROR markdups-rs: {e:#}");
        std::process::exit(2);
    }
}
