//! PAVE in Rust: a port of hmftools PAVE 2.0 (see ../pave).
//!
//! Usage: `pave-rs <the Java options>`: annotates a SAGE VCF with each variant's transcript and protein
//! impacts, its canonical and worst impact, ClinVar, PON, gnomAD, mappability and blacklist annotations
//! and PAVE's filters, writing the same VCF and tabix index as the Java PaveApplication.
//! `-java_compat` reproduces the Java faults the port corrects (REPORT.md).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod annotation;
mod compat;
mod config;
mod constants;
mod gene_cache;
mod impact;
mod input;
mod output;
mod task;
#[cfg(test)]
mod tests;
mod variant;

#[allow(unused_imports)]
use hmf_common::{debug, error, info, warn};

/// Log4j's TRACE level, which the port does not write
#[macro_export]
macro_rules! trace {
    ($($arg:tt)*) => {
        if false {
            let _ = format!($($arg)*);
        }
    };
}

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::sync::Mutex;

use hmf_common::regions::HUMAN_CHROMOSOMES;
use hmf_common::vcf::VcfFileWriter;

use crate::annotation::blacklist::Blacklistings;
use crate::annotation::clinvar::ClinvarData;
use crate::annotation::gnomad::{self, GnomadData, GnomadSource};
use crate::annotation::mappability::MappabilityData;
use crate::annotation::pon::{PonData, PonFilters};
use crate::annotation::reportability::Reportability;
use crate::config::{ConfigBuilder, PaveConfig, Parsed};
use crate::constants::*;
use crate::gene_cache::GeneDataCache;
use crate::impact::RefGenome;

pub const VERSION: &str = "2.0";

/// The reference genome as RefGenomeSource reads it; each thread opens its own reader (htslib's index is
/// not shared between threads)
pub struct FastaRefGenome {
    path: String,
}

thread_local! {
    static REF_GENOME: std::cell::RefCell<Option<(String, hmf_common::refgenome::RefGenome)>> = const { std::cell::RefCell::new(None) };
}

impl FastaRefGenome {
    pub fn open(path: &str) -> Result<FastaRefGenome> {
        // checked here, so that a thread's reader opens
        hmf_common::refgenome::RefGenome::open(path)?;
        Ok(FastaRefGenome { path: path.to_string() })
    }
}

impl RefGenome for FastaRefGenome {
    fn base_string(&self, chromosome: &str, start: i32, end: i32) -> Result<Option<String>> {
        if start < 1 || end < start {
            return Ok(None);
        }
        REF_GENOME.with(|cell| {
            let mut cell = cell.borrow_mut();
            if cell.as_ref().is_none_or(|(p, _)| p != &self.path) {
                *cell = Some((self.path.clone(), hmf_common::refgenome::RefGenome::open(&self.path)?));
            }
            let genome = &cell.as_ref().unwrap().1;
            let length = genome.chromosome_length(chromosome);
            if end > length {
                bail!("Query asks for data past end of contig. Query contig {chromosome} start:{start} stop:{end} contigLength:{length}");
            }
            let bases = genome.bases(chromosome, start, end)?;
            Ok(Some(String::from_utf8_lossy(&bases).into_owned()))
        })
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(pave_main(&argv));
}

/// PaveApplication.main: the exit status
pub fn pave_main(argv: &[String]) -> i32 {
    // the JVM options of the bioconda wrapper (pave -Xmx16g ...) are accepted and ignored
    let argv: Vec<String> = argv.iter().filter(|a| !(a.starts_with("-Xm") || a.starts_with("-Xs") || a.starts_with("-XX:"))).cloned().collect();
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e}");
            return 1;
        }
    };
    let cb = match ConfigBuilder::parse(&argv) {
        Ok(Parsed::Run(cb)) => cb,
        Ok(Parsed::Exit) => return 0,
        Err(_) => return 1,
    };
    // ConfigUtils.setLogLevel
    let level = cb.value("log_level").unwrap_or_default().to_ascii_uppercase();
    if cb.has_flag("log_debug") || level == "DEBUG" || level == "TRACE" || level == "ALL" {
        hmf_common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    info!("Pave version {VERSION}");
    let config = match PaveConfig::from(&cb) {
        Ok(c) => c,
        Err(e) => {
            error!("{e:#}");
            return 1;
        }
    };
    match run(&config) {
        Ok(()) => 0,
        Err(e) => {
            error!("{e:#}");
            1
        }
    }
}

/// Runs the jobs on up to `threads` threads, in the order given; the results in that order
fn run_parallel<T: Send>(threads: usize, jobs: Vec<Box<dyn FnOnce() -> Result<T> + Send + '_>>) -> Result<Vec<T>> {
    let n = jobs.len();
    let queue: Mutex<Vec<(usize, Box<dyn FnOnce() -> Result<T> + Send + '_>)>> = Mutex::new(jobs.into_iter().enumerate().rev().collect());
    let results: Mutex<Vec<Option<Result<T>>>> = Mutex::new((0..n).map(|_| None).collect());
    std::thread::scope(|s| {
        for _ in 0..threads.max(1).min(n.max(1)) {
            s.spawn(|| loop {
                let job = queue.lock().unwrap().pop();
                let Some((i, f)) = job else { break };
                let r = f();
                results.lock().unwrap()[i] = Some(r);
            });
        }
    });
    results.into_inner().unwrap().into_iter().map(|r| r.unwrap()).collect()
}

fn run(config: &PaveConfig) -> Result<()> {
    let started = std::time::Instant::now();

    // the PON and gnomAD annotations are skipped if SAGE has applied them already
    let input_header = input::read_header(&config.vcf_file)?;
    let has_pon_annotation = input_header.has_info_line(PON_COUNT);
    let has_gnomad_annotation = input_header.has_info_line(GNOMAD_FREQ);

    let mut valid = true;
    let gnomad_source: Option<GnomadSource> = if has_gnomad_annotation {
        None
    } else if let Some(f) = &config.gnomad_freq_file {
        Some(GnomadSource::File(f.clone()))
    } else if let Some(d) = &config.gnomad_freq_dir {
        match gnomad::find_chromosome_files(d, config.ref_gen_version) {
            Ok(files) => Some(GnomadSource::Dir(files)),
            Err(e) => {
                error!("{e:#}");
                valid = false;
                Some(GnomadSource::Dir(HashMap::new()))
            }
        }
    } else {
        None
    };
    let gnomad_enabled = match &gnomad_source {
        Some(GnomadSource::File(_)) => true,
        Some(GnomadSource::Dir(m)) => !m.is_empty(),
        None => false,
    };

    let standard_pon_file = if has_pon_annotation { None } else { config.pon_file.clone().filter(|f| !f.is_empty()) };
    let pon_filters_config = config
        .pon_filters
        .clone()
        .unwrap_or_else(|| if config.ref_gen_version.is_37() { PON_FILTERS_V37.to_string() } else { PON_FILTERS_V38.to_string() });
    let pon_filters = match PonFilters::parse(&pon_filters_config) {
        Ok(f) => f,
        Err(e) => {
            error!("invalid PON filters({}): {}", pon_filters_config, e);
            error!("invalid PON annotation data");
            bail!("invalid reference data, exiting");
        }
    };
    let artefact_pon_file = config.pon_artefact_file.clone().filter(|f| !f.is_empty());

    let blacklist = Blacklistings::load(config.blacklist_bed.as_deref(), config.blacklist_vcf.as_deref())?;

    let ref_genome = FastaRefGenome::open(&config.ref_genome).context("loading ref genome")?;

    if !valid {
        error!("invalid Gnomad frequency annotation data");
        bail!("invalid reference data, exiting");
    }

    // the chromosomes processed, in HumanChromosome order
    let chromosomes: Vec<String> = HUMAN_CHROMOSOMES
        .iter()
        .map(|c| config.ref_gen_version.versioned_chromosome(c))
        .filter(|c| !(config.specific_regions.has_filters() && config.specific_regions.exclude_chromosome(c)))
        .collect();

    info!("sample({}) processing VCF file({})", config.sample_id, config.vcf_file);
    let t_read = std::time::Instant::now();
    let inputs = input::read_variants(config, &chromosomes)?;
    debug!("input VCF read in {:.2}s", t_read.elapsed().as_secs_f64());
    let t_res = std::time::Instant::now();
    let needed = &inputs.needed;

    // the annotation resources, read in parallel, largest first
    let mut mappability: Option<MappabilityData> = None;
    let mut pon: Option<PonData> = None;
    let mut artefact_pon: Option<PonData> = None;
    let mut clinvar: Option<ClinvarData> = None;
    let mut gnomad_data = GnomadData::default();
    let mut gene_cache: Option<&'static GeneDataCache> = None;
    // with 4 threads or more the two largest files are inflated on threads of their own, beside the pool
    let threaded_reading = config.threads >= 4;
    let pool_threads = if threaded_reading { config.threads - 2 } else { config.threads };
    {
        enum Res {
            Genes(Box<GeneDataCache>),
            Map(MappabilityData),
            Pon(PonData),
            Artefact(PonData),
            Clinvar(ClinvarData),
            Gnomad(GnomadData),
        }
        let mut jobs: Vec<Box<dyn FnOnce() -> Result<Res> + Send + '_>> = Vec::new();
        // the Ensembl cache (ReferenceData's first) is loaded beside the resources
        jobs.push(Box::new(move || {
            let t = std::time::Instant::now();
            let g = GeneDataCache::load(&config.ensembl_data_dir, config.ref_gen_version, config.driver_gene_panel.as_deref(), config.only_canonical)
                .context("gene data cache loading failed, exiting")?;
            debug!("Ensembl data cache loaded in {:.2}s", t.elapsed().as_secs_f64());
            Ok(Res::Genes(Box::new(g)))
        }));
        if let Some(f) = &config.mappability_bed {
            jobs.push(Box::new(move || {
                let t = std::time::Instant::now();
                let m = annotation::mappability::load(f, needed, threaded_reading)?;
                debug!("mappability loaded in {:.2}s", t.elapsed().as_secs_f64());
                Ok(Res::Map(m))
            }));
        }
        if let Some(f) = &standard_pon_file {
            jobs.push(Box::new(move || {
                let t = std::time::Instant::now();
                let p = annotation::pon::load_pon(f, Some(needed), threaded_reading)?;
                debug!("PON loaded in {:.2}s ({} entries read)", t.elapsed().as_secs_f64(), p.entries_read);
                Ok(Res::Pon(p))
            }));
        }
        match &gnomad_source {
            Some(GnomadSource::File(f)) => {
                jobs.push(Box::new(move || {
                    let mut g = GnomadData::default();
                    let n = g.load(f, None, &|c| needed.get(c).map(|v| v.to_vec()))?;
                    info!("loaded {} Gnomad frequency records from file({})", n, f);
                    Ok(Res::Gnomad(g))
                }));
            }
            Some(GnomadSource::Dir(files)) => {
                for chr in &chromosomes {
                    // a chromosome's file is read by its task in Java; one without variants is not needed
                    if !inputs.by_chr.contains_key(chr) {
                        continue;
                    }
                    let Some(f) = files.get(chr) else { continue };
                    let chr = chr.clone();
                    jobs.push(Box::new(move || {
                        let t = std::time::Instant::now();
                        let mut g = GnomadData::default();
                        let n = g.load(f, Some(&chr), &|c| needed.get(c).map(|v| v.to_vec()))?;
                        debug!("chr({}) loaded {} Gnomad frequency records in {:.2}s", chr, n, t.elapsed().as_secs_f64());
                        Ok(Res::Gnomad(g))
                    }));
                }
            }
            None => {}
        }
        if let Some(f) = &config.clinvar_vcf {
            jobs.push(Box::new(move || {
                let t = std::time::Instant::now();
                let c = annotation::clinvar::load(f, needed)?;
                debug!("Clinvar loaded in {:.2}s", t.elapsed().as_secs_f64());
                Ok(Res::Clinvar(c))
            }));
        }
        if let Some(f) = &artefact_pon_file {
            jobs.push(Box::new(move || Ok(Res::Artefact(annotation::pon::load_pon(f, None, false)?))));
        }
        for r in run_parallel(pool_threads, jobs)? {
            match r {
                Res::Genes(g) => gene_cache = Some(Box::leak(g)),
                Res::Map(m) => mappability = Some(m),
                Res::Pon(p) => pon = Some(p),
                Res::Artefact(p) => artefact_pon = Some(p),
                Res::Clinvar(c) => clinvar = Some(c),
                Res::Gnomad(g) => gnomad_data.merge(g),
            }
        }
    }

    debug!("annotation resources loaded in {:.2}s", t_res.elapsed().as_secs_f64());
    let gene_cache = gene_cache.unwrap();
    let reportability = Reportability::new(gene_cache.driver_panel());

    // PaveApplication.initialiseVcfWriter: the header is written first
    let output_vcf = config.output_vcf();
    info!("writing VCF file({})", output_vcf);
    let header_options = output::HeaderOptions {
        write_detailed: config.write_detailed,
        pon_enabled: standard_pon_file.is_some() || artefact_pon_file.is_some(),
        gnomad_enabled,
        mappability_enabled: config.mappability_bed.is_some(),
        clinvar_enabled: config.clinvar_vcf.is_some(),
        blacklist_enabled: blacklist.enabled(),
        set_reportable: config.set_reportable,
    };
    let output_header = output::output_header(&input_header, &header_options);
    let mut writer = VcfFileWriter::create(&output_vcf, &output_header)?;

    let transcript_file = if config.write_transcript_file { Some(format!("{}{}.pave.transcript.tsv", config.output_dir, config.sample_id)) } else { None };

    // VersionInfo.write
    std::fs::write(format!("{}/pave.version", config.output_dir), format!("version={VERSION}\nbuild.date={}", env!("PAVE_BUILD_DATE")))
        .with_context(|| format!("writing the version file to {}", config.output_dir))?;


    let t_proc = std::time::Instant::now();
    let shared = task::Shared {
        config,
        gene_cache,
        ref_genome: &ref_genome,
        input_header: &inputs.header,
        output_header: &output_header,
        pon: pon.as_ref(),
        pon_filters: &pon_filters,
        artefact_pon: artefact_pon.as_ref(),
        gnomad: if gnomad_enabled { Some(&gnomad_data) } else { None },
        gnomad_apply_filter: !config.gnomad_no_filter,
        mappability: mappability.as_ref(),
        clinvar: clinvar.as_ref(),
        blacklist: &blacklist,
        reportability: &reportability,
    };

    // the chromosome tasks
    let empty: std::sync::Arc<Vec<String>> = std::sync::Arc::new(Vec::new());
    let shared_ref = &shared;
    let jobs: Vec<Box<dyn FnOnce() -> Result<task::ChromosomeOutput> + Send + '_>> = chromosomes
        .iter()
        .map(|chr| {
            let lines = inputs.by_chr.get(chr).unwrap_or(&empty);
            let count = inputs.read_counts.get(chr).copied().unwrap_or(0);
            let f: Box<dyn FnOnce() -> Result<task::ChromosomeOutput> + Send + '_> = Box::new(move || {
                let out = task::run_chromosome(shared_ref, chr, lines.as_slice())?;
                info!("chr({}) complete for {} variants", hmf_common::regions::strip_chr_prefix(chr), count);
                Ok(out)
            });
            f
        })
        .collect();
    let outputs = run_parallel(config.threads, jobs)?;
    debug!("variants processed in {:.2}s", t_proc.elapsed().as_secs_f64());
    let t_write = std::time::Instant::now();

    let mut transcript_rows: Vec<String> = Vec::new();
    for out in outputs {
        for r in &out.records {
            writer.write_line(&r.chrom, r.position, r.end, &r.line)?;
        }
        transcript_rows.extend(out.transcript_rows);
    }
    writer.close()?;

    if let Some(f) = transcript_file {
        let mut text = output::transcript_header();
        text.push('\n');
        for r in &transcript_rows {
            text.push_str(r);
            text.push('\n');
        }
        std::fs::write(&f, text).with_context(|| format!("writing {f}"))?;
    }

    debug!("output written in {:.2}s", t_write.elapsed().as_secs_f64());
    info!("Pave complete, mins({})", hmf_common::fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(())
}
