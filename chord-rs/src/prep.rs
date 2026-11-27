//! ChordDataPrep, SamplePrepTask and ChordDataWriter: the mutation context counts of each sample (96 SNV
//! bins, 30 INDEL bins, 19 SV bins), written to `<sample|cohort>.chord.mutation_contexts[.<id>].tsv`.

use anyhow::{anyhow, bail, Context, Result};
use std::io::Write;
use std::sync::OnceLock;

use hmf_common::{debug, error, info};

use crate::config::ChordConfig;
use crate::indel;
use crate::output;
use crate::refgenome::RefGenomeSource;
use crate::snv;
use crate::sv;
use crate::vcf_file;

/// MutContextCount
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MutContextCount {
    pub name: String,
    pub count: i32,
}

impl MutContextCount {
    #[cfg_attr(not(test), allow(dead_code))]
    pub fn new(name: &str, count: i32) -> MutContextCount {
        MutContextCount { name: name.to_string(), count }
    }
}

/// The reference genome of a prep: RefGenomeSource.loadRefGenome logs a failure and leaves it null, which
/// only fails once a base is read (a VCF without SNVs or INDELs never reads one)
pub struct LazyRefGenome {
    path: String,
    genome: OnceLock<std::result::Result<RefGenomeSource, String>>,
}

impl LazyRefGenome {
    pub fn new(path: &str) -> LazyRefGenome {
        LazyRefGenome { path: path.to_string(), genome: OnceLock::new() }
    }

    pub fn get(&self) -> Result<&RefGenomeSource> {
        let g = self.genome.get_or_init(|| {
            RefGenomeSource::load(&self.path).map_err(|e| {
                error!("Reference file loading failed: {e:#}");
                format!("{e:#}")
            })
        });
        g.as_ref().map_err(|e| anyhow!("java.lang.NullPointerException: the reference genome is null ({e})"))
    }
}

fn write_lines(path: &str, header: &str, lines: impl Iterator<Item = String>) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("java.io.FileNotFoundException: {path}"))?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "{header}")?;
    for l in lines {
        writeln!(w, "{l}")?;
    }
    w.flush()?;
    Ok(())
}

/// SamplePrepTask.processSample: the 145 counts of one sample, in output order
pub fn process_sample(config: &ChordConfig, sample_index: usize) -> Result<Vec<MutContextCount>> {
    let sample_id = &config.sample_ids[sample_index];
    let log_prefix = if config.is_multi_sample() {
        format!("{}/{}: sample({}): ", sample_index + 1, config.sample_ids.len(), sample_id)
    } else {
        String::new()
    };
    let show_coarse_progress = config.is_multi_sample()
        && !hmf_common::log::debug_enabled()
        && matches!(config.log_level.as_str(), "INFO" | "WARN" | "ERROR" | "FATAL" | "OFF")
        && (sample_index < 10 || sample_index % 1000 == 0);
    if show_coarse_progress {
        info!("{log_prefix}Running ChordDataPrep");
    }

    // the preps are made first: SnvPrep and IndelPrep each need the reference genome
    let Some(ref_path) = &config.ref_genome_file else {
        bail!("java.lang.IllegalStateException: Ref genome must be provided to run SnvPrep");
    };
    let ref_genome = LazyRefGenome::new(ref_path);
    let output_dir = config.output_dir_str();

    let mut counts = Vec::with_capacity(145);

    // SnvPrep, then IndelPrep: each reads the SNV/INDEL VCF
    let small_variants = (|| -> Result<Vec<hmf_common::variant_context::VcfRecord>> {
        let path = config.snv_indel_vcf_file(sample_id)?;
        vcf_file::load_variants(&path, config.include_non_pass, &log_prefix)
    })();
    let snv_counts = (|| -> Result<Vec<MutContextCount>> {
        debug!("{log_prefix}Running SnvPrep - counting trinucleotide contexts");
        let variants = small_variants.as_ref().map_err(|e| anyhow!("{e:#}"))?;
        let snvs = snv::load_snvs(variants, &log_prefix);
        let details_path = format!("{output_dir}/{sample_id}{}", snv::SNV_DETAILS_FILE_SUFFIX);
        snv::count_mutation_contexts(&snvs, &ref_genome, config.write_detailed_files, &log_prefix, |details| {
            debug!("{log_prefix}Writing SNV details to: {details_path}");
            write_lines(&details_path, &snv::SnvDetails::header(), details.iter().map(|d| d.line()))
        })
    })();
    match snv_counts {
        Ok(c) => counts.extend(c),
        Err(e) => {
            error!("{log_prefix}SnvPrep failed: {e:#}");
            bail!("java.lang.NullPointerException: Cannot invoke \"java.util.Collection.toArray()\" because \"c\" is null");
        }
    }

    let indel_counts = (|| -> Result<Vec<MutContextCount>> {
        debug!("{log_prefix}Running IndelPrep - counting microhomology and repeat contexts");
        let variants = small_variants.as_ref().map_err(|e| anyhow!("{e:#}"))?;
        let indels = indel::load_indels(variants, &log_prefix);
        let details_path = format!("{output_dir}/{sample_id}{}", indel::INDEL_DETAILS_FILE_SUFFIX);
        indel::count_mutation_contexts(&indels, &ref_genome, config.write_detailed_files, &log_prefix, |details| {
            debug!("{log_prefix}Writing indel details to: {details_path}");
            write_lines(&details_path, &indel::IndelDetails::header(), details.iter().map(|d| d.line(sample_id)))
        })
    })();
    match indel_counts {
        Ok(c) => counts.extend(c),
        Err(e) => {
            error!("{log_prefix}IndelPrep failed: {e:#}");
            bail!("java.lang.NullPointerException: Cannot invoke \"java.util.Collection.toArray()\" because \"c\" is null");
        }
    }

    let sv_counts = (|| -> Result<Vec<MutContextCount>> {
        debug!("{log_prefix}Running SvPrep - counting SVs by type and length");
        let path = config.sv_vcf_file(sample_id)?;
        let variants = vcf_file::load_variants(&path, config.include_non_pass, &log_prefix)?;
        let svs = sv::load_svs(&variants, &log_prefix)?;
        let c = sv::count_mutation_contexts(&svs, &log_prefix)?;
        if config.write_detailed_files {
            let details_path = format!("{output_dir}/{sample_id}{}", sv::SV_DETAILS_FILE_SUFFIX);
            debug!("{log_prefix}Writing indel details to: {details_path}");
            write_lines(&details_path, &sv::StructuralVariant::header(), svs.iter().map(|s| s.line()))?;
        }
        Ok(c)
    })();
    match sv_counts {
        Ok(c) => counts.extend(c),
        Err(e) => {
            error!("{log_prefix}SvPrep failed: {e:#}");
            bail!("java.lang.NullPointerException: Cannot invoke \"java.util.Collection.toArray()\" because \"c\" is null");
        }
    }
    Ok(counts)
}

/// ChordDataWriter: the header from the first sample's contexts, one line per sample
pub fn mutation_contexts_text(rows: &[(String, Vec<MutContextCount>)]) -> String {
    let mut out = String::new();
    if let Some((_, first)) = rows.first() {
        out.push_str("sample_id");
        for c in first {
            out.push('\t');
            out.push_str(&c.name);
        }
        out.push('\n');
    }
    for (id, counts) in rows {
        out.push_str(id);
        for c in counts {
            out.push('\t');
            out.push_str(&c.count.to_string());
        }
        out.push('\n');
    }
    out
}

/// ChordDataPrep.run: the counts of every sample, in sample order, with the mutation contexts file written
pub fn run(config: &ChordConfig) -> Result<Vec<(String, Vec<MutContextCount>)>> {
    info!("Starting CHORD feature extraction");
    let started = std::time::Instant::now();
    if config.sample_ids.is_empty() {
        error!("No sample ID(s) loaded");
        bail!("no samples");
    }
    let output_file = output::mutation_contexts_file(config);

    let rows: Vec<(String, Vec<MutContextCount>)> = if config.is_single_sample() {
        let counts = process_sample(config, 0).map_err(|e| {
            error!("CHORD feature extraction failed: {e:#}");
            e
        })?;
        vec![(config.sample_ids[0].clone(), counts)]
    } else {
        // TaskExecutor.executeTasks: one task per sample; at most `threads` at a time
        let n = config.sample_ids.len();
        let threads = if config.threads <= 1 { 1 } else { config.threads as usize };
        let results = run_tasks(n, threads, |i| process_sample(config, i));
        let mut rows = Vec::with_capacity(n);
        let mut failed = false;
        for (i, r) in results.into_iter().enumerate() {
            match r {
                Some(Ok(c)) => rows.push((config.sample_ids[i].clone(), c)),
                Some(Err(e)) => {
                    if !failed {
                        error!("task execution error: {e:#}");
                    }
                    failed = true;
                }
                None => failed = true,
            }
        }
        if failed {
            // Java writes the file from the samples' results and fails on the first missing one, before
            // its buffered writer has written anything for a small cohort
            info!("Writing mutation context counts to: {output_file}");
            std::fs::File::create(&output_file).ok();
            error!("CHORD feature extraction failed: java.lang.NullPointerException");
            bail!("a sample failed");
        }
        // a sample id given twice: Java's map keeps one result per id, written for each occurrence
        rows
    };

    info!("Writing mutation context counts to: {output_file}");
    let text = mutation_contexts_text(&rows);
    std::fs::write(&output_file, text).map_err(|e| {
        error!("CHORD feature extraction failed: java.io.FileNotFoundException: {output_file} ({e})");
        anyhow!("cannot write {output_file}: {e}")
    })?;
    info!("Completed CHORD feature extraction, mins({})", hmf_common::fmt::java_fixed(started.elapsed().as_secs_f64() / 60.0, 3));
    Ok(rows)
}

/// Runs `n` tasks on up to `threads` threads; with one thread they run in order and stop at the first
/// failure, as TaskExecutor.executeTasks does (None for the tasks not run)
pub fn run_tasks<T: Send>(n: usize, threads: usize, f: impl Fn(usize) -> Result<T> + Sync) -> Vec<Option<Result<T>>> {
    let mut results: Vec<Option<Result<T>>> = (0..n).map(|_| None).collect();
    if threads <= 1 {
        for (i, slot) in results.iter_mut().enumerate() {
            let r = f(i);
            let failed = r.is_err();
            *slot = Some(r);
            if failed {
                break;
            }
        }
        return results;
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<std::sync::Mutex<Option<Result<T>>>> = (0..n).map(|_| std::sync::Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..threads.min(n) {
            s.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= n {
                    break;
                }
                let r = f(i);
                *slots[i].lock().unwrap() = Some(r);
            });
        }
    });
    slots.into_iter().map(|m| m.into_inner().unwrap()).collect()
}
