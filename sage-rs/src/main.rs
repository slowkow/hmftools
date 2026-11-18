mod append;
mod bamio;
mod candidate;
mod config;
mod constants;
mod dedup;
mod evidence;
mod evidence_stage;
mod filters;
mod fmt;
mod matcher;
mod phase;
mod quality;
mod read_context;
mod refdata;
mod refseq;
mod repeat;
mod select;
mod stats;
mod sync;
mod tinc;
mod types;
mod ultima;
mod vcfout;

#[cfg(test)]
mod tests;

use anyhow::{anyhow, Result};
use candidate::*;
use config::Config;
use refdata::{RefData, RefDataPaths};
use refseq::{RefGenome, RefSequence};
use std::collections::BTreeMap;
use std::io::Write;
use std::sync::{Condvar, Mutex};
use types::*;

/// ChromosomePartition
pub(crate) fn partition(cfg: &Config, chr: &str, chr_length: i32) -> Vec<BaseRegion> {
    let threads = cfg.threads.min(4) as i32;
    let mut out = Vec::new();
    let regions: Vec<BaseRegion> =
        cfg.specific_regions.iter().filter(|(c, _, _)| c == chr).map(|(_, s, e)| BaseRegion::new(*s, *e)).collect();
    if !cfg.specific_regions.is_empty() {
        let total: i32 = regions.iter().map(|r| r.base_length()).sum();
        let dynamic = total / threads + 1;
        for region in regions {
            if region.base_length() < cfg.region_slice_size {
                out.push(region);
                continue;
            }
            let slice = dynamic.min(cfg.region_slice_size);
            let mut i = 0;
            loop {
                let start = region.start + i * slice;
                let mut end = start + slice - 1;
                if end > region.end || ((region.end - end) as f64) < 0.5 * slice as f64 {
                    end = region.end;
                }
                out.push(BaseRegion::new(start, end));
                if end >= region.end {
                    break;
                }
                i += 1;
            }
        }
        return out;
    }
    let slice = (chr_length / threads + 1).min(cfg.region_slice_size);
    let mut i = 0;
    loop {
        let start = 1 + i * slice;
        let end = (start + slice - 1).min(chr_length);
        out.push(BaseRegion::new(start, end));
        if end >= chr_length {
            break;
        }
        i += 1;
    }
    out
}

fn overlapping(regions: Option<&Vec<BaseRegion>>, r: BaseRegion) -> Vec<BaseRegion> {
    regions
        .map(|l| l.iter().filter(|x| positions_overlap(r.start, r.end, x.start, x.end)).cloned().collect())
        .unwrap_or_default()
}

/// SageCommon.setReadLength: the maximum read length over the first 1000 reads of a sample region.
fn sample_read_length(cfg: &Config, ref_data: &RefData, bam: &mut bamio::BamSlicer) -> Result<i32> {
    let mut regions: Vec<(String, BaseRegion)> = Vec::new();
    if !cfg.specific_regions.is_empty() {
        regions.extend(cfg.specific_regions.iter().map(|(c, s, e)| (c.clone(), BaseRegion::new(*s, *e))));
    } else {
        // Java iterates a HashMap keyed by chromosome here, so its choice of region is arbitrary
        let mut chrs: Vec<&String> = ref_data.panel.keys().collect();
        chrs.sort_by_key(|c| chromosome_rank(c).unwrap_or(99));
        for chr in chrs {
            for r in &ref_data.panel[chr] {
                let mut r = *r;
                if r.length() == 1 {
                    r = BaseRegion::new((r.start - 500).max(0), r.end + 500);
                }
                regions.push((chr.clone(), r));
            }
        }
    }
    for (chr, region) in regions {
        let mut count = 0;
        let mut max_len = 0;
        bam.slice_ref(&chr, &[BaseRegion::new(region.start.max(1), region.end)], false, |read| {
            if count < 1000 {
                count += 1;
                max_len = max_len.max(read.len());
            }
        })?;
        if count > 0 && max_len > 0 {
            return Ok(max_len);
        }
    }
    Ok(constants::DEFAULT_READ_LENGTH)
}

/// SageApplication: UltimaUtils.setMaxRawQual(bqrCache.maxRawQual()). BqrCache.loadBqrFiles overwrites
/// its maximum for each sample while iterating a HashMap keyed by sample id, so the value is that of
/// the sample which comes last in that map (reference ids are inserted first), and 0 with -skip_bqr.
pub(crate) fn ultima_max_raw_qual(cfg: &Config, bqr: &[quality::BqrMap]) -> i8 {
    if cfg.skip_bqr {
        return 0;
    }
    // `bqr` is ordered tumor ids then reference ids
    let tumors = cfg.tumor_ids.len();
    let mut entries: Vec<(&String, usize)> = Vec::new();
    for (id, index) in cfg.reference_ids.iter().enumerate().map(|(i, id)| (id, tumors + i)).chain(cfg.tumor_ids.iter().zip(0..)) {
        match entries.iter_mut().find(|(k, _)| *k == id) {
            Some(e) => e.1 = index,
            None => entries.push((id, index)),
        }
    }
    // the table doubles from 16 while it is more than three quarters full
    let mut capacity = 16usize;
    while entries.len() * 4 > capacity * 3 {
        capacity *= 2;
    }
    let bucket = |id: &str| {
        let h = sync::java_string_hash(id) as u32;
        ((h ^ (h >> 16)) as usize) & (capacity - 1)
    };
    // stable: entries sharing a bucket stay in insertion order
    entries.sort_by_key(|(id, _)| bucket(id));
    entries.last().map_or(0, |(_, index)| bqr[*index].max_raw_qual())
}

/// State shared by every region task.
struct Shared<'a> {
    cfg: &'a Config,
    ref_data: &'a RefData,
    bqr: &'a [quality::BqrMap],
    msi: &'a [quality::MsiJitter],
}

/// Per-thread file handles.
struct Worker {
    genome: RefGenome,
    /// one per tumor sample, then one per reference sample, as in the configuration
    tumor_bams: Vec<bamio::BamSlicer>,
    ref_bams: Vec<bamio::BamSlicer>,
}

impl Worker {
    fn open(cfg: &Config) -> Result<Worker> {
        let open_all = |paths: &[String]| -> Result<Vec<bamio::BamSlicer>> {
            paths.iter().map(|p| bamio::BamSlicer::open(p, 1)).collect()
        };
        Ok(Worker {
            genome: RefGenome::open(&cfg.ref_genome)?,
            tumor_bams: open_all(&cfg.tumor_bams)?,
            ref_bams: open_all(&cfg.reference_bams)?,
        })
    }
}

/// The parallel part of a RegionTask: candidates, evidence and soft filters.
struct RegionOutput {
    candidate_lines: Vec<String>,
    evidence_lines: Vec<String>,
    variants: Vec<filters::SageVariant>,
    phaser: phase::Phaser,
    /// candidate index -> index into `variants`
    variant_of: Vec<Option<usize>>,
    passing: Vec<bool>,
    valid: Vec<bool>,
    hotspots: Vec<SimpleVariant>,
}

fn process_region(sh: &Shared, w: &mut Worker, chr: &str, region: BaseRegion) -> Result<Option<RegionOutput>> {
    let cfg = sh.cfg;
    let ref_data = sh.ref_data;
    let panel = overlapping(ref_data.panel.get(chr), region);
    if cfg.panel_only && panel.is_empty() {
        return Ok(None);
    }
    let high_conf = overlapping(ref_data.high_confidence.get(chr), region);
    let hotspots: Vec<SimpleVariant> = ref_data
        .hotspots
        .get(chr)
        .map(|l| l.iter().filter(|h| region.contains(h.pos)).cloned().collect())
        .unwrap_or_default();

    let ref_seq = RefSequence::new(&w.genome, chr, region.start, region.end)?;
    let mut finder = CandidateFinder::new(cfg, chr, region, &ref_seq, &hotspots, &panel, ref_data);
    // CandidateStage: in panel-only mode only the panel parts of the region are read
    let slice_regions: Vec<BaseRegion> = if cfg.panel_only {
        panel.iter().map(|p| BaseRegion::new(p.start.max(region.start), p.end.min(region.end))).collect()
    } else {
        vec![region]
    };
    // CandidateStage.findCandidates: each tumor sample is read with its own RefContextCache
    w.tumor_bams[0].slice_ref(chr, &slice_regions, true, |read| finder.process_read(read))?;
    let candidates = if w.tumor_bams.len() == 1 {
        build_candidates(finder.finish(), &hotspots, &panel, &high_conf)
    } else {
        let mut sample_candidates = vec![finder.finish()];
        for bam in w.tumor_bams[1..].iter_mut() {
            let mut finder = CandidateFinder::new(cfg, chr, region, &ref_seq, &hotspots, &panel, ref_data);
            bam.slice_ref(chr, &slice_regions, true, |read| finder.process_read(read))?;
            sample_candidates.push(finder.finish());
        }
        build_candidates_of_multiple_samples(sample_candidates, &hotspots, &panel, &high_conf)
    };
    if candidates.is_empty() {
        return Ok(None);
    }
    let mut out = RegionOutput {
        candidate_lines: Vec::new(),
        evidence_lines: Vec::new(),
        variants: Vec::new(),
        phaser: phase::Phaser::new(),
        variant_of: vec![None; candidates.len()],
        passing: vec![false; candidates.len()],
        valid: vec![false; candidates.len()],
        hotspots,
    };
    if cfg.log_candidates {
        for c in &candidates {
            out.candidate_lines.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chr,
                c.position(),
                c.variant().ref_str(),
                c.variant().alt_str(),
                c.tier.name(),
                c.ctx.core_str(),
                c.ctx.read_cigar,
                c.full_match_support,
                c.core_match_support
            ));
        }
    }

    // one FragmentSync per region task, shared by the tumor and reference passes
    let mut sync = sync::FragmentSync::new();
    let tumor_count = cfg.tumor_ids.len();
    // EvidenceStage.findEvidence: phasing groups are collected for the first tumor sample only
    let mut tumor_counters: Vec<Vec<Option<evidence::ReadContextCounter>>> = Vec::with_capacity(tumor_count);
    for (s, bam) in w.tumor_bams.iter_mut().enumerate() {
        let sample_id = &cfg.tumor_ids[s];
        let mut tumor = evidence_stage::Sample {
            bam,
            bqr: &sh.bqr[s],
            msi: &sh.msi[s],
            is_reference: cfg.reference_ids.contains(sample_id),
        };
        let phaser: Option<&mut dyn evidence_stage::PhaseSink> =
            if *sample_id == cfg.tumor_ids[0] { Some(&mut out.phaser) } else { None };
        let counters =
            evidence_stage::collect_evidence(cfg, chr, &candidates, &w.genome, ref_data, &mut tumor, &mut sync, phaser)?;
        tumor_counters.push(counters.into_iter().map(Some).collect());
    }
    // ReadContextCounters.filterCandidates: a candidate is kept if any tumor sample passes the hard filters
    let keep: Vec<usize> = (0..candidates.len())
        .filter(|&i| tumor_counters.iter().any(|c| evidence_stage::passes_hard_filters(cfg, c[i].as_ref().unwrap())))
        .collect();
    if keep.is_empty() {
        // Java still reports the candidates, but never reaches the variant stage
        out.variant_of.clear();
        return Ok(Some(out));
    }
    let final_candidates: Vec<Candidate> = keep.iter().map(|&i| candidates[i].clone()).collect();
    let mut ref_counters: Vec<Vec<Option<evidence::ReadContextCounter>>> = Vec::with_capacity(w.ref_bams.len());
    for (s, bam) in w.ref_bams.iter_mut().enumerate() {
        let mut reference = evidence_stage::Sample {
            bam,
            bqr: &sh.bqr[tumor_count + s],
            msi: &sh.msi[tumor_count + s],
            is_reference: true,
        };
        let counters =
            evidence_stage::collect_evidence(cfg, chr, &final_candidates, &w.genome, ref_data, &mut reference, &mut sync, None)?;
        ref_counters.push(counters.into_iter().map(Some).collect());
    }
    if cfg.log_evidence {
        for (j, &i) in keep.iter().enumerate() {
            let c = &candidates[i];
            let mut line = format!("{}\t{}\t{}\t{}", chr, c.position(), c.variant().ref_str(), c.variant().alt_str());
            for rc in ref_counters.iter().map(|c| &c[j]).chain(tumor_counters.iter().map(|c| &c[i])) {
                line.push('\t');
                line.push_str(&vcfout::genotype(rc.as_ref().unwrap()));
            }
            out.evidence_lines.push(line);
        }
    }

    // combine reference and tumor counters into variants
    for (j, &i) in keep.iter().enumerate() {
        out.variant_of[i] = Some(j);
        let refs: Vec<evidence::ReadContextCounter> = ref_counters.iter_mut().map(|c| c[j].take().unwrap()).collect();
        let tumors: Vec<evidence::ReadContextCounter> = tumor_counters.iter_mut().map(|c| c[i].take().unwrap()).collect();
        out.variants.push(filters::SageVariant::new(final_candidates[j].clone(), refs, tumors));
    }

    filters::set_near_by_indel_status_pre_filter(&mut out.variants);
    let hla = |pos: i32| ref_data.hla_contains(chr, pos);
    let variant_filters = filters::VariantFilters { cfg, chr, hla: &hla };
    for (j, &i) in keep.iter().enumerate() {
        if !cfg.disable_soft_filter {
            variant_filters.apply_soft_filters(&mut out.variants[j]);
        }
        out.passing[i] = out.variants[j].is_passing();
        out.valid[i] = true;
    }
    filters::set_near_by_indel_status_post_filter(&mut out.variants);
    Ok(Some(out))
}

/// The sequential part of a RegionTask: phase set ids come from one global counter, so phasing,
/// dedup, the final filters and writing run in region order.
fn finalise_region(
    sh: &Shared, genome: &RefGenome, chr: &str, region: BaseRegion, mut out: RegionOutput, next_lps: &mut i32,
) -> Vec<(i32, String)> {
    let cfg = sh.cfg;
    let mut variants = std::mem::take(&mut out.variants);
    for a in out.phaser.assign_local_phase_sets(&out.passing, &out.valid, next_lps) {
        for i in a.counters {
            if let Some(j) = out.variant_of[i] {
                variants[j].tumor_counters[0].add_local_phase_set(a.lps, a.read_count, a.allocated);
            }
        }
    }
    let transcripts: Vec<refdata::TranscriptData> = sh
        .ref_data
        .transcripts
        .get(chr)
        .map(|l| l.iter().filter(|t| positions_overlap(region.start, region.end, t.trans_start, t.trans_end)).cloned().collect())
        .unwrap_or_default();
    dedup::process_variants(&mut variants, &transcripts, genome, chr);

    // finaliseResults
    let passing_sets = phase::passing_phase_sets(&variants);
    let mut final_variants: Vec<filters::SageVariant> =
        variants.into_iter().filter(|v| filters::check_final_filters(v, &passing_sets, cfg, chr)).collect();
    phase::remove_uninformative_lps(&mut final_variants, &passing_sets);
    let mut hotspot_cache = select::HotspotCache::new(&out.hotspots);
    for v in final_variants.iter_mut() {
        if v.tier() != Tier::Hotspot && v.is_passing() && hotspot_cache.is_near_hotspot(v.variant()) {
            v.near_hotspot = true;
        }
    }
    final_variants.iter().map(|v| (v.position(), vcfout::record(chr, v))).collect()
}

struct Queue {
    next: usize,
    written: usize,
    done: BTreeMap<usize, Result<Option<RegionOutput>>>,
    failed: bool,
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(|a| a.as_str()) == Some("tinc") {
        return tinc::run_standalone(&args[1..]);
    }
    if args.first().map(|a| a.as_str()) == Some("append") {
        return append::run(&args[1..]);
    }
    let mut cfg = Config::from_args(&args)?;
    let tinc_cfg = if cfg.run_tinc {
        // callerTincConfig: fail before the run rather than after it
        let exists = |p: &Option<String>| p.as_ref().is_some_and(|p| std::path::Path::new(p).exists());
        if !(exists(&cfg.gnomad_freq_file) || exists(&cfg.gnomad_freq_dir)) || !exists(&cfg.pon_file) {
            anyhow::bail!("missing PON files for TINC analyser");
        }
        Some(tinc::TincConfig {
            v38: cfg.v38,
            input_vcf: cfg.output_vcf.clone(),
            tumor_id: cfg.tumor_ids[0].clone(),
            reference_id: cfg.reference_ids[0].clone(),
            pon_file: cfg.pon_file.clone(),
            gnomad_dir: cfg.gnomad_freq_dir.clone(),
            gnomad_file: cfg.gnomad_freq_file.clone(),
            threads: cfg.threads,
            write_fit_variants: cfg.write_fit_variants,
            fit_variants_file: None,
        })
    } else {
        None
    };
    // SageCallConfig.isValid
    if cfg.reference_ids.len() != cfg.reference_bams.len() {
        anyhow::bail!("each reference sample must have matching BAM file");
    }
    if cfg.tumor_ids.len() != cfg.tumor_bams.len() {
        anyhow::bail!("Each tumor sample must have matching bam");
    }
    if cfg.tumor_ids.is_empty() {
        anyhow::bail!("At least one tumor must be supplied");
    }
    let ref_data = RefData::load(&RefDataPaths {
        ensembl_dir: cfg.ensembl_data_dir.as_deref(),
        driver_gene_panel: cfg.driver_gene_panel.as_deref(),
        panel_bed: cfg.panel_bed.as_deref(),
        hotspots: cfg.hotspots.as_deref(),
        high_confidence_bed: cfg.high_confidence_bed.as_deref(),
        v38: cfg.v38,
    })?;
    let mut main_worker = Worker::open(&cfg)?;
    if cfg.read_length == 0 {
        cfg.read_length = sample_read_length(&cfg, &ref_data, &mut main_worker.tumor_bams[0])?;
        if constants::is_ultima() {
            // SageConfig.setReadLength
            cfg.read_length = cfg.read_length.max(constants::NON_ILLUMINA_MAX_READ_LENGTH);
        }
        eprintln!("max observed read length set({})", cfg.read_length);
    }
    let cfg = cfg;

    let mut bqr = Vec::new();
    let mut msi = Vec::new();
    for sample in cfg.tumor_ids.iter().chain(cfg.reference_ids.iter()) {
        bqr.push(if cfg.skip_bqr {
            quality::BqrMap::default()
        } else {
            quality::BqrMap::load(&format!("{}/{}.redux.bqr.tsv", cfg.jitter_bqr_dir, sample))?
        });
        let dir = if cfg.skip_msi_jitter { None } else { Some(cfg.jitter_bqr_dir.as_str()) };
        msi.push(quality::MsiJitter::load(dir, sample, cfg.high_depth_mode)?);
    }
    if constants::is_ultima() {
        ultima::set_max_raw_qual(ultima_max_raw_qual(&cfg, &bqr));
    }

    let mut cand_out: Option<Box<dyn Write>> = if cfg.log_candidates {
        let path = cfg.output_vcf.replace(".vcf.gz", ".candidates.tsv");
        let mut w: Box<dyn Write> = Box::new(std::io::BufWriter::new(std::fs::File::create(path)?));
        writeln!(w, "Chromosome\tPosition\tRef\tAlt\tTier\tCoreBases\tCoreCigar\tFullMatches\tCoreMatches\tLowQualInCore")?;
        Some(w)
    } else {
        None
    };
    let mut evidence_out: Option<Box<dyn Write>> = if cfg.log_evidence {
        let path = cfg.output_vcf.replace(".vcf.gz", ".evidence.tsv");
        Some(Box::new(std::io::BufWriter::new(std::fs::File::create(path)?)))
    } else {
        None
    };

    let contigs: Vec<(String, i32)> = main_worker
        .genome
        .sequence_names()?
        .into_iter()
        .filter(|n| is_human_chromosome(n) || is_mito(n))
        .map(|n| {
            let len = main_worker.genome.chromosome_length(&n);
            (n, len)
        })
        .collect();
    let mut vcf = vcfout::VcfWriter::create(&cfg.output_vcf, &contigs, &cfg.reference_ids, &cfg.tumor_ids, cfg.run_tinc)?;

    // every region task, in the order the sequence dictionary lists the chromosomes
    let dictionary_bam = main_worker.ref_bams.first().unwrap_or(&main_worker.tumor_bams[0]);
    let mut tasks: Vec<(String, BaseRegion)> = Vec::new();
    for chr in dictionary_bam.sequence_names() {
        if !cfg.process_chromosome(&chr) {
            continue;
        }
        let chr_len = main_worker.genome.chromosome_length(&chr);
        tasks.extend(partition(&cfg, &chr, chr_len).into_iter().map(|r| (chr.clone(), r)));
    }

    let shared = Shared { cfg: &cfg, ref_data: &ref_data, bqr: &bqr, msi: &msi };
    let threads = cfg.threads.max(1);
    // workers may run this many regions ahead of the writer
    let window = threads * 4;
    let queue = Mutex::new(Queue { next: 0, written: 0, done: BTreeMap::new(), failed: false });
    let ready = Condvar::new();
    let mut next_lps = 1;

    let result: Result<()> = std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                let mut worker = Worker::open(&cfg);
                loop {
                    let index = {
                        let mut q = queue.lock().unwrap();
                        loop {
                            if q.failed || q.next >= tasks.len() {
                                return;
                            }
                            if q.next < q.written + window {
                                break;
                            }
                            q = ready.wait(q).unwrap();
                        }
                        q.next += 1;
                        q.next - 1
                    };
                    let (chr, region) = &tasks[index];
                    let res = match worker.as_mut() {
                        Ok(w) => process_region(&shared, w, chr, *region),
                        Err(e) => Err(anyhow!("failed to open inputs: {e}")),
                    };
                    queue.lock().unwrap().done.insert(index, res);
                    ready.notify_all();
                }
            });
        }

        let mut write_all = || -> Result<()> {
            for (index, (chr, region)) in tasks.iter().enumerate() {
                let res = {
                    let mut q = queue.lock().unwrap();
                    loop {
                        if let Some(r) = q.done.remove(&index) {
                            break r;
                        }
                        q = ready.wait(q).unwrap();
                    }
                };
                if let Some(out) = res? {
                    if let Some(w) = cand_out.as_mut() {
                        for l in &out.candidate_lines {
                            writeln!(w, "{l}")?;
                        }
                    }
                    if let Some(w) = evidence_out.as_mut() {
                        for l in &out.evidence_lines {
                            writeln!(w, "{l}")?;
                        }
                    }
                    if !out.variants.is_empty() {
                        let records = finalise_region(&shared, &main_worker.genome, chr, *region, out, &mut next_lps);
                        vcf.add_region(records)?;
                    }
                }
                // Java flushes the writer at the end of each chromosome
                if tasks.get(index + 1).map_or(true, |(c, _)| c != chr) {
                    vcf.flush_chromosome()?;
                }
                queue.lock().unwrap().written = index + 1;
                ready.notify_all();
            }
            Ok(())
        };
        let r = write_all();
        if r.is_err() {
            queue.lock().unwrap().failed = true;
            ready.notify_all();
        }
        r
    });
    result?;
    vcf.close()?;
    if let Some(tinc_cfg) = tinc_cfg {
        let output = if cfg.write_tinc_vcf { tinc::generate_tinc_vcf_filename(&cfg.output_vcf) } else { cfg.output_vcf.clone() };
        let mut analyser = tinc::TincAnalyser::new(tinc_cfg);
        analyser.run(true, true)?;
        analyser.write_vcf(&output)?;
    }
    Ok(())
}
