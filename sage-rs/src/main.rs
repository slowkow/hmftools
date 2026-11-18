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
mod types;
mod vcfout;

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
fn partition(cfg: &Config, chr: &str, chr_length: i32) -> Vec<BaseRegion> {
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
        bam.slice(&chr, &[BaseRegion::new(region.start.max(1), region.end)], false, |read| {
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
    tumor_bam: bamio::BamSlicer,
    ref_bam: Option<bamio::BamSlicer>,
}

impl Worker {
    fn open(cfg: &Config) -> Result<Worker> {
        Ok(Worker {
            genome: RefGenome::open(&cfg.ref_genome)?,
            tumor_bam: bamio::BamSlicer::open(&cfg.tumor_bams[0], 1)?,
            ref_bam: match cfg.reference_bams.first() {
                Some(p) => Some(bamio::BamSlicer::open(p, 1)?),
                None => None,
            },
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
    w.tumor_bam.slice(chr, &[region], true, |read| finder.process_read(&read))?;
    let candidates = build_candidates(finder.finish(), &hotspots, &panel, &high_conf);
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
    let mut tumor = evidence_stage::Sample { bam: &mut w.tumor_bam, bqr: &sh.bqr[0], msi: &sh.msi[0], is_reference: false };
    let tumor_counters = evidence_stage::collect_evidence(
        cfg, chr, &candidates, &w.genome, ref_data, &mut tumor, &mut sync, Some(&mut out.phaser),
    )?;
    let keep: Vec<usize> =
        (0..candidates.len()).filter(|&i| evidence_stage::passes_hard_filters(cfg, &tumor_counters[i])).collect();
    if keep.is_empty() {
        // Java still reports the candidates, but never reaches the variant stage
        out.variant_of.clear();
        return Ok(Some(out));
    }
    let final_candidates: Vec<Candidate> = keep.iter().map(|&i| candidates[i].clone()).collect();
    let ref_counters = match w.ref_bam.as_mut() {
        Some(bam) => {
            let mut reference = evidence_stage::Sample { bam, bqr: &sh.bqr[1], msi: &sh.msi[1], is_reference: true };
            evidence_stage::collect_evidence(cfg, chr, &final_candidates, &w.genome, ref_data, &mut reference, &mut sync, None)?
        }
        None => Vec::new(),
    };
    if cfg.log_evidence {
        for (j, &i) in keep.iter().enumerate() {
            let c = &candidates[i];
            let mut line = format!("{}\t{}\t{}\t{}", chr, c.position(), c.variant().ref_str(), c.variant().alt_str());
            if let Some(rc) = ref_counters.get(j) {
                line.push('\t');
                line.push_str(&vcfout::genotype(rc));
            }
            line.push('\t');
            line.push_str(&vcfout::genotype(&tumor_counters[i]));
            out.evidence_lines.push(line);
        }
    }

    // combine reference and tumor counters into variants
    let mut tumor_opt: Vec<Option<evidence::ReadContextCounter>> = tumor_counters.into_iter().map(Some).collect();
    let mut ref_iter = ref_counters.into_iter();
    for (j, &i) in keep.iter().enumerate() {
        out.variant_of[i] = Some(j);
        let refs: Vec<evidence::ReadContextCounter> = ref_iter.next().into_iter().collect();
        out.variants.push(filters::SageVariant::new(final_candidates[j].clone(), refs, vec![tumor_opt[i].take().unwrap()]));
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
    let mut cfg = Config::from_args(&args)?;
    if cfg.tumor_ids.len() != 1 || cfg.reference_ids.len() > 1 {
        anyhow::bail!("only one tumor and at most one reference sample are supported so far");
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
        cfg.read_length = sample_read_length(&cfg, &ref_data, &mut main_worker.tumor_bam)?;
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
        msi.push(quality::MsiJitter::load(dir, sample)?);
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
    let mut vcf = vcfout::VcfWriter::create(&cfg.output_vcf, &contigs, &cfg.reference_ids, &cfg.tumor_ids)?;

    // every region task, in the order the sequence dictionary lists the chromosomes
    let dictionary_bam = main_worker.ref_bam.as_ref().unwrap_or(&main_worker.tumor_bam);
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
    Ok(())
}
