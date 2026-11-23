//! REDUX in Rust: a port of hmftools REDUX 2.1 (see ../redux).

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod args;
mod bamin;
mod bqr;
mod capture;
mod compare;
mod compat;
mod config;
#[allow(dead_code)]
mod fmt;
mod jitter;
mod consensus;
mod dup;
mod log;
mod out;
mod partition;
mod rec;
mod refgenome;
#[allow(dead_code)]
mod sbx;
mod sbx_dup;
mod regions;
mod sorted;
mod stage1;
mod ultima;
mod umi;
#[cfg(test)]
mod umi_test;
mod unmap;
#[cfg(test)]
mod tests;

use anyhow::Result;
use config::Config;
use regions::ChrBaseRegion;
use rust_htslib::bam::{self, Read};
use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(test)]
#[path = "tests/sbx_test.rs"]
mod sbx_test;

pub const VERSION: &str = "2.1";

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    // ConfigBuilder.checkHelpAndVersion: "-version" anywhere prints the version and exits
    if argv.iter().any(|a| a == "-version") {
        println!("Redux version {VERSION}");
        return;
    }
    if argv.first().map(String::as_str) == Some("compare_bams") && argv.len() >= 3 {
        if let Err(e) = compare::run(&argv[1..]) {
            error!("{e:#}");
            std::process::exit(1);
        }
        return;
    }
    if argv.first().map(String::as_str) == Some("diff_bams") && argv.len() >= 3 {
        if let Err(e) = compare::diff(&argv[1..]) {
            error!("{e:#}");
            std::process::exit(1);
        }
        return;
    }
    // -java_compat and -java_compat_except are not options of the Java tool: taken off before the options are parsed
    let argv = match compat::take_options(&argv) {
        Ok(a) => a,
        Err(e) => {
            error!("{e:#}");
            std::process::exit(1);
        }
    };
    if let Err(e) = run(&argv) {
        error!("{e:#}");
        std::process::exit(1);
    }
}

fn run(argv: &[String]) -> Result<()> {
    let start = std::time::Instant::now();
    info!("Redux version {VERSION} (Rust port)");
    let mut config = Config::from_args(argv)?;
    log::DEBUG.store(config.log_debug, Ordering::Relaxed);
    info!("output({})", config.output_dir);
    let genome = refgenome::RefGenome::open(&config.ref_genome_file)?;
    set_read_length(&mut config, &genome)?;

    let jitter = jitter::MsJitterAnalyser::new(&config)?;
    let mut bqr = bqr::BaseQualRecalibration::new(&config, &genome)?;
    bqr.set_partitions(&create_partitions(&config, &genome));

    if config.bqr_jitter_msi_only {
        return run_capture_only(&config, &genome, jitter.as_ref(), &bqr, start);
    }
    run_full(&config, &genome, jitter.as_ref(), &bqr, start)
}

/// -bqr_jitter_msi_only: the reads of an existing REDUX BAM go to the jitter analysis and BQR only
fn run_capture_only(
    config: &Config, genome: &refgenome::RefGenome, jitter: Option<&jitter::MsJitterAnalyser>, bqr: &bqr::BaseQualRecalibration,
    start: std::time::Instant,
) -> Result<()> {
    let partitions = create_partitions(config, genome);
    let worker_count = worker_limit(config.threads).min(partitions.len()).max(1);
    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| -> Result<()> {
        let mut handles = Vec::new();
        for _ in 0..worker_count {
            handles.push(scope.spawn(|| -> Result<()> {
                let mut reader = bamin::BamReader::open(&config.bam_files, &config.ref_genome_file, 1)?;
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= partitions.len() {
                        break;
                    }
                    let mut capture = capture::ReadCapture { jitter: jitter.map(|j| j.processor()), bqr: bqr.region_reader() };
                    for region in &partitions[i] {
                        capture.bqr.initialise(region);
                        let (rs, re) = (region.start, region.end);
                        reader.slice_region(&region.chromosome, rs, re, |rec| {
                            let pos = rec.pos() as i32 + 1;
                            if pos < rs || pos > re {
                                return;
                            }
                            capture.capture(rec, &region.chromosome);
                        });
                        capture.bqr.on_region_complete();
                    }
                }
                Ok(())
            }));
        }
        for h in handles {
            h.join().unwrap()?;
        }
        Ok(())
    })?;
    info!("all partition tasks complete");
    if let Some(j) = jitter {
        info!("analysing microsatellite jitter");
        j.write_analysis_output()?;
    }
    if bqr.enabled() {
        bqr.finalise()?;
    }
    info!("Redux BQR and jitter complete, mins({:.3})", start.elapsed().as_secs_f64() / 60.0);
    Ok(())
}

fn run_full(
    config: &Config, genome: &refgenome::RefGenome, jitter: Option<&jitter::MsJitterAnalyser>, bqr: &bqr::BaseQualRecalibration,
    start: std::time::Instant,
) -> Result<()> {
    if config.umi_enabled && !config.is_illumina() {
        anyhow::bail!("UMIs with SBX or Ultima reads have not been compared with Java and are not enabled");
    }
    let method = if config.skip_duplicate_marking {
        "disabled".to_string()
    } else {
        let mut m = format!("method({})", if config.form_consensus { "consensus" } else { "max base-qual" });
        if config.drop_duplicates {
            m.push_str(" drop-duplicates");
        }
        m
    };
    info!("duplicate marking: {method}");

    let header = bam::HeaderView::from_header(&out::build_combined_header(&config.bam_files, &config.ref_genome_file)?);
    let names: Vec<String> = (0..header.target_count()).map(|t| String::from_utf8_lossy(header.tid2name(t)).to_string()).collect();
    // ReduxConstants.EXPECTED_VALID_CONTIGS: the chromosomes, MT and EBV
    let expected_valid_tid: Vec<bool> = names
        .iter()
        .map(|n| {
            regions::HUMAN_CHROMOSOMES.iter().any(|c| regions::versioned_chromosome(config.v38, c) == *n)
                || *n == if config.v38 { "chrM" } else { "MT" }
                || *n == regions::versioned_chromosome(config.v38, "EBV")
        })
        .collect();

    let mut unmapper = match &config.unmap_regions_file {
        Some(f) => unmap::ReadUnmapper::load(f, &header)?,
        None => unmap::ReadUnmapper::disabled(&header),
    };
    if config.unmap_regions_file.is_some() && config.unmap_mt {
        let mt = if config.v38 { "chrM" } else { "MT" };
        if let Some(tid) = names.iter().position(|n| n == mt) {
            if unmapper.regions(tid).is_empty() {
                unmapper.add_region(tid, unmap::UnmappingRegion { start: 1, end: 20000, max_depth: unmap::UNMAP_MIN_HIGH_DEPTH });
            }
        }
    }
    let unmap_paired_reads = unmapper.enabled() && config.is_illumina();
    if unmap_paired_reads && config.unmap_alt_decoys {
        for (tid, n) in names.iter().enumerate() {
            if !regions::is_human_chromosome(n) && !regions::is_mitochondrial(n) {
                unmapper.add_region(tid, unmap::UnmappingRegion { start: 1, end: genome.chromosome_length(n), max_depth: unmap::UNMAP_MIN_HIGH_DEPTH });
            }
        }
    }

    let partitions = create_partitions(config, genome);
    let final_bam = config.output_bam.clone().unwrap_or_else(|| {
        let mut f = format!("{}{}.redux", config.output_dir, config.sample_id);
        if let Some(id) = &config.output_id {
            f.push('.');
            f.push_str(id);
        }
        f + ".bam"
    });
    let temp_prefix = format!("{final_bam}.tmp");
    let mut temp_files: Vec<String> = Vec::new();

    // first pass: unmap reads in the problematic regions
    let mut stage1_result = None;
    let bucket_paths: Vec<String> = (0..partitions.len()).map(|i| format!("{temp_prefix}.unmapping.{i}.bam")).collect();
    let sorted_bucket_paths: Vec<String> = (0..partitions.len()).map(|i| format!("{temp_prefix}.unmapping_sorted.{i}.bam")).collect();
    let mut have_buckets = false;
    if unmap_paired_reads && config.write_bam {
        let unmap_regions = stage1::create_unmap_regions(config, &unmapper);
        if !(config.specific.has_filters() && unmap_regions.is_empty()) {
            let lookup = stage1::PartitionLookup::new(&partitions, &names);
            let mut buckets = Vec::new();
            for p in &bucket_paths {
                let mut b = out::RawBamOut::create(p, 1)?;
                b.write_header(&header)?;
                buckets.push(std::sync::Mutex::new(b));
            }
            temp_files.extend(bucket_paths.iter().cloned());
            have_buckets = true;
            let r = stage1::run(
                config, &unmapper, &unmap_regions, &expected_valid_tid, &header, &lookup, &buckets, &temp_prefix, !config.skip_fully_unmapped_reads,
            )?;
            for b in buckets {
                b.into_inner().unwrap().close()?;
            }
            // sort each partition's relocated reads; a few at a time, since a bucket is held in memory to sort it
            let next_bucket = AtomicUsize::new(0);
            std::thread::scope(|scope| -> Result<()> {
                let mut handles = Vec::new();
                for _ in 0..worker_limit(config.threads).min(4).min(bucket_paths.len()).max(1) {
                    handles.push(scope.spawn(|| -> Result<()> {
                        loop {
                            let i = next_bucket.fetch_add(1, Ordering::Relaxed);
                            if i >= bucket_paths.len() {
                                return Ok(());
                            }
                            stage1::sort_bucket(&bucket_paths[i], &sorted_bucket_paths[i], &header)?;
                            let _ = std::fs::remove_file(&bucket_paths[i]);
                        }
                    }));
                }
                for h in handles {
                    h.join().unwrap()?;
                }
                Ok(())
            })?;
            temp_files.extend(sorted_bucket_paths.iter().cloned());
            info!("readsProcessed({}) unmapped stats: {}", r.processed_reads, r.stats.paired_string());
            temp_files.extend(r.temp_files.iter().cloned());
            stage1_result = Some(r);
        }
    }

    // unpaired sequencing types have no first pass: the reads that are already unmapped are copied here, and each
    // partition adds the primary reads it unmaps (ReduxApplication: processFullyUnmappedReads)
    let unpaired_full_unmapped = !unmap_paired_reads && !config.skip_fully_unmapped_reads && config.write_bam;
    let existing_unmapped_path = format!("{temp_prefix}.full_unmapped.existing.tmp");
    let partition_unmapped_paths: Vec<String> = (0..partitions.len()).map(|i| format!("{temp_prefix}.{i}.full_unmapped.part")).collect();
    let mut existing_unmapped: i64 = 0;
    if unpaired_full_unmapped {
        use rust_htslib::bam::Read as _;
        let mut out = out::RawBamOut::create(&existing_unmapped_path, config.compress_level)?;
        // RegionUnmapper.processFullyUnmappedReads: Ultima reads lose their attributes here too, without statistics
        let mut ultima = if config.is_ultima() { Some(ultima::UltimaRoutines::new(config.ultima_drop_qual_tags)) } else { None };
        for bam_file in &config.bam_files {
            let mut r = bam::IndexedReader::from_path(bam_file)?;
            r.set_reference(&config.ref_genome_file)?;
            r.fetch(bam::FetchDefinition::Unmapped)?;
            let mut rec = bam::Record::new();
            while let Some(x) = r.read(&mut rec) {
                x?;
                if let Some(u) = &mut ultima {
                    u.strip_attributes(&mut rec, None)?;
                }
                out.write(&rec)?;
                existing_unmapped += 1;
            }
        }
        out.close()?;
        temp_files.push(existing_unmapped_path.clone());
        temp_files.extend(partition_unmapped_paths.iter().cloned());
        if existing_unmapped > 0 {
            debug!("extracted {existing_unmapped} fully-unmapped reads");
        }
    }

    // second pass: the partitions
    let part_paths: Vec<String> = (0..partitions.len()).map(|i| format!("{temp_prefix}.{i}.sorted.part")).collect();
    let next = AtomicUsize::new(0);
    let worker_count = worker_limit(config.threads).min(partitions.len()).max(1);
    let outputs: Vec<std::sync::Mutex<Option<partition::PartitionOutput>>> = (0..partitions.len()).map(|_| std::sync::Mutex::new(None)).collect();
    let cache_sync = dup::CacheSync::new(partitions.len());
    // the final BAM is written under a temporary name while the partitions run: a joiner thread appends each
    // part as soon as it and all earlier parts are complete, so little is left to copy when the last one ends
    let final_tmp = format!("{temp_prefix}.joined");
    let join_state = JoinState::new(partitions.len());
    let mut joined: Option<(out::RawBamOut, Vec<(String, u64)>)> = None;
    std::thread::scope(|scope| -> Result<()> {
        let joiner = if config.write_bam {
            let mut final_out = out::RawBamOut::create(&final_tmp, config.compress_level)?;
            final_out.write_header(&header)?;
            final_out.flush()?;
            let (join_state, part_paths) = (&join_state, &part_paths);
            Some(scope.spawn(move || -> Result<Option<(out::RawBamOut, Vec<(String, u64)>)>> {
                // each part has its own index; with the offset at which the part starts they give the final index
                let mut part_indexes = Vec::new();
                let mut part_start = final_out.compressed_offset();
                for (i, p) in part_paths.iter().enumerate() {
                    if !join_state.wait_for(i) {
                        return Ok(None);
                    }
                    part_indexes.push((format!("{p}.bai"), part_start));
                    out::append_file(&mut final_out, p)?;
                    part_start += std::fs::metadata(p)?.len() - out::BGZF_EOF_LEN;
                    if !config.keep_interim_bams {
                        let _ = std::fs::remove_file(p);
                    }
                }
                Ok(Some((final_out, part_indexes)))
            }))
        } else {
            None
        };
        let mut handles = Vec::new();
        for _ in 0..worker_count {
            handles.push(scope.spawn(|| -> Result<()> {
                let mut reader = bamin::BamReader::open(&config.bam_files, &config.ref_genome_file, 1)?;
                let mut taken = None;
                let result = (|| -> Result<()> { loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= partitions.len() {
                        break Ok(());
                    }
                    taken = Some(i);
                    reader.set_extra(if have_buckets { Some(&sorted_bucket_paths[i]) } else { None })?;
                    let out = if config.write_bam {
                        let mut o = out::RawBamOut::create(&part_paths[i], config.compress_level)?;
                        o.enable_index();
                        Some(o)
                    } else {
                        None
                    };
                    let capture = capture::ReadCapture { jitter: jitter.map(|j| j.processor()), bqr: bqr.region_reader() };
                    let full_unmapped_out =
                        if unpaired_full_unmapped { Some(out::RawBamOut::create(&partition_unmapped_paths[i], config.compress_level)?) } else { None };
                    let mut pr = partition::PartitionReader::new(config, &names, &expected_valid_tid, &unmapper, out, full_unmapped_out, capture)?;
                    pr.set_cache_sync(cache_sync.clone(), i);
                    pr.process_partition(&partitions[i], &mut reader)?;
                    pr.publish_cache_sync();
                    *outputs[i].lock().unwrap() = Some(pr.finish()?);
                    join_state.part_done(i);
                } })();
                // a failed partition must not leave later ones waiting for its count
                if let (Err(_), Some(i)) = (&result, taken) {
                    cache_sync.publish_empty(i);
                }
                if result.is_err() {
                    join_state.fail();
                }
                result
            }));
        }
        let mut result = Ok(());
        for h in handles {
            let r = h.join().unwrap();
            if result.is_ok() {
                result = r;
            }
        }
        if result.is_err() {
            join_state.fail();
        }
        if let Some(j) = joiner {
            let r = j.join().unwrap();
            if result.is_ok() {
                joined = r?;
            }
        }
        result
    })
    .inspect_err(|_| {
        let _ = std::fs::remove_file(&final_tmp);
    })?;
    info!("all partition tasks complete");

    if let Some(j) = jitter {
        info!("analysing microsatellite jitter");
        j.write_analysis_output()?;
    }
    if bqr.enabled() {
        bqr.finalise()?;
    }

    let mut stats = dup::Statistics::default();
    let mut consensus_stats = consensus::ConsensusStatistics::default();
    let mut reunmap_stats = unmap::UnmapStats::default();
    let mut unsorted_write_count = 0;
    let mut no_coor = 0u64;
    for o in &outputs {
        let o = o.lock().unwrap().take().unwrap();
        stats.merge(&o.stats);
        consensus_stats.merge(&o.consensus_stats);
        reunmap_stats.merge(&o.unmap_stats);
        unsorted_write_count += o.unsorted_write_count;
        no_coor += o.full_unmapped_written;
    }

    if config.write_bam {
        let Some((mut final_out, part_indexes)) = joined else { anyhow::bail!("the partition files were not joined") };
        if let Some(r) = &stage1_result {
            for c in &r.full_unmapped_chunks {
                out::append_file_range(&mut final_out, &c.path, c.start, c.end)?;
                no_coor += c.reads;
            }
        }
        if unpaired_full_unmapped {
            out::append_file(&mut final_out, &existing_unmapped_path)?;
            for p in &partition_unmapped_paths {
                out::append_file(&mut final_out, p)?;
            }
        }
        final_out.close()?;
        std::fs::rename(&final_tmp, &final_bam)?;
        temp_files.extend(part_paths.iter().cloned());
        let final_index = format!("{final_bam}.bai");
        if part_indexes.iter().all(|(p, _)| std::path::Path::new(p).exists()) {
            out::merge_part_indexes(&part_indexes, header.target_count() as usize, no_coor + existing_unmapped as u64, &final_index)?;
        } else {
            out::index_bam(&final_bam, worker_limit(config.threads))?;
        }
        temp_files.extend(part_indexes.into_iter().map(|(p, _)| p));
        if !config.keep_interim_bams {
            for f in &temp_files {
                let _ = std::fs::remove_file(f);
            }
        }
    }
    if unsorted_write_count > 0 {
        warn!("unsorted BAM write count({unsorted_write_count}) via sorted BAM writers");
    }
    if !config.skip_duplicate_marking {
        info!(
            "stats: totalReads({}) duplicates({}) dupGroups({}) umiGroups({})",
            stats.total_reads, stats.duplicate_reads, stats.duplicate_groups, stats.umi.umi_groups
        );
        debug!("consensus stats: {}", consensus_stats.summary());
        if config.umi.enabled {
            let c = stats.umi.merge_counts;
            debug!(
                "umi merges: duplexCollapsed({}) polyGMerged({}) polyGUnmappedGroups({}) jitterMerged({})",
                c.duplex_collapsed, c.poly_g_merged, c.poly_g_unmapped_groups, c.jitter_merged
            );
        }
    }
    if unmapper.enabled() && !unmap_paired_reads {
        reunmap_stats.existing_unmapped += existing_unmapped;
        debug!("unmapped stats: {}", reunmap_stats.unpaired_string());
    }
    if let Some(r) = &stage1_result {
        debug!("re-unmapped stats: {}", reunmap_stats.paired_string());
        if reunmap_stats.read_count != r.stats.read_count || reunmap_stats.fully_unmapped_count != r.stats.fully_unmapped_count {
            warn!("re-unmapped stats differ: {}", reunmap_stats.paired_string());
        }
    }
    if !config.skip_duplicate_marking {
        stats.write_duplicate_stats(&format!("{}{}.redux.duplicate_freq.tsv", config.output_dir, config.sample_id), config.umi_enabled)?;
    }
    if config.umi.enabled {
        stats.umi.write_position_fragments_data(&config.form_filename("umi_coord_freq"))?;
        if config.umi.base_stats {
            stats.umi.write_umi_base_diff_stats(&config.form_filename("umi_edit_distance"))?;
            stats.umi.write_umi_base_frequency_stats(&config.form_filename("umi_nucleotide_freq"))?;
        }
    }
    if config.is_ultima() {
        stats.ultima.write_stats(&config.form_filename("ultima_stats"))?;
    }
    info!("Redux complete, mins({:.3})", start.elapsed().as_secs_f64() / 60.0);
    Ok(())
}

/// The number of worker threads for `-threads n`: n, or less when the environment variable REDUX_WORKERS is set.
/// `-threads` also fixes how the genome is split into partitions, and with it the output in Java mode; this
/// limits the threads without changing the split (for comparing with a Java run made with more threads than
/// the machine should be given).
pub fn worker_limit(threads: usize) -> usize {
    match std::env::var("REDUX_WORKERS").ok().and_then(|v| v.parse::<usize>().ok()) {
        Some(limit) if limit >= 1 => threads.min(limit),
        _ => threads,
    }
}

/// Tells the joiner thread which partition files are complete.
struct JoinState {
    /// per part: complete; then whether a partition failed
    state: std::sync::Mutex<(Vec<bool>, bool)>,
    changed: std::sync::Condvar,
}

impl JoinState {
    fn new(parts: usize) -> Self {
        JoinState { state: std::sync::Mutex::new((vec![false; parts], false)), changed: std::sync::Condvar::new() }
    }
    fn part_done(&self, i: usize) {
        self.state.lock().unwrap().0[i] = true;
        self.changed.notify_all();
    }
    fn fail(&self) {
        self.state.lock().unwrap().1 = true;
        self.changed.notify_all();
    }
    /// Waits until part `i` is complete; false when a partition failed instead.
    fn wait_for(&self, i: usize) -> bool {
        let mut s = self.state.lock().unwrap();
        while !s.0[i] && !s.1 {
            s = self.changed.wait(s).unwrap();
        }
        !s.1
    }
}

/// ReduxApplication.createPartitionThreads: the regions of each partition.
pub fn create_partitions(config: &Config, genome: &refgenome::RefGenome) -> Vec<Vec<ChrBaseRegion>> {
    let mut thread_count = config.threads;
    let mut partition_thread_count = if config.write_bam && config.parallel_concat { (thread_count - 1).max(1) } else { thread_count };
    let mut partition_count = config.partition_thread_ratio * partition_thread_count;
    if thread_count > 1 && !config.specific.regions.is_empty() {
        let bases: i32 = config.specific.regions.iter().map(|r| r.base_length()).sum();
        thread_count = thread_count.min((bases as f64 / 10_000.0).round() as usize).max(1);
        partition_thread_count = thread_count;
        partition_count = partition_thread_count;
    }
    let _ = partition_thread_count;
    regions::split_regions_into_partitions(&config.specific, partition_count, config.v38, genome, config.std_chr_only)
}

pub const DEFAULT_READ_LENGTH: i32 = 151;

/// ReduxApplication.setReadLength with BamSampler: the longest of up to 1000 primary, mapped, non-duplicate reads.
fn set_read_length(config: &mut Config, genome: &refgenome::RefGenome) -> Result<()> {
    if config.read_length > 0 {
        return Ok(());
    }
    let region = match config.specific.regions.first() {
        Some(r) => r.clone(),
        None => ChrBaseRegion::new(&genome.sequence_names()[0], 100000, 200000),
    };
    let mut reader = bam::IndexedReader::from_path(&config.bam_files[0])?;
    reader.set_reference(&config.ref_genome_file)?;
    let mut count = 0;
    let mut max_len = 0;
    let (mut paired, mut mate_cigar) = (false, false);
    let mut sbx_sampler = sbx::SbxBamSampler::default();
    if let Some(tid) = reader.header().tid(region.chromosome.as_bytes()) {
        reader.fetch((tid, (region.start - 1) as i64, region.end as i64))?;
        let mut rec = bam::Record::new();
        while let Some(Ok(())) = reader.read(&mut rec) {
            if rec.flags() & (0x4 | 0x100 | 0x400 | 0x800) != 0 {
                continue;
            }
            count += 1;
            max_len = max_len.max(rec.seq_len() as i32);
            if rec.is_paired() {
                paired = true;
                mate_cigar |= rec::has_aux(&rec, b"MC");
            }
            if config.is_sbx() {
                sbx_sampler.process_sample_record(&rec);
            }
            if count >= 1000 {
                break;
            }
        }
    }
    let mut read_length = DEFAULT_READ_LENGTH;
    if count > 0 && max_len > 0 {
        read_length = max_len;
        debug!("BAM sampled max read-length({read_length})");
    } else {
        debug!("BAM read-length sampling failed, using default read length({DEFAULT_READ_LENGTH})");
    }
    if paired && !mate_cigar {
        warn!("required mate CIGAR not set, exiting");
        std::process::exit(1);
    }
    config.read_length = read_length;
    Ok(())
}
