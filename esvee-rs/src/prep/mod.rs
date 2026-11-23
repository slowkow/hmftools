//! ESVEE prep step (Java package esvee.prep): finds candidate junctions in the input BAMs and writes the reads
//! around them, with the junctions, the fragment length distribution and discordant fragment counts.
//!
//! Java runs the partitions of a chromosome on several threads that share a cache of reads spanning
//! partitions, so its output depends on which partition finishes first. Here each partition is read and its
//! junctions are found on any thread, and the stage that uses the shared cache and writes output is applied to
//! the partitions in order. The result is that of Java on one thread, whatever the thread count.

pub mod config;
pub mod depth;
pub mod discordant;
pub mod filters;
pub mod fragsize;
pub mod hotspot;
pub mod jmap;
pub mod rec;
pub mod saga;
pub mod spanning;
pub mod tracker;
pub mod types;
pub mod writer;

#[cfg(test)]
mod tests;

use anyhow::{bail, Result};
use config::*;
use depth::{DepthTracker, DEPTH_WINDOW_SIZE};
use discordant::DiscordantStats;
use filters::*;
use hotspot::HotspotCache;
use jmap::JavaHashMap;
use rec::*;
use rust_htslib::bam::{self, Read};
use saga::SagaMatcherFactory;
use spanning::SpanningReadCache;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use tracker::{JunctionTracker, TrackerConfig};
use types::*;
use writer::ResultsWriter;

/// PartitionStats
#[derive(Default, Clone)]
pub struct PartitionStats {
    pub total_reads: i64,
    pub junction_count: i32,
    pub junction_fragment_count: i32,
    pub initial_supporting_fragment_count: i32,
    pub supporting_fragment_count: i32,
    pub local_complete_groups: i32,
    pub local_incomplete_groups: i32,
    pub spanning_groups: i32,
    pub read_filter_counts: [i32; 8],
}

impl PartitionStats {
    fn add(&mut self, o: &PartitionStats) {
        self.total_reads += o.total_reads;
        self.junction_count += o.junction_count;
        self.junction_fragment_count += o.junction_fragment_count;
        self.initial_supporting_fragment_count += o.initial_supporting_fragment_count;
        self.supporting_fragment_count += o.supporting_fragment_count;
        self.local_complete_groups += o.local_complete_groups;
        self.local_incomplete_groups += o.local_incomplete_groups;
        self.spanning_groups += o.spanning_groups;
        for i in 0..8 {
            self.read_filter_counts[i] += o.read_filter_counts[i];
        }
    }
}

impl std::fmt::Display for PartitionStats {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "reads({}) junc({}) juncFrags({}) supportFrags(init={} final={}) groups(comp={} incomp={} span={})",
            self.total_reads,
            self.junction_count,
            self.junction_fragment_count,
            self.initial_supporting_fragment_count,
            self.supporting_fragment_count,
            self.local_complete_groups,
            self.local_incomplete_groups,
            self.spanning_groups
        )
    }
}

struct Partition {
    region: Region,
    last_of_chromosome: bool,
    first_of_chromosome: bool,
    chr_partitions: usize,
}

/// What the threads share without change.
struct Context {
    config: PrepConfig,
    contigs: Arc<Contigs>,
    /// per sample: the contig index of each reference of its BAM
    tid_maps: Vec<Vec<i32>>,
    hotspots: HotspotCache,
    saga: Option<SagaMatcherFactory>,
    tracker_cfg: TrackerConfig,
    poly_g_region: Region,
}

/// The first stage of a partition, done on any thread.
struct PartitionResult {
    tracker: JunctionTracker,
    stats: PartitionStats,
}

/// The second stage is applied in partition order under this lock.
struct Sequencer {
    next: usize,
    pending: BTreeMap<usize, PartitionResult>,
    /// the reads held by the partitions in `pending`
    pending_reads: usize,
    cache: SpanningReadCache,
    writer: ResultsWriter,
    chr_stats: PartitionStats,
    total_stats: PartitionStats,
    discordants: DiscordantStats,
    failed: bool,
}

/// BamSlicerFilter as prep sets it: unmapped and supplementary reads are kept, duplicates are not, and a
/// secondary read only when hard clipped
#[inline]
fn passes_slicer_filter(rec: &bam::Record) -> bool {
    let flags = rec.flags();
    if flags & FLAG_SECONDARY != 0 && !rec.raw_cigar().iter().any(|c| op(*c) == OP_H) {
        return false;
    }
    flags & FLAG_DUPLICATE == 0
}

/// PartitionSlicer.run up to the junction assignment
fn slice_partition(ctx: &Context, readers: &mut [bam::IndexedReader], region: &Region) -> Result<PartitionResult> {
    let config = &ctx.config;
    let filter = &ctx.tracker_cfg.filter;
    let chr_name = ctx.contigs.name(region.chr);
    crate::debug!("processing region({})", region.display(&ctx.contigs));
    let mut depth = DepthTracker::new(region.start, region.end, DEPTH_WINDOW_SIZE);
    let saga_matcher = ctx.saga.as_ref().map(|f| {
        let margin = 2 * filter.observed_frag_length_max;
        f.create_location_matcher(chr_name, region.start - margin, region.end + margin)
    });
    let mut tracker = JunctionTracker::new(region.clone(), ctx.tracker_cfg.clone(), ctx.hotspots.find_region_hotspots(region), saga_matcher);
    let filter_region = if region.overlaps_region(&ctx.poly_g_region) { Some(&ctx.poly_g_region) } else { None };
    let mut stats = PartitionStats::default();
    let write_reads = config.write_reads();
    let log_read_ids = !config.log_read_ids.is_empty();

    for (sample, reader) in readers.iter_mut().enumerate() {
        let Some(tid) = reader.header().tid(chr_name.as_bytes()) else {
            crate::warn!("cannot find sequence index for chromosome {chr_name} in bam header");
            continue;
        };
        reader.fetch((tid, (region.start - 1) as i64, region.end as i64))?;
        let tid_map = &ctx.tid_maps[sample];
        let mut rec = bam::Record::new();
        while let Some(r) = reader.read(&mut rec) {
            r?;
            if !passes_slicer_filter(&rec) {
                continue;
            }
            let read_start = start(&rec);
            if !region.contains_position(read_start) {
                continue;
            }
            stats.total_reads += 1;
            let read_end = end(&rec);
            if let Some(fr) = filter_region {
                if positions_overlap(read_start, read_end, fr.start, fr.end) {
                    continue;
                }
            }
            if log_read_ids && config.log_read_ids.iter().any(|id| id.as_bytes() == rec.qname()) {
                crate::debug!("specific readId({}) unmapped({})", String::from_utf8_lossy(rec.qname()), has(&rec, FLAG_UNMAPPED));
            }
            let line_tail = has_line_tail(&rec);
            let unmapped = has(&rec, FLAG_UNMAPPED);
            if ignore_read(&mut rec, line_tail, filter) {
                if !unmapped {
                    depth.process_range(read_start, read_end);
                }
                continue;
            }
            let mut read = PrepRead::new(std::mem::replace(&mut rec, bam::Record::new()), sample as u8, tid_map, &ctx.contigs);
            read.has_line_tail = line_tail;
            check_filters(&mut read, filter);
            if config.perf_debug {
                for (i, count) in stats.read_filter_counts.iter_mut().enumerate() {
                    if read.filters & (1 << i) != 0 {
                        *count += 1;
                    }
                }
            }
            if read.unfiltered() || read.filters == F_MIN_MAP_QUAL {
                // a read filtered only for its mapping quality still makes a junction
                read.set_read_type(ReadType::Junction);
                tracker.process_read(read);
            } else {
                let is_support_candidate = is_candidate_supporting_read(&read, filter);
                if is_support_candidate || write_reads {
                    if is_support_candidate {
                        read.set_read_type(ReadType::CandidateSupport);
                    }
                    tracker.process_read(read);
                }
            }
            if !unmapped {
                depth.process_range(read_start, read_end);
            }
        }
    }
    tracker.assign_fragments(&depth);
    if tracker.duplicate_groups > 100 {
        crate::debug!("region({}) marked {} supplementary duplicates", region.display(&ctx.contigs), tracker.duplicate_groups);
    }
    Ok(PartitionResult { tracker, stats })
}

/// The rest of PartitionSlicer.run: the expected reads from the earlier partitions, then writeData.
fn complete_partition(ctx: &Context, seq: &mut Sequencer, result: PartitionResult) -> Result<()> {
    let PartitionResult { mut tracker, mut stats } = result;
    let config = &ctx.config;
    let region = tracker.region.clone();
    let expected = seq.cache.get_expected_read_ids(&region);
    tracker.set_expected_reads(expected);
    tracker.mark_expected_groups();

    let junction_groups = tracker.form_unique_assigned_groups();
    let junction_text = seq.writer.junction_data_text(ctx.contigs.name(region.chr), &tracker);
    let outside_junction_rows = seq.writer.outside_junction_rows(&tracker);
    let mut write_junctions = true;
    if config.write_reads() {
        // groups with reads in other partitions go to the cache, which completes them across partitions
        let mut spanning: JavaHashMap<Box<[u8]>, (GroupId, bool, bool)> = JavaHashMap::new();
        let mut assign = |tracker: &mut JunctionTracker, gid: GroupId, from_junctions: bool, stats: &mut PartitionStats| {
            let g = &mut tracker.groups[gid as usize];
            group_set_partition_count(g, &tracker.reads, &region, config.partition_size, &ctx.contigs);
            if g.spans_partitions() && g.status != ReadGroupStatus::Expected {
                stats.spanning_groups += 1;
                group_set_state(g, &tracker.reads);
                let will_write = from_junctions && g.has_remote_junction_reads && g.status != ReadGroupStatus::Duplicate;
                spanning.insert(g.id.clone(), (gid, will_write, from_junctions));
            } else {
                group_set_state(g, &tracker.reads);
                if g.is_complete() {
                    stats.local_complete_groups += 1;
                } else {
                    stats.local_incomplete_groups += 1;
                }
            }
        };
        for &gid in &junction_groups {
            assign(&mut tracker, gid, true, &mut stats);
        }
        let remote_candidate_groups = tracker.get_remote_candidate_read_groups();
        for &gid in &remote_candidate_groups {
            assign(&mut tracker, gid, false, &mut stats);
        }
        let total_group_count = junction_groups.len() + remote_candidate_groups.len();
        let spanning_order: Vec<(GroupId, bool, bool)> = spanning.values().copied().collect();
        seq.cache.process_spanning_read_groups(&region, &spanning_order, &mut tracker)?;
        if total_group_count == 0 {
            write_junctions = false;
        } else {
            if crate::common::log::debug_enabled() {
                let expected_group_count = junction_groups.iter().filter(|&&g| tracker.groups[g as usize].status == ReadGroupStatus::Expected).count();
                crate::debug!(
                    "region({}) readGroups({} spanning={} candidate={}) expected({})",
                    region.display(&ctx.contigs),
                    total_group_count,
                    spanning_order.len(),
                    remote_candidate_groups.len(),
                    expected_group_count
                );
            }
            seq.writer.write_read_groups(&mut tracker, &junction_groups)?;
        }
    }
    if write_junctions {
        seq.writer.write_junction_text(&junction_text)?;
        seq.writer.defer_junctions(outside_junction_rows);
        stats.junction_count += (tracker.junctions.len() + tracker.outside_junctions.len()) as i32;
        let junction_fragments = junction_groups.iter().filter(|&&g| group_has_read_type(&tracker.groups[g as usize], &tracker.reads, ReadType::Junction)).count() as i32;
        stats.junction_fragment_count += junction_fragments;
        stats.supporting_fragment_count += junction_groups.len() as i32 - junction_fragments;
        stats.initial_supporting_fragment_count += tracker.initial_supporting_frags;
    }
    if stats.total_reads > 0 {
        crate::debug!("region({}) complete, stats({})", region.display(&ctx.contigs), stats);
    }
    seq.chr_stats.add(&stats);
    seq.discordants.add(&tracker.discordant_stats);
    Ok(())
}

/// ChromosomeTask.partition for every chromosome in turn
fn build_partitions(config: &PrepConfig, contigs: &Contigs) -> Vec<Partition> {
    let mut partitions: Vec<Partition> = Vec::new();
    for (i, name) in HUMAN_CHROMOSOMES.iter().enumerate() {
        let chromosome = versioned_chromosome(config.v38, name);
        if config.specific.exclude_chromosome(&chromosome) {
            continue;
        }
        let chr = contigs.id(&chromosome);
        let first = partitions.len();
        let mut add = |min_position: i32, max_position: i32| {
            let mut k = 0i64;
            loop {
                let start = min_position as i64 + k * config.partition_size as i64;
                let end = (start + config.partition_size as i64 - 1).min(max_position as i64);
                partitions.push(Partition { region: Region::new(chr, start as i32, end as i32), last_of_chromosome: false, first_of_chromosome: false, chr_partitions: 0 });
                if end >= max_position as i64 {
                    break;
                }
                k += 1;
            }
        };
        if !config.specific.regions.is_empty() {
            for r in config.specific.regions.iter().filter(|r| r.0 == chromosome) {
                add(r.1, r.2);
            }
        } else {
            add(1, chromosome_length(config.v38, i));
        }
        let count = partitions.len() - first;
        for p in &mut partitions[first..] {
            p.chr_partitions = count;
        }
        if count > 0 {
            partitions[first].first_of_chromosome = true;
            partitions.last_mut().unwrap().last_of_chromosome = true;
        } else {
            // a chromosome without partitions is still announced, as in Java
            partitions.push(Partition { region: Region::new(chr, 0, -1), last_of_chromosome: true, first_of_chromosome: true, chr_partitions: 0 });
        }
    }
    partitions
}

/// Opens an input for region queries. For CRAM, htslib is told not to make up MD and NM attributes that the
/// file does not store: htsjdk does not, and the records are written with the attributes they were read with.
pub fn open_indexed_reader(path: &str, ref_genome_file: &str) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(path)?;
    reader.set_reference(ref_genome_file)?;
    unsafe {
        rust_htslib::htslib::hts_set_opt(reader.htsfile(), rust_htslib::htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 0 as std::os::raw::c_int);
    }
    Ok(reader)
}

fn open_readers(config: &PrepConfig) -> Result<Vec<bam::IndexedReader>> {
    config.bam_files.iter().map(|f| open_indexed_reader(f, &config.ref_genome_file)).collect()
}

pub fn run(argv: &[String]) -> Result<()> {
    if crate::version_log_due() {
        crate::info!("Esvee version {}", crate::VERSION);
    }
    let mut config = PrepConfig::parse(argv)?;
    crate::info!("output({}) {}", config.output_dir, config.output_id.as_ref().map(|id| format!("outputId({id})")).unwrap_or_default());
    let start_time = std::time::Instant::now();
    std::fs::create_dir_all(&config.output_dir).ok();

    // the contigs of the inputs under one numbering, then the human chromosomes in both naming styles
    let mut contigs = Contigs::new();
    let mut tid_maps: Vec<Vec<i32>> = Vec::new();
    let mut headers: Vec<String> = Vec::new();
    for f in &config.bam_files {
        let mut reader = bam::Reader::from_path(f)?;
        reader.set_reference(&config.ref_genome_file)?;
        let header = reader.header();
        tid_maps.push(header.target_names().iter().map(|n| contigs.add(&String::from_utf8_lossy(n))).collect());
        headers.push(String::from_utf8_lossy(header.as_bytes()).into_owned());
    }
    for name in HUMAN_CHROMOSOMES {
        contigs.add(name);
        contigs.add(&format!("chr{name}"));
    }
    let contigs = Arc::new(contigs);

    let hotspots = match &config.known_hotspot_file {
        Some(f) => {
            let cache = HotspotCache::load(f, &contigs)?;
            crate::info!("loaded {} known hotspot records from file({})", cache.count, f);
            cache
        }
        None => HotspotCache::empty(),
    };
    let saga = match &config.saga_fasta_file {
        Some(f) => Some(SagaMatcherFactory::load(f)?),
        None => None,
    };

    crate::info!("running Esvee Prep for {}", config.sample_ids.join(","));

    let distribution = fragsize::FragmentSizeDistribution::run(&config)?;
    let bounds = distribution.bounds();
    crate::info!("fragment length distribution bounds(min={} max={} median={})", bounds.lower_bound, bounds.upper_bound, bounds.median);
    if bounds.is_valid() {
        let max = if config.max_fragment_length_override > 0 { config.max_fragment_length_override } else { bounds.upper_bound };
        config.filter.set_fragment_lengths(bounds.lower_bound, max);
    }

    let partitions = build_partitions(&config, &contigs);
    let poly_g_region =
        if config.v38 { Region::new(contigs.id("chr2"), 32916190, 32916630) } else { Region::new(contigs.id("2"), 33141260, 33141700) };
    let cache_header = bam::HeaderView::from_bytes(writer::header_with_sort_order(&writer::combined_header_text(&headers), "unsorted").as_bytes());
    let cache = SpanningReadCache::new(&config, contigs.clone(), Some(cache_header));
    let results_writer = ResultsWriter::new(&config, contigs.clone(), &headers)?;
    let tracker_cfg = TrackerConfig { filter: config.filter.clone(), trim_read_id: config.trim_read_id, contigs: contigs.clone() };
    let threads = config.threads;
    let ctx = Context { config, contigs: contigs.clone(), tid_maps, hotspots, saga, tracker_cfg, poly_g_region };

    let sequencer = Mutex::new(Sequencer {
        next: 0,
        pending: BTreeMap::new(),
        pending_reads: 0,
        cache,
        writer: results_writer,
        chr_stats: PartitionStats::default(),
        total_stats: PartitionStats::default(),
        discordants: DiscordantStats::default(),
        failed: false,
    });
    let turn = Condvar::new();
    let next_partition = AtomicUsize::new(0);
    // How far the first stage may run ahead of the second: a couple of partitions per thread, or further
    // while the finished partitions that wait their turn hold few reads (sparse data, where most partitions
    // are empty and one busy partition would otherwise hold up every thread).
    let window = threads * 2;
    const SPARSE_WINDOW: usize = 1024;
    const SPARSE_PENDING_READS: usize = 1_000_000;

    std::thread::scope(|scope| {
        for _ in 0..threads.min(partitions.len()).max(1) {
            scope.spawn(|| {
                let mut readers = match open_readers(&ctx.config) {
                    Ok(r) => r,
                    Err(e) => {
                        crate::error!("thread execution error: {e:#}");
                        std::process::exit(1);
                    }
                };
                loop {
                    let i = next_partition.fetch_add(1, Ordering::SeqCst);
                    if i >= partitions.len() {
                        break;
                    }
                    {
                        let mut seq = sequencer.lock().unwrap();
                        while i >= seq.next + window && !(i < seq.next + SPARSE_WINDOW && seq.pending_reads < SPARSE_PENDING_READS) && !seq.failed {
                            seq = turn.wait(seq).unwrap();
                        }
                    }
                    let part = &partitions[i];
                    let result = if part.chr_partitions == 0 {
                        Ok(PartitionResult { tracker: JunctionTracker::new(part.region.clone(), ctx.tracker_cfg.clone(), Vec::new(), None), stats: PartitionStats::default() })
                    } else {
                        slice_partition(&ctx, &mut readers, &part.region)
                    };
                    let result = match result {
                        Ok(r) => r,
                        Err(e) => {
                            crate::error!("thread execution error: {e:#}");
                            std::process::exit(1);
                        }
                    };
                    let mut seq = sequencer.lock().unwrap();
                    seq.pending_reads += result.tracker.reads.len();
                    seq.pending.insert(i, result);
                    loop {
                        let next = seq.next;
                        let Some(r) = seq.pending.remove(&next) else { break };
                        seq.pending_reads -= r.tracker.reads.len();
                        let p = &partitions[next];
                        let chr_name = ctx.contigs.name(p.region.chr);
                        if p.first_of_chromosome {
                            crate::info!("processing chromosome({chr_name})");
                            seq.chr_stats = PartitionStats::default();
                        }
                        if p.chr_partitions > 0 {
                            if let Err(e) = complete_partition(&ctx, &mut seq, r) {
                                crate::error!("thread execution error: {e:#}");
                                std::process::exit(1);
                            }
                        }
                        if p.last_of_chromosome {
                            crate::info!("chromosome({}) {} regions complete, stats: {}", chr_name, p.chr_partitions, seq.chr_stats);
                            if seq.chr_stats.read_filter_counts.iter().any(|c| *c > 0) {
                                let counts: Vec<String> = FILTER_NAMES.iter().zip(seq.chr_stats.read_filter_counts).map(|(n, c)| format!("{n}={c}")).collect();
                                crate::debug!("chromosome({}) filters({})", chr_name, counts.join(", "));
                            }
                            if ctx.config.perf_debug {
                                seq.cache.log_stats();
                            }
                            let chr_stats = seq.chr_stats.clone();
                            seq.total_stats.add(&chr_stats);
                        }
                        seq.next += 1;
                    }
                    drop(seq);
                    turn.notify_all();
                }
            });
        }
    });

    let Sequencer { mut cache, mut writer, total_stats, discordants, .. } = sequencer.into_inner().unwrap();
    let config = &ctx.config;
    if config.sample_ids.is_empty() {
        bail!("no samples");
    }
    // Java: writes the discordant statistics here, so PrepReads lacks the candidate reads claimed below
    let java_p2 = crate::common::compat::java("P2");
    let write_disc_stats = |writer: &ResultsWriter| -> Result<()> {
        std::fs::write(config.form_filename(WriteType::DiscordantStats, config.sample_id()), discordants.file_text(total_stats.total_reads, writer.written_count() as i64))?;
        Ok(())
    };
    if java_p2 {
        write_disc_stats(&writer)?;
    }

    if config.use_cache_bam {
        crate::debug!("records written before the cached candidate reads, by sample: {:?}", writer.bam_writer.sink_counts);
        cache.reset();
        cache.candidate_bam_writer.assign_candidate_reads(&mut writer, &ctx.tid_maps[0])?;
    }
    if !java_p2 {
        write_disc_stats(&writer)?;
    }
    writer.close()?;

    if config.perf_debug && (total_stats.total_reads > 10000 || start_time.elapsed().as_millis() > 10000) {
        crate::info!("final stats: {}", total_stats);
    }
    crate::debug!("{} hash buckets were held as trees, as Java does", jmap::tree_bin_events());
    crate::info!("Esvee prep complete, mins({:.3})", start_time.elapsed().as_millis() as f64 / 60000.0);
    Ok(())
}
