//! SpanningReadCache and CandidateBamWriter: what one partition leaves for the others. A partition with a
//! junction fragment records its read id against the partitions of its remote reads, where the read is then
//! "expected"; a partition with only candidate reads of a fragment parks them, in memory or in per-chromosome
//! cache BAMs, until a junction elsewhere claims them.

use super::config::*;
use super::jmap::JavaHashMap;
use super::rec::*;
use super::tracker::{JunctionTracker, ReadIdTrimmer};
use super::types::*;
use super::writer::{RawBamOut, ResultsWriter};
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

pub const SAMPLE_ID_TAG: &[u8; 2] = b"SI";

struct CachedReadGroup {
    reads: Vec<PrepRead>,
    partitions: Vec<PartKey>,
    /// the number of partition maps that hold the group
    refs: usize,
}

pub struct CandidateBamWriter {
    use_cache_bam: bool,
    no_clean_up: bool,
    file_prefix: String,
    sample_ids: Vec<String>,
    header: Option<bam::HeaderView>,
    writers: HashMap<i32, RawBamOut>,
    files: HashMap<i32, String>,
    /// by chromosome name, in the iteration order of the Java map
    chr_junction_read_ids: JavaHashMap<String, HashSet<Box<[u8]>>>,
    trimmer: ReadIdTrimmer,
    contigs: Arc<Contigs>,
}

impl CandidateBamWriter {
    fn add_junction_read_id(&mut self, remote_partitions: &[PartKey], read_id: &[u8]) {
        for p in remote_partitions {
            // chrFromChrPartition: the partition name up to its first '_'
            let name = self.contigs.name(p.chr);
            let chromosome = name.split('_').next().unwrap_or(name).to_string();
            if !self.chr_junction_read_ids.contains_key(&chromosome) {
                self.chr_junction_read_ids.insert(chromosome.clone(), HashSet::new());
            }
            self.chr_junction_read_ids.get_mut(&chromosome).unwrap().insert(read_id.into());
        }
    }

    fn write_candidate_read(&mut self, read: &PrepRead) -> Result<()> {
        if !self.use_cache_bam {
            return Ok(());
        }
        if !self.writers.contains_key(&read.chr) {
            let file = format!("{}_{}.bam", self.file_prefix, self.contigs.name(read.chr));
            let mut out = RawBamOut::create(&file, 1)?;
            if let Some(h) = &self.header {
                out.write_header(h)?;
            }
            self.writers.insert(read.chr, out);
            self.files.insert(read.chr, file);
        }
        let mut rec = read.rec.clone();
        set_attribute(&mut rec, SAMPLE_ID_TAG, Some(&encode_z(self.sample_ids[read.sample as usize].as_bytes())));
        self.writers.get_mut(&read.chr).unwrap().write(&rec)
    }

    /// assignCandidateReads: every cached read of a chromosome whose id a junction fragment registered
    /// against that chromosome is written, in the order of the cache file.
    pub fn assign_candidate_reads(&mut self, writer: &mut ResultsWriter, tid_map: &[i32]) -> Result<()> {
        if self.writers.is_empty() {
            return Ok(());
        }
        crate::info!("assigning candidate mate read to junctions");
        for (_, w) in self.writers.drain() {
            w.close()?;
        }
        let chromosomes: Vec<String> = self.chr_junction_read_ids.keys().cloned().collect();
        for chromosome in chromosomes {
            let chr = self.contigs.id(&chromosome);
            let Some(file) = self.files.get(&chr).cloned() else {
                continue;
            };
            crate::debug!("chr({}) assigning candidates from {} junction fragments", chromosome, self.chr_junction_read_ids.get(&chromosome).map_or(0, |s| s.len()));
            let mut reader = bam::Reader::from_path(&file)?;
            let names: Vec<String> = reader.header().target_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
            let mut matched = 0u64;
            loop {
                let mut rec = bam::Record::new();
                match reader.read(&mut rec) {
                    Some(r) => r?,
                    None => break,
                }
                if rec.flags() & (FLAG_DUPLICATE | FLAG_SECONDARY) != 0 {
                    continue;
                }
                let tid = rec.tid();
                if tid < 0 {
                    continue;
                }
                let read_id = self.trimmer.trim(rec.qname());
                let is_junction_read = self.chr_junction_read_ids.get(&names[tid as usize]).is_some_and(|ids| ids.contains(read_id));
                if !is_junction_read {
                    continue;
                }
                matched += 1;
                let sample = aux_str(&rec, SAMPLE_ID_TAG).and_then(|s| self.sample_ids.iter().position(|x| x.as_bytes() == s)).unwrap_or(0) as u8;
                set_attribute(&mut rec, SAMPLE_ID_TAG, None);
                let mut read = PrepRead::new(rec, sample, tid_map, &self.contigs);
                read.set_read_type(ReadType::CandidateSupport);
                writer.write_candidate_read(read)?;
            }
            crate::debug!("chr({chromosome}) matched and wrote {matched} candidate reads");
        }
        crate::info!("candidate reads assignment complete");
        if !self.no_clean_up {
            crate::debug!("deleting {} candidate BAMs", self.files.len());
            for file in self.files.values() {
                if let Err(e) = std::fs::remove_file(file) {
                    crate::error!("failed to delete candidate BAM: {e}");
                }
            }
        }
        Ok(())
    }
}

pub struct SpanningReadCache {
    partition_size: i32,
    use_cache_bam: bool,
    specific_chromosomes: Vec<String>,
    contigs: Arc<Contigs>,
    /// by partition, then read id: the cached group (an index into `cached`)
    candidate_partition_groups: HashMap<PartKey, HashMap<Box<[u8]>, usize>>,
    cached: HashMap<usize, CachedReadGroup>,
    next_handle: usize,
    junction_partition_read_ids: HashMap<PartKey, HashSet<Box<[u8]>>>,
    processed_partitions: HashSet<PartKey>,
    last_snapshot_count: i64,
    matched_candidates: i64,
    purged_candidates: i64,
    pub candidate_bam_writer: CandidateBamWriter,
}

impl SpanningReadCache {
    pub fn new(config: &PrepConfig, contigs: Arc<Contigs>, cache_header: Option<bam::HeaderView>) -> Self {
        SpanningReadCache {
            partition_size: config.partition_size,
            use_cache_bam: config.use_cache_bam,
            specific_chromosomes: config.specific.chromosomes.clone(),
            contigs: contigs.clone(),
            candidate_partition_groups: HashMap::new(),
            cached: HashMap::new(),
            next_handle: 0,
            junction_partition_read_ids: HashMap::new(),
            processed_partitions: HashSet::new(),
            last_snapshot_count: 0,
            matched_candidates: 0,
            purged_candidates: 0,
            candidate_bam_writer: CandidateBamWriter {
                use_cache_bam: config.use_cache_bam,
                no_clean_up: config.no_clean_up,
                file_prefix: if config.sample_ids.is_empty() { String::new() } else { config.form_filename(WriteType::CacheBam, config.sample_id()) },
                sample_ids: config.sample_ids.clone(),
                header: cache_header,
                writers: HashMap::new(),
                files: HashMap::new(),
                chr_junction_read_ids: JavaHashMap::new(),
                trimmer: ReadIdTrimmer::new(config.trim_read_id),
                contigs,
            },
        }
    }

    fn chr_partition(&self, chr: i32, position: i32) -> PartKey {
        PartKey::form(chr, position, self.partition_size)
    }

    /// The read ids that earlier partitions expect in this one.
    pub fn get_expected_read_ids(&self, region: &Region) -> HashSet<Box<[u8]>> {
        self.junction_partition_read_ids.get(&self.chr_partition(region.chr, region.start)).cloned().unwrap_or_default()
    }

    /// ignoreChromosome: with specific chromosomes set, a name whose part before the first '_' is not one of them
    fn ignore_chromosome(&self, name: &str) -> bool {
        if self.specific_chromosomes.is_empty() {
            return false;
        }
        let chromosome = name.split('_').next().unwrap_or(name);
        !self.specific_chromosomes.iter().any(|c| c == chromosome)
    }

    /// processSpanningReadGroups. `spanning` holds the groups in the iteration order of the Java map, each
    /// with whether the partition is about to write it.
    pub fn process_spanning_read_groups(&mut self, region: &Region, spanning: &[(GroupId, bool)], tracker: &mut JunctionTracker) -> Result<()> {
        let source = self.chr_partition(region.chr, region.start);
        self.processed_partitions.insert(source);
        let contigs = self.contigs.clone();
        for &(gid, will_write) in spanning {
            let remote_partitions = tracker.groups[gid as usize].remote_partitions_ordered(&contigs);
            let unprocessed: Vec<PartKey> = remote_partitions
                .iter()
                .copied()
                .filter(|p| !self.processed_partitions.contains(p))
                .filter(|p| self.specific_chromosomes.is_empty() || !self.ignore_chromosome(contigs.name(p.chr)))
                .collect();
            if group_conditional_on_remote_reads(&tracker.groups[gid as usize], &tracker.reads) {
                let mut k = 0;
                while k < tracker.groups[gid as usize].reads.len() {
                    let rid = tracker.groups[gid as usize].reads[k];
                    k += 1;
                    self.process_candidate_read(&unprocessed, &remote_partitions, gid, rid, will_write, tracker)?;
                }
            } else {
                self.process_junction_read(&unprocessed, &remote_partitions, gid, source, tracker);
            }
        }
        self.purge_partition(source);
        self.log_cache_count(false);
        Ok(())
    }

    fn release(&mut self, handle: usize) {
        if let Some(c) = self.cached.get_mut(&handle) {
            c.refs -= 1;
            if c.refs == 0 {
                self.cached.remove(&handle);
            }
        }
    }

    fn process_junction_read(&mut self, unprocessed: &[PartKey], remote_partitions: &[PartKey], gid: GroupId, source: PartKey, tracker: &mut JunctionTracker) {
        let first_chr = tracker.reads[tracker.groups[gid as usize].reads[0] as usize].chr;
        if self.ignore_chromosome(self.contigs.name(first_chr)) {
            return;
        }
        let id = tracker.groups[gid as usize].id.clone();
        if self.use_cache_bam {
            self.candidate_bam_writer.add_junction_read_id(remote_partitions, &id);
        } else if let Some(handle) = self.candidate_partition_groups.get_mut(&source).and_then(|m| m.remove(&id)) {
            // candidate reads found earlier for this fragment join the group
            let (reads, partitions) = {
                let c = self.cached.get_mut(&handle).unwrap();
                (std::mem::take(&mut c.reads), c.partitions.clone())
            };
            self.matched_candidates += reads.len() as i64;
            for read in reads {
                let rid = tracker.add_read(read);
                tracker.groups[gid as usize].reads.push(rid);
            }
            self.release(handle);
            for p in partitions {
                if p != source {
                    if let Some(h) = self.candidate_partition_groups.get_mut(&p).and_then(|m| m.remove(&id)) {
                        self.release(h);
                    }
                }
            }
        }
        // the read id is kept for the remote partitions still to come, where the fragment's reads are expected
        for p in unprocessed {
            self.junction_partition_read_ids.entry(*p).or_default().insert(id.clone());
        }
    }

    fn process_candidate_read(
        &mut self, unprocessed: &[PartKey], remote_partitions: &[PartKey], gid: GroupId, rid: ReadId, will_write: bool, tracker: &mut JunctionTracker,
    ) -> Result<()> {
        if unprocessed.is_empty() {
            return Ok(());
        }
        if self.ignore_chromosome(self.contigs.name(tracker.reads[rid as usize].chr)) {
            return Ok(());
        }
        if self.use_cache_bam {
            return self.candidate_bam_writer.write_candidate_read(&tracker.reads[rid as usize]);
        }
        let id = tracker.groups[gid as usize].id.clone();
        // an existing cached group for the fragment, under any of its remote partitions
        let mut cached: Option<usize> = None;
        let mut matched: Vec<PartKey> = Vec::new();
        for p in remote_partitions {
            let Some(map) = self.candidate_partition_groups.get(p) else { continue };
            if cached.is_none() {
                cached = map.get(&id).copied();
                if cached.is_some() {
                    matched.push(*p);
                }
            } else if map.contains_key(&id) {
                matched.push(*p);
            }
        }
        let handle = match cached {
            Some(h) => h,
            None => {
                let h = self.next_handle;
                self.next_handle += 1;
                self.cached.insert(h, CachedReadGroup { reads: Vec::new(), partitions: Vec::new(), refs: 0 });
                h
            }
        };
        // a read that this partition is about to write is shared with the cache in Java and found written
        // when claimed later; the cache holds a stand-in for it
        let read = if will_write { tracker.reads[rid as usize].written_stand_in() } else { tracker.reads[rid as usize].take() };
        let c = self.cached.get_mut(&handle).unwrap();
        c.reads.push(read);
        for p in remote_partitions {
            if !c.partitions.contains(p) {
                c.partitions.push(*p);
            }
        }
        let mut added = 0;
        for p in remote_partitions {
            if !unprocessed.contains(p) || matched.contains(p) {
                continue;
            }
            if self.candidate_partition_groups.entry(*p).or_default().insert(id.clone(), handle).is_none() {
                added += 1;
            }
        }
        let c = self.cached.get_mut(&handle).unwrap();
        c.refs += added;
        if c.refs == 0 {
            self.cached.remove(&handle);
        }
        Ok(())
    }

    fn purge_partition(&mut self, partition: PartKey) {
        self.junction_partition_read_ids.remove(&partition);
        let Some(groups) = self.candidate_partition_groups.remove(&partition) else { return };
        for (_, handle) in groups {
            let purged = self.cached.get(&handle).is_some_and(|c| c.partitions.iter().all(|p| self.processed_partitions.contains(p)));
            if purged {
                self.purged_candidates += 1;
            }
            self.release(handle);
        }
    }

    fn junction_read_id_count(&self) -> usize {
        self.junction_partition_read_ids.values().map(|x| x.len()).sum()
    }

    fn log_cache_count(&mut self, force_log: bool) {
        const LOG_CACHE_DIFF: i64 = 50000;
        if self.use_cache_bam {
            if force_log {
                crate::info!("spanning cache partition processed({}) junctionIds({})", self.processed_partitions.len(), self.junction_read_id_count());
            }
            return;
        }
        // groups that span several partitions are counted once for each
        let new_count: i64 = self.candidate_partition_groups.values().map(|m| m.values().map(|h| self.cached.get(h).map_or(0, |c| c.reads.len() as i64)).sum::<i64>()).sum();
        if (new_count - self.last_snapshot_count).abs() > LOG_CACHE_DIFF || force_log {
            crate::info!(
                "spanning cache partition processed({}) candidates cached({} -> {} matched={} purged={}) junctionIds({})",
                self.processed_partitions.len(),
                self.last_snapshot_count,
                new_count,
                self.matched_candidates,
                self.purged_candidates,
                self.junction_read_id_count()
            );
            self.last_snapshot_count = new_count;
        }
    }

    pub fn log_stats(&mut self) {
        self.log_cache_count(true);
    }

    pub fn reset(&mut self) {
        self.candidate_partition_groups.clear();
        self.cached.clear();
        self.junction_partition_read_ids.clear();
        self.processed_partitions.clear();
    }

    /// For tests: the read ids expected in each partition.
    #[cfg(test)]
    pub fn junction_partition_read_ids(&self) -> &HashMap<PartKey, HashSet<Box<[u8]>>> {
        &self.junction_partition_read_ids
    }

    /// getCachedReadsCount: the cached reads of all groups, or of the groups in partition maps that hold the read id
    #[cfg(test)]
    pub fn get_cached_reads_count(&self, read_id: Option<&[u8]>) -> usize {
        let mut handles: HashSet<usize> = HashSet::new();
        for map in self.candidate_partition_groups.values() {
            if read_id.is_none_or(|id| map.contains_key(id)) {
                handles.extend(map.values());
            }
        }
        handles.iter().map(|h| self.cached.get(h).map_or(0, |c| c.reads.len())).sum()
    }
}
