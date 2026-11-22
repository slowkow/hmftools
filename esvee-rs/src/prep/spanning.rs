//! SpanningReadCache and CandidateBamWriter: what one partition leaves for the others. A partition with a
//! junction fragment records its read id against the partitions of its remote reads, where the read is then
//! "expected"; a partition with only candidate reads of a fragment parks them, in memory or in per-chromosome
//! cache BAMs, until a junction elsewhere claims them.
//!
//! Java's two ways of parking differ in what gets claimed (finding P5): from the cache BAMs a read is claimed
//! when a junction fragment of its id points to the read's chromosome, and the read's filters are lost; in
//! memory a read is claimed by the partition it waits for. With P5 corrected the cache BAMs only hold the
//! records and both modes claim the same reads, by the rule of the in-memory mode (`ClaimBook`).
//!
//! The in-memory mode applies that rule as the partitions are processed. The cache BAM mode keeps no
//! bookkeeping per waiting read while it runs, since on a whole genome some ten million fragments wait at
//! once and next to none are ever claimed. Each cached record instead carries what the rule needs to know of
//! it (`CACHE_EVENT_TAG`), and a junction fragment is noted only when a candidate read of its id may be
//! waiting for its partition (a set of 32-bit hashes per partition answers that, with no false "no"). At the
//! end the cache BAMs are read, the records of the noted read ids are picked out, and the rule is applied to
//! them and the noted junction fragments in the order in which they happened; a second reading fetches the
//! claimed records. The bookkeeping of one read id does not depend on any other, so leaving out the ids that
//! no junction fragment asked for changes nothing.

use super::config::*;
use super::filters::{check_filters, has_line_tail, ReadFilterConfig};
use super::jmap::JavaHashMap;
use super::rec::*;
use super::tracker::{JunctionTracker, ReadIdTrimmer};
use super::types::*;
use super::writer::{RawBamOut, ResultsWriter};
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use crate::common::compat;

pub const SAMPLE_ID_TAG: &[u8; 2] = b"SI";

/// On a record of a cache BAM when the claims are those of the in-memory mode: an array of 32-bit numbers,
/// the event number (low, high), flags (`EVENT_CLAIMABLE`), the offset of the read id in the record's name,
/// then chromosome, partition index and "still to be processed" for each remote partition of the group, in
/// the order of the group.
pub const CACHE_EVENT_TAG: &[u8; 2] = b"ze";
const EVENT_CLAIMABLE: u32 = 1;
/// the records that a thread reading a cache BAM hands over at a time
const CACHE_BATCH: usize = 1024;
/// the batches that may wait: of the cache BAM whose turn it is, and of all of them
const HANDOVER_CURRENT: usize = 16;
const HANDOVER_AHEAD: usize = 1024;

struct CachedReadGroup<T> {
    items: Vec<T>,
    partitions: Vec<PartKey>,
    /// the number of partition maps that hold the group
    refs: usize,
}

/// The claim rule of Java's in-memory cache: a candidate read waits under the remote partitions of its group
/// that are still to come, and is claimed when one of them has a junction fragment of the same read id. `T`
/// is what waits: the read itself, or a stand-in for a record of a cache BAM.
struct ClaimBook<T> {
    /// by partition, then read id: the cached group (an index into `cached`)
    candidate_partition_groups: HashMap<PartKey, HashMap<Box<[u8]>, usize>>,
    cached: HashMap<usize, CachedReadGroup<T>>,
    next_handle: usize,
    matched_candidates: i64,
    purged_candidates: i64,
}

impl<T> ClaimBook<T> {
    fn new() -> Self {
        ClaimBook { candidate_partition_groups: HashMap::new(), cached: HashMap::new(), next_handle: 0, matched_candidates: 0, purged_candidates: 0 }
    }

    fn clear(&mut self) {
        self.candidate_partition_groups.clear();
        self.cached.clear();
    }

    fn release(&mut self, handle: usize) {
        if let Some(c) = self.cached.get_mut(&handle) {
            c.refs -= 1;
            if c.refs == 0 {
                self.cached.remove(&handle);
            }
        }
    }

    /// A junction fragment in `source`: what waits there under its read id, taken out of the book.
    fn take_matched(&mut self, source: PartKey, id: &[u8]) -> Option<Vec<T>> {
        let handle = self.candidate_partition_groups.get_mut(&source).and_then(|m| m.remove(id))?;
        let (items, partitions) = {
            let c = self.cached.get_mut(&handle).unwrap();
            (std::mem::take(&mut c.items), c.partitions.clone())
        };
        self.matched_candidates += items.len() as i64;
        self.release(handle);
        for p in partitions {
            if p != source {
                if let Some(h) = self.candidate_partition_groups.get_mut(&p).and_then(|m| m.remove(id)) {
                    self.release(h);
                }
            }
        }
        Some(items)
    }

    /// A read of a group of candidate reads with the given remote partitions, of which `unprocessed` are
    /// still to come. `item` is absent for a read that only links the group's partitions.
    fn add_candidate(&mut self, id: &[u8], remote_partitions: &[PartKey], unprocessed: &[PartKey], item: Option<T>) {
        // an existing cached group for the fragment, under any of its remote partitions
        let mut cached: Option<usize> = None;
        let mut matched: Vec<PartKey> = Vec::new();
        for p in remote_partitions {
            let Some(map) = self.candidate_partition_groups.get(p) else { continue };
            if cached.is_none() {
                cached = map.get(id).copied();
                if cached.is_some() {
                    matched.push(*p);
                }
            } else if map.contains_key(id) {
                matched.push(*p);
            }
        }
        let handle = match cached {
            Some(h) => h,
            None => {
                let h = self.next_handle;
                self.next_handle += 1;
                self.cached.insert(h, CachedReadGroup { items: Vec::new(), partitions: Vec::new(), refs: 0 });
                h
            }
        };
        let c = self.cached.get_mut(&handle).unwrap();
        c.items.extend(item);
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
            if self.candidate_partition_groups.entry(*p).or_default().insert(id.into(), handle).is_none() {
                added += 1;
            }
        }
        let c = self.cached.get_mut(&handle).unwrap();
        c.refs += added;
        if c.refs == 0 {
            self.cached.remove(&handle);
        }
    }

    /// The partition has been processed: nothing waits for it any more. `processed` is for the count of
    /// groups that no partition can claim now.
    fn purge_partition(&mut self, partition: PartKey, processed: Option<&HashSet<PartKey>>) {
        let Some(groups) = self.candidate_partition_groups.remove(&partition) else { return };
        for (_, handle) in groups {
            if let Some(processed) = processed {
                if self.cached.get(&handle).is_some_and(|c| c.partitions.iter().all(|p| processed.contains(p))) {
                    self.purged_candidates += 1;
                }
            }
            self.release(handle);
        }
    }

    /// The waiting items; groups that span several partitions are counted once for each.
    fn cached_count(&self) -> i64 {
        self.candidate_partition_groups.values().map(|m| m.values().map(|h| self.cached.get(h).map_or(0, |c| c.items.len() as i64)).sum::<i64>()).sum()
    }
}

/// A junction fragment whose partition may have candidate reads of its id waiting for it.
struct JunctionEvent {
    event: u64,
    id: Box<[u8]>,
    source: PartKey,
    /// the partition writes the fragment, so what waits for it is written too
    claim: bool,
}

/// A record of a cache BAM whose read id a junction fragment asked for.
struct CachedCandidate {
    event: u64,
    chr: i32,
    /// the record's number in its file
    ordinal: u32,
    /// the junction event that has its read id
    id_of: u32,
    /// the remote partitions of its group, each with whether it was still to be processed
    partitions: Vec<(PartKey, bool)>,
    claimable: bool,
}

/// What `in_order` has in hand.
struct Handover<T> {
    /// per item: what its thread has sent and is not yet received
    queues: Vec<VecDeque<Result<T>>>,
    finished: Vec<bool>,
    /// in all queues
    waiting: usize,
    /// the item being received
    current: usize,
}

/// Runs `work` on the items on up to `threads` threads, each item taken in turn, and hands what it sends
/// to `receive` item by item in the order of the items. The threads of later items work ahead until
/// `HANDOVER_AHEAD` batches wait in all, so that reading is not held up by the order and little is held.
fn in_order<I: Sync, T: Send>(
    items: &[I], threads: usize, work: impl Fn(&I, &dyn Fn(T)) -> Result<()> + Sync, mut receive: impl FnMut(T) -> Result<()>,
) -> Result<()> {
    let n = items.len();
    let state = Mutex::new(Handover { queues: (0..n).map(|_| VecDeque::new()).collect(), finished: vec![false; n], waiting: 0, current: 0 });
    let changed = Condvar::new();
    let next = AtomicUsize::new(0);
    let send = |i: usize, t: Result<T>| {
        let mut s = state.lock().unwrap();
        while !(if i == s.current { s.queues[i].len() < HANDOVER_CURRENT } else { s.waiting < HANDOVER_AHEAD }) {
            s = changed.wait(s).unwrap();
        }
        s.queues[i].push_back(t);
        s.waiting += 1;
        changed.notify_all();
    };
    std::thread::scope(|scope| {
        for _ in 0..threads.min(n).max(1) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= n {
                    break;
                }
                if let Err(e) = work(&items[i], &|t| send(i, Ok(t))) {
                    send(i, Err(e));
                }
                state.lock().unwrap().finished[i] = true;
                changed.notify_all();
            });
        }
        // after a failure the rest is still taken, which lets the threads run out
        let mut result = Ok(());
        for i in 0..n {
            let mut s = state.lock().unwrap();
            s.current = i;
            changed.notify_all();
            loop {
                if let Some(t) = s.queues[i].pop_front() {
                    s.waiting -= 1;
                    changed.notify_all();
                    drop(s);
                    if result.is_ok() {
                        result = t.and_then(&mut receive);
                    }
                    s = state.lock().unwrap();
                } else if s.finished[i] {
                    break;
                } else {
                    s = changed.wait(s).unwrap();
                }
            }
        }
        result
    })
}

fn id_hash(id: &[u8]) -> u32 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in id {
        h = (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3);
    }
    (h ^ (h >> 32)) as u32
}

pub struct CandidateBamWriter {
    /// the cache BAMs hold the records and the claims are those of the in-memory mode (P5 corrected)
    token_mode: bool,
    no_clean_up: bool,
    file_prefix: String,
    sample_ids: Vec<String>,
    header: Option<bam::HeaderView>,
    writers: HashMap<i32, RawBamOut>,
    files: HashMap<i32, String>,
    /// the events of the claim rule so far: cached reads and noted junction fragments (P5 corrected)
    events: u64,
    /// by partition still to be processed: the hashes of the read ids with candidate reads waiting for it
    waiting_ids: HashMap<PartKey, Vec<u32>>,
    /// `waiting_ids` of the partition being processed, sorted
    source_waiting_ids: Vec<u32>,
    junction_events: Vec<JunctionEvent>,
    /// each partition processed: the event count when it was done
    steps: Vec<(u64, PartKey)>,
    threads: usize,
    filter: ReadFilterConfig,
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

    fn cache_writer(&mut self, chr: i32) -> Result<&mut RawBamOut> {
        if !self.writers.contains_key(&chr) {
            let file = format!("{}_{}.bam", self.file_prefix, self.contigs.name(chr));
            let mut out = RawBamOut::create(&file, 1)?;
            if let Some(h) = &self.header {
                out.write_header(h)?;
            }
            self.writers.insert(chr, out);
            self.files.insert(chr, file);
        }
        Ok(self.writers.get_mut(&chr).unwrap())
    }

    /// Java's cache BAM mode: the record with its sample.
    fn write_candidate_read(&mut self, read: &PrepRead) -> Result<()> {
        let mut rec = read.rec.clone();
        set_attribute(&mut rec, SAMPLE_ID_TAG, Some(&encode_z(self.sample_ids[read.sample as usize].as_bytes())));
        self.cache_writer(read.chr)?.write(&rec)
    }

    /// P5 corrected: the record with its sample and its part in the claim rule. `id` is the read id of the
    /// group. A read that the partition is itself about to write (`stand_in`) is filed without its content:
    /// it cannot be claimed, but links the partitions of its group as the other reads do.
    fn write_waiting_read(&mut self, read: &PrepRead, id: &[u8], remote_partitions: &[PartKey], unprocessed: &[PartKey], stand_in: bool, claimable: bool) -> Result<()> {
        let event = self.note_waiting_read(id, unprocessed);
        let mut rec = if stand_in {
            let mut rec = bam::Record::new();
            rec.set(id, None, &[], &[]);
            rec.set_flags(FLAG_UNMAPPED);
            rec.set_tid(-1);
            rec.set_pos(-1);
            rec.set_mtid(-1);
            rec.set_mpos(-1);
            rec
        } else {
            read.rec.clone()
        };
        let id_offset = rec.qname().len().saturating_sub(id.len());
        if &rec.qname()[id_offset..] != id {
            anyhow::bail!("a cached read's name does not end with the read id of its group");
        }
        let numbers = 4 + 3 * remote_partitions.len();
        let mut value: Vec<u8> = Vec::with_capacity(6 + 4 * numbers);
        value.extend_from_slice(b"BI");
        value.extend_from_slice(&(numbers as u32).to_le_bytes());
        let flags = if claimable && !stand_in { EVENT_CLAIMABLE } else { 0 };
        for n in [event as u32, (event >> 32) as u32, flags, id_offset as u32] {
            value.extend_from_slice(&n.to_le_bytes());
        }
        for p in remote_partitions {
            for n in [p.chr as u32, p.part as u32, unprocessed.contains(p) as u32] {
                value.extend_from_slice(&n.to_le_bytes());
            }
        }
        // as set_attribute, for both attributes at once
        let mut aux = aux_sorted_copy(aux_bytes(&rec));
        aux_block_set(&mut aux, SAMPLE_ID_TAG, Some(&encode_z(self.sample_ids[read.sample as usize].as_bytes())));
        aux_block_set(&mut aux, CACHE_EVENT_TAG, Some(&value));
        set_aux_bytes(&mut rec, &aux);
        self.cache_writer(read.chr)?.write(&rec)
    }

    /// P5 corrected: a candidate read waits for the partitions of `unprocessed`. Returns its event number.
    fn note_waiting_read(&mut self, id: &[u8], unprocessed: &[PartKey]) -> u64 {
        let hash = id_hash(id);
        for p in unprocessed {
            let ids = self.waiting_ids.entry(*p).or_default();
            if ids.last() != Some(&hash) {
                ids.push(hash);
            }
        }
        self.events += 1;
        self.events - 1
    }

    /// P5 corrected: a partition is about to be processed.
    fn begin_partition(&mut self, source: PartKey) {
        self.source_waiting_ids = self.waiting_ids.remove(&source).unwrap_or_default();
        self.source_waiting_ids.sort_unstable();
        self.source_waiting_ids.dedup();
    }

    fn end_partition(&mut self, source: PartKey) {
        self.source_waiting_ids = Vec::new();
        self.steps.push((self.events, source));
    }

    /// P5 corrected: a junction fragment of the partition being processed. It is noted when candidate reads
    /// of its id may wait for this partition; when none can, the claim rule has nothing to do for it.
    fn note_junction_fragment(&mut self, id: &[u8], source: PartKey, claim: bool) {
        if self.source_waiting_ids.binary_search(&id_hash(id)).is_ok() {
            self.junction_events.push(JunctionEvent { event: self.events, id: id.into(), source, claim });
            self.events += 1;
        }
    }

    /// The sample of a cached record, with the attributes of the cache taken off.
    fn take_sample(sample_ids: &[String], rec: &mut bam::Record) -> u8 {
        let sample = aux_str(rec, SAMPLE_ID_TAG).and_then(|s| sample_ids.iter().position(|x| x.as_bytes() == s)).unwrap_or(0) as u8;
        set_attribute(rec, CACHE_EVENT_TAG, None);
        set_attribute(rec, SAMPLE_ID_TAG, None);
        sample
    }

    /// The records of one cache BAM whose read id is in `wanted`, sent in batches.
    fn find_cached_candidates(file: &str, chr: i32, wanted: &HashMap<&[u8], u32>, send: &dyn Fn(Vec<CachedCandidate>)) -> Result<()> {
        let mut reader = bam::Reader::from_path(file)?;
        let mut found = Vec::new();
        let mut ordinal = 0u32;
        let mut rec = bam::Record::new();
        while let Some(r) = reader.read(&mut rec) {
            r?;
            ordinal += 1;
            let Some(f) = aux_field(&rec, CACHE_EVENT_TAG).filter(|f| f.len() >= 22 && f[0] == b'B' && f[1] == b'I' && (f.len() - 22) % 12 == 0) else {
                anyhow::bail!("a record of {file} lacks the attribute of the candidate cache");
            };
            let number = |i: usize| u32::from_le_bytes([f[6 + 4 * i], f[7 + 4 * i], f[8 + 4 * i], f[9 + 4 * i]]);
            let Some(id) = rec.qname().get(number(3) as usize..) else {
                anyhow::bail!("a record of {file} has an invalid attribute of the candidate cache");
            };
            let Some(&id_of) = wanted.get(id) else { continue };
            found.push(CachedCandidate {
                event: number(0) as u64 | (number(1) as u64) << 32,
                chr,
                ordinal: ordinal - 1,
                id_of,
                partitions: (0..(f.len() - 22) / 12).map(|k| (PartKey { chr: number(4 + 3 * k) as i32, part: number(5 + 3 * k) as i32 }, number(6 + 3 * k) != 0)).collect(),
                claimable: number(2) & EVENT_CLAIMABLE != 0,
            });
            if found.len() == CACHE_BATCH {
                send(std::mem::take(&mut found));
            }
        }
        send(found);
        Ok(())
    }

    /// The records of one cache BAM with the given numbers (in ascending order), each made a read again with
    /// the filters it had when cached, sent in batches.
    fn read_claimed_records(
        file: &str, claimed: &[u32], sample_ids: &[String], tid_map: &[i32], contigs: &Contigs, filter: &ReadFilterConfig, send: &dyn Fn(Vec<PrepRead>),
    ) -> Result<()> {
        let mut reader = bam::Reader::from_path(file)?;
        let mut found = Vec::new();
        let mut next = claimed.iter().copied().peekable();
        let mut ordinal = 0u32;
        let mut rec = bam::Record::new();
        while next.peek().is_some() {
            match reader.read(&mut rec) {
                Some(r) => r?,
                None => break,
            }
            ordinal += 1;
            if next.peek() != Some(&(ordinal - 1)) {
                continue;
            }
            next.next();
            let mut rec = std::mem::replace(&mut rec, bam::Record::new());
            let sample = Self::take_sample(sample_ids, &mut rec);
            let line_tail = has_line_tail(&rec);
            let mut read = PrepRead::new(rec, sample, tid_map, contigs);
            read.has_line_tail = line_tail;
            check_filters(&mut read, filter);
            read.set_read_type(ReadType::CandidateSupport);
            found.push(read);
            if found.len() == CACHE_BATCH {
                send(std::mem::take(&mut found));
            }
        }
        send(found);
        Ok(())
    }

    /// The claim rule, applied to the noted junction fragments and the cached reads of their read ids in the
    /// order in which they happened; `steps` gives the points at which a partition was done. Returns the
    /// claimed ones of `candidates`.
    fn resolve_claims(junction_events: &[JunctionEvent], steps: &[(u64, PartKey)], candidates: &[CachedCandidate]) -> Vec<usize> {
        let mut order: Vec<(u64, bool, usize)> =
            candidates.iter().enumerate().map(|(i, c)| (c.event, true, i)).chain(junction_events.iter().enumerate().map(|(i, e)| (e.event, false, i))).collect();
        order.sort_unstable();
        let mut book: ClaimBook<usize> = ClaimBook::new();
        let mut claimed: Vec<usize> = Vec::new();
        let mut next = 0;
        for &(end, source) in steps {
            while next < order.len() && order[next].0 < end {
                let (_, is_candidate, i) = order[next];
                next += 1;
                if is_candidate {
                    let c = &candidates[i];
                    let remote_partitions: Vec<PartKey> = c.partitions.iter().map(|p| p.0).collect();
                    let unprocessed: Vec<PartKey> = c.partitions.iter().filter(|p| p.1).map(|p| p.0).collect();
                    book.add_candidate(&junction_events[c.id_of as usize].id, &remote_partitions, &unprocessed, c.claimable.then_some(i));
                } else {
                    let e = &junction_events[i];
                    if let Some(items) = book.take_matched(e.source, &e.id) {
                        if e.claim {
                            claimed.extend(items);
                        }
                    }
                }
            }
            book.purge_partition(source, None);
        }
        claimed
    }

    /// Writes the claimed records of the cache BAMs (P5 corrected), file by file in chromosome order with
    /// the file of unmapped reads last. A read gets back the filters it had when cached, so one with a
    /// poly-G soft clip is held back as it is everywhere else.
    fn write_claimed_reads(&mut self, writer: &mut ResultsWriter, tid_map: &[i32]) -> Result<()> {
        self.waiting_ids = HashMap::new();
        let junction_events = std::mem::take(&mut self.junction_events);
        let steps = std::mem::take(&mut self.steps);
        crate::debug!("{} junction fragments may have candidate reads waiting, of {} cached reads and such fragments", junction_events.len(), self.events);
        if junction_events.is_empty() {
            return Ok(());
        }
        let mut files: Vec<(i32, &str)> = self.files.iter().map(|(chr, file)| (*chr, file.as_str())).collect();
        files.sort_by_key(|&(chr, _)| (chr < 0, chr));

        // the cached records of the read ids that a junction fragment asked for
        let mut wanted: HashMap<&[u8], u32> = HashMap::with_capacity(junction_events.len());
        for (i, e) in junction_events.iter().enumerate() {
            wanted.entry(&e.id[..]).or_insert(i as u32);
        }
        let mut candidates: Vec<CachedCandidate> = Vec::new();
        in_order(&files, self.threads, |&(chr, file), send| Self::find_cached_candidates(file, chr, &wanted, send), |mut found| {
            candidates.append(&mut found);
            Ok(())
        })?;
        drop(wanted);

        let file_of: HashMap<i32, usize> = files.iter().enumerate().map(|(i, &(chr, _))| (chr, i)).collect();
        let mut claimed: Vec<Vec<u32>> = vec![Vec::new(); files.len()];
        let claimed_candidates = Self::resolve_claims(&junction_events, &steps, &candidates);
        crate::debug!("{} of {} cached reads looked at are claimed", claimed_candidates.len(), candidates.len());
        for c in claimed_candidates {
            claimed[file_of[&candidates[c].chr]].push(candidates[c].ordinal);
        }
        drop((candidates, junction_events));

        let mut claimed: Vec<(&str, Vec<u32>)> = files.iter().zip(claimed).filter(|(_, c)| !c.is_empty()).map(|(&(_, file), c)| (file, c)).collect();
        for (_, c) in &mut claimed {
            c.sort_unstable();
            c.dedup();
        }
        let (sample_ids, contigs, filter) = (&self.sample_ids[..], &*self.contigs, &self.filter);
        in_order(&claimed, self.threads, |(file, c), send| Self::read_claimed_records(file, c, sample_ids, tid_map, contigs, filter, send), |reads| {
            for read in reads {
                writer.write_candidate_read(read)?;
            }
            Ok(())
        })
    }

    /// CandidateReadMatchTask.call for one cache file: Java's rule, by read id and the record's chromosome.
    fn match_file_reads(&mut self, file: &str, writer: &mut ResultsWriter, tid_map: &[i32]) -> Result<u64> {
        let mut reader = bam::Reader::from_path(file)?;
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
            let sample = Self::take_sample(&self.sample_ids, &mut rec);
            // Java (P5): the read is made anew from the record, without its filters
            let mut read = PrepRead::new(rec, sample, tid_map, &self.contigs);
            read.set_read_type(ReadType::CandidateSupport);
            writer.write_candidate_read(read)?;
        }
        Ok(matched)
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
        if self.token_mode {
            self.write_claimed_reads(writer, tid_map)?;
        } else {
            let chromosomes: Vec<String> = self.chr_junction_read_ids.keys().cloned().collect();
            for chromosome in chromosomes {
                let chr = self.contigs.id(&chromosome);
                let Some(file) = self.files.get(&chr).cloned() else {
                    continue;
                };
                crate::debug!("chr({}) assigning candidates from {} junction fragments", chromosome, self.chr_junction_read_ids.get(&chromosome).map_or(0, |s| s.len()));
                let matched = self.match_file_reads(&file, writer, tid_map)?;
                crate::debug!("chr({chromosome}) matched and wrote {matched} candidate reads");
            }
            if !compat::java("P6") {
                // Java: files unmapped reads under chromosome "-1", which no junction read id is held against,
                // so that file is never read. An unmapped read is placed with its mate, and is claimed here
                // like the other reads: by its id and the chromosome it is placed on.
                if let Some(file) = self.files.get(&UNMAPPED_CHR).cloned() {
                    let matched = self.match_file_reads(&file, writer, tid_map)?;
                    crate::debug!("matched and wrote {matched} unmapped candidate reads");
                }
            }
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
    /// Java's cache BAM mode: claims by read id and chromosome (finding P5)
    java_cache_bam: bool,
    /// unmapped reads in the cache BAMs are never claimed (finding P6)
    java_unmapped: bool,
    specific_chromosomes: Vec<String>,
    contigs: Arc<Contigs>,
    /// the candidate reads held in memory (`-no_cache_bam`)
    book: ClaimBook<PrepRead>,
    junction_partition_read_ids: HashMap<PartKey, HashSet<Box<[u8]>>>,
    processed_partitions: HashSet<PartKey>,
    last_snapshot_count: i64,
    pub candidate_bam_writer: CandidateBamWriter,
}

impl SpanningReadCache {
    pub fn new(config: &PrepConfig, contigs: Arc<Contigs>, cache_header: Option<bam::HeaderView>) -> Self {
        let java_cache_bam = config.use_cache_bam && compat::java("P5");
        SpanningReadCache {
            partition_size: config.partition_size,
            use_cache_bam: config.use_cache_bam,
            java_cache_bam,
            java_unmapped: compat::java("P6"),
            specific_chromosomes: config.specific.chromosomes.clone(),
            contigs: contigs.clone(),
            book: ClaimBook::new(),
            junction_partition_read_ids: HashMap::new(),
            processed_partitions: HashSet::new(),
            last_snapshot_count: 0,
            candidate_bam_writer: CandidateBamWriter {
                token_mode: config.use_cache_bam && !java_cache_bam,
                no_clean_up: config.no_clean_up,
                file_prefix: if config.sample_ids.is_empty() { String::new() } else { config.form_filename(WriteType::CacheBam, config.sample_id()) },
                sample_ids: config.sample_ids.clone(),
                header: cache_header,
                writers: HashMap::new(),
                files: HashMap::new(),
                events: 0,
                waiting_ids: HashMap::new(),
                source_waiting_ids: Vec::new(),
                junction_events: Vec::new(),
                steps: Vec::new(),
                threads: config.threads,
                filter: config.filter.clone(),
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
    /// with whether the partition is about to write it while it awaits remote reads, and whether it is one of
    /// the partition's junction groups (the groups the partition writes).
    pub fn process_spanning_read_groups(&mut self, region: &Region, spanning: &[(GroupId, bool, bool)], tracker: &mut JunctionTracker) -> Result<()> {
        let source = self.chr_partition(region.chr, region.start);
        self.processed_partitions.insert(source);
        if self.candidate_bam_writer.token_mode {
            self.candidate_bam_writer.begin_partition(source);
        }
        let contigs = self.contigs.clone();
        for &(gid, will_write, from_junctions) in spanning {
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
                self.process_junction_read(&unprocessed, &remote_partitions, gid, source, from_junctions, tracker);
            }
        }
        self.purge_partition(source);
        self.log_cache_count(false);
        Ok(())
    }

    fn process_junction_read(
        &mut self, unprocessed: &[PartKey], remote_partitions: &[PartKey], gid: GroupId, source: PartKey, from_junctions: bool, tracker: &mut JunctionTracker,
    ) {
        let first_chr = tracker.reads[tracker.groups[gid as usize].reads[0] as usize].chr;
        if self.ignore_chromosome(self.contigs.name(first_chr)) {
            return;
        }
        let id = tracker.groups[gid as usize].id.clone();
        if self.java_cache_bam {
            self.candidate_bam_writer.add_junction_read_id(remote_partitions, &id);
        } else if self.candidate_bam_writer.token_mode {
            // Candidate reads found earlier for this fragment are written with the group in the in-memory
            // mode, which happens when the group is one of the partition's junction groups and not a
            // duplicate. Records in the cache BAMs are written at the end, on the same condition.
            let claim = from_junctions && tracker.groups[gid as usize].status != ReadGroupStatus::Duplicate;
            self.candidate_bam_writer.note_junction_fragment(&id, source, claim);
        } else if let Some(reads) = self.book.take_matched(source, &id) {
            // candidate reads found earlier for this fragment join the group
            for read in reads {
                let rid = tracker.add_read(read);
                tracker.groups[gid as usize].reads.push(rid);
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
        if self.java_cache_bam {
            return self.candidate_bam_writer.write_candidate_read(&tracker.reads[rid as usize]);
        }
        // a read that this partition is about to write is shared with the cache in Java and found written
        // when claimed later; the cache holds a stand-in for it
        if self.use_cache_bam {
            let read = &tracker.reads[rid as usize];
            // Java (P6): an unmapped read is filed where it is never looked for
            let claimable = !(self.java_unmapped && read.chr == UNMAPPED_CHR);
            return self.candidate_bam_writer.write_waiting_read(read, &tracker.groups[gid as usize].id, remote_partitions, unprocessed, will_write, claimable);
        }
        let read = if will_write { tracker.reads[rid as usize].written_stand_in() } else { tracker.reads[rid as usize].take() };
        self.book.add_candidate(&tracker.groups[gid as usize].id, remote_partitions, unprocessed, Some(read));
        Ok(())
    }

    fn purge_partition(&mut self, partition: PartKey) {
        self.junction_partition_read_ids.remove(&partition);
        if self.candidate_bam_writer.token_mode {
            self.candidate_bam_writer.end_partition(partition);
        } else {
            self.book.purge_partition(partition, Some(&self.processed_partitions));
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
        let new_count = self.book.cached_count();
        if (new_count - self.last_snapshot_count).abs() > LOG_CACHE_DIFF || force_log {
            crate::info!(
                "spanning cache partition processed({}) candidates cached({} -> {} matched={} purged={}) junctionIds({})",
                self.processed_partitions.len(),
                self.last_snapshot_count,
                new_count,
                self.book.matched_candidates,
                self.book.purged_candidates,
                self.junction_read_id_count()
            );
            self.last_snapshot_count = new_count;
        }
    }

    pub fn log_stats(&mut self) {
        self.log_cache_count(true);
    }

    pub fn reset(&mut self) {
        self.book.clear();
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
        for map in self.book.candidate_partition_groups.values() {
            if read_id.is_none_or(|id| map.contains_key(id)) {
                handles.extend(map.values());
            }
        }
        handles.iter().map(|h| self.book.cached.get(h).map_or(0, |c| c.items.len())).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cache BAM mode notes only the junction fragments that may have candidate reads waiting, and at the
    /// end applies the claim rule to the cached reads of their read ids alone. That must claim what the rule
    /// claims when it is applied to every read as it comes, as the in-memory mode does.
    #[test]
    fn deferred_claims_equal_immediate_claims() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut random = move |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        // partitions 0..10 of the first chromosome are processed in turn; those of the second never are
        let partition = |i: usize| if i < 10 { PartKey { chr: 0, part: i as i32 } } else { PartKey { chr: 1, part: (i - 10) as i32 } };
        let mut total_claimed = 0;
        for _ in 0..300 {
            let mut config = PrepConfig::for_tests(1000);
            config.use_cache_bam = true;
            let mut deferred = SpanningReadCache::new(&config, Arc::new(Contigs::new()), None).candidate_bam_writer;
            assert!(deferred.token_mode);
            let mut immediate: ClaimBook<usize> = ClaimBook::new();
            let mut immediate_claimed: Vec<usize> = Vec::new();
            // every cached read: its event, read id, remote partitions, whether it can be claimed
            let mut cached: Vec<(u64, Vec<u8>, Vec<(PartKey, bool)>, bool)> = Vec::new();
            for step in 0..10 {
                let source = partition(step);
                deferred.begin_partition(source);
                for id in 0..12 {
                    let read_id = format!("READ_{id}").into_bytes();
                    match random(4) {
                        0 => {
                            // a junction fragment
                            let claim = random(4) != 0;
                            deferred.note_junction_fragment(&read_id, source, claim);
                            if let Some(items) = immediate.take_matched(source, &read_id) {
                                if claim {
                                    immediate_claimed.extend(items);
                                }
                            }
                        }
                        1 | 2 => {
                            // a group of candidate reads
                            let mut remote: Vec<PartKey> = Vec::new();
                            for _ in 0..1 + random(3) {
                                let p = partition(random(12));
                                if p != source && !remote.contains(&p) {
                                    remote.push(p);
                                }
                            }
                            let unprocessed: Vec<PartKey> = remote.iter().copied().filter(|p| p.chr == 1 || p.part > step as i32).collect();
                            if unprocessed.is_empty() {
                                continue;
                            }
                            let stand_in = random(5) == 0;
                            for _ in 0..1 + random(2) {
                                let event = deferred.note_waiting_read(&read_id, &unprocessed);
                                let partitions = remote.iter().map(|p| (*p, unprocessed.contains(p))).collect();
                                immediate.add_candidate(&read_id, &remote, &unprocessed, (!stand_in).then_some(cached.len()));
                                cached.push((event, read_id.clone(), partitions, !stand_in));
                            }
                        }
                        _ => {}
                    }
                }
                deferred.end_partition(source);
                immediate.purge_partition(source, None);
            }
            // what the pass over the cache BAMs picks out: the cached reads of the noted read ids
            let junction_events = std::mem::take(&mut deferred.junction_events);
            let mut candidates: Vec<CachedCandidate> = Vec::new();
            let mut cached_index: Vec<usize> = Vec::new();
            for (i, (event, read_id, partitions, claimable)) in cached.iter().enumerate() {
                if let Some(id_of) = junction_events.iter().position(|e| &e.id[..] == &read_id[..]) {
                    candidates.push(CachedCandidate { event: *event, chr: 0, ordinal: i as u32, id_of: id_of as u32, partitions: partitions.clone(), claimable: *claimable });
                    cached_index.push(i);
                }
            }
            let deferred_claimed: Vec<usize> = CandidateBamWriter::resolve_claims(&junction_events, &deferred.steps, &candidates).into_iter().map(|c| cached_index[c]).collect();
            assert_eq!(immediate_claimed, deferred_claimed);
            total_claimed += immediate_claimed.len();
        }
        assert!(total_claimed > 1000);
    }
}
