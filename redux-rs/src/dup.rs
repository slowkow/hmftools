//! Duplicate identification: fragment coordinates, the read cache and duplicate groups
//! (Java FragmentCoords, ReadCache, DuplicateGroup, Statistics).

use crate::rec::*;
use crate::unmap::extract_supp_alignment;
use rust_htslib::bam;
use hashbrown::HashTable;

#[derive(Clone, Debug)]
pub struct FragmentCoords {
    pub pos_lower: i32,
    pub pos_upper: i32,
    pub lower_forward: bool,
    pub upper_forward: bool,
    /// forward = F1R2
    pub frag_forward: bool,
    pub read_is_lower: bool,
    /// the supplementary read's own unclipped 5' position and orientation
    pub supp: Option<(i32, bool)>,
    pub unmapped_sourced: bool,
    pub unpaired: bool,
    pub key: String,
    /// the key without the fragment orientation suffix is its first `non_oriented_len` bytes
    non_oriented_len: usize,
}

fn push_int(s: &mut String, v: i32) {
    let mut buf = [0u8; 12];
    let mut i = buf.len();
    let mut n = (v as i64).unsigned_abs();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if v < 0 {
        i -= 1;
        buf[i] = b'-';
    }
    // digits and a sign are ASCII
    s.push_str(unsafe { std::str::from_utf8_unchecked(&buf[i..]) });
}

/// chr:pos, with :R on the reverse strand
fn push_coordinate(s: &mut String, chr: &str, pos: i32, forward: bool) {
    s.push_str(chr);
    s.push(':');
    push_int(s, pos);
    if !forward {
        s.push_str(":R");
    }
}

impl FragmentCoords {
    #[allow(clippy::too_many_arguments)]
    fn new(
        chr_lower: &str, chr_upper: &str, pos_lower: i32, pos_upper: i32, lower_forward: bool, upper_forward: bool, frag_forward: bool,
        read_is_lower: bool, supp: Option<(i32, bool)>, is_unmapped: bool, key_by_frag_orient: bool, unpaired: bool,
    ) -> Self {
        // the key is built in one allocation: this runs for every read
        let mut key = String::with_capacity(chr_lower.len() + chr_upper.len() + 40);
        if unpaired {
            key.push_str(chr_lower);
            key.push(':');
            push_int(&mut key, pos_lower);
            key.push('_');
            push_int(&mut key, pos_upper);
            if if read_is_lower { !lower_forward } else { !upper_forward } {
                key.push_str("_R");
            }
            if supp.is_some() {
                key.push_str("_S");
            }
        } else {
            push_coordinate(&mut key, chr_lower, pos_lower, lower_forward);
            if pos_upper == NO_POSITION {
                if is_unmapped {
                    key.push_str("_N");
                } else if supp.is_some() {
                    key.push_str("_S");
                }
            } else {
                key.push('_');
                push_coordinate(&mut key, chr_upper, pos_upper, upper_forward);
                key.push_str(if read_is_lower { "_L" } else { "_U" });
                if supp.is_some() {
                    key.push_str("_S");
                }
            }
        }
        let non_oriented_len = key.len();
        if !unpaired && pos_upper != NO_POSITION && key_by_frag_orient && !frag_forward {
            key.push_str("_V");
        }
        FragmentCoords {
            pos_lower,
            pos_upper,
            lower_forward,
            upper_forward,
            frag_forward,
            read_is_lower,
            supp,
            unmapped_sourced: is_unmapped,
            unpaired,
            key,
            non_oriented_len,
        }
    }

    pub fn read_position(&self) -> i32 {
        match self.supp {
            Some((p, _)) => p,
            None => {
                if self.read_is_lower || self.unmapped_sourced { self.pos_lower } else { self.pos_upper }
            }
        }
    }
    pub fn read_forward(&self) -> bool {
        match self.supp {
            Some((_, f)) => f,
            None => {
                if self.read_is_lower || self.unmapped_sourced { self.lower_forward } else { self.upper_forward }
            }
        }
    }
    pub fn key_non_oriented(&self) -> &str {
        &self.key[..self.non_oriented_len]
    }

    /// FragmentCoords.fromRead. `names` are the reference names by index.
    pub fn from_read(rec: &bam::Record, names: &[String], use_fragment_orientation: bool) -> Self {
        Self::from_read_as(rec, names, use_fragment_orientation, has(rec, FLAG_SUPPLEMENTARY))
    }

    /// FragmentCoords.fromRead(read, useFragmentOrientation, isSupplementary): with `is_supplementary` false a
    /// supplementary read gives its own coordinates
    pub fn from_read_as(rec: &bam::Record, names: &[String], use_fragment_orientation: bool, is_supplementary: bool) -> Self {
        let is_paired = has(rec, FLAG_PAIRED);
        let is_unmapped = has(rec, FLAG_UNMAPPED);
        let mate_unmapped = has(rec, FLAG_MATE_UNMAPPED);
        let name = |tid: i32| -> &str { if tid >= 0 { &names[tid as usize] } else { "*" } };

        let mut mate_chr: &str = "*";
        let mut mate_forward = true;
        let mut mate_pos = NO_POSITION;
        if is_paired && !mate_unmapped {
            mate_chr = name(rec.mtid());
            mate_forward = !has(rec, FLAG_MATE_REVERSE);
            // without a mate cigar the mate is taken as fully aligned
            let fallback;
            let mate_cigar: &[u8] = match aux_str(rec, b"MC") {
                Some(mc) => mc,
                None => {
                    fallback = format!("{}M", rec.seq_len());
                    fallback.as_bytes()
                }
            };
            mate_pos = five_prime_unclipped_str(mate_start(rec), mate_cigar, mate_forward);
        }

        let mut supp = None;
        let sa_chr: std::borrow::Cow<str>;
        let read_chr: &str;
        let read_forward: bool;
        let read_pos: i32;
        if !is_paired || !is_unmapped {
            let mut chr: &str = name(rec.tid());
            let mut fwd = !has(rec, FLAG_REVERSE);
            let mut pos = five_prime_unclipped(rec);
            if is_supplementary {
                supp = Some((pos, fwd));
                let sd = extract_supp_alignment(aux_str(rec, b"SA")).expect("supplementary read without supplementary alignment data");
                sa_chr = String::from_utf8_lossy(sd.chromosome);
                chr = &sa_chr;
                fwd = sd.strand == b'+';
                pos = five_prime_unclipped_str(sd.position, sd.cigar, fwd);
                if !is_paired {
                    mate_chr = chr;
                    mate_forward = !fwd;
                    mate_pos = three_prime_unclipped_str(sd.position, sd.cigar, fwd);
                }
            } else if !is_paired {
                mate_chr = chr;
                mate_forward = !fwd;
                mate_pos = three_prime_unclipped(rec);
            }
            read_chr = chr;
            read_forward = fwd;
            read_pos = pos;
        } else {
            read_chr = mate_chr;
            read_forward = mate_forward;
            read_pos = mate_pos;
        }

        if is_unmapped || (is_paired && mate_unmapped) {
            let read_is_lower = !is_unmapped;
            return Self::new(
                read_chr, "*", read_pos, NO_POSITION, read_forward, read_forward, true, read_is_lower, supp, is_unmapped, false, !is_paired,
            );
        }

        let read_is_lower = if read_chr != mate_chr {
            rec.tid() < rec.mtid()
        } else if read_pos != mate_pos {
            read_pos < mate_pos
        } else {
            has(rec, FLAG_FIRST)
        };
        let frag_forward = !is_paired || (read_is_lower == has(rec, FLAG_FIRST));
        if read_is_lower {
            Self::new(
                read_chr, mate_chr, read_pos, mate_pos, read_forward, mate_forward, frag_forward, read_is_lower, supp, is_unmapped,
                use_fragment_orientation, !is_paired,
            )
        } else {
            Self::new(
                mate_chr, read_chr, mate_pos, read_pos, mate_forward, read_forward, frag_forward, read_is_lower, supp, is_unmapped,
                use_fragment_orientation, !is_paired,
            )
        }
    }
}

pub struct ReadInfo {
    pub read: bam::Record,
    pub coords: FragmentCoords,
}

pub struct DuplicateGroup {
    pub umi: Option<String>,
    pub coords: FragmentCoords,
    pub reads: Vec<bam::Record>,
    pub non_consensus_reads: Vec<bam::Record>,
    pub consensus_read: Option<bam::Record>,
    /// index into `reads` of the read kept as primary when no consensus is formed
    pub primary_read: Option<usize>,
    pub dual_strand: bool,
    pub poly_g_unmapped: bool,
}

impl DuplicateGroup {
    pub fn new(reads: Vec<bam::Record>, coords: FragmentCoords) -> Self {
        DuplicateGroup {
            umi: None,
            coords,
            reads,
            non_consensus_reads: Vec::new(),
            consensus_read: None,
            primary_read: None,
            dual_strand: false,
            poly_g_unmapped: false,
        }
    }
    pub fn total_read_count(&self) -> usize {
        self.reads.len() + self.non_consensus_reads.len()
    }
    /// DuplicateGroup.calculateFirstInPairCount
    pub fn first_in_pair_count(&self) -> i32 {
        let first = |r: &bam::Record| !has(r, FLAG_PAIRED) || has(r, FLAG_FIRST);
        let mut count = self.reads.iter().filter(|r| first(r)).count() + self.non_consensus_reads.iter().filter(|r| first(r)).count();
        let total = self.total_read_count();
        if !(count as f32 >= 0.5f32 * total as f32) {
            count = total - count;
        }
        count as i32
    }
}

pub struct FragmentCoordReads {
    pub duplicate_groups: Vec<DuplicateGroup>,
    pub single_reads: Vec<ReadInfo>,
}

impl FragmentCoordReads {
    pub fn total_read_count(&self) -> usize {
        self.duplicate_groups.iter().map(|g| g.total_read_count()).sum::<usize>() + self.single_reads.len()
    }
    /// FragmentCoordReads.minReadPositionStart
    pub fn min_read_position_start(&self) -> i32 {
        let mut min = self.single_reads.iter().map(|r| start(&r.read)).min().unwrap_or(-1);
        for g in &self.duplicate_groups {
            let gmin = g.reads.iter().map(start).min().unwrap_or(-1);
            if min < 0 || gmin > 0 && gmin < min {
                min = gmin;
            }
        }
        min
    }
}

/// String.hashCode over the key, then HashMap's spread of the high bits
fn java_spread_hash(key: &str) -> u32 {
    let mut h: i32 = 0;
    for &b in key.as_bytes() {
        h = h.wrapping_mul(31).wrapping_add(b as i32);
    }
    let h = h as u32;
    h ^ (h >> 16)
}

/// Spreads Java's 32-bit hash over the 64 bits the table uses
fn table_hash(hash: u32) -> u64 {
    (hash as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

struct Entry {
    coords: FragmentCoords,
    reads: Vec<bam::Record>,
    hash: u32,
    seq: u64,
}

/// ReadPositionGroup with a map that can be walked in the order of a java.util.HashMap: by bucket of the
/// table's current capacity, then by insertion. (Buckets that Java would turn into trees are not modelled.)
struct ReadPositionGroup {
    tid: i32,
    start: i32,
    end: i32,
    /// keyed by the fragment key held in each entry, hashed with Java's hash
    map: HashTable<Entry>,
    capacity: u32,
    next_seq: u64,
}

impl ReadPositionGroup {
    fn add_read(&mut self, coords: FragmentCoords, read: bam::Record) {
        let hash = java_spread_hash(&coords.key);
        if let Some(e) = self.map.find_mut(table_hash(hash), |e| e.hash == hash && e.coords.key == coords.key) {
            e.reads.push(read);
            return;
        }
        let seq = self.next_seq;
        self.next_seq += 1;
        self.map.insert_unique(table_hash(hash), Entry { coords, reads: vec![read], hash, seq }, |e| table_hash(e.hash));
        if self.map.len() as u32 > self.capacity / 4 * 3 {
            self.capacity *= 2;
        }
    }
    fn read_count(&self) -> usize {
        self.map.iter().map(|e| e.reads.len()).sum()
    }
    fn min_read_position(&self) -> i32 {
        let mut min = -1;
        for e in self.map.iter() {
            let m = e.reads.iter().map(start).min().unwrap_or(-1);
            if min < 0 || m < min {
                min = m;
            }
        }
        min
    }
}

pub const DEFAULT_GROUP_SIZE: i32 = 200;
pub const DEFAULT_MAX_SOFT_CLIP: i32 = 130;
pub const DEFAULT_POP_DISTANCE_CHECK: i32 = 100;
pub const DEFAULT_LOG_READ_COUNT_THRESHOLD: i32 = 100_000;
pub const DEFAULT_DYNAMIC_READ_COUNT_THRESHOLD: i32 = 100_000;
pub const MAX_DISTANCE_FOR_REVERSE_STRAND_DUPLICATES: i32 = 1000;

pub struct ReadCache {
    group_size: i32,
    max_soft_clip: i32,
    default_distance_check: i32,
    pop_distance_check: i32,
    use_fragment_orientation: bool,
    read_count_reduction_threshold: i32,
    check_read_count_reduction: bool,
    current_read_min_position: i32,
    current_tid: i32,
    groups: Vec<ReadPositionGroup>,
    pop_position_check: i32,
    check_size_read_count: i32,
    last_cache_read_count: i32,
    log_cache_read_count: i32,
    check_log_cache_read_count: i32,
    /// processRead calls made on this cache
    processed: u64,
    sync: Option<(std::sync::Arc<CacheSync>, usize)>,
    synced: bool,
    /// R14 as in Java: the pop distance never leaves its default
    java_fixed_pop_distance: bool,
}

/// Java keeps one ReadCache per thread and never resets its size-check counter between partitions, so
/// with one thread the counter at any read is the number of reads cached in all earlier partitions,
/// modulo the check interval. The counter decides when a cache above the read count threshold is
/// purged, which changes duplicate groups in very deep regions. To give Java's one-thread result with
/// partitions running in parallel, a partition whose cache approaches the threshold waits for the
/// earlier partitions' counts before the counter can matter.
pub struct CacheSync {
    /// per partition: reads cached, and the last cache read count if it was changed
    done: std::sync::Mutex<Vec<Option<(u64, Option<i32>)>>>,
    cv: std::sync::Condvar,
}

impl CacheSync {
    pub fn new(partitions: usize) -> std::sync::Arc<Self> {
        std::sync::Arc::new(CacheSync { done: std::sync::Mutex::new(vec![None; partitions]), cv: std::sync::Condvar::new() })
    }

    /// Marks a partition finished without a cache (also used when it failed, so that nothing waits on it)
    pub fn publish_empty(&self, index: usize) {
        let mut d = self.done.lock().unwrap();
        if d[index].is_none() {
            d[index] = Some((0, None));
        }
        self.cv.notify_all();
    }

    fn publish(&self, index: usize, reads: u64, last: Option<i32>) {
        self.done.lock().unwrap()[index] = Some((reads, last));
        self.cv.notify_all();
    }

    /// Waits for every earlier partition: the total reads they cached and the last cache read count
    fn incoming(&self, index: usize) -> (u64, i32) {
        let mut d = self.done.lock().unwrap();
        while d[..index].iter().any(|x| x.is_none()) {
            d = self.cv.wait(d).unwrap();
        }
        let total = d[..index].iter().map(|x| x.unwrap().0).sum();
        let last = d[..index].iter().rev().find_map(|x| x.unwrap().1).unwrap_or(0);
        (total, last)
    }
}

/// how often, in reads, the cache size is looked at to decide whether to wait for earlier partitions
const SYNC_CHECK_READS: u64 = 512;

impl ReadCache {
    /// Ties this cache to its partition's place in the run (see CacheSync)
    pub fn set_sync(&mut self, sync: std::sync::Arc<CacheSync>, index: usize) {
        // R21 corrected: the counters start again with each region (see evict_all), so there is nothing to wait for
        self.synced = index == 0 || !crate::compat::java("R21");
        self.sync = Some((sync, index));
    }

    /// Call when the partition is complete
    pub fn publish_sync(&mut self) {
        if let Some((sync, index)) = self.sync.take() {
            sync.publish(index, self.processed, if self.synced { Some(self.last_cache_read_count) } else { None });
        }
    }

    fn check_sync(&mut self) {
        if self.synced || self.log_cache_read_count == 0 || self.processed % SYNC_CHECK_READS != 0 {
            return;
        }
        let Some((sync, index)) = &self.sync else { return };
        // the cache grows by at most one read per call, so no size check can see the threshold before this does
        if (self.cached_read_count() as i64) < self.log_cache_read_count as i64 - SYNC_CHECK_READS as i64 - 1 {
            return;
        }
        let (total, last) = sync.incoming(*index);
        let interval = self.check_log_cache_read_count.max(1) as u64;
        self.check_size_read_count = ((total + self.processed) % interval) as i32;
        self.last_cache_read_count = last;
        self.synced = true;
    }

    pub fn new(
        group_size: i32, max_soft_clip: i32, use_fragment_orientation: bool, pop_distance_check: i32, dynamic_read_count_threshold: i32,
        log_cache_read_count: i32,
    ) -> Self {
        ReadCache {
            group_size,
            max_soft_clip,
            default_distance_check: pop_distance_check,
            pop_distance_check,
            use_fragment_orientation,
            read_count_reduction_threshold: dynamic_read_count_threshold,
            check_read_count_reduction: false,
            current_read_min_position: 0,
            current_tid: -2,
            groups: Vec::new(),
            pop_position_check: 0,
            check_size_read_count: 0,
            last_cache_read_count: 0,
            log_cache_read_count,
            check_log_cache_read_count: (log_cache_read_count as f64 * 0.1).round() as i32,
            processed: 0,
            sync: None,
            synced: true,
            java_fixed_pop_distance: crate::compat::java("R14"),
        }
    }

    pub fn illumina(use_fragment_orientation: bool) -> Self {
        Self::new(
            DEFAULT_GROUP_SIZE, DEFAULT_MAX_SOFT_CLIP, use_fragment_orientation, DEFAULT_POP_DISTANCE_CHECK,
            DEFAULT_DYNAMIC_READ_COUNT_THRESHOLD, DEFAULT_LOG_READ_COUNT_THRESHOLD,
        )
    }

    pub fn current_read_min_position(&self) -> i32 {
        self.current_read_min_position
    }

    /// ReadCache.processRead
    pub fn process_read(&mut self, read: bam::Record, names: &[String]) {
        self.check_sync();
        self.processed += 1;
        let left_clip = left_soft_clip(read.raw_cigar());
        if left_clip > self.max_soft_clip {
            crate::warn!(
                "read({}) leftClip({}) exceeds maxSoftClipLength({})",
                String::from_utf8_lossy(read.qname()), left_clip, self.max_soft_clip
            );
        }
        let coords = FragmentCoords::from_read(&read, names, self.use_fragment_orientation);
        let tid = read.tid();
        let pos = start(&read);
        let gi = self.get_or_create_group(tid, coords.read_position());
        self.groups[gi].add_read(coords, read);
        self.current_read_min_position = pos;
        if self.current_tid != tid {
            self.current_tid = tid;
            self.pop_position_check = 0;
        }
        self.check_cache_size(false);
    }

    fn get_or_create_group(&mut self, tid: i32, fragment_position: i32) -> usize {
        let mut index = 0;
        while index < self.groups.len() {
            let g = &self.groups[index];
            if g.tid != tid {
                index += 1;
                continue;
            }
            if fragment_position >= g.start && fragment_position <= g.end {
                return index;
            }
            if fragment_position < g.start {
                break;
            }
            index += 1;
        }
        let gs = self.group_size;
        // Java % and floorDiv
        let remainder = fragment_position % gs;
        let (start, end);
        if remainder == 0 {
            end = gs * fragment_position.div_euclid(gs);
            start = end - gs + 1;
        } else {
            start = gs * fragment_position.div_euclid(gs) + 1;
            end = start + gs - 1;
        }
        self.groups.insert(index, ReadPositionGroup { tid, start, end, map: HashTable::new(), capacity: 16, next_seq: 0 });
        index
    }

    /// ReadCache.popReads
    pub fn pop_reads(&mut self) -> Option<FragmentCoordReads> {
        if let Some(r) = self.pop_passed_reads() {
            return Some(r);
        }
        if self.check_read_count_reduction && self.last_cache_read_count >= self.read_count_reduction_threshold {
            return self.purge_above_read_reduction_threshold();
        }
        None
    }

    fn collect(mut popped: Vec<((usize, u32, u64), Entry)>) -> Option<FragmentCoordReads> {
        if popped.is_empty() {
            return None;
        }
        popped.sort_by_key(|p| p.0);
        let mut duplicate_groups = Vec::new();
        let mut single_reads = Vec::new();
        for (_, e) in popped {
            if e.reads.len() > 1 {
                duplicate_groups.push(DuplicateGroup::new(e.reads, e.coords));
            } else {
                let mut reads = e.reads;
                single_reads.push(ReadInfo { read: reads.pop().unwrap(), coords: e.coords });
            }
        }
        Some(FragmentCoordReads { duplicate_groups, single_reads })
    }

    fn pop_passed_reads(&mut self) -> Option<FragmentCoordReads> {
        if let Some(g) = self.groups.first() {
            if g.tid == self.current_tid && self.current_read_min_position < self.pop_position_check + self.pop_distance_check {
                return None;
            }
        }
        self.pop_position_check = (self.current_read_min_position / self.pop_distance_check) * self.pop_distance_check;
        let pop_lower = self.current_read_min_position - self.max_soft_clip;
        let pop_upper = self.current_read_min_position - 1;
        let mut popped: Vec<((usize, u32, u64), Entry)> = Vec::new();
        let mut group_index = 0;
        // position in the original list, to keep the popped fragments in list order
        let mut list_index = 0;
        while group_index < self.groups.len() {
            let current_tid = self.current_tid;
            let current_pos = self.current_read_min_position;
            let group = &mut self.groups[group_index];
            let take_all = group.tid != current_tid || group.end < pop_lower;
            let cap_mask = group.capacity - 1;
            if take_all {
                for e in group.map.drain() {
                    popped.push(((list_index, e.hash & cap_mask, e.seq), e));
                }
            } else {
                let passed = |e: &mut Entry| {
                    let read_position = e.coords.read_position();
                    if e.coords.read_forward() {
                        read_position <= pop_lower
                    } else {
                        let last = e.reads.last().unwrap();
                        !(read_position > pop_upper && start(last) + MAX_DISTANCE_FOR_REVERSE_STRAND_DUPLICATES > current_pos)
                    }
                };
                for e in group.map.extract_if(passed) {
                    popped.push(((list_index, e.hash & cap_mask, e.seq), e));
                }
            }
            list_index += 1;
            if group.map.is_empty() {
                self.groups.remove(group_index);
                continue;
            }
            group_index += 1;
        }
        Self::collect(popped)
    }

    /// ReadCache.evictAll
    pub fn evict_all(&mut self) -> Option<FragmentCoordReads> {
        let mut popped = Vec::new();
        for (gi, group) in self.groups.iter_mut().enumerate() {
            let cap_mask = group.capacity - 1;
            for e in group.map.drain() {
                popped.push(((gi, e.hash & cap_mask, e.seq), e));
            }
        }
        self.groups.clear();
        if !crate::compat::java("R21") {
            // R21: Java leaves the size-check counter and the last counted size as they are when a region ends (its
            // clear(), which resets them, is never called), so the purges of a deep region depend on the reads of
            // every region the thread took before. Corrected: each region starts as a new cache does.
            self.check_size_read_count = 0;
            self.last_cache_read_count = 0;
            self.check_read_count_reduction = false;
            self.pop_distance_check = self.default_distance_check;
        }
        self.pop_position_check = (self.current_read_min_position / self.pop_distance_check) * self.pop_distance_check;
        Self::collect(popped)
    }

    fn purge_above_read_reduction_threshold(&mut self) -> Option<FragmentCoordReads> {
        self.check_read_count_reduction = false;
        if self.groups.iter().any(|g| g.tid != self.current_tid) {
            return None;
        }
        let cached_read_count = self.cached_read_count() as i64;
        let cached_frag_count = self.cached_frag_coord_groups() as i64;
        let average_group_size = cached_read_count as f64 / cached_frag_count as f64;
        if average_group_size < 2.0 {
            return None;
        }
        let target_read_count = cached_read_count as f64 * 0.5;
        let group_size_threshold = (average_group_size * 0.5).max(2.0);
        let mut reads_popped: i64 = 0;
        let mut popped: Vec<((usize, u32, u64), Entry)> = Vec::new();
        let mut group_index = 0;
        let mut list_index = 0;
        while group_index < self.groups.len() {
            if (cached_read_count - reads_popped) as f64 <= target_read_count {
                break;
            }
            let group = &mut self.groups[group_index];
            let cap_mask = group.capacity - 1;
            let mut removed = false;
            // java.util.HashMap iteration order: by bucket of the current capacity, then by insertion
            let mut order: Vec<(u32, u64, u32, usize)> =
                group.map.iter().map(|e| (e.hash & cap_mask, e.seq, e.hash, e.reads.len())).collect();
            order.sort_unstable();
            for (_, seq, hash, len) in order {
                if (len as f64) < group_size_threshold {
                    continue;
                }
                let Ok(found) = group.map.find_entry(table_hash(hash), |e| e.seq == seq) else { unreachable!() };
                let e = found.remove().0;
                reads_popped += e.reads.len() as i64;
                popped.push(((list_index, e.hash & cap_mask, e.seq), e));
                removed = true;
                if (cached_read_count - reads_popped) as f64 <= target_read_count {
                    break;
                }
            }
            list_index += 1;
            if removed && group.map.is_empty() {
                self.groups.remove(group_index);
                continue;
            }
            group_index += 1;
        }
        if reads_popped > 0 {
            self.check_cache_size(true);
        }
        // every fragment purged here is a duplicate group, also one of a single read
        if popped.is_empty() {
            return None;
        }
        popped.sort_by_key(|p| p.0);
        Some(FragmentCoordReads {
            duplicate_groups: popped.into_iter().map(|(_, e)| DuplicateGroup::new(e.reads, e.coords)).collect(),
            single_reads: Vec::new(),
        })
    }

    pub fn min_cached_read_start(&self) -> i32 {
        self.groups.iter().map(|g| g.min_read_position()).min().unwrap_or(-1)
    }

    fn check_cache_size(&mut self, force_update: bool) {
        if self.log_cache_read_count == 0 {
            return;
        }
        if !force_update {
            self.check_size_read_count += 1;
            if self.check_size_read_count < self.check_log_cache_read_count {
                return;
            }
            self.check_size_read_count = 0;
        }
        let new_read_count = self.cached_read_count() as i32;
        if new_read_count >= self.log_cache_read_count && (new_read_count - self.last_cache_read_count).abs() >= self.log_cache_read_count {
            self.last_cache_read_count = new_read_count;
            self.check_read_count_reduction = true;
        }
        if force_update {
            self.last_cache_read_count = new_read_count;
        }
        self.check_pos_check_adjust(new_read_count);
    }

    fn check_pos_check_adjust(&mut self, cached_read_count: i32) {
        // R14: Java returns whenever the distance is the default, and nothing else changes it, so the distance is
        // never shortened. Corrected (intent uncertain): it returns only when there is also nothing to shorten.
        if self.pop_distance_check == self.default_distance_check
            && (self.java_fixed_pop_distance || cached_read_count < self.read_count_reduction_threshold)
        {
            return;
        }
        if cached_read_count < self.read_count_reduction_threshold {
            self.pop_distance_check = self.default_distance_check;
        } else {
            let multiple = ((cached_read_count as f64 / self.read_count_reduction_threshold as f64).round() as i32).min(9);
            let adjust_fraction = (10 - multiple) as f64 / 10.0;
            self.pop_distance_check = (self.default_distance_check as f64 * adjust_fraction).round() as i32;
        }
    }

    /// For tests: the pop distance, the size-check counter and the last counted size.
    #[cfg(test)]
    pub fn check_state(&self) -> (i32, i32, i32) {
        (self.pop_distance_check, self.check_size_read_count, self.last_cache_read_count)
    }

    #[cfg(test)]
    pub fn cached_read_position_groups(&self) -> usize {
        self.groups.len()
    }
    pub fn cached_read_count(&self) -> usize {
        self.groups.iter().map(|g| g.read_count()).sum()
    }
    pub fn cached_frag_coord_groups(&self) -> usize {
        self.groups.iter().map(|g| g.map.len()).sum()
    }
}

/// Statistics: duplicate group sizes, keyed by rounded size; the size written is the first one seen for the key.
#[derive(Default, Clone)]
pub struct Statistics {
    pub total_reads: i64,
    pub duplicate_reads: i64,
    pub duplicate_groups: i64,
    /// (rounded size, first size seen, count, dual strand count)
    pub frequencies: Vec<(i32, i32, i64, i64)>,
    /// Statistics.Ultima
    pub ultima: crate::ultima::UltimaStats,
    pub umi: crate::umi::UmiStatistics,
}

/// Statistics.roundFrequency
pub fn round_frequency(frequency: i32) -> i32 {
    let round = |x: f32| (x + 0.5).floor() as i32;
    if frequency <= 10 {
        frequency
    } else if frequency <= 100 {
        round(frequency as f32 / 10.0) * 10
    } else if frequency <= 1000 {
        round(frequency as f32 / 100.0) * 100
    } else {
        round(frequency as f32 / 1000.0) * 1000
    }
}

impl Statistics {
    fn add_frequency(&mut self, duplicate_count: i32, count: i64, dual: i64) {
        let rounded = round_frequency(duplicate_count);
        match self.frequencies.iter_mut().find(|f| f.0 == rounded) {
            Some(f) => {
                f.2 += count;
                f.3 += dual;
            }
            None => self.frequencies.push((rounded, duplicate_count, count, dual)),
        }
    }
    pub fn add_non_duplicate_counts(&mut self, fragment_count: usize) {
        if fragment_count > 0 {
            self.add_frequency(1, fragment_count as i64, 0);
        }
    }
    pub fn add_duplicate_group(&mut self, fragment_count: usize) {
        self.add_frequency(fragment_count as i32, 1, 0);
        self.duplicate_groups += 1;
        self.duplicate_reads += fragment_count as i64 - 1;
    }
    /// Statistics.addUmiGroup
    pub fn add_umi_group(&mut self, fragment_count: usize, has_dual_strand: bool) {
        self.add_frequency(fragment_count as i32, 1, has_dual_strand as i64);
        self.duplicate_groups += 1;
        self.duplicate_reads += fragment_count as i64 - 1;
    }
    pub fn merge(&mut self, other: &Statistics) {
        self.umi.merge(&other.umi);
        self.total_reads += other.total_reads;
        self.duplicate_reads += other.duplicate_reads;
        self.duplicate_groups += other.duplicate_groups;
        self.ultima.merge(&other.ultima);
        // Java walks the other map's values in HashMap order of the Integer keys, which matters only for
        // which first-seen size a new rounded key takes; a key has one value in `other`, so order is free
        for f in &other.frequencies {
            self.add_frequency(f.1, f.2, f.3);
        }
    }
    /// Statistics.writeDuplicateStats
    pub fn write_duplicate_stats(&self, filename: &str, include_dual_strand: bool) -> anyhow::Result<()> {
        let mut out = String::from("DuplicateReadCount\tFrequency");
        if include_dual_strand {
            out.push_str("\tDualStrandFrequency");
        }
        out.push('\n');
        let mut f = self.frequencies.clone();
        f.sort_by_key(|x| x.0);
        // R18: Java labels a row of sizes above 10 with the first size it saw for the row, which depends on the
        // order of the groups (and, with threads, on which thread took which partition). Corrected: the rounded size.
        let java_label = crate::compat::java("R18");
        for (rounded, first, count, dual) in f {
            let label = if java_label { first } else { rounded };
            out.push_str(&format!("{label}\t{count}"));
            if include_dual_strand {
                out.push_str(&format!("\t{dual}"));
            }
            out.push('\n');
        }
        std::fs::write(filename, out)?;
        Ok(())
    }
}
