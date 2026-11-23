//! UMI duplicate groups (Java UmiConfig, UmiUtils, UmiGroupBuilder, PolyGUmiMerger, JitterUmiMerger,
//! UmiStatistics, PositionFragmentCounts and the UMI branch of DuplicateGroupBuilder).
//!
//! The Java code passes groups and single reads around by reference and merges them in place. Here the fragments
//! of one call live in an arena and the lists hold indices, so that "the same object" is "the same index" and
//! reads move from one group to another without being copied.

use crate::dup::{round_frequency, DuplicateGroup, FragmentCoordReads, FragmentCoords, ReadInfo, Statistics};
use crate::rec::*;
use anyhow::{bail, Result};
use rust_htslib::bam;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::atomic::{AtomicUsize, Ordering};

pub const DEFAULT_MAX_UMI_BASE_DIFF: i32 = 1;
pub const MAX_IMBALANCED_UMI_BASE_DIFF: i32 = 4;
pub const MAX_IMBALANCED_UMI_COUNT: usize = 25;
pub const MIN_POLYG_UMI_TAIL_LENGTH: usize = 2;
pub const DEFAULT_DUPLEX_UMI_DELIM: &str = "_";
pub const SINGLE_END_JITTER_COLLAPSE_DISTANCE: i32 = 10;
pub const READ_ID_DELIM: u8 = b':';
const MAX_EDIT_DISTANCE: usize = 10;

/// String.hashCode of an ASCII string, then java.util.HashMap's spread of the high bits
fn java_spread_hash(key: &[u8]) -> u32 {
    let mut h: i32 = 0;
    for &b in key {
        h = h.wrapping_mul(31).wrapping_add(b as i32);
    }
    let h = h as u32;
    h ^ (h >> 16)
}

/// A map that can be walked in the order of a java.util.HashMap with String keys: by bucket of the table's
/// capacity, then by insertion. (Buckets that Java would turn into trees are not modelled.) The caller gives
/// the hash of the Java key, so that the key itself need not be the formatted string.
struct JavaMap<K, V> {
    index: HashMap<K, usize>,
    entries: Vec<(u32, Option<V>)>,
    capacity: u32,
    live: usize,
}

impl<K: Hash + Eq, V> JavaMap<K, V> {
    fn new() -> Self {
        JavaMap { index: HashMap::new(), entries: Vec::new(), capacity: 16, live: 0 }
    }
    fn find<Q: Hash + Eq + ?Sized>(&self, key: &Q) -> Option<usize>
    where
        K: std::borrow::Borrow<Q>,
    {
        self.index.get(key).copied()
    }
    fn insert(&mut self, key: K, hash: u32, value: V) -> usize {
        let i = self.entries.len();
        self.entries.push((hash, Some(value)));
        self.index.insert(key, i);
        self.live += 1;
        if self.live as u32 > self.capacity / 4 * 3 {
            self.capacity *= 2;
        }
        i
    }
    fn value(&self, i: usize) -> &V {
        self.entries[i].1.as_ref().unwrap()
    }
    fn value_mut(&mut self, i: usize) -> &mut V {
        self.entries[i].1.as_mut().unwrap()
    }
    /// HashMap.remove: the table keeps its capacity
    fn remove(&mut self, i: usize) {
        if self.entries[i].1.take().is_some() {
            self.live -= 1;
        }
    }
    fn is_empty(&self) -> bool {
        self.live == 0
    }
    /// entry indices in java.util.HashMap iteration order
    fn order(&self) -> Vec<usize> {
        let mask = self.capacity - 1;
        let mut v: Vec<usize> = (0..self.entries.len()).filter(|&i| self.entries[i].1.is_some()).collect();
        if v.len() > 1 {
            v.sort_by_key(|&i| (self.entries[i].0 & mask, i));
        }
        v
    }
}

/// A short formatted key on the stack, for its Java hash
struct KeyBuf {
    buf: [u8; 24],
    len: usize,
}

impl KeyBuf {
    fn new() -> Self {
        KeyBuf { buf: [0; 24], len: 0 }
    }
    fn push(mut self, c: u8) -> Self {
        self.buf[self.len] = c;
        self.len += 1;
        self
    }
    /// %d
    fn int(mut self, v: i32) -> Self {
        let mut digits = [0u8; 11];
        let mut n = 0;
        let mut x = (v as i64).abs();
        loop {
            digits[n] = b'0' + (x % 10) as u8;
            n += 1;
            x /= 10;
            if x == 0 {
                break;
            }
        }
        if v < 0 {
            self = self.push(b'-');
        }
        for i in (0..n).rev() {
            self = self.push(digits[i]);
        }
        self
    }
    fn hash(&self) -> u32 {
        java_spread_hash(&self.buf[..self.len])
    }
}

/// UmiConfig
pub struct UmiConfig {
    pub enabled: bool,
    pub duplex: bool,
    pub duplex_delim: String,
    pub base_stats: bool,
    pub uses_fixed_umi: bool,
    pub permitted_base_diff: i32,
    /// the length of the first UMI extracted, used for every read from then on; 0 until set
    umi_length: AtomicUsize,
    /// the defined UMIs in the iteration order of the Java HashSet
    defined: Vec<String>,
    defined_set: HashSet<Vec<u8>>,
}

impl UmiConfig {
    pub fn new(enabled: bool, duplex: bool, duplex_delim: &str, uses_fixed: bool, base_stats: bool) -> Self {
        UmiConfig {
            enabled,
            duplex,
            duplex_delim: duplex_delim.to_string(),
            base_stats,
            uses_fixed_umi: uses_fixed,
            permitted_base_diff: DEFAULT_MAX_UMI_BASE_DIFF,
            umi_length: AtomicUsize::new(0),
            defined: Vec::new(),
            defined_set: HashSet::new(),
        }
    }

    /// UmiConfig.from: the flags, or the settings implied by -umi_type
    pub fn from_args(args: &crate::args::Args) -> Result<Self> {
        let mut enabled = args.has_flag("umi_enabled");
        let mut duplex = args.has_flag("umi_duplex");
        let mut delim = args.value_or("umi_duplex_delim", DEFAULT_DUPLEX_UMI_DELIM).to_string();
        let mut uses_fixed = !args.has_flag("umi_variable_length");
        if let Some(umi_type) = args.value("umi_type") {
            match umi_type {
                "TWIST_DUPLEX" | "KAPA" => {
                    enabled = true;
                    duplex = true;
                    delim = "_".to_string();
                }
                "TSO500_DUPLEX" => {
                    enabled = true;
                    duplex = true;
                    delim = "+".to_string();
                }
                "MSK" => {
                    enabled = true;
                    duplex = false;
                    uses_fixed = false;
                }
                "SINGLE" => {
                    enabled = true;
                    duplex = false;
                }
                "NONE" | "OTHER_DUPLEX" => {}
                other => bail!("invalid umi_type: {other}"),
            }
        }
        let mut config = UmiConfig::new(enabled, duplex, &delim, uses_fixed, args.has_flag("umi_base_diff_stats"));
        if let Some(filename) = args.value("umi_defined_ids") {
            // Java logs a failure to read the file and carries on without defined UMIs
            match std::fs::read_to_string(filename) {
                Ok(text) => {
                    let umis: Vec<&str> = text.lines().collect();
                    crate::info!("loaded {} defined UMI IDs from {}", umis.len(), filename);
                    config.add_defined_umis(&umis);
                }
                Err(e) => crate::error!("failed to load defined UMI IDs: {e}"),
            }
        }
        Ok(config)
    }

    pub fn has_defined_umis(&self) -> bool {
        !self.defined.is_empty()
    }

    /// UmiConfig.addDefinedUmis, given the UMIs in the order they were added to the Java set. The set's iteration
    /// order (by hash bucket, then insertion) decides which defined UMI a read with a base difference matches.
    pub fn add_defined_umis(&mut self, umis: &[&str]) {
        let mut map: JavaMap<String, ()> = JavaMap::new();
        for u in self.defined.iter().map(String::as_str).chain(umis.iter().copied()) {
            if map.find(u).is_none() {
                map.insert(u.to_string(), java_spread_hash(u.as_bytes()), ());
            }
        }
        let mut keys: Vec<Option<String>> = vec![None; map.entries.len()];
        for (k, i) in map.index.drain() {
            keys[i] = Some(k);
        }
        self.defined = map.order().into_iter().map(|i| keys[i].take().unwrap()).collect();
        self.defined_set = self.defined.iter().map(|u| u.as_bytes().to_vec()).collect();
    }

    /// UmiConfig.extractUmiId: the last characters of the read name, as many as the first UMI seen had, or with
    /// variable lengths the item after the last ':'
    pub fn extract_umi_id<'a>(&self, read_id: &'a [u8]) -> &'a [u8] {
        if self.uses_fixed_umi {
            let umi_length = self.check_and_get_umi_length(read_id);
            return extract_umi_id_fixed(read_id, umi_length);
        }
        extract_umi_id_from_read_id(read_id)
    }

    /// UmiConfig.matchDefinedUmiId: the UMI itself when defined, else the first defined UMI within the permitted
    /// base difference
    pub fn match_defined_umi_id(&self, umi_id: &[u8]) -> Option<&str> {
        if self.defined_set.contains(umi_id) {
            return self.defined.iter().find(|d| d.as_bytes() == umi_id).map(String::as_str);
        }
        self.defined.iter().find(|d| !exceeds_umi_id_diff(umi_id, d.as_bytes(), self.permitted_base_diff)).map(String::as_str)
    }

    /// UmiConfig.checkAndGetUmiLength: set from the first read name whose UMI is extracted. With several threads
    /// that is whichever partition gets there first, which only matters if the UMIs differ in length.
    fn check_and_get_umi_length(&self, read_id: &[u8]) -> usize {
        let current = self.umi_length.load(Ordering::Relaxed);
        if current > 0 {
            return current;
        }
        let length = extract_umi_id_from_read_id(read_id).len();
        match self.umi_length.compare_exchange(0, length, Ordering::Relaxed, Ordering::Relaxed) {
            Ok(_) => length,
            Err(existing) => existing,
        }
    }
}

/// UmiConfig.extractUmiId(readId, umiLength): the last characters of the read name. (Java throws when the read
/// name is shorter than the UMI length; this takes the whole name.)
pub fn extract_umi_id_fixed(read_id: &[u8], umi_length: usize) -> &[u8] {
    &read_id[read_id.len().saturating_sub(umi_length)..]
}

/// UmiConfig.extractUmiIdFromReadId: the item after the last ':', or the whole name without one
pub fn extract_umi_id_from_read_id(read_id: &[u8]) -> &[u8] {
    match read_id.iter().rposition(|&c| c == READ_ID_DELIM) {
        Some(i) => &read_id[i + 1..],
        None => read_id,
    }
}

/// UmiUtils.exceedsUmiIdDiff
pub fn exceeds_umi_id_diff(first: &[u8], second: &[u8], permitted_diff: i32) -> bool {
    if first.len() != second.len() {
        return true;
    }
    let mut diffs = 0;
    for i in 0..first.len() {
        if first[i] != second[i] {
            diffs += 1;
            if diffs > permitted_diff {
                return true;
            }
        }
    }
    false
}

/// UmiUtils.calcUmiIdDiff: zero for UMIs of different lengths
pub fn calc_umi_id_diff(first: &[u8], second: &[u8]) -> usize {
    if first.len() != second.len() {
        return 0;
    }
    first.iter().zip(second).filter(|(a, b)| a != b).count()
}

/// UmiGroupBuilder.hasDuplexUmiMatch: the first half of one against the second half of the other, both ways.
/// Java's substring throws on a UMI without the delimiter, which ends the run; so does the error here.
pub fn has_duplex_umi_match(first: &str, second: &str, duplex_delim: &str, permitted_diff: i32) -> Result<bool> {
    if !crate::compat::java("R25") && (!first.contains(duplex_delim) || !second.contains(duplex_delim)) {
        // R25 corrected: a UMI without the delimiter has no two halves, so it is no duplex match
        static WARNED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if !WARNED.swap(true, Ordering::Relaxed) {
            crate::warn!("duplex UMI without the delimiter({duplex_delim}): such fragments are not collapsed as duplex; check -umi_duplex_delim");
        }
        return Ok(false);
    }
    let (a1, a2) = split_umi(first, duplex_delim)?;
    let (b1, b2) = split_umi(second, duplex_delim)?;
    Ok(!exceeds_umi_id_diff(a1, b2, permitted_diff) && !exceeds_umi_id_diff(a2, b1, permitted_diff))
}

/// UmiGroupBuilder.splitUmi. The second half starts one character after the start of the delimiter, whatever
/// the delimiter's length.
fn split_umi<'a>(duplex_umi: &'a str, duplex_delim: &str) -> Result<(&'a [u8], &'a [u8])> {
    let Some(i) = duplex_umi.find(duplex_delim) else {
        bail!("duplex UMI({duplex_umi}) has no delimiter({duplex_delim})");
    };
    let bytes = duplex_umi.as_bytes();
    // R10: Java starts the second half one character after the start of the delimiter
    let delim_length = if crate::compat::java("R10") { 1 } else { duplex_delim.len() };
    Ok((&bytes[..i], &bytes[(i + delim_length).min(bytes.len())..]))
}

/// PositionFragmentCounts
#[derive(Clone, Debug)]
pub struct PositionFragmentCounts {
    pub unique_coord_count: i32,
    pub unique_fragment_count: i32,
    pub frequency: i32,
    pub max_coord_umi_count: i32,
    pub max_umi_reads_count: i32,
    pub umi_group_details: String,
}

/// UmiStatistics.UmiGroupCounts
#[derive(Clone, Debug)]
pub struct UmiGroupCounts {
    pub duplicate_group_count: i32,
    pub read_count: i32,
    pub group_count: i32,
    pub edit_distance_frequency: [i32; MAX_EDIT_DISTANCE + 1],
}

/// How often each kind of merge between fragments happened. Not a Java statistic: logged at debug level, as a
/// check on what a run exercised.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct UmiMergeCounts {
    /// fragments merged into another by PolyGUmiMerger
    pub poly_g_merged: i64,
    /// groups of unmapped reads marked as duplicates as a whole by PolyGUmiMerger
    pub poly_g_unmapped_groups: i64,
    /// pairs of fragments of opposite orientation collapsed on their duplex UMIs
    pub duplex_collapsed: i64,
    /// fragments merged into another by JitterUmiMerger
    pub jitter_merged: i64,
}

/// UmiStatistics
#[derive(Clone, Default)]
pub struct UmiStatistics {
    pub merge_counts: UmiMergeCounts,
    pub umi_groups: i64,
    pub umi_group_frequencies: Vec<UmiGroupCounts>,
    pub position_fragments: Vec<PositionFragmentCounts>,
    /// the length of the first UMI recorded: Java sizes its table by it and ignores bases beyond it
    base_frequency_length: Option<usize>,
    /// counts of A, C, G, T, N by UMI position, for every position seen
    base_frequencies: Vec<[i32; 5]>,
}

impl UmiStatistics {
    fn position_fragment_data(&mut self, unique_coord_count: i32, unique_fragment_count: i32) -> &mut PositionFragmentCounts {
        let i = match self
            .position_fragments
            .iter()
            .position(|x| x.unique_coord_count == unique_coord_count && x.unique_fragment_count == unique_fragment_count)
        {
            Some(i) => i,
            None => {
                self.position_fragments.push(PositionFragmentCounts {
                    unique_coord_count,
                    unique_fragment_count,
                    frequency: 0,
                    max_coord_umi_count: 0,
                    max_umi_reads_count: 0,
                    umi_group_details: String::new(),
                });
                self.position_fragments.len() - 1
            }
        };
        &mut self.position_fragments[i]
    }

    /// UmiStatistics.recordFragmentPositions
    fn record_fragment_positions(
        &mut self, unique_coord_count: i32, unique_fragment_count: i32, max_coord_umi_count: i32, max_group: Option<&DuplicateGroup>,
    ) {
        let data = self.position_fragment_data(unique_coord_count, unique_fragment_count);
        data.frequency += 1;
        data.max_coord_umi_count = data.max_coord_umi_count.max(max_coord_umi_count);
        if let Some(group) = max_group {
            if group.total_read_count() as i32 > data.max_umi_reads_count {
                data.max_umi_reads_count = group.total_read_count() as i32;
                data.umi_group_details = format!("{} {}", group.coords.key_non_oriented(), String::from_utf8_lossy(group.reads[0].qname()));
            }
        }
    }

    fn umi_group_counts(&mut self, group_count: i32, read_count: usize) -> &mut UmiGroupCounts {
        let rounded = round_frequency(read_count as i32);
        let i = match self.umi_group_frequencies.iter().position(|x| x.duplicate_group_count == group_count && x.read_count == rounded) {
            Some(i) => i,
            None => {
                self.umi_group_frequencies.push(UmiGroupCounts {
                    duplicate_group_count: group_count,
                    read_count: rounded,
                    group_count: 0,
                    edit_distance_frequency: [0; MAX_EDIT_DISTANCE + 1],
                });
                self.umi_group_frequencies.len() - 1
            }
        };
        &mut self.umi_group_frequencies[i]
    }

    /// UmiStatistics.recordUmiBaseStats: the groups (and single fragments) of one orientation at one coordinate
    fn record_umi_base_stats(&mut self, config: &UmiConfig, groups: &[&DuplicateGroup]) {
        for g in groups {
            self.record_umi_base_frequencies(g.umi.as_deref().unwrap_or("").as_bytes());
        }
        // the edit distance of each read's UMI to the group's UMI, or with two groups to the other group's UMI
        let pairs: &[(usize, usize)] = match groups.len() {
            1 => &[(0, 0)],
            2 => &[(0, 1), (1, 0)],
            _ => return,
        };
        let read_count: usize = groups.iter().map(|g| g.total_read_count()).sum();
        let counts = self.umi_group_counts(groups.len() as i32, read_count);
        counts.group_count += 1;
        for &(test, reads) in pairs {
            let test_umi = groups[test].umi.as_deref().unwrap_or("").as_bytes();
            for read in groups[reads].reads.iter().chain(groups[reads].non_consensus_reads.iter()) {
                let diff = calc_umi_id_diff(config.extract_umi_id(read.qname()), test_umi);
                if diff <= MAX_EDIT_DISTANCE {
                    counts.edit_distance_frequency[diff] += 1;
                }
            }
        }
    }

    fn record_umi_base_frequencies(&mut self, umi_id: &[u8]) {
        if self.base_frequency_length.is_none() {
            self.base_frequency_length = Some(umi_id.len());
        }
        if self.base_frequencies.len() < umi_id.len() {
            self.base_frequencies.resize(umi_id.len(), [0; 5]);
        }
        for (p, base) in umi_id.iter().enumerate() {
            let index = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'N' => 4,
                _ => continue,
            };
            self.base_frequencies[p][index] += 1;
        }
    }

    /// Adds the statistics of the next partition in genome order, with the result one statistics object would
    /// have after processing both: that is what Java has with one thread. (UmiStatistics.merge, which Java uses
    /// to combine its threads, keeps the first thread's maximum group for a row instead of the largest; which
    /// partitions a Java thread processes depends on timing.)
    pub fn merge(&mut self, other: &UmiStatistics) {
        self.umi_groups += other.umi_groups;
        self.merge_counts.poly_g_merged += other.merge_counts.poly_g_merged;
        self.merge_counts.poly_g_unmapped_groups += other.merge_counts.poly_g_unmapped_groups;
        self.merge_counts.duplex_collapsed += other.merge_counts.duplex_collapsed;
        self.merge_counts.jitter_merged += other.merge_counts.jitter_merged;
        for o in &other.umi_group_frequencies {
            match self
                .umi_group_frequencies
                .iter_mut()
                .find(|x| x.duplicate_group_count == o.duplicate_group_count && x.read_count == o.read_count)
            {
                Some(m) => {
                    m.group_count += o.group_count;
                    for i in 0..=MAX_EDIT_DISTANCE {
                        m.edit_distance_frequency[i] += o.edit_distance_frequency[i];
                    }
                }
                None => self.umi_group_frequencies.push(o.clone()),
            }
        }
        for o in &other.position_fragments {
            let data = self.position_fragment_data(o.unique_coord_count, o.unique_fragment_count);
            data.frequency += o.frequency;
            data.max_coord_umi_count = data.max_coord_umi_count.max(o.max_coord_umi_count);
            if o.max_umi_reads_count > data.max_umi_reads_count {
                data.max_umi_reads_count = o.max_umi_reads_count;
                data.umi_group_details = o.umi_group_details.clone();
            }
        }
        if self.base_frequency_length.is_none() {
            self.base_frequency_length = other.base_frequency_length;
        }
        if self.base_frequencies.len() < other.base_frequencies.len() {
            self.base_frequencies.resize(other.base_frequencies.len(), [0; 5]);
        }
        for (p, counts) in other.base_frequencies.iter().enumerate() {
            for b in 0..5 {
                self.base_frequencies[p][b] += counts[b];
            }
        }
    }

    /// UmiStatistics.writePositionFragmentsData
    pub fn write_position_fragments_data(&self, filename: &str) -> Result<()> {
        let mut out = String::from("UniqueCoordsWithStartPos\tUniquePrimariesWithStartPos\tFrequency\tMaxCoordUmiCount\tMaxUmiReads\tMaxUmiDetails\n");
        for d in &self.position_fragments {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                d.unique_coord_count, d.unique_fragment_count, d.frequency, d.max_coord_umi_count, d.max_umi_reads_count, d.umi_group_details
            ));
        }
        std::fs::write(filename, out)?;
        Ok(())
    }

    /// UmiStatistics.writeUmiBaseDiffStats
    pub fn write_umi_base_diff_stats(&self, filename: &str) -> Result<()> {
        let mut out = String::from("UmiCountWithCoords\tDuplicateReadCount\tFrequency");
        for i in 0..=MAX_EDIT_DISTANCE {
            out.push_str(&format!("\tED_{i}"));
        }
        out.push('\n');
        for c in &self.umi_group_frequencies {
            out.push_str(&format!("{}\t{}\t{}", c.duplicate_group_count, c.read_count, c.group_count));
            for v in c.edit_distance_frequency {
                out.push_str(&format!("\t{v}"));
            }
            out.push('\n');
        }
        std::fs::write(filename, out)?;
        Ok(())
    }

    /// UmiStatistics.writeUmiBaseFrequencyStats: one row per position of the first UMI recorded. (Without any UMI
    /// recorded Java fails here on its unset table; this writes the header alone.)
    pub fn write_umi_base_frequency_stats(&self, filename: &str) -> Result<()> {
        let mut out = String::from("UmiPosition\tACount\tCCount\tGCount\tTCount\tNCount\n");
        for p in 0..self.base_frequency_length.unwrap_or(0) {
            out.push_str(&format!("{}", p + 1));
            for b in 0..5 {
                out.push_str(&format!("\t{}", self.base_frequencies[p][b]));
            }
            out.push('\n');
        }
        std::fs::write(filename, out)?;
        Ok(())
    }
}

/// A duplicate group or a single read: an element of the Java lists of Object
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Frag {
    Group(usize),
    Single(usize),
}

/// ReadInfo with its lazily extracted UMI, held as the offset in the read name at which it starts
struct Single {
    info: Option<ReadInfo>,
    umi_start: Option<usize>,
}

/// The groups and single reads of one call. Entries stay where they are; the lists of the algorithm hold indices.
struct Arena<'c> {
    config: &'c UmiConfig,
    groups: Vec<Option<DuplicateGroup>>,
    singles: Vec<Single>,
    counts: UmiMergeCounts,
}

fn list_remove(list: &mut Vec<usize>, id: usize) {
    if let Some(p) = list.iter().position(|&x| x == id) {
        list.remove(p);
    }
}

fn new_group(umi: Option<String>, reads: Vec<bam::Record>, coords: FragmentCoords) -> DuplicateGroup {
    let mut g = DuplicateGroup::new(reads, coords);
    g.umi = umi;
    g
}

impl<'c> Arena<'c> {
    fn new(config: &'c UmiConfig) -> Self {
        Arena { config, groups: Vec::new(), singles: Vec::new(), counts: UmiMergeCounts::default() }
    }
    fn add_group(&mut self, group: DuplicateGroup) -> usize {
        self.groups.push(Some(group));
        self.groups.len() - 1
    }
    fn add_single(&mut self, info: ReadInfo) -> usize {
        self.singles.push(Single { info: Some(info), umi_start: None });
        self.singles.len() - 1
    }
    fn group(&self, i: usize) -> &DuplicateGroup {
        self.groups[i].as_ref().unwrap()
    }
    fn group_mut(&mut self, i: usize) -> &mut DuplicateGroup {
        self.groups[i].as_mut().unwrap()
    }
    fn single_info(&self, i: usize) -> &ReadInfo {
        self.singles[i].info.as_ref().unwrap()
    }
    fn coords(&self, f: Frag) -> &FragmentCoords {
        match f {
            Frag::Group(g) => &self.group(g).coords,
            Frag::Single(s) => &self.single_info(s).coords,
        }
    }
    /// ReadInfo.getOrExtractUmi
    fn extract_single_umi(&mut self, i: usize) {
        if self.singles[i].umi_start.is_none() {
            let name = self.singles[i].info.as_ref().unwrap().read.qname();
            let start = name.len() - self.config.extract_umi_id(name).len();
            self.singles[i].umi_start = Some(start);
        }
    }
    /// the UMI of a group, or of a single read (extracting it if need be)
    fn extract_umi(&mut self, f: Frag) {
        if let Frag::Single(s) = f {
            self.extract_single_umi(s);
        }
    }
    /// the UMI after `extract_umi`
    fn umi(&self, f: Frag) -> &str {
        match f {
            Frag::Group(g) => self.group(g).umi.as_deref().unwrap_or(""),
            Frag::Single(s) => match (&self.singles[s].info, self.singles[s].umi_start) {
                (Some(info), Some(start)) => std::str::from_utf8(&info.read.qname()[start..]).unwrap_or(""),
                _ => "",
            },
        }
    }
    /// new DuplicateGroup(umi, readInfo.read(), readInfo.fragCoordinates())
    fn single_to_group(&mut self, s: usize) -> usize {
        let info = self.singles[s].info.take().unwrap();
        let umi = self.singles[s].umi_start.map(|start| String::from_utf8_lossy(&info.read.qname()[start..]).into_owned());
        self.add_group(new_group(umi, vec![info.read], info.coords))
    }
    fn take_single_read(&mut self, s: usize) -> bam::Record {
        self.singles[s].info.take().unwrap().read
    }
}

/// UmiGroupBuilder.buildUmiGroups: the reads of one set of fragment coordinates grouped by UMI, then smaller
/// groups merged into larger ones with similar UMIs
pub fn build_umi_groups(coords: &FragmentCoords, reads: Vec<bam::Record>, config: &UmiConfig) -> Vec<DuplicateGroup> {
    let mut check_defined_umis = config.has_defined_umis();
    let mut use_defined_umis = check_defined_umis;
    let mut groups: Vec<DuplicateGroup> = Vec::new();
    // group index by UMI, used once there are too many groups for a scan
    let mut index: HashMap<Vec<u8>, usize> = HashMap::new();
    const SCAN_LIMIT: usize = 12;

    for read in reads {
        let mut umi_id: &[u8] = config.extract_umi_id(read.qname());
        if check_defined_umis {
            match config.match_defined_umi_id(umi_id) {
                // once a read has no defined UMI the remaining reads keep their own
                None => {
                    use_defined_umis = false;
                    check_defined_umis = false;
                }
                Some(defined) => umi_id = defined.as_bytes(),
            }
        }
        let found = if groups.len() <= SCAN_LIMIT {
            groups.iter().position(|g| g.umi.as_deref().unwrap().as_bytes() == umi_id)
        } else {
            if index.is_empty() {
                for (i, g) in groups.iter().enumerate() {
                    index.insert(g.umi.as_deref().unwrap().as_bytes().to_vec(), i);
                }
            }
            index.get(umi_id).copied()
        };
        match found {
            Some(i) => groups[i].reads.push(read),
            None => {
                let umi = String::from_utf8_lossy(umi_id).into_owned();
                if groups.len() > SCAN_LIMIT {
                    index.insert(umi_id.to_vec(), groups.len());
                }
                groups.push(new_group(Some(umi), vec![read], coords.clone()));
            }
        }
    }

    if use_defined_umis {
        // the values of the Java HashMap keyed by UMI, in its iteration order
        let mut map: JavaMap<usize, ()> = JavaMap::new();
        for (i, g) in groups.iter().enumerate() {
            map.insert(i, java_spread_hash(g.umi.as_deref().unwrap().as_bytes()), ());
        }
        let mut slots: Vec<Option<DuplicateGroup>> = groups.into_iter().map(Some).collect();
        return map.order().into_iter().map(|i| slots[i].take().unwrap()).collect();
    }

    // by descending number of fragments, then UMI
    let mut ordered = groups;
    if ordered.len() > 1 {
        ordered.sort_by(|a, b| b.reads.len().cmp(&a.reads.len()).then_with(|| a.umi.cmp(&b.umi)));
    }
    fn umi(g: &DuplicateGroup) -> &[u8] {
        g.umi.as_deref().unwrap().as_bytes()
    }

    // the directional model: a group joins a larger (or equal) one within the permitted difference, and can then
    // draw in further groups itself
    let mut i = 0;
    while i + 1 < ordered.len() {
        let mut cluster: Vec<DuplicateGroup> = Vec::new();
        let mut j = i + 1;
        while j < ordered.len() {
            let second = &ordered[j];
            let merged = std::iter::once(&ordered[i])
                .chain(cluster.iter())
                .any(|e| e.reads.len() >= second.reads.len() && !exceeds_umi_id_diff(umi(e), umi(second), config.permitted_base_diff));
            if !merged {
                j += 1;
            } else {
                let second = ordered.remove(j);
                cluster.push(second);
                // a group just added may be close enough to one skipped earlier
                j = i + 1;
            }
        }
        for c in cluster {
            ordered[i].reads.extend(c.reads);
        }
        i += 1;
    }

    // groups whose UMIs differ by one more base than permitted
    if ordered.len() > 1 {
        merge_groups_if(&mut ordered, |first, second| !exceeds_umi_id_diff(umi(first), umi(second), config.permitted_base_diff + 1));
    }

    // up to 4 bases where one group is at least 25 times the size of the other
    let has_large_groups = ordered.iter().any(|g| g.reads.len() >= MAX_IMBALANCED_UMI_COUNT);
    if ordered.len() > 1 && has_large_groups {
        merge_groups_if(&mut ordered, |first, second| {
            let (a, b) = (first.reads.len() as f64, second.reads.len() as f64);
            let max_count_ratio = if a >= b { a / b } else { b / a };
            max_count_ratio >= MAX_IMBALANCED_UMI_COUNT as f64 && !exceeds_umi_id_diff(umi(first), umi(second), MAX_IMBALANCED_UMI_BASE_DIFF)
        });
    }
    ordered
}

/// each group takes the reads of the later groups that match it
fn merge_groups_if(ordered: &mut Vec<DuplicateGroup>, matches: impl Fn(&DuplicateGroup, &DuplicateGroup) -> bool) {
    let mut i = 0;
    while i < ordered.len() {
        let mut j = i + 1;
        while j < ordered.len() {
            if matches(&ordered[i], &ordered[j]) {
                let second = ordered.remove(j);
                ordered[i].reads.extend(second.reads);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

/// UmiGroupBuilder.CoordinateGroup: the groups and single reads at one set of coordinates, by fragment orientation
struct CoordinateGroup {
    forward: Option<Vec<Frag>>,
    reverse: Option<Vec<Frag>>,
}

/// The Java HashMap of coordinate groups by non-oriented coordinate key. The key of an entry is that of its first
/// fragment, so no key is copied; entries are found by the key's Java hash.
struct CoordinateGroupMap {
    /// (hash, first fragment, group) in insertion order
    entries: Vec<(u32, Frag, CoordinateGroup)>,
    /// hash to the latest entry with it, and from each entry to the one before with the same hash
    by_hash: HashMap<u32, usize>,
    previous_same_hash: Vec<usize>,
}

impl CoordinateGroupMap {
    const SCAN_LIMIT: usize = 16;
    const NONE: usize = usize::MAX;

    fn new() -> Self {
        CoordinateGroupMap { entries: Vec::new(), by_hash: HashMap::new(), previous_same_hash: Vec::new() }
    }

    fn get_or_create(&mut self, arena: &Arena, frag: Frag) -> &mut CoordinateGroup {
        let key = arena.coords(frag).key_non_oriented();
        let hash = java_spread_hash(key.as_bytes());
        let same_key = |entry: &(u32, Frag, CoordinateGroup)| entry.0 == hash && arena.coords(entry.1).key_non_oriented() == key;
        let found = if self.entries.len() <= Self::SCAN_LIMIT {
            self.entries.iter().position(same_key)
        } else {
            if self.by_hash.is_empty() {
                self.previous_same_hash.clear();
                for (i, e) in self.entries.iter().enumerate() {
                    self.previous_same_hash.push(self.by_hash.insert(e.0, i).unwrap_or(Self::NONE));
                }
            }
            let mut i = self.by_hash.get(&hash).copied().unwrap_or(Self::NONE);
            while i != Self::NONE && !same_key(&self.entries[i]) {
                i = self.previous_same_hash[i];
            }
            (i != Self::NONE).then_some(i)
        };
        let i = match found {
            Some(i) => i,
            None => {
                let i = self.entries.len();
                if i > Self::SCAN_LIMIT {
                    self.previous_same_hash.push(self.by_hash.insert(hash, i).unwrap_or(Self::NONE));
                }
                self.entries.push((hash, frag, CoordinateGroup { forward: None, reverse: None }));
                i
            }
        };
        &mut self.entries[i].2
    }

    /// entry indices in java.util.HashMap iteration order: by bucket of the table as grown for this many
    /// entries, then by insertion
    fn order(&self) -> Vec<usize> {
        let mut capacity: u32 = 16;
        while self.entries.len() as u32 > capacity / 4 * 3 {
            capacity *= 2;
        }
        let mask = capacity - 1;
        let mut v: Vec<usize> = (0..self.entries.len()).collect();
        if v.len() > 1 {
            v.sort_by_key(|&i| (self.entries[i].0 & mask, i));
        }
        v
    }
}

/// UmiGroupBuilder.processUmiGroups. `mergers` is true for Illumina duplex UMIs, where the poly-G and jitter
/// mergers run. Returns the UMI groups and the single fragments.
pub fn process_umi_groups(
    config: &UmiConfig, is_illumina: bool, duplicate_groups: Vec<DuplicateGroup>, single_reads: Vec<ReadInfo>, capture_stats: bool,
    stats: &mut UmiStatistics,
) -> Result<(Vec<DuplicateGroup>, Vec<ReadInfo>)> {
    let mergers = is_illumina && config.duplex;
    let duplicate_group_count = duplicate_groups.len();
    let form_coord_groups = config.base_stats || duplicate_group_count + single_reads.len() > 1;

    let mut arena = Arena::new(config);
    let mut all_umi_groups: Vec<usize> = Vec::new();
    for group in duplicate_groups {
        for umi_group in build_umi_groups(&group.coords, group.reads, config) {
            all_umi_groups.push(arena.add_group(umi_group));
        }
    }
    let mut single_fragments: Vec<usize> = single_reads.into_iter().map(|info| arena.add_single(info)).collect();

    if mergers {
        poly_g_merge_groups(&mut arena, &mut all_umi_groups, &mut single_fragments);
    }

    let mut coordinate_groups = CoordinateGroupMap::new();
    let mut coordinate_order: Vec<usize> = Vec::new();
    if form_coord_groups {
        // by descending fragment count (a stable sort)
        all_umi_groups.sort_by(|&a, &b| arena.group(b).total_read_count().cmp(&arena.group(a).total_read_count()));
        for frag in all_umi_groups.drain(..).map(Frag::Group).chain(single_fragments.iter().map(|&s| Frag::Single(s))) {
            let forward = arena.coords(frag).frag_forward;
            let coord_group = coordinate_groups.get_or_create(&arena, frag);
            let list = if forward { &mut coord_group.forward } else { &mut coord_group.reverse };
            list.get_or_insert_with(Vec::new).push(frag);
        }
        coordinate_order = coordinate_groups.order();

        if config.base_stats {
            for &i in &coordinate_order {
                let coord_group = &coordinate_groups.entries[i].2;
                capture_umi_group_stats(&arena, coord_group.forward.as_deref(), stats)?;
                capture_umi_group_stats(&arena, coord_group.reverse.as_deref(), stats)?;
            }
        }

        // collapse duplex UMIs with opposite orientations
        for &i in &coordinate_order {
            collapse_coordinate_group(&mut arena, &mut coordinate_groups.entries[i].2, &mut all_umi_groups, &mut single_fragments)?;
        }
    }

    if mergers {
        jitter_merge_groups(&mut arena, &mut all_umi_groups, &mut single_fragments)?;
    }

    // groups of one read go back to the single fragments
    let mut final_umi_groups: Vec<usize> = Vec::with_capacity(all_umi_groups.len());
    for g in all_umi_groups {
        if arena.group(g).total_read_count() == 1 {
            let mut group = arena.groups[g].take().unwrap();
            let read = group.reads.pop().unwrap();
            single_fragments.push(arena.add_single(ReadInfo { read, coords: group.coords }));
        } else {
            final_umi_groups.push(g);
        }
    }

    if capture_stats {
        // UmiGroupBuilder.captureStats
        let mut unique_coord_count = 0;
        let mut unique_fragment_count = 0;
        let mut max_coord_umi_count = 0;
        if form_coord_groups {
            for &i in &coordinate_order {
                let coord_group = &coordinate_groups.entries[i].2;
                if let Some(f) = &coord_group.forward {
                    unique_fragment_count += f.len();
                    unique_coord_count += 1;
                    max_coord_umi_count = max_coord_umi_count.max(f.len());
                }
                if let Some(r) = &coord_group.reverse {
                    unique_fragment_count += r.len();
                    if !r.is_empty() {
                        unique_coord_count += 1;
                    }
                    max_coord_umi_count = max_coord_umi_count.max(r.len());
                }
            }
        } else {
            unique_fragment_count = single_fragments.len() + final_umi_groups.len();
            unique_coord_count = duplicate_group_count + single_fragments.len();
        }
        let mut max_umi_fragment_count = 0;
        let mut max_umi_group = None;
        for &g in &final_umi_groups {
            stats.umi_groups += 1;
            let count = arena.group(g).total_read_count();
            if count > max_umi_fragment_count {
                max_umi_group = Some(g);
                max_umi_fragment_count = count;
            }
        }
        if unique_coord_count > 0 && unique_fragment_count > 0 {
            stats.record_fragment_positions(
                unique_coord_count as i32, unique_fragment_count as i32, max_coord_umi_count as i32, max_umi_group.map(|g| arena.group(g)),
            );
        }
    }

    let counts = arena.counts;
    stats.merge_counts.poly_g_merged += counts.poly_g_merged;
    stats.merge_counts.poly_g_unmapped_groups += counts.poly_g_unmapped_groups;
    stats.merge_counts.duplex_collapsed += counts.duplex_collapsed;
    stats.merge_counts.jitter_merged += counts.jitter_merged;
    let groups = final_umi_groups.into_iter().map(|g| arena.groups[g].take().unwrap()).collect();
    let singles = single_fragments.into_iter().map(|s| arena.singles[s].info.take().unwrap()).collect();
    Ok((groups, singles))
}

/// UmiGroupBuilder.captureUmiGroupStats. Java casts every element that is not a group to a read, which fails for
/// the single fragments held there (they are ReadInfo objects) and ends the run; the same is an error here.
fn capture_umi_group_stats(arena: &Arena, frag_groups: Option<&[Frag]>, stats: &mut UmiStatistics) -> Result<()> {
    let Some(frag_groups) = frag_groups else { return Ok(()) };
    if !crate::compat::java("R24") {
        // R24 corrected: a fragment without duplicates counts as a group of its one read, as the Java code intends
        // (`new DuplicateGroup(umi, read, null)`)
        let singles: Vec<Option<DuplicateGroup>> = frag_groups
            .iter()
            .map(|f| match f {
                Frag::Group(_) => None,
                Frag::Single(s) => arena.singles[*s].info.as_ref().map(|info| {
                    let umi = String::from_utf8_lossy(arena.config.extract_umi_id(info.read.qname())).into_owned();
                    new_group(Some(umi), vec![info.read.clone()], info.coords.clone())
                }),
            })
            .collect();
        let groups: Vec<&DuplicateGroup> = frag_groups
            .iter()
            .zip(&singles)
            .filter_map(|(f, single)| match f {
                Frag::Group(g) => Some(arena.group(*g)),
                Frag::Single(_) => single.as_ref(),
            })
            .collect();
        stats.record_umi_base_stats(arena.config, &groups);
        return Ok(());
    }
    let mut groups: Vec<&DuplicateGroup> = Vec::with_capacity(frag_groups.len());
    for f in frag_groups {
        match f {
            Frag::Group(g) => groups.push(arena.group(*g)),
            Frag::Single(_) => bail!(
                "umi_base_diff_stats: a fragment without duplicates cannot be used for the UMI base statistics \
                 (Java REDUX 2.1 fails here with a ClassCastException in UmiGroupBuilder.captureUmiGroupStats)"
            ),
        }
    }
    stats.record_umi_base_stats(arena.config, &groups);
    Ok(())
}

/// UmiGroupBuilder.collapseCoordinateGroup: fragments with the same coordinates and opposite orientations
/// (F1R2 and F2R1) are merged, one to one, when their duplex UMIs are complementary
fn collapse_coordinate_group(
    arena: &mut Arena, coord_group: &mut CoordinateGroup, all_umi_groups: &mut Vec<usize>, single_fragments: &mut Vec<usize>,
) -> Result<()> {
    let config = arena.config;
    let groups_of = |list: &Option<Vec<Frag>>, out: &mut Vec<usize>| {
        for f in list.iter().flatten() {
            if let Frag::Group(g) = f {
                out.push(*g);
            }
        }
    };
    if !config.duplex || coord_group.forward.is_none() || coord_group.reverse.is_none() {
        groups_of(&coord_group.forward, all_umi_groups);
        groups_of(&coord_group.reverse, all_umi_groups);
        return Ok(());
    }
    let forward = coord_group.forward.as_ref().unwrap();
    let reverse = coord_group.reverse.as_mut().unwrap();
    for &first in forward {
        arena.extract_umi(first);
        let mut first_group = match first {
            Frag::Group(g) => Some(g),
            Frag::Single(_) => None,
        };
        let mut second_index = 0;
        while second_index < reverse.len() {
            let second = reverse[second_index];
            arena.extract_umi(second);
            if !has_duplex_umi_match(arena.umi(first), arena.umi(second), &config.duplex_delim, config.permitted_base_diff)? {
                second_index += 1;
                continue;
            }
            reverse.remove(second_index);
            let target = match (first_group, first) {
                (Some(g), _) => g,
                (None, Frag::Single(s)) => {
                    let g = arena.single_to_group(s);
                    list_remove(single_fragments, s);
                    first_group = Some(g);
                    g
                }
                (None, Frag::Group(g)) => g,
            };
            match second {
                Frag::Group(g2) => {
                    let reads = std::mem::take(&mut arena.group_mut(g2).reads);
                    let non_consensus = std::mem::take(&mut arena.group_mut(g2).non_consensus_reads);
                    let t = arena.group_mut(target);
                    t.reads.extend(reads);
                    t.non_consensus_reads.extend(non_consensus);
                }
                Frag::Single(s2) => {
                    let read = arena.take_single_read(s2);
                    arena.group_mut(target).reads.push(read);
                    list_remove(single_fragments, s2);
                }
            }
            arena.group_mut(target).dual_strand = true;
            arena.counts.duplex_collapsed += 1;
            // one to one only
            break;
        }
        if let Some(g) = first_group {
            all_umi_groups.push(g);
        }
    }
    for f in reverse.iter() {
        if let Frag::Group(g) = f {
            all_umi_groups.push(*g);
        }
    }
    Ok(())
}

fn poly_g_tail_length(umi_id: &[u8]) -> usize {
    umi_id.iter().rev().take_while(|&&c| c == b'G').count()
}

/// PolyGUmiMerger.umiMatch: equal over the length of the shorter
fn umi_prefix_match(umi1: &[u8], umi2: &[u8]) -> bool {
    umi1.iter().zip(umi2).all(|(a, b)| a == b)
}

/// PolyGUmiMerger.hasUnmappedRead
fn has_unmapped_read(coords: &FragmentCoords) -> bool {
    !coords.unpaired && coords.supp.is_none() && coords.pos_upper == NO_POSITION
}

/// PolyGUmiMerger.PolyGMergeGroup
struct PolyGMergeGroup {
    /// the UMIs with their poly-G tails trimmed
    umis: Vec<String>,
    umi_has_poly_g_tail: bool,
    pos_lower: i32,
    lower_forward: bool,
    unmapped_sourced: bool,
    has_unmapped_read: bool,
    fragments: Vec<Frag>,
    merged: bool,
}

/// PolyGUmiMerger.mergeGroups (Illumina duplex UMIs): fragments at the same lower position whose UMIs differ only
/// by a tail that has degraded to Gs on the unmapped read's side are merged. Fragments with both reads mapped take
/// precedence for the consensus; groups of unmapped reads that match one are marked as duplicates as a whole.
fn poly_g_merge_groups(arena: &mut Arena, umi_groups: &mut Vec<usize>, single_fragments: &mut Vec<usize>) {
    if umi_groups.is_empty() && single_fragments.is_empty() {
        return;
    }
    // PolyGUmiMerger.hasCandidateFragments: an unmapped read and a poly-G tail, anywhere
    let mut has_unmapped = false;
    let mut has_poly_g = false;
    for &g in umi_groups.iter() {
        has_unmapped |= has_unmapped_read(&arena.group(g).coords);
        has_poly_g |= poly_g_tail_length(arena.umi(Frag::Group(g)).as_bytes()) >= MIN_POLYG_UMI_TAIL_LENGTH;
        if has_unmapped && has_poly_g {
            break;
        }
    }
    for &s in single_fragments.iter() {
        has_unmapped |= has_unmapped_read(&arena.single_info(s).coords);
        arena.extract_single_umi(s);
        has_poly_g |= poly_g_tail_length(arena.umi(Frag::Single(s)).as_bytes()) >= MIN_POLYG_UMI_TAIL_LENGTH;
        if has_unmapped && has_poly_g {
            break;
        }
    }
    if !(has_unmapped && has_poly_g) {
        return;
    }

    // PolyGUmiMerger.findCandidateFragments: by lower position, orientation, unmapped and supplementary state
    // (Java key "%d_%c_%c_%s")
    let mut group_map: JavaMap<(i32, bool, bool, bool), Vec<PolyGMergeGroup>> = JavaMap::new();
    let frags: Vec<Frag> = umi_groups.iter().map(|&g| Frag::Group(g)).chain(single_fragments.iter().map(|&s| Frag::Single(s))).collect();
    for frag in frags {
        arena.extract_umi(frag);
        let coords = arena.coords(frag);
        let key = (coords.pos_lower, coords.lower_forward, coords.unmapped_sourced, coords.supp.is_some());
        let umi = arena.umi(frag).as_bytes();
        let tail = poly_g_tail_length(umi);
        let group = PolyGMergeGroup {
            umis: vec![String::from_utf8_lossy(&umi[..umi.len() - tail]).into_owned()],
            umi_has_poly_g_tail: tail >= MIN_POLYG_UMI_TAIL_LENGTH,
            pos_lower: coords.pos_lower,
            lower_forward: coords.lower_forward,
            unmapped_sourced: coords.unmapped_sourced,
            has_unmapped_read: has_unmapped_read(coords),
            fragments: vec![frag],
            merged: false,
        };
        match group_map.find(&key) {
            Some(i) => group_map.value_mut(i).push(group),
            None => {
                let hash = KeyBuf::new()
                    .int(key.0)
                    .push(b'_')
                    .push(if key.1 { b'F' } else { b'R' })
                    .push(b'_')
                    .push(if key.2 { b'N' } else { b'M' })
                    .push(b'_')
                    .push(if key.3 { b'S' } else { b'P' })
                    .hash();
                group_map.insert(key, hash, vec![group]);
            }
        }
    }
    for i in group_map.order() {
        let groups = group_map.value(i);
        if !groups.iter().any(|g| g.umi_has_poly_g_tail) && !groups.iter().any(|g| g.has_unmapped_read) {
            group_map.remove(i);
        }
    }
    if group_map.is_empty() {
        return;
    }
    let order = group_map.order();

    let mut has_merges = false;
    for &m in &order {
        let groups = group_map.value_mut(m);
        if groups.len() == 1 {
            continue;
        }
        for i in 0..groups.len() - 1 {
            if groups[i].merged {
                continue;
            }
            for j in i + 1..groups.len() {
                // PolyGMergeGroup.canMerge: the coordinates match by the key
                let (head, tail) = groups.split_at_mut(j);
                let (group1, group2) = (&mut head[i], &mut tail[0]);
                if group2.merged {
                    continue;
                }
                let can_merge = group1.umis.iter().any(|u| group2.umis.iter().any(|x| umi_prefix_match(u.as_bytes(), x.as_bytes())));
                if can_merge {
                    group1.umis.extend(group2.umis.iter().cloned());
                    group1.fragments.extend(group2.fragments.iter().copied());
                    group2.merged = true;
                    has_merges = true;
                    arena.counts.poly_g_merged += 1;
                }
            }
        }
    }
    if !has_merges {
        return;
    }

    // remove the merged fragments and move their reads to the group they were merged into
    for &m in &order {
        let groups = group_map.value_mut(m);
        if groups.len() == 1 {
            continue;
        }
        for group in groups.iter_mut() {
            if group.merged {
                for f in &group.fragments {
                    match f {
                        Frag::Group(g) => list_remove(umi_groups, *g),
                        Frag::Single(s) => list_remove(single_fragments, *s),
                    }
                }
            } else if group.fragments.len() > 1 {
                let mut original: Option<usize> = None;
                for i in 0..group.fragments.len() {
                    let fragment = group.fragments[i];
                    let Some(target) = original else {
                        original = Some(match fragment {
                            Frag::Single(s) => {
                                let g = arena.single_to_group(s);
                                umi_groups.push(g);
                                group.fragments[0] = Frag::Group(g);
                                list_remove(single_fragments, s);
                                g
                            }
                            Frag::Group(g) => g,
                        });
                        continue;
                    };
                    let reads_to_add = match fragment {
                        Frag::Single(s) => vec![arena.take_single_read(s)],
                        Frag::Group(g) => std::mem::take(&mut arena.group_mut(g).reads),
                    };
                    let t = arena.group_mut(target);
                    if t.coords.unmapped_sourced {
                        t.non_consensus_reads.extend(reads_to_add);
                    } else {
                        t.reads.extend(reads_to_add);
                    }
                }
            }
        }
    }

    // groups of mapped reads with a mapped mate, by trimmed UMI, lower position and orientation
    let mut groups_with_mapped_mates: Vec<(String, i32, bool)> = Vec::new();
    for &m in &order {
        for group in group_map.value(m) {
            if group.merged || group.unmapped_sourced {
                continue;
            }
            if let Frag::Group(g) = group.fragments[0] {
                let duplicate_group = arena.group(g);
                if duplicate_group.reads.iter().any(|r| !has(r, FLAG_MATE_UNMAPPED)) {
                    groups_with_mapped_mates.push((group.umis[0].clone(), duplicate_group.coords.pos_lower, duplicate_group.coords.lower_forward));
                }
            }
        }
    }
    // the unmapped reads that match one of them are duplicates only
    for &m in &order {
        for group in group_map.value(m) {
            if group.merged || !group.unmapped_sourced {
                continue;
            }
            let matches = groups_with_mapped_mates.iter().any(|(trimmed, pos, forward)| {
                group.pos_lower == *pos && group.lower_forward == *forward && umi_prefix_match(trimmed.as_bytes(), group.umis[0].as_bytes())
            });
            if !matches {
                continue;
            }
            match group.fragments[0] {
                Frag::Group(g) => {
                    let duplicate_group = arena.group_mut(g);
                    if !duplicate_group.reads.iter().any(|r| has(r, FLAG_MATE_UNMAPPED)) {
                        continue;
                    }
                    duplicate_group.poly_g_unmapped = true;
                    arena.counts.poly_g_unmapped_groups += 1;
                }
                Frag::Single(s) => {
                    if !has(&arena.single_info(s).read, FLAG_UNMAPPED) {
                        continue;
                    }
                    let g = arena.single_to_group(s);
                    arena.group_mut(g).poly_g_unmapped = true;
                    arena.counts.poly_g_unmapped_groups += 1;
                    umi_groups.push(g);
                    list_remove(single_fragments, s);
                }
            }
        }
    }
}

/// JitterUmiMerger.JitterMergeGroup
struct JitterMergeGroup {
    /// the fragment the group was made from, for its UMI
    frag: Frag,
    read_is_lower: bool,
    frag_forward: bool,
    lower_forward: bool,
    upper_forward: bool,
    pos_lower: i32,
    pos_upper: i32,
    fragments: Vec<Frag>,
    merged: bool,
}

impl JitterMergeGroup {
    /// JitterMergeGroup.canMerge: the same UMI (or its duplex complement for the other orientation), without any
    /// base difference, and the other end of the fragment within the jitter distance
    fn can_merge(&self, other: &JitterMergeGroup, check_lower: bool, arena: &Arena) -> Result<bool> {
        let config = arena.config;
        let (umi, other_umi) = (arena.umi(self.frag), arena.umi(other.frag));
        if other.merged {
            return Ok(false);
        }
        // the lower and the upper read of the same fragments are kept apart
        if self.read_is_lower != other.read_is_lower {
            return Ok(false);
        }
        if self.frag_forward == other.frag_forward {
            if exceeds_umi_id_diff(umi.as_bytes(), other_umi.as_bytes(), 0) {
                return Ok(false);
            }
        } else if !has_duplex_umi_match(umi, other_umi, &config.duplex_delim, 0)? {
            return Ok(false);
        }
        Ok(if check_lower {
            self.lower_forward == other.lower_forward && (self.pos_lower - other.pos_lower).abs() <= SINGLE_END_JITTER_COLLAPSE_DISTANCE
        } else {
            self.upper_forward == other.upper_forward && (self.pos_upper - other.pos_upper).abs() <= SINGLE_END_JITTER_COLLAPSE_DISTANCE
        })
    }
}

/// JitterUmiMerger.mergeGroups (Illumina duplex UMIs): fragments that share one end exactly and whose other end
/// is within 10 bases are merged when their UMIs match exactly. The reads merged in are not used for the consensus.
fn jitter_merge_groups(arena: &mut Arena, umi_groups: &mut Vec<usize>, single_fragments: &mut Vec<usize>) -> Result<()> {
    if umi_groups.is_empty() && single_fragments.is_empty() {
        return Ok(());
    }
    // every fragment is a candidate group, listed by its lower and by its upper position (Java key "%d_%c_%c")
    let mut groups: Vec<JitterMergeGroup> = Vec::with_capacity(umi_groups.len() + single_fragments.len());
    let mut pos_groups: [JavaMap<(i32, bool, bool), Vec<usize>>; 2] = [JavaMap::new(), JavaMap::new()];
    let mut has_candidates = false;
    let frags: Vec<Frag> = umi_groups.iter().map(|&g| Frag::Group(g)).chain(single_fragments.iter().map(|&s| Frag::Single(s))).collect();
    for frag in frags {
        arena.extract_umi(frag);
        let coords = arena.coords(frag);
        let index = groups.len();
        let supp = coords.supp.is_some();
        let keys = [(coords.pos_lower, coords.lower_forward, supp), (coords.pos_upper, coords.upper_forward, supp)];
        groups.push(JitterMergeGroup {
            frag,
            read_is_lower: coords.read_is_lower,
            frag_forward: coords.frag_forward,
            lower_forward: coords.lower_forward,
            upper_forward: coords.upper_forward,
            pos_lower: coords.pos_lower,
            pos_upper: coords.pos_upper,
            fragments: vec![frag],
            merged: false,
        });
        for (g, key) in keys.into_iter().enumerate() {
            match pos_groups[g].find(&key) {
                Some(i) => {
                    has_candidates = true;
                    pos_groups[g].value_mut(i).push(index);
                }
                None => {
                    let hash =
                        KeyBuf::new().int(key.0).push(b'_').push(if key.1 { b'F' } else { b'R' }).push(b'_').push(if key.2 { b'S' } else { b'P' }).hash();
                    pos_groups[g].insert(key, hash, vec![index]);
                }
            }
        }
    }
    if !has_candidates {
        return Ok(());
    }

    let mut has_merges = false;
    for (g, map) in pos_groups.iter().enumerate() {
        let use_lower = g == 0;
        for m in map.order() {
            let list = map.value(m);
            if list.len() == 1 {
                continue;
            }
            for i in 0..list.len() - 1 {
                let i1 = list[i];
                if groups[i1].merged {
                    continue;
                }
                for &i2 in &list[i + 1..] {
                    if !groups[i1].can_merge(&groups[i2], !use_lower, arena)? {
                        continue;
                    }
                    has_merges = true;
                    arena.counts.jitter_merged += 1;
                    // JitterUmiMerger.prioritiseFirstCoords
                    let first_has_priority = groups[i1].pos_lower <= groups[i2].pos_lower && groups[i1].pos_upper <= groups[i2].pos_upper;
                    if first_has_priority {
                        let fragments = groups[i2].fragments.clone();
                        groups[i1].fragments.extend(fragments);
                        groups[i2].merged = true;
                    } else {
                        let fragments = groups[i1].fragments.clone();
                        groups[i2].fragments.extend(fragments);
                        groups[i1].merged = true;
                        break;
                    }
                }
            }
        }
    }
    if !has_merges {
        return Ok(());
    }

    // both maps hold the same groups: walk them once, by lower position
    for m in pos_groups[0].order() {
        for &gi in pos_groups[0].value(m) {
            let group = &groups[gi];
            if group.merged {
                for f in &group.fragments {
                    match f {
                        Frag::Group(g) => list_remove(umi_groups, *g),
                        Frag::Single(s) => list_remove(single_fragments, *s),
                    }
                }
            } else if group.fragments.len() > 1 {
                let mut original: Option<usize> = None;
                for &fragment in &group.fragments {
                    let Some(target) = original else {
                        original = Some(match fragment {
                            Frag::Single(s) => {
                                arena.extract_single_umi(s);
                                let g = arena.single_to_group(s);
                                umi_groups.push(g);
                                list_remove(single_fragments, s);
                                g
                            }
                            Frag::Group(g) => g,
                        });
                        continue;
                    };
                    let frag_forward = arena.coords(fragment).frag_forward;
                    let reads_to_add = match fragment {
                        Frag::Single(s) => vec![arena.take_single_read(s)],
                        Frag::Group(g) => {
                            let source = arena.group_mut(g);
                            let mut reads = std::mem::take(&mut source.reads);
                            reads.append(&mut source.non_consensus_reads);
                            reads
                        }
                    };
                    let t = arena.group_mut(target);
                    t.non_consensus_reads.extend(reads_to_add);
                    if t.coords.frag_forward != frag_forward {
                        t.dual_strand = true;
                    }
                }
            }
        }
    }
    Ok(())
}

/// DuplicateGroupBuilder.processDuplicateGroups with UMIs: the popped fragments become UMI groups and single reads
pub fn process_duplicate_groups(config: &UmiConfig, is_illumina: bool, fcr: &mut FragmentCoordReads, stats: &mut Statistics) -> Result<()> {
    let duplicate_groups = std::mem::take(&mut fcr.duplicate_groups);
    let single_reads = std::mem::take(&mut fcr.single_reads);
    let (groups, singles) = process_umi_groups(config, is_illumina, duplicate_groups, single_reads, true, &mut stats.umi)?;
    for group in &groups {
        stats.add_umi_group(group.total_read_count(), group.dual_strand);
    }
    fcr.duplicate_groups = groups;
    fcr.single_reads = singles;
    Ok(())
}

/// PolyGUmiMerger.mergeGroups on its own, for the tests
#[cfg(test)]
pub fn poly_g_merge(config: &UmiConfig, umi_groups: Vec<DuplicateGroup>, single_reads: Vec<ReadInfo>) -> (Vec<DuplicateGroup>, Vec<ReadInfo>) {
    let mut arena = Arena::new(config);
    let mut groups: Vec<usize> = umi_groups.into_iter().map(|g| arena.add_group(g)).collect();
    let mut singles: Vec<usize> = single_reads.into_iter().map(|s| arena.add_single(s)).collect();
    poly_g_merge_groups(&mut arena, &mut groups, &mut singles);
    (
        groups.into_iter().map(|g| arena.groups[g].take().unwrap()).collect(),
        singles.into_iter().map(|s| arena.singles[s].info.take().unwrap()).collect(),
    )
}

/// JitterUmiMerger.mergeGroups on its own, for the tests
#[cfg(test)]
pub fn jitter_merge(config: &UmiConfig, umi_groups: Vec<DuplicateGroup>, single_reads: Vec<ReadInfo>) -> (Vec<DuplicateGroup>, Vec<ReadInfo>) {
    let mut arena = Arena::new(config);
    let mut groups: Vec<usize> = umi_groups.into_iter().map(|g| arena.add_group(g)).collect();
    let mut singles: Vec<usize> = single_reads.into_iter().map(|s| arena.add_single(s)).collect();
    jitter_merge_groups(&mut arena, &mut groups, &mut singles).unwrap();
    (
        groups.into_iter().map(|g| arena.groups[g].take().unwrap()).collect(),
        singles.into_iter().map(|s| arena.singles[s].info.take().unwrap()).collect(),
    )
}
