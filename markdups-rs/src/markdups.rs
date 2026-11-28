//! MarkDuplicates (picard.sam.markduplicates.MarkDuplicates, Picard 3.3.0) for coordinate-sorted input.
//!
//! Pass 1 reads the input and builds, for each primary mapped read, a fragment end (5' position and
//! orientation) and, once both mates are seen, a pair end. Picard sorts all of them and walks the
//! sorted lists in chunks of equal position; the chunks of one reference only hold ends of that
//! reference (and pairs with both mates on it), so here each reference's ends are sorted and walked
//! when the input moves past it, and the pairs across references at the end. Duplicates are kept as a
//! bit per record index. Pass 2 reads the input again and writes it with the duplicate flags, the
//! DT and PG tags, and counts the metrics.

use anyhow::{Context, Result, bail};
use hashbrown::HashMap;

use crate::bai::{BaiBuilder, IndexEntry};
use crate::bam::{self, Rec, SamHeader, SortOrder, TagEdits, encode_record, written_bin};
use crate::bgzf::BgzfWriter;
use crate::compat;
use crate::input::InputSet;
use crate::metrics::{DuplicationMetrics, SetSizeHistograms, java_date_now, write_metrics};
use crate::optical::{Location, NO_LOCATION, NameParser, OpticalEnd, OpticalFinder};

pub const F: u8 = 0;
pub const R: u8 = 1;
pub const FF: u8 = 2;
pub const FR: u8 = 3;
pub const RR: u8 = 4;
pub const RF: u8 = 5;

/// ReadEnds.getOrientationByte
#[inline]
fn orientation_byte(read1_negative: bool, read2_negative: bool) -> u8 {
    match (read1_negative, read2_negative) {
        (true, true) => RR,
        (true, false) => RF,
        (false, true) => FR,
        (false, false) => FF,
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scoring {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
    Random,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tagging {
    DontTag,
    OpticalOnly,
    All,
}

pub struct Config {
    pub inputs: Vec<String>,
    pub output: String,
    pub metrics_file: String,
    pub create_index: bool,
    pub remove_duplicates: bool,
    pub remove_sequencing_duplicates: bool,
    pub tagging: Tagging,
    pub clear_dt: bool,
    pub add_pg_tag_to_reads: bool,
    pub program_record_id: Option<String>,
    pub program_group_version: Option<String>,
    pub program_group_command_line: Option<String>,
    pub program_group_name: Option<String>,
    pub comments: Vec<String>,
    pub read_name_regex: Option<String>,
    pub pixel_distance: i32,
    pub max_optical_set_size: i64,
    pub scoring: Scoring,
    pub assume_sort_order: Option<SortOrder>,
    pub compression: i32,
    pub threads: usize,
    pub stringency: crate::validate::Stringency,
    pub command_line: String,
    /// the @PG VN: GATK's version string
    pub version: String,
}

/// ReadEndsForMarkDuplicates of a pair (or of the first mate seen, while waiting for the second).
#[derive(Clone, Copy, Debug)]
struct PairEnd {
    index1: i64,
    index2: i64,
    ref1: i32,
    coord1: i32,
    ref2: i32,
    coord2: i32,
    loc: Location,
    score: i16,
    lib: i16,
    read_group: i16,
    orient: u8,
    orient_optical: i8,
}

/// ReadEndsForMarkDuplicates of a fragment (one read); its reference is the one being processed.
#[derive(Clone, Copy, Debug)]
struct FragEnd {
    index: i64,
    coord: i32,
    /// the mate's reference if paired with a mapped mate, else -1 (ReadEnds.isPaired)
    ref2: i32,
    loc: Location,
    score: i16,
    lib: i16,
    orient: u8,
}

impl FragEnd {
    #[inline]
    fn is_paired(&self) -> bool {
        self.ref2 != -1
    }
}

/// A bit per record index.
#[derive(Default)]
pub struct BitSet {
    words: Vec<u64>,
}

impl BitSet {
    #[inline]
    pub fn set(&mut self, i: i64) {
        let i = i as usize;
        let w = i >> 6;
        if w >= self.words.len() {
            self.words.resize((w + 1).max(self.words.len() * 2), 0);
        }
        self.words[w] |= 1 << (i & 63);
    }
    #[inline]
    pub fn get(&self, i: i64) -> bool {
        let i = i as usize;
        self.words.get(i >> 6).is_some_and(|w| w >> (i & 63) & 1 == 1)
    }
    pub fn count(&self) -> u64 {
        self.words.iter().map(|w| w.count_ones() as u64).sum()
    }
}

/// htsjdk's Murmur3(seed).hashUnencodedChars, over the name's bytes as Latin-1 chars.
fn murmur3_chars(seed: i32, s: &[u8]) -> i32 {
    let mix_k1 = |k: i32| (k.wrapping_mul(0xcc9e2d51u32 as i32)).rotate_left(15).wrapping_mul(0x1b873593);
    let mix_h1 = |h: i32, k: i32| (h ^ k).rotate_left(13).wrapping_mul(5).wrapping_add(0xe6546b64u32 as i32);
    let mut h1 = seed;
    let n = s.len();
    let mut i = 1;
    while i < n {
        let k1 = (s[i - 1] as i32) | ((s[i] as i32) << 16);
        h1 = mix_h1(h1, mix_k1(k1));
        i += 2;
    }
    if n & 1 == 1 {
        h1 ^= mix_k1(s[n - 1] as i32);
    }
    let mut h = h1 ^ (2 * n as i32);
    h ^= ((h as u32) >> 16) as i32;
    h = h.wrapping_mul(0x85ebca6bu32 as i32);
    h ^= ((h as u32) >> 13) as i32;
    h = h.wrapping_mul(0xc2b2ae35u32 as i32);
    h ^= ((h as u32) >> 16) as i32;
    h
}

/// DuplicateScoringStrategy.computeDuplicateScore
fn duplicate_score(rec: &Rec, ops: &[u32], scoring: Scoring) -> i16 {
    let mut score: i16 = 0;
    match scoring {
        Scoring::SumOfBaseQualities => {
            let q = rec.qual();
            // BAM's missing qualities (first byte 0xff) read as none; quality bytes are Java bytes
            let sum: i32 = if !q.is_empty() && q[0] == 0xff { 0 } else { q.iter().map(|&b| b as i8 as i32).filter(|&b| b >= 15).sum() };
            score = score.wrapping_add(sum.min(i16::MAX as i32 / 2) as i16);
        }
        Scoring::TotalMappedReferenceLength => {
            if !rec.is_unmapped() {
                score = Rec::reference_length(ops).min(i16::MAX as i32 / 2) as i16;
            }
        }
        Scoring::Random => {
            score = score.wrapping_add((murmur3_chars(1, rec.name()) & 0x3fff) as i16);
            score = score.wrapping_sub(i16::MIN / 4);
        }
    }
    if rec.flag() & bam::FLAG_QC_FAIL != 0 {
        score = score.wrapping_add(i16::MIN / 2);
    }
    score
}

/// The library names and ids (LibraryIdGenerator).
pub struct Libraries {
    /// read group id -> library name, from the header
    rg_library: HashMap<Vec<u8>, Vec<u8>>,
    /// read group id -> its position in the header's list
    rg_index: HashMap<Vec<u8>, i16>,
    n_rg: i16,
    ids: HashMap<Vec<u8>, i16>,
    next_id: i16,
}

pub const UNKNOWN_LIBRARY: &[u8] = b"Unknown Library";

impl Libraries {
    pub fn new(h: &SamHeader) -> Self {
        let mut rg_library = HashMap::new();
        let mut rg_index = HashMap::new();
        for (i, rg) in h.rg.iter().enumerate() {
            rg_library.insert(rg.id.clone(), rg.get(b"LB").map(|l| l.to_vec()).unwrap_or_else(|| UNKNOWN_LIBRARY.to_vec()));
            rg_index.entry(rg.id.clone()).or_insert(i as i16);
        }
        Libraries { rg_library, rg_index, n_rg: h.rg.len() as i16, ids: HashMap::new(), next_id: 1 }
    }
    /// LibraryIdGenerator.getLibraryName
    pub fn name<'a>(&'a self, rg: Option<&[u8]>) -> &'a [u8] {
        rg.and_then(|id| self.rg_library.get(id)).map(|v| v.as_slice()).unwrap_or(UNKNOWN_LIBRARY)
    }
    fn id(&mut self, rg: Option<&[u8]>) -> i16 {
        let name: &[u8] = rg.and_then(|id| self.rg_library.get(id)).map(|v| v.as_slice()).unwrap_or(UNKNOWN_LIBRARY);
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.next_id;
        self.next_id = self.next_id.wrapping_add(1);
        self.ids.insert(name.to_vec(), id);
        id
    }
    /// The position of the read group in the header (the number of read groups if not there; 0 without
    /// a tag), as buildReadEnds counts it for optical duplicates.
    fn read_group_index(&self, rg: Option<&[u8]>) -> i16 {
        match rg {
            None => 0,
            Some(id) => self.rg_index.get(id).copied().unwrap_or(self.n_rg),
        }
    }
}

/// What pass 1 finds.
pub struct Duplicates {
    pub dup: BitSet,
    pub optical: Option<BitSet>,
    pub num_duplicate_indices: u64,
    pub sets: SetSizeHistograms,
    /// optical duplicates by library id
    pub optical_by_library: HashMap<i16, f64>,
    pub libraries: Libraries,
    /// PG values seen on the records, in Java HashSet order (None: no PG tag)
    pub pg_seen: Vec<Option<Vec<u8>>>,
    pub records: i64,
    pub unmatched_pairs: usize,
}

struct Pass1<'c> {
    cfg: &'c Config,
    parser: NameParser,
    finder: OpticalFinder,
    track_optical_indexes: bool,
    frags: Vec<FragEnd>,
    pairs: Vec<PairEnd>,
    cross_pairs: Vec<PairEnd>,
    pending: HashMap<Vec<u8>, PairEnd>,
    current_ref: i32,
    dup: BitSet,
    optical: BitSet,
    num_duplicate_indices: u64,
    sets: SetSizeHistograms,
    optical_by_library: HashMap<i16, f64>,
    libraries: Libraries,
    key: Vec<u8>,
}

impl Pass1<'_> {
    fn add_duplicate(&mut self, index: i64) {
        self.dup.set(index);
        self.num_duplicate_indices += 1;
    }

    /// Sorts and walks the ends of the reference just finished.
    fn flush_reference(&mut self) {
        let mut pairs = std::mem::take(&mut self.pairs);
        self.process_pairs(&mut pairs);
        pairs.clear();
        self.pairs = pairs;
        let mut frags = std::mem::take(&mut self.frags);
        self.process_fragments(&mut frags);
        frags.clear();
        self.frags = frags;
    }

    /// The pair half of generateDuplicateIndexes: chunks of equal (library, read1 position, orientation,
    /// read2 position), sorted by ReadEndsMDComparator.
    fn process_pairs(&mut self, pairs: &mut [PairEnd]) {
        pairs.sort_unstable_by(|a, b| {
            (a.lib, a.ref1, a.coord1, a.orient, a.ref2, a.coord2, a.loc.tile, a.loc.x, a.loc.y, a.index1, a.index2)
                .cmp(&(b.lib, b.ref1, b.coord1, b.orient, b.ref2, b.coord2, b.loc.tile, b.loc.x, b.loc.y, b.index1, b.index2))
        });
        let mut i = 0;
        while i < pairs.len() {
            let mut j = i + 1;
            let a = pairs[i];
            while j < pairs.len() {
                let b = &pairs[j];
                if !(a.lib == b.lib && a.ref1 == b.ref1 && a.coord1 == b.coord1 && a.orient == b.orient && a.ref2 == b.ref2 && a.coord2 == b.coord2) {
                    break;
                }
                j += 1;
            }
            let chunk = &pairs[i..j];
            if chunk.len() > 1 {
                self.mark_duplicate_pairs(chunk);
            } else {
                // addSingletonToCount
                *self.sets.all_sets.entry(1).or_insert(0.0) += 1.0;
                *self.sets.non_optical_sets.entry(1).or_insert(0.0) += 1.0;
            }
            i = j;
        }
    }

    fn mark_duplicate_pairs(&mut self, list: &[PairEnd]) {
        let mut max_score: i16 = 0;
        let mut best: Option<usize> = None;
        for (k, e) in list.iter().enumerate() {
            if e.score > max_score || best.is_none() {
                max_score = e.score;
                best = Some(k);
            }
        }
        let best = best.unwrap();
        let mut optical_flags = vec![false; list.len()];
        if !self.parser.is_off() {
            self.track_optical_duplicates(list, best, &mut optical_flags);
        } else if !compat::java("MD6") {
            // Picard counts the set sizes inside its optical duplicate tracking, so without a read name
            // regex the histograms hold only the single pairs (MD6); corrected, every set counts
            *self.sets.all_sets.entry(list.len() as u64).or_insert(0.0) += 1.0;
            *self.sets.non_optical_sets.entry(list.len() as u64).or_insert(0.0) += 1.0;
        }
        for (k, e) in list.iter().enumerate() {
            if k != best {
                self.add_duplicate(e.index1);
                if e.index2 != e.index1 {
                    self.add_duplicate(e.index2);
                }
                if optical_flags[k] && self.track_optical_indexes {
                    self.optical.set(e.index1);
                    if e.index2 != e.index1 {
                        self.optical.set(e.index2);
                    }
                }
            }
        }
    }

    /// AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates: FR and RF pairs (by first and
    /// second of pair) are looked at separately.
    fn track_optical_duplicates(&mut self, list: &[PairEnd], keeper: usize, flags: &mut [bool]) {
        let has_fr = list.iter().any(|e| e.orient_optical == FR as i8);
        let has_rf = list.iter().any(|e| e.orient_optical == RF as i8);
        let mut n_optical = 0;
        let to_optical = |e: &PairEnd| OpticalEnd { read_group: e.read_group, loc: e.loc };
        if has_fr && has_rf {
            for want in [FR as i8, RF as i8] {
                let members: Vec<usize> = (0..list.len()).filter(|&k| list[k].orient_optical == want).collect();
                let sub: Vec<OpticalEnd> = members.iter().map(|&k| to_optical(&list[k])).collect();
                let kp = members.iter().position(|&k| k == keeper);
                let f = self.finder.find(&sub, kp);
                let mut count = 0;
                for (m, &k) in members.iter().enumerate() {
                    if f[m] {
                        flags[k] = true;
                        count += 1;
                    }
                }
                if count > 0 {
                    *self.optical_by_library.entry(list[members[0]].lib).or_insert(0.0) += count as f64;
                }
                n_optical += count;
            }
            if list.iter().any(|e| e.orient_optical != FR as i8 && e.orient_optical != RF as i8) {
                // Picard throws "Found an unexpected orientation"; FR/RF chunks hold only those
                panic!("Found an unexpected orientation");
            }
        } else {
            let sub: Vec<OpticalEnd> = list.iter().map(to_optical).collect();
            let f = self.finder.find(&sub, Some(keeper));
            let mut count = 0;
            for (k, &v) in f.iter().enumerate() {
                if v {
                    flags[k] = true;
                    count += 1;
                }
            }
            if count > 0 {
                *self.optical_by_library.entry(list[0].lib).or_insert(0.0) += count as f64;
            }
            n_optical = count;
        }
        // trackDuplicateCounts
        let size = list.len() as u64;
        *self.sets.all_sets.entry(size).or_insert(0.0) += 1.0;
        if size > n_optical as u64 {
            *self.sets.non_optical_sets.entry(size - n_optical as u64).or_insert(0.0) += 1.0;
        }
        if n_optical > 0 {
            *self.sets.optical_sets.entry(n_optical as u64 + 1).or_insert(0.0) += 1.0;
        }
    }

    /// The fragment half of generateDuplicateIndexes: chunks of equal (library, position, orientation);
    /// a chunk with a paired read marks all unpaired ones, else all but the best.
    fn process_fragments(&mut self, frags: &mut [FragEnd]) {
        frags.sort_unstable_by(|a, b| {
            (a.lib, a.coord, a.orient, a.ref2, a.loc.tile, a.loc.x, a.loc.y, a.index).cmp(&(b.lib, b.coord, b.orient, b.ref2, b.loc.tile, b.loc.x, b.loc.y, b.index))
        });
        let mut i = 0;
        while i < frags.len() {
            let a = frags[i];
            let mut j = i + 1;
            while j < frags.len() && frags[j].lib == a.lib && frags[j].coord == a.coord && frags[j].orient == a.orient {
                j += 1;
            }
            let chunk = &frags[i..j];
            if chunk.len() > 1 {
                let contains_pairs = chunk.iter().any(|e| e.is_paired());
                if contains_pairs {
                    for e in chunk {
                        if !e.is_paired() {
                            self.add_duplicate(e.index);
                        }
                    }
                } else {
                    let mut max_score: i16 = 0;
                    let mut best: Option<usize> = None;
                    for (k, e) in chunk.iter().enumerate() {
                        if e.score > max_score || best.is_none() {
                            max_score = e.score;
                            best = Some(k);
                        }
                    }
                    for (k, e) in chunk.iter().enumerate() {
                        if Some(k) != best {
                            self.add_duplicate(e.index);
                        }
                    }
                }
            }
            i = j;
        }
    }
}

/// A Java-HashSet-ordered key for the PG values seen (null hashes to 0).
#[derive(Clone, PartialEq, Eq, Hash)]
struct PgKey(Option<Vec<u8>>);

impl hmf_common::java_hash::JavaHash for PgKey {
    fn java_hash(&self) -> i32 {
        match &self.0 {
            None => 0,
            Some(s) => s.iter().fold(0i32, |h, &c| h.wrapping_mul(31).wrapping_add(c as i32)),
        }
    }
}

/// Checks the coordinate order as htsjdk's writer does for presorted output (SAMSortOrderChecker).
struct OrderCheck {
    last: Option<(i32, i32)>,
}

impl OrderCheck {
    fn check(&mut self, rec: &Rec, header: &SamHeader) -> Result<()> {
        let key = (rec.ref_id(), rec.pos() + 1);
        if let Some(prev) = self.last {
            let bad = match (prev.0 == -1, key.0 == -1) {
                (true, true) => false,
                (true, false) => true,
                (false, true) => false,
                (false, false) => (prev.0, prev.1) > (key.0, key.1),
            };
            if bad {
                let name = |r: i32| if r < 0 { "*".to_string() } else { String::from_utf8_lossy(&header.sq[r as usize].name).to_string() };
                bail!(
                    "Alignments added out of order in SAMFileWriterImpl.addAlignment. Sort order is coordinate. Offending records are at [{}:{}] and [{}:{}]",
                    name(prev.0),
                    prev.1,
                    name(key.0),
                    key.1
                );
            }
        }
        self.last = Some(key);
        Ok(())
    }
}

/// Pass 1: buildSortedReadEndLists and generateDuplicateIndexes.
pub fn find_duplicates(cfg: &Config, inputs: &dyn InputSet, reader_threads: usize) -> Result<Duplicates> {
    let (header, mut src) = inputs.open(reader_threads)?;
    let sort_order = header.sort_order();
    if sort_order == SortOrder::Queryname {
        bail!("markdups-rs supports coordinate-sorted input only (the header or ASSUME_SORT_ORDER says queryname)");
    }
    if sort_order != SortOrder::Coordinate {
        bail!(
            "This program requires input that are either coordinate or query sorted (according to the header, or at least ASSUME_SORT_ORDER and the content.) Found ASSUME_SORT_ORDER={} and header sortorder={}",
            cfg.assume_sort_order.map(|s| s.name()).unwrap_or("null"),
            sort_order.name()
        );
    }
    let mut p = Pass1 {
        cfg,
        parser: NameParser::new(cfg.read_name_regex.as_deref())?,
        finder: OpticalFinder { distance: cfg.pixel_distance, max_set_size: cfg.max_optical_set_size },
        track_optical_indexes: cfg.remove_sequencing_duplicates || cfg.tagging != Tagging::DontTag,
        frags: Vec::new(),
        pairs: Vec::new(),
        cross_pairs: Vec::new(),
        pending: HashMap::new(),
        current_ref: -1,
        dup: BitSet::default(),
        optical: BitSet::default(),
        num_duplicate_indices: 0,
        sets: SetSizeHistograms::default(),
        optical_by_library: HashMap::new(),
        libraries: Libraries::new(&header),
        key: Vec::new(),
    };
    let mut pg_seen: hmf_common::java_hash::JavaHashSet<PgKey> = hmf_common::java_hash::JavaHashSet::new();
    // the values already in pg_seen, to skip the set for the common case
    let mut pg_known: Vec<Option<Vec<u8>>> = Vec::new();
    let mut see_pg = |pg: Option<&[u8]>, set: &mut hmf_common::java_hash::JavaHashSet<PgKey>| {
        if !pg_known.iter().any(|k| k.as_deref() == pg) {
            pg_known.push(pg.map(|v| v.to_vec()));
            set.insert(PgKey(pg.map(|v| v.to_vec())));
        }
    };
    let mut order = OrderCheck { last: None };
    let mut buf = Vec::new();
    let mut index: i64 = 0;
    let mut done = false;
    let mut ops: Vec<u32> = Vec::new();
    while src.next(&mut buf)? {
        let rec = Rec::check(&buf).with_context(|| format!("record {index}"))?;
        let (rg, pg) = rec.rg_pg()?;
        crate::validate::check(&rec, &header, rg, index as u64 + 1, cfg.stringency)?;
        order.check(&rec, &header)?;
        if done {
            // Picard stops reading at the first unplaced read, so a PG value only on the reads after it
            // gets no @PG line and pass 2 removes it from those reads (MD7)
            if cfg.program_record_id.is_some() && !compat::java("MD7") {
                see_pg(pg, &mut pg_seen);
            }
            continue;
        }
        if cfg.program_record_id.is_some() {
            see_pg(pg, &mut pg_seen);
        }
        let ref_id = rec.ref_id();
        if ref_id != p.current_ref && ref_id >= 0 {
            if p.current_ref >= 0 {
                p.flush_reference();
            }
            p.current_ref = ref_id;
        }
        let flag = rec.flag();
        if rec.is_unmapped() {
            if ref_id == -1 {
                // the unplaced unmapped reads at the end: no reason to continue
                done = true;
                continue;
            }
        } else if !rec.is_secondary_or_supplementary() {
            ops.clear();
            ops.extend(rec.cigar_ops());
            let reverse = rec.is_reverse();
            let coord = if reverse { rec.unclipped_end(&ops) } else { rec.unclipped_start(&ops) };
            let score = duplicate_score(&rec, &ops, cfg.scoring);
            let paired_mapped = rec.is_paired() && !rec.mate_unmapped();
            let lib = p.libraries.id(rg);
            let (loc, read_group) = match p.parser.location(rec.name()) {
                Some(l) => (l, p.libraries.read_group_index(rg)),
                None => (NO_LOCATION, -1),
            };
            let frag = FragEnd { index, coord, ref2: if paired_mapped { rec.mate_ref_id() } else { -1 }, loc, score, lib, orient: if reverse { R } else { F } };
            p.frags.push(frag);
            if paired_mapped {
                // the key: read group id and read name (MD4: a read without a known read group crashes
                // Picard; MD5: the two are joined without a separator)
                let rg_id: &[u8] = match rg.filter(|id| p.libraries.rg_index.contains_key(*id)) {
                    Some(id) => id,
                    None => {
                        if compat::java("MD4") {
                            bail!(
                                "java.lang.NullPointerException: Cannot invoke \"htsjdk.samtools.SAMReadGroupRecord.getReadGroupId()\" because the return value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null"
                            );
                        }
                        rg.unwrap_or(b"")
                    }
                };
                let mut key = std::mem::take(&mut p.key);
                key.clear();
                key.extend_from_slice(&ref_id.to_le_bytes());
                key.extend_from_slice(rg_id);
                if !compat::java("MD5") {
                    // a byte that cannot be in a read group id or read name
                    key.push(0);
                }
                key.extend_from_slice(rec.name());
                match p.pending.remove(&key) {
                    None => {
                        let pe = PairEnd {
                            index1: index,
                            index2: -1,
                            ref1: ref_id,
                            coord1: coord,
                            ref2: rec.mate_ref_id(),
                            coord2: -1,
                            loc,
                            score,
                            lib,
                            read_group,
                            orient: frag.orient,
                            orient_optical: -1,
                        };
                        // stored under the mate's reference, where the mate will look for it
                        key[0..4].copy_from_slice(&rec.mate_ref_id().to_le_bytes());
                        if p.pending.insert(key.clone(), pe).is_some() {
                            bail!("Putting value into PairInfoMap that already existed. {}: read seen twice as a primary record", rec.mate_ref_id());
                        }
                    }
                    Some(mut pe) => {
                        let (mate_ref, mate_coord) = (ref_id, coord);
                        pe.orient_optical = if rec.is_first_of_pair() {
                            orientation_byte(reverse, pe.orient == R) as i8
                        } else {
                            orientation_byte(pe.orient == R, reverse) as i8
                        };
                        if mate_ref > pe.ref1 || (mate_ref == pe.ref1 && mate_coord >= pe.coord1) {
                            pe.ref2 = mate_ref;
                            pe.coord2 = mate_coord;
                            pe.index2 = index;
                            pe.orient = orientation_byte(pe.orient == R, reverse);
                            if pe.ref2 == pe.ref1 && pe.coord2 == pe.coord1 && pe.orient == RF {
                                pe.orient = FR;
                            }
                        } else {
                            pe.ref2 = pe.ref1;
                            pe.coord2 = pe.coord1;
                            pe.index2 = pe.index1;
                            pe.ref1 = mate_ref;
                            pe.coord1 = mate_coord;
                            pe.index1 = index;
                            pe.orient = orientation_byte(reverse, pe.orient == R);
                        }
                        pe.score = pe.score.wrapping_add(score);
                        if pe.ref1 == ref_id && pe.ref2 == ref_id {
                            p.pairs.push(pe);
                        } else {
                            p.cross_pairs.push(pe);
                        }
                    }
                }
                p.key = key;
            }
        }
        let _ = flag;
        index += 1;
    }
    if p.current_ref >= 0 {
        p.flush_reference();
    }
    let mut cross = std::mem::take(&mut p.cross_pairs);
    p.process_pairs(&mut cross);
    let pg_order: Vec<Option<Vec<u8>>> = pg_seen.iter().map(|k| k.0.clone()).collect();
    Ok(Duplicates {
        num_duplicate_indices: p.num_duplicate_indices,
        optical: if p.track_optical_indexes { Some(p.optical) } else { None },
        dup: p.dup,
        sets: p.sets,
        optical_by_library: p.optical_by_library,
        unmatched_pairs: p.pending.len(),
        libraries: p.libraries,
        pg_seen: pg_order,
        records: index,
    })
}

/// SAMFileHeader.PgIdGenerator
struct PgIdGenerator {
    taken: hashbrown::HashSet<Vec<u8>>,
    counter: i32,
}

impl PgIdGenerator {
    fn new(h: &SamHeader) -> Self {
        let taken: hashbrown::HashSet<Vec<u8>> = h.pg.iter().map(|p| p.id.clone()).collect();
        let counter = taken.len() as i32;
        PgIdGenerator { taken, counter }
    }
    fn non_colliding(&mut self, id: &[u8]) -> Vec<u8> {
        if self.taken.insert(id.to_vec()) {
            self.counter += 1;
            return id.to_vec();
        }
        loop {
            let mut n = id.to_vec();
            n.push(b'.');
            n.extend_from_slice(base36(self.counter).as_bytes());
            self.counter += 1;
            if self.taken.insert(n.clone()) {
                return n;
            }
        }
    }
}

/// SamFileHeaderMerger.positiveFourDigitBase36Str
pub fn base36(mut v: i32) -> String {
    if v == 0 {
        return "0".into();
    }
    let mut s = Vec::new();
    while v > 0 {
        let d = (v % 36) as u8;
        s.push(if d < 10 { b'0' + d } else { b'A' + d - 10 });
        v /= 36;
    }
    s.reverse();
    String::from_utf8(s).unwrap()
}

/// The thread split of pass 2: (reader threads, compression threads).
pub fn pass2_threads(threads: usize) -> (usize, usize) {
    if threads <= 1 {
        return (0, 0);
    }
    let workers = threads - 1;
    let c = (workers * 3).div_ceil(4).max(1);
    (workers - c, c)
}

/// The index file htsjdk's BAMFileWriter writes beside `output`.
pub fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(base) => format!("{base}.bai"),
        None => format!("{output}.bai"),
    }
}

/// MarkDuplicates.doWork after pass 1: writes the output BAM (and index) and the metrics.
pub fn write_output(cfg: &Config, inputs: &dyn InputSet, dups: &mut Duplicates, reader_threads: usize, compress_threads: usize) -> Result<()> {
    let (header, mut src) = inputs.open(reader_threads)?;
    let sort_order = header.sort_order();
    let mut out_header = header.clone();
    if cfg.assume_sort_order == Some(SortOrder::Queryname) {
        out_header.set_hd(b"GO", Some(b"query"));
        out_header.set_hd(b"SO", Some(b"unknown"));
    }
    for c in &cfg.comments {
        out_header.add_comment(c.as_bytes());
    }
    // getChainedPgIds: one @PG per PG value seen on the records, chained to it
    // few entries (one per PG value seen): a list
    let mut chained: Vec<(Option<Vec<u8>>, Vec<u8>)> = Vec::new();
    if let Some(pg_id) = &cfg.program_record_id {
        let mut ids = PgIdGenerator::new(&out_header);
        for existing in &dups.pg_seen {
            let id = ids.non_colliding(pg_id.as_bytes());
            let mut rec = bam::IdRecord { id: id.clone(), attrs: Vec::new() };
            rec.set(b"VN", Some(cfg.program_group_version.as_deref().unwrap_or(&cfg.version).as_bytes()));
            rec.set(b"CL", Some(cfg.program_group_command_line.as_deref().unwrap_or(&cfg.command_line).as_bytes()));
            if let Some(n) = &cfg.program_group_name {
                rec.set(b"PN", Some(n.as_bytes()));
            }
            if let Some(pp) = existing {
                rec.set(b"PP", Some(pp));
            }
            out_header.pg.push(rec);
            chained.push((existing.clone(), id));
        }
    }
    let create_index = cfg.create_index && out_header.sort_order() == SortOrder::Coordinate;
    let mut w = BgzfWriter::create(&cfg.output, cfg.compression, compress_threads)?;
    bam::write_header(&mut w, &out_header)?;
    let mut bai = if create_index { Some(BaiBuilder::new(out_header.sq.len())) } else { None };
    let mut pending_index: std::collections::VecDeque<(IndexEntry, u64, u64)> = std::collections::VecDeque::new();

    let mut metrics: std::collections::BTreeMap<Vec<u8>, DuplicationMetrics> = std::collections::BTreeMap::new();
    for rg in &header.rg {
        let lib = rg.get(b"LB").unwrap_or(crate::markdups::UNKNOWN_LIBRARY).to_vec();
        metrics.entry(lib.clone()).or_insert_with(|| DuplicationMetrics::new(&lib));
    }
    let set_pg = cfg.program_record_id.is_some() && cfg.add_pg_tag_to_reads;
    let mut buf = Vec::new();
    let mut out = Vec::with_capacity(1 << 16);
    let mut tag_scratch = Vec::new();
    let mut ops: Vec<u32> = Vec::new();
    let mut index: i64 = 0;
    let mut written: u64 = 0;
    while src.next(&mut buf)? {
        let rec = Rec::check(&buf)?;
        let (rg, old_pg) = rec.rg_pg()?;
        let lib = dups.libraries.name(rg);
        let m = match metrics.get_mut(lib) {
            Some(m) => m,
            None => metrics.entry(lib.to_vec()).or_insert_with(|| DuplicationMetrics::new(lib)),
        };
        m.add_read(&rec);
        let is_dup = dups.dup.get(index);
        let mut flag = rec.flag() & !bam::FLAG_DUPLICATE;
        if is_dup {
            flag |= bam::FLAG_DUPLICATE;
            m.add_duplicate(&rec);
        }
        let is_optical = dups.optical.as_ref().is_some_and(|o| o.get(index));
        index += 1;
        if cfg.remove_duplicates && is_dup {
            continue;
        }
        if cfg.remove_sequencing_duplicates && is_optical {
            continue;
        }
        let set_dt: Option<&[u8]> = if cfg.tagging != Tagging::DontTag && is_dup {
            if is_optical {
                Some(b"SQ")
            } else if cfg.tagging == Tagging::All {
                Some(b"LB")
            } else {
                None
            }
        } else {
            None
        };
        let pg_value: Option<Option<&[u8]>> = if set_pg {
            Some(chained.iter().find(|(k, _)| k.as_deref() == old_pg).map(|(_, v)| v.as_slice()))
        } else {
            None
        };
        let edits = TagEdits { clear_dt: cfg.clear_dt, set_dt, set_pg: pg_value, force: src.last_touched() };
        ops.clear();
        ops.extend(rec.cigar_ops());
        out.clear();
        encode_record(&mut out, &rec, flag, &ops, &edits, &mut tag_scratch)?;
        let start_u = w.position();
        w.write(&out)?;
        written += 1;
        if let Some(b) = &mut bai {
            let entry = IndexEntry {
                ref_id: rec.ref_id(),
                start: rec.pos() + 1,
                end: rec.alignment_end(&ops),
                bin: if rec.pos() + 1 == 0 { 0 } else { bam::indexing_bin(&rec, &ops) as u32 },
                unmapped: rec.is_unmapped(),
                chunk_start: 0,
                chunk_end: 0,
            };
            let _ = written_bin;
            pending_index.push_back((entry, start_u, w.position()));
            while let Some((e, s, t)) = pending_index.front() {
                let (Some(vs), Some(vt)) = (w.virtual_offset(*s), w.virtual_offset(*t)) else { break };
                let mut e = *e;
                e.chunk_start = vs;
                e.chunk_end = vt;
                b.add(&e)?;
                pending_index.pop_front();
            }
        }
    }
    drop(src);
    w.finish()?;
    if let Some(mut b) = bai {
        while let Some((mut e, s, t)) = pending_index.pop_front() {
            e.chunk_start = w.virtual_offset(s).context("index offset")?;
            e.chunk_end = w.virtual_offset(t).context("index offset")?;
            b.add(&e)?;
        }
        let bytes = b.finish();
        std::fs::write(index_path(&cfg.output), bytes).with_context(|| format!("cannot write {}", index_path(&cfg.output)))?;
    }
    drop(w);
    let _ = sort_order;

    // finalizeAndWriteMetrics
    let mut rows: Vec<DuplicationMetrics> = Vec::new();
    for (name, mut m) in metrics {
        m.read_pairs_examined /= 2;
        m.read_pair_duplicates /= 2;
        if let Some(id) = dups.libraries.ids.get(&name) {
            if let Some(v) = dups.optical_by_library.get(id) {
                m.read_pair_optical_duplicates = *v as i64;
            }
        }
        m.calculate_derived_fields();
        rows.push(m);
    }
    let roi = if rows.len() == 1 { rows[0].roi_histogram() } else { None };
    let text = write_metrics(&cfg.command_line, &java_date_now(), &rows, roi.as_deref(), &dups.sets);
    std::fs::write(&cfg.metrics_file, text).with_context(|| format!("cannot write {}", cfg.metrics_file))?;
    crate::log_info!("Wrote {written} records; marked {} records as duplicates", dups.num_duplicate_indices);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn murmur3_as_htsjdk() {
        // Guava's murmur3_32(seed 1).hashUnencodedChars values
        assert_eq!(murmur3_chars(1, b""), 0x514e28b7u32 as i32);
    }

    #[test]
    fn base36_ids() {
        assert_eq!(base36(0), "0");
        assert_eq!(base36(3), "3");
        assert_eq!(base36(35), "Z");
        assert_eq!(base36(36), "10");
    }

    #[test]
    fn thread_split() {
        assert_eq!(pass2_threads(1), (0, 0));
        assert_eq!(pass2_threads(2), (0, 1));
        assert_eq!(pass2_threads(3), (0, 2));
        assert_eq!(pass2_threads(8), (1, 6));
    }

    #[test]
    fn index_names() {
        assert_eq!(index_path("d/x.md.bam"), "d/x.md.bai");
        assert_eq!(index_path("d/x.out"), "d/x.out.bai");
    }
}
