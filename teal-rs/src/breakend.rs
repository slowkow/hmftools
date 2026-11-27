//! BreakEndApp.findBreakEnds: SampleBreakEndFinder, TelbamBreakEndFinder, CandidateBreakEndFinder,
//! BreakEndSupportCounter, BreakEndEvidenceFilter and BreakEndTsvWriter.
//!
//! Telomeric break ends are soft-clip sites where the clipped part is telomeric and the aligned part is not.
//! Java's TelomericBreakEnd is a mutable object that the tumor and the germline analysis share: the germline
//! candidates start from the tumor's break ends, and merging close candidates moves the position of the
//! object it keeps. The port keeps the break ends in an arena (`Keys`) so that this sharing is reproduced
//! (TE4); corrected, the germline analysis gets copies.

use anyhow::{Context, Result};
use std::io::Write;
use std::sync::Arc;

use crate::aligner::{find_telomere_segment, matches_telomere};
use crate::read::{Read, ReadGroup};
use crate::telbam_reader::{read_telbam, GenomeRegion, TelbamGroups};
use crate::telo::{is_poly_gc, reverse_complement};

/// CandidateBreakEndFinder: BREAKEND_MERGE_DISTANCE
const BREAKEND_MERGE_DISTANCE: i32 = 2;
/// BreakEndSupportCounter
const MAX_DISTANCE_FROM_BREAKEND: i32 = 500;
const SPLIT_POSITION_TOLERANCE: i32 = 5;
const MIN_TELOMERE_MATCH_LENGTH: usize = 12;

/// BreakEndParams: the options and fixed values
#[derive(Clone, Debug)]
pub struct BreakEndParams {
    pub mark_duplicate_distance: i32,
    /// not applied: TelomereMatcher ignores its thresholds (TE3)
    #[allow(dead_code)]
    pub aligned_segment_telomere_match_threshold: f64,
    pub max_aligned_poly_base_threshold: f64,
    pub v37: bool,
    pub excluded_regions: Vec<GenomeRegion>,
    pub included_regions: Option<Vec<GenomeRegion>>,
    pub sample_id: String,
    pub tumor_telbam: String,
    pub germline_telbam: String,
    pub output_file: String,
}

impl BreakEndParams {
    pub fn new(v37: bool) -> BreakEndParams {
        BreakEndParams {
            mark_duplicate_distance: 60,
            aligned_segment_telomere_match_threshold: 0.8,
            max_aligned_poly_base_threshold: 0.85,
            v37,
            excluded_regions: crate::telbam_reader::excluded_regions(v37),
            included_regions: None,
            sample_id: String::new(),
            tumor_telbam: String::new(),
            germline_telbam: String::new(),
            output_file: String::new(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BreakEndType {
    // potential new telomere
    RightGTelomeric,
    LeftCTelomeric,
    // potential join with telomere
    RightCTelomeric,
    LeftGTelomeric,
}

impl BreakEndType {
    pub fn is_right_telomeric(self) -> bool {
        matches!(self, BreakEndType::RightGTelomeric | BreakEndType::RightCTelomeric)
    }
    pub fn is_g_telomeric(self) -> bool {
        matches!(self, BreakEndType::LeftGTelomeric | BreakEndType::RightGTelomeric)
    }
}

/// TelomericBreakEnd
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Key {
    pub ty: BreakEndType,
    pub chromosome: Arc<str>,
    pub position: i32,
    pub is_duplicate: bool,
}

pub type KeyId = usize;

/// The TelomericBreakEnd objects
#[derive(Default)]
pub struct Keys(pub Vec<Key>);

impl Keys {
    pub fn add(&mut self, k: Key) -> KeyId {
        self.0.push(k);
        self.0.len() - 1
    }
    /// TelomericBreakEnd.compareTo: type, chromosome (as text), position
    pub fn compare(&self, a: KeyId, b: KeyId) -> std::cmp::Ordering {
        let (x, y) = (&self.0[a], &self.0[b]);
        x.ty.cmp(&y.ty).then_with(|| x.chromosome.cmp(&y.chromosome)).then(x.position.cmp(&y.position))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignedReadType {
    /// split reads with the split section telomeric
    SplitReadTelomeric,
    /// reads across the break end with no split
    NotSplitRead,
    /// split read with the split section not telomeric
    SplitReadNotTelomeric,
    DiscordantPair,
}

impl AlignedReadType {
    pub fn is_split_read(self) -> bool {
        matches!(self, AlignedReadType::SplitReadTelomeric | AlignedReadType::SplitReadNotTelomeric)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PairedReadType {
    DiscordantPairTelomeric,
    DiscordantPairNotTelomeric,
    NotDiscordant,
}

use AlignedReadType as A;
use PairedReadType as P;

/// Fragment.supportFragTypes
pub const SUPPORT_FRAG_TYPES: [(AlignedReadType, PairedReadType); 5] = [
    (A::SplitReadTelomeric, P::DiscordantPairTelomeric),
    (A::SplitReadTelomeric, P::DiscordantPairNotTelomeric),
    (A::SplitReadTelomeric, P::NotDiscordant),
    (A::SplitReadNotTelomeric, P::DiscordantPairTelomeric),
    (A::DiscordantPair, P::DiscordantPairTelomeric),
];

/// BreakEndEvidenceFilter.supportSplitReadTypes
const SUPPORT_SPLIT_READ_TYPES: [(AlignedReadType, PairedReadType); 4] = [
    (A::SplitReadTelomeric, P::DiscordantPairTelomeric),
    (A::SplitReadTelomeric, P::DiscordantPairNotTelomeric),
    (A::SplitReadTelomeric, P::NotDiscordant),
    (A::SplitReadNotTelomeric, P::DiscordantPairTelomeric),
];

/// BreakEndEvidenceFilter.supportDiscordantPairTypes
const SUPPORT_DISCORDANT_PAIR_TYPES: [(AlignedReadType, PairedReadType); 3] = [
    (A::DiscordantPair, P::DiscordantPairTelomeric),
    (A::SplitReadTelomeric, P::DiscordantPairTelomeric),
    (A::SplitReadNotTelomeric, P::DiscordantPairTelomeric),
];

/// Fragment: a read group's evidence for a break end
#[derive(Clone, Debug)]
pub struct Fragment<'a> {
    pub group: usize,
    pub aligned_type: AlignedReadType,
    pub paired_type: PairedReadType,
    pub aligned_read: &'a Read,
    pub paired_read: &'a Read,
}

impl Fragment<'_> {
    fn types(&self) -> (AlignedReadType, PairedReadType) {
        (self.aligned_type, self.paired_type)
    }
}

/// TelomericBreakEndSupport
#[derive(Clone, Debug)]
pub struct Support<'a> {
    pub key: KeyId,
    pub fragments: Vec<Fragment<'a>>,
    /// the length of longestTelomereSegment
    pub longest_telomere_segment: Option<usize>,
    pub distance_to_telomere: Option<i32>,
    pub longest_split_read_align_length: i32,
}

impl Support<'_> {
    pub fn fragment_count(&self, f: impl Fn(&Fragment) -> bool) -> i32 {
        self.fragments.iter().filter(|x| f(x)).count() as i32
    }
    pub fn count_of(&self, a: AlignedReadType, p: PairedReadType) -> i32 {
        self.fragment_count(|f| f.aligned_type == a && f.paired_type == p)
    }
    pub fn support_fragment_count(&self) -> i32 {
        self.fragment_count(|f| SUPPORT_FRAG_TYPES.contains(&f.types()))
    }
    pub fn num_split_reads(&self) -> i32 {
        self.fragment_count(|f| f.aligned_type.is_split_read())
    }
    pub fn sum_mapq(&self, f: impl Fn(&Fragment) -> bool) -> i32 {
        self.fragments.iter().filter(|x| f(x)).map(|x| x.aligned_read.mapq as i32).sum()
    }
}

// ---- CandidateBreakEndFinder

#[derive(Clone, Copy, Debug)]
struct Candidate {
    key: KeyId,
    split_reads: usize,
}

/// CandidateBreakEndFinder.isInExcludedBaseRegion
pub fn is_in_excluded_base_region(params: &BreakEndParams, chromosome: &str, start: i32, end: i32) -> bool {
    let rank = hmf_common::regions::chromosome_rank;
    // LINC00486, where one end maps in the poly-G region
    let (lc, ls, le) = if params.v37 { ("2", 33141260, 33141700) } else { ("chr2", 32916190, 32916630) };
    let contains = |p: i32| p >= ls && p <= le;
    if (contains(start) || contains(end)) && rank(chromosome) == rank(lc) {
        return true;
    }
    params.excluded_regions.iter().any(|r| start >= r.start && start <= r.end && end >= r.start && end <= r.end && rank(chromosome) == rank(&r.chromosome))
}

/// CandidateBreakEndFinder.findCandidateBreakEnd: the type and position of a split read's break end
pub fn find_candidate_break_end(params: &BreakEndParams, r: &Read, boundary_zone: usize) -> Option<(BreakEndType, i32)> {
    if is_in_excluded_base_region(params, &r.chrom, r.start, r.end()) {
        return None;
    }
    let s = r.read_string();
    let left = r.left_soft_clip();
    let right = r.right_soft_clip();
    let left_bases = &s[..left];
    let right_bases = &s[s.len() - right..];
    let aligned = &s[left..s.len() - right];
    // poly-A/C/G/T aligned part
    let mut counts = [0usize; 4];
    for &b in aligned {
        match b {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let poly_ratio = *counts.iter().max().unwrap() as f64 / aligned.len() as f64;
    if poly_ratio > params.max_aligned_poly_base_threshold {
        return None;
    }
    if aligned.len() <= MIN_TELOMERE_MATCH_LENGTH {
        return None;
    }
    if is_poly_gc(aligned) {
        return None;
    }
    let mut tbe: Option<(BreakEndType, i32)> = None;
    if left_bases.len() >= MIN_TELOMERE_MATCH_LENGTH + boundary_zone {
        let m = &left_bases[..left_bases.len() - boundary_zone];
        if matches_telomere(m, MIN_TELOMERE_MATCH_LENGTH, true) {
            tbe = Some((BreakEndType::LeftGTelomeric, r.start));
        } else if matches_telomere(m, MIN_TELOMERE_MATCH_LENGTH, false) {
            tbe = Some((BreakEndType::LeftCTelomeric, r.start));
        }
    }
    if tbe.is_none() && right_bases.len() >= MIN_TELOMERE_MATCH_LENGTH + boundary_zone {
        let m = &right_bases[boundary_zone..];
        if matches_telomere(m, MIN_TELOMERE_MATCH_LENGTH, true) {
            tbe = Some((BreakEndType::RightGTelomeric, r.end()));
        } else if matches_telomere(m, MIN_TELOMERE_MATCH_LENGTH, false) {
            tbe = Some((BreakEndType::RightCTelomeric, r.end()));
        }
    }
    let (ty, pos) = tbe?;
    // the aligned part must not be telomeric of the same kind
    if matches_telomere(aligned, MIN_TELOMERE_MATCH_LENGTH, ty.is_g_telomeric()) {
        return None;
    }
    Some((ty, pos))
}

/// CandidateBreakEndFinder: the candidate break ends of the groups, consolidated
fn candidate_break_ends(params: &BreakEndParams, groups: &TelbamGroups, keys: &mut Keys, extra: &[KeyId]) -> Vec<KeyId> {
    let mut candidates: Vec<Candidate> = extra.iter().map(|&k| Candidate { key: k, split_reads: 0 }).collect();
    for rg in groups.iter() {
        for r in rg.all_reads() {
            if r.unmapped() {
                continue;
            }
            if let Some((ty, position)) = find_candidate_break_end(params, r, 1) {
                let key = keys.add(Key { ty, chromosome: r.chrom.clone(), position, is_duplicate: false });
                candidates.push(Candidate { key, split_reads: 1 });
            }
        }
    }
    // consolidatePotentialBreakEnds
    hmf_common::java_sort::java_sort(&mut candidates, |a, b| keys.compare(a.key, b.key));
    let same = |keys: &Keys, a: KeyId, b: KeyId| {
        let (x, y) = (&keys.0[a], &keys.0[b]);
        x.ty == y.ty && x.chromosome == y.chromosome
    };
    let mut combined: Vec<Candidate> = Vec::new();
    for tbe in candidates {
        if let Some(m) = combined.last_mut() {
            if same(keys, m.key, tbe.key) && keys.0[m.key].position == keys.0[tbe.key].position {
                m.split_reads += tbe.split_reads;
                continue;
            }
        }
        combined.push(tbe);
    }
    let mut out: Vec<Candidate> = Vec::new();
    let mut max_split = 0usize;
    for tbe in combined {
        if let Some(cur) = out.last_mut() {
            if same(keys, cur.key, tbe.key) && (keys.0[cur.key].position - keys.0[tbe.key].position).abs() <= BREAKEND_MERGE_DISTANCE {
                if tbe.split_reads > max_split {
                    // the position with the most split reads
                    let p = keys.0[tbe.key].position;
                    keys.0[cur.key].position = p;
                    max_split = tbe.split_reads;
                }
                cur.split_reads += tbe.split_reads;
                continue;
            }
        }
        max_split = tbe.split_reads;
        out.push(tbe);
    }
    out.into_iter().map(|c| c.key).collect()
}

// ---- BreakEndSupportCounter

/// findAlignedRead: the read aligned next to the break end and facing it
fn find_aligned_read<'a>(rg: &'a ReadGroup, key: &Key) -> Option<&'a Read> {
    let right = key.ty.is_right_telomeric();
    let mut aligned: Option<&Read> = None;
    for r in rg.all_reads() {
        if r.unmapped() || r.chrom != key.chromosome {
            continue;
        }
        if right {
            let e = r.end();
            if e <= key.position && e > key.position - MAX_DISTANCE_FROM_BREAKEND && aligned.is_none_or(|a| e > a.end()) {
                aligned = Some(r);
            }
        } else if r.start >= key.position && r.start < key.position + MAX_DISTANCE_FROM_BREAKEND && aligned.is_none_or(|a| r.start < a.start) {
            aligned = Some(r);
        }
    }
    let a = aligned?;
    if right && a.neg() && a.end() < key.position {
        return None;
    }
    if !right && !a.neg() && a.start > key.position {
        return None;
    }
    Some(a)
}

/// The part of a read a telomere segment is looked for in
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Part {
    LeftClip,
    RightClip,
    /// the bases as stored
    Stored,
    /// their reverse complement
    Reverse,
}

/// TelomereMatcher results by read, part and G or C richness: the counter looks at the same reads for every
/// break end near them
#[derive(Default)]
struct SegmentCache {
    map: std::collections::HashMap<(usize, Part, bool), Option<usize>, crate::prep::FxBuild>,
}

impl SegmentCache {
    fn bases(read: &Read, part: Part) -> Option<Vec<u8>> {
        match part {
            Part::LeftClip => read.left_soft_clip_bases().map(|b| b.to_vec()),
            Part::RightClip => read.right_soft_clip_bases().map(|b| b.to_vec()),
            Part::Stored => Some(read.read_string().to_vec()),
            Part::Reverse => Some(reverse_complement(read.read_string())),
        }
    }
    /// findTelomereSegment: the length of the matched segment
    fn segment(&mut self, read: &Read, part: Part, g_rich: bool) -> Option<usize> {
        let key = (read as *const Read as usize, part, g_rich);
        if let Some(v) = self.map.get(&key) {
            return *v;
        }
        let v = Self::bases(read, part).and_then(|b| find_telomere_segment(&b, g_rich)).map(|m| m.matched_sequence.len());
        self.map.insert(key, v);
        v
    }
    /// matchesTelomere: GGG (or CCC) somewhere and a segment of 12 matched bases
    fn matches(&mut self, read: &Read, part: Part, g_rich: bool) -> bool {
        let triple: &[u8] = if g_rich { b"GGG" } else { b"CCC" };
        let has = match Self::bases(read, part) {
            Some(b) => b.windows(3).any(|w| w == triple),
            None => false,
        };
        has && self.segment(read, part, g_rich).is_some()
    }
}

fn classify_aligned_read(key: &Key, a: &Read, cache: &mut SegmentCache) -> Option<AlignedReadType> {
    let right = key.ty.is_right_telomeric();
    let clip = if right && (a.end() - key.position).abs() <= SPLIT_POSITION_TOLERANCE && a.right_soft_clip() > 0 {
        Some(Part::RightClip)
    } else if !right && (a.start - key.position).abs() <= SPLIT_POSITION_TOLERANCE && a.left_soft_clip() > 0 {
        Some(Part::LeftClip)
    } else {
        None
    };
    if let Some(c) = clip {
        return Some(if cache.matches(a, c, key.ty.is_g_telomeric()) { A::SplitReadTelomeric } else { A::SplitReadNotTelomeric });
    }
    if a.start + SPLIT_POSITION_TOLERANCE < key.position && a.end() - SPLIT_POSITION_TOLERANCE > key.position {
        return Some(A::NotSplitRead);
    }
    if (right && a.end() < key.position && !a.neg()) || (!right && a.start > key.position && a.neg()) {
        return Some(A::DiscordantPair);
    }
    None
}

/// The paired read as it points away from the aligned read
fn oriented_part(aligned: &Read, paired: &Read) -> Part {
    if aligned.neg() == paired.neg() { Part::Reverse } else { Part::Stored }
}

fn classify_paired_read(key: &Key, paired: &Read, aligned: &Read, cache: &mut SegmentCache) -> PairedReadType {
    if key.ty.is_right_telomeric() == aligned.neg() {
        return P::NotDiscordant;
    }
    if cache.matches(paired, oriented_part(aligned, paired), key.ty.is_g_telomeric()) {
        P::DiscordantPairTelomeric
    } else {
        P::DiscordantPairNotTelomeric
    }
}

fn break_end_fragment<'a>(gi: usize, rg: &'a ReadGroup, key: &Key, cache: &mut SegmentCache) -> Option<Fragment<'a>> {
    let aligned = find_aligned_read(rg, key)?;
    if aligned.duplicate() {
        return None;
    }
    let paired = if aligned.first_in_pair() { rg.second_of_pair() } else { rg.first_of_pair() }?;
    let aligned_type = classify_aligned_read(key, aligned, cache)?;
    let paired_type = classify_paired_read(key, paired, aligned, cache);
    Some(Fragment { group: gi, aligned_type, paired_type, aligned_read: aligned, paired_read: paired })
}

/// GRCh37 and GRCh38 chromosome lengths (hmf-common refgenome/lengths.*.tsv)
fn chromosome_length(chromosome: &str, v37: bool) -> Option<i32> {
    // HumanChromosome.fromString: a "chr" prefix (any case) is dropped; 1-22, X and Y only
    let bare = if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") { &chromosome[3..] } else { chromosome };
    hmf_common::ref_coords::chromosome_length(bare, v37)
}

/// countBreakEndSupports, for the groups with a read near each break end
fn count_break_end_supports<'a>(params: &BreakEndParams, break_ends: &[KeyId], keys: &Keys, groups: &'a [&'a ReadGroup]) -> Vec<Support<'a>> {
    // an index of the mapped reads by contig and by start and end: only groups with a read within
    // MAX_DISTANCE_FROM_BREAKEND of a break end can give it a fragment
    let mut by_end: Vec<(&str, i32, usize)> = Vec::new();
    let mut by_start: Vec<(&str, i32, usize)> = Vec::new();
    for (gi, rg) in groups.iter().enumerate() {
        for r in rg.all_reads() {
            if r.unmapped() {
                continue;
            }
            by_end.push((&r.chrom, r.end(), gi));
            by_start.push((&r.chrom, r.start, gi));
        }
    }
    by_end.sort_unstable();
    by_start.sort_unstable();
    let range = |v: &[(&str, i32, usize)], chr: &str, lo: i32, hi: i32| -> Vec<usize> {
        let from = v.partition_point(|x| (x.0, x.1) < (chr, lo));
        let to = v.partition_point(|x| (x.0, x.1) <= (chr, hi));
        let mut g: Vec<usize> = v[from..to].iter().map(|x| x.2).collect();
        g.sort_unstable();
        g.dedup();
        g
    };
    let mut cache = SegmentCache::default();
    let mut supports: Vec<Support> = Vec::with_capacity(break_ends.len());
    for &k in break_ends {
        let key = &keys.0[k];
        let candidates = if key.ty.is_right_telomeric() {
            range(&by_end, &key.chromosome, key.position - MAX_DISTANCE_FROM_BREAKEND + 1, key.position)
        } else {
            range(&by_start, &key.chromosome, key.position, key.position + MAX_DISTANCE_FROM_BREAKEND - 1)
        };
        let fragments: Vec<Fragment> = candidates.into_iter().filter_map(|gi| break_end_fragment(gi, groups[gi], key, &mut cache)).collect();
        supports.push(Support { key: k, fragments, longest_telomere_segment: None, distance_to_telomere: None, longest_split_read_align_length: 0 });
    }
    for s in supports.iter_mut() {
        let key = &keys.0[s.key];
        // populateLongestTelomereSegment: the clips and the paired reads of discordant pairs
        let g_rich = key.ty.is_g_telomeric();
        let mut longest: Option<usize> = None;
        for f in &s.fragments {
            let clip = if key.ty.is_right_telomeric() { Part::RightClip } else { Part::LeftClip };
            let clipped = if clip == Part::RightClip { f.aligned_read.right_soft_clip() > 0 } else { f.aligned_read.left_soft_clip() > 0 };
            if clipped {
                longest = longest.max(cache.segment(f.aligned_read, clip, g_rich));
            }
            if f.paired_type != P::NotDiscordant {
                // the paired read as it is stored, not as classifyPairedRead orients it (TE5: kept, its intent
                // is not clear; oriented_part would give classifyPairedRead's view)
                longest = longest.max(cache.segment(f.paired_read, Part::Stored, g_rich));
            }
        }
        s.longest_telomere_segment = longest;
        // populateLongestSplitReadAlign
        for f in &s.fragments {
            if f.aligned_type.is_split_read() {
                let r = f.aligned_read;
                let len = r.read_length() as i32 - r.left_soft_clip() as i32 - r.right_soft_clip() as i32;
                s.longest_split_read_align_length = s.longest_split_read_align_length.max(len);
            }
        }
        // populateDistanceFromTelomere
        match chromosome_length(&key.chromosome, params.v37) {
            Some(length) => s.distance_to_telomere = Some(key.position.min(length - key.position)),
            None => crate::debug!("cannot find chromosome length for chromosome: {}", key.chromosome),
        }
    }
    supports
}

/// TelbamBreakEndFinder.run: the break ends of one telbam with at least one fragment
fn telbam_break_ends<'a>(params: &BreakEndParams, groups: &'a TelbamGroups, group_refs: &'a [&'a ReadGroup], keys: &mut Keys, extra: &[KeyId]) -> Vec<Support<'a>> {
    let t = std::time::Instant::now();
    let potential = candidate_break_ends(params, groups, keys, extra);
    crate::debug!("PerfStats: {} candidate break ends {:.1}s", potential.len(), t.elapsed().as_secs_f64());
    let t = std::time::Instant::now();
    let supports = count_break_end_supports(params, &potential, keys, group_refs);
    crate::debug!("PerfStats: break end supports {:.1}s", t.elapsed().as_secs_f64());
    supports.into_iter().filter(|s| !s.fragments.is_empty()).collect()
}

/// TelomericBreakEndEvidence (indices into the tumor and germline supports)
#[derive(Clone, Copy, Debug)]
struct Evidence {
    tumor: Option<usize>,
    germline: Option<usize>,
}

/// BreakEndEvidenceFilter.getFilterReasons
fn filter_reasons(tumor: Option<&Support>, germline: Option<&Support>, keys: &Keys) -> Vec<&'static str> {
    let mut f = Vec::new();
    match tumor {
        Some(t) => {
            if keys.0[t.key].is_duplicate {
                f.push("duplicate");
            }
            if t.longest_telomere_segment.is_none_or(|l| l < 20) {
                f.push("minTelomericLength");
            }
            if t.longest_split_read_align_length < 50 {
                f.push("minAnchorLength");
            }
            if t.fragment_count(|x| SUPPORT_SPLIT_READ_TYPES.contains(&x.types())) < 3 {
                f.push("minSplitReads");
            }
            if t.fragment_count(|x| SUPPORT_DISCORDANT_PAIR_TYPES.contains(&x.types())) < 1 {
                f.push("minDiscordantPairs");
            }
            if t.sum_mapq(|x| SUPPORT_FRAG_TYPES.contains(&x.types())) < 300 {
                f.push("minTumorMAPQ");
            }
        }
        None => f.push("notInTumor"),
    }
    if let Some(g) = germline {
        let tumor_count = tumor.map(|t| t.support_fragment_count()).unwrap_or(0);
        if g.support_fragment_count() as f64 > 5f64.min(tumor_count as f64 * 0.02) {
            f.push("maxGermlineSupport");
        }
    }
    if f.is_empty() {
        f.push("PASS");
    }
    f
}

/// Kotlin data class equality of two TelomericBreakEndSupport objects of one telbam
fn supports_equal(a: &Support, b: &Support, keys: &Keys) -> bool {
    keys.0[a.key] == keys.0[b.key]
        && a.longest_telomere_segment == b.longest_telomere_segment
        && a.distance_to_telomere == b.distance_to_telomere
        && a.longest_split_read_align_length == b.longest_split_read_align_length
        && a.fragments.len() == b.fragments.len()
        && a.fragments.iter().zip(&b.fragments).all(|(x, y)| {
            x.group == y.group && x.types() == y.types() && std::ptr::eq(x.aligned_read, y.aligned_read) && std::ptr::eq(x.paired_read, y.paired_read)
        })
}

/// SampleBreakEndFinder.run
pub fn find_break_ends(params: &BreakEndParams) -> Result<()> {
    let included = params.included_regions.as_deref();
    let tumor_groups = read_telbam(&params.tumor_telbam, &params.excluded_regions, included)?;
    let germline_groups = read_telbam(&params.germline_telbam, &params.excluded_regions, included)?;
    let tumor_refs: Vec<&ReadGroup> = tumor_groups.iter().collect();
    let germline_refs: Vec<&ReadGroup> = germline_groups.iter().collect();
    let mut keys = Keys::default();

    crate::info!("reading telbam file: {}", params.tumor_telbam);
    let tumor = telbam_break_ends(params, &tumor_groups, &tumor_refs, &mut keys, &[]);
    // the germline analysis looks at the tumor's break ends too: Java hands over the objects (TE4)
    let extra: Vec<KeyId> = if crate::compat::java("TE4") {
        tumor.iter().map(|s| s.key).collect()
    } else {
        tumor.iter().map(|s| { let k = keys.0[s.key].clone(); keys.add(k) }).collect()
    };
    crate::info!("reading telbam file: {}", params.germline_telbam);
    let germline = telbam_break_ends(params, &germline_groups, &germline_refs, &mut keys, &extra);

    // matchTumorGermlineBreakEnds
    let mut t_order: Vec<usize> = (0..tumor.len()).collect();
    let mut g_order: Vec<usize> = (0..germline.len()).collect();
    hmf_common::java_sort::java_sort(&mut t_order, |&a, &b| keys.compare(tumor[a].key, tumor[b].key));
    hmf_common::java_sort::java_sort(&mut g_order, |&a, &b| keys.compare(germline[a].key, germline[b].key));
    let merge = params.mark_duplicate_distance;
    let mut evidence: Vec<Evidence> = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < t_order.len() && j < g_order.len() {
        let (t, g) = (t_order[i], g_order[j]);
        i += 1;
        j += 1;
        let (tk, gk) = (&keys.0[tumor[t].key], &keys.0[germline[g].key]);
        if tk.ty == gk.ty && tk.chromosome == gk.chromosome && (tk.position - gk.position).abs() <= merge {
            evidence.push(Evidence { tumor: Some(t), germline: Some(g) });
        } else if keys.compare(tumor[t].key, germline[g].key) == std::cmp::Ordering::Greater {
            evidence.push(Evidence { tumor: None, germline: Some(g) });
            i -= 1;
        } else {
            evidence.push(Evidence { tumor: Some(t), germline: None });
            j -= 1;
        }
    }
    for &g in &g_order[j..] {
        evidence.push(Evidence { tumor: None, germline: Some(g) });
    }
    for &t in &t_order[i..] {
        evidence.push(Evidence { tumor: Some(t), germline: None });
    }
    let main_key = |e: &Evidence| -> KeyId { e.tumor.map(|t| tumor[t].key).or(e.germline.map(|g| germline[g].key)).unwrap() };
    // by chromosome rank (ContigComparator), then position
    hmf_common::java_sort::java_sort(&mut evidence, |a, b| {
        let (ka, kb) = (&keys.0[main_key(a)], &keys.0[main_key(b)]);
        hmf_common::regions::chromosome_rank(&ka.chromosome).cmp(&hmf_common::regions::chromosome_rank(&kb.chromosome)).then(ka.position.cmp(&kb.position))
    });

    // markDuplicateBreakEnds
    {
        // (in tumor, support index)
        let mut dups: Vec<(bool, usize)> = Vec::new();
        let mut current: Option<KeyId> = None;
        let support = |c: &(bool, usize)| -> &Support { if c.0 { &tumor[c.1] } else { &germline[c.1] } };
        let handle = |dups: &[(bool, usize)], keys: &mut Keys| {
            if dups.len() <= 1 {
                return;
            }
            let in_tumor: Vec<&(bool, usize)> = dups.iter().filter(|c| c.0).collect();
            let main_candidates: Vec<&(bool, usize)> = if in_tumor.is_empty() { dups.iter().collect() } else { in_tumor };
            // maxByOrNull: the first of the highest
            let mut main = main_candidates[0];
            for c in &main_candidates[1..] {
                if support(c).num_split_reads() > support(main).num_split_reads() {
                    main = c;
                }
            }
            let main = *main;
            // one by one, as Java marks them: a marked key can change a later comparison
            for c in dups {
                let equal = c.0 == main.0 && (c.1 == main.1 || supports_equal(support(c), support(&main), keys));
                if !equal {
                    keys.0[support(c).key].is_duplicate = true;
                }
            }
        };
        for e in &evidence {
            let main = if let Some(t) = e.tumor { (true, t) } else { (false, e.germline.unwrap()) };
            let mk = support(&main).key;
            let continues = current.is_some_and(|c| {
                let (ck, ek) = (&keys.0[c], &keys.0[mk]);
                ck.ty == ek.ty && ck.chromosome == ek.chromosome && (ck.position - ek.position).abs() <= params.mark_duplicate_distance
            });
            if continues {
                dups.push(main);
            } else {
                handle(&dups, &mut keys);
                dups.clear();
                dups.push(main);
                current = Some(mk);
            }
        }
        handle(&dups, &mut keys);
    }

    write_break_ends(&params.output_file, &evidence, &tumor, &germline, &keys)
}

fn column_name(a: AlignedReadType, p: PairedReadType, prefix: &str) -> String {
    if a == A::DiscordantPair && p == P::DiscordantPairTelomeric {
        return format!("{prefix}DPTelNoSR");
    }
    let at = match a {
        A::SplitReadTelomeric => "SRTel",
        A::NotSplitRead => "NotSplit",
        A::SplitReadNotTelomeric => "SRNotTel",
        A::DiscordantPair => "DP",
    };
    let pt = match p {
        P::DiscordantPairTelomeric => "DPTel",
        P::DiscordantPairNotTelomeric => "DPNotTel",
        P::NotDiscordant => "NoDP",
    };
    format!("{prefix}{at}{pt}")
}

/// BreakEndTsvWriter.writeBreakEnds: a tablesaw table written as gzipped TSV
fn write_break_ends(path: &str, evidence: &[Evidence], tumor: &[Support], germline: &[Support], keys: &Keys) -> Result<()> {
    let mut cols: Vec<String> = [
        "chromosome",
        "position",
        "orientation",
        "cOrGRich",
        "filter",
        "inTumor",
        "inGermline",
        "distanceToTelomere",
        "maxTelomericLength",
        "maxAnchorLength",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for prefix in ["tumor", "germline"] {
        for (a, p) in SUPPORT_FRAG_TYPES {
            cols.push(column_name(a, p, prefix));
        }
        cols.push(format!("{prefix}TotalSupport"));
        cols.push(format!("{prefix}MAPQ"));
    }
    let mut text = cols.join("\t");
    text.push('\n');
    let data = |s: Option<&Support>, v: &mut Vec<String>| {
        for (a, p) in SUPPORT_FRAG_TYPES {
            v.push(s.map(|x| x.count_of(a, p)).unwrap_or(0).to_string());
        }
        v.push(s.map(|x| x.support_fragment_count()).unwrap_or(0).to_string());
        v.push(s.map(|x| x.sum_mapq(|f| SUPPORT_FRAG_TYPES.contains(&f.types()))).unwrap_or(0).to_string());
    };
    for e in evidence {
        let t = e.tumor.map(|i| &tumor[i]);
        let g = e.germline.map(|i| &germline[i]);
        let main = t.or(g).unwrap();
        let key = &keys.0[main.key];
        let mut v: Vec<String> = vec![
            key.chromosome.to_string(),
            key.position.to_string(),
            (if key.ty.is_right_telomeric() { "1" } else { "-1" }).to_string(),
            (if matches!(key.ty, BreakEndType::RightGTelomeric | BreakEndType::LeftCTelomeric) { "G" } else { "C" }).to_string(),
            filter_reasons(t, g, keys).join(","),
            t.is_some().to_string(),
            g.is_some().to_string(),
            main.distance_to_telomere.map(|d| d.to_string()).unwrap_or_default(),
            t.and_then(|x| x.longest_telomere_segment).unwrap_or(0).to_string(),
            t.map(|x| x.longest_split_read_align_length).unwrap_or(0).to_string(),
        ];
        data(t, &mut v);
        data(g, &mut v);
        text.push_str(&v.join("\t"));
        text.push('\n');
    }
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let mut gz = flate2::write::GzEncoder::new(std::io::BufWriter::new(file), flate2::Compression::default());
    gz.write_all(text.as_bytes()).with_context(|| format!("writing {path}"))?;
    gz.finish().with_context(|| format!("writing {path}"))?.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::{FLAG_FIRST, FLAG_PAIRED, FLAG_REVERSE};

    fn params() -> BreakEndParams {
        let mut p = BreakEndParams::new(true);
        p.mark_duplicate_distance = 10;
        p.aligned_segment_telomere_match_threshold = 0.8;
        p.max_aligned_poly_base_threshold = 0.9;
        p
    }

    /// TelomericSplitReadAnalyserTest.testFindTelomericBreakEnd (boundary zone 2)
    #[test]
    fn find_telomeric_break_end() {
        let p = params();
        let base = Read::test("TestRead", "1", 1000, FLAG_PAIRED | FLAG_FIRST).with_mate("2", 200);
        let r = base.clone().with_cigar("20S30M").with_seq("TAACCCTAACCCTAACCCAGCTAGCTAGCTACGTTGCAACCTGACGAA");
        assert_eq!(find_candidate_break_end(&p, &r, 2), Some((BreakEndType::LeftCTelomeric, 1000)));
        let r = base.clone().with_cigar("30M20S").with_seq("AGCTAGCTAGCTACGTTGCAACCTGACGAATAACCCTAACCCTAACCC");
        assert_eq!(find_candidate_break_end(&p, &r, 2), Some((BreakEndType::RightCTelomeric, 1029)));
        let r = base.with_cigar("20S30M").with_seq("TTAGGGTTAGGGTTAGGGAGCTAGCTAGCTACGTTGCAACCTGACGAA");
        assert_eq!(find_candidate_break_end(&p, &r, 2), Some((BreakEndType::LeftGTelomeric, 1000)));
    }

    /// BreakEndSupportCounterTest.testBreakEndSupport1
    #[test]
    fn break_end_support() {
        let mut rg = ReadGroup::new(b"ReadGroup");
        rg.reads.push(
            Read::test("ReadGroup", "X", 100, FLAG_PAIRED | FLAG_FIRST | FLAG_REVERSE).with_mate("2", 200).with_seq("TTAGGGTTAGGGCTGAAGCTTGACCG").with_cigar("12S14M"),
        );
        assert!(!rg.is_complete().unwrap());
        rg.reads.push(Read::test("ReadGroup", "2", 200, FLAG_PAIRED).with_mate("X", 100).with_seq("TTAGGGTTAGGGTTA"));
        assert!(rg.is_complete().unwrap());
        assert!(rg.invariant());
        let mut keys = Keys::default();
        let k = keys.add(Key { ty: BreakEndType::LeftGTelomeric, chromosome: Arc::from("X"), position: 100, is_duplicate: false });
        let groups = vec![&rg];
        let s = count_break_end_supports(&params(), &[k], &keys, &groups);
        assert_eq!(s[0].fragments.len(), 1);
        let f = &s[0].fragments[0];
        assert!(std::ptr::eq(f.aligned_read, rg.first_of_pair().unwrap()));
        assert!(std::ptr::eq(f.paired_read, rg.second_of_pair().unwrap()));
        assert_eq!(f.aligned_type, A::SplitReadTelomeric);
        assert_eq!(f.paired_type, P::DiscordantPairTelomeric);
    }

    #[test]
    fn columns() {
        let names: Vec<String> = SUPPORT_FRAG_TYPES.iter().map(|&(a, p)| column_name(a, p, "tumor")).collect();
        assert_eq!(names, vec!["tumorSRTelDPTel", "tumorSRTelDPNotTel", "tumorSRTelNoDP", "tumorSRNotTelDPTel", "tumorDPTelNoSR"]);
    }
}
