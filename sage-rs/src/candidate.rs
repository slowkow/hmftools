//! Candidate stage: RefContextConsumer, EvictingArray, RefContext/AltContext, RegionBlock, Candidates.

use crate::config::Config;
use crate::constants::*;
use crate::matcher::*;
use crate::read_context::*;
use crate::refdata::RefData;
use crate::refseq::RefSequence;
use crate::select::*;
use crate::types::*;
use std::sync::Arc as Rc;

/// NumberEvents.calcAdjustedNumMutations
pub fn calc_adjusted_num_mutations(read: &Read, ref_seq: &RefSequence) -> i32 {
    let nm = read.nm.unwrap_or_else(|| calculate_nm(read, ref_seq));
    let mut additional = 0;
    for el in &read.cigar {
        if el.op.is_indel() {
            additional += el.len - 1;
        }
    }
    nm - additional
}

/// Fallback for a missing NM tag: mismatches in aligned blocks plus inserted and deleted bases.
fn calculate_nm(read: &Read, ref_seq: &RefSequence) -> i32 {
    let mut nm = 0;
    let mut read_index = 0usize;
    let mut ref_pos = read.start;
    for el in &read.cigar {
        match el.op {
            Op::M | Op::Eq | Op::X => {
                for i in 0..el.len {
                    let ri = ref_seq.index(ref_pos + i);
                    if ri >= 0 && (ri as usize) < ref_seq.bases.len() {
                        let r = ref_seq.bases[ri as usize].to_ascii_uppercase();
                        let b = read.bases[read_index + i as usize].to_ascii_uppercase();
                        if r != b {
                            nm += 1;
                        }
                    }
                }
                read_index += el.len as usize;
                ref_pos += el.len;
            }
            Op::I => {
                nm += el.len;
                read_index += el.len as usize;
            }
            Op::D => {
                nm += el.len;
                ref_pos += el.len;
            }
            Op::N => ref_pos += el.len,
            Op::S => read_index += el.len as usize,
            _ => {}
        }
    }
    nm
}

/// NumberEvents.calcSoftClipAdjustment
pub fn calc_soft_clip_adjustment(sc_len: i32) -> f64 {
    if sc_len > 0 {
        1f64.max(sc_len as f64 / SC_READ_EVENTS_FACTOR)
    } else {
        0.0
    }
}

/// NumberEvents.calcWithMnvRaw
pub fn calc_with_mnv_raw(events: i32, ref_: &[u8], alt: &[u8]) -> i32 {
    let diff = alt.iter().zip(ref_.iter()).filter(|(a, r)| a != r).count() as i32;
    events - diff + 1
}

/// QualityCalculator.calcEventPenalty
pub fn calc_event_penalty(num_events: f64, read_length: i32, penalty: f64) -> i32 {
    let nm_proportion = jround(0f64.max(num_events - 1.0)) as f64 / read_length as f64;
    jround(nm_proportion * penalty) as i32
}

pub struct ReadContextCandidate {
    pub ctx: Rc<VariantReadContext>,
    matcher: Matcher,
    pub full_match: i32,
    pub core_match: i32,
    pub min_events: i32,
}

struct AltContext {
    variant: SimpleVariant,
    candidates: Vec<ReadContextCandidate>,
    above_min_alt_support: bool,
    unique_read_ids: Vec<String>,
    has_sufficient_qual_core_support: bool,
    /// Ultima: mCheckCoreQualSupport, decided on the first read (None until then)
    ultima_check_core_qual_support: Option<bool>,
    /// Ultima: mHighQualCoreSupport, the reads without low-qual bases in the core
    high_qual_core_support: i32,
}

/// IlluminaUtils.lowQualInVariantCore
fn low_qual_in_variant_core(var_read_index: i32, ctx: &VariantReadContext, read: &Read) -> bool {
    let quals = &read.quals;
    if ctx.variant.is_indel() {
        if ctx.max_repeat.as_ref().map_or(false, |r| r.count >= 15) {
            return false;
        }
        let s = (var_read_index - ctx.left_core_length()).max(0);
        let e = (var_read_index + ctx.right_core_length()).min(quals.len() as i32 - 1);
        (s..=e).any(|i| is_low_base_qual(quals[i as usize]))
    } else {
        (var_read_index..=var_read_index + ctx.variant.alt_len() - 1).any(|i| is_low_base_qual(quals[i as usize]))
    }
}

impl AltContext {
    fn new(variant: SimpleVariant) -> Self {
        AltContext {
            variant,
            candidates: Vec::new(),
            above_min_alt_support: false,
            unique_read_ids: Vec::new(),
            // Illumina and Ultima check for high-quality core support
            has_sufficient_qual_core_support: !is_illumina()
                && !(is_ultima() && crate::ultima::candidate_min_high_bq_threshold() > 0),
            ultima_check_core_qual_support: if is_ultima() && crate::ultima::candidate_min_high_bq_threshold() <= 0 {
                Some(false)
            } else {
                None
            },
            high_qual_core_support: 0,
        }
    }

    fn check_qual_status(&mut self, ctx: &VariantReadContext, var_read_index: i32, read: &Read) {
        if self.has_sufficient_qual_core_support {
            return;
        }
        if !low_qual_in_variant_core(var_read_index, ctx, read) {
            self.has_sufficient_qual_core_support = true;
        }
    }

    fn add_read_context(
        &mut self, num_events: i32, read: &Read, var_read_index: i32, ref_seq: &RefSequence, min_alt_support: i32,
    ) {
        let mut full_match_count = 0;
        let mut has_full_match = false;
        let ultima = is_ultima();
        // Ultima: the read's low-qual indices, only while they are still needed (null in Java otherwise)
        let mut low_qual_indices: Option<Vec<i32>> = None;
        let mut ultima_check = false;
        if ultima {
            ultima_check = match self.ultima_check_core_qual_support {
                Some(check) => check,
                None => {
                    // whether this variant ought to check low-qual bases for this and subsequent reads.
                    // An out of range read index throws in Java, which ends the run.
                    let check = !crate::ultima::has_long_homopolymer_in_context(&self.variant, read, var_read_index)
                        .expect("variant read index within the read");
                    self.ultima_check_core_qual_support = Some(check);
                    self.has_sufficient_qual_core_support = !check;
                    check
                }
            };
            if ultima_check && !self.has_sufficient_qual_core_support {
                low_qual_indices = Some(crate::ultima::extract_low_qual_indices(read));
            }
        }
        // Illumina: mCheckCoreQualSupport flips to false exactly when sufficient support is found
        for i in 0..self.candidates.len() {
            let info = {
                let c = &self.candidates[i];
                c.matcher.determine(&c.ctx, &read.bases, None, var_read_index, true)
            };
            match info.kind {
                MatchType::Full => {
                    let c = &mut self.candidates[i];
                    c.full_match += 1;
                    c.min_events = c.min_events.min(num_events);
                    has_full_match = true;
                }
                MatchType::Core => {
                    let c = &mut self.candidates[i];
                    c.core_match += 1;
                    full_match_count += c.full_match;
                }
                _ => {}
            }
            if ultima {
                if ultima_check && matches!(info.kind, MatchType::Full | MatchType::Core) {
                    if let Some(low_qual_indices) = &low_qual_indices {
                        if !crate::ultima::low_qual_in_read_context_core(var_read_index, low_qual_indices, &self.candidates[i].ctx)
                        {
                            self.high_qual_core_support += 1;
                        }
                    }
                }
            } else if !self.has_sufficient_qual_core_support && matches!(info.kind, MatchType::Full | MatchType::Core) {
                let ctx = self.candidates[i].ctx.clone();
                self.check_qual_status(&ctx, var_read_index, read);
            }
        }

        if !has_full_match {
            // Ultima: skip creating a candidate if this variant checks low-quals in the core and finds them
            let skip_context_creation = ultima
                && ultima_check
                && crate::ultima::low_qual_in_read_core(&self.variant, var_read_index, low_qual_indices.as_deref());
            let mut new_ctx = None;
            if !skip_context_creation {
                let Some(ctx) = create_context(&self.variant, read, var_read_index, ref_seq, DEFAULT_FLANK_LENGTH) else {
                    return;
                };
                new_ctx = Some(ctx);
            }
            let mut has_valid_core = new_ctx.is_some();
            if let Some(ctx) = &new_ctx {
                // now the core has been established, check again
                if let Some(low_qual_indices) = &low_qual_indices {
                    if crate::ultima::low_qual_in_read_context_core(var_read_index, low_qual_indices, ctx) {
                        has_valid_core = false;
                    } else {
                        self.high_qual_core_support += 1;
                    }
                } else if !ultima {
                    self.check_qual_status(ctx, var_read_index, read);
                }
            }
            if let (true, Some(ctx)) = (has_valid_core, new_ctx) {
                let matcher = Matcher::new(&ctx, false, false);
                self.candidates.push(ReadContextCandidate {
                    ctx: Rc::new(ctx),
                    matcher,
                    full_match: 0,
                    core_match: full_match_count,
                    min_events: num_events,
                });
            }
        }

        if !self.above_min_alt_support && (self.unique_read_ids.len() as i32) < min_alt_support {
            if !self.unique_read_ids.contains(&read.name) {
                self.unique_read_ids.push(read.name.clone());
            }
            if self.unique_read_ids.len() as i32 >= min_alt_support {
                self.above_min_alt_support = true;
                self.unique_read_ids.clear();
            }
        }

        if ultima && !self.has_sufficient_qual_core_support {
            self.has_sufficient_qual_core_support =
                self.high_qual_core_support >= crate::ultima::candidate_min_high_bq_threshold();
        }
    }

    fn has_min_alt_support(&self) -> bool {
        self.above_min_alt_support && self.has_sufficient_qual_core_support
    }

    /// selectCandidates: returns (candidate, second candidate)
    fn select_candidates(mut self) -> (Option<ReadContextCandidate>, Option<ReadContextCandidate>) {
        if !self.has_sufficient_qual_core_support || self.candidates.is_empty() {
            return (None, None);
        }
        // stable, as Collections.sort
        self.candidates.sort_by(|a, b| {
            b.full_match.cmp(&a.full_match).then(b.core_match.cmp(&a.core_match)).then(a.min_events.cmp(&b.min_events))
        });
        let mut second_index = None;
        if self.candidates.len() > 1 {
            let top_core = self.candidates[0].ctx.core().to_vec();
            let threshold = self.candidates[0].full_match as f64 * MIN_SECOND_CANDIDATE_FULL_READS_PERC;
            for (i, c) in self.candidates.iter().enumerate() {
                if c.full_match < MIN_SECOND_CANDIDATE_FULL_READS || (c.full_match as f64) < threshold {
                    break;
                }
                let core = c.ctx.core();
                if java_index_of(core, &top_core, 0) >= 0 || java_index_of(&top_core, core, 0) >= 0 {
                    continue;
                }
                second_index = Some(i);
                break;
            }
        }
        let second = second_index.map(|i| self.candidates.remove(i));
        let first = self.candidates.swap_remove(0);
        (Some(first), second)
    }
}

pub(crate) struct RefContext {
    position: i32,
    alts: Vec<AltContext>,
}

/// EvictingArray: a ring buffer of RefContext by position.
pub(crate) struct EvictingArray {
    elements: Vec<Option<RefContext>>,
    min_position: i32,
    min_position_index: i32,
    read_length_buffer: i32,
    capacity: i32,
}

impl EvictingArray {
    fn new(max_read_length: i32) -> Self {
        let capacity = max_read_length * 2;
        let max_expected_del_perc = 80.0 / DEFAULT_READ_LENGTH as f64;
        let buffer = jround(max_read_length as f64 * (1.0 + max_expected_del_perc)) as i32;
        EvictingArray {
            elements: (0..capacity).map(|_| None).collect(),
            min_position: 0,
            min_position_index: 0,
            read_length_buffer: buffer,
            capacity,
        }
    }
    fn calc_index(&self, distance: i32) -> i32 {
        if self.min_position_index + distance < self.capacity {
            self.min_position_index + distance
        } else {
            distance + self.min_position_index - self.capacity
        }
    }
    fn reset_min_position(&mut self, position: i32) {
        self.min_position_index = 0;
        self.min_position = 1.max(position - self.read_length_buffer);
    }
    /// checkFlush; evicted contexts are appended to `evicted`.
    pub(crate) fn check_flush(&mut self, position: i32, evicted: &mut Vec<RefContext>) {
        let flush_count;
        if position > 0 {
            if self.min_position == 0 {
                self.reset_min_position(position);
                return;
            }
            let distance = position - self.min_position;
            if distance < self.capacity {
                return;
            }
            let new_min = position - self.read_length_buffer;
            flush_count = (new_min - self.min_position).min(self.capacity);
        } else {
            flush_count = self.capacity;
        }
        for _ in 0..flush_count.min(self.capacity) {
            if let Some(el) = self.elements[self.min_position_index as usize].take() {
                evicted.push(el);
            }
            self.min_position += 1;
            if self.min_position_index + 1 >= self.capacity {
                self.min_position_index = 0;
            } else {
                self.min_position_index += 1;
            }
        }
        if flush_count >= self.capacity {
            self.reset_min_position(position);
        }
    }
    /// getOrCreateRefContext: false when the position precedes the window (Java returns null).
    pub(crate) fn ensure(&mut self, position: i32, evicted: &mut Vec<RefContext>) -> bool {
        if self.min_position > 0 && position < self.min_position {
            eprintln!("WARN ignoring read with position({}) before prior position({})", position, self.min_position);
            return false;
        }
        self.check_flush(position, evicted);
        let index = self.calc_index(position - self.min_position);
        assert!(index >= 0 && index < self.capacity, "invalid evicting array index");
        let slot = &mut self.elements[index as usize];
        if slot.is_none() {
            *slot = Some(RefContext { position, alts: Vec::new() });
        }
        true
    }
    /// Look up a context created earlier while processing the same read; None if since evicted.
    fn get_mut(&mut self, position: i32) -> Option<&mut RefContext> {
        if position < self.min_position || position - self.min_position >= self.capacity {
            return None;
        }
        let index = self.calc_index(position - self.min_position);
        self.elements[index as usize].as_mut().filter(|c| c.position == position)
    }
}

/// Test-only constructor and accessors (EvictingArray's @VisibleForTesting methods in Java).
#[cfg(test)]
impl EvictingArray {
    pub(crate) fn with_capacity(capacity: i32, read_length_buffer: i32) -> Self {
        EvictingArray {
            elements: (0..capacity).map(|_| None).collect(),
            min_position: 0,
            min_position_index: 0,
            read_length_buffer,
            capacity,
        }
    }
    pub(crate) fn capacity(&self) -> i32 {
        self.capacity
    }
    pub(crate) fn read_length_buffer(&self) -> i32 {
        self.read_length_buffer
    }
    pub(crate) fn min_position(&self) -> i32 {
        self.min_position
    }
    pub(crate) fn item_count(&self) -> usize {
        self.elements.iter().filter(|x| x.is_some()).count()
    }
}

pub(crate) struct RegionBlock {
    pub(crate) start: i32,
    pub(crate) end: i32,
    hotspot_positions: Vec<i32>,
    pub(crate) apply_event_penalty: bool,
    pub(crate) depth_limit: i32,
    pub(crate) depth: i32,
}

impl RegionBlock {
    fn contains(&self, pos: i32) -> bool {
        pos >= self.start && pos <= self.end
    }
    pub(crate) fn limit_reached(&self) -> bool {
        self.depth >= self.depth_limit
    }
    fn is_hotspot(&self, pos: i32) -> bool {
        self.hotspot_positions.contains(&pos)
    }
    pub(crate) fn covers_hotspot(&self, s: i32, e: i32) -> bool {
        self.hotspot_positions.iter().any(|&p| position_within(p, s, e))
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_region_blocks(
    chr: &str, region: BaseRegion, panel: &[BaseRegion], hotspots: &[SimpleVariant], ref_data: &RefData,
    max_depth_panel: i32, max_depth_non_panel: i32,
) -> Vec<RegionBlock> {
    let is_mito = is_mito(chr);
    let block_size = REGION_BLOCK_SIZE;
    let mut selector = PanelSelector::new(panel);
    let mut blocks = Vec::new();
    let mut block_pos = region.start;
    while block_pos <= region.end {
        let next_start = block_pos + block_size;
        let block_end = if region.end - next_start < block_size / 2 { region.end } else { next_start - 1 };
        let hotspot_positions: Vec<i32> = {
            let mut v: Vec<i32> =
                hotspots.iter().filter(|h| position_within(h.pos, block_pos, block_end)).map(|h| h.pos).collect();
            v.dedup();
            v
        };
        let status = selector.panel_status(block_pos, block_end);
        let depth_limit =
            if !is_mito && status == PanelStatus::OutsidePanel { max_depth_non_panel } else { max_depth_panel };
        blocks.push(RegionBlock {
            start: block_pos,
            end: block_end,
            hotspot_positions,
            apply_event_penalty: !ref_data.hla_overlaps(chr, block_pos, block_end),
            depth_limit,
            depth: 0,
        });
        if block_end == region.end {
            break;
        }
        block_pos += block_size;
    }
    blocks
}

struct PendingAlt {
    position: i32,
    ref_: Vec<u8>,
    alt: Vec<u8>,
    num_events: i32,
    read_index: i32,
}

struct ReadInfo {
    num_events: i32,
    aligned_length: i32,
    soft_clip_length: i32,
    exceeds_quality: bool,
    exceeds_sc_adjusted_quality: bool,
    alts: Vec<PendingAlt>,
}

pub struct CandidateFinder<'a> {
    cfg: &'a Config,
    bounds: BaseRegion,
    ref_seq: &'a RefSequence,
    array: EvictingArray,
    blocks: Vec<RegionBlock>,
    current: usize,
    next: Option<usize>,
    pub read_count: i32,
    hotspot_cache: HotspotCache<'a>,
    saved: Vec<ReadContextCandidate>,
    evicted: Vec<RefContext>,
}

/// RefContextConsumer.processSoftClip (static): returns (ref, alt) of the implied insert.
pub fn soft_clip_insert(
    read_start: i32, read_end: i32, read_bases: &[u8], sc_length: i32, sc_read_index: i32, ref_seq: &RefSequence,
    on_left: bool,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let test_len = SC_INSERT_REF_TEST_LENGTH;
    if on_left {
        let prev_ref_pos = read_start - 1;
        let ref_index_start = prev_ref_pos - ref_seq.start - test_len + 1;
        let ref_index_end = ref_index_start + test_len;
        if ref_index_start < 0 || ref_index_end > ref_seq.bases.len() as i32 {
            return None;
        }
        let required = &ref_seq.bases[ref_index_start as usize..ref_index_end as usize];
        let sc_bases = &read_bases[..sc_length as usize];
        let max_index = sc_length - test_len - SC_INSERT_MIN_LENGTH;
        let mut match_index = -1;
        let mut next = java_index_of(sc_bases, required, 0);
        while next >= 0 {
            if next > max_index {
                break;
            }
            match_index = next;
            next = java_index_of(sc_bases, required, next + 1);
        }
        if match_index <= 0 {
            return None;
        }
        let implied = match_index + test_len - 1;
        let alt_length = sc_length - implied - 1;
        Some((
            read_bases[implied as usize..=implied as usize].to_vec(),
            read_bases[implied as usize..=(implied + alt_length) as usize].to_vec(),
        ))
    } else {
        let ref_index_start = read_end + 1 - ref_seq.start;
        let ref_index_end = ref_index_start + test_len;
        if ref_index_start < 0 || ref_index_end > ref_seq.bases.len() as i32 {
            return None;
        }
        let required = &ref_seq.bases[ref_index_start as usize..ref_index_end as usize];
        let sc_bases = &read_bases[sc_read_index as usize..];
        let match_index = java_index_of(sc_bases, required, 0);
        if match_index <= 0 || match_index < SC_INSERT_MIN_LENGTH {
            return None;
        }
        let implied = sc_read_index - 1;
        Some((
            read_bases[implied as usize..=implied as usize].to_vec(),
            read_bases[implied as usize..=(implied + match_index) as usize].to_vec(),
        ))
    }
}

/// CigarUtils.getPositionFromReadIndex: (refPosition, readIndexShift)
pub fn position_from_read_index(
    alignment_start: i32, cigar: &[CigarEl], read_index: i32, last_if_inserted: bool, infer_from_soft_clip: bool,
) -> Option<(i32, i32)> {
    if read_index < 0 {
        return None;
    }
    let mut ref_pos = alignment_start;
    let mut index = 0;
    for el in cigar {
        if el.op == Op::S && ref_pos == alignment_start && infer_from_soft_clip {
            ref_pos -= el.len;
        }
        let consumes_ref = el.op.consumes_ref() || (el.op == Op::S && infer_from_soft_clip);
        if el.op.consumes_read() && index + el.len - 1 >= read_index {
            let mut shift = 0;
            if consumes_ref {
                ref_pos += read_index - index;
            } else {
                if !last_if_inserted {
                    return None;
                }
                ref_pos -= 1;
                shift = index - read_index - 1;
            }
            if ref_pos <= 0 {
                return None;
            }
            return Some((ref_pos, shift));
        }
        if consumes_ref {
            ref_pos += el.len;
        }
        if el.op.consumes_read() {
            index += el.len;
        }
    }
    None
}

/// The index of the first position at which two equally long sequences differ.
pub(crate) fn first_mismatch(a: &[u8], b: &[u8]) -> Option<usize> {
    let mut offset = 0;
    let mut a_words = a.chunks_exact(8);
    let mut b_words = b.chunks_exact(8);
    for (x, y) in (&mut a_words).zip(&mut b_words) {
        let diff = u64::from_le_bytes(x.try_into().unwrap()) ^ u64::from_le_bytes(y.try_into().unwrap());
        if diff != 0 {
            return Some(offset + (diff.trailing_zeros() / 8) as usize);
        }
        offset += 8;
    }
    a_words.remainder().iter().zip(b_words.remainder()).position(|(x, y)| x != y).map(|i| offset + i)
}

pub(crate) fn mnv_length(read_index: i32, ref_index: i32, read_bases: &[u8], ref_bases: &[u8]) -> i32 {
    let diff = |i: i32| {
        ref_index + i < ref_bases.len() as i32
            && read_index + i < read_bases.len() as i32
            && ref_bases[(ref_index + i) as usize] != read_bases[(read_index + i) as usize]
    };
    if diff(2) {
        3
    } else if diff(1) {
        2
    } else {
        1
    }
}

impl<'a> CandidateFinder<'a> {
    pub fn new(
        cfg: &'a Config, chr: &str, bounds: BaseRegion, ref_seq: &'a RefSequence, hotspots: &'a [SimpleVariant],
        panel: &'a [BaseRegion], ref_data: &RefData,
    ) -> Self {
        let blocks =
            build_region_blocks(chr, bounds, panel, hotspots, ref_data, cfg.max_read_depth_panel, cfg.max_read_depth);
        CandidateFinder {
            cfg,
            bounds,
            ref_seq,
            array: EvictingArray::new(cfg.read_length.max(DEFAULT_READ_LENGTH)),
            blocks,
            current: 0,
            next: None,
            read_count: 0,
            hotspot_cache: HotspotCache::new(hotspots),
            saved: Vec::new(),
            evicted: Vec::new(),
        }
    }

    fn update_blocks(&mut self, read_start: i32) {
        if self.blocks[self.current].contains(read_start) {
            return;
        }
        while read_start > self.blocks[self.current].end && self.current < self.blocks.len() {
            self.current += 1;
        }
        if self.current + 1 < self.blocks.len() {
            self.next = Some(self.current + 1);
        }
    }
    fn next_block(&self) -> Option<&RegionBlock> {
        self.next.map(|i| &self.blocks[i])
    }
    fn reached_depth_limit(&self, pos: i32) -> bool {
        let cur = &self.blocks[self.current];
        if cur.contains(pos) {
            return cur.limit_reached();
        }
        match self.next_block() {
            Some(n) if n.contains(pos) => n.limit_reached(),
            _ => false,
        }
    }
    fn is_hotspot_position(&self, pos: i32) -> bool {
        self.blocks[self.current].is_hotspot(pos) || self.next_block().map_or(false, |n| n.is_hotspot(pos))
    }

    fn ref_context(&mut self, pos: i32) -> bool {
        let ok = self.array.ensure(pos, &mut self.evicted);
        self.process_evicted();
        ok
    }

    fn process_evicted(&mut self) {
        if self.evicted.is_empty() {
            return;
        }
        let evicted = std::mem::take(&mut self.evicted);
        for ctx in evicted {
            for alt in ctx.alts {
                // hasValidDnaBases checks the ref only
                if !alt.variant.ref_.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    continue;
                }
                if !(self.hotspot_cache.is_hotspot(&alt.variant) || alt.has_min_alt_support()) {
                    continue;
                }
                let (first, second) = alt.select_candidates();
                self.saved.extend(first);
                self.saved.extend(second);
            }
        }
    }

    pub fn process_read(&mut self, read: &Read) {
        let read_start = read.start;
        let read_end = read.end;
        if !positions_overlap(read_start, read_end, self.bounds.start, self.bounds.end) {
            return;
        }
        self.update_blocks(read_start);
        self.read_count += 1;
        if self.reached_depth_limit(read_start) && self.reached_depth_limit(read_end) {
            return;
        }

        let mut info = ReadInfo {
            num_events: 0,
            aligned_length: 0,
            soft_clip_length: 0,
            exceeds_quality: false,
            exceeds_sc_adjusted_quality: false,
            alts: Vec::new(),
        };
        for el in &read.cigar {
            match el.op {
                Op::M => info.aligned_length += el.len,
                Op::S => info.soft_clip_length += el.len,
                _ => {}
            }
        }
        if !read.supplementary() {
            info.num_events = calc_adjusted_num_mutations(read, self.ref_seq);
        }

        let sc_events = calc_soft_clip_adjustment(info.soft_clip_length) as i32;

        // applyMapQualEventPenalty
        let cur = &self.blocks[self.current];
        let mut apply_penalty = true;
        if cur.contains(read_start) && !cur.apply_event_penalty {
            apply_penalty = false;
        }
        if let Some(n) = self.next_block() {
            if n.contains(read_end) && !n.apply_event_penalty {
                apply_penalty = false;
            }
        }
        let adjusted_map_qual = if !apply_penalty {
            read.mapq
        } else {
            let event_penalty = calc_event_penalty(info.num_events as f64, read.len(), self.cfg.read_events_penalty);
            let improper =
                if read.improper_pair() || read.supplementary() { self.cfg.improper_pair_penalty } else { 0 };
            read.mapq - self.cfg.fixed_map_qual_penalty - event_penalty - improper
        };
        info.exceeds_quality = adjusted_map_qual > 0;

        let covers_hotspot = self.blocks[self.current].covers_hotspot(read_start, read_end)
            || self.next_block().map_or(false, |n| n.covers_hotspot(read_start, read_end));
        if !info.exceeds_quality && !covers_hotspot {
            return;
        }

        // updateRegionBlockDepth
        if self.blocks[self.current].contains(read_start) {
            self.blocks[self.current].depth += 1;
        }
        if let Some(n) = self.next {
            if self.blocks[n].contains(read_end) {
                self.blocks[n].depth += 1;
            }
        }

        let sc_proportion = sc_events as f64 / read.len() as f64;
        let sc_adjusted = jround(adjusted_map_qual as f64 - sc_proportion * self.cfg.read_events_penalty) as i32;
        info.exceeds_sc_adjusted_quality = sc_adjusted > 0;

        // ignoreSoftClipAdapter
        let ignore_sc_adapter = sc_events > 0 && {
            let frag_len = read.insert_size.abs();
            if frag_len >= read.len() + MIN_INSERT_ALIGNMENT_OVERLAP {
                false
            } else {
                (frag_len - info.aligned_length).abs() < MIN_INSERT_ALIGNMENT_OVERLAP
            }
        };

        // CigarHandler.traverseCigar
        let secondary = read.secondary_or_supplementary();
        let n = read.cigar.len();
        let mut read_index = 0;
        let mut ref_base = read.start;
        for (i, el) in read.cigar.iter().enumerate() {
            match el.op {
                Op::H | Op::P => {}
                Op::S => {
                    if !ignore_sc_adapter {
                        if i == 0 {
                            self.process_soft_clip(read, el.len, 0, &mut info, true);
                        } else if i == n - 1 {
                            self.process_soft_clip(read, el.len, read_index, &mut info, false);
                        }
                    }
                    read_index += el.len;
                }
                Op::N => ref_base += el.len,
                Op::D => {
                    if !secondary {
                        self.process_del(el.len, read, read_index - 1, ref_base - 1, &mut info);
                    }
                    ref_base += el.len;
                }
                Op::I => {
                    if !secondary {
                        let ri = if read_index == 0 { 0 } else { read_index - 1 };
                        self.process_insert(el.len, read, ri, ref_base - 1, &mut info);
                    }
                    read_index += el.len;
                }
                Op::M | Op::Eq | Op::X => {
                    if !secondary && (info.exceeds_sc_adjusted_quality || covers_hotspot) {
                        self.process_alignment(read, read_index, ref_base, el.len, &mut info);
                    }
                    read_index += el.len;
                    ref_base += el.len;
                }
            }
        }

        let min_alt_support = self.cfg.hard_min_tumor_raw_alt_support;
        let ref_seq = self.ref_seq;
        for alt in std::mem::take(&mut info.alts) {
            let Some(ctx) = self.array.get_mut(alt.position) else { continue };
            let idx = match ctx.alts.iter().position(|a| a.variant.ref_ == alt.ref_ && a.variant.alt == alt.alt) {
                Some(i) => i,
                None => {
                    ctx.alts.push(AltContext::new(SimpleVariant { pos: alt.position, ref_: alt.ref_, alt: alt.alt }));
                    ctx.alts.len() - 1
                }
            };
            ctx.alts[idx].add_read_context(alt.num_events, read, alt.read_index, ref_seq, min_alt_support);
        }
    }

    fn process_insert(&mut self, len: i32, read: &Read, read_index: i32, ref_position: i32, info: &mut ReadInfo) {
        if !self.bounds.contains(ref_position) {
            return;
        }
        let exceeds =
            if len as f64 <= SC_READ_EVENTS_FACTOR { info.exceeds_sc_adjusted_quality } else { info.exceeds_quality };
        if !exceeds && !self.is_hotspot_position(ref_position) {
            return;
        }
        let ref_index = self.ref_seq.index(ref_position) as usize;
        let ref_ = self.ref_seq.bases[ref_index..ref_index + 1].to_vec();
        let alt = read.bases[read_index as usize..(read_index + len + 1) as usize].to_vec();
        if alt.contains(&b'N') {
            return;
        }
        if self.ref_context(ref_position) {
            info.alts.push(PendingAlt { position: ref_position, ref_, alt, num_events: info.num_events, read_index });
        }
    }

    fn process_del(&mut self, len: i32, read: &Read, read_index: i32, ref_position: i32, info: &mut ReadInfo) {
        if !self.bounds.contains(ref_position) {
            return;
        }
        let exceeds =
            if len as f64 <= SC_READ_EVENTS_FACTOR { info.exceeds_sc_adjusted_quality } else { info.exceeds_quality };
        if !exceeds && !self.is_hotspot_position(ref_position) {
            return;
        }
        let ref_index = self.ref_seq.index(ref_position) as usize;
        let ref_ = self.ref_seq.bases[ref_index..ref_index + len as usize + 1].to_vec();
        let alt = read.bases[read_index as usize..read_index as usize + 1].to_vec();
        if self.ref_context(ref_position) {
            info.alts.push(PendingAlt { position: ref_position, ref_, alt, num_events: info.num_events, read_index });
        }
    }

    fn process_alignment(&mut self, read: &Read, read_start_index: i32, ref_pos_start: i32, len: i32, info: &mut ReadInfo) {
        let ref_index = self.ref_seq.index(ref_pos_start);
        let ref_bases = &self.ref_seq.bases;
        // only the part of the element within the region bounds is looked at
        let first = (self.bounds.start - ref_pos_start).max(0);
        let end = (self.bounds.end - ref_pos_start + 1).min(len);
        let mut i = first - 1;
        loop {
            i += 1;
            if i >= end {
                break;
            }
            if info.exceeds_quality {
                // nothing happens at a base matching the reference, so go straight to the next mismatch
                let ref_part = &ref_bases[(ref_index + i) as usize..(ref_index + end) as usize];
                let read_part = &read.bases[(read_start_index + i) as usize..(read_start_index + end) as usize];
                match first_mismatch(ref_part, read_part) {
                    Some(offset) => i += offset as i32,
                    None => break,
                }
            }
            let ref_position = ref_pos_start + i;
            let read_base_index = read_start_index + i;
            let ref_base_index = ref_index + i;
            if !info.exceeds_quality && !self.is_hotspot_position(ref_position) {
                continue;
            }
            let ref_byte = ref_bases[ref_base_index as usize];
            let read_byte = read.bases[read_base_index as usize];
            if read_byte == ref_byte {
                continue;
            }
            if !self.ref_context(ref_position) {
                continue;
            }
            info.alts.push(PendingAlt {
                position: ref_position,
                ref_: vec![ref_byte],
                alt: vec![read_byte],
                num_events: info.num_events,
                read_index: read_base_index,
            });
            let max_len = mnv_length(read_base_index, ref_base_index, &read.bases, ref_bases);
            let mut next_read_index = i;
            for mnv_len in 2..=max_len {
                next_read_index += 1;
                if next_read_index >= len {
                    break;
                }
                let mnv_ref = &ref_bases[ref_base_index as usize..(ref_base_index + mnv_len) as usize];
                let mnv_alt = &read.bases[read_base_index as usize..(read_base_index + mnv_len) as usize];
                if mnv_ref[mnv_len as usize - 1] != mnv_alt[mnv_len as usize - 1] {
                    info.alts.push(PendingAlt {
                        position: ref_position,
                        ref_: mnv_ref.to_vec(),
                        alt: mnv_alt.to_vec(),
                        num_events: calc_with_mnv_raw(info.num_events, mnv_ref, mnv_alt),
                        read_index: read_base_index,
                    });
                }
            }
        }
    }

    fn process_soft_clip(&mut self, read: &Read, sc_length: i32, sc_read_index: i32, info: &mut ReadInfo, on_left: bool) {
        if !info.exceeds_quality {
            return;
        }
        if sc_length < SC_INSERT_REF_TEST_LENGTH + 1 {
            return;
        }
        let Some((alt_ref, alt_alt)) =
            soft_clip_insert(read.start, read.end, &read.bases, sc_length, sc_read_index, self.ref_seq, on_left)
        else {
            return;
        };
        let (mut ref_position, mut read_index) = if on_left {
            (read.start - 1, sc_length - alt_alt.len() as i32)
        } else {
            (read.end, read.len() - sc_length - 1)
        };
        if !self.bounds.contains(ref_position) {
            return;
        }
        if !(read_index >= DEFAULT_FLANK_LENGTH && read_index < read.len() - DEFAULT_FLANK_LENGTH) {
            return;
        }
        let mut variant = SimpleVariant { pos: ref_position, ref_: alt_ref, alt: alt_alt };
        if variant.is_insert() && !on_left {
            let shift = find_left_homology_shift(&variant, self.ref_seq, &read.bases, read_index);
            if shift > 0 {
                let new_read_index = read_index - shift;
                let Some((pos, idx_shift)) = position_from_read_index(read.start, &read.cigar, new_read_index, true, true)
                else {
                    return;
                };
                ref_position = pos;
                read_index = new_read_index + idx_shift;
                let Some(new_ref) = self.ref_seq.position_bases(ref_position, ref_position) else { return };
                let mut new_alt = new_ref.to_vec();
                new_alt.extend_from_slice(&read.bases[(read_index + 1) as usize..=(read_index + variant.alt_len() - 1) as usize]);
                variant = SimpleVariant { pos: ref_position, ref_: new_ref.to_vec(), alt: new_alt };
            }
        }
        if self.ref_context(variant.pos) {
            info.alts.push(PendingAlt {
                position: variant.pos,
                ref_: variant.ref_,
                alt: variant.alt,
                num_events: info.num_events,
                read_index,
            });
        }
    }

    /// RefContextCache.altCandidates: flush everything and sort by support.
    pub fn finish(mut self) -> Vec<ReadContextCandidate> {
        self.array.check_flush(-1, &mut self.evicted);
        self.process_evicted();
        self.saved.sort_by(|a, b| {
            b.full_match.cmp(&a.full_match).then(b.core_match.cmp(&a.core_match)).then(a.min_events.cmp(&b.min_events))
        });
        self.saved
    }
}

#[derive(Clone)]
pub struct Candidate {
    pub tier: Tier,
    pub ctx: Rc<VariantReadContext>,
    pub min_events: i32,
    pub full_match_support: i32,
    pub core_match_support: i32,
}

impl Candidate {
    pub fn variant(&self) -> &SimpleVariant {
        &self.ctx.variant
    }
    pub fn position(&self) -> i32 {
        self.ctx.variant.pos
    }
}

/// Candidates.addSingleSample + candidates()
pub fn build_candidates(
    alt_candidates: Vec<ReadContextCandidate>, hotspots: &[SimpleVariant], panel: &[BaseRegion], high_conf: &[BaseRegion],
) -> Vec<Candidate> {
    let mut tiers = TierSelector::new(hotspots, panel, high_conf);
    let mut list: Vec<Candidate> = Vec::new();
    for alt in alt_candidates {
        let tier = tiers.tier(&alt.ctx.variant);
        let candidate = Candidate {
            tier,
            ctx: alt.ctx,
            min_events: alt.min_events,
            full_match_support: alt.full_match,
            core_match_support: alt.core_match,
        };
        let mut index = 0;
        while index < list.len() {
            if list[index].variant().compare_full(candidate.variant()) == std::cmp::Ordering::Greater {
                break;
            }
            index += 1;
        }
        list.insert(index, candidate);
    }
    list
}

/// Candidates.addOfMultipleSamples + candidates(): the candidates of several tumor samples. Those of a
/// variant which share a core are merged, keeping the read context with the most full matches.
pub fn build_candidates_of_multiple_samples(
    sample_alt_candidates: Vec<Vec<ReadContextCandidate>>, hotspots: &[SimpleVariant], panel: &[BaseRegion],
    high_conf: &[BaseRegion],
) -> Vec<Candidate> {
    // Java keeps a HashMap of variant -> candidates and sorts its values by variant afterwards; the sort
    // is stable and only candidates of one variant compare equal, so the map's order does not matter
    let mut index_of: std::collections::HashMap<SimpleVariant, usize> = std::collections::HashMap::new();
    let mut by_variant: Vec<Vec<Candidate>> = Vec::new();
    for alt_candidates in sample_alt_candidates {
        let mut tiers = TierSelector::new(hotspots, panel, high_conf);
        for alt in alt_candidates {
            let index = *index_of.entry(alt.ctx.variant.clone()).or_insert_with(|| {
                by_variant.push(Vec::new());
                by_variant.len() - 1
            });
            let candidates = &mut by_variant[index];
            let tier = tiers.tier(&alt.ctx.variant);
            match candidates.iter_mut().find(|c| c.ctx.core_str() == alt.ctx.core_str()) {
                Some(matching) => {
                    // Candidate.updateAltCandidate
                    if alt.full_match > matching.full_match_support {
                        matching.full_match_support = alt.full_match;
                        matching.ctx = alt.ctx;
                        matching.min_events = matching.min_events.min(alt.min_events);
                    }
                }
                None => candidates.push(Candidate {
                    tier,
                    ctx: alt.ctx,
                    min_events: alt.min_events,
                    full_match_support: alt.full_match,
                    core_match_support: alt.core_match,
                }),
            }
        }
    }
    let mut list: Vec<Candidate> = by_variant.into_iter().flatten().collect();
    list.sort_by(|a, b| a.variant().compare_full(b.variant()));
    list
}
