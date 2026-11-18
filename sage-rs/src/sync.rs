//! Fragment sync: overlapping mates are merged into a single consensus read (FragmentSync,
//! CombinedSyncData, CigarState).

use crate::constants::*;
use crate::evidence::FragmentData;
use crate::types::*;
use std::collections::HashMap;

pub const SYNC_COMBINED: usize = 0;
pub const SYNC_INVERSION: usize = 1;
pub const SYNC_NO_OVERLAP: usize = 2;
pub const SYNC_CIGAR_MISMATCH: usize = 3;
pub const SYNC_EXCEPTION: usize = 5;

pub(crate) enum Outcome {
    Combined(Read),
    NoOverlap,
    Inversion,
    CigarMismatch,
    /// stands in for any exception thrown by the Java code, which drops both reads
    Exception,
}

struct CigarState<'a> {
    elements: &'a [CigarEl],
    index: usize,
    base_index: i32,
    exhausted: bool,
}

impl<'a> CigarState<'a> {
    fn new(elements: &'a [CigarEl]) -> Self {
        CigarState { elements, index: 0, base_index: 0, exhausted: false }
    }
    fn current(&self) -> Option<CigarEl> {
        if self.exhausted {
            None
        } else {
            Some(self.elements[self.index])
        }
    }
    fn op(&self) -> Option<Op> {
        self.current().map(|e| e.op)
    }
    fn move_next(&mut self) -> bool {
        let Some(cur) = self.current() else { return false };
        if self.base_index < cur.len - 1 {
            self.base_index += 1;
            return false;
        }
        if self.index == self.elements.len() - 1 {
            self.exhausted = true;
        } else {
            self.index += 1;
            self.base_index = 0;
        }
        true
    }
    /// Err when Java would index past the last element
    fn move_to_position(&mut self, effective_start: i32, required_start: i32) -> Result<(), ()> {
        let mut diff = required_start - effective_start;
        while diff > 0 {
            let cur = self.elements[self.index];
            let min_bases = cur.len.min(diff);
            if cur.op.consumes_ref() || cur.op == Op::S {
                diff -= min_bases;
            }
            if min_bases < cur.len {
                self.base_index += min_bases;
                break;
            } else {
                self.index += 1;
                if self.index >= self.elements.len() {
                    return Err(());
                }
                self.base_index = 0;
            }
            if diff <= 0 {
                break;
            }
        }
        Ok(())
    }
}

fn soft_clip_start(c: &[CigarEl]) -> i32 {
    if c[0].op == Op::S {
        c[0].len
    } else {
        0
    }
}
fn soft_clip_end(c: &[CigarEl]) -> i32 {
    let l = c[c.len() - 1];
    if l.op == Op::S {
        l.len
    } else {
        0
    }
}

fn moves_effective_position(op: Op) -> bool {
    op.consumes_ref() || op == Op::S
}

/// BaseQualAdjustment.adjustBaseQual with the standard quals { 1, 11, 25, 37 }
fn adjust_base_qual(q: f64) -> u8 {
    const STD: [i32; 4] = [1, 11, 25, 37];
    if q < STD[1] as f64 - 1.5 {
        return STD[0] as u8;
    }
    for i in (1..STD.len()).rev() {
        if q > STD[i] as f64 {
            return STD[i] as u8;
        }
        if STD[i] as f64 - q <= 1.5 {
            return STD[i] as u8;
        }
    }
    STD[0] as u8
}

pub(crate) fn form_fragment_read(first: &Read, second: &Read) -> Outcome {
    if !positions_overlap(first.start, first.end, second.start, second.end) {
        return Outcome::NoOverlap;
    }
    if first.negative_strand() == second.negative_strand() {
        return Outcome::Inversion;
    }
    if first.cigar.is_empty() || second.cigar.is_empty() {
        return Outcome::Exception;
    }

    // buildCombinedCigar
    let first_sc_start = soft_clip_start(&first.cigar);
    let first_sc_end = soft_clip_end(&first.cigar);
    let second_sc_start = soft_clip_start(&second.cigar);
    let second_sc_end = soft_clip_end(&second.cigar);
    let first_eff_start = first.start - first_sc_start;
    let second_eff_start = second.start - second_sc_start;
    let first_eff_end = first.end + first_sc_end;
    let second_eff_end = second.end + second_sc_end;
    let (fragment_start, fragment_end, combined_eff_start, combined_eff_end) = if !first.negative_strand() {
        (first.start, second.end, first_eff_start, second_eff_end)
    } else {
        (second.start, first.end, second_eff_start, first_eff_end)
    };

    let mut first_cigar = CigarState::new(&first.cigar);
    let mut second_cigar = CigarState::new(&second.cigar);
    if first_eff_start < combined_eff_start {
        if first_cigar.move_to_position(first_eff_start, combined_eff_start).is_err() {
            return Outcome::Exception;
        }
    } else if second_eff_start < combined_eff_start
        && second_cigar.move_to_position(second_eff_start, combined_eff_start).is_err()
    {
        return Outcome::Exception;
    }

    let mut combined: Vec<CigarEl> = Vec::new();
    let mut el_len = 0;
    let mut el_op: Option<Op> = None;
    let mut first_active = false;
    let mut second_active = false;
    let mut current_pos = combined_eff_start;
    while current_pos <= combined_eff_end {
        let mut first_change = false;
        let mut second_change = false;
        if current_pos >= first_eff_start {
            if !first_active {
                first_active = true;
                first_change = true;
            } else if !first_cigar.exhausted {
                first_change = first_cigar.move_next();
            }
        }
        if current_pos >= second_eff_start {
            if !second_active {
                second_active = true;
                second_change = true;
            } else if !second_cigar.exhausted {
                second_change = second_cigar.move_next();
            }
        }
        if first_change || second_change {
            let first_op = if first_active { first_cigar.op() } else { None };
            let second_op = if second_active { second_cigar.op() } else { None };
            let new_op = match (first_op, second_op) {
                (Some(f), Some(s)) => {
                    if f != s {
                        let ignorable = matches!(f, Op::M | Op::S) && matches!(s, Op::M | Op::S);
                        if !ignorable {
                            return Outcome::CigarMismatch;
                        }
                        if position_within(current_pos, fragment_start, fragment_end) {
                            Some(Op::M)
                        } else {
                            Some(Op::S)
                        }
                    } else {
                        Some(f)
                    }
                }
                (Some(f), None) => Some(f),
                (None, s) => s,
            };
            if el_len > 0 && el_op != new_op {
                match el_op {
                    Some(op) => combined.push(CigarEl::new(el_len, op)),
                    None => return Outcome::Exception,
                }
                el_len = 0;
            }
            el_op = new_op;
        }
        el_len += 1;
        if el_op == Some(Op::I) {
            current_pos -= 1;
        }
        current_pos += 1;
    }
    match el_op {
        Some(op) => combined.push(CigarEl::new(el_len, op)),
        None => return Outcome::Exception,
    }

    // setBaseAndQualities
    let combined_len: i32 = combined.iter().filter(|e| e.op.consumes_read()).map(|e| e.len).sum();
    let mut bases = vec![0u8; combined_len as usize];
    let mut quals = vec![0u8; combined_len as usize];
    let find_start = |read: &Read, effective_start: i32| -> i32 {
        if effective_start > combined_eff_start {
            return -1;
        }
        if effective_start == combined_eff_start {
            return 0;
        }
        let mut read_position = effective_start;
        let mut read_index = 0;
        for el in &read.cigar {
            let mut position_move = 0;
            if moves_effective_position(el.op) {
                position_move = (combined_eff_start - read_position).min(el.len);
                read_position += position_move;
            }
            if el.op.consumes_read() {
                read_index += if position_move == 0 { el.len } else { position_move };
            }
            if read_position >= combined_eff_start {
                break;
            }
        }
        read_index
    };
    let mut first_index = find_start(first, first_eff_start);
    let mut second_index = find_start(second, second_eff_start);
    let first_len = first.len();
    let second_len = second.len();
    let mut state = CigarState::new(&combined);
    let mut current_position = combined_eff_start;
    let mut ci = 0;
    while ci < combined_len {
        let Some(op) = state.op() else { return Outcome::Exception };
        let move_read = op.consumes_read();
        let move_ref = moves_effective_position(op);
        if move_read {
            let first_ok = first_index >= 0 && first_index < first_len;
            let second_ok = second_index >= 0 && second_index < second_len;
            if first_ok && second_ok {
                let (fi, si) = (first_index as usize, second_index as usize);
                let first_in_sc = first_index < first_sc_start || first_index >= first_len - first_sc_end;
                let second_in_sc = second_index < second_sc_start || second_index >= second_len - second_sc_end;
                let mut fq = first.quals[fi] as i32;
                let mut sq = second.quals[si] as i32;
                if first_in_sc && !second_in_sc {
                    fq = MATCHING_BASE_QUALITY as i32;
                } else if !first_in_sc && second_in_sc {
                    sq = MATCHING_BASE_QUALITY as i32;
                }
                if first.bases[fi] == second.bases[si] {
                    bases[ci as usize] = first.bases[fi];
                    quals[ci as usize] = fq.max(sq) as u8;
                } else {
                    // getCombinedBaseAndQual for differing bases
                    let (base, qual) = if fq > sq {
                        (first.bases[fi], (BASE_QUAL_MINIMUM as i32).max(fq - sq))
                    } else {
                        (second.bases[si], (BASE_QUAL_MINIMUM as i32).max(sq - fq))
                    };
                    bases[ci as usize] = base;
                    quals[ci as usize] = adjust_base_qual(qual as f64);
                }
            } else if first_ok {
                bases[ci as usize] = first.bases[first_index as usize];
                quals[ci as usize] = first.quals[first_index as usize];
            } else if second_ok {
                bases[ci as usize] = second.bases[second_index as usize];
                quals[ci as usize] = second.quals[second_index as usize];
            }
        }
        if move_ref {
            current_position += 1;
        }
        state.move_next();
        if move_read {
            ci += 1;
        }
        if first_index < 0 {
            if current_position == first_eff_start {
                first_index = 0;
            }
        } else if move_read {
            first_index += 1;
        }
        if second_index < 0 {
            if current_position == second_eff_start {
                second_index = 0;
            }
        } else if move_read {
            second_index += 1;
        }
    }

    // buildSyncedRead: flags and tags are taken from the first read
    let ref_len: i32 = combined.iter().filter(|e| e.op.consumes_ref()).map(|e| e.len).sum();
    Outcome::Combined(Read {
        name: first.name.clone(),
        flags: first.flags,
        tid: first.tid,
        start: fragment_start,
        end: fragment_start + ref_len - 1,
        mapq: first.mapq.max(second.mapq),
        cigar: combined,
        bases,
        quals,
        mate_tid: first.mate_tid,
        mate_start: fragment_start,
        insert_size: first.insert_size.abs(),
        nm: first.nm,
        mate_cigar: first.mate_cigar.clone(),
        consensus: first.consensus,
        has_consensus_tag: first.has_consensus_tag,
        duplex_index: first.duplex_index,
    })
}

fn java_string_hash(s: &str) -> i32 {
    let mut h: i32 = 0;
    for c in s.encode_utf16() {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

/// A multiply-rotate hasher for read names (as rustc's FxHasher). The cache is only ever looked up by
/// name, and is ordered explicitly when it is emptied, so the hash function does not affect output.
#[derive(Default, Clone, Copy)]
pub struct FastHasher(u64);

impl std::hash::Hasher for FastHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        const K: u64 = 0x517c_c1b7_2722_0a95;
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            let v = u64::from_le_bytes(c.try_into().unwrap());
            self.0 = (self.0.rotate_left(5) ^ v).wrapping_mul(K);
        }
        for b in chunks.remainder() {
            self.0 = (self.0.rotate_left(5) ^ *b as u64).wrapping_mul(K);
        }
    }
}

type FastHash = std::hash::BuildHasherDefault<FastHasher>;

pub struct FragmentSync {
    cached: HashMap<String, (Read, u64), FastHash>,
    seq: u64,
    /// the Java HashMap never shrinks, so its capacity follows the peak size over its lifetime
    peak: usize,
    pub counts: [i32; 6],
}

pub type SyncOutput = Vec<(Read, Option<FragmentData>)>;

impl FragmentSync {
    pub fn new() -> Self {
        FragmentSync { cached: HashMap::default(), seq: 0, peak: 0, counts: [0; 6] }
    }

    /// handleOverlappingReads: Some(read) when the caller should process the read itself; otherwise
    /// the read was consumed and any reads to process have been appended to `out`.
    pub fn handle(&mut self, read: Read, out: &mut SyncOutput) -> Option<Read> {
        if !read.paired() || read.flag(FLAG_MATE_UNMAPPED) {
            return Some(read);
        }
        if let Some((other, other_seq)) = self.cached.remove(&read.name) {
            match form_fragment_read(&other, &read) {
                Outcome::Combined(combined) => {
                    self.counts[SYNC_COMBINED] += 1;
                    out.push((combined, Some(FragmentData { first: other, second: read })));
                }
                Outcome::CigarMismatch => {
                    self.counts[SYNC_CIGAR_MISMATCH] += 1;
                    let indel_len = |r: &Read| -> i32 { r.cigar.iter().filter(|e| e.op.is_indel()).map(|e| e.len).sum() };
                    if indel_len(&read) > indel_len(&other) {
                        out.push((read, None));
                    } else {
                        out.push((other, None));
                    }
                }
                Outcome::NoOverlap => {
                    self.counts[SYNC_NO_OVERLAP] += 1;
                    out.push((other, None));
                    out.push((read, None));
                }
                Outcome::Inversion => {
                    self.counts[SYNC_INVERSION] += 1;
                    out.push((other, None));
                    out.push((read, None));
                }
                Outcome::Exception => {
                    // Java throws before removing the cached read, so it stays cached and the new read is dropped
                    self.counts[SYNC_EXCEPTION] += 1;
                    self.cached.insert(other.name.clone(), (other, other_seq));
                }
            }
            return None;
        }
        if read.tid != read.mate_tid {
            return Some(read);
        }
        if read.mate_start < read.start {
            return Some(read);
        }
        if !positions_overlap(read.start, read.end, read.mate_start, read.mate_start + read.len()) {
            self.counts[SYNC_NO_OVERLAP] += 1;
            return Some(read);
        }
        self.seq += 1;
        self.cached.insert(read.name.clone(), (read, self.seq));
        self.peak = self.peak.max(self.cached.len());
        None
    }

    /// emptyCachedReads: unpaired leftovers in Java HashMap iteration order.
    pub fn empty_cached_reads(&mut self) -> Vec<Read> {
        if self.cached.is_empty() {
            return Vec::new();
        }
        let mut capacity = 16usize;
        while self.peak as f64 > capacity as f64 * 0.75 {
            capacity *= 2;
        }
        let mut items: Vec<(usize, u64, Read)> = self
            .cached
            .drain()
            .map(|(name, (read, seq))| {
                let h = java_string_hash(&name);
                let spread = (h ^ ((h as u32) >> 16) as i32) as u32 as usize;
                (spread & (capacity - 1), seq, read)
            })
            .collect();
        items.sort_by_key(|(bucket, seq, _)| (*bucket, *seq));
        items.into_iter().map(|(_, _, r)| r).collect()
    }
}
