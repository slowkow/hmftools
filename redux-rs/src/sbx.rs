//! Roche SBX reads (Java consensus/SbxRoutines, SbxDuplexIndelBuilder, SbxDuplexIndel, SbxBamSampler and
//! hmf-common sequencing/SbxBamUtils): the changes to a read before duplicate marking, the consensus base rule
//! and the final base qualities of a read that is written as a non-duplicate.

use crate::consensus::{BaseQualPair, RefBases, INVALID_POSITION};
use crate::rec::*;
use crate::unmap::extract_supp_alignment;
use rust_htslib::bam;
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, Ordering};

pub const RAW_DUPLEX_QUAL: u8 = 39;
pub const RAW_SIMPLEX_QUAL: u8 = 22;
pub const RAW_DUPLEX_MISMATCH_QUAL: u8 = 5;

pub const SBX_SIMPLEX_QUAL: u8 = 27;
pub const SBX_DUPLEX_QUAL: u8 = 40;
pub const SBX_DUPLEX_MISMATCH_QUAL: u8 = 1;
pub const SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH: u8 = 15;
pub const SBX_DUPLEX_ADJACENT_1_QUAL: u8 = 20;
pub const SBX_DUPLEX_ADJACENT_2_3_QUAL: u8 = 25;

/// internal base qualities that keep a simplex and a duplex position without consensus apart until finaliseRead
pub const DUPLEX_NO_CONSENSUS_QUAL: u8 = 3;
pub const SIMPLEX_NO_CONSENSUS_QUAL: u8 = 2;

pub const SBX_READ_CACHE_GROUP_SIZE: i32 = 1000;
pub const SBX_READ_CACHE_MAX_SOFT_CLIP: i32 = SBX_READ_CACHE_GROUP_SIZE - 30;
pub const SBX_READ_CACHE_LOG_READ_COUNT_THRESHOLD: i32 = 5000;
pub const SBX_CONSENSUS_MAX_DEPTH: usize = 20;
pub const DEFAULT_SBX_MAX_DUPLICATE_DISTANCE: i32 = 2;
pub const SBX_CONSENSUS_BASE_THRESHOLD: f64 = 0.5;

const BWA_MATCH_SCORE: i64 = 1;
const BWA_MISMATCH_PENALTY: i64 = 4;
const BWA_GAP_OPEN_PENALTY: i64 = 6;
const BWA_GAP_EXTEND_PENALTY: i64 = 1;

/// SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL: set once from the sampled reads before the partitions start
static HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL: AtomicBool = AtomicBool::new(true);

pub fn homopolymer_5_prime_low_base_qual() -> bool {
    HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL.load(Ordering::Relaxed)
}
pub fn set_homopolymer_5_prime_low_base_qual(value: bool) {
    HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL.store(value, Ordering::Relaxed);
}

/// A Java exception inside the SBX routines: the callers log it and carry on with the read as it is.
#[derive(Debug)]
pub struct SbxError(pub String);

type SbxResult<T> = Result<T, SbxError>;

fn err<T>(msg: &str) -> SbxResult<T> {
    Err(SbxError(msg.to_string()))
}

// ---------------------------------------------------------------------------------------------------------
// record access

const BASE_CHARS: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

#[inline]
fn base_code(base: u8) -> u8 {
    match base {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        _ => 15,
    }
}

#[inline]
fn seq_offset(rec: &bam::Record) -> usize {
    let c = &rec.inner.core;
    c.l_qname as usize + c.n_cigar as usize * 4
}

/// The base qualities, writable in place.
pub fn quals_mut(rec: &mut bam::Record) -> &mut [u8] {
    let len = rec.inner.core.l_qseq as usize;
    let off = seq_offset(rec) + ((len + 1) >> 1);
    unsafe { std::slice::from_raw_parts_mut(rec.inner.data.add(off), len) }
}

#[inline]
fn base_at(rec: &bam::Record, i: usize) -> u8 {
    let p = unsafe { *rec.inner.data.add(seq_offset(rec) + (i >> 1)) };
    BASE_CHARS[(if i & 1 == 0 { p >> 4 } else { p & 0xf }) as usize]
}

#[inline]
fn set_base(rec: &mut bam::Record, i: usize, base: u8) {
    let code = base_code(base);
    unsafe {
        let p = rec.inner.data.add(seq_offset(rec) + (i >> 1));
        *p = if i & 1 == 0 { (*p & 0x0f) | (code << 4) } else { (*p & 0xf0) | code };
    }
}

fn decode_bases(rec: &bam::Record, out: &mut Vec<u8>) {
    let len = rec.inner.core.l_qseq as usize;
    out.clear();
    out.reserve(len);
    let packed = unsafe { std::slice::from_raw_parts(rec.inner.data.add(seq_offset(rec)), (len + 1) >> 1) };
    for i in 0..len {
        let p = packed[i >> 1];
        out.push(BASE_CHARS[(if i & 1 == 0 { p >> 4 } else { p & 0xf }) as usize]);
    }
}

unsafe extern "C" {
    fn realloc(p: *mut std::ffi::c_void, size: usize) -> *mut std::ffi::c_void;
}

/// Replaces the cigar and, when given, the bases and qualities; the name and the attributes stay.
/// `buf` is scratch space.
fn replace_alignment_data(rec: &mut bam::Record, cigar: &[u32], seq_qual: Option<(&[u8], &[u8])>, buf: &mut Vec<u8>) {
    let (l_qname, old_seq_off, old_aux_off, l_data) = {
        let c = &rec.inner.core;
        let l_qname = c.l_qname as usize;
        let old_len = c.l_qseq as usize;
        let seq_off = l_qname + c.n_cigar as usize * 4;
        (l_qname, seq_off, seq_off + ((old_len + 1) >> 1) + old_len, rec.inner.l_data as usize)
    };
    buf.clear();
    {
        let data = unsafe { std::slice::from_raw_parts(rec.inner.data, l_data) };
        buf.extend_from_slice(&data[..l_qname]);
        for c in cigar {
            buf.extend_from_slice(&c.to_le_bytes());
        }
        match seq_qual {
            Some((seq, qual)) => {
                let mut i = 0;
                while i + 1 < seq.len() {
                    buf.push(base_code(seq[i]) << 4 | base_code(seq[i + 1]));
                    i += 2;
                }
                if i < seq.len() {
                    buf.push(base_code(seq[i]) << 4);
                }
                buf.extend_from_slice(qual);
            }
            None => buf.extend_from_slice(&data[old_seq_off..old_aux_off]),
        }
        buf.extend_from_slice(&data[old_aux_off..]);
    }
    let b = &mut rec.inner;
    unsafe {
        if buf.len() > b.m_data as usize {
            let p = realloc(b.data as *mut std::ffi::c_void, buf.len()) as *mut u8;
            assert!(!p.is_null());
            b.data = p;
            b.m_data = buf.len() as u32;
        }
        std::ptr::copy_nonoverlapping(buf.as_ptr(), b.data, buf.len());
    }
    b.l_data = buf.len() as i32;
    b.core.n_cigar = cigar.len() as u32;
    if let Some((seq, _)) = seq_qual {
        b.core.l_qseq = seq.len() as i32;
    }
    // SAMRecord.setCigar clears the indexing bin, which is then computed when the record is written
    fix_bin(rec);
}

/// Length in bytes of the attribute that starts at `aux[0]` (tag, type, value).
fn aux_field_len(aux: &[u8]) -> usize {
    let fixed = |t: u8| match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        _ => 0,
    };
    match aux[2] {
        b'Z' | b'H' => 3 + aux[3..].iter().position(|&c| c == 0).unwrap() + 1,
        b'B' => 8 + u32::from_le_bytes([aux[4], aux[5], aux[6], aux[7]]) as usize * fixed(aux[3]),
        t => 3 + fixed(t),
    }
}

/// The integer value of an attribute field, when it is an integer.
fn aux_field_int(f: &[u8]) -> Option<i64> {
    Some(match f[2] {
        b'c' => f[3] as i8 as i64,
        b'C' => f[3] as i64,
        b's' => i16::from_le_bytes([f[3], f[4]]) as i64,
        b'S' => u16::from_le_bytes([f[3], f[4]]) as i64,
        b'i' => i32::from_le_bytes([f[3], f[4], f[5], f[6]]) as i64,
        b'I' => u32::from_le_bytes([f[3], f[4], f[5], f[6]]) as i64,
        _ => return None,
    })
}

/// `rec::aux_sorted_copy` into a reused buffer: the attributes sorted by binary tag, each integer in the
/// smallest type that holds it.
fn aux_sorted_into(aux: &[u8], out: &mut Vec<u8>) {
    const MAX_FIELDS: usize = 48;
    let mut fields = [(0u16, 0u32, 0u32); MAX_FIELDS];
    let mut count = 0;
    let mut i = 0;
    let mut plain = true;
    while i + 3 <= aux.len() {
        let n = aux_field_len(&aux[i..]);
        let key = (aux[i + 1] as u16) << 8 | aux[i] as u16;
        // insertion sort; a repeated tag or too many attributes go the general way
        let mut at = count;
        while at > 0 && fields[at - 1].0 > key {
            at -= 1;
        }
        if count == MAX_FIELDS || (at > 0 && fields[at - 1].0 == key) {
            plain = false;
            break;
        }
        fields.copy_within(at..count, at + 1);
        fields[at] = (key, i as u32, n as u32);
        count += 1;
        i += n;
    }
    out.clear();
    if !plain {
        out.extend_from_slice(&aux_sorted_copy(aux));
        return;
    }
    for &(_, s, n) in &fields[..count] {
        let f = &aux[s as usize..(s + n) as usize];
        match aux_field_int(f) {
            Some(v) => {
                out.extend_from_slice(&f[..2]);
                // BinaryTagCodec.getIntegerType
                if v > i32::MAX as i64 {
                    out.push(b'I');
                    out.extend_from_slice(&(v as u32).to_le_bytes());
                } else if (-128..=127).contains(&v) {
                    out.push(b'c');
                    out.push(v as i8 as u8);
                } else if (0..=255).contains(&v) {
                    out.push(b'C');
                    out.push(v as u8);
                } else if (-32768..=32767).contains(&v) {
                    out.push(b's');
                    out.extend_from_slice(&(v as i16).to_le_bytes());
                } else if (0..=65535).contains(&v) {
                    out.push(b'S');
                    out.extend_from_slice(&(v as u16).to_le_bytes());
                } else {
                    out.push(b'i');
                    out.extend_from_slice(&(v as i32).to_le_bytes());
                }
            }
            None => out.extend_from_slice(f),
        }
    }
}

/// `rec::decode_attributes` without allocation.
fn decode_attributes_with(rec: &mut bam::Record, buf: &mut Vec<u8>) {
    aux_sorted_into(aux_bytes(rec), buf);
    set_aux_bytes(rec, buf);
}

/// The attributes that finaliseRead reads, found in one pass.
#[derive(Default)]
struct FinaliseAttributes {
    /// 0 absent or NONE, 1 SINGLE, 2 DUAL, 3 any other value
    consensus_type: u8,
    duplex_base_index: Option<i64>,
    has_num_mutations: bool,
    num_mutations: Option<i64>,
    alignment_score: Option<i64>,
}

fn finalise_attributes(aux: &[u8]) -> FinaliseAttributes {
    let mut a = FinaliseAttributes::default();
    let mut i = 0;
    while i + 3 <= aux.len() {
        let n = aux_field_len(&aux[i..]);
        let f = &aux[i..i + n];
        match &f[..2] {
            b"UT" => {
                a.consensus_type = if f[2] != b'Z' {
                    3
                } else {
                    match &f[3..n - 1] {
                        b"NONE" => 0,
                        b"SINGLE" => 1,
                        b"DUAL" => 2,
                        _ => 3,
                    }
                }
            }
            b"YX" => a.duplex_base_index = aux_field_int(f),
            b"NM" => {
                a.has_num_mutations = true;
                a.num_mutations = aux_field_int(f);
            }
            b"AS" => a.alignment_score = aux_field_int(f),
            _ => {}
        }
        i += n;
    }
    a
}

/// Whether htsjdk would write the record from its decoded fields: a setter for the name, bases, qualities, cigar
/// or an attribute has run on it. Kept in the record's `id`, which is not part of the BAM record.
#[inline]
pub fn is_decoded(rec: &bam::Record) -> bool {
    rec.inner.id == DECODED_ID
}
#[inline]
pub fn set_decoded(rec: &mut bam::Record, decoded: bool) {
    rec.inner.id = if decoded { DECODED_ID } else { 0 };
}
const DECODED_ID: u64 = 0x5342_5844;

/// The fields that SAMRecord.equals compares before it decodes both records.
pub fn fixed_fields_equal(a: &bam::Record, b: &bam::Record) -> bool {
    a.pos() == b.pos()
        && a.flags() == b.flags()
        && a.insert_size() == b.insert_size()
        && a.mapq() == b.mapq()
        && a.mpos() == b.mpos()
        && a.mtid() == b.mtid()
        && a.tid() == b.tid()
}

/// BAMRecord.eagerDecode, which SAMRecord.equals runs on two records with the same fixed fields (Java compares
/// reads by value in List.remove and List.contains): from then on the record is written from its decoded fields.
pub fn eager_decode(rec: &mut bam::Record) {
    if is_decoded(rec) {
        return;
    }
    // a quality change made in memory is now written
    if has_pending_quals(rec) {
        remove_pending_quals(rec);
    }
    decode_attributes(rec);
    set_decoded(rec, true);
}

/// List.remove(read) for each of the reads at `remove` in turn: the read is compared with each read before it
/// that is still in the list.
pub fn eager_decode_on_remove(reads: &mut [bam::Record], remove: &[usize]) {
    let mut current: Vec<usize> = (0..reads.len()).collect();
    for &x in remove {
        let Some(at) = current.iter().position(|&e| e == x) else { continue };
        for &e in &current[..at] {
            if fixed_fields_equal(&reads[x], &reads[e]) {
                eager_decode(&mut reads[x]);
                eager_decode(&mut reads[e]);
            }
        }
        current.remove(at);
    }
}

/// ConsensusReads.buildFromRead gives the new record the base and quality arrays of the read it copies, so what
/// finaliseRead then does to the bases and qualities of the consensus read is done to that read as well. It
/// reaches the BAM when the read is written from its decoded fields.
pub fn share_finalised_bases(consensus: &bam::Record, source: &mut bam::Record) {
    if !is_decoded(source) {
        return;
    }
    let len = consensus.inner.core.l_qseq as usize;
    if len != source.inner.core.l_qseq as usize {
        return;
    }
    let n = ((len + 1) >> 1) + len;
    unsafe {
        std::ptr::copy_nonoverlapping(consensus.inner.data.add(seq_offset(consensus)), source.inner.data.add(seq_offset(source)), n);
    }
}

/// The attribute that keeps base qualities which Java changed in memory only; see `prep_process_read`.
/// Its binary tag sorts first, so it stays the first attribute through `aux_sorted_copy` and `aux_block_set`.
const PENDING_QUALS_TAG: [u8; 2] = [0, 0];

#[inline]
fn has_pending_quals(rec: &bam::Record) -> bool {
    let aux = aux_bytes(rec);
    aux.len() >= 8 && aux[0] == 0 && aux[1] == 0 && aux[2] == b'B'
}

fn pending_quals_len(aux: &[u8]) -> usize {
    8 + u32::from_le_bytes([aux[4], aux[5], aux[6], aux[7]]) as usize
}

fn remove_pending_quals(rec: &mut bam::Record) {
    let aux = aux_bytes(rec);
    let n = pending_quals_len(aux);
    let rest = aux[n..].to_vec();
    set_aux_bytes(rec, &rest);
}

/// A read written without finaliseRead (a duplicate or a secondary read): base qualities that Java changed in
/// the decoded array of a record it then wrote from its stored bytes go back to their stored values.
pub fn restore_stored_quals(rec: &mut bam::Record) {
    if !has_pending_quals(rec) {
        return;
    }
    let aux = aux_bytes(rec);
    let n = pending_quals_len(aux);
    let stored = aux[8..n].to_vec();
    let rest = aux[n..].to_vec();
    let quals = quals_mut(rec);
    if quals.len() == stored.len() {
        quals.copy_from_slice(&stored);
    }
    set_aux_bytes(rec, &rest);
}

// ---------------------------------------------------------------------------------------------------------
// YC tag

/// The read indices of the duplex indels of a YC tag.
pub enum DuplexIndelIndices {
    /// Java null: the tag has no closing '+' or no duplex part
    Invalid,
    Found,
}

/// Integer.parseInt of a run of digits
fn parse_digits(s: &[u8]) -> SbxResult<i32> {
    let mut v: i64 = 0;
    if s.is_empty() {
        return err("java.lang.NumberFormatException");
    }
    for &c in s {
        if !c.is_ascii_digit() {
            return err("java.lang.NumberFormatException");
        }
        v = v * 10 + (c - b'0') as i64;
        if v > i32::MAX as i64 {
            return err("java.lang.NumberFormatException");
        }
    }
    Ok(v as i32)
}

/// Integer.parseInt: an optional sign, then digits
fn java_parse_int(s: &[u8]) -> SbxResult<i32> {
    match s.first() {
        Some(b'-') if s.len() > 1 => {
            let mut v: i64 = 0;
            for &c in &s[1..] {
                if !c.is_ascii_digit() {
                    return err("java.lang.NumberFormatException");
                }
                v = v * 10 + (c - b'0') as i64;
                if v > i32::MAX as i64 + 1 {
                    return err("java.lang.NumberFormatException");
                }
            }
            Ok((-v) as i32)
        }
        Some(b'+') if s.len() > 1 => parse_digits(&s[1..]),
        _ => parse_digits(s),
    }
}

/// SbxBamUtils.getDuplexIndelIndices: the indices go to `out`.
pub fn get_duplex_indel_indices(yc_tag: &[u8], out: &mut Vec<i32>) -> SbxResult<DuplexIndelIndices> {
    out.clear();
    let last_delim = match yc_tag.iter().rposition(|&c| c == b'+') {
        Some(i) if i > 0 => i,
        _ => return Ok(DuplexIndelIndices::Invalid),
    };
    let yc = &yc_tag[..last_delim];
    let mut simplex_head_length = 0;
    let duplex_region: &[u8];
    if yc.first() != Some(&b'+') {
        // String.split drops trailing empty items
        let mut components: Vec<&[u8]> = yc.split(|&c| c == b'+').collect();
        while components.last().is_some_and(|c| c.is_empty()) {
            components.pop();
        }
        if components.len() < 2 {
            return Ok(DuplexIndelIndices::Invalid);
        }
        simplex_head_length = java_parse_int(components[0])?;
        duplex_region = components[1];
    } else {
        duplex_region = &yc[1..];
    }
    let mut read_index = simplex_head_length;
    let mut i = 0;
    while i < duplex_region.len() {
        let c = duplex_region[i];
        if c.is_ascii_digit() {
            let mut j = i;
            while j < duplex_region.len() && duplex_region[j].is_ascii_digit() {
                j += 1;
            }
            read_index = read_index.wrapping_add(parse_digits(&duplex_region[i..j])?);
            i = j;
            continue;
        }
        i += 1;
        // the indel codes; the others are single base mismatches
        if matches!(c, b'I' | b'L' | b'P' | b'Q' | b'J' | b'O' | b'X' | b'Z') {
            out.push(read_index);
        }
        read_index = read_index.wrapping_add(1);
    }
    Ok(DuplexIndelIndices::Found)
}

/// SbxBamUtils.isHomopolymerLowBaseQualAtStart: for the first homopolymer whose two end qualities differ and
/// one of them is a duplex mismatch quality, whether that one is at the start
pub fn is_homopolymer_low_base_qual_at_start(bases: &[u8], quals: &[u8]) -> Option<bool> {
    if bases.is_empty() {
        return None;
    }
    let mut previous_base = bases[0];
    let mut previous_qual = quals[0] as i8;
    let mut hp_start_qual: i8 = 0;
    let mut in_homopolymer = false;
    for i in 1..bases.len() {
        let base = bases[i];
        let qual = quals[i] as i8;
        if in_homopolymer {
            if base != previous_base {
                in_homopolymer = false;
                if hp_start_qual != previous_qual {
                    if hp_start_qual <= RAW_DUPLEX_MISMATCH_QUAL as i8 {
                        return Some(true);
                    } else if previous_qual <= RAW_DUPLEX_MISMATCH_QUAL as i8 {
                        return Some(false);
                    }
                }
            }
        } else if base == previous_base {
            in_homopolymer = true;
            hp_start_qual = previous_qual;
        }
        previous_qual = qual;
        previous_base = base;
    }
    None
}

/// SbxBamSampler: the first sampled read that tells sets which end of a homopolymer has the low quality.
#[derive(Default)]
pub struct SbxBamSampler {
    direction_set: bool,
}

impl SbxBamSampler {
    pub fn process_sample_record(&mut self, rec: &bam::Record) {
        if self.direction_set {
            return;
        }
        let bases = rec.seq().as_bytes();
        let Some(at_start) = is_homopolymer_low_base_qual_at_start(&bases, rec.qual()) else { return };
        self.direction_set = true;
        set_homopolymer_5_prime_low_base_qual(at_start == !has(rec, FLAG_REVERSE));
        crate::debug!("SBX homopolymer low base-qual at 5' end({})", homopolymer_5_prime_low_base_qual());
    }
}

// ---------------------------------------------------------------------------------------------------------
// duplex indels

/// SbxDuplexIndel
#[derive(Debug, Clone)]
pub struct SbxDuplexIndel {
    pub duplex_indel_index_start: i32,
    pub duplex_indel_index_end: i32,
    pub repeat_bases: Vec<u8>,
    pub first_read_insert_index: i32,
    pub low_base_qual_indices: Vec<i32>,
    pub deleted_indel_index_start: i32,
    pub deleted_indel_index_end: i32,
}

impl SbxDuplexIndel {
    fn within_bounds(&self, read_index: i32) -> bool {
        read_index >= self.first_read_insert_index && read_index <= self.duplex_indel_index_end
    }
    fn within_delete_bounds(&self, read_index: i32) -> bool {
        read_index >= self.deleted_indel_index_start && read_index <= self.deleted_indel_index_end
    }
    fn is_low_qual_base(&self, read_index: i32) -> bool {
        self.low_base_qual_indices.contains(&read_index)
    }
}

fn cigar_from_str(cigar: &[u8]) -> SbxResult<Vec<u32>> {
    // CigarUtils.cigarElementsFromStr
    if cigar == b"*" {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    let mut len: i32 = 0;
    for &c in cigar {
        if c.is_ascii_digit() {
            len = len.wrapping_mul(10).wrapping_add((c - b'0') as i32);
        } else {
            let op = match c {
                b'M' => OP_M,
                b'I' => OP_I,
                b'D' => OP_D,
                b'N' => OP_N,
                b'S' => OP_S,
                b'H' => OP_H,
                b'P' => OP_P,
                b'=' => OP_EQ,
                b'X' => OP_X,
                _ => return err("java.lang.IllegalArgumentException: Unrecognized CigarOperator"),
            };
            if len < 0 {
                return err("java.lang.IllegalArgumentException: Cigar element being constructed with negative length");
            }
            out.push((len as u32) << 4 | op);
            len = 0;
        }
    }
    Ok(out)
}

#[inline]
fn el(len: i32, op: u32) -> u32 {
    (len as u32) << 4 | op
}

/// SbxDuplexIndelBuilder.trimCigar
pub fn trim_cigar(elements: &mut Vec<u32>, length: i32, from_start: bool) -> SbxResult<()> {
    let mut remaining = length;
    let mut index: i64 = if from_start { 0 } else { elements.len() as i64 - 1 };
    while remaining > 0 {
        if index < 0 || index as usize >= elements.len() {
            return err("java.lang.IndexOutOfBoundsException");
        }
        let e = elements[index as usize];
        if oplen(e) > remaining {
            elements[index as usize] = el(oplen(e) - remaining, op(e));
            return Ok(());
        }
        elements.remove(index as usize);
        remaining -= oplen(e);
        if !from_start {
            index -= 1;
        }
    }
    Ok(())
}

/// SbxDuplexIndelBuilder.checkSupplementaryCigar: the read-base elements of the supplementary alignment's
/// cigar, made to cover exactly the read's soft clip
pub fn check_supplementary_cigar(supp_elements: &[u32], soft_clip_length: i32, is_left_clipped: bool) -> SbxResult<Vec<u32>> {
    if supp_elements.is_empty() {
        return err("java.lang.IndexOutOfBoundsException");
    }
    let count = supp_elements.len();
    // left clipped is about the read, so the other side of the supplementary cigar is trimmed
    let trim_left_side = !is_left_clipped;
    let (matching_soft_clip, rest) =
        if trim_left_side { (supp_elements[0], &supp_elements[1..]) } else { (supp_elements[count - 1], &supp_elements[..count - 1]) };
    let mut new_elements: Vec<u32> = rest.iter().copied().filter(|e| consumes_read(op(*e))).collect();
    let read_base_length: i32 = new_elements.iter().map(|e| oplen(*e)).sum();
    if read_base_length == soft_clip_length {
        return Ok(new_elements);
    }
    let diff = read_base_length - soft_clip_length;
    if diff > 0 {
        trim_cigar(&mut new_elements, diff, trim_left_side)?;
    } else {
        let extra = el(diff.abs(), op(matching_soft_clip));
        if trim_left_side {
            new_elements.insert(0, extra);
        } else {
            new_elements.push(extra);
        }
    }
    Ok(new_elements)
}

#[derive(Clone, Copy)]
struct ReadBaseInfo {
    index: i32,
    base: u8,
    cigar_index: i32,
    cigar_op: u32,
    index_in_cigar: i32,
    valid: bool,
}

/// SbxDuplexIndelBuilder: the duplex indels of a read from the indices of its duplex mismatches, walked over the
/// read's cigar with the soft clip replaced by the supplementary alignment's elements
struct DuplexIndelBuilder<'a> {
    bases: &'a [u8],
    read_length: i32,
    negative_strand: bool,
    adjusted_cigar: Vec<u32>,
    duplex_indels: Vec<SbxDuplexIndel>,
}

impl<'a> DuplexIndelBuilder<'a> {
    fn base(&self, index: i32) -> SbxResult<u8> {
        if index < 0 || index >= self.bases.len() as i32 {
            return err("java.lang.ArrayIndexOutOfBoundsException");
        }
        Ok(self.bases[index as usize])
    }
    fn element(&self, cigar_index: i32) -> SbxResult<u32> {
        if cigar_index < 0 || cigar_index as usize >= self.adjusted_cigar.len() {
            return err("java.lang.IndexOutOfBoundsException");
        }
        Ok(self.adjusted_cigar[cigar_index as usize])
    }

    fn build(
        bases: &'a [u8], cigar: &[u32], negative_strand: bool, supp: Option<(&[u8], bool)>, duplex_indel_indices: &[i32],
    ) -> SbxResult<Vec<SbxDuplexIndel>> {
        let mut b = DuplexIndelBuilder {
            bases,
            read_length: bases.len() as i32,
            negative_strand,
            adjusted_cigar: Vec::new(),
            duplex_indels: Vec::with_capacity(duplex_indel_indices.len()),
        };
        match supp {
            Some((supp_cigar, supp_forward)) => b.build_adjusted_supplementary_cigar(cigar, supp_cigar, supp_forward)?,
            None => b.adjusted_cigar.extend_from_slice(cigar),
        }
        let mut i = 0;
        while i < duplex_indel_indices.len() {
            let duplex_indel_index_start = duplex_indel_indices[i];
            let mut duplex_indel_index_end = duplex_indel_index_start;
            let read_base_info = b.get_read_base_info(duplex_indel_index_start)?;
            let mut next_read_base_info: Option<ReadBaseInfo> = None;
            for &next_duplex_indel_index in &duplex_indel_indices[i + 1..] {
                if next_duplex_indel_index > duplex_indel_index_end + 1 {
                    break;
                }
                if next_read_base_info.is_none() {
                    next_read_base_info = Some(read_base_info.ok_or_else(|| SbxError("java.lang.NullPointerException".to_string()))?);
                }
                let next = next_read_base_info.as_mut().unwrap();
                b.move_next(next)?;
                if !consumes_read(next.cigar_op) {
                    break;
                }
                duplex_indel_index_end = next_duplex_indel_index;
            }
            b.process_duplex_indel(read_base_info, duplex_indel_index_start, duplex_indel_index_end)?;
            i += (duplex_indel_index_end - duplex_indel_index_start + 1) as usize;
        }
        Ok(b.duplex_indels)
    }

    fn build_adjusted_supplementary_cigar(&mut self, cigar: &[u32], supp_cigar: &[u8], supp_forward: bool) -> SbxResult<()> {
        if cigar.is_empty() {
            return err("java.lang.IndexOutOfBoundsException");
        }
        let left_soft_clip_length = left_soft_clip(cigar);
        let last_element = cigar.len() - 1;
        let right_soft_clip_length = right_soft_clip(cigar);
        let mut supp_elements = cigar_from_str(supp_cigar)?;
        if supp_forward == self.negative_strand {
            supp_elements.reverse();
        }
        if left_soft_clip_length > right_soft_clip_length {
            let supp = check_supplementary_cigar(&supp_elements, left_soft_clip_length, true)?;
            self.adjusted_cigar.extend_from_slice(&supp);
            self.adjusted_cigar.extend_from_slice(&cigar[1..]);
        } else {
            // also taken without any soft clip, when it drops the read's last element (as Java does)
            self.adjusted_cigar.extend_from_slice(&cigar[..last_element]);
            let supp = check_supplementary_cigar(&supp_elements, right_soft_clip_length, false)?;
            self.adjusted_cigar.extend_from_slice(&supp);
        }
        Ok(())
    }

    fn get_read_base_info(&self, target_read_index: i32) -> SbxResult<Option<ReadBaseInfo>> {
        let mut read_index = 0;
        for (cigar_index, &e) in self.adjusted_cigar.iter().enumerate() {
            if !consumes_read(op(e)) {
                continue;
            }
            if read_index + oplen(e) > target_read_index {
                return Ok(Some(ReadBaseInfo {
                    index: target_read_index,
                    base: self.base(target_read_index)?,
                    cigar_index: cigar_index as i32,
                    cigar_op: op(e),
                    index_in_cigar: target_read_index - read_index,
                    valid: true,
                }));
            }
            read_index += oplen(e);
        }
        Ok(None)
    }

    fn move_to(&self, info: &mut ReadBaseInfo, target_read_index: i32) -> SbxResult<()> {
        while target_read_index > info.index {
            self.move_next(info)?;
        }
        while target_read_index < info.index {
            self.move_previous(info)?;
        }
        Ok(())
    }

    fn move_previous(&self, info: &mut ReadBaseInfo) -> SbxResult<()> {
        info.index -= 1;
        if info.index < 0 {
            info.valid = false;
            return Ok(());
        }
        info.base = self.base(info.index)?;
        self.element(info.cigar_index)?;
        info.index_in_cigar -= 1;
        if info.index_in_cigar < 0 {
            info.cigar_index -= 1;
            if info.cigar_index < 0 {
                info.valid = false;
                return Ok(());
            }
            let e = self.element(info.cigar_index)?;
            info.cigar_op = op(e);
            info.index_in_cigar = oplen(e) - 1;
        }
        Ok(())
    }

    /// Steps one read base on. An element that has no read bases is walked as if it had (as Java does).
    fn move_next(&self, info: &mut ReadBaseInfo) -> SbxResult<()> {
        info.index += 1;
        if info.index >= self.read_length {
            info.valid = false;
            return Ok(());
        }
        info.base = self.base(info.index)?;
        let e = self.element(info.cigar_index)?;
        info.index_in_cigar += 1;
        if info.index_in_cigar >= oplen(e) {
            info.cigar_index += 1;
            let e = self.element(info.cigar_index)?;
            info.cigar_op = op(e);
            info.index_in_cigar = 0;
        }
        Ok(())
    }

    fn can_shift_by_length(&self, initial: &ReadBaseInfo, repeat_bases: &[u8], existing_shift: i32, new_shift: i32) -> SbxResult<bool> {
        let target_start_index = initial.index + existing_shift + new_shift;
        let mut info = *initial;
        self.move_to(&mut info, target_start_index)?;
        if !info.valid {
            return Ok(false);
        }
        for &repeat_base in repeat_bases {
            if !info.valid || !(info.cigar_op == OP_M || info.cigar_op == OP_I) || repeat_base != info.base {
                return Ok(false);
            }
            self.move_next(&mut info)?;
        }
        Ok(true)
    }

    /// SbxDuplexIndelBuilder.findRepeatShiftLength: how far the repeat can move in one direction and still
    /// match the read, in whole repeats and then a part of one
    fn find_repeat_shift_length(&self, initial: &ReadBaseInfo, repeat_bases: &[u8], shift_up: bool) -> SbxResult<i32> {
        let mut shift_bases = 0;
        let repeat_length = repeat_bases.len() as i32;
        loop {
            let shift_count = if shift_up { repeat_length } else { -repeat_length };
            if !self.can_shift_by_length(initial, repeat_bases, shift_bases, shift_count)? {
                break;
            }
            shift_bases += shift_count;
        }
        if repeat_length > 1 {
            let is_homopolymer = repeat_bases.iter().all(|b| *b == repeat_bases[0]);
            let mut max_partial_shift = 0;
            for i in 1..repeat_length {
                let partial_shift = if shift_up { i } else { -i };
                if !self.can_shift_by_length(initial, repeat_bases, shift_bases, partial_shift)? {
                    if is_homopolymer {
                        break;
                    }
                } else {
                    max_partial_shift = partial_shift;
                }
            }
            shift_bases += max_partial_shift;
        }
        Ok(shift_bases)
    }

    fn process_duplex_indel(
        &mut self, read_base_info: Option<ReadBaseInfo>, mut duplex_indel_index_start: i32, mut duplex_indel_index_end: i32,
    ) -> SbxResult<()> {
        let duplex_mismatch_length = duplex_indel_index_end - duplex_indel_index_start + 1;
        let mut repeat_bases = Vec::with_capacity(duplex_mismatch_length as usize);
        for i in 0..duplex_mismatch_length {
            repeat_bases.push(self.base(duplex_indel_index_start + i)?);
        }
        let Some(mut info) = read_base_info else { return err("java.lang.NullPointerException") };

        let right_alignment_duplex_indel = homopolymer_5_prime_low_base_qual() == !self.negative_strand;
        if right_alignment_duplex_indel {
            let right_shift = self.find_repeat_shift_length(&info, &repeat_bases, true)?;
            if right_shift > 0 {
                duplex_indel_index_start += right_shift;
                duplex_indel_index_end += right_shift;
                self.move_to(&mut info, duplex_indel_index_start)?;
            }
        }

        // the inserted bases at the start of the duplex mismatch
        let mut duplex_mismatch_insert_count = 0;
        while info.index <= duplex_indel_index_end {
            if info.cigar_op == OP_I {
                duplex_mismatch_insert_count += 1;
            } else {
                break;
            }
            self.move_next(&mut info)?;
        }
        let duplex_mismatch_in_insert = duplex_mismatch_insert_count > 0;
        let mut first_read_insert_index = if duplex_mismatch_in_insert { duplex_indel_index_start } else { -1 };

        // back through the repeat to find its first inserted base
        self.move_to(&mut info, duplex_indel_index_start)?;
        let insert_shift_count = self.find_repeat_shift_length(&info, &repeat_bases, false)?;
        for _ in 0..insert_shift_count.abs() {
            self.move_previous(&mut info)?;
            if !info.valid {
                break;
            }
            if info.cigar_op == OP_I {
                first_read_insert_index = info.index;
            }
        }
        // the repeat is either inside the inserted bases or runs on past the insert without a gap
        if first_read_insert_index < 0 || (insert_shift_count == 0 && !duplex_mismatch_in_insert) {
            return Ok(());
        }
        let insert_repeat_shift = insert_shift_count.abs();
        let total_inserted_bases = insert_repeat_shift + duplex_mismatch_insert_count;
        let trim_length = total_inserted_bases.min(duplex_mismatch_length);
        let mut low_base_qual_count = duplex_mismatch_length;
        let trimmed_repeat_base_length = insert_repeat_shift + duplex_mismatch_length - trim_length;
        // enough low quality bases to be symmetric around the repeat
        if low_base_qual_count < trimmed_repeat_base_length && duplex_mismatch_length % 2 == 1 {
            low_base_qual_count += 1;
        }
        let mut trim_count = 0;
        let mut deleted_indel_index_start = -1;
        let mut deleted_indel_index_end = -1;
        self.move_to(&mut info, first_read_insert_index)?;
        if !info.valid {
            return Ok(());
        }
        let mut low_qual_indices: Vec<i32> = Vec::with_capacity(low_base_qual_count as usize);
        let mut remaining_low_qual_count = low_base_qual_count;
        let mut added_low_qual_at_start_count = 0;
        for i in first_read_insert_index..=duplex_indel_index_end {
            if info.cigar_op == OP_I && trim_count < trim_length {
                trim_count += 1;
                if deleted_indel_index_start < 0 {
                    deleted_indel_index_start = i;
                }
                deleted_indel_index_end = i;
            } else if remaining_low_qual_count > 0 {
                if !low_qual_indices.contains(&i) {
                    low_qual_indices.push(i);
                }
                remaining_low_qual_count -= 1;
                added_low_qual_at_start_count += 1;
                if remaining_low_qual_count > 0 {
                    let other_low_qual_index = duplex_indel_index_end - (added_low_qual_at_start_count - 1);
                    if !low_qual_indices.contains(&other_low_qual_index) {
                        low_qual_indices.push(other_low_qual_index);
                    }
                    remaining_low_qual_count -= 1;
                }
            }
            if trim_count == trim_length && remaining_low_qual_count == 0 {
                break;
            }
            self.move_next(&mut info)?;
            if !info.valid {
                break;
            }
        }
        if deleted_indel_index_start < 0 || deleted_indel_index_end < deleted_indel_index_start {
            return Ok(());
        }
        self.duplex_indels.push(SbxDuplexIndel {
            duplex_indel_index_start,
            duplex_indel_index_end,
            repeat_bases,
            first_read_insert_index,
            low_base_qual_indices: low_qual_indices,
            deleted_indel_index_start,
            deleted_indel_index_end,
        });
        Ok(())
    }
}

/// The duplex indels of a read (SbxDuplexIndelBuilder over the read and its duplex indel indices).
pub fn build_duplex_indels(rec: &bam::Record, duplex_indel_indices: &[i32]) -> SbxResult<Vec<SbxDuplexIndel>> {
    let bases = rec.seq().as_bytes();
    let supp = extract_supp_alignment(aux_str(rec, b"SA")).map(|s| (s.cigar, s.strand == b'+'));
    DuplexIndelBuilder::build(&bases, rec.raw_cigar(), has(rec, FLAG_REVERSE), supp, duplex_indel_indices)
}

// ---------------------------------------------------------------------------------------------------------
// cigar corrections

/// SbxRoutines.hasInvalidCigar: an insert before the first or after the last aligned section
pub fn has_invalid_cigar(cigar: &[u32]) -> bool {
    if op(cigar[0]) == OP_I {
        return true;
    }
    if cigar.len() > 1 {
        let last_index = cigar.len() - 1;
        if op(cigar[last_index]) == OP_I {
            return true;
        }
        if op(cigar[0]) == OP_S && op(cigar[1]) == OP_I {
            return true;
        }
        if op(cigar[last_index - 1]) == OP_I && op(cigar[last_index]) == OP_S {
            return true;
        }
    }
    false
}

/// SbxRoutines.correctInvalidCigar: the inserts at the ends become soft clips
pub fn correct_invalid_cigar(cigar: &mut Vec<u32>) {
    let mut last_index = cigar.len() as i64 - 1;
    let mut i: i64 = 0;
    while (i as usize) < cigar.len() {
        let element = cigar[i as usize];
        if op(element) == OP_I {
            if i == 0 || i == last_index {
                cigar[i as usize] = el(oplen(element), OP_S);
            } else if i == 1 && op(cigar[0]) == OP_S {
                cigar.remove(i as usize);
                last_index -= 1;
                cigar[0] = el(oplen(cigar[0]) + oplen(element), OP_S);
                continue;
            } else if i == last_index - 1 && op(cigar[last_index as usize]) == OP_S {
                cigar.remove(i as usize);
                last_index -= 1;
                let last = cigar[last_index as usize];
                cigar[last_index as usize] = el(oplen(last) + oplen(element), OP_S);
                continue;
            }
        }
        i += 1;
    }
}

/// CigarUtils.checkLeftAlignment: inserts that can move left through a repeat in the aligned bases before them
/// are moved; true when the cigar changed
pub fn check_left_alignment(cigar: &mut Vec<u32>, read_bases: &[u8]) -> bool {
    let mut modified = false;
    let mut convert_first_insert = false;
    let mut read_index: i32 = 0;
    let mut i = 0;
    while i + 1 < cigar.len() {
        let element = cigar[i];
        let next_element = cigar[i + 1];
        let prev_element = if i > 0 { Some(cigar[i - 1]) } else { None };
        let mut check_left_shift = match prev_element {
            Some(p) => op(element) == OP_I && op(p) == OP_M && op(next_element) == OP_M,
            None => false,
        };
        if check_left_shift {
            let p = prev_element.unwrap();
            // no shift back through the whole previous element, except as a conversion to a soft clip
            if oplen(element) > oplen(p) || (oplen(element) == oplen(p) && i != 1) {
                check_left_shift = false;
            }
        }
        if check_left_shift {
            let p = prev_element.unwrap();
            let insert_len = oplen(element) as usize;
            if read_index < 0 || read_index as usize + insert_len > read_bases.len() {
                return false;
            }
            let insert = &read_bases[read_index as usize..read_index as usize + insert_len];
            let is_homopolymer = insert.iter().all(|b| *b == insert[0]);
            let repeat: &[u8] = if is_homopolymer && insert_len > 1 { &insert[..1] } else { insert };
            let repeat_len = repeat.len() as i32;
            let mut prior_index = read_index;
            let mut matched_repeat_count = 0;
            loop {
                prior_index -= repeat_len;
                if prior_index < 0 {
                    break;
                }
                if repeat != &read_bases[prior_index as usize..(prior_index + repeat_len) as usize] {
                    break;
                }
                matched_repeat_count += 1;
                if (matched_repeat_count + 1) * repeat_len > oplen(p) {
                    break;
                }
            }
            if matched_repeat_count > 0 {
                let mut alignment_shift = matched_repeat_count * repeat_len;
                if alignment_shift == oplen(p) {
                    if i == 1 {
                        // the first element goes and the insert becomes a soft clip
                        convert_first_insert = true;
                    } else {
                        alignment_shift = (matched_repeat_count - 1) * repeat_len;
                    }
                }
                cigar[i - 1] = el((oplen(p) - alignment_shift).max(1), op(p));
                cigar[i + 1] = el(oplen(next_element) + alignment_shift, op(next_element));
                read_index -= alignment_shift;
                modified = true;
            }
        }
        if consumes_read(op(element)) {
            read_index += oplen(element);
        }
        i += 1;
    }
    if convert_first_insert {
        cigar.remove(0);
        cigar[0] = el(oplen(cigar[0]), OP_S);
    }
    modified
}

/// One alignment of an SA attribute with all its fields, as SupplementaryReadData holds it.
struct SuppAlignment {
    chromosome: Vec<u8>,
    position: i32,
    strand: u8,
    cigar: Vec<u8>,
    map_quality: i32,
    nm: i32,
}

impl SuppAlignment {
    /// SupplementaryReadData.asSamTag
    fn write_sam_tag(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.chromosome);
        out.extend_from_slice(format!(",{},", self.position).as_bytes());
        out.push(self.strand);
        out.push(b',');
        out.extend_from_slice(&self.cigar);
        out.extend_from_slice(format!(",{},{}", self.map_quality, self.nm).as_bytes());
    }
}

/// SupplementaryReadData.extractAlignments with every field: Ok(None) when absent or malformed
fn extract_full_supp_alignments(sa: Option<&[u8]>) -> SbxResult<Option<Vec<SuppAlignment>>> {
    let Some(sa) = sa else { return Ok(None) };
    if sa.is_empty() {
        return Ok(None);
    }
    let mut parts: Vec<&[u8]> = sa.split(|&c| c == b';').collect();
    if sa.last() == Some(&b';') {
        parts.pop();
    }
    let mut out = Vec::with_capacity(parts.len());
    for p in parts {
        let mut items: Vec<&[u8]> = p.split(|&c| c == b',').collect();
        while items.last().is_some_and(|i| i.is_empty()) {
            items.pop();
        }
        if items.len() != 6 {
            return Ok(None);
        }
        let Some(&strand) = items[2].first() else { return err("java.lang.StringIndexOutOfBoundsException") };
        out.push(SuppAlignment {
            chromosome: items[0].to_vec(),
            position: java_parse_int(items[1])?,
            strand,
            cigar: items[3].to_vec(),
            map_quality: java_parse_int(items[4])?,
            nm: java_parse_int(items[5])?,
        });
    }
    Ok(Some(out))
}

/// SbxRoutines.correctSupplementaryData: the longer soft clip of the supplementary alignment's cigar loses the
/// bases stripped from this read
fn correct_supplementary_cigar(cigar: &[u8], soft_clip_reduction: i32) -> SbxResult<Vec<u8>> {
    let mut elements = cigar_from_str(cigar)?;
    if elements.is_empty() {
        return err("java.lang.IndexOutOfBoundsException");
    }
    let last_index = elements.len() - 1;
    let left = left_soft_clip(&elements);
    let right = right_soft_clip(&elements);
    let target = if left > right {
        Some(0)
    } else if right > left {
        Some(last_index)
    } else {
        None
    };
    if let Some(t) = target {
        let new_len = oplen(elements[t]) - soft_clip_reduction;
        if new_len < 0 {
            return err("java.lang.IllegalArgumentException: Cigar element being constructed with negative length");
        }
        elements[t] = el(new_len, OP_S);
    }
    Ok(cigar_string(&elements).into_bytes())
}

// ---------------------------------------------------------------------------------------------------------
// pre-processing

#[derive(Default)]
struct Scratch {
    buf: Vec<u8>,
    indices: Vec<i32>,
    bases: Vec<u8>,
    mismatch_indices: Vec<i32>,
    mismatch_ref_bases: Vec<(i32, u8)>,
}

thread_local! {
    static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
}

/// What stripDuplexIndels did to the record.
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum StripOutcome {
    Unchanged,
    /// a Java setter ran, so htsjdk writes the record from its decoded fields
    Modified,
}

fn fatal(msg: String) -> ! {
    crate::error!("{msg}");
    std::process::exit(1);
}

fn read_to_string(rec: &bam::Record) -> String {
    format!(
        "id({}) tid({}) start({}) cigar({}) flags({})",
        String::from_utf8_lossy(rec.qname()), rec.tid(), start(rec), cigar_string(rec.raw_cigar()), rec.flags()
    )
}

/// SbxRoutines.stripDuplexIndels: bases that only one strand of a duplex read inserted are removed where they
/// are in an insert of the alignment, and the bases left around them get the duplex mismatch quality
pub fn strip_duplex_indels(rec: &mut bam::Record) -> SbxResult<StripOutcome> {
    let mut decoded = false;
    SCRATCH.with(|s| strip_duplex_indels_with(rec, &mut s.borrow_mut(), &mut decoded))?;
    Ok(if decoded { StripOutcome::Modified } else { StripOutcome::Unchanged })
}

/// `decoded` is set as soon as a setter has run on the record.
fn strip_duplex_indels_with(rec: &mut bam::Record, scratch: &mut Scratch, decoded: &mut bool) -> SbxResult<()> {
    if has(rec, FLAG_UNMAPPED) {
        return Ok(());
    }
    let Some(yc_tag) = aux_str(rec, b"YC") else {
        return Err(SbxError(format!("java.lang.IllegalArgumentException: read missing YC tag: {}", read_to_string(rec))));
    };
    if rec.raw_cigar().is_empty() {
        return err("java.lang.IndexOutOfBoundsException");
    }
    let found = get_duplex_indel_indices(yc_tag, &mut scratch.indices);

    // inserts at the ends become soft clips
    if has_invalid_cigar(rec.raw_cigar()) {
        let mut fixed = rec.raw_cigar().to_vec();
        correct_invalid_cigar(&mut fixed);
        replace_alignment_data(rec, &fixed, None, &mut scratch.buf);
        decode_attributes_with(rec, &mut scratch.buf);
        *decoded = true;
    }
    match found? {
        DuplexIndelIndices::Invalid => fatal(format!(
            "invalid SBX ycTag({}) in read: {}",
            String::from_utf8_lossy(aux_str(rec, b"YC").unwrap_or_default()), read_to_string(rec)
        )),
        DuplexIndelIndices::Found => {}
    }
    {
        let cigar = rec.raw_cigar();
        if op(cigar[0]) == OP_H || op(cigar[cigar.len() - 1]) == OP_H {
            fatal(format!("hard-clipped reads not supported in SBX: {}", read_to_string(rec)));
        }
    }
    if scratch.indices.is_empty() {
        return Ok(());
    }
    decode_bases(rec, &mut scratch.bases);
    let old_base_length = scratch.bases.len();
    if has(rec, FLAG_REVERSE) {
        // SbxRoutines.reverseDuplexIndelIndices
        scratch.indices.reverse();
        for i in scratch.indices.iter_mut() {
            *i = old_base_length as i32 - *i - 1;
        }
    }
    let supp = extract_supp_alignment(aux_str(rec, b"SA")).map(|s| (s.cigar, s.strand == b'+'));
    let duplex_indels = DuplexIndelBuilder::build(&scratch.bases, rec.raw_cigar(), has(rec, FLAG_REVERSE), supp, &scratch.indices)?;
    if duplex_indels.is_empty() {
        return Ok(());
    }

    // the deleted indices, which several duplex indels can share
    let mut net_stripped: Vec<i32> = Vec::new();
    for d in &duplex_indels {
        for i in d.deleted_indel_index_start..=d.deleted_indel_index_end {
            if !net_stripped.contains(&i) {
                net_stripped.push(i);
            }
        }
    }
    let new_base_length = old_base_length as i64 - net_stripped.len() as i64;
    if new_base_length < 0 {
        return err("java.lang.NegativeArraySizeException");
    }
    let new_base_length = new_base_length as usize;
    let mut read_bases = vec![0u8; new_base_length];
    let mut read_quals = vec![0u8; new_base_length];
    let old_cigar: Vec<u32> = rec.raw_cigar().to_vec();
    let mut new_cigar: Vec<u32> = Vec::with_capacity(old_cigar.len());
    {
        let old_quals = rec.qual();
        let old_bases = &scratch.bases;
        let mut min_duplex_indel_index = 0usize;
        let mut indel_trimmed_count = 0i32;
        let mut old_cigar_index = 0usize;
        let mut old_cigar_element_index = 0i32;
        let mut old_element = old_cigar[0];
        let mut current_op: Option<u32> = None;
        let mut current_len = 0i32;
        let mut old_insert_count = 0i64;
        let mut old_inserted_bases = 0i64;
        let mut new_read_index = 0usize;
        let mut old_read_index = 0usize;
        while old_read_index < old_base_length {
            if old_cigar_element_index >= oplen(old_element) {
                old_cigar_element_index = 0;
                old_cigar_index += 1;
                if old_cigar_index >= old_cigar.len() {
                    return err("java.lang.IndexOutOfBoundsException");
                }
                old_element = old_cigar[old_cigar_index];
                if op(old_element) == OP_I {
                    old_insert_count += 1;
                    old_inserted_bases += oplen(old_element) as i64;
                }
            }
            let old_op = op(old_element);
            // a run of bases before the next duplex indel is copied as it is; the walk below would do the same
            // base by base
            {
                let element_remaining = oplen(old_element) - old_cigar_element_index;
                let mut run = element_remaining;
                if consumes_read(old_op) {
                    let effect_read_index = old_read_index as i32;
                    run = run.min((old_base_length - old_read_index) as i32);
                    for d in &duplex_indels[min_duplex_indel_index..] {
                        if effect_read_index > d.duplex_indel_index_end {
                            continue;
                        }
                        run = run.min(d.first_read_insert_index - effect_read_index);
                        if run <= 1 {
                            break;
                        }
                    }
                    if run > 1 && new_read_index + run as usize > new_base_length {
                        run = 0;
                    }
                }
                if run > 1 && new_read_index as i32 + indel_trimmed_count == old_read_index as i32 {
                    if current_op != Some(old_op) {
                        if current_len > 0 {
                            new_cigar.push(el(current_len, current_op.unwrap()));
                        }
                        current_op = Some(old_op);
                        current_len = run;
                    } else {
                        current_len += run;
                    }
                    if consumes_read(old_op) {
                        let n = run as usize;
                        read_bases[new_read_index..new_read_index + n].copy_from_slice(&old_bases[old_read_index..old_read_index + n]);
                        read_quals[new_read_index..new_read_index + n].copy_from_slice(&old_quals[old_read_index..old_read_index + n]);
                        new_read_index += n;
                        old_read_index += n;
                    }
                    old_cigar_element_index += run;
                    continue;
                }
            }
            let mut is_stripped_indel_base = false;
            let mut is_low_qual_base = false;
            let mut within_duplex_indel_bounds = false;
            // the index in the old read
            let effect_read_index = new_read_index as i32 + indel_trimmed_count;
            if consumes_read(old_op) {
                let mut first_active_index: i64 = -1;
                for (i, d) in duplex_indels.iter().enumerate().skip(min_duplex_indel_index) {
                    if effect_read_index > d.duplex_indel_index_end {
                        continue;
                    }
                    if first_active_index < 0 {
                        first_active_index = i as i64;
                    }
                    if !d.within_bounds(effect_read_index) {
                        continue;
                    }
                    within_duplex_indel_bounds = true;
                    if d.within_delete_bounds(effect_read_index) {
                        is_stripped_indel_base = true;
                        break;
                    }
                    if d.is_low_qual_base(effect_read_index) {
                        is_low_qual_base = true;
                    }
                }
                if first_active_index > 0 {
                    min_duplex_indel_index = first_active_index as usize;
                }
            }
            if is_stripped_indel_base {
                indel_trimmed_count += 1;
                is_low_qual_base = false;
            } else if current_op != Some(old_op) {
                if current_len > 0 {
                    new_cigar.push(el(current_len, current_op.unwrap()));
                }
                current_op = Some(old_op);
                current_len = 1;
            } else {
                current_len += 1;
            }
            if consumes_read(old_op) {
                if !is_stripped_indel_base {
                    if new_read_index >= new_base_length || old_read_index >= old_base_length {
                        break;
                    }
                    read_bases[new_read_index] = old_bases[old_read_index];
                    let existing_qual = old_quals[old_read_index];
                    read_quals[new_read_index] = if is_low_qual_base {
                        SBX_DUPLEX_MISMATCH_QUAL
                    } else if existing_qual == RAW_DUPLEX_MISMATCH_QUAL && within_duplex_indel_bounds {
                        RAW_DUPLEX_QUAL
                    } else {
                        existing_qual
                    };
                    new_read_index += 1;
                }
                old_read_index += 1;
            }
            old_cigar_element_index += 1;
        }
        let Some(last_op) = current_op else { return err("java.lang.NullPointerException") };
        new_cigar.push(el(current_len, last_op));

        correct_invalid_cigar(&mut new_cigar);
        check_left_alignment(&mut new_cigar, &read_bases);

        let mut new_insert_count = 0i64;
        let mut new_inserted_bases = 0i64;
        for &e in &new_cigar {
            if op(e) == OP_I {
                new_insert_count += 1;
                new_inserted_bases += oplen(e) as i64;
            }
        }
        let nm_diff = new_inserted_bases - old_inserted_bases;
        *decoded = true;
        let mut aux = aux_sorted_copy(aux_bytes(rec));
        if nm_diff != 0 {
            if let Some(old_nm) = aux_int(rec, b"NM") {
                aux_block_set(&mut aux, b"NM", Some(&encode_int((old_nm + nm_diff) as i32 as i64)));
            }
        }
        let old_insert_alignment_score = -BWA_GAP_OPEN_PENALTY * old_insert_count - BWA_GAP_EXTEND_PENALTY * old_inserted_bases;
        let new_insert_alignment_score = -BWA_GAP_OPEN_PENALTY * new_insert_count - BWA_GAP_EXTEND_PENALTY * new_inserted_bases;
        let alignment_score_diff = new_insert_alignment_score - old_insert_alignment_score;
        if alignment_score_diff != 0 {
            if let Some(old_as) = aux_int(rec, b"AS") {
                aux_block_set(&mut aux, b"AS", Some(&encode_int((old_as + alignment_score_diff) as i32 as i64)));
            }
        }
        set_aux_bytes(rec, &aux);
    }
    replace_alignment_data(rec, &new_cigar, Some((&read_bases, &read_quals)), &mut scratch.buf);

    // the first supplementary alignment's soft clip covers this read's aligned and inserted bases, so it loses
    // the stripped bases too
    if let Some(mut supp_list) = extract_full_supp_alignments(aux_str(rec, b"SA"))? {
        let read_base_length = |cigar: &[u32]| -> i32 { cigar.iter().filter(|e| matches!(op(**e), OP_M | OP_I)).map(|e| oplen(*e)).sum() };
        let old_read_base_length = read_base_length(&old_cigar);
        let new_read_base_length = read_base_length(&new_cigar);
        if old_read_base_length > new_read_base_length {
            if supp_list.is_empty() {
                return err("java.lang.IndexOutOfBoundsException");
            }
            supp_list[0].cigar = correct_supplementary_cigar(&supp_list[0].cigar, old_read_base_length - new_read_base_length)?;
            let mut tag = Vec::new();
            for (i, s) in supp_list.iter().enumerate() {
                if i > 0 {
                    tag.push(b';');
                }
                s.write_sam_tag(&mut tag);
            }
            set_attribute(rec, b"SA", Some(&encode_z(&tag)));
        }
    }
    Ok(())
}

/// SbxRoutines.ensureMirroredMismatchRepeatQuals: where a duplex mismatch quality is at one end of a repeat in an
/// aligned section, the other end of the repeat gets it too
pub fn ensure_mirrored_mismatch_repeat_quals(read_index_start: i32, element_length: i32, bases: &[u8], base_quals: &mut [u8]) -> SbxResult<bool> {
    let low = |q: u8| (q as i8) <= SBX_DUPLEX_MISMATCH_QUAL as i8;
    let mut changed = false;
    let read_index_end = read_index_start + element_length - 1;
    let mut read_index = read_index_start;
    while read_index <= read_index_end {
        if read_index as usize >= bases.len() {
            break;
        }
        if !low(base_quals[read_index as usize]) {
            read_index += 1;
            continue;
        }
        let mut repeat_length = 1;
        let mut is_homopolymer = true;
        let repeat_base = bases[read_index as usize];
        for i in read_index + 1..=read_index_end {
            if i as usize >= base_quals.len() {
                return err("java.lang.ArrayIndexOutOfBoundsException");
            }
            if !low(base_quals[i as usize]) {
                break;
            }
            if bases[i as usize] != repeat_base {
                is_homopolymer = false;
            }
            repeat_length += 1;
        }
        if is_homopolymer {
            repeat_length = 1;
        }
        let repeat = &bases[read_index as usize..(read_index + repeat_length) as usize];

        // the furthest whole repeat in either direction
        let mut mirrored_index_start = -1;
        let mut max_repeat_shift = 0;
        for search_up in [true, false] {
            let step = if search_up { repeat_length } else { -repeat_length };
            let mut repeat_start = read_index + step;
            let mut last_repeat_start_index = read_index;
            while repeat_start >= 0 && repeat_start + repeat_length - 1 <= read_index_end {
                if repeat_start as usize + repeat_length as usize > bases.len() {
                    return err("java.lang.ArrayIndexOutOfBoundsException");
                }
                if repeat != &bases[repeat_start as usize..(repeat_start + repeat_length) as usize] {
                    break;
                }
                last_repeat_start_index = repeat_start;
                repeat_start += step;
            }
            if last_repeat_start_index != read_index {
                let repeat_shift = (last_repeat_start_index - read_index).abs();
                if repeat_shift > max_repeat_shift {
                    max_repeat_shift = repeat_shift;
                    mirrored_index_start = last_repeat_start_index;
                }
            }
        }
        if mirrored_index_start >= 0 {
            for i in 0..repeat_length {
                let q = &mut base_quals[(mirrored_index_start + i) as usize];
                if *q != SBX_DUPLEX_MISMATCH_QUAL {
                    *q = SBX_DUPLEX_MISMATCH_QUAL;
                    changed = true;
                }
            }
        }
        if mirrored_index_start > read_index {
            read_index = mirrored_index_start + repeat_length - 1;
        }
        read_index += 1;
    }
    Ok(changed)
}

/// SbxRoutines.prepProcessRead. `decoded` says that a Java setter has already run on the record (the unmapper
/// removed its supplementary data).
///
/// Java changes the mirrored repeat qualities in the array that htsjdk decoded, without a setter. htsjdk writes
/// a record that never had a setter called from its stored bytes, so such a change reaches the BAM only when the
/// record is decoded for another reason or is later finalised. For a record that is not, the stored qualities
/// are kept in a leading private attribute until the record is either finalised (the change stays) or written
/// (`restore_stored_quals`).
pub fn prep_process_read(rec: &mut bam::Record, decoded: bool) {
    let mut decoded = decoded;
    let result = SCRATCH.with(|s| {
        let scratch = &mut *s.borrow_mut();
        strip_duplex_indels_with(rec, scratch, &mut decoded)?;
        // nothing to mirror without a duplex mismatch quality
        if !rec.qual().iter().any(|q| (*q as i8) <= SBX_DUPLEX_MISMATCH_QUAL as i8) {
            return Ok(());
        }
        decode_bases(rec, &mut scratch.bases);
        scratch.buf.clear();
        scratch.buf.extend_from_slice(rec.qual());
        let cigar_len = rec.raw_cigar().len();
        let mut read_index = 0;
        let mut changed = false;
        let mut result = Ok(());
        for i in 0..cigar_len {
            let e = rec.raw_cigar()[i];
            if op(e) == OP_M {
                match ensure_mirrored_mismatch_repeat_quals(read_index, oplen(e), &scratch.bases, quals_mut(rec)) {
                    Ok(c) => changed |= c,
                    Err(e) => {
                        changed = true;
                        result = Err(e);
                        break;
                    }
                }
            }
            if consumes_read(op(e)) {
                read_index += oplen(e);
            }
        }
        if changed && !decoded && rec.qual() != &scratch.buf[..] {
            let mut aux = Vec::with_capacity(8 + scratch.buf.len() + aux_bytes(rec).len());
            aux.extend_from_slice(&PENDING_QUALS_TAG);
            aux.extend_from_slice(b"BC");
            aux.extend_from_slice(&(scratch.buf.len() as u32).to_le_bytes());
            aux.extend_from_slice(&scratch.buf);
            aux.extend_from_slice(aux_bytes(rec));
            set_aux_bytes(rec, &aux);
        }
        result
    });
    set_decoded(rec, decoded);
    if let Err(e) = result {
        crate::error!("pre-process read({}) error: {}", read_to_string(rec), e.0);
    }
}

// ---------------------------------------------------------------------------------------------------------
// consensus

const DNA_BASE_BYTES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const DNA_N_BYTE: u8 = b'N';

fn find_most_common_base_count(base_counts: &[i32; 4]) -> i32 {
    base_counts.iter().copied().max().unwrap().max(0)
}

fn find_most_common_base(base_counts: &[i32; 4], ref_base: u8, max_count: i32) -> u8 {
    let mut max_base = DNA_N_BYTE;
    for b in 0..4 {
        if base_counts[b] == max_count && (DNA_BASE_BYTES[b] == ref_base || max_base == DNA_N_BYTE) {
            max_base = DNA_BASE_BYTES[b];
        }
    }
    max_base
}

/// SbxRoutines.determineBaseAndQual: the consensus of the duplex bases when there are any, else of the simplex
/// bases, needs more than half of the reads (those with any other quality count against it); without it the
/// reference base is taken with a quality that marks the position
pub fn determine_base_and_qual(location_bases: &[u8], location_quals: &[i32], chromosome: &str, position: i32, refs: &mut RefBases) -> BaseQualPair {
    let mut simplex_counts = [0i32; 4];
    let mut duplex_counts = [0i32; 4];
    let mut low_qual_count = 0;
    let mut simplex_count = 0;
    let mut duplex_count = 0;
    for i in 0..location_bases.len() {
        let qual = location_quals[i];
        let counts = if qual == RAW_SIMPLEX_QUAL as i32 {
            simplex_count += 1;
            &mut simplex_counts
        } else if qual == RAW_DUPLEX_QUAL as i32 {
            duplex_count += 1;
            &mut duplex_counts
        } else {
            low_qual_count += 1;
            continue;
        };
        if let Some(b) = DNA_BASE_BYTES.iter().position(|b| *b == location_bases[i]) {
            counts[b] += 1;
        }
    }
    let ref_base = if position != INVALID_POSITION { refs.ref_base(chromosome, position) } else { DNA_N_BYTE };
    if simplex_count == 0 && duplex_count == 0 {
        return BaseQualPair { base: ref_base, qual: SBX_DUPLEX_MISMATCH_QUAL as i32 };
    }
    if duplex_count == 0 {
        let max_base_count = find_most_common_base_count(&simplex_counts);
        let max_base = find_most_common_base(&simplex_counts, ref_base, max_base_count);
        let total_read_count = simplex_count + low_qual_count;
        return if max_base_count as f64 > SBX_CONSENSUS_BASE_THRESHOLD * total_read_count as f64 {
            BaseQualPair { base: max_base, qual: RAW_SIMPLEX_QUAL as i32 }
        } else {
            BaseQualPair { base: ref_base, qual: SIMPLEX_NO_CONSENSUS_QUAL as i32 }
        };
    }
    let max_base_count = find_most_common_base_count(&duplex_counts);
    let max_base = find_most_common_base(&duplex_counts, ref_base, max_base_count);
    if max_base_count as f64 > SBX_CONSENSUS_BASE_THRESHOLD * (duplex_count + low_qual_count) as f64 {
        BaseQualPair { base: max_base, qual: RAW_DUPLEX_QUAL as i32 }
    } else {
        BaseQualPair { base: ref_base, qual: DUPLEX_NO_CONSENSUS_QUAL as i32 }
    }
}

/// SbxRoutines.findMaxDuplexBaseIndex: the first index from the 5' end (of the first read) at which any read
/// has a quality other than simplex; -1 when there is none. `reads` gives each read's qualities.
pub fn find_max_duplex_base_index(first_is_negative_strand: bool, reads: &[&[u8]]) -> i32 {
    let from_start = !first_is_negative_strand;
    let max_read_index = reads.iter().map(|q| q.len()).max().unwrap_or(0) as i32;
    let mut index = if from_start { 0 } else { max_read_index - 1 };
    while index >= 0 && index < max_read_index {
        for quals in reads {
            if index as usize >= quals.len() {
                continue;
            }
            if quals[index as usize] != RAW_SIMPLEX_QUAL {
                return index;
            }
        }
        index += if from_start { 1 } else { -1 };
    }
    -1
}

// ---------------------------------------------------------------------------------------------------------
// finalising

const ADJACENT_INDICES: [i32; 6] = [-3, -2, -1, 1, 2, 3];

/// SbxRoutines.finaliseRead: the raw and internal base qualities become the final ones, bases with a duplex
/// mismatch quality in aligned sections become the reference base (with NM and AS following), X elements
/// become M, and the bases next to a duplex mismatch get lower qualities.
pub fn finalise_read(refs: &mut RefBases, names: &[String], rec: &mut bam::Record) {
    SCRATCH.with(|s| finalise_read_with(refs, names, rec, &mut s.borrow_mut()))
}

fn finalise_read_with(refs: &mut RefBases, names: &[String], rec: &mut bam::Record, scratch: &mut Scratch) {
    // a quality change that Java made in memory only stays, since the qualities are set here
    if has_pending_quals(rec) {
        remove_pending_quals(rec);
    }
    let read_length = rec.inner.core.l_qseq as usize;
    let last_read_index = read_length as i32 - 1;
    let unmapped = has(rec, FLAG_UNMAPPED);
    let negative_strand = has(rec, FLAG_REVERSE);

    // Java reads the reference over the whole alignment before changing the record, and gives the read up when
    // that fails (after having changed the qualities in memory, which is not reproduced)
    if !unmapped {
        let chromosome_length = if rec.tid() >= 0 { refs.chromosome_length(&names[rec.tid() as usize]) } else { 0 };
        let read_end = end(rec);
        if rec.tid() < 0 || start(rec) < 1 || read_end > chromosome_length || read_end < start(rec) {
            crate::error!("post-process read error: {}", read_to_string(rec));
            return;
        }
    }

    let attributes = finalise_attributes(aux_bytes(rec));
    if attributes.consensus_type == 3 {
        // ConsensusType.valueOf fails
        crate::error!("post-process read error: {}", read_to_string(rec));
        return;
    }
    let first_duplex_base_index: Option<i32> = if attributes.consensus_type == 0 {
        Some(find_max_duplex_base_index(negative_strand, &[rec.qual()]))
    } else {
        attributes.duplex_base_index.map(|v| v as i32)
    };
    let mut duplex_region_start = -1;
    let mut duplex_region_end = -1;
    let mut set_dual = false;
    if let Some(first) = first_duplex_base_index {
        if first >= 0 {
            // a read with any duplex base is DUAL; the duplex base index then tells the bases apart
            set_dual = attributes.consensus_type == 1;
            if negative_strand {
                duplex_region_start = 0;
                duplex_region_end = first;
            } else {
                duplex_region_start = first;
                duplex_region_end = last_read_index;
            }
        }
    }

    // the final qualities, noting the duplex mismatches whose neighbours are adjusted below
    scratch.mismatch_indices.clear();
    {
        let quals = quals_mut(rec);
        let has_region = duplex_region_start >= 0;
        for i in 0..read_length {
            let qual = quals[i];
            if qual == DUPLEX_NO_CONSENSUS_QUAL || (qual as i8) <= SBX_DUPLEX_MISMATCH_QUAL as i8 {
                scratch.mismatch_indices.push(i as i32);
            }
            let in_region = has_region && i as i32 >= duplex_region_start && i as i32 <= duplex_region_end;
            quals[i] = match qual {
                // a duplex quality outside the duplex region, or a simplex one inside it, came from a consensus
                RAW_DUPLEX_QUAL => {
                    if has_region && !in_region { SBX_DUPLEX_MISMATCH_QUAL } else { SBX_DUPLEX_QUAL }
                }
                RAW_SIMPLEX_QUAL => {
                    if in_region { SBX_DUPLEX_MISMATCH_QUAL } else { SBX_SIMPLEX_QUAL }
                }
                _ => SBX_DUPLEX_MISMATCH_QUAL,
            };
        }
    }

    let has_num_events_attribute = attributes.has_num_mutations;
    let mut new_nm: Option<i64> = None;
    let mut new_as: Option<i64> = None;
    scratch.mismatch_ref_bases.clear();

    if !unmapped {
        let chromosome = &names[rec.tid() as usize];
        let read_position_start = start(rec);
        let mut ref_pos = read_position_start;
        let mut read_index: usize = 0;
        let mut num_events: i64 = 0;
        let mut nm_diff: i64 = 0;
        let mut alignment_score_diff: i64 = 0;
        let mut requires_cigar_update = false;
        let cigar_count = rec.raw_cigar().len();
        // the qualities are not changed in this walk; the bases are
        let final_quals: &[u8] = unsafe { std::slice::from_raw_parts(rec.qual().as_ptr(), read_length) };
        for i in 0..cigar_count {
            let e = rec.raw_cigar()[i];
            let (mut o, n) = (op(e), oplen(e));
            if o == OP_X {
                requires_cigar_update = true;
                o = OP_M;
            }
            if o != OP_M {
                if consumes_read(o) {
                    read_index += n as usize;
                }
                if consumes_ref(o) {
                    ref_pos += n;
                }
                if !has_num_events_attribute && (o == OP_I || o == OP_D) {
                    num_events += n as i64;
                }
                continue;
            }
            for _ in 0..n {
                if read_index >= read_length {
                    break;
                }
                let low = final_quals[read_index] as i8 <= SBX_DUPLEX_MISMATCH_QUAL as i8;
                if has_num_events_attribute && !low {
                    read_index += 1;
                    ref_pos += 1;
                    continue;
                }
                let read_base = base_at(rec, read_index);
                let ref_base = refs.ref_base(chromosome, ref_pos);
                if !has_num_events_attribute && ref_base != read_base {
                    num_events += 1;
                }
                if low && ref_pos >= 1 {
                    scratch.mismatch_ref_bases.push((read_index as i32, ref_base));
                    if ref_base != read_base {
                        set_base(rec, read_index, ref_base);
                        nm_diff -= 1;
                        alignment_score_diff += BWA_MISMATCH_PENALTY + BWA_MATCH_SCORE;
                    }
                }
                read_index += 1;
                ref_pos += 1;
            }
        }
        if has_num_events_attribute {
            if nm_diff != 0 {
                let num_events = attributes.num_mutations.unwrap_or(0);
                new_nm = Some((num_events + nm_diff).max(0));
            }
        } else {
            new_nm = Some(num_events);
        }
        if alignment_score_diff != 0 {
            if let Some(old_as) = attributes.alignment_score {
                new_as = Some(old_as + alignment_score_diff);
            }
        }
        if requires_cigar_update {
            let mut cigar: Vec<u32> = rec.raw_cigar().iter().map(|e| if op(*e) == OP_X { el(oplen(*e), OP_M) } else { *e }).collect();
            let mut index = 0;
            while index + 1 < cigar.len() {
                let mut element = cigar[index];
                let next_index = index + 1;
                while next_index < cigar.len() {
                    let next_element = cigar[next_index];
                    if op(element) == OP_M && op(next_element) == OP_M {
                        element = el(oplen(element) + oplen(next_element), OP_M);
                        cigar[index] = element;
                        // List.remove(Object): the first element equal to the next one goes, which is an
                        // earlier element when one has the same length and operator
                        let first_equal = cigar.iter().position(|c| *c == next_element).unwrap();
                        cigar.remove(first_equal);
                    } else {
                        break;
                    }
                }
                index += 1;
            }
            replace_alignment_data(rec, &cigar, None, &mut scratch.buf);
        }
    }

    if !scratch.mismatch_indices.is_empty() {
        let has_region = duplex_region_start >= 0;
        for mi in 0..scratch.mismatch_indices.len() {
            let read_index = scratch.mismatch_indices[mi];
            let mismatch_ref_base =
                scratch.mismatch_ref_bases.binary_search_by_key(&read_index, |x| x.0).ok().map(|j| scratch.mismatch_ref_bases[j].1);
            for adjacent in ADJACENT_INDICES {
                let mut adjust_qual = if adjacent.abs() == 1 { SBX_DUPLEX_ADJACENT_1_QUAL } else { SBX_DUPLEX_ADJACENT_2_3_QUAL };
                let adjust_index = read_index + adjacent;
                if adjust_index < 0 || adjust_index > last_read_index {
                    continue;
                }
                // the base next to the mismatch that is the reference base at the mismatch gets a lower quality
                if adjust_qual == SBX_DUPLEX_ADJACENT_1_QUAL {
                    if let Some(ref_base) = mismatch_ref_base {
                        if base_at(rec, adjust_index as usize) == ref_base {
                            adjust_qual = SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH;
                        }
                    }
                }
                if has_region && (adjust_index < duplex_region_start || adjust_index > duplex_region_end) {
                    continue;
                }
                let q = &mut quals_mut(rec)[adjust_index as usize];
                if (*q as i8) > adjust_qual as i8 && *q != SBX_SIMPLEX_QUAL {
                    *q = adjust_qual;
                }
            }
        }
    }

    // setBaseQualities and setAttribute: the record is written in decoded form
    let aux = &mut scratch.buf;
    aux_sorted_into(aux_bytes(rec), aux);
    if set_dual {
        aux_block_set(aux, b"UT", Some(&encode_z(b"DUAL")));
    }
    if let Some(nm) = new_nm {
        aux_block_set(aux, b"NM", Some(&encode_int(nm as i32 as i64)));
    }
    if let Some(score) = new_as {
        aux_block_set(aux, b"AS", Some(&encode_int(score as i32 as i64)));
    }
    set_aux_bytes(rec, aux);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorted_attributes_match_general_routine() {
        // file order with integers in wider types than needed, as aligners write them
        let mut aux = Vec::new();
        aux.extend_from_slice(b"YCZ332+35H90+\0");
        aux.extend_from_slice(b"PGZMarkDuplicates\0");
        aux.extend_from_slice(b"ASi");
        aux.extend_from_slice(&453i32.to_le_bytes());
        aux.extend_from_slice(b"XSS");
        aux.extend_from_slice(&241u16.to_le_bytes());
        aux.extend_from_slice(b"MDZ325A132\0");
        aux.extend_from_slice(b"NMC\x01");
        aux.extend_from_slice(b"xbBs");
        aux.extend_from_slice(&2u32.to_le_bytes());
        aux.extend_from_slice(&[1, 0, 2, 0]);
        aux.extend_from_slice(b"xnc\xff");
        aux.extend_from_slice(b"xls");
        aux.extend_from_slice(&(-300i16).to_le_bytes());
        aux.extend_from_slice(b"xuI");
        aux.extend_from_slice(&4_000_000_000u32.to_le_bytes());
        let mut out = Vec::new();
        aux_sorted_into(&aux, &mut out);
        assert_eq!(aux_sorted_copy(&aux), out);
        // sorted input is unchanged
        let sorted = out.clone();
        aux_sorted_into(&sorted, &mut out);
        assert_eq!(sorted, out);
        // a repeated tag takes the general routine
        let mut repeated = aux.clone();
        repeated.extend_from_slice(b"NMC\x05");
        aux_sorted_into(&repeated, &mut out);
        assert_eq!(aux_sorted_copy(&repeated), out);
        let attributes = finalise_attributes(&aux);
        assert_eq!(0, attributes.consensus_type);
        assert_eq!(Some(1), attributes.num_mutations);
        assert_eq!(Some(453), attributes.alignment_score);
        assert_eq!(None, attributes.duplex_base_index);
    }
}
