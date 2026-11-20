//! Helpers over htslib records that follow htsjdk's SAMRecord semantics where REDUX depends on them.

use rust_htslib::bam::{self, record::Aux};
use rust_htslib::htslib;

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST: u16 = 0x40;
pub const FLAG_SECOND: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

pub const OP_M: u32 = 0;
pub const OP_I: u32 = 1;
pub const OP_D: u32 = 2;
pub const OP_N: u32 = 3;
pub const OP_S: u32 = 4;
pub const OP_H: u32 = 5;
pub const OP_P: u32 = 6;
pub const OP_EQ: u32 = 7;
pub const OP_X: u32 = 8;

pub const NO_POSITION: i32 = 0;

#[inline]
pub fn has(rec: &bam::Record, flag: u16) -> bool {
    rec.flags() & flag != 0
}

/// 1-based alignment start (0 when unplaced)
#[inline]
pub fn start(rec: &bam::Record) -> i32 {
    rec.pos() as i32 + 1
}

#[inline]
pub fn mate_start(rec: &bam::Record) -> i32 {
    rec.mpos() as i32 + 1
}

#[inline]
pub fn op(c: u32) -> u32 {
    c & 0xf
}
#[inline]
pub fn oplen(c: u32) -> i32 {
    (c >> 4) as i32
}
#[inline]
pub fn consumes_read(op: u32) -> bool {
    matches!(op, OP_M | OP_I | OP_S | OP_EQ | OP_X)
}
#[inline]
pub fn consumes_ref(op: u32) -> bool {
    matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X)
}

pub fn ref_length(cigar: &[u32]) -> i32 {
    cigar.iter().filter(|c| consumes_ref(op(**c))).map(|c| oplen(*c)).sum()
}

/// SAMRecord.getAlignmentEnd: 0 for an unmapped read, start - 1 for a mapped read without a cigar
#[inline]
pub fn end(rec: &bam::Record) -> i32 {
    if has(rec, FLAG_UNMAPPED) {
        return 0;
    }
    start(rec) + ref_length(rec.raw_cigar()) - 1
}

pub fn left_soft_clip(cigar: &[u32]) -> i32 {
    match cigar.first() {
        Some(c) if op(*c) == OP_S => oplen(*c),
        _ => 0,
    }
}
pub fn right_soft_clip(cigar: &[u32]) -> i32 {
    match cigar.last() {
        Some(c) if op(*c) == OP_S => oplen(*c),
        _ => 0,
    }
}

/// SamRecordUtils.getUnclippedPosition
pub fn unclipped_position(rec: &bam::Record, is_start: bool) -> i32 {
    if is_start { start(rec) - left_soft_clip(rec.raw_cigar()) } else { end(rec) + right_soft_clip(rec.raw_cigar()) }
}

const OP_CHARS: &[u8; 9] = b"MIDNSHP=X";

/// SAMRecord.getCigarString: "*" without a cigar
pub fn cigar_string(cigar: &[u32]) -> String {
    if cigar.is_empty() {
        return "*".to_string();
    }
    let mut s = String::with_capacity(cigar.len() * 4);
    for c in cigar {
        s.push_str(&oplen(*c).to_string());
        s.push(OP_CHARS[op(*c) as usize] as char);
    }
    s
}

/// SamRecordUtils.getReadBoundaryPosition over a cigar string
pub fn read_boundary_position(read_start: i32, cigar: &[u8], get_read_start: bool, include_soft_clipped: bool) -> i32 {
    if get_read_start && !include_soft_clipped {
        return read_start;
    }
    let mut current = read_start;
    let mut len: i32 = 0;
    for &c in cigar {
        if matches!(c, b'D' | b'M' | b'S' | b'N') {
            if get_read_start {
                return if c == b'S' { read_start - len } else { read_start };
            }
            if c == b'S' && read_start == current {
                // left clip, ignored for the read end
            } else if c == b'S' && !include_soft_clipped {
                break;
            } else {
                current += len;
            }
            len = 0;
            continue;
        }
        let digit = c as i32 - b'0' as i32;
        if (0..=9).contains(&digit) {
            len = len.wrapping_mul(10).wrapping_add(digit);
        } else {
            len = 0;
        }
    }
    current - 1
}

pub fn five_prime_unclipped_str(read_start: i32, cigar: &[u8], forward: bool) -> i32 {
    read_boundary_position(read_start, cigar, forward, true)
}
pub fn three_prime_unclipped_str(read_start: i32, cigar: &[u8], forward: bool) -> i32 {
    read_boundary_position(read_start, cigar, !forward, true)
}
pub fn five_prime_unclipped(rec: &bam::Record) -> i32 {
    unclipped_position(rec, !has(rec, FLAG_REVERSE))
}
pub fn three_prime_unclipped(rec: &bam::Record) -> i32 {
    unclipped_position(rec, has(rec, FLAG_REVERSE))
}

/// CigarUtils.calcCigarAlignedLength over a cigar string
pub fn calc_cigar_aligned_length(cigar: &[u8]) -> i32 {
    let mut total = 0;
    let mut len: i32 = 0;
    for &c in cigar {
        if matches!(c, b'D' | b'M' | b'N') {
            total += len;
            len = 0;
            continue;
        }
        let digit = c as i32 - b'0' as i32;
        if (0..=9).contains(&digit) {
            len = len * 10 + digit;
        } else {
            len = 0;
        }
    }
    total
}

/// A string (Z) attribute.
pub fn aux_str<'a>(rec: &'a bam::Record, tag: &[u8; 2]) -> Option<&'a [u8]> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.as_bytes()),
        _ => None,
    }
}

pub fn aux_int(rec: &bam::Record, tag: &[u8; 2]) -> Option<i64> {
    match rec.aux(tag) {
        Ok(Aux::I8(v)) => Some(v as i64),
        Ok(Aux::U8(v)) => Some(v as i64),
        Ok(Aux::I16(v)) => Some(v as i64),
        Ok(Aux::U16(v)) => Some(v as i64),
        Ok(Aux::I32(v)) => Some(v as i64),
        Ok(Aux::U32(v)) => Some(v as i64),
        _ => None,
    }
}

pub fn has_aux(rec: &bam::Record, tag: &[u8; 2]) -> bool {
    rec.aux(tag).is_ok()
}

fn aux_offset(b: &htslib::bam1_t) -> usize {
    let c = &b.core;
    c.l_qname as usize + c.n_cigar as usize * 4 + ((c.l_qseq as usize + 1) >> 1) + c.l_qseq as usize
}

/// The raw attribute block.
pub fn aux_bytes(rec: &bam::Record) -> &[u8] {
    let b = &rec.inner;
    let off = aux_offset(b);
    unsafe { std::slice::from_raw_parts(b.data.add(off), b.l_data as usize - off) }
}

unsafe extern "C" {
    fn realloc(p: *mut std::ffi::c_void, size: usize) -> *mut std::ffi::c_void;
}

/// Replaces the whole attribute block.
pub fn set_aux_bytes(rec: &mut bam::Record, aux: &[u8]) {
    let b = &mut rec.inner;
    let off = aux_offset(b);
    let new_len = off + aux.len();
    unsafe {
        if new_len > b.m_data as usize {
            let p = realloc(b.data as *mut std::ffi::c_void, new_len) as *mut u8;
            assert!(!p.is_null());
            b.data = p;
            b.m_data = new_len as u32;
        }
        std::ptr::copy_nonoverlapping(aux.as_ptr(), b.data.add(off), aux.len());
    }
    b.l_data = new_len as i32;
}

/// Length in bytes of the attribute that starts at `aux[0]` (tag, type, value).
fn aux_field_len(aux: &[u8]) -> usize {
    let t = aux[2];
    let fixed = |t: u8| match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        _ => 0,
    };
    match t {
        b'Z' | b'H' => 3 + aux[3..].iter().position(|&c| c == 0).unwrap() + 1,
        b'B' => {
            let n = u32::from_le_bytes([aux[4], aux[5], aux[6], aux[7]]) as usize;
            8 + n * fixed(aux[3])
        }
        _ => 3 + fixed(t),
    }
}

/// The attributes as (start, end) spans of the block.
pub fn aux_spans(aux: &[u8]) -> Vec<(usize, usize)> {
    let mut v = Vec::new();
    let mut i = 0;
    while i + 3 <= aux.len() {
        let n = aux_field_len(&aux[i..]);
        v.push((i, i + n));
        i += n;
    }
    v
}

/// SAMTag.makeBinaryTag: the second character is the high byte.
#[inline]
fn java_tag(tag: &[u8]) -> u16 {
    (tag[1] as u16) << 8 | tag[0] as u16
}

/// An attribute value encoded as type byte plus payload.
pub fn encode_z(value: &[u8]) -> Vec<u8> {
    let mut v = Vec::with_capacity(value.len() + 2);
    v.push(b'Z');
    v.extend_from_slice(value);
    v.push(0);
    v
}

/// htsjdk BinaryTagCodec.getIntegerType: the smallest type that holds the value
pub fn encode_int(val: i64) -> Vec<u8> {
    if val > i32::MAX as i64 {
        let mut v = vec![b'I'];
        v.extend_from_slice(&(val as u32).to_le_bytes());
        return v;
    }
    if val >= 0 {
        if val <= 127 {
            vec![b'c', val as u8]
        } else if val <= 255 {
            vec![b'C', val as u8]
        } else if val <= 32767 {
            let mut v = vec![b's'];
            v.extend_from_slice(&(val as i16).to_le_bytes());
            v
        } else if val <= 65535 {
            let mut v = vec![b'S'];
            v.extend_from_slice(&(val as u16).to_le_bytes());
            v
        } else {
            let mut v = vec![b'i'];
            v.extend_from_slice(&(val as i32).to_le_bytes());
            v
        }
    } else if val >= -128 {
        vec![b'c', val as i8 as u8]
    } else if val >= -32768 {
        let mut v = vec![b's'];
        v.extend_from_slice(&(val as i16).to_le_bytes());
        v
    } else {
        let mut v = vec![b'i'];
        v.extend_from_slice(&(val as i32).to_le_bytes());
        v
    }
}

/// SAMRecord.setAttribute on an attribute block, as SAMBinaryTagAndValue.insert does it: the new attribute goes
/// before the first one with a greater binary tag, replaces one with the same tag, or goes at the end. `None`
/// removes the attribute.
pub fn aux_block_set(aux: &mut Vec<u8>, tag: &[u8; 2], value: Option<&[u8]>) {
    let key = java_tag(tag);
    let mut i = 0;
    let mut at: Option<(usize, usize)> = None;
    while i + 3 <= aux.len() {
        let n = aux_field_len(&aux[i..]);
        let k = java_tag(&aux[i..i + 2]);
        match value {
            None => {
                if k == key {
                    aux.drain(i..i + n);
                    return;
                }
            }
            Some(_) => {
                if key < k {
                    at = Some((i, i));
                    break;
                } else if key == k {
                    at = Some((i, i + n));
                    break;
                }
            }
        }
        i += n;
    }
    let Some(value) = value else { return };
    let (s, e) = at.unwrap_or((aux.len(), aux.len()));
    let mut field = Vec::with_capacity(2 + value.len());
    field.extend_from_slice(tag);
    field.extend_from_slice(value);
    aux.splice(s..e, field);
}

/// The integer value of an attribute field (tag, type, value), when it is an integer.
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

/// The attribute block as htsjdk writes it once it has decoded the attributes: sorted by binary tag (a repeated
/// tag keeps its last value) and with each integer in the smallest type that holds it. htsjdk writes a record's
/// stored bytes unchanged until its name, bases, qualities, cigar or attributes are set; from then on it
/// writes the decoded form.
pub fn aux_sorted_copy(aux: &[u8]) -> Vec<u8> {
    let mut seen: Vec<(u16, (usize, usize))> = Vec::with_capacity(16);
    let mut i = 0;
    let mut sorted = true;
    while i + 3 <= aux.len() {
        let n = aux_field_len(&aux[i..]);
        let k = java_tag(&aux[i..i + 2]);
        match seen.last() {
            Some(last) if k <= last.0 => {
                sorted = false;
                match seen.binary_search_by(|e| e.0.cmp(&k)) {
                    Ok(j) => seen[j].1 = (i, i + n),
                    Err(j) => seen.insert(j, (k, (i, i + n))),
                }
            }
            _ => seen.push((k, (i, i + n))),
        }
        i += n;
    }
    let _ = sorted;
    let mut out = Vec::with_capacity(aux.len());
    for (_, (s, e)) in seen {
        let f = &aux[s..e];
        match aux_field_int(f) {
            Some(v) => {
                out.extend_from_slice(&f[..2]);
                out.extend_from_slice(&encode_int(v));
            }
            None => out.extend_from_slice(f),
        }
    }
    out
}

/// Puts a record's attributes in htsjdk's decoded form; see `aux_sorted_copy`.
pub fn decode_attributes(rec: &mut bam::Record) {
    let aux = aux_sorted_copy(aux_bytes(rec));
    set_aux_bytes(rec, &aux);
}

/// SAMRecord.setAttribute on a record read from a BAM: the attributes take their decoded form and the new value
/// its sorted place. Removing an attribute that is absent changes nothing (the Java callers test for it first).
pub fn set_attribute(rec: &mut bam::Record, tag: &[u8; 2], value: Option<&[u8]>) {
    if value.is_none() && !has_aux(rec, tag) {
        return;
    }
    let mut aux = aux_sorted_copy(aux_bytes(rec));
    aux_block_set(&mut aux, tag, value);
    set_aux_bytes(rec, &aux);
}

/// Sets the BAM bin from the position and cigar, as writers expect after the alignment changes.
pub fn fix_bin(rec: &mut bam::Record) {
    let beg = rec.pos();
    let mut end = if has(rec, FLAG_UNMAPPED) { beg + 1 } else { beg + ref_length(rec.raw_cigar()) as i64 };
    if end <= beg {
        end = beg + 1;
    }
    rec.inner.core.bin = reg2bin(beg, end);
}

fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    if beg < 0 {
        return 4680;
    }
    if beg >> 14 == end >> 14 {
        return (((1 << 15) - 1) / 7 + (beg >> 14)) as u16;
    }
    if beg >> 17 == end >> 17 {
        return (((1 << 12) - 1) / 7 + (beg >> 17)) as u16;
    }
    if beg >> 20 == end >> 20 {
        return (((1 << 9) - 1) / 7 + (beg >> 20)) as u16;
    }
    if beg >> 23 == end >> 23 {
        return (((1 << 6) - 1) / 7 + (beg >> 23)) as u16;
    }
    if beg >> 26 == end >> 26 {
        return (((1 << 3) - 1) / 7 + (beg >> 26)) as u16;
    }
    0
}

/// Sets or clears one flag.
#[inline]
pub fn set_flag(rec: &mut bam::Record, flag: u16, on: bool) {
    let f = rec.flags();
    rec.set_flags(if on { f | flag } else { f & !flag });
}

/// Removes the cigar (SAMRecord.setCigarString("*")).
pub fn clear_cigar(rec: &mut bam::Record) {
    let b = &mut rec.inner;
    let n = b.core.n_cigar as usize * 4;
    if n == 0 {
        return;
    }
    let off = b.core.l_qname as usize;
    unsafe {
        std::ptr::copy(b.data.add(off + n), b.data.add(off), b.l_data as usize - off - n);
    }
    b.l_data -= n as i32;
    b.core.n_cigar = 0;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_positions() {
        assert_eq!(five_prime_unclipped_str(100, b"10S90M", true), 90);
        assert_eq!(five_prime_unclipped_str(100, b"90M10S", false), 199);
        assert_eq!(five_prime_unclipped_str(100, b"5H10S50M", true), 90);
        assert_eq!(three_prime_unclipped_str(100, b"10S90M", true), 189);
        assert_eq!(read_boundary_position(100, b"10S50M5D40M3S", false, false), 194);
    }

    #[test]
    fn attribute_insert_order() {
        // file order NM, MD, AS
        let mut aux = Vec::new();
        for (t, v) in [(b"NM", encode_int(1)), (b"MD", encode_z(b"100")), (b"AS", encode_int(100))] {
            aux.extend_from_slice(t);
            aux.extend_from_slice(&v);
        }
        let mut aux = aux_sorted_copy(&aux);
        aux_block_set(&mut aux, b"UM", Some(&encode_z(b"chr1:5")));
        let tags: Vec<&[u8]> = aux_spans(&aux).iter().map(|(s, _)| &aux[*s..*s + 2]).collect();
        // by binary tag: MD 0x444D, NM 0x4D4E, UM 0x4D55, AS 0x5341
        assert_eq!(tags, vec![&b"MD"[..], b"NM", b"UM", b"AS"]);
        aux_block_set(&mut aux, b"MD", None);
        assert_eq!(aux_spans(&aux).len(), 3);
    }
}
