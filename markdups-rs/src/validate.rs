//! htsjdk's checks of every record read under VALIDATION_STRINGENCY STRICT (MarkDuplicates' default):
//! SAMRecord.isValid (the first error) as BAMFileReader's iterator runs it, then the CIGAR checks
//! (Cigar.isValid, SAMUtils.validateCigar) that BAMRecord.getCigar runs when it first decodes the
//! CIGAR of a mapped read. Under STRICT an error stops the program as it stops GATK; LENIENT warns.

use anyhow::{Result, bail};

use crate::bam::{CIGAR_D, CIGAR_EQ, CIGAR_H, CIGAR_M, CIGAR_N, CIGAR_S, CIGAR_X, Rec, SamHeader};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Stringency {
    Strict,
    Lenient,
    Silent,
}

impl Stringency {
    pub fn parse(s: &str) -> Stringency {
        match s {
            "LENIENT" => Stringency::Lenient,
            "SILENT" => Stringency::Silent,
            _ => Stringency::Strict,
        }
    }
}

const CIGAR_I: u32 = 1;
const CIGAR_P: u32 = 6;

fn is_real(op: u32) -> bool {
    matches!(op, CIGAR_M | CIGAR_I | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X)
}
fn is_indel(op: u32) -> bool {
    matches!(op, CIGAR_I | CIGAR_D)
}

/// One validation error: its type and message.
type Error = (&'static str, String);

/// SAMRecord.isValid(firstOnly = true)
fn record_error(rec: &Rec, h: &SamHeader, rg: Option<&[u8]>) -> Option<Error> {
    let flag = rec.flag();
    let has_ref = |r: i32| r != -1;
    let alignment_start = rec.pos() + 1;
    let mate_start = rec.mate_pos() + 1;
    if !rec.is_paired() {
        let checks: [(u16, &str, &str); 5] = [
            (0x2, "INVALID_FLAG_PROPER_PAIR", "Proper pair flag should not be set for unpaired read."),
            (0x8, "INVALID_FLAG_MATE_UNMAPPED", "Mate unmapped flag should not be set for unpaired read."),
            (0x20, "INVALID_FLAG_MATE_NEG_STRAND", "Mate negative strand flag should not be set for unpaired read."),
            (0x40, "INVALID_FLAG_FIRST_OF_PAIR", "First of pair flag should not be set for unpaired read."),
            (0x80, "INVALID_FLAG_SECOND_OF_PAIR", "Second of pair flag should not be set for unpaired read."),
        ];
        for (bit, ty, msg) in checks {
            if flag & bit != 0 {
                return Some((ty, msg.into()));
            }
        }
        if rec.mate_ref_id() != -1 {
            return Some(("INVALID_MATE_REF_INDEX", "MRNM should not be set for unpaired read.".into()));
        }
    } else {
        if let Some(e) = position_error(rec.mate_ref_id(), mate_start, true, h) {
            return Some(e);
        }
        if !has_ref(rec.mate_ref_id()) && !rec.mate_unmapped() {
            return Some(("INVALID_FLAG_MATE_UNMAPPED", "Mapped mate should have mate reference name".into()));
        }
        if flag & (0x40 | 0x80) == 0 {
            return Some(("PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND", "Paired read should be marked as first of pair or second of pair.".into()));
        }
    }
    if rec.is_unmapped() {
        if flag & 0x100 != 0 {
            return Some(("INVALID_FLAG_NOT_PRIM_ALIGNMENT", "Secondary alignment flag should not be set for unmapped read.".into()));
        }
        if flag & 0x800 != 0 {
            return Some(("INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT", "Supplementary alignment flag should not be set for unmapped read.".into()));
        }
        if rec.mapq() != 0 {
            return Some(("INVALID_MAPPING_QUALITY", "MAPQ should be 0 for unmapped read.".into()));
        }
    } else {
        if rec.n_cigar() == 0 {
            return Some(("INVALID_CIGAR", "CIGAR should have > zero elements for mapped read.".into()));
        }
        if h.sq.is_empty() {
            return Some(("MISSING_SEQUENCE_DICTIONARY", "Empty sequence dictionary.".into()));
        }
        if !has_ref(rec.ref_id()) {
            return Some(("INVALID_FLAG_READ_UNMAPPED", "Mapped read should have valid reference name".into()));
        }
    }
    if let Some(rg) = rg {
        if !h.rg.iter().any(|r| r.id == rg) {
            return Some(("READ_GROUP_NOT_FOUND", format!("RG ID on SAMRecord not found in header: {}", String::from_utf8_lossy(rg))));
        }
    }
    if let Some(e) = position_error(rec.ref_id(), alignment_start, false, h) {
        return Some(e);
    }
    let q = rec.qual();
    let null_quals = rec.l_seq() == 0 || q[0] == 0xff;
    if !has_ref(rec.mate_ref_id()) && mate_start != 0 {
        return Some(("INVALID_UNALIGNED_MATE_START", format!("The unaligned mate start position is {alignment_start}, should be 0")));
    }
    let _ = null_quals;
    None
}

/// The checks of isValid after the CIGAR is decoded (which runs the CIGAR checks first)
fn length_error(rec: &Rec) -> Option<Error> {
    let ops = rec.cigar_ops();
    let cigar_read_len: u32 = ops.iter().filter(|&&o| matches!(o & 0xf, CIGAR_M | CIGAR_I | CIGAR_S | CIGAR_EQ | CIGAR_X)).map(|&o| o >> 4).sum();
    if rec.l_seq() != 0 && cigar_read_len != 0 && cigar_read_len as usize != rec.l_seq() {
        return Some(("MISMATCH_CIGAR_SEQ_LENGTH", format!("CIGAR covers {cigar_read_len} bases but the sequence is {} read bases ", rec.l_seq())));
    }
    None
}

/// SAMValidationError.Type severities that are not ERROR
fn severity(ty: &str) -> &'static str {
    if matches!(ty, "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND" | "ADJACENT_INDEL_IN_CIGAR") { "WARNING" } else { "ERROR" }
}

/// isValidReferenceIndexAndPosition
fn position_error(ref_id: i32, start: i32, mate: bool, h: &SamHeader) -> Option<Error> {
    let msg = |m: &str| if mate { format!("Mate {}{}", m[..1].to_lowercase(), &m[1..]) } else { m.to_string() };
    if ref_id == -1 {
        if start != 0 {
            return Some(("INVALID_ALIGNMENT_START", msg("Alignment start should be 0 because reference name = *.")));
        }
    } else {
        if start == 0 {
            return Some(("INVALID_ALIGNMENT_START", msg("Alignment start should != 0 because reference name != *.")));
        }
        if !h.sq.is_empty() {
            match h.sq.get(ref_id as usize) {
                None => return Some(("INVALID_REFERENCE_INDEX", msg("Reference sequence not found in sequence dictionary."))),
                Some(s) => {
                    if start > s.len {
                        return Some((
                            "INVALID_ALIGNMENT_START",
                            msg(&format!("Alignment start ({start}) must be <= reference sequence length ({}) on reference {}", s.len, String::from_utf8_lossy(&s.name))),
                        ));
                    }
                }
            }
        }
    }
    None
}

/// Cigar.isValid and the off-the-end check of SAMUtils.validateCigar (the first error).
fn cigar_error(rec: &Rec, h: &SamHeader) -> Option<Error> {
    if rec.is_unmapped() {
        return None;
    }
    let ops = rec.cigar_ops();
    if ops.is_empty() {
        return None;
    }
    let n = ops.len();
    let mut real = false;
    for (i, &o) in ops.iter().enumerate() {
        let op = o & 0xf;
        if o >> 4 == 0 {
            return Some(("INVALID_CIGAR", "CIGAR element with zero length".into()));
        }
        if op == CIGAR_H {
            if i != 0 && i != n - 1 {
                return Some(("INVALID_CIGAR", "Hard clipping operator not at start or end of CIGAR".into()));
            }
        } else if op == CIGAR_S {
            if i == 0 || i == n - 1 {
            } else if i == 1 {
                if !(n == 3 && ops[2] & 0xf == CIGAR_H) && ops[0] & 0xf != CIGAR_H {
                    return Some(("INVALID_CIGAR", "Soft clipping CIGAR operator can only be inside of hard clipping operator".into()));
                }
            } else if i == n - 2 {
                if ops[n - 1] & 0xf != CIGAR_H {
                    return Some(("INVALID_CIGAR", "Soft clipping CIGAR operator can only be inside of hard clipping operator".into()));
                }
            } else {
                return Some(("INVALID_CIGAR", "Soft clipping CIGAR operator can at start or end of read, or be inside of hard clipping operator".into()));
            }
        } else if is_real(op) {
            real = true;
            if is_indel(op) {
                for &p in &ops[i + 1..] {
                    let next = p & 0xf;
                    if (is_real(next) && !is_indel(next)) || next == CIGAR_P {
                        break;
                    }
                    if is_indel(next) && op == next {
                        let name = if op == CIGAR_I { "I" } else { "D" };
                        return Some(("ADJACENT_INDEL_IN_CIGAR", format!("No M or N operator between pair of {name} operators in CIGAR")));
                    }
                }
            }
        } else if op == CIGAR_P {
            if i == n - 1 {
                return Some(("INVALID_CIGAR", "Padding operator not valid at end of CIGAR".into()));
            } else if i > 0 && (!is_real(ops[i - 1] & 0xf) || !is_real(ops[i + 1] & 0xf)) {
                return Some(("INVALID_CIGAR", "Padding operator not between real operators in CIGAR".into()));
            }
        }
    }
    if !real {
        return Some(("INVALID_CIGAR", "No real operator (M|I|D|N) in CIGAR".into()));
    }
    // alignment blocks (M, =, X) past the end of the reference
    if rec.ref_id() >= 0 {
        if let Some(s) = h.sq.get(rec.ref_id() as usize) {
            let mut ref_pos = rec.pos() + 1;
            for &o in &ops {
                let (op, len) = (o & 0xf, (o >> 4) as i32);
                match op {
                    CIGAR_M | CIGAR_EQ | CIGAR_X => {
                        if ref_pos + len - 1 > s.len {
                            return Some(("CIGAR_MAPS_OFF_REFERENCE", "Read CIGAR M operator maps off end of reference".into()));
                        }
                        ref_pos += len;
                    }
                    CIGAR_D | CIGAR_N => ref_pos += len,
                    _ => {}
                }
            }
        }
    }
    None
}

/// Validates a record as htsjdk reads it; `number` is its 1-based number in the file.
pub fn check(rec: &Rec, h: &SamHeader, rg: Option<&[u8]>, number: u64, s: Stringency) -> Result<()> {
    if s == Stringency::Silent {
        return Ok(());
    }
    let name = NameDisplay(rec.name());
    let report = |text: String| -> Result<()> {
        if s == Stringency::Strict {
            bail!("htsjdk.samtools.SAMFormatException: SAM validation error: {text}");
        }
        eprintln!("Ignoring SAM validation error: {text}");
        Ok(())
    };
    // isValid stops at its first error; the CIGAR is decoded (and checked, without a record number) at
    // its CIGAR length check
    if let Some((ty, msg)) = record_error(rec, h, rg) {
        report(format!("{}::{ty}:Record {number}, Read name {name}, {msg}", severity(ty)))?;
        if let Some((ty, msg)) = cigar_error(rec, h) {
            report(format!("{}::{ty}:Read name {name}, {msg}", severity(ty)))?;
        }
        return Ok(());
    }
    if let Some((ty, msg)) = cigar_error(rec, h) {
        report(format!("{}::{ty}:Read name {name}, {msg}", severity(ty)))?;
    }
    if let Some((ty, msg)) = length_error(rec) {
        report(format!("{}::{ty}:Record {number}, Read name {name}, {msg}", severity(ty)))?;
    }
    Ok(())
}

/// A read name printed only when an error is reported
struct NameDisplay<'a>(&'a [u8]);

impl std::fmt::Display for NameDisplay<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(self.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::SqRecord;

    fn header() -> SamHeader {
        let mut h = SamHeader::default();
        h.sq.push(SqRecord { name: b"chr1".to_vec(), len: 1000, attrs: Vec::new() });
        h.rg.push(crate::bam::IdRecord { id: b"rg1".to_vec(), attrs: Vec::new() });
        h
    }

    fn rec(flag: u16, tid: i32, pos: i32, mtid: i32, mpos: i32, cigar: &[u32], mapq: u8, aux: &[u8]) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(&tid.to_le_bytes());
        b.extend_from_slice(&pos.to_le_bytes());
        b.push(2);
        b.push(mapq);
        b.extend_from_slice(&0u16.to_le_bytes());
        b.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
        b.extend_from_slice(&flag.to_le_bytes());
        let l: u32 = cigar.iter().filter(|&&o| matches!(o & 0xf, 0 | 1 | 4)).map(|&o| o >> 4).sum();
        b.extend_from_slice(&(l as i32).to_le_bytes());
        b.extend_from_slice(&mtid.to_le_bytes());
        b.extend_from_slice(&mpos.to_le_bytes());
        b.extend_from_slice(&0i32.to_le_bytes());
        b.extend_from_slice(b"r\0");
        for c in cigar {
            b.extend_from_slice(&c.to_le_bytes());
        }
        b.extend(std::iter::repeat_n(0x11u8, (l as usize).div_ceil(2)));
        b.extend(std::iter::repeat_n(30u8, l as usize));
        b.extend_from_slice(aux);
        b
    }

    #[test]
    fn strict_errors_as_htsjdk() {
        let h = header();
        let ok = rec(0, 0, 10, -1, -1, &[50 << 4], 60, b"RGZrg1\0");
        assert!(check(&Rec(&ok), &h, Rec(&ok).rg_pg().unwrap().0, 1, Stringency::Strict).is_ok());
        let bad_rg = rec(0, 0, 10, -1, -1, &[50 << 4], 60, b"RGZrgZ\0");
        let e = check(&Rec(&bad_rg), &h, Rec(&bad_rg).rg_pg().unwrap().0, 40, Stringency::Strict).unwrap_err().to_string();
        assert!(e.ends_with("ERROR::READ_GROUP_NOT_FOUND:Record 40, Read name r, RG ID on SAMRecord not found in header: rgZ"), "{e}");
        assert!(check(&Rec(&bad_rg), &h, Rec(&bad_rg).rg_pg().unwrap().0, 40, Stringency::Lenient).is_ok());
        let unmapped_mapq = rec(4, 0, 10, -1, -1, &[], 60, b"");
        assert!(check(&Rec(&unmapped_mapq), &h, Rec(&unmapped_mapq).rg_pg().unwrap().0, 1, Stringency::Strict).unwrap_err().to_string().contains("MAPQ should be 0 for unmapped read."));
        let off_end = rec(0, 0, 990, -1, -1, &[50 << 4], 60, b"");
        assert!(check(&Rec(&off_end), &h, Rec(&off_end).rg_pg().unwrap().0, 1, Stringency::Strict).unwrap_err().to_string().contains("CIGAR_MAPS_OFF_REFERENCE"));
        let mid_soft = rec(0, 0, 10, -1, -1, &[10 << 4, 5 << 4 | 4, 10 << 4], 60, b"");
        assert!(check(&Rec(&mid_soft), &h, Rec(&mid_soft).rg_pg().unwrap().0, 1, Stringency::Strict).is_err());
        let unpaired_mate_flag = rec(0x8, 0, 10, -1, -1, &[50 << 4], 60, b"");
        assert!(check(&Rec(&unpaired_mate_flag), &h, Rec(&unpaired_mate_flag).rg_pg().unwrap().0, 1, Stringency::Strict).is_err());
    }
}
