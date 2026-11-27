//! A BAM record as the metrics code reads it: the fixed fields, the name, the CIGAR and qualities, and
//! the tags it needs (CR, UT, SA, MC) found in one pass over the aux data. Nothing is allocated per read.

use rust_htslib::bam;

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST: u16 = 0x40;
pub const FLAG_SECOND: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_QC_FAIL: u16 = 0x200;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

pub const OP_M: u32 = 0;
pub const OP_I: u32 = 1;
pub const OP_D: u32 = 2;
pub const OP_N: u32 = 3;
pub const OP_S: u32 = 4;
pub const OP_EQ: u32 = 7;
pub const OP_X: u32 = 8;

/// ConsensusType from the UT tag
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConsensusType {
    None,
    Single,
    Dual,
}

/// A view of one record
pub struct Rec<'a> {
    pub tid: i32,
    /// 1-based alignment start (htsjdk getAlignmentStart)
    pub start: i32,
    pub flags: u16,
    pub mapq: u8,
    pub mtid: i32,
    /// 1-based mate start
    pub mate_start: i32,
    pub tlen: i32,
    pub name: &'a [u8],
    pub cigar: &'a [u32],
    pub seq_len: usize,
    pub quals: &'a [u8],
    aux: &'a [u8],
    /// htsjdk getAlignmentEnd: start + reference length - 1, 0 for an unmapped read
    pub end: i32,
    pub has_cr: bool,
    ut: Option<&'a [u8]>,
    sa: Option<&'a [u8]>,
    mc: Option<&'a [u8]>,
}

/// The reference bases a CIGAR covers (M, D, N, = and X), as htsjdk's getReferenceLength
#[inline]
pub fn reference_length(cigar: &[u32]) -> i32 {
    let mut len = 0i32;
    for &c in cigar {
        let op = c & 0xf;
        if op == OP_M || op == OP_D || op == OP_N || op == OP_EQ || op == OP_X {
            len += (c >> 4) as i32;
        }
    }
    len
}

impl<'a> Rec<'a> {
    pub fn of(record: &'a bam::Record) -> Rec<'a> {
        let b = record.inner();
        let c = &b.core;
        // SAFETY: htslib guarantees `data` holds `l_data` bytes for as long as the record is borrowed
        let data: &'a [u8] = if b.l_data > 0 { unsafe { std::slice::from_raw_parts(b.data, b.l_data as usize) } } else { &[] };
        let l_qname = c.l_qname as usize;
        let name_len = l_qname.saturating_sub(1 + c.l_extranul as usize);
        let n_cigar = c.n_cigar as usize;
        let cigar_bytes = &data[l_qname..l_qname + n_cigar * 4];
        // SAFETY: htslib aligns the CIGAR to 4 bytes (l_extranul pads the name); the slice is n_cigar u32s
        let cigar: &'a [u32] = unsafe { std::slice::from_raw_parts(cigar_bytes.as_ptr() as *const u32, n_cigar) };
        let seq_len = c.l_qseq.max(0) as usize;
        let seq_off = l_qname + n_cigar * 4;
        let qual_off = seq_off + seq_len.div_ceil(2);
        let aux_off = qual_off + seq_len;
        let quals = &data[qual_off..aux_off];
        let aux = &data[aux_off.min(data.len())..];
        let flags = c.flag;
        let start = c.pos as i32 + 1;
        let end = if flags & FLAG_UNMAPPED != 0 { 0 } else { start + reference_length(cigar) - 1 };
        let mut r = Rec {
            tid: c.tid,
            start,
            flags,
            mapq: c.qual,
            mtid: c.mtid,
            mate_start: c.mpos as i32 + 1,
            tlen: c.isize_ as i32,
            name: &data[..name_len],
            cigar,
            seq_len,
            quals,
            aux,
            end,
            has_cr: false,
            ut: None,
            sa: None,
            mc: None,
        };
        r.scan_aux();
        r
    }

    /// one pass over the aux data for CR, UT, SA and MC
    fn scan_aux(&mut self) {
        let a = self.aux;
        let mut i = 0usize;
        while i + 3 <= a.len() {
            let (t0, t1, ty) = (a[i], a[i + 1], a[i + 2]);
            let v = i + 3;
            let next = match ty {
                b'A' | b'c' | b'C' => v + 1,
                b's' | b'S' => v + 2,
                b'i' | b'I' | b'f' => v + 4,
                b'd' => v + 8,
                b'Z' | b'H' => {
                    let j = match memchr::memchr(0, &a[v.min(a.len())..]) {
                        Some(k) => v + k,
                        None => a.len(),
                    };
                    let value = &a[v..j];
                    match (t0, t1) {
                        (b'U', b'T') => self.ut = Some(value),
                        (b'S', b'A') => self.sa = Some(value),
                        (b'M', b'C') => self.mc = Some(value),
                        _ => {}
                    }
                    j + 1
                }
                b'B' => {
                    if v + 5 > a.len() {
                        break;
                    }
                    let sub = a[v];
                    let n = u32::from_le_bytes([a[v + 1], a[v + 2], a[v + 3], a[v + 4]]) as usize;
                    let size = match sub {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        _ => 4,
                    };
                    v + 5 + n * size
                }
                _ => break,
            };
            if t0 == b'C' && t1 == b'R' {
                self.has_cr = true;
            }
            i = next;
        }
    }

    #[inline]
    pub fn flag(&self, f: u16) -> bool {
        self.flags & f != 0
    }
    #[inline]
    pub fn paired(&self) -> bool {
        self.flag(FLAG_PAIRED)
    }
    #[inline]
    pub fn unmapped(&self) -> bool {
        self.flag(FLAG_UNMAPPED)
    }
    /// getMateUnmappedFlag (only asked of paired reads)
    #[inline]
    pub fn mate_unmapped(&self) -> bool {
        self.flag(FLAG_MATE_UNMAPPED)
    }
    #[inline]
    pub fn reverse(&self) -> bool {
        self.flag(FLAG_REVERSE)
    }
    #[inline]
    pub fn mate_reverse(&self) -> bool {
        self.flag(FLAG_MATE_REVERSE)
    }
    #[inline]
    pub fn duplicate(&self) -> bool {
        self.flag(FLAG_DUPLICATE)
    }
    #[inline]
    pub fn secondary(&self) -> bool {
        self.flag(FLAG_SECONDARY)
    }
    #[inline]
    pub fn supplementary(&self) -> bool {
        self.flag(FLAG_SUPPLEMENTARY)
    }
    #[inline]
    pub fn first_of_pair(&self) -> bool {
        self.flag(FLAG_FIRST)
    }
    #[inline]
    pub fn second_of_pair(&self) -> bool {
        self.flag(FLAG_SECOND)
    }
    #[inline]
    pub fn qc_fail(&self) -> bool {
        self.flag(FLAG_QC_FAIL)
    }

    /// SamRecordUtils.extractConsensusType: from the UT tag; Java's valueOf throws for other values
    pub fn consensus_type(&self) -> Result<ConsensusType, String> {
        match self.ut {
            None => Ok(ConsensusType::None),
            Some(b"NONE") => Ok(ConsensusType::None),
            Some(b"SINGLE") => Ok(ConsensusType::Single),
            Some(b"DUAL") => Ok(ConsensusType::Dual),
            Some(v) => Err(format!("No enum constant ConsensusType.{}", String::from_utf8_lossy(v))),
        }
    }

    /// SvUtils.isDiscordant with the default fragment length of 1000
    #[inline]
    pub fn is_discordant(&self) -> bool {
        if !self.paired() {
            return false;
        }
        if self.tid != self.mtid {
            return true;
        }
        if self.reverse() == self.mate_reverse() {
            return true;
        }
        let fragment = self.tlen.unsigned_abs();
        fragment == 0 || fragment >= 1000
    }

    /// SamRecordUtils.getMateAlignmentEnd: from the MC tag, 0 when absent or "*"
    pub fn mate_alignment_end(&self) -> i32 {
        match self.mc {
            None => 0,
            Some(b"*") => 0,
            Some(c) => alignment_end_from_cigar(self.mate_start, c),
        }
    }

    /// The first alignment of the SA tag as SupplementaryReadData.extractAlignment reads it: (chromosome,
    /// position); None where Java gets null, an error where it throws
    pub fn first_supplementary(&self) -> Option<Result<(&'a [u8], i32), String>> {
        let sa = self.sa?;
        if sa.is_empty() {
            return None;
        }
        let first: &[u8] = if sa.contains(&b',') {
            match sa.iter().position(|&b| b == b';') {
                Some(p) => &sa[..p],
                None => sa,
            }
        } else {
            // the alternative delimitation of Java's tests: items separated by ';'
            return parse_supplementary(sa, b';');
        };
        parse_supplementary(first, b',')
    }

    /// The bases of the read's M operators (htsjdk CigarOperator.M only)
    #[cfg(test)]
    pub fn m_bases(&self) -> i64 {
        let mut n = 0i64;
        for &c in self.cigar {
            if c & 0xf == OP_M {
                n += (c >> 4) as i64;
            }
        }
        n
    }

    #[inline]
    pub fn has_op(&self, op: u32) -> bool {
        self.cigar.iter().any(|&c| c & 0xf == op)
    }
}

/// SupplementaryReadData.fromAlignment: exactly 6 fields once split as Java's String.split does
fn parse_supplementary(s: &[u8], delim: u8) -> Option<Result<(&[u8], i32), String>> {
    let mut fields: Vec<&[u8]> = s.split(|&b| b == delim).collect();
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    if fields.len() != 6 {
        return None;
    }
    let num = |f: &[u8]| std::str::from_utf8(f).ok().and_then(hmf_common::java_text::parse_int);
    let Some(pos) = num(fields[1]) else {
        return Some(Err(format!("NumberFormatException: For input string: \"{}\"", String::from_utf8_lossy(fields[1]))));
    };
    if fields[2].is_empty() {
        return Some(Err("StringIndexOutOfBoundsException".into()));
    }
    if num(fields[4]).is_none() || num(fields[5]).is_none() {
        return Some(Err("NumberFormatException in the SA tag".into()));
    }
    Some(Ok((fields[0], pos)))
}

/// CigarUtils.getReadBoundaryPosition(readStart, cigar, getReadStart=false, includeSoftClipped=false):
/// the alignment end from a CIGAR string
pub fn alignment_end_from_cigar(read_start: i32, cigar: &[u8]) -> i32 {
    let mut current = read_start as i64;
    let mut len: i64 = 0;
    for &c in cigar {
        let add = matches!(c, b'D' | b'M' | b'S' | b'N');
        if add {
            if c == b'S' && read_start as i64 == current {
                // the left clip
            } else if c == b'S' {
                break;
            } else {
                current += len;
            }
            len = 0;
            continue;
        }
        if c.is_ascii_digit() {
            len = len * 10 + (c - b'0') as i64;
        } else {
            len = 0;
        }
    }
    (current - 1) as i32
}

#[cfg(test)]
pub mod test_records {
    //! Records for tests, as SamRecordTestUtils.createSamRecord makes them
    use super::*;
    use rust_htslib::bam::record::{Aux, CigarString};

    pub struct Spec<'s> {
        pub name: &'s str,
        pub tid: i32,
        pub start: i32,
        pub bases_len: usize,
        pub cigar: &'s str,
        /// -1 for an unmapped mate ("*")
        pub mtid: i32,
        pub mate_start: i32,
        pub is_reversed: bool,
        pub supplementary: bool,
        /// (chromosome, position, cigar) of an SA tag
        pub supp: Option<(&'s str, i32, &'s str)>,
        pub mate_cigar: Option<&'s str>,
    }

    impl Default for Spec<'_> {
        fn default() -> Self {
            Spec {
                name: "READ",
                tid: 0,
                start: 1,
                bases_len: 10,
                cigar: "10M",
                mtid: 0,
                mate_start: 1,
                is_reversed: false,
                supplementary: false,
                supp: None,
                mate_cigar: None,
            }
        }
    }

    /// SamRecordTestUtils.createSamRecord: paired and first of pair, mapq 60, base quals 37, a proper
    /// pair when the mate is mapped (else mate unmapped at the read's own start), insert size the
    /// distance between the starts
    pub fn make(s: &Spec) -> bam::Record {
        let mut r = bam::Record::new();
        let cigar = CigarString::try_from(s.cigar.as_bytes()).unwrap();
        let bases = vec![b'A'; s.bases_len];
        r.set(s.name.as_bytes(), Some(&cigar), &bases, &vec![37u8; s.bases_len]);
        r.set_tid(s.tid);
        r.set_pos((s.start - 1) as i64);
        r.set_mapq(60);
        let mut flags = FLAG_PAIRED | FLAG_FIRST;
        if s.is_reversed {
            flags |= FLAG_REVERSE;
        }
        if s.mtid >= 0 {
            flags |= FLAG_PROPER_PAIR;
            r.set_mtid(s.mtid);
            r.set_mpos((s.mate_start - 1) as i64);
        } else {
            flags |= FLAG_MATE_UNMAPPED;
            r.set_mtid(-1);
            r.set_mpos((s.start - 1) as i64);
        }
        if s.supplementary {
            flags |= FLAG_SUPPLEMENTARY;
        }
        r.set_flags(flags);
        let tlen = if s.mtid == s.tid { (s.start - s.mate_start).abs() } else { 0 };
        r.set_insert_size(tlen as i64);
        if let Some(mc) = s.mate_cigar {
            r.push_aux(b"MC", Aux::String(mc)).unwrap();
        }
        if let Some((chr, pos, cig)) = s.supp {
            r.push_aux(b"SA", Aux::String(&format!("{chr},{pos},+,{cig},60,0;"))).unwrap();
        }
        r
    }

    /// SamRecordUtils.addConsensusReadAttribute
    pub fn add_consensus(r: &mut bam::Record, read_count: i32, first_in_pair: i32, consensus_type: &str) {
        r.push_aux(b"CR", Aux::String(&format!("{read_count};{first_in_pair}"))).unwrap();
        r.push_aux(b"UT", Aux::String(consensus_type)).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::test_records::*;
    use super::*;

    #[test]
    fn view_and_tags() {
        let mut r = make(&Spec {
            name: "READ_001",
            start: 10,
            bases_len: 45,
            cigar: "5S10M10D10M5I10M5S",
            mate_start: 30,
            supp: Some(("1", 40, "5S10M10N10M5S")),
            mate_cigar: Some("5S10M10D10M5I10M5S"),
            ..Default::default()
        });
        let v = Rec::of(&r);
        assert_eq!((v.start, v.end, v.mate_start), (10, 49, 30));
        assert_eq!(v.name, b"READ_001");
        assert_eq!(v.quals.len(), 45);
        assert!(!v.has_cr && v.paired() && v.first_of_pair());
        assert_eq!(v.first_supplementary().unwrap().unwrap(), (&b"1"[..], 40));
        assert_eq!(v.mate_alignment_end(), 30 + 40 - 1);
        assert_eq!(v.m_bases(), 30);
        add_consensus(&mut r, 2, 1, "DUAL");
        let v = Rec::of(&r);
        assert!(v.has_cr);
        assert_eq!(v.consensus_type().unwrap(), ConsensusType::Dual);
    }

    #[test]
    fn cigar_ends() {
        assert_eq!(alignment_end_from_cigar(100, b"10M"), 109);
        assert_eq!(alignment_end_from_cigar(100, b"5S10M5S"), 109);
        assert_eq!(alignment_end_from_cigar(100, b"5H10M2I3D4M"), 116);
        assert_eq!(alignment_end_from_cigar(100, b"10M100N10M"), 219);
    }

    #[test]
    fn supplementary_parsing() {
        assert!(parse_supplementary(b"chr1,100,+,10M", b',').is_none());
        assert!(parse_supplementary(b"chr1,x,+,10M,60,0", b',').unwrap().is_err());
        assert_eq!(parse_supplementary(b"chr2,5,-,10M,0,1", b',').unwrap().unwrap(), (&b"chr2"[..], 5));
    }
}
