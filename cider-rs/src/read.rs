//! A BAM record as CIDER uses it (htsjdk SAMRecord): the fields, with htsjdk's conventions (1-based
//! alignment start, alignment end 0 for unmapped reads, "*" for no reference), the CIGAR, bases and
//! qualities, the mate CIGAR (MC), and the record's BAM encoding for cider.bam.

use std::sync::Arc;

/// htsjdk CigarOperator, in BAM code order (M I D N S H P = X)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CigarOp {
    M,
    I,
    D,
    N,
    S,
    H,
    P,
    Eq,
    X,
}

impl CigarOp {
    pub fn from_bam(code: u32) -> CigarOp {
        match code {
            0 => CigarOp::M,
            1 => CigarOp::I,
            2 => CigarOp::D,
            3 => CigarOp::N,
            4 => CigarOp::S,
            5 => CigarOp::H,
            6 => CigarOp::P,
            7 => CigarOp::Eq,
            _ => CigarOp::X,
        }
    }
    pub fn from_char(c: u8) -> Option<CigarOp> {
        Some(match c {
            b'M' => CigarOp::M,
            b'I' => CigarOp::I,
            b'D' => CigarOp::D,
            b'N' => CigarOp::N,
            b'S' => CigarOp::S,
            b'H' => CigarOp::H,
            b'P' => CigarOp::P,
            b'=' => CigarOp::Eq,
            b'X' => CigarOp::X,
            _ => return None,
        })
    }
    pub fn to_char(self) -> char {
        match self {
            CigarOp::M => 'M',
            CigarOp::I => 'I',
            CigarOp::D => 'D',
            CigarOp::N => 'N',
            CigarOp::S => 'S',
            CigarOp::H => 'H',
            CigarOp::P => 'P',
            CigarOp::Eq => '=',
            CigarOp::X => 'X',
        }
    }
    pub fn consumes_read_bases(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::I | CigarOp::S | CigarOp::Eq | CigarOp::X)
    }
    pub fn consumes_reference_bases(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::D | CigarOp::N | CigarOp::Eq | CigarOp::X)
    }
    pub fn is_clipping(self) -> bool {
        matches!(self, CigarOp::S | CigarOp::H)
    }
    pub fn is_alignment(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::Eq | CigarOp::X)
    }
    pub fn is_indel(self) -> bool {
        matches!(self, CigarOp::I | CigarOp::D)
    }
}

pub type Cigar = Vec<(CigarOp, u32)>;

/// CigarUtils.cigarElementsFromStr (and TextCigarCodec): "*" or empty is no elements
pub fn cigar_from_str(s: &str) -> Cigar {
    let mut out = Vec::new();
    let mut n: u32 = 0;
    for &b in s.as_bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + (b - b'0') as u32;
        } else if let Some(op) = CigarOp::from_char(b) {
            out.push((op, n));
            n = 0;
        }
    }
    out
}

pub fn cigar_to_string(c: &[(CigarOp, u32)]) -> String {
    if c.is_empty() {
        return "*".into();
    }
    let mut s = String::new();
    for (op, len) in c {
        s.push_str(&len.to_string());
        s.push(op.to_char());
    }
    s
}

pub fn reference_length(c: &[(CigarOp, u32)]) -> i32 {
    c.iter().filter(|(op, _)| op.consumes_reference_bases()).map(|(_, l)| *l as i32).sum()
}

pub fn left_soft_clip(c: &[(CigarOp, u32)]) -> i32 {
    match c.first() {
        Some((CigarOp::S, l)) => *l as i32,
        _ => 0,
    }
}

pub fn right_soft_clip(c: &[(CigarOp, u32)]) -> i32 {
    match c.last() {
        Some((CigarOp::S, l)) => *l as i32,
        _ => 0,
    }
}

/// htsjdk AlignmentBlock: read start (1-based), reference start (1-based), length
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HtsAlignmentBlock {
    pub read_start: i32,
    pub reference_start: i32,
    pub length: i32,
}

/// SAMUtils.getAlignmentBlocks: one block per M, = or X element
pub fn alignment_blocks(c: &[(CigarOp, u32)], alignment_start: i32) -> Vec<HtsAlignmentBlock> {
    let mut blocks = Vec::new();
    let mut read_base = 1;
    let mut ref_base = alignment_start;
    for &(op, len) in c {
        let len = len as i32;
        match op {
            CigarOp::H | CigarOp::P => {}
            CigarOp::S | CigarOp::I => read_base += len,
            CigarOp::N | CigarOp::D => ref_base += len,
            CigarOp::M | CigarOp::Eq | CigarOp::X => {
                blocks.push(HtsAlignmentBlock { read_start: read_base, reference_start: ref_base, length: len });
                read_base += len;
                ref_base += len;
            }
        }
    }
    blocks
}

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST: u16 = 0x40;
#[allow(dead_code)]
pub const FLAG_DUPLICATE: u16 = 0x400;

/// A record. Positions follow htsjdk: `alignment_start` 1-based (0 when there is none).
#[derive(Debug)]
#[allow(dead_code)]
pub struct CiderRead {
    pub name: String,
    pub flags: u16,
    /// reference index, -1 for none
    pub tid: i32,
    pub contig: Arc<str>,
    pub alignment_start: i32,
    pub mapq: u8,
    pub mate_tid: i32,
    pub mate_contig: Arc<str>,
    pub mate_alignment_start: i32,
    pub insert_size: i32,
    pub cigar: Cigar,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    /// the MC tag
    pub mate_cigar: Option<String>,
    /// the trim CiderUtils.setAdapterDnaTrim stores in the record (left, right)
    pub adapter_trim: (i32, i32),
    /// position in the input processing order (for a stable identity)
    pub order: usize,
}

impl CiderRead {
    pub fn paired(&self) -> bool {
        self.flags & FLAG_PAIRED != 0
    }
    pub fn unmapped(&self) -> bool {
        self.flags & FLAG_UNMAPPED != 0
    }
    pub fn mate_unmapped_flag(&self) -> bool {
        self.flags & FLAG_MATE_UNMAPPED != 0
    }
    pub fn negative_strand(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }
    pub fn mate_negative_strand(&self) -> bool {
        self.flags & FLAG_MATE_REVERSE != 0
    }
    pub fn first_of_pair_flag(&self) -> bool {
        self.flags & FLAG_FIRST != 0
    }
    /// SamRecordUtils.firstInPair: unpaired reads count as first
    pub fn first_in_pair(&self) -> bool {
        !self.paired() || self.first_of_pair_flag()
    }
    /// SamRecordUtils.mateUnmapped
    pub fn mate_unmapped(&self) -> bool {
        self.paired() && self.mate_unmapped_flag()
    }
#[allow(dead_code)]
    pub fn duplicate(&self) -> bool {
        self.flags & FLAG_DUPLICATE != 0
    }
    pub fn read_length(&self) -> i32 {
        self.bases.len() as i32
    }
    /// getReadString: "*" when there are no bases
    pub fn read_string(&self) -> &[u8] {
        if self.bases.is_empty() { b"*" } else { &self.bases }
    }
    /// getAlignmentEnd: 0 for an unmapped read
    pub fn alignment_end(&self) -> i32 {
        if self.unmapped() { 0 } else { self.alignment_start + reference_length(&self.cigar) - 1 }
    }
    pub fn cigar_string(&self) -> String {
        cigar_to_string(&self.cigar)
    }
    pub fn alignment_blocks(&self) -> Vec<HtsAlignmentBlock> {
        alignment_blocks(&self.cigar, self.alignment_start)
    }
    /// SAMRecord.toString
    pub fn to_hts_string(&self) -> String {
        let mut s = self.name.clone();
        if self.paired() {
            s.push_str(if self.first_of_pair_flag() { " 1/2" } else { " 2/2" });
        }
        s.push_str(&format!(" {}b", self.read_length()));
        if self.unmapped() {
            s.push_str(" unmapped read.");
        } else {
            s.push_str(&format!(" aligned to {}:{}-{}.", self.contig, self.alignment_start, self.alignment_end()));
        }
        s
    }
}

/// A read shared between the candidate lists, the layouts and the output
pub type ReadRef = Arc<CiderRead>;

#[cfg(test)]
pub mod test_util {
    use super::*;

    /// a record for unit tests (htsjdk `SAMRecord(null)` with its defaults: unmapped fields, no flags)
    pub fn test_read() -> CiderRead {
        CiderRead {
            name: "read".into(),
            flags: 0,
            tid: -1,
            contig: Arc::from("*"),
            alignment_start: 0,
            mapq: 255,
            mate_tid: -1,
            mate_contig: Arc::from("*"),
            mate_alignment_start: 0,
            insert_size: 0,
            cigar: Vec::new(),
            bases: Vec::new(),
            quals: Vec::new(),
            mate_cigar: None,
            adapter_trim: (0, 0),
            order: 0,
        }
    }

    /// SAMUtils.fastqToPhred
    pub fn phred(s: &str) -> Vec<u8> {
        s.bytes().map(|b| b - 33).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_and_blocks() {
        let c = cigar_from_str("20S50M7000N20M1000N30M80S");
        assert_eq!(cigar_to_string(&c), "20S50M7000N20M1000N30M80S");
        assert_eq!(reference_length(&c), 50 + 7000 + 20 + 1000 + 30);
        let b = alignment_blocks(&c, 100000);
        assert_eq!(b.len(), 3);
        assert_eq!(b[0], HtsAlignmentBlock { read_start: 21, reference_start: 100000, length: 50 });
        assert_eq!(b[1], HtsAlignmentBlock { read_start: 71, reference_start: 107050, length: 20 });
        assert_eq!(b[2], HtsAlignmentBlock { read_start: 91, reference_start: 108070, length: 30 });
        assert_eq!(left_soft_clip(&c), 20);
        assert_eq!(right_soft_clip(&c), 80);
        assert!(cigar_from_str("*").is_empty());
    }
}
