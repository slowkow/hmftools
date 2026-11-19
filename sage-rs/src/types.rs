//! Core value types shared by every stage.

use std::cmp::Ordering;

pub const SE_START: usize = 0;
pub const SE_END: usize = 1;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
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

impl Op {
    pub fn consumes_read(self) -> bool {
        matches!(self, Op::M | Op::I | Op::S | Op::Eq | Op::X)
    }
    pub fn consumes_ref(self) -> bool {
        matches!(self, Op::M | Op::D | Op::N | Op::Eq | Op::X)
    }
    pub fn is_indel(self) -> bool {
        matches!(self, Op::I | Op::D)
    }
    pub fn ch(self) -> char {
        match self {
            Op::M => 'M',
            Op::I => 'I',
            Op::D => 'D',
            Op::N => 'N',
            Op::S => 'S',
            Op::H => 'H',
            Op::P => 'P',
            Op::Eq => '=',
            Op::X => 'X',
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CigarEl {
    pub len: i32,
    pub op: Op,
}

impl CigarEl {
    pub fn new(len: i32, op: Op) -> Self {
        CigarEl { len, op }
    }
}

pub fn cigar_str(cigar: &[CigarEl]) -> String {
    let mut s = String::new();
    for el in cigar {
        s.push_str(&el.len.to_string());
        s.push(el.op.ch());
    }
    s
}

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST: u16 = 0x40;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConsensusType {
    None = 0,
    Single = 1,
    Dual = 2,
}

/// An owned alignment record: the subset of SAMRecord that SAGE reads.
#[derive(Clone, Debug)]
pub struct Read {
    pub name: String,
    pub flags: u16,
    pub tid: i32,
    pub start: i32,
    pub end: i32,
    pub mapq: i32,
    pub cigar: Vec<CigarEl>,
    pub bases: Vec<u8>,
    pub quals: Vec<u8>,
    pub mate_tid: i32,
    pub mate_start: i32,
    pub insert_size: i32,
    pub nm: Option<i32>,
    pub mate_cigar: Option<String>,
    pub consensus: ConsensusType,
    pub has_consensus_tag: bool,
    /// SBX: the YX tag, the read index where the duplex region starts; -1 if absent
    pub duplex_index: i32,
    /// Ultima: the tp tag, the homopolymer length change each base quality refers to; empty if absent
    pub tp: Vec<i8>,
    /// Ultima: the t0 tag as quals (the tag string less the phred offset); empty if absent
    pub t0: Vec<i8>,
    /// Ultima: the UQ tag set by REDUX, "count=index,start-end,..." of the low-qual read indices
    pub low_qual_tag: Option<String>,
    /// the CR tag set by REDUX on a consensus read; only kept when the visualiser is on
    pub consensus_read_tag: Option<String>,
}

impl Read {
    pub fn flag(&self, f: u16) -> bool {
        self.flags & f != 0
    }
    pub fn paired(&self) -> bool {
        self.flag(FLAG_PAIRED)
    }
    pub fn negative_strand(&self) -> bool {
        self.flag(FLAG_REVERSE)
    }
    pub fn supplementary(&self) -> bool {
        self.flag(FLAG_SUPPLEMENTARY)
    }
    pub fn secondary_or_supplementary(&self) -> bool {
        self.flag(FLAG_SECONDARY) || self.flag(FLAG_SUPPLEMENTARY)
    }
    pub fn improper_pair(&self) -> bool {
        self.paired() && !self.flag(FLAG_PROPER)
    }
    pub fn left_soft_clip(&self) -> i32 {
        match self.cigar.first() {
            Some(el) if el.op == Op::S => el.len,
            _ => 0,
        }
    }
    pub fn right_soft_clip(&self) -> i32 {
        match self.cigar.last() {
            Some(el) if el.op == Op::S => el.len,
            _ => 0,
        }
    }
    pub fn len(&self) -> i32 {
        self.bases.len() as i32
    }
    pub fn has_op(&self, op: Op) -> bool {
        self.cigar.iter().any(|x| x.op == op)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BaseRegion {
    pub start: i32,
    pub end: i32,
}

impl BaseRegion {
    pub fn new(start: i32, end: i32) -> Self {
        BaseRegion { start, end }
    }
    pub fn contains(&self, pos: i32) -> bool {
        pos >= self.start && pos <= self.end
    }
    pub fn length(&self) -> i32 {
        self.end - self.start
    }
    pub fn base_length(&self) -> i32 {
        self.end - self.start + 1
    }
}

pub fn positions_overlap(s1: i32, e1: i32, s2: i32, e2: i32) -> bool {
    !(s1 > e2 || e1 < s2)
}
pub fn position_within(p: i32, s: i32, e: i32) -> bool {
    p >= s && p <= e
}
pub fn positions_within(inner_s: i32, inner_e: i32, outer_s: i32, outer_e: i32) -> bool {
    inner_s <= inner_e && inner_s >= outer_s && inner_e <= outer_e
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Tier {
    Hotspot,
    Panel,
    HighConfidence,
    LowConfidence,
}

impl Tier {
    pub fn name(self) -> &'static str {
        match self {
            Tier::Hotspot => "HOTSPOT",
            Tier::Panel => "PANEL",
            Tier::HighConfidence => "HIGH_CONFIDENCE",
            Tier::LowConfidence => "LOW_CONFIDENCE",
        }
    }
}

/// A variant on the chromosome currently being processed (chromosome is implicit).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SimpleVariant {
    pub pos: i32,
    pub ref_: Vec<u8>,
    pub alt: Vec<u8>,
}

impl SimpleVariant {
    pub fn new(pos: i32, ref_: &[u8], alt: &[u8]) -> Self {
        SimpleVariant { pos, ref_: ref_.to_vec(), alt: alt.to_vec() }
    }
    pub fn ref_len(&self) -> i32 {
        self.ref_.len() as i32
    }
    pub fn alt_len(&self) -> i32 {
        self.alt.len() as i32
    }
    pub fn is_indel(&self) -> bool {
        self.ref_.len() != self.alt.len()
    }
    pub fn is_snv(&self) -> bool {
        self.ref_.len() == 1 && self.alt.len() == 1
    }
    pub fn is_mnv(&self) -> bool {
        self.ref_.len() == self.alt.len() && self.ref_.len() > 1
    }
    pub fn is_delete(&self) -> bool {
        self.ref_.len() > self.alt.len()
    }
    pub fn is_insert(&self) -> bool {
        self.alt.len() > self.ref_.len()
    }
    pub fn indel_length(&self) -> i32 {
        self.alt_len() - self.ref_len()
    }
    pub fn indel_length_abs(&self) -> i32 {
        self.indel_length().abs()
    }
    /// SimpleVariant.end()
    pub fn end(&self) -> i32 {
        self.pos + self.ref_len() - 1
    }
    /// SimpleVariant.positionEnd()
    pub fn position_end(&self) -> i32 {
        if self.is_snv() {
            self.pos
        } else if self.is_insert() {
            self.pos + 1
        } else {
            self.pos + self.ref_len() - 1
        }
    }
    pub fn is_long_insert(&self) -> bool {
        self.is_insert() && self.indel_length() >= crate::constants::LONG_INSERT_LENGTH
    }
    pub fn ref_str(&self) -> &str {
        std::str::from_utf8(&self.ref_).unwrap()
    }
    pub fn alt_str(&self) -> &str {
        std::str::from_utf8(&self.alt).unwrap()
    }
    /// Candidates.SimpleVariantComparator
    pub fn compare_full(&self, o: &SimpleVariant) -> Ordering {
        self.pos
            .cmp(&o.pos)
            .then_with(|| self.ref_.len().max(self.alt.len()).cmp(&o.ref_.len().max(o.alt.len())))
            .then_with(|| self.ref_.cmp(&o.ref_))
            .then_with(|| self.alt.cmp(&o.alt))
    }
}

/// HumanChromosome.chromosomeRank: 1-22, X=23, Y=24, MT=25; None for other contigs.
pub fn chromosome_rank(chr: &str) -> Option<i32> {
    let t = chr.strip_prefix("chr").unwrap_or(chr);
    if t.eq_ignore_ascii_case("X") {
        Some(23)
    } else if t.eq_ignore_ascii_case("Y") {
        Some(24)
    } else if t.eq_ignore_ascii_case("MT") || t.eq_ignore_ascii_case("M") {
        Some(25)
    } else {
        t.parse::<i32>().ok().filter(|x| (1..=22).contains(x))
    }
}

pub fn is_human_chromosome(chr: &str) -> bool {
    matches!(chromosome_rank(chr), Some(r) if r <= 24)
}
pub fn is_mito(chr: &str) -> bool {
    chromosome_rank(chr) == Some(25)
}
