//! A BAM record as TEAL looks at it (the htsjdk SAMRecord accessors it uses), and ReadGroup: the records
//! of one fragment (one read name), primary reads and supplementary ones apart.

use anyhow::{bail, Result};
use rust_htslib::bam;
use std::sync::Arc;

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_FIRST: u16 = 0x40;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";

/// The 4-bit base codes of BAM
const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

#[derive(Clone, Debug)]
pub struct Read {
    pub name: Box<[u8]>,
    pub flags: u16,
    pub tid: i32,
    /// getReferenceName ("*" when unplaced)
    pub chrom: Arc<str>,
    /// getAlignmentStart (1-based; 0 when unplaced)
    pub start: i32,
    /// the reference bases the CIGAR covers
    pub ref_len: i32,
    /// raw CIGAR operations (length << 4 | op)
    pub cigar: Vec<u32>,
    pub mapq: u8,
    pub mtid: i32,
    pub mate_chrom: Arc<str>,
    /// getMateAlignmentStart (0 when none)
    pub mate_start: i32,
    pub isize: i32,
    /// the SA attribute
    pub sa: Option<Box<str>>,
    /// the bases as text (only when decoded; see `seq_len`)
    pub seq: Vec<u8>,
    pub seq_len: usize,
    /// the record itself, for writing it out
    pub raw: Option<Box<bam::Record>>,
}

/// The contig names of a header, shared by the reads
pub struct Contigs {
    pub names: Vec<Arc<str>>,
    pub star: Arc<str>,
}

impl Contigs {
    pub fn from_header(h: &bam::HeaderView) -> Contigs {
        let names = (0..h.target_count()).map(|t| Arc::from(String::from_utf8_lossy(h.tid2name(t)).as_ref())).collect();
        Contigs { names, star: Arc::from("*") }
    }
    pub fn name(&self, tid: i32) -> Arc<str> {
        if tid >= 0 && (tid as usize) < self.names.len() { self.names[tid as usize].clone() } else { self.star.clone() }
    }
    pub fn tid(&self, name: &str) -> Option<i32> {
        self.names.iter().position(|n| &**n == name).map(|i| i as i32)
    }
}

/// The 4-bit packed sequence as text
pub fn unpack_seq(packed: &[u8], len: usize) -> Vec<u8> {
    (0..len).map(|i| BASES[((packed[i / 2] >> if i % 2 == 0 { 4 } else { 0 }) & 15) as usize]).collect()
}

pub fn cigar_ref_len(cigar: &[u32]) -> i32 {
    cigar.iter().filter(|&&c| matches!(c & 15, 0 | 2 | 3 | 7 | 8)).map(|&c| (c >> 4) as i32).sum()
}

impl Read {
    /// A Read of a BAM record. `decode_seq` decodes the bases as text; `keep_raw` keeps the record.
    pub fn from_record(rec: &bam::Record, contigs: &Contigs, decode_seq: bool, keep_raw: bool) -> Read {
        let cigar: Vec<u32> = rec.raw_cigar().to_vec();
        let sa = match rec.aux(b"SA") {
            Ok(bam::record::Aux::String(s)) => Some(Box::from(s)),
            _ => None,
        };
        let seq_len = rec.seq_len();
        Read {
            name: Box::from(rec.qname()),
            flags: rec.flags(),
            tid: rec.tid(),
            chrom: contigs.name(rec.tid()),
            start: rec.pos() as i32 + 1,
            ref_len: cigar_ref_len(&cigar),
            cigar,
            mapq: rec.mapq(),
            mtid: rec.mtid(),
            mate_chrom: contigs.name(rec.mtid()),
            mate_start: rec.mpos() as i32 + 1,
            isize: rec.insert_size() as i32,
            sa,
            seq: if decode_seq { unpack_seq(rec.seq().encoded, seq_len) } else { Vec::new() },
            seq_len,
            raw: if keep_raw { Some(Box::new(rec.clone())) } else { None },
        }
    }

    pub fn paired(&self) -> bool {
        self.flags & FLAG_PAIRED != 0
    }
    pub fn unmapped(&self) -> bool {
        self.flags & FLAG_UNMAPPED != 0
    }
    pub fn neg(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }
    pub fn first_of_pair_flag(&self) -> bool {
        self.flags & FLAG_FIRST != 0
    }
    /// SamRecordUtils.firstInPair
    pub fn first_in_pair(&self) -> bool {
        !self.paired() || self.first_of_pair_flag()
    }
    pub fn secondary(&self) -> bool {
        self.flags & FLAG_SECONDARY != 0
    }
    pub fn supplementary(&self) -> bool {
        self.flags & FLAG_SUPPLEMENTARY != 0
    }
    pub fn secondary_or_supplementary(&self) -> bool {
        self.flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0
    }
    pub fn duplicate(&self) -> bool {
        self.flags & FLAG_DUPLICATE != 0
    }
    /// getAlignmentEnd: 0 for an unmapped read
    pub fn end(&self) -> i32 {
        if self.unmapped() { 0 } else { self.start + self.ref_len - 1 }
    }
    /// getReadString: "*" when there are no bases
    pub fn read_string(&self) -> &[u8] {
        if self.seq.is_empty() { b"*" } else { &self.seq }
    }
    /// getReadLength
    pub fn read_length(&self) -> usize {
        self.seq_len
    }
    /// getCigarString ("*" for no CIGAR)
    pub fn cigar_string(&self) -> String {
        if self.cigar.is_empty() {
            return "*".to_string();
        }
        let mut s = String::new();
        for &c in &self.cigar {
            s.push_str(&(c >> 4).to_string());
            s.push(CIGAR_OPS.get((c & 15) as usize).copied().unwrap_or(b'?') as char);
        }
        s
    }
    /// CigarUtils.leftSoftClipLength
    pub fn left_soft_clip(&self) -> usize {
        match self.cigar.first() {
            Some(&c) if c & 15 == 4 => (c >> 4) as usize,
            _ => 0,
        }
    }
    /// CigarUtils.rightSoftClipLength
    pub fn right_soft_clip(&self) -> usize {
        match self.cigar.last() {
            Some(&c) if c & 15 == 4 => (c >> 4) as usize,
            _ => 0,
        }
    }
    /// CigarUtils.leftSoftClipBases (None when not clipped)
    pub fn left_soft_clip_bases(&self) -> Option<&[u8]> {
        let n = self.left_soft_clip();
        if n == 0 { None } else { Some(&self.read_string()[..n]) }
    }
    /// CigarUtils.rightSoftClipBases (None when not clipped)
    pub fn right_soft_clip_bases(&self) -> Option<&[u8]> {
        let n = self.right_soft_clip();
        let s = self.read_string();
        if n == 0 { None } else { Some(&s[s.len() - n..]) }
    }
    /// Sets the flags, of the kept record too
    pub fn set_flags(&mut self, flags: u16) {
        self.flags = flags;
        if let Some(r) = self.raw.as_mut() {
            r.set_flags(flags);
        }
    }

    /// A read for tests (no record kept)
    #[cfg(test)]
    pub fn test(name: &str, chrom: &str, start: i32, flags: u16) -> Read {
        Read {
            name: Box::from(name.as_bytes()),
            flags,
            tid: 0,
            chrom: Arc::from(chrom),
            start,
            ref_len: 0,
            cigar: Vec::new(),
            mapq: 0,
            mtid: -1,
            mate_chrom: Arc::from("*"),
            mate_start: 0,
            isize: 0,
            sa: None,
            seq: Vec::new(),
            seq_len: 0,
            raw: None,
        }
    }
    #[cfg(test)]
    pub fn with_cigar(mut self, cigar: &str) -> Read {
        self.cigar = parse_cigar(cigar);
        self.ref_len = cigar_ref_len(&self.cigar);
        self
    }
    #[cfg(test)]
    pub fn with_seq(mut self, seq: &str) -> Read {
        self.seq = seq.as_bytes().to_vec();
        self.seq_len = seq.len();
        self
    }
    #[cfg(test)]
    pub fn with_mate(mut self, chrom: &str, start: i32) -> Read {
        self.mate_chrom = Arc::from(chrom);
        self.mate_start = start;
        self.mtid = 1;
        self
    }
    #[cfg(test)]
    pub fn with_sa(mut self, sa: &str) -> Read {
        self.sa = Some(Box::from(sa));
        self
    }
}

#[cfg(test)]
pub fn parse_cigar(s: &str) -> Vec<u32> {
    let mut out = Vec::new();
    let mut n = 0u32;
    for c in s.bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else {
            let op = CIGAR_OPS.iter().position(|&o| o == c).unwrap() as u32;
            out.push((n << 4) | op);
            n = 0;
        }
    }
    out
}

/// ReadGroup.SupplementaryAlignment: one entry of an SA attribute
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SuppAlignment {
    pub first_of_pair: bool,
    pub chromosome: String,
    pub negative_strand: bool,
    pub position: i32,
    pub cigar: String,
}

impl SuppAlignment {
    /// isMatch: hard and soft clips are equivalent
    pub fn is_match(&self, read: &Read) -> bool {
        read.first_in_pair() == self.first_of_pair
            && read.start == self.position
            && read.neg() == self.negative_strand
            && &*read.chrom == self.chromosome.as_str()
            && read.cigar_string().replace('H', "S") == self.cigar.replace('H', "S")
    }
}

/// Integer.parseInt
fn parse_java_int(s: &str) -> Result<i32> {
    let digits = s.strip_prefix('-').or_else(|| s.strip_prefix('+')).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        bail!("java.lang.NumberFormatException: For input string: \"{s}\"");
    }
    s.parse::<i32>().map_err(|_| anyhow::anyhow!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

/// ReadGroup.suppAlignmentPositions: `rname,pos,strand,CIGAR,mapQ,NM;` entries; those with fewer than five
/// fields are skipped
pub fn supp_alignment_positions(first_of_pair: bool, sa: &str) -> Result<Vec<SuppAlignment>> {
    let mut out = Vec::new();
    for si in sa.split(';') {
        let items: Vec<&str> = si.split(',').collect();
        if items.len() < 5 {
            continue;
        }
        out.push(SuppAlignment {
            first_of_pair,
            chromosome: items[0].to_string(),
            position: parse_java_int(items[1])?,
            negative_strand: items[2] == "-",
            cigar: items[3].to_string(),
        });
    }
    Ok(out)
}

/// ReadGroup
#[derive(Clone, Debug, Default)]
pub struct ReadGroup {
    pub name: Box<[u8]>,
    pub reads: Vec<Read>,
    pub supps: Vec<Read>,
}

impl ReadGroup {
    pub fn new(name: &[u8]) -> ReadGroup {
        ReadGroup { name: Box::from(name), reads: Vec::new(), supps: Vec::new() }
    }

    /// isPairedRead: the paired flag of the first read, else of the first supplementary read
    pub fn is_paired_read(&self) -> bool {
        self.reads.first().or(self.supps.first()).map(|r| r.paired()).unwrap_or(false)
    }
    pub fn first_of_pair(&self) -> Option<&Read> {
        self.reads.iter().find(|r| r.first_in_pair())
    }
    pub fn second_of_pair(&self) -> Option<&Read> {
        self.reads.iter().find(|r| !r.first_in_pair())
    }
    /// allReads: the reads, then the supplementary reads
    pub fn all_reads(&self) -> impl Iterator<Item = &Read> {
        self.reads.iter().chain(self.supps.iter())
    }
    pub fn all_reads_len(&self) -> usize {
        self.reads.len() + self.supps.len()
    }

    /// isComplete: None when complete, else why not (as Java logs it)
    pub fn incomplete_reason(&self) -> Result<Option<String>> {
        if self.reads.is_empty() {
            return Ok(Some("Read is empty".to_string()));
        }
        if self.reads[0].paired() && self.reads.len() != 2 {
            return Ok(Some(format!("{} missing mate pair", describe(&self.reads[0]))));
        }
        for read in &self.reads {
            if let Some(sa) = &read.sa {
                for s in supp_alignment_positions(read.first_in_pair(), sa)? {
                    if !self.supps.iter().any(|r| s.is_match(r)) {
                        return Ok(Some(format!("{} Missing supplementary read: aligned to {}:{}", describe(&self.reads[0]), s.chromosome, s.position)));
                    }
                }
            }
        }
        Ok(None)
    }

    pub fn is_complete(&self) -> Result<bool> {
        Ok(self.incomplete_reason()?.is_none())
    }

    /// contains: a read of the same kind, mate, position, strand, contig and CIGAR
    pub fn contains(&self, read: &Read) -> bool {
        if *read.name != *self.name {
            return false;
        }
        let list = if read.secondary_or_supplementary() { &self.supps } else { &self.reads };
        list.iter().any(|x| {
            x.first_in_pair() == read.first_in_pair() && x.start == read.start && x.neg() == read.neg() && x.chrom == read.chrom && x.cigar == read.cigar
        })
    }

    /// acceptRead: secondary alignments are taken as supplementary ones (their flags changed)
    pub fn accept_read(&mut self, mut read: Read) -> bool {
        if self.contains(&read) {
            return false;
        }
        if read.secondary_or_supplementary() {
            if read.secondary() {
                read.set_flags((read.flags | FLAG_SUPPLEMENTARY) & !FLAG_SECONDARY);
            }
            self.supps.push(read);
        } else {
            self.reads.push(read);
        }
        true
    }

    /// invariant
    pub fn invariant(&self) -> bool {
        if self.reads.len() > 2 {
            return false;
        }
        if !self.reads.iter().all(|x| *x.name == *self.name) || !self.supps.iter().all(|x| *x.name == *self.name) {
            return false;
        }
        if self.reads.iter().any(|x| x.supplementary()) {
            return false;
        }
        self.supps.iter().all(|x| x.supplementary())
    }

    /// findMissingReadBaseRegions: the mate's position when one read of a pair is missing, and the positions
    /// of missing supplementary reads (contig, position)
    pub fn find_missing_read_base_regions(&self) -> Result<Vec<(String, i32)>> {
        let mut out = Vec::new();
        if self.reads.len() == 1 {
            let read = &self.reads[0];
            if read.paired() && read.mtid >= 0 {
                if read.mate_chrom.is_empty() || read.mate_start <= 0 {
                    crate::warn!("read({}) invalid mate reference, mate ref index({})", describe(read), read.mtid);
                } else {
                    out.push((read.mate_chrom.to_string(), read.mate_start));
                }
            }
        }
        for read in &self.reads {
            if let Some(sa) = &read.sa {
                for s in supp_alignment_positions(read.first_in_pair(), sa)? {
                    if !self.supps.iter().any(|r| s.is_match(r)) {
                        out.push((s.chromosome.clone(), s.position));
                    }
                }
            }
        }
        Ok(out)
    }
}

/// SAMRecord.toString
pub fn describe(r: &Read) -> String {
    let name = String::from_utf8_lossy(&r.name);
    let which = if r.paired() { if r.first_of_pair_flag() { " 1/2" } else { " 2/2" } } else { "" };
    if r.unmapped() {
        format!("{name}{which} {}b unmapped read.", r.seq_len)
    } else {
        format!("{name}{which} {}b aligned to {}:{}-{}.", r.seq_len, r.chrom, r.start, r.end())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ReadGroupTest.testSuppAlignmentPositions
    #[test]
    fn supp_positions() {
        let sas = supp_alignment_positions(true, "16,59998,+,64S38M49S,55,0;").unwrap();
        assert_eq!(sas.len(), 1);
        assert_eq!(sas[0].chromosome, "16");
        assert_eq!(sas[0].position, 59998);
        let sas = supp_alignment_positions(true, "16,59998,+,64S38M49S,55,0;15,42243201,+,108S35M8S,9,0;,").unwrap();
        assert_eq!(sas.len(), 2);
        assert_eq!(sas[0].chromosome, "16");
        assert_eq!(sas[0].position, 59998);
        assert_eq!(sas[1].chromosome, "15");
        assert_eq!(sas[1].position, 42243201);
    }

    /// ReadGroupTest.testReadGroupIsComplete
    #[test]
    fn group_complete() {
        let mut rg = ReadGroup::new(b"ReadGroup");
        rg.reads.push(Read::test("ReadGroup", "1", 100, FLAG_PAIRED | FLAG_FIRST).with_mate("2", 200));
        assert!(!rg.is_complete().unwrap());
        rg.reads.push(Read::test("ReadGroup", "2", 200, FLAG_PAIRED).with_mate("1", 100));
        assert!(rg.is_complete().unwrap());
        assert!(rg.invariant());
    }

    /// ReadGroupTest.testReadGroupIsCompleteWithSupplementary
    #[test]
    fn group_complete_with_supplementary() {
        let (c1, s1, s2) = ("64S38M49S", "13S43M26S", "66M43S");
        let mut rg = ReadGroup::new(b"ReadGroup");
        rg.reads.push(Read::test("ReadGroup", "1", 100, FLAG_PAIRED | FLAG_FIRST | FLAG_REVERSE).with_mate("2", 200));
        assert!(!rg.is_complete().unwrap());
        rg.reads.push(
            Read::test("ReadGroup", "2", 200, FLAG_PAIRED).with_mate("1", 100).with_cigar(c1).with_sa(&format!("3,300,+,{s1},55,0;4,400,-,{s2},55,0;,")),
        );
        assert!(!rg.is_complete().unwrap());
        rg.supps.push(
            Read::test("ReadGroup", "3", 300, FLAG_PAIRED | FLAG_SUPPLEMENTARY).with_mate("1", 100).with_cigar(s1).with_sa(&format!("2,200,+,{c1},55,0;")),
        );
        // one supplementary read still missing
        assert!(!rg.is_complete().unwrap());
        rg.supps.push(
            Read::test("ReadGroup", "4", 400, FLAG_PAIRED | FLAG_SUPPLEMENTARY | FLAG_REVERSE)
                .with_mate("1", 100)
                .with_cigar(s2)
                .with_sa(&format!("2,200,+,{c1},55,0;")),
        );
        assert!(rg.is_complete().unwrap());
        assert!(rg.invariant());
        assert_eq!(rg.find_missing_read_base_regions().unwrap(), vec![]);
    }

    #[test]
    fn accept_and_missing_regions() {
        let mut rg = ReadGroup::new(b"x");
        let r1 = Read::test("x", "chr1", 100, FLAG_PAIRED | FLAG_FIRST).with_cigar("100M").with_mate("chr2", 500).with_sa("chr3,700,-,50S50M,60,0;");
        assert!(rg.accept_read(r1.clone()));
        assert!(!rg.accept_read(r1));
        assert_eq!(rg.find_missing_read_base_regions().unwrap(), vec![("chr2".to_string(), 500), ("chr3".to_string(), 700)]);
        // a secondary alignment is kept as a supplementary one, and matches the SA entry (hard clip as soft)
        let sec = Read::test("x", "chr3", 700, FLAG_PAIRED | FLAG_FIRST | FLAG_SECONDARY | FLAG_REVERSE).with_cigar("50H50M");
        assert!(rg.accept_read(sec));
        assert!(rg.supps[0].supplementary() && !rg.supps[0].secondary());
        assert_eq!(rg.find_missing_read_base_regions().unwrap(), vec![("chr2".to_string(), 500)]);
    }
}
