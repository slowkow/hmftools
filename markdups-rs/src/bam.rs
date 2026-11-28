//! BAM headers and records as htsjdk 4.1.3 reads and writes them.
//!
//! Records are kept as their bytes (from refID on, without the length prefix). htsjdk decodes the
//! optional fields of a record into a list sorted by tag (the tag's second character first: it reads
//! tags as `(c1 << 8) | c0` shorts and inserts in order) and writes integers in the smallest type that
//! holds the value; Picard sets a tag on every record (it clears DT and sets PG), so every record it
//! writes goes through that. `encode_record` does the same.

use anyhow::{Context, Result, bail};

use crate::bgzf::{BgzfReader, BgzfWriter};

pub const FLAG_PAIRED: u16 = 0x1;
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

// ---------------------------------------------------------------------------------------------------
// header

/// One header line's tags in their order (htsjdk keeps a LinkedHashMap per header record).
pub type Attrs = Vec<(Vec<u8>, Vec<u8>)>;

fn attr_get<'a>(attrs: &'a Attrs, key: &[u8]) -> Option<&'a [u8]> {
    attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_slice())
}

/// Sets a tag as AbstractSAMHeaderRecord.setAttribute: in place if present, else at the end; `None`
/// removes it.
fn attr_set(attrs: &mut Attrs, key: &[u8], value: Option<&[u8]>) {
    match value {
        Some(v) => {
            if let Some(e) = attrs.iter_mut().find(|(k, _)| k == key) {
                e.1 = v.to_vec();
            } else {
                attrs.push((key.to_vec(), v.to_vec()));
            }
        }
        None => attrs.retain(|(k, _)| k != key),
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SqRecord {
    pub name: Vec<u8>,
    pub len: i32,
    pub attrs: Attrs,
}

/// An @RG or @PG line: its ID and the other tags.
#[derive(Clone, Debug, PartialEq)]
pub struct IdRecord {
    pub id: Vec<u8>,
    pub attrs: Attrs,
}

impl IdRecord {
    pub fn get(&self, key: &[u8]) -> Option<&[u8]> {
        attr_get(&self.attrs, key)
    }
    pub fn set(&mut self, key: &[u8], value: Option<&[u8]>) {
        attr_set(&mut self.attrs, key, value)
    }
}

/// SAMFileHeader: @HD tags (VN first), sequences, read groups, programs, comments.
#[derive(Clone, Debug, PartialEq)]
pub struct SamHeader {
    pub hd: Attrs,
    pub sq: Vec<SqRecord>,
    pub rg: Vec<IdRecord>,
    pub pg: Vec<IdRecord>,
    pub co: Vec<Vec<u8>>,
    /// whether the text had an @HD line
    pub has_hd: bool,
}

/// SAMFileHeader.CURRENT_VERSION
pub const CURRENT_VERSION: &[u8] = b"1.6";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortOrder {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
    Unknown,
}

impl SortOrder {
    pub fn parse(s: &[u8]) -> Option<SortOrder> {
        Some(match s {
            b"unsorted" => SortOrder::Unsorted,
            b"queryname" => SortOrder::Queryname,
            b"coordinate" => SortOrder::Coordinate,
            b"duplicate" => SortOrder::Duplicate,
            b"unknown" => SortOrder::Unknown,
            _ => return None,
        })
    }
    pub fn name(self) -> &'static str {
        match self {
            SortOrder::Unsorted => "unsorted",
            SortOrder::Queryname => "queryname",
            SortOrder::Coordinate => "coordinate",
            SortOrder::Duplicate => "duplicate",
            SortOrder::Unknown => "unknown",
        }
    }
}

fn header_error(msg: &str, line: &[u8]) -> anyhow::Error {
    anyhow::anyhow!("Error parsing SAM header. {msg}. Line:\n{}", String::from_utf8_lossy(line))
}

impl Default for SamHeader {
    fn default() -> Self {
        SamHeader { hd: vec![(b"VN".to_vec(), CURRENT_VERSION.to_vec())], sq: Vec::new(), rg: Vec::new(), pg: Vec::new(), co: Vec::new(), has_hd: false }
    }
}

impl SamHeader {
    /// SAMTextHeaderCodec.decode with STRICT validation: lines are read while they start with '@'.
    pub fn parse(text: &[u8]) -> Result<SamHeader> {
        let mut h = SamHeader::default();
        for raw in text.split(|&b| b == b'\n') {
            let line = raw.strip_suffix(b"\r").unwrap_or(raw);
            if line.first() != Some(&b'@') {
                break;
            }
            let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
            let kind = &fields[0][1..];
            if kind == b"CO" {
                h.co.push(line.to_vec());
                continue;
            }
            if !matches!(kind, b"HD" | b"SQ" | b"RG" | b"PG") {
                return Err(header_error("Unrecognized header record type", line));
            }
            let mut kv: Attrs = Vec::new();
            for f in &fields[1..] {
                let Some(colon) = f.iter().position(|&b| b == b':') else {
                    return Err(header_error(&format!("Problem parsing @{} key:value pair", String::from_utf8_lossy(kind)), line));
                };
                let (k, v) = (&f[..colon], &f[colon + 1..]);
                if let Some(old) = attr_get(&kv, k) {
                    if old != v {
                        return Err(header_error("Problem parsing key:value pair: tag multiply defined", line));
                    }
                    continue;
                }
                if k.len() != 2 {
                    return Err(header_error("Problem parsing key:value pair: key is not two characters", line));
                }
                kv.push((k.to_vec(), v.to_vec()));
            }
            let take = |kv: &mut Attrs, key: &[u8]| -> Option<Vec<u8>> {
                let i = kv.iter().position(|(k, _)| k == key)?;
                Some(kv.remove(i).1)
            };
            match kind {
                b"HD" => {
                    h.has_hd = true;
                    if attr_get(&kv, b"VN").is_none() {
                        return Err(header_error("@HD line missing VN tag", line));
                    }
                    if let Some(so) = attr_get(&kv, b"SO") {
                        if SortOrder::parse(so).is_none() {
                            return Err(header_error("@HD line has non-conforming SO tag value", line));
                        }
                    }
                    for (k, v) in kv {
                        attr_set(&mut h.hd, &k, Some(&v));
                    }
                }
                b"SQ" => {
                    let Some(name) = take(&mut kv, b"SN") else { return Err(header_error("@SQ line missing SN tag", line)) };
                    let Some(len) = take(&mut kv, b"LN") else { return Err(header_error("@SQ line missing LN tag", line)) };
                    let len: i32 = std::str::from_utf8(&len).ok().and_then(|s| s.parse().ok()).ok_or_else(|| header_error("bad LN", line))?;
                    h.sq.push(SqRecord { name, len, attrs: kv });
                }
                b"RG" => {
                    let Some(id) = take(&mut kv, b"ID") else { return Err(header_error("@RG line missing ID tag", line)) };
                    if attr_get(&kv, b"SM").is_none() {
                        return Err(header_error("@RG line missing SM tag", line));
                    }
                    if attr_get(&kv, b"DT").is_some() {
                        // htsjdk rewrites the date in its own form; not reproduced
                        bail!("@RG lines with a DT tag are not supported");
                    }
                    if let Some(pi) = attr_get(&kv, b"PI") {
                        if std::str::from_utf8(pi).ok().and_then(|s| s.parse::<i32>().ok()).is_none() {
                            return Err(header_error("PI is not numeric", line));
                        }
                    }
                    h.rg.push(IdRecord { id, attrs: kv });
                }
                _ => {
                    let Some(id) = take(&mut kv, b"ID") else { return Err(header_error("@PG line missing ID tag", line)) };
                    if h.pg.iter().any(|p| p.id == id) {
                        bail!("Program record with group id {} already exists in SAMFileHeader!", String::from_utf8_lossy(&id));
                    }
                    h.pg.push(IdRecord { id, attrs: kv });
                }
            }
        }
        Ok(h)
    }

    pub fn sort_order(&self) -> SortOrder {
        attr_get(&self.hd, b"SO").and_then(SortOrder::parse).unwrap_or(SortOrder::Unsorted)
    }

    pub fn set_hd(&mut self, key: &[u8], value: Option<&[u8]>) {
        attr_set(&mut self.hd, key, value)
    }

    pub fn read_group(&self, id: &[u8]) -> Option<&IdRecord> {
        self.rg.iter().find(|r| r.id == id)
    }

    /// SAMTextHeaderCodec.encode: @HD with the current version first, then @SQ, @RG, @PG and @CO.
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(b"@HD\tVN:");
        out.extend_from_slice(CURRENT_VERSION);
        for (k, v) in &self.hd {
            if k != b"VN" {
                push_tag(&mut out, k, v);
            }
        }
        out.push(b'\n');
        for s in &self.sq {
            out.extend_from_slice(b"@SQ\tSN:");
            out.extend_from_slice(&s.name);
            out.extend_from_slice(b"\tLN:");
            out.extend_from_slice(s.len.to_string().as_bytes());
            for (k, v) in &s.attrs {
                push_tag(&mut out, k, v);
            }
            out.push(b'\n');
        }
        for (kind, recs) in [(b"@RG", &self.rg), (b"@PG", &self.pg)] {
            for r in recs {
                out.extend_from_slice(kind);
                out.extend_from_slice(b"\tID:");
                out.extend_from_slice(&r.id);
                for (k, v) in &r.attrs {
                    push_tag(&mut out, k, v);
                }
                out.push(b'\n');
            }
        }
        for c in &self.co {
            out.extend_from_slice(c);
            out.push(b'\n');
        }
        out
    }

    /// SAMFileHeader.addComment
    pub fn add_comment(&mut self, comment: &[u8]) {
        if comment.starts_with(b"@CO\t") {
            self.co.push(comment.to_vec());
        } else {
            let mut c = b"@CO\t".to_vec();
            c.extend_from_slice(comment);
            self.co.push(c);
        }
    }
}

fn push_tag(out: &mut Vec<u8>, k: &[u8], v: &[u8]) {
    out.push(b'\t');
    out.extend_from_slice(k);
    out.push(b':');
    out.extend_from_slice(v);
}

/// A BAM file's header: the parsed text and the binary reference list.
pub struct BamHeader {
    pub sam: SamHeader,
    pub refs: Vec<(Vec<u8>, u32)>,
}

/// Reads the BAM header (magic, text, references) as BAMFileReader.readHeader does.
pub fn read_header(r: &mut BgzfReader) -> Result<BamHeader> {
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic).with_context(|| format!("{}: not a BAM file", r.path()))?;
    if &magic != b"BAM\x01" {
        bail!("{}: Invalid BAM file header (only BAM input is supported)", r.path());
    }
    let mut b4 = [0u8; 4];
    r.read_exact(&mut b4)?;
    let l_text = i32::from_le_bytes(b4);
    if l_text < 0 {
        bail!("{}: bad BAM header length", r.path());
    }
    let mut text = vec![0u8; l_text as usize];
    r.read_exact(&mut text)?;
    let mut sam = SamHeader::parse(&text).with_context(|| r.path().to_string())?;
    r.read_exact(&mut b4)?;
    let n_ref = i32::from_le_bytes(b4);
    if n_ref < 0 {
        bail!("{}: bad reference count", r.path());
    }
    let mut refs = Vec::with_capacity(n_ref as usize);
    for _ in 0..n_ref {
        r.read_exact(&mut b4)?;
        let l_name = i32::from_le_bytes(b4);
        if l_name <= 0 {
            bail!("{}: bad reference name length", r.path());
        }
        let mut name = vec![0u8; l_name as usize];
        r.read_exact(&mut name)?;
        name.pop();
        r.read_exact(&mut b4)?;
        refs.push((name, u32::from_le_bytes(b4)));
    }
    if !sam.sq.is_empty() {
        if sam.sq.len() != refs.len() {
            bail!(
                "{}: Number of sequences in text header ({}) != number of sequences in binary header ({})",
                r.path(),
                sam.sq.len(),
                refs.len()
            );
        }
        for (i, (s, (n, l))) in sam.sq.iter().zip(refs.iter()).enumerate() {
            if &s.name != n {
                bail!("{}: For sequence {i}, text and binary have different names", r.path());
            }
            if s.len as u32 != *l {
                bail!("{}: For sequence {i}, text and binary have different lengths", r.path());
            }
        }
    } else {
        sam.sq = refs.iter().map(|(n, l)| SqRecord { name: n.clone(), len: *l as i32, attrs: Vec::new() }).collect();
    }
    Ok(BamHeader { sam, refs })
}

/// Writes the BAM header as BAMFileWriter.writeHeader: magic, the text, the references from the
/// sequence dictionary.
pub fn write_header(w: &mut BgzfWriter, sam: &SamHeader) -> Result<()> {
    let text = sam.encode();
    let mut out = Vec::with_capacity(text.len() + 64);
    out.extend_from_slice(b"BAM\x01");
    out.extend_from_slice(&(text.len() as i32).to_le_bytes());
    out.extend_from_slice(&text);
    out.extend_from_slice(&(sam.sq.len() as i32).to_le_bytes());
    for s in &sam.sq {
        out.extend_from_slice(&((s.name.len() + 1) as i32).to_le_bytes());
        out.extend_from_slice(&s.name);
        out.push(0);
        out.extend_from_slice(&s.len.to_le_bytes());
    }
    w.write(&out)
}

// ---------------------------------------------------------------------------------------------------
// records

#[inline]
fn i32_at(b: &[u8], at: usize) -> i32 {
    i32::from_le_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])
}
#[inline]
fn u16_at(b: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([b[at], b[at + 1]])
}

/// A record's bytes, from refID on.
#[derive(Clone, Copy)]
pub struct Rec<'a>(pub &'a [u8]);

/// CIGAR operations (BAM codes)
pub const CIGAR_M: u32 = 0;
pub const CIGAR_D: u32 = 2;
pub const CIGAR_N: u32 = 3;
pub const CIGAR_S: u32 = 4;
pub const CIGAR_H: u32 = 5;
pub const CIGAR_EQ: u32 = 7;
pub const CIGAR_X: u32 = 8;

impl<'a> Rec<'a> {
    /// Checks the lengths of the variable fields.
    pub fn check(b: &'a [u8]) -> Result<Rec<'a>> {
        if b.len() < 32 {
            bail!("truncated BAM record");
        }
        let r = Rec(b);
        let need = 32 + r.l_read_name() + r.n_cigar() * 4 + r.l_seq().div_ceil(2) + r.l_seq();
        if r.l_read_name() == 0 || b.len() < need || (r.l_seq() as i64) < 0 {
            bail!("truncated or malformed BAM record");
        }
        Ok(r)
    }
    #[inline]
    pub fn ref_id(&self) -> i32 {
        i32_at(self.0, 0)
    }
    /// 0-based position
    #[inline]
    pub fn pos(&self) -> i32 {
        i32_at(self.0, 4)
    }
    #[inline]
    pub fn l_read_name(&self) -> usize {
        self.0[8] as usize
    }
    #[inline]
    pub fn mapq(&self) -> u8 {
        self.0[9]
    }
    #[inline]
    pub fn n_cigar(&self) -> usize {
        u16_at(self.0, 12) as usize
    }
    #[inline]
    pub fn flag(&self) -> u16 {
        u16_at(self.0, 14)
    }
    #[inline]
    pub fn l_seq(&self) -> usize {
        i32_at(self.0, 16) as usize
    }
    #[inline]
    pub fn mate_ref_id(&self) -> i32 {
        i32_at(self.0, 20)
    }
    #[inline]
    pub fn mate_pos(&self) -> i32 {
        i32_at(self.0, 24)
    }
    /// The read name without its terminating NUL (htsjdk keeps l_read_name - 1 bytes).
    #[inline]
    pub fn name(&self) -> &'a [u8] {
        &self.0[32..32 + self.l_read_name() - 1]
    }
    #[inline]
    fn cigar_at(&self) -> usize {
        32 + self.l_read_name()
    }
    #[inline]
    pub fn cigar_raw(&self) -> &'a [u8] {
        let at = self.cigar_at();
        &self.0[at..at + self.n_cigar() * 4]
    }
    #[inline]
    fn seq_at(&self) -> usize {
        self.cigar_at() + self.n_cigar() * 4
    }
    #[inline]
    fn qual_at(&self) -> usize {
        self.seq_at() + self.l_seq().div_ceil(2)
    }
    #[inline]
    pub fn qual(&self) -> &'a [u8] {
        let at = self.qual_at();
        &self.0[at..at + self.l_seq()]
    }
    #[inline]
    fn aux_at(&self) -> usize {
        self.qual_at() + self.l_seq()
    }
    #[inline]
    pub fn aux(&self) -> &'a [u8] {
        &self.0[self.aux_at()..]
    }
    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.flag() & FLAG_UNMAPPED != 0
    }
    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.flag() & FLAG_REVERSE != 0
    }
    #[inline]
    pub fn is_paired(&self) -> bool {
        self.flag() & FLAG_PAIRED != 0
    }
    /// SAMRecord.getMateUnmappedFlag (only meaningful for paired reads; htsjdk throws otherwise, the
    /// callers check paired first)
    #[inline]
    pub fn mate_unmapped(&self) -> bool {
        self.flag() & FLAG_MATE_UNMAPPED != 0
    }
    #[inline]
    pub fn is_secondary_or_supplementary(&self) -> bool {
        self.flag() & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0
    }
    #[inline]
    pub fn is_first_of_pair(&self) -> bool {
        self.flag() & FLAG_FIRST != 0
    }

    /// The CIGAR htsjdk uses: the CG tag's when the stored one is the long-CIGAR placeholder
    /// (`<l_seq>S<ref length>N`).
    pub fn cigar_ops(&self) -> Vec<u32> {
        let raw = self.cigar_raw();
        let ops: Vec<u32> = raw.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
        if ops.len() == 2 && ops[0] & 0xf == CIGAR_S && (ops[0] >> 4) as usize == self.l_seq() && ops[1] & 0xf == CIGAR_N {
            if let Some(TagValue::Array(b'I', data)) = find_tag(self.aux(), *b"CG").map(|(_, v)| v) {
                return data.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
            }
        }
        ops
    }

    /// Cigar.getReferenceLength: M, D, N, =, X
    pub fn reference_length(ops: &[u32]) -> i32 {
        ops.iter().filter(|&&o| matches!(o & 0xf, CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X)).map(|&o| (o >> 4) as i32).sum()
    }

    /// SAMRecord.getAlignmentEnd: 0 for an unmapped read
    pub fn alignment_end(&self, ops: &[u32]) -> i32 {
        if self.is_unmapped() { 0 } else { self.pos() + 1 + Self::reference_length(ops) - 1 }
    }

    /// SAMRecord.getUnclippedStart
    pub fn unclipped_start(&self, ops: &[u32]) -> i32 {
        let mut s = self.pos() + 1;
        for &o in ops {
            if matches!(o & 0xf, CIGAR_S | CIGAR_H) {
                s -= (o >> 4) as i32;
            } else {
                break;
            }
        }
        s
    }

    /// SAMRecord.getUnclippedEnd
    pub fn unclipped_end(&self, ops: &[u32]) -> i32 {
        let mut e = self.alignment_end(ops);
        for &o in ops.iter().rev() {
            if matches!(o & 0xf, CIGAR_S | CIGAR_H) {
                e += (o >> 4) as i32;
            } else {
                break;
            }
        }
        e
    }

    /// The RG and PG values in one scan of the optional fields (the last of a repeated tag); an error if
    /// either is not a string (htsjdk's String casts).
    pub fn rg_pg(&self) -> Result<(Option<&'a [u8]>, Option<&'a [u8]>)> {
        let aux = self.aux();
        let (mut rg, mut pg) = (None, None);
        let mut at = 0;
        while at < aux.len() {
            let (t, v, next) = parse_tag(aux, at)?;
            if t == *b"RG" || t == *b"PG" {
                let TagValue::Str(s) = v else { bail!("Value for tag {} is not a String", String::from_utf8_lossy(&t)) };
                if t == *b"RG" {
                    rg = Some(s);
                } else {
                    pg = Some(s);
                }
            }
            at = next;
        }
        Ok((rg, pg))
    }

    /// The string value of tag `tag` (SAMRecord.getStringAttribute); an error if it is not a string.
    pub fn string_tag(&self, tag: [u8; 2]) -> Result<Option<&'a [u8]>> {
        match find_tag(self.aux(), tag) {
            None => Ok(None),
            Some((_, TagValue::Str(s))) => Ok(Some(s)),
            Some(_) => bail!("Value for tag {} is not a String", String::from_utf8_lossy(&tag)),
        }
    }
}

/// A decoded optional field value, borrowing the record's bytes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TagValue<'a> {
    Char(u8),
    /// an integer of any of the types c C s S i I
    Int(i64),
    /// the 4 bytes of a float
    Float([u8; 4]),
    Str(&'a [u8]),
    /// an 'H' hex string, as stored
    Hex(&'a [u8]),
    /// B array: element type (c C s S i I f) and the raw little-endian elements
    Array(u8, &'a [u8]),
}

/// Parses the next optional field at `at`: (tag, value, next offset).
pub fn parse_tag(aux: &[u8], at: usize) -> Result<([u8; 2], TagValue<'_>, usize)> {
    if at + 3 > aux.len() {
        bail!("truncated optional field");
    }
    let tag = [aux[at], aux[at + 1]];
    let ty = aux[at + 2];
    let p = at + 3;
    let need = |n: usize| -> Result<()> { if p + n > aux.len() { bail!("truncated optional field") } else { Ok(()) } };
    Ok(match ty {
        b'A' => {
            need(1)?;
            (tag, TagValue::Char(aux[p]), p + 1)
        }
        b'c' => {
            need(1)?;
            (tag, TagValue::Int(aux[p] as i8 as i64), p + 1)
        }
        b'C' => {
            need(1)?;
            (tag, TagValue::Int(aux[p] as i64), p + 1)
        }
        b's' => {
            need(2)?;
            (tag, TagValue::Int(i16::from_le_bytes([aux[p], aux[p + 1]]) as i64), p + 2)
        }
        b'S' => {
            need(2)?;
            (tag, TagValue::Int(u16::from_le_bytes([aux[p], aux[p + 1]]) as i64), p + 2)
        }
        b'i' => {
            need(4)?;
            (tag, TagValue::Int(i32_at(aux, p) as i64), p + 4)
        }
        b'I' => {
            need(4)?;
            (tag, TagValue::Int(u32::from_le_bytes([aux[p], aux[p + 1], aux[p + 2], aux[p + 3]]) as i64), p + 4)
        }
        b'f' => {
            need(4)?;
            (tag, TagValue::Float([aux[p], aux[p + 1], aux[p + 2], aux[p + 3]]), p + 4)
        }
        b'Z' | b'H' => {
            let Some(end) = aux[p..].iter().position(|&b| b == 0) else { bail!("unterminated string field") };
            let s = &aux[p..p + end];
            (tag, if ty == b'Z' { TagValue::Str(s) } else { TagValue::Hex(s) }, p + end + 1)
        }
        b'B' => {
            need(5)?;
            let sub = aux[p];
            let n = u32::from_le_bytes([aux[p + 1], aux[p + 2], aux[p + 3], aux[p + 4]]) as usize;
            let size = match sub {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => bail!("Unrecognized tag array type: {}", sub as char),
            };
            let q = p + 5;
            if q + n * size > aux.len() {
                bail!("truncated optional field");
            }
            (tag, TagValue::Array(sub, &aux[q..q + n * size]), q + n * size)
        }
        _ => bail!("Unrecognized tag type: {}", ty as char),
    })
}

/// The value of tag `tag` as htsjdk sees it (the last one, if a record repeats a tag).
pub fn find_tag(aux: &[u8], tag: [u8; 2]) -> Option<([u8; 2], TagValue<'_>)> {
    let mut at = 0;
    let mut found = None;
    while at < aux.len() {
        let Ok((t, v, next)) = parse_tag(aux, at) else { return found };
        if t == tag {
            found = Some((t, v));
        }
        at = next;
    }
    found
}

/// htsjdk's short for a tag: second character in the high byte (SAMTag.makeBinaryTag), compared signed.
#[inline]
pub fn tag_key(tag: [u8; 2]) -> i16 {
    ((tag[1] as u16) << 8 | tag[0] as u16) as i16
}

/// BinaryTagCodec.getIntegerType
fn integer_type(v: i64) -> Result<u8> {
    Ok(if v > u32::MAX as i64 {
        bail!("Integer attribute value too large to be encoded in BAM")
    } else if v > i32::MAX as i64 {
        b'I'
    } else if v > u16::MAX as i64 {
        b'i'
    } else if v > i16::MAX as i64 {
        b'S'
    } else if v > u8::MAX as i64 {
        b's'
    } else if v > i8::MAX as i64 {
        b'C'
    } else if v >= i8::MIN as i64 {
        b'c'
    } else if v >= i16::MIN as i64 {
        b's'
    } else if v >= i32::MIN as i64 {
        b'i'
    } else {
        bail!("Integer attribute value too negative to be encoded in BAM")
    })
}

/// Writes one tag as BinaryTagCodec.writeTag would for the value htsjdk decoded.
fn write_tag(out: &mut Vec<u8>, tag: [u8; 2], v: &TagValue) -> Result<()> {
    out.extend_from_slice(&tag);
    match *v {
        TagValue::Char(c) => {
            out.push(b'A');
            out.push(c);
        }
        TagValue::Int(i) => {
            let ty = integer_type(i)?;
            out.push(ty);
            match ty {
                b'c' | b'C' => out.push(i as u8),
                b's' | b'S' => out.extend_from_slice(&(i as u16).to_le_bytes()),
                _ => out.extend_from_slice(&(i as u32).to_le_bytes()),
            }
        }
        TagValue::Float(f) => {
            out.push(b'f');
            out.extend_from_slice(&f);
        }
        TagValue::Str(s) => {
            out.push(b'Z');
            out.extend_from_slice(s);
            out.push(0);
        }
        TagValue::Hex(h) => {
            // decoded to a byte[] (StringUtil.hexStringToBytes) and written as a signed byte array
            if h.len() % 2 != 0 {
                bail!("Hex string has odd length");
            }
            let bytes: Vec<u8> = h
                .chunks_exact(2)
                .map(|c| {
                    let d = |x: u8| match x {
                        b'0'..=b'9' => Ok(x - b'0'),
                        b'a'..=b'f' => Ok(x - b'a' + 10),
                        b'A'..=b'F' => Ok(x - b'A' + 10),
                        _ => Err(anyhow::anyhow!("bad hex digit")),
                    };
                    Ok(d(c[0])? << 4 | d(c[1])?)
                })
                .collect::<Result<_>>()?;
            out.push(b'B');
            out.push(b'c');
            out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
            out.extend_from_slice(&bytes);
        }
        TagValue::Array(sub, data) => {
            let size = match sub {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                _ => 4,
            };
            out.push(b'B');
            out.push(sub);
            out.extend_from_slice(&((data.len() / size) as u32).to_le_bytes());
            out.extend_from_slice(data);
        }
    }
    Ok(())
}

/// How Picard changes a record's optional fields before writing it.
#[derive(Default, Clone, Copy)]
pub struct TagEdits<'a> {
    /// remove DT (CLEAR_DT)
    pub clear_dt: bool,
    /// set DT to this value
    pub set_dt: Option<&'a [u8]>,
    /// set PG to this value (`Some(None)`: remove it, as setAttribute(PG, null) does)
    pub set_pg: Option<Option<&'a [u8]>>,
    /// a tag was set before (by the merging of several inputs)
    pub force: bool,
}

impl TagEdits<'_> {
    /// Whether htsjdk would write the record's variable-length part afresh (any setAttribute call)
    pub fn touches(&self) -> bool {
        self.clear_dt || self.set_dt.is_some() || self.set_pg.is_some() || self.force
    }
}

/// GenomicIndexUtil.regionToBin (0-based, half-open)
pub fn region_to_bin(beg: i32, end: i32) -> i32 {
    let end = end - 1;
    if beg >> 14 == end >> 14 {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if beg >> 17 == end >> 17 {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if beg >> 20 == end >> 20 {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if beg >> 23 == end >> 23 {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if beg >> 26 == end >> 26 {
        return ((1 << 3) - 1) / 7 + (beg >> 26);
    }
    0
}

/// SAMRecord.computeIndexingBin
pub fn indexing_bin(rec: &Rec, ops: &[u32]) -> i32 {
    let start = rec.pos();
    let mut end = rec.alignment_end(ops);
    if end <= 0 {
        end = start + 1;
    }
    region_to_bin(start, end)
}

/// The bin htsjdk's BAMRecordCodec writes: 0 for a record without a position or with an end beyond
/// what the bin field can index (GenomicIndexUtil.BIN_GENOMIC_SPAN).
pub fn written_bin(rec: &Rec, ops: &[u32]) -> u16 {
    if rec.pos() + 1 == 0 {
        return 0;
    }
    if rec.alignment_end(ops) as i64 > (1i64 << 29) {
        return 0;
    }
    indexing_bin(rec, ops) as u16
}

/// Appends the record as htsjdk's BAMRecordCodec.encode writes it (length prefix included), with flag
/// `flag` and the tag edits. `ops` are the record's CIGAR operations (for the bin).
pub fn encode_record(out: &mut Vec<u8>, rec: &Rec, flag: u16, ops: &[u32], edits: &TagEdits, tags: &mut Vec<(i16, [u8; 2], usize, usize)>) -> Result<()> {
    let b = rec.0;
    let start = out.len();
    out.extend_from_slice(&[0, 0, 0, 0]);
    out.extend_from_slice(&b[0..10]);
    out.extend_from_slice(&written_bin(rec, ops).to_le_bytes());
    out.extend_from_slice(&b[12..14]);
    out.extend_from_slice(&flag.to_le_bytes());
    out.extend_from_slice(&b[16..32]);
    if !edits.touches() {
        // the variable-length part is written as read
        out.extend_from_slice(&b[32..]);
    } else {
        let seq_at = rec.seq_at();
        let qual_at = rec.qual_at();
        let l_seq = rec.l_seq();
        out.extend_from_slice(&b[32..seq_at]);
        // bases: decoded and encoded again, which zeroes the unused half of the last byte
        out.extend_from_slice(&b[seq_at..qual_at]);
        if l_seq % 2 == 1 {
            let last = out.len() - 1;
            out[last] &= 0xf0;
        }
        let q = &b[qual_at..qual_at + l_seq];
        if l_seq > 0 && q[0] == 0xff {
            out.extend(std::iter::repeat_n(0xffu8, l_seq));
        } else {
            out.extend_from_slice(q);
        }
        // optional fields: sorted by htsjdk's tag short, the last of repeated tags kept
        let aux = rec.aux();
        tags.clear();
        let mut at = 0;
        while at < aux.len() {
            let (t, _, next) = parse_tag(aux, at)?;
            tags.push((tag_key(t), t, at, next));
            at = next;
        }
        // stable: file order among equal tags; of a run of equal tags the last one stays
        tags.sort_by_key(|e| e.0);
        let mut w = 0;
        for r in 0..tags.len() {
            if r + 1 < tags.len() && tags[r + 1].0 == tags[r].0 {
                continue;
            }
            tags[w] = tags[r];
            w += 1;
        }
        tags.truncate(w);
        const NEW_PG: usize = usize::MAX;
        const NEW_DT: usize = usize::MAX - 1;
        let pg_key = tag_key(*b"PG");
        let dt_key = tag_key(*b"DT");
        let insert = |tags: &mut Vec<(i16, [u8; 2], usize, usize)>, e: (i16, [u8; 2], usize, usize)| {
            let i = tags.partition_point(|x| x.0 < e.0);
            tags.insert(i, e);
        };
        if edits.clear_dt || edits.set_dt.is_some() {
            tags.retain(|e| e.0 != dt_key);
        }
        if edits.set_dt.is_some() {
            insert(tags, (dt_key, *b"DT", NEW_DT, 0));
        }
        if let Some(pg) = edits.set_pg {
            tags.retain(|e| e.0 != pg_key);
            if pg.is_some() {
                insert(tags, (pg_key, *b"PG", NEW_PG, 0));
            }
        }
        for &(_, t, from, _) in tags.iter() {
            match from {
                NEW_PG => write_tag(out, t, &TagValue::Str(edits.set_pg.unwrap().unwrap()))?,
                NEW_DT => write_tag(out, t, &TagValue::Str(edits.set_dt.unwrap()))?,
                _ => {
                    let (_, v, _) = parse_tag(aux, from)?;
                    write_tag(out, t, &v)?;
                }
            }
        }
    }
    let len = (out.len() - start - 4) as u32;
    out[start..start + 4].copy_from_slice(&len.to_le_bytes());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a record with the given flag, position, cigar (ops), name, quality and raw aux bytes
    pub fn make_record(tid: i32, pos: i32, flag: u16, name: &[u8], cigar: &[u32], seq_len: usize, qual: u8, aux: &[u8]) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(&tid.to_le_bytes());
        b.extend_from_slice(&pos.to_le_bytes());
        b.push((name.len() + 1) as u8);
        b.push(60);
        b.extend_from_slice(&0u16.to_le_bytes());
        b.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
        b.extend_from_slice(&flag.to_le_bytes());
        b.extend_from_slice(&(seq_len as i32).to_le_bytes());
        b.extend_from_slice(&(-1i32).to_le_bytes());
        b.extend_from_slice(&(-1i32).to_le_bytes());
        b.extend_from_slice(&0i32.to_le_bytes());
        b.extend_from_slice(name);
        b.push(0);
        for &c in cigar {
            b.extend_from_slice(&c.to_le_bytes());
        }
        b.extend(std::iter::repeat_n(0x12u8, seq_len.div_ceil(2)));
        b.extend(std::iter::repeat_n(qual, seq_len));
        b.extend_from_slice(aux);
        b
    }

    #[test]
    fn clipping_and_ends() {
        // 5S20M100N30M7H at 0-based 999
        let ops = [5 << 4 | CIGAR_S, 20 << 4 | CIGAR_M, 100 << 4 | CIGAR_N, 30 << 4 | CIGAR_M, 7 << 4 | CIGAR_H];
        let b = make_record(0, 999, 0, b"r", &ops, 55, 30, b"");
        let r = Rec::check(&b).unwrap();
        assert_eq!(r.unclipped_start(&ops), 1000 - 5);
        assert_eq!(r.alignment_end(&ops), 1000 + 150 - 1);
        assert_eq!(r.unclipped_end(&ops), 1149 + 7);
    }

    #[test]
    fn tags_are_sorted_narrowed_and_edited() {
        // NH:C:1 HI:i:1 AS:i:300 nM:i:0 RG:Z:x (STAR's order) and a DT to clear
        let mut aux = Vec::new();
        aux.extend_from_slice(b"NHC\x01");
        aux.extend_from_slice(b"HIi\x01\x00\x00\x00");
        aux.extend_from_slice(b"ASi\x2c\x01\x00\x00");
        aux.extend_from_slice(b"nMi\x00\x00\x00\x00");
        aux.extend_from_slice(b"RGZx\x00");
        aux.extend_from_slice(b"DTZLB\x00");
        let b = make_record(0, 10, 0, b"r1", &[5 << 4], 5, 30, &aux);
        let r = Rec::check(&b).unwrap();
        let mut out = Vec::new();
        let mut scratch = Vec::new();
        let edits = TagEdits { clear_dt: true, set_dt: None, set_pg: Some(Some(b"MarkDuplicates")), force: false };
        encode_record(&mut out, &r, 0, &[5 << 4], &edits, &mut scratch).unwrap();
        let rec2 = Rec(&out[4..]);
        let mut names = Vec::new();
        let mut at = 0;
        let aux2 = rec2.aux();
        while at < aux2.len() {
            let (t, v, next) = parse_tag(aux2, at).unwrap();
            names.push((String::from_utf8(t.to_vec()).unwrap(), aux2[at + 2] as char, v));
            at = next;
        }
        let order: Vec<&str> = names.iter().map(|n| n.0.as_str()).collect();
        assert_eq!(order, ["PG", "RG", "NH", "HI", "nM", "AS"]);
        // 1 and 0 as signed bytes, 300 as a signed short
        assert_eq!(names[2].1, 'c');
        assert_eq!(names[4].1, 'c');
        assert_eq!(names[5].1, 's');
        assert_eq!(names[5].2, TagValue::Int(300));
    }

    #[test]
    fn untouched_records_keep_their_fields() {
        let aux = b"NHC\x01RGZx\x00";
        let b = make_record(0, 10, 0, b"r1", &[5 << 4], 5, 30, aux);
        let r = Rec::check(&b).unwrap();
        let mut out = Vec::new();
        let mut scratch = Vec::new();
        encode_record(&mut out, &r, 0x400, &[5 << 4], &TagEdits::default(), &mut scratch).unwrap();
        assert_eq!(&out[4 + 32..], &b[32..]);
        assert_eq!(Rec(&out[4..]).flag(), 0x400);
        assert_eq!(u16_at(&out[4..], 10), region_to_bin(10, 15) as u16);
    }

    #[test]
    fn repeated_tags_keep_the_last_value() {
        let aux = b"XAi\x05\x00\x00\x00XAi\x07\x00\x00\x00";
        let b = make_record(0, 10, 0, b"r1", &[5 << 4], 5, 30, aux);
        let r = Rec::check(&b).unwrap();
        let mut out = Vec::new();
        let mut scratch = Vec::new();
        let edits = TagEdits { clear_dt: true, ..Default::default() };
        encode_record(&mut out, &r, 0, &[5 << 4], &edits, &mut scratch).unwrap();
        assert_eq!(Rec(&out[4..]).aux(), b"XAc\x07");
    }

    #[test]
    fn missing_qualities_and_odd_bases() {
        let mut b = make_record(0, 10, 0, b"r1", &[3 << 4], 3, 0xff, b"");
        // a stray low nibble in the last base byte, and qualities 0xff then others
        let seq_at = 32 + 3 + 4;
        b[seq_at + 1] = 0x4f;
        b[seq_at + 2 + 1] = 20;
        let r = Rec::check(&b).unwrap();
        let mut out = Vec::new();
        let mut scratch = Vec::new();
        encode_record(&mut out, &r, 0, &[3 << 4], &TagEdits { clear_dt: true, ..Default::default() }, &mut scratch).unwrap();
        let r2 = Rec(&out[4..]);
        assert_eq!(r2.qual(), &[0xff, 0xff, 0xff]);
        assert_eq!(out[4 + seq_at + 1], 0x40);
    }

    #[test]
    fn header_round_trip_as_htsjdk() {
        let text = b"@HD\tSO:coordinate\tVN:1.4\n@SQ\tSN:chr1\tLN:100\tM5:abc\n@PG\tID:STAR\tPN:STAR\n@RG\tID:a\tSM:s\tLB:l\n@CO\tuser command line: x\n";
        let h = SamHeader::parse(text).unwrap();
        assert_eq!(h.sort_order(), SortOrder::Coordinate);
        assert_eq!(
            String::from_utf8(h.encode()).unwrap(),
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\tM5:abc\n@RG\tID:a\tSM:s\tLB:l\n@PG\tID:STAR\tPN:STAR\n@CO\tuser command line: x\n"
        );
        assert!(SamHeader::parse(b"@XY\tID:1\n").is_err());
        assert!(SamHeader::parse(b"@RG\tID:1\n").is_err());
    }

    #[test]
    fn bins() {
        assert_eq!(region_to_bin(0, 1), 4681);
        assert_eq!(region_to_bin(16383, 16385), 585);
        assert_eq!(region_to_bin(-1, 0), 4680);
    }
}
