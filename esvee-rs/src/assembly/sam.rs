//! SAM record model with the parts of htsjdk's SAMRecord, Cigar and hmf-common CigarUtils / SamRecordUtils /
//! SupplementaryReadData that the assembly step uses.

use std::fmt::Write as _;
use std::sync::Mutex;

pub const NO_CHROMOSOME_NAME: &str = "*";
pub const NO_CIGAR: &str = "*";
pub const NO_POSITION: i32 = 0;
pub const INVALID_READ_INDEX: i32 = -1;

pub const MATE_CIGAR_ATTRIBUTE: &str = "MC";
pub const SUPPLEMENTARY_ATTRIBUTE: &str = "SA";
pub const ALIGNMENT_SCORE_ATTRIBUTE: &str = "AS";
pub const NUM_MUTATONS_ATTRIBUTE: &str = "NM";

pub const FLAG_READ_PAIRED: i32 = 1;
pub const FLAG_PROPER_PAIR: i32 = 2;
pub const FLAG_READ_UNMAPPED: i32 = 4;
pub const FLAG_MATE_UNMAPPED: i32 = 8;
pub const FLAG_READ_REVERSE_STRAND: i32 = 16;
pub const FLAG_MATE_REVERSE_STRAND: i32 = 32;
pub const FLAG_FIRST_OF_PAIR: i32 = 64;
pub const FLAG_SECOND_OF_PAIR: i32 = 128;
pub const FLAG_SECONDARY: i32 = 256;
pub const FLAG_FAILS_QC: i32 = 512;
pub const FLAG_DUPLICATE: i32 = 1024;
pub const FLAG_SUPPLEMENTARY: i32 = 2048;

/// Returns one shared `&'static str` per distinct name (chromosome names are few and used everywhere).
pub fn intern(name: &str) -> &'static str {
    static TABLE: Mutex<Vec<&'static str>> = Mutex::new(Vec::new());
    let mut table = TABLE.lock().unwrap();
    if let Some(s) = table.iter().find(|s| **s == name) {
        return s;
    }
    let s: &'static str = Box::leak(name.to_string().into_boxed_str());
    table.push(s);
    s
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CigarOp {
    M,
    I,
    D,
    N,
    S,
    H,
    P,
    EQ,
    X,
}

impl CigarOp {
    #[inline]
    pub fn consumes_read_bases(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::I | CigarOp::S | CigarOp::EQ | CigarOp::X)
    }
    #[inline]
    pub fn consumes_reference_bases(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::D | CigarOp::N | CigarOp::EQ | CigarOp::X)
    }
    #[inline]
    pub fn is_indel(self) -> bool {
        matches!(self, CigarOp::I | CigarOp::D)
    }
    #[inline]
    pub fn is_clipping(self) -> bool {
        matches!(self, CigarOp::S | CigarOp::H)
    }
    pub fn as_char(self) -> char {
        match self {
            CigarOp::M => 'M',
            CigarOp::I => 'I',
            CigarOp::D => 'D',
            CigarOp::N => 'N',
            CigarOp::S => 'S',
            CigarOp::H => 'H',
            CigarOp::P => 'P',
            CigarOp::EQ => '=',
            CigarOp::X => 'X',
        }
    }
    pub fn from_char(c: char) -> CigarOp {
        match c {
            'M' => CigarOp::M,
            'I' => CigarOp::I,
            'D' => CigarOp::D,
            'N' => CigarOp::N,
            'S' => CigarOp::S,
            'H' => CigarOp::H,
            'P' => CigarOp::P,
            '=' => CigarOp::EQ,
            'X' => CigarOp::X,
            _ => panic!("Unrecognized CigarOperator: {c}"),
        }
    }
    pub fn from_bam(op: u32) -> CigarOp {
        match op {
            0 => CigarOp::M,
            1 => CigarOp::I,
            2 => CigarOp::D,
            3 => CigarOp::N,
            4 => CigarOp::S,
            5 => CigarOp::H,
            6 => CigarOp::P,
            7 => CigarOp::EQ,
            _ => CigarOp::X,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct CigarElement {
    pub len: i32,
    pub op: CigarOp,
}

impl CigarElement {
    pub fn new(len: i32, op: CigarOp) -> Self {
        CigarElement { len, op }
    }
}

pub fn cigar_elements_from_str(cigar: &str) -> Vec<CigarElement> {
    if cigar == NO_CIGAR {
        return Vec::new();
    }
    let mut elements = Vec::new();
    let mut length = 0i32;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            length = length * 10 + (c as i32 - '0' as i32);
        } else {
            elements.push(CigarElement::new(length, CigarOp::from_char(c)));
            length = 0;
        }
    }
    elements
}

pub fn cigar_elements_to_str(elements: &[CigarElement]) -> String {
    let mut s = String::new();
    for e in elements {
        let _ = write!(s, "{}{}", e.len, e.op.as_char());
    }
    s
}

/// htsjdk Cigar.toString / SAMRecord.getCigarString: "*" when empty
pub fn cigar_string(elements: &[CigarElement]) -> String {
    if elements.is_empty() { NO_CIGAR.to_string() } else { cigar_elements_to_str(elements) }
}

pub fn cigar_base_length(elements: &[CigarElement]) -> i32 {
    elements.iter().filter(|x| x.op.consumes_read_bases()).map(|x| x.len).sum()
}

pub fn cigar_aligned_length(elements: &[CigarElement]) -> i32 {
    elements.iter().filter(|x| x.op.consumes_reference_bases()).map(|x| x.len).sum()
}

pub fn calc_cigar_aligned_length(cigar: &str) -> i32 {
    let mut base_length = 0;
    let mut current = 0;
    for c in cigar.chars() {
        if c == 'D' || c == 'M' || c == 'N' {
            base_length += current;
            current = 0;
            continue;
        }
        if c.is_ascii_digit() {
            current = current * 10 + (c as i32 - '0' as i32);
        } else {
            current = 0;
        }
    }
    base_length
}

pub fn left_soft_clip_length(elements: &[CigarElement]) -> i32 {
    match elements.first() {
        Some(e) if e.op == CigarOp::S => e.len,
        _ => 0,
    }
}

pub fn right_soft_clip_length(elements: &[CigarElement]) -> i32 {
    match elements.last() {
        Some(e) if e.op == CigarOp::S => e.len,
        _ => 0,
    }
}

pub fn left_clip_length(elements: &[CigarElement]) -> i32 {
    match elements.first() {
        Some(e) if e.op.is_clipping() => e.len,
        _ => 0,
    }
}

pub fn right_clip_length(elements: &[CigarElement]) -> i32 {
    match elements.last() {
        Some(e) if e.op.is_clipping() => e.len,
        _ => 0,
    }
}

/// CigarUtils.getReadIndexFromPosition(alignmentStart, elements, position, gapMode, inferFromSoftClip)
pub fn get_read_index_from_position_gap(alignment_start: i32, elements: &[CigarElement], position: i32, gap_mode: i32, infer_from_soft_clip: bool) -> i32 {
    if position < alignment_start && !infer_from_soft_clip {
        return INVALID_READ_INDEX;
    }
    let mut ref_position = alignment_start;
    let mut index = 0;
    for element in elements {
        if element.op == CigarOp::S && ref_position == alignment_start && infer_from_soft_clip {
            ref_position -= element.len;
        }
        let consumes_ref = element.op.consumes_reference_bases() || (element.op == CigarOp::S && infer_from_soft_clip);
        if consumes_ref && ref_position + element.len - 1 >= position {
            if element.op.consumes_read_bases() {
                index += position - ref_position;
            } else if gap_mode < 0 {
                index -= 1;
            } else if gap_mode == 0 {
                return INVALID_READ_INDEX;
            }
            return index;
        }
        if consumes_ref {
            ref_position += element.len;
        }
        if element.op.consumes_read_bases() {
            index += element.len;
        }
    }
    INVALID_READ_INDEX
}

pub fn get_read_index_from_position_ext(alignment_start: i32, elements: &[CigarElement], position: i32, last_if_gapped: bool, infer_from_soft_clip: bool) -> i32 {
    get_read_index_from_position_gap(alignment_start, elements, position, if last_if_gapped { -1 } else { 0 }, infer_from_soft_clip)
}

pub fn get_read_index_from_position(alignment_start: i32, elements: &[CigarElement], position: i32) -> i32 {
    get_read_index_from_position_ext(alignment_start, elements, position, false, false)
}

pub fn collapse_cigar_ops(ops: &[CigarOp]) -> Vec<CigarElement> {
    let mut elems = Vec::new();
    let mut current: Option<CigarOp> = None;
    let mut current_length = 0;
    for &op in ops {
        match current {
            None => {
                current = Some(op);
                current_length = 1;
            }
            Some(c) if c == op => current_length += 1,
            Some(c) => {
                elems.push(CigarElement::new(current_length, c));
                current = Some(op);
                current_length = 1;
            }
        }
    }
    if current_length > 0 {
        if let Some(c) = current {
            elems.push(CigarElement::new(current_length, c));
        }
    }
    elems
}

pub fn has_valid_cigar(elements: &[CigarElement]) -> bool {
    if elements.is_empty() {
        return false;
    }
    let count = elements.len();
    if count > 1 {
        let mut last = elements[0].op;
        if last != CigarOp::S && last != CigarOp::M {
            return false;
        }
        for (i, e) in elements.iter().enumerate().skip(1) {
            let op = e.op;
            if op == last {
                return false;
            }
            if last == CigarOp::S && op != CigarOp::M {
                return false;
            }
            if i == count - 1 && op != CigarOp::M && op != CigarOp::S {
                return false;
            }
            last = op;
        }
    } else if elements[0].op != CigarOp::M {
        return false;
    }
    true
}

/// CigarUtils.getReadBoundaryPosition
pub fn get_read_boundary_position(read_start: i32, cigar: &str, get_read_start: bool, include_soft_clipped: bool) -> i32 {
    if get_read_start && !include_soft_clipped {
        return read_start;
    }
    let mut current_position = read_start;
    let mut element_length = 0;
    for c in cigar.chars() {
        let is_add_item = c == 'D' || c == 'M' || c == 'S' || c == 'N';
        if is_add_item {
            if get_read_start {
                return if c == 'S' { read_start - element_length } else { read_start };
            }
            if c == 'S' && read_start == current_position {
                // ignore the left clip
            } else if c == 'S' && !include_soft_clipped {
                break;
            } else {
                current_position += element_length;
            }
            element_length = 0;
            continue;
        }
        if c.is_ascii_digit() {
            element_length = element_length * 10 + (c as i32 - '0' as i32);
        } else {
            element_length = 0;
        }
    }
    current_position - 1
}

pub fn get_alignment_end_from_cigar(read_start: i32, cigar: &str) -> i32 {
    get_read_boundary_position(read_start, cigar, false, false)
}

pub fn get_five_prime_unclipped_position(read_start: i32, cigar: &str, forward_strand: bool) -> i32 {
    get_read_boundary_position(read_start, cigar, forward_strand, true)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SupplementaryReadData {
    pub chromosome: &'static str,
    pub position: i32,
    pub strand: char,
    pub cigar: String,
    pub map_quality: i32,
    pub nm: i32,
}

pub const SUPP_POS_STRAND: char = '+';
pub const SUPP_NEG_STRAND: char = '-';

impl SupplementaryReadData {
    pub fn new(chromosome: &str, position: i32, strand: char, cigar: &str, map_quality: i32, nm: i32) -> Self {
        SupplementaryReadData { chromosome: intern(chromosome), position, strand, cigar: cigar.to_string(), map_quality, nm }
    }
    pub fn is_forward_orient(&self) -> bool {
        self.strand == SUPP_POS_STRAND
    }
    pub fn orientation(&self) -> i8 {
        if self.strand == SUPP_POS_STRAND { 1 } else { -1 }
    }
    pub fn from_alignment(alignment: &str, delim: char) -> Option<Self> {
        // String.split drops trailing empty strings
        let mut items: Vec<&str> = alignment.split(delim).collect();
        while items.last().is_some_and(|x| x.is_empty()) {
            items.pop();
        }
        if items.len() != 6 {
            return None;
        }
        Some(SupplementaryReadData::new(
            items[0],
            items[1].parse().expect("supplementary position"),
            items[2].chars().next().expect("supplementary strand"),
            items[3],
            items[4].parse().expect("supplementary map qual"),
            items[5].parse().expect("supplementary NM"),
        ))
    }
    /// SupplementaryReadData.extractAlignment(String): the first alignment in the tag
    pub fn extract_alignment(supp_data: Option<&str>) -> Option<Self> {
        let supp_data = supp_data?;
        if supp_data.is_empty() {
            return None;
        }
        if supp_data.contains(',') {
            let first = supp_data.split(';').next().unwrap();
            Self::from_alignment(first, ',')
        } else {
            Self::from_alignment(supp_data, ';')
        }
    }
    /// SupplementaryReadData.extractAlignments(String)
    pub fn extract_alignments(supp_data: Option<&str>) -> Option<Vec<Self>> {
        let supp_data = supp_data?;
        if supp_data.is_empty() {
            return None;
        }
        if supp_data.contains(',') {
            let alignments: Vec<&str> = supp_data.split(';').collect();
            let end_index = if supp_data.ends_with(';') { alignments.len() as isize - 2 } else { alignments.len() as isize - 1 };
            let mut output = Vec::new();
            let mut i = 0isize;
            while i <= end_index {
                output.push(Self::from_alignment(alignments[i as usize], ',')?);
                i += 1;
            }
            Some(output)
        } else {
            Some(vec![Self::from_alignment(supp_data, ';').expect("supplementary alignment")])
        }
    }
    pub fn as_sam_tag(&self) -> String {
        format!("{},{},{},{},{},{}", self.chromosome, self.position, self.strand, self.cigar, self.map_quality, self.nm)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum AttrValue {
    Int(i64),
    Str(String),
    Float(f32),
    /// B arrays and anything else, kept as the raw BAM value bytes after the type character
    Other(u8, Vec<u8>),
}

/// The fields of an htsjdk SAMRecord that the assembly step reads or changes.
#[derive(Clone, Debug)]
pub struct SamRecord {
    pub read_name: String,
    pub flags: i32,
    /// "*" when the record has no reference
    pub reference_name: &'static str,
    /// 1-based, 0 when unset
    pub alignment_start: i32,
    pub mapping_quality: i32,
    pub cigar: Vec<CigarElement>,
    pub mate_reference_name: &'static str,
    pub mate_alignment_start: i32,
    pub inferred_insert_size: i32,
    pub read_bases: Vec<u8>,
    pub base_qualities: Vec<u8>,
    /// raw BAM auxiliary data
    pub aux: Vec<u8>,
    /// the header attribute `sampleId` that the Java BamReader sets on each file's header
    pub sample_index: i32,
}

impl SamRecord {
    pub fn empty() -> Self {
        SamRecord {
            read_name: String::new(),
            flags: 0,
            reference_name: NO_CHROMOSOME_NAME,
            alignment_start: 0,
            mapping_quality: 0,
            cigar: Vec::new(),
            mate_reference_name: NO_CHROMOSOME_NAME,
            mate_alignment_start: 0,
            inferred_insert_size: 0,
            read_bases: Vec::new(),
            base_qualities: Vec::new(),
            aux: Vec::new(),
            sample_index: 0,
        }
    }
    #[inline]
    pub fn flag(&self, f: i32) -> bool {
        self.flags & f != 0
    }
    pub fn set_flag(&mut self, f: i32, on: bool) {
        if on {
            self.flags |= f;
        } else {
            self.flags &= !f;
        }
    }
    pub fn read_unmapped(&self) -> bool {
        self.flag(FLAG_READ_UNMAPPED)
    }
    pub fn read_paired(&self) -> bool {
        self.flag(FLAG_READ_PAIRED)
    }
    pub fn mate_unmapped(&self) -> bool {
        self.flag(FLAG_MATE_UNMAPPED)
    }
    pub fn read_negative_strand(&self) -> bool {
        self.flag(FLAG_READ_REVERSE_STRAND)
    }
    pub fn mate_negative_strand(&self) -> bool {
        self.flag(FLAG_MATE_REVERSE_STRAND)
    }
    pub fn first_of_pair(&self) -> bool {
        self.flag(FLAG_FIRST_OF_PAIR)
    }
    pub fn second_of_pair(&self) -> bool {
        self.flag(FLAG_SECOND_OF_PAIR)
    }
    pub fn supplementary(&self) -> bool {
        self.flag(FLAG_SUPPLEMENTARY)
    }
    pub fn secondary(&self) -> bool {
        self.flag(FLAG_SECONDARY)
    }
    pub fn duplicate(&self) -> bool {
        self.flag(FLAG_DUPLICATE)
    }
    pub fn cigar_string(&self) -> String {
        cigar_string(&self.cigar)
    }
    /// SAMRecord.getAlignmentEnd: start + reference length - 1, or 0 for an unmapped read
    pub fn alignment_end(&self) -> i32 {
        if self.read_unmapped() {
            return 0;
        }
        self.alignment_start + cigar_aligned_length(&self.cigar) - 1
    }
    pub fn unclipped_start(&self) -> i32 {
        let mut pos = self.alignment_start;
        for e in &self.cigar {
            if e.op.is_clipping() {
                pos -= e.len;
            } else {
                break;
            }
        }
        pos
    }
    pub fn unclipped_end(&self) -> i32 {
        let mut pos = self.alignment_end();
        for e in self.cigar.iter().rev() {
            if e.op.is_clipping() {
                pos += e.len;
            } else {
                break;
            }
        }
        pos
    }

    /// Finds an auxiliary tag; returns (type, offset of value, length of value) in `aux`.
    fn find_aux(&self, tag: &str) -> Option<(u8, usize, usize)> {
        let t = tag.as_bytes();
        let a = &self.aux;
        let mut i = 0;
        while i + 3 <= a.len() {
            let ty = a[i + 2];
            let start = i + 3;
            let len = match ty {
                b'A' | b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                b'Z' | b'H' => {
                    let mut j = start;
                    while j < a.len() && a[j] != 0 {
                        j += 1;
                    }
                    j - start + 1
                }
                b'B' => {
                    let sub = a[start];
                    let n = u32::from_le_bytes([a[start + 1], a[start + 2], a[start + 3], a[start + 4]]) as usize;
                    let size = match sub {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        _ => 4,
                    };
                    5 + n * size
                }
                _ => return None,
            };
            if a[i] == t[0] && a[i + 1] == t[1] {
                return Some((ty, start, len));
            }
            i = start + len;
        }
        None
    }
    pub fn has_attribute(&self, tag: &str) -> bool {
        self.find_aux(tag).is_some()
    }
    pub fn get_attribute(&self, tag: &str) -> Option<AttrValue> {
        let (ty, s, len) = self.find_aux(tag)?;
        let a = &self.aux;
        Some(match ty {
            b'A' => AttrValue::Str((a[s] as char).to_string()),
            b'c' => AttrValue::Int(a[s] as i8 as i64),
            b'C' => AttrValue::Int(a[s] as i64),
            b's' => AttrValue::Int(i16::from_le_bytes([a[s], a[s + 1]]) as i64),
            b'S' => AttrValue::Int(u16::from_le_bytes([a[s], a[s + 1]]) as i64),
            b'i' => AttrValue::Int(i32::from_le_bytes([a[s], a[s + 1], a[s + 2], a[s + 3]]) as i64),
            b'I' => AttrValue::Int(u32::from_le_bytes([a[s], a[s + 1], a[s + 2], a[s + 3]]) as i64),
            b'f' => AttrValue::Float(f32::from_le_bytes([a[s], a[s + 1], a[s + 2], a[s + 3]])),
            b'Z' | b'H' => AttrValue::Str(String::from_utf8_lossy(&a[s..s + len - 1]).into_owned()),
            _ => AttrValue::Other(ty, a[s..s + len].to_vec()),
        })
    }
    pub fn get_string_attribute(&self, tag: &str) -> Option<String> {
        match self.get_attribute(tag) {
            Some(AttrValue::Str(s)) => Some(s),
            _ => None,
        }
    }
    pub fn get_integer_attribute(&self, tag: &str) -> Option<i32> {
        match self.get_attribute(tag) {
            Some(AttrValue::Int(v)) => Some(v as i32),
            _ => None,
        }
    }
    /// raw value bytes of a B array (sub-type, elements) or any other tag
    pub fn get_raw_attribute(&self, tag: &str) -> Option<(u8, &[u8])> {
        let (ty, s, len) = self.find_aux(tag)?;
        Some((ty, &self.aux[s..s + len]))
    }
    fn remove_attribute(&mut self, tag: &str) {
        if let Some((_, s, len)) = self.find_aux(tag) {
            self.aux.drain(s - 3..s + len);
        }
    }
    pub fn set_string_attribute(&mut self, tag: &str, value: &str) {
        self.remove_attribute(tag);
        self.aux.extend_from_slice(tag.as_bytes());
        self.aux.push(b'Z');
        self.aux.extend_from_slice(value.as_bytes());
        self.aux.push(0);
    }
    pub fn set_integer_attribute(&mut self, tag: &str, value: i32) {
        self.remove_attribute(tag);
        self.aux.extend_from_slice(tag.as_bytes());
        self.aux.push(b'i');
        self.aux.extend_from_slice(&value.to_le_bytes());
    }
}

/// SamRecordUtils.inferredInsertSize
pub fn inferred_insert_size(record: &SamRecord) -> i32 {
    if record.read_paired() {
        return record.inferred_insert_size;
    }
    record.alignment_end() - record.alignment_start + 1
}

pub fn inferred_insert_size_abs(record: &SamRecord) -> i32 {
    inferred_insert_size(record).abs()
}

/// SamRecordUtils.getMateAlignmentEnd
pub fn get_mate_alignment_end(record: &SamRecord) -> i32 {
    match record.get_string_attribute(MATE_CIGAR_ATTRIBUTE) {
        Some(mc) if mc != NO_CIGAR => get_alignment_end_from_cigar(record.mate_alignment_start, &mc),
        _ => NO_POSITION,
    }
}
