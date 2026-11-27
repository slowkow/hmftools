//! Read, MappedCoords and the per-read maps, as the `common` package has them.
//!
//! A read's exon regions (Java `Map<RegionReadData,RegionMatchType>`, keyed by identity) iterate as a HashMap
//! of constant identity hashes does (`-XX:hashCode=2`, the defined Java reference): insertion order up to ten
//! entries, then the order of the bucket's red-black tree, whose keys are Comparable by region start
//! (`RegionMap`). Its transcript classifications (`Map<Integer,TransMatchType>`) iterate in Java's Integer
//! order (`IntJavaMap`).

use crate::common::*;
use hmf_common::java_id_hash::{IdHashMap, IdKey};
use hmf_common::regions::BaseRegion;

// BAM cigar operation codes
pub const OP_M: u8 = 0;
pub const OP_I: u8 = 1;
pub const OP_D: u8 = 2;
pub const OP_N: u8 = 3;
pub const OP_S: u8 = 4;
pub const OP_H: u8 = 5;
pub const OP_P: u8 = 6;
pub const OP_EQ: u8 = 7;
pub const OP_X: u8 = 8;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CigarElem {
    pub op: u8,
    pub len: u32,
}

pub fn op_char(op: u8) -> char {
    match op {
        OP_M => 'M',
        OP_I => 'I',
        OP_D => 'D',
        OP_N => 'N',
        OP_S => 'S',
        OP_H => 'H',
        OP_P => 'P',
        OP_EQ => '=',
        _ => 'X',
    }
}

pub fn consumes_ref(op: u8) -> bool {
    matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X)
}
pub fn consumes_read(op: u8) -> bool {
    matches!(op, OP_M | OP_I | OP_S | OP_EQ | OP_X)
}

pub fn cigar_to_string(elems: &[CigarElem]) -> String {
    if elems.is_empty() {
        return "*".to_string();
    }
    let mut s = String::with_capacity(elems.len() * 4);
    for e in elems {
        s.push_str(&e.len.to_string());
        s.push(op_char(e.op));
    }
    s
}

/// TextCigarCodec.decode: None for an invalid cigar
pub fn parse_cigar(s: &str) -> Option<Vec<CigarElem>> {
    if s == "*" {
        return Some(Vec::new());
    }
    let mut out = Vec::new();
    let mut n: u32 = 0;
    let mut have = false;
    for c in s.bytes() {
        if c.is_ascii_digit() {
            n = n.checked_mul(10)?.checked_add((c - b'0') as u32)?;
            have = true;
        } else {
            if !have {
                return None;
            }
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
                _ => return None,
            };
            out.push(CigarElem { op, len: n });
            n = 0;
            have = false;
        }
    }
    if have {
        return None;
    }
    Some(out)
}

pub fn reference_length(elems: &[CigarElem]) -> i32 {
    elems.iter().filter(|e| consumes_ref(e.op)).map(|e| e.len as i32).sum()
}

/// SupplementaryReadData: the first alignment of an SA tag
#[derive(Clone, Debug, PartialEq)]
pub struct SuppData {
    pub chromosome: String,
    pub position: i32,
    pub strand: char,
    pub cigar: String,
    pub map_quality: i32,
    pub nm: i32,
}

impl SuppData {
    /// SupplementaryReadData.extractAlignment(String)
    pub fn extract(sa: &str) -> Option<SuppData> {
        if sa.is_empty() {
            return None;
        }
        if sa.contains(',') {
            let first = sa.split(';').next().unwrap_or("");
            Self::from_alignment(first, ',')
        } else {
            Self::from_alignment(sa, ';')
        }
    }
    fn from_alignment(s: &str, delim: char) -> Option<SuppData> {
        // String.split drops trailing empty strings
        let items = hmf_common::java_text::split(s, delim);
        if items.len() != 6 {
            return None;
        }
        Some(SuppData {
            chromosome: items[0].to_string(),
            position: items[1].parse().ok()?,
            strand: items[2].chars().next()?,
            cigar: items[3].to_string(),
            map_quality: items[4].parse().ok()?,
            nm: items[5].parse().ok()?,
        })
    }
}

/// Read.AltAlignment: an alternate locus from the XA tag
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct AltAlignment {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub spliced: bool,
}

fn parse_alt_loci(xa: Option<&str>) -> Option<Vec<AltAlignment>> {
    let xa = xa?;
    if xa.is_empty() {
        return None;
    }
    let mut out = Vec::new();
    for entry in hmf_common::java_text::split(xa, ';') {
        if entry.is_empty() {
            continue;
        }
        let fields = hmf_common::java_text::split(entry, ',');
        if fields.len() < 2 {
            continue;
        }
        let Some(p) = hmf_common::java_text::parse_int(fields[1]) else { continue };
        let position = p.wrapping_abs();
        let mut ref_len = 1;
        let mut spliced = false;
        if let Some(c) = fields.get(2) {
            if !c.is_empty() {
                // TextCigarCodec.decode throws on an invalid cigar: Java does not catch that
                let elems = parse_cigar(c).unwrap_or_default();
                ref_len = reference_length(&elems).max(1);
                spliced = c.contains('N');
            }
        }
        out.push(AltAlignment { chromosome: fields[0].to_string(), start: position, end: position + ref_len - 1, spliced });
    }
    if out.is_empty() { None } else { Some(out) }
}

/// MappedCoords: the read's aligned sections, with any inferred sections added from soft-clipped bases
#[derive(Clone, Debug, Default)]
pub struct MappedCoords {
    original_count: usize,
    pub alignments: Vec<BaseRegion>,
    inferred_added: Option<[bool; 2]>,
    soft_clip_matched: Option<[i32; 2]>,
}

pub const INVALID_INDEX: i32 = -1;

impl MappedCoords {
    pub fn build(elems: &[CigarElem], pos_start: i32) -> MappedCoords {
        let mut alignments: Vec<BaseRegion> = Vec::with_capacity(2);
        let mut have_current = false;
        let mut pos_offset: i32 = 0;
        let mut continue_region = false;
        for e in elems {
            match e.op {
                OP_S => continue,
                OP_I => {
                    continue_region = true;
                    continue;
                }
                OP_D => {
                    pos_offset += e.len as i32;
                    continue_region = true;
                }
                OP_N => {
                    pos_offset += e.len as i32;
                    continue_region = false;
                }
                OP_M => {
                    let s = pos_start + pos_offset;
                    let end = s + e.len as i32 - 1;
                    if continue_region && have_current {
                        alignments.last_mut().unwrap().end = end;
                    } else {
                        alignments.push(BaseRegion::new(s, end));
                        have_current = true;
                    }
                    pos_offset += e.len as i32;
                    continue_region = false;
                }
                _ => {}
            }
        }
        MappedCoords { original_count: alignments.len(), alignments, inferred_added: None, soft_clip_matched: None }
    }
    pub fn alignment_count(&self) -> usize {
        self.alignments.len()
    }
    pub fn original_alignment_count(&self) -> usize {
        self.original_count
    }
    /// alignmentsWithoutInferred
    pub fn without_inferred(&self) -> Vec<BaseRegion> {
        match self.inferred_added {
            None => self.alignments.clone(),
            Some(inf) => {
                let start = if inf[0] { 1 } else { 0 };
                self.alignments[start..start + self.original_count].to_vec()
            }
        }
    }
    /// index of lowestAlignment(false)
    pub fn lowest_index(&self, include_inferred: bool) -> usize {
        if include_inferred || !self.inferred(SE_START) { 0 } else { 1 }
    }
    pub fn highest_index(&self, include_inferred: bool) -> usize {
        let n = self.alignments.len();
        if include_inferred || !self.inferred(SE_END) { n - 1 } else { n - 2 }
    }
    pub fn find_region_index(&self, start: i32, end: i32) -> i32 {
        for (i, a) in self.alignments.iter().enumerate() {
            if positions_overlap(a.start, a.end, start, end) {
                return i as i32;
            }
        }
        INVALID_INDEX
    }
    pub fn inferred(&self, se: usize) -> bool {
        self.inferred_added.map(|x| x[se]).unwrap_or(false)
    }
    pub fn add_inferred_region(&mut self, is_lower: bool, pos_start: i32, pos_end: i32) {
        let inf = self.inferred_added.get_or_insert([false, false]);
        if is_lower {
            if !inf[0] {
                inf[0] = true;
                self.alignments.insert(0, BaseRegion::new(pos_start, pos_end));
            } else {
                let e = &mut self.alignments[0];
                e.start = e.start.min(pos_start);
            }
        } else if !inf[1] {
            inf[1] = true;
            self.alignments.push(BaseRegion::new(pos_start, pos_end));
        } else {
            let e = self.alignments.last_mut().unwrap();
            e.end = e.end.max(pos_end);
        }
    }
    pub fn alignments_overlap(&self, start: i32, end: i32) -> bool {
        self.alignments.iter().any(|x| positions_overlap(start, end, x.start, x.end))
    }
    pub fn alignments_within(&self, start: i32, end: i32) -> bool {
        self.alignments.iter().any(|x| positions_within(x.start, x.end, start, end))
    }
    pub fn coords_boundary(&self, se: usize) -> i32 {
        if se == SE_START { self.alignments[0].start } else { self.alignments[self.alignments.len() - 1].end }
    }
    pub fn add_soft_clip_region_matched(&mut self, on_start: bool, count: i32) {
        let m = self.soft_clip_matched.get_or_insert([0, 0]);
        m[if on_start { 0 } else { 1 }] += count;
    }
    pub fn is_soft_clip_region_matched(&self, se: usize) -> bool {
        self.soft_clip_matched.map(|m| m[se] > 0).unwrap_or(false)
    }
}

/// A region key: the index of the RegionReadData in its gene collection, compared by identity, and its start
/// for RegionReadData.compareTo
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RKey {
    pub idx: u32,
    pub start: i32,
}

impl IdKey for RKey {
    fn id_addr(&self) -> usize {
        self.idx as usize
    }
    fn compare_key(&self, other: &Self) -> i32 {
        self.start.wrapping_sub(other.start)
    }
}

/// a read's regions and their match types, in Java's iteration order (see the module comment)
#[derive(Clone, Default)]
pub struct RegionMap {
    small: Vec<(RKey, RegionMatchType)>,
    big: Option<Box<IdHashMap<RKey, RegionMatchType>>>,
}

/// a HashMap's single bucket is a list up to ten entries; the eleventh makes it a tree
const LIST_MAX: usize = 10;

#[allow(dead_code)]
impl RegionMap {
    pub fn is_empty(&self) -> bool {
        match &self.big {
            Some(b) => b.is_empty(),
            None => self.small.is_empty(),
        }
    }
    pub fn len(&self) -> usize {
        match &self.big {
            Some(b) => b.len(),
            None => self.small.len(),
        }
    }
    pub fn get(&self, idx: usize) -> Option<RegionMatchType> {
        match &self.big {
            Some(b) => b.iter().find(|(k, _)| k.idx as usize == idx).map(|(_, v)| *v),
            None => self.small.iter().find(|(k, _)| k.idx as usize == idx).map(|(_, v)| *v),
        }
    }
    pub fn contains(&self, idx: usize) -> bool {
        self.get(idx).is_some()
    }
    pub fn put(&mut self, key: RKey, value: RegionMatchType) {
        if let Some(b) = &mut self.big {
            b.insert(key, value);
            return;
        }
        if let Some(e) = self.small.iter_mut().find(|(k, _)| k.idx == key.idx) {
            e.1 = value;
            return;
        }
        if self.small.len() < LIST_MAX {
            self.small.push((key, value));
            return;
        }
        let mut b: IdHashMap<RKey, RegionMatchType> = IdHashMap::new();
        for (k, v) in self.small.drain(..) {
            b.insert(k, v);
        }
        b.insert(key, value);
        self.big = Some(Box::new(b));
    }
    /// entries in Java iteration order
    pub fn entries(&self) -> Vec<(usize, RegionMatchType)> {
        match &self.big {
            Some(b) => b.iter().map(|(k, v)| (k.idx as usize, *v)).collect(),
            None => self.small.iter().map(|(k, v)| (k.idx as usize, *v)).collect(),
        }
    }
    pub fn keys(&self) -> Vec<usize> {
        self.entries().into_iter().map(|(k, _)| k).collect()
    }
    pub fn values(&self) -> Vec<RegionMatchType> {
        self.entries().into_iter().map(|(_, v)| v).collect()
    }
}

#[inline]
fn spread(h: i32) -> u32 {
    let h = h as u32;
    h ^ (h >> 16)
}

/// a HashMap<Integer,V> (or HashSet<Integer> with V = ()) with Java's iteration order, for small maps
#[derive(Clone, Debug)]
pub struct IntJavaMap<V> {
    entries: Vec<(i32, V, u32)>,
    next_seq: u32,
    capacity: usize,
    initial: usize,
}

impl<V> Default for IntJavaMap<V> {
    fn default() -> Self {
        Self::new()
    }
}

#[allow(dead_code)]
impl<V> IntJavaMap<V> {
    pub fn new() -> Self {
        IntJavaMap { entries: Vec::new(), next_seq: 0, capacity: 0, initial: 16 }
    }
    pub fn with_capacity(cap: usize) -> Self {
        IntJavaMap { entries: Vec::new(), next_seq: 0, capacity: 0, initial: cap.max(1).next_power_of_two() }
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn get(&self, k: i32) -> Option<&V> {
        self.entries.iter().find(|e| e.0 == k).map(|e| &e.1)
    }
    pub fn get_mut(&mut self, k: i32) -> Option<&mut V> {
        self.entries.iter_mut().find(|e| e.0 == k).map(|e| &mut e.1)
    }
    pub fn contains_key(&self, k: i32) -> bool {
        self.entries.iter().any(|e| e.0 == k)
    }
    pub fn insert(&mut self, k: i32, v: V) {
        if let Some(e) = self.entries.iter_mut().find(|e| e.0 == k) {
            e.1 = v;
            return;
        }
        if self.capacity == 0 {
            self.capacity = self.initial;
        }
        self.entries.push((k, v, self.next_seq));
        self.next_seq += 1;
        while self.entries.len() > self.capacity * 3 / 4 {
            self.capacity *= 2;
        }
    }
    pub fn remove(&mut self, k: i32) -> Option<V> {
        let i = self.entries.iter().position(|e| e.0 == k)?;
        Some(self.entries.remove(i).1)
    }
    fn order(&self) -> Vec<usize> {
        let mask = (self.capacity.max(1) - 1) as u32;
        let mut idx: Vec<usize> = (0..self.entries.len()).collect();
        idx.sort_unstable_by_key(|i| (spread(self.entries[*i].0) & mask, self.entries[*i].2));
        idx
    }
    /// keys in Java iteration order
    pub fn keys(&self) -> Vec<i32> {
        self.order().into_iter().map(|i| self.entries[i].0).collect()
    }
    pub fn iter(&self) -> Vec<(i32, &V)> {
        self.order().into_iter().map(|i| (self.entries[i].0, &self.entries[i].1)).collect()
    }
    pub fn values_mut_unordered(&mut self) -> impl Iterator<Item = (i32, &mut V)> {
        self.entries.iter_mut().map(|e| (e.0, &mut e.1))
    }
    pub fn values_unordered(&self) -> impl Iterator<Item = &V> {
        self.entries.iter().map(|e| &e.1)
    }
}

/// a BAM record as ISOFOX uses it
#[derive(Clone)]
pub struct Read {
    pub id: Box<str>,
    pub chromosome: &'static str,
    pub pos_start: i32,
    pub pos_end: i32,
    pub flags: u16,
    pub mate_chromosome: &'static str,
    pub mate_pos: i32,
    pub insert_size: i32,
    pub mapq: u8,
    pub cigar: Vec<CigarElem>,
    /// the cigar before adapter trimming, if trimmed
    orig_cigar: Option<Vec<CigarElem>>,
    pub bases: Vec<u8>,
    /// the trimmed bases: an offset and length into `bases`
    base_view: Option<(usize, usize)>,
    pub unclipped_start: i32,
    pub unclipped_end: i32,
    pub has_split: bool,
    pub gene_collections: [i32; 2],
    pub is_genic: [bool; 2],
    pub mapped_coords: MappedCoords,
    pub supp_data: Option<SuppData>,
    pub has_inter_gene_split: bool,
    pub alt_loci: Option<Vec<AltAlignment>>,
    /// the CR (consensus read) tag
    pub consensus: Option<Box<str>>,
    pub junction_positions: Option<[i32; 2]>,
    pub mapped_regions: RegionMap,
    pub trans_class: IntJavaMap<TransMatchType>,
    /// the INTRON transcript-exon references (getReadTransExonRefs, which only ever holds INTRON)
    pub intron_trans_refs: Option<Vec<TransExonRef>>,
}

pub const FLAG_PAIRED: u16 = 0x1;
#[allow(dead_code)]
pub const FLAG_PROPER: u16 = 0x2;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST: u16 = 0x40;
#[allow(dead_code)]
pub const FLAG_SECOND: u16 = 0x80;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[allow(clippy::too_many_arguments)]
impl Read {
    pub fn new(
        id: &str,
        chromosome: &'static str,
        pos_start: i32,
        flags: u16,
        cigar: Vec<CigarElem>,
        bases: Vec<u8>,
        mate_chromosome: &'static str,
        mate_pos: i32,
        insert_size: i32,
        mapq: u8,
        sa: Option<&str>,
        xa: Option<&str>,
        cr: Option<&str>,
    ) -> Read {
        let pos_end = pos_start + reference_length(&cigar) - 1;
        let has_split = cigar.iter().any(|e| e.op == OP_N);
        let mapped_coords = MappedCoords::build(&cigar, pos_start);
        let mut r = Read {
            id: id.into(),
            chromosome,
            pos_start,
            pos_end,
            flags,
            mate_chromosome,
            mate_pos,
            insert_size,
            mapq,
            cigar,
            orig_cigar: None,
            bases,
            base_view: None,
            unclipped_start: 0,
            unclipped_end: 0,
            has_split,
            gene_collections: [NO_GENE_ID, NO_GENE_ID],
            is_genic: [false, false],
            mapped_coords,
            supp_data: sa.and_then(SuppData::extract),
            has_inter_gene_split: false,
            alt_loci: parse_alt_loci(xa),
            consensus: cr.map(|s| s.into()),
            junction_positions: None,
            mapped_regions: RegionMap::default(),
            trans_class: IntJavaMap::new(),
            intron_trans_refs: None,
        };
        r.set_boundaries();
        r
    }

    fn set_boundaries(&mut self) {
        self.unclipped_start = self.pos_start;
        if let Some(f) = self.cigar.first() {
            if f.op == OP_S {
                self.unclipped_start -= f.len as i32;
            }
        }
        self.unclipped_end = self.pos_end;
        if self.cigar.len() >= 2 {
            let l = self.cigar[self.cigar.len() - 1];
            if l.op == OP_S {
                self.unclipped_end += l.len as i32;
            }
        }
    }

    pub fn is_read_paired(&self) -> bool {
        self.flags & FLAG_PAIRED != 0
    }
    pub fn is_read_reversed(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }
    pub fn is_first_of_pair(&self) -> bool {
        !self.is_read_paired() || self.flags & FLAG_FIRST != 0
    }
    pub fn is_duplicate(&self) -> bool {
        self.flags & FLAG_DUPLICATE != 0
    }
    pub fn is_supplementary(&self) -> bool {
        self.flags & FLAG_SUPPLEMENTARY != 0
    }
    pub fn is_mate_neg_strand(&self) -> bool {
        self.is_read_paired() && self.flags & FLAG_MATE_REVERSE != 0
    }
    pub fn is_mate_unmapped(&self) -> bool {
        self.is_read_paired() && self.flags & FLAG_MATE_UNMAPPED != 0
    }
    pub fn is_translocation(&self) -> bool {
        self.is_read_paired() && self.chromosome != self.mate_chromosome
    }
    pub fn is_inversion(&self) -> bool {
        self.is_read_paired() && self.is_read_reversed() == self.is_mate_neg_strand()
    }
    pub fn orient(&self) -> i8 {
        if !self.is_read_reversed() { ORIENT_FWD } else { ORIENT_REV }
    }
    pub fn has_supp_alignment(&self) -> bool {
        self.supp_data.is_some()
    }
    pub fn is_chimeric(&self) -> bool {
        self.is_translocation() || self.is_inversion() || self.is_supplementary() || self.supp_data.is_some()
    }
    /// SamRecordUtils.inferredInsertSize
    pub fn fragment_insert_size(&self) -> i32 {
        if self.is_read_paired() { self.insert_size } else { self.pos_end - self.pos_start + 1 }
    }
    pub fn num_loci(&self) -> usize {
        self.alt_loci.as_ref().map(|a| 1 + a.len()).unwrap_or(1)
    }
    pub fn is_multi_mapped(&self) -> bool {
        self.num_loci() > 1
    }
    pub fn read_bases(&self) -> &[u8] {
        match self.base_view {
            Some((o, l)) => &self.bases[o..o + l],
            None => &self.bases,
        }
    }
    pub fn base_length(&self) -> i32 {
        self.read_bases().len() as i32
    }
    pub fn cigar_str(&self) -> String {
        cigar_to_string(&self.cigar)
    }
    pub fn is_left_clipped(&self) -> bool {
        self.unclipped_start != self.pos_start
    }
    pub fn is_right_clipped(&self) -> bool {
        self.unclipped_end != self.pos_end
    }
    pub fn contains_soft_clipping(&self) -> bool {
        self.is_left_clipped() || self.is_right_clipped()
    }
    pub fn left_clip_length(&self) -> i32 {
        (self.pos_start - self.unclipped_start).max(0)
    }
    pub fn right_clip_length(&self) -> i32 {
        (self.unclipped_end - self.pos_end).max(0)
    }
    #[allow(dead_code)]
    pub fn longest_soft_clip(&self) -> i32 {
        self.left_clip_length().max(self.right_clip_length())
    }
    pub fn longest_soft_clip_is_left(&self) -> Option<bool> {
        let (l, r) = (self.left_clip_length(), self.right_clip_length());
        if l == 0 && r == 0 {
            return None;
        }
        Some(l >= r)
    }
    pub fn is_soft_clipped_no_region_match(&self, se: usize) -> bool {
        if self.mapped_coords.is_soft_clip_region_matched(se) {
            return false;
        }
        if se == SE_START { self.is_left_clipped() } else { self.is_right_clipped() }
    }
    pub fn likely_adapter_soft_clipping(&self) -> bool {
        self.fragment_insert_size() < self.base_length()
    }
    pub fn set_gene_collection(&mut self, se: usize, gc: i32, is_genic: bool) {
        self.gene_collections[se] = gc;
        self.is_genic[se] = is_genic;
    }
    pub fn within_gene_collection(&self) -> bool {
        self.is_genic[0] && self.is_genic[1]
    }
    pub fn overlaps_gene_collection(&self) -> bool {
        self.is_genic[0] || self.is_genic[1]
    }
    pub fn fully_non_genic(&self) -> bool {
        !self.is_genic[0] && !self.is_genic[1] && self.gene_collections[0] != NO_GENE_ID && self.gene_collections[1] != NO_GENE_ID
    }
    pub fn matches(&self, other_id: &str, other_flags: u16) -> bool {
        &*self.id == other_id && self.flags == other_flags
    }
    pub fn spans_gene_collections(&self) -> bool {
        self.gene_collections[0] != self.gene_collections[1]
    }
    pub fn coords_boundary(&self, se: usize) -> i32 {
        self.mapped_coords.coords_boundary(se)
    }
    pub fn mapped_region_count(&self) -> usize {
        self.mapped_coords.original_alignment_count()
    }
    pub fn trans_classification(&self, trans_id: i32) -> TransMatchType {
        self.trans_class.get(trans_id).copied().unwrap_or(TransMatchType::Unknown)
    }
    pub fn set_junction_position(&mut self, se: usize, pos: i32) {
        let j = self.junction_positions.get_or_insert([0, 0]);
        j[se] = pos;
    }

    /// trimAdapterSoftClipBases
    pub fn trim_adapter_soft_clip_bases(&mut self, trim_length: i32) {
        let forward = !self.is_read_reversed();
        let sc_len = if forward { self.right_clip_length() } else { self.left_clip_length() };
        let trim = trim_length.min(sc_len);
        if trim == 0 {
            return;
        }
        if self.orig_cigar.is_none() {
            self.orig_cigar = Some(self.cigar.clone());
        }
        let n = self.bases.len();
        if forward {
            self.base_view = Some((0, n - trim as usize));
            self.unclipped_end -= trim;
            let last = self.cigar.len() - 1;
            if trim < sc_len {
                self.cigar[last] = CigarElem { op: OP_S, len: (sc_len - trim) as u32 };
            } else {
                self.cigar.remove(last);
            }
        } else {
            self.base_view = Some((trim as usize, n - trim as usize));
            self.unclipped_start += trim;
            if trim < sc_len {
                self.cigar[0] = CigarElem { op: OP_S, len: (sc_len - trim) as u32 };
            } else {
                self.cigar.remove(0);
            }
        }
    }
}

/// CigarUtils.getReadIndexFromPosition(alignmentStart, cigarElements, position) with gaps invalid
pub fn read_index_from_position(alignment_start: i32, elems: &[CigarElem], position: i32) -> i32 {
    if position < alignment_start {
        return -1;
    }
    let mut ref_pos = alignment_start;
    let mut index = 0;
    for e in elems {
        let consumes = consumes_ref(e.op);
        if consumes && ref_pos + e.len as i32 - 1 >= position {
            if consumes_read(e.op) {
                index += position - ref_pos;
            } else {
                return -1;
            }
            return index;
        }
        if consumes {
            ref_pos += e.len as i32;
        }
        if consumes_read(e.op) {
            index += e.len as i32;
        }
    }
    -1
}

/// ReadUtils.trimAdapterBases: trims soft-clipped bases beyond the other read's 5' end
pub fn trim_adapter_bases(read1: &mut Read, read2: &mut Read) {
    if read1.orient() == read2.orient() {
        return;
    }
    if !positions_overlap(read1.pos_start, read1.pos_end, read2.pos_start, read2.pos_end) {
        return;
    }
    let mut ref_start = 0;
    'outer: for c1 in &read1.mapped_coords.alignments {
        for c2 in &read2.mapped_coords.alignments {
            if positions_overlap(c1.start, c1.end, c2.start, c2.end) {
                ref_start = c1.start.max(c2.start);
                break 'outer;
            }
        }
    }
    if ref_start == 0 {
        return;
    }
    let i1 = read_index_from_position(read1.pos_start, &read1.cigar, ref_start);
    let i2 = read_index_from_position(read2.pos_start, &read2.cigar, ref_start);
    let upper1 = read1.base_length() - i1 - 1;
    let upper2 = read2.base_length() - i2 - 1;
    let (t1, t2) = if !read1.is_read_reversed() {
        ((upper1 - upper2).max(0), (i2 - i1).max(0))
    } else {
        ((i1 - i2).max(0), (upper2 - upper1).max(0))
    };
    read1.trim_adapter_soft_clip_bases(t1);
    read2.trim_adapter_soft_clip_bases(t2);
}

/// ReadUtils.consensusDuplicateCount from the CR tag
pub fn consensus_duplicate_count(cr: Option<&str>) -> i32 {
    let Some(info) = cr else { return 0 };
    let first = info.split(';').next().unwrap_or("");
    // Integer.parseInt throws on a bad value: Java does not catch it
    let dup: i32 = first.parse().unwrap_or(0);
    (dup - 1).max(1)
}

/// ReadUtils.clippedSide: (side, length, isSoft); never null
pub fn clipped_side(read: &Read) -> (usize, i32) {
    let l = read.left_clip_length();
    let r = read.right_clip_length();
    if l > 0 && r > 0 {
        if l >= r { (SE_START, l) } else { (SE_END, r) }
    } else if l > 0 {
        (SE_START, l)
    } else {
        (SE_END, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cig(s: &str) -> Vec<CigarElem> {
        parse_cigar(s).unwrap()
    }

    #[test]
    fn mapped_coords() {
        let mc = MappedCoords::build(&cig("10S20M100N30M5D10M2I5M"), 100);
        // a deletion and an insertion continue the section
        assert_eq!(mc.alignments, vec![BaseRegion::new(100, 119), BaseRegion::new(220, 269)]);
        let r = Read::new("r", "1", 100, 0x41, cig("10S20M100N30M5D10M2I5M5S"), vec![b'A'; 82], "1", 300, 0, 60, None, None, None);
        assert_eq!(r.pos_end, 269);
        assert_eq!(r.unclipped_start, 90);
        assert_eq!(r.unclipped_end, 274);
        assert!(r.has_split);
    }

    #[test]
    fn supp_data_parse() {
        let s = SuppData::extract("chr2,33141317,+,94S57M,5,0;chr3,10,-,10M,0,0;").unwrap();
        assert_eq!(s.chromosome, "chr2");
        assert_eq!(s.position, 33141317);
        assert_eq!(s.cigar, "94S57M");
        assert!(SuppData::extract("chr2,1,+,10M,5").is_none());
    }

    #[test]
    fn int_java_map_order() {
        let mut m: IntJavaMap<()> = IntJavaMap::new();
        for k in [200001, 17, 1, 33] {
            m.insert(k, ());
        }
        // spread(200001) = 200001 ^ 3 = 200002, bucket 2; 17 -> 1, 1 -> 1, 33 -> 1
        assert_eq!(m.keys(), vec![17, 1, 33, 200001]);
    }

    #[test]
    fn region_map_list_order() {
        let mut m = RegionMap::default();
        for (i, s) in [(3u32, 30), (1, 10), (2, 20)] {
            m.put(RKey { idx: i, start: s }, RegionMatchType::WithinExon);
        }
        m.put(RKey { idx: 1, start: 10 }, RegionMatchType::ExonBoundary);
        assert_eq!(m.keys(), vec![3, 1, 2]);
        assert_eq!(m.get(1), Some(RegionMatchType::ExonBoundary));
    }

    #[test]
    fn read_index() {
        let c = cig("5S20M100N30M");
        assert_eq!(read_index_from_position(100, &c, 100), 5);
        assert_eq!(read_index_from_position(100, &c, 119), 24);
        assert_eq!(read_index_from_position(100, &c, 150), -1);
        assert_eq!(read_index_from_position(100, &c, 220), 25);
    }

    #[test]
    fn consensus_counts() {
        assert_eq!(consensus_duplicate_count(None), 0);
        assert_eq!(consensus_duplicate_count(Some("2;x;y")), 1);
        assert_eq!(consensus_duplicate_count(Some("5;x")), 4);
    }
}
