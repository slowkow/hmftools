//! Reads, read groups and junctions of a partition (Java package esvee.prep.types).
//!
//! The Java objects refer to each other freely; here a partition keeps its reads, groups and junctions in
//! vectors and the references are indices into them (`ReadId`, `GroupId`, `JuncId`).

use super::jmap::{set_order, string_hash};
use super::rec::*;
use rust_htslib::bam;
use std::collections::HashMap;

pub type ReadId = u32;
pub type GroupId = u32;
pub type JuncId = u32;

#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum Orientation {
    Forward,
    Reverse,
}

impl Orientation {
    pub fn is_forward(self) -> bool {
        self == Orientation::Forward
    }
    pub fn is_reverse(self) -> bool {
        self == Orientation::Reverse
    }
    pub fn as_byte(self) -> i32 {
        if self.is_forward() { 1 } else { -1 }
    }
    pub fn opposite(self) -> Orientation {
        if self.is_forward() { Orientation::Reverse } else { Orientation::Forward }
    }
    /// Orientation.fromChar
    pub fn from_char(c: u8) -> Option<Orientation> {
        match c {
            b'+' => Some(Orientation::Forward),
            b'-' => Some(Orientation::Reverse),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ReadType {
    NoSupport,
    CandidateSupport,
    Support,
    ExactSupport,
    Junction,
    Expected,
    #[allow(dead_code)]
    Recovered,
}

impl ReadType {
    pub fn rank(self) -> i32 {
        match self {
            ReadType::Junction => 7,
            ReadType::ExactSupport => 6,
            ReadType::Support => 5,
            ReadType::CandidateSupport => 4,
            _ => 0,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            ReadType::NoSupport => "NO_SUPPORT",
            ReadType::CandidateSupport => "CANDIDATE_SUPPORT",
            ReadType::Support => "SUPPORT",
            ReadType::ExactSupport => "EXACT_SUPPORT",
            ReadType::Junction => "JUNCTION",
            ReadType::Expected => "EXPECTED",
            ReadType::Recovered => "RECOVERED",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ReadGroupStatus {
    Unset,
    Incomplete,
    Expected,
    Paired,
    Supplementary,
    Complete,
    Duplicate,
}

impl ReadGroupStatus {
    pub fn name(self) -> &'static str {
        match self {
            ReadGroupStatus::Unset => "UNSET",
            ReadGroupStatus::Incomplete => "INCOMPLETE",
            ReadGroupStatus::Expected => "EXPECTED",
            ReadGroupStatus::Paired => "PAIRED",
            ReadGroupStatus::Supplementary => "SUPPLEMENTARY",
            ReadGroupStatus::Complete => "COMPLETE",
            ReadGroupStatus::Duplicate => "DUPLICATE",
        }
    }
}

/// ReadFilterType flags
pub const F_MIN_ALIGN_MATCH: i32 = 1;
pub const F_MIN_MAP_QUAL: i32 = 2;
pub const F_INSERT_MAP_OVERLAP: i32 = 4;
pub const F_SOFT_CLIP_LENGTH: i32 = 8;
pub const F_SOFT_CLIP_BASE_QUAL: i32 = 16;
pub const F_BREAK_IN_REPEAT: i32 = 32;
pub const F_POLY_G_SC: i32 = 64;
pub const F_SOFT_CLIP_LOW_BASE_QUAL: i32 = 128;
pub const FILTER_NAMES: [&str; 8] = [
    "MIN_ALIGN_MATCH",
    "MIN_MAP_QUAL",
    "INSERT_MAP_OVERLAP",
    "SOFT_CLIP_LENGTH",
    "SOFT_CLIP_BASE_QUAL",
    "BREAK_IN_REPEAT",
    "POLY_G_SC",
    "SOFT_CLIP_LOW_BASE_QUAL",
];

/// PrepRead.UNMAPPED_CHR ("-1"): the chromosome of an unmapped read or of a missing mate
pub const UNMAPPED_CHR: i32 = -1;
/// a contig name that is in none of the input BAM headers
pub const UNKNOWN_CHR: i32 = -2;

fn strip_chr(contig: &str) -> &str {
    contig.strip_prefix("chr").unwrap_or(contig)
}

/// HumanChromosome.contains: 1 to 22 as a number, X or Y, after a "chr" prefix
pub fn is_human_chromosome(contig: &str) -> bool {
    let t = strip_chr(contig);
    if !t.is_empty() && t.bytes().all(|c| c.is_ascii_digit()) {
        return t.parse::<i64>().is_ok_and(|n| (1..=22).contains(&n));
    }
    t == "X" || t == "Y"
}

/// HumanChromosome.chromosomeRank: 1-22, X 23, Y 24, MT 25, any other contig 26
pub fn chromosome_rank(chr: &str) -> i32 {
    let bare = strip_chr(chr);
    if bare.eq_ignore_ascii_case("X") {
        23
    } else if bare.eq_ignore_ascii_case("Y") {
        24
    } else if bare.eq_ignore_ascii_case("MT") || bare.eq_ignore_ascii_case("M") {
        25
    } else {
        bare.parse::<i32>().unwrap_or(26)
    }
}

/// The contigs of the input BAMs under one numbering. The Java code carries chromosome names as strings;
/// here a chromosome is an index into this table, with `UNMAPPED_CHR` for Java's "-1".
#[derive(Default, Debug)]
pub struct Contigs {
    pub names: Vec<String>,
    pub human: Vec<bool>,
    index: HashMap<String, i32>,
}

impl Contigs {
    pub fn new() -> Self {
        Self::default()
    }
    /// The index of a name, added when new.
    pub fn add(&mut self, name: &str) -> i32 {
        if let Some(&i) = self.index.get(name) {
            return i;
        }
        let i = self.names.len() as i32;
        self.names.push(name.to_string());
        self.human.push(is_human_chromosome(name));
        self.index.insert(name.to_string(), i);
        i
    }
    pub fn id(&self, name: &str) -> i32 {
        self.index.get(name).copied().unwrap_or(UNKNOWN_CHR)
    }
    /// The name as the Java code has it: "-1" for an unmapped read.
    pub fn name(&self, id: i32) -> &str {
        if id >= 0 { &self.names[id as usize] } else { "-1" }
    }
    pub fn is_human(&self, id: i32) -> bool {
        id >= 0 && self.human[id as usize]
    }
}

/// A chromosome region, 1-based and inclusive.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    pub chr: i32,
    pub start: i32,
    pub end: i32,
}

impl Region {
    pub fn new(chr: i32, start: i32, end: i32) -> Self {
        Region { chr, start, end }
    }
    pub fn contains_position(&self, pos: i32) -> bool {
        pos >= self.start && pos <= self.end
    }
    pub fn contains(&self, chr: i32, pos: i32) -> bool {
        self.chr == chr && self.contains_position(pos)
    }
    pub fn overlaps(&self, chr: i32, start: i32, end: i32) -> bool {
        self.chr == chr && positions_overlap(self.start, self.end, start, end)
    }
    pub fn overlaps_region(&self, other: &Region) -> bool {
        self.overlaps(other.chr, other.start, other.end)
    }
    pub fn display(&self, contigs: &Contigs) -> String {
        format!("{}:{}-{}", contigs.name(self.chr), self.start, self.end)
    }
}

#[inline]
pub fn positions_overlap(start1: i32, end1: i32, start2: i32, end2: i32) -> bool {
    !(start1 > end2 || end1 < start2)
}

#[inline]
pub fn position_within(pos: i32, start: i32, end: i32) -> bool {
    pos >= start && pos <= end
}

/// SupplementaryReadData: the first alignment of the SA attribute
#[derive(Clone, Debug, PartialEq)]
pub struct SuppData {
    pub chr_name: String,
    pub chr: i32,
    pub position: i32,
    pub strand: u8,
    pub cigar: String,
    pub map_quality: i32,
    pub nm: i32,
}

impl SuppData {
    /// SupplementaryReadData.extractAlignment
    pub fn extract(sa: Option<&[u8]>, contigs: &Contigs) -> Option<SuppData> {
        let sa = std::str::from_utf8(sa?).ok()?;
        if sa.is_empty() {
            return None;
        }
        if sa.contains(',') {
            let first = sa.split(';').next().unwrap_or("");
            Self::from_alignment(first, ',', contigs)
        } else {
            // the delimiter of the Java unit tests
            Self::from_alignment(sa, ';', contigs)
        }
    }

    fn from_alignment(alignment: &str, delim: char, contigs: &Contigs) -> Option<SuppData> {
        // String.split drops trailing empty items
        let mut items: Vec<&str> = alignment.split(delim).collect();
        while items.last().is_some_and(|s| s.is_empty()) {
            items.pop();
        }
        if items.len() != 6 {
            return None;
        }
        Some(SuppData {
            chr_name: items[0].to_string(),
            chr: contigs.id(items[0]),
            position: items[1].parse().ok()?,
            strand: *items[2].as_bytes().first()?,
            cigar: items[3].to_string(),
            map_quality: items[4].parse().ok()?,
            nm: items[5].parse().ok()?,
        })
    }

    pub fn orient(&self) -> Orientation {
        if self.strand == b'+' { Orientation::Forward } else { Orientation::Reverse }
    }

    pub fn as_delim_str(&self) -> String {
        format!("{};{};{};{};{};{}", self.chr_name, self.position, self.strand as char, self.cigar, self.map_quality, self.nm)
    }
}

/// IndelCoords
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndelCoords {
    pub pos_start: i32,
    pub pos_end: i32,
    pub length: i32,
}

impl IndelCoords {
    pub fn is_insert(&self) -> bool {
        self.pos_end == self.pos_start + 1
    }
}

/// IndelCoords.findIndelCoords: the longest insert or delete of the cigar
pub fn find_indel_coords(read_start: i32, cigar: &[u32], min_indel_length: i32) -> Option<IndelCoords> {
    let mut max_indel_length = 0;
    let mut indel_start = read_start - 1;
    let mut indel_end = 0;
    let mut ref_position = read_start;
    for &c in cigar {
        let (o, len) = (op(c), oplen(c));
        if (o == OP_I || o == OP_D) && len > max_indel_length {
            max_indel_length = len;
            indel_start = ref_position - 1;
            indel_end = if o == OP_D { indel_start + len + 1 } else { indel_start + 1 };
        }
        if consumes_ref(o) {
            ref_position += len;
        }
    }
    if indel_end <= indel_start || max_indel_length < min_indel_length {
        return None;
    }
    Some(IndelCoords { pos_start: indel_start, pos_end: indel_end, length: max_indel_length })
}

pub const MIN_INDEL_SUPPORT_LENGTH: i32 = 3;

/// PrepRead: a BAM record with the values prep works from.
pub struct PrepRead {
    pub rec: bam::Record,
    /// index of the sample (input BAM) the record came from; Java keeps it in an attribute of the record
    pub sample: u8,
    pub chr: i32,
    pub start: i32,
    pub end: i32,
    pub unclipped_start: i32,
    pub unclipped_end: i32,
    pub mate_chr: i32,
    pub mate_pos_start: i32,
    pub fragment_insert_size: i32,
    pub supp: Option<Box<SuppData>>,
    /// hasSuppAlignment: a supplementary alignment on a human chromosome
    pub supp_is_human: bool,
    pub aligned_base_length: i32,
    pub sc_left: i32,
    pub sc_right: i32,
    pub max_indel_length: i32,
    pub has_line_tail: bool,
    pub filters: i32,
    pub read_type: ReadType,
    pub written: bool,
    /// the read bases as letters
    pub bases: Vec<u8>,
}

/// SamRecordUtils.inferredInsertSize, absolute
pub fn inferred_insert_size_abs(rec: &bam::Record) -> i32 {
    if has(rec, FLAG_PAIRED) { (rec.insert_size() as i32).wrapping_abs() } else { (end(rec) - start(rec) + 1).wrapping_abs() }
}

/// SamRecordUtils.mateUnmapped
#[inline]
pub fn mate_unmapped(rec: &bam::Record) -> bool {
    rec.flags() & (FLAG_PAIRED | FLAG_MATE_UNMAPPED) == (FLAG_PAIRED | FLAG_MATE_UNMAPPED)
}

impl PrepRead {
    /// `tid_map` gives the contig index of each reference of the record's BAM.
    pub fn new(rec: bam::Record, sample: u8, tid_map: &[i32], contigs: &Contigs) -> PrepRead {
        let cigar = rec.raw_cigar();
        let (mut aligned, mut sc_left, mut sc_right, mut max_indel) = (0, 0, 0, 0);
        for (i, &c) in cigar.iter().enumerate() {
            let len = oplen(c);
            match op(c) {
                OP_M => aligned += len,
                OP_S => {
                    if i == 0 {
                        sc_left = len
                    } else {
                        sc_right = len
                    }
                }
                OP_D | OP_I => max_indel = max_indel.max(len),
                _ => {}
            }
        }
        let lookup = |tid: i32| if tid >= 0 && (tid as usize) < tid_map.len() { tid_map[tid as usize] } else { UNMAPPED_CHR };
        let (chr, s, e, us, ue);
        if !has(&rec, FLAG_UNMAPPED) {
            chr = lookup(rec.tid());
            s = start(&rec);
            e = end(&rec);
            us = s - sc_left;
            ue = e + sc_right;
        } else {
            (chr, s, e, us, ue) = (UNMAPPED_CHR, 0, 0, 0, 0);
        }
        let (mate_chr, mate_pos_start) =
            if !mate_unmapped(&rec) && mate_start(&rec) > 0 { (lookup(rec.mtid()), mate_start(&rec)) } else { (UNMAPPED_CHR, 0) };
        let fragment_insert_size = inferred_insert_size_abs(&rec);
        let supp = SuppData::extract(aux_str(&rec, b"SA"), contigs).map(Box::new);
        let bases = seq_bytes(&rec);
        let supp_is_human = supp.as_ref().is_some_and(|s| is_human_chromosome(&s.chr_name));
        PrepRead {
            rec,
            sample,
            chr,
            start: s,
            end: e,
            unclipped_start: us,
            unclipped_end: ue,
            mate_chr,
            mate_pos_start,
            fragment_insert_size,
            supp,
            supp_is_human,
            aligned_base_length: aligned,
            sc_left,
            sc_right,
            max_indel_length: max_indel,
            has_line_tail: false,
            filters: 0,
            read_type: ReadType::NoSupport,
            written: false,
            bases,
        }
    }

    /// A copy without the record, marked as written: stands in for a read that another holder has written.
    pub fn written_stand_in(&self) -> PrepRead {
        PrepRead {
            rec: bam::Record::new(),
            sample: self.sample,
            chr: self.chr,
            start: self.start,
            end: self.end,
            unclipped_start: self.unclipped_start,
            unclipped_end: self.unclipped_end,
            mate_chr: self.mate_chr,
            mate_pos_start: self.mate_pos_start,
            fragment_insert_size: self.fragment_insert_size,
            supp: None,
            supp_is_human: self.supp_is_human,
            aligned_base_length: self.aligned_base_length,
            sc_left: self.sc_left,
            sc_right: self.sc_right,
            max_indel_length: self.max_indel_length,
            has_line_tail: self.has_line_tail,
            filters: self.filters,
            read_type: self.read_type,
            written: true,
            bases: Vec::new(),
        }
    }

    /// Takes the read out of its slot, leaving an empty one.
    pub fn take(&mut self) -> PrepRead {
        let mut empty = self.written_stand_in();
        empty.written = self.written;
        std::mem::replace(self, empty)
    }

    pub fn id(&self) -> &[u8] {
        self.rec.qname()
    }
    pub fn is_left_clipped(&self) -> bool {
        self.sc_left > 0
    }
    pub fn is_right_clipped(&self) -> bool {
        self.sc_right > 0
    }
    pub fn orientation(&self) -> Orientation {
        if !self.is_read_reversed() { Orientation::Forward } else { Orientation::Reverse }
    }
    pub fn mate_orientation(&self) -> Orientation {
        if !has(&self.rec, FLAG_MATE_REVERSE) { Orientation::Forward } else { Orientation::Reverse }
    }
    pub fn flags(&self) -> u16 {
        self.rec.flags()
    }
    pub fn cigar(&self) -> &[u32] {
        self.rec.raw_cigar()
    }
    pub fn is_read_reversed(&self) -> bool {
        has(&self.rec, FLAG_REVERSE)
    }
    /// SamRecordUtils.firstInPair: an unpaired read counts as first
    pub fn is_first_of_pair(&self) -> bool {
        !has(&self.rec, FLAG_PAIRED) || has(&self.rec, FLAG_FIRST)
    }
    pub fn is_paired(&self) -> bool {
        has(&self.rec, FLAG_PAIRED)
    }
    pub fn is_supplementary(&self) -> bool {
        has(&self.rec, FLAG_SUPPLEMENTARY)
    }
    pub fn is_unmapped(&self) -> bool {
        has(&self.rec, FLAG_UNMAPPED)
    }
    pub fn has_mate(&self) -> bool {
        self.mate_pos_start > 0
    }
    pub fn is_mate_unmapped(&self) -> bool {
        has(&self.rec, FLAG_MATE_UNMAPPED)
    }
    pub fn is_mate_mapped(&self) -> bool {
        self.is_paired() && !self.is_mate_unmapped()
    }
    pub fn has_supp_alignment(&self) -> bool {
        self.supp_is_human
    }
    pub fn base_qualities(&self) -> &[u8] {
        self.rec.qual()
    }
    pub fn add_filter(&mut self, f: i32) {
        self.filters |= f;
    }
    pub fn has_filter(&self, f: i32) -> bool {
        self.filters & f != 0
    }
    pub fn unfiltered(&self) -> bool {
        self.filters == 0
    }
    pub fn set_read_type(&mut self, t: ReadType) {
        self.read_type = t;
    }
    /// setReadType(type, true): keeps the higher ranked type
    pub fn set_read_type_ranked(&mut self, t: ReadType) {
        if t.rank() > self.read_type.rank() {
            self.read_type = t;
        }
    }
    pub fn map_quality(&self) -> i32 {
        self.rec.mapq() as i32
    }
    pub fn indel_coords(&self) -> Option<IndelCoords> {
        if self.max_indel_length >= MIN_INDEL_SUPPORT_LENGTH {
            find_indel_coords(self.start, self.cigar(), MIN_INDEL_SUPPORT_LENGTH)
        } else {
            None
        }
    }
}

/// A chromosome partition as the spanning cache names it: "<chromosome>_<index>". Partitions start at
/// 1 + index * size; Java takes the index as position / size (finding P4).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct PartKey {
    pub chr: i32,
    pub part: i32,
}

impl PartKey {
    /// SpanningReadCache.formChromosomePartition
    pub fn form(chr: i32, position: i32, partition_size: i32) -> PartKey {
        if crate::common::compat::java("P4") {
            // Java: the last base of a partition is named as the next partition
            PartKey { chr, part: position / partition_size }
        } else {
            PartKey { chr, part: (position - 1).max(0) / partition_size }
        }
    }
    pub fn name(&self, contigs: &Contigs) -> String {
        format!("{}_{}", contigs.name(self.chr), self.part)
    }
}

/// ReadGroup: the reads of one read id within a partition.
pub struct ReadGroup {
    pub id: Box<[u8]>,
    pub reads: Vec<ReadId>,
    pub status: ReadGroupStatus,
    /// a Java HashSet; `remote_partitions_ordered` gives its iteration order
    pub remote_partitions: Vec<PartKey>,
    pub expected_read_count: i32,
    pub junction_positions: Vec<(i32, Orientation)>,
    pub has_remote_junction_reads: bool,
    /// member of the tracker's set of remote candidate groups
    pub remote_candidate: bool,
    /// member of the tracker's list of expected groups
    pub expected: bool,
}

impl ReadGroup {
    pub fn new(id: Box<[u8]>, read: ReadId) -> ReadGroup {
        ReadGroup {
            id,
            reads: vec![read],
            status: ReadGroupStatus::Unset,
            remote_partitions: Vec::new(),
            expected_read_count: 0,
            junction_positions: Vec::new(),
            has_remote_junction_reads: false,
            remote_candidate: false,
            expected: false,
        }
    }
    pub fn size(&self) -> usize {
        self.reads.len()
    }
    pub fn is_complete(&self) -> bool {
        self.status == ReadGroupStatus::Complete
    }
    pub fn spans_partitions(&self) -> bool {
        !self.remote_partitions.is_empty()
    }
    pub fn has_junction_position(&self, position: i32, orient: Orientation) -> bool {
        self.junction_positions.iter().any(|&(p, o)| p == position && o == orient)
    }
    pub fn has_junction_positions(&self) -> bool {
        !self.junction_positions.is_empty()
    }
    pub fn add_junction_position(&mut self, position: i32, orient: Orientation) {
        if !self.has_junction_position(position, orient) {
            self.junction_positions.push((position, orient));
        }
    }
    pub fn junction_positions_str(&self) -> String {
        self.junction_positions.iter().map(|(p, o)| format!("{}:{}", p, o.as_byte())).collect::<Vec<_>>().join(";")
    }
    pub fn add_remote_partition(&mut self, key: PartKey) {
        if !self.remote_partitions.contains(&key) {
            self.remote_partitions.push(key);
        }
    }
    /// The remote partitions in the iteration order of the Java HashSet of their names.
    pub fn remote_partitions_ordered(&self, contigs: &Contigs) -> Vec<PartKey> {
        if self.remote_partitions.len() <= 1 {
            return self.remote_partitions.clone();
        }
        let hashes: Vec<i32> = self.remote_partitions.iter().map(|k| string_hash(k.name(contigs).as_bytes())).collect();
        set_order(&hashes).into_iter().map(|i| self.remote_partitions[i]).collect()
    }
}

/// The group functions that look at the reads.
pub fn group_is_simple_complete(g: &ReadGroup, reads: &[PrepRead]) -> bool {
    if g.reads.len() == 1 {
        let r = &reads[g.reads[0] as usize];
        return !r.has_mate() && !r.has_supp_alignment();
    }
    g.reads.len() == 2 && g.reads.iter().all(|&i| !reads[i as usize].has_supp_alignment() && !reads[i as usize].is_supplementary())
}

pub fn group_all_no_support(g: &ReadGroup, reads: &[PrepRead]) -> bool {
    g.reads.iter().all(|&i| reads[i as usize].read_type == ReadType::NoSupport)
}

pub fn group_has_read_type(g: &ReadGroup, reads: &[PrepRead], t: ReadType) -> bool {
    g.reads.iter().any(|&i| reads[i as usize].read_type == t)
}

/// conditionalOnRemoteReads: a group of candidate, unsupporting or unmapped reads, or of supporting reads with
/// low-quality soft clips, is only kept when its remote reads support a junction
pub fn group_conditional_on_remote_reads(g: &ReadGroup, reads: &[PrepRead]) -> bool {
    g.reads.iter().all(|&i| {
        let r = &reads[i as usize];
        r.read_type == ReadType::CandidateSupport
            || r.read_type == ReadType::NoSupport
            || r.is_unmapped()
            || (r.read_type == ReadType::Support && r.has_filter(F_SOFT_CLIP_LOW_BASE_QUAL))
    })
}

fn group_has_read_mate(g: &ReadGroup, reads: &[PrepRead], read: ReadId) -> bool {
    let r = &reads[read as usize];
    if !r.has_mate() {
        return false;
    }
    for &oi in &g.reads {
        if oi == read {
            continue;
        }
        let o = &reads[oi as usize];
        if r.is_first_of_pair() == o.is_first_of_pair() {
            continue;
        }
        if r.is_supplementary() != o.is_supplementary() {
            continue;
        }
        if r.mate_chr == o.chr && r.mate_pos_start == o.start {
            return true;
        }
    }
    false
}

/// ReadGroup.setGroupState()
pub fn group_set_state(g: &mut ReadGroup, reads: &[PrepRead]) {
    if group_conditional_on_remote_reads(g, reads) {
        if g.status != ReadGroupStatus::Expected {
            g.status = if g.reads.len() == 2 { ReadGroupStatus::Paired } else { ReadGroupStatus::Supplementary };
        }
        g.expected_read_count = g.reads.len() as i32;
        return;
    }
    let (mut first_has_supp, mut second_has_supp, mut non_supp_paired) = (false, false, false);
    for &i in &g.reads {
        let r = &reads[i as usize];
        if !r.is_supplementary() && !non_supp_paired {
            non_supp_paired = group_has_read_mate(g, reads, i);
        }
        if r.has_supp_alignment() || r.is_supplementary() {
            if r.is_first_of_pair() {
                first_has_supp = true;
            } else {
                second_has_supp = true;
            }
        }
    }
    let supp_count = first_has_supp as i32 + second_has_supp as i32;
    let main_read_count = if reads[g.reads[0] as usize].has_mate() { 2 } else { 1 };
    g.expected_read_count = main_read_count + supp_count;
    g.status = if g.reads.len() as i32 >= g.expected_read_count {
        ReadGroupStatus::Complete
    } else if non_supp_paired {
        ReadGroupStatus::Paired
    } else {
        ReadGroupStatus::Incomplete
    };
}

/// ReadGroup.setPartitionCount: the partitions of mates and supplementaries outside the region
pub fn group_set_partition_count(g: &mut ReadGroup, reads: &[PrepRead], region: &Region, partition_size: i32, contigs: &Contigs) {
    for k in 0..g.reads.len() {
        let r = &reads[g.reads[k] as usize];
        if r.is_unmapped() {
            continue;
        }
        if r.has_supp_alignment() {
            let s = r.supp.as_ref().unwrap();
            if !region.contains(s.chr, s.position) {
                g.add_remote_partition(PartKey::form(s.chr, s.position, partition_size));
            }
        }
        if r.has_mate() && contigs.is_human(r.mate_chr) && !region.contains(r.mate_chr, r.mate_pos_start) {
            g.add_remote_partition(PartKey::form(r.mate_chr, r.mate_pos_start, partition_size));
        }
    }
}

/// RemoteJunction
#[derive(Clone, Debug)]
pub struct RemoteJunction {
    pub chr: i32,
    pub position: i32,
    pub orient: Orientation,
    pub fragments: i32,
}

/// JunctionData
pub struct JunctionData {
    pub position: i32,
    pub orient: Orientation,
    pub junction_groups: Vec<GroupId>,
    pub supporting_groups: Vec<GroupId>,
    pub exact_support_groups: Vec<GroupId>,
    pub remote_junctions: Vec<RemoteJunction>,
    /// the reads by type: junction, support, exact support
    pub junction_reads: Vec<ReadId>,
    pub support_reads: Vec<ReadId>,
    pub exact_support_reads: Vec<ReadId>,
    pub top_junction_read: Option<ReadId>,
    pub internal_indel: bool,
    pub discordant_group: bool,
    pub hotspot: bool,
    pub depth: i32,
    pub linked_indel: Option<JuncId>,
    pub max_valid_soft_clip_length: i32,
    pub saga_match: Option<String>,
    /// for a junction at the remote region of a discordant group that lies outside the partition: its
    /// chromosome
    pub outside_chr: Option<i32>,
}

impl JunctionData {
    pub fn new(position: i32, orient: Orientation, read: Option<ReadId>) -> JunctionData {
        JunctionData {
            position,
            orient,
            junction_groups: Vec::new(),
            supporting_groups: Vec::new(),
            exact_support_groups: Vec::new(),
            remote_junctions: Vec::new(),
            junction_reads: Vec::new(),
            support_reads: Vec::new(),
            exact_support_reads: Vec::new(),
            top_junction_read: read,
            internal_indel: false,
            discordant_group: false,
            hotspot: false,
            depth: 0,
            linked_indel: None,
            max_valid_soft_clip_length: 0,
            saga_match: None,
            outside_chr: None,
        }
    }
    pub fn is_forward(&self) -> bool {
        self.orient.is_forward()
    }
    pub fn is_reverse(&self) -> bool {
        self.orient.is_reverse()
    }
    pub fn junction_fragment_count(&self) -> i32 {
        self.junction_groups.len() as i32
    }
    pub fn supporting_fragment_count(&self) -> i32 {
        self.supporting_groups.len() as i32
    }
    pub fn exact_support_fragment_count(&self) -> i32 {
        self.exact_support_groups.len() as i32
    }
    /// addReadType: only junction, support and exact support reads are kept
    pub fn add_read_type(&mut self, read: ReadId, t: ReadType) {
        match t {
            ReadType::Junction => self.junction_reads.push(read),
            ReadType::Support => self.support_reads.push(read),
            ReadType::ExactSupport => self.exact_support_reads.push(read),
            _ => {}
        }
    }
    pub fn add_remote_junction(&mut self, rj: RemoteJunction) {
        if let Some(m) = self.remote_junctions.iter_mut().find(|x| x.chr == rj.chr && x.position == rj.position && x.orient == rj.orient) {
            m.fragments += 1;
            return;
        }
        self.remote_junctions.push(rj);
    }
}
