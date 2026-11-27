//! The telbam: TelbamApp.processBam, i.e. BamProcessor, PartitionReader and BamRecordWriter.
//!
//! Every record of the BAM is tested for telomeric content (TTAGGGTTAGGG or CCCTAACCCTAA). A telomeric
//! record opens a read group under its name; later records of the same name join it, until the group has
//! both reads of the pair and every supplementary alignment their SA attributes list; then it is written.
//! Groups still incomplete after the pass look for their missing reads at the mate and SA positions,
//! repeatedly, while that finds new reads. The incomplete groups left are written at the end.
//!
//! Java reads the BAM in partitions: the unmapped reads without a position first, then 10 Mb regions of
//! each contig in header order, with `threads - 1` reader threads (at least one) feeding one writer.
//! With one reader thread (`-threads` 1 or 2) the records reach the writer in exactly this order: the
//! unplaced unmapped records, then every placed record in file order. The port reads them in that order
//! with one pass over the file (its BGZF decompression on `threads - 1` helper threads), so it gives
//! Java's single-reader result for every thread count (TE6). The missing-read passes read the unplaced
//! records and the consolidated regions in Java's order (chromosome rank, then position).
//!
//! The telbam is written as htsjdk's writer does with `presorted = false`: with the input header, and
//! for a coordinate-sorted header the records sorted with SAMRecordCoordinateComparator, stable over the
//! order they were added (group by group, as groups complete, then the incomplete groups in Java HashMap
//! order).

use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use crate::read::{describe, Contigs, Read, ReadGroup};
use crate::telo;

/// TealConstants.DEFAULT_PARTITION_SIZE
const PARTITION_SIZE: i64 = 10_000_000;
/// BamProcessor.REGION_CONSOLIDATION_DISTANCE
const REGION_CONSOLIDATION_DISTANCE: i32 = 100_000;
/// The file buffer of the reader: few large reads keep a disk shared with other jobs streaming
const READ_BUFFER: i32 = 4 << 20;

/// FxHash (rustc's): fast hashing of read names
#[derive(Default, Clone, Copy)]
pub struct FxHasher {
    hash: u64,
}

const FX_SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

impl FxHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.hash = (self.hash.rotate_left(5) ^ w).wrapping_mul(FX_SEED);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut b = bytes;
        while b.len() >= 8 {
            self.add(u64::from_le_bytes(b[..8].try_into().unwrap()));
            b = &b[8..];
        }
        if b.len() >= 4 {
            self.add(u32::from_le_bytes(b[..4].try_into().unwrap()) as u64);
            b = &b[4..];
        }
        for &x in b {
            self.add(x as u64);
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

pub type FxBuild = BuildHasherDefault<FxHasher>;

/// One of TelbamApp's parameters
pub struct TelbamParams<'a> {
    pub bam_file: &'a str,
    pub ref_genome: Option<&'a str>,
    pub telbam_file: &'a str,
    pub threads: i32,
}

struct Entry {
    /// insertion order (Java's order of entries within a hash bucket)
    seq: u64,
    group: ReadGroup,
}

/// BamRecordWriter
struct RecordWriter {
    contigs: Contigs,
    groups: HashMap<Box<[u8]>, Entry, FxBuild>,
    /// how many group names fall in each slot of a cheap key of the name: most records are not in a group,
    /// and a zero here says so without hashing the whole name
    name_filter: Vec<u32>,
    next_seq: u64,
    /// the most groups held at once: it sets the size of Java's hash table
    max_groups: usize,
    processing_mate_regions: bool,
    num_accepted: u64,
    num_completed: u64,
    /// the records written, in the order they were added
    out: Vec<Read>,
    /// CRAM: the attributes in the order htsjdk's CRAM reader gives them
    htsjdk_tag_order: bool,
}

impl RecordWriter {
    /// PartitionReader: consensus reads are skipped; a record goes to the writer when it is telomeric or
    /// its name belongs to an incomplete group. Returns whether it went to the writer.
    #[inline]
    fn offer(&mut self, rec: &bam::Record) -> Result<bool> {
        let seq_len = rec.seq_len();
        let telo = telo::has_telomeric_content_packed(rec.seq().encoded, seq_len);
        if !telo && !self.may_have_group(rec.qname()) {
            return Ok(false);
        }
        if rec.aux(b"CR").is_ok() {
            return Ok(false);
        }
        self.process_read_record(rec, telo)?;
        Ok(true)
    }

    /// `offer` for a raw record of the main pass: the record is only built when it goes to the writer
    #[inline]
    fn offer_raw(&mut self, raw: &[u8], v: &crate::bgzf::RawView, rec: &mut bam::Record) -> Result<bool> {
        let telo = telo::has_telomeric_content_packed(v.seq, v.seq_len);
        if !telo && !self.may_have_group(v.name) {
            return Ok(false);
        }
        crate::bgzf::record_from_raw(raw, rec);
        if rec.aux(b"CR").is_ok() {
            return Ok(false);
        }
        self.process_read_record(rec, telo)?;
        Ok(true)
    }

    /// processReadRecord
    fn process_read_record(&mut self, rec: &bam::Record, has_telomere_content: bool) -> Result<()> {
        let name = rec.qname();
        if !self.groups.contains_key(name) {
            if self.processing_mate_regions {
                // no new group in the mate regions: it was written before
                return Ok(());
            }
            if !has_telomere_content {
                return Ok(());
            }
            let seq = self.next_seq;
            self.next_seq += 1;
            self.groups.insert(Box::from(name), Entry { seq, group: ReadGroup::new(name) });
            let slot = name_slot(name);
            self.name_filter[slot] += 1;
            self.max_groups = self.max_groups.max(self.groups.len());
        }
        let mut read = Read::from_record(rec, &self.contigs, false, true);
        if self.htsjdk_tag_order {
            sort_tags_as_htsjdk(read.raw.as_mut().unwrap());
        }
        let entry = self.groups.get_mut(name).unwrap();
        if !entry.group.contains(&read) && entry.group.accept_read(read) {
            self.num_accepted += 1;
        }
        debug_assert!(entry.group.invariant() || entry.group.reads.len() > 2);
        if entry.group.is_complete()? {
            let e = self.groups.remove(name).unwrap();
            self.name_filter[name_slot(name)] -= 1;
            self.write_group(e.group);
            self.num_completed += 1;
        }
        Ok(())
    }

    fn write_group(&mut self, group: ReadGroup) {
        self.out.extend(group.reads);
        self.out.extend(group.supps);
    }

    /// The incomplete groups in Java HashMap order: by bucket of the name's String.hashCode in the table
    /// Java has grown to, then in insertion order
    fn groups_in_java_order(&self) -> Vec<&Entry> {
        let mut capacity = 16usize;
        while self.max_groups > capacity * 3 / 4 {
            capacity *= 2;
        }
        let mask = (capacity - 1) as u32;
        let mut entries: Vec<&Entry> = self.groups.values().collect();
        entries.sort_by_key(|e| {
            let h = java_bytes_hash(&e.group.name) as u32;
            ((h ^ (h >> 16)) & mask, e.seq)
        });
        entries
    }
}

const NAME_FILTER_BITS: u32 = 20;

/// The slot of a name in the name filter: its last eight bytes (they hold the varying part of
/// sequencer read names) and its length
#[inline]
fn name_slot(name: &[u8]) -> usize {
    let n = name.len();
    let w = if n >= 8 {
        u64::from_le_bytes(name[n - 8..].try_into().unwrap())
    } else {
        let mut b = [0u8; 8];
        b[..n].copy_from_slice(name);
        u64::from_le_bytes(b)
    };
    ((w ^ (n as u64)).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - NAME_FILTER_BITS)) as usize
}

impl RecordWriter {
    /// Whether a group may have this name (false: certainly not)
    #[inline]
    fn may_have_group(&self, name: &[u8]) -> bool {
        self.name_filter[name_slot(name)] != 0 && self.groups.contains_key(name)
    }
}

/// htsjdk's CRAM reader sets each attribute of a record with setAttribute, which inserts it into the list
/// sorted by its binary code (SAMTag.makeBinaryTag: second character << 8 | first), and writes integers in
/// the smallest type: the attributes of a record read from CRAM come out that way (from BAM their bytes
/// are kept). Checked on the HG008 SBX CRAMs: the telbams' uncompressed bytes equal Java's.
pub fn sort_tags_as_htsjdk(rec: &mut bam::Record) {
    // and a record without a reference gets bin 0 (htsjdk writes 0 for those whatever their stored bin;
    // `write_bam` does the same for records from a BAM)
    if rec.tid() < 0 {
        rec.set_bin(0);
    }
    let core = rec.inner().core;
    let l_data = rec.inner().l_data as usize;
    let data = unsafe { std::slice::from_raw_parts(rec.inner().data, l_data) }.to_vec();
    let aux_off = core.l_qname as usize + 4 * core.n_cigar as usize + (core.l_qseq as usize).div_ceil(2) + core.l_qseq as usize;
    if aux_off >= data.len() {
        return;
    }
    let mut fields: Vec<(u16, Vec<u8>)> = Vec::new();
    let mut p = aux_off;
    while p + 3 <= data.len() {
        let size = match data[p + 2] {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'Z' | b'H' => match data[p + 3..].iter().position(|&b| b == 0) {
                Some(n) => n + 1,
                None => return,
            },
            b'B' => {
                if p + 8 > data.len() {
                    return;
                }
                let n = u32::from_le_bytes(data[p + 4..p + 8].try_into().unwrap()) as usize;
                let each = match data[p + 3] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    _ => 4,
                };
                5 + n * each
            }
            _ => return,
        };
        let end = p + 3 + size;
        if end > data.len() {
            return;
        }
        fields.push((((data[p + 1] as u16) << 8) | data[p] as u16, retype_integer(&data[p..end])));
        p = end;
    }
    fields.sort_by_key(|f| f.0);
    let mut out = data[..aux_off].to_vec();
    for (_, f) in &fields {
        out.extend_from_slice(f);
    }
    rec.set_data(&out);
}

/// htsjdk writes an integer attribute it holds as a number with the smallest type that fits
/// (BinaryTagCodec.getIntegerType): `tag type value` of a c/C/s/S/i/I field in that type; others unchanged
fn retype_integer(field: &[u8]) -> Vec<u8> {
    let v: i64 = match field[2] {
        b'c' => field[3] as i8 as i64,
        b'C' => field[3] as i64,
        b's' => i16::from_le_bytes([field[3], field[4]]) as i64,
        b'S' => u16::from_le_bytes([field[3], field[4]]) as i64,
        b'i' => i32::from_le_bytes(field[3..7].try_into().unwrap()) as i64,
        b'I' => u32::from_le_bytes(field[3..7].try_into().unwrap()) as i64,
        _ => return field.to_vec(),
    };
    let mut out = field[..2].to_vec();
    if v > i32::MAX as i64 {
        out.push(b'I');
        out.extend_from_slice(&(v as u32).to_le_bytes());
    } else if v > u16::MAX as i64 {
        out.push(b'i');
        out.extend_from_slice(&(v as i32).to_le_bytes());
    } else if v > i16::MAX as i64 {
        out.push(b'S');
        out.extend_from_slice(&(v as u16).to_le_bytes());
    } else if v > u8::MAX as i64 {
        out.push(b's');
        out.extend_from_slice(&(v as i16).to_le_bytes());
    } else if v > i8::MAX as i64 {
        out.push(b'C');
        out.push(v as u8);
    } else if v >= i8::MIN as i64 {
        out.push(b'c');
        out.push(v as i8 as u8);
    } else if v >= i16::MIN as i64 {
        out.push(b's');
        out.extend_from_slice(&(v as i16).to_le_bytes());
    } else {
        out.push(b'i');
        out.extend_from_slice(&(v as i32).to_le_bytes());
    }
    out
}

/// String.hashCode of an ASCII (Latin-1) name
fn java_bytes_hash(b: &[u8]) -> i32 {
    let mut h: i32 = 0;
    for &c in b {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

pub fn open_reader(path: &str, ref_genome: Option<&str>, threads: i32) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening {path} (an indexed BAM or CRAM is required)"))?;
    unsafe {
        use rust_htslib::htslib;
        htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_HTS_OPT_BLOCK_SIZE, READ_BUFFER as std::os::raw::c_int);
    }
    if let Some(r) = ref_genome.filter(|r| !r.is_empty()) {
        if is_cram(path) {
            reader.set_reference(r).with_context(|| format!("setting CRAM reference {r}"))?;
        }
    }
    if is_cram(path) {
        // htsjdk's CRAM reader does not make MD and NM up; htslib does unless told not to
        unsafe {
            use rust_htslib::htslib;
            htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 0 as std::os::raw::c_int);
        }
    }
    if threads > 1 {
        reader.set_threads((threads - 1) as usize).context("setting reader threads")?;
    }
    Ok(reader)
}

fn is_cram(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".cram")
}

/// A missing-read region (ChrBaseRegion)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
}

/// BamProcessor.getMissingReadRegions without the unmapped partition: the missing positions sorted as
/// ChrBaseRegion.compareTo sorts them and merged when within 100 kb
pub fn consolidate_regions(mut regions: Vec<Region>) -> Vec<Region> {
    regions.sort_by(|a, b| {
        if a.chromosome == b.chromosome { a.start.cmp(&b.start) } else { hmf_common::regions::compare_chromosomes(&a.chromosome, &b.chromosome) }
    });
    let mut compact: Vec<Region> = Vec::new();
    for br in regions {
        if let Some(last) = compact.last_mut() {
            if last.chromosome == br.chromosome && last.end + REGION_CONSOLIDATION_DISTANCE >= br.start {
                if br.end > last.end {
                    last.end = br.end;
                }
                continue;
            }
        }
        compact.push(br);
    }
    compact
}

/// htsjdk's SAMRecordCoordinateComparator
pub fn coordinate_compare(a: &Read, b: &Read) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    // fileOrderCompare: unplaced records last
    let c = match (a.tid == -1, b.tid == -1) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        _ => a.tid.cmp(&b.tid).then(a.start.cmp(&b.start)),
    };
    if c != Ordering::Equal {
        return c;
    }
    if a.neg() != b.neg() {
        return if a.neg() { Ordering::Greater } else { Ordering::Less };
    }
    a.name
        .cmp(&b.name)
        .then(a.flags.cmp(&b.flags))
        .then(a.mapq.cmp(&b.mapq))
        .then(a.mtid.cmp(&b.mtid))
        .then(a.mate_start.cmp(&b.mate_start))
        .then(a.isize.cmp(&b.isize))
}

/// The partition (Java's 10 Mb region) a position falls in, for the logs
fn partition_label(contig: &str, index: i64, length: i64) -> String {
    let start = (index * PARTITION_SIZE).max(1);
    let end = ((index + 1) * PARTITION_SIZE - 1).min(length);
    format!("{contig}:{}-{}", group_thousands(start), group_thousands(end))
}

fn group_thousands(v: i64) -> String {
    let s = v.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// TelbamApp.processBam / BamProcessor.processBam
pub fn process_bam(params: &TelbamParams) -> Result<()> {
    crate::info!("processing bam: {}", params.bam_file);
    let started = std::time::Instant::now();
    let mut reader = open_reader(params.bam_file, params.ref_genome, params.threads)?;
    // htsjdk's writer creates the file at the start (an empty file is left when the run fails)
    std::fs::File::create(params.telbam_file).with_context(|| format!("creating {}", params.telbam_file))?;
    let header_view = reader.header().clone();
    let contigs = Contigs::from_header(&header_view);
    // htsjdk writes the input header with its own version, VN:1.6, whatever the input's
    let text = hmf_common::sam_header::htsjdk_header_version(&String::from_utf8_lossy(header_view.as_bytes()));
    let header = bam::Header::from_template(&bam::HeaderView::from_bytes(text.as_bytes()));
    let coordinate_sorted = sort_order(&header_view) == Some("coordinate".to_string());

    let mut w = RecordWriter {
        contigs,
        groups: HashMap::default(),
        name_filter: vec![0u32; 1 << NAME_FILTER_BITS],
        next_seq: 0,
        max_groups: 0,
        processing_mate_regions: false,
        num_accepted: 0,
        num_completed: 0,
        out: Vec::new(),
        htsjdk_tag_order: is_cram(params.bam_file),
    };
    let mut rec = bam::Record::new();
    let mut tasks = 0usize;

    // the unplaced unmapped reads first
    let n = read_unmapped(params, &mut reader, &mut rec, &mut w, None)?;
    crate::info!("processed partition(unmapped) num reads({n})");
    tasks += 1;

    // then every placed record, in file order (Java: the 10 Mb partitions of each contig in header order)
    let lengths: Vec<i64> = (0..header_view.target_count()).map(|t| header_view.target_len(t).unwrap_or(0) as i64).collect();
    let mut cur: Option<(i32, i64)> = None;
    let mut count = 0u64;
    let log_partition = |w: &RecordWriter, tid: i32, index: i64, count: u64| {
        let name = w.contigs.name(tid);
        crate::debug!("processed partition({}) num reads({count})", partition_label(&name, index, lengths[tid as usize]));
    };
    let mut track = |w: &RecordWriter, tid: i32, pos: i64, offered: bool| {
        let index = pos.max(0) / PARTITION_SIZE;
        if cur != Some((tid, index)) {
            if let Some((t, i)) = cur {
                log_partition(w, t, i, count);
            }
            cur = Some((tid, index));
            count = 0;
        }
        if offered {
            count += 1;
        }
    };
    if is_cram(params.bam_file) {
        reader.fetch(bam::FetchDefinition::All).context("reading the CRAM from its start")?;
        while let Some(r) = reader.read(&mut rec) {
            r.with_context(|| format!("reading {}", params.bam_file))?;
            let tid = rec.tid();
            if tid < 0 {
                break;
            }
            let offered = w.offer(&rec)?;
            track(&w, tid, rec.pos(), offered);
        }
    } else {
        // BGZF inflated on every thread, the records parsed in place
        let mut stream = crate::bgzf::BgzfStream::open(params.bam_file, params.threads.max(1) as usize)?;
        stream.skip_header().with_context(|| format!("reading {}", params.bam_file))?;
        while let Some(raw) = stream.next_record().with_context(|| format!("reading {}", params.bam_file))? {
            let v = crate::bgzf::view(raw)?;
            if v.tid < 0 {
                break;
            }
            let offered = w.offer_raw(raw, &v, &mut rec)?;
            track(&w, v.tid, v.pos as i64, offered);
        }
    }
    if let Some((t, i)) = cur {
        log_partition(&w, t, i, count);
    }
    tasks += lengths.iter().map(|l| ((l + PARTITION_SIZE) / PARTITION_SIZE) as usize).sum::<usize>();
    crate::info!("{tasks} bam reader tasks finished");
    crate::info!("initial BAM file processing complete");
    crate::debug!("PerfStats: initial pass {:.1}s", started.elapsed().as_secs_f64());

    // the missing reads, until no new read is found
    let mut unplaced_cache: Option<Vec<bam::Record>> = None;
    while !w.groups.is_empty() {
        let accepted = w.num_accepted;
        let mut missing: Vec<Region> = Vec::new();
        for e in w.groups.values() {
            for (chromosome, pos) in e.group.find_missing_read_base_regions()? {
                missing.push(Region { chromosome, start: pos, end: pos });
            }
        }
        let regions = consolidate_regions(missing);
        crate::info!("processing {} missing read regions", regions.len() + 1);
        // Java sets this after submitting the tasks, which may then start without it (TE1); here it is set
        // first
        w.processing_mate_regions = true;
        match unplaced_cache.as_ref() {
            // the later passes: the unplaced records of the names that had a group at the first
            Some(cache) => {
                for r in cache {
                    w.offer(r)?;
                }
            }
            None => {
                let mut cache = Vec::new();
                read_unmapped(params, &mut reader, &mut rec, &mut w, Some(&mut cache))?;
                crate::debug!("unplaced records of incomplete groups: {}", cache.len());
                unplaced_cache = Some(cache);
            }
        }
        // the regions are read on every thread; no group opens in this pass, so only records of the names
        // open now can matter: they are kept per region and offered in Java's region order
        let open: std::collections::HashSet<Box<[u8]>, FxBuild> = w.groups.keys().cloned().collect();
        let tids: Vec<Option<i32>> = regions.iter().map(|r| w.contigs.tid(&r.chromosome)).collect();
        // TE11: a mate or SA contig the header does not have: htsjdk's query fails on reference index -1
        if let Some(i) = tids.iter().position(|t| t.is_none()) {
            if crate::compat::java("TE11") {
                bail!("java.util.concurrent.ExecutionException: java.lang.IllegalArgumentException: Invalid reference index -1 (contig {} of a missing read is not in the BAM header)", regions[i].chromosome);
            }
            crate::warn!("missing read positions on contig {} which is not in the BAM header: skipped", regions[i].chromosome);
        }
        for found in read_regions(params, &regions, &tids, &open)? {
            for r in &found {
                w.offer(r)?;
            }
        }
        crate::info!("{} bam reader tasks finished", regions.len() + 1);
        if w.num_accepted == accepted {
            break;
        }
    }
    finish(&mut w, header, coordinate_sorted, params)?;
    crate::debug!("PerfStats: telbam {} {:.1}s", params.telbam_file, started.elapsed().as_secs_f64());
    Ok(())
}

/// The unplaced unmapped reads (Java's unmapped partition): with the parallel reader from where the BAI says
/// they start, or with htslib's iterator (CRAM, no BAI). With `cache`, the records of names that have a group
/// now are kept: no new group comes in the missing-read passes, so the later passes need only those.
fn read_unmapped(
    params: &TelbamParams,
    reader: &mut bam::IndexedReader,
    rec: &mut bam::Record,
    w: &mut RecordWriter,
    mut cache: Option<&mut Vec<bam::Record>>,
) -> Result<u64> {
    let mut n = 0u64;
    let start = if is_cram(params.bam_file) { None } else { crate::bgzf::unplaced_start(params.bam_file) };
    if let Some(voffset) = start {
        let mut stream = crate::bgzf::BgzfStream::open_at(params.bam_file, voffset, params.threads.max(1) as usize)?;
        if voffset == 0 {
            stream.skip_header()?;
        }
        while let Some(raw) = stream.next_record().with_context(|| format!("reading {}", params.bam_file))? {
            let v = crate::bgzf::view(raw)?;
            if v.tid >= 0 {
                continue;
            }
            if let Some(c) = cache.as_deref_mut() {
                if w.may_have_group(v.name) {
                    let mut r = bam::Record::new();
                    crate::bgzf::record_from_raw(raw, &mut r);
                    c.push(r);
                }
            }
            if w.offer_raw(raw, &v, rec)? {
                n += 1;
            }
        }
        return Ok(n);
    }
    // an index without unplaced reads gives no iterator
    if reader.fetch(bam::FetchDefinition::Unmapped).is_err() {
        return Ok(0);
    }
    while let Some(r) = reader.read(rec) {
        r.context("reading the unmapped reads")?;
        if rec.tid() >= 0 {
            continue;
        }
        if let Some(c) = cache.as_deref_mut() {
            if w.may_have_group(rec.qname()) {
                c.push(rec.clone());
            }
        }
        if w.offer(rec)? {
            n += 1;
        }
    }
    Ok(n)
}

/// The records of the regions (start within the region) whose names are in `open`, per region in file order;
/// the regions are shared out to `threads` workers, each with its own reader
fn read_regions(params: &TelbamParams, regions: &[Region], tids: &[Option<i32>], open: &std::collections::HashSet<Box<[u8]>, FxBuild>) -> Result<Vec<Vec<bam::Record>>> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let next = AtomicUsize::new(0);
    let workers = (params.threads.max(1) as usize).min(regions.len().max(1));
    let results: Vec<Result<Vec<(usize, Vec<bam::Record>)>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                let next = &next;
                scope.spawn(move || -> Result<Vec<(usize, Vec<bam::Record>)>> {
                    let mut reader = open_reader(params.bam_file, params.ref_genome, 1)?;
                    let mut rec = bam::Record::new();
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(region) = regions.get(i) else { break };
                        let mut found = Vec::new();
                        if let Some(tid) = tids[i] {
                            reader.fetch((tid, (region.start - 1).max(0) as i64, region.end as i64)).with_context(|| format!("querying {}", params.bam_file))?;
                            while let Some(r) = reader.read(&mut rec) {
                                r.with_context(|| format!("reading {}", params.bam_file))?;
                                let start = rec.pos() as i32 + 1;
                                if start < region.start || start > region.end || !open.contains(rec.qname()) {
                                    continue;
                                }
                                found.push(rec.clone());
                            }
                        }
                        out.push((i, found));
                    }
                    Ok(out)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut by_region: Vec<Vec<bam::Record>> = (0..regions.len()).map(|_| Vec::new()).collect();
    for r in results {
        for (i, found) in r? {
            by_region[i] = found;
        }
    }
    Ok(by_region)
}

/// The header's SO
fn sort_order(h: &bam::HeaderView) -> Option<String> {
    let text = String::from_utf8_lossy(h.as_bytes()).to_string();
    let hd = text.lines().find(|l| l.starts_with("@HD"))?.to_string();
    hd.split('\t').find_map(|f| f.strip_prefix("SO:").map(|s| s.to_string()))
}

/// BamRecordWriter.finish: the incomplete groups are written too; the records sorted and written
fn finish(w: &mut RecordWriter, header: bam::Header, coordinate_sorted: bool, params: &TelbamParams) -> Result<()> {
    let incomplete: Vec<ReadGroup> = w.groups_in_java_order().into_iter().map(|e| e.group.clone()).collect();
    let num_incomplete = incomplete.len();
    for rg in &incomplete {
        let reason = rg.incomplete_reason()?;
        if let Some(reason) = &reason {
            crate::warn!("{reason}");
        }
        crate::warn!("incomplete read group: id={}, is complete={}", String::from_utf8_lossy(&rg.name), reason.is_none());
        for r in &rg.reads {
            crate::warn!(
                "record({}) cigar({}) neg strand({}) suppl flag({}) suppl attr({})",
                describe(r),
                r.cigar_string(),
                r.neg(),
                r.supplementary(),
                r.sa.as_deref().unwrap_or("null")
            );
        }
        for r in &rg.supps {
            crate::warn!(
                "supplementary({}) cigar({}) neg strand({}) suppl flag({}) suppl attr({})",
                describe(r),
                r.cigar_string(),
                r.neg(),
                r.supplementary(),
                r.sa.as_deref().unwrap_or("null")
            );
        }
    }
    for rg in incomplete {
        w.write_group(rg);
    }
    w.groups.clear();
    w.name_filter.iter_mut().for_each(|c| *c = 0);
    let mut out = std::mem::take(&mut w.out);
    if coordinate_sorted {
        out.sort_by(coordinate_compare);
    } else if sort_order_is_queryname(&header) {
        bail!("queryname-sorted input is not supported");
    }
    write_bam(params.telbam_file, &header, &out, params.threads)?;
    crate::info!("wrote {} read groups, complete({}), incomplete({num_incomplete})", w.num_completed as usize + num_incomplete, w.num_completed);
    Ok(())
}

fn sort_order_is_queryname(header: &bam::Header) -> bool {
    String::from_utf8_lossy(&header.to_bytes()).lines().any(|l| l.starts_with("@HD") && l.contains("\tSO:queryname"))
}

pub fn write_bam(path: &str, header: &bam::Header, records: &[Read], threads: i32) -> Result<()> {
    let mut writer = bam::Writer::from_path(path, header, bam::Format::Bam).with_context(|| format!("creating {path}"))?;
    // htsjdk's default compression level
    writer.set_compression_level(bam::CompressionLevel::Level(5)).context("setting the compression level")?;
    if threads > 1 {
        writer.set_threads((threads - 1) as usize).context("setting writer threads")?;
    }
    for r in records {
        let rec = r.raw.as_ref().expect("a record to write");
        // htsjdk's BAMRecordCodec writes bin 0 for a record without a reference, whatever its stored bin
        // (htslib-written inputs have 4680, reg2bin(-1, 0)); placed records keep theirs
        if rec.tid() < 0 && rec.inner().core.bin != 0 {
            let mut c = rec.clone();
            c.set_bin(0);
            writer.write(&c).with_context(|| format!("writing {path}"))?;
        } else {
            writer.write(rec).with_context(|| format!("writing {path}"))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions_consolidate() {
        let r = |c: &str, p: i32| Region { chromosome: c.to_string(), start: p, end: p };
        let out = consolidate_regions(vec![r("chr2", 500), r("chr1", 300_000), r("chr10", 5), r("chr1", 100), r("chr1", 100_100), r("chrX", 1), r("chr1_alt", 7)]);
        assert_eq!(
            out,
            vec![
                Region { chromosome: "chr1".into(), start: 100, end: 100_100 },
                Region { chromosome: "chr1".into(), start: 300_000, end: 300_000 },
                r("chr2", 500),
                r("chr10", 5),
                r("chrX", 1),
                r("chr1_alt", 7),
            ]
        );
    }

    #[test]
    fn htsjdk_integer_types() {
        let f = |t: u8, v: &[u8]| [&b"AS"[..], &[t], v].concat();
        assert_eq!(retype_integer(&f(b'I', &240u32.to_le_bytes())), f(b'C', &[240]));
        assert_eq!(retype_integer(&f(b'C', &[60])), f(b'c', &[60]));
        assert_eq!(retype_integer(&f(b'i', &300i32.to_le_bytes())), f(b's', &300i16.to_le_bytes()));
        assert_eq!(retype_integer(&f(b'i', &40000i32.to_le_bytes())), f(b'S', &40000u16.to_le_bytes()));
        assert_eq!(retype_integer(&f(b'i', &70000i32.to_le_bytes())), f(b'i', &70000i32.to_le_bytes()));
        assert_eq!(retype_integer(&f(b'i', &(-5i32).to_le_bytes())), f(b'c', &[(-5i8) as u8]));
        assert_eq!(retype_integer(&f(b'i', &(-300i32).to_le_bytes())), f(b's', &(-300i16).to_le_bytes()));
        assert_eq!(retype_integer(&f(b'Z', b"x\0")), f(b'Z', b"x\0"));
    }

    #[test]
    fn thousands() {
        assert_eq!(group_thousands(1), "1");
        assert_eq!(group_thousands(9_999_999), "9,999,999");
        assert_eq!(partition_label("chr1", 0, 248_956_422), "chr1:1-9,999,999");
        assert_eq!(partition_label("chr21", 4, 46_709_983), "chr21:40,000,000-46,709,983");
    }
}
