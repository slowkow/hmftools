//! ResultsWriter and BamWriter: the junction and read TSVs and the prep BAM of each sample.
//!
//! Java writes each sample's records to an unsorted BAM as they come and then has samtools sort and index it.
//! Here the records are sorted in the process, in the order `samtools sort` gives: by reference (unplaced
//! last), position and strand, records that tie staying in the order they were written. Runs of records are
//! spilled to temporary files once they pass a size limit and merged at the end. Without a BAM tool on the
//! command line Java leaves the unsorted BAM, and so does this.

use super::config::*;
use super::rec::*;
use super::tracker::JunctionTracker;
use super::types::*;
use anyhow::{bail, Result};
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;
use std::ffi::CString;
use std::io::Write;
use std::sync::Arc;

/// A BGZF stream of BAM records, with or without the BAM header.
pub struct RawBamOut {
    fp: *mut htslib::BGZF,
    pub path: String,
}

unsafe impl Send for RawBamOut {}

impl RawBamOut {
    pub fn create(path: &str, level: u32) -> Result<Self> {
        let cpath = CString::new(path)?;
        let mode = CString::new(format!("w{level}"))?;
        let fp = unsafe { htslib::bgzf_open(cpath.as_ptr(), mode.as_ptr()) };
        if fp.is_null() {
            bail!("failed to open {path} for writing");
        }
        Ok(RawBamOut { fp, path: path.to_string() })
    }
    /// Compresses on `threads` extra threads.
    pub fn set_threads(&mut self, threads: usize) {
        if threads > 1 {
            unsafe { htslib::bgzf_mt(self.fp, threads as i32, 256) };
        }
    }
    pub fn write_header(&mut self, header: &bam::HeaderView) -> Result<()> {
        if unsafe { htslib::bam_hdr_write(self.fp, header.inner_ptr()) } < 0 {
            bail!("failed to write the BAM header to {}", self.path);
        }
        Ok(())
    }
    #[inline]
    pub fn write(&mut self, rec: &bam::Record) -> Result<()> {
        if unsafe { htslib::bam_write1(self.fp, &rec.inner) } < 0 {
            bail!("failed to write a record to {}", self.path);
        }
        Ok(())
    }
    pub fn close(mut self) -> Result<()> {
        let fp = std::mem::replace(&mut self.fp, std::ptr::null_mut());
        if unsafe { htslib::bgzf_close(fp) } < 0 {
            bail!("failed to close {}", self.path);
        }
        Ok(())
    }
}

impl Drop for RawBamOut {
    fn drop(&mut self) {
        if !self.fp.is_null() {
            unsafe { htslib::bgzf_close(self.fp) };
        }
    }
}

/// Builds the .bai index of a coordinate-sorted BAM.
pub fn index_bam(path: &str, threads: usize) -> Result<()> {
    let cpath = CString::new(path)?;
    let r = unsafe { htslib::sam_index_build3(cpath.as_ptr(), std::ptr::null(), 0, threads as i32) };
    if r < 0 {
        bail!("failed to index {path}");
    }
    Ok(())
}

/// The header text with the given sort order in its @HD line, as htsjdk and then samtools leave it.
pub fn header_with_sort_order(text: &str, sort_order: &str) -> String {
    let mut out = String::with_capacity(text.len() + 32);
    let mut has_hd = false;
    for line in text.lines() {
        if line.starts_with("@HD") {
            has_hd = true;
            let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            match fields.iter_mut().find(|f| f.starts_with("SO:")) {
                Some(f) => *f = format!("SO:{sort_order}"),
                None => fields.push(format!("SO:{sort_order}")),
            }
            out.push_str(&fields.join("\t"));
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }
    if !has_hd {
        out = format!("@HD\tVN:1.6\tSO:{sort_order}\n{out}");
    }
    out
}

/// BamMerger.buildCombinedHeader, for the cache BAMs: the first input's header with the read groups of the
/// other inputs that it lacks, and the first program record of each other input under the id `<id>.<index>`.
pub fn combined_header_text(headers: &[String]) -> String {
    let id_of = |line: &str| line.split('\t').find_map(|f| f.strip_prefix("ID:")).map(str::to_string);
    let mut lines: Vec<String> = headers[0].lines().map(str::to_string).collect();
    for (i, text) in headers.iter().enumerate().skip(1) {
        let mut first_pg = true;
        for line in text.lines() {
            if line.starts_with("@RG") {
                let id = id_of(line);
                if !lines.iter().any(|l| l.starts_with("@RG") && id_of(l) == id) {
                    let at = lines.iter().rposition(|l| l.starts_with("@RG")).map_or(lines.len(), |p| p + 1);
                    lines.insert(at, line.to_string());
                }
            } else if line.starts_with("@PG") && first_pg {
                first_pg = false;
                let Some(id) = id_of(line) else { continue };
                let new_id = format!("{id}.{i}");
                if !lines.iter().any(|l| l.starts_with("@PG") && id_of(l).as_deref() == Some(&new_id)) {
                    let fields: Vec<String> = line.split('\t').map(|f| if f.starts_with("ID:") { format!("ID:{new_id}") } else { f.to_string() }).collect();
                    let at = lines.iter().rposition(|l| l.starts_with("@PG")).map_or(lines.len(), |p| p + 1);
                    lines.insert(at, fields.join("\t"));
                }
            }
        }
    }
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

/// The key of `samtools sort`: reference as unsigned (unplaced last), position + 1, reverse strand.
#[inline]
fn sort_key(rec: &bam::Record) -> (u32, u64, bool) {
    (rec.tid() as u32, (rec.pos() + 1) as u64, rec.flags() & FLAG_REVERSE != 0)
}

/// the size in bytes of held records at which a run is written out
const RUN_BYTES: usize = 768 << 20;

struct SortedSink {
    path: String,
    header: bam::HeaderView,
    records: Vec<bam::Record>,
    bytes: usize,
    runs: Vec<String>,
    run_bytes: usize,
}

impl SortedSink {
    fn push(&mut self, rec: bam::Record) -> Result<()> {
        self.bytes += rec.inner.l_data as usize + 96;
        self.records.push(rec);
        if self.bytes >= self.run_bytes {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> Result<()> {
        let path = format!("{}.tmp.{:04}.bam", self.path, self.runs.len());
        self.records.sort_by_key(sort_key);
        let mut out = RawBamOut::create(&path, 1)?;
        out.set_threads(4);
        out.write_header(&self.header)?;
        for rec in self.records.drain(..) {
            out.write(&rec)?;
        }
        out.close()?;
        self.bytes = 0;
        self.runs.push(path);
        Ok(())
    }

    fn finish(mut self, threads: usize) -> Result<()> {
        let mut out = RawBamOut::create(&self.path, 6)?;
        out.set_threads(threads);
        out.write_header(&self.header)?;
        self.records.sort_by_key(sort_key);
        if self.runs.is_empty() {
            for rec in &self.records {
                out.write(rec)?;
            }
        } else {
            // the runs in the order written, then the records still held: equal keys come out in that order
            let mut readers: Vec<bam::Reader> = Vec::new();
            for run in &self.runs {
                readers.push(bam::Reader::from_path(run)?);
            }
            let n = readers.len();
            let mut heads: Vec<Option<bam::Record>> = Vec::with_capacity(n + 1);
            for r in readers.iter_mut() {
                let mut rec = bam::Record::new();
                heads.push(match r.read(&mut rec) {
                    Some(x) => {
                        x?;
                        Some(rec)
                    }
                    None => None,
                });
            }
            let mut held = std::mem::take(&mut self.records).into_iter();
            heads.push(held.next());
            loop {
                let mut best: Option<usize> = None;
                for (i, h) in heads.iter().enumerate() {
                    if let Some(rec) = h {
                        if best.is_none_or(|b| sort_key(rec) < sort_key(heads[b].as_ref().unwrap())) {
                            best = Some(i);
                        }
                    }
                }
                let Some(b) = best else { break };
                out.write(heads[b].as_ref().unwrap())?;
                if b < n {
                    let rec = heads[b].as_mut().unwrap();
                    match readers[b].read(rec) {
                        Some(x) => x?,
                        None => heads[b] = None,
                    }
                } else {
                    heads[b] = held.next();
                }
            }
            for run in &self.runs {
                let _ = std::fs::remove_file(run);
            }
        }
        out.close()?;
        index_bam(&self.path, threads)
    }
}

enum Sink {
    Sorted(SortedSink),
    Unsorted(RawBamOut),
}

/// BamWriter: one output per sample
pub struct BamWriter {
    sinks: Vec<Sink>,
    pub record_write_count: u64,
    /// records written per output
    pub sink_counts: Vec<u64>,
    threads: usize,
}

impl BamWriter {
    /// `headers` holds the header text of each sample's input BAM.
    pub fn new(config: &PrepConfig, headers: &[String]) -> Result<Self> {
        let mut sinks = Vec::new();
        if config.has_write_type(WriteType::PrepBam) {
            // the tools the Java step knows: samtools sorts and indexes, sambamba sorts
            let sorted = config.bam_tool_path.is_some();
            for (i, sample) in config.sample_ids.iter().enumerate() {
                if sorted {
                    let text = header_with_sort_order(&headers[i], "coordinate");
                    sinks.push(Sink::Sorted(SortedSink {
                        path: config.form_filename(WriteType::PrepBam, sample),
                        header: bam::HeaderView::from_bytes(text.as_bytes()),
                        records: Vec::new(),
                        bytes: 0,
                        runs: Vec::new(),
                        run_bytes: std::env::var("ESVEE_PREP_SORT_RUN_MB").ok().and_then(|v| v.parse::<usize>().ok()).map_or(RUN_BYTES, |mb| mb << 20),
                    }));
                } else {
                    let text = header_with_sort_order(&headers[i], "unsorted");
                    let mut out = RawBamOut::create(&config.form_filename(WriteType::UnsortedBam, sample), 6)?;
                    out.set_threads(config.threads.min(4));
                    out.write_header(&bam::HeaderView::from_bytes(text.as_bytes()))?;
                    sinks.push(Sink::Unsorted(out));
                }
            }
        }
        let sink_counts = vec![0; sinks.len()];
        Ok(BamWriter { sinks, record_write_count: 0, sink_counts, threads: config.threads })
    }

    pub fn is_active(&self) -> bool {
        !self.sinks.is_empty()
    }

    /// writeRecord: the record goes out in the form htsjdk writes once an attribute has been set on it, with
    /// its attributes in tag order and each integer in the smallest type.
    pub fn write_record(&mut self, sample: u8, mut rec: bam::Record) -> Result<()> {
        if self.sinks.is_empty() {
            return Ok(());
        }
        self.record_write_count += 1;
        decode_attributes(&mut rec);
        let index = if self.sinks.len() == 1 { 0 } else { sample as usize };
        self.sink_counts[index] += 1;
        match &mut self.sinks[index] {
            Sink::Sorted(s) => s.push(rec),
            Sink::Unsorted(out) => out.write(&rec),
        }
    }

    pub fn close(self) -> Result<()> {
        if self.sinks.is_empty() {
            return Ok(());
        }
        crate::info!("{} records written to BAM", self.record_write_count);
        for sink in self.sinks {
            match sink {
                Sink::Sorted(s) => {
                    crate::debug!("writing sorted BAM: {}", s.path);
                    s.finish(self.threads)?
                }
                Sink::Unsorted(out) => out.close()?,
            }
        }
        Ok(())
    }
}

/// A junction row that belongs elsewhere in the junction file than where its partition writes.
pub struct DeferredJunction {
    /// HumanChromosome.chromosomeRank of the chromosome, the order of the file
    pub rank: i32,
    pub chromosome: String,
    pub position: i32,
    pub orient: i32,
    pub text: Vec<u8>,
}

/// The lines of a junction file with the deferred rows put in: each after the rows of its chromosome up to
/// its position. A deferred row for a site (chromosome, position, orientation) that the file or an earlier
/// deferred row already has is left out.
pub fn merge_junction_lines(lines: impl Iterator<Item = Vec<u8>>, mut deferred: Vec<DeferredJunction>, out: &mut impl Write) -> std::io::Result<()> {
    deferred.sort_by(|a, b| (a.rank, &a.chromosome, a.position).cmp(&(b.rank, &b.chromosome, b.position)));
    let mut deferred = deferred.into_iter().peekable();
    // the site last written and the orientations written for it
    let mut site: Option<(i32, String, i32)> = None;
    let mut orients: Vec<i32> = Vec::new();
    let mut emit = |key: (i32, String, i32), orient: i32, text: &[u8], is_deferred: bool, out: &mut dyn Write| -> std::io::Result<()> {
        if site.as_ref() == Some(&key) {
            if is_deferred && orients.contains(&orient) {
                return Ok(());
            }
        } else {
            site = Some(key);
            orients.clear();
        }
        orients.push(orient);
        out.write_all(text)
    };
    for (n, mut line) in lines.enumerate() {
        if !line.ends_with(b"\n") {
            line.push(b'\n');
        }
        if n == 0 {
            out.write_all(&line)?;
            continue;
        }
        let mut fields = line.split(|&c| c == b'\t');
        let chromosome = String::from_utf8_lossy(fields.next().unwrap_or(b"")).into_owned();
        let number = |f: Option<&[u8]>| f.and_then(|x| std::str::from_utf8(x).ok()).and_then(|x| x.trim().parse::<i32>().ok()).unwrap_or(0);
        let position = number(fields.next());
        let orient = number(fields.next());
        let key = (chromosome_rank(&chromosome), chromosome, position);
        while deferred.peek().is_some_and(|d| (d.rank, &d.chromosome, d.position) < (key.0, &key.1, key.2)) {
            let d = deferred.next().unwrap();
            emit((d.rank, d.chromosome, d.position), d.orient, &d.text, true, out)?;
        }
        emit(key, orient, &line, false, out)?;
    }
    for d in deferred {
        emit((d.rank, d.chromosome, d.position), d.orient, &d.text, true, out)?;
    }
    Ok(())
}

fn merge_deferred_junctions(path: &str, deferred: Vec<DeferredJunction>) -> Result<()> {
    use std::io::BufRead;
    let tmp = format!("{path}.tmp");
    {
        let reader = std::io::BufReader::new(std::fs::File::open(path)?);
        let mut out = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
        let mut failed: Option<std::io::Error> = None;
        let lines = reader.split(b'\n').map_while(|l| match l {
            Ok(l) => Some(l),
            Err(e) => {
                failed = Some(e);
                None
            }
        });
        merge_junction_lines(lines, deferred, &mut out)?;
        if let Some(e) = failed {
            return Err(e.into());
        }
        out.flush()?;
    }
    std::fs::rename(&tmp, path)?;
    Ok(())
}

pub struct ResultsWriter {
    junction_file: String,
    deferred_junctions: Vec<DeferredJunction>,
    junction_writer: Option<std::io::BufWriter<std::fs::File>>,
    read_writer: Option<std::io::BufWriter<std::fs::File>>,
    pub bam_writer: BamWriter,
    sample_ids: Vec<String>,
    contigs: Arc<Contigs>,
    saga: bool,
}

impl ResultsWriter {
    pub fn new(config: &PrepConfig, contigs: Arc<Contigs>, headers: &[String]) -> Result<Self> {
        let mut junction_writer = None;
        let junction_file = config.form_filename(WriteType::PrepJunction, config.sample_id());
        if config.has_write_type(WriteType::PrepJunction) {
            let mut w = std::io::BufWriter::new(std::fs::File::create(&junction_file)?);
            let mut fields = vec![
                "Chromosome", "Position", "Orientation", "JunctionFrags", "ExactSupportFrags", "OtherSupportFrags", "RemoteFrags", "Indel", "Hotspot",
                "LowMapQualFrags", "InitialReadId", "Depth", "MaxValidSoftClip",
            ];
            if config.saga_fasta_file.is_some() {
                fields.push("SagaMatch");
            }
            writeln!(w, "{}", fields.join("\t"))?;
            junction_writer = Some(w);
        }
        let mut read_writer = None;
        if config.has_write_type(WriteType::PrepRead) {
            let mut w = std::io::BufWriter::new(std::fs::File::create(config.form_filename(WriteType::PrepRead, config.sample_id()))?);
            writeln!(
                w,
                "ReadId\tSampleId\tGroupCount\tExpectedCount\tGroupStatus\tHasExternal\tReadType\tChromosome\tPosStart\tPosEnd\tCigar\tFragLength\tMateChr\tMatePosStart\tMapQual\tSuppData\tFlags\tFilters\tFirstInPair\tReadReversed\tMateReversed\tUnmapped\tMateUnmapped\tSupplementary\tJunctionPositions"
            )?;
            read_writer = Some(w);
        }
        Ok(ResultsWriter {
            junction_file,
            deferred_junctions: Vec::new(),
            junction_writer,
            read_writer,
            bam_writer: BamWriter::new(config, headers)?,
            sample_ids: config.sample_ids.clone(),
            contigs,
            saga: config.saga_fasta_file.is_some(),
        })
    }

    pub fn written_count(&self) -> u64 {
        self.bam_writer.record_write_count
    }

    #[allow(clippy::too_many_arguments)]
    fn write_read_data(&mut self, read: &PrepRead, read_count: usize, expected_read_count: i32, status: ReadGroupStatus, spans_partitions: bool, junction_positions: &str) -> Result<()> {
        let Some(w) = &mut self.read_writer else { return Ok(()) };
        let flags = read.flags();
        let fields: Vec<String> = vec![
            String::from_utf8_lossy(read.id()).into_owned(),
            self.sample_ids.get(read.sample as usize).cloned().unwrap_or_default(),
            read_count.to_string(),
            expected_read_count.to_string(),
            status.name().to_string(),
            spans_partitions.to_string(),
            read.read_type.name().to_string(),
            self.contigs.name(read.chr).to_string(),
            read.start.to_string(),
            read.end.to_string(),
            cigar_string(read.cigar()),
            read.fragment_insert_size.to_string(),
            self.contigs.name(read.mate_chr).to_string(),
            read.mate_pos_start.to_string(),
            read.map_quality().to_string(),
            read.supp.as_ref().map_or("N/A".to_string(), |s| s.as_delim_str()),
            flags.to_string(),
            read.filters.to_string(),
            read.is_first_of_pair().to_string(),
            read.is_read_reversed().to_string(),
            (flags & FLAG_MATE_REVERSE != 0).to_string(),
            (flags & FLAG_UNMAPPED != 0).to_string(),
            (read.has_mate() && read.is_mate_unmapped()).to_string(),
            (flags & FLAG_SUPPLEMENTARY != 0).to_string(),
            junction_positions.to_string(),
        ];
        writeln!(w, "{}", fields.join("\t"))?;
        Ok(())
    }

    /// writeReadGroups: the groups that are not held back for their remote reads, less reads with poly-G
    /// soft clips and reads already written
    pub fn write_read_groups(&mut self, tracker: &mut JunctionTracker, groups: &[GroupId]) -> Result<()> {
        for &gid in groups {
            let g = &tracker.groups[gid as usize];
            if group_conditional_on_remote_reads(g, &tracker.reads) && !g.has_remote_junction_reads {
                continue;
            }
            if g.status == ReadGroupStatus::Duplicate {
                continue;
            }
            if self.read_writer.is_some() {
                let junction_positions = g.junction_positions_str();
                let (size, expected, status, spans) = (g.size(), g.expected_read_count, g.status, g.spans_partitions());
                for k in 0..size {
                    let rid = tracker.groups[gid as usize].reads[k];
                    if tracker.reads[rid as usize].written {
                        continue;
                    }
                    self.write_read_data(&tracker.reads[rid as usize], size, expected, status, spans, &junction_positions)?;
                }
            }
            for k in 0..tracker.groups[gid as usize].reads.len() {
                let rid = tracker.groups[gid as usize].reads[k];
                let read = &mut tracker.reads[rid as usize];
                if !read.has_filter(F_POLY_G_SC) && !read.written && self.bam_writer.is_active() {
                    let rec = std::mem::replace(&mut read.rec, bam::Record::new());
                    self.bam_writer.write_record(read.sample, rec)?;
                }
                read.written = true;
            }
        }
        Ok(())
    }

    /// A candidate read claimed from the cache: a group of its own, expected, with remote junction reads
    pub fn write_candidate_read(&mut self, read: PrepRead) -> Result<()> {
        self.write_read_data(&read, 1, 0, ReadGroupStatus::Expected, false, "")?;
        if !read.has_filter(F_POLY_G_SC) && self.bam_writer.is_active() {
            self.bam_writer.write_record(read.sample, read.rec)?;
        }
        Ok(())
    }

    /// The junction file rows of a partition (writeJunctionData). Made before the partition's reads are
    /// written, since that hands their records on.
    pub fn junction_data_text(&self, chromosome: &str, tracker: &JunctionTracker) -> Vec<u8> {
        let mut w: Vec<u8> = Vec::new();
        if self.junction_writer.is_none() {
            return w;
        }
        for &jid in &tracker.junctions {
            self.junction_row(&mut w, chromosome, tracker, jid);
        }
        w
    }

    /// The rows of the junctions that lie outside the partition that found them (finding P8, corrected): each
    /// under its own chromosome, to be put in its place when the file is closed.
    pub fn outside_junction_rows(&self, tracker: &JunctionTracker) -> Vec<DeferredJunction> {
        let mut rows = Vec::new();
        if self.junction_writer.is_none() {
            return rows;
        }
        for &jid in &tracker.outside_junctions {
            let j = &tracker.jstore[jid as usize];
            let chromosome = self.contigs.name(j.outside_chr.unwrap_or(tracker.region.chr)).to_string();
            let mut text = Vec::new();
            self.junction_row(&mut text, &chromosome, tracker, jid);
            rows.push(DeferredJunction { rank: chromosome_rank(&chromosome), chromosome, position: j.position, orient: j.orient.as_byte(), text });
        }
        rows
    }

    fn junction_row(&self, w: &mut Vec<u8>, chromosome: &str, tracker: &JunctionTracker, jid: JuncId) {
        {
            let j = &tracker.jstore[jid as usize];
            let mut low_map_qual_frags = 0;
            let mut disc_remote_frags = 0;
            if !j.discordant_group {
                let expect_left_clipped = j.orient.is_reverse();
                for &rid in &j.junction_reads {
                    let read = &tracker.reads[rid as usize];
                    // the read can also be at another junction
                    let supports_junction = (expect_left_clipped && read.start == j.position && read.is_left_clipped())
                        || (!expect_left_clipped && read.end == j.position && read.is_right_clipped());
                    if supports_junction && read.has_filter(F_MIN_MAP_QUAL) {
                        low_map_qual_frags += 1;
                    }
                }
                low_map_qual_frags += j.exact_support_reads.iter().filter(|&&rid| tracker.reads[rid as usize].has_filter(F_MIN_MAP_QUAL)).count();
            } else {
                disc_remote_frags = j.remote_junctions.first().map_or(0, |r| r.fragments);
            }
            let _ = write!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                chromosome,
                j.position,
                j.orient.as_byte(),
                j.junction_fragment_count(),
                j.exact_support_fragment_count(),
                j.supporting_fragment_count(),
                disc_remote_frags,
                j.internal_indel,
                j.hotspot,
                low_map_qual_frags
            );
            match j.top_junction_read {
                Some(rid) => w.extend_from_slice(tracker.reads[rid as usize].id()),
                None => w.extend_from_slice(b"EXISTING"),
            }
            let _ = write!(w, "\t{}\t{}", j.depth, j.max_valid_soft_clip_length);
            if self.saga {
                let _ = write!(w, "\t{}", j.saga_match.as_deref().unwrap_or(""));
            }
            w.push(b'\n');
        }
    }

    pub fn defer_junctions(&mut self, rows: Vec<DeferredJunction>) {
        self.deferred_junctions.extend(rows);
    }

    pub fn write_junction_text(&mut self, text: &[u8]) -> Result<()> {
        if let Some(w) = &mut self.junction_writer {
            w.write_all(text)?;
        }
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        if let Some(mut w) = self.read_writer.take() {
            w.flush()?;
        }
        if let Some(mut w) = self.junction_writer.take() {
            w.flush()?;
            drop(w);
            if !self.deferred_junctions.is_empty() {
                merge_deferred_junctions(&self.junction_file, std::mem::take(&mut self.deferred_junctions))?;
            }
        }
        self.bam_writer.close()
    }
}
