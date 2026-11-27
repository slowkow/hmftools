//! BamReadCounter and ReadDepthAccumulator: the bases of the aligned blocks of every read that passes the
//! filters, and how many of them are G or C, summed per 1,000-base window of each human chromosome.
//!
//! As in Java, each chromosome is read in partitions from position 1 (10 Mb here, 100 Mb in Java), each partition by one task
//! with its own reader, counting only the bases inside the partition. Partitions start at window starts,
//! so every window belongs to one partition: a task writes its own windows and the counts do not depend on
//! the number of threads or the order of the tasks. The tasks of the tumor and the reference BAM share one
//! pool of `threads` workers, as Java's executor does. Records are read into one reused buffer per worker;
//! the GC bases are counted on the packed 4-bit sequence.

use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::chromosome::Chr;
use crate::constants::{PARTITION_SIZE, WINDOW_SIZE};

const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// The read counts of one chromosome: bases and G/C bases per window (the partial window at the end of the
/// chromosome is left out, as Java does).
pub struct ChrCounts {
    pub chr: Chr,
    /// the contig name in the BAM header
    pub contig: String,
    pub bases: Vec<i32>,
    pub gc: Vec<i32>,
}

/// The counts of one BAM: its human chromosomes in header order.
pub struct SampleCounts {
    pub chromosomes: Vec<ChrCounts>,
}

impl SampleCounts {
    /// The number of DepthReadings (one per window)
    pub fn window_count(&self) -> usize {
        self.chromosomes.iter().map(|c| c.bases.len()).sum()
    }
}

pub struct ReadFilter {
    pub min_mapping_quality: i32,
    pub include_duplicates: bool,
}

/// What one BAM needs: its path, the reference for a CRAM, and the chromosomes to count.
pub struct BamJob<'a> {
    pub path: &'a str,
    pub chromosomes: Vec<(Chr, String, i32)>,
}

/// BamReadCounter.loadChromosomes: the human chromosomes of the header, in its order, without those the
/// specific chromosomes leave out. With `-specific_chr` or `-specific_regions` only the chromosomes are
/// taken (whole); the regions do not limit the reading.
pub fn load_chromosomes(path: &str, ref_genome: Option<&str>, specific: &hmf_common::regions::SpecificRegions) -> Result<Vec<(Chr, String, i32)>> {
    let reader = open_reader(path, ref_genome)?;
    let header = reader.header();
    let mut out = Vec::new();
    for tid in 0..header.target_count() {
        let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
        let Some(chr) = Chr::from_contig(&name) else { continue };
        if specific.has_filters() && !specific.include_chromosome(&name) {
            continue;
        }
        let length = header.target_len(tid).unwrap_or(0) as i32;
        out.push((chr, name, length));
    }
    Ok(out)
}

/// The file buffer of each reader: few large reads keep a disk shared with other jobs streaming
const READ_BUFFER: i32 = 4 << 20;

pub fn open_reader(path: &str, ref_genome: Option<&str>) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening {path}"))?;
    unsafe {
        use rust_htslib::htslib;
        htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_HTS_OPT_BLOCK_SIZE, READ_BUFFER as std::os::raw::c_int);
    }
    if let Some(r) = ref_genome.filter(|r| !r.is_empty() && path.ends_with(".cram")) {
        reader.set_reference(r).with_context(|| format!("setting CRAM reference {r}"))?;
    }
    Ok(reader)
}

struct Task {
    job: usize,
    chromosome: usize,
    start: i32,
    end: i32,
}

/// Counts the reads of every job's BAM: returns one SampleCounts per job, in job order.
pub fn count_reads(jobs: &[BamJob], ref_genome: Option<&str>, filter: &ReadFilter, threads: usize) -> Result<Vec<SampleCounts>> {
    // the window arrays of each chromosome of each job; tasks get disjoint slices of them
    let mut results: Vec<SampleCounts> = jobs
        .iter()
        .map(|j| SampleCounts {
            chromosomes: j
                .chromosomes
                .iter()
                .map(|(chr, name, length)| {
                    let windows = (length / WINDOW_SIZE) as usize;
                    ChrCounts { chr: *chr, contig: name.clone(), bases: vec![0; windows], gc: vec![0; windows] }
                })
                .collect(),
        })
        .collect();
    // BamReadCounter.partitionGenome
    let mut tasks: Vec<Task> = Vec::new();
    for (ji, job) in jobs.iter().enumerate() {
        crate::info!("calculating read depths from {}", job.path);
        for (ci, (_, _, length)) in job.chromosomes.iter().enumerate() {
            let mut start = 1;
            while start < *length {
                let end = (start + PARTITION_SIZE - 1).min(*length);
                tasks.push(Task { job: ji, chromosome: ci, start, end });
                start += PARTITION_SIZE;
            }
        }
    }
    // the longest first, to balance the workers
    tasks.sort_by_key(|t| std::cmp::Reverse(t.end - t.start));

    // the slices of the window arrays, one per task: partitions start at window boundaries
    type Slice<'a> = Option<(&'a mut [i32], &'a mut [i32])>;
    type Slices<'a> = Vec<Slice<'a>>;
    let mut slices: Slices = Vec::with_capacity(tasks.len());
    {
        let mut by_key: std::collections::HashMap<(usize, usize, i32), usize> = std::collections::HashMap::new();
        for (i, t) in tasks.iter().enumerate() {
            by_key.insert((t.job, t.chromosome, t.start), i);
            slices.push(None);
        }
        for (ji, sample) in results.iter_mut().enumerate() {
            for (ci, counts) in sample.chromosomes.iter_mut().enumerate() {
                let mut bases: &mut [i32] = &mut counts.bases;
                let mut gc: &mut [i32] = &mut counts.gc;
                let mut start = 1;
                let length = jobs[ji].chromosomes[ci].2;
                while start < length {
                    let window_count = (PARTITION_SIZE / WINDOW_SIZE) as usize;
                    let take = window_count.min(bases.len());
                    let (b, rest_b) = std::mem::take(&mut bases).split_at_mut(take);
                    let (g, rest_g) = std::mem::take(&mut gc).split_at_mut(take);
                    bases = rest_b;
                    gc = rest_g;
                    slices[by_key[&(ji, ci, start)]] = Some((b, g));
                    start += PARTITION_SIZE;
                }
            }
        }
    }
    let work: Vec<Mutex<Slice>> = slices.into_iter().map(Mutex::new).collect();
    let next = AtomicUsize::new(0);
    let threads = threads.max(1).min(tasks.len().max(1));
    let outcome: Vec<Result<()>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let (tasks, work, next) = (&tasks, &work, &next);
                scope.spawn(move || -> Result<()> {
                    let mut readers: Vec<Option<bam::IndexedReader>> = (0..jobs.len()).map(|_| None).collect();
                    let mut record = bam::Record::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(task) = tasks.get(i) else { break };
                        let (bases, gc) = work[i].lock().unwrap().take().unwrap();
                        let job = &jobs[task.job];
                        if readers[task.job].is_none() {
                            let mut r = open_reader(job.path, ref_genome)?;
                            if job.path.ends_with(".cram") {
                                set_cram_fields(&mut r, filter.include_duplicates);
                            }
                            readers[task.job] = Some(r);
                        }
                        let reader = readers[task.job].as_mut().unwrap();
                        let contig = &job.chromosomes[task.chromosome].1;
                        crate::debug!("region({contig}:{}-{}) accumulating read depth", task.start, task.end);
                        count_partition(reader, &mut record, contig, task.start, task.end, filter, bases, gc).with_context(|| format!("reading {}", job.path))?;
                        crate::debug!("region({contig}:{}-{}) complete", task.start, task.end);
                    }
                    Ok(())
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for r in outcome {
        r?;
    }
    Ok(results)
}

/// For CRAM: decode only the fields the counting reads (no names or qualities).
fn set_cram_fields(reader: &mut bam::IndexedReader, aux: bool) {
    use rust_htslib::htslib;
    let mut fields = htslib::sam_fields_SAM_FLAG
        | htslib::sam_fields_SAM_RNAME
        | htslib::sam_fields_SAM_POS
        | htslib::sam_fields_SAM_MAPQ
        | htslib::sam_fields_SAM_CIGAR
        | htslib::sam_fields_SAM_SEQ;
    if aux {
        fields |= htslib::sam_fields_SAM_AUX;
    }
    unsafe {
        htslib::hts_set_opt(reader.htsfile(), htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields as std::os::raw::c_int);
    }
}

/// The number of G and C bases (codes 2 and 4 of =ACMGRSVTWYHKDBN) in each byte of two packed bases
static GC_PAIR: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let (hi, lo) = (i >> 4, i & 15);
        t[i] = ((hi == 2 || hi == 4) as u8) + ((lo == 2 || lo == 4) as u8);
        i += 1;
    }
    t
};

/// G or C bases among read bases [from, to) of a packed sequence (high nibble first)
#[inline]
fn gc_bases(packed: &[u8], from: usize, to: usize) -> u32 {
    let mut n = 0u32;
    let mut i = from;
    if i < to && i % 2 == 1 {
        let c = packed[i / 2] & 15;
        n += (c == 2 || c == 4) as u32;
        i += 1;
    }
    if i >= to {
        return n;
    }
    // whole bytes, eight at a time: a nibble is C (2) or G (4) when it xors one of them to zero
    const LOW: u64 = 0x1111_1111_1111_1111;
    let mut b = i / 2;
    let end_byte = to / 2;
    while b + 8 <= end_byte {
        let x = u64::from_le_bytes(packed[b..b + 8].try_into().unwrap());
        let nonzero = |y: u64| (y | (y >> 1) | (y >> 2) | (y >> 3)) & LOW;
        let c = nonzero(x ^ (LOW * 2));
        let g = nonzero(x ^ (LOW * 4));
        n += (!(c & g) & LOW).count_ones();
        b += 8;
    }
    while b < end_byte {
        n += GC_PAIR[packed[b] as usize] as u32;
        b += 1;
    }
    i = b * 2;
    if i < to {
        let c = packed[i / 2] >> 4;
        n += (c == 2 || c == 4) as u32;
    }
    n
}

/// One partition [start, end] of a chromosome: BamSlicer(minMappingQuality, keepDuplicates, false, false)
/// over the reads overlapping it, then processRead and accumulateAlignmentBlock for each.
#[allow(clippy::too_many_arguments)]
fn count_partition(
    reader: &mut bam::IndexedReader,
    record: &mut bam::Record,
    contig: &str,
    start: i32,
    end: i32,
    filter: &ReadFilter,
    bases: &mut [i32],
    gc: &mut [i32],
) -> Result<()> {
    let Some(tid) = reader.header().tid(contig.as_bytes()) else {
        crate::warn!("cannot find sequence index for chromosome {contig} in bam header");
        return Ok(());
    };
    reader.fetch((tid, (start - 1) as i64, end as i64))?;
    // the first window of the slice
    let first_window = ((start - 1) / WINDOW_SIZE) as i64;
    let windows = bases.len() as i64;
    while let Some(r) = reader.read(record) {
        r?;
        let flags = record.flags();
        if (record.mapq() as i32) < filter.min_mapping_quality || flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            continue;
        }
        if filter.include_duplicates {
            // only the raw reads: consensus reads (attribute CR) are left out
            if record.aux(b"CR").is_ok() {
                continue;
            }
        } else if flags & FLAG_DUPLICATE != 0 {
            continue;
        }
        let packed = record.seq().encoded;
        let read_len = record.seq_len();
        let mut ref_pos = record.pos() as i32 + 1;
        let mut read_pos = 0usize;
        for &op in record.raw_cigar() {
            let len = (op >> 4) as i32;
            match op & 15 {
                // M, =, X: an alignment block
                0 | 7 | 8 => {
                    let genome_start = ref_pos.max(start);
                    let block_end = (ref_pos + len).min(end + 1);
                    let length = block_end - genome_start;
                    if length > 0 {
                        let read_start = read_pos + (genome_start - ref_pos) as usize;
                        if read_start + length as usize > read_len {
                            bail!("java.lang.ArrayIndexOutOfBoundsException: read bases shorter than its alignment");
                        }
                        add_block(bases, gc, first_window, windows, genome_start, length, packed, read_start);
                    }
                    ref_pos += len;
                    read_pos += len as usize;
                }
                // I, S: read bases only
                1 | 4 => read_pos += len as usize,
                // D, N: reference only
                2 | 3 => ref_pos += len,
                // H, P
                _ => {}
            }
        }
    }
    Ok(())
}

/// ReadDepthAccumulator.addReadAlignmentToCounts for one block within the partition (int counters, as Java's,
/// which wrap on overflow)
#[allow(clippy::too_many_arguments)]
#[inline]
fn add_block(bases: &mut [i32], gc: &mut [i32], first_window: i64, windows: i64, genome_start: i32, length: i32, packed: &[u8], read_start: usize) {
    let block_end = genome_start + length; // exclusive
    let mut window = ((genome_start - 1) / WINDOW_SIZE) as i64;
    loop {
        let window_start = (window * WINDOW_SIZE as i64 + 1) as i32;
        if window_start >= block_end {
            break;
        }
        let local = window - first_window;
        if local >= windows {
            // the partial window at the end of the chromosome
            break;
        }
        let from = genome_start.max(window_start);
        let to = block_end.min(window_start + WINDOW_SIZE);
        let n = to - from;
        let offset = read_start + (from - genome_start) as usize;
        bases[local as usize] = bases[local as usize].wrapping_add(n);
        gc[local as usize] = gc[local as usize].wrapping_add(gc_bases(packed, offset, offset + n as usize) as i32);
        window += 1;
    }
}

/// add_block for a block of unpacked bases at `genome_start` (1-based), the window arrays from window 0
#[cfg(test)]
pub fn add_block_for_test(bases: &mut [i32], gc: &mut [i32], genome_start: i32, length: i32, read: &[u8]) {
    let code = |b: u8| -> u8 { b"=ACMGRSVTWYHKDBN".iter().position(|&c| c == b).unwrap() as u8 };
    let packed: Vec<u8> = read.chunks(2).map(|p| (code(p[0]) << 4) | if p.len() > 1 { code(p[1]) } else { 0 }).collect();
    let windows = bases.len() as i64;
    add_block(bases, gc, 0, windows, genome_start, length, &packed, 0);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack(seq: &[u8]) -> Vec<u8> {
        let code = |b: u8| -> u8 { b"=ACMGRSVTWYHKDBN".iter().position(|&c| c == b).unwrap() as u8 };
        seq.chunks(2).map(|p| (code(p[0]) << 4) | if p.len() > 1 { code(p[1]) } else { 0 }).collect()
    }

    #[test]
    fn gc_count_on_packed_bases() {
        let seq: Vec<u8> = b"GCATNGGCCATGGCCGCGNAGTTCGCGCGAACGTGGGCCAAATCGCGTAGGCAGCCGGC=MRSVTWYHKDBN".to_vec();
        let p = pack(&seq);
        for from in 0..seq.len() {
            for to in from..=seq.len() {
                let expected = seq[from..to].iter().filter(|&&b| b == b'G' || b == b'C').count() as u32;
                assert_eq!(gc_bases(&p, from, to), expected, "{from}..{to}");
            }
        }
    }

    /// ReadDepthAccumulatorTest-style: a block across windows, a block past the last whole window
    #[test]
    fn blocks_split_at_windows() {
        let mut bases = vec![0i32; 3];
        let mut gc = vec![0i32; 3];
        let seq = vec![b'G'; 2500];
        let p = pack(&seq);
        // 1-based 901 to 2100: 100 bases in window 0, 1000 in window 1, 100 in window 2
        add_block(&mut bases, &mut gc, 0, 3, 901, 1200, &p, 0);
        assert_eq!(bases, vec![100, 1000, 100]);
        assert_eq!(gc, vec![100, 1000, 100]);
        // 2901 to 3200: window 3 does not exist
        add_block(&mut bases, &mut gc, 0, 3, 2901, 300, &p, 0);
        assert_eq!(bases, vec![100, 1000, 200]);
    }
}
