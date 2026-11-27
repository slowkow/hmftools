//! BamEvidenceReader and BamReaderThread: the sites of each chromosome are split into regions where they
//! lie closer than the minimum gap, and worker threads, each with its own reader, take the regions from a
//! queue and count the reads of each region at its sites. A region's sites are counted by its own task
//! only, so the counts do not depend on the number of threads or the order the tasks run in.
//!
//! Unlike Java, a queue item is a run of consecutive regions, and a BAM's item is read with one index
//! query over its span, each read given to the regions it overlaps (see `read_regions`). Every region
//! still sees exactly the reads of its own query, in the same order, so the counts are unchanged; the
//! reading becomes sequential, which matters most for a BAM read through a bucket mount, where each
//! query is a new ranged request.

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read as BamRead};
use std::sync::Mutex;

use crate::config::AmberConfig;
use crate::constants::{BAM_MIN_GAP_START, CRAM_MIN_GAP_START};
use crate::evidence::*;
use crate::sites::{versioned_chromosome, ByChr};

struct Task<'a> {
    chromosome: String,
    start: i32,
    end: i32,
    positions: &'a mut [Evidence],
}

fn is_cram(path: &str) -> bool {
    path.ends_with(".cram")
}

/// The region tasks of BamEvidenceReader.populateTaskQueue: (chromosome index, first and last site index,
/// region start and end) per task, the sites of a chromosome in their list order.
pub fn region_tasks(lists: &ByChr<Evidence>, min_gap: i32) -> Vec<(usize, usize, usize, i32, i32)> {
    let mut tasks = Vec::new();
    for (ci, positions) in lists.chromosomes() {
        let mut first = 0usize;
        let (mut start, mut end) = (positions[0].pos, positions[0].pos);
        for (i, p) in positions.iter().enumerate().skip(1) {
            if end + min_gap < p.pos {
                tasks.push((ci, first, i, start, end));
                first = i;
                start = p.pos;
                end = p.pos;
            } else {
                end = end.max(p.pos);
            }
        }
        tasks.push((ci, first, positions.len(), start, end));
    }
    tasks
}

/// Opens a BAM or CRAM for indexed reading, with the reference for a CRAM.
pub fn open_reader(path: &str, ref_genome: Option<&str>) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(path).with_context(|| format!("opening {path}"))?;
    if is_cram(path) {
        // without a reference htslib would look for one through REF_PATH, which by default downloads the
        // sequences from the EBI; htsjdk (Java) fails instead, and so does this
        let Some(r) = ref_genome else { anyhow::bail!("CRAM {path} needs -ref_genome") };
        reader.set_reference(r).with_context(|| format!("setting CRAM reference {r}"))?;
    }
    let block_kb = FILE_BLOCK_KB.load(std::sync::atomic::Ordering::Relaxed);
    if block_kb > 0 {
        // htslib reads the file in blocks of the file system's preferred size (often 32 KB or less); with two
        // readers each in its own part of a large file, bigger blocks mean far fewer requests (as sage-rs)
        // SAFETY: the option takes one int argument
        let ret = unsafe {
            rust_htslib::htslib::hts_set_opt(
                reader.htsfile(),
                rust_htslib::htslib::hts_fmt_option_HTS_OPT_BLOCK_SIZE,
                (block_kb << 10) as std::os::raw::c_int,
            )
        };
        if ret != 0 {
            anyhow::bail!("setting the file block size for {path}");
        }
    }
    Ok(reader)
}

/// The block size (KB) htslib reads BAM/CRAM files in; 0 for htslib's default. The environment variable
/// AMBER_FILE_BLOCK_KB sets it (not a Java option).
pub static FILE_BLOCK_KB: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(1024);

/// Grouping of the region tasks into queue items (not Java options; the environment variables
/// AMBER_STREAM_GAP_KB and AMBER_STREAM_CHUNK_MB set them). An item holds consecutive tasks of one
/// chromosome over at most STREAM_CHUNK_BP bases, with no gap above STREAM_GAP_BP between them. A gap
/// of 0 gives each task its own item and query, as Java reads.
pub static STREAM_GAP_BP: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(256_000);
pub static STREAM_CHUNK_BP: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(10_000_000);

/// The queue items: ranges of `tasks` (as from region_tasks, in order), a new item where the chromosome
/// changes, where the gap from the previous task's end is above `gap`, or where the item would span more
/// than `chunk` bases.
pub fn group_tasks(tasks: &[(usize, usize, usize, i32, i32)], gap: i32, chunk: i32) -> Vec<std::ops::Range<usize>> {
    let mut groups = Vec::new();
    let mut first = 0usize;
    for i in 1..tasks.len() {
        let (head, prev, t) = (&tasks[first], &tasks[i - 1], &tasks[i]);
        if t.0 != head.0 || t.3 - prev.4 > gap || t.4 as i64 - head.3 as i64 + 1 > chunk as i64 {
            groups.push(first..i);
            first = i;
        }
    }
    if !tasks.is_empty() {
        groups.push(first..tasks.len());
    }
    groups
}

/// Reads the regions of `tasks` (consecutive and disjoint, in position order) with one index query over
/// their span, and gives each read to every region it overlaps by htslib's own test for a query
/// [start - 1, end): read start < end and bam_endpos > start - 1. A query's reads are the BAM's reads
/// passing that test, in file order, so each region sees exactly the reads, in the same order, that a
/// query over its own region returns. Returns the reads counted (a read once per region given to).
fn read_regions(
    reader: &mut bam::IndexedReader,
    tid: u32,
    tasks: &mut [Task],
    record: &mut bam::Record,
    checker: &EvidenceChecker,
    scratch: &mut Vec<CigarElem>,
    bam_path: &str,
) -> Result<u64> {
    let (start, end) = (tasks[0].start, tasks[tasks.len() - 1].end);
    reader.fetch((tid, (start - 1) as i64, end as i64))?;
    let mut regions: Vec<RegionTask> = (0..tasks.len()).map(|_| RegionTask::new()).collect();
    let mut first = 0usize; // the first region not ended before the reads still to come
    let mut open = tasks.len(); // regions not complete
    let mut reads = 0u64;
    while let Some(r) = reader.read(record) {
        r.with_context(|| format!("reading {bam_path}"))?;
        if record.flags() & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_DUPLICATE) != 0 {
            continue;
        }
        let read_start = record.pos();
        // SAFETY: a record just read is a valid bam1_t
        let read_end = unsafe { rust_htslib::htslib::bam_endpos(record.inner()) };
        // the reads come in start order, so a region ending at or before this start has had all of its reads
        while first < tasks.len() && (tasks[first].end as i64) <= read_start {
            first += 1;
        }
        if first == tasks.len() {
            break;
        }
        let mut view: Option<RecView> = None;
        for t in first..tasks.len() {
            if tasks[t].start as i64 > read_end {
                break;
            }
            if regions[t].complete {
                continue;
            }
            reads += 1;
            let v = view.get_or_insert_with(|| RecView::of(record));
            regions[t].process_record(tasks[t].positions, v, checker, scratch);
            if regions[t].complete {
                open -= 1;
            }
        }
        if open == 0 {
            break;
        }
    }
    Ok(reads)
}

/// BamEvidenceReader.processBam: counts the reads of the BAM at every site of `lists`.
pub fn process_bam(config: &AmberConfig, bam_path: &str, lists: &mut ByChr<Evidence>) -> Result<()> {
    let started = std::time::Instant::now();
    let limit_slice_regions = is_cram(bam_path);
    let min_gap = if config.position_gap > 0 {
        config.position_gap
    } else if limit_slice_regions {
        CRAM_MIN_GAP_START
    } else {
        BAM_MIN_GAP_START
    };
    let task_list = region_tasks(lists, min_gap);
    if hmf_common::log::debug_enabled() {
        // the first task that is longer or has more sites than the one before
        let mut max: Option<&(usize, usize, usize, i32, i32)> = None;
        for t in &task_list {
            if max.is_none_or(|m| t.4 - t.3 + 1 > m.4 - m.3 + 1 || t.2 - t.1 > m.2 - m.1) {
                max = Some(t);
            }
        }
        if let Some(m) = max {
            crate::debug!(
                "split {} sites across {} regions, max region({}:{}-{} size={} length={}) minGap({})",
                lists.len(),
                task_list.len(),
                versioned_chromosome(config.v38, m.0),
                m.3,
                m.4,
                m.2 - m.1,
                m.4 - m.3 + 1,
                min_gap
            );
        }
    }

    // split each chromosome's sites into the tasks' disjoint slices
    let mut tasks: Vec<Task> = Vec::with_capacity(task_list.len());
    let mut by_chr: Vec<Vec<&(usize, usize, usize, i32, i32)>> = vec![Vec::new(); 25];
    for t in &task_list {
        by_chr[t.0].push(t);
    }
    for (ci, list) in lists.lists.iter_mut().enumerate() {
        // a chromosome's tasks cover its sites in consecutive ranges
        let mut rest: &mut [Evidence] = list.as_mut_slice();
        for t in &by_chr[ci] {
            let (positions, tail) = std::mem::take(&mut rest).split_at_mut(t.2 - t.1);
            rest = tail;
            tasks.push(Task { chromosome: versioned_chromosome(config.v38, ci), start: t.3, end: t.4, positions });
        }
    }
    let task_count = tasks.len();
    // `tasks` is in chromosome-index order, as `task_list`, so the groups index both
    let gap = STREAM_GAP_BP.load(std::sync::atomic::Ordering::Relaxed);
    let chunk = STREAM_CHUNK_BP.load(std::sync::atomic::Ordering::Relaxed);
    let groups = group_tasks(&task_list, gap, chunk);
    // a CRAM query filters with CRAM's own alignment end, so each CRAM region keeps its own query; an item
    // still gives one reader consecutive regions
    let query_per_task = limit_slice_regions || gap <= 0;
    let mut items: Vec<Vec<Task>> = Vec::with_capacity(groups.len());
    let mut it = tasks.into_iter();
    for g in &groups {
        items.push(it.by_ref().take(g.len()).collect());
    }
    let item_count = items.len();
    items.reverse(); // popped from the back, so taken in order
    let queue = Mutex::new(items);
    let checker = EvidenceChecker {
        min_map_quality: config.min_mapping_quality,
        min_base_quality: config.min_base_quality,
        ultima: config.is_ultima(),
    };
    let threads = config.reader_threads().min(item_count.max(1));
    let ref_genome = config.ref_genome_file.as_deref();
    let results: Vec<Result<u64>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let queue = &queue;
                let checker = &checker;
                scope.spawn(move || -> Result<u64> {
                    let mut reader = open_reader(bam_path, ref_genome)?;
                    let mut record = bam::Record::new();
                    let mut scratch: Vec<CigarElem> = Vec::with_capacity(16);
                    let mut reads = 0u64;
                    let mut warned: Vec<String> = Vec::new();
                    loop {
                        let item = queue.lock().unwrap().pop();
                        let Some(mut item) = item else { break };
                        let chromosome = item[0].chromosome.clone();
                        let Some(tid) = reader.header().tid(chromosome.as_bytes()) else {
                            if !warned.contains(&chromosome) {
                                crate::warn!("cannot find sequence index for chromosome {} in bam header", chromosome);
                                warned.push(chromosome);
                            }
                            continue;
                        };
                        if query_per_task {
                            for task in item.iter_mut() {
                                let one = std::slice::from_mut(task);
                                reads += read_regions(&mut reader, tid, one, &mut record, checker, &mut scratch, bam_path)?;
                            }
                        } else {
                            reads += read_regions(&mut reader, tid, &mut item, &mut record, checker, &mut scratch, bam_path)?;
                        }
                    }
                    Ok(reads)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut reads = 0;
    for r in results {
        reads += r?;
    }
    if config.perf_debug {
        crate::info!(
            "PerfStats: Read: {} regions, {} reads, {:.3} seconds",
            task_count,
            reads,
            started.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::AmberBase;

    #[test]
    fn tasks_split_on_gaps() {
        let mut lists: ByChr<Evidence> = ByChr::new();
        for p in [100, 200, 5000, 9001, 20000] {
            lists.put(1, Evidence::new("chr1", p, AmberBase::A, AmberBase::C));
        }
        lists.put(3, Evidence::new("chr3", 7, AmberBase::A, AmberBase::C));
        let t = region_tasks(&lists, 4000);
        // 200 + 4000 < 5000 starts a new region, 5000 + 4000 < 9001 another
        assert_eq!(t, vec![(1, 0, 2, 100, 200), (1, 2, 3, 5000, 5000), (1, 3, 4, 9001, 9001), (1, 4, 5, 20000, 20000), (3, 0, 1, 7, 7)]);
        assert_eq!(region_tasks(&lists, 5000).len(), 3);
    }

    #[test]
    fn groups_of_tasks() {
        // (chromosome, first, last, start, end)
        let t = vec![(1, 0, 2, 100, 200), (1, 2, 3, 5000, 5000), (1, 3, 4, 9001, 9001), (1, 4, 5, 20000, 20000), (3, 0, 1, 7, 7)];
        // a new item at the chromosome change only
        assert_eq!(group_tasks(&t, 256_000, 10_000_000), vec![0..4, 4..5]);
        // gap 0: one item per task (Java's reading)
        assert_eq!(group_tasks(&t, 0, 10_000_000), vec![0..1, 1..2, 2..3, 3..4, 4..5]);
        // the gap from 9001 to 20000 (10999) is above 10000
        assert_eq!(group_tasks(&t, 10_000, 10_000_000), vec![0..3, 3..4, 4..5]);
        // an item spans at most 8000 bases: 100-5000 (4901) fits, 100-9001 (8902) does not
        assert_eq!(group_tasks(&t, 256_000, 8_000), vec![0..2, 2..3, 3..4, 4..5]);
        assert!(group_tasks(&[], 1, 1).is_empty());
    }
}
