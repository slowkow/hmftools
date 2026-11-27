//! Reading a BAM from its start in file order, with the decompression spread over all threads.
//!
//! htslib's multi-threaded BGZF reader leaves the reading thread to wait while its pool inflates; with two
//! threads, one inflates and the other mostly waits. Here the file is cut into batches of whole BGZF
//! blocks (about 1 MB compressed). Any thread claims the next batch, inflates it (libdeflate) and checks
//! its CRC32s; the main thread consumes the batches in file order and, rather than wait, inflates a batch
//! itself when the next one is not ready. Records are parsed in place from the inflated bytes; one that
//! crosses a batch boundary is put together in a small buffer.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::BTreeMap;
use std::sync::{Arc, Condvar, Mutex};

/// Compressed bytes per batch
const BATCH_BYTES: usize = 1 << 20;

/// Read-ahead of the file (not Java options; the environment variables TEAL_READ_STREAMS and
/// TEAL_READ_CHUNK_MB set them): the number of file handles reading ahead, and the bytes each reads
/// in one go.
const READ_STREAMS: usize = 4;
const READ_CHUNK_MB: usize = 32;

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name).ok().and_then(|v| v.parse::<usize>().ok()).filter(|&v| v > 0).unwrap_or(default)
}

struct ChunkState {
    done: BTreeMap<u64, std::io::Result<Vec<u8>>>,
    claimed: u64,
    consumed: u64,
    stop: bool,
    free: Vec<Vec<u8>>,
}

struct ChunkShared {
    state: Mutex<ChunkState>,
    cv: Condvar,
}

/// The file's bytes from an offset to its end, in order, read ahead by several threads, each with its own
/// file handle and each chunk one sequential read. Through a bucket mount (gcsfuse) one handle reading in
/// order gives one stream's throughput (about 160 MB/s on Google Batch); the chunks of several handles
/// add up, while the bytes still reach the BGZF parser in file order, so nothing else changes.
struct ChunkSource {
    shared: Arc<ChunkShared>,
    threads: Vec<std::thread::JoinHandle<()>>,
    cur: Vec<u8>,
    pos: usize,
    next: u64,
    chunks: u64,
}

impl ChunkSource {
    fn open(path: &str, start: u64) -> Result<ChunkSource> {
        let len = std::fs::metadata(path).with_context(|| format!("opening {path}"))?.len();
        let streams = env_usize("TEAL_READ_STREAMS", READ_STREAMS);
        let chunk = (env_usize("TEAL_READ_CHUNK_MB", READ_CHUNK_MB) << 20) as u64;
        let chunks = len.saturating_sub(start).div_ceil(chunk);
        let max_ahead = streams as u64 + 2;
        let shared = Arc::new(ChunkShared {
            state: Mutex::new(ChunkState { done: BTreeMap::new(), claimed: 0, consumed: 0, stop: false, free: Vec::new() }),
            cv: Condvar::new(),
        });
        let mut threads = Vec::new();
        for _ in 0..streams.min(chunks.max(1) as usize) {
            let mut file = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
            let s = shared.clone();
            threads.push(std::thread::spawn(move || {
                use std::io::{Read, Seek};
                let mut at: Option<u64> = None;
                loop {
                    let (k, mut buf) = {
                        let mut st = s.state.lock().unwrap();
                        while !st.stop && st.claimed < chunks && st.claimed >= st.consumed + max_ahead {
                            st = s.cv.wait(st).unwrap();
                        }
                        if st.stop || st.claimed >= chunks {
                            return;
                        }
                        let k = st.claimed;
                        st.claimed += 1;
                        (k, st.free.pop().unwrap_or_default())
                    };
                    let off = start + k * chunk;
                    let n = chunk.min(len - off) as usize;
                    buf.resize(n, 0);
                    let r = (|| {
                        if at != Some(off) {
                            file.seek(std::io::SeekFrom::Start(off))?;
                        }
                        file.read_exact(&mut buf)?;
                        Ok(())
                    })();
                    at = if r.is_ok() { Some(off + n as u64) } else { None };
                    let mut st = s.state.lock().unwrap();
                    st.done.insert(k, r.map(|_| buf));
                    drop(st);
                    s.cv.notify_all();
                }
            }));
        }
        Ok(ChunkSource { shared, threads, cur: Vec::new(), pos: 0, next: 0, chunks })
    }

    /// Up to `buf.len()` bytes into `buf`; 0 at the end of the file
    fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        if self.pos == self.cur.len() {
            if self.next == self.chunks {
                return Ok(0);
            }
            let chunk = {
                let mut st = self.shared.state.lock().unwrap();
                loop {
                    if let Some(c) = st.done.remove(&self.next) {
                        break c;
                    }
                    st = self.shared.cv.wait(st).unwrap();
                }
            };
            let chunk = chunk.context("reading the BAM")?;
            let old = std::mem::replace(&mut self.cur, chunk);
            self.pos = 0;
            self.next += 1;
            let mut st = self.shared.state.lock().unwrap();
            st.consumed = self.next;
            st.free.push(old);
            drop(st);
            self.shared.cv.notify_all();
        }
        let n = buf.len().min(self.cur.len() - self.pos);
        buf[..n].copy_from_slice(&self.cur[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl Drop for ChunkSource {
    fn drop(&mut self) {
        self.shared.state.lock().unwrap().stop = true;
        self.shared.cv.notify_all();
        for t in self.threads.drain(..) {
            let _ = t.join();
        }
    }
}

/// The whole BGZF blocks of the file, a batch at a time
struct BlockReader {
    file: ChunkSource,
    /// emptied compressed buffers, for reuse
    free: Vec<Vec<u8>>,
    /// bytes read but not yet part of a whole block
    carry: Vec<u8>,
    eof: bool,
    next_id: u64,
}

impl BlockReader {
    /// The next batch: its id and its compressed blocks; None at the end of the file
    fn next_batch(&mut self) -> Result<Option<(u64, Vec<u8>)>> {
        let mut buf = self.free.pop().unwrap_or_default();
        buf.clear();
        buf.reserve(BATCH_BYTES + self.carry.len() + (64 << 10));
        buf.extend_from_slice(&self.carry);
        self.carry.clear();
        while buf.len() < BATCH_BYTES && !self.eof {
            let old = buf.len();
            buf.resize(old + BATCH_BYTES, 0);
            let n = self.file.read(&mut buf[old..])?;
            buf.truncate(old + n);
            if n == 0 {
                self.eof = true;
            }
        }
        if buf.is_empty() {
            return Ok(None);
        }
        // whole blocks only
        let mut end = 0usize;
        while end + 18 <= buf.len() {
            let size = block_size(&buf[end..])?;
            if end + size > buf.len() {
                break;
            }
            end += size;
        }
        if end == 0 {
            if self.eof {
                bail!("truncated BGZF block at the end of the file");
            }
            // a block larger than the batch: read on
            self.carry = buf;
            return self.next_batch();
        }
        self.carry.extend_from_slice(&buf[end..]);
        buf.truncate(end);
        let id = self.next_id;
        self.next_id += 1;
        Ok(Some((id, buf)))
    }
}

/// The size of the BGZF block at the start of `b` (from its BC extra field)
fn block_size(b: &[u8]) -> Result<usize> {
    if b.len() < 18 || b[0] != 31 || b[1] != 139 || b[2] != 8 || b[3] & 4 == 0 {
        bail!("not a BGZF block (invalid gzip header)");
    }
    let xlen = u16::from_le_bytes([b[10], b[11]]) as usize;
    let mut p = 12;
    while p + 4 <= 12 + xlen && p + 4 <= b.len() {
        let slen = u16::from_le_bytes([b[p + 2], b[p + 3]]) as usize;
        if b[p] == b'B' && b[p + 1] == b'C' && slen == 2 && p + 6 <= b.len() {
            return Ok(u16::from_le_bytes([b[p + 4], b[p + 5]]) as usize + 1);
        }
        p += 4 + slen;
    }
    bail!("not a BGZF block (no BC field)")
}

struct Inflater(*mut libdeflate_sys::libdeflate_decompressor);

impl Inflater {
    fn new() -> Inflater {
        Inflater(unsafe { libdeflate_sys::libdeflate_alloc_decompressor() })
    }
    /// Inflates every block of a batch into `out`, checking sizes and CRC32s
    fn inflate_batch(&mut self, comp: &[u8], out: &mut Vec<u8>) -> Result<()> {
        // the inflated size of the batch, from the blocks' ISIZE fields
        let mut total = 0usize;
        let mut p = 0usize;
        while p < comp.len() {
            let size = block_size(&comp[p..])?;
            if size < 26 || p + size > comp.len() {
                bail!("corrupt BGZF block (size)");
            }
            total += u32::from_le_bytes(comp[p + size - 4..p + size].try_into().unwrap()) as usize;
            p += size;
        }
        out.clear();
        out.reserve(total);
        let dst = out.spare_capacity_mut();
        let mut at = 0usize;
        let mut p = 0usize;
        while p < comp.len() {
            let size = block_size(&comp[p..])?;
            let block = &comp[p..p + size];
            let xlen = u16::from_le_bytes([block[10], block[11]]) as usize;
            let cdata = &block[12 + xlen..size - 8];
            let crc = u32::from_le_bytes(block[size - 8..size - 4].try_into().unwrap());
            let isize = u32::from_le_bytes(block[size - 4..].try_into().unwrap()) as usize;
            let target = dst[at..at + isize].as_mut_ptr() as *mut u8;
            if isize > 0 {
                let mut actual = 0usize;
                let r = unsafe {
                    libdeflate_sys::libdeflate_deflate_decompress(self.0, cdata.as_ptr() as *const _, cdata.len(), target as *mut _, isize, &mut actual)
                };
                if r != libdeflate_sys::libdeflate_result_LIBDEFLATE_SUCCESS || actual != isize {
                    bail!("corrupt BGZF block (inflate failed)");
                }
            }
            // the bytes are written now
            let c = unsafe { libdeflate_sys::libdeflate_crc32(0, target as *const _, isize) };
            if c != crc {
                bail!("corrupt BGZF block (CRC32 mismatch)");
            }
            at += isize;
            p += size;
        }
        // every byte up to `total` was written by libdeflate
        unsafe { out.set_len(total) };
        Ok(())
    }
}

impl Drop for Inflater {
    fn drop(&mut self) {
        unsafe { libdeflate_sys::libdeflate_free_decompressor(self.0) }
    }
}

// the decompressor is used by one thread at a time
unsafe impl Send for Inflater {}

struct State {
    done: BTreeMap<u64, Result<Vec<u8>, String>>,
    /// batches handed out so far, and consumed
    issued: u64,
    consumed: u64,
    eof: bool,
    stop: bool,
}

struct Shared {
    reader: Mutex<BlockReader>,
    /// emptied inflated buffers, for reuse
    free: Mutex<Vec<Vec<u8>>>,
    state: Mutex<State>,
    cv: Condvar,
    max_ahead: u64,
}

impl Shared {
    /// Claims, inflates and stores the next batch; false when there is none to claim now
    fn work(&self, inflater: &mut Inflater) -> bool {
        {
            let st = self.state.lock().unwrap();
            if st.stop || st.eof || st.issued - st.consumed >= self.max_ahead {
                return false;
            }
        }
        let batch = {
            let mut r = self.reader.lock().unwrap();
            let b = r.next_batch();
            if let Ok(Some((id, _))) = &b {
                let mut st = self.state.lock().unwrap();
                st.issued = st.issued.max(id + 1);
            }
            b
        };
        match batch {
            Ok(None) => {
                self.state.lock().unwrap().eof = true;
                self.cv.notify_all();
                false
            }
            Err(e) => {
                let mut st = self.state.lock().unwrap();
                let id = st.issued;
                st.issued += 1;
                st.done.insert(id, Err(format!("{e:#}")));
                st.eof = true;
                drop(st);
                self.cv.notify_all();
                false
            }
            Ok(Some((id, comp))) => {
                let mut out = self.free.lock().unwrap().pop().unwrap_or_default();
                let r = inflater.inflate_batch(&comp, &mut out).map(|_| out).map_err(|e| format!("{e:#}"));
                self.reader.lock().unwrap().free.push(comp);
                self.state.lock().unwrap().done.insert(id, r);
                self.cv.notify_all();
                true
            }
        }
    }
}

/// The inflated bytes of a BGZF file, in order
pub struct BgzfStream {
    shared: Arc<Shared>,
    workers: Vec<std::thread::JoinHandle<()>>,
    inflater: Inflater,
    cur: Vec<u8>,
    pos: usize,
    tmp: Vec<u8>,
    finished: bool,
}

impl BgzfStream {
    /// `threads` threads in all inflate: this one and threads - 1 helpers
    pub fn open(path: &str, threads: usize) -> Result<BgzfStream> {
        Self::open_at(path, 0, threads)
    }

    /// From a virtual offset (the compressed offset of a block << 16 | the offset within it)
    pub fn open_at(path: &str, voffset: u64, threads: usize) -> Result<BgzfStream> {
        let file = ChunkSource::open(path, voffset >> 16)?;
        let threads = threads.max(1);
        let shared = Arc::new(Shared {
            reader: Mutex::new(BlockReader { file, free: Vec::new(), carry: Vec::new(), eof: false, next_id: 0 }),
            free: Mutex::new(Vec::new()),
            state: Mutex::new(State { done: BTreeMap::new(), issued: 0, consumed: 0, eof: false, stop: false }),
            cv: Condvar::new(),
            max_ahead: 2 * threads as u64 + 2,
        });
        let workers = (1..threads)
            .map(|_| {
                let s = shared.clone();
                std::thread::spawn(move || {
                    let mut inflater = Inflater::new();
                    loop {
                        if s.work(&mut inflater) {
                            continue;
                        }
                        let st = s.state.lock().unwrap();
                        if st.stop || st.eof {
                            break;
                        }
                        // wait for the reader to consume
                        if st.issued - st.consumed >= s.max_ahead {
                            let _st = s.cv.wait(st).unwrap();
                        }
                    }
                })
            })
            .collect();
        let mut stream = BgzfStream { shared, workers, inflater: Inflater::new(), cur: Vec::new(), pos: 0, tmp: Vec::new(), finished: false };
        let within = (voffset & 0xFFFF) as usize;
        if within > 0 {
            stream.get(within)?.ok_or_else(|| anyhow!("virtual offset past the end of {path}"))?;
        }
        Ok(stream)
    }

    /// The next batch in order, None at the end
    fn next_batch(&mut self) -> Result<Option<Vec<u8>>> {
        if self.finished {
            return Ok(None);
        }
        loop {
            {
                let mut st = self.shared.state.lock().unwrap();
                let want = st.consumed;
                if let Some(b) = st.done.remove(&want) {
                    st.consumed += 1;
                    drop(st);
                    self.shared.cv.notify_all();
                    return b.map(Some).map_err(|e| anyhow!(e));
                }
                if st.eof && st.issued == st.consumed {
                    self.finished = true;
                    return Ok(None);
                }
            }
            // help rather than wait
            if !self.shared.work(&mut self.inflater) {
                let st = self.shared.state.lock().unwrap();
                let want = st.consumed;
                if !st.done.contains_key(&want) && !(st.eof && st.issued == st.consumed) {
                    let _st = self.shared.cv.wait(st).unwrap();
                }
            }
        }
    }

    /// The next `n` bytes, contiguous; None at the end of the stream (an error when it ends within them)
    pub fn get(&mut self, n: usize) -> Result<Option<&[u8]>> {
        if self.cur.len() - self.pos >= n {
            let s = self.pos;
            self.pos += n;
            return Ok(Some(&self.cur[s..s + n]));
        }
        self.tmp.clear();
        self.tmp.extend_from_slice(&self.cur[self.pos..]);
        self.pos = self.cur.len();
        while self.tmp.len() < n {
            match self.next_batch()? {
                Some(b) => {
                    let old = std::mem::replace(&mut self.cur, b);
                    self.shared.free.lock().unwrap().push(old);
                    let take = (n - self.tmp.len()).min(self.cur.len());
                    self.tmp.extend_from_slice(&self.cur[..take]);
                    self.pos = take;
                }
                None => {
                    if self.tmp.is_empty() {
                        return Ok(None);
                    }
                    bail!("truncated BAM: the file ends within a record");
                }
            }
        }
        Ok(Some(&self.tmp))
    }

    fn get_i32(&mut self) -> Result<Option<i32>> {
        Ok(self.get(4)?.map(|b| i32::from_le_bytes(b.try_into().unwrap())))
    }

    /// Skips the BAM header (magic, text and reference list)
    pub fn skip_header(&mut self) -> Result<()> {
        let magic = self.get(4)?.ok_or_else(|| anyhow!("empty BAM"))?;
        if magic != b"BAM\x01" {
            bail!("not a BAM file");
        }
        let l_text = self.get_i32()?.ok_or_else(|| anyhow!("truncated BAM header"))?;
        self.get(l_text as usize)?.ok_or_else(|| anyhow!("truncated BAM header"))?;
        let n_ref = self.get_i32()?.ok_or_else(|| anyhow!("truncated BAM header"))?;
        for _ in 0..n_ref {
            let l_name = self.get_i32()?.ok_or_else(|| anyhow!("truncated BAM header"))?;
            self.get(l_name as usize + 4)?.ok_or_else(|| anyhow!("truncated BAM header"))?;
        }
        Ok(())
    }

    /// The next record (its bytes after block_size), None at the end
    pub fn next_record(&mut self) -> Result<Option<&[u8]>> {
        let Some(size) = self.get_i32()? else { return Ok(None) };
        if size < 32 {
            bail!("invalid BAM record (block size {size})");
        }
        match self.get(size as usize)? {
            Some(r) => Ok(Some(r)),
            None => bail!("truncated BAM: the file ends within a record"),
        }
    }
}

impl Drop for BgzfStream {
    fn drop(&mut self) {
        self.shared.state.lock().unwrap().stop = true;
        self.shared.cv.notify_all();
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
    }
}

#[inline]
pub fn le_i32(b: &[u8], at: usize) -> i32 {
    i32::from_le_bytes(b[at..at + 4].try_into().unwrap())
}

#[inline]
pub fn le_u16(b: &[u8], at: usize) -> u16 {
    u16::from_le_bytes([b[at], b[at + 1]])
}

/// The parts of a raw record the test of every record needs
pub struct RawView<'a> {
    pub tid: i32,
    pub pos: i32,
    pub name: &'a [u8],
    pub seq: &'a [u8],
    pub seq_len: usize,
}

#[inline]
pub fn view(raw: &[u8]) -> Result<RawView<'_>> {
    let l_read_name = raw[8] as usize;
    let n_cigar = le_u16(raw, 12) as usize;
    let l_seq = le_i32(raw, 16).max(0) as usize;
    let seq_off = 32 + l_read_name + 4 * n_cigar;
    let seq_end = seq_off + l_seq.div_ceil(2);
    if l_read_name == 0 || seq_end + l_seq > raw.len() {
        bail!("invalid BAM record");
    }
    Ok(RawView { tid: le_i32(raw, 0), pos: le_i32(raw, 4), name: &raw[32..32 + l_read_name - 1], seq: &raw[seq_off..seq_end], seq_len: l_seq })
}

/// A rust-htslib record of the raw bytes, as htslib's bam_read1 sets it up (the name padded with NULs to
/// align the CIGAR)
pub fn record_from_raw(raw: &[u8], rec: &mut rust_htslib::bam::Record) {
    let l_read_name = raw[8] as usize;
    let extranul = (4 - l_read_name % 4) % 4;
    let mut data = Vec::with_capacity(raw.len() - 32 + extranul);
    data.extend_from_slice(&raw[32..32 + l_read_name]);
    data.extend(std::iter::repeat_n(0u8, extranul));
    data.extend_from_slice(&raw[32 + l_read_name..]);
    rec.set_data(&data);
    let c = &mut rec.inner_mut().core;
    c.tid = le_i32(raw, 0);
    c.pos = le_i32(raw, 4) as i64;
    c.l_qname = (l_read_name + extranul) as u16;
    c.l_extranul = extranul as u8;
    c.qual = raw[9];
    c.bin = le_u16(raw, 10);
    c.n_cigar = le_u16(raw, 12) as u32;
    c.flag = le_u16(raw, 14);
    c.l_qseq = le_i32(raw, 16);
    c.mtid = le_i32(raw, 20);
    c.mpos = le_i32(raw, 24) as i64;
    c.isize_ = le_i32(raw, 28) as i64;
}

/// The virtual offset where a coordinate-sorted BAM's unplaced reads start, from its BAI index: after the
/// last record of every reference (the end offset in each reference's pseudo-bin), as htslib's
/// HTS_IDX_NOCOOR iterator finds it. None without a BAI (a CSI index, for one) or when it cannot be read.
pub fn unplaced_start(bam: &str) -> Option<u64> {
    let candidates = [format!("{bam}.bai"), bam.strip_suffix(".bam").map(|b| format!("{b}.bai")).unwrap_or_default()];
    let path = candidates.iter().find(|p| !p.is_empty() && std::path::Path::new(p).exists())?;
    let b = std::fs::read(path).ok()?;
    if b.len() < 8 || &b[..4] != b"BAI\x01" {
        return None;
    }
    let u32_at = |p: usize| -> Option<u32> { b.get(p..p + 4).map(|x| u32::from_le_bytes(x.try_into().unwrap())) };
    let u64_at = |p: usize| -> Option<u64> { b.get(p..p + 8).map(|x| u64::from_le_bytes(x.try_into().unwrap())) };
    let n_ref = u32_at(4)? as usize;
    let mut p = 8usize;
    let mut start: Option<u64> = None;
    for _ in 0..n_ref {
        let n_bin = u32_at(p)? as usize;
        p += 4;
        for _ in 0..n_bin {
            let bin = u32_at(p)?;
            let n_chunk = u32_at(p + 4)? as usize;
            if bin == 37450 && n_chunk >= 1 {
                // the reference's first and last record offsets
                let end = u64_at(p + 8 + 8)?;
                start = Some(start.map_or(end, |s| s.max(end)));
            }
            p += 8 + 16 * n_chunk;
        }
        let n_intv = u32_at(p)? as usize;
        p += 4 + 8 * n_intv;
    }
    // no mapped reads: from the first record (the caller skips the header)
    Some(start.unwrap_or(0))
}
