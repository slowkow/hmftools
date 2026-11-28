//! BGZF files through libdeflate: a reader that inflates blocks on worker threads and hands out the
//! uncompressed stream in order, and a writer that deflates full blocks on worker threads, writes them
//! in order and knows the virtual offset of every byte it was given (for the BAM index).

use anyhow::{Context, Result, bail};
use libdeflate_sys::*;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

/// Uncompressed bytes per block written (htslib's BGZF_BLOCK_SIZE): even stored without compression a
/// block of this size fits the 64 KiB limit of a BGZF block.
pub const BLOCK_DATA: usize = 0xff00;

/// The empty block that ends a BGZF file.
pub const EOF_BLOCK: [u8; 28] = [31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0];

const HEADER_LEN: usize = 18;
const MAX_BLOCK: usize = 65536;

pub fn crc32(data: &[u8]) -> u32 {
    unsafe { libdeflate_crc32(0, data.as_ptr() as *const _, data.len()) }
}

struct Inflater(*mut libdeflate_decompressor);
unsafe impl Send for Inflater {}

impl Inflater {
    fn new() -> Self {
        let p = unsafe { libdeflate_alloc_decompressor() };
        assert!(!p.is_null(), "libdeflate: out of memory");
        Inflater(p)
    }
    /// Inflates raw deflate data into exactly `out.len()` bytes.
    fn inflate(&mut self, cdata: &[u8], out: &mut [u8]) -> Result<()> {
        let mut n = 0usize;
        let r = unsafe {
            libdeflate_deflate_decompress(self.0, cdata.as_ptr() as *const _, cdata.len(), out.as_mut_ptr() as *mut _, out.len(), &mut n)
        };
        if r != libdeflate_result_LIBDEFLATE_SUCCESS || n != out.len() {
            bail!("corrupt BGZF block (inflate failed)");
        }
        Ok(())
    }
}

impl Drop for Inflater {
    fn drop(&mut self) {
        unsafe { libdeflate_free_decompressor(self.0) }
    }
}

struct Deflater {
    p: *mut libdeflate_compressor,
    stored: *mut libdeflate_compressor,
}
unsafe impl Send for Deflater {}

impl Deflater {
    fn new(level: i32) -> Self {
        let level = level.clamp(0, 12);
        let p = unsafe { libdeflate_alloc_compressor(level) };
        let stored = unsafe { libdeflate_alloc_compressor(0) };
        assert!(!p.is_null() && !stored.is_null(), "libdeflate: out of memory");
        Deflater { p, stored }
    }
    /// Appends one whole BGZF block holding `data` (at most BLOCK_DATA bytes) to `out`.
    fn block(&mut self, data: &[u8], out: &mut Vec<u8>) {
        debug_assert!(data.len() <= BLOCK_DATA);
        let start = out.len();
        out.resize(start + MAX_BLOCK, 0);
        let avail = MAX_BLOCK - HEADER_LEN - 8;
        let dst = &mut out[start + HEADER_LEN..start + HEADER_LEN + avail];
        let mut n = unsafe { libdeflate_deflate_compress(self.p, data.as_ptr() as *const _, data.len(), dst.as_mut_ptr() as *mut _, avail) };
        if n == 0 {
            // does not fit compressed: store it
            n = unsafe { libdeflate_deflate_compress(self.stored, data.as_ptr() as *const _, data.len(), dst.as_mut_ptr() as *mut _, avail) };
            assert!(n > 0, "a stored BGZF block must fit");
        }
        let total = HEADER_LEN + n + 8;
        let h = &mut out[start..start + HEADER_LEN];
        h.copy_from_slice(&[31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0]);
        h[16..18].copy_from_slice(&((total - 1) as u16).to_le_bytes());
        let crc = crc32(data);
        let t = start + HEADER_LEN + n;
        out[t..t + 4].copy_from_slice(&crc.to_le_bytes());
        out[t + 4..t + 8].copy_from_slice(&(data.len() as u32).to_le_bytes());
        out.truncate(start + total);
    }
}

impl Drop for Deflater {
    fn drop(&mut self) {
        unsafe {
            libdeflate_free_compressor(self.p);
            libdeflate_free_compressor(self.stored);
        }
    }
}

/// Reads the next BGZF block's bytes (header included) from `r` into `buf`; false at the end of the file.
fn read_raw_block(r: &mut impl Read, buf: &mut Vec<u8>, path: &str) -> Result<bool> {
    let mut h = [0u8; HEADER_LEN];
    let mut got = 0;
    while got < HEADER_LEN {
        let n = r.read(&mut h[got..])?;
        if n == 0 {
            if got == 0 {
                return Ok(false);
            }
            bail!("{path}: truncated BGZF block header");
        }
        got += n;
    }
    if h[0] != 31 || h[1] != 139 || h[2] != 8 || h[3] & 4 == 0 {
        bail!("{path}: not a BGZF file (BAM input is required)");
    }
    let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
    // the BC subfield is normally the only one; read the extra field to find it
    let mut extra = vec![0u8; xlen];
    extra[..xlen.min(6)].copy_from_slice(&h[12..12 + xlen.min(6)]);
    if xlen > 6 {
        r.read_exact(&mut extra[6..]).with_context(|| format!("{path}: truncated BGZF block"))?;
    }
    let mut bsize = None;
    let mut i = 0;
    while i + 4 <= xlen {
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if extra[i] == 66 && extra[i + 1] == 67 && slen == 2 && i + 6 <= xlen {
            bsize = Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize + 1);
        }
        i += 4 + slen;
    }
    let Some(bsize) = bsize else { bail!("{path}: BGZF block without its size (BC) field") };
    let header_len = 12 + xlen;
    if bsize < header_len + 8 {
        bail!("{path}: bad BGZF block size");
    }
    buf.clear();
    buf.extend_from_slice(&h[..12]);
    buf.extend_from_slice(&extra);
    let rest = bsize - header_len;
    let at = buf.len();
    buf.resize(at + rest, 0);
    r.read_exact(&mut buf[at..]).with_context(|| format!("{path}: truncated BGZF block"))?;
    Ok(true)
}

/// Inflates one raw block (as read by `read_raw_block`), appending its data to `out`.
fn inflate_block(inf: &mut Inflater, raw: &[u8], out: &mut Vec<u8>, path: &str) -> Result<()> {
    let xlen = u16::from_le_bytes([raw[10], raw[11]]) as usize;
    let header_len = 12 + xlen;
    let n = raw.len();
    let cdata = &raw[header_len..n - 8];
    let crc = u32::from_le_bytes(raw[n - 8..n - 4].try_into().unwrap());
    let isize = u32::from_le_bytes(raw[n - 4..].try_into().unwrap()) as usize;
    if isize > MAX_BLOCK {
        bail!("{path}: bad BGZF block (uncompressed size {isize})");
    }
    let at = out.len();
    out.resize(at + isize, 0);
    inf.inflate(cdata, &mut out[at..]).with_context(|| path.to_string())?;
    if crc32(&out[at..]) != crc {
        bail!("{path}: BGZF block CRC mismatch (corrupt file)");
    }
    Ok(())
}

/// Compressed bytes read per job when inflating on threads.
const JOB_COMPRESSED: usize = 1 << 20;

type ChunkResult = (u64, Result<Vec<u8>>);

enum Source {
    Inline { file: BufReader<File>, inf: Inflater, raw: Vec<u8>, done: bool },
    Threads { rx: Receiver<ChunkResult>, pending: BTreeMap<u64, Result<Vec<u8>>>, next: u64, handles: Vec<JoinHandle<()>>, done: bool },
}

/// The uncompressed stream of a BGZF file.
pub struct BgzfReader {
    path: String,
    src: Source,
    buf: Vec<u8>,
    pos: usize,
}

impl BgzfReader {
    /// Opens `path`; with `threads` > 0 the blocks are inflated on that many worker threads (and one
    /// thread reads the file).
    pub fn open(path: &str, threads: usize) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("cannot open {path}"))?;
        let src = if threads == 0 {
            Source::Inline { file: BufReader::with_capacity(1 << 20, file), inf: Inflater::new(), raw: Vec::new(), done: false }
        } else {
            let (job_tx, job_rx) = sync_channel::<(u64, Vec<Vec<u8>>)>(threads * 2);
            let (res_tx, res_rx) = sync_channel::<ChunkResult>(threads * 2 + 2);
            let mut handles = Vec::new();
            let p = path.to_string();
            let reader_res_tx = res_tx.clone();
            handles.push(std::thread::spawn(move || {
                let mut file = BufReader::with_capacity(1 << 20, file);
                let mut seq = 0u64;
                loop {
                    let mut blocks = Vec::new();
                    let mut size = 0;
                    let mut end = false;
                    while size < JOB_COMPRESSED {
                        let mut raw = Vec::new();
                        match read_raw_block(&mut file, &mut raw, &p) {
                            Ok(true) => {
                                size += raw.len();
                                blocks.push(raw);
                            }
                            Ok(false) => {
                                end = true;
                                break;
                            }
                            Err(e) => {
                                let _ = reader_res_tx.send((seq, Err(e)));
                                return;
                            }
                        }
                    }
                    if !blocks.is_empty() && job_tx.send((seq, blocks)).is_err() {
                        return;
                    }
                    seq += 1;
                    if end {
                        return;
                    }
                }
            }));
            let job_rx = Arc::new(Mutex::new(job_rx));
            for _ in 0..threads {
                let job_rx = job_rx.clone();
                let res_tx = res_tx.clone();
                let p = path.to_string();
                handles.push(std::thread::spawn(move || {
                    let mut inf = Inflater::new();
                    loop {
                        let job = { job_rx.lock().unwrap().recv() };
                        let Ok((seq, blocks)) = job else { return };
                        let total: usize = blocks.iter().map(|b| u32::from_le_bytes(b[b.len() - 4..].try_into().unwrap()) as usize).sum();
                        let mut out = Vec::with_capacity(total);
                        let mut r = Ok(());
                        for b in &blocks {
                            r = inflate_block(&mut inf, b, &mut out, &p);
                            if r.is_err() {
                                break;
                            }
                        }
                        if res_tx.send((seq, r.map(|_| out))).is_err() {
                            return;
                        }
                    }
                }));
            }
            drop(res_tx);
            Source::Threads { rx: res_rx, pending: BTreeMap::new(), next: 0, handles, done: false }
        };
        Ok(BgzfReader { path: path.to_string(), src, buf: Vec::new(), pos: 0 })
    }

    /// Replaces the current buffer with the next stretch of the stream; false at its end.
    fn refill(&mut self) -> Result<bool> {
        self.buf.clear();
        self.pos = 0;
        match &mut self.src {
            Source::Inline { file, inf, raw, done } => {
                while !*done && self.buf.len() < (1 << 20) {
                    if !read_raw_block(file, raw, &self.path)? {
                        *done = true;
                        break;
                    }
                    inflate_block(inf, raw, &mut self.buf, &self.path)?;
                }
                Ok(!self.buf.is_empty())
            }
            Source::Threads { rx, pending, next, done, .. } => loop {
                if let Some(r) = pending.remove(next) {
                    *next += 1;
                    self.buf = r?;
                    if self.buf.is_empty() {
                        continue;
                    }
                    return Ok(true);
                }
                if *done {
                    if let Some((_, r)) = pending.pop_first() {
                        r?;
                    }
                    return Ok(false);
                }
                match rx.recv() {
                    Ok((seq, r)) => {
                        if let Err(e) = r {
                            return Err(e);
                        }
                        pending.insert(seq, r);
                    }
                    Err(_) => *done = true,
                }
            },
        }
    }

    /// Fills `out` from the stream; false if the stream ended before the first byte (true end of
    /// data), an error if it ended part way.
    pub fn read_exact_or_eof(&mut self, out: &mut [u8]) -> Result<bool> {
        let mut got = 0;
        while got < out.len() {
            if self.pos == self.buf.len() && !self.refill()? {
                if got == 0 {
                    return Ok(false);
                }
                bail!("{}: unexpected end of the BAM data", self.path);
            }
            let n = (out.len() - got).min(self.buf.len() - self.pos);
            out[got..got + n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
            self.pos += n;
            got += n;
        }
        Ok(true)
    }

    pub fn read_exact(&mut self, out: &mut [u8]) -> Result<()> {
        if !self.read_exact_or_eof(out)? {
            bail!("{}: unexpected end of the BAM data", self.path);
        }
        Ok(())
    }

    /// Reads the next length-prefixed BAM record into `rec` (without its 4-byte length); false at the
    /// end of the file.
    pub fn read_record(&mut self, rec: &mut Vec<u8>) -> Result<bool> {
        let mut len = [0u8; 4];
        if !self.read_exact_or_eof(&mut len)? {
            return Ok(false);
        }
        let n = u32::from_le_bytes(len) as usize;
        if n < 32 {
            bail!("{}: bad BAM record length {n}", self.path);
        }
        rec.clear();
        // fast path: the record lies within the current buffer
        if self.buf.len() - self.pos >= n {
            rec.extend_from_slice(&self.buf[self.pos..self.pos + n]);
            self.pos += n;
            return Ok(true);
        }
        rec.resize(n, 0);
        self.read_exact(rec)?;
        Ok(true)
    }

    pub fn path(&self) -> &str {
        &self.path
    }
}

impl Drop for BgzfReader {
    fn drop(&mut self) {
        if let Source::Threads { rx, handles, .. } = &mut self.src {
            // unblock the workers: drop the receiver, then wait for them
            let (_tx, dummy) = sync_channel(1);
            let old = std::mem::replace(rx, dummy);
            drop(old);
            for h in handles.drain(..) {
                let _ = h.join();
            }
        }
    }
}

type BlockJob = (u64, Vec<u8>);
type BlockDone = (u64, Vec<u8>, Vec<u8>);

struct Pool {
    job_tx: Option<SyncSender<BlockJob>>,
    res_rx: Receiver<BlockDone>,
    handles: Vec<JoinHandle<()>>,
    workers: usize,
}

/// A BGZF file being written.
pub struct BgzfWriter {
    path: String,
    out: BufWriter<File>,
    deflater: Deflater,
    pool: Option<Pool>,
    block: Vec<u8>,
    spare: Vec<Vec<u8>>,
    /// uncompressed bytes taken so far
    taken: u64,
    next_seq: u64,
    written_seq: u64,
    pending: BTreeMap<u64, (Vec<u8>, Vec<u8>)>,
    /// file address of every block written, and of the next one
    block_addr: Vec<u64>,
    file_pos: u64,
    compressed: Vec<u8>,
}

impl BgzfWriter {
    /// Creates `path`; `level` is the deflate level (0 to 12), `threads` the number of worker threads
    /// (0: deflate on the calling thread).
    pub fn create(path: &str, level: i32, threads: usize) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
        let pool = if threads == 0 {
            None
        } else {
            let (job_tx, job_rx) = sync_channel::<BlockJob>(threads * 4);
            let (res_tx, res_rx) = sync_channel::<BlockDone>(threads * 4 + 4);
            let job_rx = Arc::new(Mutex::new(job_rx));
            let mut handles = Vec::new();
            for _ in 0..threads {
                let job_rx = job_rx.clone();
                let res_tx = res_tx.clone();
                handles.push(std::thread::spawn(move || {
                    let mut d = Deflater::new(level);
                    loop {
                        let job = { job_rx.lock().unwrap().recv() };
                        let Ok((seq, data)) = job else { return };
                        let mut out = Vec::with_capacity(MAX_BLOCK);
                        d.block(&data, &mut out);
                        if res_tx.send((seq, out, data)).is_err() {
                            return;
                        }
                    }
                }));
            }
            Some(Pool { job_tx: Some(job_tx), res_rx, handles, workers: threads })
        };
        Ok(BgzfWriter {
            path: path.to_string(),
            out: BufWriter::with_capacity(4 << 20, file),
            deflater: Deflater::new(level),
            pool,
            block: Vec::with_capacity(BLOCK_DATA),
            spare: Vec::new(),
            taken: 0,
            next_seq: 0,
            written_seq: 0,
            pending: BTreeMap::new(),
            block_addr: vec![0],
            file_pos: 0,
            compressed: Vec::with_capacity(MAX_BLOCK),
        })
    }

    /// Uncompressed bytes taken so far: the stream offset of the next byte.
    #[inline]
    pub fn position(&self) -> u64 {
        self.taken
    }

    /// The virtual offset of stream offset `u`, if the address of its block is known yet. Blocks hold
    /// BLOCK_DATA bytes each but the last, so an offset at the end of a full block is the start of the
    /// next one (as htsjdk's writer has it).
    #[inline]
    pub fn virtual_offset(&self, u: u64) -> Option<u64> {
        let k = (u / BLOCK_DATA as u64) as usize;
        let off = u % BLOCK_DATA as u64;
        self.block_addr.get(k).map(|a| (a << 16) | off)
    }

    pub fn write(&mut self, mut data: &[u8]) -> Result<()> {
        self.taken += data.len() as u64;
        while !data.is_empty() {
            let n = (BLOCK_DATA - self.block.len()).min(data.len());
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.block.len() == BLOCK_DATA {
                self.submit()?;
            }
        }
        Ok(())
    }

    fn submit(&mut self) -> Result<()> {
        let fresh = self.spare.pop().unwrap_or_else(|| Vec::with_capacity(BLOCK_DATA));
        let data = std::mem::replace(&mut self.block, fresh);
        let seq = self.next_seq;
        self.next_seq += 1;
        if let Some(pool) = &self.pool {
            pool.job_tx.as_ref().unwrap().send((seq, data)).map_err(|_| anyhow::anyhow!("{}: compression thread failed", self.path))?;
            let limit = pool.workers as u64 * 3 + 2;
            while self.next_seq - self.written_seq > limit {
                self.collect(true)?;
            }
            self.collect(false)?;
        } else {
            let mut c = std::mem::take(&mut self.compressed);
            c.clear();
            self.deflater.block(&data, &mut c);
            self.put(&c)?;
            self.compressed = c;
            let mut data = data;
            data.clear();
            self.spare.push(data);
            self.written_seq += 1;
        }
        Ok(())
    }

    /// Takes finished blocks from the workers (waiting for one if `wait`) and writes those next in line.
    fn collect(&mut self, wait: bool) -> Result<()> {
        let Some(pool) = &self.pool else { return Ok(()) };
        if wait {
            let (seq, c, d) = pool.res_rx.recv().map_err(|_| anyhow::anyhow!("{}: compression thread failed", self.path))?;
            self.pending.insert(seq, (c, d));
        }
        while let Ok((seq, c, d)) = pool.res_rx.try_recv() {
            self.pending.insert(seq, (c, d));
        }
        while let Some((c, mut d)) = self.pending.remove(&self.written_seq) {
            self.put(&c)?;
            d.clear();
            self.spare.push(d);
            self.written_seq += 1;
        }
        Ok(())
    }

    fn put(&mut self, block: &[u8]) -> Result<()> {
        self.out.write_all(block).with_context(|| format!("cannot write {}", self.path))?;
        self.file_pos += block.len() as u64;
        self.block_addr.push(self.file_pos);
        Ok(())
    }

    /// Writes the last block and the end-of-file block. `on_block` sees nothing new: all virtual
    /// offsets are known afterwards.
    pub fn finish(&mut self) -> Result<()> {
        if !self.block.is_empty() {
            self.submit()?;
        }
        if let Some(pool) = &mut self.pool {
            pool.job_tx = None;
        }
        while self.written_seq < self.next_seq {
            self.collect(true)?;
        }
        if let Some(pool) = &mut self.pool {
            for h in pool.handles.drain(..) {
                let _ = h.join();
            }
        }
        self.out.write_all(&EOF_BLOCK)?;
        self.out.flush().with_context(|| format!("cannot write {}", self.path))?;
        Ok(())
    }
}

impl Drop for BgzfWriter {
    fn drop(&mut self) {
        if let Some(pool) = &mut self.pool {
            pool.job_tx = None;
            while pool.res_rx.recv().is_ok() {}
            for h in pool.handles.drain(..) {
                let _ = h.join();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(threads_w: usize, threads_r: usize, n: usize) {
        let dir = std::env::temp_dir().join(format!("mdrs_bgzf_{}_{threads_w}_{threads_r}_{n}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.bgzf").to_string_lossy().to_string();
        let data: Vec<u8> = (0..n).map(|i| ((i * 7919) % 251) as u8).collect();
        let mut w = BgzfWriter::create(&path, 2, threads_w).unwrap();
        for chunk in data.chunks(1000) {
            w.write(chunk).unwrap();
        }
        w.finish().unwrap();
        // every offset maps to a virtual offset after finishing
        assert!(w.virtual_offset(n as u64).is_some());
        drop(w);
        let mut r = BgzfReader::open(&path, threads_r).unwrap();
        let mut back = vec![0u8; n];
        r.read_exact(&mut back).unwrap();
        assert_eq!(back, data);
        let mut one = [0u8; 1];
        assert!(!r.read_exact_or_eof(&mut one).unwrap());
        drop(r);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn roundtrip_inline_and_threads() {
        roundtrip(0, 0, 0);
        roundtrip(0, 0, 5);
        roundtrip(0, 0, BLOCK_DATA);
        roundtrip(2, 0, BLOCK_DATA * 3 + 17);
        roundtrip(0, 2, BLOCK_DATA * 40 + 1);
        roundtrip(3, 3, 3_000_000);
    }

    #[test]
    fn incompressible_blocks_are_stored() {
        let mut d = Deflater::new(9);
        let mut x = 12345u64;
        let data: Vec<u8> = (0..BLOCK_DATA)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                x as u8
            })
            .collect();
        let mut out = Vec::new();
        d.block(&data, &mut out);
        assert!(out.len() <= MAX_BLOCK);
        let mut back = Vec::new();
        inflate_block(&mut Inflater::new(), &out, &mut back, "t").unwrap();
        assert_eq!(back, data);
    }
}
