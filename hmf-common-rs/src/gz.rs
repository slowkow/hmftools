//! Reading gzip and BGZF text files fast: a streaming inflater (zlib-rs, about twice as fast as the
//! default backends) that reads member after member as GZIPInputStream does, and a line reader that
//! hands out each line as a byte slice without copying it.
//!
//! The resource files the wave B tools read (PAVE's PON, gnomAD and mappability files, ClinVar) are
//! hundreds of megabytes of gzip each; inflating them is most of their reading time.

use std::io::Read;

/// A gzip file as a stream of its decompressed bytes, member after member (as GZIPInputStream reads
/// it). Moved here from esvee-rs.
pub struct GzReader<R: Read> {
    inner: R,
    input: Vec<u8>,
    pos: usize,
    len: usize,
    eof: bool,
    inflater: Option<zlib_rs::Inflate>,
}

impl<R: Read> GzReader<R> {
    pub fn new(inner: R) -> GzReader<R> {
        GzReader { inner, input: vec![0; 1 << 17], pos: 0, len: 0, eof: false, inflater: None }
    }

    fn fill(&mut self) -> std::io::Result<()> {
        if self.pos < self.len {
            self.input.copy_within(self.pos..self.len, 0);
        }
        self.len -= self.pos;
        self.pos = 0;
        while !self.eof && self.len < self.input.len() {
            match self.inner.read(&mut self.input[self.len..]) {
                Ok(0) => self.eof = true,
                Ok(n) => self.len += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
}

impl<R: Read> Read for GzReader<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            if self.len - self.pos < 1 << 15 && !self.eof {
                self.fill()?;
            }
            if self.inflater.is_none() {
                // the next member, if what follows is one
                let rest = &self.input[self.pos..self.len];
                if rest.len() < 2 || rest[0] != 0x1f || rest[1] != 0x8b {
                    return Ok(0);
                }
                self.inflater = Some(zlib_rs::Inflate::new(true, 31));
            }
            let inflater = self.inflater.as_mut().unwrap();
            let (before_in, before_out) = (inflater.total_in(), inflater.total_out());
            let status = inflater
                .decompress(&self.input[self.pos..self.len], out, zlib_rs::InflateFlush::NoFlush)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("invalid gzip data: {e:?}")))?;
            let consumed = (inflater.total_in() - before_in) as usize;
            let produced = (inflater.total_out() - before_out) as usize;
            self.pos += consumed;
            if status == zlib_rs::Status::StreamEnd {
                self.inflater = None;
            } else if consumed == 0 && produced == 0 && self.eof {
                return Err(std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Unexpected end of ZLIB input stream"));
            }
            if produced > 0 {
                return Ok(produced);
            }
        }
    }
}

/// A file as a byte stream: `.gz` and `.bgz` files (by extension, as the Java readers decide) are
/// decompressed.
pub fn open(path: &str) -> std::io::Result<Box<dyn Read + Send>> {
    let file = std::fs::File::open(path)?;
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        Ok(Box::new(GzReader::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// A gzip file inflated on a thread of its own, the reader taking the inflated chunks as they come: for a
/// large single-member gzip, whose inflating cannot be split, this overlaps it with the parsing.
pub fn open_threaded(path: &str) -> std::io::Result<Box<dyn Read + Send>> {
    if !(path.ends_with(".gz") || path.ends_with(".bgz")) {
        return open(path);
    }
    let file = std::fs::File::open(path)?;
    let (tx, rx) = std::sync::mpsc::sync_channel::<std::io::Result<Vec<u8>>>(4);
    std::thread::spawn(move || {
        let mut r = GzReader::new(file);
        loop {
            let mut buf = vec![0u8; 1 << 20];
            let mut filled = 0;
            while filled < buf.len() {
                match r.read(&mut buf[filled..]) {
                    Ok(0) => break,
                    Ok(n) => filled += n,
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                }
            }
            if filled == 0 {
                return;
            }
            buf.truncate(filled);
            if tx.send(Ok(buf)).is_err() {
                return;
            }
        }
    });
    Ok(Box::new(ChannelReader { rx, current: Vec::new(), pos: 0 }))
}

struct ChannelReader {
    rx: std::sync::mpsc::Receiver<std::io::Result<Vec<u8>>>,
    current: Vec<u8>,
    pos: usize,
}

impl Read for ChannelReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        while self.pos == self.current.len() {
            match self.rx.recv() {
                Ok(Ok(chunk)) => {
                    self.current = chunk;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Ok(0),
            }
        }
        let n = out.len().min(self.current.len() - self.pos);
        out[..n].copy_from_slice(&self.current[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Lines of a stream as byte slices, split at '\n' with a trailing '\r' removed (BufferedReader.readLine
/// for files with Unix or Windows line ends).
pub struct LineReader {
    inner: Box<dyn Read + Send>,
    buf: Vec<u8>,
    start: usize,
    end: usize,
    eof: bool,
}

impl LineReader {
    pub fn new(inner: Box<dyn Read + Send>) -> LineReader {
        LineReader { inner, buf: vec![0; 1 << 20], start: 0, end: 0, eof: false }
    }

    pub fn open(path: &str) -> std::io::Result<LineReader> {
        Ok(LineReader::new(open(path)?))
    }

    /// As `open`, a gzip file inflated on a thread of its own (`open_threaded`)
    pub fn open_threaded(path: &str) -> std::io::Result<LineReader> {
        Ok(LineReader::new(open_threaded(path)?))
    }

    /// The next line, None at the end of the stream
    pub fn next_line(&mut self) -> std::io::Result<Option<&[u8]>> {
        loop {
            if let Some(i) = memchr::memchr(b'\n', &self.buf[self.start..self.end]) {
                let line_start = self.start;
                let mut line_end = self.start + i;
                self.start = line_end + 1;
                if line_end > line_start && self.buf[line_end - 1] == b'\r' {
                    line_end -= 1;
                }
                return Ok(Some(&self.buf[line_start..line_end]));
            }
            if self.eof {
                if self.start < self.end {
                    let line_start = self.start;
                    let mut line_end = self.end;
                    self.start = self.end;
                    if self.buf[line_end - 1] == b'\r' {
                        line_end -= 1;
                    }
                    return Ok(Some(&self.buf[line_start..line_end]));
                }
                return Ok(None);
            }
            // keep the partial line and read more
            if self.start > 0 {
                self.buf.copy_within(self.start..self.end, 0);
                self.end -= self.start;
                self.start = 0;
            }
            if self.end == self.buf.len() {
                let n = self.buf.len() * 2;
                self.buf.resize(n, 0);
            }
            match self.inner.read(&mut self.buf[self.end..]) {
                Ok(0) => self.eof = true,
                Ok(n) => self.end += n,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn members_and_lines() {
        let dir = std::env::temp_dir().join(format!("hmf_common_gz_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.gz");
        // two gzip members, as bgzip writes them
        let mut bytes = Vec::new();
        for part in ["a\tb\r\nc", "d\n\nlast"] {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(5));
            enc.write_all(part.as_bytes()).unwrap();
            bytes.extend(enc.finish().unwrap());
        }
        std::fs::write(&path, bytes).unwrap();
        let mut r = LineReader::open(path.to_str().unwrap()).unwrap();
        let mut lines = Vec::new();
        while let Some(l) = r.next_line().unwrap() {
            lines.push(String::from_utf8(l.to_vec()).unwrap());
        }
        assert_eq!(lines, vec!["a\tb", "cd", "", "last"]);
        let mut r = LineReader::open_threaded(path.to_str().unwrap()).unwrap();
        let mut threaded = Vec::new();
        while let Some(l) = r.next_line().unwrap() {
            threaded.push(String::from_utf8(l.to_vec()).unwrap());
        }
        assert_eq!(threaded, lines);
        std::fs::remove_dir_all(&dir).ok();
    }
}
