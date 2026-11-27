//! hmf-common `utils.file.DelimFileReader` and the writing side of `FileWriterUtils.createBufferedWriter`
//! / `DelimFileWriter` (added 2026-11-27 for qsee).
//!
//! - `DelimFileReader`: a tab-separated file (gzip when the name ends in `.gz`) with its header as the
//!   column names; a repeated column name is an error; a value `null` or `NULL` is an error unless asked
//!   for with the `_or_null` getters; numbers parse as Integer.parseInt, Long.parseLong and
//!   Double.parseDouble.
//! - `JavaWriter`: a file written as FileWriterUtils.createBufferedWriter writes it: plain, or for a name
//!   ending in `.gz` a GZIPOutputStream (zlib level 6, header with no time stamp and OS 255), whose
//!   compressed bytes equal Java's.
//! - `DelimRow`: a row in the column order of the header, as DelimFileWriter.Row; an unset value is
//!   written as `null`; `write_delim_file` writes a header and rows as DelimFileWriter.write does.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::HashMap;
use std::io::Write;

use crate::java_text::{parse_double, parse_int, parse_long, read_lines, split_all};

/// DelimFileReader on a whole file: the header and the remaining lines
pub struct DelimFileReader {
    columns: Vec<String>,
    index: HashMap<String, usize>,
    lines: Vec<String>,
}

impl DelimFileReader {
    /// new DelimFileReader(filename): reads the file and its header (an empty file is an error, as the
    /// NullPointerException of Java's split on a missing first line; a file that cannot be read is an
    /// UncheckedIOException)
    pub fn open(path: &str) -> Result<DelimFileReader> {
        let mut lines = read_lines(path).with_context(|| format!("java.io.UncheckedIOException: {path}"))?;
        if lines.is_empty() {
            bail!("java.lang.NullPointerException: {path} has no header line");
        }
        let header = lines.remove(0);
        let r = Self::from_header(&header, lines)?;
        Ok(r)
    }

    /// a reader on lines already read: `header` gives the column names
    pub fn from_header(header: &str, lines: Vec<String>) -> Result<DelimFileReader> {
        let columns: Vec<String> = split_all(header, '\t').into_iter().map(String::from).collect();
        let mut index = HashMap::new();
        for (i, c) in columns.iter().enumerate() {
            if index.insert(c.clone(), i).is_some() {
                bail!("duplicate column name: {c}");
            }
        }
        Ok(DelimFileReader { columns, index, lines })
    }

    pub fn column_names(&self) -> &[String] {
        &self.columns
    }
    pub fn has_column(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }
    /// getColumnIndex: None for a column the file does not have
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    /// the rows, each split at every tab (String.split(delim, -1))
    pub fn rows(&self) -> impl Iterator<Item = DelimRowValues<'_>> {
        self.lines.iter().map(move |l| DelimRowValues { index: &self.index, values: split_all(l, '\t') })
    }
}

/// DelimFileReader.Row
pub struct DelimRowValues<'a> {
    index: &'a HashMap<String, usize>,
    values: Vec<&'a str>,
}

impl<'a> DelimRowValues<'a> {
    fn column(&self, name: &str) -> Result<usize> {
        self.index.get(name).copied().ok_or_else(|| anyhow!("java.util.NoSuchElementException: column: {name} not found"))
    }
    fn raw(&self, i: usize) -> Result<&'a str> {
        self.values.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", self.values.len()))
    }
    fn value(&self, i: usize, allow_null: bool) -> Result<Option<&'a str>> {
        let v = self.raw(i)?;
        if v == "null" || v == "NULL" {
            if allow_null {
                return Ok(None);
            }
            bail!("java.util.NoSuchElementException: column {i} is null");
        }
        Ok(Some(v))
    }
    /// Row.get(String)
    pub fn get(&self, name: &str) -> Result<&'a str> {
        Ok(self.value(self.column(name)?, false)?.unwrap())
    }
    /// Row.get(int)
    pub fn get_at(&self, i: usize) -> Result<&'a str> {
        Ok(self.value(i, false)?.unwrap())
    }
    pub fn get_or_null(&self, name: &str) -> Result<Option<&'a str>> {
        self.value(self.column(name)?, true)
    }
    pub fn get_int(&self, name: &str) -> Result<i32> {
        let v = self.get(name)?;
        parse_int(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    pub fn get_long(&self, name: &str) -> Result<i64> {
        let v = self.get(name)?;
        parse_long(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    pub fn get_double(&self, name: &str) -> Result<f64> {
        let v = self.get(name)?;
        parse_double(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
    pub fn get_double_at(&self, i: usize) -> Result<f64> {
        let v = self.get_at(i)?;
        parse_double(v).ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{v}\""))
    }
}

/// A file opened as FileWriterUtils.createBufferedWriter opens it
pub enum JavaWriter {
    Plain(std::io::BufWriter<std::fs::File>),
    Gzip(std::io::BufWriter<flate2::write::GzEncoder<std::fs::File>>),
}

impl JavaWriter {
    /// a gzip stream for a name ending in `.gz`, else a plain file
    pub fn create(path: &str) -> Result<JavaWriter> {
        let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
        Ok(if path.ends_with(".gz") {
            // java.util.zip.GZIPOutputStream: Deflater.DEFAULT_COMPRESSION (zlib's level 6)
            JavaWriter::Gzip(std::io::BufWriter::with_capacity(1 << 16, flate2::write::GzEncoder::new(file, flate2::Compression::new(6))))
        } else {
            JavaWriter::Plain(std::io::BufWriter::with_capacity(1 << 16, file))
        })
    }

    /// flushes and, for gzip, writes the trailer
    pub fn finish(self) -> Result<()> {
        match self {
            JavaWriter::Plain(mut w) => w.flush()?,
            JavaWriter::Gzip(w) => {
                let enc = w.into_inner().map_err(|e| anyhow!("writing: {}", e.error()))?;
                enc.finish()?.flush()?;
            }
        }
        Ok(())
    }
}

impl Write for JavaWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            JavaWriter::Plain(w) => w.write(buf),
            JavaWriter::Gzip(w) => w.write(buf),
        }
    }
    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        match self {
            JavaWriter::Plain(w) => w.write_all(buf),
            JavaWriter::Gzip(w) => w.write_all(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            JavaWriter::Plain(w) => w.flush(),
            JavaWriter::Gzip(w) => w.flush(),
        }
    }
}

/// DelimFileWriter.Row: values by column name, in header order; a value never set is written as "null"
pub struct DelimRow<'a> {
    columns: &'a [String],
    values: Vec<Option<String>>,
}

impl<'a> DelimRow<'a> {
    pub fn new(columns: &'a [String]) -> DelimRow<'a> {
        DelimRow { columns, values: vec![None; columns.len()] }
    }
    /// Row.set(key, value); an unknown column is an error, as Java's IllegalArgumentException
    pub fn set(&mut self, key: &str, value: impl Into<String>) -> Result<()> {
        let i = self.columns.iter().position(|c| c == key).ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: invalid column: {key}"))?;
        self.values[i] = Some(value.into());
        Ok(())
    }
    /// Row.set by the column's index (the position in the header)
    pub fn set_at(&mut self, index: usize, value: impl Into<String>) {
        self.values[index] = Some(value.into());
    }
    /// Row.set(key, (String) null): written as "null"
    pub fn set_null(&mut self, key: &str) -> Result<()> {
        let i = self.columns.iter().position(|c| c == key).ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: invalid column: {key}"))?;
        self.values[i] = None;
        Ok(())
    }
    /// the row joined by tabs, ended by a newline
    pub fn line(&self) -> String {
        let mut s = String::new();
        for (i, v) in self.values.iter().enumerate() {
            if i > 0 {
                s.push('\t');
            }
            s.push_str(v.as_deref().unwrap_or("null"));
        }
        s.push('\n');
        s
    }
}

/// DelimFileWriter.write(filename, columns, objects, encoder): the header, then a line per object; a
/// repeated column name is an error (RuntimeException "duplicate column")
pub fn write_delim_file<T>(path: &str, columns: &[String], objects: &[T], mut encode: impl FnMut(&T, &mut DelimRow) -> Result<()>) -> Result<()> {
    for (i, c) in columns.iter().enumerate() {
        if columns[..i].contains(c) {
            bail!("duplicate column: {c}");
        }
    }
    let mut w = JavaWriter::create(path)?;
    w.write_all(columns.join("\t").as_bytes())?;
    w.write_all(b"\n")?;
    // as Java's try-with-resources: a row that fails ends the writing, and the file is closed with the
    // rows before it
    let mut result = Ok(());
    for o in objects {
        let mut row = DelimRow::new(columns);
        if let Err(e) = encode(o, &mut row) {
            result = Err(e);
            break;
        }
        w.write_all(row.line().as_bytes())?;
    }
    w.finish()?;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> String {
        let dir = std::env::temp_dir().join(format!("hmf_common_delim_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name).to_str().unwrap().to_string()
    }

    #[test]
    fn reader_like_java() {
        let p = tmp("a.tsv");
        std::fs::write(&p, "a\tb\tc\n1\tnull\t2.5\n-3\tx\t\n").unwrap();
        let r = DelimFileReader::open(&p).unwrap();
        assert_eq!(r.column_names(), ["a", "b", "c"]);
        let rows: Vec<_> = r.rows().collect();
        assert_eq!(rows[0].get_int("a").unwrap(), 1);
        assert!(rows[0].get("b").is_err());
        assert_eq!(rows[0].get_or_null("b").unwrap(), None);
        assert_eq!(rows[0].get_double("c").unwrap(), 2.5);
        assert_eq!(rows[1].get("c").unwrap(), "");
        assert!(rows[1].get_double("c").is_err());
        assert!(rows[1].get("d").is_err());
        std::fs::write(&p, "a\ta\n").unwrap();
        assert!(DelimFileReader::open(&p).is_err());
        std::fs::write(&p, "").unwrap();
        assert!(DelimFileReader::open(&p).is_err());
    }

    #[test]
    fn writer_gzip_and_plain() {
        let cols: Vec<String> = ["x", "y"].iter().map(|s| s.to_string()).collect();
        let p = tmp("w.tsv.gz");
        write_delim_file(&p, &cols, &[1, 2], |v, row| {
            row.set("x", v.to_string())?;
            if *v == 1 { row.set("y", "one") } else { Ok(()) }
        })
        .unwrap();
        let bytes = std::fs::read(&p).unwrap();
        // the GZIPOutputStream header: no time stamp, no extra flags, OS 255
        assert_eq!(&bytes[..10], &[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff]);
        assert_eq!(read_lines(&p).unwrap(), vec!["x\ty", "1\tone", "2\tnull"]);
        let q = tmp("w.tsv");
        write_delim_file(&q, &cols, &[3], |v, row| row.set("y", v.to_string())).unwrap();
        assert_eq!(std::fs::read_to_string(&q).unwrap(), "x\ty\nnull\t3\n");
        let dup: Vec<String> = ["x", "x"].iter().map(|s| s.to_string()).collect();
        assert!(write_delim_file(&q, &dup, &[3], |_, _| Ok(())).is_err());
    }
}
