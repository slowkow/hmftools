//! commons-csv 1.9 CSVPrinter as the Kotlin tools (CIDER) write their TSV files: a delimiter, '\n' as the
//! record separator, the quote character '"' and QuoteMode MINIMAL (the default), and an optional null
//! string.
//!
//! MINIMAL quotes a value that is empty and first in its record, that starts with a character up to '#'
//! (space, '!', '"', '#' and control characters), that contains the delimiter, a quote, CR or LF, or that
//! ends with a character up to ' '. Quotes inside a quoted value are doubled. A null value is written as
//! the null string (or nothing), never quoted.

use std::io::Write;

pub struct CsvPrinter<W: Write> {
    out: W,
    delimiter: u8,
    null_string: Option<String>,
    new_record: bool,
}

impl<W: Write> CsvPrinter<W> {
    /// A printer with `delimiter` and no null string (a null is written as nothing)
    pub fn new(out: W, delimiter: u8) -> Self {
        CsvPrinter { out, delimiter, null_string: None, new_record: true }
    }

    /// CSVFormat.setNullString
    pub fn with_null_string(mut self, s: &str) -> Self {
        self.null_string = Some(s.to_string());
        self
    }

    /// CSVPrinter.print(value), for a value already turned into its string (Kotlin/Java toString)
    pub fn print(&mut self, value: &str) -> std::io::Result<()> {
        if !self.new_record {
            self.out.write_all(&[self.delimiter])?;
        }
        let first = self.new_record;
        self.new_record = false;
        let b = value.as_bytes();
        let quote = if b.is_empty() {
            first
        } else if b[0] <= b'#' {
            true
        } else {
            b.iter().any(|&c| c == b'\n' || c == b'\r' || c == b'"' || c == self.delimiter) || b[b.len() - 1] <= b' '
        };
        if !quote {
            return self.out.write_all(b);
        }
        self.out.write_all(b"\"")?;
        let mut start = 0;
        for (i, &c) in b.iter().enumerate() {
            if c == b'"' {
                self.out.write_all(&b[start..i])?;
                self.out.write_all(b"\"")?;
                start = i;
            }
        }
        self.out.write_all(&b[start..])?;
        self.out.write_all(b"\"")
    }

    /// CSVPrinter.print(null)
    pub fn print_null(&mut self) -> std::io::Result<()> {
        if !self.new_record {
            self.out.write_all(&[self.delimiter])?;
        }
        self.new_record = false;
        if let Some(s) = &self.null_string {
            self.out.write_all(s.as_bytes())?;
        }
        Ok(())
    }

    /// print(value) or print(null)
    pub fn print_opt<S: AsRef<str>>(&mut self, value: Option<S>) -> std::io::Result<()> {
        match value {
            Some(v) => self.print(v.as_ref()),
            None => self.print_null(),
        }
    }

    /// CSVPrinter.println
    pub fn println(&mut self) -> std::io::Result<()> {
        self.new_record = true;
        self.out.write_all(b"\n")
    }

    /// printRecord(values): the header line of CSVFormat.setHeader
    pub fn print_record(&mut self, values: &[&str]) -> std::io::Result<()> {
        for v in values {
            self.print(v)?;
        }
        self.println()
    }

    pub fn get_mut(&mut self) -> &mut W {
        &mut self.out
    }

    pub fn into_inner(self) -> W {
        self.out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(values: &[Option<&str>], null: Option<&str>) -> String {
        let mut p = CsvPrinter::new(Vec::new(), b'\t');
        if let Some(n) = null {
            p = p.with_null_string(n);
        }
        for v in values {
            p.print_opt(*v).unwrap();
        }
        p.println().unwrap();
        String::from_utf8(p.into_inner()).unwrap()
    }

    #[test]
    fn minimal_quoting() {
        assert_eq!(line(&[Some("a"), Some(""), Some("b")], None), "a\t\tb\n");
        assert_eq!(line(&[Some(""), Some("b")], None), "\"\"\tb\n");
        assert_eq!(line(&[None, Some("b")], None), "\tb\n");
        assert_eq!(line(&[None, None], Some("null")), "null\tnull\n");
        assert_eq!(line(&[Some("#12"), Some("-1"), Some("x y"), Some("x ")], None), "\"#12\"\t-1\tx y\t\"x \"\n");
        assert_eq!(line(&[Some("a\"b"), Some("a\tb")], None), "\"a\"\"b\"\t\"a\tb\"\n");
        assert_eq!(line(&[Some("$1"), Some("!")], None), "$1\t\"!\"\n");
    }
}
