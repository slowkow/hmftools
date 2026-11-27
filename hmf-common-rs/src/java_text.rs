//! Java's text handling where the hmf-common file readers depend on it: `Integer.parseInt`,
//! `Long.parseLong`, `Double.parseDouble`, `Boolean.parseBoolean`, `String.split`, `String.compareTo`
//! and `FileReaderUtils.createFieldsIndexMap`.

use std::cmp::Ordering;
use std::collections::HashMap;

/// Integer.parseInt (ASCII digits, an optional sign); None where Java throws NumberFormatException
pub fn parse_int(s: &str) -> Option<i32> {
    if s.is_empty() || s == "+" || s == "-" {
        return None;
    }
    s.parse::<i32>().ok()
}

/// Long.parseLong
pub fn parse_long(s: &str) -> Option<i64> {
    if s.is_empty() || s == "+" || s == "-" {
        return None;
    }
    s.parse::<i64>().ok()
}

/// Byte.parseByte
pub fn parse_byte(s: &str) -> Option<i8> {
    parse_int(s).and_then(|v| i8::try_from(v).ok())
}

/// Double.parseDouble for decimal input: white space (chars up to ' ') around the number is ignored, a
/// trailing d, D, f or F is accepted, and "NaN" and "Infinity" with an optional sign; hexadecimal floats
/// are not handled.
pub fn parse_double(s: &str) -> Option<f64> {
    let t = s.trim_matches(|c: char| c <= ' ');
    let (sign, body) = match t.as_bytes().first() {
        Some(b'-') => (-1.0, &t[1..]),
        Some(b'+') => (1.0, &t[1..]),
        _ => (1.0, t),
    };
    match body {
        "NaN" => return Some(f64::NAN),
        "Infinity" => return Some(sign * f64::INFINITY),
        _ => {}
    }
    let body = body.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(body);
    if body.is_empty() || !body.bytes().any(|b| b.is_ascii_digit()) {
        return None;
    }
    if !body.bytes().all(|b| b.is_ascii_digit() || b == b'.' || b == b'e' || b == b'E' || b == b'+' || b == b'-') {
        return None;
    }
    body.parse::<f64>().ok().map(|v| sign * v)
}

/// Boolean.parseBoolean: true for "true" in any case, false for anything else
pub fn parse_bool(s: &str) -> bool {
    s.eq_ignore_ascii_case("true")
}

/// String.split(delim) with a literal one-character delimiter and no limit: trailing empty strings are
/// removed, and an empty input gives one empty string.
pub fn split(s: &str, delim: char) -> Vec<&str> {
    if s.is_empty() {
        return vec![""];
    }
    let mut parts: Vec<&str> = s.split(delim).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// String.split(delim, -1): every field, trailing empty ones included
pub fn split_all(s: &str, delim: char) -> Vec<&str> {
    s.split(delim).collect()
}

/// String.compareTo: by UTF-16 code units
pub fn compare(a: &str, b: &str) -> Ordering {
    if a.is_ascii() && b.is_ascii() {
        return a.cmp(b);
    }
    a.encode_utf16().cmp(b.encode_utf16())
}

/// FileReaderUtils.createFieldsIndexMap: column name to index; a repeated name keeps its last index
pub fn fields_index_map(header: &str, delim: char) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    for (i, f) in split_all(header, delim).into_iter().enumerate() {
        m.insert(f.to_string(), i);
    }
    m
}

/// Lines of a text file as Files.readAllLines and BufferedReader.readLine give them: split at "\n",
/// "\r\n" or "\r", without a last empty line; a `.gz` file is decompressed (htslib's BGZF reader reads
/// plain gzip as well).
pub fn read_lines(path: &str) -> std::io::Result<Vec<String>> {
    let bytes = if path.ends_with(".gz") {
        use std::io::Read as _;
        let mut r = rust_htslib::bgzf::Reader::from_path(path).map_err(|e| std::io::Error::other(format!("{path}: {e}")))?;
        let mut out = Vec::new();
        r.read_to_end(&mut out)?;
        out
    } else {
        std::fs::read(path)?
    };
    let text = match String::from_utf8(bytes) {
        Ok(s) => s,
        Err(e) => String::from_utf8_lossy(e.as_bytes()).into_owned(),
    };
    Ok(lines_of(&text))
}

/// Each line of a text file in turn, split as BufferedReader.readLine splits them (at "\n", "\r\n" or
/// "\r"), without holding the file in memory; a `.gz` file is decompressed. Stops at the first error `f`
/// returns.
pub fn for_each_line<F: FnMut(&str) -> anyhow::Result<()>>(path: &str, mut f: F) -> anyhow::Result<()> {
    use std::io::BufRead as _;
    let raw: Box<dyn std::io::Read> = if path.ends_with(".gz") {
        Box::new(rust_htslib::bgzf::Reader::from_path(path).map_err(|e| anyhow::anyhow!("{path}: {e}"))?)
    } else {
        Box::new(std::fs::File::open(path).map_err(|e| anyhow::anyhow!("{path}: {e}"))?)
    };
    let mut reader = std::io::BufReader::with_capacity(1 << 20, raw);
    let mut buf: Vec<u8> = Vec::with_capacity(512);
    loop {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf).map_err(|e| anyhow::anyhow!("{path}: {e}"))?;
        if n == 0 {
            return Ok(());
        }
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        if buf.last() == Some(&b'\r') {
            buf.pop();
        }
        let text = String::from_utf8_lossy(&buf);
        if text.contains('\r') {
            for part in text.split('\r') {
                f(part)?;
            }
        } else {
            f(&text)?;
        }
    }
}

/// The first line of a text file (None for an empty file)
pub fn first_line(path: &str) -> anyhow::Result<Option<String>> {
    let mut first: Option<String> = None;
    let stop = "first line read";
    match for_each_line(path, |line| {
        first = Some(line.to_string());
        Err(anyhow::anyhow!(stop))
    }) {
        Ok(()) => Ok(first),
        Err(e) if first.is_some() && e.to_string() == stop => Ok(first),
        Err(e) => Err(e),
    }
}

/// The lines of a text as BufferedReader.readLine splits them
pub fn lines_of(text: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let b = text.as_bytes();
    let mut start = 0;
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'\n' || b[i] == b'\r' {
            lines.push(text[start..i].to_string());
            if b[i] == b'\r' && i + 1 < b.len() && b[i + 1] == b'\n' {
                i += 1;
            }
            start = i + 1;
        }
        i += 1;
    }
    if start < b.len() {
        lines.push(text[start..].to_string());
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers() {
        assert_eq!(parse_int("-12"), Some(-12));
        assert_eq!(parse_int("+7"), Some(7));
        assert_eq!(parse_int(" 7"), None);
        assert_eq!(parse_int("2147483648"), None);
        assert_eq!(parse_byte("-1"), Some(-1));
        assert_eq!(parse_byte("200"), None);
        assert_eq!(parse_double(" 1.5d "), Some(1.5));
        assert_eq!(parse_double("-Infinity"), Some(f64::NEG_INFINITY));
        assert!(parse_double("NaN").unwrap().is_nan());
        assert_eq!(parse_double("inf"), None);
        assert_eq!(parse_double("1e3"), Some(1000.0));
        assert_eq!(parse_double(""), None);
        assert!(parse_bool("TRUE") && !parse_bool("yes"));
    }

    #[test]
    fn splitting() {
        assert_eq!(split("a,b,,", ','), vec!["a", "b"]);
        assert_eq!(split(",,", ','), Vec::<&str>::new());
        assert_eq!(split("", ','), vec![""]);
        assert_eq!(split(",a", ','), vec!["", "a"]);
        assert_eq!(split_all("a,b,,", ','), vec!["a", "b", "", ""]);
        assert_eq!(lines_of("a\r\nb\n\nc"), vec!["a", "b", "", "c"]);
        assert_eq!(lines_of("a\n"), vec!["a"]);
        assert_eq!(compare("B", "a"), Ordering::Less);
        let p = std::env::temp_dir().join(format!("hmf_common_lines_{}", std::process::id()));
        std::fs::write(&p, "a\r\nb\n\nc\rd\n").unwrap();
        let mut seen = Vec::new();
        for_each_line(p.to_str().unwrap(), |l| {
            seen.push(l.to_string());
            Ok(())
        })
        .unwrap();
        assert_eq!(seen, read_lines(p.to_str().unwrap()).unwrap());
        assert_eq!(seen, vec!["a", "b", "", "c", "d"]);
        std::fs::remove_file(&p).ok();
    }
}
