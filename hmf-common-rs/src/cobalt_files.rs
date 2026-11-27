//! hmf-common `cobalt` files as qsee reads them (added 2026-11-27): CobaltGcMedianFile
//! (`.cobalt.gc.median.tsv`), GcMedianReadDepth.

use anyhow::{anyhow, Context, Result};

use crate::java_text::{parse_double, parse_int, read_lines, split};

/// GcMedianReadDepth.NO_READ_DEPTH_VALUE
pub const NO_READ_DEPTH_VALUE: f64 = -1.0;

/// GcMedianReadDepth: the sample mean and median and the median per GC bucket. Java keeps the buckets in
/// a HashMap; here they are in file order, a repeated bucket keeping its last value (as HashMap.put).
#[derive(Clone, Debug, PartialEq)]
pub struct GcMedianReadDepth {
    pub mean: f64,
    pub median: f64,
    pub median_per_bucket: Vec<(i32, f64)>,
}

/// CobaltGcMedianFile.generateFilename: basePath + File.separator + sample + ".cobalt.gc.median.tsv"
pub fn cobalt_gc_median_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.cobalt.gc.median.tsv")
}

/// CobaltGcMedianFile.read
pub fn read_gc_median_file(path: &str) -> Result<GcMedianReadDepth> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    gc_median_from_lines(&lines).with_context(|| format!("reading {path}"))
}

/// CobaltGcMedianFile.fromLines: line 1 is skipped, line 2 has the mean and median (when present),
/// line 3 is skipped, and every later line of exactly two fields (String.split, trailing empty fields
/// dropped) is a bucket and its median
pub fn gc_median_from_lines(lines: &[String]) -> Result<GcMedianReadDepth> {
    let bad = |s: &str| anyhow!("java.lang.NumberFormatException: For input string: \"{s}\"");
    let (mut mean, mut median) = (0.0, 0.0);
    let mut i = 1;
    if lines.len() > i {
        let v = split(&lines[i], '\t');
        i += 1;
        let m = v.first().copied().unwrap_or("");
        mean = parse_double(m).ok_or_else(|| bad(m))?;
        let d = v.get(1).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index 1 out of bounds for length 1"))?;
        median = parse_double(d).ok_or_else(|| bad(d))?;
    }
    i += 1;
    let mut buckets: Vec<(i32, f64)> = Vec::new();
    while i < lines.len() {
        let v = split(&lines[i], '\t');
        if v.len() == 2 {
            let b = parse_int(v[0]).ok_or_else(|| bad(v[0]))?;
            let m = parse_double(v[1]).ok_or_else(|| bad(v[1]))?;
            match buckets.iter_mut().find(|(k, _)| *k == b) {
                Some(e) => e.1 = m,
                None => buckets.push((b, m)),
            }
        }
        i += 1;
    }
    Ok(GcMedianReadDepth { mean, median, median_per_bucket: buckets })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_median_lines() {
        let lines: Vec<String> = ["#sampleMean\tsampleMedian", "73.76\t62.55", "#gcBucket\tmedian", "24\t60.63", "25\t59.02\t", "x", "24\t61.00"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let g = gc_median_from_lines(&lines).unwrap();
        assert_eq!((g.mean, g.median), (73.76, 62.55));
        assert_eq!(g.median_per_bucket, vec![(24, 61.0), (25, 59.02)]);
        let empty = gc_median_from_lines(&["#sampleMean\tsampleMedian".to_string()]).unwrap();
        assert_eq!((empty.mean, empty.median, empty.median_per_bucket.len()), (0.0, 0.0, 0));
    }
}
