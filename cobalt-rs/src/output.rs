//! The output files: CobaltRatioFile (`<sample>.cobalt.ratio.tsv.gz`), CobaltGcMedianFile
//! (`.cobalt.gc.median.tsv`), CobaltMedianRatioFile (`.cobalt.ratio.median.tsv`) and the version file. The
//! .pcf files are written by hmf_common::segmentation.

use anyhow::{Context, Result};
use std::io::Write;

use crate::calc::{CollatedChr, GcMedianReadDepth, MedianRatio};
use hmf_common::fmt::{decimal_format, java_fixed};

/// DecimalFormat("0.0000"): four decimals, rounded as DecimalFormat("#.####") rounds
pub fn decimal_format_4(v: f64) -> String {
    let mut s = decimal_format(v, 4);
    if !v.is_finite() {
        return s;
    }
    let decimals = match s.find('.') {
        Some(i) => s.len() - i - 1,
        None => {
            s.push('.');
            0
        }
    };
    for _ in decimals..4 {
        s.push('0');
    }
    s
}

/// CobaltRatioFile.write: the rows sorted by chromosome and position, DecimalFormat("#.####"), gzipped
/// (zlib level 6, as java.util.zip writes it)
pub fn write_ratios(path: &str, v38: bool, rows: &[CollatedChr]) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let gz = flate2::write::GzEncoder::new(file, flate2::Compression::new(6));
    let mut w = std::io::BufWriter::with_capacity(1 << 17, gz);
    w.write_all(b"chromosome\tposition\treferenceReadDepth\ttumorReadDepth\treferenceGCRatio\ttumorGCRatio\treferenceGCDiploidRatio\treferenceGCContent\ttumorGCContent\n")?;
    let mut line = String::with_capacity(128);
    for c in rows {
        let chromosome = c.chr.versioned(v38);
        for i in 0..c.len() {
            let r = c.row(i);
            use std::fmt::Write as _;
            line.clear();
            let _ = write!(line, "{chromosome}\t{}\t", r.position);
            for (i, v) in [
                r.reference_read_depth,
                r.tumor_read_depth,
                r.reference_gc_ratio,
                r.tumor_gc_ratio,
                r.reference_gc_diploid_ratio,
                r.reference_gc_content,
                r.tumor_gc_content,
            ]
            .into_iter()
            .enumerate()
            {
                if i > 0 {
                    line.push('\t');
                }
                push_decimal(&mut line, v);
            }
            line.push('\n');
            w.write_all(line.as_bytes())?;
        }
    }
    let gz = w.into_inner().map_err(|e| anyhow::anyhow!("writing {path}: {}", e.error()))?;
    gz.finish().with_context(|| format!("writing {path}"))?;
    Ok(())
}

/// DecimalFormat("#.####"), with the common values written directly
#[inline]
fn push_decimal(out: &mut String, v: f64) {
    if v == -1.0 {
        out.push_str("-1");
    } else if v == 0.0 && v.is_sign_positive() {
        out.push('0');
    } else if v == 1.0 {
        out.push('1');
    } else {
        out.push_str(&decimal_format(v, 4));
    }
}

/// CobaltGcMedianFile.write
pub fn write_gc_median(path: &str, stats: &GcMedianReadDepth) -> Result<()> {
    let mut text = String::new();
    text.push_str("#sampleMean\tsampleMedian\n");
    text.push_str(&format!("{}\t{}\n", java_fixed(stats.mean, 2), java_fixed(stats.median, 2)));
    text.push_str("#gcBucket\tmedian\n");
    for (bucket, depth) in &stats.bucket_medians {
        if *depth > 0.0 {
            text.push_str(&format!("{bucket}\t{}\n", java_fixed(*depth, 2)));
        }
    }
    std::fs::write(path, text).with_context(|| format!("writing {path}"))
}

/// CobaltMedianRatioFile.write
pub fn write_median_ratios(path: &str, ratios: &[MedianRatio]) -> Result<()> {
    let mut text = String::from("chromosome\tmedianRatio\tcount\n");
    for r in ratios {
        text.push_str(&format!("{}\t{}\t{}\n", r.chromosome, decimal_format_4(r.median_ratio), r.count));
    }
    std::fs::write(path, text).with_context(|| format!("writing {path}"))
}

/// VersionInfo.write: the version resource as is (Maven fills in the version and build time)
pub fn write_version(dir: &str) -> Result<()> {
    let text = format!("version={}\nbuild.date={}", crate::VERSION, env!("COBALT_BUILD_DATE"));
    let path = format!("{dir}/cobalt.version");
    std::fs::write(&path, text).with_context(|| format!("writing {path}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats() {
        assert_eq!(decimal_format_4(1.0), "1.0000");
        assert_eq!(decimal_format_4(0.98925), "0.9892");
        assert_eq!(decimal_format_4(f64::NAN), "NaN");
        let mut s = String::new();
        for v in [-1.0, 0.0, 1.0, 0.79324, 94.609] {
            push_decimal(&mut s, v);
            s.push(' ');
        }
        assert_eq!(s, "-1 0 1 0.7932 94.609 ");
    }
}
