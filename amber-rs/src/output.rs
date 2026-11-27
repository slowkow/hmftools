//! The output files, written as the Java writers write them: AmberBAFFile, AmberQCFile,
//! TumorContaminationFile, RegionOfHomozygosityFile, PositionEvidenceFile and the version file.

use anyhow::{Context, Result};
use std::io::Write;

use crate::constants::{doubles, HOMOZYGOUS_REGION_LONG_SIZE, HOMOZYGOUS_REGION_MAX_HET_RATIO};
use crate::contamination::TumorContamination;
use crate::roh::RegionOfHomozygosity;
use crate::types::{AmberBaf, PositionEvidence};

/// java.text.DecimalFormat("0.0000"): the rounding of hmf_common::fmt::decimal_format (Java's DigitList
/// over the old FloatingDecimal digits, half-even), padded to four decimals as the pattern's minimum
pub fn decimal_format_4(v: f64) -> String {
    let mut s = hmf_common::fmt::decimal_format(v, 4);
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

/// The file under a text writer: plain, or gzipped as java.util.zip.GZIPOutputStream writes it (zlib
/// level 6, no flushes), so that the compressed bytes are Java's too
pub enum Sink {
    Plain(std::fs::File),
    Gz(flate2::write::GzEncoder<std::fs::File>),
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Sink::Plain(f) => f.write(buf),
            Sink::Gz(g) => g.write(buf),
        }
    }
    /// a flush of the gzip stream would add an empty block (zlib Z_SYNC_FLUSH), which Java does not write
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Sink::Plain(f) => f.flush(),
            Sink::Gz(_) => Ok(()),
        }
    }
}

pub type TextWriter = std::io::BufWriter<Sink>;

/// A writer that gzips when the name ends with .gz (FileWriterUtils.createBufferedWriter)
pub fn create_writer(path: &str) -> Result<TextWriter> {
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let sink =
        if path.ends_with(".gz") { Sink::Gz(flate2::write::GzEncoder::new(file, flate2::Compression::new(6))) } else { Sink::Plain(file) };
    Ok(std::io::BufWriter::with_capacity(1 << 16, sink))
}

fn finish(w: TextWriter, path: &str) -> Result<()> {
    let sink = w.into_inner().map_err(|e| anyhow::anyhow!("writing {path}: {}", e.error()))?;
    match sink {
        Sink::Plain(mut f) => f.flush().with_context(|| format!("writing {path}"))?,
        Sink::Gz(g) => {
            g.finish().with_context(|| format!("writing {path}"))?;
        }
    }
    Ok(())
}

pub fn baf_filename(dir: &str, sample: &str) -> String {
    format!("{dir}{sample}.amber.baf.tsv.gz")
}

/// AmberBAFFile.write
pub fn write_baf(path: &str, bafs: &[AmberBaf]) -> Result<()> {
    let mut w = create_writer(path)?;
    writeln!(w, "chromosome\tposition\ttumorBAF\ttumorModifiedBAF\ttumorDepth\tnormalBAF\tnormalModifiedBAF\tnormalDepth\tlohRecovered")?;
    for b in bafs {
        let normal = if b.normal_baf.is_finite() { decimal_format_4(b.normal_baf) } else { "0".into() };
        let normal_mod = if b.normal_modified_baf().is_finite() { decimal_format_4(b.normal_modified_baf()) } else { "0".into() };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            b.chromosome,
            b.position,
            decimal_format_4(b.tumor_baf),
            decimal_format_4(b.tumor_modified_baf()),
            b.tumor_depth,
            normal,
            normal_mod,
            b.normal_depth,
            b.loh_recovered
        )?;
    }
    finish(w, path)
}

/// AmberQC.status
fn qc_status(contamination: f64) -> &'static str {
    if doubles::greater_than(contamination, 0.1) {
        "FAIL"
    } else if doubles::greater_than(contamination, 0.0) {
        "WARN"
    } else {
        "PASS"
    }
}

/// AmberQCFile.write
pub fn write_qc(path: &str, contamination: f64, consanguinity_proportion: f64, uniparental_disomy: Option<&str>) -> Result<()> {
    let text = format!(
        "QCStatus\t{}\nContamination\t{}\nConsanguinityProportion\t{}\nUniparentalDisomy\t{}\n",
        qc_status(contamination),
        decimal_format_4(contamination),
        decimal_format_4(consanguinity_proportion),
        uniparental_disomy.unwrap_or("NONE")
    );
    std::fs::write(path, text).with_context(|| format!("writing {path}"))
}

/// TumorContaminationFile.write, the sites in the order given
pub fn write_contamination(path: &str, sites: &[TumorContamination]) -> Result<()> {
    let mut w = create_writer(path)?;
    writeln!(w, "chromosome\tposition\tref\talt\tnormalDepth\tnormalRefSupport\tnormalAltSupport\ttumorDepth\ttumorRefSupport\ttumorAltSupport")?;
    for s in sites {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.chr,
            s.pos,
            s.ref_.as_str(),
            s.alt.as_str(),
            s.normal.read_depth,
            s.normal.ref_support,
            s.normal.alt_support,
            s.tumor.read_depth,
            s.tumor.ref_support,
            s.tumor.alt_support
        )?;
    }
    finish(w, path)
}

/// RegionOfHomozygosityFile.softFilter
fn roh_filter(r: &RegionOfHomozygosity) -> String {
    let mut filters = Vec::new();
    if r.length() < HOMOZYGOUS_REGION_LONG_SIZE {
        filters.push("minLength");
    }
    let het_ratio = r.num_heterozygous as f64 / r.snp_count() as f64;
    if het_ratio > HOMOZYGOUS_REGION_MAX_HET_RATIO {
        filters.push("maxHetProportion");
    }
    if filters.is_empty() { "PASS".into() } else { filters.join(";") }
}

/// RegionOfHomozygosityFile.write
pub fn write_homozygous_regions(path: &str, regions: &[RegionOfHomozygosity]) -> Result<()> {
    let mut w = create_writer(path)?;
    writeln!(w, "chromosome\tstartPosition\tendPosition\tsnpCount\thomCount\thetCount\tfilter")?;
    for r in regions {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.chromosome_name(),
            r.start,
            r.end,
            r.snp_count(),
            r.num_homozygous,
            r.num_heterozygous,
            roh_filter(r)
        )?;
    }
    finish(w, path)
}

/// VersionInfo.write: the amber.version resource of the build
pub fn write_version(dir: &str) -> Result<()> {
    let text = format!("version={}\nbuild.date={}", crate::VERSION, env!("AMBER_BUILD_DATE"));
    let path = format!("{dir}/amber.version");
    std::fs::write(&path, text).with_context(|| format!("writing {path}"))
}

pub fn tumor_data_filename(dir: &str, sample: &str) -> String {
    format!("{dir}{sample}.amber.tumor.raw.tsv.gz")
}

/// PositionEvidenceFile.write: the LOH columns when the first row has LOH data
pub fn write_tumor_data(path: &str, rows: &[&PositionEvidence]) -> Result<()> {
    let mut w = create_writer(path)?;
    write!(w, "Chromosome\tPosition\tRef\tAlt\tReadDepth\tIndelCount\tRefCount\tAltCount\tBaseQualFiltered\tMapQualFiltered\tSeqTechFiltered")?;
    if rows.first().is_some_and(|r| r.loh_data.is_some()) {
        write!(w, "\tValidBafHigh\tValidBafMed\tValidBafLow\tValidRateHigh\tValidRateMed\tValidRateLow\tLohCandidate\tLohRecovered")?;
    }
    writeln!(w)?;
    for p in rows {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p.chromosome,
            p.position,
            p.ref_.as_str(),
            p.alt.as_str(),
            p.read_depth,
            p.indel_count,
            p.ref_support,
            p.alt_support,
            p.base_qual_filtered,
            p.map_qual_filtered,
            p.seq_tech_filtered
        )?;
        if let Some(l) = &p.loh_data {
            use hmf_common::fmt::java_fixed;
            write!(
                w,
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                l.valid_baf_high,
                l.valid_baf_medium,
                l.valid_baf_low,
                java_fixed(l.valid_rate_high, 2),
                java_fixed(l.valid_rate_medium, 2),
                java_fixed(l.valid_rate_low, 2),
                l.loh_candidate,
                l.loh_recovered
            )?;
        }
        writeln!(w)?;
    }
    finish(w, path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimal_format() {
        assert_eq!(decimal_format_4(0.5), "0.5000");
        assert_eq!(decimal_format_4(-1.0), "-1.0000");
        assert_eq!(decimal_format_4(2.0), "2.0000");
        // values from Java's DecimalFormat("0.0000")
        assert_eq!(decimal_format_4(0.12345), "0.1235"); // binary value above the tie
        assert_eq!(decimal_format_4(0.5 + 0.00005), "0.5000"); // 0.50004999999999999449
        assert_eq!(decimal_format_4(1.0 / 160.0), "0.0063");
        assert_eq!(decimal_format_4(3.0 / 160.0), "0.0187");
        assert_eq!(decimal_format_4(0.03125), "0.0312"); // exact tie: to even
        assert_eq!(decimal_format_4(0.96875), "0.9688");
        assert_eq!(decimal_format_4(0.99995), "1.0000");
        assert_eq!(decimal_format_4(-0.0), "-0.0000");
        assert_eq!(decimal_format_4(0.00005), "0.0000"); // Java's FloatingDecimal rounding, not the exact value
        assert_eq!(decimal_format_4(0.00015), "0.0001");
        assert_eq!(decimal_format_4(1.23455), "1.2346");
        assert_eq!(decimal_format_4(0.0000499999), "0.0000");
        assert_eq!(decimal_format_4(f64::NAN), "NaN");
        assert_eq!(decimal_format_4(f64::INFINITY), "\u{221e}");
    }

    #[test]
    fn qc_status_levels() {
        assert_eq!(qc_status(0.0), "PASS");
        assert_eq!(qc_status(0.02), "WARN");
        assert_eq!(qc_status(0.2), "FAIL");
    }
}
