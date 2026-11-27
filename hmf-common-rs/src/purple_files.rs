//! hmf-common `purple` and `driver` files as LINX reads (and writes) them: PurpleCopyNumberFile
//! (`.purple.cnv.somatic.tsv`), PurityContextFile (`.purple.purity.tsv` with `.purple.qc`; the fields LINX
//! uses), GeneCopyNumberFile (`.purple.cnv.gene.tsv`) and DriverCatalogFile (`.driver.catalog.*.tsv`).
//! Enum-valued fields are kept as their names.

use anyhow::{anyhow, bail, Context, Result};

use crate::fmt::decimal_format_fixed;
use crate::java_text::{fields_index_map, parse_bool, parse_double, parse_int, read_lines};

fn check_add_dir_separator(path: &str) -> String {
    if path.is_empty() || path.ends_with('/') { path.to_string() } else { format!("{path}/") }
}

fn col(f: &std::collections::HashMap<String, usize>, name: &str, file: &str) -> Result<usize> {
    f.get(name).copied().ok_or_else(|| anyhow!("{file}: missing column {name}"))
}

fn int(v: &str) -> Result<i32> {
    parse_int(v).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{v}\""))
}
fn dbl(v: &str) -> Result<f64> {
    parse_double(v).ok_or_else(|| anyhow!("NumberFormatException: For input string: \"{v}\""))
}

/// PurpleCopyNumber
#[derive(Clone, Debug)]
pub struct PurpleCopyNumber {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub average_tumor_copy_number: f64,
    pub baf_count: i32,
    pub average_observed_baf: f64,
    pub average_actual_baf: f64,
    pub segment_start_support: String,
    pub segment_end_support: String,
    pub method: String,
    pub depth_window_count: i32,
    pub gc_content: f64,
    pub min_start: i32,
    pub max_start: i32,
}

impl PurpleCopyNumber {
    pub fn minor_allele_copy_number(&self) -> f64 {
        (1.0 - self.average_actual_baf) * self.average_tumor_copy_number
    }
    pub fn major_allele_copy_number(&self) -> f64 {
        self.average_tumor_copy_number - self.minor_allele_copy_number()
    }
}

/// PurpleCopyNumberFile.generateFilenameForReading
pub fn purple_copy_number_filename(base_path: &str, sample: &str) -> String {
    let f = format!("{base_path}/{sample}.purple.cnv.somatic.tsv");
    if std::path::Path::new(&f).exists() { f } else { format!("{base_path}/{sample}.purple.cnv") }
}

/// PurpleCopyNumberFile.read
pub fn read_purple_copy_numbers(path: &str) -> Result<Vec<PurpleCopyNumber>> {
    let lines = read_lines(path).with_context(|| format!("failed to read {path}"))?;
    let Some(header) = lines.first() else { bail!("{path}: empty file") };
    let f = fields_index_map(header, '\t');
    let c = |n: &str| col(&f, n, path);
    let (chr, start, end, cn, baf_count, obs_baf, baf) =
        (c("chromosome")?, c("start")?, c("end")?, c("copyNumber")?, c("bafCount")?, c("observedBAF")?, c("baf")?);
    let (seg_start, seg_end, method, dwc, gc, min_start, max_start) = (
        c("segmentStartSupport")?,
        c("segmentEndSupport")?,
        c("method")?,
        c("depthWindowCount")?,
        c("gcContent")?,
        c("minStart")?,
        c("maxStart")?,
    );
    let mut out = Vec::with_capacity(lines.len());
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = line.split('\t').collect();
        let g = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("{path}: short line: {line}"));
        out.push(PurpleCopyNumber {
            chromosome: g(chr)?.to_string(),
            start: int(g(start)?)?,
            end: int(g(end)?)?,
            average_tumor_copy_number: dbl(g(cn)?)?,
            baf_count: int(g(baf_count)?)?,
            average_observed_baf: dbl(g(obs_baf)?)?,
            average_actual_baf: dbl(g(baf)?)?,
            segment_start_support: g(seg_start)?.to_string(),
            segment_end_support: g(seg_end)?.to_string(),
            method: g(method)?.to_string(),
            depth_window_count: int(g(dwc)?)?,
            gc_content: dbl(g(gc)?)?,
            min_start: int(g(min_start)?)?,
            max_start: int(g(max_start)?)?,
        });
    }
    Ok(out)
}

/// The parts of PurityContext that LINX uses
#[derive(Clone, Debug)]
pub struct PurityContext {
    pub purity: f64,
    pub ploidy: f64,
    /// Gender's name: MALE, FEMALE, MALE_KLINEFELTER
    pub gender: String,
    pub whole_genome_duplication: bool,
}

/// PurityContextFile.read(basePath, sample): the QC file must exist and parse as well
pub fn read_purity_context(base_path: &str, sample: &str) -> Result<PurityContext> {
    let base = check_add_dir_separator(base_path);
    let qc_file = format!("{base}{sample}.purple.qc");
    let purity_file = format!("{base}{sample}.purple.purity.tsv");
    read_lines(&qc_file).with_context(|| format!("java.nio.file.NoSuchFileException: {qc_file}"))?;
    let lines = read_lines(&purity_file).with_context(|| format!("java.nio.file.NoSuchFileException: {purity_file}"))?;
    if lines.len() < 2 {
        bail!("{purity_file}: no purity line");
    }
    let f = fields_index_map(&lines[0], '\t');
    let v: Vec<&str> = lines[1].split('\t').collect();
    let g = |n: &str| -> Result<&str> {
        let i = col(&f, n, &purity_file)?;
        v.get(i).copied().ok_or_else(|| anyhow!("{purity_file}: short line"))
    };
    Ok(PurityContext {
        purity: dbl(g("purity")?)?,
        ploidy: dbl(g("ploidy")?)?,
        gender: g("gender")?.to_string(),
        whole_genome_duplication: parse_bool(g("wholeGenomeDuplication")?),
    })
}

/// GeneCopyNumber
#[derive(Clone, Debug)]
pub struct GeneCopyNumber {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub gene_name: String,
    pub trans_name: String,
    pub is_canonical: bool,
    pub chromosome_band: String,
    pub max_copy_number: f64,
    pub min_copy_number: f64,
    pub min_minor_allele_copy_number: f64,
    pub somatic_regions: i32,
    pub min_regions: i32,
    pub min_region_start: i32,
    pub min_region_end: i32,
    pub depth_window_count: i32,
    pub gc_content: f64,
    pub min_region_start_support: String,
    pub min_region_end_support: String,
    pub min_region_method: String,
    pub relative_min_copy_number: f64,
    pub reported_status: String,
    pub driver_type: String,
}

impl GeneCopyNumber {
    pub fn min_region_bases(&self) -> i32 {
        self.min_region_end - self.min_region_start + 1
    }
}

/// GeneCopyNumberFile.generateFilename
pub fn gene_copy_number_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.purple.cnv.gene.tsv", check_add_dir_separator(base_path))
}

/// GeneCopyNumberFile.read
pub fn read_gene_copy_numbers(path: &str) -> Result<Vec<GeneCopyNumber>> {
    let lines = read_lines(path).with_context(|| format!("failed to read {path}"))?;
    let Some(header) = lines.first() else { bail!("{path}: empty file") };
    let f = fields_index_map(header, '\t');
    let c = |n: &str| col(&f, n, path);
    let idx = [
        c("chromosome")?,
        c("start")?,
        c("end")?,
        c("gene")?,
        c("minCopyNumber")?,
        c("maxCopyNumber")?,
        c("somaticRegions")?,
        c("transcriptId")?,
        c("isCanonical")?,
        c("chromosomeBand")?,
        c("minRegions")?,
        c("minRegionStart")?,
        c("minRegionEnd")?,
        c("minRegionStartSupport")?,
        c("minRegionEndSupport")?,
        c("minRegionMethod")?,
        c("minMinorAlleleCopyNumber")?,
        c("depthWindowCount")?,
    ];
    let gc_i = f.get("gcContent").copied();
    let reported_i = f.get("reportedStatus").copied();
    let driver_i = f.get("driverType").copied();
    let rel_i = f.get("relativeMinCopyNumber").copied();
    let mut out = Vec::with_capacity(lines.len());
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = line.split('\t').collect();
        let g = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("{path}: short line: {line}"));
        out.push(GeneCopyNumber {
            chromosome: g(idx[0])?.to_string(),
            start: int(g(idx[1])?)?,
            end: int(g(idx[2])?)?,
            gene_name: g(idx[3])?.to_string(),
            min_copy_number: dbl(g(idx[4])?)?,
            max_copy_number: dbl(g(idx[5])?)?,
            somatic_regions: int(g(idx[6])?)?,
            trans_name: g(idx[7])?.to_string(),
            is_canonical: parse_bool(g(idx[8])?),
            chromosome_band: g(idx[9])?.to_string(),
            min_regions: int(g(idx[10])?)?,
            min_region_start: int(g(idx[11])?)?,
            min_region_end: int(g(idx[12])?)?,
            min_region_start_support: g(idx[13])?.to_string(),
            min_region_end_support: g(idx[14])?.to_string(),
            min_region_method: g(idx[15])?.to_string(),
            min_minor_allele_copy_number: dbl(g(idx[16])?)?,
            depth_window_count: int(g(idx[17])?)?,
            gc_content: match gc_i {
                Some(i) => dbl(g(i)?)?,
                None => 0.0,
            },
            relative_min_copy_number: match rel_i {
                Some(i) => dbl(g(i)?)?,
                None => -1.0,
            },
            reported_status: match reported_i {
                Some(i) => g(i)?.to_string(),
                None => "NONE".to_string(),
            },
            driver_type: match driver_i {
                Some(i) => g(i)?.to_string(),
                None => "UNKNOWN".to_string(),
            },
        });
    }
    Ok(out)
}

/// DriverCatalog
#[derive(Clone, Debug, PartialEq)]
pub struct DriverCatalog {
    pub chromosome: String,
    pub chromosome_band: String,
    pub gene: String,
    pub transcript: String,
    pub is_canonical: bool,
    /// DriverType name
    pub driver: String,
    /// DriverCategory name
    pub category: String,
    pub likelihood_method: String,
    pub reported_status: String,
    pub driver_likelihood: f64,
    pub missense: i32,
    pub nonsense: i32,
    pub splice: i32,
    pub inframe: i32,
    pub frameshift: i32,
    pub biallelic: bool,
    pub min_copy_number: f64,
    pub max_copy_number: f64,
}

const DRIVER_COLUMNS: [&str; 18] = [
    "chromosome",
    "chromosomeBand",
    "gene",
    "transcript",
    "isCanonical",
    "driver",
    "category",
    "likelihoodMethod",
    "reportedStatus",
    "driverLikelihood",
    "missense",
    "nonsense",
    "splice",
    "inframe",
    "frameshift",
    "biallelic",
    "minCopyNumber",
    "maxCopyNumber",
];

/// The DriverType names (checkConvertType maps any other to UNKNOWN)
pub const DRIVER_TYPES: [&str; 14] = [
    "AMP",
    "PARTIAL_AMP",
    "DEL",
    "HET_DEL",
    "LOH",
    "HOM_DUP_DISRUPTION",
    "HOM_DEL_DISRUPTION",
    "DISRUPTION",
    "MUTATION",
    "GERMLINE_MUTATION",
    "GERMLINE_DELETION",
    "GERMLINE_AMP",
    "GERMLINE_DISRUPTION",
    "UNKNOWN",
];

/// DriverCatalogFile.generateSomaticFilename
pub fn driver_catalog_somatic_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.purple.driver.catalog.somatic.tsv")
}

/// DriverCatalogFile.read
pub fn read_driver_catalog(path: &str) -> Result<Vec<DriverCatalog>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let Some(header) = lines.first() else { bail!("{path}: empty file") };
    let f = fields_index_map(header, '\t');
    let reported_i = f.get("reportedStatus").copied();
    let mut out = Vec::new();
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = line.split('\t').collect();
        let g = |n: &str| -> Result<&str> {
            let i = col(&f, n, path)?;
            v.get(i).copied().ok_or_else(|| anyhow!("{path}: short line: {line}"))
        };
        let driver = g("driver")?;
        out.push(DriverCatalog {
            chromosome: g("chromosome")?.to_string(),
            chromosome_band: g("chromosomeBand")?.to_string(),
            gene: g("gene")?.to_string(),
            transcript: g("transcript")?.to_string(),
            is_canonical: parse_bool(g("isCanonical")?),
            driver: if DRIVER_TYPES.contains(&driver) { driver.to_string() } else { "UNKNOWN".to_string() },
            category: g("category")?.to_string(),
            likelihood_method: g("likelihoodMethod")?.to_string(),
            reported_status: match reported_i {
                Some(i) => v.get(i).copied().unwrap_or("").to_string(),
                None => "REPORTED".to_string(),
            },
            driver_likelihood: dbl(g("driverLikelihood")?)?,
            missense: int(g("missense")?)?,
            nonsense: int(g("nonsense")?)?,
            splice: int(g("splice")?)?,
            inframe: int(g("inframe")?)?,
            frameshift: int(g("frameshift")?)?,
            biallelic: parse_bool(g("biallelic")?),
            min_copy_number: dbl(g("minCopyNumber")?)?,
            max_copy_number: dbl(g("maxCopyNumber")?)?,
        });
    }
    Ok(out)
}

/// DriverCatalogFile.toLines
pub fn driver_catalog_lines(catalog: &[DriverCatalog]) -> Vec<String> {
    let mut lines = vec![DRIVER_COLUMNS.join("\t")];
    for d in catalog {
        lines.push(
            [
                d.chromosome.clone(),
                d.chromosome_band.clone(),
                d.gene.clone(),
                d.transcript.clone(),
                d.is_canonical.to_string(),
                d.driver.clone(),
                d.category.clone(),
                d.likelihood_method.clone(),
                d.reported_status.clone(),
                decimal_format_fixed(d.driver_likelihood, 4),
                d.missense.to_string(),
                d.nonsense.to_string(),
                d.splice.to_string(),
                d.inframe.to_string(),
                d.frameshift.to_string(),
                d.biallelic.to_string(),
                decimal_format_fixed(d.min_copy_number, 4),
                decimal_format_fixed(d.max_copy_number, 4),
            ]
            .join("\t"),
        );
    }
    lines
}

/// Files.write(path, lines): each line followed by '\n'
pub fn write_lines(path: &str, lines: &[String]) -> Result<()> {
    let mut s = String::with_capacity(lines.iter().map(|l| l.len() + 1).sum());
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    std::fs::write(path, s).with_context(|| format!("failed to write {path}"))
}

/// DriverCatalogFile.write
pub fn write_driver_catalog(path: &str, catalog: &[DriverCatalog]) -> Result<()> {
    write_lines(path, &driver_catalog_lines(catalog))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn driver_catalog_round_trip() {
        let d = DriverCatalog {
            chromosome: "chr1".into(),
            chromosome_band: "p36".into(),
            gene: "GENE1".into(),
            transcript: "ENST1".into(),
            is_canonical: true,
            driver: "DEL".into(),
            category: "TSG".into(),
            likelihood_method: "DEL".into(),
            reported_status: "REPORTED".into(),
            driver_likelihood: 1.0,
            missense: 0,
            nonsense: 0,
            splice: 0,
            inframe: 0,
            frameshift: 0,
            biallelic: false,
            min_copy_number: 0.12344,
            max_copy_number: 2.5,
        };
        let lines = driver_catalog_lines(std::slice::from_ref(&d));
        assert_eq!(lines[1], "chr1\tp36\tGENE1\tENST1\ttrue\tDEL\tTSG\tDEL\tREPORTED\t1.0000\t0\t0\t0\t0\t0\tfalse\t0.1234\t2.5000");
        let dir = std::env::temp_dir().join(format!("hmf_common_driver_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("d.tsv");
        write_driver_catalog(p.to_str().unwrap(), &[d.clone()]).unwrap();
        let back = read_driver_catalog(p.to_str().unwrap()).unwrap();
        assert_eq!(back[0].gene, "GENE1");
        assert_eq!(back[0].min_copy_number, 0.1234);
        std::fs::remove_dir_all(&dir).ok();
    }
}

// ---------------------------------------------------------------------------------------------------------
// PurpleQCFile (added 2026-11-27 for Virus Interpreter)

/// PurpleQCStatus, in enum order
pub const PURPLE_QC_STATUSES: [&str; 10] = [
    "PASS",
    "WARN_DELETED_GENES",
    "WARN_HIGH_COPY_NUMBER_NOISE",
    "WARN_GENDER_MISMATCH",
    "WARN_LOW_PURITY",
    "WARN_TINC",
    "WARN_CONTAMINATION",
    "FAIL_CONTAMINATION",
    "FAIL_NO_TUMOR",
    "FAIL_TINC",
];

const FITTED_PURITY_METHODS: [&str; 4] = ["NORMAL", "HIGHLY_DIPLOID", "SOMATIC", "NO_TUMOR"];
const GENDERS: [&str; 2] = ["MALE", "FEMALE"];
const GERMLINE_ABERRATIONS: [&str; 10] =
    ["NONE", "MOSAIC_X", "KLINEFELTER", "XYY", "TRISOMY_X", "TRISOMY_13", "TRISOMY_15", "TRISOMY_18", "TRISOMY_21", "TERTRASOMY_9"];

/// PurpleQC as PurpleQCFile.read gives it; enum values kept as their names
#[derive(Clone, Debug, PartialEq)]
pub struct PurpleQc {
    /// the QC statuses (a set in Java), each a PurpleQCStatus name, in file order without repeats
    pub status: Vec<String>,
    pub method: String,
    pub copy_number_segments: i32,
    pub unsupported_copy_number_segments: i32,
    pub purity: f64,
    pub amber_gender: String,
    pub cobalt_gender: String,
    pub deleted_genes: i32,
    pub contamination: f64,
    pub germline_aberrations: Vec<String>,
    pub amber_mean_depth: i32,
    pub loh_percent: f64,
    pub tinc_level: f64,
}

impl PurpleQc {
    pub fn has_status(&self, status: &str) -> bool {
        self.status.iter().any(|s| s == status)
    }
}

/// PurpleQCFile.generateFilename: basePath + File.separator + sample + ".purple.qc"
pub fn purple_qc_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.purple.qc")
}

fn enum_value(value: &str, names: &[&str], enum_name: &str) -> Result<String> {
    if names.contains(&value) {
        Ok(value.to_string())
    } else {
        bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.{enum_name}.{value}")
    }
}

/// PurpleQCFile.read: each field from the first line whose first tab-separated field is its name (Java's
/// FileReaderUtils.getValue), else its default; errors where Java throws (an unknown status, method, gender
/// or aberration name, a missing gender line, a number that does not parse)
pub fn read_purple_qc(path: &str) -> Result<PurpleQc> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let get = |field: &str, default: Option<&str>| -> Result<Option<String>> {
        for line in &lines {
            let values: Vec<&str> = line.split('\t').collect();
            if values[0] == field {
                return values
                    .get(1)
                    .map(|v| Some(v.to_string()))
                    .ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: {path}: no value for {field}"));
            }
        }
        Ok(default.map(String::from))
    };
    let num_int = |field: &str, default: &str| -> Result<i32> { int(&get(field, Some(default))?.unwrap()) };
    let num_dbl = |field: &str, default: &str| -> Result<f64> { dbl(&get(field, Some(default))?.unwrap()) };
    let gender = |field: &str| -> Result<String> {
        match get(field, None)? {
            None => bail!("java.lang.NullPointerException: Name is null ({path}: no {field})"),
            Some(v) => enum_value(&v, &GENDERS, "purple.Gender"),
        }
    };

    // PurpleQCStatus.fromString, with the old FAIL_ names read as their WARN_ statuses
    let mut status: Vec<String> = Vec::new();
    for s in get("QCStatus", Some(""))?.unwrap().split(',') {
        let s = match s {
            "FAIL_SEGMENT" => "WARN_HIGH_COPY_NUMBER_NOISE",
            "FAIL_GENDER" => "WARN_GENDER_MISMATCH",
            "FAIL_DELETED_GENES" => "WARN_DELETED_GENES",
            other => other,
        };
        let s = enum_value(s, &PURPLE_QC_STATUSES, "purple.PurpleQCStatus")?;
        if !status.contains(&s) {
            status.push(s);
        }
    }
    let method = enum_value(&get("Method", Some("NORMAL"))?.unwrap(), &FITTED_PURITY_METHODS, "purple.FittedPurityMethod")?;
    let copy_number_segments = num_int("CopyNumberSegments", "0")?;
    let unsupported_copy_number_segments = num_int("UnsupportedCopyNumberSegments", "0")?;
    let purity = num_dbl("Purity", "0")?;
    let amber_gender = gender("AmberGender")?;
    let cobalt_gender = gender("CobaltGender")?;
    let deleted_genes = num_int("DeletedGenes", "0")?;
    let contamination = num_dbl("Contamination", "0")?;
    let mut germline_aberrations = Vec::new();
    for a in get("GermlineAberrations", Some("NONE"))?.unwrap().split(',') {
        let a = enum_value(a, &GERMLINE_ABERRATIONS, "genome.chromosome.GermlineAberration")?;
        if !germline_aberrations.contains(&a) {
            germline_aberrations.push(a);
        }
    }
    Ok(PurpleQc {
        status,
        method,
        copy_number_segments,
        unsupported_copy_number_segments,
        purity,
        amber_gender,
        cobalt_gender,
        deleted_genes,
        contamination,
        germline_aberrations,
        amber_mean_depth: num_int("AmberMeanDepth", "0")?,
        loh_percent: num_dbl("LohPercent", "-1")?,
        tinc_level: num_dbl("TincLevel", "0")?,
    })
}

#[cfg(test)]
mod purple_qc_tests {
    use super::*;

    #[test]
    fn reads_the_qc_file() {
        let dir = std::env::temp_dir().join(format!("purple_qc_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("S.purple.qc");
        std::fs::write(
            &path,
            "QCStatus\tWARN_LOW_PURITY,FAIL_SEGMENT\nMethod\tNO_TUMOR\nCopyNumberSegments\t12\nPurity\t0.1500\nAmberGender\tFEMALE\nCobaltGender\tMALE\nGermlineAberrations\tNONE\n",
        )
        .unwrap();
        let qc = read_purple_qc(path.to_str().unwrap()).unwrap();
        assert_eq!(qc.status, vec!["WARN_LOW_PURITY", "WARN_HIGH_COPY_NUMBER_NOISE"]);
        assert_eq!(qc.method, "NO_TUMOR");
        assert_eq!(qc.copy_number_segments, 12);
        assert_eq!(qc.loh_percent, -1.0);
        assert!(qc.has_status("WARN_LOW_PURITY") && !qc.has_status("PASS"));
        std::fs::write(&path, "QCStatus\tPASS\nAmberGender\tFEMALE\n").unwrap();
        assert!(read_purple_qc(path.to_str().unwrap()).is_err());
        std::fs::write(&path, "QCStatus\tNOT_A_STATUS\nAmberGender\tFEMALE\nCobaltGender\tFEMALE\n").unwrap();
        assert!(read_purple_qc(path.to_str().unwrap()).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

// ---------------------------------------------------------------------------------------------------------
// PurplePurity (added 2026-11-27 for qsee): every column of `.purple.purity.tsv`

const RUN_MODES: [&str; 3] = ["TUMOR_GERMLINE", "TUMOR", "GERMLINE"];
const MICROSATELLITE_STATUSES: [&str; 3] = ["MSI", "MSS", "UNKNOWN"];
const TUMOR_MUTATIONAL_STATUSES: [&str; 3] = ["HIGH", "LOW", "UNKNOWN"];

/// PurplePurity, as PurplePurity.read gives it; enum values kept as their names
#[derive(Clone, Debug, PartialEq)]
pub struct PurplePurity {
    pub purity: f64,
    pub norm_factor: f64,
    pub score: f64,
    pub diploid_proportion: f64,
    pub ploidy: f64,
    /// Gender: MALE or FEMALE
    pub sex: String,
    /// FittedPurityMethod (the `status` column)
    pub fit_method: String,
    pub polyclonal_proportion: f64,
    pub min_purity: f64,
    pub max_purity: f64,
    pub min_ploidy: f64,
    pub max_ploidy: f64,
    pub min_diploid_proportion: f64,
    pub max_diploid_proportion: f64,
    pub somatic_penalty: f64,
    pub whole_genome_duplication: bool,
    pub ms_indels_per_mb: f64,
    pub ms_status: String,
    pub tml: i32,
    pub tml_status: String,
    pub tmb_per_mb: f64,
    pub tmb_status: String,
    pub sv_tumor_mutational_burden: i32,
    pub run_mode: String,
    pub targeted: bool,
}

/// PurplePurity.generateFilename: basePath + File.separator + sample + ".purple.purity.tsv" (no check for
/// a separator already there, so an empty base path gives a path at the file system root, as in Java)
pub fn purple_purity_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.purple.purity.tsv")
}

/// PurplePurity.read: the header and the first values line; a missing column, a short line, a number
/// that does not parse or an unknown enum name is an error, as Java's NullPointerException,
/// ArrayIndexOutOfBoundsException, NumberFormatException and IllegalArgumentException
pub fn read_purple_purity(path: &str) -> Result<PurplePurity> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    purple_purity_from_lines(&lines).with_context(|| format!("reading {path}"))
}

/// PurplePurity.read on lines already read
pub fn purple_purity_from_lines(lines: &[String]) -> Result<PurplePurity> {
    if lines.len() < 2 {
        bail!("java.lang.IndexOutOfBoundsException: no values line");
    }
    let f = fields_index_map(&lines[0], '\t');
    let v: Vec<&str> = lines[1].split('\t').collect();
    let g = |n: &str| -> Result<&str> {
        let i = f.get(n).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {n}"))?;
        v.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: no value for {n}"))
    };
    let d = |n: &str| -> Result<f64> { dbl(g(n)?) };
    let e = |n: &str, names: &[&str], enum_name: &str| -> Result<String> { enum_value(g(n)?, names, enum_name) };
    Ok(PurplePurity {
        purity: d("purity")?,
        norm_factor: d("normFactor")?,
        score: d("score")?,
        diploid_proportion: d("diploidProportion")?,
        ploidy: d("ploidy")?,
        sex: e("gender", &GENDERS, "purple.Gender")?,
        fit_method: e("status", &FITTED_PURITY_METHODS, "purple.FittedPurityMethod")?,
        polyclonal_proportion: d("polyclonalProportion")?,
        min_purity: d("minPurity")?,
        max_purity: d("maxPurity")?,
        min_ploidy: d("minPloidy")?,
        max_ploidy: d("maxPloidy")?,
        min_diploid_proportion: d("minDiploidProportion")?,
        max_diploid_proportion: d("maxDiploidProportion")?,
        somatic_penalty: d("somaticPenalty")?,
        whole_genome_duplication: parse_bool(g("wholeGenomeDuplication")?),
        ms_indels_per_mb: d("msIndelsPerMb")?,
        ms_status: e("msStatus", &MICROSATELLITE_STATUSES, "purple.MicrosatelliteStatus")?,
        tml: int(g("tml")?)?,
        tml_status: e("tmlStatus", &TUMOR_MUTATIONAL_STATUSES, "purple.TumorMutationalStatus")?,
        tmb_per_mb: d("tmbPerMb")?,
        tmb_status: e("tmbStatus", &TUMOR_MUTATIONAL_STATUSES, "purple.TumorMutationalStatus")?,
        sv_tumor_mutational_burden: int(g("svTumorMutationalBurden")?)?,
        run_mode: e("runMode", &RUN_MODES, "purple.RunMode")?,
        targeted: parse_bool(g("targeted")?),
    })
}

#[cfg(test)]
mod purple_purity_tests {
    use super::*;

    #[test]
    fn reads_every_column() {
        let header = "purity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\truntimeMode";
        let header = header.replace("runtimeMode", "runMode\ttargeted");
        let values = "0.5000\t1.0\t0.6\t0.3\t3.1000\tMALE\tSOMATIC\t0.01\t0.45\t0.55\t3.0\t3.2\t0.2\t0.4\t0.1\ttrue\t0.25\tMSS\t120\tLOW\t5.5\tLOW\t77\tTUMOR\tfalse";
        let p = purple_purity_from_lines(&[header.clone(), values.to_string()]).unwrap();
        assert_eq!((p.purity, p.ploidy, p.sex.as_str(), p.fit_method.as_str()), (0.5, 3.1, "MALE", "SOMATIC"));
        assert_eq!((p.tml, p.sv_tumor_mutational_burden, p.whole_genome_duplication, p.targeted), (120, 77, true, false));
        assert_eq!((p.ms_indels_per_mb, p.tmb_per_mb, p.run_mode.as_str()), (0.25, 5.5, "TUMOR"));
        assert!(purple_purity_from_lines(&[header.clone(), values.replace("MALE", "OTHER")]).is_err());
        assert!(purple_purity_from_lines(&[header.replace("tmbPerMb", "x"), values.to_string()]).is_err());
        assert!(purple_purity_from_lines(&[header]).is_err());
        assert_eq!(purple_purity_filename("", "S"), "/S.purple.purity.tsv");
    }
}

/// PurityContextFile.readWithQC: the QC file, then the purity file
pub fn read_purity_context_with_qc(qc_file: &str, purity_file: &str) -> Result<(PurpleQc, PurplePurity)> {
    let qc = read_purple_qc(qc_file)?;
    let purity = read_purple_purity(purity_file)?;
    Ok((qc, purity))
}

/// DriverCatalogFile.generateGermlineFilename
pub fn driver_catalog_germline_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.purple.driver.catalog.germline.tsv")
}

/// ChrArmCopyNumber (chromosome as HumanChromosome's name: "1".."22", "X", "Y"; arm "P" or "Q")
#[derive(Clone, Debug, PartialEq)]
pub struct ChrArmCopyNumber {
    pub chromosome: String,
    pub arm: String,
    pub mean_copy_number: f64,
    pub median_copy_number: f64,
    pub min_copy_number: f64,
    pub max_copy_number: f64,
}

/// HumanChromosome.fromString: the enum name without a "chr" prefix (any case); Err where valueOf throws
pub fn human_chromosome(chromosome: &str) -> Result<String> {
    let name = if chromosome.len() >= 3 && chromosome[..3].eq_ignore_ascii_case("chr") { &chromosome[3..] } else { chromosome };
    if crate::ref_coords::chromosome_length(name, false).is_some() {
        Ok(name.to_string())
    } else {
        bail!("java.lang.IllegalArgumentException: No enum constant HumanChromosome._{name}")
    }
}

/// ChrArmCopyNumbersFile.generateFilename
pub fn chr_arm_copy_number_filename(base_path: &str, sample: &str) -> String {
    format!("{}{sample}.purple.chromosome_arm.tsv", check_add_dir_separator(base_path))
}

/// ChrArmCopyNumbersFile.read
pub fn read_chr_arm_copy_numbers(path: &str) -> Result<Vec<ChrArmCopyNumber>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let Some(header) = lines.first() else { bail!("{path}: IndexOutOfBoundsException: Index 0 out of bounds for length 0") };
    let f = fields_index_map(header, '\t');
    let c = |n: &str| col(&f, n, path);
    let (chr, arm, mean, median, min, max) =
        (c("chromosome")?, c("arm")?, c("meanCopyNumber")?, c("medianCopyNumber")?, c("minCopyNumber")?, c("maxCopyNumber")?);
    let mut out = Vec::new();
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = line.split('\t').collect();
        let g = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("{path}: short line: {line}"));
        let arm_value = g(arm)?;
        if arm_value != "P" && arm_value != "Q" {
            bail!("java.lang.IllegalArgumentException: No enum constant segmentation.Arm.{arm_value}");
        }
        out.push(ChrArmCopyNumber {
            chromosome: human_chromosome(g(chr)?)?,
            arm: arm_value.to_string(),
            mean_copy_number: dbl(g(mean)?)?,
            median_copy_number: dbl(g(median)?)?,
            min_copy_number: dbl(g(min)?)?,
            max_copy_number: dbl(g(max)?)?,
        });
    }
    Ok(out)
}

const GERMLINE_DETECTION_METHODS: [&str; 2] = ["SEGMENT", "EXON_COVERAGE"];
/// purple.GermlineStatus names
pub const GERMLINE_STATUSES: [&str; 9] =
    ["HOM_DELETION", "HET_DELETION", "AMPLIFICATION", "NOISE", "DIPLOID", "LIKELY_DIPLOID", "EXCLUDED", "CENTROMETIC", "UNKNOWN"];
const REPORTED_STATUSES: [&str; 3] = ["NONE", "NOT_REPORTED", "REPORTED"];

/// GermlineAmpDel; enum values kept as their names
#[derive(Clone, Debug, PartialEq)]
pub struct GermlineAmpDel {
    pub gene_name: String,
    pub transcript: String,
    pub chromosome: String,
    pub chromosome_band: String,
    pub region_start: i32,
    pub region_end: i32,
    pub depth_window_count: i32,
    pub exon_start: i32,
    pub exon_end: i32,
    pub is_partial: bool,
    pub detection_method: String,
    pub normal_status: String,
    pub tumor_status: String,
    pub germline_copy_number: f64,
    pub tumor_copy_number: f64,
    pub filter: String,
    pub cohort_frequency: i32,
    pub reported: String,
}

/// GermlineAmpDel.generateFilename
pub fn germline_amp_del_filename(base_path: &str, sample: &str) -> String {
    format!("{base_path}/{sample}.purple.germline_amp_del.tsv")
}

/// GermlineAmpDel.read (with the pre-v3.0 `reported` column and the optional transcript and isPartial)
pub fn read_germline_amp_dels(path: &str) -> Result<Vec<GermlineAmpDel>> {
    let lines = read_lines(path).with_context(|| format!("java.nio.file.NoSuchFileException: {path}"))?;
    let Some(header) = lines.first() else { bail!("{path}: IndexOutOfBoundsException: Index 0 out of bounds for length 0") };
    let f = fields_index_map(header, '\t');
    let reported_status_i = f.get("reportedStatus").copied();
    let reported_i = f.get("reported").copied();
    let transcript_i = f.get("transcript").copied();
    let partial_i = f.get("isPartial").copied();
    let mut out = Vec::new();
    for line in lines.iter().skip(1) {
        let v: Vec<&str> = line.split('\t').collect();
        let at = |i: usize| v.get(i).copied().ok_or_else(|| anyhow!("{path}: ArrayIndexOutOfBoundsException: Index {i}"));
        let g = |n: &str| -> Result<&str> { at(col(&f, n, path)?) };
        let reported = match reported_status_i {
            Some(i) => enum_value(at(i)?, &REPORTED_STATUSES, "purple.ReportedStatus")?,
            None => {
                let i = reported_i.ok_or_else(|| anyhow!("{path}: NullPointerException: no reported column"))?;
                if parse_bool(at(i)?) { "REPORTED".to_string() } else { "NONE".to_string() }
            }
        };
        out.push(GermlineAmpDel {
            gene_name: g("gene")?.to_string(),
            transcript: match transcript_i {
                Some(_) => g("transcript")?.to_string(),
                None => String::new(),
            },
            chromosome: g("chromosome")?.to_string(),
            chromosome_band: g("chromosomeBand")?.to_string(),
            region_start: int(g("regionStart")?)?,
            region_end: int(g("regionEnd")?)?,
            depth_window_count: int(g("depthWindowCount")?)?,
            exon_start: int(g("exonStart")?)?,
            exon_end: int(g("exonEnd")?)?,
            is_partial: match partial_i {
                Some(i) => parse_bool(at(i)?),
                None => false,
            },
            detection_method: enum_value(g("detectionMethod")?, &GERMLINE_DETECTION_METHODS, "purple.GermlineDetectionMethod")?,
            normal_status: enum_value(g("germlineStatus")?, &GERMLINE_STATUSES, "purple.GermlineStatus")?,
            tumor_status: enum_value(g("tumorStatus")?, &GERMLINE_STATUSES, "purple.GermlineStatus")?,
            germline_copy_number: dbl(g("germlineCopyNumber")?)?,
            tumor_copy_number: dbl(g("tumorCopyNumber")?)?,
            filter: g("filter")?.to_string(),
            cohort_frequency: int(g("cohortFrequency")?)?,
            reported,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod orange_addition_tests {
    use super::*;

    #[test]
    fn chromosome_names() {
        assert_eq!(human_chromosome("chr1").unwrap(), "1");
        assert_eq!(human_chromosome("CHRX").unwrap(), "X");
        assert_eq!(human_chromosome("22").unwrap(), "22");
        assert!(human_chromosome("chrM").is_err());
        assert!(human_chromosome("x").is_err());
    }

    #[test]
    fn purity_and_arms() {
        let dir = std::env::temp_dir().join(format!("purple_orange_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("S.purple.purity.tsv");
        std::fs::write(&p, "purity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\trunMode\ttargeted\n\
0.5000\t1.0000\t0.2000\t0.9000\t2.1000\tFEMALE\tNORMAL\t0.0000\t0.4500\t0.5500\t2.0000\t2.2000\t0.8000\t1.0000\t0.0000\tfalse\t0.1000\tMSS\t50\tLOW\t2.5000\tLOW\t30\tTUMOR_GERMLINE\tfalse\n").unwrap();
        let pp = read_purple_purity(p.to_str().unwrap()).unwrap();
        assert_eq!(pp.ploidy, 2.1);
        assert_eq!(pp.tml, 50);
        assert_eq!(pp.run_mode, "TUMOR_GERMLINE");
        let a = dir.join("S.purple.chromosome_arm.tsv");
        std::fs::write(&a, "chromosome\tarm\tmeanCopyNumber\tmedianCopyNumber\tminCopyNumber\tmaxCopyNumber\nchr1\tP\t2.0\t2.0\t1.5\t2.5\nchrX\tQ\t1.0\t1.0\t1.0\t1.0\n").unwrap();
        let arms = read_chr_arm_copy_numbers(a.to_str().unwrap()).unwrap();
        assert_eq!(arms[1].chromosome, "X");
        assert_eq!(arms[0].arm, "P");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
